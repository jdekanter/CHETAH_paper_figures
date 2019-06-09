data.dir <- "DEFINE"
output.dir <- "DEFINE"
source.dir <- "DEFINE"

library(CHETAH)
library(ggplot2)

## ----------------------------------------------- Head-Neck load + prepare data
## Load
load(paste0(data.dir, "HeadNeck_seurat.Rdata"))
load(paste0(data.dir,"CHETAH_Tumor_reference.Rdata"))
## Remove ribosomal:
ribo <- read.table(paste0(data.dir, "ribosomal.txt"), header=F, sep = '\t')
Celldata <- lapply(names(Celldata), function(x) Celldata[[x]][!rownames(Celldata[[x]]) %in% ribo[,1], ])
## Extract
input_hn_ct <- as.vector(seurat@meta.data$celltype)
names(input_hn_ct) <- rownames(seurat@meta.data)
input_hn <- as.matrix(seurat@raw.data)
tsne <- seurat@dr$tsne@cell.embeddings
## Prepare Ref
Celldata_hn <- lapply(Celldata, function(x) x[ ,!(colnames(x) %in% colnames(input_hn))])
Celldata_hn <- Celldata_hn[unlist(lapply(Celldata_hn, ncol)) >= 10]
## Prepare Input
input_hn_ct[input_hn_ct == 'Fibroblast'] <- 'CAF' ## Give same names
input_hn <- input_hn[ ,names(input_hn_ct)]
tsne <- tsne[names(input_hn_ct), ]
rm( seurat)
gc()

## ----------------------------------------------- Melanoma load + prepare data
## Load
load(paste0(data.dir, "Melanoma_seurat.Rdata"))
## Extract 
input_mel_ct <- as.vector(seurat@meta.data$celltype)
names(input_mel_ct) <- rownames(seurat@meta.data)
input_mel <- as.matrix(seurat@raw.data)
## Prepare Ref
Celldata_mel <- lapply(Celldata, function(x) x[ ,!(colnames(x) %in% colnames(input_mel))])
Celldata_mel <- Celldata[unlist(lapply(Celldata_mel, ncol)) >= 10]
tsne <- seurat@dr$tsne@cell.embeddings

input_mel_ct <- input_mel_ct[!(input_mel_ct %in% c('Unknown'))] ## Delete unclassified
input_mel <- input_mel[ ,names(input_mel_ct)]
tsne <- tsne[names(input_mel_ct), ]
rm(Celldata, seurat)
gc()

## ----------------------------------------------- CHETAH
## Group data
refs <- list("hn" = Celldata_hn,"mel" =  Celldata_mel)
input <- list("hn" = input_hn, "mel" = input_mel)
paper <- list("hn" = input_hn_ct, "mel" = input_mel_ct)
excl <- list('hn' = c('Mycote', 'Mast', 'Tumor', 'Resting Fibroblast', 'Myofibroblast'),
             'mel' = 'Tumor')
cor_meths <- c("spearman", "pearson", "kendall", "cosine")
nmbr_genes <- c(10, 25, 50, 100, 200, 500, 1000, 2000, 5000)
fld_threshs <- c(1.2, 1.5, 1.8, 2, 2.5, 3, 4, 5)

## Run CHETAH with the different parameters
data <- lapply(c('hn', 'mel'), function(dataset) {
    message(dataset, "---------------")
    ## Run correlation measures
    cors <- lapply(cor_meths, function(cor_meth) {
        message(cor_meth)
        CHETAHclassifier(input = input[[dataset]], 
                         ref_cells = refs[[dataset]],
                         cor_method = cor_meth)
    })
    names(cors) <- cor_meths
    ## Run with different # of genes
    ngenes <- lapply(nmbr_genes, function(number) {
        message(number)
        CHETAHclassifier(input = input[[dataset]], 
                         ref_cells = Celldata_hn,
                         n_genes = number)
    })
    names(ngenes) <- nmbr_genes
    ## Run with different fld-changes
    flds <- lapply(fld_threshs, function(fld) {
        message(fld)
        CHETAHclassifier(input = input[[dataset]], 
                         ref_cells = refs[[dataset]],
                         fix_ngenes = FALSE,
                         fc_thresh = fld)
    })
    names(flds) <- fld_threshs
    list("cors" = cors, 
         "ngenes" = ngenes,
         'flds' = flds)
})
names(data) <- c('hn', 'mel')
save(data, paper, excl , file = paste0(data.dir, "Parameter_Sweep.Rdata"))

## -------------------------------------------------------------------- Calculate performance
cor_meths <- c("spearman", "pearson", "kendall", "cosine")
nmbr_genes <- c(10, 25, 50, 100, 200, 500, 1000, 2000, 5000)
fld_threshs <- c(1.2, 1.5, 1.8, 2, 2.5, 3, 4, 5)
CalcPerf <- function(ref, type, excl_ct, info) {
    # % of no.ref unclassified
    norefcells <- type[ref %in% excl_ct]
    noref <- sum(grepl('Node|Unassigned', norefcells))
    noref_cl <- sum(!grepl('Node|Unassigned', norefcells))
    #% of yes.ref correct
    yesref <- type[!(ref %in% excl_ct)]
    yrref <- ref[names(yesref)]
    corr <- sum(yesref == yrref)
    #% of yes.ref incorrect{
    incorr <- sum(yesref != yrref & !grepl("Node|Unassigned", yesref))
    ## % of yes.ref in Node1 (new Node0 == Unassigned)
    yesrefi <- yesref[grepl('Node|Unassigned', yesref)]
    yrref <- ref[names(yesrefi)]
    incunass <- sum(grepl('Unassigned', yesrefi))
    ## % of yes.ref in correct branch
    perc_gb <- 0 ## percentage in good branch,
    for (i in unique(yrref)) {
        Nodes <- paste0('Node', which(unlist(lapply(info$nodetypes, function (x) i %in% names(x)))) - 1)
        Nodes <- Nodes[!grepl("Node0", Nodes)]
        perc_gb <- (perc_gb + sum(yesrefi[yrref == i] %in% Nodes))
    }
    perc_gb <- perc_gb
    ## % of yes.ref in incorrect branch
    perc_wb <- length(yesref) - incunass - perc_gb - corr - incorr
    out <- c(perc_wb + perc_gb, incunass, corr, incorr, noref, noref_cl)/length(ref)
    names(out) <- c('intermediate', 'unassigned', 'correct', 
                    'incorrect', 'correctly unclassified', 'erroneously classified')
    return(out)
}

perf <- do.call(rbind, lapply(c('mel', 'hn'), function(dataset) {
    out2 <- lapply(names(data[[dataset]]), function(dt) {
        out1 <- lapply(names(data[[dataset]][[dt]]), function(type) {
            chetah <- data[[dataset]][[dt]][[type]]
            perf <- CalcPerf(ref = paper[[dataset]], 
                             type = chetah$classification, 
                             excl_ct = excl[[dataset]], 
                             info = chetah)
            data.frame(perf_type = names(perf),
                       perf = unname(perf),
                       dataset = dataset,
                       comparison = dt,
                       selected = type,
                       stringsAsFactors = FALSE)
        })
        do.call(rbind, out1)
    })
    do.call(rbind, out2)
}))

clrs <- c(`correct`=              '#267c21',
         `correctly unclassified`='#40e038',
         `intermediate`          ='#BFEDBD',
         `unassigned`            ='#C1C1C1',
         `erroneously classified`='#d60202',
         `incorrect`             ='#a30914',
         `incorrect branch` =     '#f2b5b5')

## Set the order
perf$selected <- factor(perf$selected, 
                        levels = c(fld_threshs, nmbr_genes, c('pearson', 'spearman', 'kendall', 'cosine')),
                        labels = c(fld_threshs, nmbr_genes, c('Pearson', 'Spearman', 'Kendall', 'cosine')))
perf$dataset <- factor(perf$dataset, labels = c("Head-Neck", "Melanoma"))
perf$perf_type <- factor(perf$perf_type, levels = c('incorrect', 'erroneously classified',
                                                    'unassigned', 'intermediate', 
                                                    'correctly unclassified', 'correct'))

## Plot
PlotBars <- function(perf) {
    ggplot(data = perf, aes(x = selected, y = perf, fill = perf_type)) +
        geom_bar(position = "stack" ,stat = "identity") +
        facet_wrap(~ dataset) +
        theme_minimal() +
        scale_fill_manual(values = clrs, name = '') +
        theme(axis.text.x = element_text(angle = 60, hjust = 1),
              axis.ticks.x = element_blank()) +
        labs(x = 'method', y = 'percentage of cells') +
        scale_y_continuous(labels = c(0, 25, 50, 75, 100))
}

pl_1 <- PlotBars(perf[perf$comparison == "cors", ])
pl_2 <- PlotBars(perf[perf$comparison == "ngenes", ])
pl_3 <- PlotBars(perf[perf$comparison == "flds", ])
plts <- cowplot::plot_grid(pl_1, pl_3, pl_2, ncol = 1)
ggsave(plts, filename = paste0(output.dir, "FigureS1.pdf"), width = 12, height = 10)
