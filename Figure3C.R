args <- commandArgs(trailingOnly = TRUE)
cut_off <- as.numeric(args[1])
cat(cut_off, fill = T)

## Run script from command line with as first argument: 5, 10, 30, 50, 100, 250, 500, 1000, 2500, 5000, 10000
## Then run 2nd part of script once.

## Load data
data.dir <- "DEFINE"
output.dir <- "DEFINE"
source.dir <- "DEFINE"

load(file = paste0(data.dir, "CITE_final.Rdata")) ## CITE, ctcite
load(file = paste0(data.dir, "68K_pbmc.Rdata")) ## pbmc10x
load(file = paste0(data.dir, "pbmc_10x_clusters.Rdata")) ## ct10x, tsne10x

library(Matrix)
library(CHETAH)

## Remove ribosomal:
ribo <- read.table(paste0(data.dir, "ribosomal.txt"), header=F, sep = '\t')
pbmc10x <- pbmc10x[!rownames(pbmc10x) %in% ribo[,1], ]

## Subsample 10x pbmc ref
ref <- lapply(unique(ct10x), function(type) {
  dt <- pbmc10x[ ,ct10x == type]
  if (ncol(dt) > cut_off) {
    dt[ , sample(ncol(dt), cut_off)]
  } else  dt
})
cat(unlist(lapply(ref, ncol)), fill = T)
names(ref) <- unique(ct10x)
ref <- ref[!names(ref) %in% c("remove")]

## Classify using CHETAH and save
chetah <- CHETAHclassifier(input = CITE@raw.data,
                           ref_cells = ref)
saveRDS(chetah, file = paste0(data.dir, 'CHETAH_cutoff_', cut_off, '.rds'))

## -------------------------------------------- Figure script
if (FALSE) {
    library(ggplot2)
    load(file = paste0(data.dir, "CITE_final.Rdata")) ## CITE, ctcite
    load(file = paste0(data.dir, "68K_pbmc.Rdata")) ## pbmc10x
    load(file = paste0(data.dir, "pbmc_10x_clusters.Rdata")) ## ct10x, tsne10x
    source(paste0(source.dir, 'Calculate_performance.R'))
    
    ## load output from each run and calculate performance
    files <- list.files(data.dir, pattern = "CHETAH_cutoff_.*.rds")
    chs <- sapply(files, function(file) {
        ch <- readRDS(paste0(data.dir, file))
        CalcPerf(ref = ctcite, type = ch$classification, 
                 excl_ct = c('Cycling', 'Erythrocyte', 'doublets'), info = ch)
    })
    
    ## Prepare data
    colnames(chs) <- gsub("CHETAH_cutoff_", "", gsub(".rds", "", colnames(chs)))
    chs <- reshape2::melt(chs)
    colnames(chs) <- c('variable', 'cut_off', 'value')
    chs$variable <- factor(chs$variable, levels = c('incorrect', 'erroneously classified', 
                                                'unassigned', 'intermediate', 
                                                'correctly unclassified', 'correct'))
    chs$cut_off <- factor(chs$cut_off)
    
    clrs <- c(`correct`             ='#267c21',
            `correctly unclassified`='#40e038',
            `intermediate`          ='#BFEDBD',
            `unassigned`            ='#C1C1C1',
            `erroneously classified`='#d60202',
            `incorrect`             ='#a30914',
            `incorrect branch`      ='#f2b5b5')
    
    ## Plot
    lbs_x <- levels(chs$cut_off)
    lbs_x[length(lbs_x)] <- "all cells"
    ggplot(data = chs, aes(x = cut_off, y = value, fill = variable)) +
        geom_bar(position = "stack" ,stat = "identity") +
        theme_minimal() +
        scale_fill_manual(values = clrs, name = '') +
        theme(axis.text.x = element_text(angle = 60, hjust = 1),
              axis.ticks.x = element_blank()) +
        labs(x = 'max. cells per reference type', y = 'percentage of cells') +
        scale_y_continuous(labels = c(0, 25, 50, 75, 100)) +
        scale_x_discrete(labels = lbs_x)
    ggsave(filename = paste0(out.dir, "Figure3C.pdf", width = 10, height = 8))
}