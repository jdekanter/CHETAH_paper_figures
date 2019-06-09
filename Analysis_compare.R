## Parameters + libraries
rm(list=ls())

data.dir <- "DEFINE"
source.dir <- "DEFINE"
output.dir <- "DEFINE"

library(SingleR)
library(scmap)
library(SingleCellExperiment)
library(tictoc)
library(irr)
library(reshape2)
library(CHETAH)
library(igraph)
library(xgboost)

## Define the comparing function ---------------------------
RunAllMethods <- function(ref, ref_ct, input, input_ct, 
                          noreftypes, origin, tsne, ch_only = F, 
                          ch_thresh = 0.1, return_raw = FALSE, 
                          sgl_ft = F, Node0 = F) {
    ## _ct = celltype
    cat(origin, '\n')
    time <- list()
    if (!sgl_ft) {
        ## Run CHETAH
        tic()
        chetah <-  CHETAHclassifier(ref_cells = ref,
                                    ref_types = ref_ct,
                                    input = input)
        chetah <- Classify(chetah, ch_thresh)
        tm <- toc()
        time$ch <- tm$toc - tm$tic
    } else {
        refSlr <- MeanRef(ref, ref_types = ref_ct)
        refSlr <-  list(name = "Ref",
                        data = refSlr,
                        types = colnames(refSlr),
                        main_types = colnames(refSlr))
        tic()
        refSlr$de.genes = CreateVariableGeneSet(refSlr$data, refSlr$types, 200) ## All standard settings
        refSlr$de.genes.main = CreateVariableGeneSet(refSlr$data, refSlr$main_types, 300)
        singler = CreateSinglerObject(input, annot = NULL, "Data", min.genes = 500,
                                      technology = "Sequencing", species = "Human", citation = "",
                                      ref.list = list(refSlr), normalize.gene.length = F, variable.genes = "de",
                                      fine.tune = T, do.signatures = F, clusters = NULL)
        singler <- singler$singler[[1]]$SingleR.single$labels[ ,1]
        cells_noref <- names(input_ct)[input_ct %in% noreftypes]
        cells_wref <- names(input_ct)[!(input_ct %in% noreftypes)]
        s <- c(0, 
               sum(input_ct[cells_wref] == singler[cells_wref]),
               length(cells_wref) - sum(input_ct[cells_wref] == singler[cells_wref]),
               0, 
               length(cells_noref))
        names(s) <- c('uncl', 'correct', 'incor', 'FN', 'FP')
        return(s)
    }
    if (!ch_only) {
        ## Run SingleR
        refSlr <- MeanRef(ref, ref_types = ref_ct)
        refSlr <-  list(name = "Ref",
                        data = refSlr,
                        types = colnames(refSlr),
                        main_types = colnames(refSlr))
        tic()
        refSlr$de.genes = CreateVariableGeneSet(refSlr$data, refSlr$types, 200) ## All standard settings
        refSlr$de.genes.main = CreateVariableGeneSet(refSlr$data, refSlr$main_types, 300)
        singler = CreateSinglerObject(input, annot = NULL, "Data", min.genes = 500,
                                      technology = "Sequencing", species = "Human", citation = "",
                                      ref.list = list(refSlr), normalize.gene.length = F, variable.genes = "de",
                                      fine.tune = T, do.signatures = F, clusters = NULL)
        singler <- singler$singler[[1]]$SingleR.single$labels[ ,1]
        tm <- toc()
        time$s <- tm$toc - tm$tic
        
        ## Run scmap
        # Prepare the data
        if (origin %in% c('HN', 'Mel')) {
            ref_raw <- 2 ^ ref - 1
            sce <- SingleCellExperiment(assays = list(counts = ref_raw), colData = data.frame('cell_type1' = ref_ct))
            logcounts(sce) <- ref
            inp <- 2 ^ input - 1
            inp <- SingleCellExperiment(assays = list(counts = inp))
            logcounts(inp) <- input
            rm(ref_raw)
        } else {
            ref_log <- log2(ref + 1)
            sce <- SingleCellExperiment(assays = list(counts = ref), colData = data.frame('cell_type1' = ref_ct))
            logcounts(sce) <- ref_log
            inplog <- log2(input + 1)
            inp <- SingleCellExperiment(assays = list(counts = input))
            logcounts(inp) <- inplog
            rm(inplog, ref_log)
        }
        rowData(sce)$feature_symbol <- rownames(ref)
        rowData(inp)$feature_symbol <- rownames(input)
        
        tic()
        sce <- selectFeatures(sce, suppress_plot = T)
        sce <- indexCluster(sce)
        ## Run scmap_Cluster
        scmapCluster_results <- scmapCluster(
            projection = inp,
            index_list = list(
                data1 = metadata(sce)$scmap_cluster_index
            )
        )
        sccluster <- as.vector(scmapCluster_results$scmap_cluster_labs)
        names(sccluster) <- colnames(input)
        tm <- toc()
        time$mapc <- tm$toc - tm$tic
        
        # Run scmap_Cell
        tic()
        sce <- selectFeatures(sce, suppress_plot = T)
        sce <- indexCell(sce)
        scmapCell_results <- scmapCell(
            inp,
            list(
                data2 = metadata(sce)$scmap_cell_index
            )
        )
        scmapCell_clusters <- scmapCell2Cluster(
            scmapCell_results,
            list(
                as.character(colData(sce)$cell_type1)
            )
        )
        sccell <- as.vector(scmapCell_clusters$scmap_cluster_labs)
        names(sccell) <- colnames(input)
        tm <- toc()
        time$map <- tm$toc - tm$tic
        
        ## Run CaSTLe
        BREAKS = c(-1, 0, 1, 6, Inf)
        nFeatures <- 100
        cell_type1 <- as.factor(ref_ct)
        tic()
        # 2. Unify sets, excluding low expressed genes
        source_n_cells_counts = apply(ref, 1, function(x) { sum(x > 0) } )
        target_n_cells_counts = apply(input, 1, function(x) { sum(x > 0) } )
        common_genes = intersect( rownames(ref)[source_n_cells_counts>10],
                                  rownames(input)[target_n_cells_counts>10]
        )
        remove(source_n_cells_counts, target_n_cells_counts)
        refc = t(ref[rownames(ref) %in% common_genes, ])
        inputc = t(input[rownames(input) %in% common_genes, ])
        ds = rbind(refc[,common_genes], inputc[,common_genes])
        isSource = c(rep(TRUE,nrow(refc)), rep(FALSE,nrow(inputc)))
        rm(inputc, refc)
        
        # 3. Highest mean in both source and target
        topFeaturesAvg = colnames(ds)[order(apply(ds, 2, mean), decreasing = T)]
        
        # for each cell - what is the most probable classification?
        L = length(levels(cell_type1))
        targetClassification = as.data.frame(matrix(rep(0,L*sum(!isSource)), nrow=L), row.names = levels(cell_type1))
        
        
        # iterate over all source cell types
        for (cellType in levels(cell_type1)) {
            
            inSourceCellType = as.factor(ifelse(cell_type1 == cellType, cellType, paste0("NOT",cellType)))
            
            # 4. Highest mutual information in source
            topFeaturesMi = names(sort(apply(ds[isSource,],2,function(x) { compare(cut(x,breaks=BREAKS),inSourceCellType,method = "nmi") }), decreasing = T))
            
            # 5. Top n genes that appear in both mi and avg
            selectedFeatures = union(head(topFeaturesAvg, nFeatures) , head(topFeaturesMi, nFeatures) )
            
            # 6. remove correlated features
            tmp = cor(ds[,selectedFeatures], method = "pearson")
            tmp[!lower.tri(tmp)] = 0
            selectedFeatures = selectedFeatures[apply(tmp,2,function(x) any(x < 0.9))]
            remove(tmp)
            
            # 7,8. Convert data from continous to binned dummy vars
            # break datasets to bins
            dsBins = apply(ds[, selectedFeatures], 2, cut, breaks= BREAKS)
            # use only bins with more than one value
            nUniq = apply(dsBins, 2, function(x) { length(unique(x)) })
            # convert to dummy vars
            ds0 = model.matrix(~ . , as.data.frame(dsBins[,nUniq>1]))
            remove(dsBins, nUniq)
            
            cat(paste0("<h2>Classifier for ",cellType,"</h2>"))
            
            inTypeSource = cell_type1 == cellType
            # 9. Classify
            xg=xgboost(data=ds0[isSource,] ,
                       label=inTypeSource,
                       objective="binary:logistic",
                       eta=0.7 , nthread=1, nround=20, verbose=0,
                       gamma=0.001, max_depth=5, min_child_weight=10)
            
            # 10. Predict
            inTypeProb = predict(xg, ds0[!isSource, ])
            
            targetClassification[cellType,] = inTypeProb
        }
        
        # use the targetClassification values to determine the predicted cell type for each cell in the target dataset
        # see "EvaluateCaSTLePerCellType.R" for an example
        castle <- apply(targetClassification, 2, which.max)
        max <- apply(targetClassification, 2, max)
        castle <- rownames(targetClassification)[castle]
        castle[max <= 0.5] <- 'unassigned'
        names(castle) <- colnames(input)
        tm <- toc()
        time$ca <- tm$toc-tm$tic
        time <- time[c('s', 'map', 'mapc', 'ca', 'ch')]
    }
    
    if (return_raw) {
        output <- list(castle, chetah, sccell, sccluster, singler)
        return(output)
    }
    ## ------------------------------------------- Quantify
    uncl <- list() ## unclassified
    cor <- list() ## correct
    incor <- list() ## incorrect
    FP <- list()
    FN <- list() ## False Negative == not classified, while reference was available
    cells_noref <- names(input_ct)[input_ct %in% noreftypes]
    cells_wref <- names(input_ct)[!(input_ct %in% noreftypes)]
    
    if (!ch_only) {
        ## SingleR
        uncl$s <- 0 ## SingleR classifies everything ## s == SingleR
        cor$s <- sum(input_ct[cells_wref] == singler[cells_wref])
        incor$s <- length(cells_wref) - cor$s
        FN$s <- 0
        FP$s <- length(cells_noref)
        
        ## scmap
        uncl$map <- sum(sccell[cells_noref] == 'unassigned')
        uncl$mapc <- sum(sccluster[cells_noref] == 'unassigned')
        uncl$ca <- sum(castle[cells_noref] == 'unassigned')
        FP$map <- length(cells_noref) - uncl$map
        FP$mapc <- length(cells_noref) - uncl$mapc
        FP$ca <- length(cells_noref) - uncl$ca
        
        # scmap_Cell
        map <- sccell[cells_wref]
        FN$map <- sum(map == 'unassigned')
        map <- map[map != 'unassigned']
        mapr <- input_ct[names(map)] ## mapr = scmap reference
        cor$map <- sum(map == mapr)
        incor$map <- sum(map != mapr)
        
        # scmap_Cluster
        map <- sccluster[cells_wref]
        FN$mapc <- sum(map == 'unassigned')
        map <- map[map != 'unassigned']
        mapr <- input_ct[names(map)]
        cor$mapc <- sum(map == mapr)
        incor$mapc <- sum(map != mapr)
        
        ## CaSTLe
        cast <- castle[cells_wref]
        FN$ca <- sum(cast == 'unassigned')
        cast <- cast[cast != 'unassigned']
        cr <- input_ct[names(cast)] ## cr = castle reference
        cor$ca <- sum(cast == cr)
        incor$ca <- sum(cast != cr)
        
    }
    ## CHETAH
    uncl$ch <- sum(grepl('Node', chetah[cells_noref]))
    FP$ch <- length(cells_noref) - uncl$ch
    ch <- chetah[cells_wref]
    if (Node0) {
        FN$ch <- sum(grepl("Node1", ch))
    } else {
        FN$ch <- sum(grepl("Node", ch))
    }
    ch <- ch[!grepl('Node', ch )] ## All types with 'Node' in it are intermediate types
    chr <- input_ct[names(ch)]
    cor$ch <- sum(ch == chr)
    incor$ch <- sum(ch != chr)
    
    ## ---------------------------
    
    output <- list(cor, incor, uncl, FN, FP)
    return(output)
} ## RunAllMethods

## -------------------------------- Run for HN data
## Load data
load(paste0(data.dir, "HeadNeck_seurat.Rdata"))
load(paste0(data.dir,"CHETAH_Tumor_reference.Rdata")) ## See https://figshare.com/s/aaf026376912366f81b6 
input_ct <- as.vector(seurat@meta.data$celltype)
names(input_ct) <- rownames(seurat@meta.data)
input <- as.matrix(seurat@raw.data)
Celldata <- lapply(Celldata, function(x) x[ ,!(colnames(x) %in% colnames(input))]); dim(Celldata[[1]])
Celldata <- Celldata[unlist(lapply(Celldata, ncol)) >= 10]
tsne <- seurat@dr$tsne@cell.embeddings
## Remove ribosomal:
ribo <- read.table(paste0(data.dir, "ribosomal.txt"), header=F, sep = '\t')
Celldata <- lapply(names(Celldata), function(x) Celldata[[x]][!rownames(Celldata[[x]]) %in% ribo[,1], ])

## Make and save CHETAH classification:
info <- CHETAHclassifier(input = input,
                         ref_cells = Celldata)
save(info, file = paste0(data.dir, "HN_output.Rdata"))

## Prepare reference
ref <- do.call(cbind, Celldata)
ref_ct <- vector()
for(i in 1:length(Celldata)) ref_ct <- c(ref_ct, rep(names(Celldata)[i], ncol(Celldata[[i]])))
names(ref_ct) <- colnames(ref)

input <- input[ ,names(input_ct)]
tsne <- tsne[names(input_ct), ]
rm(Celldata, seurat)
gc()

hn_output <- RunAllMethods(ref, ref_ct, input, input_ct, noreftypes = c('Mast', 'Tumor', 'Myocyte', 'Resting Fibroblast', 'Myofibroblast'),
                           origin = 'HN', tsne, ch_thresh = 0.1, sgl_ft = T)
hn_output_cho <- RunAllMethods(ref, ref_ct, input, input_ct, noreftypes = c('Mast', 'Tumor', 'Myocyte', 'Resting Fibroblast', 'Myofibroblast'),
                           origin = 'HN', tsne, ch_thresh = 0, sgl_ft = F, ch_only = T)
hn_output_node0 <- RunAllMethods(ref, ref_ct, input, input_ct, noreftypes = c('Mast', 'Tumor', 'Myocyte', 'Resting Fibroblast', 'Myofibroblast'),
                               origin = 'HN', tsne, ch_thresh = 0.1, sgl_ft = F, ch_only = T, Node0 = T)
rm(input_ct, input, ref, ref_ct, tsne, i)

## -------------------------------- Run for Mel data
## Load data
load(paste0(data.dir, "Melanoma_seurat.Rdata"))
load(paste0(data.dir,"CHETAH_Tumor_reference.Rdata"))
input_ct <- as.vector(seurat@meta.data$celltype)
names(input_ct) <- rownames(seurat@meta.data)
input <- as.matrix(seurat@raw.data)
Celldata <- lapply(Celldata, function(x) x[ ,!(colnames(x) %in% colnames(input))]); dim(Celldata[[1]])
Celldata <- Celldata[unlist(lapply(Celldata, ncol)) >= 10]
tsne <- seurat@dr$tsne@cell.embeddings
## Remove ribosomal:
ribo <- read.table(paste0(data.dir, "ribosomal.txt"), header=F, sep = '\t')
Celldata <- lapply(names(Celldata), function(x) Celldata[[x]][!rownames(Celldata[[x]]) %in% ribo[,1], ])

## Make and save CHETAH classification:
info <- CHETAHclassifier(input = input,
                         ref_cells = Celldata)
save(info, file = paste0(data.dir, "Mel_output.Rdata"))

## Prepare reference
ref <- do.call(cbind, Celldata)
ref_ct <- vector()
for(i in 1:length(Celldata)) ref_ct <- c(ref_ct, rep(names(Celldata)[i], ncol(Celldata[[i]])))
names(ref_ct) <- colnames(ref)

input_ct <- input_ct[!(input_ct %in% c('Unknown'))] ## Delete cells that were unclassified by the authors: those can never match
input <- input[ ,names(input_ct)]
tsne <- tsne[names(input_ct), ]
rm(Celldata, seurat)
gc()

mel_output <- RunAllMethods(ref, ref_ct, input, input_ct, noreftypes = c('Tumor'),
                            origin = 'Mel', tsne, ch_thresh = 0.1, sgl_ft = T)
mel_output_cho <- RunAllMethods(ref, ref_ct, input, input_ct, noreftypes = c('Tumor'),
                            origin = 'Mel', tsne, ch_thresh = 0, sgl_ft = T, ch_only = T)
mel_output_node0 <- RunAllMethods(ref, ref_ct, input, input_ct, noreftypes = c('Tumor'),
                            origin = 'Mel', tsne, ch_thresh = 0.1, sgl_ft = T, ch_only = T, Node0 = T)
rm(input_ct, input, ref, ref_ct, tsne, i)

## -------------------------------- Run for Ovarian
load(paste0(data.dir, "Ovarian_seurat.Rdata"))
load(paste0(data.dir,"CHETAH_Tumor_reference.Rdata"))
input_ct <- as.vector(seurat@meta.data$celltype)
names(input_ct) <- rownames(seurat@meta.data)
input <- as.matrix(seurat@raw.data)
tsne <- seurat@dr$tsne@cell.embeddings
## Remove ribosomal:
ribo <- read.table(paste0(data.dir, "ribosomal.txt"), header=F, sep = '\t')
Celldata <- lapply(names(Celldata), function(x) Celldata[[x]][!rownames(Celldata[[x]]) %in% ribo[,1], ])

## Make and save CHETAH classification:
info <- CHETAHclassifier(input = input,
                         ref_cells = Celldata)
save(info, file = paste0(data.dir, "OV_output.Rdata"))

## Prepare reference
ref <- do.call(cbind, Celldata)
ref_ct <- vector()
for(i in 1:length(Celldata)) ref_ct <- c(ref_ct, rep(names(Celldata)[i], ncol(Celldata[[i]])))
names(ref_ct) <- colnames(ref)

input_ct <- input_ct[input_ct != 'Unknown'] ## Delete Unclassified cells
input <- input[ ,names(input_ct)]
tsne <- tsne[names(input_ct), ]
rm(Celldata, seurat)
gc()

ov_output <- RunAllMethods(ref, ref_ct, input, input_ct, noreftypes = c('Tumor'),
                           origin = 'Ov', tsne, ch_thresh = 0.1, sgl_ft = T)
ov_output_cho <- RunAllMethods(ref, ref_ct, input, input_ct, noreftypes = c('Tumor'),
                           origin = 'Ov', tsne, ch_thresh = 0, sgl_ft = T, ch_only = T)
ov_output_node0 <- RunAllMethods(ref, ref_ct, input, input_ct, noreftypes = c('Tumor'),
                           origin = 'Ov', tsne, ch_thresh = 0.1, sgl_ft = T, ch_only = T, Node0 = T)
rm(input_ct, input, ref, ref_ct, tsne, i)
## -------------------------------- Run for Pancreas datasets
## Load data
load(paste0(data.dir, "Pancreas1.Rdata")); baron <- seurat
load(paste0(data.dir, "Pancreas2.Rdata")); muraro <- seurat; rm(seurat)
tsneb <- baron@dr$tsne@cell.embeddings
tsnem <- muraro@dr$tsne@cell.embeddings

## Normalize input + remove RP genes
inputm <- as.matrix(muraro@raw.data)
inputb <- as.matrix(baron@raw.data)
inputm <- t(t(inputm)/colSums(inputm) * 10000)
inputb <- t(t(inputb)/colSums(inputb) * 10000)
## Remove ribosomal:
ribo <- read.table(paste0(data.dir, "ribosomal.txt"), header=F, sep = '\t')
inputb <- inputb[!rownames(inputb) %in% ribo[,1], ]
inputm <- inputm[!rownames(inputm) %in% ribo[,1], ]

## Select cells
mcelt <- mcelt[colnames(muraro@raw.data)]
mcelt[mcelt == 'mesenchymal'] <- 'activated_stellate'
mcelt <- mcelt[mcelt != c('unclear')] ## unclassified
inputm <- inputm[ ,names(mcelt)]
mcelt[mcelt == 'pp'] <- 'gamma' ## for comparability
bcel <- bcelt
mcel <- mcelt

## remove types with low number of cells
mcelt <- mcelt[mcelt %in% names(table(mcelt))[table(mcelt) >= 10]]
bcelt <- bcelt[bcelt %in% names(table(bcelt))[table(bcelt) >= 10]]

## Run all three
muraro_output <- RunAllMethods(inputb[ ,names(bcelt)], bcelt, inputm, mcel,
                              noreftypes = c('') ,
                              origin = 'muraro', tsnem, ch_thresh = 0.1, sgl_ft = T)

baron_output <- RunAllMethods(inputm[ ,names(mcelt)], mcelt, inputb, bcel,
                               noreftypes = c('epsilon', 'macrophage', 'mast', 'schwann', 'T cell', 'quiescent_stellate'),
                               origin = 'baron', tsneb, ch_thresh = 0.1, sgl_ft = T)

muraro_output_cho <- RunAllMethods(inputb[ ,names(bcelt)], bcelt, inputm, mcel,
                               noreftypes = c('') ,
                               origin = 'muraro', tsnem, ch_thresh = 0, sgl_ft = T)

baron_output_cho <- RunAllMethods(inputm[ ,names(mcelt)], mcelt, inputb, bcel,
                              noreftypes = c('epsilon', 'macrophage', 'mast', 'schwann', 'T cell', 'quiescent_stellate'),
                              origin = 'baron', tsneb, ch_thresh = 0, sgl_ft = T)

muraro_output_node0 <- RunAllMethods(inputb[ ,names(bcelt)], bcelt, inputm, mcel,
                               noreftypes = c('') ,
                               origin = 'muraro', tsnem, ch_thresh = 0.1, sgl_ft = T, ch_only = T)

baron_output_node0 <- RunAllMethods(inputm[ ,names(mcelt)], mcelt, inputb, bcel,
                              noreftypes = c('epsilon', 'macrophage', 'mast', 'schwann', 'T cell', 'quiescent_stellate'),
                              origin = 'baron', tsneb, ch_thresh = 0.1, sgl_ft = T, ch_only = T, Node0 = T)
rm(bcel, mcel, bcelt, mcelt, muraro, baron, tsneb, tsnem, inputb, inputm)
gc()
## --------------------------------------------- Combine results and save
ch_only <- list(hn_output, mel_output, ov_output, baron_output, muraro_output)
zero_ch <- list(hn_output, mel_output, ov_output, baron_output, muraro_output)
comparison <- list(hn_output, mel_output, ov_output, baron_output, muraro_output)
comparison <- lapply(outp, function(x) {
    matrix <- do.call(rbind, lapply(x, unlist))
    rownames(matrix) <- c("correct", 'incorrect', 'correctly unclassified', 'intermediate', 'erroneously classified')
    colnames(matrix) <- c('s', 'map', 'mapc', 'ca', 'ch')
})
save(zero_ch, file = paste0(output.dir, 'CompareThreeMethods_CH_Node0.Rdata'))
save(ch_only, file = paste0(output.dir, 'CompareThreeMethods_CH_0thresh.Rdata'))
save(comparison, file = paste0(output.dir, 'CompareThreeMethods.Rdata'))