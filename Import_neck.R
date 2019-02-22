## Load libraries
library(Seurat)

data.dir <- "DEFINE"

## Load the matrix
## Download at: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103322 
data <- read.table(paste0(data.dir, "HNSCC_all_data.txt"), 
                   row.names = NULL, header = TRUE, sep = "\t"); dim(data)

## Extract metadata
meta <- as.data.frame(t(data[1:5,-1]))
colnames(meta) <- c('enzyme', 'lymph', 'malignant', 'benign', 'celltype')
meta$patient <- gsub("_.*", "", colnames(data)[-1])
levels(meta$celltype) <- c(' ', 'Tumor', 'B cell', 'Dendritic', 'Endothelial', 'Fibroblast', 
                           'Macrophage', 'Mast', 'Myocyte', 'T cell')
## --- T cell subtypes have been added later, by classifying a separate T cell seurat -- ##
Tcellsubtypes <- readRDS("HeadNeck_T_cell_subtypes.rds")
meta$celltype[names(Tcellsubtypes)] <- Tcellsubtypes

## Delete duplicated rownames
data <- data[-c(1:5), ]
data <- data[!duplicated(data[,1]), ]
rownames(data) <- data[ ,1]
data <- data[ ,-c(1)]

## Delete cells without information
data <- data[ ,!(meta$malignant == 0 & meta$benign == 0)]
data <- data[ ,meta$celltype != '']
meta <- meta[colnames(data), ]

## Change back to numeric
data <- as.matrix(data)
class(data) <- 'numeric'
data <- Matrix(data)

## Create the Seurat Object
seurat <- CreateSeuratObject(raw.data = data)
seurat <- NormalizeData(seurat, normalization.method = "LogNormalize", scale.factor = 10000)
seurat <- ScaleData(seurat)
seurat <- FindVariableGenes(seurat, x.low.cutoff = 0.0125, x.high.cutoff = 4, y.cutoff = 0.5); length(seurat@var.genes)

seurat <- RunPCA(object = seurat, pc.genes = seurat@var.genes, do.print = F, pcs.print = 1:5, genes.print = 5)
seurat <- ProjectPCA(object = seurat, do.print = F)

PCElbowPlot(object = seurat)
PCHeatmap(object = seurat, pc.use = 1:20, cells.use = 500, do.balanced = TRUE, label.columns = FALSE)

seurat <- RunTSNE(object = seurat, dims.use = 1:20, do.fast = TRUE)
TSNEPlot(object = seurat)

seurat <- AddMetaData(seurat, metadata = meta)

save(seurat, file = paste0(data.dir, 'HeadNeck_seurat.Rdata'))
