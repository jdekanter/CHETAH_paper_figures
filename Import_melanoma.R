## Load libraries
library(Seurat)

data.dir <- "DEFINE"

## Load the matrix
## Download at: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE72056 
data <- read.table(paste0(data.dir, "GSE72056_melanoma_single_cell_revised_v2.txt"), 
                    row.names = NULL, header = TRUE, sep = "\t"); dim(data)

## Extract metadata
meta <- as.data.frame(t(data[1:3,-1]))
colnames(meta) <- c('tumor', 'malignant', 'celltype')
meta$celltype <- as.factor(meta$celltype)
levels(meta$celltype) <- c('Tumor', 'T cell', 'B cell', 'Macrophage', 'Endothelial', 'CAF', 'NK', 'Unknown')
meta$celltype[meta$celltype == 'Tumor' & meta$malignant == 1] <- 'Unknown'
## --- T cell subtypes have been added later, by classifying a separate T cell seurat -- ##
## Add them by:
## melanoma_ct <- readRDS(paste0(data.dir, "Melanoma_celltypes.rds"))
## meta$celltype <- melanoma_ct

## Delete duplicated rownames
data <- data[-c(1:3), ]
data <- data[!duplicated(data[,1]), ]
rownames(data) <- data[ ,1]
data <- data[ ,-c(1)]

## Delete 'unresolved' cells
data <- data[ ,meta[ ,'malignant'] != 0]

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

save(seurat, file = paste0(data.dir, 'Melanoma_seurat.Rdata'))


