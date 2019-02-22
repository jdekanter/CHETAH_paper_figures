## Load libraries
library(Seurat)

data.dir <- "DEFINE"

## Load the matrix
## Download matlab files at: https://figshare.com/s/711d3fb2bd3288c8483 
data <- read.csv(paste0(data.dir, "ascites_scRNAseq_data_data_matrix.csv"), header = FALSE)
column <- read.csv(paste0(data.dir, "ascites_scRNAseq_data_column_information.csv"), sep = ";")
row <- read.csv(paste0(data.dir, "ascites_scRNAseq_data_row_information.csv"), stringsAsFactors = FALSE)

## Integrate the names
data <- as.matrix(data)
column <- colnames(column)
column <- column[-1] ## delete 'Samples'
row <- as.vector(as.matrix(row))
rownames(data) <- row
colnames(data) <- column
data <- data[!(duplicated(rownames(data))), ]

## Clear
rm(row, column)
gc()

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

## Load Cell types
meta <- read.csv(paste0(data.dir, 'Samples_celltypes.csv'), header = T, sep = ';')
meta <- meta[grepl('ascites', meta[ ,1]), ]
celltypes <- as.factor(meta[ ,3])
names(celltypes) <- colnames(seurat@raw.data)
levels(celltypes) <- c('Unknown', 'CD4 T cell', 'CD8 T cell', 'reg. T cell', 
                       'B cell', 'Macrophage', 'Dendritic',
                       'NK', 'CAF', 'Tumor')
seurat <- AddMetaData(seurat, metadata = celltypes, 'celltypes')

save(seurat, file = paste0(data.dir, 'Ovarian_seurat.Rdata'))
