## Load libraries
library(Seurat)

data.dir <- "DEFINE"

## Load Baron ------------------------------
## https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84133 
celltypes <- vector()
for(i in 1:4) {
  mat <- read.table(paste0(data.dir, "GSM22307", 56 +i, "_human", i, "_umifm_counts.csv"), 
                    header = T, row.names = 1, sep = ',')
  celltypes <- c(celltypes, as.vector(mat$assigned_cluster))
  mat <- mat[ ,-c(1,2)]
  mat <- t(mat)
  if(i == 1) data <- mat else data <- cbind(data, mat)
}
names(celltypes) <- colnames(data)
celltypes[celltypes == 'ductal'] <- 'duct'
celltypes[celltypes == 't_cell'] <- 'T cell'
rm(mat)

## Make the seurat object
seurat <- CreateSeuratObject(raw.data = data)
seurat <- Normalizedata(seurat, normalization.method = "LogNormalize", scale.factor = 10000)
seurat <- Scaledata(seurat)
seurat <- FindVariableGenes(seurat, x.low.cutoff = 0.0125, x.high.cutoff = 4, y.cutoff = 0.5); length(seurat@var.genes)

seurat <- RunPCA(object = seurat, pc.genes = seurat@var.genes, do.print = T, pcs.print = 1:5, genes.print = 5)
seurat <- ProjectPCA(object = seurat, do.print = F)

PCElbowPlot(object = seurat)
PCHeatmap(object = seurat, pc.use = 1:20, cells.use = 500, do.balanced = TRUE, label.columns = FALSE)

seurat <- RunTSNE(object = seurat, dims.use = 1:20, do.fast = TRUE)
seurat <- AddMetadata(seurat, metadata = as.data.frame(celltypes))
TSNEPlot(object = seurat, group.by = "celltypes")

save(seurat, bcelt, file = paste0(data.dir, "Pancreas1.Rdata"))

## Muraro -------------------------------------------------------------------
## Read data + extract gene symbols
## Download at: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE85241
data <- read.table(paste0(data.dir, "GSE85241_cellsystems_dataset_4donors_updated.csv"),
                  header = T, row.names = 1)
rwns <- gsub('__.*', '', rownames(data))
keep <- !duplicated(rwns)
data <- data[keep, ]
rownames(data) <- rwns[keep]

## Celltypes
meta <- read.table(paste0(data.dir, "cell_type_annotation_Cels2016.csv"),
                    header = T, sep = '\t', stringsAsFactors = F)
mcelt <- meta[,1] ## muraro cell types
names(mcelt) <- rownames(meta);rm(meta)
data <- data[ ,intersect(names(mcelt), colnames(data))]

seurat <- CreateSeuratObject(raw.data = data)
seurat <- NormalizeData(seurat, normalization.method = "LogNormalize", scale.factor = 10000)
seurat <- ScaleData(seurat)
seurat <- FindVariableGenes(seurat, x.low.cutoff = 0.0125, x.high.cutoff = 4, y.cutoff = 0.5); length(seurat@var.genes)

seurat <- RunPCA(object = seurat, pc.genes = seurat@var.genes, do.print = T, pcs.print = 1:5, genes.print = 5)
seurat <- ProjectPCA(object = seurat, do.print = F)

PCElbowPlot(object = seurat)
PCHeatmap(object = seurat, pc.use = 1:20, cells.use = 500, do.balanced = TRUE, label.columns = FALSE)

seurat <- RunTSNE(object = seurat, dims.use = 1:20, do.fast = TRUE)
TSNEPlot(object = seurat)

patient <- colnames(seurat@raw.data)
patient <- gsub('\\..*', '', patient)
names(patient) <- colnames(seurat@raw.data)
seurat <- AddMetaData(seurat, patient, 'patient')
seurat <- AddMetaData(seurat, mcelt, 'celltypes')
TSNEPlot(seurat, group.by = 'celltypes')

save(seurat, mcelt, file = paste0(data.dir, "Pancreas2.Rdata"))
