## Load libraries
library(Seurat)
data.dir <- "DEFINE"

## Load data
## http://s3-us-west-2.amazonaws.com/10x.files/samples/cell/pbmc68k_rds/pbmc68k_data.rds 
## More information at: https://github.com/10XGenomics/single-cell-3prime-paper/tree/master/pbmc68k_analysis 
pbmc <- readRDS(paste0(data.dir, "pbmc68k_data.rds"))
pbmc <- pbmc$all_data$`17820`$hg19
data <- t(pbmc$mat)
colnames(data) <- pbmc$barcodes
rownames(data) <- pbmc$gene_symbols

## Normalize & save
data <- t(as.matrix(data)/colSums(data))
data <- t(data) * 1000
data <- log2(data + 1)
data <- Matrix(data)
save(pbmc10x, file = paste0(data.dir, "pbmc_10x_counts.Rdata"))

## Make Seurat 
data <- data[!duplicated(rownames(data)), ]
pbmc <- CreateSeuratObject(raw.data = data)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 100000)
pbmc <- ScaleData(pbmc)
pbmc <- FindVariableGenes(pbmc, x.low.cutoff = 0.0125, x.high.cutoff = 4, y.cutoff = 0.5); length(pbmc@var.genes)

pbmc <- RunPCA(object = pbmc, pc.genes = pbmc@var.genes, do.print = F, pcs.print = 1:5, genes.print = 5)
pbmc <- ProjectPCA(object = pbmc, do.print = F)

pbmc <- FindClusters(object = pbmc, dims.use = 1:20, resolution = 1, print.output = 0, save.SNN = TRUE)
pbmc <- RunTSNE(object = pbmc, dims.use = 1:20, do.fast = TRUE)  ## Not needed for figures

## Assign cell labels and save
labels <- c(
  '0' = "CD4+ active",
  '1' = "CD4+ T cell", 
  '2' = "CD4+ T cell",
  '3' = "CD8+ GZMK",
  '4' = "CD8+ naive",
  '5' = "CD8+ GNLY",
  '6' = "NK",
  '7' = "B cell",
  '8' = "CD14+ Macro",
  '9' = "CD16+ Macro",
  '10' = "CD4+ memory",
  '11' = "NK/CD8",
  '12' = "cDC",
  '13' = "pDC", 
  '14' = "Megak.",
  '15' = "remove",
  '16' = "remove",
  '17' = "remove",
  '18' = "remove",
  '19' = "DC CLEC9A",
  '20' = "remove",
  '21' = "remove",
  '22' = "HSC"
)
ct10x <- as.vector(pbmc@ident)
names(ct10x) <- names(pbmc@ident)
for(i in 1:length(labels)) ct10x[ct10x == names(labels[i])] <- labels[i]
tsne10x <- seurat@dr$tsne@cell.embeddings
save(ct10x,tsne10x, file = paste0(data.dir, "pbmc_10x_clusters.Rdata"))
## Don't save seurat: not needed for figures