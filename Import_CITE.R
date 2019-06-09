library(Seurat)
data.dir <- "DEFINE"

## Load data 
## https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE100866 
ADT <- as.matrix(read.csv(paste0(data.dir, "CBMC_8K_13AB_10X-ADT_umi.csv"), 
                          row.names = 1, header = TRUE, sep = ",")); dim(ADT)
rownames(ADT) <- paste0("CITE_", rownames(ADT))
RNA <- as.matrix(read.table(paste0(data.dir, "CBMC_8K_13AB_10X-RNA_umi.csv"), 
                            row.names = 1, header = TRUE, sep = ",")); dim(RNA)

## Filtering of cells and genes: 
# % of Human genes (>90% is human)
mouse <- apply(RNA[grepl("MOUSE_", rownames(RNA)), ], 2, sum)
human <- apply(RNA[grepl("HUMAN_", rownames(RNA)), ], 2, sum)
hum_perc <- human/(human + mouse)
keep_hum <- hum_perc > 0.9
# > 500 genes per cell
RNA <- RNA[grepl("HUMAN_", rownames(RNA)), ]
ngenes <- apply(RNA, 2, function(x) sum(x != 0))
keep_ngl <- ngenes > 500
# not more 3 sd above means ngenes
trans <- log10(ngenes + 1)
cut_ng <- mean(trans) + 3 * sd(trans)
keep_ngh <- trans < cut_ng

RNA <- RNA[ ,keep_ngl & keep_ngh & keep_hum]
rownames(RNA) <- gsub("HUMAN_", "", rownames(RNA))

## Plot
CITE <- CreateSeuratObject(raw.data = RNA)
CITE <- NormalizeData(CITE, normalization.method = "LogNormalize", scale.factor = 10000)
CITE <- ScaleData(CITE)
CITE <- FindVariableGenes(CITE, x.low.cutoff = 0.0125, x.high.cutoff = 4, y.cutoff = 0.5); length(CITE@var.genes)

CITE <- RunPCA(object = CITE, pc.genes = CITE@var.genes, do.print = F, pcs.print = 1:5, genes.print = 5)
CITE <- ProjectPCA(object = CITE, do.print = F)

CITE <- FindClusters(object = CITE, dims.use = 1:20, save.SNN = TRUE, resolution = 0.4)
CITE <- RunTSNE(object = CITE, dims.use = 1:20, do.fast = TRUE)
TSNEPlot(object = CITE)

## Add Assay CITE data
ADT <- ADT[ ,colnames(CITE@raw.data)]
CITE <- SetAssayData(CITE, assay.type = "CITE", slot = "raw.data", new.data = ADT)
CITE <- NormalizeData(CITE, assay.type = "CITE", normalization.method = "genesCLR")
CITE <- ScaleData(CITE, assay.type = "CITE", display.progress = FALSE)

labels <- c(
  '5' = "C16+ Macro",
  '1' = "C14+ Macro", 
  '0' = "CD4+ T cell",
  '4' = "CD8+ T cell", 
  '11' = "pDC",
  '10' = "cDC",
  '13' = "Ery",
  '8' = "Precursor",
  '2' = "NK",
  '9' = "Megakaryocyte/Platelets",
  '7' = "Mixed",
  '12' = "Cycling",
  '3' = "B cell",
  '6' = "C14+ Macro"
)
ctcite <- as.vector(CITE@ident)
names(ctcite) <- names(CITE@ident)
for(i in seq_len(length(labels))) ctcite[ctcite == names(labels)[i]] <- labels[i]
CITE <- AddMetaData(CITE, ctcite, 'celltypes')
save(CITE, ctcite, file = paste0(data.dir, "CITE_seurat.Rdata"))
