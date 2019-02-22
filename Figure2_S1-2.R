## Libraries
library(ggplot2)
library(reshape2)
library(ggdendro)
library(bioDist)
library(scales)
library(dendextend)
library(CHETAH)

data.dir <- "DEFINE"
output.dir <- "DEFINE"
source.dir <- "DEFINE"

## Select colors
col <- c('seagreen4', 'brown', 'darkolivegreen3', 'cyan3',
         'green3', 'blue', 'yellow1',
         'navy', 'orange', 'red2','mediumspringgreen', 'lightyellow3',
         'purple', 'turquoise1', 'gray80', 'gold3', 'goldenrod2', 'darkorange1', 'saddlebrown', 'deeppink3')
names(col) <- c('CD8 T cell', 'Tumor', 'CD4 T cell', 'Dendritic', 'reg. T cell', 'B cell', 'Resting Fibroblast',
                'Endothelial', 'Monocyte', 'Mast', 'T cell', 'Splits', 'NK',
                'Myocyte', 'Unknown', 'CAF', 'Myofibroblast', 'Macrophage', 'Plasma', 'pDC')
gray <- paste0("gray", rev(seq(20, 76, 4)))
names(gray) <- c(paste0("Node", 1:15))
col <- c(col, gray)

## ------------------------------------------ Melanoma
load(paste0(data.dir, "Melanoma_seurat.Rdata"))
load(paste0(data.dir, "Mel_output.Rdata")) ## Output of CHETAH, when running as in script "RunAnalysis_Figure3.R"
info <- mel_output; rm(mel_output)
type <- info$classification
tsne <- seurat@dr$tsne@cell.embeddings
M_tree <- PlotTree(info, col = col, plot_limits = c(-0.5, 0.1))
ggsave(plot=M_tree, file = paste0(output.dir, 'Mela_tree.png'), height = 6, width = 8.5) ## For Figure S2

## Make paper type
paper <- as.vector(seurat@meta.data$celltype)
names(paper) <- rownames(seurat@meta.data)
type <- type[names(paper)]

## Make dataframe
data <- data.frame(as.data.frame(tsne), "type" = type, "paper" = paper)
rm(seurat)

## ------------------------------------------ For Ovarian
load(paste0(data.dir, "Ovarian_seurat.Rdata"))
paper <- as.vector(seurat@meta.data$celltype)
names(paper) <- rownames(seurat@meta.data)

load(paste0(data.dir, "OV_output.Rdata")) ## Output of CHETAH, when running as in script "RunAnalysis_Figure3.R"
info <- ov_output; rm(ov_output)
type <- info$classification
tsne <- seurat@dr$tsne@cell.embeddings
O_tree <- PlotTree(info, col = col)
ggsave(plot = O_tree, file = paste0(output.dir, 'Ov_tree.png'), height = 6, width = 8.5)

## Make paper type
type <- type[names(paper)]
tsne <- tsne[names(paper), ]

data2 <- data.frame(as.data.frame(tsne), "type" = type, "paper" = paper)
rm(seurat)

## ----------------  Make the plot  ----------------------------- ##
plot1 <- ggplot(data, aes(x=tSNE_1, y=tSNE_2)) +
  geom_point(aes(color=type), size = 0.8) +
  scale_color_manual(values = col) +
  labs(x = "", y = "tSNE 2") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        legend.position = "none",
        plot.margin = unit(c(3,3,3,3), "pt"))

plot2 <- ggplot(data, aes(x=tSNE_1, y=tSNE_2)) +
  geom_point(aes(color=paper), size = 0.8) +
  scale_color_manual(values = col) +
  labs(x = "", y = "") +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        legend.position = "none",
        plot.margin = unit(c(3,3,3,3), "pt"))

plot3 <- ggplot(data2, aes(x=tSNE_1, y=tSNE_2)) +
  labs(x = "tSNE 1", y = "tSNE 2") +
  geom_point(aes(color=type), size = 0.8) +
  scale_color_manual(values = col) +
  theme(legend.position = "none",
        plot.margin = unit(c(3,3,3,3), "pt"))

plot4 <- ggplot(data2, aes(x=tSNE_1, y=tSNE_2)) +
  geom_point(aes(color=paper), size = 0.8) +
  scale_color_manual(values = col) +
  labs(x = "tSNE 1", y = "") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        legend.position = "none",
        plot.margin = unit(c(3,3,3,3), "pt"))

## Create óne legend
legend <- data.frame("type" = as.factor(sort(union(union(unique(data$paper), unique(data$type)),
                                        union(unique(data2$paper), unique(data2$type))))))
legend$V1 <- 1:nrow(legend)
legend$V2 <- 1:nrow(legend)
col_leg <- col[levels(legend$type)]
nnodes <- col_leg[which(!(grepl('Node', names(col_leg))))]
nodes <- col_leg[which(grepl('Node', names(col_leg)))]
nodes <- nodes[order(as.numeric(gsub("[^\\d]+", "", names(nodes), perl=TRUE)))]
col_leg <- c(nnodes, nodes)
legend$type <- factor(c(names(nnodes), names(nodes)), levels = c(names(nnodes), names(nodes)))

plot_legend <- ggplot(data = legend, aes(x = V1, y = V2)) +
  geom_point(aes(color = type)) +
  scale_color_manual(values=col_leg, name="celltypes", labels = names(col_leg)) +
  theme(legend.text = element_text(size=14),
        legend.title = element_text(size=14)) +
  guides(colour = guide_legend(override.aes = list(size = 8), ncol = 1, title = "cell types"))
library(cowplot)
leg <- get_legend(plot_legend)

## Final Plot
final <- plot_grid(plot1, plot2, plot3, plot4, ncol = 2)
final <- plot_grid(final, leg, ncol = 2, rel_widths = c(1, 0.2))

ggsave(paste0(output.dir, "OV_Mel_CHETAH.pdf"), width = 178, height = 140, units = 'mm', scale = 2) ## png
rm(final, plot_legend, legend, nnodes, nodes, col_leg)
## --------------------------------------------------------------------------------------- Classification of HN
load(paste0(data.dir, "HeadNeck_seurat.Rdata"))
load(paste0(data.dir, "HN_output.Rdata")) ## Output of CHETAH, when running as in script "RunAnalysis_Figure3.R"
info <- hn_output; rm(hn_output)
type          <- info$classification
tsne <- seurat@dr$tsne@cell.embeddings
H_tree <- PlotTree(info, col = col)
ggsave(plot=H_tree, file = paste0(output.dir, 'HN_tree.png'), height = 6, width = 8.5)

## Make paper type
paper <- as.vector(seurat@meta.data$celltype)
names(paper) <- rownames(seurat@meta.data)
type <- type[names(paper)]

## Make dataframe
data3 <- data.frame(as.data.frame(tsne), "type" = type, "paper" = paper)
rm(seurat, Celldata)

## ----------------  Make the plot  ----------------------------- ##
plot1 <- ggplot(data3, aes(x=tSNE_1, y=tSNE_2)) +
  geom_point(aes(color=type), size = 0.8) +
  labs(x = "tSNE 1", y = "tSNE 2") +
  scale_color_manual(values = col) +
  theme(legend.position = "none",
        plot.margin = unit(c(3,3,3,3), "pt"))

plot2 <- ggplot(data3, aes(x=tSNE_1, y=tSNE_2)) +
  geom_point(aes(color=paper), size = 0.8) +
  scale_color_manual(values = col) +
  labs(x = "tSNE 1", y = "") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        legend.position = "none",
        plot.margin = unit(c(3,3,3,3), "pt"))

## Create óne legend
legend <- data.frame("type" = as.factor(sort(union(unique(data3$paper), unique(data3$type)))))
legend$V1 <- 1:nrow(legend)
legend$V2 <- 1:nrow(legend)
col_leg <- col[levels(legend$type)]
nnodes <- col_leg[which(!(grepl('Node', names(col_leg))))]
nodes <- col_leg[which(grepl('Node', names(col_leg)))]
nodes <- nodes[order(as.numeric(gsub("[^\\d]+", "", names(nodes), perl=TRUE)))]
col_leg <- c(nnodes, nodes)
legend$type <- factor(c(names(nnodes), names(nodes)), levels = c(names(nnodes), names(nodes)))

plot_legend <- ggplot(data = legend, aes(x = V1, y = V2)) +
  geom_point(aes(color = type)) +
  scale_color_manual(values=col_leg, name="celltypes", labels = names(col_leg)) +
  theme(legend.text = element_text(size=14),
        legend.title = element_text(size=14)) +
  guides(colour = guide_legend(override.aes = list(size = 8), ncol = 2, title = "cell types"))
library(cowplot)
leg <- get_legend(plot_legend)

## Final plot
final <- plot_grid(plot1, plot2, ncol = 2)
final <- plot_grid(final, leg, ncol = 2, rel_widths = c(1, 0.4))

ggsave(paste0(output.dir, "HeadNeck_CHETAH.pdf"), width = 208, height = 70, units = 'mm', scale = 2) ## png