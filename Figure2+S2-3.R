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
         'purple', 'turquoise1', 'gray80', 'gold3', 'goldenrod2', 'darkorange1', 'saddlebrown', 'deeppink3', 'black')
names(col) <- c('CD8 T cell', 'Tumor', 'CD4 T cell', 'Dendritic', 'reg. T cell', 'B cell', 'Resting Fibroblast',
                'Endothelial', 'Monocyte', 'Mast', 'T cell', 'Splits', 'NK',
                'Myocyte', 'Unknown', 'CAF', 'Myofibroblast', 'Macrophage', 'Plasma', 'pDC', 'Unassigned')
gray <- paste0("gray", rev(seq(20, 76, 4)))
names(gray) <- c(paste0("Node", 1:15))
col <- c(col, gray)

## ------------------------------------------ Melanoma
load(paste0(data.dir, "Melanoma_seurat.Rdata"))
load(paste0(data.dir, "Mel_output.Rdata")) ## Output of CHETAH, when running "Analysis_compare.R"
type <- info$classification
tsne <- seurat@dr$tsne@cell.embeddings
M_tree <- PlotTree(info, col = col, plot_limits = c(-0.5, 0.1))
ggsave(plot=M_tree, file = paste0(output.dir, 'Mela_tree.png'), height = 6, width = 8.5) ## For Figure S3

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

load(paste0(data.dir, "OV_output.Rdata")) ## Output of CHETAH, when running "Analysis_compare.R"
type <- info$classification
tsne <- seurat@dr$tsne@cell.embeddings
O_tree <- PlotTree(info, col = col)
ggsave(plot = O_tree, file = paste0(output.dir, 'Ov_tree.png'), height = 6, width = 8.5) ## For Figure S3

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

ggsave(paste0(output.dir, "Figure2"), width = 178, height = 140, units = 'mm', scale = 2) ## png
rm(final, plot_legend, legend, nnodes, nodes, col_leg)
## --------------------------------------------------------------------------------------- Classification of HN
load(paste0(data.dir, "HeadNeck_seurat.Rdata"))
load(paste0(data.dir, "HN_output.Rdata")) ## Output of CHETAH, when running "Analysis_compare.R"
type          <- info$classification
tsne <- seurat@dr$tsne@cell.embeddings
H_tree <- PlotTree(info, col = col)
ggsave(plot=H_tree, file = paste0(output.dir, 'HN_tree.png'), height = 6, width = 8.5) ## For Figure S3

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

ggsave(paste0(output.dir, "FigureS2.pdf"), width = 208, height = 70, units = 'mm', scale = 2) ## png

## ----------------------------------------- Make the river plots
library(riverplot)

## River plot function
PlotRiver <- function(data, exclude = NULL, tumor, rows, cols, y) {
    ## Discriminate cols from rows with a space
    rows = paste0(rows, " ")
    
    ## exclude data
    if (!is.null(exclude)) {
        keep <- !data$paper %in% exclude
        data <- data[keep, ]
    }
    
    ## width of flows + give correct names
    dt <- sapply(unique(data$paper), function(pap_t) {
        sapply(unique(data$type), function(ch_t) {
            sum(data$type == ch_t & data$paper == pap_t)
        })
    })*10
    colnames(dt) <- unique(data$paper)
    rownames(dt) <- paste0(unique(data$type), " ")
    dt <- reshape2::melt(dt)
    colnames(dt) <- c("N1", "N2", "Value")
    dt$N1 <- as.vector(dt$N1)
    dt$N2 <- as.vector(dt$N2)
    
    ## Set order + merge + delete empty flows
    dt <- lapply(rows, function(row) {
        sub <- dt[dt$N1 == row, ]
        sub <- sub[match(cols, sub$N2), ]
        sub
    })
    dt <- do.call(rbind, dt)
    dt <- dt[dt$Value != 0, ]
    
    ## Make node information
    nodes = data.frame(ID = c(rows, cols), 
                       x = rep(1:2, c(length(unique(dt$N1)), length(unique(dt$N2)))),
                       y = y,
                       stringsAsFactors = FALSE)
    if(is.na(y[1])) nodes <- nodes[ ,c(1,2)]
    rownames(nodes) = nodes$ID
    
    ## Make style list
    styles <- lapply(rownames(nodes), function(mn) {
        list(srt = "0", 
             col = rgb(0.8, 0.8, 0.8, 0.7), 
             edgecol = "col",
             textcex = 1.5,
             nodestyle = "invisible")
    })
    names(styles) <- rownames(nodes)
    
    ## Make final riverplot
    rp <- list(nodes = nodes, edges = dt, styles = styles)
    class(rp) <- c(class(rp), "riverplot")
    png(file = paste0(output.dir, "river_", tumor, ".png"), width = 8, height = 10, units = "in", res = 800)
        riverplot(rp, plot_area = 1, direction = "rl")
    dev.off()
}
## --------------------------------- Melanoma
rows_mel <- c("B cell", 'Plasma', 'Dendritic', 'reg. T cell', "CD4 T cell", 
              "CD8 T cell", "Intermediate", 
              'NK','Macrophage', 'Endothelial', 'CAF', 
              'Mast', 'Myofibroblast', 'Unassigned')
cols_mel <- c("B cell", 'Unknown', 'reg. T cell', "CD4 T cell", "CD8 T cell",  
              "Tumor",'T cell', 'NK',
              'Macrophage', 'Endothelial', 'CAF')
mel_y <- c(1:7, 10, 14, 17, 19:21, 23, 4:8 ,10, 14, 17, 19, 21, 23)
PlotRiver(data = data, tumor = "Melanoma", rows = rev(rows_mel), cols = rev(cols_mel), y = mel_y)

## --------------------------------- Head-Neck
rows_hn <- c("Intermediate", 'Unassigned', "B cell", 'Plasma', "CD4 T cell", 
             "CD8 T cell", 'reg. T cell', 
             'Macrophage', 'Dendritic', 'Endothelial', 'NK','CAF')
cols_hn <- c("Tumor", 'Myocyte', 'Mast', "B cell", "CD4 T cell", "CD8 T cell", 'reg. T cell',
             'Macrophage', 'Dendritic', 'Endothelial', 'Myofibroblast', 'CAF',
             "Resting Fibroblast")
hn_y <- c(2:12, 14, 1:12, 14.4)
PlotRiver(data = data3, exclude = "Unknown", tumor = "HeadNeck", 
          rows = rev(rows_hn), cols = rev(cols_hn), y = hn_y)

## --------------------------------- Ovarian
rows_ov <- c('Macrophage', "Unassigned", 'Intermediate', 
             "CD4 T cell", "CD8 T cell", 
             'NK', "B cell", 'Plasma', 'Dendritic', 'CAF')
cols_ov <- c('Macrophage', "Tumor", "Unknown", "CD4 T cell", "CD8 T cell", 
             'reg. T cell', "NK","B cell", 
             'Dendritic', 'CAF')
ov_y <- c(1:7, 8.5, 10, 12, 1:7, 8.5, 10,12)
PlotRiver(data = data2, tumor = "Ovarian", rows = rev(rows_ov), 
          cols = rev(cols_ov), y = NA)
