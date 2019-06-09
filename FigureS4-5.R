data.dir <- "DEFINE"
source.dir <- "DEFINE"
output.dir <- "DEFINE"

## --------------------------
## Source Utils
## --------------------------
library(CHETAH)
library(cowplot)
library(Seurat)
library(cowplot)
library(reshape2)
options(stringsAsFactors = F)
## --------------------------
## Plots
## --------------------------
boxpl <- function(seurat, type, gns, ncol =2, sub = NULL, col = NULL) {
  genes <- as.data.frame(as.matrix(seurat@raw.data[gns, ]))
  genes <- genes[ ,names(type)]
  genes <- cbind(t(genes), as.data.frame(type))
  if(!is.null(sub)) genes <- genes[genes$type %in% sub, ]
  genes <- melt(genes)

  plot <- ggplot(genes, aes(x=type, y = value)) +
    geom_jitter(aes(color = type), size = 1) +
    geom_boxplot(aes(color = type), fill = rgb(1,1,1,0.5), outlier.colour = "NA")  +
    facet_wrap(~variable, ncol = ncol) +
    theme_cowplot() +
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 10),
          axis.title.x = element_blank(),
          strip.background =element_rect(fill="white")) +
    guides(color=FALSE) +
    labs(y = 'expression')
  if (!is.null(col)) plot <- plot + scale_color_manual(values = col)
  plot
}
## Color
col <- c('seagreen4', 'brown', 'darkolivegreen3', 'cyan3',
         'green3', 'blue', 'yellow1',
         'navy', 'orange', 'red2','mediumspringgreen', 'lightyellow3',
         'purple', 'gray80', 'lemonchiffon3', 'gold3', 'goldenrod2', 'darkorange1', 'saddlebrown', 'deeppink3', 'black')
names(col) <- c('CD8 T cell', 'Tumor', 'CD4 T cell', 'Dendritic', 'reg. T cell', 'B cell', 'Resting Fibroblast',
                'Endothelial', 'Monocyte', 'Mast', 'T cell', 'Splits', 'NK',
                'Myocyte', 'Unknown', 'CAF', 'Myofibroblast', 'Macrophage', 'Plasma', 'pDC', 'Unassigned')
gray <- paste0("gray", rev(seq(20, 76, 4)))
names(gray) <- paste0("Node", 1:15)
col <- c(col, gray)

## ---------------- HN
# Neck tumor
load(paste0(data.dir, "HeadNeck_seurat.Rdata"));seurat2 <- seurat
load(paste0(data.dir, "HN_output.Rdata")) ## From "Analysis_compare.R"
info2 <- info; rm(info)
type <- info2$classification
alltype <- unique(type)
no_node <- alltype[!grepl("Node", alltype)]

plot1 <- boxpl(seurat2, type, c('MZB1', 'CD79A', 'MS4A1', 'CD3D'), ##
               col = col, ncol = 4, sub = c('Plasma'))
plot7 <- boxpl(seurat2, type, c('MZB1', 'CD79A', 'MS4A1', 'CD3D'), ##
               col = col, ncol = 4, sub = c('B cell'))
plot2 <- boxpl(seurat2, type, c('KLRF1', 'NCAM1'), ##
               col = col, ncol = 3, sub = c('NK', 'CD8 T cell', 'CD4 T cell'))

tsne1 <- PlotTSNE(info2$classification, seurat2@dr$tsne@cell.embeddings, col =  col, pt.size = 0.7) +
  theme(legend.position = 'none')
plots <- gridExtra::grid.arrange(tsne1, plot2, plot1, plot7,
                                 layout_matrix = matrix(c(1,2,1,2,1,2,3,4,3,4), ncol = 2, byrow = T))

ggsave(plots, filename = paste0(output.dir, "FigureS5.pdf"),
       width = 10.5, height = 7.4, dpi = 800)

## ---------------- Mel
load(paste0(data.dir, "Melanoma_seurat.Rdata"))
load(paste0(data.dir, "Mel_output.Rdata")) ## From "Analysis_compare.R"
tsne <- seurat@dr$tsne@cell.embeddings
type <- info$classification

## Mast boxplots
plot3 <- boxpl(seurat, type, c('TPSAB1', 'TPSB2'), col = col,
               sub = unique(type)[!(grepl("Node", unique(type)))], ncol = 3)

## CAF/MyoFibro
plot4 <- boxpl(seurat, type, c('DCN', 'LUM', 'FAP', 'ACTA2', 'MYL9', 'CDH6'),
               ncol = 6, col = col,sub = c('CAF', 'Myofibroblast', 'Split10'))

## Macro's
plot5 <- boxpl(seurat, type, c('CD14', 'TLR2', 'CCR7', 'FLT3'),
                 ncol = 4, sub = c("Macrophage", 'Dendritic'),
                 col = col)

## Plasma vs. B
plot6 <- boxpl(seurat, type, c('MZB1', 'CD79A', 'MS4A1'),
               ncol = 4, col = col,sub = c("Plasma", "B cell"))

tsne2 <- PlotTSNE(info$classification, seurat@dr$tsne@cell.embeddings, col =  col, pt.size = 0.7) +
  theme(legend.position = 'none') +
  annotate(geom = "rect",
           xmin = -20,
           xmax = 15,
           ymin = 5,
           ymax = 30,
           color = 'black', fill = rgb(0,0,0,0))

tsne3 <- PlotTSNE(info$classification, seurat@dr$tsne@cell.embeddings, col =  col, pt.size = 1) +
  theme(legend.position = 'none') + ylim(c(5, 30)) + xlim(c(-20,15))

plots <- gridExtra::grid.arrange(tsne2, tsne3, plot6, plot3, plot4, plot5,
                                 layout_matrix = matrix(c(1,2,1,2,1,2,3,4,3,4,5,6,5,6), ncol = 2, byrow = T))


ggsave(plots, filename = paste0(output.dir, "FigureS4.pdf"),
       width = 10.5, height = 11, dpi = 800)