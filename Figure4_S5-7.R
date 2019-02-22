## Libraries
rm(list=ls())

data.dir <- "DEFINE"
output.dir <- "DEFINE"

library(CHETAH)
library(scales)
library(SDMTools)

## Load data
load(paste0(data.dir, "Pancreas1.Rdata")); baron <- seurat
load(paste0(data.dir, "Pancreas2.Rdata")); muraro <- seurat; rm(seurat)
tsneb <- baron@dr$tsne@cell.embeddings
tsnem <- muraro@dr$tsne@cell.embeddings

## Normalize
inputm <- as.matrix(muraro@raw.data)
inputb <- as.matrix(baron@raw.data)
inputm <- t(t(inputm)/colSums(inputm) * 10000)
inputb <- t(t(inputb)/colSums(inputb) * 10000)

## Select cells
mcelt <- mcelt[colnames(muraro@raw.data)]
mcelt[mcelt == 'mesenchymal'] <- 'activated stellate'
bcelt[bcelt == 'activated_stellate'] <- 'activated stellate'
bcelt[bcelt == 'quiescent_stellate'] <-'quiescent stellate'
mcelt <- mcelt[mcelt != c('unclear')] ## unclassified
inputm <- inputm[ ,names(mcelt)]
mcelt[mcelt == 'pp'] <- 'gamma' ## for comparability
bcel <- bcelt
mcel <- mcelt

## remove types with low number of cells
mcelt <- mcelt[mcelt %in% names(table(mcelt))[table(mcelt) >= 10]]
bcelt <- bcelt[bcelt %in% names(table(bcelt))[table(bcelt) >= 10]]

cb <-  CHETAHclassifier(ref_cells = inputm[ ,names(mcelt)],
                        ref_types = mcelt,
                        input = inputb)

cm <-  CHETAHclassifier(ref_cells = inputb[ ,names(bcelt)],
                        ref_types = bcelt,
                        input = inputm)

## Pictures Paper
tps <- union(unique(bcel), unique(mcel))
colors <- c ('blue', 'gold', 'cyan3', 'navy',
             'forestgreen', 'orange', 'darkolivegreen3',
             'brown', 'green', 'purple','deepskyblue', 'cyan',
             'orangered3', 'coral', 'yellow3', "black",
             'yellow1', 'darkorchid1', 'darksalmon', 'darkseagreen1',
             'darkslategrey', 'deeppink4', 'green2', 'lemonchiffon1',
             'lightcyan', 'midnightblue', 'maroon1', 'orange3', 'palegreen',
             'palevioletred1', 'peru', 'seagreen1', 'red3', 'snow2',
             'steelblue1', 'turquoise')
names(colors) <- tps
gray <- paste0('gray', rev(seq(22, 88, 4)))
names(gray) <- paste0('Node', 1:17)
colors <- c(colors, gray)
colors <- colors[!is.na(names(colors))]

gradient_c1 <- c(gplots::colorpanel(n = 20, low = '#120f70', mid = '#1170dd', high = '#53cfd6'),
          gplots::colorpanel(n = 20, low = '#53cfd6', high = '#fbffb4'),
          gplots::colorpanel(n = 20, low = '#fbffb4', mid = "#ffad66", high = '#d60000'))

## ----------------------- Muraro figures ----------------------------------------------------------
M_class <- PlotTSNE(cm$classification, tsnem, col = colors) +
  labs(color='Cell types') + guides(colour = guide_legend(override.aes = list(size=6), ncol = 1))
M_score <- PlotTSNE(cm$prof_scores[[7]][ ,'duct'], tsnem, col = gradient_c1, limits = range(cm$prof_scores[[7]][ ,'duct'])*1.02) +
  labs(color='Profile score')
M_tree <- PlotTree(cm, col = colors, plot_limits = c(-0.25, 0.05)) + ggtitle('')
ggsave(plot=M_class, file = paste0(output.dir, 'Pancreas_clas.pdf'), height = 6, width = 8.5)
ggsave(plot=M_score, file = paste0(output.dir, 'Pancreas_score.pdf'), height = 6, width = 8.1) ## different dimensions for different legends
ggsave(plot=M_tree, file = paste0(output.dir, 'P1_tree.pdf'), height = 10, width = 10)

## -------------------------------------------- Heatmap
gradient_c2 <- c(gplots::colorpanel(n = 20, low = '#120f70', mid = '#56a5e2', high = '#f8ffba'),
            gplots::colorpanel(n = 10, low = '#f8ffba', high = '#ffdd7a'),
            gplots::colorpanel(n = 20, low = '#ffdd7a', mid = "#f77838", high = '#d60000'))

## Select the genes that correlate highly w the confidence score
genes <- cm$genes[[7]]$duct
classified <- mcel[mcel %in% c('acinar')]
score <- sort(cm$prof_scores[[7]][ ,'duct'], decreasing = T)
score <- score[names(score) %in% names(classified)]
data <- as.matrix(muraro@data[names(genes),names(score)])
corrs <- cor(t(data), score, method = 'spearman')
data2 <- rbind(data[corrs > 0.5 , ], data[corrs < -0.5 , ])

## correlation with transcriptcount)
cor(score, colSums(muraro@raw.data)[names(score)], method = 'spearman')

## Score silhouette
dat <- data.frame('score' = score + 1, 'cell' = 1:length(score))
M_scoresilh <- ggplot(dat, aes(x = cell, y = score)) + geom_area(fill = 'black') +
  scale_y_continuous(labels=c("0" = "-1",
                              "2" = "1"), breaks = c(0,2)) +
  theme(axis.line =  element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank())

ggsave(plot=M_scoresilh, file = paste0(output.dir, 'Pancreas1_silh.pdf'), height = 0.7, width = 8.5)

## Heatmap
dev.off()
pdf(file = paste0(output.dir, 'Heatmap_P1.pdf'), width = 9.84252, height = 7.086)
library(SDMTools)
heatmap(as.matrix(data2), Rowv = NA, Colv = NA, col = gradient_c2,
        labCol = NA)
pnts = cbind(x = c(0.9, 0.94, 0.94, 0.9), y = c(0.73, 0.73, 0.48, 0.48))
legend.gradient(pnts,
                cols = gradient_c2,
                limits = c(round(min(data), 2), round(max(data),2)),
                title = "Norm. expr."
)
dev.off()
## ----------------- Now for Baron ----------------------------------------------------------------------
## Pictures Paper

B_class <- PlotTSNE(cb$classification, tsneb, col = colors) +
  labs(color='Cell types') + guides(colour = guide_legend(override.aes = list(size=6), ncol = 1))
B_score <- PlotTSNE(cb$prof_scores[[6]][ ,'duct'], tsneb, col = gradient_c1, limits = range(cb$prof_scores[[6]][ ,'duct'])*1.02) +
  labs(color='Profile score')
B_tree <- PlotTree(cb, col = colors, plot_limits = c(-0.25, 0.08)) + ggtitle('')
ggsave(plot=B_class, file = paste0(output.dir, 'Pancreas2_clas.pdf'), height = 6, width = 8.5)
ggsave(plot=B_score, file = paste0(output.dir, 'Pancreas2_score.pdf'), height = 6, width = 8.1)
ggsave(plot=B_tree, file = paste0(output.dir, 'P2_tree.pdf'), height = 8, width = 10)

## ------------------------------------------------ Heatmap
## Select the genes that correlate highly w confidence score
genes <- cb$genes[[6]]$duct
classified <- readRDS(paste0(data.dir, 'Gradient_cells_P2.rds'))
score <- sort(cb$prof_scores[[6]][ ,'duct'], decreasing = T)
score <- score[names(score) %in% classified]
data <- as.matrix(baron@data[names(genes),names(score)])
corrs <- cor(t(data), score, method = 'spearman')
corrs[is.na(corrs)] <- 0
data2 <- rbind(data[corrs > 0.5 , ], data[corrs < -0.5 , ])

## correlation with transcriptcount)
cor(score, colSums(baron@raw.data)[names(score)], method = 'spearman')

## score silhouette
dat <- data.frame('score' = score + 1, 'cell' = 1:length(score))
B_scoresilh <- ggplot(dat, aes(x = cell, y = score)) + geom_area(fill = 'black') +
  scale_y_continuous(labels=c("0" = "-1",
                              "2" = "1"), breaks = c(0,2)) +
  theme(axis.line =  element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank())

ggsave(plot=B_scoresilh, file = paste0(output.dir, 'Pancreas2_silh.pdf'), height = 0.7, width = 8.5)

## Heatmap
dev.off()
pdf(file = paste0(output.dir, 'Heatmap_P2.pdf'), width = 9.84252, height = 7.086)
library(SDMTools)
heatmap(as.matrix(data2), Rowv = NA, Colv = NA, col = gradient_c2,
        labCol = NA)
pnts = cbind(x = c(0.9, 0.94, 0.94, 0.9), y = c(0.73, 0.73, 0.48, 0.48))
legend.gradient(pnts,
                cols = gradient_c2,
                limits = c(round(min(data), 2), round(max(data),2)),
                title = "Norm. expr."
)
dev.off()

## ------------- Marker gene plots
Markerplts <- function (data, tsne) {
  plts <- list()
    for (i in 1:ncol(data)) {
      plts[[i]] <- PlotTSNE(data[ ,i], tsne, return = T,
                            col = gplots::colorpanel(n = 50, low = 'lightgray', high = 'red')) +
        labs(color='Norm. exp') +
        theme(axis.text = element_blank(),
              axis.ticks = element_blank(),
              axis.line = element_blank())
    }
  return(plts)
}
p_m <- Markerplts(as.data.frame(as.matrix(t(muraro@data[c('RGS5', 'PDGFRA'), ]))), tsnem)
p_b <- Markerplts(as.data.frame(as.matrix(t(baron@data[c('RGS5', 'PDGFRA'), ]))), tsneb)

plts1 <- cowplot::plot_grid(plotlist = c(p_m, p_b))
ggsave(plot=plts1, file = paste0(output.dir, 'Stellate_markers.pdf'), height = 12, width = 15.8)