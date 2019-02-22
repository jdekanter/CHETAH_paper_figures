data.dir <- "DEFINE"

library(reshape2)
library(ggplot2)

load(paste0(data.dir, 'CompareThreeMethods.Rdata'))
load(paste0(data.dir, 'CompareThreeMethods_CH_Node0.Rdata'))
load(paste0(data.dir, 'CompareThreeMethods_CH_0thresh.Rdata')) ## Fourth value contains the #cells of interm in Node0

## Calculate percentage
c <- lapply(1:5, function(x) {
    comparison[[x]] <- cbind.data.frame(comparison[[x]], zero_ch[[x]])
    comparison[[x]] <- rbind.data.frame(comparison[[x]], rep(0, ncol(comparison[[x]])))
    rownames(comparison[[x]]) <- c('correct', 'incorrect', 'correctly unclassified', 'unassigned', 'erroneously classified', 'intermediate')
    comparison[[x]]['intermediate', 'ch'] <- comparison[[x]]['unassigned', 'ch'] - node0[[x]][4, ]
    comparison[[x]]['unassigned', 'ch'] <- node0[[x]][4, ]
    comparison[[x]] <- comparison[[x]]/sum(comparison[[x]][ ,1])*100
    comparison[[x]] <- cbind.data.frame(comparison[[x]], rep(c('Head-Neck', 'Melanoma', 'Ovarian', 'Pancreas1', 'Pancreas2')[x], 6), rownames(comparison[[x]]))
    colnames(comparison[[x]]) <- c('SingleR', 'scmap_cell', 'scmap_cl.', 'CaSTLe', 'CHETAH', 'CHETAH_0', 'dataset', 'type')
    melt(comparison[[x]], id.vars = c('dataset', 'type'))
})

## Order correctly
c <- do.call(rbind, c)
c$dataset <- factor(c$dataset, levels = c('Head-Neck', 'Melanoma', 'Ovarian', 'Pancreas1', 'Pancreas2'))
c$variable <- factor(c$variable, levels = c('CHETAH', 'CHETAH_0', 'SingleR', 'CaSTLe', 'scmap_cell', 'scmap_cl.'))
c$type <- factor(c$type, levels = c('incorrect', 'erroneously classified', 'unassigned', 'intermediate', 'correctly unclassified', 'correct'))

## Colors
colors <- c(`correct`=               '#267c21',
            `correctly unclassified`='#40e038',
            `intermediate`          ='#BFEDBD',
            `unassigned`            ='#C1C1C1',
            `erroneously classified`='#d60202',
            `incorrect`             ='#a30914')

## Plot
ggplot(c, aes(x=variable, y=value, fill = type)) +
    theme_minimal(base_size = 14) +
    geom_col() +
    facet_wrap(~dataset, ncol = 5) +
    scale_fill_manual(values = colors, name = '') +
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 14),
          axis.ticks.x = element_blank(),
          axis.title = element_text(size = 14),
          strip.text = element_text(size = 14)) +
    labs(x = 'method', y = 'percentage of cells')

ggsave(file = 'C:\\Users\\Jurrian\\Documents\\Internship_Celseq\\Figures\\Paper_figures\\Compare_three\\NewFig4_allSplit.pdf',
       height = 6, width = 12)
