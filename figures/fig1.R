rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

require(ggplot2)
require(ggrepel)
require(cowplot)
require(gridExtra)
require(latex2exp)
library(cowplot)
library(ggpubr)
library(survminer)
library(survival)
library(readxl)

source("../copy_number_analysis_organoids/helper_functions.R")

exposures <- readRDS("../copy_number_analysis_organoids/robjects/exposures.RDS")
dendrograminputclr <- readRDS("../copy_number_analysis_organoids/robjects/dendrograminputclr.RDS")
ticks_in_1e = T
if(ticks_in_1e){
  heatmapinputclr <- readRDS("../copy_number_analysis_organoids/robjects/heatmapinputclr_with_ticks.RDS")
}else{
  heatmapinputclr <- readRDS("../copy_number_analysis_organoids/robjects/heatmapinputclr.RDS")
}

ascites_organoid_exposures <- readRDS("../copy_number_analysis_organoids/robjects/ascites_organoid_exposures.RDS")

size_axis_text <- 8
size_text <- 12
size_legend <- 8
size_text_table_surv <- 3.5
# size_legend_key <- 8

## survival 2
a <- createBarplot(exposures, remove_labels = FALSE, order_labels = gsub('Sample ', 'PDO', names(sort(exposures[,1])))) + 
  theme_bw()+scale_fill_brewer(palette="Dark2")+labs(y='Copy number\nsignature activity')+
  guides(fill='none')+
  # ggtitle('Exposures for the organoids')+labs(x='')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title.x=element_blank())+

  # guides(fill=F)+
  labs(x=NULL)+
  geom_bracket(data = head(melt(exposures)),
               xmin = c("PDO3"), xmax = c("PDO9"),
               y.position = c(1.1), label = c(""),
               tip.length = 0.05
  )+
  geom_bracket(data = head(melt(exposures)),
               xmin = c("PDO7"), xmax = c("PDO8"),
               y.position = c(1.1), label = c(""),
               tip.length = 0.05
  )+
  geom_bracket(data = head(melt(exposures)),
               xmin = c("PDO5"), xmax = c("PDO6"),
               y.position = c(1.15), label = c(""),
               tip.length = 0.09
  )+theme(legend.position = "bottom")

ascites_organoid_exposures$label=paste0(gsub( " .*", "", ascites_organoid_exposures$OVO4), ' (', ascites_organoid_exposures$org, ')')

ascites_organoid_exposures$label <- factor(ascites_organoid_exposures$label,
                                           levels=ascites_organoid_exposures$label[gtools::mixedorder(unique(ascites_organoid_exposures$org))])
b <- ggplot(ascites_organoid_exposures, aes(x=group, y=value, fill=variable))+
  geom_bar(stat="identity")+facet_wrap(.~label, nrow=2)+
  scale_fill_brewer(palette="Dark2")+
  theme(legend.position = "bottom")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.title.x=element_blank())+
  guides(fill='none')+
  labs(y='Copy number signature activity')

e1 <- dendrograminputclr#+geom_label_repel(label.size = NA)
e2 <- heatmapinputclr+  theme(axis.title.x=element_text(), axis.title.y=element_text(angle=90))+
  labs(x='Public tumour datasets (TCGA, PCAWG and BriTROC) and organoids', y='Copy number signature activity')
# e <- plot_grid(e1, e2, nrow = 2, labels = c('e'), rel_widths = c(3,5))
if(ticks_in_1e){
  e <- plot_grid(plot_grid(plot.new(), e1, rel_widths=c(0, 5)),
                 e2, nrow = 2, labels = c('e'))  ## adding different space in dendrogram
}else{
  e <- plot_grid(plot_grid(plot.new(), e1, rel_widths=c(0.1, 5)),
                 e2, nrow = 2, labels = c('e'))
}

e1_v2 <- e1; e1_v2$layers[[2]] = NULL; e1_v2$scales$scales[[2]] <- NULL; e1_v2$theme$plot.margin[3] = unit(0, "cm"); e1_v2 <- e1_v2+scale_y_continuous(expand=c(0,0))
e1_v2
e2_v2 <- e2+geom_label_repel(fill='white', col='black',nudge_x = 0, nudge_y = -5,
                             aes(y=0,label=ifelse(grepl('PDO', Var2) & Var1 == 's1',
                                                  as.character(Var2), NA)))+
  scale_y_continuous(expand=c(0,0))+lims(y=c(-0.6, 1))
e2_v2$theme$plot.margin[1] = unit(0, "cm")
e1_v2$theme$plot.margin[2] = unit(0, "cm")
e1_v2$theme$plot.margin[3] = unit(0, "cm")
e1_v2$theme$plot.margin[4] = unit(0, "cm")

e_v2 <- plot_grid(plot_grid(plot.new(), e1_v2, plot.new(), nrow=1, rel_widths=c(0.085, 5, 0.00)),
                  plot_grid(plot.new(), e2_v2, plot.new(), nrow=1, rel_widths=c(0.13, 5, 0.315)),
                  nrow = 2, labels = c('c'), rel_heights = c(2,5))

pdf("fig1.pdf", height = 10, width = 10, onefile = F)
plot_grid(plot_grid(a, labels='a'),
          plot_grid(plot.new(), b, labels='b', ncol=1, rel_heights = c(0.0, 0.9), vjust = -.1),
                    e_v2, label_size = 12, ncol=1, scale = 0.9, rel_heights = c(.4,0.5,0.8))
dev.off()
