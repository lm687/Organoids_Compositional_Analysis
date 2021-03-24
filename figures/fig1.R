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

source("../copy_number_analysis_organoids/helper_functions.R")

exposures <- readRDS("../copy_number_analysis_organoids/robjects/exposures.RDS")
dendrograminputclr <- readRDS("../copy_number_analysis_organoids/robjects/dendrograminputclr.RDS")
heatmapinputclr <- readRDS("../copy_number_analysis_organoids/robjects/heatmapinputclr.RDS")
rank_nsegs <- readRDS("../copy_number_analysis_organoids/robjects/rank_nsegments.RDS")
rank_ploidy <- readRDS("../copy_number_analysis_organoids/robjects/rank_ploidy.RDS")
surv_org <- readRDS("../survival_analysis/robjects/km_as_one.RDS")
survival <- read.csv("../survival_analysis/data/OrganoidSurvival.csv")

size_axis_text <- 8
size_text <- 10
size_legend <- 8
size_text_table_surv <- 3.5
# size_legend_key <- 8

a <- ggsurvplot(surv_org, conf=TRUE,legend.labs = c("Ascites", "Solid", "xenograft"), pval = TRUE, risk.table = TRUE, risk.table.fontsize = size_text_table_surv, pval.size=size_text_table_surv)
b <- rank_nsegs+labs(y='Number of segments')+guides(fill=F)
  # theme(axis.text.x = element_text(size = size_axis_text), axis.text = element_text(size=size_text), legend.text = element_text(size = size_legend))#, legend.key.size = size_legend_key)
c <- rank_ploidy+labs(y='Ploidy')+guides(fill=F)
d <- createBarplot(exposures, remove_labels = FALSE, order_labels = gsub('Sample ', 'PDO', names(sort(exposures[,1])))) + 
  scale_fill_brewer(palette="Dark2")+labs(y='Exposure')+
  # ggtitle('Exposures for the organoids')+labs(x='')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  guides(fill=F)+labs(x=NULL)+
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
  )
e1 <- dendrograminputclr#+geom_label_repel(label.size = NA)
e2 <- heatmapinputclr
e <- plot_grid(e1, e2, nrow = 2, labels = c('e'))

pdf("fig1.pdf", height = 13, width = 9)
plot_grid(plot_grid(plot_grid(ggplotGrob(a$plot), ggplotGrob(a$table), nrow=2, rel_heights = c(7,3.5)),
                    d, labels = c('a', 'b')),
          plot_grid(b, c, labels = c('c', 'd')), e, label_size = 12, ncol=1, scale = 0.9, rel_heights = c(1,.9,1.1))
dev.off()

