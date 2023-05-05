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

source("../../copy_number_analysis_organoids/helper_functions.R")

rank_nsegs0 <- readRDS("../../copy_number_analysis_organoids/robjects/rank_nsegments_df.RDS")
rank_ploidy0 <- readRDS("../../copy_number_analysis_organoids/robjects/rank_ploidy_df.RDS")

size_axis_text <- 8
size_text <- 12
size_legend <- 8
size_text_table_surv <- 3.5

rank_nsegs <- plot_rank(rank_nsegs0, size_labels = 3)
rank_ploidy <- plot_rank(rank_ploidy0, nudge_scalar=0.02, size_labels = 3)
rank_nsegs$coordinates$limits$y <-  c(-600, 1500)
rank_ploidy$coordinates$limits$y <-  c(-10, 15)

b <- rank_nsegs+
  theme(axis.title.x=element_text(), text = element_text(size=size_text))+
  labs(y='Number of segments', x='Rank')+guides(fill=F)#+
c <- rank_ploidy+
  theme(axis.title.x=element_text(), text = element_text(size=size_text))+
  labs(y='Ploidy', x='Rank')+guides(fill=F)

b
c

nsegments_pcawg_plot <- readRDS("../../copy_number_analysis_organoids/robjects/nsegments_pcawg_plot.RDS")

pdf("suppl8.pdf", height = 12, width = 10)
plot_grid(plot_grid(b,c, nrow=1, labels=c('a', 'b')),nsegments_pcawg_plot, nrow=2, rel_heights = c(1,2), labels=c('', 'c'))
dev.off()

