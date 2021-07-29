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
heatmapinputclr <- readRDS("../copy_number_analysis_organoids/robjects/heatmapinputclr.RDS")
# rank_nsegs <- readRDS("../copy_number_analysis_organoids/robjects/rank_nsegments.RDS")
# rank_ploidy <- readRDS("../copy_number_analysis_organoids/robjects/rank_ploidy.RDS")
rank_nsegs0 <- readRDS("../copy_number_analysis_organoids/robjects/rank_nsegments_df.RDS")
rank_ploidy0 <- readRDS("../copy_number_analysis_organoids/robjects/rank_ploidy_df.RDS")
# surv_org <- readRDS("../survival_analysis/robjects/km_as_one.RDS")
survival <- read.csv("../survival_analysis/data/OrganoidSurvival.csv")
OrganoidTODO <- read_excel("../survival_analysis/data/OrganoidCulturesSurvival_091117OR.xlsx")


size_axis_text <- 8
size_text <- 12
size_legend <- 8
size_text_table_surv <- 3.5
# size_legend_key <- 8

## survival 2
km.as.one <- survfit(formula=Surv(CulturedTime, Death==1) ~Tissue, data = OrganoidTODO, conf.type = "log-log")
a <- ggsurvplot(km.as.one, conf=TRUE, risk.table = TRUE, risk.table.height = 0.3,
                pval = TRUE,
                risk.table.fontsize = size_text_table_surv,
                pval.size=size_text_table_surv, xlab="Time (days)",
                legend.labs=c("Ascites", "Solid", "Xenograft"))
coxmodel <- coxph(formula = Surv(CulturedTime, Death == 1) ~ Tissue, data = OrganoidTODO)
a$plot <- a$plot +
  ggplot2::annotate(
    "text",
    x = Inf, y = Inf,
    vjust = 3.2, hjust = 1,
    label = paste0(c('HR Solid/Ascites: ', 'HR Xenograph/Ascites: '), round(coxmodel$coefficients, 2), collapse = "\n"),
    size = 3.5
  )
a+text(cbind(x=c(0,0), y=c(0,0.2), label=paste0(c('HR Solid/Ascites: ', 'HR Xenograph/Ascites: '), round(coxmodel$coefficients, 2))))

# coxmodel <- coxph(formula = Surv(CulturedTime, Death == 1) ~ Tissue, data = survival)
# a_old <- ggsurvplot(surv_org, conf=TRUE,legend.labs = c("Ascites", "Solid", "Xenograft"), pval = TRUE,
#                 risk.table = TRUE, risk.table.fontsize = size_text_table_surv, 
#                 pval.size=size_text_table_surv, xlab="Time (days)")
# a$plot <- a$plot +
#   ggplot2::annotate(
#     "text",
#     x = Inf, y = Inf,
#     vjust = 1, hjust = 1,
#     label = paste0(c('HR Ascites/Solid: ', 'HR Ascites/Xenograph: '), round(coxmodel$coefficients, 2), collapse = "\n"),
#     size = 3.5
#   )
# a+text(cbind(x=c(0,0), y=c(0,0.4), label=paste0(c('HR Ascites/Solid: ', 'HR Ascites/Xenograph: '), round(coxmodel$coefficients, 2))))
rank_nsegs <- plot_rank(rank_nsegs0, size_labels = 3)
rank_ploidy <- plot_rank(rank_ploidy0, nudge_scalar=0.02, size_labels = 3)
rank_nsegs$coordinates$limits$y <-  c(-600, 1500)
rank_ploidy$coordinates$limits$y <-  c(-10, 15)

b <- rank_nsegs+
  theme(axis.title.x=element_text(), text = element_text(size=size_text))+
  labs(y='Number of segments', x='Rank')+guides(fill=F)#+
  # geom_label_repel(nudge_y=400*1)+lims(x=c())
  # theme(axis.text.x = element_text(size = size_axis_text), axis.text = element_text(size=size_text), legend.text = element_text(size = size_legend))#, legend.key.size = size_legend_key)
c <- rank_ploidy+
  theme(axis.title.x=element_text(), text = element_text(size=size_text))+
  labs(y='Ploidy', x='Rank')+guides(fill=F)
d <- createBarplot(exposures, remove_labels = FALSE, order_labels = gsub('Sample ', 'PDO', names(sort(exposures[,1])))) + 
  scale_fill_brewer(palette="Dark2")+labs(y='Copy number signature activity')+
  # ggtitle('Exposures for the organoids')+labs(x='')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
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
e1 <- dendrograminputclr#+geom_label_repel(label.size = NA)
e2 <- heatmapinputclr+  theme(axis.title.x=element_text(), axis.title.y=element_text(angle=90))+
  labs(x='Public tumour datasets (TCGA, PCAWG and BriTROC) and organoids', y='Copy number signature activity')
# e <- plot_grid(e1, e2, nrow = 2, labels = c('e'), rel_widths = c(3,5))
e <- plot_grid(plot_grid(plot.new(), e1, rel_widths=c(0.1, 5)),
               e2, nrow = 2, labels = c('e'))

e1_v2 <- e1; e1_v2$layers[[2]] = NULL; e1_v2$scales$scales[[2]] <- NULL; e1_v2$theme$plot.margin[3] = unit(0, "cm"); e1_v2 <- e1_v2+scale_y_continuous(expand=c(0,0))
e1_v2
e2_v2 <- e2+geom_label_repel(fill='white', col='black',nudge_x = 0, nudge_y = -5,
                       aes(y=0,label=ifelse(grepl('PDO', Var2) & Var1 == 's1',
                                        as.character(Var2), NA)))+
  scale_y_continuous(expand=c(0,0))+lims(y=c(-1.6, 1))
e2_v2$theme$plot.margin[1] = unit(0, "cm")
e_v2 <- plot_grid(plot_grid(plot.new(), e1_v2, rel_widths=c(0.1, 5), rel_heights = c(3,5)),
                  e2_v2, nrow = 2, labels = c('e'))

# pdf("fig1.pdf", height = 13, width = 9)
# plot_grid(plot_grid(plot_grid(ggplotGrob(a$plot), ggplotGrob(a$table), nrow=2, rel_heights = c(7,3.5)),
#                     d, labels = c('a', 'b')),
#           plot_grid(b, c, labels = c('c', 'd')), e, label_size = 12, ncol=1, scale = 0.9, rel_heights = c(1,.9,1.1))
# dev.off()

pdf("fig1.pdf", height = 13, width = 9)
plot_grid(plot_grid(plot_grid(ggplotGrob(a$plot), ggplotGrob(a$table), nrow=2, rel_heights = c(7,3.5)),
                    d, labels = c('a', 'b')),
          plot_grid(b, c, labels = c('c', 'd')), e_v2, label_size = 12, ncol=1, scale = 0.9, rel_heights = c(1,.9,1.1))
dev.off()

e_v2_b <- plot_grid(e2+geom_label_repel(fill='white', col='black',nudge_x = 0, nudge_y = -5,
                                        aes(y=0,label=ifelse(grepl('PDO', Var2) & Var1 == 's1',
                                                             as.character(Var2), NA)))+
                      scale_y_continuous(expand=c(0,0))+lims(y=c(-0.3, 1)))

pdf("fig1_clustering_exposures.pdf", height = 4, width = 9)
e_v2
dev.off()

e_v2_b$layers[[1]]$position
pdf("fig1_clustering_exposures.pdf", height = 4, width = 9)
e_v2_b
dev.off()

## changing the ratio of the dendrogram and barplot
e1_v2_c <- e1; e1_v2_c$layers[[2]] = NULL; e1_v2_c$scales$scales[[2]] <- NULL
e1_v2_c$theme$plot.margin[3] = unit(0, "cm")
e1_v2_c <- e1_v2_c+scale_y_continuous(expand=c(0,0))
e1_v2_c
e2_v2_c <- e2+geom_label_repel(fill='white', col='black',nudge_x = 0, nudge_y = -.8,
                             aes(y=0,label=ifelse(grepl('PDO', Var2) & Var1 == 's1',
                                                  as.character(Var2), NA)))+
  scale_y_continuous(expand=c(0,0))+lims(y=c(-0.6, 1))
e2_v2_c$theme$plot.margin[1] = unit(0, "cm")
e_v2_c <- plot_grid(plot_grid(plot.new(), e1_v2_c, rel_widths=c(0.1, 5), rel_heights = c(3,5)),
                    e2_v2_c, nrow = 2, labels = c('e'), rel_heights = c(1.5,5))
e_v2_c


pdf("fig1.pdf", height = 13, width = 9)
plot_grid(plot_grid(plot_grid(ggplotGrob(a$plot), ggplotGrob(a$table), nrow=2, rel_heights = c(7,3.5)),
                    d, labels = c('a', 'b')),
          plot_grid(b, c, labels = c('c', 'd')), e_v2_c, label_size = 12, ncol=1, scale = 0.9, rel_heights = c(1,.9,1.1))
dev.off()
