rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(cowplot)
library(ggplot2)
library(ggdendro)
library(reshape2)
library(dplyr)
library(grid)
library(dendextend)
source("../scDNAseq-Organoids/code/helper.R")

dendextend_clades_PDO2 <- readRDS("../scDNAseq-Organoids/robjects/dendextend_clades_PDO2.RDS")
dendextend_clades_PDO3 <- readRDS("../scDNAseq-Organoids/robjects/dendextend_clades_PDO3.RDS")
dendextend_clades_PDO6 <- readRDS("../scDNAseq-Organoids/robjects/dendextend_clades_PDO6.RDS")

PDO2scDNA <- readRDS("../scDNAseq-Organoids/robjects/fig2_subclonal_hclust_nosexchromPDO2.RDS")
PDO3scDNA <- readRDS("../scDNAseq-Organoids/robjects/fig2_subclonal_hclust_nosexchromPDO3.RDS")
PDO6scDNA <- readRDS("../scDNAseq-Organoids/robjects/fig2_subclonal_hclust_nosexchromPDO6.RDS")

PDO2scDNAmat <- readRDS("../scDNAseq-Organoids/robjects/fig2_subclonal_hclust_nosexchromPDO2_2.RDS")
PDO3scDNAmat <- readRDS("../scDNAseq-Organoids/robjects/fig2_subclonal_hclust_nosexchromPDO3_2.RDS")
PDO6scDNAmat <- readRDS("../scDNAseq-Organoids/robjects/fig2_subclonal_hclust_nosexchromPDO6_2.RDS")

a <- ggplotify::as.grob(PDO2scDNA)
b <- cowplot::ggdraw()+draw_image("../copy_number_analysis_organoids/data/absolute_profiles/PDO2_annotated_2.pdf", scale = 1)
c <- ggplotify::as.grob(PDO3scDNA)
d <- cowplot::ggdraw()+draw_image("../copy_number_analysis_organoids/data/absolute_profiles/PDO3_annotated_2.pdf", scale = 1)
e <- ggplotify::as.grob(PDO6scDNA)
f <- cowplot::ggdraw()+draw_image("../copy_number_analysis_organoids/data/absolute_profiles/PDO6_annotated_2.pdf", scale = 1)
g <- cowplot::ggdraw()+draw_image("fig2_chrom_annotation_crop_nosexchrom.png", scale = 1)

rect <- rectGrob(
  x = unit(1.5, "in"),
  y = unit(.835, "npc"),
  width = unit(.04, "in"),
  height = unit(2.25, "in"),
  hjust = 0, vjust = 1,
  # gp = gpar(fill = "skyblue2", col='black', alpha = 0.9)
  gp = gpar(col='black', alpha = 0.9)
)

rect2 <- rectGrob(
  x = unit(1.30, "in"),
  y = unit(.835, "npc"),
  width = unit(.06, "in"),
  height = unit(2.25, "in"),
  hjust = 0, vjust = 1,
  # gp = gpar(fill = "skyblue2", col='black', alpha = 0.9)
  gp = gpar(col='black', alpha = 0.9)
)

rect3 <- rectGrob(
  x = unit(.95, "in"),
  y = unit(.835, "npc"),
  width = unit(.06, "in"),
  height = unit(2.25, "in"),
  hjust = 0, vjust = 1,
  # gp = gpar(fill = "skyblue2", col='black', alpha = 0.9)
  gp = gpar(col='black', alpha = 0.9)
)

rect4 <- rectGrob(
  x = unit(1.38, "in"),
  y = unit(.885, "npc"),
  width = unit(.06, "in"),
  height = unit(2.05, "in"),
  hjust = 0, vjust = 1,
  # gp = gpar(fill = "skyblue2", col='black', alpha = 0.9)
  gp = gpar(col='red', alpha = 0.9)
)

rect5 <- rectGrob(
  x = unit(1.11, "in"),
  y = unit(.885, "npc"),
  width = unit(.06, "in"),
  height = unit(2.05, "in"),
  hjust = 0, vjust = 1,
  # gp = gpar(fill = "skyblue2", col='black', alpha = 0.9)
  gp = gpar(col='red', alpha = 0.9)
)

rect6 <- rectGrob(
  x = unit(.65, "in"),
  y = unit(.885, "npc"),
  width = unit(.06, "in"),
  height = unit(2.05, "in"),
  hjust = 0, vjust = 1,
  # gp = gpar(fill = "skyblue2", col='black', alpha = 0.9)
  gp = gpar(col='red', alpha = 0.9)
)

aplus <- plot_grid(a)#+draw_grob(rect)+
  # draw_line(
  #   x = c(.22, .27),
  #   y = c(-0.05, 0.0),
  #   color = "black", size = .5,
  #   arrow = arrow(length = unit(0.02, "npc"))
  # )+draw_line(
  #   x = c(.47, .51),
  #   y = c(-0.05, 0.0),
  #   color = "red", size = .5,
  #   arrow = arrow(length = unit(0.03, "npc"))
  # )+
  # draw_line(
  #   x = c(.35, .40),
  #   y = c(-0.05, 0.0),
  #   color = "black", size = .5,
  #   arrow = arrow(length = unit(0.02, "npc"))
  # )+
  # draw_line(
  #   x = c(.26, .31),
  #   y = c(-0.05, 0.0),
  #   color = "black", size = .5,
  #   arrow = arrow(length = unit(0.02, "npc"))
  # )+
  # draw_line(
  #   x = c(.77, .82),
  #   y = c(-0.05, 0.0),
  #   color = "black", size = .5,
  #   arrow = arrow(length = unit(0.02, "npc"))
  # )

cplus <- plot_grid(c)#+draw_grob(rect2)+
  # draw_line(
  #   x = c(.22, .27),
  #   y = c(-0.05, 0.0),
  #   color = "black", size = .5,
  #   arrow = arrow(length = unit(0.02, "npc"))
  # )+draw_line(
  #   x = c(.39, .44),
  #   y = c(-0.05, 0.0),
  #   color = "red", size = .5,
  #   arrow = arrow(length = unit(0.03, "npc"))
  # )+
  # draw_line(
  #   x = c(.35, .40),
  #   y = c(-0.05, 0.0),
  #   color = "black", size = .5,
  #   arrow = arrow(length = unit(0.02, "npc"))
  # )+
  # draw_line(
  #   x = c(.26, .31),
  #   y = c(-0.05, 0.0),
  #   color = "black", size = .5,
  #   arrow = arrow(length = unit(0.02, "npc"))
  # )+
  # draw_line(
  #   x = c(.77, .82),
  #   y = c(-0.05, 0.0),
  #   color = "black", size = .5,
  #   arrow = arrow(length = unit(0.02, "npc"))
  # )+
  # draw_line(
  #   x = c(.44, .49),
  #   y = c(-0.05, 0.0),
  #   color = "black", size = .5,
  #   arrow = arrow(length = unit(0.02, "npc"))
  # )+
  # draw_line(
  #   x = c(.61, .66),
  #   y = c(-0.05, 0.0),
  #   color = "black", size = .5,
  #   arrow = arrow(length = unit(0.02, "npc"))
  # )

eplus <- plot_grid(e)#+draw_grob(rect3)+
  # draw_line(
  #   x = c(.22, .27),
  #   y = c(-0.05, 0.0),
  #   color = "black", size = .5,
  #   arrow = arrow(length = unit(0.02, "npc"))
  # )+draw_line(
  #   x = c(.28, .33),
  #   y = c(-0.05, 0.0),
  #   color = "red", size = .5,
  #   arrow = arrow(length = unit(0.03, "npc"))
  # )+
  # draw_line(
  #   x = c(.44, .49),
  #   y = c(-0.05, 0.0),
  #   color = "black", size = .5,
  #   arrow = arrow(length = unit(0.02, "npc"))
  # )+
  # draw_line(
  #   x = c(.61, .66),
  #   y = c(-0.05, 0.0),
  #   color = "black", size = .5,
  #   arrow = arrow(length = unit(0.02, "npc"))
  # )

give_ph <- function(pdo){
  require(pheatmap)

  scDNAmat <- get(paste0(pdo, 'scDNAmat'))
  mat_breaks <- c(-Inf, 0:6, seq(8, 14, by=2)) - 0.01 ## adding a small number because otherwise the binning is done wrong
  col_list <- c("#2670af", "#2670af", "#8daecf", "#eaeaea", "#ffd2b6", "#f5b788", "#f39d5f", "#f17e37", "#ee1200", "#b60000", "#7a0e04", "#401004", "#000000")
  annotation_chroms = data.frame(row.names = colnames(scDNAmat),
                                 chrom=clean_chrom(gsub("\\..*", "",
                                                        colnames(scDNAmat))))
  sexchrom_bool = annotation_chroms$chrom %in% c('X', 'Y')
  annotation_chroms_nosexchrom <- data.frame(row.names=rownames(annotation_chroms)[!sexchrom_bool],
                                             chrom=annotation_chroms$chrom[!sexchrom_bool])
  
  scDNA_pdo <- get(paste0(pdo, 'scDNA'))
  ph_nosexchrom = pheatmap::pheatmap(scDNAmat[rev(scDNA_pdo$tree_row$order),]+0.02, show_colnames = FALSE, show_rownames = FALSE,
                                     color             = col_list, cluster_rows = FALSE,
                                     breaks            = mat_breaks, cluster_cols = FALSE,
                                     annotation_col = annotation_chroms_nosexchrom, annotation_legend=F)
  return(ph_nosexchrom)
}

ph_nosexchroms <- lapply(c("PDO2", 'PDO3', 'PDO6'), give_ph)
names(ph_nosexchroms) <- c("PDO2", 'PDO3', 'PDO6')

labels(dendextend_clades_PDO2) <- NULL
labels(dendextend_clades_PDO3) <- NULL
labels(dendextend_clades_PDO6) <- NULL

dendropdo2 <- lapply(list('dendextend_clades_PDO2',
                          'dendextend_clades_PDO3',
                          'dendextend_clades_PDO6'),
                     function(i_dend_name){
   i_dend <- get(i_dend_name)
   if(i_dend_name == 'dendextend_clades_PDO2'){
     add_colour <- scale_color_manual(values=c( 'red', 'sienna3', 'grey', 'grey','grey'))
   }else if(i_dend_name == 'dendextend_clades_PDO3'){
     add_colour <- scale_color_manual(values=c( 'red', 'grey', 'grey', 'sienna3','grey',  'seagreen4', 'grey','grey'))
   }else if(i_dend_name == 'dendextend_clades_PDO6'){
     add_colour <- scale_color_manual(values=c( 'red', 'grey', 'grey', 'grey', 'sienna3',  'seagreen4', 'grey','grey'))
   }
   
   ggplot(as.ggdend(i_dend%>% set("branches_lwd", 0.5)),
          horiz = TRUE, theme = NULL)+
     theme_bw()+
     scale_x_continuous(expand = c(0, 0)) +
     scale_y_continuous(expand = c(0, 0)) +
     # annotate(label=paste0('n=', nrow(get(paste0(gsub('dendextend_clades_', '',
     #                                               i_dend_name), 'scDNAmat')))),
     #       x = c(.28), geom='text',
     #       y = c(20))+
     add_colour+
     theme(axis.line=element_blank(),axis.text.x=element_blank(),
           axis.text.y=element_blank(),axis.ticks=element_blank(),
           axis.title.x=element_blank(),
           axis.title.y=element_blank(),legend.position="none",
           panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
           panel.grid.minor=element_blank(),plot.background=element_blank())
   })

dendropdo2; names(dendropdo2) <- c("PDO2", 'PDO3', 'PDO6')

list_pdo <- c("PDO2", 'PDO3', 'PDO6')

heatmaps_trees <- function(pdo){
  if(pdo == 'PDO2'){
    extra=0.2
    extra1dendro = 0.7
  }else if(pdo == 'PDO3'){
    extra=0.0
    extra1dendro = 0.7
  }else if(pdo == 'PDO6'){
    extra=0.0
    extra1dendro = 0.65
  }
  cowplot::plot_grid(plot_grid(ph_nosexchroms[[pdo]]$gtable,
                               plot.new(), ncol=1, rel_heights=c(10, 0.2)),
                     plot_grid(#plot.new()+
                       ggplot()+ theme_void()+
                         annotate(label=paste0('n=', nrow(get(paste0(pdo, 'scDNAmat')))),
                                  x = c(0), geom='text',
                                  y = c(0)),
                               dendropdo2[[pdo]],
                               plot.new(), ncol=1,
                               rel_heights = c(extra1dendro, 9, 0.0)), rel_widths = c(0.8, 0.4))
}

heatmaps_trees_pdo <- lapply(list_pdo, heatmaps_trees)
names(heatmaps_trees_pdo) <- list_pdo


pdf("fig4clades.pdf", height = 9, width = 7)
plot_grid(plot_grid(b,#+draw_grob(rect4),
                    # aplus,
                    heatmaps_trees_pdo$PDO2,
                    nrow=1, rel_widths = c(4,4,1.6)),
          plot_grid(d,#+draw_grob(rect5),
                    # cplus, 
                    heatmaps_trees_pdo$PDO3,
                    nrow=1, rel_widths = c(4,4,1.6)),
          plot_grid(f,#+draw_grob(rect6),
                    # eplus,
                    heatmaps_trees_pdo$PDO6,
                    nrow=1, rel_widths = c(4,4,1.6)),
          g, nrow=4, labels = c('PDO2', 'PDO3', 'PDO6'), 
          rel_heights = c(4,4,4,1), rel_widths = c(2,2))
dev.off()

