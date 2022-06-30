rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
setwd("../")
set.seed(234)

library(ggplot2)
# library(ggplotify)
# library(reshape2)
# library(pheatmap)
library(cowplot)
# library(dplyr)
# library(readxl)
# library(ComplexHeatmap)

exposures_org <- readRDS("../copy_number_analysis_organoids/robjects/exposures.RDS")

org_list <- rownames(exposures_org)#c('PDO1', 'PDO2', 'PDO3', 'PDO4')
plots_small_exposures=lapply(org_list, function(k){
  ggplot(cbind.data.frame(names=names(exposures_org[k,]), exposures=exposures_org[k,]), aes(y=exposures, x=1, fill=names))+geom_bar(stat = "identity")+ coord_flip()+
    theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),legend.position="none",
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank())+scale_fill_brewer(palette="Dark2")+labs(y='Exposure')
}); names(plots_small_exposures) = org_list

# a <- cowplot::ggdraw()+draw_image("~/Desktop/fig4a.png", scale = 1)
a <- plot_grid(plot_grid(cowplot::ggdraw()+draw_image("../copy_number_analysis_organoids/data/absolute_profiles/PDO1_annotated_1.pdf", scale = 1),
                         plots_small_exposures$PDO1, rel_heights=c(5,1), nrow=2),
               plot_grid(cowplot::ggdraw()+draw_image("../copy_number_analysis_organoids/data/absolute_profiles/PDO2_annotated_1.pdf", scale = 1),
                         plots_small_exposures$PDO2, rel_heights=c(5,1), nrow=2),
               plot_grid(cowplot::ggdraw()+draw_image("../copy_number_analysis_organoids/data/absolute_profiles/PDO3_annotated_1.pdf", scale = 1),
                         plots_small_exposures$PDO3, rel_heights=c(5,1), nrow=2),
               plot_grid(cowplot::ggdraw()+draw_image("../copy_number_analysis_organoids/data/absolute_profiles/PDO4_annotated_1.pdf", scale = 1),
                         plots_small_exposures$PDO4, rel_heights=c(5,1), nrow=2),
               plot_grid(cowplot::ggdraw()+draw_image("../copy_number_analysis_organoids/data/absolute_profiles/119127org_anno.pdf", scale = 1),
                         plots_small_exposures$PDO5, rel_heights=c(5,1), nrow=2),
               plot_grid(cowplot::ggdraw()+draw_image("../copy_number_analysis_organoids/data/absolute_profiles/119148org_anno.pdf", scale = 1),
                         plots_small_exposures$PDO6, rel_heights=c(5,1), nrow=2),
               plot_grid(cowplot::ggdraw()+draw_image("../copy_number_analysis_organoids/data/absolute_profiles/32077org_anno.pdf", scale = 1),
                         plots_small_exposures$PDO7, rel_heights=c(5,1), nrow=2),
               plot_grid(cowplot::ggdraw()+draw_image("../copy_number_analysis_organoids/data/absolute_profiles/54059org_anno.pdf", scale = 1),
                         plots_small_exposures$PDO8, rel_heights=c(5,1), nrow=2),
               plot_grid(cowplot::ggdraw()+draw_image("../copy_number_analysis_organoids/data/absolute_profiles/119058org_anno.pdf", scale = 1),
                         plots_small_exposures$PDO9, rel_heights=c(5,1), nrow=2),
               plot_grid(cowplot::ggdraw()+draw_image("../copy_number_analysis_organoids/data/absolute_profiles/119178org_anno.pdf", scale = 1),
                         plots_small_exposures$PDO10, rel_heights=c(5,1), nrow=2),
               plot_grid(cowplot::ggdraw()+draw_image("../copy_number_analysis_organoids/data/absolute_profiles/54327org_anno.pdf", scale = 1),
                         plots_small_exposures$PDO11, rel_heights=c(5,1), nrow=2),
               plot_grid(cowplot::ggdraw()+draw_image("../copy_number_analysis_organoids/data/absolute_profiles/151761org_anno.pdf", scale = 1),
                         plots_small_exposures$PDO12, rel_heights=c(5,1), nrow=2),
               plot_grid(cowplot::ggdraw()+draw_image("../copy_number_analysis_organoids/data/absolute_profiles/54276org_anno.pdf", scale = 1),
                         plots_small_exposures$PDO13, rel_heights=c(5,1), nrow=2),
               plot_grid(cowplot::ggdraw()+draw_image("../copy_number_analysis_organoids/data/absolute_profiles/54288org_anno.pdf", scale = 1),
                         plots_small_exposures$PDO14, rel_heights=c(5,1), nrow=2),
               plot_grid(cowplot::ggdraw()+draw_image("../copy_number_analysis_organoids/data/absolute_profiles/118947org_anno.pdf", scale = 1),
                         plots_small_exposures$PDO15, rel_heights=c(5,1), nrow=2),
               plot_grid(cowplot::ggdraw()+draw_image("../copy_number_analysis_organoids/data/absolute_profiles/119025org_anno.pdf", scale = 1),
                         plots_small_exposures$PDO16, rel_heights=c(5,1), nrow=2),
               plot_grid(cowplot::ggdraw()+draw_image("../copy_number_analysis_organoids/data/absolute_profiles/50495org_anno.pdf", scale = 1),
                         plots_small_exposures$PDO17, rel_heights=c(5,1), nrow=2),
               plot_grid(cowplot::ggdraw()+draw_image("../copy_number_analysis_organoids/data/absolute_profiles/32070org_anno.pdf", scale = 1),
                         plots_small_exposures$PDO18, rel_heights=c(5,1), nrow=2), ncol=4)


pdf("supplementary_figures/suppl7.pdf", height = 7, width = 4)
a
dev.off()
