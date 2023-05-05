# S10A	/Users/morril01/Documents/PhD/other_repos/Vias_Brenton/RNASeq_and_CN/figures_GRCh37/fig3d_chromosomes.pdf
# S10B	/Users/morril01/Documents/PhD/other_repos/Vias_Brenton/RNASeq_and_CN/figures_GRCh37/fig3d_chromosomes_highlyvarsubset.pdf

rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(cowplot)
library(ggpubr)


a <- plot_grid(plot_grid(cowplot::ggdraw()+draw_image("/Users/morril01/Documents/PhD/other_repos/Vias_Brenton/RNASeq_and_CN/figures_GRCh37/fig3d_chromosomes.pdf", scale = 1),
                         rel_heights=c(5,0), nrow=2, labels = "a", label_size = 14),
               plot_grid(cowplot::ggdraw()+draw_image("/Users/morril01/Documents/PhD/other_repos/Vias_Brenton/RNASeq_and_CN/figures_GRCh37/fig3d_chromosomes_highlyvarsubset.pdf", scale = 1),
                         rel_heights=c(5,0), nrow=2, labels = "b", label_size = 14), nrow=2, rel_heights = c(2.2, 2))

plot_grid(a)

pdf(paste0("suppl10.pdf"), width = 20*0.394, height = 30*0.394)
a
dev.off()

system(" mutool draw -w 901 -o suppl10.ng suppl10.pdf")


