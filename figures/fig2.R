rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(cowplot)
library(ggplot2)
library(ggdendro)
library(reshape2)

absCN <- readRDS("../scDNAseq-Organoids/robjects/fig2_absCN.RDS")
PDO2scDNA <- readRDS("../scDNAseq-Organoids/robjects/fig2_subclonal_hclustPDO2.RDS")
PDO3scDNA <- readRDS("../scDNAseq-Organoids/robjects/fig2_subclonal_hclustPDO3.RDS")
PDO6scDNA <- readRDS("../scDNAseq-Organoids/robjects/fig2_subclonal_hclustPDO6.RDS")
PDO2abs <- readRDS("../scDNAseq-Organoids/robjects/fig2_subclonal_hclustPDO2_2.RDS")
PDO3abs <- readRDS("../scDNAseq-Organoids/robjects/fig2_subclonal_hclustPDO3_2.RDS")
PDO6abs <- readRDS("../scDNAseq-Organoids/robjects/fig2_subclonal_hclustPDO6_2.RDS")
colours <- readRDS("../scDNAseq-Organoids/robjects/fig2_colours.RDS")
mat_breaks <- colours$mat_breaks
col_list <- colours$col_list

# give_plot = function(ph, org_it, absCN_it){
#   gA <- ggdendrogram(as.dendrogram(ph$tree_row), no.margin = TRUE)+
#     theme(axis.line=element_blank(),axis.text.x=element_blank(),
#           axis.text.y=element_blank(),axis.ticks=element_blank(),
#           plot.margin = unit(c(1, 0, 0, 0), "cm"),
#           # axis.title.x=element_blank(),
#           axis.title.y=element_blank(),legend.position="none")+scale_x_continuous(expand = c(0,0))+
#     theme(axis.line=element_blank(),axis.text.x=element_blank(),
#           axis.text.y=element_blank(),axis.ticks=element_blank(),
#           axis.title.x=element_blank(),
#           axis.title.y=element_blank(),legend.position="none",
#           panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
#           panel.grid.minor=element_blank(),plot.background=element_blank())
#   gB <- ggplot(melt(t(absCN_it)), aes(y=Var1, x=factor(Var2, levels=ph$tree_row$order), fill=value))+geom_tile()+
#     # scale_colour_steps(breaks=mat_breaks,value=col_list)
#     scale_fill_gradientn(colours=col_list, breaks=mat_breaks)+
#     theme(axis.line=element_blank(),axis.text.x=element_blank(),
#           axis.text.y=element_blank(),axis.ticks=element_blank(),
#           axis.title.x=element_blank(),
#           axis.title.y=element_blank(),legend.position="none",
#           panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
#           panel.grid.minor=element_blank(),plot.background=element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm"))
#   # maxWidth = grid::unit.pmax(gA$widths[2:3], gB$widths[2:3])
#   # gA$widths[2:3] <- as.list(maxWidth)
#   # gB$widths[2:3] <- as.list(maxWidth)
#   # grid.arrange(gA,
#   #              gB, nrow=2, top=org_it)
#   plot_grid(gA, gB, nrow=2)
# }

a <- ggplotify::as.grob(PDO2scDNA)#give_plot(PDO2scDNA,"PDO2", PDO2abs)
b <- cowplot::ggdraw()+draw_image("~/Desktop/fig2_1.png", scale = 1)
c <- ggplotify::as.grob(PDO3scDNA)#give_plot(PDO3scDNA, "PDO3", PDO3abs)
d <- cowplot::ggdraw()+draw_image("~/Desktop/fig2_2.png", scale = 1)
e <- ggplotify::as.grob(PDO6scDNA)#give_plot(PDO6scDNA, "PDO6", PDO6abs)
f <- cowplot::ggdraw()+draw_image("~/Desktop/fig2_3.png", scale = 1)

# pdf("fig2.pdf", height = 9, width = 9)
# plot_grid(plot_grid(a, b, labels = c('a', 'b')),
#           plot_grid(c, ggplotify::as.grob(d), labels = c('c', 'd')),
#           label_size = 12, ncol=1, scale = 0.9, rel_heights = c(1.5,2))
# dev.off()

pdf("fig2.pdf", height = 9, width = 9)
plot_grid(a, b, c, d, e, f, nrow=3, labels = c('a', '', 'b', '', 'c', ''))
dev.off()

# dim(PDO2abs)
# dim(PDO3abs)
# dim(PDO6abs)
# 
# PDO3abs_matched = PDO3abs[,match(colnames(PDO2abs), colnames(PDO3abs))]
# PDO6abs_matched = PDO6abs[,match(colnames(PDO2abs), colnames(PDO6abs))]
# dim(PDO2abs)
# dim(PDO3abs_matched)
# dim(PDO6abs_matched)
# 
# pairs(cbind(colMeans(PDO2abs),
#             colMeans(PDO3abs_matched),
#             colMeans(PDO6abs_matched)))
# abline(coef = c(0,1), lty='dashed', col='blue')
# 
# plot(colMeans(PDO2abs),
#      colMeans(PDO3abs_matched))
# 
# ggplot(cbind.data.frame(PDO2=colMeans(PDO2abs),
#                  PDO3=colMeans(PDO3abs_matched),
#                  PDO6=colMeans(PDO6abs_matched),
#                  chrom=gsub("\\..*", "", colnames(PDO2abs))),
#        aes(x=PDO2, y=PDO3))+geom_point()+facet_wrap(.~chrom)
