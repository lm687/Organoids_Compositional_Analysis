rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(cowplot)
library(ggplot2)
library(ggdendro)
library(reshape2)
library(dplyr)

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
b <- cowplot::ggdraw()+draw_image("../copy_number_analysis_organoids/data/absolute_profiles/PDO2.pdf", scale = 1)
c <- ggplotify::as.grob(PDO3scDNA)#give_plot(PDO3scDNA, "PDO3", PDO3abs)
d <- cowplot::ggdraw()+draw_image("../copy_number_analysis_organoids/data/absolute_profiles/PDO3.pdf", scale = 1)
e <- ggplotify::as.grob(PDO6scDNA)#give_plot(PDO6scDNA, "PDO6", PDO6abs)
f <- cowplot::ggdraw()+draw_image("../copy_number_analysis_organoids/data/absolute_profiles/PDO6.pdf", scale = 1)
g <- cowplot::ggdraw()+draw_image("fig2_chrom_annotation_crop.pdf", scale = 1)

# pdf("fig2.pdf", height = 9, width = 9)
# plot_grid(plot_grid(a, b, labels = c('a', 'b')),
#           plot_grid(c, ggplotify::as.grob(d), labels = c('c', 'd')),
#           label_size = 12, ncol=1, scale = 0.9, rel_heights = c(1.5,2))
# dev.off()

# pdf("fig2.pdf", height = 9, width = 9)
# plot_grid(a, b, c, d, e, f, nrow=4, labels = c('a', '', 'b', '', 'c', ''))
# dev.off()

pdf("fig2.pdf", height = 9, width = 9)
plot_grid(plot_grid(a, b), plot_grid(c, d), plot_grid(e, f), g, nrow=4, labels = c('a', '', 'b', '', 'c', ''), rel_heights = c(4,4,4,1))
dev.off()

## add a specific region for each organoid
specific_regions <- readRDS("../scDNAseq-Organoids/robjects/fig2_regions_subclonalb.RDS")

add_extra_rows_step = function(df){
  rbind.data.frame(c(df[1,1], -Inf, df[1,2], NA, NA, NA, df[1,7]),
                   df,
                   c(df[nrow(df),1], df[nrow(df),3], Inf, NA, NA, NA, df[nrow(df),7]))
}

ggplot(specific_regions[[1]][[1]][[1]]# %>% filter(event_confidence > 5)# %>%
       # filter(`names_cells[i]` == 10)
)+
  geom_rect(aes(xmin=start, xmax=end, ymin=0, ymax=CNval,group=cell), col='black',
            fill=NA, alpha=0.1, size=0.1)+
  geom_rect(data=as(specific_regions[[1]][[1]][[2]], "data.frame"),
            aes(xmin=start, xmax=end, ymin=0, ymax=CNval), col='red', fill=NA, alpha=0.01)+
  ggtitle(specific_regions[[1]][[2]][[1]]$title)

ggplot(specific_regions[[2]][[1]][[1]] %>% filter(event_confidence > 5, width >1) %>%
        filter(cell < 10)
)+
  # geom_rect(aes(xmin=start, xmax=end, ymin=0, ymax=CNval,group=cell), col='black',
  #           fill=NA, alpha=0.1, size=0.1)+
  # geom_rect(data=as(specific_regions[[2]][[1]][[2]], "data.frame"),
  #           aes(xmin=start, xmax=end, ymin=0, ymax=CNval), col='red', fill=NA, alpha=0.01)+
  geom_segment(aes(x=start, xend=end, y=CNval, yend=CNval,group=cell),fill='black', 
               size=0.1)+
  geom_step(aes(x=start, y=CNval,group=cell),fill='black', 
               size=0.1)+
  geom_step(data=as(specific_regions[[2]][[1]][[2]], "data.frame"),
            aes(x=start, y=CNval), col='red', fill=NA)+
  ggtitle(specific_regions[[2]][[2]][[1]]$title)

pdo2specific <- ggplot(specific_regions[[2]][[1]][[1]] %>% filter(event_confidence > 5, width >1) %>%
         group_by(cell) %>% 
         mutate(mean_CN = mean(CNval)))+
  geom_rect(aes(xmin=start, xmax=end, ymin=0, ymax=CNval,group=cell),
            fill='black', alpha=0.008, size=0.001)+
  geom_step(data=as(specific_regions[[2]][[1]][[2]], "data.frame"),
            aes(x=start, y=CNval), col='red', fill=NA)+
  ggtitle(specific_regions[[2]][[2]][[1]]$title)+
  lims(y=c(0,6))
pdo2specific2 <- ggplot(specific_regions[[2]][[1]][[1]] %>%
                          filter(event_confidence > 5, width >1), aes(x=CNval))+
  geom_density()+coord_flip()+lims(x=c(0,6))

# pdo6specific <- ggplot(specific_regions[[5]][[1]][[1]] %>%
#          filter(event_confidence > 10, width >1, end < 1.6e7))+
#   geom_rect(aes(xmin=start, xmax=end, ymin=0, ymax=CNval,group=cell),
#             fill='black', fill=NA,alpha=0.002, size=.1)+
#   geom_step(data=add_extra_rows_step(as(specific_regions[[5]][[1]][[2]], "data.frame")),
#             aes(x=start, y=CNval), col='red', fill=NA)+
#   ggtitle(specific_regions[[5]][[2]][[1]]$title)+
#   lims(x=c(min(specific_regions[[5]][[1]][[1]]$start), 1.6e7), y=c(0,4))
# # pdo6specific2 <- ggplot(specific_regions[[5]][[1]][[1]] %>%
#                          filter(event_confidence > 10, width >1, end < 1.6e7), aes(x=CNval))+
#   geom_density()+coord_flip()

ggplot(specific_regions[[7]][[1]][[1]])+
  geom_rect(aes(xmin=start, xmax=end, ymin=0, ymax=CNval,group=cell),
            fill='black', fill=NA,alpha=0.9, size=.1)+
  geom_step(data=add_extra_rows_step(as(specific_regions[[7]][[1]][[2]], "data.frame")),
            aes(x=start, y=CNval), col='red', fill=NA)+
  ggtitle(specific_regions[[7]][[2]][[1]]$title)

pdo6specific <- ggplot(specific_regions[[8]][[1]][[1]] %>%
         filter(event_confidence > 10, width >1))+
  geom_rect(aes(xmin=start, xmax=end, ymin=0, ymax=CNval,group=cell),
            fill='black', fill=NA,alpha=0.01, size=.1)+
  geom_step(data=add_extra_rows_step(as(specific_regions[[8]][[1]][[2]], "data.frame")),
            aes(x=start, y=CNval), col='red', fill=NA)+
  ggtitle(specific_regions[[8]][[2]][[1]]$title)+lims(y=c(0,4))
pdo6specific2 <- ggplot(specific_regions[[8]][[1]][[1]] %>%
                         filter(event_confidence > 10, width >1), aes(x=CNval))+
  geom_density()+coord_flip()+lims(y=c(0,4))

pdo3specific <- ggplot(specific_regions[[6]][[1]][[1]] %>%
         filter(event_confidence > 10, width>1))+
  geom_rect(aes(xmin=start, xmax=end, ymin=0, ymax=CNval,group=cell),
            fill='black', size=1, alpha=0.02)+
  geom_step(data=add_extra_rows_step(as(specific_regions[[6]][[1]][[2]], "data.frame")),
            aes(x=start, y=CNval), col='red', fill=NA)+
  ggtitle(specific_regions[[6]][[2]][[1]]$title)+lims(y=c(0,3))
pdo3specific
pdo3specific2 <- ggplot(specific_regions[[6]][[1]][[1]] %>%
                         filter(event_confidence > 10, width>1), aes(x=CNval))+
  geom_density()+coord_flip()+lims(y=c(0,3))
  
  
ggplot(specific_regions[[5]][[1]][[1]] %>% filter(event_confidence > 5, width >1),
       aes(x=CNval))+geom_density()
ggplot(specific_regions[[1]][[1]][[1]] %>% filter(event_confidence > 5, width >1),
       aes(x=CNval))+geom_density()

ggplot(specific_regions[[4]][[1]][[1]])+
  geom_rect(aes(xmin=start, xmax=end, ymin=0, ymax=CNval,group=cell),
            fill='black', alpha=0.2, size=0.1)+
  geom_step(data=as(specific_regions[[4]][[1]][[2]], "data.frame"),
            aes(x=start, y=CNval), col='red', fill=NA)+
  ggtitle(specific_regions[[4]][[2]][[1]]$title)

# outliers_scDNA = readRDS("../scDNAseq-Organoids/robjects/scDNA_outliers.RDS")
# .tree = PDO2scDNA$tree_row
# .tree = hclust(dist(PDO2abs))
# .tree$labels[cutree(.tree, k = 2) == 1] = 'Y'
# plot(.tree)

ggplot(specific_regions[[1]][[1]][[1]] %>%
         # mutate(clone=cutree(PDO2scDNA$tree_row, k = 2)[match(cell-1, 1:length(PDO2scDNA$tree_row$order))]) %>%
         filter(event_confidence > 5, width >1)# %>%
       # filter(`names_cells[i]` == 10)
)+
  # geom_segment(aes(x=start, xend=end, y=CNval, yend=CNval,group=cell),fill='black', 
  #              size=0.1)+
  geom_step(aes(x=start, y=CNval,group=cell), alpha=0.2,
            size=1)+
  geom_step(data=add_extra_rows_step(as(specific_regions[[1]][[1]][[2]], "data.frame")),
            aes(x=start, y=CNval), col='red', fill=NA)+
  ggtitle(specific_regions[[1]][[2]][[1]]$title)#+facet_wrap(.~clone)


ggplot(specific_regions[[3]][[1]][[1]] #%>% filter(event_confidence > 5)# %>%
       # filter(`names_cells[i]` == 10)
)+
  # geom_rect(aes(xmin=start, xmax=end, ymin=0, ymax=CNval,group=cell,
  #               fill=CNval),
  # geom_rect(aes(xmin=start, xmax=end, ymin=0, ymax=CNval,group=cell),fill='black', 
  #           alpha=0.3,
  #           size=0.1)+
  geom_segment(aes(x=start, xend=end, y=CNval, yend=CNval,group=cell),fill='black', 
            alpha=0.8,
            size=0.1)+
  geom_step(data=as(specific_regions[[3]][[1]][[2]], "data.frame"),
            aes(x=start, y=CNval), col='red', fill=NA)+
  ggtitle(specific_regions[[3]][[2]][[1]]$title)+
  lims(x=c(74280000, 81920000))

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

pdf("fig2_v2.pdf", height = 9, width = 11)
plot_grid(plot_grid(a, b, pdo2specific, pdo2specific2+ggtitle('')+labs(x='', y=''), nrow=1, rel_widths = c(4,4,4,1.6)),
          plot_grid(c, d, pdo3specific, pdo3specific2+ggtitle('')+labs(x='', y=''), nrow=1, rel_widths = c(4,4,4,1.6)),
          plot_grid(e, f, pdo6specific, pdo6specific2+ggtitle('')+labs(x='', y=''), nrow=1, rel_widths = c(4,4,4,1.6)),
          g, nrow=4, labels = c('a', '', 'b', '', 'c', ''), 
          rel_heights = c(4,4,4,1))
dev.off()
