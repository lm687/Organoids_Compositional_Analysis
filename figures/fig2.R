rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(cowplot)
library(ggplot2)
library(ggdendro)
library(reshape2)
library(dplyr)
library(grid)

absCN <- readRDS("../scDNAseq-Organoids/robjects/fig2_absCN.RDS")
PDO2scDNA <- readRDS("../scDNAseq-Organoids/robjects/fig2_subclonal_hclust_nosexchromPDO2.RDS")
PDO3scDNA <- readRDS("../scDNAseq-Organoids/robjects/fig2_subclonal_hclust_nosexchromPDO3.RDS")
PDO6scDNA <- readRDS("../scDNAseq-Organoids/robjects/fig2_subclonal_hclust_nosexchromPDO6.RDS")
PDO2abs <- readRDS("../scDNAseq-Organoids/robjects/fig2_subclonal_hclust_nosexchromPDO2_2.RDS")
PDO3abs <- readRDS("../scDNAseq-Organoids/robjects/fig2_subclonal_hclust_nosexchromPDO3_2.RDS")
PDO6abs <- readRDS("../scDNAseq-Organoids/robjects/fig2_subclonal_hclust_nosexchromPDO6_2.RDS")
colours <- readRDS("../scDNAseq-Organoids/robjects/fig2_colours.RDS")
mat_breaks <- colours$mat_breaks
col_list <- colours$col_list

# number of cells
dim(PDO2abs)
dim(PDO3abs)
dim(PDO6abs)


PDO2scDNA$gtable$grobs[[1]]$label = ""
PDO2scDNA$gtable$grobs[[1]]$label
PDO3scDNA$gtable$grobs[[1]]$label = ""
PDO6scDNA$gtable$grobs[[1]]$label = ""

a <- ggplotify::as.grob(PDO2scDNA)#give_plot(PDO2scDNA,"PDO2", PDO2abs)
b <- cowplot::ggdraw()+draw_image("../copy_number_analysis_organoids/data/absolute_profiles/PDO2_annotated_2.pdf", scale = 1)
c <- ggplotify::as.grob(PDO3scDNA)#give_plot(PDO3scDNA, "PDO3", PDO3abs)
d <- cowplot::ggdraw()+draw_image("../copy_number_analysis_organoids/data/absolute_profiles/PDO3_annotated_2.pdf", scale = 1)
e <- ggplotify::as.grob(PDO6scDNA)#give_plot(PDO6scDNA, "PDO6", PDO6abs)
f <- cowplot::ggdraw()+draw_image("../copy_number_analysis_organoids/data/absolute_profiles/PDO6_annotated_2.pdf", scale = 1)
g <- cowplot::ggdraw()+draw_image("fig2_chrom_annotation_crop_nosexchrom.png", scale = 1)


pdf("fig2.pdf", height = 9, width = 9)
plot_grid(plot_grid(a, b), plot_grid(c, d), plot_grid(e, f), g, nrow=4, labels = c('a', '', 'b', '', 'c', ''), rel_heights = c(4,4,4,1))
dev.off()

## add a specific region for each organoid
# specific_regions <- readRDS("../scDNAseq-Organoids/robjects/fig2_regions_subclonalb.RDS")
# 
# add_extra_rows_step = function(df){
#   rbind.data.frame(c(df[1,1], -Inf, df[1,2], NA, NA, NA, df[1,7]),
#                    df,
#                    c(df[nrow(df),1], df[nrow(df),3], Inf, NA, NA, NA, df[nrow(df),7]))
# }
# 
# 
# pdo2specific <- ggplot(specific_regions[[2]][[1]][[1]] %>% filter(event_confidence > 5, width >1) %>%
#          group_by(cell) %>% 
#          mutate(mean_CN = mean(CNval)))+
#   geom_rect(aes(xmin=start, xmax=end, ymin=0, ymax=CNval,group=cell),
#             fill='black', alpha=0.008, size=0.001)+
#   geom_step(data=as(specific_regions[[2]][[1]][[2]], "data.frame"),
#             aes(x=start, y=CNval), col='red', fill=NA)+
#   # ggtitle(specific_regions[[2]][[2]][[1]]$title)+
#   ggtitle("Chromosome 7")+
#   lims(y=c(0,6))+labs(x='Chromosome location', y='Absolute tumour CN')
# pdo2specific2 <- ggplot(specific_regions[[2]][[1]][[1]] %>%
#                           filter(event_confidence > 5, width >1), aes(x=CNval))+
#   geom_density()+coord_flip()+lims(x=c(0,6))+labs(x='Absolute tumour CN', y='Density')
# 
# pdo6specific <- ggplot(specific_regions[[8]][[1]][[1]] %>%
#          filter(event_confidence > 10, width >1))+
#   geom_rect(aes(xmin=start, xmax=end, ymin=0, ymax=CNval,group=cell),
#             fill='black', fill=NA,alpha=0.01, size=.1)+
#   geom_step(data=add_extra_rows_step(as(specific_regions[[8]][[1]][[2]], "data.frame")),
#             aes(x=start, y=CNval), col='red', fill=NA)+
#   ggtitle(specific_regions[[8]][[2]][[1]]$title)+lims(y=c(0,4))+
#   labs(x='Chromosome location', y='Absolute tumour CN')+
#   scale_x_continuous(breaks = (seq(from = min(specific_regions[[8]][[1]][[1]]$start),
#                                    to = max(specific_regions[[8]][[1]][[1]]$end), length.out = 3)),
#                      labels = function(x) format(x, scientific = TRUE))
# 
# pdo6specific2 <- ggplot(specific_regions[[8]][[1]][[1]] %>%
#                          filter(event_confidence > 10, width >1), aes(x=CNval))+
#   geom_density()+coord_flip()+lims(y=c(0,4))+labs(x='Absolute tumour CN', y='Density')
# 
# pdo3specific <- ggplot(specific_regions[[6]][[1]][[1]] %>%
#          filter(event_confidence > 10, width>1))+
#   geom_rect(aes(xmin=start, xmax=end, ymin=0, ymax=CNval,group=cell),
#             fill='black', size=1, alpha=0.02)+
#   geom_step(data=add_extra_rows_step(as(specific_regions[[6]][[1]][[2]], "data.frame")),
#             aes(x=start, y=CNval), col='red', fill=NA)+
#   # ggtitle(specific_regions[[6]][[2]][[1]]$title)+
#   ggtitle('Start of chromosome 5')+
#   lims(y=c(0,3))+
#   labs(x='Chromosome location', y='Absolute tumour CN')+
#   scale_x_continuous(breaks = (seq(from = min(specific_regions[[6]][[1]][[1]]$start),
#                                    to = max(specific_regions[[6]][[1]][[1]]$end), length.out = 3)),
#                      labels = function(x) format(x, scientific = TRUE))
# 
# pdo3specificb <- ggplot(specific_regions[[6]][[1]][[1]] %>%
#                          filter(event_confidence > 10, width>1))+
#   geom_point(aes(x=mean(c(start,end)), y=CNval),
#             fill='black', size=1, alpha=0.02)+
#   geom_step(data=add_extra_rows_step(as(specific_regions[[6]][[1]][[2]], "data.frame")),
#             aes(x=start, y=CNval), col='red', fill=NA)+
#   # ggtitle(specific_regions[[6]][[2]][[1]]$title)+
#   ggtitle('Start of chromosome 5')+
#   lims(y=c(0,3))+
#   labs(x='Chromosome location', y='Absolute tumour CN')+
#   scale_x_continuous(breaks = (seq(from = min(specific_regions[[6]][[1]][[1]]$start),
#                                    to = max(specific_regions[[6]][[1]][[1]]$end), length.out = 3)),
#                      labels = function(x) format(x, scientific = TRUE))
# 
# pdo3specific2 <- ggplot(specific_regions[[6]][[1]][[1]] %>%
#                          filter(event_confidence > 10, width>1), aes(x=CNval))+
#   geom_density()+coord_flip()+lims(y=c(0,3))+labs(x='Absolute tumour CN', y='Density')

hist(specific_regions[[6]][[1]][[1]]$CNval, xlim = c(0, 4))
table(specific_regions[[6]][[1]][[1]]$CNval == 1)
table(droplevels(specific_regions[[6]][[1]][[1]]$seqnames))
specific_regions[[6]][[1]][[1]]$start[1]
specific_regions[[6]][[1]][[1]]$end[1]

pdo2specific2
PDO6abs

plot(density(PDO2abs[,grepl('X7.',colnames(PDO6abs))][,39]))
plot(density(PDO2abs[,"X7.5.120.001.10.240.000"]))
pdo2specific2 <- ggplot(melt(PDO2abs[,"X7.5.120.001.10.240.000"]), aes(x=value))+
  geom_density()+
  coord_flip()+lims(x=c(0,4))+labs(x='Absolute tumour CN', y='Density')

table(round(PDO2abs[,"X7.5.120.001.10.240.000"]))/nrow(PDO2abs)


plot(density(PDO3abs[, "X5.56.320.001.61.440.000"    ]))
pdo3specific2 <- ggplot(melt(PDO3abs[,"X5.56.320.001.61.440.000"]), aes(x=value))+
  geom_density()+
  coord_flip()+lims(x=c(0,4))+labs(x='Absolute tumour CN', y='Density')
table(round(PDO3abs[,"X5.56.320.001.61.440.000"]))/nrow(PDO3abs)

plot(density(PDO6abs[,grepl('X2.',colnames(PDO6abs))][,39]))
pdo6specific2 <- ggplot(melt(PDO6abs[,"X2.199.680.001.204.800.000"]), aes(x=value))+
  geom_density()+
  coord_flip()+lims(x=c(0,4))+labs(x='Absolute tumour CN', y='Density')
table(round(PDO6abs[,"X2.199.680.001.204.800.000"]))/nrow(PDO6abs)


pdf("fig2_v2.pdf", height = 9, width = 11)
plot_grid(plot_grid(b, a, pdo2specific, pdo2specific2+ggtitle(''), nrow=1, rel_widths = c(4,4,4,1.6)),
          plot_grid(d, c, pdo3specific, pdo3specific2+ggtitle(''), nrow=1, rel_widths = c(4,4,4,1.6)),
          plot_grid(f, e, pdo6specific, pdo6specific2+ggtitle(''), nrow=1, rel_widths = c(4,4,4,1.6)),
          g, nrow=4, labels = c('PDO2', 'PDO3', 'PDO4'), 
          rel_heights = c(4,4,4,1))
dev.off()


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

pdf("fig2_v3.pdf", height = 9, width = 7)
plot_grid(plot_grid(b+draw_grob(rect4), plot_grid(a)+draw_grob(rect)+
                      draw_line(
                        x = c(.47, .51),
                        y = c(-0.05, 0.0),
                        color = "red", size = .5,
                        arrow = arrow(length = unit(0.03, "npc"))
                      ),
                      # draw_image(img , x = 1, y = 1, scale = .9),
               pdo2specific2+ggtitle(''), nrow=1, rel_widths = c(4,4,1.6)),
          plot_grid(d+draw_grob(rect5), plot_grid(c)+draw_grob(rect2)+
                      draw_line(
                        x = c(.41, .44),
                        y = c(-0.05, 0.0),
                        color = "red", size = .5,
                        arrow = arrow(length = unit(0.03, "npc"))
                      ), pdo3specific2+ggtitle(''), nrow=1, rel_widths = c(4,4,1.6)),
          plot_grid(f+draw_grob(rect6), plot_grid(e)+draw_grob(rect3)+
                      draw_line(
                        x = c(.31, .33),
                        y = c(-0.05, 0.0),
                        color = "red", size = .5,
                        arrow = arrow(length = unit(0.03, "npc"))
                      ), pdo6specific2+ggtitle(''), nrow=1, rel_widths = c(4,4,1.6)),
          g, nrow=4, labels = c('PDO2', 'PDO3', 'PDO6'), 
          rel_heights = c(4,4,4,1))
dev.off()

