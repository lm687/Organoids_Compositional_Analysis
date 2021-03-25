rm(list = ls())

library(reshape2)
library(ggplot2)
library(ggdendro)
library(pheatmap)
library(readxl)
library(gridExtra)
library(ggplotify)
library(cowplot)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

renaming <- readxl::read_excel("../../RNASeq_DE_resistant_sensitive/files/PDOnameProperSample_sWGS_RNAseq.xlsx")

organoid_list = c('118976org', '119148orgb', '23868org')

absCN = lapply(paste0("../data/absCN_clean_", organoid_list, ".RDS"), readRDS)
names(absCN) = organoid_list
names(absCN) = gsub("orgb", "org", names(absCN))
names(absCN) = renaming$PDO[match(names(absCN), renaming$ID)]
saveRDS(absCN, "../robjects/fig2_absCN.RDS")
organoid_list = renaming$PDO[match(gsub("orgb", "org", organoid_list), renaming$ID)]


# image(absCN[[1]])
# 
# plot(hclust(dist(absCN[[org_it]])))

# org_it = '23868org'

## remove the outliers
outliers = list()
outliers$`PDO3` = c(12, 10, 11, 7, 5, 6, 4, 3, 1, 2, 8, 9, 158)
outliers$`PDO6` = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 23, 377, 378, 379, 380, 381)
outliers$`PDO2` = c(1,2,3, 4, 5, 6, 7, 8)

# plot(hclust(dist(absCN[[org_it]][-outliers[[org_it]],])))

clean_chrom = function(i){
  sapply(i, function(j){
    if(j %in% c('X', 'Y')){
      j
    }else{
      ## autosomal
      substr(j, 2, 1000)
    }
  })
}

for(org_it in organoid_list){
  absCN[[org_it]] = absCN[[org_it]][-outliers[[org_it]],]
  
  absCN[[org_it]] = absCN[[org_it]][,rowSums(apply(absCN[[org_it]], 1, is.na)) == 0]
  
  ## select only top variable areas of the genome
  # absCN[[org_it]] = absCN[[org_it]][,order(apply(absCN[[org_it]] , 2, var), decreasing = T)[1:1000]]
  ## doesn't work too well
  
  absCN[[org_it]][absCN[[org_it]] > 14] = 14
  
  mat_breaks <- c(-Inf, 0:6, seq(8, 14, by=2)) #seq(min(org_clean, na.rm = T), max(org_clean, na.rm = T), length.out = 10)
  col_list <- c("#2670af", "#8daecf", "#eaeaea", "#ffd2b6", "#f5b788", "#f39d5f", "#f17e37", "#ee1200", "#b60000", "#7a0e04", "#401004", "#000000")
  annotation_chroms = data.frame(row.names = colnames(absCN[[org_it]]),
                                 chrom=clean_chrom(gsub("\\..*", "", colnames(absCN[[org_it]]))))
  
  ph = pheatmap::pheatmap(absCN[[org_it]], show_colnames = FALSE, show_rownames = FALSE,
                     color             = col_list,
                     breaks            = mat_breaks, cluster_cols = FALSE,
                     annotation_col = annotation_chroms, annotation_legend=F, main=org_it)
  ph
  saveRDS(list(mat_breaks=mat_breaks,col_list=col_list), paste0("../robjects/fig2_colours.RDS"))
  saveRDS(absCN[[org_it]], paste0("../robjects/fig2_subclonal_hclust", org_it, "_2.RDS"))
  saveRDS(ph, paste0("../robjects/fig2_subclonal_hclust", org_it, ".RDS"))
  # ph$gtable$grobs[[1]]$gp <- gpar(lwd = 5)
  # ph$gtable$grobs[[2]]$gp <- gpar(col = 'blue')
  
  plot(as.dendrogram(ph$tree_row))
  pdf(paste0("../plots/subclonal_hclust", org_it, ".pdf"))
  # quartz()
  # print(grid.arrange(ggdendrogram(as.dendrogram(ph$tree_row), no.margin = TRUE)+
  #                theme(axis.line=element_blank(),axis.text.x=element_blank(),
  #                      axis.text.y=element_blank(),axis.ticks=element_blank(),
  #                      # axis.title.x=element_blank(),
  #                      axis.title.y=element_blank(),legend.position="none")+scale_x_continuous(expand = c(0,0)),
  #              ggplot(melt(t(absCN[[org_it]])), aes(y=Var1, x=factor(Var2, levels=ph$tree_row$order), fill=value))+geom_tile()+
  #   # scale_colour_steps(breaks=mat_breaks,value=col_list)
  #     scale_fill_gradientn(colours=col_list, breaks=mat_breaks)+
  #   theme(axis.line=element_blank(),axis.text.x=element_blank(),
  #         axis.text.y=element_blank(),axis.ticks=element_blank(),
  #         axis.title.x=element_blank(),
  #         axis.title.y=element_blank(),legend.position="none",
  #         panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
  #         panel.grid.minor=element_blank(),plot.background=element_blank()), nrow=2, top=org_it))
  
  ph
  # dendro <- ggdendrogram(as.dendrogram(ph$tree_row), no.margin = TRUE)+
  #                  theme(axis.line=element_blank(),axis.text.x=element_blank(),
  #                        axis.text.y=element_blank(),axis.ticks=element_blank(),
  #                        # axis.title.x=element_blank(),
  #                        axis.title.y=element_blank(),legend.position="none")+scale_x_continuous(expand = c(0,0))+
  #   coord_flip()
  # cowplot::plot_grid(ggplotify::as.grob(ph), dendro)
  # dev.off()
  
}

# quartz()
ggplot(melt(t(absCN[[org_it]])), aes(y=Var1, x=factor(Var2, levels=ph$tree_row$order), fill=value))+
  geom_tile()+
  # scale_colour_steps(breaks=mat_breaks,value=col_list)
  scale_fill_gradientn(colours=col_list, breaks=mat_breaks)
# 
# 
