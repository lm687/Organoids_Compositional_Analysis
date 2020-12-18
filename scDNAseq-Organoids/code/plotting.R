## plotting the 10X data

rm(list=ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(viridis)
library(grid)

offset_x_vec = c('118976org'=45, '23868org'=45, '119148orgb' = 45)
scaling_x_vec = c('118976org'=0.88, '23868org'=0.88, '119148orgb' = 0.88)
position_chrom_labels_vec = c('118976org'=0.93, '23868org'=0.93, '119148orgb' = 0.95)
cellheight_vec = c('118976org'=2.5, '23868org'=4.8, '119148orgb' = 1.1)
cellwidth_vec = c('118976org'=1.6, '23868org'=1.6, '119148orgb' = 1.6)

for(org in c('118976org', '23868org', '119148orgb')){
  org_tab = read.csv(paste0("../UID-", org, ".csv"))
  org_tab = apply(org_tab, 2, as.numeric)
  
  ## printing the cells that have been removed
  print(org)
  print(org_tab[org_tab[,'num_noisy'] != 0,'node_id'])
  
  org_clean = org_tab[org_tab[,1:4][,'num_noisy'] == 0,-(1:4)]
  mat_breaks <- c(-Inf, 0:6, seq(8, 14, by=2)) #seq(min(org_clean, na.rm = T), max(org_clean, na.rm = T), length.out = 10)
  labels_chroms0 = sapply(colnames(org_clean), function(i) strsplit(i, '[.]')[[1]][1])
  labels_chroms = unique(labels_chroms0)
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
  
  offset_x = offset_x_vec[org]
  scaling_x = scaling_x_vec[org]
  labels_chroms[1:22] = substr(labels_chroms[1:22], 2, 1000)
  labels_chroms = cbind.data.frame(labels_chroms,
                        idx_first=sapply(unique(sapply(colnames(org_clean), function(i) strsplit(i, '[.]')[[1]][1])), function(i) scaling_x*(offset_x + mean(min(which(sapply(colnames(org_clean), function(i) strsplit(i, '[.]')[[1]][1]) == i)),
                                                                                                                                               max(which(sapply(colnames(org_clean), function(i) strsplit(i, '[.]')[[1]][1]) == i))))))
  labels_chroms$idx_first_norm =labels_chroms$idx_first/(ncol(org_clean)-4)
  annotation_chroms = data.frame(row.names = colnames(org_clean), chrom=clean_chrom(labels_chroms0))
  
  pdf(paste0("../plots/basicplot_", org, ".pdf"), width=15)
  pheatmap::pheatmap(org_clean, cluster_rows = FALSE, cluster_cols = FALSE, show_colnames = FALSE,
                     color             = c("#2670af", "#8daecf", "#eaeaea", "#ffd2b6", "#f5b788", "#f39d5f", "#f17e37", "#ee1200", "#b60000", "#7a0e04", "#401004", "#000000"),
                     breaks            = mat_breaks, cellheight=cellheight_vec[org],cellwidth=cellwidth_vec[org], annotation_col = annotation_chroms, annotation_legend=F)
  grid::grid.text(labels_chroms$labels_chroms,
                  x=labels_chroms$idx_first_norm,y=position_chrom_labels_vec[org], gp=gpar(fontsize=10))
  dev.off()
}