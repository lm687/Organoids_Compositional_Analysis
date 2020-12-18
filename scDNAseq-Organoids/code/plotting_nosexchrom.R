## plotting the 10X data

rm(list=ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(viridis)
library(grid)
library(pheatmap)
library(gridExtra)

offset_x_vec = c('118976org'=70, '23868org'=65, '119148orgb' = 65)
scaling_x_vec = c('118976org'=0.815, '23868org'=0.82, '119148orgb' = 0.82)
position_chrom_labels_vec = c('118976org'=0.93, '23868org'=0.93, '119148orgb' = 0.95)
cellheight_vec = c('118976org'=2.5, '23868org'=4.8, '119148orgb' = 1.1)
cellwidth_vec = c('118976org'=1.6, '23868org'=1.6, '119148orgb' = 1.6)

renaming = c('118976org'='PDO3', '23868org'='PDO2', '119148orgb'='PDO6')

vec_colours = c("#2670af", "#8daecf", "#eaeaea", "#ffd2b6", "#f5b788", "#f39d5f", "#f17e37", "#ee1200", "#b60000", "#7a0e04", "#401004", "#000000")

x = list()
x_basic = list()
idx = 0
for(org in c('118976org', '23868org', '119148orgb')){
  idx = idx+1
  org_tab = read.csv(paste0("../UID-", org, ".csv"))
  org_tab = apply(org_tab, 2, as.numeric)
  
  ## printing the cells that have been removed
  print(org)
  print(org_tab[org_tab[,'num_noisy'] != 0,'node_id'])
  
  org_clean = org_tab[org_tab[,1:4][,'num_noisy'] == 0,-(1:4)]
  mat_breaks <- c(-Inf, 0:6, seq(8, 14, by=2)) #seq(min(org_clean, na.rm = T), max(org_clean, na.rm = T), length.out = 10)
  labels_chroms0 = sapply(colnames(org_clean), function(i) strsplit(i, '[.]')[[1]][1])
  keep = !(labels_chroms0 %in% c('X', 'Y'))
  org_clean = org_clean[,keep]
  labels_chroms0 = labels_chroms0[keep]
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
  
  x[[idx]] = pheatmap(org_clean, cluster_rows = FALSE, cluster_cols = FALSE, show_colnames = FALSE,
                         color             = c("#2670af", "#8daecf", "#eaeaea", "#ffd2b6", "#f5b788", "#f39d5f", "#f17e37", "#ee1200", "#b60000", "#7a0e04", "#401004", "#000000"),
                         breaks            = mat_breaks, cellheight=cellheight_vec[org],cellwidth=cellwidth_vec[org], annotation_col = annotation_chroms, annotation_legend=F)
  if(org == '23868org'){
    legend_bool=T
  }else{
    legend_bool =F
  }
  x_basic[[idx]] = pheatmap(org_clean, cluster_rows = FALSE, cluster_cols = FALSE, show_colnames = FALSE,
                      color             = vec_colours,
                      breaks            = mat_breaks, annotation_col = annotation_chroms, annotation_legend=F,
                      legend = legend_bool, cellwidth=cellwidth_vec[org], main=renaming[org], annotation_names_col = FALSE)
  pdf(paste0("../plots/basicplot_", org, "_nosexchrom.pdf"), width=15)
  print(x[[idx]])
  grid::grid.text(labels_chroms$labels_chroms,
                  x=labels_chroms$idx_first_norm,y=position_chrom_labels_vec[org], gp=gpar(fontsize=10))
  dev.off()
}

sapply(x_basic, function(i) i[[4]])

grab_legend = x_basic[[2]][[4]]
x_basic_mod = x_basic
x_basic_mod[[2]][[4]] = x_basic_mod[[2]][[4]][,1:4]
grab_legend_mod = grab_legend[,5]


scaling_x_all = 0.81
offset_x_all = 42
labels_chroms = unique(labels_chroms0)
labels_chroms[1:22] = substr(labels_chroms[1:22], 2, 1000)
labels_chroms = cbind.data.frame(labels_chroms,
           idx_first=sapply(unique(sapply(colnames(org_clean), function(i) strsplit(i, '[.]')[[1]][1])),
                  function(i) scaling_x_all*(offset_x_all + mean(min(which(sapply(colnames(org_clean),
                      function(i) strsplit(i, '[.]')[[1]][1]) == i)),
                               max(which(sapply(colnames(org_clean), function(i) strsplit(i, '[.]')[[1]][1]) == i))))))
labels_chroms$idx_first_norm =labels_chroms$idx_first/(ncol(org_clean)-4)

pdf(paste0("../plots/basicplot_all_plots_nosexchrom.pdf"), width=15)
grid.arrange(grobs=list(x_basic_mod[[1]][[4]], rectGrob(gp=gpar(col=NA)),
                        x_basic_mod[[2]][[4]],  grab_legend_mod,
                        x_basic_mod[[3]][[4]], rectGrob(gp=gpar(col=NA))), ncol=2, widths = c(8/9, 1/9))
grid::grid.text(labels_chroms$labels_chroms,
                x=labels_chroms$idx_first_norm,y=0.9505, gp=gpar(fontsize=10))
grid::grid.text(labels_chroms$labels_chroms,
                x=labels_chroms$idx_first_norm,y=0.618, gp=gpar(fontsize=10))
grid::grid.text(labels_chroms$labels_chroms,
                x=labels_chroms$idx_first_norm,y=0.286, gp=gpar(fontsize=10))
dev.off()





