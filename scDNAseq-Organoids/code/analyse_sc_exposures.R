## read

rm(list=ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

source("../../copy_number_analysis_organoids/helper_functions.R")
library(gridExtra)
library(dplyr)

orgs = c('23868org', '119148orgb', '118976org')
# orgs = c('23868org', '118976org')

exposures = lapply(paste0("../data/signature_object_sc_cache_",  orgs, ".RDS"), readRDS)
absCN = lapply(paste0("../data/absCN_clean_",  orgs, ".RDS"), readRDS)
names(exposures) = names(absCN)= orgs

pdf("../plots/sc_exposures.pdf", width = 15)
do.call('grid.arrange', list(grobs=lapply(orgs, function(i){
  createBarplot(t(exposures[[i]]$sigs), remove_labels = T)+ggtitle(i)+theme(legend.position = "bottom")
}), nrow=1))
dev.off()


## plotting the 10X data
## checking if the samples without S3 look weird

dim(exposures$`23868org`$sigs)
dim(absCN$`23868org`)

sum(exposures$`23868org`$sigs['s3',] == 0)
sum(exposures$`119148orgb`$sigs['s3',] == 0)
sum(exposures$`118976org`$sigs['s3',] == 0)

mat_breaks <- c(-Inf, 0:6, seq(8, 14, by=2)) #seq(min(org_clean, na.rm = T), max(org_clean, na.rm = T), length.out = 10)

for(org in c('118976org', '23868org', '119148orgb')){
  org_clean_it = absCN[[org]]
  annotation_s3 = data.frame(zero_s3=as.factor(exposures[[org]]$sigs['s3',] == 0))
  rownames(org_clean_it) = colnames(exposures[[org]]$sigs)
  pdf(paste0("../plots/basicplot_s3_", org, ".pdf"), width=15)
  (pheatmap::pheatmap(org_clean_it, cluster_rows = FALSE, cluster_cols = FALSE, show_colnames = FALSE,
                     color             = c("#2670af", "#8daecf", "#eaeaea", "#ffd2b6", "#f5b788", "#f39d5f", "#f17e37", "#ee1200", "#b60000", "#7a0e04", "#401004", "#000000"),
                     breaks            = mat_breaks,
                     annotation_row = annotation_s3,
                     annotation_legend=T, show_rownames = FALSE))
  dev.off()
  
}

