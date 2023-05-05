#' Checking if there are normal cells in the scDNA data

rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

organoid_list = c('118976org', '119148orgb', '23868org')

absCN = lapply(paste0("../data/absCN_clean_", organoid_list, ".RDS"), readRDS)
names(absCN) = organoid_list
names(absCN) = gsub("orgb", "org", names(absCN))

rowmeans_all <- lapply(absCN, rowMeans, na.rm = T)
rowvar_all <- lapply(absCN, function(i) apply(i, 1, 'var', na.rm = T) )

par(mfrow=c(1,3))
sapply(rowmeans_all, function(i) plot(density(i)))
sapply(1:3, function(i) plot(rowmeans_all[[i]], rowvar_all[[i]]))
## no plot with cells at x=2 and variance=0

## cells closest to diploid (i.e. those we are interested in)
give_lowest_CN <- function(i) pheatmap::pheatmap(t(apply(absCN[[i]][order(rowmeans_all[[i]])[1:10],] < 2.2, 1, as.numeric)), cluster_cols = F)
## cells furthest from diploid
give_highest_CN <- function(i) pheatmap::pheatmap(t(apply(absCN[[i]][order(rowmeans_all[[i]],decreasing=T)[1:10],] < 2.2, 1, as.numeric)), cluster_cols = F)

## in red, diploid regions
give_lowest_CN(1) ## only one possible normal cell, but even that has CN around 2, but not even everywhere
give_lowest_CN(2) ## no normal cell
give_lowest_CN(3) ## no normal cell

# give_highest_CN(1)
