## plotting the 10X data

rm(list=ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

source("../../../../../other_repos/britroc-cnsignatures-bfb69cd72c50/main_functions.R")
source("../../../../../other_repos/britroc-cnsignatures-bfb69cd72c50/helper_functions.R")

library(grid)

for(org in c('23868org', '119148orgb', '118976org')){
  org_tab = read.csv(paste0("../UID-", org, ".csv"), check.names = FALSE)
  org_tab = apply(org_tab, 2, as.numeric)
  keep_idx = org_tab[,1:4][,'num_noisy'] == 0
  org_clean = org_tab[keep_idx,-(1:4)]
  org_tab = org_tab[keep_idx,]
  rownames(org_clean) = paste0('scSample_', org_tab[,1])
  
  ## remove segments which don't have values
  org_clean = org_clean[,!(rowSums(apply(org_clean, 1, is.na)) == nrow(org_clean))]
  
  .colnames_pos = colnames(org_clean)
  .colnames_pos2 = lapply(.colnames_pos, function(i) strsplit(i, split = ":")[[1]])
  chrom = sapply(.colnames_pos2, function(i) i[1])
  .colnames_pos3 = lapply(.colnames_pos2, function(i) strsplit(i[2], '-')[[1]])
  starts = sapply(.colnames_pos3, function(i) as.numeric(gsub(",", "", i[1])) )
  ends = sapply(.colnames_pos3, function(i) as.numeric(gsub(",", "", i[2])) )
  
  gr = as(paste0('chr', chrom, ':', starts, '-', ends), "GRanges")
  ## now do it by 
  
  segs_sc = lapply(rownames(org_clean), function(i){
    a = cbind.data.frame(chromosome = chrom, start=starts, end=ends, segVal = org_clean[i,], sample=i)
    a = a[!is.na(a$segVal),]
  })
  names(segs_sc) = rownames(org_clean)
  sapply(segs_sc, function(k) table(is.na(k$segVal)))
  # $`JBLAB-4998`
  # chromosome     start       end      segVal     sample
  # 26837          1    840001  16110000  3.12972788 JBLAB-4998
  # 26838          1  16110001  16770000  4.27279458 JBLAB-4998
  # 26839          1  16770001  54660000  3.07608674 JBLAB-4998
  # 26840          1  54660001  56220000  4.19125637 JBLAB-4998
  # 26841          1  56220001 110130000  2.96476171 JBLAB-4998
  # 26842          1 110130001 111180000  4.07866743 JBLAB-4998
  # 26843          1 111180001 119880000  3.04572234 JBLAB-4998
  
  cnfeatures = extractCopynumberFeatures(segs_sc)
  SxCmatrix = generateSampleByComponentMatrix(cnfeatures)
  sigs = quantifySignatures(SxCmatrix)
  saveRDS(list(cnfeatures=cnfeatures, SxCmatrix=SxCmatrix, sigs=sigs), file = paste0("../data/signature_object_sc_cache_", org, ".RDS"))
}
