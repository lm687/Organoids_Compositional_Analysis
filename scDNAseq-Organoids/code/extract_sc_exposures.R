## plotting the 10X data

rm(list=ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

source("../../../../../other_repos/britroc-cnsignatures-bfb69cd72c50/main_functions.R")
source("../../../../../other_repos/britroc-cnsignatures-bfb69cd72c50/helper_functions.R")

library(grid)
library(GenomicRanges)
library(parallel)

chromlen = readRDS("../../copy_number_analysis_organoids/data/chrlen.RDS")
rownames(chromlen) = chromlen[,1]
chroms_sorted = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", 
                  "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY")
chromlen = chromlen[chroms_sorted,]
## Create 30KB bins
bins30KB = do.call('rbind', apply(chromlen[,1:2], 1, function(i){.sq = seq(from = 1, to = as.numeric(i[2]), by = 30000); return(cbind(rep(i[1], length(.sq)),.sq,c(.sq[-1]-1,chromlen[chromlen$V1 == as.character(i[1]),2])))} ))
gr = as(paste0(bins30KB[,1], ':', bins30KB[,2], '-', bins30KB[,3]), "GRanges")
# table(width(gr))
# seqnames(gr)
seqlevels(gr)
seqlengths(gr) = chromlen[seqlevels(gr),2]

for(org in c('23868org', '119148orgb', '118976org')){
  org_tab = data.frame(read.table(paste0("../UID-", org, ".bed"), check.names = FALSE))
  org_tab$V2 = as.numeric(org_tab$V2)
  org_tab$V3 = as.numeric(org_tab$V3)
  org_tab$V4 = as.numeric(org_tab$V4)
  org_tab$V5 = as.numeric(org_tab$V5)
  
  # org_tab = apply(org_tab, 2, as.numeric)
  org_clean = org_tab

  ## remove segments which don't have values
  org_clean = org_clean[,!(rowSums(apply(org_clean, 1, is.na)) == nrow(org_clean))]

  ## find each sample in 
  segs_sc = cbind.data.frame(chrom=org_clean[,1], start=org_clean[,2], end=org_clean[,3], segval=org_clean[,5], sample=paste0('sample', org_clean[,4]))
  samples_name = unique(segs_sc$sample)
  segs_sc_list = lapply(unique(segs_sc$sample), function(i){
    .x = segs_sc[segs_sc$sample == i,]
    # paste0(.x[,1], ':', .x[,2], '-', .x[,3])
    data.frame(seqnames=paste0('chr', .x[,'chrom']), start=.x[,'start'], end=.x[,'end'], segVal=.x[,'segval'])
  })
  names(segs_sc_list) = samples_name

  grs_samples = mclapply(samples_name, function(idx_sample){
    cat('Sample:', idx_sample, '\n')
    gr_org = as(segs_sc_list[[idx_sample]], "GRanges")
    full_GR <- c(gr,gr_org)
    table(width(gr))
    table(width(gr_org)) ## some huge segments
    gr_org[order(width(gr_org), decreasing = T)[1:10]]
    disjoint_gr <- GenomicRanges::disjoin(full_GR, with.revmap=TRUE, ignore.strand=TRUE)
    seqlengths(disjoint_gr) <- chromlen[seqlevels(disjoint_gr),2]
    disjoint_gr = disjoint_gr[width(disjoint_gr)>1]
    
    disjoint_gr_revmap_first = sapply(disjoint_gr$revmap, function(i) i[1])
    disjoint_gr_revmap_second = sapply(disjoint_gr$revmap, function(i) i[2])
    
    idx_gr = disjoint_gr_revmap_second-length(gr)
    disjoint_gr$CN_val = gr_org$segVal[idx_gr]
    disjoint_gr$CN_val[is.na(disjoint_gr$CN_val)] = 2
    disjoint_gr$revmap = disjoint_gr_revmap_first

    CN_RleList = mcolAsRleList(disjoint_gr, "CN_val")
    averageCN = binnedAverage(bins = gr, numvar = CN_RleList, varname = "CN_bin_averaged")
    
    return(averageCN=averageCN)
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
