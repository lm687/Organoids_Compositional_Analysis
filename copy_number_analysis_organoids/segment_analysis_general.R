rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

source("../../../other_repos/britroc-cnsignatures-bfb69cd72c50/main_functions.R")
source("../../../other_repos/britroc-cnsignatures-bfb69cd72c50/helper_functions.R")

pcawg_CN_features = readRDS("data/pcawg_CN_features.rds")
tcga_CN_features = readRDS("data/tcga_CN_features.rds")

BriTROC_absolute_copynumber = readRDS("../../cnsignatures/manuscript_Rmarkdown/data/BriTROC_absolute_copynumber.rds")
BriTROC2_CN_features = readRDS("data/6_TCGA_Signatures_on_BRITROC/0_BRITROC_absolute_CN.rds")

organoids_absolute_copynumber = readRDS("data/organoid_absolute_CN.rds")
# sampleNames(organoids_absolute_copynumber) = names_orgs$`new name`[match(gsub("org", "", sampleNames(organoids_absolute_copynumber)), names_orgs$`old name`)]
organoids_CN_features = extractCopynumberFeatures(organoids_absolute_copynumber)

BriTROC_CN_features = readRDS("data/BriTROC_CN_features.rds")

pcawg_CN_features
tcga_CN_features
BriTROC_absolute_copynumber
BriTROC2_CN_features


sum_segs_per_chrom <- BriTROC2_CN_features %>% group_by(sample) %>% summarise(num_segs=sapply(sort(unique(BriTROC2_CN_features$chromosome)),
              function(chrom) c(sum(chromosome == chrom))))
samples_britroc2 <- sum_segs_per_chrom$sample
sum_segs_per_chrom <- matrix(sum_segs_per_chrom$num_segs, nrow=length(unique(BriTROC2_CN_features$chromosome)))
rownames(sum_segs_per_chrom) <- sort(unique(BriTROC2_CN_features$chromosome))
colnames(sum_segs_per_chrom) <- unique(samples_britroc2)
# sum_segs_per_chrom <- dcast(sum_segs_per_chrom, num_segs~sample)
# rownames(sum_segs_per_chrom) <- sum_segs_per_chrom$sample; sum_segs_per_chrom$sample <- NULL
sum_segs_per_chrom_pca <- prcomp(t(sum_segs_per_chrom))
plot(sum_segs_per_chrom_pca$x[,1:2]) ## no clustering at all
