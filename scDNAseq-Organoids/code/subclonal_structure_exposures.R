##----------------------------------------------------------------------------
rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("../../../../other_repos/britroc-cnsignatures-bfb69cd72c50/main_functions.R")
source("../../../../other_repos/britroc-cnsignatures-bfb69cd72c50/helper_functions.R")
source("../../copy_number_analysis_organoids/helper_functions.R")

##----------------------------------------------------------------------------

##----------------------------------------------------------------------------
get_exposures <- function(segs){
  features <- extractCopynumberFeatures(segs)
  # features <- extractCopynumberFeatures(data.frame(segs))
  # features_rounded <- extractCopynumberFeatures(segs_rounded)
  SxC <- generateSampleByComponentMatrix(CN_features = features)
  # SxC_rounded <- generateSampleByComponentMatrix(CN_features = features_rounded)
  sigs_ascites <- t(quantifySignatures(SxC, feat_sig_mat))
  # sigs_ascites_rounded <- t(quantifySignatures(SxC_rounded, feat_sig_mat))
}
##----------------------------------------------------------------------------

#-------------------------------------------------------------------------
component_parameters = readRDS("../../../britroc-cnsignatures-bfb69cd72c50/data/component_parameters.rds")
feat_sig_mat = readRDS("../../../britroc-cnsignatures-bfb69cd72c50/data/feat_sig_mat.rds")
sig_data = readRDS("../../copy_number_analysis_organoids/data/sig_data_unorm.RDS")
sig_data = cbind(sweep(sig_data[,1:7], 1, rowSums(sig_data[,1:7]), '/'),
                 sig_data[,8:ncol(sig_data)])
sig_data <- as.matrix(sig_data[,1:7])
#-------------------------------------------------------------------------

##----------------------------------------------------------------------------
get_exposures <- function(segs){
  features <- extractCopynumberFeatures(segs)
  # features <- extractCopynumberFeatures(data.frame(segs))
  # features_rounded <- extractCopynumberFeatures(segs_rounded)
  SxC <- generateSampleByComponentMatrix(CN_features = features, all_components = component_parameters)
  # SxC_rounded <- generateSampleByComponentMatrix(CN_features = features_rounded)
  sigs_ascites <- t(quantifySignatures(SxC, feat_sig_mat))
  # sigs_ascites_rounded <- t(quantifySignatures(SxC_rounded, feat_sig_mat))
}
##----------------------------------------------------------------------------

##----------------------------------------------------------------------------

### Using unified segments (common bins in all clades)
#### reading unified segments
# PDO2_unified <- readRDS("../../copy_number_analysis_organoids/data/clade_analysis_signatures/unifyCladesSegments_PDO2.rds") ## DO NOT USE
# PDO3_unified <- readRDS("../../copy_number_analysis_organoids/data/clade_analysis_signatures/unifyCladesSegments_PDO3.rds") ## DO NOT USE
# PDO6_unified <- readRDS("../../copy_number_analysis_organoids/data/clade_analysis_signatures/unifyCladesSegments_PDO6.rds") ## DO NOT USE

PDO2_unified <- readRDS("../../copy_number_analysis_organoids/data/clade_analysis_signatures/sudobulk_SegTabClades_PDO2.rds")
PDO3_unified <- readRDS("../../copy_number_analysis_organoids/data/clade_analysis_signatures/sudobulk_SegTabClades_PDO3.rds")
PDO6_unified <- readRDS("../../copy_number_analysis_organoids/data/clade_analysis_signatures/sudobulk_SegTabClades_PDO6.rds")

## PDO2

rename_segval <- function(i){
  colnames(i)[4] <- 'segVal'
  i
}

exposures_orgs <- readRDS("../../copy_number_analysis_organoids/robjects/exposures.RDS")

# PDO2 <- get_exposures(list(CladeA=rename_segval(data.frame(PDO2_unified[,c(1:4)])),
#                                        CladeB=rename_segval(data.frame(PDO2_unified[,c(1:3, 5)]))))
# PDO3 <- get_exposures(list(CladeA=rename_segval(data.frame(PDO3_unified[,c(1:4)])),
#                            CladeB=rename_segval(data.frame(PDO3_unified[,c(1:3, 5)])),
#                            CladeC=rename_segval(data.frame(PDO3_unified[,c(1:3, 6)]))))
# PDO6 <- get_exposures(list(CladeA=rename_segval(data.frame(PDO6_unified[,c(1:4)])),
#                            CladeB=rename_segval(data.frame(PDO6_unified[,c(1:3, 5)])),
#                            CladeC=rename_segval(data.frame(PDO6_unified[,c(1:3, 6)]))))

PDO2 <- get_exposures(PDO2_unified)
PDO3 <- get_exposures(PDO3_unified)
PDO6 <- get_exposures(PDO3_unified)

PDO2
PDO3
PDO6

createBarplot(rbind(PDO2=exposures_orgs['PDO2',], PDO2))
ggsave("../plots/PDO2_subclonal_exposures.pdf", height = 3)

createBarplot(rbind(PDO3=exposures_orgs['PDO3',], PDO3))
ggsave("../plots/PDO3_subclonal_exposures.pdf", height = 3)

createBarplot(rbind(PDO6=exposures_orgs['PDO6',], PDO6))
ggsave("../plots/PDO6_subclonal_exposures.pdf", height = 3)

##----------------------------------------------------------------------------

##----------------------------------------------------------------------------
table_absCNscDNA_PDO2_cladeA_1 <- readRDS("../robjects/table_absCNscDNA_PDO2_cladeA_1.RDS")
table_absCNscDNA_PDO2_cladeB_4 <- readRDS("../robjects/table_absCNscDNA_PDO2_cladeB_4.RDS")
table_absCNscDNA_PDO3_cladeA_3 <- readRDS("../robjects/table_absCNscDNA_PDO3_cladeA_3.RDS")
table_absCNscDNA_PDO3_cladeB_4 <- readRDS("../robjects/table_absCNscDNA_PDO3_cladeB_4.RDS")
table_absCNscDNA_PDO3_cladeC_6 <- readRDS("../robjects/table_absCNscDNA_PDO3_cladeC_6.RDS")
table_absCNscDNA_PDO6_cladeA_1 <- readRDS("../robjects/table_absCNscDNA_PDO6_cladeA_1.RDS")
table_absCNscDNA_PDO6_cladeB_2 <- readRDS("../robjects/table_absCNscDNA_PDO6_cladeB_2.RDS")
table_absCNscDNA_PDO6_cladeC_3 <- readRDS("../robjects/table_absCNscDNA_PDO6_cladeC_3.RDS")
##----------------------------------------------------------------------------

##----------------------------------------------------------------------------
### Using an alternative to unified segments

###' - Using 30kb bins
###' - the segments don't necessarily have to have the same length

## end and start positions are repeated. subtract 1 from end to avoid disjoint segments of width=1
end(table_absCNscDNA_PDO2_cladeA_1) <- end(table_absCNscDNA_PDO2_cladeA_1)-1
disjoin_PDO2_cladeA <- GenomicRanges::disjoin(table_absCNscDNA_PDO2_cladeA_1, with.revmap=T)
disjoin_PDO2_cladeA

min(table(width(disjoin_PDO2_cladeA)))
table_absCNscDNA_PDO2_cladeA_1[which(width(disjoin_PDO2_cladeA) == 1),]

## get the average of CN for each of the bins
disjoin_PDO2_cladeA$segval = sapply(disjoin_PDO2_cladeA$revmap, function(i){
  sum(table_absCNscDNA_PDO2_cladeA_1[i]$segVal)/length(i)
})
table(width(disjoin_PDO2_cladeA))

## outstanding problem: merge adjacent segments of very similar CN
### for each chromosome arm, get the changepoints
.pelt <- changepoint::PELT(matrix(disjoin_PDO2_cladeA[1:100,]$segval))
plot(disjoin_PDO2_cladeA[1:100,]$segval)

library(QDNAseq)
bins <- getBinAnnotations(binSize=30)

## outstanding problem: the minimum width is 20, i.e. these are 20kb bins
min(width(disjoin_PDO2_cladeA))
new("QDNAseqReadCounts")
?QDNAseqReadCounts
binReadCounts(disjoin_PDO2_cladeA)
binReadCounts(bins = annotatedDataFrameFrom(as(data.frame(start=start(disjoin_PDO2_cladeA),
      end=end(disjoin_PDO2_cladeA),
      segval=disjoin_PDO2_cladeA$segval), 'matrix'), byrow=T))

example_binanno <- getBinAnnotations(15) ## example
mat_bins <- as(data.frame(start=start(disjoin_PDO2_cladeA),
                          end=end(disjoin_PDO2_cladeA),
                          segval=disjoin_PDO2_cladeA$segval,
                          use=T), 'matrix')
rownames(mat_bins) <- paste0(seqnames(disjoin_PDO2_cladeA), ':', start(disjoin_PDO2_cladeA), '-', end(disjoin_PDO2_cladeA))
mat_binsanno <- annotatedDataFrameFrom(mat_bins, byrow = T)
mat_binsanno@dimLabels = c("rowNames",    "columnNames")
mat_binsanno
example_binanno

bins <- getBinAnnotations(15)
readCounts <- binReadCounts(bins)

binReadCounts(example_binanno)

## create directly a QDNAseqReadCounts object??
aaa <- new(Class = "QDNAseqReadCounts", bins=bins, phenodata=data.frame(100),
           counts=mat_bins)

# annotatedDataFrameFrom(as(data.frame(start=start(disjoin_PDO2_cladeA),
#                                     end=end(disjoin_PDO2_cladeA),
#                                     segval=disjoin_PDO2_cladeA$segval), 'matrix'), byrow=T)

##----------------------------------------------------------------------------

##----------------------------------------------------------------------------
