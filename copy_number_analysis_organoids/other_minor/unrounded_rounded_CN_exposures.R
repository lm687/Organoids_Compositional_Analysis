#-------------------------------------------------------------------------

rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(devtools)
library(QDNAseqmod)
source("../../../../other_repos/britroc-cnsignatures-bfb69cd72c50/main_functions.R")
source("../../../../other_repos/britroc-cnsignatures-bfb69cd72c50/helper_functions.R")
source("../../../../other_repos/britroc-1/code/models/helper/functions.R")
source("../helper_functions.R")
#-------------------------------------------------------------------------


#-------------------------------------------------------------------------
read_ICGC <- function(){
  ICGC_absolute_copynumber_AU = readRDS("../data/CN_Calls_ABSOLUTE_PCAWG/OV-AU.segments.raw.rds")
  ICGC_absolute_copynumber_US = readRDS("../data/CN_Calls_ABSOLUTE_PCAWG/OV-US.segments.raw.rds")
  ICGC_absolute_copynumber_AU = ICGC_absolute_copynumber_AU[,c('sample', 'chr', 'startpos', 'endpos', 'segVal')]
  ICGC_absolute_copynumber_US = ICGC_absolute_copynumber_US[,c('sample', 'chr', 'startpos', 'endpos', 'segVal')]
  
  segtables_ICGC_absolute_copynumber_AU = lapply(sort(unique(ICGC_absolute_copynumber_AU$sample)),
                                                 function(samplename)
                                                   ICGC_absolute_copynumber_AU[ICGC_absolute_copynumber_AU$sample == samplename,])
  segtables_ICGC_absolute_copynumber_AU = lapply(segtables_ICGC_absolute_copynumber_AU, function(i) { colnames(i)[colnames(i) == "chr"] = "chromosome";
  colnames(i)[colnames(i) == "endpos"] = "end";
  return(i) } )
  names(segtables_ICGC_absolute_copynumber_AU) = unique(ICGC_absolute_copynumber_AU$sample)
  
  segtables_ICGC_absolute_copynumber_US = lapply(sort(unique(ICGC_absolute_copynumber_US$sample)),
                                                 function(samplename) ICGC_absolute_copynumber_US[ICGC_absolute_copynumber_US$sample == samplename,])
  segtables_ICGC_absolute_copynumber_US = lapply(segtables_ICGC_absolute_copynumber_US, function(i) { colnames(i)[colnames(i) == "chr"] = "chromosome";
  colnames(i)[colnames(i) == "endpos"] = "end";
  return(i) } )
  names(segtables_ICGC_absolute_copynumber_US) = unique(ICGC_absolute_copynumber_US$sample)
  
  ## for ICGC, remove the samples row and put it in the rows
  segtables_ICGC_absolute_copynumber_US = lapply(segtables_ICGC_absolute_copynumber_US, function(i){
    rownames(i) = i$samples
    i = i[,-1]
    i})
  segtables_ICGC_absolute_copynumber_AU = lapply(segtables_ICGC_absolute_copynumber_AU, function(i){
    rownames(i) = i$samples
    i = i[,-1]
    i})
  return(list(segtables_ICGC_absolute_copynumber_US=segtables_ICGC_absolute_copynumber_US,
              segtables_ICGC_absolute_copynumber_AU=segtables_ICGC_absolute_copynumber_AU))
}

ICGC_absolute_copynumber <- read_ICGC()
BriTROC_absolute_copynumber = readRDS("../../../cnsignatures/manuscript_Rmarkdown/data/BriTROC_absolute_copynumber.rds")
BriTROC_absolute_copynumber_segtable <- getSegTable(BriTROC_absolute_copynumber)
table(BriTROC_absolute_copynumber_segtable$chromosome)
table(is.na(BriTROC_absolute_copynumber_segtable$chromosome)) ## the NAs are already in the original data??
organoids_absolute_copynumber = readRDS("../data/organoid_absolute_CN.rds")
organoids_absolute_copynumber_segtable <- getSegTable(organoids_absolute_copynumber)
table(is.na(organoids_absolute_copynumber_segtable$chromosome))
TCGA_absolute_copynumber = readRDS("../data/combined.ascat.segments.filt.rds")

segs_ascites <- readRDS("../data/20220511BH_ascites_absoluteCN_bestfit.rds")
segs_ascites_segtable <- getSegTable(segs_ascites)
table(is.na(segs_ascites_segtable$chromosome))
segs_ascites_segtable$sample <- NA#rep(colnames(segs_ascites@assayData$copynumber), each=nrow(segs_ascites@assayData$copynumber))


segs_rounded <- segs
segs_rounded@assayData <- list(copynumber=round(segs_rounded@assayData$copynumber),
                               segmented=round(segs_rounded@assayData$segmented))
# assign(segs_rounded@assayData$copynumber, round(segs_rounded@assayData$copynumber))

## BriTROC_absolute_copynumber ## ::::: QDNAseqCopyNumbers
## ICGC_absolute_copynumber ## ::::::: segment table
## organoids_absolute_copynumber ## ::::: QDNAseqCopyNumbers
## segs_ascites ## ::::: QDNAseqCopyNumbers
## TCGA_absolute_copynumber ## ::::::: segment table

get_rounded <- function(qdnaseq_obj){
  qdnaseq_obj@assayData <- list(copynumber=round(qdnaseq_obj@assayData$copynumber),
                                 segmented=round(qdnaseq_obj@assayData$segmented))
  qdnaseq_obj
}

# c(segs_ascites, organoids_absolute_copynumber)
# get_rounded(segs_ascites)
# segs <- rbind(segs_ascites_segtable, data.frame(TCGA_absolute_copynumber))
# segs_rounded <- segs
# segs_rounded$segVal <- round(as.numeric(segs_rounded$segVal))
#-------------------------------------------------------------------------

#-------------------------------------------------------------------------
component_parameters = readRDS("../../../britroc-cnsignatures-bfb69cd72c50/data/component_parameters.rds")
feat_sig_mat = readRDS("../../../britroc-cnsignatures-bfb69cd72c50/data/feat_sig_mat.rds")
sig_data = readRDS("../data/sig_data_unorm.RDS")
sig_data = cbind(sweep(sig_data[,1:7], 1, rowSums(sig_data[,1:7]), '/'),
                 sig_data[,8:ncol(sig_data)])
sig_data <- as.matrix(sig_data[,1:7])
#-------------------------------------------------------------------------

#-------------------------------------------------------------------------
get_exposures <- function(segs){
  features <- extractCopynumberFeatures(segs)
  # features <- extractCopynumberFeatures(data.frame(segs))
  # features_rounded <- extractCopynumberFeatures(segs_rounded)
  SxC <- generateSampleByComponentMatrix(CN_features = features)
  # SxC_rounded <- generateSampleByComponentMatrix(CN_features = features_rounded)
  sigs_ascites <- t(quantifySignatures(SxC, feat_sig_mat))
  # sigs_ascites_rounded <- t(quantifySignatures(SxC_rounded, feat_sig_mat))
}

ascites_exposures <- list(unrounded=get_exposures(segs_ascites),
                          rounded=get_exposures(get_rounded(segs_ascites)))
saveRDS(ascites_exposures, "../robjects/rounded_unrounded_ascites_exposures.RDS")
BriTROC_exposures <- list(unrounded=get_exposures(BriTROC_absolute_copynumber),
                          rounded=get_exposures(get_rounded(BriTROC_absolute_copynumber)))
saveRDS(BriTROC_exposures, "../robjects/rounded_unrounded_BriTROC_exposures.RDS")

#-------------------------------------------------------------------------

#-------------------------------------------------------------------------
compare_matrices(sigs_ascites, sigs_ascites_rounded, facets = T, ncol=7)
imputation_value <- 0.001
ggplot(melt(as(compositions::alr(sigs_ascites[,c(2:7, 1)]+imputation_value), 'matrix') - as(compositions::alr(sigs_ascites_rounded[,c(2:7, 1)]+imputation_value), 'matrix')),
       aes(x=paste0(Var2, '/s1'), y=value, group=interaction(Var2, Var1), col=Var2))+geom_boxplot()+theme_bw()
## s5 is bimodal: in some cases there are no changes, in others extreme cases
#-------------------------------------------------------------------------

#-------------------------------------------------------------------------
give_alr_subtraction <- function(exp){
  as(compositions::alr(exp$unrounded[,c(2:7, 1)]+imputation_value), 'matrix') - as(compositions::alr(exp$rounded[,c(2:7, 1)]+imputation_value), 'matrix')
}

alr_subtractions <- lapply(list(BriTROC=BriTROC_exposures, ascites=ascites_exposures), give_alr_subtraction)

alr_subtractions

ggplot(melt(alr_subtractions), aes(x=interaction(L1, Var2, paste0(Var2, '/s1')),#x=paste0(Var2, '/s1'),
                                   y=value, group=interaction(L1, Var2), col=Var2, shape=L1))+
  geom_violin()+geom_point(alpha=0.2, width = 0.1)+theme_bw()+ggtitle('Difference between rounded and unrounded CN')+
  labs(shape='Cohort', col='Signature', y='ALR wrt s1', x='LogR of signatures')+
  scale_fill_brewer(palette="Dark2")+
  annotate("text", x = 7.5, y = 10, label = "Higher in unrounded")+
  annotate("text", x = 7.5, y = -6, label = "Higher in rounded")+#+facet_wrap(.~Var2,scales = "free_x", nrow=1)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave("../../copy_number_analysis_organoids/figures/unrounded_rounded_copynumber_exposures.pdf", height = 5)
#-------------------------------------------------------------------------

pairs(alr_subtractions$BriTROC)
plot(density(alr_subtractions$BriTROC[,5]))
