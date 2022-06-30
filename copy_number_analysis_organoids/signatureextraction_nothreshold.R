## re-extract CN signatures from the Nature Genetics paper

rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

source("../../../other_repos/britroc-cnsignatures-bfb69cd72c50/main_functions.R")
source("../../../other_repos/britroc-cnsignatures-bfb69cd72c50/helper_functions.R")
source("helper_functions.R")

##--------------------------------------------------------------------------------------------------------------##

component_parameters = readRDS("../../britroc-cnsignatures-bfb69cd72c50/data/component_parameters.rds")
feat_sig_mat = readRDS("../../britroc-cnsignatures-bfb69cd72c50/data/feat_sig_mat.rds")
sig_data = readRDS("data/sig_data_unorm.RDS")
sig_data = cbind(sweep(sig_data[,1:7], 1, rowSums(sig_data[,1:7]), '/'),
                 sig_data[,8:ncol(sig_data)])
sig_data <- as.matrix(sig_data[,1:7])
##--------------------------------------------------------------------------------------------------------------##
## get segment tables

pcawg_CN_features = readRDS("data/pcawg_CN_features.rds")
tcga_CN_features = readRDS("data/tcga_CN_features.rds")

BriTROC_absolute_copynumber = readRDS("../../cnsignatures/manuscript_Rmarkdown/data/BriTROC_absolute_copynumber.rds")
# BriTROC2_CN_features = readRDS("data/6_TCGA_Signatures_on_BRITROC/0_BRITROC_absolute_CN.rds")

organoids_absolute_copynumber = readRDS("data/organoid_absolute_CN.rds")
sampleNames(organoids_absolute_copynumber) = names_orgs$`new name`[match(gsub("org", "", sampleNames(organoids_absolute_copynumber)), names_orgs$`old name`)]
organoids_CN_features = extractCopynumberFeatures(organoids_absolute_copynumber)

BriTROC_CN_features = readRDS("data/BriTROC_CN_features.rds")

segtables_BriTROC_absolute_copynumber = lapply(sampleNames(BriTROC_absolute_copynumber), function(samplename) getSegTable(BriTROC_absolute_copynumber[,samplename]))
names(segtables_BriTROC_absolute_copynumber) = sampleNames(BriTROC_absolute_copynumber)
segtables_organoids_absolute_copynumber = lapply(sampleNames(organoids_absolute_copynumber), function(samplename) getSegTable(organoids_absolute_copynumber[,samplename]))
names(segtables_organoids_absolute_copynumber) = sampleNames(organoids_absolute_copynumber)

TCGA_absolute_copynumber = readRDS("data/combined.ascat.segments.filt.rds")

## we only want the ovarian ones
summary_ascat = read.table("data/summary.ascatTCGA.penalty70.txt", header = TRUE, stringsAsFactors = FALSE)
table(pass=summary_ascat$pass,
      dCIN=summary_ascat$dCIN)
tcga_samples = summary_ascat$name[which(summary_ascat$dCIN & (summary_ascat$cancer_type == 'OV'))]
## select the TCGA samples which are in the subset in which we are interested (the ones that passed QC and that are only OV)
TCGA_absolute_copynumber = TCGA_absolute_copynumber[TCGA_absolute_copynumber$sample %in% tcga_samples,]
segtables_TCGA_absolute_copynumber = lapply(as.character(sort(unique(TCGA_absolute_copynumber$sample))), function(samplename) TCGA_absolute_copynumber[TCGA_absolute_copynumber$sample == samplename,])
segtables_TCGA_absolute_copynumber = lapply(segtables_TCGA_absolute_copynumber, function(i) i[,1:(ncol(i)-1)]) ## removing sample name
segtables_TCGA_absolute_copynumber = lapply(segtables_TCGA_absolute_copynumber, function(i) data.frame(i))
names(segtables_TCGA_absolute_copynumber) = as.character(sort(unique(TCGA_absolute_copynumber$sample)))

##--------------------------------------------------------------------------------------------------------------##

## except for the original TCGA and PCAWG, we need to extract features
head(tcga_CN_features[[1]])
head(segtables_BriTROC_absolute_copynumber[[1]])

## get features
features_BriTROC <- extractCopynumberFeatures(segtables_BriTROC_absolute_copynumber)
features_TCGA <- extractCopynumberFeatures(segtables_TCGA_absolute_copynumber)
features_organoids <- extractCopynumberFeatures(segtables_organoids_absolute_copynumber)

##--------------------------------------------------------------------------------------------------------------##

## SxC matrices
length(unique(pcawg_CN_features[[1]]$ID))
SxC_pcawg <- generateSampleByComponentMatrix(CN_features = pcawg_CN_features)

##-------------

length(unique(tcga_CN_features[[1]]$ID))
split_oldTCGA <- split(unique(tcga_CN_features[[1]]$ID), f= factor(1:length(unique(tcga_CN_features[[1]]$ID)) %/% 120))
sapply(split_oldTCGA, length)
tcga_CN_featuresSPLIT <- lapply(1:length(split_oldTCGA), function(split_it){
  .x <- lapply(1:length(tcga_CN_features), function(i){
    tcga_CN_features[[i]][tcga_CN_features[[i]]$ID %in% split_oldTCGA[[split_it]],]
  })
  names(.x) <- names(tcga_CN_features)
  .x
})

SxC_oldTCGAPART1 <- generateSampleByComponentMatrix(CN_features = tcga_CN_featuresSPLIT[[1]])
SxC_oldTCGAPART2 <- generateSampleByComponentMatrix(CN_features = tcga_CN_featuresSPLIT[[2]])
SxC_oldTCGAPART3 <- generateSampleByComponentMatrix(CN_features = tcga_CN_featuresSPLIT[[3]])
SxC_oldTCGAPART4 <- generateSampleByComponentMatrix(CN_features = tcga_CN_featuresSPLIT[[4]])

SxC_oldTCGA <- rbind(SxC_oldTCGAPART1, SxC_oldTCGAPART2, SxC_oldTCGAPART3, SxC_oldTCGAPART4)
##-------------

SxC_BriTROC <- generateSampleByComponentMatrix(CN_features = features_BriTROC)

##-------------

SxC_organoids <- generateSampleByComponentMatrix(features_organoids)

##-------------
# split_TCGA <- split(unique(segtables_TCGA_absolute_copynumber[[1]]$ID), f= factor(1:length(unique(segtables_TCGA_absolute_copynumber[[1]]$ID)) %/% 120))
# tcga2_CN_featuresSPLIT <- lapply(1:length(split_oldTCGA), function(split_it){
#   .x <- lapply(1:length(segtables_TCGA_absolute_copynumber), function(i){
#     segtables_TCGA_absolute_copynumber[[i]][segtables_TCGA_absolute_copynumber[[i]]$ID %in% split_oldTCGA[[split_it]],]
#   })
#   names(.x) <- names(segtables_TCGA_absolute_copynumber)
#   .x
# })
SxC_TCGA <- generateSampleByComponentMatrix(features_TCGA)

dim(SxC_pcawg)

saveRDS(SxC_pcawg, file = "robjects/no_threshold_signature_extraction/SxC_pcawg.RDS")
saveRDS(SxC_oldTCGA, file = "robjects/no_threshold_signature_extraction/SxC_oldTCGA.RDS")
saveRDS(SxC_BriTROC, file = "robjects/no_threshold_signature_extraction/SxC_BriTROC.RDS")
saveRDS(SxC_organoids, file = "robjects/no_threshold_signature_extraction/SxC_organoids.RDS")
saveRDS(SxC_TCGA, file = "robjects/no_threshold_signature_extraction/SxC_TCGA.RDS")

##--------------------------------------------------------------------------------------------------------------##
## extract signatures
sigs_oldTCGA <- lapply(list(SxC_oldTCGAPART1,
                            SxC_oldTCGAPART2,
                            SxC_oldTCGAPART3,
                            SxC_oldTCGAPART4), function(i){
                              quantifySignaturesLM(sample_by_component = i, component_by_signature = feat_sig_mat,
                                                   sig_thresh=0)
                            })
sigs_oldTCGA <- t(do.call('cbind', sigs_oldTCGA))
sigs_pcawg <- t(quantifySignaturesLM(SxC_pcawg, feat_sig_mat, sig_thresh=0))
sigs_BriTROC <- t(quantifySignaturesLM(SxC_BriTROC, feat_sig_mat, sig_thresh=0))
sigs_organoids <- t(quantifySignaturesLM(SxC_organoids, feat_sig_mat, sig_thresh=0))
sigs_TCGA <- t(quantifySignaturesLM(SxC_TCGA, feat_sig_mat, sig_thresh=0))


saveRDS(sigs_pcawg, file = "robjects/no_threshold_signature_extraction/sigs_pcawg.RDS")
saveRDS(sigs_oldTCGA, file = "robjects/no_threshold_signature_extraction/sigs_oldTCGA.RDS")
saveRDS(sigs_BriTROC, file = "robjects/no_threshold_signature_extraction/sigs_BriTROC.RDS")
saveRDS(sigs_organoids, file = "robjects/no_threshold_signature_extraction/sigs_organoids.RDS")
saveRDS(sigs_TCGA, file = "robjects/no_threshold_signature_extraction/sigs_TCGA.RDS")

createBarplot(sigs_oldTCGA)

##--------------------------------------------------------------------------------------------------------------##
## compare to the other signatures

NAtozero <- function(m){m[is.na(m)] <- 0; m}
NAtoInf <- function(m){m[is.na(m)] <- Inf; m}
org<- as(readRDS("data/organoid_exposures.rds"), 'matrix')

df_TCGA_comparison <- data.frame(sigs_oldTCGA=as.vector(sigs_oldTCGA),
                                 sigs_TCGA=as.vector(NAtoInf(sigs_TCGA[match(substr(rownames(sigs_oldTCGA), 1, 12),
                                                                             rownames(sigs_TCGA)),])))
df_TCGA_comparison$sig = rep(1:7, each=nrow(sigs_oldTCGA))

ggplot(df_TCGA_comparison, aes(x=sigs_oldTCGA, y=sigs_TCGA, col=sig))+geom_point()+facet_wrap(.~sig)


ggplot(data.frame(pcawgNT=as.vector(sigs_pcawg),
                  natgen_pcawg=as.vector(sig_data[match(rownames(sigs_pcawg), rownames(sig_data)),]),
                  sig=rep(1:7, each=nrow(sigs_pcawg))), aes(x=natgen_pcawg, y=pcawgNT, col=sig))+
  geom_point()+facet_wrap(.~sig)+theme_bw()+geom_abline(slope = 1, intercept = 0, lty='dashed')


ggplot(data.frame(BriTROCNT=as.vector(sigs_BriTROC),
                  natgen_BriTROC=as.vector(sig_data[match(rownames(sigs_BriTROC), rownames(sig_data)),]),
                  sig=rep(1:7, each=nrow(sigs_BriTROC))), aes(x=natgen_BriTROC, y=BriTROCNT, col=sig))+
  geom_point()+facet_wrap(.~sig)+theme_bw()+geom_abline(slope = 1, intercept = 0, lty='dashed')


ggplot(data.frame(organoidsNT=as.vector(sigs_organoids),
                  natgen_organoids=as.vector(org[match(rownames(sigs_organoids), rownames(org)),]),
                  sig=rep(1:7, each=nrow(sigs_organoids))), aes(x=natgen_organoids, y=organoidsNT, col=sig))+
  geom_point()+facet_wrap(.~sig)+theme_bw()+geom_abline(slope = 1, intercept = 0, lty='dashed')


ggplot(data.frame(TCGANT=as.vector(sigs_TCGA),
                  natgen_TCGA=as.vector(sig_data[match(rownames(sigs_TCGA), substr(rownames(sig_data), 1, 12)),]),
                  sig=rep(1:7, each=nrow(sigs_TCGA))), aes(x=natgen_TCGA, y=TCGANT, col=sig))+
  geom_point()+facet_wrap(.~sig)+theme_bw()+geom_abline(slope = 1, intercept = 0, lty='dashed')
