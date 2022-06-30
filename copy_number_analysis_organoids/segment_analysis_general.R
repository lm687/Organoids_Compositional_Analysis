pcawg_CN_features = readRDS("data/pcawg_CN_features.rds")
tcga_CN_features = readRDS("data/tcga_CN_features.rds")

BriTROC_absolute_copynumber = readRDS("../../cnsignatures/manuscript_Rmarkdown/data/BriTROC_absolute_copynumber.rds")
BriTROC2_CN_features = readRDS("data/6_TCGA_Signatures_on_BRITROC/0_BRITROC_absolute_CN.rds")

organoids_absolute_copynumber = readRDS("data/organoid_absolute_CN.rds")
sampleNames(organoids_absolute_copynumber) = names_orgs$`new name`[match(gsub("org", "", sampleNames(organoids_absolute_copynumber)), names_orgs$`old name`)]
organoids_CN_features = extractCopynumberFeatures(organoids_absolute_copynumber)

BriTROC_CN_features = readRDS("data/BriTROC_CN_features.rds")