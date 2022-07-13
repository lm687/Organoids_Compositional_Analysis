rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(devtools)
library(gridExtra)
# install('~/software/QDNAseqmod/') ## I have had to modify DESCRIPTION and added "mod" to its name
library(QDNAseqmod)
library(ggplot2)
library(reshape2)

# segs <- readRDS("data/20220511BH_ascites_absoluteCN_bestfit.rds")
## for one signature, get the updated version of the ascites, and then save them all in a file with all used optimal fits
segs <- readRDS("data/20220629BH_ascites_absoluteCN_bestfit.rds")
length(sampleNames(segs))
segs@phenoData@data$name

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

# pcawg_CN_features = readRDS("data/pcawg_CN_features.rds")
# tcga_CN_features = readRDS("data/tcga_CN_features.rds")
# BriTROC_CN_features = readRDS("data/BriTROC_CN_features.rds")

recompute_exposures <- T
if(recompute_exposures){
  features <- extractCopynumberFeatures(segs)
  SxC <- generateSampleByComponentMatrix(CN_features = features)
  # sigs_ascites <- t(quantifySignaturesLM(SxC, feat_sig_mat, sig_thresh=0))
  sigs_ascites <- t(quantifySignatures(SxC, feat_sig_mat))
  rownames(sigs_ascites) <- gsub("_", "-", rownames(sigs_ascites))
  
  saveRDS(sigs_ascites, "robjects/sigs_ascites.RDS")
}else{
  sigs_ascites <- readRDS("robjects/sigs_ascites.RDS")
}

## For one of the organoids (PDO11), use the previous ascites exposures that Geoff computed
## this is the sample in which of the ascites was a vial with multiple cell populations, not just the one corresponding to the organoids
previous_ascites <- readRDS("../copy_number_analysis_organoids/robjects/fig4_ascites.RDS")
sigs_ascites['14369.A004',] = unlist(previous_ascites[which(previous_ascites$sample == 'PDO11') - 1,-c(1,9,10)])

## read the exposures from the organoids
exposures_orgs <- readRDS("../copy_number_analysis_organoids/robjects/exposures.RDS")

rownames(sigs_ascites) %in% segs@phenoData@data$name
segs@phenoData@data$name %in% rownames(sigs_ascites)

pdf("figures/barplot_all_ascites.pdf", width=6, height=4)
createBarplot(sigs_ascites, angle_rotation_axis = 45)
dev.off()

rownames(sigs_ascites)[rownames(sigs_ascites) == "16421.D705tp-D501tp"] =  "16421.D705-D501"
rownames(sigs_ascites)[rownames(sigs_ascites) == "16421.D705tp-D503tp"] = "16421.D705-D503"

previous_match <- F
if(previous_match){
  ascc <- readxl::read_excel("data/AscitesSLXforOrganoidProject.xlsx")
  
  ascc
  sigs_ascites
  
  gsub("*.[.]", "", rownames(sigs_ascites))
  
  # organoidsMatched <- ascc$Derived_organoid[match(gsub(".*[.]","", rownames(sigs_ascites)), ascc$sWGS_barcode)]
  
  ###' BELOW:::: using the previous name, which is the correct match of ascites and PDO
  ###' i.e. the first three columns are the correct match, but the names for the
  ###' ascites are matched to the PDO number in INCORRECT_PREVIOUS_ORGANOID_NAME
  organoidsMatched <- ascc$INCORRECT_PREVIOUS_ORGANOID_NAME[match(rownames(sigs_ascites), paste0(ascc$sWGS_SLX, '.', ascc$sWGS_barcode))]
  
  ###### BELOW::: it leads to an incorrect match
  # organoidsMatched <- ascc$Derived_organoid[match(rownames(sigs_ascites), paste0(ascc$sWGS_SLX, '.', ascc$sWGS_barcode))]
  
  add_to_figs <- ''
}else{
  ## using my final matching
  add_to_figs <- '_latest'
  ascc <- readxl::read_excel("data/matching_ascites_samples_Lena.xlsx")
  stopifnot(max(table(paste0(ascc$sWGS_SLX, ascc$sWGS_barcode))) == 1)
  stopifnot(max(table(ascc$Derived_organoid)) == 1)
  organoidsMatched <- ascc$Derived_organoid[match(rownames(sigs_ascites), paste0(ascc$sWGS_SLX, '.', ascc$sWGS_barcode))]
  
}
  
rownames(sigs_ascites)[is.na(organoidsMatched)]

sigs_ascites <- sigs_ascites[!is.na(organoidsMatched),]
organoidsMatched <- organoidsMatched[!is.na(organoidsMatched)]
exposures_orgs[organoidsMatched,]

exp2 <- exposures_orgs[organoidsMatched,]
rownames(exp2)[duplicated(rownames(exp2))] <- paste0(rownames(exp2)[duplicated(rownames(exp2))], '_2')

sigs_ascites
pdf(paste0("figures/ascites_all_matched", add_to_figs, ".pdf"), width=6, height=7)
grid.arrange(createBarplot(sigs_ascites, angle_rotation_axis = 45),
             createBarplot(exp2, angle_rotation_axis = 45))
dev.off()


sigs_ascites

make_unique <- function(i){
  for(j in unique(i)){
    if(sum(i == j) > 1){
      i[i == j] = paste0((i[i == j]), '_', 1:(sum(i == j)))
    }
  }
  i
}

make_unique(ascc$OV04[match(organoidsMatched, ascc$Derived_organoid)])

df_ascites_orgs <- rbind.data.frame(cbind.data.frame(sigs_ascites, group='ascites',
      OVO4=make_unique(ascc$OV04[match(organoidsMatched, ascc$Derived_organoid)])),
      cbind.data.frame(exposures_orgs[organoidsMatched,], group='organoids', OVO4=make_unique(ascc$OV04[match(organoidsMatched, ascc$Derived_organoid)])))
df_ascites_orgs$org = organoidsMatched
df_ascites_orgs <- df_ascites_orgs[!(df_ascites_orgs[,'OVO4'] == '466_1'),] ## remove this pair due to bad quality (underpowered)

df_ascites_orgs$OVO4 <- gsub("_1", " (1)", df_ascites_orgs$OVO4)
df_ascites_orgs$OVO4 <- gsub("_2", " (2)", df_ascites_orgs$OVO4)

ggplot(melt(df_ascites_orgs,
          id.vars = c('group', 'OVO4', 'org')), aes(x=group, y=value, fill=variable))+
  geom_bar(stat="identity")+facet_wrap(.~OVO4, ncol=6)+
  scale_fill_brewer(palette="Dark2")+
  theme(legend.position = "bottom")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave(paste0("figures/barplot_all_ascites_2", add_to_figs, ".pdf"), width=6, height=5)

saveRDS(melt(df_ascites_orgs,
             id.vars = c('group', 'OVO4', 'org')), file = "robjects/ascites_organoid_exposures.RDS")

df_ascites_orgs_melt <- melt(df_ascites_orgs,
     id.vars = c('group', 'OVO4'))
for(ovo4_it in unique(df_ascites_orgs$OVO4)){
  ggplot(df_ascites_orgs_melt[df_ascites_orgs_melt$OVO4 == ovo4_it,],
         aes(x=group, y=value, fill=variable))+
    geom_bar(stat="identity")+facet_wrap(.~OVO4, ncol=6)+
    scale_fill_brewer(palette="Dark2")+
    theme(legend.position = "bottom")+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+labs(fill='', x='')+guides(fill='none')
  ggsave(paste0("figures/barplot_all_ascites_per_ovo4/", ovo4_it, add_to_figs, ".pdf"), width=2, height=3)
}

add_to_rownames <- function(i,j){
  rownames(i) <- paste0( rownames(i), '_', j)
  i
}

### ----- Previous ascites ------
ascites_geoff <- previous_ascites
# ascites_geoff = readRDS("../copy_number_analysis_organoids/robjects/fig4_ascites.RDS")
rownames(ascites_geoff) <- ascites_geoff$sample
ascites_geoff$sample <- NULL; ascites_geoff$bool_ascites <- NULL; ascites_geoff$sample_paired <- NULL
rownames(ascites_geoff)[c(T,F)] <- paste0(rownames(ascites_geoff)[c(T,F)], '(i.e. ', rownames(ascites_geoff)[c(F,T)], ')')
### ------------------------------


rbind_all_exposures <- rbind(add_to_rownames(sigs_ascites, paste0(ascc$OV04[match(organoidsMatched, ascc$Derived_organoid)],
                                                                  '_', organoidsMatched)),
      add_to_rownames(exposures_orgs[organoidsMatched,], ascc$OV04[match(organoidsMatched, ascc$Derived_organoid)]),
      add_to_rownames(ascites_geoff, 'geoff'))

png(paste0("~/Desktop/heatmap_organoids_2", add_to_figs, ".png"), height = 8, width = 6, res = 300, units = "in")
ComplexHeatmap::Heatmap(rbind_all_exposures)
dev.off()

