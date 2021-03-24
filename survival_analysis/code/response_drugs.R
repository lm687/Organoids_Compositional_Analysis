rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

## Code from Carolin Sauer 2020
library(tidyverse)
library(ggpubr)
library(RColorBrewer)
library(ComplexHeatmap)
library(reshape2)

#----------------------------------------------------------------------------------#

AUC_all <- read_csv("../data/20200419-AUC-organoids.csv")
AUC <- read_csv("../data/20200419-AUC-organoids_PFI.csv")
head(melt(AUC_all[,-(2:7)], id.vars="X1"))
melt_AUC_all <- melt(AUC_all[,-(2:7)], id.vars="X1")

melt_AUC_all$PDO = AUC$new[match(melt_AUC_all$X1, AUC$id)]
melt_AUC_all$PFI = AUC$PFI[match(melt_AUC_all$X1, AUC$id)]

AUC_all$PDO = AUC$new[match(AUC_all$X1, AUC$id)]
AUC_all$PFI = AUC$PFI[match(AUC_all$X1, AUC$id)]

ggplot(melt_AUC_all, aes(x=factor(X1), y=variable, fill=value))+geom_tile()

AUC_all_df = data.frame(AUC_all)
rownames(AUC_all_df) = AUC_all_df$PD
colnames(AUC_all_df) = gsub("auc_ll5.", "", colnames(AUC_all_df))

AUC_all_df = AUC_all_df %>% select(-Elescamol)

saveRDS(AUC_all_df, "../robjects/AUC_all_df.RDS")

pdf("../data/AUC_all_drugs.pdf")
pheatmap(data.frame(AUC_all_df[,-(c(1:7, 21:22))]), annotation_row = AUC_all_df %>% select(PFI),
         annotation_colors = list(Resistant='red', Sensitive='purple'))
dev.off()
