rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

## Code from Carolin Sauer 2020
library(tidyverse)
library(ggpubr)
library(RColorBrewer)
library(ComplexHeatmap)
library(reshape2)

renaming <- readxl::read_excel("../../RNASeq_DE_resistant_sensitive/files/PDOnameProperSample_sWGS_RNAseq.xlsx")

#----------------------------------------------------------------------------------#

AUC_all <- read_csv("../data/20200419-AUC-organoids.csv")
AUC <- read_csv("../data/20200419-AUC-organoids_PFI.csv")
colnames(AUC_all)[colnames(AUC_all) == "...1"] = "X1"
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


#### drugs - auc of several drugs

# AUC_all_df = readRDS("../survival_analysis/robjects/AUC_all_df.RDS")
AUC_all_df_clean <- AUC_all_df[,-c(1:7)]
## remove APR.246
AUC_all_df_clean$APR.246 <- NULL

AUC_all_df_clean$PFI <- NULL
AUC_all_df_clean <- (melt(AUC_all_df_clean, id.vars='PDO'))

PDS_PDO_summary <- AUC_all_df_clean %>% group_by(PDO) %>% summarise(stddevAUC=sd(value),
                                                                    meanAUC=mean(value))
n <- 60
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unique(unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))))


## for a second, change from po to the original ID
AUC_all_df_clean$PDO = renaming$ID[match(AUC_all_df_clean$PDO, renaming$PDO)]
PDS_PDO_summary$PDO = renaming$ID[match(PDS_PDO_summary$PDO, renaming$PDO)]


AUC_all_df_clean$PDO = factor(AUC_all_df_clean$PDO,
                              levels=PDS_PDO_summary$PDO[order(PDS_PDO_summary$meanAUC)])
PDS_PDO_summary$PDO = factor(PDS_PDO_summary$PDO,
                             levels=PDS_PDO_summary$PDO[order(PDS_PDO_summary$meanAUC)])

left_group <- sort(PDS_PDO_summary$PDO[which(PDS_PDO_summary$meanAUC < mean(PDS_PDO_summary$meanAUC))])
rightmost_PDO_of_left <- left_group[length(left_group)]
right_group <- sort(PDS_PDO_summary$PDO[which(PDS_PDO_summary$meanAUC > mean(PDS_PDO_summary$meanAUC))])
leftmost_PDO_of_right <- right_group[1]

coef_std <- 3
ggplot()+
  theme_bw()+
  geom_line(data = AUC_all_df_clean, aes(y=value, x=PDO, group=variable, col=variable))+
  geom_point(data = AUC_all_df_clean, aes(y=value, x=PDO, group=variable, col=variable), shape=5)+
  geom_line(data = PDS_PDO_summary, aes(x=PDO, y=meanAUC, group=1),
            col='black', lty=1, size=1.1)+
  geom_line(data = PDS_PDO_summary, aes(x=PDO, y=stddevAUC*coef_std, group=1),
            col='blue', lty='dashed', size=1.1)+
  scale_color_manual(name = "Organoid",values = col_vector)+
  scale_y_continuous(
    
    # Features of the first axis
    name = "AUC",
    
    # Add a second axis and specify its features
    sec.axis = sec_axis(~./coef_std, name="Standard deviation")
  )+
  theme(
    axis.title.y = element_text(color = 'black', size=13),
    axis.title.y.right = element_text(color = 'blue', size=10),
    axis.line.y.right = element_line(color = "blue"),
    axis.ticks.y.right = element_line(color = "blue"),
    axis.text.y.right = element_text(color = "blue")
  )+
  geom_vline(xintercept = leftmost_PDO_of_right, lty='dashed')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("../../survival_analysis/figures/drugs_AUC_all_lineplot.pdf", width = 10, height = 4)

## now only using 4 drugs
AUC_all_df_clean_4drugs <- AUC_all_df_clean %>% dplyr::filter(variable %in% c('Paclitaxel',
                                                                                     'Oxaliplatin',
                                                                                     'Doxorubicin',
                                                                                     'Gemcitabine' ))
# AUC_all_df_clean_4drugs$PDO = paste0(AUC_all_df_clean_4drugs$PDO, ' - ',
#                                      renaming$PDO[match(AUC_all_df_clean_4drugs$PDO, renaming$ID)])
AUC_all_df_clean_4drugs$PDO = renaming$PDO[match(AUC_all_df_clean_4drugs$PDO, renaming$ID)]
PDS_PDO_summary_4drugs <- AUC_all_df_clean_4drugs %>% group_by(PDO) %>% summarise(stddevAUC=sd(value),
                                                                    meanAUC=mean(value))

AUC_all_df_clean_4drugs$PDO = factor(AUC_all_df_clean_4drugs$PDO,
                              levels=PDS_PDO_summary_4drugs$PDO[order(PDS_PDO_summary_4drugs$meanAUC)])
PDS_PDO_summary_4drugs$PDO = factor(PDS_PDO_summary_4drugs$PDO,
                             levels=PDS_PDO_summary_4drugs$PDO[order(PDS_PDO_summary_4drugs$meanAUC)])

left_group <- sort(PDS_PDO_summary_4drugs$PDO[which(PDS_PDO_summary_4drugs$meanAUC < median(PDS_PDO_summary_4drugs$meanAUC))])
rightmost_PDO_of_left <- left_group[length(left_group)]
right_group <- sort(PDS_PDO_summary_4drugs$PDO[which(PDS_PDO_summary_4drugs$meanAUC > median(PDS_PDO_summary_4drugs$meanAUC))])
leftmost_PDO_of_right <- right_group[1]

AUC_all_df_clean_4drugs

ggplot()+
  theme_bw()+
  geom_line(data = AUC_all_df_clean_4drugs, aes(y=value, x=PDO, group=variable, col=variable))+
  geom_point(data = AUC_all_df_clean_4drugs, aes(y=value, x=PDO, group=variable, col=variable), shape=5)+
  geom_line(data = PDS_PDO_summary_4drugs, aes(x=PDO, y=meanAUC, group=1),
            col='black', lty=1, size=1.1)+
  geom_line(data = PDS_PDO_summary_4drugs, aes(x=PDO, y=stddevAUC*coef_std, group=1),
            col='blue', lty='dashed', size=1.1)+
  scale_color_manual(name = "Organoid",values = col_vector)+
  scale_y_continuous(
    
    # Features of the first axis
    name = "AUC",
    
    # Add a second axis and specify its features
    sec.axis = sec_axis(~./coef_std, name="Standard deviation")
  )+
  theme(
    axis.title.y = element_text(color = 'black', size=13),
    axis.title.y.right = element_text(color = 'blue', size=10),
    axis.line.y.right = element_line(color = "blue"),
    axis.ticks.y.right = element_line(color = "blue"),
    axis.text.y.right = element_text(color = "blue")
  )+
  geom_vline(xintercept = as.numeric(leftmost_PDO_of_right)-0.5, lty='dashed')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("../../survival_analysis/figures/drugs_AUC_all_lineplot_4drugs.pdf", width = 10, height = 4)

