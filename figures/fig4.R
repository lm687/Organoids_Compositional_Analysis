rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

set.seed(234)

library(ggplot2)
library(ggplotify)
library(reshape2)
library(pheatmap)
library(cowplot)
library(dplyr)
library(readxl)
library(ComplexHeatmap)

renaming <- readxl::read_excel("../RNASeq_DE_resistant_sensitive/files/PDOnameProperSample_sWGS_RNAseq.xlsx")
renaming_patients <- read.table("../RNASeq_and_CN/20191218_ViasM_BJ_orgaBrs/Input/samples.csv", sep = ";", header = T)
renaming_ascites <- readxl::read_excel("../copy_number_analysis_organoids/data/Book1_ascites.xlsx")
ascites = readRDS("../copy_number_analysis_organoids/robjects/fig4_ascites.RDS")
AUC_all_df = readRDS("../survival_analysis/robjects/AUC_all_df.RDS")
PDS_PDO = readRDS("../survival_analysis/robjects/fig4_PDS_PDO.RDS")
exposures_org <- readRDS("../copy_number_analysis_organoids/robjects/exposures.RDS")

ascites_df <- melt(ascites, id.vars=c('sample', 'bool_ascites', 'sample_paired'))

renaming_ascites
ascites_df$patients = renaming_patients$OV04_number[match(renaming$ID[match(ascites_df$sample, renaming$PDO)],
                                    renaming_patients$organoid_name)]
ascites_df$patients[c(T,F)] = ascites_df$patients[c(F,T)]

b <- ggplot(ascites_df, aes(x=bool_ascites, y=value, fill=variable))+geom_bar(stat = "identity")+
  facet_wrap(.~patients, scales = "free_x", nrow=2)+
  scale_fill_brewer(palette="Dark2", name = NULL)+theme(legend.position="bottom")+labs(x=NULL, y=NULL)+
  guides(fill=guide_legend(nrow=1))


PDS_PDO$sample <- renaming$PDO[match(as.character(PDS_PDO$sample), gsub("org", "", renaming$ID))]
PDS_PDO$name = paste0(PDS_PDO$sample, ' (',
                      renaming_patients$OV04_number[match(renaming$ID[match(PDS_PDO$sample, renaming$PDO)],
                                                          renaming_patients$organoid_name)], ')')

renaming_patients
c <- ggplot(data = PDS_PDO, aes(auc_ll5.sph, auc_ll5.org, colour = sample))+
  geom_point() +
  labs(x="Patient-derived spheroids", y="Patient-derived organoids")+
  theme_bw() +
  geom_smooth(method=lm, se=FALSE)+facet_wrap(.~name, ncol=1)+theme(legend.position = "bottom")

# d <- pheatmap(data.frame(AUC_all_df[,-(c(1:7, 21:22))]), annotation_row = AUC_all_df %>% select(PFI),
#               annotation_colors = list(Resistant='red', Sensitive='purple'))

drugs_name <- read.table("additional_files/drugs_nomenclature.txt", header = T)
drugs_name
rename_drugs <- drugs_name$name[match(colnames(AUC_all_df), drugs_name$company_name)]
rename_drugs[is.na(rename_drugs)] = colnames(AUC_all_df)[is.na(rename_drugs)]
colnames(AUC_all_df) <- rename_drugs

d <- ComplexHeatmap::Heatmap(as(AUC_all_df[,-(c(1:7, 21:22))], 'matrix'),
                             left_annotation = (ComplexHeatmap::rowAnnotation(df =  AUC_all_df %>% select(PFI))),
                             heatmap_legend_param = list(direction = "horizontal",heatmap_legend_side = "bottom"))

# d <- draw(d, heatmap_legend_side = "bottom")


d2 <- ggplot(droplevels(melt(exposures_org[rownames(AUC_all_df),])),
       aes(x=value,y = droplevels(factor(Var1, levels=rev(rownames(AUC_all_df)[row_order(d)]))), fill=Var2))+
  geom_bar(stat = "identity")+
  scale_fill_brewer(palette="Dark2")+labs(y='Exposure')+
  guides(fill=FALSE)+theme_bw()+ theme(plot.margin = unit(c(2.2,0,8,0), "lines"))+
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

# d$gtable$layout = d$gtable$layout[!(d$gtable$layout$name == "legend"),]
# d$gtable$grobs[[9]] = NULL


org_list <- c('PDO1', 'PDO2', 'PDO3', 'PDO4')
plots_small_exposures=lapply(org_list, function(k){
  ggplot(cbind.data.frame(names=names(exposures_org[k,]), exposures=exposures_org[k,]), aes(y=exposures, x=1, fill=names))+geom_bar(stat = "identity")+ coord_flip()+
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
  axis.text.y=element_blank(),axis.ticks=element_blank(),
  axis.title.x=element_blank(),
  axis.title.y=element_blank(),legend.position="none",
  panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
  panel.grid.minor=element_blank(),plot.background=element_blank())+scale_fill_brewer(palette="Dark2")+labs(y='Exposure')
}); names(plots_small_exposures) = org_list
  
# a <- cowplot::ggdraw()+draw_image("~/Desktop/fig4a.png", scale = 1)
a <- plot_grid(plot_grid(cowplot::ggdraw()+draw_image("../copy_number_analysis_organoids/data/absolute_profiles/PDO1_annotated.pdf", scale = 1),
                         plots_small_exposures$PDO1, rel_heights=c(5,1), nrow=2),
               plot_grid(cowplot::ggdraw()+draw_image("../copy_number_analysis_organoids/data/absolute_profiles/PDO2_annotated.pdf", scale = 1),
                         plots_small_exposures$PDO2, rel_heights=c(5,1), nrow=2),
               plot_grid(cowplot::ggdraw()+draw_image("../copy_number_analysis_organoids/data/absolute_profiles/PDO3_annotated.pdf", scale = 1),
                         plots_small_exposures$PDO3, rel_heights=c(5,1), nrow=2),
               plot_grid(cowplot::ggdraw()+draw_image("../copy_number_analysis_organoids/data/absolute_profiles/PDO4_annotated.pdf", scale = 1),
                         plots_small_exposures$PDO4, rel_heights=c(5,1), nrow=2))

pdf("fig4.pdf", height = 9, width = 9)
plot_grid(plot_grid(a, b, labels = c('a', 'b')),
          plot_grid(c, plot_grid( grid.grabExpr(draw(d, heatmap_legend_side = "bottom")),
                                  d2, rel_widths = c(3,1)), labels = c('c', 'd'), rel_widths = c(3,4)),
          label_size = 12, ncol=1, scale = 0.9, rel_heights = c(1.5,2))
dev.off()




