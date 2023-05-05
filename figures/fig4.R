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
# renaming_ascites <- readxl::read_excel("../copy_number_analysis_organoids/data/Book1_ascites.xlsx")
renaming_ascites <- readxl::read_excel("../copy_number_analysis_organoids/data/matching_ascites_samples_Lena.xlsx")
ascites = readRDS("../copy_number_analysis_organoids/robjects/fig4_ascites.RDS")
AUC_all_df = readRDS("../survival_analysis/robjects/AUC_all_df.RDS")
PDS_PDO = readRDS("../survival_analysis/robjects/fig4_PDS_PDO.RDS")
exposures_org <- readRDS("../copy_number_analysis_organoids/robjects/exposures.RDS")
fgseaResTidy <- readRDS("../RNASeq_DE_resistant_sensitive/objects/fig4_fgseaResTidy.RDS")

ascites_df <- melt(ascites, id.vars=c('sample', 'bool_ascites', 'sample_paired'))

renaming_ascites
# ascites_df$patients = renaming_patients$OV04_number[match(renaming$ID[match(ascites_df$sample, renaming$PDO)],
#                                     renaming_patients$organoid_name)]
ascites_df$patients = renaming_patients$OV04_number[match(renaming$ID[match(ascites_df$sample, renaming$PDO)],
                                                          renaming_patients$organoid_name)]
ascites_df$patients[c(T,F)] = ascites_df$patients[c(F,T)]

b <- ggplot(ascites_df, aes(x=bool_ascites, y=value, fill=variable))+geom_bar(stat = "identity")+
  facet_wrap(.~patients, scales = "free_x", nrow=2)+
  scale_fill_brewer(palette="Dark2", name = NULL)+theme(legend.position="bottom")+labs(x=NULL, y='Copy number signature activity')+
  guides(fill=guide_legend(nrow=1))


PDS_PDO$sample <- renaming$PDO[match(as.character(PDS_PDO$sample), gsub("org", "", renaming$ID))]
PDS_PDO$name = paste0(PDS_PDO$sample, ' (',
                      renaming_patients$OV04_number[match(renaming$ID[match(PDS_PDO$sample, renaming$PDO)],
                                                          renaming_patients$organoid_name)], ')')

renaming_patients
c <- ggplot(data = PDS_PDO, aes(auc_ll5.sph, auc_ll5.org, colour = sample))+
  geom_point() +
  labs(x="AUC (Patient-derived ascites)", y="AUC (Patient-derived organoids)")+
  theme_bw() +
  geom_smooth(method=lm, se=FALSE)+facet_wrap(.~name, ncol=1)+theme(legend.position = "bottom")+
  guides(col=FALSE)

# d <- pheatmap(data.frame(AUC_all_df[,-(c(1:7, 21:22))]), annotation_row = AUC_all_df %>% select(PFI),
#               annotation_colors = list(Resistant='red', Sensitive='purple'))

drugs_name <- read.table("additional_files/drugs_nomenclature.txt", header = T)
drugs_name
rename_drugs <- drugs_name$name[match(colnames(AUC_all_df), drugs_name$company_name)]
rename_drugs[is.na(rename_drugs)] = colnames(AUC_all_df)[is.na(rename_drugs)]
colnames(AUC_all_df) <- rename_drugs
AUC_all_df <- AUC_all_df[order(AUC_all_df$PFI),]

AUC_all_df$PFI = ifelse(AUC_all_df$PFI == 'Sensitive', yes = 'Platinum sensitive', no = "Platinum resistant")
annotation_left <- ComplexHeatmap::rowAnnotation(df =  AUC_all_df %>% dplyr::select(PFI), show_annotation_name =F)
# annotation_left@anno_list$PFI@show_legend=T
# annotation_left@anno_list$PFI@name_param$show=F
annotation_left@anno_list$PFI@color_mapping@name = 'PFI'
d <- ComplexHeatmap::Heatmap(as(AUC_all_df[,-(c(1:7, 21:22))], 'matrix'),
                             left_annotation = annotation_left,
                             heatmap_legend_param = list(direction = "horizontal",heatmap_legend_side = "bottom"),
                             name='AUC', row_names_gp = grid::gpar(fontsize = 10),
                             column_names_gp = grid::gpar(fontsize = 10), cluster_rows = F)
d
AUC_all_df <- AUC_all_df[,colnames(AUC_all_df) %in% c("Paclitaxel", "Oxaliplatin",
                                                      "Doxorubicin", "Gemcitabine")]
# d_v2 <- ComplexHeatmap::Heatmap(as(AUC_all_df[,-(c(1:7, 21:22))], 'matrix'),
d_v2 <- ComplexHeatmap::Heatmap(as(AUC_all_df, 'matrix'),
                                heatmap_legend_param = list(direction = "horizontal",
                                                         heatmap_legend_side = "bottom"),
                             name='AUC', row_names_gp = grid::gpar(fontsize = 10),
                             column_names_gp = grid::gpar(fontsize = 10), cluster_rows = T)
d_v2
# d@left_annotation@param
# # d@name = ""
# d@row_names_param$anno@show_name
# d@left_annotation$PFI
# d@layout
# d <- draw(d, heatmap_legend_side = "bottom")


d2 <- ggplot(droplevels(melt(exposures_org[rownames(AUC_all_df),])),
       aes(x=value,y = droplevels(factor(Var1, levels=rev(rownames(AUC_all_df)[row_order(d)]))), fill=Var2))+
  geom_bar(stat = "identity")+
  scale_fill_brewer(palette="Dark2")+labs(y='Exposure')+
  guides(fill=FALSE)+theme_bw()+ theme(plot.margin = unit(c(2.35,0,4.9,0), "lines"))+
  theme(axis.line=element_blank(),axis.text.x=element_text(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_text(),
        axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())+
  labs(x='Copy number\nsignature activity')

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
a <- plot_grid(plot_grid(cowplot::ggdraw()+draw_image("../copy_number_analysis_organoids/data/absolute_profiles/PDO1_annotated_1.pdf", scale = 1),
                         plots_small_exposures$PDO1, rel_heights=c(5,1), nrow=2),
               plot_grid(cowplot::ggdraw()+draw_image("../copy_number_analysis_organoids/data/absolute_profiles/PDO2_annotated_1.pdf", scale = 1),
                         plots_small_exposures$PDO2, rel_heights=c(5,1), nrow=2),
               plot_grid(cowplot::ggdraw()+draw_image("../copy_number_analysis_organoids/data/absolute_profiles/PDO3_annotated_1.pdf", scale = 1),
                         plots_small_exposures$PDO3, rel_heights=c(5,1), nrow=2),
               plot_grid(cowplot::ggdraw()+draw_image("../copy_number_analysis_organoids/data/absolute_profiles/PDO4_annotated_1.pdf", scale = 1),
                         plots_small_exposures$PDO4, rel_heights=c(5,1), nrow=2))

# pdf("fig4.pdf", height = 9, width = 9)
# plot_grid(plot_grid(a, b, labels = c('a', 'b'), rel_widths = c(2,2)),
#           plot_grid(c, plot_grid( grid.grabExpr(draw(d, heatmap_legend_side = "bottom")),
#                                   d2, rel_widths = c(3,1)), labels = c('c', 'd'), rel_widths = c(3,4)),
#           label_size = 12, ncol=1, scale = 0.9, rel_heights = c(2,2))
# dev.off()

# pdf("fig4.pdf", height = 9, width = 9)
# plot_grid(plot_grid(a, b, labels = c('a', 'b'), rel_widths = c(2,2)),
#           plot_grid(c, plot_grid( grid.grabExpr(draw(d_v2, heatmap_legend_side = "bottom")),
#                                   d2, rel_widths = c(3,1)), labels = c('c', 'd'), rel_widths = c(3,4)),
#           label_size = 12, ncol=1, scale = 0.9, rel_heights = c(2,2))
# dev.off()


fgseaResTidy$pathway <- stringr::str_to_title(tolower(gsub("HALLMARK ", "", gsub("_", " ", fgseaResTidy$pathway))))
fgseaResTidy$P <- stringr::str_to_title(tolower(gsub("HALLMARK ", "", gsub("_", " ", fgseaResTidy$pathway))))
fgseaResTidy$Pathway <- stringr::str_to_title(tolower(gsub("HALLMARK ", "", gsub("_", " ", fgseaResTidy$pathway))))
## change the names of genes
fgseaResTidy$Pathway <- gsub('Myc', 'MYC', fgseaResTidy$Pathway)
fgseaResTidy$Pathway <- gsub('E2f', 'E2F', fgseaResTidy$Pathway)
fgseaResTidy$Pathway <- gsub('G2m', 'G2M', fgseaResTidy$Pathway)
fgseaResTidy$Pathway <- gsub('Mtorc1', 'MTORC1', fgseaResTidy$Pathway)
fgseaResTidy$Pathway <- gsub('Uv Response', 'UV Response', fgseaResTidy$Pathway)
fgseaResTidy$Pathway <- gsub('Il6', 'IL6', fgseaResTidy$Pathway)
fgseaResTidy$Pathway <- gsub('Jak', 'JAK', fgseaResTidy$Pathway)
fgseaResTidy$Pathway <- gsub('Stat3', 'STAT3', fgseaResTidy$Pathway)
fgseaResTidy$Pathway <- gsub('Dna', 'DNA', fgseaResTidy$Pathway)
fgseaResTidy$Pathway <- gsub('Notch', 'NOTCH', fgseaResTidy$Pathway)
fgseaResTidy$Pathway <- gsub('Tnfa', 'TNFA', fgseaResTidy$Pathway)
fgseaResTidy$Pathway <- gsub('Nfkb', 'NFKB', fgseaResTidy$Pathway)
fgseaResTidy$Pathway <- gsub('Pi3k Akt Mtor', 'PI3K AKT MTOR', fgseaResTidy$Pathway)
fgseaResTidy$Pathway <- gsub('Wnt', 'WNT', fgseaResTidy$Pathway)
fgseaResTidy$Pathway <- gsub('Il2 Stat5', 'IL2 STAT5', fgseaResTidy$Pathway)
fgseaResTidy$Pathway <- gsub('Tgf', 'TGF', fgseaResTidy$Pathway)
fgseaResTidy$Pathway <- gsub('Kras', 'KRAS', fgseaResTidy$Pathway)
fgseaResTidy$Pathway <- gsub('Response Dn', 'Response DN', fgseaResTidy$Pathway)

d3 <- fgseaResTidy %>% 
  dplyr::filter(padj < 0.05) %>% 
  ggplot(aes(x = reorder(Pathway, NES), y = NES)) +
  geom_col(aes(fill=padj)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized GSEA Score"
       # title="Hallmark pathways NES"
       ) + 
  theme_minimal(7) +
  theme(legend.key.size = unit(0.1, "in"), legend.key.width = unit(0.35, "in"),
        legend.position = 'bottom', text = element_text(size=13), legend.text = element_text(size=10))

pdf("fig4_v2.pdf", height = 9, width = 9)
plot_grid(plot_grid(a, b, labels = c('a', 'b'), rel_widths = c(2,2)),
          plot_grid(c, plot_grid( grid.grabExpr(draw(d_v2, heatmap_legend_side = "bottom")),
                                  d3, rel_widths = c(.55,1),
                                  labels = c('d', 'e')), labels = c('c'), rel_widths = c(1.8,4)),
          label_size = 12, ncol=1, scale = 0.9, rel_heights = c(3,2.5))
dev.off()

pdf("fig4_d.pdf", height = 4, width = 4)
d_v2
dev.off()

## only keeping a subset of drugs for fig 4c
c_subset <- ggplot(data = PDS_PDO %>% dplyr::filter(labdrug %in% c('Doxorubicin', 'Gemcitabine',
                                                      'Oxaliplatin', 'Paclitaxel')),
                   aes(auc_ll5.sph, auc_ll5.org, colour = sample))+
  geom_point() +
  labs(x="AUC (Patient-derived ascites)", y="AUC (Patient-derived organoids)")+
  theme_bw() +
  geom_smooth(method=lm, se=FALSE)+facet_wrap(.~name, ncol=1)+theme(legend.position = "bottom")+
  guides(col=FALSE)

pdf("fig4_v3.pdf", height = 9, width = 9)
plot_grid(plot_grid(a, b, labels = c('a', 'b'), rel_widths = c(2,2)),
          plot_grid(c_subset, plot_grid( grid.grabExpr(draw(d_v2, heatmap_legend_side = "bottom")),
                                  d3, rel_widths = c(.55,1),
                                  labels = c('d', 'e')), labels = c('c'), rel_widths = c(1.8,4)),
          label_size = 12, ncol=1, scale = 0.9, rel_heights = c(3,2.5))
dev.off()

### Recreating a figure that I thought we already had as it was in the manuscript - maybe the one in the manuscript had been manually edited
pdf("fig4_v4.pdf", height = 4, width = 9)
plot_grid(plot_grid(c, plot_grid( grid.grabExpr(draw(d_v2, heatmap_legend_side = "bottom")),
                                         d3, rel_widths = c(.55,1),
                                         labels = c('b', 'c')), labels = c('a'), rel_widths = c(1.8,4)),
          label_size = 12, ncol=1, scale = 0.9, rel_heights = c(3,2.5))
dev.off()
