# S9a	/Users/morril01/Documents/PhD/other_repos/Vias_Brenton/RNASeq_and_CN/figures_GRCh37/genes_of_interest_CN_log2.pdf
# S9b	/Users/morril01/Documents/PhD/other_repos/Vias_Brenton/RNASeq_and_CN/figures_GRCh37/most_highly_amplified.pdf

rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(cowplot)
library(ggpubr)
add <- ''

cowplot::ggdraw()+draw_image("../../copy_number_analysis_organoids/data/absolute_profiles/PDO1_GM_anno.pdf", scale = 1)
bbb <- cowplot::ggdraw()+draw_image("../../RNASeq_and_CN/figures_GRCh37/genes_of_interest_CN_log2.pdf", scale = 1)

a <- plot_grid(plot_grid(cowplot::ggdraw()+draw_image("../../RNASeq_and_CN/figures_GRCh37/genes_of_interest_CN_log2.pdf", scale = 1),
                         rel_heights=c(5,0), nrow=2, labels = "a", label_size = 7),
               plot_grid(cowplot::ggdraw()+draw_image("../../RNASeq_and_CN/figures_GRCh37/most_highly_amplified.pdf", scale = 1),
                         rel_heights=c(5,0), nrow=2, labels = "b", label_size = 7), nrow=2)

pdf(paste0("suppl9", add, ".pdf"), height = 8, width = 4.5)
a
dev.off()

########################################################################################

joint_counts_CN0 <- readRDS("../../RNASeq_and_CN/20191218_ViasM_BJ_orgaBrs/output/output_GRCh37_with_14_orgs/joint_counts_CN0_added_cols.RDS")
subset_genes_of_interest <- readRDS("../../RNASeq_and_CN/20191218_ViasM_BJ_orgaBrs/output/output_GRCh37_with_14_orgs/subset_genes_of_interest.RDS")
col_vector <- readRDS("../../RNASeq_and_CN/20191218_ViasM_BJ_orgaBrs/output/output_GRCh37_with_14_orgs/col_vector.RDS")

a <- ggplot(joint_counts_CN0 %>% dplyr::filter(CN.gene_name %in% subset_genes_of_interest),
       aes(x=CN.gene_name, y=log2(CN.value)))+geom_violin()+
  geom_jitter(aes(col=PDO))+
  facet_wrap(.~CN.gene_name, scales = "free", nrow=3)+
  scale_colour_manual(values = col_vector)+
  geom_hline(yintercept = log2(2), lty='dashed', col='blue')+theme_bw()+
  theme(legend.position="bottom", legend.spacing.x = unit(.15, 'cm'),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  labs(x='Gene name', y='log2(CN value)')+
  guides(col=guide_legend(nrow=2,byrow=TRUE, title = 'PDO'),)

give_order_colsums <- function(f){
  unlist(f[order(f[,2], decreasing = T),1])
}
joint_counts_CN0_MHA <- readRDS("../../RNASeq_and_CN/20191218_ViasM_BJ_orgaBrs/output/output_GRCh37_with_14_orgs/joint_counts_CN0_MHA.RDS")
most_highly_amplified_genes_list <- readRDS("../../RNASeq_and_CN/20191218_ViasM_BJ_orgaBrs/output/output_GRCh37_with_14_orgs/most_highly_amplified_genes_list.RDS")
b <- ggplot(joint_counts_CN0_MHA,
       aes(y=log2(CN.value), x=factor(CN.gene_name, levels=rev(most_highly_amplified_genes_list)),
           fill=factor(chrom, levels=sort(as.numeric(unique(chrom))))))+
  geom_bar(stat = "identity")+coord_flip()+
  facet_wrap(.~factor(gsub('PDO', '', PDO), levels = gsub('PDO', '', give_order_colsums(joint_counts_CN0_MHA %>% 
                                                          dplyr::group_by(PDO) %>% dplyr::summarise(sum(CN.value))))), nrow=1)+
  geom_hline(yintercept = log(2), col='blue')+
  labs(x='Gene name', y='CN (log2)', fill='Chromosome')+
  theme(legend.position="bottom")

pdf(paste0("suppl8_v2", add, ".pdf"), height = 7.5, width = 8)
cowplot::plot_grid(a,b, nrow=2)
dev.off()
