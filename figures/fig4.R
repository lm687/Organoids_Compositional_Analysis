rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

require(ggplot2)
require(ggrepel)
require(cowplot)
require(gridExtra)
require(latex2exp)

subset_genes_of_interest = c('MYC', 'CCNE1', 'PIK3CA', 'TERT', 'KRAS', 'PTEN', 'RB1', 'AKT1',
                             'AKT2', 'PARP1', 'PARP2', 'ATM', 'ATR', 'WEE1', 'TOP1', 'TUBB1',
                             'AKT3', 'CCND1', 'CCND2', 'CCND3', 'CDKN2A', 'CDKN2B', 'MECOM', 'CDK12')

df_gene_characteristics <- readRDS("../RNASeq_and_CN/20191218_ViasM_BJ_orgaBrs/output/fig4_df_gene_characteristics.RDS")
df_average_bottomCN <- readRDS("../RNASeq_and_CN/20191218_ViasM_BJ_orgaBrs/output/fig4_df_average_bottomCN.RDS")
pca_with_gsva_annotation_NC <- readRDS("../RNASeq_DE_resistant_sensitive/objects/fig4_pca_with_gsva_annotation_NC.RDS")
df_colmeans_deseqcounts_correlation_tcga_org <- readRDS("../RNASeq_DE_resistant_sensitive/objects/fig4_df_colmeans_deseqcounts_correlation_tcga_org.RDS")

a <- ggplot(df_colmeans_deseqcounts_correlation_tcga_org,
       aes(x=means_tcga, y=means_org, col=TME))+geom_point()+
  scale_x_continuous(trans = "log2")+scale_y_continuous(trans = "log2")+
  geom_abline(slope = 1, intercept = 0, lty='dashed')+facet_wrap(.~TME)+
  theme(legend.position = "bottom")+#ggtitle('Comparison of DESeq counts between TCGA\nand organoid samples')+
  theme_bw()+labs(x='Count means for TCGA', y='Count means for organoids')+theme(legend.position = "bottom")

b <- ggplot(pca_with_gsva_annotation_NC, aes(x=PC1, y=PC2, col=factor(BRCA1), label=labels))+
  geom_point()+
  geom_label_repel()+theme_bw()+theme(legend.position = "bottom")+ theme(legend.title = element_blank())#+
  ggtitle('Number of mutations in BRCA (TAMSeq)')

c <- ggplot(droplevels(df_average_bottomCN),
            aes(x=Gene, y=average_comparison_CN_DESeq, label=as.character(label)))+
  geom_point()+#geom_label_repel(max.overlaps = 30, segment.size = 0.01)+
  # geom_label()+
  geom_text_repel(force = .01, direction = "y", nudge_x = 0, nudge_y = .1)+
  theme(axis.text.x=element_blank(), axis.line.x.bottom =element_blank(), 
        panel.grid.minor = element_line(size = 0.1, colour = "black"))+
  geom_hline(yintercept = 0.5, lty='dashed')

d <- ggplot(df_gene_characteristics, aes(x=df_average_bottomCN.average_comparison_CN_DESeq, y=r2_normCNnormDESeq,
                                    label=ifelse( (df_average_bottomCN.average_comparison_CN_DESeq >= .8) &
                                                    (r2_normCNnormDESeq > 0.5),
                                                  yes = Gene, no = NA )))+
  geom_boxplot(aes(group=df_average_bottomCN.average_comparison_CN_DESeq),
               width=0.5/length(unique(df_gene_characteristics$df_average_bottomCN.average_comparison_CN_DESeq)),
               col='blue')+
  geom_point()+geom_label_repel()+theme_bw()+labs(x='Averaged higher CN and higher GE', y=TeX('R^2 between CN and GE'))

e <- ggplot(df_gene_characteristics, aes(x=df_average_bottomCN.average_comparison_CN_DESeq, y=r2_normCNnormDESeq,
                                    label=ifelse( (Gene %in% subset_genes_of_interest) | ((df_average_bottomCN.average_comparison_CN_DESeq >= .8) &
                                                    (r2_normCNnormDESeq > 0.5)),
                                                  yes = Gene, no = NA )))+
  geom_boxplot(aes(group=df_average_bottomCN.average_comparison_CN_DESeq),
               width=0.5/length(unique(df_gene_characteristics$df_average_bottomCN.average_comparison_CN_DESeq)),
               col='blue')+
  geom_point()+geom_label_repel(aes(col=factor(ifelse( test = Gene %in% subset_genes_of_interest, yes = 'Gene of interest', no='High CN/GE correlation'))))+
  theme_bw()+labs(x='Averaged higher CN and higher GE', y=TeX('R^2 between CN and GE'))+
  theme(legend.title = element_blank(), legend.position = "bottom")

# cowplot::

pdf("fig4.pdf", width = 7, height = 7)
grid.arrange(a, b, c, d)
dev.off()

pdf("fig4_v2.pdf", width = 12, height = 4)
grid.arrange(a, b, e, nrow=1)
dev.off()
