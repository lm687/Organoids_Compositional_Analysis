rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

require(ggplot2)
require(ggrepel)
require(cowplot)
require(gridExtra)
require(latex2exp)
require(scales)
require(pheatmap)
require(dplyr)
require(RColorBrewer)

subset_genes_of_interest = c('MYC', 'CCNE1', 'PIK3CA', 'TERT', 'KRAS', 'PTEN', 'RB1', 'AKT1',
                             'AKT2', 'PARP1', 'PARP2', 'ATM', 'ATR', 'WEE1', 'TOP1', 'TUBB1',
                             'AKT3', 'CCND1', 'CCND2', 'CCND3', 'CDKN2A', 'CDKN2B', 'MECOM', 'CDK12')

df_gene_characteristics <- readRDS("../RNASeq_and_CN/20191218_ViasM_BJ_orgaBrs/output/output_GRCh37/fig3_df_gene_characteristics.RDS")
df_average_bottomCN <- readRDS("../RNASeq_and_CN/20191218_ViasM_BJ_orgaBrs/output/output_GRCh37/fig3_df_average_bottomCN.RDS")
df_average_bottomCN <- df_average_bottomCN[!is.na(df_average_bottomCN$Gene),]
select_genes_mostvar <- readRDS("../RNASeq_and_CN/20191218_ViasM_BJ_orgaBrs/output/output_GRCh37/joint_counts_CN_subset.RDS")
df_average_bottomCN = df_average_bottomCN[df_average_bottomCN$Gene %in% unique(select_genes_mostvar$CN.gene_name),]
pca_with_gsva_annotation_NC <- readRDS("../RNASeq_DE_resistant_sensitive/objects/fig4_pca_with_gsva_annotation_NC.RDS")
pca_with_gsva_annotation_NC_prcomp <- readRDS("../RNASeq_DE_resistant_sensitive/objects/fig4_pca_with_gsva_annotation_NC_prcomp.RDS")
df_colmeans_deseqcounts_correlation_tcga_org <- readRDS("../RNASeq_DE_resistant_sensitive/objects/fig4_df_colmeans_deseqcounts_correlation_tcga_org.RDS")
ssgsea_repair <- readRDS("../RNASeq_DE_resistant_sensitive/objects/fig3_ssgsea_repair.RDS")

.x <- df_colmeans_deseqcounts_correlation_tcga_org[df_colmeans_deseqcounts_correlation_tcga_org$TME == "TME",]
.x$TME <- "Other"
df_colmeans_deseqcounts_correlation_tcga_org = rbind(df_colmeans_deseqcounts_correlation_tcga_org, .x)
df_colmeans_deseqcounts_correlation_tcga_org$TME = sapply(df_colmeans_deseqcounts_correlation_tcga_org$TME,
   function(i) ifelse(i == 'TME', yes = 'Consensus TME', no = 'All genes'))
# a <- ggplot(df_colmeans_deseqcounts_correlation_tcga_org,
#        aes(x=log(means_tcga), y=log(means_org), col=TME))+geom_point()+
#   # scale_x_continuous(trans = "log2")+#scale_y_continuous(trans = "log2")+
#   geom_abline(slope = 1, intercept = 0, lty='dashed')+facet_wrap(.~TME)+
#   theme(legend.position = "bottom")+#ggtitle('Comparison of DESeq counts between TCGA\nand organoid samples')+
#   theme_bw()+labs(x='Count means for TCGA (log2)', y='Count means for organoids (log2)')+theme(legend.position = "bottom")
df_colmeans_deseqcounts_correlation_tcga_org$xidentity=c(NA,rep(c(.01, max(df_colmeans_deseqcounts_correlation_tcga_org$means_tcga, na.rm=T)),
                                                           nrow(df_colmeans_deseqcounts_correlation_tcga_org)/2))
df_colmeans_deseqcounts_correlation_tcga_org$yidentity=c(NA,rep(c(.01, max(df_colmeans_deseqcounts_correlation_tcga_org$means_tcga, na.rm=T)),
                                                                nrow(df_colmeans_deseqcounts_correlation_tcga_org)/2))
# df_colmeans_deseqcounts_correlation_tcga_org$yidentity=c(NA,rep(c(.01, max(df_colmeans_deseqcounts_correlation_tcga_org$means_org, na.rm=T)),
#                                                            nrow(df_colmeans_deseqcounts_correlation_tcga_org)/2))
a <- ggplot(df_colmeans_deseqcounts_correlation_tcga_org,
            aes(x=(means_tcga), y=(means_org), col=TME))+
  geom_point(alpha=0.2)+
  # geom_abline(slope = (1), intercept = 0, lty='dashed')
  facet_wrap(.~TME)+
  scale_x_continuous(trans = "log2",
                     labels = scales::label_scientific(digits=0))+
  # scale_y_continuous(trans = "log2")+
  scale_y_continuous(trans = "log2",
                     labels = scales::label_scientific(digits=0))+
  geom_line(aes(x=xidentity, y=yidentity), lty='dashed', col='black')+
  theme(legend.position = "bottom")+
  theme_bw()+labs(x='Count means for TCGA (log2)', y='Count means for organoids (log2)')+
  theme(legend.position = "bottom",
        plot.margin = unit(c(0.3,1.3,0.3,.3), "lines"))+
  guides(col="none")

# > rownames(pca_with_gsva_annotation_NC)
# [1] "PDO6"  "PDO5"  "PDO15" "PDO11" "PDO3"  "PDO10" "PDO17" "PDO2"  "PDO9"  "PDO8"  "PDO1"  "PDO12"
# [13] "PDO13" "PDO7"  "PDO4" 
pca_with_gsva_annotation_NC$germline_BRCA = c()
sum_pca <- summary(pca_with_gsva_annotation_NC_prcomp)
varexplained <- round(sum_pca$importance[2,1:2]*100, digits = 1)
# b <- ggplot(pca_with_gsva_annotation_NC, aes(x=PC1, y=PC2, col=factor(BRCA1), label=labels))+
b <- ggplot(pca_with_gsva_annotation_NC, aes(x=PC1, y=PC2, label=labels))+
  geom_point()+
  geom_text_repel()+theme_bw()+theme(legend.position = "bottom")+ theme(legend.title = element_blank())+
  labs(x=paste0('PC1 (', varexplained[1], '%)'),
       y=paste0('PC2 (', varexplained[2], '%)')) #+
  ggtitle('Number of mutations in BRCA (TAMSeq)')

# b0 <- ggplot(pca_with_gsva_annotation_NC, aes(x=PC1, y=PC2, col=factor(BRCA1), label=labels))+
  b0 <- ggplot(pca_with_gsva_annotation_NC, aes(x=PC1, y=PC2, label=labels))+
    geom_point()+
    labs(x=paste0('PC1 (', varexplained[1], '%)'),
         y=paste0('PC2 (', varexplained[2], '%)'))+
  theme_bw()+theme(legend.position = "bottom")+ theme(legend.title = element_blank())

df_average_bottomCN[df_average_bottomCN$Gene == 'ERBB2','label'] = 'ERBB2'
df_average_bottomCN$label[!(df_average_bottomCN$label %in% c('MYC', 'KRAS', 'CCNE1', 'PIK3CA', 'AKT2'))] = NA
c <- ggplot(droplevels(df_average_bottomCN),
            aes(x=Gene, y=average_comparison_CN_DESeq, label=as.character(label)))+
  geom_point()+#geom_label_repel(max.overlaps = 30, segment.size = 0.01)+
  # geom_label()+
  # geom_text(force = .01, direction = "y", nudge_x = 0, nudge_y = .1)+
  geom_text_repel(force = .8, size=3)+
  theme(axis.text.x=element_blank(), axis.line.x.bottom =element_blank(), 
        panel.grid.minor = element_line(size = 0.1, colour = "black"))+
  geom_hline(yintercept = 0.5, lty='dashed')+labs(x='Ranked genes', y='Averaged higher CN and higher GE')

d <- ggplot(df_gene_characteristics[,!duplicated(colnames(df_gene_characteristics))], aes(x=df_average_bottomCN.average_comparison_CN_DESeq, y=r2_normCNnormDESeq,
                                    label=ifelse( (df_average_bottomCN.average_comparison_CN_DESeq >= .91) &
                                                    (r2_normCNnormDESeq > 0.91),
                                                  yes = Gene, no = NA )))+
  geom_boxplot(aes(group=df_average_bottomCN.average_comparison_CN_DESeq),
               width=0.5/length(unique(df_gene_characteristics$df_average_bottomCN.average_comparison_CN_DESeq)),
               col='blue')+
  geom_point()+geom_label_repel(max.overlaps=100, size=3)+theme_bw()+
  labs(x='Averaged higher CN and higher GE', y=TeX('R^2 between CN and GE'))+coord_flip()

e <- ggplot(df_gene_characteristics[,!duplicated(colnames(df_gene_characteristics))], aes(x=df_average_bottomCN.average_comparison_CN_DESeq, y=r2_normCNnormDESeq,
                                    label=ifelse( (Gene %in% subset_genes_of_interest) | ((df_average_bottomCN.average_comparison_CN_DESeq >= .8) &
                                                    (r2_normCNnormDESeq > 0.5)),
                                                  yes = Gene, no = NA )))+
  geom_boxplot(aes(group=df_average_bottomCN.average_comparison_CN_DESeq),
               width=0.5/length(unique(df_gene_characteristics$df_average_bottomCN.average_comparison_CN_DESeq)),
               col='blue')+
  geom_point()+
  geom_label_repel(aes(col=factor(ifelse( test = Gene %in% subset_genes_of_interest,
                                          yes = 'Gene of interest', no='High CN/GE correlation'))),
                   size=3)+
  theme_bw()+labs(x='Averaged higher CN and higher GE', y=TeX('R^2 between CN and GE'))+
  theme(legend.title = element_blank(), legend.position = "bottom")

# cowplot::

pdf("fig3.pdf", width = 7, height = 7)
grid.arrange(a, b, c, d)
dev.off()

pdf("fig3_v2.pdf", width = 12, height = 4)
grid.arrange(a, b, e, nrow=1)
dev.off()

df_gene_characteristics$df_average_bottomCN.average_comparison_CN_DESeq = df_gene_characteristics$df_average_bottomCN.average_comparison_CN_DESeq*8

df_average_bottomCN_summary <- df_gene_characteristics %>%
  group_by(df_average_bottomCN.average_comparison_CN_DESeq) %>%
  summarise(median=median(r2_normCNnormDESeq, na.rm = T))
df_average_bottomCN_summary$Gene = NA
c <- ggplot(df_gene_characteristics, aes(x=df_average_bottomCN.average_comparison_CN_DESeq, y=r2_normCNnormDESeq,
                                         label=ifelse(Gene %in% subset_genes_of_interest,
                                                      yes = Gene, no = NA )))+
  geom_jitter(alpha=0.2, width = 0.03*8)+theme_bw()+
  geom_boxplot(aes(group=df_average_bottomCN.average_comparison_CN_DESeq),
               width=0.5*8/length(unique(df_gene_characteristics$df_average_bottomCN.average_comparison_CN_DESeq)),
               col='white', outlier.shape = NA, fill=NA)+
  geom_line(data=df_average_bottomCN_summary,
            aes(x=df_average_bottomCN.average_comparison_CN_DESeq, y=median,
                group=1), col='blue', size=2)+
  geom_label_repel()+
  labs(x='Averaged higher CN and higher GE', y=TeX('R^2 between CN and GE'))

subset_genes_of_interest = c('MYC', 'KRAS', 'CCNE1', 'PIK3CA', 'AKT2')

conlyboxplot <- ggplot(df_gene_characteristics, aes(x=df_average_bottomCN.average_comparison_CN_DESeq, y=r2_normCNnormDESeq,
                                         label=ifelse(Gene %in% subset_genes_of_interest,
                                                      yes = Gene, no = NA )))+
  geom_boxplot(aes(group=df_average_bottomCN.average_comparison_CN_DESeq),
               width=0.5*8/length(unique(df_gene_characteristics$df_average_bottomCN.average_comparison_CN_DESeq)),
               col='grey', outlier.shape = NA, fill=NA)+
  geom_line(data=df_average_bottomCN_summary,
            aes(x=df_average_bottomCN.average_comparison_CN_DESeq, y=median,
                group=1), col='blue', size=2)+
  geom_label_repel()+theme_bw()+
  labs(x='Averaged higher CN and higher GE', y=TeX('R^2 between CN and GE'))

d <- ggplot(df_gene_characteristics,
            aes(x=cor_normCNnormDESeq, y=averageCN,
                # label=ifelse( (log(averageCN) > 1) & (cor_normCNnormDESeq > 0.8), yes = Gene, no = NA  )))+
                label=ifelse( ((log(averageCN) > 1.85) & (cor_normCNnormDESeq > 0.5)), yes = Gene, no = NA  )))+
  geom_point(alpha=0.2)+
  geom_density_2d()+
  scale_y_continuous(trans = "log2")+
  geom_label_repel(max.overlaps = 100)+theme_bw()+
  labs(x='Correlation between CN and GE', y='Average copy number')
d

rownames(ssgsea_repair) <- stringr::str_to_title(tolower(gsub('_', ' ', gsub('KEGG_', '', rownames(ssgsea_repair)))))
rownames(ssgsea_repair) <- gsub('Erbb', 'ERBB', rownames(ssgsea_repair))

nbreaks <- 11

e <- pheatmap(ssgsea_repair, legend = F,
              breaks = seq(min(ssgsea_repair), max(ssgsea_repair), length.out = nbreaks),
              color = rev(brewer.pal(nbreaks+1, "PRGn")),
              legend_breaks=seq(min(ssgsea_repair), max(ssgsea_repair), length.out = nbreaks))

elegend <- pheatmap(ssgsea_repair, legend = T,
                    breaks = seq(min(ssgsea_repair), max(ssgsea_repair), length.out = nbreaks),
                    color = rev(brewer.pal(nbreaks+1, "PRGn")),
                    legend_breaks=seq(min(ssgsea_repair), max(ssgsea_repair), length.out = nbreaks))
elegend
name_rect <- names(elegend$gtable$grobs[[6]]$children)[grep('rect', names(elegend$gtable$grobs[[6]]$children))]
name_text <- names(elegend$gtable$grobs[[6]]$children)[grep('text', names(elegend$gtable$grobs[[6]]$children))]

plot_grid(elegend$gtable$grobs[[6]])
plot_grid(elegend$gtable$grobs[[6]]$children[[name_rect]])
plot_grid(elegend$gtable$grobs[[6]]$children[[name_text]])

plt <- elegend$gtable$grobs[[6]]$children[[name_rect]]
plt$y <- plt$x
plt$x <- elegend$gtable$grobs[[6]]$children[[name_rect]]$y
plt$width <- plt$height
plt$height <- elegend$gtable$grobs[[6]]$children[[name_rect]]$width
plot_grid(plt)

plt2 <- elegend$gtable$grobs[[6]]$children[[name_text]]
plt2$y <- plt2$x
plt2$x <- elegend$gtable$grobs[[6]]$children[[name_text]]$y - unit(10, "bigpts")
plt2$width <- plt2$height
plt2$height <- elegend$gtable$grobs[[6]]$children[[name_text]]$width
plot_grid(plt2)

plot_grid(plt)
plot_grid(plt2)

elegend$gtable$grobs[[6]]$children[[name_rect]] = plt
elegend$gtable$grobs[[6]]$children[[name_text]] = plt2
elegend$gtable$widths[1] <- unit(0, 'bigpts')
elegend$gtable$heights[2] <- unit(0, 'bigpts')
elegend$gtable$grobs[[1]] <- NULL; elegend$gtable$grobs[[1]] <- NULL; elegend$gtable$grobs[[1]] <- NULL; elegend$gtable$grobs[[1]] <- NULL; elegend$gtable$grobs[[1]] <- NULL

elegend$gtable$layout <- elegend$gtable$layout[6,]
elegend$gtable$layout[4] = 2
elegend$gtable$widths[3] <- unit(0, 'bigpts')
elegend$gtable$widths[5] <- unit(0, 'bigpts')
elegend$gtable$heights[4] <- unit(2, 'bigpts')
elegend$gtable$heights[5] <- unit(0, 'bigpts')
elegend$gtable$widths[1] <- unit(0, 'bigpts')
elegend$gtable$widths[3] <- unit(0, 'bigpts')
# elegend$gtable$widths[4] <- unit(500, 'bigpts')
elegend$gtable$grobs[[1]]$children[[name_text]]$gp$fontsize = 8
elegend$gtable$layout
# elegend$gtable$layout[1] = 2
# elegend$gtable$layout[2] = 5
# elegend$gtable$layout[5] = 8
elegend$gtable$vp$width = unit(2, "npc")
elegend$gtable$vp$x = unit(0.3, "npc")

elegend$gtable$grobs[[1]]$children[[name_text]]$label[c(T,F)] = ""
elegend$gtable$grobs[[1]]$children[[name_text]]$label = substr(elegend$gtable$grobs[[1]]$children[[name_text]]$label, 1, 4)

plot_grid(elegend$gtable)

elegend$gtable$layout
elegend$gtable$widths

e_composite <- plot_grid(e[[4]],
          elegend$gtable, nrow = 2, rel_heights = c(4, 0.8))
e_composite

# pdf("fig3_v4.pdf", width = 10, height = 6.5)
# plot_grid(plot_grid(a, b0, b, rel_widths=c(3,2,2), ncol=3, labels=c('a', 'b')), plot_grid(c, d, e_composite, labels=c('c', 'd', 'e'),
#                                                                                           ncol=3, rel_widths=c(1.8,1.6,2.5)),
#           nrow=2, rel_heights = c(1.8, 2))
# dev.off()

pdf("fig3_v4.pdf", width = 10, height = 6.5)
plot_grid(plot_grid(a, b0, b, rel_widths=c(3,2,2), ncol=3, labels=c('a', 'b')), plot_grid(d, c, e_composite, labels=c('c', 'd', 'e'), ncol=3, rel_widths=c(1.6,1.8,2.5)), nrow=2, rel_heights = c(1.8, 2))
dev.off()

pdf("fig3_v5.pdf", width = 10, height = 6.5)
plot_grid(plot_grid(a, b0, b, rel_widths=c(3,2,2), ncol=3, labels=c('a', 'b')), plot_grid(d, conlyboxplot, e_composite, labels=c('c', 'd', 'e'), ncol=3, rel_widths=c(1.6,1.8,2.5)), nrow=2, rel_heights = c(1.8, 2))
dev.off()

### Recreating a figure that I thought we already had as it was in the manuscript - maybe the one in the manuscript had been manually edited
pdf("fig3_v6.pdf", width = 10, height = 6.5)
plot_grid(plot_grid(a, b, d, rel_widths=c(3,2,2), ncol=3, labels=c('a', 'b', 'c')), plot_grid(c, plot_grid(plot.new(), e_composite, plot.new(), nrow=1, rel_widths=c(1,4,1)), labels=c('d', 'e'),
                                                                                          ncol=2, rel_widths=c(1.8,2.6)),
          nrow=2, rel_heights = c(1.8, 2))
dev.off()

pdf("fig3_v5.pdf", width = 10, height = 6.5)
plot_grid(plot_grid(a, b0, b, rel_widths=c(3,2,2), ncol=3, labels=c('a', 'b')), plot_grid(d, conlyboxplot, e_composite, labels=c('c', 'd', 'e'), ncol=3, rel_widths=c(1.6,1.8,2.5)), nrow=2, rel_heights = c(1.8, 2))
dev.off()

