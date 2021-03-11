##------------------------------------------------------------------------------------------------------------##
rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

require(ggplot2)
require(ggrepel)
require(dplyr)
require(reshape2)
require(jcolors)
# require(biomaRt)
##------------------------------------------------------------------------------------------------------------##

##------------------------------------------------------------------------------------------------------------##
##' 3' bias samples to remove
remove_samples = rbind(c('PDO14', 'PDO16', 'PDO18'), c('54288org', '119025org', '32070org'))
##------------------------------------------------------------------------------------------------------------##

##------------------------------------------------------------------------------------------------------------##
joint_counts_CN = readRDS("../output/joint_counts_CN_TPM_20210301210511.RDS")
joint_counts_CN = joint_counts_CN %>% filter( !(counts.Var1 %in% remove_samples[2,]))
topvariable = read.table("../../top_variable.txt", comment.char = "#")

## now add DESeq counts
deseq_counts = read.table("../../../RNASeq_DE_resistant_sensitive/files/counts_norm.csv",
                          sep=',', header = T)

deseq_counts = (melt(deseq_counts))
colnames(deseq_counts) = c('Gene', 'DESEq.sample', 'DESeq.count')

## Gene name conversion
# mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
# ensembl <- useEnsembl(biomart = "genes")
# ensembl <- useDataset(dataset = "hsapiens_gene_ensembl", mart = ensembl)
# 
# t2g <- getBM(
#   attributes = c('ensembl_gene_id', 'external_gene_name'),
#   values = deseq_counts$Gene,
#   filter = 'ensembl_gene_id',
#   mart = ensembl, useCache = FALSE)
# saveRDS(t2g, file = "~/Desktop/t2g.RDS")
t2g = readRDS("~/Desktop/t2g.RDS")

genename_ensmid_matched = t2g[match(deseq_counts$Gene, t2g$ensembl_gene_id),'external_gene_name']
deseq_counts$Gene = genename_ensmid_matched

renaming = readxl::read_xlsx("/Users/morril01/Documents/PhD/other_repos/b_tape/Vias_Brenton/RNASeq_DE_resistant_sensitive/files/PDOnameProperSample_sWGS_RNAseq.xlsx")

deseq_counts$DESEq.sample = renaming$ID[match(gsub("[.]", "-", deseq_counts$DESEq.sample), renaming$sampleNameRNAseq)]

joint_counts_CN = cbind.data.frame(joint_counts_CN, deseq_counts[match(paste0(joint_counts_CN$CN.L1, joint_counts_CN$CN.gene_name), paste0(deseq_counts$DESEq.sample, deseq_counts$Gene)),])
joint_counts_CN[4,]

# saveRDS(joint_counts_CN, file = "../output/joint_counts_CN_TPM_DESeq2.RDS")

# corDfAll = read.csv("../GexVsCnGwDeseq2/corStatTable.csv")
corDfAll = readRDS("../output/corDfAll.RDS")
joint_counts_CN_TPM_subset = readRDS("../output/joint_counts_CN_TPM_subset.RDS")

length(table(corDfAll$GeneName))
length(table(joint_counts_CN_TPM_subset$nearestGeneCN.gene)) ## ?? we're losing genes?
length(table(joint_counts_CN_TPM_DESeq2_subset$nearestGeneCN.gene)) ## ?? we're losing genes?
length(table(topvariable))

##------------------------------------------------------------------------------------------------------------##

##------------------------------------------------------------------------------------------------------------##
## General plots
ggplot(joint_counts_CN %>% filter(CN.L1 == '118947org' ), aes(x=CN.value, y=counts.value))+geom_point()

## per organoid separately in facets
ggplot(joint_counts_CN, aes(x=CN.value, y=counts.value))+geom_point()+facet_wrap(.~CN.L1, scales = "free")+theme_bw()
ggsave("../../figures/joint_counts_CN_TPM_all.pdf", width=10, height=10)

## same, log scale
ggplot(joint_counts_CN, aes(x=CN.value, y=counts.value))+geom_point()+facet_wrap(.~CN.L1, scales = "free")+theme_bw()+
  scale_x_continuous(trans = "log2")+scale_y_continuous(trans = "log2")
ggsave("../../figures/joint_counts_CN_TPM_all_loglog.pdf", width=10, height=8)

## only most variable genes (when it comes to neighbouring CN)
ggplot(joint_counts_CN %>% filter(CN.gene_name %in% topvariable$V1), aes(x=CN.value, y=counts.value, col=CN.gene_name))+
  geom_point()

pdf("~/Desktop/Vias/figures/joint_counts_CN_TPM_topvar.pdf")
for(i in topvariable$V1){
  print(ggplot(joint_counts_CN %>% filter(CN.gene_name == i), aes(x=CN.value, y=counts.value, label=CN.L1))+
          geom_point()+geom_smooth()+ggtitle(i)+geom_label_repel()+theme_bw())
}
dev.off()
##------------------------------------------------------------------------------------------------------------##

##------------------------------------------------------------------------------------------------------------##

ggplot(joint_counts_CN, aes(x=CN.value, y=DESeq.count))+geom_point()+facet_wrap(.~CN.L1, scales = "free")+theme_bw()
ggsave("../../figures/joint_counts_CN_DESeq_all.pdf", width=10, height=10)

ggplot(joint_counts_CN, aes(x=CN.value, y=DESeq.count))+geom_point()+facet_wrap(.~CN.L1, scales = "free")+theme_bw()+
  scale_x_continuous(trans = "log2")+scale_y_continuous(trans = "log2")
ggsave("../../figures/joint_counts_CN_DESeq_all_loglog.pdf", width=10, height=8)

pdf("../../figures/joint_counts_CN_DESeq_topvar.pdf")
for(i in topvariable$V1){
  print(ggplot(joint_counts_CN %>% filter(CN.gene_name == i), aes(x=CN.value, y=DESeq.count, label=CN.L1))+
          geom_point()+geom_smooth()+ggtitle(i)+geom_label_repel()+theme_bw())
}
dev.off()

## raw correlation between CN and deseq count
plot(joint_counts_CN$counts.value, joint_counts_CN$DESeq.count)

##------------------------------------------------------------------------------------------------------------##

## Looking at genes in specific
ggplot(joint_counts_CN %>% filter(CN.gene_name == "CCNE1"), aes(x=CN.value, y=counts.value, label=CN.L1))+
  geom_point()+geom_smooth()+ggtitle("CCNE1")+geom_label_repel()+theme_bw()


ggplot(joint_counts_CN %>% filter(CN.gene_name == "CDK12"), aes(x=CN.value, y=counts.value, label=CN.L1))+
  geom_point()+geom_smooth()+ggtitle("CDK12")+geom_label_repel()+theme_bw()


##------------------------------------------------------------------------------------------------------------##

## only genes included in Stephane's analysis
## length(table(joint_counts_CN_TPM_subset$nearestGeneCN.gene))
ggplot(joint_counts_CN %>% filter(CN.gene_name %in% joint_counts_CN_TPM_subset$nearestGeneCN.gene),
       aes(x=CN.value, y=DESeq.count))+geom_point()+facet_wrap(.~CN.L1, scales = "free")+theme_bw()+
  scale_x_continuous(trans = "log2")+scale_y_continuous(trans = "log2")+
  geom_abline(slope = 1, intercept = 0, col='blue', lty='dashed')

ggplot(joint_counts_CN %>% filter(DESeq.count > 0) %>%
         filter(CN.gene_name %in% joint_counts_CN_TPM_subset$nearestGeneCN.gene),
       aes(x=CN.value, y=DESeq.count))+geom_point()+theme_bw()+
  scale_x_continuous(trans = "log2")+scale_y_continuous(trans = "log2")+
  geom_abline(slope = 1, intercept = 0, col='blue', lty='dashed')+geom_smooth()


ggplot(joint_counts_CN %>% filter(DESeq.count > 0) %>%
         filter(CN.gene_name %in% joint_counts_CN_TPM_subset$nearestGeneCN.gene),
       aes(x=CN.value, y=DESeq.count))+geom_point()+theme_bw()+
  scale_x_continuous(trans = "log2")+scale_y_continuous(trans = "log2")+
  geom_abline(slope = 1, intercept = 0, col='blue', lty='dashed')+geom_smooth()

joint_counts_CN = joint_counts_CN %>% group_by(counts.Var2) %>% mutate(var_deseq=var(DESeq.count))

## only highly variable expressed genes
thres_var_GE = quantile(joint_counts_CN$var_deseq, na.rm = T, prob=0.9)
high_GEvar=joint_counts_CN %>% filter(DESeq.count > 0) %>%
  filter(var_deseq > thres_var_GE)
ggplot(high_GEvar,
       aes(x=CN.value, y=DESeq.count))+geom_point()+theme_bw()+
  scale_x_continuous(trans = "log2")+scale_y_continuous(trans = "log2")+
  geom_abline(slope = 1, intercept = 0, col='blue', lty='dashed')+geom_smooth()


## Link this to my weighted copy number for all genes

joint_counts_CN_TPM_DESeq2 = readRDS("../output/joint_counts_CN_TPM_DESeq2.RDS")
joint_counts_CN_TPM_DESeq2 = joint_counts_CN_TPM_DESeq2 %>% filter( !(counts.Var1 %in% remove_samples[2,]))

head(corDfAll)

corDfAll_df = (melt(cbind.data.frame(gene=rownames(corDfAll),
                                     corDfAll[,grepl('segVal_', colnames(corDfAll))])))
corDfAll_df$gene = sapply(corDfAll_df$gene, function(i) strsplit(i, '::')[[1]][2])
corDfAll_df$variable = gsub("segVal_", "", corDfAll_df$variable)
head(corDfAll_df)

joint_counts_CN_TPM_DESeq2_subset = cbind.data.frame(nearestGeneCN=corDfAll_df,
                                                     joint_counts_CN_TPM_DESeq2[match(paste0(corDfAll_df$gene, corDfAll_df$variable), paste0(joint_counts_CN_TPM_DESeq2$counts.Var2, joint_counts_CN_TPM_DESeq2$counts.Var1)),])

# saveRDS(joint_counts_CN_TPM_DESeq2_subset, "../output/joint_counts_CN_TPM_DESeq2_subset.RDS")
dim(joint_counts_CN_TPM_DESeq2_subset)
dim(joint_counts_CN_TPM_DESeq2_subset)

plot(joint_counts_CN_TPM_DESeq2_subset$CN.value,
     joint_counts_CN_TPM_DESeq2_subset$nearestGeneCN.value)

ggplot(joint_counts_CN_TPM_DESeq2_subset, aes(x=CN.value, y=nearestGeneCN.value))+geom_point(alpha=0.2)+
  labs(x='Weighted CN value of gene (Lena)', y='CN value of gene highly variable region nearby')+
  scale_x_continuous(trans = "log2")+scale_y_continuous(trans = "log2")
ggsave("../../figures/scatterplot_stephane_lena.pdf", width = 7, height = 5)

## removing one outlier
joint_counts_CN_TPM_DESeq2_subset = joint_counts_CN_TPM_DESeq2_subset[-which.max(joint_counts_CN_TPM_DESeq2_subset$counts.value),]

ggplot(joint_counts_CN_TPM_DESeq2_subset, aes(x=CN.value, y=counts.value))+geom_point()+
  scale_y_continuous(trans='log2')+geom_smooth()
ggplot(joint_counts_CN_TPM_DESeq2_subset %>% filter(counts.value > 0), aes(x=CN.value, y=counts.value))+geom_point()+
  scale_y_continuous(trans='log2')+geom_smooth()+lims(x=c(0,10))

## all genes
ggplot(joint_counts_CN_TPM_DESeq2, aes(x=CN.value, y=counts.value))+geom_point()+
  scale_y_continuous(trans='log2')+geom_smooth()

ggplot(joint_counts_CN_TPM_DESeq2, aes(x=CN.value, y=counts.value))+geom_point(alpha=0.2)+
  scale_y_continuous(trans='log2')+geom_smooth(col='red')+lims(x=c(0,20))+theme_bw()+
  labs(x='Weighted CN in gene', y='TPM (log2)')
ggsave("../../figures/scatterplot_CNweighted_TPM_all.pdf", width=8)


ggplot(joint_counts_CN_TPM_DESeq2, aes(x=CN.value, y=DESeq.count))+geom_point(alpha=0.2)+
  scale_y_continuous(trans='log2')+geom_smooth(col='red')+lims(x=c(0,20))+theme_bw()+
  labs(x='Weighted CN in gene', y='DEseq2 counts (log2)')
ggsave("../../figures/scatterplot_CNweighted_DEseq2counts_all.pdf", width=8)

## deseq counts and tmp seem to be quite similar
ggplot(joint_counts_CN_TPM_DESeq2, aes(x=counts.value, y=DESeq.count))+geom_point(alpha=0.2)

## genes differentially abundant between normal tissue and cancer tissue
## colour these genes
load("../RnaSeqPip/DEAnalysis/deObject_SampleGroup.RData")
results_deseq = DESeq2::results(SampleGroup)

top_genes_normal_tumour = rownames(results_deseq)[order(results_deseq$padj, decreasing = F)[1:400]]
t2g  = readRDS("~/Desktop/t2g.RDS")
top_genes_normal_tumour = cbind.data.frame(ensembl=top_genes_normal_tumour, gene=t2g$external_gene_name[match(top_genes_normal_tumour, t2g$ensembl_gene_id)])

joint_counts_CN_TPM_DESeq2$topDiff_norm_tissue = ifelse(joint_counts_CN_TPM_DESeq2$CN.gene_name %in% top_genes_normal_tumour$gene,
                                                        yes = 'DE', no = 'nonDE' )

ggplot(joint_counts_CN_TPM_DESeq2, aes(x=CN.value, y=DESeq.count, col=topDiff_norm_tissue))+geom_point(alpha=0.2)+
  scale_y_continuous(trans='log2')+geom_smooth(col='red')+lims(x=c(0,20))+theme_bw()+
  labs(x='Weighted CN in gene', y='DEseq2 counts (log2)')+theme(legend.position = "bottom")
ggsave("../../figures/scatterplot_CNweighted_DEseq2counts_all_colourDE.pdf", width=4, height=4)

ggplot(joint_counts_CN_TPM_DESeq2[joint_counts_CN_TPM_DESeq2$topDiff_norm_tissue == 'DE',], aes(x=CN.value, y=DESeq.count, col=topDiff_norm_tissue))+geom_point(alpha=0.2)+
  scale_y_continuous(trans='log2')+geom_smooth(col='red')+lims(x=c(0,20))+theme_bw()+
  labs(x='Weighted CN in gene', y='DEseq2 counts (log2)')+theme(legend.position = "bottom")
ggsave("../../figures/scatterplot_CNweighted_DEseq2counts_onlyDE.pdf", width=4, height=4)

ggplot(joint_counts_CN_TPM_DESeq2[joint_counts_CN_TPM_DESeq2$counts.value > 0,],
       aes(x=CN.value, y=counts.value))+geom_point(alpha=0.2)+
  scale_y_continuous(trans='log2')+geom_smooth(col='red')+lims(x=c(0,20))+theme_bw()+
  labs(x='Weighted CN in gene', y='TPM (log2)')+theme(legend.position = "bottom")
ggsave("../../figures/scatterplot_CNweighted_TPM_nonzero.pdf", width=8)

ggplot(joint_counts_CN_TPM_DESeq2_subset[joint_counts_CN_TPM_DESeq2_subset$counts.value > 0,],
       aes(x=nearestGeneCN.value, y=counts.value))+geom_point(alpha=0.2)+
  scale_y_continuous(trans='log2')+geom_smooth(col='red')+lims(x=c(0,20))+theme_bw()+
  labs(x='Weighted CN in gene', y='TPM (log2)')+theme(legend.position = "bottom")
ggsave("../../figures/scatterplot_CNclosest_TPM_nonzero.pdf", width=8)

## coding
ensembl <- useMart("ensembl",
                   dataset = "hsapiens_gene_ensembl")
coding_genes = getBM(attributes=c("ensembl_gene_id","external_gene_name","description"),
                     filters='biotype', values=c('protein_coding'), mart=ensembl)
coding_genes$external_gene_name

length(coding_genes$external_gene_name)
length(unique(joint_counts_CN_TPM_DESeq2$CN.gene_name))
dim(joint_counts_CN_TPM_DESeq2)
dim(joint_counts_CN_TPM_DESeq2[joint_counts_CN_TPM_DESeq2$CN.gene_name %in% coding_genes$external_gene_name,])
## they already are coding, the vast majority

## for each gene, find the rank of both ploidy and RNASeq
head(joint_counts_CN_TPM_DESeq2)

## group by gene
joint_counts_CN_TPM_DESeq2_ranked = joint_counts_CN_TPM_DESeq2 %>% group_by(CN.gene_name) %>%
  mutate(rank_weighted_CN = order(CN.value), rank_DESeq = order(DESeq.count), rank_TPM = order(counts.value))
joint_counts_CN_TPM_DESeq2_normalised = joint_counts_CN_TPM_DESeq2 %>% group_by(CN.gene_name) %>%
  mutate(scaled_centered_weighted_CN = scale(CN.value), scaled_centered_DESeq = scale(DESeq.count),
         scaled_centered_TPM = scale(counts.value))
joint_counts_CN_TPM_DESeq2_normalised_excludingnormal = joint_counts_CN_TPM_DESeq2 %>% filter(CN.value != 2) %>%
  group_by(CN.gene_name) %>%
  mutate(scaled_centered_weighted_CN = scale(CN.value), scaled_centered_DESeq = scale(DESeq.count),
         scaled_centered_TPM = scale(counts.value))

example_centering = joint_counts_CN_TPM_DESeq2_normalised %>% filter(counts.Var2 == 'NOC2L') %>% dplyr::select(scaled_centered_weighted_CN)
mean(example_centering$scaled_centered_weighted_CN); sd(example_centering$scaled_centered_weighted_CN)
table(joint_counts_CN_TPM_DESeq2_ranked$rank_weighted_CN)

## there are genes which are duplicated - remove them
joint_counts_CN_TPM_DESeq2_ranked = joint_counts_CN_TPM_DESeq2_ranked[!(joint_counts_CN_TPM_DESeq2_ranked$CN.gene_name %in% joint_counts_CN_TPM_DESeq2_ranked$CN.gene_name[joint_counts_CN_TPM_DESeq2_ranked$rank_weighted_CN > length(unique(joint_counts_CN_TPM_DESeq2_ranked$counts.Var1))]),]
table(joint_counts_CN_TPM_DESeq2_ranked$rank_weighted_CN)

plot(joint_counts_CN_TPM_DESeq2_ranked$rank_weighted_CN, joint_counts_CN_TPM_DESeq2_ranked$rank_DESeq)

ggplot(joint_counts_CN_TPM_DESeq2_ranked, aes(rank_weighted_CN, rank_DESeq)) + geom_point()
ggplot(joint_counts_CN_TPM_DESeq2_ranked, aes(x = rank_weighted_CN, y = rank_DESeq)) + geom_tile()

ggplot((melt(table(joint_counts_CN_TPM_DESeq2_ranked$rank_weighted_CN, joint_counts_CN_TPM_DESeq2_ranked$rank_DESeq)/nrow(joint_counts_CN_TPM_DESeq2_ranked))),
       aes(x=Var1, y=Var2, fill=value))+geom_raster()+scale_fill_jcolors_contin("pal2", bias = 1.75) + theme_bw()+
  labs(x='Organoid rank for CN (weighted)', y='Organoid rank for gene expression (DESeq2 counts)')
ggsave("../../figures/rankplot_CNweighted_deseq.pdf", width=8)

df <- expand.grid(x = sample(x = 1:100, size = 50), y = sample(1:100, size = 50))
# default is compatible with geom_tile()
ggplot(df, aes(x, y)) + geom_raster()

## blob, but clear trend
ggplot(joint_counts_CN_TPM_DESeq2_normalised %>% filter(CN.value != 2), aes(x=scaled_centered_weighted_CN, y=scaled_centered_DESeq))+
  geom_point(alpha=0.01)+geom_smooth()

ggplot(joint_counts_CN_TPM_DESeq2_normalised %>% filter(CN.value != 2), aes(x=scaled_centered_weighted_CN, y=scaled_centered_DESeq))+
  geom_point(alpha=0.01)+geom_smooth()+facet_wrap(.~(CN.value < 2), scales = "free")

ggplot(joint_counts_CN_TPM_DESeq2_normalised %>% filter(CN.value > 2), aes(x=scaled_centered_weighted_CN, y=scaled_centered_DESeq))+
  geom_point(alpha=0.01)+geom_smooth()
ggsave("../../figures/scatterplot_normCNweighted_normDESeq.pdf", width=8)

ggplot(joint_counts_CN_TPM_DESeq2_normalised_excludingnormal %>% filter(CN.value > 2), aes(x=scaled_centered_weighted_CN, y=scaled_centered_DESeq))+
  geom_point(alpha=0.01)+geom_smooth()+facet_wrap(.~counts.Var1)
ggsave("../../figures/scatterplot_normCNweighted_normDESeq_perorg.pdf", width=8)

## same but normalising only with nonnormal segments
ggplot(joint_counts_CN_TPM_DESeq2_normalised_excludingnormal, aes(x=scaled_centered_weighted_CN, y=scaled_centered_DESeq))+
  geom_point(alpha=0.01)+geom_smooth()
ggsave("../../figures/scatterplot_normCNweighted_normDESeq_nonormal.pdf", width=8)

ggplot(joint_counts_CN_TPM_DESeq2_normalised_excludingnormal, aes(x=scaled_centered_weighted_CN, y=scaled_centered_DESeq))+
  geom_point(alpha=0.01)+geom_smooth()+facet_wrap(.~(CN.value < 2), scales = "free")

## only nonnormal segments, and only with Stephane's genes
ggplot(joint_counts_CN_TPM_DESeq2_normalised_excludingnormal %>% filter(counts.Var2 %in% unlist(topvariable)),
       aes(x=scaled_centered_weighted_CN, y=scaled_centered_DESeq, col=CN.value < 2))+
  geom_point(alpha=0.1)+geom_smooth()+facet_wrap(.~(CN.value < 2), scales = "free")
ggsave("../../figures/scatterplot_normCNweighted_normDESeq_nonormal_highlyvarregions.pdf", width=8)

## not good
ggplot(joint_counts_CN_TPM_DESeq2_normalised %>% filter(CN.value != 2), aes(x=scaled_centered_weighted_CN, y=DESeq.count))+
  geom_point(alpha=0.2)+geom_smooth()+scale_y_continuous(trans = "log2")

## not great
ggplot(joint_counts_CN_TPM_DESeq2_normalised %>% filter(CN.value != 2), aes(x=CN.value, y=scaled_centered_DESeq))+
  geom_point(alpha=0.2)+geom_smooth()+lims(x=c(0,10))

## good trend
ggplot(joint_counts_CN_TPM_DESeq2_normalised  %>% filter(CN.value != 2), aes(x=scaled_centered_weighted_CN, y=scaled_centered_TPM))+
  geom_point(alpha=0.2)+geom_smooth()

## not great, but good trend
ggplot(joint_counts_CN_TPM_DESeq2_normalised_excludingnormal, aes(x=scaled_centered_weighted_CN, y=scaled_centered_TPM))+
  geom_point(alpha=0.2)+geom_smooth()

