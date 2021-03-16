##------------------------------------------------------------------------------------------------------------##
rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

require(ggplot2)
require(ggrepel)
require(dplyr)
require(reshape2)
require(jcolors)
require(biomaRt)
require(gridExtra)
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

# saveRDS(joint_counts_CN, file = "../output/joint_counts_CN.RDS")

# corDfAll = read.csv("../GexVsCnGwDeseq2/corStatTable.csv")
corDfAll = readRDS("../output/corDfAll.RDS")
# joint_counts_CN_TPM_subset = readRDS("../output/joint_counts_CN_TPM_subset.RDS")

## add TPM
joint_counts_CN_TPM_subset = joint_counts_CN %>% filter( !(counts.Var1 %in% remove_samples[2,]))
joint_counts_CN_TPM_subset = cbind(joint_counts_CN_TPM_subset,
                                   counts=joint_counts_CN[match(paste0(joint_counts_CN_TPM_subset$nearestGeneCN.variable, joint_counts_CN_TPM_subset$counts.Var2),
      paste0(joint_counts_CN$counts.Var1, joint_counts_CN$counts.Var2)),])

plot(joint_counts_CN_TPM_subset$counts.value, joint_counts_CN_TPM_subset$counts.DESeq.count) ## we have added DESeq2
plot(joint_counts_CN_TPM_subset$counts.value, joint_counts_CN_TPM_subset$counts.counts.value)

length(table(corDfAll$GeneName))
length(table(joint_counts_CN_TPM_subset$nearestGeneCN.gene)) ## ?? we're losing genes?
# length(table(joint_counts_CN_subset$nearestGeneCN.gene)) ## ?? we're losing genes?
length(table(topvariable))

## only highly variable expressed genes
joint_counts_CN = joint_counts_CN %>% group_by(counts.Var2) %>% mutate(var_deseq=var(DESeq.count))
thres_var_GE = quantile(joint_counts_CN$var_deseq, na.rm = T, prob=0.9)
high_GEvar=joint_counts_CN %>% filter(DESeq.count > 0) %>%
  filter(var_deseq > thres_var_GE)

## Link this to my weighted copy number for all genes
# joint_counts_CN = readRDS("../output/joint_counts_CN.RDS")
# joint_counts_CN = joint_counts_CN %>% filter( !(counts.Var1 %in% remove_samples[2,]))

dim(joint_counts_CN)

# x = joint_counts_CN[match(paste0(joint_counts_CN$counts.Var1, joint_counts_CN$counts.Var2),
#       paste0(joint_counts_CN$counts.Var1, joint_counts_CN$counts.Var2)),]
# plot(x$counts.value, joint_counts_CN$counts.value)
# 
# head(joint_counts_CN[is.na(x$counts.value),])
# head(x[is.na(x$counts.value),])

head(corDfAll)

corDfAll_df = (melt(cbind.data.frame(gene=rownames(corDfAll),
                                     corDfAll[,grepl('segVal_', colnames(corDfAll))])))
corDfAll_df$gene = sapply(corDfAll_df$gene, function(i) strsplit(i, '::')[[1]][2])
corDfAll_df$variable = gsub("segVal_", "", corDfAll_df$variable)
head(corDfAll_df)

joint_counts_CN_subset = cbind.data.frame(nearestGeneCN=corDfAll_df,
   joint_counts_CN[match(paste0(corDfAll_df$gene, corDfAll_df$variable),
                         paste0(joint_counts_CN$counts.Var2, joint_counts_CN$counts.Var1)),])

# saveRDS(joint_counts_CN_subset, "../output/joint_counts_CN_subset.RDS")


## genes differentially abundant between normal tissue and cancer tissue
## colour these genes
load("../RnaSeqPip/DEAnalysis/deObject_SampleGroup.RData")
results_deseq = DESeq2::results(SampleGroup)

top_genes_normal_tumour = rownames(results_deseq)[order(results_deseq$padj, decreasing = F)[1:400]]
top_genes_normal_tumour = cbind.data.frame(ensembl=top_genes_normal_tumour, gene=t2g$external_gene_name[match(top_genes_normal_tumour, t2g$ensembl_gene_id)])

joint_counts_CN$topDiff_norm_tissue = ifelse(joint_counts_CN$CN.gene_name %in% top_genes_normal_tumour$gene,
                                             yes = 'DE', no = 'nonDE' )
## coding
# ensembl <- useMart("ensembl",
#                    dataset = "hsapiens_gene_ensembl")
# coding_genes = getBM(attributes=c("ensembl_gene_id","external_gene_name","description"),
#                      filters='biotype', values=c('protein_coding'), mart=ensembl)
# saveRDS(coding_genes, "~/Desktop/coding_genes.RDS")
coding_genes.RDS <- readRDS("~/Desktop/coding_genes.RDS")
coding_genes$external_gene_name

## Removing outliers from joint_counts_CN_subset
joint_counts_CN_subset = joint_counts_CN_subset[!is.na(joint_counts_CN_subset$nearestGeneCN.value),]

probs_outlier=0.999
plot(sort(log(joint_counts_CN_subset$nearestGeneCN.value)))
joint_counts_CN_subset$flag_outlier_nearestGeneCN.value = ifelse(joint_counts_CN_subset$nearestGeneCN.value > quantile(joint_counts_CN_subset$nearestGeneCN.value, probs = probs_outlier),
                                                                yes = TRUE, no = FALSE)

outliers = joint_counts_CN_subset[joint_counts_CN_subset$flag_outlier_nearestGeneCN.value,]
table(outliers$nearestGeneCN.gene)
table(outliers$counts.Var1)

pdf("../../figures/joint_counts_CN_subset_remove_outliers.pdf")
plot(sort(log(joint_counts_CN_subset$nearestGeneCN.value)),
     col=as.factor(joint_counts_CN_subset$flag_outlier_nearestGeneCN.value[order(log(joint_counts_CN_subset$nearestGeneCN.value))]),
     main=paste0('Removing ', probs_outlier, ' percentile'))
dev.off()

subset_genes_of_interest = c('MYC', 'CCNE1', 'PIK3CA', 'TERT', 'KRAS', 'PTEN', 'RB1', 'AKT1',
                             'AKT2', 'PARP1', 'PARP2', 'ATM', 'ATR', 'WEE1', 'TOP1', 'TUBB1',
                             'AKT3', 'CCND1', 'CCND2', 'CCND3', 'CDKN2A', 'CDKN2B', 'MECOM', 'CDK12')

joint_counts_CN_subset = joint_counts_CN_subset[!joint_counts_CN_subset$flag_outlier_nearestGeneCN.value,]
# plot(sort(log(joint_counts_CN$CN.value)))
joint_counts_CN_subset = joint_counts_CN_subset %>%
  group_by(CN.gene_name) %>%
  mutate(scaled_centered_weighted_CN = exp(scale(log(CN.value))),
         scaled_centered_weighted_CN_mean = (scale((CN.value))),
         scaled_centered_DESeq = scale(DESeq.count),
         scaled_centered_TPM = scale(counts.value),
         scaled_centered_nearestCN = exp(scale(log(nearestGeneCN.value))),
         scaled_centered_nearestCN_mean = (scale((nearestGeneCN.value)))
  )
joint_counts_CN_subset = joint_counts_CN_subset[!(joint_counts_CN_subset$nearestGeneCN.variable %in% remove_samples[2,]),]
joint_counts_CN_subset$labels = joint_counts_CN_subset$CN.gene_name;  joint_counts_CN_subset$labels[!(joint_counts_CN_subset$labels %in% subset_genes_of_interest)] = NA

joint_counts_CN = joint_counts_CN %>% group_by(CN.gene_name) %>%
                                     mutate(rank_weighted_CN = order(CN.value), rank_DESeq = order(DESeq.count), rank_TPM = order(counts.value))
joint_counts_CN = joint_counts_CN %>% group_by(CN.gene_name) %>%
                                     mutate(scaled_centered_weighted_CN = scale(CN.value),
                                            scaled_centered_DESeq = scale(DESeq.count),
                                            scaled_centered_TPM = scale(counts.value))
joint_counts_CN_normalised_excludingnormal = joint_counts_CN %>% filter(CN.value != 2) %>%
  group_by(CN.gene_name) %>%
  mutate(scaled_centered_weighted_CN = scale(CN.value),
         scaled_centered_DESeq = scale(DESeq.count),
         scaled_centered_TPM = scale(counts.value))

joint_counts_CN_normalised_highlyvar = joint_counts_CN_subset %>%
  group_by(CN.gene_name) %>%
  mutate(scaled_neighbour_CN=scale(nearestGeneCN.value),
         scaled_centered_weighted_CN = scale(CN.value),
         scaled_centered_DESeq = scale(DESeq.count),
         scaled_centered_TPM = scale(counts.value))
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

pdf("../../figures/joint_counts_CN_TPM_topvar.pdf")
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

ggplot(high_GEvar,
       aes(x=CN.value, y=DESeq.count))+geom_point()+theme_bw()+
  scale_x_continuous(trans = "log2")+scale_y_continuous(trans = "log2")+
  geom_abline(slope = 1, intercept = 0, col='blue', lty='dashed')+geom_smooth()

ggplot(joint_counts_CN_subset, aes(x=CN.value, y=nearestGeneCN.value,
       col=nearestGeneCN.gene %in% subset_genes_of_interest))+geom_point(alpha=0.2)+
  labs(x='Weighted CN value of gene (Lena)', y='CN value of gene highly variable region nearby')+
  facet_wrap(.~(nearestGeneCN.gene %in% subset_genes_of_interest))+
  scale_x_continuous(trans = "log2")+scale_y_continuous(trans = "log2")+guides(col=F)
ggsave("../../figures/scatterplot_stephane_lena.pdf", width = 7, height = 5)

## CN vs count for genes in highly variable areas
ggplot(joint_counts_CN_subset, aes(x=CN.value, y=counts.value))+geom_point()+
  scale_y_continuous(trans='log2')+geom_smooth()
## nearest segment vs count
ggplot(joint_counts_CN_subset %>% filter(counts.value > 0), aes(x=nearestGeneCN.value, y=counts.value))+geom_point()+
  scale_y_continuous(trans='log2')+geom_smooth()+lims(x=c(0,10))
## nearest segment vs count, excluding normal
ggplot()+
  geom_point(data = joint_counts_CN_subset %>% filter(counts.value > 1, abs(nearestGeneCN.value - 2) > 0.1) %>%
               dplyr::select(-labels),
             aes(x=nearestGeneCN.value, y=counts.value), alpha=0.2)+
  geom_point(data = joint_counts_CN_subset %>% filter(counts.value > 0, abs(nearestGeneCN.value - 2) > 0.1) %>% filter(!is.na(labels)),
             aes(x=nearestGeneCN.value, y=counts.value, label=labels, col=labels))+
  scale_y_continuous(trans='log2')+geom_smooth()+lims(x=c(0,10))+geom_label_repel()+facet_wrap(.~labels)+
  guides(col=FALSE)+theme_bw()
ggsave("../../figures/selected_genes_nearestgene_TPM_nonorm.png", width = 7, height = 5)

ggplot()+
  geom_point(data = joint_counts_CN_subset %>% filter(counts.value > 1, abs(nearestGeneCN.value - 2) > 0.1) %>%
               dplyr::select(-labels),
             aes(x=CN.value, y=counts.value), alpha=0.2)+
  geom_point(data = joint_counts_CN_subset %>% filter(counts.value > 0, abs(CN.value - 2) > 0.1) %>% filter(!is.na(labels)),
             aes(x=CN.value, y=counts.value, label=labels, col=labels))+
  scale_y_continuous(trans='log2')+geom_smooth()+lims(x=c(0,10))+geom_label_repel()+facet_wrap(.~labels)+
  guides(col=FALSE)+theme_bw()
ggsave("../../figures/selected_genes_weightedCN_TPM_nonorm.png", width = 7, height = 5)

ggplot()+
  geom_point(data = joint_counts_CN_subset %>% filter(counts.value > 0, abs(nearestGeneCN.value - 2) > 0.1) %>%
               dplyr::select(-labels),
             aes(x=scaled_centered_nearestCN, y=scaled_centered_TPM), alpha=0.2)+
  geom_point(data = joint_counts_CN_subset %>% filter(counts.value > 0, abs(nearestGeneCN.value - 2) > 0.1) %>% filter(!is.na(labels)),
             aes(x=scaled_centered_nearestCN, y=scaled_centered_TPM, label=labels, col=labels))+
  geom_smooth()+geom_label_repel()+facet_wrap(.~labels)+
  guides(col=FALSE)+theme_bw()

ggplot()+
  geom_point(data = joint_counts_CN_subset %>% filter(counts.value > 1, abs(nearestGeneCN.value - 2) > 0.1) %>%
               dplyr::select(-labels),
             aes(x=scaled_centered_weighted_CN, y=scaled_centered_DESeq), alpha=0.2)+
  geom_point(data = joint_counts_CN_subset %>% filter(counts.value > 0, abs(nearestGeneCN.value - 2) > 0.1) %>% filter(!is.na(labels)),
             aes(x=scaled_centered_weighted_CN, y=scaled_centered_DESeq, label=labels, col=labels))+
  geom_smooth()+geom_label_repel()+facet_wrap(.~labels)+
  guides(col=FALSE)+theme_bw()+scale_x_continuous(trans="log2")
ggsave("../../figures/selected_genes_scaled_weightedCN_scaledDESeq_nonorm.png", width = 7, height = 5)

## it doesn't look centered to me? e.g. MYC doesn't seem to be centered for CN
# (joint_counts_CN_subset %>% filter(nearestGeneCN.gene == 'MYC') %>% dplyr::select(scaled_centered_weighted_CN))[,2] %>% unlist  %>% mean
# (joint_counts_CN_subset %>% filter(nearestGeneCN.gene == 'MYC') %>% dplyr::select(scaled_centered_weighted_CN))[,2] %>% unlist  %>% sd
## it is. there is just a bunch of observations in the same CN (this is what happens if you don't remove outliers)

ggplot(joint_counts_CN_subset %>% filter(counts.value > 0), aes(x=CN.value, y=counts.value))+geom_point()+
  scale_y_continuous(trans='log2')+geom_smooth()+lims(x=c(0,10))

ggplot(data = joint_counts_CN_subset %>% filter(counts.value > 0, abs(nearestGeneCN.value - 2) > 0.1) %>%
               filter(nearestGeneCN.gene == 'MYC') %>%
               dplyr::select(-labels),
             aes(x=scaled_centered_weighted_CN, y=scaled_centered_DESeq, col=log(CN.value)))+
  geom_point()
ggplot(data = joint_counts_CN_subset %>% filter(nearestGeneCN.gene == 'MYC') %>%
         dplyr::select(-labels),
       aes(x=nearestGeneCN.value, y=counts.value))+
  geom_point()+ggtitle('MYC CN and GE')


ggplot(data = joint_counts_CN_subset,
       aes(x=scaled_centered_weighted_CN, y=scaled_centered_DESeq, label=labels, col=labels))+
  geom_point(alpha=1)+scale_x_continuous(trans="log2")+geom_smooth()+
  # geom_label()+
  facet_wrap(.~labels, scales="free", nrow=2)

##------------------------------------------------------------------------------------------------------------##

corr_centered =  joint_counts_CN[,!duplicated(colnames(joint_counts_CN))] %>%
  group_by(CN.gene_name) %>% summarise(cor(scaled_centered_weighted_CN, scaled_centered_DESeq))
r2_centered =  joint_counts_CN %>% group_by(CN.gene_name) %>% summarise( if(all(is.na(scaled_centered_DESeq)) | all(is.na(scaled_centered_weighted_CN))){NA}else{ try(summary(lm( scaled_centered_DESeq ~ scaled_centered_weighted_CN))$r.squared)} )
slope_centered =  joint_counts_CN %>% group_by(CN.gene_name) %>% summarise( if(all(is.na(scaled_centered_DESeq)) | all(is.na(scaled_centered_weighted_CN))){NA}else{ try(summary(lm( scaled_centered_DESeq ~ scaled_centered_weighted_CN))$coefficients[2])} )
var_CN =  joint_counts_CN %>% group_by(CN.gene_name) %>% summarise( var(CN.value) )
var_DESeq2 =  joint_counts_CN %>% group_by(CN.gene_name) %>% summarise( var(DESeq.count) )
colnames(r2_centered) = c('Gene', 'r2')
top_cor = corr_centered$CN.gene_name[order(corr_centered$`cor(scaled_centered_weighted_CN, scaled_centered_DESeq)`, decreasing = T)[1:30]]
top_r2 = r2_centered$Gene[order(r2_centered$r2, decreasing = T)[1:30]]

df_gene_characteristics = cbind.data.frame(cor_normCNnormDESeq=corr_centered$`cor(scaled_centered_weighted_CN, scaled_centered_DESeq)`,
                                           var_DESeq=var_DESeq2$`var(DESeq.count)`,
                                           r2_normCNnormDESeq=r2_centered$r2,
                                           var_CN=var_CN$`var(CN.value)`,
                                           slope_centered=slope_centered,
                                           Gene=as.character(corr_centered$CN.gene_name))

##------------------------------------------------------------------------------------------------------------##


ggplot(data = joint_counts_CN_subset %>% filter(CN.gene_name %in% top_cor),
       aes(x=scaled_centered_weighted_CN, y=scaled_centered_DESeq, col=CN.gene_name))+
  geom_point(alpha=1)+scale_x_continuous(trans="log2")+geom_smooth()+
  facet_wrap(.~CN.gene_name, scales="free", nrow=2)
ggplot(data = joint_counts_CN_subset %>% filter(CN.gene_name %in% top_r2),
       aes(x=scaled_centered_weighted_CN, y=scaled_centered_DESeq, col=CN.gene_name))+
  geom_point(alpha=1)+scale_x_continuous(trans="log2")+geom_smooth()+
  facet_wrap(.~CN.gene_name, scales="free", nrow=2)

plot(corr_centered$`cor(scaled_centered_weighted_CN, scaled_centered_DESeq)`,
     r2_centered$r2)
plot(corr_centered$`cor(scaled_centered_weighted_CN, scaled_centered_DESeq)`,
     (var_CN$`var(CN.value)`))
plot(corr_centered$`cor(scaled_centered_weighted_CN, scaled_centered_DESeq)`,
     log(var_DESeq2$`var(DESeq.count)`))

ggplot(df_gene_characteristics, aes(x=cor_normCNnormDESeq, y=var_DESeq,
                                    label=as.character(ifelse(cor_normCNnormDESeq>0.8, Gene, NA))))+
  geom_point()+geom_label_repel()+
  scale_y_continuous(trans = "log2")

ggplot(data = joint_counts_CN_subset,
       aes(x=scaled_centered_weighted_CN, y=scaled_centered_DESeq, label=labels, col=labels))+
  geom_point(alpha=1)+scale_x_continuous(trans="log2")+geom_smooth()+
  # geom_label()+
  facet_wrap(.~labels, scales="free", nrow=2)
ggsave("../../figures/selected_genes_scaled_weightedCN_scaledDESeq_2.png", width = 7, height = 5)

ggplot(data = joint_counts_CN_subset,
       aes(x=scaled_centered_weighted_CN_mean, y=scaled_centered_DESeq, label=labels, col=labels))+
  geom_point(alpha=1)+geom_smooth()+
  # geom_label()+
  facet_wrap(.~labels, scales="free", nrow=2)


ggplot(data = joint_counts_CN %>% filter(counts.Var2 %in% subset_genes_of_interest),
       aes(x=scaled_centered_weighted_CN, y=scaled_centered_DESeq, label=counts.Var2))+
  geom_point(alpha=1)+scale_x_continuous(trans="log2")+geom_smooth()+
  # geom_label()+
  facet_wrap(.~counts.Var2, scales="free", nrow=2)

ggplot(data = joint_counts_CN_subset %>% filter(nearestGeneCN.gene == 'PARP1') %>%
         dplyr::select(-labels),
       aes(x=nearestGeneCN.value, y=scaled_centered_DESeq))+
  geom_point()
ggplot(data = joint_counts_CN_subset %>% filter(nearestGeneCN.gene == 'PARP1'),
       aes(x=scaled_centered_weighted_CN, y=scaled_centered_DESeq, label=labels, col=labels))+
  geom_point(alpha=1)+scale_x_continuous(trans="log2")+geom_smooth()+
  # geom_label()+
  facet_wrap(.~labels, scales="free")


# ggplot(data = joint_counts_CN_subset %>% filter(nearestGeneCN.gene %in% subset_genes_of_interest) %>%
#          dplyr::select(-labels),
#        aes(x=nearestGeneCN.value, y=counts.value, nearestGeneCN.gene=nearestGeneCN.gene, label=nearestGeneCN.variable))+
#   facet_wrap(.~nearestGeneCN.gene, scales = "free")+
#   geom_text_repel(size=2.5)+
#   geom_point()+ggtitle('genes of interest CN and GE')
# ggsave("../../figures/outliers_annotated_orgs.png", width = 7, height = 5)
# 
# logCN = joint_counts_CN_subset %>% group_by(nearestGeneCN.variable) %>% summarise(mean=mean(log(nearestGeneCN.value)))
# ggplot(joint_counts_CN_subset, aes(x=log(nearestGeneCN.value)))+geom_density()+
#   facet_wrap(.~factor(nearestGeneCN.variable, levels=logCN$nearestGeneCN.variable[order(logCN$mean)]), nrow=2)+
#   geom_vline(xintercept = mean(log(joint_counts_CN_subset$nearestGeneCN.value)), col='red')
# ggsave("../../figures/mean_CN_segment_per_org.pdf", width = 15, height = 5)

joint_counts_CN_subset %>% filter(nearestGeneCN.gene == 'MYC')

#-----------------------------------------------------------------------------------------#
## all genes
ggplot(joint_counts_CN, aes(x=CN.value, y=counts.value))+geom_point()+
  scale_y_continuous(trans='log2')+geom_smooth()

ggplot(joint_counts_CN, aes(x=CN.value, y=counts.value))+geom_point(alpha=0.2)+
  scale_y_continuous(trans='log2')+geom_smooth(col='red')+lims(x=c(0,20))+theme_bw()+
  labs(x='Weighted CN in gene', y='TPM (log2)')
ggsave("../../figures/scatterplot_CNweighted_TPM_all.pdf", width=8)


ggplot(joint_counts_CN, aes(x=CN.value, y=DESeq.count))+geom_point(alpha=0.2)+
  scale_y_continuous(trans='log2')+geom_smooth(col='red')+lims(x=c(0,20))+theme_bw()+
  labs(x='Weighted CN in gene', y='DEseq2 counts (log2)')
ggsave("../../figures/scatterplot_CNweighted_DEseq2counts_all.pdf", width=8)

## deseq counts and tmp seem to be quite similar
ggplot(joint_counts_CN, aes(x=counts.value, y=DESeq.count))+geom_point(alpha=0.2)


ggplot(joint_counts_CN, aes(x=CN.value, y=DESeq.count, col=topDiff_norm_tissue))+geom_point(alpha=0.2)+
  scale_y_continuous(trans='log2')+geom_smooth(col='red')+lims(x=c(0,20))+theme_bw()+
  labs(x='Weighted CN in gene', y='DEseq2 counts (log2)')+theme(legend.position = "bottom")
ggsave("../../figures/scatterplot_CNweighted_DEseq2counts_all_colourDE.pdf", width=4, height=4)

ggplot(joint_counts_CN[joint_counts_CN$topDiff_norm_tissue == 'DE',], aes(x=CN.value, y=DESeq.count, col=topDiff_norm_tissue))+geom_point(alpha=0.2)+
  scale_y_continuous(trans='log2')+geom_smooth(col='red')+lims(x=c(0,20))+theme_bw()+
  labs(x='Weighted CN in gene', y='DEseq2 counts (log2)')+theme(legend.position = "bottom")
ggsave("../../figures/scatterplot_CNweighted_DEseq2counts_onlyDE.pdf", width=4, height=4)

ggplot(joint_counts_CN[joint_counts_CN$counts.value > 0,],
       aes(x=CN.value, y=counts.value))+geom_point(alpha=0.2)+
  scale_y_continuous(trans='log2')+geom_smooth(col='red')+lims(x=c(0,20))+theme_bw()+
  labs(x='Weighted CN in gene', y='TPM (log2)')+theme(legend.position = "bottom")
ggsave("../../figures/scatterplot_CNweighted_TPM_nonzero.pdf", width=8)

ggplot(joint_counts_CN_subset[joint_counts_CN_subset$counts.value > 0,],
       aes(x=nearestGeneCN.value, y=counts.value))+geom_point(alpha=0.2)+
  scale_y_continuous(trans='log2')+geom_smooth(col='red')+lims(x=c(0,20))+theme_bw()+
  labs(x='Weighted CN in gene', y='TPM (log2)')+theme(legend.position = "bottom")
ggsave("../../figures/scatterplot_CNclosest_TPM_nonzero.pdf", width=8)

#-----------------------------------------------------------------------------------------------#

## for each gene, find the rank of both ploidy and RNASeq
## group by gene

ggplot(joint_counts_CN_normalised_highlyvar, aes(x=scaled_centered_weighted_CN,
                                                 y=scaled_centered_DESeq))+
  geom_point(alpha=0.9)+
  geom_smooth()+
  facet_wrap(.~factor(counts.Var1,
    levels=c( "119058org", "119148org", "151723org", "151761org", "151773org", "23868org",  "32077org",
              "50495org", "54059org",  "54276org", '118976org', '119127org', '54327org', '119178org','118947org' )),
    nrow=2, scales = "free")
ggplot(joint_counts_CN_normalised_highlyvar, aes(x=scaled_centered_weighted_CN,
                                                 y=scaled_centered_DESeq))+
  geom_point(alpha=0.9)+
  geom_smooth()

ggplot(joint_counts_CN_normalised_highlyvar, aes(x=scaled_neighbour_CN,
                                                 y=scaled_centered_DESeq))+
  geom_point(alpha=0.9)+
  geom_smooth()+
  facet_wrap(.~factor(counts.Var1,
                      levels=c( "119058org", "119148org", "151723org", "151761org", "151773org", "23868org",  "32077org",
                                "50495org", "54059org",  "54276org", '118976org', '119127org', '54327org', '119178org','118947org' )),
             nrow=2, scales = "free")

ggplot(joint_counts_CN_normalised_highlyvar, aes(x=nearestGeneCN.value,
                                                 y=counts.counts.value))+
  geom_point(alpha=0.9)+
  geom_smooth()+
  facet_wrap(.~factor(counts.Var1,
                      levels=c( "119058org", "119148org", "151723org", "151761org", "151773org", "23868org",  "32077org",
                                "50495org", "54059org",  "54276org", '118976org', '119127org', '54327org', '119178org','118947org' )),
             nrow=2, scales = "free")

# ggsave("../../figures/scatterplot_normCNweighted_normDESeq_perorg2.png", width=8, height=2.5)

#-----------------------------------------------------------------------------------------------#

example_centering = joint_counts_CN_normalised %>% filter(counts.Var2 == 'NOC2L') %>% dplyr::select(scaled_centered_weighted_CN)
mean(example_centering$scaled_centered_weighted_CN); sd(example_centering$scaled_centered_weighted_CN)
table(joint_counts_CN_ranked$rank_weighted_CN)

## there are genes which are duplicated - remove them
joint_counts_CN_ranked = joint_counts_CN %>% dplyr::select(c("counts.Var1", "counts.Var2", "counts.value", "rank_weighted_CN", "rank_DESeq", "rank_TPM"))
joint_counts_CN_ranked = joint_counts_CN_ranked[!(joint_counts_CN_ranked$CN.gene_name %in% joint_counts_CN_ranked$CN.gene_name[joint_counts_CN_ranked$rank_weighted_CN > length(unique(joint_counts_CN_ranked$counts.Var1))]),]
table(joint_counts_CN_ranked$rank_weighted_CN)

plot(joint_counts_CN_ranked$rank_weighted_CN, joint_counts_CN_ranked$rank_DESeq)

ggplot(joint_counts_CN_ranked, aes(rank_weighted_CN, rank_DESeq)) + geom_point()
ggplot(joint_counts_CN_ranked, aes(x = rank_weighted_CN, y = rank_DESeq)) + geom_tile()

ggplot((melt(table(joint_counts_CN_ranked$rank_weighted_CN, joint_counts_CN_ranked$rank_DESeq)/nrow(joint_counts_CN_ranked))),
       aes(x=Var1, y=Var2, fill=value))+geom_raster()+scale_fill_jcolors_contin("pal2", bias = 1.75) + theme_bw()+
  labs(x='Organoid rank for CN (weighted)', y='Organoid rank for gene expression (DESeq2 counts)')
ggsave("../../figures/rankplot_CNweighted_deseq.pdf", width=8)

joint_counts_CN_ranked_subsetorgs = joint_counts_CN_ranked %>% filter(! (counts.Var1 %in% c('118976org', '119127org', '54327org', '119178org', '118947org')))
ggplot((melt(table(joint_counts_CN_ranked_subsetorgs$rank_weighted_CN, joint_counts_CN_ranked_subsetorgs$rank_DESeq)/nrow(joint_counts_CN_ranked_subsetorgs))),
       aes(x=Var1, y=Var2, fill=value))+geom_raster()+scale_fill_jcolors_contin("pal2", bias = 1.75) + theme_bw()+
  labs(x='Organoid rank for CN (weighted)', y='Organoid rank for gene expression (DESeq2 counts)')

## is the rank of the excluded organoids any different from the rank of the included organoids?
grid.arrange(ggplot(joint_counts_CN_ranked, aes(x=(counts.Var1 %in% c('118976org', '119127org', '54327org', '119178org', '118947org')),
                                      y= rank_weighted_CN))+geom_violin(),
             ggplot(joint_counts_CN_ranked, aes(x=(counts.Var1 %in% c('118976org', '119127org', '54327org', '119178org', '118947org')),
                                   y= rank_DESeq))+geom_violin(), nrow=1)

df <- expand.grid(x = sample(x = 1:100, size = 50), y = sample(1:100, size = 50))
# default is compatible with geom_tile()
ggplot(df, aes(x, y)) + geom_raster()

## blob, but clear trend
ggplot(joint_counts_CN_normalised %>% filter(CN.value != 2), aes(x=scaled_centered_weighted_CN, y=scaled_centered_DESeq))+
  geom_point(alpha=0.01)+geom_smooth()

ggplot(joint_counts_CN_normalised %>% filter(CN.value != 2), aes(x=scaled_centered_weighted_CN, y=scaled_centered_DESeq))+
  geom_point(alpha=0.01)+geom_smooth()+facet_wrap(.~(CN.value < 2), scales = "free")

ggplot(joint_counts_CN_normalised %>% filter(CN.value > 2), aes(x=scaled_centered_weighted_CN, y=scaled_centered_DESeq))+
  geom_point(alpha=0.01)+geom_smooth()
ggsave("../../figures/scatterplot_normCNweighted_normDESeq.pdf", width=8)

ggplot(joint_counts_CN_normalised_excludingnormal %>% filter(CN.value > 2), aes(x=scaled_centered_weighted_CN, y=scaled_centered_DESeq))+
  geom_point(alpha=0.01)+geom_smooth()+facet_wrap(.~counts.Var1)
ggsave("../../figures/scatterplot_normCNweighted_normDESeq_perorg.pdf", width=8)

ggplot(joint_counts_CN_normalised_excludingnormal %>% filter(CN.value > 2), aes(x=scaled_centered_weighted_CN,
                                                                                y=scaled_centered_DESeq))+
  geom_point(alpha=0.01)+geom_smooth()+facet_wrap(.~factor(counts.Var1,
   levels=c( "119058org", "119148org", "151723org", "151761org", "151773org", "23868org",  "32077org",  "50495org", 
             "54059org",  "54276org", '118976org', '119127org', '54327org', '119178org', '118947org' )), nrow=2)
ggsave("../../figures/scatterplot_normCNweighted_normDESeq_perorg2.png", width=8, height=2.5)

ggplot(joint_counts_CN_normalised_excludingnormal %>% filter(CN.value > 2) %>%
         filter(counts.Var1 %in% c( "119058org", "119148org", "151723org", "151761org", "151773org", "23868org",  "32077org",  "50495org", 
                         "54059org",  "54276org")),
                 aes(x=scaled_centered_weighted_CN, y=scaled_centered_DESeq, col=counts.Var1))+
  geom_point(alpha=0.01)+geom_smooth()
ggsave("../../figures/scatterplot_normCNweighted_normDESeq_perorg2_subsetorgs.png", width=5, height=4)

ggplot(joint_counts_CN_normalised_excludingnormal %>% filter(CN.value > 2) ,
       aes(x=counts.Var1, y=scaled_centered_DESeq))+
  geom_boxplot(alpha=0.01)+
    facet_wrap(.~(counts.Var1 %in% c( "119058org", "119148org", "151723org", "151761org", "151773org", "23868org",  "32077org",  "50495org", 
                                     "54059org",  "54276org")), scales = "free_x")


## same when we include normal segments
ggplot(joint_counts_CN_normalised, aes(x=scaled_centered_weighted_CN,
                                                                                y=scaled_centered_DESeq))+
  geom_point(alpha=0.01)+geom_smooth()+facet_wrap(.~factor(counts.Var1,
                                                           levels=c( "119058org", "119148org", "151723org", "151761org", "151773org", "23868org",  "32077org",  "50495org", 
                                                                     "54059org",  "54276org", '118976org', '119127org', '54327org', '119178org', '118947org' )), nrow=2)
## same when we only include the very variable segments?
ggplot(joint_counts_CN_normalised %>% filter(CN.value > 2), aes(x=scaled_centered_weighted_CN,
                                                                y=scaled_centered_DESeq))+
  geom_point(alpha=0.01)+geom_smooth()+facet_wrap(.~factor(counts.Var1,
                                                           levels=c( "119058org", "119148org", "151723org", "151761org", "151773org", "23868org",  "32077org",  "50495org", 
                                                                     "54059org",  "54276org", '118976org', '119127org', '54327org', '119178org', '118947org' )), nrow=2)


## same but normalising only with nonnormal segments
ggplot(joint_counts_CN_normalised_excludingnormal, aes(x=scaled_centered_weighted_CN, y=scaled_centered_DESeq))+
  geom_point(alpha=0.002)+geom_smooth()
ggsave("../../figures/scatterplot_normCNweighted_normDESeq_nonormal.pdf", width=8)
ggsave("../../figures/scatterplot_normCNweighted_normDESeq_nonormal.png", width=4)

ggplot(joint_counts_CN_normalised_excludingnormal, aes(x=scaled_centered_weighted_CN, y=scaled_centered_TPM))+
  geom_point(alpha=0.002)+geom_smooth()
ggsave("../../figures/scatterplot_normCNweighted_normTPM_nonormal.png", width=4)

ggplot(joint_counts_CN_normalised_excludingnormal, aes(x=scaled_centered_weighted_CN, y=scaled_centered_DESeq))+
  geom_point(alpha=0.01)+geom_smooth()+facet_wrap(.~(CN.value < 2), scales = "free")

ggplot(joint_counts_CN_normalised_excludingnormal %>% filter(CN.value > 2), aes(x=scaled_centered_weighted_CN, y=scaled_centered_DESeq))+
  geom_point(alpha=0.01)+geom_smooth()
ggsave("../../figures/scatterplot_normCNweighted_normDESeq_nonormal_onlyamp.png", width=4)

## only nonnormal segments, and only with Stephane's genes
ggplot(joint_counts_CN_normalised_excludingnormal %>% filter(counts.Var2 %in% unlist(topvariable)),
       aes(x=scaled_centered_weighted_CN, y=scaled_centered_DESeq, col=CN.value < 2))+
  geom_point(alpha=0.1)+geom_smooth()+facet_wrap(.~(CN.value < 2), scales = "free")
ggsave("../../figures/scatterplot_normCNweighted_normDESeq_nonormal_highlyvarregions.pdf", width=8)

ggplot(joint_counts_CN_normalised %>% filter(CN.value != 2), aes(x=scaled_centered_weighted_CN, y=DESeq.count))+
  geom_point(alpha=0.2)+geom_smooth()+scale_y_continuous(trans = "log2")

## not great
ggplot(joint_counts_CN_normalised %>% filter(CN.value != 2), aes(x=CN.value, y=scaled_centered_DESeq))+
  geom_point(alpha=0.2)+geom_smooth()+lims(x=c(0,10))

## good trend
ggplot(joint_counts_CN_normalised  %>% filter(CN.value != 2), aes(x=scaled_centered_weighted_CN, y=scaled_centered_TPM))+
  geom_point(alpha=0.2)+geom_smooth()

## not great, but good trend
ggplot(joint_counts_CN_normalised_excludingnormal, aes(x=scaled_centered_weighted_CN, y=scaled_centered_TPM))+
  geom_point(alpha=0.2)+geom_smooth()

# save.image("../output/analyse_joing_counts_CN.RData")


##------------------------------------------------------------------------------------------------------------##
## if a gene has a copy number greater than two, is its expression higher than in normal cells?
df_gene_characteristics = cbind(df_gene_characteristics,
 joint_counts_CN %>% group_by(CN.gene_name) %>%
  summarise(averageCN=mean(CN.value),
            sdCN=sd(CN.value),
            averageDESeq2=mean(DESeq.count),
            sdDESeq2=sd(DESeq.count)
            ))

plot(df_gene_characteristics$averageCN, df_gene_characteristics$averageDESeq2)

df_gene_characteristics_normalsamples = deseq_counts %>%
  filter(DESEq.sample %in% c('fal01', 'fal07', 'fal08', 'fal09', 'fal10')) %>%
  group_by(Gene) %>% 
  summarise(averageDESeq2=mean(DESeq.count),
            sdDESeq2=sd(DESeq.count)
  )

df_gene_characteristics_normalsamples = df_gene_characteristics_normalsamples[match(df_gene_characteristics$Gene,
      df_gene_characteristics_normalsamples$Gene),]

df_gene_characteristics = cbind.data.frame(df_gene_characteristics,
                                           normalsamples=df_gene_characteristics_normalsamples)
df_gene_characteristics = df_gene_characteristics[!is.na(df_gene_characteristics$Gene),]
##' for each gene and each organoid sample, check if the cn is greater than 2.
##' Then check if its GE value is greater than the average GE for normal samples

df_gene_characteristics = as.data.frame(df_gene_characteristics)
DESeq2_increase_withCN = apply(joint_counts_CN, 1, function(i){
  gene_match <- match(i[['Gene']], df_gene_characteristics$Gene)
  if( (sum(gene_match, na.rm = T) == 0)){
    return(NA)
  }else{
    return(as.numeric(i[['DESeq.count']]) > df_gene_characteristics$averageDESeq2[gene_match])
  }
})

table(DESeq2_increase_withCN)
joint_counts_CN$DESeq2_increase_withCN  = (DESeq2_increase_withCN)

df_gene_characteristics[df_gene_characteristics$Gene == "AAMDC",'averageDESeq2']
joint_counts_CN[joint_counts_CN$counts.Var2 == "AAMDC",]$DESeq.count
joint_counts_CN[joint_counts_CN$counts.Var2 == "AAMDC",]$DESeq2_increase_withCN

table((joint_counts_CN$DESeq2_increase_withCN))
length(joint_counts_CN$DESeq2_increase_withCN)

# DESeq2_increase_withCN[sapply(DESeq2_increase_withCN, length) > 1]
# table(sapply(DESeq2_increase_withCN, length))
# 
# table(deseq_counts$DESEq.sample)

table(joint_counts_CN$DESeq2_increase_withCN)
table( (joint_counts_CN %>% filter(CN.value > 2) %>% dplyr::select(DESeq2_increase_withCN))[,2] )

mean_coord_CN_DESeq2 = joint_counts_CN %>% dplyr::group_by(Gene) %>% dplyr::summarise(mean_coord_CN_DESeq2=mean(DESeq2_increase_withCN, na.rm=T))
# mean_coord_CN_DESeq2 = mean_coord_CN_DESeq2[!is.na(mean_coord_CN_DESeq2$mean_coord_CN_DESeq2),]
ntop_coord = 100
plot(sort(mean_coord_CN_DESeq2$mean_coord_CN_DESeq2))
plot(sort(mean_coord_CN_DESeq2$mean_coord_CN_DESeq2, decreasing = T)[1:ntop_coord])

# mean_coord_CN_DESeq2$Gene = factor(mean_coord_CN_DESeq2$Gene, levels=mean_coord_CN_DESeq2$Gene[order(mean_coord_CN_DESeq2$mean_coord_CN_DESeq2, decreasing = T)])
mean_coord_CN_DESeq2 = mean_coord_CN_DESeq2[order(mean_coord_CN_DESeq2$mean_coord_CN_DESeq2, decreasing = T),]
mean_coord_CN_DESeq2$rank = factor(1:nrow(mean_coord_CN_DESeq2))
ggplot(droplevels(mean_coord_CN_DESeq2[1:ntop_coord,]), aes(label=Gene, x=rank,
                                 y=mean_coord_CN_DESeq2))+
  geom_point()+geom_label_repel()

ggplot((mean_coord_CN_DESeq2), aes(label=ifelse(test = Gene %in% subset_genes_of_interest, yes = Gene, no=NA),
                                             x=rank, y=mean_coord_CN_DESeq2))+
  geom_point()+theme_classic()+
  theme(axis.text.x=element_blank())+geom_label_repel()
ggsave("../../figures/coordinated_CN_DESeq_goi.pdf", width = 8, height = 8)

dim(mean_coord_CN_DESeq2)
dim(df_gene_characteristics) 

plot(mean_coord_CN_DESeq2$mean_coord_CN_DESeq2, df_gene_characteristics$r2_normCNnormDESeq[match(mean_coord_CN_DESeq2$Gene, df_gene_characteristics$Gene)])
plot(df_gene_characteristics$cor_normCNnormDESeq, df_gene_characteristics$r2_normCNnormDESeq)
plot(df_gene_characteristics$cor_normCNnormDESeq, log(df_gene_characteristics$averageCN))

ggplot(df_gene_characteristics,
  aes(x=cor_normCNnormDESeq, y=averageCN,
                                    # label=ifelse( (log(averageCN) > 1) & (cor_normCNnormDESeq > 0.8), yes = Gene, no = NA  )))+
      label=ifelse( (log(averageCN) > 1.5) & (cor_normCNnormDESeq > 0.5), yes = Gene, no = NA  )))+
  geom_point()+scale_y_continuous(trans = "log2")+geom_label_repel(size=2, max.overlaps = 30)+theme_bw()
ggsave("../../figures/coordinated_CN_DESeq_and_CNvalue.pdf", width = 8, height = 8)

## example of a 'good' result
mean_coord_CN_DESeq2 %>% filter(Gene == 'TOP1')
df_gene_characteristics %>% filter(Gene == 'TOP1') %>% dplyr::select(normalsamples.averageDESeq2)
joint_counts_CN%>% filter(counts.Var2 == 'TOP1') %>% dplyr::select(DESeq.count)
  
ggplot(joint_counts_CN_subset %>% filter(counts.Var2 == 'FBL'),
       aes(x=scaled_centered_weighted_CN, y=scaled_centered_DESeq, label=counts.Var1), alpha=0.2)+
  geom_smooth()+geom_label_repel()+facet_wrap(.~labels)+
  guides(col=FALSE)+theme_bw()

ggplot(joint_counts_CN_subset %>% filter(counts.Var2 == 'MYC'),
       aes(x=scaled_centered_weighted_CN, y=scaled_centered_DESeq, label=counts.Var1), alpha=0.2)+
  geom_smooth()+geom_label_repel()+facet_wrap(.~labels)+geom_point()+
  guides(col=FALSE)+theme_bw()
ggplot(joint_counts_CN_subset %>% filter(counts.Var2 == 'MYC'),
       aes(x=CN.value, y=DESeq.count, label=counts.Var1), alpha=0.2)+
  geom_smooth()+geom_label_repel()+facet_wrap(.~labels)+geom_point()+
  guides(col=FALSE)+theme_bw()
ggplot(joint_counts_CN_subset %>% filter(counts.Var2 == 'KRT6B'),
       aes(x=CN.value, y=DESeq.count, label=counts.Var1), alpha=0.2)+
  geom_smooth()+geom_label_repel()+geom_point()+
  guides(col=FALSE)+theme_bw()+ggtitle('KRT6Bs')
ggsave("../../figures/example_KRT6B.pdf", width = 4, height = 4)
ggplot(joint_counts_CN_subset %>% filter(counts.Var2 == 'AKT2'),
       aes(x=CN.value, y=DESeq.count, label=counts.Var1), alpha=0.2)+
  geom_smooth()+geom_label_repel()+geom_point()+
  guides(col=FALSE)+theme_bw()+ggtitle('AKT2')
ggsave("../../figures/example_AKT2.pdf", width = 4, height = 4)


plot(df_gene_characteristics$slope_centered, df_gene_characteristics$r2_normCNnormDESeq)
head(df_gene_characteristics$slope_centered)
head(df_gene_characteristics$r2_normCNnormDESeq)

##' for each gene, average the bottom three organoids (sorted by CN) and check if
##'  the expression of the remaining is higher than for those

average_comparison_CN_DESeq = sapply(unique(joint_counts_CN_subset$counts.Var2), function(i){
  .x = joint_counts_CN[joint_counts_CN$counts.Var2 == i,]
  if(any(is.na(.x$DESeq.count))){
    return(NA)
  }else{
    .averg_DESeq = mean(unlist(.x[order(.x$CN.value, decreasing = F)[1:3],'DESeq.count']))
    return(.x[order(.x$CN.value, decreasing = F)[-(1:3)],'DESeq.count'] > .averg_DESeq)
  }
})

## in general we would expect the value to be positive

length(unique(joint_counts_CN_subset$counts.Var2))
average_comparison_CN_DESeq[[1]]

sum(is.na(average_comparison_CN_DESeq))
df_average_bottomCN <- cbind.data.frame(Gene=unique(joint_counts_CN_subset$counts.Var2),
                 average_comparison_CN_DESeq=sapply(average_comparison_CN_DESeq, mean),
                 label=ifelse(unique(joint_counts_CN_subset$counts.Var2) %in% subset_genes_of_interest,
                              yes = as.character(unique(joint_counts_CN_subset$counts.Var2)), no = NA))
df_average_bottomCN$Gene <- factor(df_average_bottomCN$Gene, levels=df_average_bottomCN$Gene[order(df_average_bottomCN$average_comparison_CN_DESeq, na.last = T)])
ggplot(droplevels(df_average_bottomCN),
       aes(x=Gene, y=average_comparison_CN_DESeq, label=as.character(label)))+
  geom_point()+#geom_label_repel(max.overlaps = 30, segment.size = 0.01)+
  # geom_label()+
  geom_text_repel(force = .01, direction = "y", nudge_x = 0, nudge_y = .1)+
  theme(axis.text.x=element_blank(), axis.line.x.bottom =element_blank(), 
        panel.grid.minor = element_line(size = 0.1, colour = "black"))+
  geom_hline(yintercept = 0.5, lty='dashed')
ggsave("../../figures/average_bottomCN_DESeq.pdf", width = 4, height = 4)

## For each gene, plot the CN

median_CN = aggregate(data = joint_counts_CN %>% filter(CN.gene_name %in% topvariable$V1),
          CN.value ~ CN.gene_name, FUN = median)

ggplot(joint_counts_CN %>% filter(CN.gene_name %in% topvariable$V1),
       aes(x=factor(CN.gene_name, levels=median_CN$CN.gene_name[order(median_CN$CN.value)]),
           y=CN.value, group=CN.gene_name))+geom_violin()
ggsave("../../figures/CN_violinplots.pdf", width = 9, height = 4)

median_CN_goi = aggregate(data = joint_counts_CN %>% filter(counts.Var2 %in% subset_genes_of_interest),
                      CN.value ~ counts.Var2, FUN = median)
ggplot(joint_counts_CN %>% filter(counts.Var2 %in% subset_genes_of_interest),
       aes(x=factor(counts.Var2, levels=median_CN_goi$counts.Var2[order(median_CN_goi$CN.value)]),
           y=CN.value, group=counts.Var2, col=(scaled_centered_DESeq)))+geom_violin()+
  geom_point()+scale_y_continuous(trans = "log2")+
  geom_line(aes(group=counts.Var1), col='black', alpha=0.1)+
  geom_line(data = median_CN_goi, aes(x=factor(counts.Var2, levels=median_CN_goi$counts.Var2[order(median_CN_goi$CN.value)]),
                                      y=CN.value, group=1), col='red')+
  theme(legend.position = "bottom")+labs(x='Gene, sorted by median weighted CN', y='Weighted CN (log scale)')+
  theme(axis.text.x=element_text(angle = 90, hjust = 0))+
  scale_color_jcolors_contin("pal2", reverse = TRUE, bias = 2.25)
ggsave("../../figures/CN_violinplots_goi.pdf", width = 6, height = 6)

library(jcolors)
ggplot(joint_counts_CN %>% filter(counts.Var2 %in% subset_genes_of_interest),
       aes(x=factor(counts.Var2, levels=median_CN_goi$counts.Var2[order(median_CN_goi$CN.value)]),
           y=scaled_centered_DESeq, col=(scaled_centered_weighted_CN)))+geom_violin()+
  geom_point(size=3)+
  geom_line(aes(group=counts.Var1), col='black', alpha=0.1)+
  theme(legend.position = "bottom")+labs(x='Gene, sorted by median weighted CN', y='Weighted CN (log scale)')+
  scale_color_jcolors_contin("pal2", reverse = TRUE, bias = 2.25)+
  theme(axis.text.x=element_text(angle = 90, hjust = 0))
ggsave("../../figures/CN_violinplots_goi_2.pdf", width = 6, height = 6)

ggplot(joint_counts_CN %>% filter(counts.Var2 %in% c('AKT1', 'AKT3')),
       aes(x=factor(counts.Var2, levels=median_CN_goi$counts.Var2[order(median_CN_goi$CN.value)]),
           y=scaled_centered_DESeq, col=(scaled_centered_weighted_CN)))+geom_violin()+
  geom_point(size=3)+
  geom_line(aes(group=counts.Var1), col='black', alpha=0.1)+
  theme(legend.position = "bottom")+labs(x='Gene, sorted by median weighted CN', y='Weighted CN (log scale)')+
  scale_color_jcolors_contin("pal2", reverse = TRUE, bias = 2.25)+
  theme(axis.text.x=element_text(angle = 90, hjust = 0))

df_gene_characteristics = cbind.data.frame(df_gene_characteristics,
                                           df_average_bottomCN = df_average_bottomCN[match(df_gene_characteristics$Gene,
                                                                                           df_average_bottomCN$Gene),])


plot(df_gene_characteristics$df_average_bottomCN.average_comparison_CN_DESeq,
df_gene_characteristics$r2_normCNnormDESeq)

ggplot(df_gene_characteristics, aes(x=df_average_bottomCN.average_comparison_CN_DESeq, y=r2_normCNnormDESeq,
                                    label=ifelse(Gene %in% subset_genes_of_interest,
                                                 yes = Gene, no = NA )))+
  geom_boxplot(aes(group=df_average_bottomCN.average_comparison_CN_DESeq),
               width=0.5/length(unique(df_gene_characteristics$df_average_bottomCN.average_comparison_CN_DESeq)),
               col='blue')+
  geom_point()+geom_label_repel()
ggsave("../../figures/r2normCNnormDESeq_vs_averagebottomCN.pdf", width = 6, height = 6)

ggplot(df_gene_characteristics, aes(x=df_average_bottomCN.average_comparison_CN_DESeq, y=r2_normCNnormDESeq,
                                    label=ifelse( (df_average_bottomCN.average_comparison_CN_DESeq >= .8) &
                                                    (r2_normCNnormDESeq > 0.5),
                                                 yes = Gene, no = NA )))+
  geom_boxplot(aes(group=df_average_bottomCN.average_comparison_CN_DESeq),
               width=0.5/length(unique(df_gene_characteristics$df_average_bottomCN.average_comparison_CN_DESeq)),
               col='blue')+
  geom_point()+geom_label_repel()
ggsave("../../figures/r2normCNnormDESeq_vs_averagebottomCN_best.pdf", width = 6, height = 6)

plot(joint_counts_CN_subset$scaled_centered_weighted_CN_mean,
     log(joint_counts_CN_subset$scaled_centered_weighted_CN))
ggplot(joint_counts_CN_subset, aes(x=scaled_centered_weighted_CN_mean, y=scaled_centered_weighted_CN, col=scaled_centered_DESeq))+
  geom_point()+scale_y_continuous(trans = "log")

ggplot(joint_counts_CN_subset, aes(x=scaled_centered_weighted_CN_mean, y=scaled_centered_DESeq))+
  geom_point()+geom_smooth()#+scale_y_continuous(trans = "log")

