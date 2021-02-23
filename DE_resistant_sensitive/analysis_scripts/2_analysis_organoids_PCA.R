#--------------------------------------------------------------------------------#
rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(ggplot2)
library(ggrepel)
library(reshape2)
library(readxl)
require(jcolors)
library(biomaRt)
require(GSVA)
require(GSVAdata)
library(lsa) ## for cosine similarity
source("other_scripts/functions.R")

#--------------------------------------------------------------------------------#

#--------------------------------------------------------------------------------#
TPM = read.csv("../files/20191218_ViasM_BJ_orgaBrs_tpm.csv", stringsAsFactors = FALSE)
normalised_counts = read.csv("../files/counts_norm.csv", stringsAsFactors = FALSE)
rownames(normalised_counts) = normalised_counts[,1]
normalised_counts = normalised_counts[,-1]
TPM = TPM[!duplicated(TPM$gene_name),]
# TPM = TPM[!(is.na(TPM)),]

## Remove mitocondrial genes (outliers) and other weird chromosomes
TPM = TPM[TPM$seq_name %in% c(as.character(1:22), 'X', 'Y'),]

## Remove NA genes
TPM = TPM[!(is.na(TPM$gene_id)),]

# TPM[duplicated(TPM$gene_name),]

TPM_data = TPM[,!(colnames(TPM) %in% c('gene_id', 'gene_name', 'entrezid', 'gene_biotype',
                                       'gene_seq_start', 'gene_seq_end', 'seq_name',
                                       'seq_strand', 'symbol'))]
rownames(TPM_data) = TPM$gene_id

#--------------------------------------------------------------------------------#

#--------------------------------------------------------------------------------#
# Renaming
renaming1 = read_excel("../files/PDOnameProperSample_sWGS_RNAseq.xlsx")
 
colnames(TPM_data)
renaming1 = renaming1[match(colnames(TPM_data), gsub('-', '', renaming1$sampleNameRNAseq)),]
renaming_NC = renaming1[match( gsub('[.]', '', colnames(normalised_counts)), gsub('-', '', renaming1$sampleNameRNAseq)),]

colnames(TPM_data)[!is.na(renaming1$PDO)] = renaming1$PDO[!is.na(renaming1$PDO)]
colnames(normalised_counts)[!is.na(renaming_NC$PDO)] = renaming_NC$PDO[!is.na(renaming_NC$PDO)]
colnames(TPM_data)
colnames(normalised_counts)

## Remove normal samples
TPM_data = TPM_data[,!grepl('FT', colnames(TPM_data))]
normalised_counts = normalised_counts[,!grepl('FT', colnames(normalised_counts))]
normalised_counts = normalised_counts[rowSums(normalised_counts)>0,]
#--------------------------------------------------------------------------------#

#--------------------------------------------------------------------------------#
# DESeq2::normalizationFactors(TPM_data)
# 
# library(DESeq2)
# DESeq2::DESeqDataSet(SummarizedExperiment(assays=round(t(TPM_data*1e6))))
prcomp_normalisedcounts = prcomp(t(normalised_counts), scale. = TRUE, center = TRUE)
# prcomp_normalisedcounts = prcomp(t(normalised_counts[,!grepl('FT', colnames(normalised_counts))]), scale. = FALSE, center = TRUE)
# prcomp_normalisedcounts_notscaled = prcomp(t(normalised_counts), scale. = FALSE, center = TRUE)
plot(prcomp_normalisedcounts$x[,1:2])

size_points = 2
prcomp_res = prcomp(t(TPM_data), scale. = FALSE, center = TRUE)
labels = colnames(TPM_data)
eigs <- prcomp_res$sdev^2
ggplot(data=cbind.data.frame(prcomp_res$x[,1:2], names=labels),
       aes(x=PC1, y=PC2, label=names))+
  geom_point(size=size_points)+
  geom_label_repel()+
  labs(x=paste0('PC1 (', round(100*eigs[1]/sum(eigs), 2), '%)' ),
       y=paste0('PC2 (', round(100*eigs[2]/sum(eigs), 2), '%)' ))+ggtitle('PCA')
ggsave("../figures/PCA_RNASeq/PCA_TPM.pdf")

all(rownames(prcomp_res$rotation) == TPM$gene_id)
## good, it just takes a bit of time
# ggplot(cbind.data.frame(prcomp_res$rotation[,1:2], col=factor(TPM$seq_name)), aes(x=PC1, y=PC2))+
#   geom_point(data=cbind.data.frame(prcomp_res$rotation[,1:2]), aes(x=PC1, y=PC2), alpha=0.2)+
#   geom_point(col='blue')+facet_wrap(.~col)
# ggsave("../figures/PCA_RNASeq/PCA_TPM_per_chrom.png")

## Very time consuming
## GO term analysis on clustered loadings
# hclust_on_pc1_and_2 = hclust(dist(prcomp_res$rotation[,1:2]))
# cutree_on_pc1_and_2 = cutree(hclust_on_pc1_and_2, k = 100)
# cutree_on_pc1_and_2_k5 = cutree(hclust_on_pc1_and_2, k = 5)

# table(cutree_on_pc1_and_2)
# 
# plot(prcomp_res$rotation[,1:2], col=cutree_on_pc1_and_2)
#--------------------------------------------------------------------------------#

#--------------------------------------------------------------------------------#
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- TPM$gene_id

t2g <- getBM(
  attributes = c('ensembl_gene_id', 'entrezgene_id', 'external_gene_name'),
  values = TPM$gene_id,
  filter = 'ensembl_gene_id',
  mart = mart, useCache = FALSE)

entrezgene_id_matched = t2g[match(TPM$gene_id, t2g$ensembl_gene_id),'entrezgene_id']
entrezgene_id_matched_normalisedcounts = t2g[match(rownames(normalised_counts), t2g$ensembl_gene_id),'entrezgene_id']
naming_id_matched_normalisedcounts = t2g[match(rownames(normalised_counts), t2g$ensembl_gene_id),]

#--------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------#
data(c2BroadSets) ## from GSVAdata

TPM_data_entrez = TPM_data[!is.na(entrezgene_id_matched),]
entrezgene_id_matched = entrezgene_id_matched[!is.na(entrezgene_id_matched)]
dim(TPM_data_entrez)
length(entrezgene_id_matched)
TPM_data_entrez = TPM_data_entrez[!duplicated(entrezgene_id_matched),]
entrezgene_id_matched = entrezgene_id_matched[!duplicated(entrezgene_id_matched)]
rownames(TPM_data_entrez) = entrezgene_id_matched
dim(TPM_data_entrez)
length(entrezgene_id_matched)

# library(lsa)
# hclust_on_pc1_and_2_NC_cosim_top = hclust(dist(prcomp_normalisedcounts$rotation[,1:2]))
# cosine_dist= outer(1:nrow(prcomp_res$rotation), 1:nrow(prcomp_res$rotation),
#        function(i,j) cosine(prcomp_res$rotation[,1:2][i,], prcomp_res$rotation[,1:2][j,]))
# plot(prcomp_normalisedcounts$rotation[,1:2], col=cutree_on_pc1_and_2_NC)

results_Gsva <- gsva(expr = as(TPM_data_entrez, 'matrix'), gset.idx.list = c2BroadSets, min.sz=10, max.sz=500, verbose=TRUE)

dim(results_Gsva)

results_Gsva[,match(colnames(results_Gsva), rownames(prcomp_res$x))]
pca_with_gsva_annotation = cbind.data.frame(prcomp_res$x, t(results_Gsva[,match(colnames(results_Gsva), rownames(prcomp_res$x))]),
                                            labels=rownames(prcomp_res$x))
pca_with_gsva_annotation = cbind.data.frame(prcomp_res$x, t(results_Gsva[,match(colnames(results_Gsva), rownames(prcomp_res$x))]),
                                            labels=rownames(prcomp_res$x))

ggplot(pca_with_gsva_annotation, aes(x=PC1, y=PC2, col=KAUFFMANN_DNA_REPAIR_GENES, label=labels))+
  geom_point()+
  geom_label_repel()+
  scale_color_jcolors_contin("pal3", reverse = TRUE, bias = 2.25)+labs(col = "Threshold", shape="Amount")+
  ggtitle('KAUFFMANN_DNA_REPAIR_GENES')
ggsave("../figures/PCA_RNASeq/PCA_TPM_KAUFFMANN_DNA_REPAIR_GENES.png", width = 6.5, height = 5)

ggplot(pca_with_gsva_annotation, aes(x=PC1, y=PC2, col=OUELLET_CULTURED_OVARIAN_CANCER_INVASIVE_VS_LMP_UP, label=labels))+
  geom_point()+
  geom_label_repel()+
  scale_color_jcolors_contin("pal3", reverse = TRUE, bias = 2.25)+labs(col = "Threshold", shape="Amount")+
  ggtitle('OUELLET_CULTURED_OVARIAN_CANCER_INVASIVE_VS_LMP_UP')
ggsave("../figures/PCA_RNASeq/PCA_TPM_OUELLET_CULTURED_OVARIAN_CANCER_INVASIVE_VS_LMP_UP.png", width = 6.5, height = 5)

colnames(pca_with_gsva_annotation)[grepl('OVARIAN', colnames(pca_with_gsva_annotation))]

ggplot(pca_with_gsva_annotation, aes(x=PC1, y=PC2, col=BONOME_OVARIAN_CANCER_POOR_SURVIVAL_DN, label=labels))+
  geom_point()+
  geom_label_repel()+
  scale_color_jcolors_contin("pal3", reverse = TRUE, bias = 2.25)+labs(col = "Threshold", shape="Amount")+
  ggtitle('')
ggsave("../figures/PCA_RNASeq/PCA_TPM_BONOME_OVARIAN_CANCER_POOR_SURVIVAL_DN.png", width = 6.5, height = 5)
#--------------------------------------------------------------------------------#

#--------------------------------------------------------------------------------#
## same with normalised counts
TPM_data_entrez_NC = normalised_counts[!is.na(entrezgene_id_matched_normalisedcounts),]
entrezgene_id_matched_normalisedcounts = entrezgene_id_matched_normalisedcounts[!is.na(entrezgene_id_matched_normalisedcounts)]
TPM_data_entrez_NC = TPM_data_entrez_NC[!duplicated(entrezgene_id_matched_normalisedcounts),]
entrezgene_id_matched_normalisedcounts = entrezgene_id_matched_normalisedcounts[!duplicated(entrezgene_id_matched_normalisedcounts)]
rownames(TPM_data_entrez_NC) = entrezgene_id_matched_normalisedcounts
results_Gsva_counts <- gsva(expr = as(TPM_data_entrez_NC, 'matrix'), gset.idx.list = c2BroadSets, min.sz=10, max.sz=500, verbose=TRUE)
pca_with_gsva_annotation_NC = cbind.data.frame(prcomp_normalisedcounts$x, t(results_Gsva_counts[,match(colnames(results_Gsva_counts), rownames(prcomp_normalisedcounts$x))]),
                                               labels=rownames(prcomp_normalisedcounts$x))

ggplot(pca_with_gsva_annotation_NC, aes(x=PC1, y=PC2, col=KAUFFMANN_DNA_REPAIR_GENES, label=labels))+
  geom_point()+
  geom_label_repel()+
  scale_color_jcolors_contin("pal3", reverse = TRUE)+labs(col = "Threshold", shape="Amount")+
  ggtitle('KAUFFMANN_DNA_REPAIR_GENES')
ggsave("../figures/PCA_RNASeq/PCA_counts_KAUFFMANN_DNA_REPAIR_GENES.png", width = 6.5, height = 5)

ggplot(pca_with_gsva_annotation_NC, aes(x=PC1, y=PC2, col=OUELLET_CULTURED_OVARIAN_CANCER_INVASIVE_VS_LMP_UP, label=labels))+
  geom_point()+
  geom_label_repel()+
  scale_color_jcolors_contin("pal3", reverse = TRUE, bias = 2.25)+labs(col = "Threshold", shape="Amount")+
  ggtitle('OUELLET_CULTURED_OVARIAN_CANCER_INVASIVE_VS_LMP_UP')
ggsave("../figures/PCA_RNASeq/PCA_counts_OUELLET_CULTURED_OVARIAN_CANCER_INVASIVE_VS_LMP_UP.png", width = 6.5, height = 5)

results_Gsva_counts_only_variable = results_Gsva_counts[order(apply(results_Gsva_counts, 1, var), decreasing = T)[1:100],]
pdf("../figures/PCA_RNASeq/heatmap_gsva_counts.pdf", height = 20, width = 18)
pheatmap::pheatmap(results_Gsva_counts_only_variable)
dev.off()
pdf("../figures/PCA_RNASeq/heatmap_gsva_counts_all.pdf", height = 80, width = 30)
pheatmap::pheatmap(results_Gsva_counts)
dev.off()
pdf("../figures/PCA_RNASeq/heatmap_gsva_counts_ovarian.pdf", height = 6, width = 15)
pheatmap::pheatmap(results_Gsva_counts[grepl('OVAR', rownames(results_Gsva_counts)),])
dev.off()

head(melt(results_Gsva))
tail(sort(apply(results_Gsva, 1, var)), n=20)

# Very time consuming
#hclust_on_pc1_and_2_NC = hclust(dist(prcomp_normalisedcounts$rotation[,1:2]))
# cutree_on_pc1_and_2_NC = cutree(hclust_on_pc1_and_2_NC, k = 10)
# table(cutree_on_pc1_and_2_NC)

## get only the most variable genes
highlyVar = names(sort(apply(normalised_counts, 1, var), decreasing = T)[1:500])
prcomp_normalisedcounts_highlyVar = prcomp_normalisedcounts
prcomp_normalisedcounts_highlyVar$rotation = prcomp_normalisedcounts_highlyVar$rotation[rownames(prcomp_normalisedcounts_highlyVar$rotation) %in% highlyVar,]
size_prcomp_normalisedcounts_highlyVar = apply(prcomp_normalisedcounts_highlyVar$rotation[,1:2], 1, function(i) sqrt(sum(i**2)))
top_size_size_prcomp_normalisedcounts_highlyVar = order(size_prcomp_normalisedcounts_highlyVar, decreasing = T)[1:250]
top_size_size_prcomp_normalisedcounts_highlyVar2 = order(size_prcomp_normalisedcounts_highlyVar, decreasing = T)[1:25]
ggplot()+
  geom_segment(data=cbind.data.frame(names=t2g$external_gene_name[match(rownames(prcomp_normalisedcounts_highlyVar$rotation[top_size_size_prcomp_normalisedcounts_highlyVar2,1:2]), t2g$ensembl_gene_id)],
                                     prcomp_normalisedcounts_highlyVar$rotation[top_size_size_prcomp_normalisedcounts_highlyVar2,1:2]),
               aes(x=0, xend=PC1, y=0, yend=PC2, label=names))+
  geom_label_repel(data=cbind.data.frame(names=t2g$external_gene_name[match(rownames(prcomp_normalisedcounts_highlyVar$rotation[top_size_size_prcomp_normalisedcounts_highlyVar2,1:2]), t2g$ensembl_gene_id)],
                                         prcomp_normalisedcounts_highlyVar$rotation[top_size_size_prcomp_normalisedcounts_highlyVar2,1:2]),
                   aes(x=PC1, y=PC2, label=names))

## See the genes with the highest loadings, even if they are not the most variable
size_prcomp_normalisedcounts = apply(prcomp_normalisedcounts$rotation[,1:2], 1, function(i) sqrt(sum(i**2)))
top_size_size_prcomp_normalisedcounts = order(size_prcomp_normalisedcounts, decreasing = T)[1:250]
ggplot()+
  geom_segment(data=cbind.data.frame(names=t2g$external_gene_name[match(rownames(prcomp_normalisedcounts$rotation[top_size_size_prcomp_normalisedcounts,1:2]), t2g$ensembl_gene_id)],
                                     prcomp_normalisedcounts$rotation[top_size_size_prcomp_normalisedcounts,1:2]),
               aes(x=0, xend=PC1, y=0, yend=PC2, label=names))
  # geom_label_repel(data=cbind.data.frame(names=t2g$external_gene_name[match(rownames(prcomp_normalisedcounts$rotation[top_size_size_prcomp_normalisedcounts,1:2]), t2g$ensembl_gene_id)],
  #                                        prcomp_normalisedcounts$rotation[top_size_size_prcomp_normalisedcounts,1:2]),
  #                  aes(x=PC1, y=PC2, label=names))
ggplot()+
  geom_segment(data=cbind.data.frame(names=t2g$external_gene_name[match(rownames(prcomp_normalisedcounts$rotation[,1:2]), t2g$ensembl_gene_id)],
                                     prcomp_normalisedcounts$rotation[,1:2]),
               aes(x=0, xend=PC1, y=0, yend=PC2, label=names), alpha=0.02)+
  geom_point(data=(cbind.data.frame(label=rownames(prcomp_normalisedcounts$x), prcomp_normalisedcounts$x[,1:2])),
                   aes(x=PC1*1e-4, y=PC2*1e-4), col='red', size=3)+
  geom_label_repel(data=(cbind.data.frame(label=rownames(prcomp_normalisedcounts$x), prcomp_normalisedcounts$x[,1:2])),
             aes(x=PC1*1e-4, y=PC2*1e-4, label=label), col='red', size=3)
ggsave("../figures/PCA_RNASeq/PCA_counts_all_loadings.png", width = 8, height = 8)

## cosine similarity of PC1, PC2
cosine_sim_NC= outer(1:nrow(prcomp_normalisedcounts_highlyVar$rotation), 1:nrow(prcomp_normalisedcounts_highlyVar$rotation),
                      Vectorize(function(i,j) cosine(prcomp_normalisedcounts_highlyVar$rotation[,1:2][i,], prcomp_normalisedcounts_highlyVar$rotation[,1:2][j,])))
cosine_dist_NC =  acos(cosine_sim_NC)/pi
hclust_on_pc1_and_2_NC_cos = hclust(as.dist(cosine_dist_NC))
cutree_on_pc1_and_2_NC_cos = cutree(hclust_on_pc1_and_2_NC_cos, k = 10)
cutree_on_pc1_and_2_NC_cos = cutree(hclust_on_pc1_and_2_NC_cos, k = 10)

ggplot(cbind.data.frame(prcomp_normalisedcounts_highlyVar$rotation[,1:2], col=cutree_on_pc1_and_2_NC_cos),
       aes(x=0, xend=PC1, y=0, yend=PC2, col=factor(col)))+
  geom_segment()

## plot only the genes which are most important for the first two PCs
dim(prcomp_normalisedcounts$rotation)
## compute the length of the vectors in PCs 1 and 2
size_prcomp_normalisedcounts = apply(prcomp_normalisedcounts$rotation[,1:2], 1, function(i) sqrt(sum(i**2)))
# size_prcomp_normalisedcounts_notscaled = apply(prcomp_normalisedcounts_notscaled$rotation[,1:2], 1, function(i) sqrt(sum(i**2)))
top_size_prcomp_normalisedcounts = order(size_prcomp_normalisedcounts, decreasing = T)[1:1500]
top_size_prcomp_normalisedcounts2 = order(size_prcomp_normalisedcounts, decreasing = T)[1:700]

## the top genes are all in the same quadrant
eigs_normalisedcounts <- prcomp_normalisedcounts$sdev^2
ggplot()+
  geom_segment(data=cbind.data.frame(names=t2g$external_gene_name[match(rownames(prcomp_normalisedcounts$rotation[top_size_prcomp_normalisedcounts2,1:2]), t2g$ensembl_gene_id)],
                                     prcomp_normalisedcounts$rotation[top_size_prcomp_normalisedcounts2,1:2]),
               aes(x=0, xend=PC1, y=0, yend=PC2, label=names))+
  geom_label_repel(data=cbind.data.frame(names=t2g$external_gene_name[match(rownames(prcomp_normalisedcounts$rotation[top_size_prcomp_normalisedcounts2,1:2]), t2g$ensembl_gene_id)],
                                         prcomp_normalisedcounts$rotation[top_size_prcomp_normalisedcounts2,1:2]),
                   aes(x=PC1, y=PC2, label=names))+
  geom_label_repel(data=(cbind.data.frame(label=rownames(prcomp_normalisedcounts$x), prcomp_normalisedcounts$x[,1:2])),
                   # aes(x=PC1*5e-5, y=PC2*5e-5, label=label), col='red', size=3)
                   aes(x=PC1*5e-7, y=PC2*5e-7, label=label), col='red', size=3)+
  labs(x=paste0('PC1 (', round(100*eigs_normalisedcounts[1]/sum(eigs_normalisedcounts), 2), '%)' ),
       y=paste0('PC2 (', round(100*eigs_normalisedcounts[2]/sum(eigs_normalisedcounts), 2), '%)' ))
ggsave("../figures/PCA_RNASeq/PCA_counts_mostvar_annotated.png", width = 8, height = 8)

## see which highly variable genes are part of the Tier 1 Cancer Gene Census (CGC) list
CGC_tier1 = read.table("../files/Census_allFri Feb 12 09 56 06 2021.tsv", sep = '\t', header = T)

ggplot(cbind.data.frame(prcomp_normalisedcounts$rotation[top_size_prcomp_normalisedcounts,1:2],
                        t2g[match(rownames(prcomp_normalisedcounts$rotation[top_size_prcomp_normalisedcounts,]),
                                        t2g$ensembl_gene_id),c('external_gene_name', 'entrezgene_id')]))+
  geom_segment(aes(x=0, xend=PC1, y=0, yend=PC2))+
  geom_label_repel(aes(x=PC1, y=PC2, label=ifelse(test = entrezgene_id %in% CGC_tier1$Entrez.GeneId,
                                                            yes = external_gene_name, no=NA)))+
  geom_label_repel(data=(cbind.data.frame(label=rownames(prcomp_normalisedcounts$x), prcomp_normalisedcounts$x[,1:2])),
             # aes(x=PC1*5e-5, y=PC2*5e-5, label=label), col='red', size=3)
             aes(x=PC1*5e-7, y=PC2*5e-7, label=label), col='red', size=3)
ggsave("../figures/PCA_RNASeq/PCA_counts_mostvar_CGC_tier1.png", width = 8, height = 8)

ggplot(cbind.data.frame(prcomp_normalisedcounts$rotation[top_size_prcomp_normalisedcounts2,1:2],
                        t2g[match(rownames(prcomp_normalisedcounts$rotation[top_size_prcomp_normalisedcounts2,]),
                                  t2g$ensembl_gene_id),c('external_gene_name', 'entrezgene_id')]))+
  geom_segment(aes(x=0, xend=PC1, y=0, yend=PC2))+
  geom_label_repel(aes(x=PC1, y=PC2, label=ifelse(test = entrezgene_id %in% CGC_tier1$Entrez.GeneId,
                                                  yes = external_gene_name, no=NA)))
ggsave("../figures/PCA_RNASeq/PCA_counts_mostvar_CGC_tier1_subset2.png", width = 8, height = 8)

## plot only loadings of these tier 1 genes
dim(naming_id_matched_normalisedcounts)
ggplot(cbind.data.frame(prcomp_normalisedcounts$rotation[naming_id_matched_normalisedcounts$entrezgene_id %in% CGC_tier1$Entrez.GeneId,1:2],
                        label=naming_id_matched_normalisedcounts$external_gene_name[naming_id_matched_normalisedcounts$entrezgene_id %in% CGC_tier1$Entrez.GeneId]))+
  geom_segment(aes(x=0, xend=PC1, y=0, yend=PC2))+
  geom_label_repel(aes(x=PC1, y=PC2, label=label))+
  geom_point(data=(cbind.data.frame(label=rownames(prcomp_normalisedcounts$x), prcomp_normalisedcounts$x[,1:2])),
                   # aes(x=PC1*5e-5, y=PC2*5e-5, label=label), col='red', size=3)
                   aes(x=PC1*1e-4, y=PC2*1e-4, label=label), col='red', size=3, max.overlaps = 18)+
  geom_label_repel(data=(cbind.data.frame(label=rownames(prcomp_normalisedcounts$x), prcomp_normalisedcounts$x[,1:2])),
                   # aes(x=PC1*5e-5, y=PC2*5e-5, label=label), col='red', size=3)
                   aes(x=PC1*1e-4, y=PC2*1e-4, label=label), col='red', size=3, max.overlaps = 18)
ggsave("../figures/PCA_RNASeq/PCA_counts_only_CGC_tier1.png", width = 8, height = 8)

## not used
goterm_invid = function(i){
  .x = i
  names(.x) = t2g[match(names(.x), t2g$ensembl_gene_id),'entrezgene_id']
  .x = as(.x[!is.na(names(.x))], 'matrix')
  colnames(.x) = 'Sample'
  res_gsva = gsva(expr = .x, gset.idx.list = c2BroadSets, min.sz=10, max.sz=500, verbose=TRUE)
  }

# library(topGO)
# topGO::runTest(object = rownames(prcomp_normalisedcounts_highlyVar$rotation)[cutree_on_pc1_and_2_NC_cos==1])
# topGOdata(rownames(prcomp_normalisedcounts_highlyVar$rotation)[cutree_on_pc1_and_2_NC_cos==1])
# sampleGOdata <- new("topGOdata",
#                     description = "Simple session", ontology = "BP",
#                     allGenes = rownames(prcomp_normalisedcounts_highlyVar$rotation),
#                     geneSel = rownames(prcomp_normalisedcounts_highlyVar$rotation)[cutree_on_pc1_and_2_NC_cos==1],
#                     nodeSize = 10,
#                     annot = annFUN.db, affyLib = affyLib)
# allGenes_fake0 = rownames(prcomp_normalisedcounts_highlyVar$rotation)
# allGenes_fake = rep(1, length(allGenes_fake0)); names(allGenes_fake) = allGenes_fake0
# topGenes_fake0 = rownames(prcomp_normalisedcounts_highlyVar$rotation)
# topGenes_fake = rep(1, length(topGenes_fake0)); names(allGenes_fake) = topGenes_fake0
# select_which_genes = function() rownames(prcomp_normalisedcounts_highlyVar$rotation)[cutree_on_pc1_and_2_NC_cos==1]
# sampleGOdata <- new("topGOdata",
#                     description = "Simple session", ontology = "BP",
#                     allGenes = allGenes_fake,
#                     geneSel = select_which_genes,
#                     nodeSize = 10,
#                     annot = annFUN.db, affyLib = affyLib)

## in the top 1500 most variable genes we have three clusters (judging from the PCA) and cancer genes in each cluster
## takes some time as it's 1500 genes
cosine_sim_NC_top = outer(1:nrow(prcomp_normalisedcounts$rotation[top_size_prcomp_normalisedcounts2,]),
                          1:nrow(prcomp_normalisedcounts$rotation[top_size_prcomp_normalisedcounts2,]),
                     Vectorize(function(i,j) cosine(prcomp_normalisedcounts$rotation[top_size_prcomp_normalisedcounts2,][i,],
                                                    prcomp_normalisedcounts$rotation[top_size_prcomp_normalisedcounts2,][j,])))
cosine_sim_NC_top =  acos(cosine_sim_NC_top)/pi
hclust_on_pc1_and_2_NC_cos_top = hclust(as.dist(cosine_sim_NC_top))

cutree_on_pc1_and_2_NC_cos_top = cutree(hclust_on_pc1_and_2_NC_cos, k = 10)

dim(prcomp_normalisedcounts$rotation[top_size_prcomp_normalisedcounts,])


library(GOSim)
goterm_per_cat = mclapply(unique(cutree_on_pc1_and_2_NC_cos), function(i){
  require(topGO)
  .x = rownames(prcomp_normalisedcounts_highlyVar$rotation)[cutree_on_pc1_and_2_NC_cos==i]
  .xentrez = t2g[match(.x, t2g$ensembl_gene_id),'entrezgene_id']
  .xentrez = .xentrez[!is.na(.xentrez)]
  
  .xres_gosimenrich = GOSim::GOenrichment(genesOfInterest = as.character(.xentrez), allgenes = as.character(entrezgene_id_matched))
  .xres_gosim = GOSim::getGOInfo(.xentrez)
  # .xres_gosimenrich$GOTerms[.xres_gosimenrich$GOTerms$go_id == names(sort(.xres_gosimenrich$p.values)[1]),'Term']
  .xres_gosimenrich
})
names(goterm_per_cat) = unique(cutree_on_pc1_and_2_NC_cos)
data_goterms_top = cbind.data.frame(x=t(sapply(unique(cutree_on_pc1_and_2_NC_cos), function(i) colMeans(subset(prcomp_normalisedcounts_highlyVar$rotation[,1:2], cutree_on_pc1_and_2_NC_cos == i)))),
                 go=sapply(goterm_per_cat, function(j){
                   j$GOTerms[j$GOTerms$go_id == names(sort(j$p.values)[1]),'Term']
                 }))
# data_goterms_top$go = paste0(data_goterms_top$go, 1:nrow(data_goterms_top))
data_goterms_top$go = sapply(as.character(data_goterms_top$go), give_linebreaks, max_line = 16)

ggplot()+
  geom_segment(data = cbind.data.frame(prcomp_normalisedcounts_highlyVar$rotation[,1:2], col=cutree_on_pc1_and_2_NC_cos),
               aes(x=0, xend=PC1, y=0, yend=PC2, col=factor(col)))+
  # xlim(-0.2, 0.7)+
  geom_label_repel(data=(cbind.data.frame(label=rownames(prcomp_normalisedcounts$x), prcomp_normalisedcounts$x[,1:2])),
                   # aes(x=PC1*5e-5, y=PC2*5e-5, label=label), col='red', size=3)
                   aes(x=PC1*1e-4, y=PC2*1e-4, label=label), col='red', size=3, max.overlaps = 18)+
  geom_label_repel(data = data_goterms_top,
                   aes(x=x.PC1*2, y=x.PC2*2, label=go))+labs(x='PC1', y='PC2')+
  ggtitle('Go terms for clustered loadings')
ggsave("../figures/PCA_RNASeq/grouped_loadings_GO.pdf", width = 10, height = 10)


GOinfo_NC = getGOInfo(naming_id_matched_normalisedcounts$entrezgene_id)
GOinfo_NC1_goterms = GOinfo_NC[1,][(match(naming_id_matched_normalisedcounts$entrezgene_id, names(GOinfo_NC[1,])))]

go_term_list =  c('GO:0097194'= 'execution phase of apoptosis',
                  'GO:0070227'= ' lymphocyte apoptotic process',
                  'GO:0005125' = 'cytokine activity', 
                  #'GO:0006357'= 'regulation of transcription by RNA polymerase II',
                  #'GO:0045944' = 'positive regulation of transcription by RNA polymerase II',
                  'GO:0006915'= 'apoptotic process',
                  'GO:0030154' = 'cell differentiation',
                  'GO:0043312' = 'neutrophil degranulation',
                  'GO:0008284' = 'positive regulation of cell population proliferation',
                  'GO:0043066' = 'negative regulation of apoptotic process'
)
# GO:0007186 G protein-coupled receptor signaling pathway
GOinfo_NC_annotated = lapply(GOinfo_NC1_goterms, function(i) go_term_list[match(names(go_term_list), i)] )
GOinfo_NC_annotated = sapply(GOinfo_NC_annotated, function(i) i[!is.na(i)][1])
table(unlist(GOinfo_NC_annotated))
head(sort(table(unlist(GOinfo_NC1_goterms)), d=T), n=20)

ggplot(cbind.data.frame(prcomp_normalisedcounts$rotation[,1:2],
                        naming_id_matched_normalisedcounts$ensembl_gene_id,
                        go_term=GOinfo_NC_annotated))+
  geom_segment(aes(x=0, xend=PC1, y=0, yend=PC2, col=go_term))+
  geom_label_repel(data=(cbind.data.frame(label=rownames(prcomp_normalisedcounts$x), prcomp_normalisedcounts$x[,1:2])),
                   # aes(x=PC1*5e-5, y=PC2*5e-5, label=label), col='red', size=3)
                   aes(x=PC1*5e-7, y=PC2*5e-7, label=label), col='red', size=3)+facet_wrap(.~go_term)
ggsave("../figures/PCA_RNASeq/PCA_counts_loadings_goterms.png", width = 12, height = 8)

## genes interesting for ovarian cancer
subset_genes = c('MYC', 'CCNE1', 'PIK3CA', 'TERT', 'KRAS', 'PTEN', 'RB1', 'AKT1',
                 'AKT2', 'PARP1', 'PARP2', 'ATM', 'ATR', 'WEE1', 'TOP1', 'TUBB1')
t2g$ensembl_gene_id[match(subset_genes, t2g$external_gene_name)]
ggplot(cbind.data.frame(prcomp_normalisedcounts$rotation[t2g$ensembl_gene_id[match(subset_genes, t2g$external_gene_name)],1:2],
                        label=subset_genes))+
  geom_segment(aes(x=0, xend=PC1, y=0, yend=PC2))+
  geom_label_repel(aes(x=PC1, y=PC2, label=label))+
  geom_point(data=(cbind.data.frame(label=rownames(prcomp_normalisedcounts$x), prcomp_normalisedcounts$x[,1:2])),
             # aes(x=PC1*5e-5, y=PC2*5e-5, label=label), col='red', size=3)
             aes(x=PC1*1e-4, y=PC2*1e-4, label=label), col='red', size=3, max.overlaps = 18)+
  geom_label_repel(data=(cbind.data.frame(label=rownames(prcomp_normalisedcounts$x), prcomp_normalisedcounts$x[,1:2])),
                   # aes(x=PC1*5e-5, y=PC2*5e-5, label=label), col='red', size=3)
                   aes(x=PC1*1e-4, y=PC2*1e-4, label=label), col='red', size=3, max.overlaps = 18)
ggsave("../figures/PCA_RNASeq/PCA_counts_genes_of_interest.png", width = 8, height = 8)

## look at the weird categories
# lapply(unique(cutree_on_pc1_and_2)[table(cutree_on_pc1_and_2) < 5], function(i){
#   # sapply(cutree_on_pc1_and_2[cutree_on_pc1_and_2 == i], function(j) TPM$gene_name[TPM$gene_id == j] )
#   TPM[match(names(cutree_on_pc1_and_2[cutree_on_pc1_and_2 == i]), TPM$gene_id), c('gene_name', 'gene_id')]
# })
# 
# lapply(unique(cutree_on_pc1_and_2_k5)[table(cutree_on_pc1_and_2_k5) < 5], function(i){
#   TPM[match(names(cutree_on_pc1_and_2_k5[cutree_on_pc1_and_2_k5 == i]), TPM$gene_id), c('gene_name', 'gene_id')]
# })

# table(cutree_on_pc1_and_2)

names(cutree_on_pc1_and_2[cutree_on_pc1_and_2 == 2])

df = cbind.data.frame(prcomp_res$x[,1:2], names=labels)
max_loadings = 10
selected_genes = sort(apply(prcomp_res$rotation[,1:2], 1, function(i) sqrt(sum(i**2))), decreasing = T)[1:max_loadings]
subset_loadings = prcomp_res$rotation[,1:2][rownames(prcomp_res$rotation) %in% names(selected_genes),]
TPM[match(rownames(subset_loadings), TPM$gene_id),1:3]
ggplot(data=df, aes(x=PC1, y=PC2, label=names))+
  geom_point(size=size_points)+
  geom_label_repel()+
  geom_segment(data=cbind.data.frame(subset_loadings,
                                     names=rownames(subset_loadings)),
               aes(x=0, y=0, xend=PC1*1e6, yend=PC2*1e6), col='red',
               arrow = arrow(length = unit(0.03, "npc")))+
  labs(x=paste0('PC1 (', round(100*eigs[1]/sum(eigs), 2), '%)' ),
       y=paste0('PC2 (', round(100*eigs[2]/sum(eigs), 2), '%)' ))+ggtitle('PCA with loadings')+
  theme_bw()
ggsave("../figures/PCA_RNASeq/PCA_TPM_loadings.pdf", width = 6.5, height = 5)
#--------------------------------------------------------------------------------#

#--------------------------------------------------------------------------------#
colnames(TPM_data)
rownames(TPM_data)

any(duplicated(rownames(TPM_data)))
any(duplicated(colnames(TPM_data)))

TPM_data = TPM_data[-which(duplicated(TPM_data)),]
#--------------------------------------------------------------------------------#

#--------------------------------------------------------------------------------#
library(M3C)
perplex = 2
tsne_result = tsne(TPM_data, perplex = perplex)
ggplot(data.frame(tsne_result$data, label=gsub('JBLAB', '', colnames(TPM_data)),
                  TPM_p53=as.numeric(TPM_data[rownames(TPM_data) == 'ENSG00000141510',])),
       aes(x=X1, y=X2, col=TPM_p53, label=label))+
  geom_point()+
  geom_label_repel(size=3)+ggtitle(paste0('TSNE with perplexity=',perplex ))
ggsave(paste0("../figures/PCA_RNASeq/TSNE_TPM_perplexity", perplex, ".pdf"), width = 6.5, height = 5)
#--------------------------------------------------------------------------------#

#--------------------------------------------------------------------------------#
save.image("../objects/image_RNASeq_organoids.RData")
#--------------------------------------------------------------------------------#

#--------------------------------------------------------------------------------#
## Without PDO17
prcomp_normalisedcounts_noPDO17 = prcomp(t(normalised_counts[,!grepl('PDO17', colnames(normalised_counts))]), scale. = FALSE, center = TRUE)
eigs_prcomp_normalisedcounts_noPDO17 <- prcomp_normalisedcounts_noPDO17$sdev^2
ggplot(data=cbind.data.frame(prcomp_normalisedcounts_noPDO17$x[,1:2], names=rownames(prcomp_normalisedcounts_noPDO17$x)),
       aes(x=PC1, y=PC2, label=names))+
  geom_point(size=size_points)+
  geom_label_repel()+
  labs(x=paste0('PC1 (', round(100*eigs_prcomp_normalisedcounts_noPDO17[1]/sum(eigs_prcomp_normalisedcounts_noPDO17), 2), '%)' ),
       y=paste0('PC2 (', round(100*eigs_prcomp_normalisedcounts_noPDO17[2]/sum(eigs_prcomp_normalisedcounts_noPDO17), 2), '%)' ))+ggtitle('PCA')
ggsave("../figures/PCA_RNASeq/PCA_counts_wihout_PDO17.pdf", width = 6, height = 6)
#--------------------------------------------------------------------------------#


highlyVar_TPM = names(sort(apply(TPM_data, 1, var), decreasing = T)[1:500])
prcomp_TPM_highlyVar = prcomp_res
prcomp_TPM_highlyVar$rotation = prcomp_TPM_highlyVar$rotation[rownames(prcomp_TPM_highlyVar$rotation) %in% highlyVar_TPM,]
size_prcomp_TPM_highlyVar = apply(prcomp_TPM_highlyVar$rotation[,1:2], 1, function(i) sqrt(sum(i**2)))
top_size_size_prcomp_TPM_highlyVar = order(size_prcomp_TPM_highlyVar, decreasing = T)[1:250]
ggplot()+
  geom_segment(data=cbind.data.frame(names=t2g$external_gene_name[match(rownames(prcomp_TPM_highlyVar$rotation[top_size_size_prcomp_TPM_highlyVar,1:2]),
                                                                        t2g$ensembl_gene_id)],
                                     prcomp_TPM_highlyVar$rotation[top_size_size_prcomp_TPM_highlyVar,1:2]),
               aes(x=0, xend=PC1, y=0, yend=PC2, label=names))+
  geom_label_repel(data=cbind.data.frame(names=t2g$external_gene_name[match(rownames(prcomp_TPM_highlyVar$rotation[top_size_size_prcomp_TPM_highlyVar,1:2]),
                                                                            t2g$ensembl_gene_id)],
                                         prcomp_TPM_highlyVar$rotation[top_size_size_prcomp_TPM_highlyVar,1:2]),
                   aes(x=PC1, y=PC2, label=names))+
  geom_label_repel(data=(cbind.data.frame(label=rownames(prcomp_TPM_highlyVar$x), prcomp_TPM_highlyVar$x[,1:2])),
                   # aes(x=PC1*5e-5, y=PC2*5e-5, label=label), col='red', size=3)
                   aes(x=PC1*1e-7, y=PC2*1e-7, label=label), col='red', size=3)+
  labs(x=paste0('PC1 (', round(100*eigs[1]/sum(eigs), 2), '%)' ),
       y=paste0('PC2 (', round(100*eigs[2]/sum(eigs), 2), '%)' ))+ggtitle('PCA')
ggsave("../figures/PCA_RNASeq/PCA_TPM_loadings_highlyVar.pdf", height = 10, width = 10)

library(DESeq2)
DESeq2::plotPCA(object = normalised_counts)


apply(log(normalised_counts), 1, mean)

exp(sweep(log(normalised_counts), 1, rowSums(log(normalised_counts)), '/'))

log(normalised_counts)

exp(sweep(log(normalised_counts), 1, rowSums(log(normalised_counts)), ''))

a = matrix(rpois(n = 100, 30), nrow=20)
## geometric means
geom_means = apply(a, 1, function(i) exp(mean(log(i))) )
exp(sweep(log(a), 1, rowMeans(log(a)), '/' ))

normalised_counts_no_zero = normalised_counts
normalised_counts_no_zero = normalised_counts_no_zero[rowSums(normalised_counts_no_zero)>0,]
# normalised_counts_no_zero = normalised_counts_no_zero[! apply(normalised_counts_no_zero, 1, function(i) sum(i > 0) <3 ),]
x_prcomp_normalisedcounts = prcomp(t(normalised_counts_no_zero), scale. = TRUE, center = TRUE)




