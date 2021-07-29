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
require(pheatmap)
library(ReactomePA)
library(DESeq2)
library(lsa) ## for cosine similarity
source("other_scripts/functions.R")

gene_to_ensembl <- function(genename){
  t2g_GRCh38$ensembl_gene_id[match(genename, t2g_GRCh38$external_gene_name)]
}

#--------------------------------------------------------------------------------#

#--------------------------------------------------------------------------------#
raw_counts0 = read.csv("../../RNASeq_and_CN/20191218_ViasM_BJ_orgaBrs/RnaSeqPip/counts/counts_raw_subsetno3pbias.csv", stringsAsFactors = FALSE)
rownames(raw_counts0) = raw_counts0[,1]; raw_counts0 <- raw_counts0[,-1]
normalised_counts = read.csv("../files/counts_norm.csv", stringsAsFactors = FALSE)
rownames(normalised_counts) = normalised_counts[,1]
normalised_counts = normalised_counts[,-1]
#--------------------------------------------------------------------------------#

#--------------------------------------------------------------------------------#
# Renaming
renaming1 = read_excel("../files/PDOnameProperSample_sWGS_RNAseq.xlsx")

renaming_NC = renaming1[match( gsub('[.]', '', colnames(normalised_counts)), gsub('-', '', renaming1$sampleNameRNAseq)),]

colnames(normalised_counts)[!is.na(renaming_NC$PDO)] = renaming_NC$PDO[!is.na(renaming_NC$PDO)]
colnames(normalised_counts)

## Remove normal samples
normalised_counts = normalised_counts[,!grepl('FT', colnames(normalised_counts))]
normalised_counts = normalised_counts[rowSums(normalised_counts)>0,]
#--------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------#
## 3' bias
bias3prime = read_xlsx("../files/PDOnameProperSample_sWGS_RNAseq_3bias.xlsx")
#--------------------------------------------------------------------------------#

#--------------------------------------------------------------------------------#
prcomp_normalisedcounts = prcomp(t(normalised_counts), scale. = TRUE, center = TRUE)
eigs_normalisedcounts <- prcomp_normalisedcounts$sdev^2
round(eigs_normalisedcounts/sum(eigs_normalisedcounts)*100)

ggplot(cbind.data.frame(prcomp_normalisedcounts$x,
                        col=bias3prime$`3'Bias`[match(rownames(prcomp_normalisedcounts$x),bias3prime$PDO)],
                        label=rownames(prcomp_normalisedcounts$x)),
       aes(x=PC1, y=PC2, col=col, label=label))+
  geom_label_repel()+
  geom_point()
ggsave("../figures/PCA_RNASeq/3primebias_PCA_counts.pdf", width = 6.5, height = 5)
ggplot(cbind.data.frame(prcomp_normalisedcounts$x,
                        col=bias3prime$`3'Bias`[match(rownames(prcomp_normalisedcounts$x),bias3prime$PDO)],
                        label=rownames(prcomp_normalisedcounts$x),
                        pc1=prcomp_normalisedcounts$x[,1]),
       aes(x=factor(label, levels=bias3prime$PDO[order(bias3prime$`3'Bias`)]),
           y=col, fill=pc1))+
  geom_bar(stat = "identity")
ggsave("../figures/PCA_RNASeq/3primebias_barplot_counts.pdf", width = 9.5, height = 5)
rm(prcomp_normalisedcounts)
#--------------------------------------------------------------------------------#

#--------------------------------------------------------------------------------#
## Re-normalise counts

colnames(raw_counts0) <- renaming1$PDO[match(gsub('[.]', '-', colnames(raw_counts0)), renaming1$sampleNameRNAseq)]
raw_counts0 = raw_counts0[,!(colnames(raw_counts0) %in% c('PDO14', 'PDO16', 'PDO18',
                                                                            'PDO13', 'PDO4', 'PDO9',
                                                                            'PDO17'))]
## in CN/GE script we have already re-derived the normalised DESeq2 counts with non-3' bias samples (PDO and FT)
renormalised_counts <- readRDS("../../RNASeq_and_CN/20191218_ViasM_BJ_orgaBrs/output/fig3_renormalised_counts_obj_11orgs.RDS")
renormalised_counts <- DESeq2::counts(renormalised_counts, normalized=T)
## keeping only PDO
renormalised_counts <- renormalised_counts[,grepl('PDO', colnames(renormalised_counts))]
ncol(renormalised_counts) ## we have 11 organoids
# renormalised_counts <- raw_counts0
# renormalised_counts <- DESeqDataSetFromMatrix(countData = renormalised_counts,
#                               colData = cbind.data.frame(Sample=colnames(renormalised_counts),
#                                                          Group=1),
#                               design = ~ 1)
# renormalised_counts <- estimateSizeFactors(renormalised_counts)
# renormalised_counts <- counts(renormalised_counts, normalized=T)
normalised_counts <- renormalised_counts

#--------------------------------------------------------------------------------#

#--------------------------------------------------------------------------------#
## Removing samples with 3' bias

normalised_counts = normalised_counts[,!(colnames(normalised_counts) %in% c('PDO14', 'PDO16', 'PDO18',
                                                                            'PDO13', 'PDO4', 'PDO9',
                                                                            'PDO17'))]
normalised_counts = normalised_counts[rowSums(normalised_counts)>0,]
prcomp_normalisedcounts = prcomp(t(normalised_counts), scale. = TRUE, center = TRUE)
#--------------------------------------------------------------------------------#

#--------------------------------------------------------------------------------#
hclust_rnaseq <- (hclust(dist(t(normalised_counts))))
saveRDS(hclust_rnaseq, "../objects/hclust_rnaseq.RDS")
#--------------------------------------------------------------------------------#

#--------------------------------------------------------------------------------#
## TAMSeq
TAMSeq = read_xlsx("../files/SupplementaryTable1.xlsx")
TAMSeq[,'org'] = renaming1$PDO[match(TAMSeq$name, renaming1$ID)]
TAMSeq = dcast(TAMSeq[,c(8, 20)], org~`Symbol (Gene ID)`)
TAMSeq = TAMSeq[match(colnames(normalised_counts), TAMSeq$org),]
#--------------------------------------------------------------------------------#

#--------------------------------------------------------------------------------#
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
# genes <- TPM$gene_id

t2g <- getBM(
  attributes = c('ensembl_gene_id', 'entrezgene_id', 'external_gene_name', 'chromosome_name'),
  values = rownames(normalised_counts),
  filter = 'ensembl_gene_id',
  mart = mart, useCache = FALSE)

entrezgene_id_matched_normalisedcounts = t2g[match(rownames(normalised_counts), t2g$ensembl_gene_id),'entrezgene_id']
naming_id_matched_normalisedcounts = t2g[match(rownames(normalised_counts), t2g$ensembl_gene_id),]

t2g_GRCh38 <- readRDS("../../copy_number_analysis_organoids/robjects/t2g2.RDS")
coding <- readRDS("../../copy_number_analysis_organoids/robjects/coding_genes.RDS")

#--------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------#
data(c2BroadSets) ## from GSVAdata

#--------------------------------------------------------------------------------#
## same with normalised counts
TPM_data_entrez_NC = normalised_counts[!is.na(entrezgene_id_matched_normalisedcounts),]
entrezgene_id_matched_normalisedcounts = entrezgene_id_matched_normalisedcounts[!is.na(entrezgene_id_matched_normalisedcounts)]
TPM_data_entrez_NC = TPM_data_entrez_NC[!duplicated(entrezgene_id_matched_normalisedcounts),]
entrezgene_id_matched_normalisedcounts = entrezgene_id_matched_normalisedcounts[!duplicated(entrezgene_id_matched_normalisedcounts)]
rownames(TPM_data_entrez_NC) = entrezgene_id_matched_normalisedcounts
results_Gsva_counts <- gsva(expr = as(TPM_data_entrez_NC, 'matrix'), gset.idx.list = c2BroadSets, min.sz=10, max.sz=500, verbose=TRUE)
pca_with_gsva_annotation_NC = cbind.data.frame(prcomp_normalisedcounts$x, t(results_Gsva_counts[,match(colnames(results_Gsva_counts), rownames(prcomp_normalisedcounts$x))]),
                                               labels=rownames(prcomp_normalisedcounts$x),
                                               TAMSeq)

ggplot(pca_with_gsva_annotation_NC, aes(x=PC1, y=PC2, col=KAUFFMANN_DNA_REPAIR_GENES, label=labels))+
  geom_point()+
  geom_label_repel()+
  scale_color_jcolors_contin("pal3", reverse = TRUE)+labs(col = "Threshold", shape="Amount")+
  ggtitle('KAUFFMANN_DNA_REPAIR_GENES')
ggsave("../figures/PCA_RNASeq/PCA_counts_subset_KAUFFMANN_DNA_REPAIR_GENES.png", width = 6.5, height = 5)

saveRDS(pca_with_gsva_annotation_NC, "../objects/fig4_pca_with_gsva_annotation_NC.RDS")
ggplot(pca_with_gsva_annotation_NC, aes(x=PC1, y=PC2, col=OUELLET_CULTURED_OVARIAN_CANCER_INVASIVE_VS_LMP_UP, label=labels))+
  geom_point()+
  geom_label_repel()+
  scale_color_jcolors_contin("pal3", reverse = TRUE, bias = 2.25)+labs(col = "Threshold", shape="Amount")+
  ggtitle('OUELLET_CULTURED_OVARIAN_CANCER_INVASIVE_VS_LMP_UP')
ggsave("../figures/PCA_RNASeq/PCA_counts_subset_OUELLET_CULTURED_OVARIAN_CANCER_INVASIVE_VS_LMP_UP.png", width = 6.5, height = 5)

ggplot(pca_with_gsva_annotation_NC, aes(x=PC1, y=PC2, col=factor(BRCA1), label=labels))+
  geom_point()+
  geom_label_repel()+
  ggtitle('Number of mutations in BRCA (TAMSeq)')
ggsave("../figures/PCA_RNASeq/PCA_counts_subset_BRCA_TAMSeq.png", width = 6.5, height = 5)

results_Gsva_counts_only_variable = results_Gsva_counts[order(apply(results_Gsva_counts, 1, var), decreasing = T)[1:100],]
pdf("../figures/PCA_RNASeq/heatmap_gsva_counts_subset.pdf", height = 20, width = 18)
pheatmap::pheatmap(results_Gsva_counts_only_variable)
dev.off()
pdf("../figures/PCA_RNASeq/heatmap_gsva_subset_counts_all.pdf", height = 80, width = 30)
pheatmap::pheatmap(results_Gsva_counts)
dev.off()
pdf("../figures/PCA_RNASeq/heatmap_gsva_counts_subset_ovarian.pdf", height = 6, width = 15)
pheatmap::pheatmap(results_Gsva_counts[grepl('OVAR', rownames(results_Gsva_counts)),])
dev.off()


#--------------------------------------------------------------------------------#
# Very time consuming
## get only the most variable genes
highlyVar = names(sort(apply(normalised_counts, 1, var), decreasing = T)[1:500])
prcomp_normalisedcounts_highlyVar = prcomp_normalisedcounts
prcomp_normalisedcounts_highlyVar$rotation = prcomp_normalisedcounts_highlyVar$rotation[rownames(prcomp_normalisedcounts_highlyVar$rotation) %in% highlyVar,]
size_prcomp_normalisedcounts_highlyVar = apply(prcomp_normalisedcounts_highlyVar$rotation[,1:2], 1, function(i) sqrt(sum(i**2)))
top_size_size_prcomp_normalisedcounts_highlyVar = order(size_prcomp_normalisedcounts_highlyVar, decreasing = T)[1:250]
top_size_size_prcomp_normalisedcounts_highlyVar2 = order(size_prcomp_normalisedcounts_highlyVar, decreasing = T)[1:25]
# ggplot()+
#   geom_segment(data=cbind.data.frame(names=t2g$external_gene_name[match(rownames(prcomp_normalisedcounts_highlyVar$rotation[top_size_size_prcomp_normalisedcounts_highlyVar2,1:2]), t2g$ensembl_gene_id)],
#                                      prcomp_normalisedcounts_highlyVar$rotation[top_size_size_prcomp_normalisedcounts_highlyVar2,1:2]),
#                aes(x=0, xend=PC1, y=0, yend=PC2, label=names))+
#   geom_label_repel(data=cbind.data.frame(names=t2g$external_gene_name[match(rownames(prcomp_normalisedcounts_highlyVar$rotation[top_size_size_prcomp_normalisedcounts_highlyVar2,1:2]), t2g$ensembl_gene_id)],
#                                          prcomp_normalisedcounts_highlyVar$rotation[top_size_size_prcomp_normalisedcounts_highlyVar2,1:2]),
#                    aes(x=PC1, y=PC2, label=names))

## See the genes with the highest loadings, even if they are not the most variable
size_prcomp_normalisedcounts = apply(prcomp_normalisedcounts$rotation[,1:2], 1, function(i) sqrt(sum(i**2)))
top_size_size_prcomp_normalisedcounts = order(size_prcomp_normalisedcounts, decreasing = T)[1:250]
ggplot()+
  geom_segment(data=cbind.data.frame(names=t2g$external_gene_name[match(rownames(prcomp_normalisedcounts$rotation[top_size_size_prcomp_normalisedcounts,1:2]), t2g$ensembl_gene_id)],
                                     prcomp_normalisedcounts$rotation[top_size_size_prcomp_normalisedcounts,1:2]),
               aes(x=0, xend=PC1, y=0, yend=PC2, label=names), alpha=0.3)+
  geom_point(data=(cbind.data.frame(label=rownames(prcomp_normalisedcounts$x), prcomp_normalisedcounts$x[,1:2])),
           aes(x=PC1*1e-4, y=PC2*1e-4), col='red', size=3, alpha=0.7)+
  geom_label_repel(data=(cbind.data.frame(label=rownames(prcomp_normalisedcounts$x), prcomp_normalisedcounts$x[,1:2])),
                   aes(x=PC1*1e-4, y=PC2*1e-4, label=label), col='red', size=3)
ggsave("../figures/PCA_RNASeq/PCA_counts_subset_all_loadings_top_size.png", width = 8, height = 8)

ggplot()+
  geom_segment(data=cbind.data.frame(names=t2g$external_gene_name[match(rownames(prcomp_normalisedcounts$rotation[,1:2]), t2g$ensembl_gene_id)],
                                     prcomp_normalisedcounts$rotation[,1:2]),
               aes(x=0, xend=PC1, y=0, yend=PC2, label=names), alpha=0.02)+
  geom_point(data=(cbind.data.frame(label=rownames(prcomp_normalisedcounts$x), prcomp_normalisedcounts$x[,1:2])),
             aes(x=PC1*1e-4, y=PC2*1e-4), col='red', size=3)+
  geom_label_repel(data=(cbind.data.frame(label=rownames(prcomp_normalisedcounts$x), prcomp_normalisedcounts$x[,1:2])),
                   aes(x=PC1*1e-4, y=PC2*1e-4, label=label), col='red', size=3)
ggsave("../figures/PCA_RNASeq/PCA_counts_subset_all_loadings.png", width = 8, height = 8)

## cosine similarity of PC1, PC2
cosine_sim_NC= outer(1:nrow(prcomp_normalisedcounts_highlyVar$rotation), 1:nrow(prcomp_normalisedcounts_highlyVar$rotation),
                     Vectorize(function(i,j) cosine(prcomp_normalisedcounts_highlyVar$rotation[,1:2][i,], prcomp_normalisedcounts_highlyVar$rotation[,1:2][j,])))
cosine_dist_NC =  acos(cosine_sim_NC)/pi
hclust_on_pc1_and_2_NC_cos = hclust(as.dist(cosine_dist_NC))
cutree_on_pc1_and_2_NC_cos = cutree(hclust_on_pc1_and_2_NC_cos, k = 4)
cutree_on_pc1_and_2_NC_cos = cutree(hclust_on_pc1_and_2_NC_cos, k = 10)


ggplot(cbind.data.frame(prcomp_normalisedcounts_highlyVar$rotation[,1:2]),
       aes(x=0, xend=PC1, y=0, yend=PC2))+
  geom_segment()

ggplot(cbind.data.frame(prcomp_normalisedcounts_highlyVar$rotation[,1:2], col=cutree_on_pc1_and_2_NC_cos),
       aes(x=0, xend=PC1, y=0, yend=PC2, col=factor(col)))+
  geom_segment()

## normalise to same length
normalise_magnitude = function(i){
  i/sqrt(sum(i**2))
}

prcomp_normalisedcounts_highlyVar_norm_rotation = t(apply(prcomp_normalisedcounts_highlyVar$rotation[,1:2], 1, normalise_magnitude))
ggplot(cbind.data.frame(prcomp_normalisedcounts_highlyVar_norm_rotation[,1:2]),
       aes(x=0, xend=PC1, y=0, yend=PC2))+
  geom_segment(alpha=0.2)

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
ggsave("../figures/PCA_RNASeq/PCA_counts_subset_mostvar_annotated.png", width = 8, height = 8)

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
ggsave("../figures/PCA_RNASeq/PCA_counts_subset_mostvar_CGC_tier1.png", width = 8, height = 8)

ggplot(cbind.data.frame(prcomp_normalisedcounts$rotation[top_size_prcomp_normalisedcounts2,1:2],
                        t2g[match(rownames(prcomp_normalisedcounts$rotation[top_size_prcomp_normalisedcounts2,]),
                                  t2g$ensembl_gene_id),c('external_gene_name', 'entrezgene_id')]))+
  geom_segment(aes(x=0, xend=PC1, y=0, yend=PC2))+
  geom_label_repel(aes(x=PC1, y=PC2, label=ifelse(test = entrezgene_id %in% CGC_tier1$Entrez.GeneId,
                                                  yes = external_gene_name, no=NA)))
ggsave("../figures/PCA_RNASeq/PCA_counts_subset_mostvar_CGC_tier1_subset2.png", width = 8, height = 8)

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
ggsave("../figures/PCA_RNASeq/PCA_counts_subset_only_CGC_tier1.png", width = 8, height = 8)

## not used
goterm_invid = function(i){
  .x = i
  names(.x) = t2g[match(names(.x), t2g$ensembl_gene_id),'entrezgene_id']
  .x = as(.x[!is.na(names(.x))], 'matrix')
  colnames(.x) = 'Sample'
  res_gsva = gsva(expr = .x, gset.idx.list = c2BroadSets, min.sz=10, max.sz=500, verbose=TRUE)
  return(res_gsva)
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
give_goterm_list_of_genes <- function(i){
  require(topGO)
  .x = rownames(prcomp_normalisedcounts_highlyVar$rotation)[cutree_on_pc1_and_2_NC_cos==i]
  .xentrez = t2g[match(.x, t2g$ensembl_gene_id),'entrezgene_id']
  .xentrez = .xentrez[!is.na(.xentrez)]
  
  .xres_gosimenrich = GOSim::GOenrichment(genesOfInterest = as.character(.xentrez), allgenes = as.character(entrezgene_id_matched_normalisedcounts))
  .xres_gosim = GOSim::getGOInfo(.xentrez)
  # .xres_gosimenrich$GOTerms[.xres_gosimenrich$GOTerms$go_id == names(sort(.xres_gosimenrich$p.values)[1]),'Term']
  .xres_gosimenrich
}
goterm_per_cat = mclapply(unique(cutree_on_pc1_and_2_NC_cos), give_goterm_list_of_genes)
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
ggsave("../figures/PCA_RNASeq/grouped_loadings_GO_subset.pdf", width = 10, height = 10)

lapply(goterm_per_cat, function(j) j$GOTerms$Term[1:20])

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
ggsave("../figures/PCA_RNASeq/PCA_counts_subset_loadings_goterms.png", width = 12, height = 8)

## Looking at other particular go terms
GOinfo_NC_annotated_matched <- GOinfo_NC1_goterms[match(rownames(prcomp_normalisedcounts$rotation), naming_id_matched_normalisedcounts$external_gene_name)]
dim(prcomp_normalisedcounts$rotation)
length(GOinfo_NC_annotated_matched)

proliferation_0008283 <- sapply(GOinfo_NC_annotated_matched, function(i) "GO:0008283" %in% i)
HRD_0000278 <- sapply(GOinfo_NC_annotated_matched, function(i) "GO:0006303" %in% i)
apoptosisp53_0000278 <- sapply(GOinfo_NC_annotated_matched, function(i) "GO:0072332" %in% i)
hypoxia_0000278 <- sapply(GOinfo_NC_annotated_matched, function(i) "GO:0001666" %in% i)
cellcycle_0000278 <- sapply(GOinfo_NC_annotated_matched, function(i) "GO:0000278" %in% i)
apoptosis_0006915 <- sapply(GOinfo_NC_annotated_matched, function(i) "GO:0006915" %in% i)
apoptosis_0006915 <- sapply(GOinfo_NC_annotated_matched, function(i) "GO:0009653" %in% i)
table(apoptosis_0006915)

df_goterm <- cbind.data.frame(gene=rownames(prcomp_normalisedcounts$rotation[,1:2]),
                 prcomp_normalisedcounts$rotation[,1:2],
                 go_term=apoptosis_0006915)[apoptosis_0006915,]
ggplot()+
  geom_segment(data = cbind.data.frame(gene=rownames(prcomp_normalisedcounts$rotation[,1:2]),
                                       prcomp_normalisedcounts$rotation[,1:2]),
               aes(x=0, xend=PC1, y=0, yend=PC2), alpha=0.3)+
  geom_segment(data = df_goterm, aes(x=0, xend=PC1, y=0, yend=PC2, col=go_term))+
  geom_label_repel(data = df_goterm, aes(x=PC1, y=PC2, label=gene), size=3)

## genes interesting for ovarian cancer
subset_genes = c('MYC', 'CCNE1', 'PIK3CA', 'TERT', 'KRAS', 'PTEN', 'RB1', 'AKT1',
                 'AKT2', 'PARP1', 'PARP2', 'ATM', 'ATR', 'WEE1', 'TOP1', 'TUBB1')
t2g$ensembl_gene_id[match(subset_genes, t2g$external_gene_name)]
ggplot(cbind.data.frame(prcomp_normalisedcounts$rotation[t2g$ensembl_gene_id[match(subset_genes, t2g$external_gene_name)],1:2],
                        label=subset_genes))+
  geom_segment(aes(x=0, xend=PC1, y=0, yend=PC2),  arrow = arrow(length = unit(0.02, "npc")))+
  geom_label_repel(aes(x=PC1, y=PC2, label=label))+
  geom_point(data=(cbind.data.frame(label=rownames(prcomp_normalisedcounts$x), prcomp_normalisedcounts$x[,1:2])),
             # aes(x=PC1*5e-5, y=PC2*5e-5, label=label), col='red', size=3)
             aes(x=PC1*1e-4, y=PC2*1e-4, label=label), col='red', size=3, max.overlaps = 18)+
  geom_label_repel(data=(cbind.data.frame(label=rownames(prcomp_normalisedcounts$x), prcomp_normalisedcounts$x[,1:2])),
                   # aes(x=PC1*5e-5, y=PC2*5e-5, label=label), col='red', size=3)
                   aes(x=PC1*1e-4, y=PC2*1e-4, label=label), col='red', size=3, max.overlaps = 18)
ggsave("../figures/PCA_RNASeq/PCA_counts_subset_genes_of_interest.png", width = 8, height = 8)
#--------------------------------------------------------------------------------#

#--------------------------------------------------------------------------------#

#--------------------------------------------------------------------------------#
save.image("../objects/image_RNASeq_subset_organoids.RData")
#--------------------------------------------------------------------------------#


extreme_vals_pca <- (prcomp_normalisedcounts$rotation[(abs(prcomp_normalisedcounts$rotation[,1]) > 0.01) &
                                        (abs(prcomp_normalisedcounts$rotation[,2]) < 0.001),])
extreme_vals_pca <- cbind.data.frame(extreme_vals_pca[,1:2], gene=t2g$external_gene_name[match(rownames(extreme_vals_pca), t2g$ensembl_gene_id)])
head(extreme_vals_pca)

give_goterm_list_of_genes_general <- function(i, allgenes){
  require(topGO)
  .x = i
  .xentrez = t2g[match(.x, t2g$ensembl_gene_id),'entrezgene_id']
  .xentrez = .xentrez[!is.na(.xentrez)]
  
  .xres_gosimenrich = GOSim::GOenrichment(genesOfInterest = as.character(.xentrez), allgenes = allgenes)
  .xres_gosim = GOSim::getGOInfo(.xentrez)
  .xres_gosimenrich
}

vector_and_names <- function(i){names(i)=i; return(i)}
remove_na <- function(i) i[!is.na(i)]

## extreme values in PC1
extreme_vals_pca_PC1 <- (prcomp_normalisedcounts$rotation[(abs(prcomp_normalisedcounts$rotation[,1]) > 0.01) &
                                                            (abs(prcomp_normalisedcounts$rotation[,2]) < 0.001),])
extreme_vals_pca_PC1 <- cbind.data.frame(extreme_vals_pca_PC1[,1:2], gene=t2g$external_gene_name[match(rownames(extreme_vals_pca_PC1), t2g$ensembl_gene_id)])
plot(extreme_vals_pca_PC1[,1:2])
extreme_vals_pca_PC1
extreme_vals_pca_PC1_GO = give_goterm_list_of_genes_general(i = vector_and_names(gene_to_ensembl(rownames(extreme_vals_pca_PC1)[extreme_vals_pca_PC1$PC1 > 0])),
                                                            allgenes=(as.character(t2g$entrezgene_id)))

extreme_vals_pca_PC1_GO_2 = give_goterm_list_of_genes_general(i = vector_and_names(gene_to_ensembl(rownames(extreme_vals_pca_PC1)[extreme_vals_pca_PC2$PC1 < 0])),
                                                              allgenes=as.character(t2g$entrezgene_id))

sapply(1:10, function(k) extreme_vals_pca_PC1_GO$GOTerms[extreme_vals_pca_PC1_GO$GOTerms$go_id == names(sort(extreme_vals_pca_PC1_GO$p.values)[k]),'Term'])
sapply(1:10, function(k) extreme_vals_pca_PC1_GO_2$GOTerms[extreme_vals_pca_PC1_GO_2$GOTerms$go_id == names(sort(extreme_vals_pca_PC1_GO_2$p.values)[k]),'Term'])


## extreme values in PC2
extreme_vals_pca_PC2 <- (prcomp_normalisedcounts$rotation[(abs(prcomp_normalisedcounts$rotation[,2]) > 0.01) &
                                                        (abs(prcomp_normalisedcounts$rotation[,1]) < 0.001),])
extreme_vals_pca_PC2 <- cbind.data.frame(extreme_vals_pca_PC2[,1:2], gene=t2g$external_gene_name[match(rownames(extreme_vals_pca_PC2), t2g$ensembl_gene_id)])
plot(extreme_vals_pca_PC2[,1:2])
extreme_vals_pca_PC2
extreme_vals_pca_PC2_GO = give_goterm_list_of_genes_general(i = vector_and_names(rownames(extreme_vals_pca_PC2)[extreme_vals_pca_PC2$PC2 > 0]),
                                                            allgenes=as.character(t2g$entrezgene_id))

extreme_vals_pca_PC2_GO_2 = give_goterm_list_of_genes_general(i = vector_and_names(rownames(extreme_vals_pca_PC2)[extreme_vals_pca_PC2$PC2 < 0]),
                                                            allgenes=as.character(t2g$entrezgene_id))

sapply(1:10, function(k) extreme_vals_pca_PC2_GO$GOTerms[extreme_vals_pca_PC2_GO$GOTerms$go_id == names(sort(extreme_vals_pca_PC2_GO$p.values)[k]),'Term'])
sapply(1:10, function(k) extreme_vals_pca_PC2_GO_2$GOTerms[extreme_vals_pca_PC2_GO_2$GOTerms$go_id == names(sort(extreme_vals_pca_PC2_GO_2$p.values)[k]),'Term'])

## Summary
## along the x axis, we've got on the one hand (separation of PDO5 and PDO6 to all others)
## negative regulation of mitosis, on the right
## apoptosis, receptor-mediated signalling, progesterone
sapply(1:10, function(k) extreme_vals_pca_PC1_GO$GOTerms[extreme_vals_pca_PC1_GO$GOTerms$go_id == names(sort(extreme_vals_pca_PC1_GO$p.values)[k]),'Term'])
sapply(1:10, function(k) extreme_vals_pca_PC1_GO_2$GOTerms[extreme_vals_pca_PC1_GO_2$GOTerms$go_id == names(sort(extreme_vals_pca_PC1_GO_2$p.values)[k]),'Term'])
## along the y axis (most variation of organoids, with PDO2 and PDO10 in an extreme)
## negative regulation of cell cycle pathways and signalling on upwards
## strange GO terms downwards related perhaps to mitosis
sapply(1:10, function(k) extreme_vals_pca_PC2_GO$GOTerms[extreme_vals_pca_PC2_GO$GOTerms$go_id == names(sort(extreme_vals_pca_PC2_GO$p.values)[k]),'Term'])
sapply(1:10, function(k) extreme_vals_pca_PC2_GO_2$GOTerms[extreme_vals_pca_PC2_GO_2$GOTerms$go_id == names(sort(extreme_vals_pca_PC2_GO_2$p.values)[k]),'Term'])


normalised_counts2 = normalised_counts
rownames(normalised_counts2) = make.names(t2g$external_gene_name[match(rownames(normalised_counts2), t2g$ensembl_gene_id)], unique = T)
ggplot(cbind.data.frame(PC=prcomp_normalisedcounts$x[,1:2],
                        # exp=unlist(normalised_counts2['CCNE1',])),
                        # exp=unlist(normalised_counts2['PARP1',])),
                        # exp=unlist(normalised_counts2['MYC',])),
                        # exp=unlist(normalised_counts2['CDK12',])),
                        # exp=unlist(normalised_counts2['BRCA1',])),
                        exp=unlist(normalised_counts2['BRCA2',])),
       aes(x=PC.PC1, y=PC.PC2, col=exp))+
  geom_point()

       
rownames(normalised_counts2)[grepl('MT', rownames(normalised_counts2))]
normalised_counts2_nomt = (normalised_counts2)[!grepl('MT[.]', rownames(normalised_counts2)),]

prcomp_normalisedcounts = prcomp(t(normalised_counts2_nomt), scale. = TRUE, center = TRUE)
plot(prcomp_normalisedcounts$x[,1:2])

#---------------------------------------------------------------------------------#
## ssgsea

subset_df <- function(i, col){
  x <- i[,col]
  names(x) <- rownames(i)
  x
}

# t2g <- readRDS("../../copy_number_analysis_organoids/robjects/t2g2.RDS")
gsva(subset_df(normalised_counts, 'PDO6'), c2BroadSets, method = 'ssGSEA')

rename_rows <- function(i){
  .entrez <- t2g$entrezgene_id[match(rownames(i), t2g$ensembl_gene_id)]
  .dup <- duplicated(.entrez) | is.na(.entrez)
  i <- i[!.dup,]
  .entrez <- .entrez[!.dup]
  rownames(i) <- .entrez
  i
}

ssgsea <- gsva(as(rename_rows(normalised_counts), 'matrix'), c2BroadSets, method = 'ssgsea')

ssgsea_repair <- ssgsea[c('KEGG_HOMOLOGOUS_RECOMBINATION', 'KEGG_MISMATCH_REPAIR',
         'KEGG_BASE_EXCISION_REPAIR', 'KEGG_NUCLEOTIDE_EXCISION_REPAIR',
         'KEGG_NON_HOMOLOGOUS_END_JOINING',
         'KEGG_ERBB_SIGNALING_PATHWAY'),]
saveRDS(ssgsea_repair, "../objects/fig3_ssgsea_repair.RDS")


pdf("../figures/other_DE/pathways_heatmap_orgs.pdf", height = 4, width = 7)
pheatmap(ssgsea_repair)
dev.off()

ggplot(cbind.data.frame(prcomp_normalisedcounts$x[,1:2], PDO=colnames(normalised_counts),
                        GAPDH=log(unlist(normalised_counts[gene_to_ensembl('POLE'),]))),
       aes(x=PC1, y=PC2, label=PDO, col=GAPDH))+geom_point()+geom_label_repel()+
  theme_bw()

ggplot(cbind.data.frame(prcomp_normalisedcounts$x[,1:2], PDO=colnames(normalised_counts),
                        GAPDH=log(unlist(normalised_counts[gene_to_ensembl('POLQ'),]))),
       aes(x=PC1, y=PC2, label=PDO, col=GAPDH))+geom_point()+geom_label_repel()+
  theme_bw()

ggplot(cbind.data.frame(prcomp_normalisedcounts$x[,1:2], PDO=colnames(normalised_counts),
                        GAPDH=log(unlist(normalised_counts[gene_to_ensembl('ERBB2'),]))),
       aes(x=PC1, y=PC2, label=PDO, col=GAPDH))+geom_point()+geom_label_repel()+
  theme_bw()

##-----
## Using TPM
## Remove mitocondrial genes
# t2g_only_autosomal <- t2g[t2g$chromosome_name %in% c(1:22),]
TPM_counts <- sweep(raw_counts0, 2, colSums(raw_counts0), '/')*1e6
raw_counts <- raw_counts0[rownames(raw_counts0) %in% t2g_only_autosomal$ensembl_gene_id,]
TPM_counts <- sweep(raw_counts, 2, colSums(raw_counts), '/')*1e6

prcomp_TPM = prcomp(t(TPM_counts), scale. = FALSE, center = TRUE)
eigs_TPM <- prcomp_TPM$sdev^2
ggplot(cbind.data.frame(prcomp_TPM$x[,1:2], PDO=colnames(TPM_counts)),
       aes(x=PC1, y=PC2, label=PDO))+geom_point()+geom_label_repel()+theme_bw()
ggsave("../figures/PCA_RNASeq/PCA_counts_subset_TPM.png", width = 6.5, height = 5)

ggplot(cbind.data.frame(prcomp_TPM$rotation[t2g$ensembl_gene_id[match(subset_genes, t2g$external_gene_name)],1:2],
                        label=subset_genes))+
  geom_segment(aes(x=0, xend=PC1, y=0, yend=PC2),  arrow = arrow(length = unit(0.02, "npc")))+
  geom_label_repel(aes(x=PC1, y=PC2, label=label))+
  geom_point(data=(cbind.data.frame(label=rownames(prcomp_TPM$x), prcomp_TPM$x[,1:2])),
             # aes(x=PC1*5e-5, y=PC2*5e-5, label=label), col='red', size=3)
             aes(x=PC1*1e-6, y=PC2*1e-6, label=label), col='red', size=3, max.overlaps = 18)+
  geom_label_repel(data=(cbind.data.frame(label=rownames(prcomp_TPM$x), prcomp_TPM$x[,1:2])),
                   # aes(x=PC1*5e-5, y=PC2*5e-5, label=label), col='red', size=3)
                   aes(x=PC1*1e-6, y=PC2*1e-6, label=label), col='red', size=3, max.overlaps = 18)
ggsave("../figures/PCA_RNASeq/PCA_counts_subset_genes_of_interest_TPM.png", width = 8, height = 8)

# ## GAPDH
ggplot(cbind.data.frame(prcomp_TPM$x[,1:2], PDO=colnames(TPM_counts),
                        GAPDH=unlist(TPM_counts['ENSG00000111640',])),
       aes(x=PC1, y=PC2, label=PDO, col=GAPDH))+geom_point()+geom_label_repel()+
  theme_bw()

ggplot(cbind.data.frame(prcomp_TPM$x[,1:2], PDO=colnames(TPM_counts),
                        GAPDH=unlist(TPM_counts[gene_to_ensembl('CDK4'),])),
       aes(x=PC1, y=PC2, label=PDO, col=GAPDH))+geom_point()+geom_label_repel()+
  theme_bw()

ggplot(cbind.data.frame(prcomp_TPM$x[,1:2], PDO=colnames(TPM_counts),
                        GAPDH=log(unlist(TPM_counts[gene_to_ensembl('MYC'),]))),
       aes(x=PC1, y=PC2, label=PDO, col=GAPDH))+geom_point()+geom_label_repel()+
  theme_bw()

ggplot(cbind.data.frame(prcomp_TPM$x[,1:2], PDO=colnames(TPM_counts),
                        GAPDH=unlist(TPM_counts[gene_to_ensembl('MSI2'),])),
       aes(x=PC1, y=PC2, label=PDO, col=GAPDH))+geom_point()+geom_label_repel()+
  theme_bw()

ggplot(cbind.data.frame(prcomp_TPM$x[,1:2], PDO=colnames(TPM_counts),
                        GAPDH=unlist(TPM_counts[gene_to_ensembl('NPM1'),])),
       aes(x=PC1, y=PC2, label=PDO, col=GAPDH))+geom_point()+geom_label_repel()+
  theme_bw()

ggplot(cbind.data.frame(prcomp_TPM$x[,1:2], PDO=colnames(TPM_counts),
                        GAPDH=unlist(TPM_counts[gene_to_ensembl('BRCA1'),])),
       aes(x=PC1, y=PC2, label=PDO, col=GAPDH))+geom_point()+geom_label_repel()+
  theme_bw()
 
# ## remove all the ones with relatively high 3' bias: PDO13, PDO4, PDO9,
# ## (plus the ones already removed: PDO17 PDO14 PDO16 PDO18)
# TPMcounts_subset3pbias_2 <- TPM_counts[,!(colnames(TPM_counts) %in% c('PDO13', 'PDO4', 'PDO9',
#                                           'PDO17', 'PDO14', 'PDO16',
#                                           'PDO18'))]
# prcomp_TPM_subset3pbias_2 = prcomp(t(TPMcounts_subset3pbias_2),
#                                    scale. = FALSE, center = TRUE)
# ggplot(cbind.data.frame(prcomp_TPM_subset3pbias_2$x[,1:2],
#                         PDO=colnames(TPMcounts_subset3pbias_2)),
#        aes(x=PC1, y=PC2, label=PDO))+geom_point()+geom_label_repel()+theme_bw()


hclust_rnaseq <- (hclust(dist(t(TPM_counts))))
saveRDS(hclust_rnaseq, "../objects/hclust_rnaseq_TPM.RDS")

## check genes for proliferation
VEGFA_pathway <- ReactomePA::viewPathway("VEGFA-VEGFR2 Pathway")
VEGFA_pathway2 <- ReactomePA::viewPathway("Signaling by VEGF")
# insulin_pathway <- ReactomePA::viewPathway("Regulation of Insulin-like Growth Factor (IGF) transport and uptake by Insulin-like Growth Factor Binding Proteins (IGFBPs)")
# insulin_pathway <- ReactomePA::viewPathway("Metabolism of proteins Pathway", readable = TRUE)

# ReactomePA::enrichPathway()
library(pathfindR)

ggplot(cbind.data.frame(prcomp_TPM$rotation[,1:2],
                        Gene=t2g_GRCh38$external_gene_name[match(rownames(prcomp_TPM$rotation), t2g_GRCh38$ensembl_gene_id)]),
       aes(x=PC1, y=PC2, label=ifelse(PC1 > 0.1, Gene, NA)))+
  geom_point()+geom_label_repel()

df_TPM_PCA_loadings <- cbind.data.frame(prcomp_TPM$rotation[,1:2],
                 Gene=t2g_GRCh38$external_gene_name[match(rownames(prcomp_TPM$rotation),
                                                          t2g_GRCh38$ensembl_gene_id)])

df_TPM_PCA_loadings$label = ifelse(df_TPM_PCA_loadings$Gene %in% VEGFA_pathway2$data$name,
                                   df_TPM_PCA_loadings$Gene, NA)

df_TPM_PCA_loadings$PACHER_TARGETS_OF_IGF1_AND_IGF2_UP = ifelse(df_TPM_PCA_loadings$Gene %in% t2g_GRCh38$external_gene_name[match(c2BroadSets[['PACHER_TARGETS_OF_IGF1_AND_IGF2_UP']]@geneIds, t2g_GRCh38$entrezgene_id)],
                                                                df_TPM_PCA_loadings$Gene, NA)

ggplot(df_TPM_PCA_loadings,
       aes(x=PC1, y=PC2, label=label))+
  geom_point()+geom_label_repel(max.overlaps = 30)+facet_wrap(.~is.na(label))

ggplot(df_TPM_PCA_loadings,
       aes(x=PC1, y=PC2, label=PACHER_TARGETS_OF_IGF1_AND_IGF2_UP))+
  geom_point()+geom_label_repel(max.overlaps = 30)+facet_wrap(.~is.na(PACHER_TARGETS_OF_IGF1_AND_IGF2_UP))

### select genes that explain only PC1
magnitude_genes <- apply(prcomp_TPM$rotation[,1:2], 1, function(i) sqrt(sum(i**2)))
loadingsPC1_genes_percentage <- abs(prcomp_TPM$rotation[,1])/magnitude_genes

which_genes_PC1 <- (( magnitude_genes >= quantile(magnitude_genes, 0.75) ) & (abs(loadingsPC1_genes_percentage) >= quantile(abs(loadingsPC1_genes_percentage), 0.9, na.rm=T)))

df_TPM_PCA_loadings$important_PC1_genes = which_genes_PC1
ggplot(df_TPM_PCA_loadings,
       aes(x=PC1, y=PC2, label=PACHER_TARGETS_OF_IGF1_AND_IGF2_UP))+
  geom_point()+geom_label_repel(max.overlaps = 30)+facet_wrap(.~important_PC1_genes)

hist(magnitude_genes)
plot(density(log(magnitude_genes)))

PC1_genes_enrichmentpathway <- ReactomePA::enrichPathway(gene = t2g$entrezgene_id[match(df_TPM_PCA_loadings$Gene[df_TPM_PCA_loadings$important_PC1_genes], t2g$external_gene_name)])
PC1_genes_enrichmentpathway@result[1,]
ReactomePA::gsePathway( t2g$entrezgene_id[match(df_TPM_PCA_loadings$Gene[df_TPM_PCA_loadings$important_PC1_genes], t2g$external_gene_name)])

### Analysing the counts
raw_counts0_sorted <- (apply(raw_counts0, 2, sort))
ggplot(melt(raw_counts0_sorted+1), aes(x=Var1, col=Var2, y=value))+geom_step()+scale_y_continuous(trans = 'log')
ggsave("../figures/PCA_RNASeq/raw_counts_persample_sorted.pdf", width = 8, height = 8)

## get loadings in PC1
size_prcomp_normalisedcounts_TPM = apply(prcomp_TPM$rotation[,1:2], 1, function(i) sqrt(sum(i**2)))
# size_prcomp_normalisedcounts_notscaled = apply(prcomp_normalisedcounts_notscaled$rotation[,1:2], 1, function(i) sqrt(sum(i**2)))
top_size_prcomp_normalisedcounts_TPM = order(size_prcomp_normalisedcounts_TPM, decreasing = T)[1:1500]
top_size_prcomp_normalisedcounts2_TPM = order(size_prcomp_normalisedcounts_TPM, decreasing = T)[1:700]

ggplot()+
  geom_segment(data=cbind.data.frame(names=t2g$external_gene_name[match(rownames(prcomp_TPM$rotation[top_size_prcomp_normalisedcounts_TPM,]), t2g$ensembl_gene_id)],
                                     prcomp_TPM$rotation[top_size_prcomp_normalisedcounts_TPM,1:2]),
               aes(x=0, xend=PC1, y=0, yend=PC2, label=names))+
  geom_label_repel(data=cbind.data.frame(names=t2g$external_gene_name[match(rownames(prcomp_TPM$rotation[top_size_prcomp_normalisedcounts_TPM,1:2]), t2g$ensembl_gene_id)],
                                         prcomp_TPM$rotation[top_size_prcomp_normalisedcounts_TPM,1:2]),
                   aes(x=PC1, y=PC2, label=names))+
  geom_label_repel(data=(cbind.data.frame(label=rownames(prcomp_TPM$x), prcomp_TPM$x[,1:2])),
                   aes(x=PC1*5e-7, y=PC2*5e-7, label=label), col='red', size=3)+
  geom_point(data=(cbind.data.frame(label=rownames(prcomp_TPM$x),
                                    # col=unlist(TPM_counts[gene_to_ensembl('IGF2'),]),
                                    col=unlist(TPM_counts[gene_to_ensembl('FTH1'),]),
                                    prcomp_TPM$x[,1:2])),
             aes(x=PC1*1e-4, y=PC2*1e-4, label=label, col=col), size=3, max.overlaps = 18)+
  geom_label_repel(data=(cbind.data.frame(label=rownames(prcomp_TPM$x), prcomp_TPM$x[,1:2])),
                   aes(x=PC1*1e-4, y=PC2*1e-4, label=label), col='red', size=3, max.overlaps = 18)+
  labs(x=paste0('PC1 (', round(100*eigs_TPM[1]/sum(eigs_TPM), 2), '%)' ),
       y=paste0('PC2 (', round(100*eigs_TPM[2]/sum(eigs_TPM), 2), '%)' ))

cosine_sim_NC_top_TPM = outer(1:nrow(prcomp_TPM$rotation[top_size_prcomp_normalisedcounts2_TPM,]),
                          1:nrow(prcomp_TPM$rotation[top_size_prcomp_normalisedcounts2_TPM,]),
                          Vectorize(function(i,j) cosine(prcomp_TPM$rotation[top_size_prcomp_normalisedcounts2_TPM,][i,],
                                                         prcomp_TPM$rotation[top_size_prcomp_normalisedcounts2_TPM,][j,])))
cosine_sim_NC_top_TPM =  acos(cosine_sim_NC_top_TPM)/pi
hclust_on_pc1_and_2_NC_cos_top_TPM = hclust(as.dist(cosine_sim_NC_top_TPM))

cutree_on_pc1_and_2_NC_cos_top_TPM = cutree(hclust_on_pc1_and_2_NC_cos_top_TPM, k = 4)

goterm_per_cat_TPM = mclapply(unique(cutree_on_pc1_and_2_NC_cos_top_TPM), give_goterm_list_of_genes)
names(goterm_per_cat_TPM) = unique(cutree_on_pc1_and_2_NC_cos_top_TPM)

data_goterms_top_TPM = cbind.data.frame(x=t(sapply(unique(cutree_on_pc1_and_2_NC_cos_top_TPM), function(i) colMeans(subset(prcomp_TPM$rotation[,1:2], cutree_on_pc1_and_2_NC_cos_top_TPM == i)))),
                                    go=sapply(goterm_per_cat_TPM, function(j){
                                      j$GOTerms[j$GOTerms$go_id == names(sort(j$p.values)[1]),'Term']
                                    }))
data_goterms_top_TPM$go = sapply(as.character(data_goterms_top_TPM$go), give_linebreaks, max_line = 16)

data_goterms_top_TPM$normPC = t(apply(data_goterms_top_TPM[,c('x.PC1', 'x.PC2')], 1, function(i) i/(sqrt(sum(i**2)))))

data_goterms_top_TPM$x.PC1_mod = data_goterms_top_TPM$x.PC1*2000
data_goterms_top_TPM$x.PC1_mod[6] = 0.3
data_goterms_top_TPM$x.PC2_mod = data_goterms_top_TPM$x.PC2*2000

prcomp_TPM_df <- cbind.data.frame(prcomp_TPM$rotation[top_size_prcomp_normalisedcounts2_TPM,1:2],
                 col=cutree_on_pc1_and_2_NC_cos_top_TPM)
prcomp_TPM_df = cbind(prcomp_TPM_df, t(apply(prcomp_TPM_df[,c('PC1', 'PC2')], 1, function(i) i/(sqrt(sum(i**2))))))
colnames(prcomp_TPM_df)[((ncol(prcomp_TPM_df)-1)):ncol(prcomp_TPM_df)] = c('PC1norm', 'PC2norm')

ggplot()+
  geom_segment(data = prcomp_TPM_df,
               aes(x=0, xend=PC1, y=0, yend=PC2, col=factor(col)))+
  # geom_label_repel(data=(cbind.data.frame(label=rownames(prcomp_TPM$x), prcomp_TPM$x[,1:2])),
  #                  aes(x=PC1*1e-5, y=PC2*1e-5, label=label), col='red', size=3, max.overlaps = 18)+
  # geom_label_repel(data = data_goterms_top_TPM,
  #                  aes(x=x.PC1_mod, y=x.PC2_mod, label=go))+labs(x='PC1', y='PC2')+
  ggtitle('Go terms for clustered loadings')

ggplot()+
  geom_segment(data = prcomp_TPM_df,
               # aes(x=0, xend=PC1, y=0, yend=PC2, col=factor(col)))+
               aes(x=0, xend=PC1norm, y=0, yend=PC2norm, col=factor(col)))+
  xlim(-0.2, 0.7)+
  geom_label_repel(data=(cbind.data.frame(label=rownames(prcomp_TPM$x), prcomp_TPM$x[,1:2])),
                   aes(x=PC1*1e-5, y=PC2*1e-5, label=label), col='red', size=3, max.overlaps = 18)+
  geom_label_repel(data = data_goterms_top_TPM,
                   aes(x=x.PC1_mod, y=x.PC2_mod, label=go))+labs(x='PC1', y='PC2')+
  ggtitle('Go terms for clustered loadings')
ggsave("../figures/PCA_RNASeq/grouped_loadings_GO_subset_TPM.pdf", width = 10, height = 10)


###
## if we only select codingf genes, do we get a better correlation?
pdo7_pdo8 <- cbind.data.frame(PDO7=raw_counts0$PDO7,
                 PDO8=raw_counts0$PDO8,
                 coding=(rownames(raw_counts0) %in% coding$ensembl_gene_id))
pdo7_pdo8_coding_TPM <- sweep(pdo7_pdo8[pdo7_pdo8$coding,1:2], 2, colSums(pdo7_pdo8[pdo7_pdo8$coding,1:2]), '/')*1e6
plot(log(pdo7_pdo8_coding_TPM))

pdo7_pdo8 <- cbind(pdo7_pdo8, sweep(pdo7_pdo8[,1:2], 2, colSums(pdo7_pdo8[,1:2]), '/')*1e6)
colnames(pdo7_pdo8)[(ncol(pdo7_pdo8)-1):ncol(pdo7_pdo8)] = paste0(colnames(pdo7_pdo8)[(ncol(pdo7_pdo8)-1):ncol(pdo7_pdo8)], '_TPM')
colnames(pdo7_pdo8)

ggplot(pdo7_pdo8, aes(x=PDO7_TPM, y=PDO8_TPM, col=coding))+geom_point(alpha=0.1)+
  geom_density_2d(alpha = 0.5)+
  scale_x_continuous(trans = "log2")+scale_y_continuous(trans = "log2")

pdo7_pdo8[!pdo7_pdo8$coding,]


