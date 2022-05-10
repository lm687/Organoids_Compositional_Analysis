
rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(rnaseqRpkg)
library(biomaRt)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(DESeq2)
library(reshape2)
library(ConsensusTME)
library(readxl)

renaming1 = read_excel("../files/PDOnameProperSample_sWGS_RNAseq.xlsx")

## Load file output from DESeq2 for TCGA samples, which has been created with 1_run_DE.sh
load("../objects/deaObjectFile")
deObj = `~response`

deseq_obj_11orgs <- readRDS("../../RNASeq_and_CN/20191218_ViasM_BJ_orgaBrs/output/fig3_renormalised_counts_obj_11orgs.RDS")

##---------------------------------------------------------------------------------------------------------------##
## Run differential expression with DESeq
# In any case, the contrast argument of the function results takes a character vector of length three: the name of the variable, the name of the factor level for the numerator of the log2 ratio, and the name of the factor level for the denominator. The contrast argument can also take other forms, as described in the help page for results and below
# results <- DESeq2::results(deObj, c("response", "complete_remission_or_response", "progressive_disease"),
#                   alpha = 0.05, format = "DataFrame")
# saveRDS(results, file = "../objects/resultsDESeq_TCGA.RDS")
results <- readRDS("../objects/resultsDESeq_TCGA.RDS")
rownames_short = sapply(rownames(results), function(i) strsplit(i, '[.]')[[1]][1])

## Re-name
# mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
# gene_conversion <- getBM(attributes = c("external_gene_name", "ensembl_gene_id", "entrezgene_id"),
#                  filters = "ensembl_gene_id", values = rownames_short,
#                  mart = mart)
# gene_conversion = gene_conversion[match(rownames_short, gene_conversion$ensembl_gene_id),]
# saveRDS(object = gene_conversion, file = "~/Desktop/t2g2.RDS")
t2g = readRDS("../../copy_number_analysis_organoids/robjects/t2g2.RDS")
coding_genes = readRDS("../../copy_number_analysis_organoids/robjects/coding_genes.RDS")
# t2g = readRDS("~/Desktop/t2g.RDS")
gene_conversion = t2g[match(rownames_short, t2g$ensembl_gene_id),]
dim(gene_conversion)
dim(results)

results = cbind(gene_conversion, results)
results = results[order(results$log2FoldChange),]

highDA = results[1:20,]
lowDA = highDA[nrow(highDA):(nrow(highDA)-20),]

ggplot(cbind.data.frame(log2fc= results$log2FoldChange, pval=-log10(results$padj)))+
  geom_point(aes(x=log2fc, y=pval))

vlcano = function (res) {
  df <- data.frame(labls = res$external_gene_name, lfc = res$log2FoldChange, pv = -log10(res$padj), 
                   isDE = ifelse(is.na(res$padj), FALSE, res$padj < 0.05))
  colNonSig <- "gray32"
  colSig <- "red3"
  colOut <- "blue"
  dotShape <- 20
  upShape <- 24
  outShape <- 19
  xlim <- c(-1, 1) * quantile(abs(df$lfc[is.finite(df$lfc)]), 
                              probs = 0.995) * 1.1
  ylim <- c(0, quantile(df$pv[is.finite(df$pv)], probs = 0.995) * 
              1.1)
  df$vy <- pmin(df$pv, ylim[2])
  df$vx <- pmax(xlim[1], pmin(xlim[2], df$lfc))
  df$outY <- df$pv > ylim[2]
  df$outX <- df$lfc < xlim[1] | df$lfc > xlim[2]
  df$outXY <- df$outX & df$outY
  dfNonSig <- df[!df$isDE, ]
  dfNonSig$vy[is.na(dfNonSig$vy)] = 0
  dfSigInlier <- df[df$isDE & !df$outX & !df$outY, ]
  dfSigOutlier <- df[df$isDE & (df$outX | df$outY), ]
  dfSigOutlier$shape <- ifelse(dfSigOutlier$outXY, outShape, 
                               ifelse(dfSigOutlier$outY, upShape, dotShape))
    p <- ggplot() +
      geom_point(data = dfNonSig, aes(x = .data$vx,   y = .data$vy),
                 colour = colNonSig, pch = dotShape) + 
      geom_point(data = dfSigInlier, aes(x = .data$vx, y = .data$vy), 
               colour = colSig, shape = dotShape) +
      ggrepel::geom_label_repel(data=dfSigInlier, aes(x=  .data$vx, y = .data$vy, label=labls))+
      geom_point(data = dfSigOutlier,  aes(x = .data$vx, y = .data$vy),
               colour = colOut, shape = dfSigOutlier$shape) + 
      ggrepel::geom_label_repel(data=dfSigOutlier, aes(x=  .data$vx, y = .data$vy, label=labls))+
      theme(legend.position = "none") + ggtitle("Volcano Plot") + 
      labs(x = "Log of Fold Change", y = "-log10 of FDR")+
      ylim(c(0, 1.8))
  return(p)
}
vlcano(results)
ggsave("../figures/Sensitive_resistant_figures/volcano_plot.pdf", width = 20, height = 20)


## save results in a different format
for(flename in c('differential_all', 'differential_sig')){
  dif_all = read.table(paste0("../objects/", flename), sep = ",", header = TRUE, stringsAsFactors = FALSE)
  # if(flename == "differential_all"){
    dif_all = cbind(gene_conversion$external_gene_name, dif_all[match(gene_conversion$ensembl_gene_id, sapply(dif_all$X, function(i) strsplit(i, '[.]')[[1]][1])  ),])
    dif_all = dif_all[!is.na(dif_all$X),]
  # }else if(flename == "differential_sig"){
  #   dif_all = cbind(gene_conversion$external_gene_name, dif_all[match(gene_conversion$ensembl_gene_id, sapply(dif_all$X, function(i) strsplit(i, '[.]')[[1]][1])
  #   ),])
  # }
  write.table(x = dif_all, file = paste0("../files/", flename, ".csv"), sep = ",", quote = FALSE)
  
}
##---------------------------------------------------------------------------------------------------------------##


##---------------------------------------------------------------------------------------------------------------##
## Here I am comparing the TCGA DE and the organoid DE for sensitive vs resistant

DE_results_TCGA = read.table("../objects/differential_all", sep = ',', header = T)

## Load the organoids DE
load("../files/deObject_SampleGroup_sensitive_vs_resistant.RData")
SampleGroup_14orgs <- SampleGroup
rm(SampleGroup)
## e.g. JBLAB-19920 is PDO18 and has been excluded based on 3' bias
# "JBLAB-19920" %in% colnames(SampleGroup)
# length(colnames(SampleGroup)) ## 3' bias have been excluded. No!!! one has been left
# ## but I am not sure why there are 14 instead of 15
# "JBLAB-19939" %in% colnames(SampleGroup) ## the fourth sample with highest 3' bias

## re-run DESeq2 with only 11 organoids
## use the previous object to get the group of sensitive/resistant
coldata_SampleGroup_14orgs <- colData(SampleGroup_14orgs)
coldata_SampleGroup_14orgs$PDO = renaming1$PDO[match(coldata_SampleGroup_14orgs$SampleName,
                                                     renaming1$sampleNameRNAseq)]
dim(DESeq2::counts(deseq_obj_11orgs))
deseq_obj_11orgs_normalFT <- deseq_obj_11orgs
## remove FP as they are not sensitive/resistant, and remove PDO15 as well, because we don't
## have drug results for it
deseq_obj_11orgs <- deseq_obj_11orgs[,grepl('PDO', colnames(deseq_obj_11orgs))]
deseq_obj_11orgs <- deseq_obj_11orgs[,colnames(deseq_obj_11orgs) != 'PDO15']
deseq_obj_11orgs$responseGroup = coldata_SampleGroup_14orgs$responseGroup[match(rownames(colData(deseq_obj_11orgs)), coldata_SampleGroup_14orgs$PDO)]
deseq_obj_11orgs_sensitivity <- DESeq2::DESeqDataSetFromMatrix(DESeq2::counts(deseq_obj_11orgs),
                                                               colData=colData(deseq_obj_11orgs),
                                                               design = ~responseGroup)
rm(deseq_obj_11orgs)
deseq_obj_11orgs_sensitivity <- estimateSizeFactors(deseq_obj_11orgs_sensitivity)
deseq_obj_11orgs_sensitivity <- estimateDispersions(deseq_obj_11orgs_sensitivity,fitType="local")
deseq_obj_11orgs_sensitivity <- DESeq2::DESeq(deseq_obj_11orgs_sensitivity)
deseq_obj_11orgs_results_sensitivity <- DESeq2::results(deseq_obj_11orgs_sensitivity)
saveRDS(deseq_obj_11orgs_results_sensitivity, "../objects/deseq_obj_11orgs_results_sensitivity.RDS")
deseq_obj_11orgs <- deseq_obj_11orgs_sensitivity ## below we are studying sensitivity vs specificity

DE_results_org <- deseq_obj_11orgs
# remove_na = function(i)i[!is.na(i)]

# SampleGroup <- SampleGroup[,-remove_na(match(c('JBLAB-19920', 'JBLAB-19925', 'JBLAB-19936'), SampleGroup$SampleName))]
# SampleGroup = DESeq2::DESeq(SampleGroup)
# DE_results_org = DESeq2::results(SampleGroup)
# DE_results_org = DESeq2::results(deseq_obj_11orgs)
# saveRDS(DE_results_org, file = "../objects/resultsDESeq_org_11orgs.RDS")
# DE_results_org <- readRDS("../objects/resultsDESeq_org_11orgs.RDS") ## this is fallopian tube vs normal
# saveRDS(SampleGroup, file = "../objects/SampleGroup_org_no3primebias.RDS")
# SampleGroup <- readRDS("../objects/SampleGroup_org_no3primebias.RDS")
SampleGroup <- deseq_obj_11orgs
DE_results_org

plot(DE_results_TCGA$log2FoldChange, -log(DE_results_TCGA$padj), cex=.1, pch=19)

ggplot(droplevels(DE_results_TCGA), aes(x=log2FoldChange, y=-log(padj)))+geom_point(size=0.1)+
  lims(y=c(0, 1.5), x=c(-2.5,2.5))

all_DE = cbind.data.frame(TCGA=DE_results_TCGA[match(rownames(DE_results_org), gsub("[.].*", "", DE_results_TCGA$X)),],
                          org=DE_results_org, gene=rownames(DE_results_org))
all_DE$GeneName=t2g$external_gene_name[match(gsub("[.].*", "", all_DE$TCGA.X), t2g$ensembl_gene_id)]
head(all_DE)

plot(all_DE$TCGA.log2FoldChange, all_DE$org.log2FoldChange)
plot(all_DE$TCGA.padj, all_DE$org.padj)

tier1 = read.table("../files/Census_allFri Feb 12 09 56 06 2021.tsv", sep = "\t", header = TRUE)
tier1$ensembl = sapply(tier1$Synonyms, function(i) gsub("[.].*", "", strsplit(i, ',')[[1]][2]))

ggplot(all_DE, aes(x=TCGA.log2FoldChange, y=org.log2FoldChange))+geom_point(alpha=0.2)
ggsave("../figures/Sensitive_resistant_figures/tcga_org_comparison_log2FC.pdf")

ggplot(all_DE %>% filter(abs(TCGA.log2FoldChange) > quantile(abs(all_DE$TCGA.log2FoldChange), probs = c(.9), na.rm = T)),
       aes(x=TCGA.log2FoldChange, y=org.log2FoldChange, label=GeneName))+geom_point(alpha=0.2)+geom_label_repel()
ggsave("../figures/Sensitive_resistant_figures/tcga_org_comparison_log2FC_onlyhighDE.pdf")

ggplot(all_DE %>% filter(gene %in% tier1$ensembl), aes(x=TCGA.log2FoldChange, y=org.log2FoldChange))+geom_point(alpha=0.2)
ggsave("../figures/Sensitive_resistant_figures/tcga_org_comparison_log2FC_onlytier1.pdf")

ggplot(all_DE %>% filter( (TCGA.padj  < 0.001) | (org.padj  < 0.001) ),
       aes(x=TCGA.log2FoldChange, y=org.log2FoldChange, label=paste0(signif(TCGA.padj, 3), ' - ', signif(org.padj, 3))))+
  geom_point(alpha=0.2)+geom_label_repel()
ggsave("../figures/Sensitive_resistant_figures/tcga_org_comparison_log2FC_onlylowpval.pdf")

ggplot(all_DE, aes(x=TCGA.padj, y=org.padj))+geom_point(alpha=0.2)
ggsave("../figures/Sensitive_resistant_figures/tcga_org_comparison_pvals.pdf")
##---------------------------------------------------------------------------------------------------------------##

##---------------------------------------------------------------------------------------------------------------##
## not look at differential abundance, but simply at gene expression and whether it correlates
## counts for organoids
# counts_DESeq_org = read.table("../files/counts_norm.csv", sep = ",", header=T, stringsAsFactors = FALSE)
# rownames(counts_DESeq_org) = counts_DESeq_org[,1]; counts_DESeq_org = counts_DESeq_org[,-1]
## Remove the organoids that have a 3' bias, and the normal samples
## PDO14 (JBLAB-19925) PDO16 (JBLAB-19936) PDO18 (JBLAB-19920) and normal (JBLAB-19950, JBLAB-19953, JBLAB-19952
## JBLAB-19954, JBLAB-19955, JBLAB-19951)
# counts_DESeq_org = counts_DESeq_org[,!(colnames(counts_DESeq_org) %in% c('JBLAB.19925', 'JBLAB.19936', 'JBLAB.19920', 'JBLAB.19950',
#                                                                          'JBLAB.19953', 'JBLAB.19952',
#                                                                          'JBLAB.19954', 'JBLAB.19955', 'JBLAB.19951'))]

counts_DESeq_org <- DESeq2::counts(deseq_obj_11orgs, normalized=T)
## select only PDO
counts_DESeq_org <- counts_DESeq_org[,grepl('PDO', colnames(counts_DESeq_org))]
apply(counts_DESeq_org, 2, function(i) as.numeric(i))
colnames(counts_DESeq_org) ## we have more than 18 because we also have normal samples

list_to_df = function(a){
  .x = apply(a, 2, function(i) as.numeric(i))
  rownames(.x) = rownames(a); colnames(.x) = colnames(a)
  return(.x)
}
counts_DESeq_org = list_to_df(counts_DESeq_org)

## counts for TCGA
counts_DESeq_TCGA = DESeq2::counts(deObj, normalized=TRUE)
counts_DESeq_TCGA_raw = DESeq2::counts(deObj, normalized=FALSE)
rownames(counts_DESeq_TCGA) = gsub("[.].*", "", rownames(counts_DESeq_TCGA))
rownames(counts_DESeq_TCGA) == rownames(counts_DESeq_org)

## why the difference?
length(rownames(counts_DESeq_org))
length(rownames(counts_DESeq_TCGA))

counts_DESeq_TCGA = counts_DESeq_TCGA[match(rownames(counts_DESeq_org), rownames(counts_DESeq_TCGA)),]
length(rownames(counts_DESeq_org))
length(rownames(counts_DESeq_TCGA))

rownames(counts_DESeq_TCGA) = rownames(counts_DESeq_org) = make.names(gene_conversion$external_gene_name[match(rownames(counts_DESeq_org), gene_conversion$ensembl_gene_id)], unique = T)

subset_genes = c('MYC', 'CCNE1', 'PIK3CA', 'TERT', 'KRAS', 'PTEN', 'RB1', 'AKT1',
                 'AKT2', 'PARP1', 'PARP2', 'ATM', 'ATR', 'WEE1', 'TOP1', 'TUBB1',
                 'AKT3', 'CCND1', 'CCND2', 'CCND3', 'CDKN2A', 'CDKN2B', 'MECOM', 'CDK12')
df_all_genes <- melt(list(TCGA=counts_DESeq_TCGA,
                             org=counts_DESeq_org))
df_subset_genes <- melt(list(TCGA=counts_DESeq_TCGA[rownames(counts_DESeq_TCGA) %in% subset_genes,],
                             org=counts_DESeq_org[rownames(counts_DESeq_org) %in% subset_genes,]))

# ggplot(df_all_genes, aes(y=value, x=L1))+geom_boxplot()+geom_jitter()

means_tcga <- rowMeans(counts_DESeq_TCGA)
means_org <- rowMeans(counts_DESeq_org)

pdf("../figures/Sensitive_resistant_figures/colmeans_deseqcounts_correlation_tcga_org.pdf")
plot(log(means_tcga), log(means_org), xlab='DESeq normalised counts from TCGA',ylab='DESeq normalised counts from organoids')
abline(coef=c(0,1), lty='dashed', col='blue')
dev.off()

# require(GSVA)
# require(GSVAdata)
# data("c2BroadSets")
# c2BroadSets$NAKAMURA_CANCER_MICROENVIRONMENT_UP
# table(sapply(c2BroadSets, length))

ConsensusTMB_OV = read.table("https://raw.githubusercontent.com/cansysbio/ConsensusTME/master/Consensus_Signatures/OV_Consensus_Signatures.txt", header = T)
genes_TMB_OV = unique(unlist(ConsensusTMB_OV))
colours_TME_tcga_org_cor = match(names(means_tcga), genes_TMB_OV)
colours_TME_tcga_org_cor[!is.na(colours_TME_tcga_org_cor)] = 'TME'
colours_TME_tcga_org_cor[is.na(colours_TME_tcga_org_cor)] = 'Other'
df_colmeans_deseqcounts_correlation_tcga_org <- cbind.data.frame(means_tcga, means_org, TME=colours_TME_tcga_org_cor)
df_colmeans_deseqcounts_correlation_tcga_org = df_colmeans_deseqcounts_correlation_tcga_org[colSums(apply(df_colmeans_deseqcounts_correlation_tcga_org, 1, is.na)) == 0,]
saveRDS(df_colmeans_deseqcounts_correlation_tcga_org, "../objects/fig4_df_colmeans_deseqcounts_correlation_tcga_org.RDS")
df_colmeans_deseqcounts_correlation_tcga_org$xidentity=rep(c(.01, max(means_tcga, na.rm=T)), nrow(df_colmeans_deseqcounts_correlation_tcga_org)/2)
df_colmeans_deseqcounts_correlation_tcga_org$yidentity=rep(c(.01, max(means_org, na.rm=T)), nrow(df_colmeans_deseqcounts_correlation_tcga_org)/2)
ggplot(df_colmeans_deseqcounts_correlation_tcga_org,
       aes(x=(means_tcga), y=(means_org), col=TME))+geom_point()+
  # geom_abline(slope = (1), intercept = 0, lty='dashed')
  facet_wrap(.~TME)+
  scale_x_continuous(trans = "log2")+scale_y_continuous(trans = "log2")+
  geom_line(aes(x=xidentity, y=yidentity), lty='dashed', col='black')+
  theme(legend.position = "bottom")+ggtitle('Comparison of DESeq counts between TCGA\nand organoid samples')+
  theme_bw()+labs(x='DESEq count means for TCGA', y='DESEq count means for organoid samples')
ggsave("../figures/Sensitive_resistant_figures/colmeans_deseqcounts_correlation_tcga_org_TME.pdf", width = 6, height = 4)
ggsave("../figures/Sensitive_resistant_figures/colmeans_deseqcounts_correlation_tcga_org_TME.png", width = 6, height = 4)

keep_pos = (df_colmeans_deseqcounts_correlation_tcga_org$means_tcga > 0) & (df_colmeans_deseqcounts_correlation_tcga_org$means_org > 0)
lm_colmeans_deseqcounts_correlation_tcga_org = lm(log(means_org) ~ log(means_tcga),
  data = df_colmeans_deseqcounts_correlation_tcga_org[keep_pos,])
quantile_lm_cor = quantile(abs(lm_colmeans_deseqcounts_correlation_tcga_org$effects), 0.9)
lower_diag = lm_colmeans_deseqcounts_correlation_tcga_org$residuals < -quantile_lm_cor
upper_diag = lm_colmeans_deseqcounts_correlation_tcga_org$residuals > quantile_lm_cor

df_colmeans_deseqcounts_correlation_tcga_org[keep_pos,'cat'] = paste0(lower_diag, upper_diag)
ggplot(df_colmeans_deseqcounts_correlation_tcga_org,
       aes(x=means_tcga, y=means_org, col=cat))+geom_point()+
  scale_x_continuous(trans = "log2")+scale_y_continuous(trans = "log2")+
  geom_abline(slope = lm_colmeans_deseqcounts_correlation_tcga_org$coefficients[2],
              intercept = lm_colmeans_deseqcounts_correlation_tcga_org$coefficients[1], lty='dashed')+
  theme(legend.position = "bottom")+ggtitle('Comparison of DESeq counts between TCGA\nand organoid samples')+
  theme_bw()+labs(x='DESEq count means for TCGA', y='DESEq count means for organoid samples')

## do GE on this
names(lower_diag[lower_diag])
names(upper_diag[upper_diag])

counts_DESeq_TCGA_upperdiag = counts_DESeq_TCGA[match(names(upper_diag[upper_diag]), rownames(counts_DESeq_TCGA)),]
clust_counts_DESeq_TCGA_upperdiag = hclust(dist(log(counts_DESeq_TCGA_upperdiag + 0.001)))
clust_counts_DESeq_TCGA_upperdiag_cutree = cutree(tree = clust_counts_DESeq_TCGA_upperdiag, k = 10)
pheatmap::pheatmap(log(counts_DESeq_TCGA_upperdiag+0.001), show_rownames = F, show_colnames = F)
counts_DESeq_TCGA_lowerdiag = counts_DESeq_TCGA[match(names(lower_diag[lower_diag]), rownames(counts_DESeq_TCGA)),]
clust_counts_DESeq_TCGA_lowerdiag = hclust(dist(log(counts_DESeq_TCGA_lowerdiag + 0.001)))
clust_counts_DESeq_TCGA_lowerdiag_cutree = cutree(tree = clust_counts_DESeq_TCGA_lowerdiag, k = 10)
pheatmap::pheatmap(log(counts_DESeq_TCGA_lowerdiag+0.001), show_rownames = F, show_colnames = F)


cutreecats = c(paste0('lowerdiag', clust_counts_DESeq_TCGA_lowerdiag_cutree),
  paste0('upperdiag', clust_counts_DESeq_TCGA_upperdiag_cutree))
df_colmeans_deseqcounts_correlation_tcga_org$cutreecats<- cutreecats[match(rownames(df_colmeans_deseqcounts_correlation_tcga_org),
    c(names(clust_counts_DESeq_TCGA_lowerdiag_cutree), names(clust_counts_DESeq_TCGA_upperdiag_cutree)))]

ggplot(df_colmeans_deseqcounts_correlation_tcga_org,
       aes(x=means_tcga, y=means_org, col=cutreecats))+geom_point()+
  scale_x_continuous(trans = "log2")+scale_y_continuous(trans = "log2")+
  geom_abline(slope = lm_colmeans_deseqcounts_correlation_tcga_org$coefficients[2],
              intercept = lm_colmeans_deseqcounts_correlation_tcga_org$coefficients[1], lty='dashed')+
  theme(legend.position = "bottom")+ggtitle('Comparison of DESeq counts between TCGA\nand organoid samples')+
  theme_bw()+labs(x='DESEq count means for TCGA', y='DESEq count means for organoid samples')+
  facet_wrap(.~cutreecats)

lower_diag[lower_diag]

t2g_matched_entrez = t2g$entrezgene_id[match(rownames(counts_DESeq_org), t2g$external_gene_name)]
library(GOSim)
require(topGO)
list_categories_cor = list()

cutreecats_list = split(x = rownames(df_colmeans_deseqcounts_correlation_tcga_org),
      f = factor(df_colmeans_deseqcounts_correlation_tcga_org$cutreecats))
names(cutreecats_list) = levels(factor(df_colmeans_deseqcounts_correlation_tcga_org$cutreecats))
# goterm_per_cutreecat = mclapply(cutreecats_list, function(i){
#   .x = i
#   .xentrez = t2g[match(.x, t2g$external_gene_name),'entrezgene_id']
#   .xentrez = .xentrez[!is.na(.xentrez)]
#   .xres_gosimenrich = GOSim::GOenrichment(genesOfInterest = as.character(.xentrez),
#                                         allgenes = as.character(t2g_matched_entrez[!is.na(t2g_matched_entrez)]))
#   .xres_gosim = GOSim::getGOInfo(.xentrez)
#   return(.xres_gosim)
# })
# saveRDS(goterm_per_cutreecat, "../objects/goterm_per_cutreecat.RDS")
goterm_per_cutreecat <- readRDS("../objects/goterm_per_cutreecat.RDS")

cbind.data.frame(sapply(cutreecats_list, length), sapply(goterm_per_cutreecat, typeof))
lapply(goterm_per_cutreecat, function(j) try(sapply(1:10, function(i) j[,i]$Term[1])))

short_go = c('Signalling, immune response', NA, 'Hypoxia +', NA, 'Cell proliferation, immune', NA, 'Metabolism +', NA, NA, NA, 'Metabolism +', NA, 'Chemotaxis +', NA, NA, NA, 'Apoptosis, DNA', NA, NA, NA)
df_colmeans_deseqcounts_correlation_tcga_org = cbind.data.frame(df_colmeans_deseqcounts_correlation_tcga_org,
                                                                short_go=short_go[match(df_colmeans_deseqcounts_correlation_tcga_org$cutreecats, names(cutreecats_list))])
ggplot(droplevels(df_colmeans_deseqcounts_correlation_tcga_org[!is.na(df_colmeans_deseqcounts_correlation_tcga_org$short_go),]),
       aes(x=means_tcga, y=means_org, col=cutreecats))+geom_point(alpha=0.4)+
  scale_y_continuous(trans = "log2")+
  geom_abline(slope = lm_colmeans_deseqcounts_correlation_tcga_org$coefficients[2],
              intercept = lm_colmeans_deseqcounts_correlation_tcga_org$coefficients[1], lty='dashed')+
  theme(legend.position = "bottom")+ggtitle('Comparison of DESeq counts between TCGA and organoid samples')+
  theme_bw()+labs(x='DESEq count means for TCGA', y='DESEq count means for organoid samples')+
  facet_wrap(.~short_go)+guides(col=FALSE)+scale_x_continuous(trans = "log2")
ggsave("../figures/Sensitive_resistant_figures/colmeans_deseqcounts_correlation_tcga_org_GOterms.png", width = 6.5, height = 5)


ggplot()+
  geom_violin(data=df_subset_genes %>% filter(L1 == 'TCGA'), aes(x=factor(Var1, levels=names(sort(means_tcga[subset_genes]))),
                                                                 y=value, col=L1))+
  geom_point(data=df_subset_genes %>% filter(L1 == 'org'), aes(x=factor(Var1, levels=names(sort(means_tcga[subset_genes]))),
                                                               y=value, col=L1))+
  scale_y_continuous(trans='log10')+theme_bw()+labs(x='Genes', y='Normalised DESeq2 counts')
ggsave("../figures/Sensitive_resistant_figures/deseqcounts_correlation_tcga_org_selected_genes.pdf", width = 14)

## I should do the same for TPMs
counts_DESeq_TCGA_raw
counts_DESeq_org_raw = read.table("../files/20191218_ViasM_BJ_orgaBrs_tpm.csv", sep = ",", header=T, stringsAsFactors = FALSE)
rownames(counts_DESeq_org_raw) = counts_DESeq_org_raw$gene_id; counts_DESeq_org_raw = counts_DESeq_org_raw[,-c(1:9)]
counts_DESeq_org_raw = list_to_df(counts_DESeq_org_raw)
## Remove the organoids that have a 3' bias, and the normal samples
## PDO14 (JBLAB-19925) PDO16 (JBLAB-19936) PDO18 (JBLAB-19920) and normal (JBLAB-19950, JBLAB-19953, JBLAB-19952
## JBLAB-19954, JBLAB-19955, JBLAB-19951)
counts_DESeq_org_raw = counts_DESeq_org_raw[,!(colnames(counts_DESeq_org_raw) %in% c('JBLAB19925', 'JBLAB19936', 'JBLAB19920', 'JBLAB19950',
                                                                         'JBLAB19953', 'JBLAB19952',
                                                                         'JBLAB19954', 'JBLAB19955', 'JBLAB19951'))]

counts_DESeq_TCGA_TPM = sweep(counts_DESeq_TCGA_raw, 2, colSums(counts_DESeq_TCGA_raw), '/')
counts_DESeq_org_TPM = sweep(counts_DESeq_org_raw, 2, colSums(counts_DESeq_org_raw), '/')

rownames(counts_DESeq_TCGA_TPM) = gsub("[.].*", "", rownames(counts_DESeq_TCGA_TPM))
counts_DESeq_TCGA_TPM = counts_DESeq_TCGA_TPM[match(rownames(counts_DESeq_org_TPM), rownames(counts_DESeq_TCGA_TPM)),]
dim(counts_DESeq_TCGA_TPM)[1] == dim(counts_DESeq_org_TPM)[1]
rownames(counts_DESeq_TCGA_TPM) = rownames(counts_DESeq_org_TPM) = make.names(gene_conversion$external_gene_name[match(rownames(counts_DESeq_org_TPM), gene_conversion$ensembl_gene_id)], unique = T)
df_all_genes_TPM <- melt(list(TCGA=counts_DESeq_TCGA_TPM,
                          org=counts_DESeq_org_TPM))
df_subset_genes_TPM <- df_all_genes_TPM %>% filter(Var1 %in% subset_genes)
means_tcga_TPM <- rowMeans(counts_DESeq_TCGA_TPM)
means_org_TPM <- rowMeans(counts_DESeq_org_TPM)

pvals_KS_bool = sapply(names(sort(means_tcga_TPM[subset_genes])), function(gene_it){
  ks.test(x = df_subset_genes_TPM %>% filter(L1 == 'TCGA', Var1 == gene_it) %>% select(value) %>% unlist,
          y = df_subset_genes_TPM %>% filter(L1 == 'org', Var1 == gene_it) %>% select(value) %>% unlist)$p.value
})*length(subset_genes) < 0.05

ggplot()+
  geom_violin(data=df_subset_genes_TPM %>% filter(L1 == 'TCGA'), aes(x=factor(Var1, levels=names(sort(means_tcga_TPM[subset_genes]))),
                                                                 y=value))+
  geom_point(data=df_subset_genes_TPM %>% filter(L1 == 'org'), aes(x=factor(Var1, levels=names(sort(means_tcga_TPM[subset_genes]))),
                                                               y=value, col=L1))+
  geom_text(data=cbind.data.frame(gene=names(pvals_KS_bool)[pvals_KS_bool],
                                   y=means_tcga_TPM[names(pvals_KS_bool)[pvals_KS_bool]]), aes(x=gene, y=1e-3, label='*'), size=12)+
  scale_y_continuous(trans='log10')+theme_bw()+labs(x='Genes', y='TPM')
ggsave("../figures/Sensitive_resistant_figures/TPM_correlation_tcga_org_selected_genes.pdf", width = 14)


pdf("../figures/Sensitive_resistant_figures/colmeans_TPM_correlation_tcga_org.pdf")
plot(log(means_tcga_TPM), log(means_org_TPM), xlab='TPM from TCGA',ylab='TPM from organoids')
abline(coef=c(0,1), lty='dashed', col='blue')
dev.off()

cor(x = means_tcga_TPM[!(is.na(means_tcga_TPM) | is.na(means_org_TPM))], y = means_org_TPM[!(is.na(means_tcga_TPM) | is.na(means_org_TPM))])

length(table(droplevels(df_subset_genes_TPM[df_subset_genes_TPM$L1 == 'org','Var2'])))

##---------------------------------------------------------------------------------------------------------------##
## Check that we have not included normal TCGA samples
colnames(counts_DESeq_TCGA_raw)
table(deObj$response)
deObj$response

## also to confirm:
query_GDC <- readRDS("../objects/query")
dim(query_GDC$results[[1]])
table(query_GDC$results[[1]]$sample_type)

##---------------------------------------------------------------------------------------------------------------##
## PCA of TCGA
pca_tcga <- prcomp(t(counts_DESeq_TCGA[!apply(counts_DESeq_TCGA, 1, function(i) any(is.na(i)) ),]))
genes_pca <- c('PARP2', 'TERT', 'ATR', 'AKT1', 'KRAS', 'TOP1',
               'PIK3CA', 'ATM', 'PARP1', 'CCNE1', 'PTEN')
plot(pca_tcga$x)
pca_tcga2 <- pca_tcga$rotation[match(genes_pca, rownames(pca_tcga$rotation)),][,1:2]
ggplot(cbind.data.frame(pca_tcga2, gene=rownames(pca_tcga2)),
       aes(label=gene))+geom_segment(aes(x=0, y=0, xend=PC1, yend=PC2),
                                     arrow = arrow(length = unit(0.02, "npc")))+
  geom_label_repel(aes(x=PC1, y=PC2))

##---------------------------------------------------------------------------------------------------------------##
## Clustering of TCGA
# hclust_genes <- hclust(dist(counts_DESeq_TCGA))
# clashes
library(umap)
counts_DESeq_org_nonorm <- counts_DESeq_org[! (colnames(counts_DESeq_org) %in% c('JBLAB.19950',
                                                                                 'JBLAB.19953',
                                                                                 'JBLAB.19952',
                                                                                 'JBLAB.19954',
                                                                                 'JBLAB.19955',
                                                                                 'JBLAB.19951'))]
joint_counts <- cbind(counts_DESeq_TCGA[match(rownames(counts_DESeq_org_nonorm), rownames(counts_DESeq_TCGA)),],
                      counts_DESeq_org_nonorm)
joint_counts <- joint_counts[!apply(joint_counts, 1, function(i) any(is.na(i))),]
umap_2_with_orgs <- umap(t(joint_counts))
prcomp_with_orgs <- prcomp(t(joint_counts[,]), center = T, scale. = T)

umap_TCGA <- umap(counts_DESeq_TCGA)
set.seed(1325)
umap_2 <- umap(t(counts_DESeq_TCGA))

annotation_umap <- cbind.data.frame(Gene=rownames(counts_DESeq_TCGA))
annotation_umap$replication <- sapply(coding_genes$description, grepl, pattern = 'replication' )[match(rownames(counts_DESeq_TCGA), coding_genes$ensembl_gene_id)]
annotation_umap$consensusTMB <- (coding_genes$external_gene_name[match(rownames(counts_DESeq_TCGA), coding_genes$ensembl_gene_id)] %in% unique(unlist(ConsensusTMB_OV)))[match(rownames(counts_DESeq_TCGA), coding_genes$ensembl_gene_id)]
plot_umap(umap_TCGA, factor(annotation_umap$replication, levels=c('FALSE', 'TRUE')))
plot_umap(umap_TCGA, factor(annotation_umap$consensusTMB))

## from file 0_preparing_files.R
clinical_tcga <- read.table("../files/clinical.cases_selection.2021-05-12/clinical.tsv", sep = "\t", fill = T, header = T)
matched_clinical <- clinical_tcga[match(query_GDC$results[[1]]$cases.submitter_id[match(colnames(counts_DESeq_TCGA),
                                                                                        query_GDC$results[[1]]$id)], clinical_tcga$case_submitter_id),]
wgd=read.table("../../../../CDA_in_Cancer/data/TCGA_WGD/Haase2019_TCGA.giScores.wgd.txt", header = T)
matched_wgs <- wgd$wgd[match(query_GDC$results[[1]]$cases.submitter_id[match(colnames(counts_DESeq_TCGA), query_GDC$results[[1]]$id)],
                                  wgd$patient)]
exposures_TCGA <- readRDS("../../copy_number_analysis_organoids/data/Export-matrix_OV_Sigs_on_TCGA-OV_12112019.rds")
#consensusTME
change_rownames <- function(i){
  rownames(i) = t2g$external_gene_name[match(gsub("\\..*","",rownames(i)), t2g$ensembl_gene_id)]
  i
}
matched_consensusTME <- ConsensusTME::consensusTMEAnalysis(change_rownames(counts_DESeq_TCGA_raw), cancerType = "OV")

## clinical ov data
TCGA_genes <- readRDS("../../../../other_repos/cnsigs_Initial_submission/survival_analysis/from_ruben/survival_models/TCGA_OVBRCAonly_Exposures_and_BRCA_Status_plusGene.rds")
matched_exposures <- exposures_TCGA[match(substring(query_GDC$results[[1]]$cases[match(colnames(counts_DESeq_TCGA), query_GDC$results[[1]]$id)], 1, 12),
                                          rownames(exposures_TCGA)),]

matched_group_WGD <- cutree(cluster_fig1, k=2)[match(substring(query_GDC$results[[1]]$cases[match(colnames(counts_DESeq_TCGA), query_GDC$results[[1]]$id)], 1, 12),
                                          names(cutree(cluster_fig1, k=2)))]
matched_group_genes <- TCGA_genes[match(substring(query_GDC$results[[1]]$cases[match(colnames(counts_DESeq_TCGA), query_GDC$results[[1]]$id)], 1, 12),
                                        TCGA_genes$Sample)]
TCGA_genes
matched_queryGDC <- query_GDC$results[[1]][match(colnames(counts_DESeq_TCGA),
                                                 query_GDC$results[[1]]$id),]
table(matched_queryGDC$type) ## only gene expression, as expected
ggplot(data.frame(umap_TCGA$layout, mean_counts=log(rowMeans(counts_DESeq_TCGA))),
       aes(x=X1, y=X2, col=mean_counts))+geom_point()

df_umap_2 <- data.frame(umap_2$layout, mean_exprs=log(colMeans(counts_DESeq_TCGA)),
                        matched_clinical,WGD=matched_wgs,matched_queryGDC,
                        matched_exposures,
                        group_clr=matched_group_WGD,
                        genes=matched_group_genes,
                        consensusTME=t(matched_consensusTME))
colnames(df_umap_2)


## there seems to be a very clear split by the second umap component
ggplot(df_umap_2,
       aes(x=X1, y=X2, col=type))+geom_point() ## 
ggplot(df_umap_2,
       aes(x=X1, y=X2, col=mean_exprs))+geom_point() ## mostly seems to be due to the average expression
ggplot(df_umap_2,
       aes(x=X1, y=X2, col=vital_status))+geom_point() ## not separated by vital status
ggplot(df_umap_2,
       aes(x=X1, y=X2, col=WGD))+geom_point() ## not separated by WGD
ggplot(df_umap_2,
       aes(x=X1, y=X2, col=age_at_index))+geom_point() ## not separated by age
ggplot(df_umap_2,
       aes(x=X1, y=X2, col=as.numeric(days_to_birth)))+geom_point() ## again, not separate by age
ggplot(df_umap_2,
       aes(x=X1, y=X2, col=as.numeric(days_to_death)))+geom_point() ## not separated by life expectancy
ggplot(df_umap_2,
       aes(x=X1, y=X2, col=race))+geom_point() ## not separated by race
ggplot(df_umap_2,
       aes(x=X1, y=X2, col=race))+geom_point() ## not separated by race
ggplot(df_umap_2,
       aes(x=X1, y=X2, col=group_clr))+geom_point() ## not separated by WGD (clustering by clr)
ggplot(df_umap_2,
       aes(x=X1, y=X2, col=genes.Status))+geom_point() ## not separated by BRCA
ggplot(df_umap_2,
       aes(x=X1, y=X2, col=consensusTME.Immune_Score))+geom_point() ## 
ggplot(df_umap_2,
       aes(x=X1, y=X2, col=consensusTME.B_cells))+geom_point() ## 
ggplot(df_umap_2,
       aes(x=X1, y=X2, col=consensusTME.Fibroblasts))+geom_point() ## 


ggplot(cbind(df_umap_2, gene=counts_DESeq_TCGA['ENSG00000136997',]),
       aes(x=X1, y=X2, col=gene))+geom_point() ## not separated by MYC

ggplot(cbind(df_umap_2, gene=counts_DESeq_TCGA['ENSG00000175054',]),
       aes(x=X1, y=X2, col=gene))+geom_point() ## not separated by ATR
ggplot(cbind(df_umap_2, gene=counts_DESeq_TCGA['ENSG00000095585',]),
       aes(x=X1, y=X2, col=gene))+geom_point() ## not separated by BLNK
to_ens <- function(i){
  t2g$ensembl_gene_id[match(i, t2g$external_gene_name)]
}
ggplot(cbind(df_umap_2, gene=counts_DESeq_TCGA[to_ens('ESR1'),]),
       aes(x=X1, y=X2, col=gene))+geom_point() ## not separated by ER
ggplot(cbind(df_umap_2, gene=counts_DESeq_TCGA[to_ens('PGR'),]),
       aes(x=X1, y=X2, col=gene))+geom_point() ## not separated by PR
ggplot(cbind(df_umap_2, gene=counts_DESeq_TCGA[to_ens('RECQL'),]),
       aes(x=X1, y=X2, col=log(gene)))+geom_point() ## not separated by RECQL
ggplot(cbind(df_umap_2, gene=counts_DESeq_TCGA[to_ens('RNVU1-7'),]),
       aes(x=X1, y=X2, col=log(gene)))+geom_point()
ggplot(cbind(df_umap_2, gene=counts_DESeq_TCGA[to_ens('RNVU1-18'),]),
       aes(x=X1, y=X2, col=log(gene)))+geom_point()
ggplot(cbind(df_umap_2, gene=counts_DESeq_TCGA[to_ens('AL355075.4'),]),
       aes(x=X1, y=X2, col=log(gene)))+geom_point()
ggplot(cbind(df_umap_2, gene=counts_DESeq_TCGA[to_ens('AL355075.4'),]),
       aes(x=X1, y=X2, col=log(gene)))+geom_point()

ggplot(cbind(df_umap_2),
       aes(x=X1, y=X2, col=s1))+geom_point() ## not separated by s1
ggplot(cbind(df_umap_2),
       aes(x=X1, y=X2, col=s3))+geom_point() ## not separated by s3
ggplot(cbind(df_umap_2),
       aes(x=X1, y=X2, col=s4))+geom_point() ## not separated by s4

ggplot(data.frame(umap_2_with_orgs$layout, mean_counts=log(colMeans(joint_counts)),
                  group =c(rep('TCGA', ncol(counts_DESeq_TCGA)),
                           rep('Organoid', ncol(counts_DESeq_org_nonorm))),
                  label=c(rep(NA, ncol(counts_DESeq_TCGA)),
                          colnames(counts_DESeq_org_nonorm))),
       aes(x=X1, y=X2,
           # col=mean_counts
           col=group,
           label=label
       ))+geom_point()+geom_label_repel()

ggplot(data.frame(prcomp_with_orgs$x[,1:2], mean_counts=log(colMeans(joint_counts)),
           group =c(rep('TCGA', ncol(counts_DESeq_TCGA)),
                    rep('Organoid', ncol(counts_DESeq_org_nonorm))),
           label=c(rep(NA, ncol(counts_DESeq_TCGA)),
                   colnames(counts_DESeq_org_nonorm))), aes(x=PC1, y=PC2, col=group,
                                                     label=label))+
  geom_point()+geom_label()

dist_patients <- dist(df_umap_2)
pheatmap::pheatmap(dist_patients)

grouping_patients <- split(rownames(df_umap_2), f = df_umap_2$X2 < -2)

groups_umap <- cbind.data.frame(row.names=unlist(grouping_patients),
                 grouping_umap=rep(c('Group1','Group2'), sapply(grouping_patients, length)))
raw_counts_for_umap <- DESeq2::counts(`~response`)[,match(unlist(grouping_patients), colnames(DESeq2::counts(`~response`)))]
mat <- DESeqDataSetFromMatrix(countData=raw_counts_for_umap,
                              colData=groups_umap,
                              design=~grouping_umap)
mat <- estimateSizeFactors(mat)
mat <- suppressMessages(estimateDispersions(mat,fitType="local"))
results_DE_umap0 <- nbinomWaldTest(mat)
results_DE_umap <- DESeq2::results(results_DE_umap0, c("grouping_umap", "Group1", "Group2"),
                                      alpha = 0.05, format = "DataFrame")
results_DE_umap
results_DE_umap = cbind(t2g[match(gsub("\\..*","",rownames(raw_counts_for_umap)), t2g$ensembl_gene_id),], results_DE_umap)

vlcano(results_DE_umap)
table(results_DE_umap$padj < 0.00000001)
head(results_DE_umap[order(results_DE_umap$padj, decreasing = F),])
top_umap_DE <- change_rownames((DESeq2::counts(mat, normalized=F)[which(results_DE_umap$padj < 0.000000001),]))
pheatmap::pheatmap(top_umap_DE,
                   annotation_col = groups_umap, cluster_cols = T, scale = "row",
                   main = 'Genes with adjusted p-value < 0.05', show_colnames = FALSE)
top_umap_DE2 <- change_rownames(head(DESeq2::counts(mat, normalized=T)[order(results_DE_umap$padj),], n=35))
pheatmap::pheatmap(log(top_umap_DE2+0.001),
                   annotation_col = groups_umap, cluster_cols = F, scale = "row",
                   main = 'Genes with adjusted p-value < 0.05', show_colnames = FALSE)

pc1_loadings <- prcomp_with_orgs$rotation[,1]
pc12_loadings <- cbind.data.frame(PC1_loading=prcomp_with_orgs$rotation[,1],
                                  PC2_loading=prcomp_with_orgs$rotation[,2],
                                  ensembl=names(pc1_loadings),
                                  gene=t2g$external_gene_name[match(names(pc1_loadings),
                                                                    t2g$ensembl_gene_id)])
pc12_loadings$mostly_y = (abs(pc12_loadings$PC1_loading) > 0.011)
pc12_loadings$mostly_y_v2 = (abs(pc12_loadings$PC1_loading) > 0.012)

give_goterm_list_of_genes_general <- function(i, allgenes, nomenclature_subset='ensemble', nomeclature_all='entrez'){
  require(topGO)
  .x = i
  if(nomenclature_subset == 'ensemble'){
    .xentrez = t2g[match(.x, t2g$ensembl_gene_id),'entrezgene_id']
    .xentrez = .xentrez[!is.na(.xentrez)]
  }
  
  if(nomenclature_subset == 'ensemble'){
    allgenes = t2g[match(allgenes, t2g$ensembl_gene_id),'entrezgene_id']
    allgenes = as.character(allgenes[!is.na(allgenes)])
  }
  
  .xres_gosimenrich = GOSim::GOenrichment(genesOfInterest = as.character(.xentrez), allgenes = allgenes)
  .xres_gosim = GOSim::getGOInfo(.xentrez)
  .xres_gosimenrich
}

pc1_topgo <- give_goterm_list_of_genes_general(i = pc12_loadings[pc12_loadings$mostly_y_v2,]$ensembl,
                                               allgenes = pc12_loadings$ensembl,
                                               nomenclature_subset ='ensemble',
                                               nomeclature_all = 'ensemble')

ggplot(pc12_loadings)+geom_segment(aes(x = 0, xend=PC1_loading,
                                       y=0, yend=PC2_loading,
                                       col=mostly_y))


plot(log(rowMeans(counts_DESeq_org_nonorm)),
     log(rowMeans(counts_DESeq_TCGA[match(rownames(counts_DESeq_org_nonorm), rownames(counts_DESeq_TCGA)),])))

changes_average_genes <- cbind.data.frame(logmean_org=log(rowMeans(counts_DESeq_org_nonorm)),
                                          logmean_TCGA=log(rowMeans(counts_DESeq_TCGA[match(rownames(counts_DESeq_org_nonorm), rownames(counts_DESeq_TCGA)),])),
                                          ensemble=rownames(counts_DESeq_org_nonorm))
changes_average_genes$gene=t2g$external_gene_name[match(changes_average_genes$ensemble, t2g$ensembl_gene_id)]
changes_average_genes$difference = changes_average_genes$logmean_org - changes_average_genes$logmean_TCGA

changes_average_genes[duplicated(changes_average_genes$gene),]

ggplot(changes_average_genes, aes(x=factor(ensemble, levels=changes_average_genes$ensemble[order(changes_average_genes$difference)]),
                                  y=difference))+
  geom_point()

## DE between WGD and non-WGD (for organoids)
## but this is including the three samples with 3' bias
## also removing PDO17	xPDO17	JBLAB-19916
cluster_fig1 <- readRDS("../../copy_number_analysis_organoids/robjects/dendrograminputclr_tree.RDS")
t(t(cutree(cluster_fig1, k=2)[grepl('PDO', cluster_fig1$labels)]))
group_WGD <- cutree(cluster_fig1, k=2)[grepl('PDO', cluster_fig1$labels)]
group_WGD
## run 1_run_DE for analysis. This is the output file:

load("../objects/deaObjectFile_organoids_cluster")
DE_cluster = `~cluster`
results_DE_cluster <- DESeq2::results(DE_cluster, c("cluster", "1", "2"),
                                      alpha = 0.05, format = "DataFrame")
results_DE_cluster
results_DE_cluster = cbind(t2g[match(rownames(DESeq2::counts(DE_cluster)), t2g$ensembl_gene_id),], results_DE_cluster)

pdf("../figures/other_DE/DE_clusterWGD_org_volcano.pdf", height = 8, width = 8)
vlcano(results_DE_cluster)#+lims(y=c(0, 0.0001))
dev.off()


counts_DE_cluster <- DESeq2::counts(DE_cluster)
counts_DE_cluster[which(results_DE_cluster$padj < 0.05),]
results_DE_cluster[which(results_DE_cluster$padj < 0.05),]
hist(results_DE_cluster$pvalue)  

renaming <- readxl::read_excel("../files/PDOnameProperSample_sWGS_RNAseq.xlsx")
groups_cluster <- read.table("files/sample_sheet_organoids_cluster.csv",sep = ',', header = T)
colnames(counts_DE_cluster) <- renaming$PDO[match(colnames(counts_DE_cluster), renaming$sampleNameRNAseq)]
rownames(counts_DE_cluster) <- t2g$external_gene_name[match(rownames(counts_DE_cluster), t2g$ensembl_gene_id)]
groups_cluster$sample <- renaming$PDO[match(groups_cluster$sample, renaming$sampleNameRNAseq)]
rownames(groups_cluster) = groups_cluster$sample; groups_cluster$sample <- NULL

counts_DE_cluster <- counts_DE_cluster[,match(rownames(groups_cluster)[order(groups_cluster$cluster)], colnames(counts_DE_cluster))]
pdf("../figures/other_DE/DE_clusterWGD_org_highlyDE.pdf")
pheatmap::pheatmap(counts_DE_cluster[which(results_DE_cluster$padj < 0.05),],
                   annotation_col = groups_cluster, cluster_cols = F, scale = "row",
                   main = 'Genes with adjusted p-value < 0.05')
dev.off()


