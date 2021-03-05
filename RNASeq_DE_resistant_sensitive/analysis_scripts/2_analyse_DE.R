
rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(rnaseqRpkg)
library(biomaRt)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(DESeq2)
library(reshape2)

## Load file output from DESeq2, which has been created with 1_run_DE.sh
load("../objects/deaObjectFile")
deObj = `~response`

##---------------------------------------------------------------------------------------------------------------##
## Run differential expression with DESeq
# In any case, the contrast argument of the function results takes a character vector of length three: the name of the variable, the name of the factor level for the numerator of the log2 ratio, and the name of the factor level for the denominator. The contrast argument can also take other forms, as described in the help page for results and below
results <- DESeq2::results(deObj, c("response", "complete_remission_or_response", "progressive_disease"), 
                  alpha = 0.05, format = "DataFrame")
rownames_short = sapply(rownames(results), function(i) strsplit(i, '[.]')[[1]][1])

## Re-name
# mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
# gene_conversion <- getBM(attributes = c("external_gene_name", "ensembl_gene_id"),
#                  filters = "ensembl_gene_id", values = rownames_short,
#                  mart = mart)
# gene_conversion = gene_conversion[match(rownames_short, gene_conversion$ensembl_gene_id),]
t2g = readRDS("~/Desktop/t2g.RDS")
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
## e.g. JBLAB-19920 is PDO18 and has been excluded based on 3' bias
"JBLAB-19920" %in% colnames(SampleGroup)
length(colnames(SampleGroup)) ## 3' bias have been excluded
## but I am not sure why there are 14 instead of 15
"JBLAB-19939" %in% colnames(SampleGroup) ## the fourth sample with highest 3' bias

DE_results_org = DESeq2::results(SampleGroup)
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
counts_DESeq_org = read.table("../files/counts_norm.csv", sep = ",", header=T, stringsAsFactors = FALSE)
apply(counts_DESeq_org, 2, function(i) as.numeric(i))
rownames(counts_DESeq_org) = counts_DESeq_org[,1]; counts_DESeq_org = counts_DESeq_org[,-1]
colnames(counts_DESeq_org) ## we have more than 18 because we also have normal samples
## Remove the organoids that have a 3' bias, and the normal samples
## PDO14 (JBLAB-19925) PDO16 (JBLAB-19936) PDO18 (JBLAB-19920) and normal (JBLAB-19950, JBLAB-19953, JBLAB-19952
## JBLAB-19954, JBLAB-19955, JBLAB-19951)
counts_DESeq_org = counts_DESeq_org[,!(colnames(counts_DESeq_org) %in% c('JBLAB.19925', 'JBLAB.19936', 'JBLAB.19920', 'JBLAB.19950',
                                                                         'JBLAB.19953', 'JBLAB.19952',
                                                                         'JBLAB.19954', 'JBLAB.19955', 'JBLAB.19951'))]
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

ggplot(df_all_genes, aes(y=value, x=L1))+geom_boxplot()+geom_jitter()

means_tcga <- rowMeans(counts_DESeq_TCGA)
means_org <- rowMeans(counts_DESeq_org)

pdf("../figures/Sensitive_resistant_figures/colmeans_deseqcounts_correlation_tcga_org.pdf")
plot(log(means_tcga), log(means_org), xlab='DESeq normalised counts from TCGA',ylab='DESeq normalised counts from organoids')
abline(coef=c(0,1), lty='dashed', col='blue')
dev.off()

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



