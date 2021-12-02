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
library(jcolors)
library(parallel)
library(latex2exp)
library(GSVA)
library(GSVAdata)
library(DESeq2)
##------------------------------------------------------------------------------------------------------------##

## Note: any plots containing nearestGeneCN/cordfAll is wrong, as it had been computed with the incorrect organoid renaming

include_14_orgs = FALSE
include_11_orgs = TRUE

##------------------------------------------------------------------------------------------------------------##
##' 3' bias samples to remove
renaming = readxl::read_xlsx("../../../RNASeq_DE_resistant_sensitive/files/PDOnameProperSample_sWGS_RNAseq.xlsx")
if(include_14_orgs){
  remove_samples = (c('PDO14', 'PDO16', 'PDO18'))
}else if(include_11_orgs){
  remove_samples = (c('PDO14', 'PDO16', 'PDO18',
                           'PDO13', 'PDO4', 'PDO9',
                           'PDO17', 'FT10', 'FT7', 'FT7-bis'))
}
remove_samples_with_normal_fal = cbind.data.frame(PDO=remove_samples,
                                  ID=renaming$ID[match(remove_samples, renaming$PDO)],
                                  RNASeq=renaming$sampleNameRNAseq[match(remove_samples, renaming$PDO)])
remove_samples = remove_samples_with_normal_fal[grepl('PDO', remove_samples_with_normal_fal$PDO),]
##------------------------------------------------------------------------------------------------------------##

##------------------------------------------------------------------------------------------------------------##
joint_counts_CN0 = readRDS("../output/output_GRCh37_with_14_orgs/joint_counts_CN_TPM_20210506104217.RDS")
joint_counts_CN = joint_counts_CN0
joint_counts_CN = joint_counts_CN %>% filter( !(counts.Var1 %in% remove_samples[2,]))
topvariable = read.table("../../top_variable.txt", comment.char = "#")

## now add DESeq counts
## recompute them, depending on which organoids we are including

if(include_14_orgs){
  deseq_counts0 = read.table("../../../RNASeq_DE_resistant_sensitive/files/counts_norm.csv",
                            sep=',', header = T)
  deseq_counts = deseq_counts0
  # deseq_counts_renormalised: without the 3'- biased samples
  deseq_counts_renormalised0 <- read.table("../../../RNASeq_and_CN/20191218_ViasM_BJ_orgaBrs/RnaSeqPip/counts/counts_raw.csv", sep = ",", header = T)
  deseq_counts_renormalised0 <- deseq_counts_renormalised0[,!(renaming$PDO[match(colnames(deseq_counts_renormalised0), gsub("-", ".", renaming$sampleNameRNAseq))] %in% remove_samples[1,])]
  rownames(deseq_counts_renormalised0) <- deseq_counts_renormalised0[,1]; deseq_counts_renormalised0[,1] <- NULL
  deseq_counts_renormalised0 <- DESeq2::DESeqDataSetFromMatrix(deseq_counts_renormalised0, colData=data.frame(org_bool=grepl('PDO', renaming$PDO[match(colnames(deseq_counts_renormalised0), gsub("-", ".", renaming$sampleNameRNAseq))])), design = ~org_bool)
  deseq_counts_renormalised0 <- estimateSizeFactors(deseq_counts_renormalised0)
  deseq_counts_renormalised0 <- estimateDispersions(deseq_counts_renormalised0,fitType="local")
  deseq_counts_renormalised0 <- DESeq2::counts(deseq_counts_renormalised0, normalized=T)
}else if(include_11_orgs){
  raw_counts0 = read.csv("../../../RNASeq_and_CN/20191218_ViasM_BJ_orgaBrs/RnaSeqPip/counts/counts_raw.csv", stringsAsFactors = FALSE)
  rownames(raw_counts0) = raw_counts0[,1]; raw_counts0 <- raw_counts0[,-1]
  colnames(raw_counts0) <- renaming$PDO[match(gsub('[.]', '-', colnames(raw_counts0)),
                                              renaming$sampleNameRNAseq)]
  raw_counts0 = raw_counts0[,!(colnames(raw_counts0) %in% remove_samples_with_normal_fal$PDO)]
  ## here we normalising with only 3 FT samples without 3' bias
  renormalised_counts_obj <- DESeq2::DESeqDataSetFromMatrix(countData = raw_counts0,
                                                colData = cbind.data.frame(Sample=colnames(raw_counts0),
                                                                           Group=grepl('PDO', colnames(raw_counts0))),
                                                design = ~ Group)
  renormalised_counts_obj <- DESeq2::estimateSizeFactors(renormalised_counts_obj)
  saveRDS(renormalised_counts_obj, "../../20191218_ViasM_BJ_orgaBrs/output/fig3_renormalised_counts_obj_11orgs.RDS")
  renormalised_counts <- DESeq2::counts(renormalised_counts_obj, normalized=T)
  deseq_counts_renormalised0 <- renormalised_counts
  deseq_counts <- deseq_counts_renormalised0
}else{
  stop('set include_11_orgs or include_14_orgs to T')
}
deseq_counts = (melt(deseq_counts))
colnames(deseq_counts) = c('Gene', 'DESEq.sample', 'DESeq.count')

## Gene name conversion
# mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl", host="grch37.ensembl.org", ))
# # # ensembl <- useEnsembl(biomart = "genes")
# ensembl <- useDataset(dataset = "hsapiens_gene_ensembl", mart = mart)
# 
# t2g <- getBM(
#   attributes = c('ensembl_gene_id', 'external_gene_name'),
#   values = deseq_counts$Gene,
#   filter = 'ensembl_gene_id',
#   mart = ensembl, useCache = FALSE)
# saveRDS(t2g, file = "~/Desktop/t2g_GRCh37.RDS")
t2g_37 = readRDS("../../../copy_number_analysis_organoids/robjects/t2g_GRCh37.RDS")
t2g_38 = readRDS("../../../copy_number_analysis_organoids/robjects/t2g2.RDS")
EnsDb.Hsapiens.v87 <- ensembldb::EnsDb("../../../RNASeq_and_CN/20191218_ViasM_BJ_orgaBrs/Data/Homo_sapiens.GRCh37.87.sqlite")
ag <- ensembldb::genes(EnsDb.Hsapiens.v87, filter=list(AnnotationFilter::GeneBiotypeFilter('protein_coding')), return.type="DataFrame")

genename_ensmid_matched = t2g_38[match(deseq_counts$Gene, t2g_38$ensembl_gene_id),'external_gene_name']
deseq_counts$Gene = genename_ensmid_matched
deseq_counts = deseq_counts[!is.na(deseq_counts$Gene),]

if(include_14_orgs){
  deseq_counts$DESEq.sample = renaming$ID[match(gsub("[.]", "-", deseq_counts$DESEq.sample), renaming$sampleNameRNAseq)]
}else if(include_11_orgs){
  deseq_counts$DESEq.sample = renaming$ID[match(gsub("[.]", "-", deseq_counts$DESEq.sample), renaming$PDO)]
}

joint_counts_CN = cbind.data.frame(joint_counts_CN, deseq_counts[match(paste0(joint_counts_CN$counts.Var1, joint_counts_CN$CN.gene_name), paste0(deseq_counts$DESEq.sample, deseq_counts$Gene)),])

plot(joint_counts_CN$counts.value, joint_counts_CN$DESeq.count)

# saveRDS(joint_counts_CN, file = "../output/joint_counts_CN.RDS")

corDfAll_do_not_use_values = read.csv("../GexVsCnGwDeseq2/corStatTable.csv")
# corDfAll_do_not_use_values = readRDS("../output/corDfAll.RDS") ## this is likely to be with incorrect PDO

joint_counts_CN$PDO = renaming$PDO[match(joint_counts_CN$counts.Var1, renaming$ID)]
## remove samples with 3' bias
joint_counts_CN = joint_counts_CN[!(joint_counts_CN$PDO %in% remove_samples$PDO),]
# which_all_na <- apply(apply(joint_counts_CN, 1, is.na), 2, all)

# table(is.na(joint_counts_CN_TPM_subset$counts.value))
# table(is.na(joint_counts_CN_TPM_subset$counts.DESeq.count))

# length(table(corDfAll$GeneName))
# length(table(joint_counts_CN_TPM_subset$nearestGeneCN.gene)) ## ?? we're losing genes?
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

# joint_counts_CN_subset = cbind.data.frame(nearestGeneCN=corDfAll_df,
#    joint_counts_CN[match(paste0(corDfAll_df$gene, corDfAll_df$variable),
#                          paste0(joint_counts_CN$counts.Var2, joint_counts_CN$counts.Var1)),])

joint_counts_CN_subset = joint_counts_CN %>% filter(CN.gene_name %in% corDfAll_do_not_use_values$GeneName)

saveRDS(joint_counts_CN_subset, "../output/output_GRCh37/joint_counts_CN_subset.RDS")
# joint_counts_CN_TPM_subset = readRDS("../output/joint_counts_CN_TPM_subset.RDS")
# joint_counts_CN_TPM_subset = joint_counts_CN %>% filter( !(counts.Var1 %in% remove_samples[2,]))
# joint_counts_CN_TPM_subset = cbind(joint_counts_CN_TPM_subset,
#                                    counts=joint_counts_CN[match(paste0(joint_counts_CN_TPM_subset$nearestGeneCN.variable, joint_counts_CN_TPM_subset$counts.Var2),
#       paste0(joint_counts_CN$counts.Var1, joint_counts_CN$counts.Var2)),])
# 
# plot(joint_counts_CN_TPM_subset$counts.value, joint_counts_CN_TPM_subset$counts.DESeq.count) ## we have added DESeq2
# plot(joint_counts_CN_TPM_subset$counts.value, joint_counts_CN_TPM_subset$counts.counts.value)

## coding
# ensembl <- useMart("ensembl",
#                    dataset = "hsapiens_gene_ensembl")
# coding_genes = getBM(attributes=c("ensembl_gene_id","external_gene_name","description"),
#                      filters='biotype', values=c('protein_coding'), mart=ensembl)
# saveRDS(coding_genes, "~/Desktop/coding_genes.RDS")
coding_genes <- readRDS("../../../copy_number_analysis_organoids/robjects/coding_genes.RDS")
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

pdf("../../figures_GRCh37/joint_counts_CN_subset_remove_outliers.pdf")
plot(sort(log(joint_counts_CN_subset$nearestGeneCN.value)),
     col=as.factor(joint_counts_CN_subset$flag_outlier_nearestGeneCN.value[order(log(joint_counts_CN_subset$nearestGeneCN.value))]),
     main=paste0('Removing ', probs_outlier, ' percentile'))
dev.off()

subset_genes_of_interest = c('MYC', 'CCNE1', 'PIK3CA', 'TERT', 'KRAS', 'PTEN', 'RB1', 'AKT1',
                             'AKT2', 'PARP1', 'PARP2', 'ATM', 'ATR', 'WEE1', 'TOP1', 'TUBB1',
                             'AKT3', 'CCND1', 'CCND2', 'CCND3', 'CDKN2A', 'CDKN2B', 'MECOM', 'CDK12')

# joint_counts_CN_subset = joint_counts_CN_subset[!joint_counts_CN_subset$flag_outlier_nearestGeneCN.value,]
# plot(sort(log(joint_counts_CN$CN.value)))
joint_counts_CN_subset = joint_counts_CN_subset %>%
  group_by(CN.gene_name) %>%
  mutate(scaled_centered_weighted_CN = exp(scale(log(CN.value))),
         scaled_centered_weighted_CN_mean = (scale((CN.value))),
         scaled_centered_DESeq = scale(DESeq.count),
         scaled_centered_TPM = scale(counts.value)
         # scaled_centered_nearestCN = exp(scale(log(nearestGeneCN.value))),
         # scaled_centered_nearestCN_mean = (scale((nearestGeneCN.value)))
  )
# joint_counts_CN_subset = joint_counts_CN_subset[!(joint_counts_CN_subset$nearestGeneCN.variable %in% remove_samples[2,]),]
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

## genes differentially abundant between normal tissue and cancer tissue
## colour these genes
## this is including samples with 3' bias, however
# load("../RnaSeqPip/DEAnalysis/deObject_SampleGroup.RData")
# results_deseq = DESeq2::results(SampleGroup)
## using the newly computed DESeq object, only with non-biased samples
renormalised_counts_obj = DESeq2::DESeq(renormalised_counts_obj)
results_deseq = DESeq2::results(renormalised_counts_obj)

top_genes_normal_tumour = rownames(results_deseq)[order(results_deseq$padj, decreasing = F)[1:400]]
top_genes_normal_tumour = cbind.data.frame(ensembl=top_genes_normal_tumour, gene=t2g_38$external_gene_name[match(top_genes_normal_tumour, t2g_38$ensembl_gene_id)])

joint_counts_CN$topDiff_norm_tissue = ifelse(joint_counts_CN$CN.gene_name %in% top_genes_normal_tumour$gene,
                                             yes = 'DE', no = 'nonDE' )
##------------------------------------------------------------------------------------------------------------##
corr_centered =  joint_counts_CN[,!duplicated(colnames(joint_counts_CN))] %>%
  group_by(CN.gene_name) %>% summarise(cor(scaled_centered_weighted_CN, scaled_centered_DESeq))
r2_centered =  joint_counts_CN %>% group_by(CN.gene_name) %>% summarise( if(all(is.na(scaled_centered_DESeq)) | all(is.na(scaled_centered_weighted_CN))){NA}else{ try(summary(lm( scaled_centered_DESeq ~ scaled_centered_weighted_CN))$r.squared)} )
r2_raw_DESeq2 =  joint_counts_CN %>% group_by(CN.gene_name) %>% summarise( if(all(is.na(DESeq.count)) | all(is.na(CN.value))){NA}else{ try(summary(lm( DESeq.count ~ CN.value))$r.squared)} )
slope_centered =  joint_counts_CN %>% group_by(CN.gene_name) %>% summarise( if(all(is.na(scaled_centered_DESeq)) | all(is.na(scaled_centered_weighted_CN))){NA}else{ try(summary(lm( scaled_centered_DESeq ~ scaled_centered_weighted_CN))$coefficients[2])} )
var_CN =  joint_counts_CN %>% group_by(CN.gene_name) %>% summarise( var(CN.value) )
var_DESeq2 =  joint_counts_CN %>% group_by(CN.gene_name) %>% summarise( var(DESeq.count) )
colnames(r2_centered) = c('Gene', 'r2')
colnames(r2_raw_DESeq2) = c('Gene', 'r2')
top_cor = corr_centered$CN.gene_name[order(corr_centered$`cor(scaled_centered_weighted_CN, scaled_centered_DESeq)`, decreasing = T)[1:30]]
top_r2 = r2_centered$Gene[order(r2_centered$r2, decreasing = T)[1:30]]

plot(r2_centered$r2, r2_raw_DESeq2$r2)


##------------------------------------------------------------------------------------------------------------##

logCN = joint_counts_CN_subset %>% group_by(PDO) %>% summarise(mean=mean(log(CN.value)))

# example_centering = joint_counts_CN_normalised %>% filter(counts.Var2 == 'NOC2L') %>% dplyr::select(scaled_centered_weighted_CN)
# mean(example_centering$scaled_centered_weighted_CN); sd(example_centering$scaled_centered_weighted_CN)
# table(joint_counts_CN_ranked$rank_weighted_CN)

## there are genes which are duplicated - remove them
joint_counts_CN_ranked = joint_counts_CN %>% dplyr::select(c("counts.Var1", "counts.Var2", "counts.value", "rank_weighted_CN", "rank_DESeq", "rank_TPM"))
joint_counts_CN_ranked = joint_counts_CN_ranked[!(joint_counts_CN_ranked$CN.gene_name %in% joint_counts_CN_ranked$CN.gene_name[joint_counts_CN_ranked$rank_weighted_CN > length(unique(joint_counts_CN_ranked$counts.Var1))]),]
table(joint_counts_CN_ranked$rank_weighted_CN)

##------------------------------------------------------------------------------------------------------------##

joint_counts_CN_ranked_subsetorgs = joint_counts_CN_ranked %>% filter(! (counts.Var1 %in% c('118976org', '119127org', '54327org', '119178org', '118947org')))
## if a gene has a copy number greater than two, is its expression higher than in normal cells?
df_gene_characteristics = cbind.data.frame(cor_normCNnormDESeq=corr_centered$`cor(scaled_centered_weighted_CN, scaled_centered_DESeq)`,
                                             var_DESeq=var_DESeq2$`var(DESeq.count)`,
                                             r2_normCNnormDESeq=r2_centered$r2,
                                             r2_raw_DESeq2=r2_raw_DESeq2,
                                             var_CN=var_CN$`var(CN.value)`,
                                             slope_centered=slope_centered,
                                             Gene=as.character(corr_centered$CN.gene_name),
                                           joint_counts_CN %>% group_by(CN.gene_name) %>%
                                             summarise(averageCN=mean(CN.value),
                                            sdCN=sd(CN.value),
                                            averageDESeq2=mean(DESeq.count),
                                            sdDESeq2=sd(DESeq.count)
                                  ))

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
##------------------------------------------------------------------------------------------------------------##

mean_coord_CN_DESeq2 = joint_counts_CN %>% dplyr::group_by(Gene) %>% dplyr::summarise(mean_coord_CN_DESeq2=mean(DESeq2_increase_withCN, na.rm=T))
mean_coord_CN_DESeq2 = mean_coord_CN_DESeq2[order(mean_coord_CN_DESeq2$mean_coord_CN_DESeq2, decreasing = T),]
mean_coord_CN_DESeq2$rank = factor(1:nrow(mean_coord_CN_DESeq2))

saveRDS(mean_coord_CN_DESeq2, "../output/output_GRCh37/fig3_mean_coord_CN_DESeq2.RDS")

##' for each gene, average the bottom three organoids (sorted by CN) and check if
##'  the expression of the remaining is higher than for those

if(include_14_orgs){
  num_orgs <- 15
}else if(include_11_orgs){
  num_orgs <- ncol(deseq_counts_renormalised0[,grepl('PDO', colnames(deseq_counts_renormalised0))])
}else{
  stop()
}

joint_counts_CN_cleanDESeq = joint_counts_CN[!(is.na(joint_counts_CN$Gene)),]
give_average_comparison_CN_DESeq <- function(i, num_bottom_orgs=3){
  .x = joint_counts_CN_cleanDESeq[joint_counts_CN_cleanDESeq$counts.Var2 == i,]
  if(!( (!any(is.na(.x$DESeq.count))) & (length(.x$counts.Var1) == num_orgs) )){
    return(NA)
  }else{
    .averg_DESeq = mean(unlist(.x[order(.x$CN.value, decreasing = F)[1:num_bottom_orgs],'DESeq.count']))
    return(.x[order(.x$CN.value, decreasing = F)[-(1:num_bottom_orgs)],'DESeq.count'] > .averg_DESeq)
  }
}

give_average_comparison_CN_TPM <- function(i, num_bottom_orgs=3){
  .x = joint_counts_CN_cleanDESeq[joint_counts_CN_cleanDESeq$counts.Var2 == i,]
  if(!( (!any(is.na(.x$DESeq.count))) & (length(.x$counts.Var1) == 15) )){
    return(NA)
  }else{
    .averg_TPM = mean(unlist(.x[order(.x$CN.value, decreasing = F)[1:num_bottom_orgs],'counts.value']))
    return(.x[order(.x$CN.value, decreasing = F)[-(1:num_bottom_orgs)],'counts.value'] > .averg_TPM)
  }
}


average_comparison_CN_DESeq = mclapply(unique(joint_counts_CN$counts.Var2),
                                     give_average_comparison_CN_DESeq)
average_comparison_CN_DESeq_2 = mclapply(unique(joint_counts_CN$counts.Var2),
                                       give_average_comparison_CN_DESeq, num_bottom_orgs=2)
average_comparison_CN_DESeq_4 = mclapply(unique(joint_counts_CN$counts.Var2),
                                         give_average_comparison_CN_DESeq, num_bottom_orgs=4)
average_comparison_CN_TPM = mclapply(unique(joint_counts_CN$counts.Var2),
                                       give_average_comparison_CN_TPM)
## in general we would expect the value to be positive

length(unique(joint_counts_CN_subset$counts.Var2))
average_comparison_CN_DESeq[[1]]

sum(is.na(average_comparison_CN_DESeq))
df_average_bottomCN <- cbind.data.frame(Gene=unique(joint_counts_CN$counts.Var2),
                                        average_comparison_CN_DESeq=sapply(average_comparison_CN_DESeq, mean),
                                        average_comparison_CN_DESeq2=sapply(average_comparison_CN_DESeq_2, mean),
                                        average_comparison_CN_DESeq4=sapply(average_comparison_CN_DESeq_4, mean),
                                        average_comparison_CN_TPM = sapply(average_comparison_CN_TPM, mean),
                                        label=ifelse(unique(joint_counts_CN$counts.Var2) %in% subset_genes_of_interest,
                                                     yes = as.character(unique(joint_counts_CN$counts.Var2)), no = NA))
df_average_bottomCN$Gene <- factor(df_average_bottomCN$Gene, levels=df_average_bottomCN$Gene[order(df_average_bottomCN$average_comparison_CN_DESeq, na.last = T)])
saveRDS(df_average_bottomCN, "../output/output_GRCh37/fig3_df_average_bottomCN.RDS")
# df_average_bottomCN <- readRDS("../output/output_GRCh37/fig3_df_average_bottomCN.RDS")

## For each gene, plot the CN

pairs(df_average_bottomCN[,c('average_comparison_CN_DESeq', 'average_comparison_CN_DESeq2', 'average_comparison_CN_DESeq4')])

median_CN = aggregate(data = joint_counts_CN %>% filter(CN.gene_name %in% topvariable$V1),
                      CN.value ~ CN.gene_name, FUN = median)
median_CN_goi = aggregate(data = joint_counts_CN %>% filter(counts.Var2 %in% subset_genes_of_interest),
                          CN.value ~ counts.Var2, FUN = median)

df_gene_characteristics = cbind.data.frame(as.data.frame(df_gene_characteristics),
                                           df_average_bottomCN = data.frame(df_average_bottomCN[match(df_gene_characteristics$Gene,
                                                                                           df_average_bottomCN$Gene),])
                                           )
# saveRDS(df_gene_characteristics, "../output/fig3_df_gene_characteristics.RDS")
saveRDS(df_gene_characteristics, "../output/output_GRCh37/fig3_df_gene_characteristics.RDS")
# df_gene_characteristics <- readRDS("../output/output_GRCh37/fig3_df_gene_characteristics.RDS")


##------------------------------------------------------------------------------------------------------------##
## General plots
# ggplot(joint_counts_CN %>% filter(CN.L1 == '118947org' ), aes(x=CN.value, y=counts.value))+geom_point()

## per organoid separately in facets
ggplot(joint_counts_CN, aes(x=CN.value, y=counts.value))+geom_point()+facet_wrap(.~CN.L1, scales = "free")+theme_bw()
ggsave("../../figures_GRCh37/joint_counts_CN_TPM_all.pdf", width=10, height=10)

## same, log scale
ggplot(joint_counts_CN, aes(x=CN.value, y=counts.value))+geom_point()+facet_wrap(.~CN.L1, scales = "free")+theme_bw()+
  scale_x_continuous(trans = "log2")+scale_y_continuous(trans = "log2")
ggsave("../../figures_GRCh37//joint_counts_CN_TPM_all_loglog.pdf", width=10, height=8)

## only most variable genes (when it comes to neighbouring CN)
# ggplot(joint_counts_CN %>% filter(CN.gene_name %in% topvariable$V1), aes(x=CN.value, y=counts.value, col=CN.gene_name))+
  # geom_point()

# pdf("../../figures/joint_counts_CN_TPM_topvar.pdf")
# for(i in topvariable$V1){
#   print(ggplot(joint_counts_CN %>% filter(CN.gene_name == i), aes(x=CN.value, y=counts.value, label=CN.L1))+
#           geom_point()+geom_smooth()+ggtitle(i)+geom_label_repel()+theme_bw())
# }
# dev.off()
##------------------------------------------------------------------------------------------------------------##

##------------------------------------------------------------------------------------------------------------##

ggplot(joint_counts_CN, aes(x=CN.value, y=DESeq.count))+geom_point()+facet_wrap(.~CN.L1, scales = "free")+theme_bw()
ggsave("../../figures_GRCh37/joint_counts_CN_DESeq_all.pdf", width=10, height=10)

ggplot(joint_counts_CN, aes(x=CN.value, y=DESeq.count))+geom_point()+facet_wrap(.~CN.L1, scales = "free")+theme_bw()+
  scale_x_continuous(trans = "log2")+scale_y_continuous(trans = "log2")
ggsave("../../figures_GRCh37//joint_counts_CN_DESeq_all_loglog.pdf", width=10, height=8)

# pdf("../../figures/joint_counts_CN_DESeq_topvar.pdf")
# for(i in topvariable$V1){
#   print(ggplot(joint_counts_CN %>% filter(CN.gene_name == i), aes(x=CN.value, y=DESeq.count, label=CN.L1))+
#           geom_point()+geom_smooth()+ggtitle(i)+geom_label_repel()+theme_bw())
# }
# dev.off()

## raw correlation between CN and deseq count
plot(joint_counts_CN$counts.value, joint_counts_CN$DESeq.count)

##------------------------------------------------------------------------------------------------------------##

## Looking at genes in specific
# ggplot(joint_counts_CN %>% filter(CN.gene_name == "CCNE1"), aes(x=CN.value, y=counts.value, label=CN.L1))+
#   geom_point()+geom_smooth()+ggtitle("CCNE1")+geom_label_repel()+theme_bw()


# ggplot(joint_counts_CN %>% filter(CN.gene_name == "CDK12"), aes(x=CN.value, y=counts.value, label=CN.L1))+
#   geom_point()+geom_smooth()+ggtitle("CDK12")+geom_label_repel()+theme_bw()


##------------------------------------------------------------------------------------------------------------##

## only genes included in Stephane's analysis
## length(table(joint_counts_CN_TPM_subset$nearestGeneCN.gene))
# ggplot(joint_counts_CN %>% filter(CN.gene_name %in% joint_counts_CN_TPM_subset$nearestGeneCN.gene),
#        aes(x=CN.value, y=DESeq.count))+geom_point()+facet_wrap(.~CN.L1, scales = "free")+theme_bw()+
#   scale_x_continuous(trans = "log2")+scale_y_continuous(trans = "log2")+
#   geom_abline(slope = 1, intercept = 0, col='blue', lty='dashed')

# ggplot(joint_counts_CN %>% filter(DESeq.count > 0) %>%
#          filter(CN.gene_name %in% joint_counts_CN_TPM_subset$nearestGeneCN.gene),
#        aes(x=CN.value, y=DESeq.count))+geom_point()+theme_bw()+
#   scale_x_continuous(trans = "log2")+scale_y_continuous(trans = "log2")+
#   geom_abline(slope = 1, intercept = 0, col='blue', lty='dashed')+geom_smooth()


# ggplot(joint_counts_CN %>% filter(DESeq.count > 0) %>%
#          filter(CN.gene_name %in% joint_counts_CN_TPM_subset$nearestGeneCN.gene),
#        aes(x=CN.value, y=DESeq.count))+geom_point()+theme_bw()+
#   scale_x_continuous(trans = "log2")+scale_y_continuous(trans = "log2")+
#   geom_abline(slope = 1, intercept = 0, col='blue', lty='dashed')+geom_smooth()

# ggplot(high_GEvar,
#        aes(x=CN.value, y=DESeq.count))+geom_point()+theme_bw()+
#   scale_x_continuous(trans = "log2")+scale_y_continuous(trans = "log2")+
#   geom_abline(slope = 1, intercept = 0, col='blue', lty='dashed')+geom_smooth()

ggplot(joint_counts_CN_subset, aes(x=CN.value, y=nearestGeneCN.value,
       col=nearestGeneCN.gene %in% subset_genes_of_interest))+geom_point(alpha=0.2)+
  labs(x='Weighted CN value of gene (Lena)', y='CN value of gene highly variable region nearby')+
  facet_wrap(.~(nearestGeneCN.gene %in% subset_genes_of_interest))+
  scale_x_continuous(trans = "log2")+scale_y_continuous(trans = "log2")+guides(col=F)
ggsave("../../figures_GRCh37/scatterplot_stephane_lena.pdf", width = 7, height = 5)

## CN vs count for genes in highly variable areas
# ggplot(joint_counts_CN_subset, aes(x=CN.value, y=counts.value))+geom_point()+
#   scale_y_continuous(trans='log2')+geom_smooth()
## nearest segment vs count
# ggplot(joint_counts_CN_subset %>% filter(counts.value > 0), aes(x=nearestGeneCN.value, y=counts.value))+geom_point()+
#   scale_y_continuous(trans='log2')+geom_smooth()+lims(x=c(0,10))
## nearest segment vs count, excluding normal
ggplot()+
  geom_point(data = joint_counts_CN_subset %>% filter(counts.value > 1, abs(nearestGeneCN.value - 2) > 0.1) %>%
               dplyr::select(-labels),
             aes(x=nearestGeneCN.value, y=counts.value), alpha=0.2)+
  geom_point(data = joint_counts_CN_subset %>% filter(counts.value > 0, abs(nearestGeneCN.value - 2) > 0.1) %>% filter(!is.na(labels)),
             aes(x=nearestGeneCN.value, y=counts.value, label=labels, col=labels))+
  scale_y_continuous(trans='log2')+geom_smooth()+lims(x=c(0,10))+geom_label_repel()+facet_wrap(.~labels)+
  guides(col=FALSE)+theme_bw()
ggsave("../../figures_GRCh37/selected_genes_nearestgene_TPM_nonorm.png", width = 7, height = 5)

ggplot()+
  geom_point(data = joint_counts_CN_subset %>% filter(counts.value > 1, abs(CN.value - 2) > 0.1) %>%
               dplyr::select(-labels),
             aes(x=CN.value, y=counts.value), alpha=0.2)+
  geom_point(data = joint_counts_CN_subset %>% filter(counts.value > 0, abs(CN.value - 2) > 0.1) %>% filter(!is.na(labels)),
             aes(x=CN.value, y=counts.value, label=labels, col=labels))+
  scale_y_continuous(trans='log2')+geom_smooth()+lims(x=c(0,10))+geom_label_repel()+facet_wrap(.~labels)+
  guides(col=FALSE)+theme_bw()
ggsave("../../figures_GRCh37/selected_genes_weightedCN_TPM_nonorm.png", width = 7, height = 5)

# ggplot()+
#   geom_point(data = joint_counts_CN_subset %>% filter(counts.value > 0, abs(nearestGeneCN.value - 2) > 0.1) %>%
#                dplyr::select(-labels),
#              aes(x=scaled_centered_nearestCN, y=scaled_centered_TPM), alpha=0.2)+
#   geom_point(data = joint_counts_CN_subset %>% filter(counts.value > 0, abs(nearestGeneCN.value - 2) > 0.1) %>% filter(!is.na(labels)),
#              aes(x=scaled_centered_nearestCN, y=scaled_centered_TPM, label=labels, col=labels))+
#   geom_smooth()+geom_label_repel()+facet_wrap(.~labels)+
#   guides(col=FALSE)+theme_bw()

ggplot()+
  geom_point(data = joint_counts_CN_subset %>% filter(counts.value > 1, abs(CN.value - 2) > 0.1) %>%
               dplyr::select(-labels),
             aes(x=scaled_centered_weighted_CN, y=scaled_centered_DESeq), alpha=0.2)+
  geom_point(data = joint_counts_CN_subset %>% filter(counts.value > 0, abs(CN.value - 2) > 0.1) %>% filter(!is.na(labels)),
             aes(x=scaled_centered_weighted_CN, y=scaled_centered_DESeq, label=labels, col=labels))+
  geom_smooth()+geom_label_repel()+facet_wrap(.~labels)+
  guides(col=FALSE)+theme_bw()+scale_x_continuous(trans="log2")
ggsave("../../figures_GRCh37/selected_genes_scaled_weightedCN_scaledDESeq_nonorm.png", width = 7, height = 5)

## it doesn't look centered to me? e.g. MYC doesn't seem to be centered for CN
# (joint_counts_CN_subset %>% filter(nearestGeneCN.gene == 'MYC') %>% dplyr::select(scaled_centered_weighted_CN))[,2] %>% unlist  %>% mean
# (joint_counts_CN_subset %>% filter(nearestGeneCN.gene == 'MYC') %>% dplyr::select(scaled_centered_weighted_CN))[,2] %>% unlist  %>% sd
## it is. there is just a bunch of observations in the same CN (this is what happens if you don't remove outliers)

# ggplot(joint_counts_CN_subset %>% filter(counts.value > 0), aes(x=CN.value, y=counts.value))+geom_point()+
#   scale_y_continuous(trans='log2')+geom_smooth()+lims(x=c(0,10))

# ggplot(data = joint_counts_CN_subset %>% filter(counts.value > 0, abs(nearestGeneCN.value - 2) > 0.1) %>%
#                filter(nearestGeneCN.gene == 'MYC') %>%
#                dplyr::select(-labels),
#              aes(x=scaled_centered_weighted_CN, y=scaled_centered_DESeq, col=log(CN.value)))+
#   geom_point()
# ggplot(data = joint_counts_CN_subset %>% filter(nearestGeneCN.gene == 'MYC') %>%
#          dplyr::select(-labels),
#        aes(x=nearestGeneCN.value, y=counts.value))+
#   geom_point()+ggtitle('MYC CN and GE')


# ggplot(data = joint_counts_CN_subset,
#        aes(x=scaled_centered_weighted_CN, y=scaled_centered_DESeq, label=labels, col=labels))+
#   geom_point(alpha=1)+scale_x_continuous(trans="log2")+geom_smooth()+
#   # geom_label()+
#   facet_wrap(.~labels, scales="free", nrow=2)

##------------------------------------------------------------------------------------------------------------##


# ggplot(data = joint_counts_CN_subset %>% filter(CN.gene_name %in% top_cor),
#        aes(x=scaled_centered_weighted_CN, y=scaled_centered_DESeq, col=CN.gene_name))+
#   geom_point(alpha=1)+scale_x_continuous(trans="log2")+geom_smooth()+
#   facet_wrap(.~CN.gene_name, scales="free", nrow=2)
# ggplot(data = joint_counts_CN_subset %>% filter(CN.gene_name %in% top_r2),
#        aes(x=scaled_centered_weighted_CN, y=scaled_centered_DESeq, col=CN.gene_name))+
#   geom_point(alpha=1)+scale_x_continuous(trans="log2")+geom_smooth()+
#   facet_wrap(.~CN.gene_name, scales="free", nrow=2)

plot(corr_centered$`cor(scaled_centered_weighted_CN, scaled_centered_DESeq)`,
     r2_centered$r2)
plot(corr_centered$`cor(scaled_centered_weighted_CN, scaled_centered_DESeq)`,
     (var_CN$`var(CN.value)`))
plot(corr_centered$`cor(scaled_centered_weighted_CN, scaled_centered_DESeq)`,
     log(var_DESeq2$`var(DESeq.count)`))

# ggplot(df_gene_characteristics, aes(x=cor_normCNnormDESeq, y=var_DESeq,
#                                     label=as.character(ifelse(cor_normCNnormDESeq>0.8, Gene, NA))))+
#   geom_point()+geom_label_repel()+
#   scale_y_continuous(trans = "log2")

ggplot(data = joint_counts_CN_subset,
       aes(x=scaled_centered_weighted_CN, y=scaled_centered_DESeq, label=labels, col=labels))+
  geom_point(alpha=1)+scale_x_continuous(trans="log2")+geom_smooth()+
  # geom_label()+
  facet_wrap(.~labels, scales="free", nrow=2)
ggsave("../../figures_GRCh37/selected_genes_scaled_weightedCN_scaledDESeq_2.png", width = 7, height = 5)

# ggplot(data = joint_counts_CN_subset,
#        aes(x=scaled_centered_weighted_CN_mean, y=scaled_centered_DESeq, label=labels, col=labels))+
#   geom_point(alpha=1)+geom_smooth()+
#   # geom_label()+
#   facet_wrap(.~labels, scales="free", nrow=2)


# ggplot(data = joint_counts_CN %>% filter(counts.Var2 %in% subset_genes_of_interest),
#        aes(x=scaled_centered_weighted_CN, y=scaled_centered_DESeq, label=counts.Var2))+
#   geom_point(alpha=1)+scale_x_continuous(trans="log2")+geom_smooth()+
#   # geom_label()+
#   facet_wrap(.~counts.Var2, scales="free", nrow=2)

# ggplot(data = joint_counts_CN_subset %>% filter(nearestGeneCN.gene == 'PARP1') %>%
#          dplyr::select(-labels),
#        aes(x=nearestGeneCN.value, y=scaled_centered_DESeq))+
#   geom_point()
# ggplot(data = joint_counts_CN_subset %>% filter(nearestGeneCN.gene == 'PARP1'),
#        aes(x=scaled_centered_weighted_CN, y=scaled_centered_DESeq, label=labels, col=labels))+
#   geom_point(alpha=1)+scale_x_continuous(trans="log2")+geom_smooth()+
#   # geom_label()+
#   facet_wrap(.~labels, scales="free")


ggplot(data = joint_counts_CN_subset %>% filter(nearestGeneCN.gene %in% subset_genes_of_interest) %>%
         dplyr::select(-labels),
       aes(x=nearestGeneCN.value, y=counts.value, nearestGeneCN.gene=nearestGeneCN.gene, label=PDO))+
  facet_wrap(.~nearestGeneCN.gene, scales = "free")+
  geom_text_repel(size=2.5)+
  geom_point()+ggtitle('genes of interest CN and GE')
ggsave("../../figures_GRCh37/outliers_annotated_orgs.png", width = 7, height = 5)

ggplot(joint_counts_CN_subset, aes(x=log(CN.value)))+geom_density()+
  facet_wrap(.~factor(PDO, levels=logCN$PDO[order(logCN$mean)]), nrow=2)+
  geom_vline(xintercept = mean(log(joint_counts_CN_subset$CN.value)), col='red')
ggsave("../../figures_GRCh37/mean_CN_segment_per_org.pdf", width = 15, height = 5)

#-----------------------------------------------------------------------------------------#
## all genes
# ggplot(joint_counts_CN, aes(x=CN.value, y=counts.value))+geom_point()+
#   scale_y_continuous(trans='log2')+geom_smooth()

ggplot(joint_counts_CN, aes(x=CN.value, y=counts.value))+geom_point(alpha=0.2)+
  scale_y_continuous(trans='log2')+geom_smooth(col='red')+lims(x=c(0,20))+theme_bw()+
  labs(x='Weighted CN in gene', y='TPM (log2)')
ggsave("../../figures_GRCh37/scatterplot_CNweighted_TPM_all.pdf", width=8)

ggplot(joint_counts_CN, aes(x=CN.value, y=DESeq.count))+geom_point(alpha=0.2)+
  scale_y_continuous(trans='log2')+geom_smooth(col='red')+lims(x=c(0,20))+theme_bw()+
  labs(x='Weighted CN in gene', y='DEseq2 counts (log2)')
ggsave("../../figures_GRCh37/scatterplot_CNweighted_DEseq2counts_all.pdf", width=8)

## deseq counts and tmp seem to be quite similar
# ggplot(joint_counts_CN, aes(x=counts.value, y=DESeq.count))+geom_point(alpha=0.2)


ggplot(joint_counts_CN, aes(x=CN.value, y=DESeq.count, col=topDiff_norm_tissue))+geom_point(alpha=0.2)+
  scale_y_continuous(trans='log2')+geom_smooth(col='red')+lims(x=c(0,20))+theme_bw()+
  labs(x='Weighted CN in gene', y='DEseq2 counts (log2)')+theme(legend.position = "bottom")
ggsave("../../figures_GRCh37/scatterplot_CNweighted_DEseq2counts_all_colourDE.pdf", width=4, height=4)

ggplot(joint_counts_CN[joint_counts_CN$topDiff_norm_tissue == 'DE',], aes(x=CN.value, y=DESeq.count, col=topDiff_norm_tissue))+geom_point(alpha=0.2)+
  scale_y_continuous(trans='log2')+geom_smooth(col='red')+lims(x=c(0,20))+theme_bw()+
  labs(x='Weighted CN in gene', y='DEseq2 counts (log2)')+theme(legend.position = "bottom")
ggsave("../../figures_GRCh37/scatterplot_CNweighted_DEseq2counts_onlyDE.pdf", width=4, height=4)

ggplot(joint_counts_CN[joint_counts_CN$counts.value > 0,],
       aes(x=CN.value, y=counts.value))+geom_point(alpha=0.2)+
  scale_y_continuous(trans='log2')+geom_smooth(col='red')+lims(x=c(0,20))+theme_bw()+
  labs(x='Weighted CN in gene', y='TPM (log2)')+theme(legend.position = "bottom")
ggsave("../../figures_GRCh37/scatterplot_CNweighted_TPM_nonzero.pdf", width=8)

ggplot(joint_counts_CN_subset[joint_counts_CN_subset$counts.value > 0,],
       aes(x=nearestGeneCN.value, y=counts.value))+geom_point(alpha=0.2)+
  scale_y_continuous(trans='log2')+geom_smooth(col='red')+lims(x=c(0,20))+theme_bw()+
  labs(x='Weighted CN in gene', y='TPM (log2)')+theme(legend.position = "bottom")
ggsave("../../figures_GRCh37/scatterplot_CNclosest_TPM_nonzero.pdf", width=8)

#-----------------------------------------------------------------------------------------------#

## for each gene, find the rank of both ploidy and RNASeq
## group by gene

# ggplot(joint_counts_CN_normalised_highlyvar, aes(x=scaled_centered_weighted_CN,
#                                                  y=scaled_centered_DESeq))+
#   geom_point(alpha=0.9)+
#   geom_smooth()+
#   facet_wrap(.~factor(counts.Var1,
#     levels=c( "119058org", "119148org", "151723org", "151761org", "151773org", "23868org",  "32077org",
#               "50495org", "54059org",  "54276org", '118976org', '119127org', '54327org', '119178org','118947org' )),
#     nrow=2, scales = "free")
# ggplot(joint_counts_CN_normalised_highlyvar, aes(x=scaled_centered_weighted_CN,
#                                                  y=scaled_centered_DESeq))+
#   geom_point(alpha=0.9)+
#   geom_smooth()

# ggplot(joint_counts_CN_normalised_highlyvar, aes(x=scaled_neighbour_CN,
#                                                  y=scaled_centered_DESeq))+
#   geom_point(alpha=0.9)+
#   geom_smooth()+
#   facet_wrap(.~factor(counts.Var1,
#                       levels=c( "119058org", "119148org", "151723org", "151761org", "151773org", "23868org",  "32077org",
#                                 "50495org", "54059org",  "54276org", '118976org', '119127org', '54327org', '119178org','118947org' )),
#              nrow=2, scales = "free")

# ggplot(joint_counts_CN_normalised_highlyvar, aes(x=nearestGeneCN.value,
#                                                  y=counts.value))+
#   geom_point(alpha=0.9)+
#   geom_smooth()+
#   facet_wrap(.~factor(counts.Var1,
#                       levels=c( "119058org", "119148org", "151723org", "151761org", "151773org", "23868org",  "32077org",
#                                 "50495org", "54059org",  "54276org", '118976org', '119127org', '54327org', '119178org','118947org' )),
#              nrow=2, scales = "free")
# ggsave("../../figures_GRCh37/scatterplot_normCNweighted_normDESeq_perorg2.png", width=8, height=2.5)

#-----------------------------------------------------------------------------------------------#

plot(joint_counts_CN_ranked$rank_weighted_CN, joint_counts_CN_ranked$rank_DESeq)

# ggplot(joint_counts_CN_ranked, aes(rank_weighted_CN, rank_DESeq)) + geom_point()
# ggplot(joint_counts_CN_ranked, aes(x = rank_weighted_CN, y = rank_DESeq)) + geom_tile()

ggplot((melt(table(joint_counts_CN_ranked$rank_weighted_CN, joint_counts_CN_ranked$rank_DESeq)/nrow(joint_counts_CN_ranked))),
       aes(x=Var1, y=Var2, fill=value))+geom_raster()+scale_fill_jcolors_contin("pal2", bias = 1.75) + theme_bw()+
  labs(x='Organoid rank for CN (weighted)', y='Organoid rank for gene expression (DESeq2 counts)')
ggsave("../../figures_GRCh37/rankplot_CNweighted_deseq.pdf", width=8)
ggsave("../../figures_GRCh37/rankplot_CNweighted_deseq.png", width=8)

# ggplot((melt(table(joint_counts_CN_ranked_subsetorgs$rank_weighted_CN, joint_counts_CN_ranked_subsetorgs$rank_DESeq)/nrow(joint_counts_CN_ranked_subsetorgs))),
#        aes(x=Var1, y=Var2, fill=value))+geom_raster()+scale_fill_jcolors_contin("pal2", bias = 1.75) + theme_bw()+
#   labs(x='Organoid rank for CN (weighted)', y='Organoid rank for gene expression (DESeq2 counts)')

## is the rank of the excluded organoids any different from the rank of the included organoids?
grid.arrange(ggplot(joint_counts_CN_ranked, aes(x=(counts.Var1 %in% c('118976org', '119127org', '54327org', '119178org', '118947org')),
                                      y= rank_weighted_CN))+geom_violin(),
             ggplot(joint_counts_CN_ranked, aes(x=(counts.Var1 %in% c('118976org', '119127org', '54327org', '119178org', '118947org')),
                                   y= rank_DESeq))+geom_violin(), nrow=1)

df <- expand.grid(x = sample(x = 1:100, size = 50), y = sample(1:100, size = 50))
# default is compatible with geom_tile()
# ggplot(df, aes(x, y)) + geom_raster()

## blob, but clear trend
# ggplot(joint_counts_CN_normalised %>% filter(CN.value != 2), aes(x=scaled_centered_weighted_CN, y=scaled_centered_DESeq))+
#   geom_point(alpha=0.01)+geom_smooth()

# ggplot(joint_counts_CN_normalised %>% filter(CN.value != 2), aes(x=scaled_centered_weighted_CN, y=scaled_centered_DESeq))+
#   geom_point(alpha=0.01)+geom_smooth()+facet_wrap(.~(CN.value < 2), scales = "free")

ggplot(joint_counts_CN_normalised %>% filter(CN.value > 2), aes(x=scaled_centered_weighted_CN, y=scaled_centered_DESeq))+
  geom_point(alpha=0.01)+geom_smooth()
ggsave("../../figures_GRCh37/scatterplot_normCNweighted_normDESeq.pdf", width=8)

ggplot(joint_counts_CN_normalised_excludingnormal %>% filter(CN.value > 2), aes(x=scaled_centered_weighted_CN, y=scaled_centered_DESeq))+
  geom_point(alpha=0.01)+geom_smooth()+facet_wrap(.~PDO)
ggsave("../../figures_GRCh37/scatterplot_normCNweighted_normDESeq_perorg.pdf", width=8)

ggplot(joint_counts_CN_normalised_excludingnormal %>% filter(CN.value > 2), aes(x=scaled_centered_weighted_CN,
                                                                                y=scaled_centered_DESeq))+
  geom_point(alpha=0.01)+geom_smooth()+facet_wrap(.~factor(PDO,
   levels=renaming$PDO[match(c( "119058org", "119148org", "151723org", "151761org", "151773org", "23868org",  "32077org",  "50495org", 
                                "54059org",  "54276org", '118976org', '119127org', '54327org', '119178org', '118947org' ), renaming$ID)]), nrow=2)
ggsave("../../figures_GRCh37/scatterplot_normCNweighted_normDESeq_perorg2.png", width=8, height=2.5)

ggplot(joint_counts_CN_normalised_excludingnormal %>% filter(CN.value > 2) %>%
         filter(counts.Var1 %in% c( "119058org", "119148org", "151723org", "151761org", "151773org", "23868org",  "32077org",  "50495org", 
                         "54059org",  "54276org")),
                 aes(x=scaled_centered_weighted_CN, y=scaled_centered_DESeq, col=PDO))+
  geom_point(alpha=0.01)+geom_smooth()
ggsave("../../figures_GRCh37/scatterplot_normCNweighted_normDESeq_perorg2_subsetorgs.png", width=5, height=4)

ggplot(joint_counts_CN_normalised_excludingnormal %>% filter(CN.value > 2) ,
       aes(x=counts.Var1, y=scaled_centered_DESeq))+
  geom_boxplot(alpha=0.01)+
    facet_wrap(.~(counts.Var1 %in% c( "119058org", "119148org", "151723org", "151761org", "151773org", "23868org",  "32077org",  "50495org", 
                                     "54059org",  "54276org")), scales = "free_x")


## same when we include normal segments
# ggplot(joint_counts_CN_normalised, aes(x=scaled_centered_weighted_CN,
#                                                                                 y=scaled_centered_DESeq))+
#   geom_point(alpha=0.01)+geom_smooth()+facet_wrap(.~factor(counts.Var1,
#                                                            levels=c( "119058org", "119148org", "151723org", "151761org", "151773org", "23868org",  "32077org",  "50495org", 
#                                                                      "54059org",  "54276org", '118976org', '119127org', '54327org', '119178org', '118947org' )), nrow=2)
## same when we only include the very variable segments?
# ggplot(joint_counts_CN_normalised %>% filter(CN.value > 2), aes(x=scaled_centered_weighted_CN,
#                                                                 y=scaled_centered_DESeq))+
#   geom_point(alpha=0.01)+geom_smooth()+facet_wrap(.~factor(counts.Var1,
#                                                            levels=c( "119058org", "119148org", "151723org", "151761org", "151773org", "23868org",  "32077org",  "50495org", 
#                                                                      "54059org",  "54276org", '118976org', '119127org', '54327org', '119178org', '118947org' )), nrow=2)


## same but normalising only with nonnormal segments
ggplot(joint_counts_CN_normalised_excludingnormal, aes(x=scaled_centered_weighted_CN, y=scaled_centered_DESeq))+
  geom_point(alpha=0.002)+geom_smooth()
ggsave("../../figures_GRCh37/scatterplot_normCNweighted_normDESeq_nonormal.pdf", width=8)
ggsave("../../figures_GRCh37/scatterplot_normCNweighted_normDESeq_nonormal.png", width=4)

ggplot(joint_counts_CN_normalised_excludingnormal, aes(x=scaled_centered_weighted_CN, y=scaled_centered_TPM))+
  geom_point(alpha=0.002)+geom_smooth()
ggsave("../../figures_GRCh37/scatterplot_normCNweighted_normTPM_nonormal.png", width=4)

# ggplot(joint_counts_CN_normalised_excludingnormal, aes(x=scaled_centered_weighted_CN, y=scaled_centered_DESeq))+
#   geom_point(alpha=0.01)+geom_smooth()+facet_wrap(.~(CN.value < 2), scales = "free")

ggplot(joint_counts_CN_normalised_excludingnormal %>% filter(CN.value > 2), aes(x=scaled_centered_weighted_CN, y=scaled_centered_DESeq))+
  geom_point(alpha=0.01)+geom_smooth()
ggsave("../../figures_GRCh37/scatterplot_normCNweighted_normDESeq_nonormal_onlyamp.png", width=4)

## only nonnormal segments, and only with Stephane's genes
ggplot(joint_counts_CN_normalised_excludingnormal %>% filter(counts.Var2 %in% unlist(topvariable)),
       aes(x=scaled_centered_weighted_CN, y=scaled_centered_DESeq, col=CN.value < 2))+
  geom_point(alpha=0.1)+geom_smooth()+facet_wrap(.~(CN.value < 2), scales = "free")
ggsave("../../figures_GRCh37/scatterplot_normCNweighted_normDESeq_nonormal_highlyvarregions.pdf", width=8)

ggplot(joint_counts_CN_normalised %>% filter(CN.value != 2), aes(x=scaled_centered_weighted_CN, y=DESeq.count))+
  geom_point(alpha=0.2)+geom_smooth()+scale_y_continuous(trans = "log2")

## not great
# ggplot(joint_counts_CN_normalised %>% filter(CN.value != 2), aes(x=CN.value, y=scaled_centered_DESeq))+
#   geom_point(alpha=0.2)+geom_smooth()+lims(x=c(0,10))

## good trend
# ggplot(joint_counts_CN_normalised  %>% filter(CN.value != 2), aes(x=scaled_centered_weighted_CN, y=scaled_centered_TPM))+
#   geom_point(alpha=0.2)+geom_smooth()

## not great, but good trend
# ggplot(joint_counts_CN_normalised_excludingnormal, aes(x=scaled_centered_weighted_CN, y=scaled_centered_TPM))+
#   geom_point(alpha=0.2)+geom_smooth()

# save.image("../output/analyse_joing_counts_CN.RData")


##------------------------------------------------------------------------------------------------------------##
plot(df_gene_characteristics$averageCN, df_gene_characteristics$averageDESeq2)

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

ntop_coord = 100
plot(sort(mean_coord_CN_DESeq2$mean_coord_CN_DESeq2))
plot(sort(mean_coord_CN_DESeq2$mean_coord_CN_DESeq2, decreasing = T)[1:ntop_coord])

# ggplot(droplevels(mean_coord_CN_DESeq2[1:ntop_coord,]), aes(label=Gene, x=rank,
#                                  y=mean_coord_CN_DESeq2))+
#   geom_point()+geom_label_repel()

ggplot((mean_coord_CN_DESeq2), aes(label=ifelse(test = Gene %in% subset_genes_of_interest, yes = Gene, no=NA),
                                             x=rank, y=mean_coord_CN_DESeq2))+
  geom_point()+theme_classic()+
  theme(axis.text.x=element_blank())+geom_label_repel()
ggsave("../../figures_GRCh37/coordinated_CN_DESeq_goi.pdf", width = 8, height = 8)

dim(mean_coord_CN_DESeq2)
dim(df_gene_characteristics) 

plot(mean_coord_CN_DESeq2$mean_coord_CN_DESeq2, df_gene_characteristics$r2_normCNnormDESeq[match(mean_coord_CN_DESeq2$Gene, df_gene_characteristics$Gene)])
plot(df_gene_characteristics$cor_normCNnormDESeq, df_gene_characteristics$r2_normCNnormDESeq)
plot(df_gene_characteristics$cor_normCNnormDESeq, log(df_gene_characteristics$averageCN))

ggplot(df_gene_characteristics,
  aes(x=cor_normCNnormDESeq, y=averageCN,
                                    # label=ifelse( (log(averageCN) > 1) & (cor_normCNnormDESeq > 0.8), yes = Gene, no = NA  )))+
      label=ifelse( (log(averageCN) > 1.5) & (cor_normCNnormDESeq > 0.5), yes = Gene, no = NA  )))+
  geom_point()+scale_y_continuous(trans = "log2")+geom_label_repel(size=3, max.overlaps = 30)+theme_bw()
ggsave("../../figures_GRCh37/coordinated_CN_DESeq_and_CNvalue.pdf", width = 8, height = 8)

## example of a 'good' result
mean_coord_CN_DESeq2 %>% filter(Gene == 'TOP1')
df_gene_characteristics %>% filter(Gene == 'TOP1') %>% dplyr::select(normalsamples.averageDESeq2)
joint_counts_CN%>% filter(counts.Var2 == 'TOP1') %>% dplyr::select(DESeq.count)
  
# ggplot(joint_counts_CN_subset %>% filter(counts.Var2 == 'FBL'),
#        aes(x=scaled_centered_weighted_CN, y=scaled_centered_DESeq, label=counts.Var1), alpha=0.2)+
#   geom_smooth()+geom_label_repel()+facet_wrap(.~labels)+
#   guides(col=FALSE)+theme_bw()

# ggplot(joint_counts_CN_subset %>% filter(counts.Var2 == 'MYC'),
#        aes(x=scaled_centered_weighted_CN, y=scaled_centered_DESeq, label=counts.Var1), alpha=0.2)+
#   geom_smooth()+geom_label_repel()+facet_wrap(.~labels)+geom_point()+
#   guides(col=FALSE)+theme_bw()
# ggplot(joint_counts_CN_subset %>% filter(counts.Var2 == 'MYC'),
#        aes(x=CN.value, y=DESeq.count, label=counts.Var1), alpha=0.2)+
#   geom_smooth()+geom_label_repel()+facet_wrap(.~labels)+geom_point()+
#   guides(col=FALSE)+theme_bw()
ggplot(joint_counts_CN_subset %>% filter(counts.Var2 == 'KRT6B'),
       aes(x=CN.value, y=DESeq.count, label=PDO), alpha=0.2)+
  geom_smooth()+geom_label_repel()+geom_point()+
  guides(col=FALSE)+theme_bw()+ggtitle('KRT6Bs')
ggsave("../../figures_GRCh37/example_KRT6B.pdf", width = 4, height = 4)
ggplot(joint_counts_CN_subset %>% filter(counts.Var2 == 'AKT2'),
       aes(x=CN.value, y=DESeq.count, label=PDO), alpha=0.2)+
  geom_smooth()+geom_label_repel()+geom_point()+
  guides(col=FALSE)+theme_bw()+ggtitle('AKT2')
ggsave("../../figures_GRCh37/example_AKT2.pdf", width = 4, height = 4)


plot(df_gene_characteristics$slope_centered, df_gene_characteristics$r2_normCNnormDESeq)
head(df_gene_characteristics$slope_centered)
head(df_gene_characteristics$r2_normCNnormDESeq)


ggplot(droplevels(df_average_bottomCN),
       aes(x=Gene, y=average_comparison_CN_DESeq, label=as.character(label)))+
  geom_point()+#geom_label_repel(max.overlaps = 30, segment.size = 0.01)+
  # geom_text_repel(force = .01, direction = "y", nudge_x = 0, nudge_y = .1)+
  geom_text_repel(direction = "y", nudge_x=-1, size=3)+
  theme(axis.text.x=element_blank(), axis.line.x.bottom =element_blank(), 
        panel.grid.minor = element_line(size = 0.1, colour = "black"))+
  geom_hline(yintercept = 0.5, lty='dashed')
ggsave("../../figures_GRCh37/average_bottomCN_DESeq.pdf", width = 6, height = 6)

ggplot(droplevels(df_average_bottomCN),
       aes(x=factor(Gene, levels=Gene[order(df_average_bottomCN$average_comparison_CN_TPM)]),
           y=average_comparison_CN_TPM, label=as.character(label)))+
  geom_point()+#geom_label_repel(max.overlaps = 30, segment.size = 0.01)+
  # geom_text_repel(force = .01, direction = "y", nudge_x = 0, nudge_y = .1)+
  geom_text_repel(direction = "y", nudge_x=-1, size=3)+
  theme(axis.text.x=element_blank(), axis.line.x.bottom =element_blank(), 
        panel.grid.minor = element_line(size = 0.1, colour = "black"))+
  geom_hline(yintercept = 0.5, lty='dashed')+labs(x='Gene')
ggsave("../../figures_GRCh37/average_bottomCN_DESeq_TPM.pdf", width = 6, height = 6)

ggplot(droplevels(df_average_bottomCN),
       aes(x=Gene, y=average_comparison_CN_TPM, label=as.character(label)))+
  geom_point()+#geom_label_repel(max.overlaps = 30, segment.size = 0.01)+
  # geom_text_repel(force = .01, direction = "y", nudge_x = 0, nudge_y = .1)+
  geom_text_repel(direction = "y", nudge_x=-1, size=3)+
  theme(axis.text.x=element_blank(), axis.line.x.bottom =element_blank(), 
        panel.grid.minor = element_line(size = 0.1, colour = "black"))+
  geom_hline(yintercept = 0.5, lty='dashed')
ggsave("../../figures_GRCh37/average_bottomCN_TPM.pdf", width = 6, height = 6)

# ggplot(droplevels(df_gene_characteristics),
#        aes(x=r2_normCNnormDESeq, col=df_average_bottomCN.average_comparison_CN_DESeq,
#            group=df_average_bottomCN.average_comparison_CN_DESeq))+
#   geom_density()+#geom_label_repel(max.overlaps = 30, segment.size = 0.01)+
#   # geom_text_repel(force = .01, direction = "y", nudge_x = 0, nudge_y = .1)+
#   theme(axis.text.x=element_blank(), axis.line.x.bottom =element_blank(), 
#         panel.grid.minor = element_line(size = 0.1, colour = "black"))+
#   scale_x_continuous(trans = "log2")+
#   geom_hline(yintercept = 0.5, lty='dashed')


ggplot(joint_counts_CN %>% filter(CN.gene_name %in% topvariable$V1),
       aes(x=factor(CN.gene_name, levels=median_CN$CN.gene_name[order(median_CN$CN.value)]),
           y=CN.value, group=CN.gene_name))+geom_violin()
ggsave("../../figures_GRCh37/CN_violinplots.pdf", width = 9, height = 4)

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
ggsave("../../figures_GRCh37/CN_violinplots_goi.pdf", width = 6, height = 6)

ggplot(joint_counts_CN %>% filter(counts.Var2 %in% subset_genes_of_interest),
       aes(x=factor(counts.Var2, levels=median_CN_goi$counts.Var2[order(median_CN_goi$CN.value)]),
           y=scaled_centered_DESeq, col=(scaled_centered_weighted_CN)))+geom_violin()+
  geom_point(size=3)+
  geom_line(aes(group=counts.Var1), col='black', alpha=0.1)+
  theme(legend.position = "bottom")+labs(x='Gene', y='Scaled and centered DESeq2 counts')+
  scale_color_jcolors_contin("pal2", reverse = TRUE, bias = 2.25)+
  theme(axis.text.x=element_text(angle = 90, hjust = 0))
ggsave("../../figures_GRCh37/CN_violinplots_goi_2.pdf", width = 6, height = 6)

# ggplot(joint_counts_CN %>% filter(counts.Var2 %in% c('AKT1', 'AKT3')),
#        aes(x=factor(counts.Var2, levels=median_CN_goi$counts.Var2[order(median_CN_goi$CN.value)]),
#            y=scaled_centered_DESeq, col=(scaled_centered_weighted_CN)))+geom_violin()+
#   geom_point(size=3)+
#   geom_line(aes(group=counts.Var1), col='black', alpha=0.1)+
#   theme(legend.position = "bottom")+labs(x='Gene, sorted by median weighted CN', y='Weighted CN (log scale)')+
#   scale_color_jcolors_contin("pal2", reverse = TRUE, bias = 2.25)+
#   theme(axis.text.x=element_text(angle = 90, hjust = 0))


plot(df_gene_characteristics$df_average_bottomCN.average_comparison_CN_DESeq,
df_gene_characteristics$r2_normCNnormDESeq)

ggplot(df_gene_characteristics, aes(x=df_average_bottomCN.average_comparison_CN_DESeq, y=r2_normCNnormDESeq,
                                    label=ifelse(Gene %in% subset_genes_of_interest,
                                                 yes = Gene, no = NA )))+
  geom_boxplot(aes(group=df_average_bottomCN.average_comparison_CN_DESeq),
               width=0.5/length(unique(df_gene_characteristics$df_average_bottomCN.average_comparison_CN_DESeq)),
               col='blue')+
  geom_point()+geom_label_repel()
ggsave("../../figures_GRCh37/r2normCNnormDESeq_vs_averagebottomCN.pdf", width = 6, height = 6)

ggplot(df_gene_characteristics[,!duplicated(colnames(df_gene_characteristics))],
       aes(x=df_average_bottomCN.average_comparison_CN_TPM, y=r2_normCNnormDESeq,
           label=ifelse(Gene %in% subset_genes_of_interest,
                        yes = Gene, no = NA )))+
  geom_boxplot(aes(group=df_average_bottomCN.average_comparison_CN_TPM),
               width=0.5/length(unique(df_gene_characteristics$df_average_bottomCN.average_comparison_CN_TPM)),
               col='blue')+
  geom_point()+geom_label_repel()
ggsave("../../figures_GRCh37/r2normCNnormDESeq_vs_averagebottomCN_TPM.pdf", width = 6, height = 6)

ggplot(df_gene_characteristics, aes(x=df_average_bottomCN.average_comparison_CN_DESeq, y=r2_normCNnormDESeq,
                                    label=ifelse( (df_average_bottomCN.average_comparison_CN_DESeq >= .8) &
                                                    (r2_normCNnormDESeq > 0.5),
                                                 yes = Gene, no = NA )))+
  geom_boxplot(aes(group=df_average_bottomCN.average_comparison_CN_DESeq),
               width=0.5/length(unique(df_gene_characteristics$df_average_bottomCN.average_comparison_CN_DESeq)),
               col='blue')+
  geom_point()+geom_label_repel()
ggsave("../../figures_GRCh37/r2normCNnormDESeq_vs_averagebottomCN_best.pdf", width = 6, height = 6)

ggplot((melt(df_gene_characteristics, measure.vars = c('df_average_bottomCN.average_comparison_CN_DESeq',
                                                                   'df_average_bottomCN.average_comparison_CN_DESeq2',
                                                                   'df_average_bottomCN.average_comparison_CN_DESeq4'))),
       aes(y=r2_normCNnormDESeq, x=value, col=log2(averageDESeq2)))+
  geom_point()+
  geom_line(aes(group=r2_raw_DESeq2.Gene), alpha=0.1)+
  scale_y_continuous(trans = "log2")+
  geom_smooth()
ggsave("../../figures_GRCh37/r2normCNnormDESeq_vs_averagebottomCN_2_3_4_bottom.pdf", width = 6, height = 6)

df_gene_characteristics$chrom = ag$seq_name[match(df_gene_characteristics$Gene, ag$gene_name)]

ggplot((melt(df_gene_characteristics, measure.vars = c('df_average_bottomCN.average_comparison_CN_DESeq',
                                                       'df_average_bottomCN.average_comparison_CN_DESeq2',
                                                       'df_average_bottomCN.average_comparison_CN_DESeq4'))),
       aes(y=r2_normCNnormDESeq, x=value, col=log2(averageDESeq2)))+
  geom_point()+
  geom_line(aes(group=r2_raw_DESeq2.Gene), alpha=0.1)+
  scale_y_continuous(trans = "log2")+
  facet_wrap(.~chrom)
ggsave("../../figures_GRCh37/r2normCNnormDESeq_vs_averagebottomCN_2_3_4_bottom_bychrom.pdf", width = 14, height = 14)

median_r2_by_chrom <- df_gene_characteristics[,!duplicated(t(df_gene_characteristics))] %>% group_by(chrom) %>% summarize(median_r2=median(r2_normCNnormDESeq, na.rm = T))

ggplot(df_gene_characteristics[,!duplicated(t(df_gene_characteristics))],
       aes(x=factor(chrom,levels=median_r2_by_chrom$chrom[order(median_r2_by_chrom$median_r2)]),
           y=r2_normCNnormDESeq, col=chrom))+geom_jitter(alpha=0.1)+
  geom_boxplot()+
  # geom_violin()+
  labs(x='Chromosome (sorted)', y='R^2 between GE and CN in organoids')
  # scale_y_continuous(trans = "log2")
ggsave("../../figures_GRCh37/r2normCNnormDESeq_sorted_by_chrom.pdf", width = 6.5, height = 5)

grid.arrange(ggplot(df_gene_characteristics, aes(x=df_average_bottomCN.average_comparison_CN_DESeq,
                                    y=df_average_bottomCN.average_comparison_CN_DESeq2))+
  geom_tile(alpha=0.01),
ggplot(df_gene_characteristics, aes(x=df_average_bottomCN.average_comparison_CN_DESeq,
                                    y=df_average_bottomCN.average_comparison_CN_DESeq4))+
  geom_tile(alpha=0.01))

ggplot(df_gene_characteristics)+
  geom_bin2d(aes(x=df_average_bottomCN.average_comparison_CN_DESeq,
                   y=df_average_bottomCN.average_comparison_CN_DESeq2))+
  scale_fill_continuous(type = "viridis")


ggplot(df_gene_characteristics, aes(x=df_average_bottomCN.average_comparison_CN_DESeq, y=r2_raw_DESeq2.r2,
                                    label=ifelse( (df_average_bottomCN.average_comparison_CN_DESeq >= .8) &
                                                    (r2_raw_DESeq2.r2 > 0.5),
                                                  yes = Gene, no = NA )))+
  geom_boxplot(aes(group=df_average_bottomCN.average_comparison_CN_DESeq),
               width=0.5/length(unique(df_gene_characteristics$df_average_bottomCN.average_comparison_CN_DESeq)),
               col='blue')+
  geom_point()+geom_label_repel()
ggsave("../../figures_GRCh37/r2CNDESeq_vs_averagebottomCN.pdf", width = 6, height = 6)

plot(joint_counts_CN_subset$scaled_centered_weighted_CN_mean,
     log(joint_counts_CN_subset$scaled_centered_weighted_CN))
# ggplot(joint_counts_CN_subset, aes(x=scaled_centered_weighted_CN_mean, y=scaled_centered_weighted_CN, col=scaled_centered_DESeq))+
#   geom_point()+scale_y_continuous(trans = "log")

# ggplot(joint_counts_CN_subset, aes(x=scaled_centered_weighted_CN_mean, y=scaled_centered_DESeq))+
#   geom_point()+geom_smooth()#+scale_y_continuous(trans = "log")

ggplot(df_gene_characteristics,
       aes(x=cor_normCNnormDESeq, y=averageCN,
           # label=ifelse( (log(averageCN) > 1) & (cor_normCNnormDESeq > 0.8), yes = Gene, no = NA  )))+
           label=ifelse( (log(averageCN) > 1.5) & (cor_normCNnormDESeq > 0.5), yes = Gene, no = NA  )))+
  geom_point(alpha=0.2)+
  geom_density_2d()+
  scale_y_continuous(trans = "log2")+
  geom_label_repel(max.overlaps = 30)+theme_bw()+
  labs(x='Correlation between CN and GE', y='Average copy number')
df_gene_characteristics[df_gene_characteristics$Gene == 'MYC','cor_normCNnormDESeq']
MYCdf_gene_characteristics <- df_gene_characteristics[df_gene_characteristics$Gene == 'MYC',]
MYCdf_gene_characteristics
MYCcorr_centered =  joint_counts_CN[,!duplicated(colnames(joint_counts_CN))] %>%
  filter(Gene == 'MYC') %>% summarise(cor(scaled_centered_weighted_CN, scaled_centered_DESeq))

plot(remove_na(joint_counts_CN$scaled_centered_weighted_CN[joint_counts_CN$Gene == 'MYC',]),
     remove_na(joint_counts_CN$scaled_centered_DESeq[joint_counts_CN$Gene == 'MYC',]))
cor(remove_na(joint_counts_CN$scaled_centered_weighted_CN[joint_counts_CN$Gene == 'MYC',]),
     remove_na(joint_counts_CN$scaled_centered_DESeq[joint_counts_CN$Gene == 'MYC',]))

TP53 <- joint_counts_CN %>% filter(CN.gene_name == 'TP53')
ggplot(TP53, aes(x=log2(CN.value), y=DESeq.count, label=PDO))+geom_point()+
  geom_label_repel(alpha=0.2)+
  ggtitle('TP53')+labs(x='Absolute CN (log2)', y='DESeq2 counts')
ggsave("../../figures_GRCh37/TP53.pdf", height = 4, width = 4)

ZWINT <- joint_counts_CN %>% filter(CN.gene_name == 'ZWINT')
ggplot(ZWINT, aes(x=log2(CN.value), y=DESeq.count, label=PDO))+geom_point()+
  geom_label_repel(alpha=0.2)+
  ggtitle('ZWINT')+labs(x='Absolute CN (log2)', y='DESeq2 counts')
ggsave("../../figures_GRCh37/ZWINT.pdf", height = 4, width = 4)

ERBB2 <- joint_counts_CN %>% filter(CN.gene_name == 'ERBB2')
ggplot(ERBB2, aes(x=log2(CN.value), y=DESeq.count, label=PDO))+geom_point()+
  geom_label_repel(alpha=0.2)+
  ggtitle('ERBB2')+labs(x='Absolute CN (log2)', y='DESeq2 counts')
ggsave("../../figures_GRCh37/ERBB2.pdf", height = 4, width = 4)

ggplot( joint_counts_CN %>% filter(CN.gene_name == 'AKT2'), aes(x=log2(CN.value), y=DESeq.count, label=PDO))+geom_point()+
  geom_label_repel(alpha=0.2)+
  ggtitle('AKT2')+labs(x='Absolute CN (log2)', y='DESeq2 counts')
ggsave("../../figures_GRCh37/AKT2.pdf", height = 4, width = 4)

ggplot( joint_counts_CN %>% filter(CN.gene_name == 'WEE1'), aes(x=log2(CN.value), y=DESeq.count, label=PDO))+geom_point()+
  geom_label_repel(alpha=0.2)+
  ggtitle('WEE1')+labs(x='Absolute CN (log2)', y='DESeq2 counts')
ggsave("../../figures_GRCh37/WEE1.pdf", height = 4, width = 4)

ggplot( joint_counts_CN %>% filter(CN.gene_name == 'ATR'), aes(x=log2(CN.value), y=DESeq.count, label=PDO))+geom_point()+
  geom_label_repel(alpha=0.2)+
  ggtitle('ATR')+labs(x='Absolute CN (log2)', y='DESeq2 counts')
ggsave("../../figures_GRCh37/ATR.pdf", height = 4, width = 4)

ggplot( joint_counts_CN %>% filter(CN.gene_name == 'CCNE1'), aes(x=log2(CN.value), y=DESeq.count, label=PDO))+geom_point()+
  geom_label_repel(alpha=0.2)+
  ggtitle('CCNE1')+labs(x='Absolute CN (log2)', y='DESeq2 counts')
ggsave("../../figures_GRCh37/CCNE1.pdf", height = 4, width = 4)

ccne1 <- joint_counts_CN %>% filter(CN.gene_name == 'CCNE1')
brcamut <- c(PDO4='BRCAmut', PDO5='BRCAmut', PDO6='BRCAmut', PDO7='BRCAmut', PDO8='BRCAmut',
             PDO10='BRCAmut', PDO13='BRCAmut', PDO17='BRCAmut')
ccne1$brcamut = brcamut[match(ccne1$PDO, names(brcamut))]
ccne1$brcamut[is.na(ccne1$brcamut)] <- 'WT'
ggplot( ccne1, aes(x=log2(CN.value), shape=brcamut, y=DESeq.count, col=PDO, label=PDO))+geom_point()+
  # geom_label_repel(alpha=0.5, col='black')+
  scale_colour_manual(values = col_vector)+
  ggtitle('CCNE1')+labs(x='Absolute CN (log2)', y='DESeq2 counts')+theme_bw()
ggsave("../../figures_GRCh37/CCNE1_brca_mut.pdf", height = 4, width = 4.5)


ggplot( joint_counts_CN %>% filter(CN.gene_name == 'MYC'), aes(x=log2(CN.value), y=DESeq.count, label=PDO))+geom_point()+
  geom_label_repel(alpha=0.2)+
  ggtitle('MYC')+labs(x='Absolute CN (log2)', y='DESeq2 counts')
ggsave("../../figures_GRCh37/MYC.pdf", height = 4, width = 4)

ggplot( joint_counts_CN %>% filter(CN.gene_name == 'MAD1L1'), aes(x=log2(CN.value), y=DESeq.count, label=PDO))+geom_point()+
  geom_label_repel(alpha=0.2)+
  ggtitle('MAD1L1')+labs(x='Absolute CN (log2)', y='DESeq2 counts')
ggsave("../../figures_GRCh37/MAD1L1.pdf", height = 4, width = 4)

ggplot( joint_counts_CN %>% filter(CN.gene_name == 'MNT'), aes(x=log2(CN.value), y=DESeq.count, label=PDO))+geom_point()+
  geom_label_repel(alpha=0.2)+
  ggtitle('MNT')+labs(x='Absolute CN (log2)', y='DESeq2 counts')
ggsave("../../figures_GRCh37/MNT.pdf", height = 4, width = 4)

ggplot( joint_counts_CN %>% filter(CN.gene_name == 'PARP1'), aes(x=log2(CN.value), y=DESeq.count, label=PDO))+geom_point()+
  geom_label_repel(alpha=0.2)+
  ggtitle('PARP1')+labs(x='Absolute CN (log2)', y='DESeq2 counts')
ggsave("../../figures_GRCh37/PARP1.pdf", height = 4, width = 4)

ggplot( joint_counts_CN %>% filter(CN.gene_name == 'PARP2'), aes(x=log2(CN.value), y=DESeq.count, label=PDO))+geom_point()+
  geom_label_repel(alpha=0.2)+
  ggtitle('PARP2')+labs(x='Absolute CN (log2)', y='DESeq2 counts')
ggsave("../../figures_GRCh37/PARP2.pdf", height = 4, width = 4)

ggplot( joint_counts_CN %>% filter(CN.gene_name == 'PTEN'), aes(x=log2(CN.value), y=DESeq.count, label=PDO))+geom_point()+
  geom_label_repel(alpha=0.2)+
  ggtitle('PTEN')+labs(x='Absolute CN (log2)', y='DESeq2 counts')
ggsave("../../figures_GRCh37/PTEN.pdf", height = 4, width = 4)

ggplot( joint_counts_CN %>% filter(CN.gene_name == 'ATM'), aes(x=log2(CN.value), y=DESeq.count, label=PDO))+geom_point()+
  geom_label_repel(alpha=0.2)+
  ggtitle('ATM')+labs(x='Absolute CN (log2)', y='DESeq2 counts')
ggsave("../../figures_GRCh37/ATM.pdf", height = 4, width = 4)


ggplot( joint_counts_CN %>% filter(CN.gene_name == 'MTOR'), aes(x=log2(CN.value), y=DESeq.count, label=PDO))+geom_point()+
  geom_label_repel(alpha=0.2)+
  ggtitle('MTOR')+labs(x='Absolute CN (log2)', y='DESeq2 counts')
ggsave("../../figures_GRCh37/MTOR.pdf", height = 4, width = 4)

ggplot( joint_counts_CN %>% filter(CN.gene_name == 'PIK3CA'), aes(x=log2(CN.value), y=DESeq.count, label=PDO))+geom_point()+
  geom_label_repel(alpha=0.2)+
  ggtitle('PIK3CA')+labs(x='Absolute CN (log2)', y='DESeq2 counts')
ggsave("../../figures_GRCh37/PIK3CA.pdf", height = 4, width = 4)

ggplot( joint_counts_CN %>% filter(CN.gene_name == 'PIK3CB'), aes(x=log2(CN.value), y=DESeq.count, label=PDO))+geom_point()+
  geom_label_repel(alpha=0.2)+
  ggtitle('PIK3CB')+labs(x='Absolute CN (log2)', y='DESeq2 counts')
ggsave("../../figures_GRCh37/PIK3CB.pdf", height = 4, width = 4)

ggplot( joint_counts_CN %>% filter(CN.gene_name == 'PIK3CD'), aes(x=log2(CN.value), y=DESeq.count, label=PDO))+geom_point()+
  geom_label_repel(alpha=0.2)+
  ggtitle('PIK3CD')+labs(x='Absolute CN (log2)', y='DESeq2 counts')
ggsave("../../figures_GRCh37/PIK3CD.pdf", height = 4, width = 4)

exposures <- readRDS("../../../copy_number_analysis_organoids/robjects/exposures.RDS")
pik3 <- joint_counts_CN %>% filter(grepl('PIK3', CN.gene_name))
pik3$s4 <- exposures[match(pik3$PDO, rownames(exposures)),'s4']
ggplot( pik3, aes(x=log2(CN.value), y=DESeq.count, label=PDO, col=s4))+geom_point()+
  geom_label_repel(alpha=0.2, size=2)+
  ggtitle('PI3K+')+labs(x='Absolute CN (log2)', y='DESeq2 counts')+facet_wrap(.~CN.gene_name, scales = "free")
ggsave("../../figures_GRCh37/PI3K.pdf", height = 8, width = 8)

ggplot( joint_counts_CN %>% filter(CN.gene_name == 'CDK12'), aes(x=log2(CN.value), y=DESeq.count, label=PDO))+geom_point()+
  geom_label_repel(alpha=0.2)+
  ggtitle('CDK12')+labs(x='Absolute CN (log2)', y='DESeq2 counts')
ggsave("../../figures_GRCh37/CDK12.pdf", height = 4, width = 4)


xls_out <- df_gene_characteristics[,c('df_average_bottomCN.average_comparison_CN_DESeq', 'r2_normCNnormDESeq', 'Gene', 'averageCN', 'averageDESeq2')]
xls_out <- xls_out[!is.na(xls_out$df_average_bottomCN.average_comparison_CN_DESeq),]
xls_out <- xls_out[order(xls_out$df_average_bottomCN.average_comparison_CN_DESeq),]
write.csv(xls_out, file = "../../../RNASeq_and_CN/20191218_ViasM_BJ_orgaBrs/output/output_GRCh37/df_fig3cd.csv")

xls_out2 <- xls_out[order((xls_out$df_average_bottomCN.average_comparison_CN_DESeq + xls_out$r2_normCNnormDESeq), decreasing = T),]
write.csv(xls_out2, file = "../../../RNASeq_and_CN/20191218_ViasM_BJ_orgaBrs/output/output_GRCh37/df_fig3cd_sort2.csv")
# xls_out2 <- read.csv(file = "../../../RNASeq_and_CN/20191218_ViasM_BJ_orgaBrs/output/output_GRCh37/df_fig3cd_sort2.csv", header = T)

xls_out2$chrom = ag$seq_name[match(xls_out2$Gene, ag$gene_name)]
df_gene_characteristics$chrom = ag$seq_name[match(df_gene_characteristics$CN.gene_name, ag$gene_name)]

xls_out2_highlyvar <- xls_out2 %>% filter(Gene %in% (unique(joint_counts_CN_subset$CN.gene_name)))

# ggplot(xls_out2, aes(x=log(averageCN), y=log(averageDESeq2),
#                      col=df_average_bottomCN.average_comparison_CN_DESeq))+
#          geom_point()

# ggplot(xls_out2[!(df_gene_characteristics$var_CN == 0),], aes(x=log(averageCN), y=log(averageDESeq2),
#                      col=df_average_bottomCN.average_comparison_CN_DESeq))+
#   geom_point()


# ggplot(xls_out2, aes(y=df_average_bottomCN.average_comparison_CN_DESeq, x=r2_normCNnormDESeq,
#                      col=averageCN))+
#   geom_point()

ggplot(xls_out2, aes(y=df_average_bottomCN.average_comparison_CN_DESeq, x=r2_normCNnormDESeq,
                     col=averageCN))+
  geom_point()+facet_wrap(.~factor(chrom, levels=1:22))
ggsave("../../figures_GRCh37/r2normCNnormDESeq_vs_averagebottomCN_perchrom.pdf", width = 8)

ggplot(xls_out2_highlyvar, aes(y=df_average_bottomCN.average_comparison_CN_DESeq, x=r2_normCNnormDESeq,
                     col=averageCN))+
  geom_point()+facet_wrap(.~factor(chrom, levels=1:22))
ggsave("../../figures_GRCh37/r2normCNnormDESeq_vs_averagebottomCN_perchrom_subsethighlyvar.pdf", width = 8)

ggplot(df_gene_characteristics, aes(x=df_average_bottomCN.average_comparison_CN_DESeq, y=r2_normCNnormDESeq,
                                    label=ifelse( (df_average_bottomCN.average_comparison_CN_DESeq >= .8) &
                                                    (r2_normCNnormDESeq > 0.5),
                                                  yes = Gene, no = NA )))+
  geom_point()+geom_label_repel(max.overlaps=10, size=3)+theme_bw()+labs(x='Averaged higher CN and higher GE', y=TeX('R^2 between CN and GE'))+coord_flip()+
  lims(x=c(0.5, 1), y=c(0.5, 1))+facet_wrap(.~factor(chrom, levels=c(1:22)), drop = T)
ggsave("../../figures_GRCh37/fig3d_chromosomes.pdf", width = 10, height = 10)

# df_gene_characteristics$label_top20 = ifelse( df_gene_characteristics$Gene %in% xls_out2_highlyvar$Gene[1:20] | ((df_gene_characteristics$df_average_bottomCN.average_comparison_CN_DESeq >= .8) &
#                                                                                            (df_gene_characteristics$r2_normCNnormDESeq > 0.5)),
#                                               yes = df_gene_characteristics$Gene, no = NA )
df_gene_characteristics$label_top20 = ifelse( df_gene_characteristics$Gene %in% xls_out2_highlyvar$Gene[1:20] | !(df_gene_characteristics$chrom %in% c(8,12)),
                                              yes = df_gene_characteristics$Gene, no = NA )
df_gene_characteristics$label_top20_col = ifelse( df_gene_characteristics$Gene %in% xls_out2_highlyvar$Gene[1:20], T, F)

min_x <- 0.5
min_y <- 0.5
ggplot(df_gene_characteristics %>% filter(Gene %in% (unique(joint_counts_CN_subset$CN.gene_name))) %>%
         filter(df_average_bottomCN.average_comparison_CN_DESeq > min_x,
                r2_normCNnormDESeq > min_y),
       aes(x=df_average_bottomCN.average_comparison_CN_DESeq,
                                    y=r2_normCNnormDESeq))+
  geom_point()+
  geom_label_repel(aes(label=label_top20, col=label_top20_col), size=3)+
  theme_bw()+labs(y='Averaged higher CN and higher GE', x=TeX('R^2 between CN and GE'))+coord_flip()+
  lims(x=c(min_x, 1), y=c(min_y, 1))+facet_wrap(.~drop(factor(chrom, levels=c(1:22))), drop = T)+
  guides(col=FALSE)
ggsave("../../figures_GRCh37/fig3d_chromosomes_highlyvarsubset.pdf", width = 12, height = 10)


# pdf("../../figures_GRCh37/scatterplot_top20genes.pdf", width = 5, height = 5)
for(top20gene in xls_out2_highlyvar$Gene[1:20]){
  pdf(paste0("../../figures_GRCh37/scatterplot_top20genes/scatterplot_top20genes_highlyvarsubset_",
             top20gene, ".pdf"), width = 5, height = 5)
  print(ggplot(joint_counts_CN %>% filter(CN.gene_name == top20gene),
         aes(x=CN.value, y=DESeq.count, label=PDO))+geom_point()+geom_label_repel(alpha=0.2)+
    ggtitle(top20gene))
  dev.off()
}

ggplot(joint_counts_CN %>% filter(CN.gene_name %in% xls_out2_highlyvar$Gene[1:20]),
       aes(x=CN.value, y=DESeq.count, label=PDO))+geom_point()+geom_label_repel(alpha=0.2)+
  facet_wrap(.~factor(Gene, levels=xls_out2_highlyvar$Gene[1:20]), scales = "free")+
  theme(strip.text = element_text(size=14))
ggsave("../../figures_GRCh37/scatterplot_alltop20genes.pdf", width = 13, height = 13)

#setwd("~/Documents/PhD/other_repos/Vias_Brenton/RNASeq_and_CN/20191218_ViasM_BJ_orgaBrs/Scripts/")
# old <- readRDS("../output/output_GRCh38_incorrect/joint_counts_CN_TPM_20210301210511.RDS")
# new <- readRDS("../output/output_GRCh37/joint_counts_CN_TPM_20210414211439.RDS")
# 
# old_bot3 <- readRDS("../output/output_GRCh38_incorrect/fig3_df_average_bottomCN.RDS")
# dim(old_bot3)
# 
# require(scales)
# mtch <- match(paste0(new$counts.Var1, new$counts.Var2), paste0(old$counts.Var1, old$counts.Var2))
# 
# pdf("~/Desktop/comparison_CNvalues_genes_intactgenenames_loglog.pdf")
# plot(log(old$CN.value[mtch]),
#      log(new$CN.value), pch=19, col=alpha('black', 0.05))
# dev.off()

## find most highly amplified genes
## use the original joint_counts_CN, as we want to keep all samples
joint_counts_CN0$chrom = ag$seq_name[match(joint_counts_CN0$CN.gene_name, ag$gene_name)]

most_highly_amplified = joint_counts_CN0 %>% group_by(CN.gene_name) %>% summarise(mean_CN_org= mean(CN.value))
most_highly_amplified <- most_highly_amplified[order(most_highly_amplified$mean_CN_org, decreasing = T),]

most_highly_amplified_genes_list <- most_highly_amplified$CN.gene_name[1:20]
joint_counts_CN0$PDO = renaming$PDO[match(joint_counts_CN0$CN.L1, renaming$ID)] ## because PDO was created with RNASeq data, so if there's none we didn't have the PDO name
joint_counts_CN0_MHA = joint_counts_CN0 %>% filter(CN.gene_name %in% most_highly_amplified_genes_list)
give_order_colsums <- function(f){
  unlist(f[order(f[,2], decreasing = T),1])
}

ggplot(joint_counts_CN0_MHA,
       aes(y=log2(CN.value), x=factor(CN.gene_name, levels=rev(most_highly_amplified_genes_list)),
           fill=factor(chrom, levels=sort(as.numeric(unique(chrom))))))+
  geom_bar(stat = "identity")+coord_flip()+
  facet_wrap(.~factor(PDO2, levels = give_order_colsums(joint_counts_CN0_MHA %>% group_by(PDO2) %>% summarise(sum(CN.value)))), nrow=1)+
  geom_hline(yintercept = log(2), col='blue')+
  labs(x='Gene name', y='CN (log2)', fill='Chromosome')
ggsave("../../figures_GRCh37/most_highly_amplified.pdf", height = 7, width = 10)

# xxx <- joint_counts_CN %>% filter(CN.gene_name %in% most_highly_amplified_genes_list)
# table(xxx$PDO2)
# table(joint_counts_CN$CN.L1)

# > table(xxx$PDO2)
# 
# PDO1 PDO10 PDO11 PDO12 PDO13 PDO14 PDO15 PDO16 PDO17 PDO18  PDO2  PDO3  PDO4  PDO5  PDO6  PDO7  PDO8 
# 20    20    20    20    20     1    20     1    20     1    20    20    20    20    20    20    20 
# PDO9 
# 20 
# > table(joint_counts_CN$CN.L1)
# 
# 118947org 118976org 119025org 119058org 119127org 119148org 119178org 151723org 151761org 151773org 
# 19430     19430      1028     19430     19430     19430     19430     19430     19430     19430 
# 23868org  32070org  32077org  50495org  54059org  54276org  54288org  54327org 
# 19430      1028     19430     19430     19430     19430      1028     19430 

table(joint_counts_CN$CN.L1)
table(joint_counts_CN$counts.Var1)

View(joint_counts_CN %>% filter(CN.L1 == '119025org'))

library(RColorBrewer)
n <- length(unique(joint_counts_CN$PDO))
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

ggplot(joint_counts_CN0 %>% filter(CN.gene_name %in% subset_genes_of_interest),
       aes(x=CN.gene_name, y=CN.value))+geom_violin()+
  geom_jitter(aes(col=PDO))+
  # geom_label(aes(label=PDO))+
  facet_wrap(.~CN.gene_name, scales = "free")+
  scale_colour_manual(values = col_vector)+
  geom_hline(yintercept = 2, lty='dashed', col='blue')+theme_bw()
ggsave("../../figures_GRCh37/genes_of_interest_CN.pdf", width = 7)

ggplot(joint_counts_CN0 %>% filter(CN.gene_name %in% subset_genes_of_interest),
       aes(x=CN.gene_name, y=log2(CN.value)))+geom_violin()+
  geom_jitter(aes(col=PDO))+
  # geom_label(aes(label=PDO))+
  facet_wrap(.~CN.gene_name, scales = "free")+
  scale_colour_manual(values = col_vector)+
  geom_hline(yintercept = log2(2), lty='dashed', col='blue')+theme_bw()
  # scale_y_continuous(trans = "log2")
ggsave("../../figures_GRCh37/genes_of_interest_CN_log2.pdf", width = 7)

GOI_summary <- joint_counts_CN0 %>% filter(CN.gene_name %in% subset_genes_of_interest) %>%
  dplyr::group_by(CN.gene_name) %>% dplyr::summarise(median_cn=median(log2(CN.value)))
ggplot(joint_counts_CN0 %>% filter(CN.gene_name %in% subset_genes_of_interest),
       aes(x=factor(CN.gene_name, levels=as.character(GOI_summary$CN.gene_name[order(GOI_summary$median_cn)])),
           y=log2(CN.value)))+geom_boxplot()+
  geom_jitter(aes(col=PDO), alpha=0.2)+
  scale_colour_manual(values = col_vector)+
  geom_hline(yintercept = log2(2), lty='dashed', col='blue')+theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  labs(x='Gene', y='Absolute CN (log2)', col='')
ggsave("../../figures_GRCh37/genes_of_interest_CN_log2_boxplot.pdf", width = 5.25, height = 5)

ggplot(joint_counts_CN0 %>% filter(CN.gene_name %in% subset_genes_of_interest),
       aes(x=factor(CN.gene_name, levels=as.character(GOI_summary$CN.gene_name[order(GOI_summary$median_cn)])),
           y=log2(CN.value)))+geom_violin()+
  geom_jitter(alpha=0.2)+
  scale_colour_manual(values = col_vector)+
  geom_hline(yintercept = log2(2), lty='dashed', col='blue')+theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  labs(x='Gene', y='Absolute CN (log2)')
ggsave("../../figures_GRCh37/genes_of_interest_CN_log2_violin.pdf", width = 5.25, height = 5)


ggplot(joint_counts_CN0 %>% filter(CN.gene_name %in% subset_genes_of_interest),
       aes(x=factor(CN.gene_name, levels=as.character(GOI_summary$CN.gene_name[order(GOI_summary$median_cn)])),
           y=log2(CN.value), fill=PDO))+geom_bar(stat='identity')+
  scale_fill_manual(values = col_vector)+
  facet_wrap(.~PDO, ncol=4)+
  geom_hline(yintercept = log2(2), lty='dashed', col='blue')+theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  labs(x='Gene', y='Absolute CN (log2)', col='')+guides(fill=F)
ggsave("../../figures_GRCh37/genes_of_interest_CN_log2_barplot.pdf", width = 15, height = 11)


joint_counts_CN0$PDO <- factor(joint_counts_CN0$PDO, levels=gtools::mixedsort(unique(joint_counts_CN0$PDO)))
CN_genes_of_interest_mat <- (dcast(joint_counts_CN0 %>% filter(CN.gene_name %in% subset_genes_of_interest) %>%
             dplyr::select(c('counts.Var2', 'CN.value', 'PDO')), CN.gene_name~PDO, value.var='CN.value'))
write.table(CN_genes_of_interest_mat,
              file = "../../tables_GRCh37/table_CN_genes_of_interest.csv", sep = ",",
              row.names = F, col.names = T)

rownames(CN_genes_of_interest_mat) <- CN_genes_of_interest_mat$CN.gene_name
remove_cols_with_no_variance <- function(i)  t(i[,-1])[,!(apply(t(i[,-1]), 2, var) == 0)]

pheatmap::pheatmap(cor(remove_cols_with_no_variance((CN_genes_of_interest_mat[,-1]))),
                   cluster_cols = T, cluster_rows = T)
cor_CN_genes_of_interest <- cor(remove_cols_with_no_variance((CN_genes_of_interest_mat[,-1])))
library(grid)

pdf("../../figures_GRCh37/correlation_CN_genes_of_interest.pdf", width = 8, height = 7.5)
ComplexHeatmap::Heatmap(cor_CN_genes_of_interest,
                        cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
                          grid.text(round(cor_CN_genes_of_interest[i, j], 1), x, y)
                        })
dev.off()

# mat_breaks <- c( 0:6, seq(8, 14, by=2)) - 0.01 ## adding a small number because otherwise the binning is done wrong
# breaks <- mat_breaks
# colours_vec <- c("#2670af", "#2670af", "#8daecf", "#eaeaea", "#ffd2b6", "#f5b788", "#f39d5f", "#f17e37", "#ee1200", "#b60000", "#7a0e04", "#401004", "#000000")
# ggplot(joint_counts_CN0 %>% filter(CN.gene_name %in% subset_genes_of_interest),
#        aes(x=counts.Var2, y=PDO, fill=CN.value))+geom_tile()+
#   # scale_fill_jcolors_contin("pal3", reverse = T, bias = 1.0)
#   # scale_fill_gradient2(low = 1,mid = 2, high = 3, midpoint = 2)
#   # scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = 1, breaks = mat_breaks)
#   # scale_fill_gradientn(colours = c("red","white","blue"),
#   #                        breaks = breaks, labels = format(breaks))
#   scale_fill_steps2(low = 0, mid = 2, high = 100, midpoint = 2, )
# ggsave("../../figures_GRCh37/genes_of_interest_CN_tile.pdf", width = 7)

ggplot(joint_counts_CN0 %>% filter(CN.gene_name %in% subset_genes_of_interest),
       aes(x=counts.Var2, y=PDO,
           # fill=paste0((CN.value==0), (CN.value>0 & CN.value<=1),
           #                                  (CN.value>1 & CN.value<=2), CN.value>2), size=10*log(CN.value),
           fill=cut(round(CN.value, 1), breaks = c(0, 1.9, 2.1, 3, 5, 50, 150)),
           label=round(CN.value, 1)))+
  geom_tile(alpha=0.2)+
  scale_fill_manual(values=c('blue', 'white', 'yellow', 'orange',  'red', 'purple'))+theme_bw()+
  geom_text(col='black', size=3)+guides(fill=F, size=F)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+labs(x='Genes', y='PDO')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave("../../figures_GRCh37/genes_of_interest_CN_tile.pdf", width = 7, height = 5)

ggplot(joint_counts_CN0 %>% filter(CN.gene_name %in% subset_genes_of_interest),
       aes(x=CN.gene_name, y=CN.value))+geom_violin()+
  geom_jitter(aes(col=PDO))+
  # geom_label(aes(label=PDO))+
  facet_wrap(.~CN.gene_name, scales = "free_x")+
  scale_colour_manual(values = col_vector)+
  scale_y_continuous(trans = "log2")+
  geom_hline(yintercept = 2, lty='dashed', col='blue')+theme_bw()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
ggsave("../../figures_GRCh37/genes_of_interest_CN_v2.pdf", width = 7, height = 7)

ggplot(joint_counts_CN0 %>% filter(CN.gene_name %in% subset_genes_of_interest),
       aes(x=CN.gene_name, y=CN.value))+
  geom_jitter(aes(col=PDO))+
  geom_text_repel(aes(label=gsub('PDO', '', PDO), col=PDO),
                  # position = position_jitterdodge(),
                  segment.color = NA, max.overlaps = 30, size=2)+
  # geom_label(aes(label=PDO))+
  facet_wrap(.~CN.gene_name, scales = "free_x")+
  scale_colour_manual(values = col_vector)+
  scale_y_continuous(trans = "log2")+
  geom_hline(yintercept = 2, lty='dashed', col='blue')+theme_bw()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
ggsave("../../figures_GRCh37/genes_of_interest_CN_v3.pdf", width = 7, height = 7)

#'First, a group of genes including PARP2, TERT ATR and AKT1 show an inverse genomic 
#'expression pattern than KRAS and TOP1. Secondly, and orthogonal to the previous genes,
#' PIK3CA and ATM are expressed with inverse relation to PARP1, CCNE1 and PTEN. 
#' Organoids derived from the same patient do not cluster together in regard to the
#'transcriptome. 

genes_rnaseq <- joint_counts_CN %>% filter(CN.gene_name %in% c('PARP2', 'TERT', 'ATR', 'AKT1', 'KRAS', 'TOP1',
                                               'PIK3CA', 'ATM', 'PARP1', 'CCNE1', 'PTEN'))

ggplot(genes_rnaseq, aes(x=Gene, y=DESeq.count))+geom_bar(stat = "identity")+facet_wrap(.~PDO)

genes_rnaseq2 <- dcast(genes_rnaseq[,c('Gene', 'DESeq.count', 'PDO')], PDO~Gene, value.var = 'DESeq.count')
rownames(genes_rnaseq2) <- genes_rnaseq2[,1]; genes_rnaseq2 <- genes_rnaseq2[,-1]
pairs(genes_rnaseq2, pch=18)

ggplot(joint_counts_CN %>% filter(CN.gene_name == "LRP1B"), aes(x=PDO, y=DESeq.count, group=1))+geom_line()+
  scale_y_continuous(trans = "log2")

#-----------------------------------------------------------------------------------------------#

## Link RNASeq to drug response

joint_counts_CN$pairs_pdo[joint_counts_CN$PDO %in% c('PDO3', 'PDO9')] = 'PDO3-PDO9'
joint_counts_CN$pairs_pdo[joint_counts_CN$PDO %in% c('PDO5', 'PDO6')] = 'PDO5-PDO6'
joint_counts_CN$pairs_pdo[joint_counts_CN$PDO %in% c('PDO7', 'PDO8')] = 'PDO7-PDO8'

gene_subsets_iterate_list = list(AZD8186_Vistusertib_Capivasertib_AZD8835=c('PI3', 'AKT1', 'MTOR'),
                                 Olaparib=c('PARP1','PARP2'),
                                 AZD0156='ATM',
                                 Ceralasertib='ATR',
                                 Adavosertib='WEE1',
                                 Doxorubicin='TOP2A',
                                 Paclitaxel=NA,
                                 Gemcitabine=c('RRM1', 'RRM2B', 'RRM2B', 'DCK', 'CDA', 'Slc29a1', 'Slc29a2', 'Slc28a1', 'Slc28a3'),
                                 Oxaliplatin=NA,
                                 APR246='TP53')

pdf("../../figures_GRCh37/geneexpression_and_drugsusceptibility.pdf", width = 8, height = 4)
for(gene_subsets_iterate_it in 1:length(gene_subsets_iterate_list)){
  gene_subsets_iterate = gene_subsets_iterate_list[[gene_subsets_iterate_it]]
  if(is.na(gene_subsets_iterate)){
    print(ggplot()+ggtitle(names(gene_subsets_iterate_list)[gene_subsets_iterate_it]))
  }else{
    print(ggplot(joint_counts_CN %>% filter(CN.gene_name %in% gene_subsets_iterate) %>%
             filter(PDO %in% c('PDO3', 'PDO9', 'PDO5', 'PDO6', 'PDO7', 'PDO8')),
           aes(x=interaction(PDO, CN.gene_name), y=DESeq.count, group=CN.gene_name, col=CN.gene_name))+
      geom_point()+geom_line()+facet_wrap(.~pairs_pdo, scales = "free_x")+
      theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))+
      ggtitle(names(gene_subsets_iterate_list)[gene_subsets_iterate_it])+
      labs(x=''))
  }
}
dev.off()

joint_counts_CN_dcast <- dcast(joint_counts_CN[!duplicated(joint_counts_CN),] %>%
                                 dplyr::select(PDO, DESeq.count, CN.gene_name) %>%
                                 filter(!is.na(DESeq.count) & !is.na(CN.gene_name)) %>%
                                 distinct(),
      value.var = "DESeq.count", formula = CN.gene_name~PDO)

joint_counts_CN_dcast_TPM <- dcast(joint_counts_CN[!duplicated(joint_counts_CN),] %>%
                                 dplyr::select(PDO, counts.value, CN.gene_name) %>%
                                 filter(!is.na(counts.value) & !is.na(CN.gene_name)) %>%
                                 distinct(),
                               value.var = "counts.value", formula = CN.gene_name~PDO)

rownames(joint_counts_CN_dcast) = joint_counts_CN_dcast[,1]; joint_counts_CN_dcast[,1] <- NULL
rownames(joint_counts_CN_dcast_TPM) = joint_counts_CN_dcast_TPM[,1]; joint_counts_CN_dcast_TPM[,1] <- NULL

plot(joint_counts_CN_dcast$PDO3, joint_counts_CN_dcast$PDO9)
plot(joint_counts_CN_dcast$PDO5, joint_counts_CN_dcast$PDO6)
plot(joint_counts_CN_dcast$PDO7, joint_counts_CN_dcast$PDO8)
plot((joint_counts_CN_dcast$PDO6), (joint_counts_CN_dcast$PDO10))
# pairs(joint_counts_CN_dcast)

cor_rnaseq <- cor(joint_counts_CN_dcast %>% filter(rownames(joint_counts_CN_dcast) %in% coding_genes$external_gene_name))
cor_rnaseq_TPM <- cor(joint_counts_CN_dcast_TPM)
pdf("../../figures_GRCh37/heatmap_correlation_rnaseq.pdf")
print(pheatmap::pheatmap(cor_rnaseq))
dev.off()
## makes you think that the values might have been reshuffled e.g. PDO10 (119178org) and PDO6 seem to be very similar

# pheatmap::pheatmap(cor_rnaseq_TPM)

plot(joint_counts_CN$DESeq.count, joint_counts_CN$counts.value)

dcast_2 <- dcast(deseq_counts, Gene~DESEq.sample, value.var = "DESeq.count", fun.aggregate = mean)
plot(dcast_2$`119148org`, dcast_2$`119178org`)

deseq_counts0 <- read.table("../../../RNASeq_DE_resistant_sensitive/files/counts_norm.csv",
                                           sep=',', header = T)
renaming[renaming$PDO %in% c('PDO6', 'PDO10'),]
plot(deseq_counts0$JBLAB.19902, deseq_counts0$JBLAB.19907)

## working on the hypothesis that JBLAB.19907 is actually PDO5, and not PDO10:
JBLAB.19907

## get the gene expression from PDO10 and the copy number from PDO5, to see if there is a good correlation

pdo_current_name <- 'PDO10'
pdo_of_interest <- 'PDO8'
pdo5_CN <- joint_counts_CN[joint_counts_CN$PDO == pdo_of_interest,c('CN.gene_name', 'scaled_centered_weighted_CN')] %>% filter(!is.na(CN.gene_name))
possiblepdo5_GE <- joint_counts_CN[joint_counts_CN$PDO == pdo_current_name,c('CN.gene_name', 'scaled_centered_DESeq')] %>% filter(!is.na(CN.gene_name))
currentpdo5_GE <- joint_counts_CN[joint_counts_CN$PDO == pdo_of_interest,c('CN.gene_name', 'scaled_centered_DESeq')] %>% filter(!is.na(CN.gene_name))
grid.arrange(ggplot(cbind.data.frame(cn=pdo5_CN$scaled_centered_weighted_CN,
                        ge=unlist(possiblepdo5_GE[match(unlist(pdo5_CN[,'CN.gene_name']),
      unlist(possiblepdo5_GE[,'CN.gene_name'])),'scaled_centered_DESeq'])), aes(x=cn, y=ge))+geom_point()+geom_smooth()+
  ggtitle(paste0('Suggested ', pdo_of_interest, '(current ', pdo_current_name, ')')),
ggplot(cbind.data.frame(cn=pdo5_CN$scaled_centered_weighted_CN,
                        ge=unlist(currentpdo5_GE[match(unlist(pdo5_CN[,'CN.gene_name']),
                                unlist(currentpdo5_GE[,'CN.gene_name'])),'scaled_centered_DESeq'])),
       aes(x=cn, y=ge))+geom_point()+geom_smooth()+
  ggtitle(paste0('Current ', pdo_of_interest)), nrow=1)

##' for each RNASeq vector, and for each CN vector, match the names of any two organoids and compute the
##' correlation of GE of organoid 1 with the CN of organoid 2

orgs <- unique(joint_counts_CN0$PDO)[!is.na(unique(joint_counts_CN0$PDO))]
joint_counts_CN0 = cbind.data.frame(joint_counts_CN0,
                                    deseq_counts[match(paste0(joint_counts_CN0$counts.Var1, joint_counts_CN0$CN.gene_name), paste0(deseq_counts$DESEq.sample, deseq_counts$Gene)),])
joint_counts_CN0 = joint_counts_CN0 %>%
  group_by(CN.gene_name) %>%
  mutate(scaled_centered_weighted_CN = exp(scale(log(CN.value))),
         scaled_centered_weighted_CN_mean = (scale((CN.value))),
         scaled_centered_DESeq = scale(DESeq.count),
         scaled_centered_TPM = scale(counts.value)
  )

cors_GE_CN <- outer(orgs, orgs, Vectorize(function(pdo_A, pdo_B){
  cat(pdo_A, '\t')
  cat(pdo_B, '\n')
  .x1 <- joint_counts_CN0 %>% filter(PDO == pdo_A) %>% filter(!is.na(scaled_centered_weighted_CN))
  .x2 <- joint_counts_CN0 %>% filter(PDO == pdo_B)
  .x2 <- .x2[match(paste0(.x1$CN.gene_name), paste0(.x2$CN.gene_name)),] %>% filter(!is.na(scaled_centered_DESeq))
  .x1 <- .x1[match(paste0(.x2$CN.gene_name), paste0(.x1$CN.gene_name)),]
  cor(.x1$scaled_centered_weighted_CN, .x2$scaled_centered_DESeq)
}))
rownames(cors_GE_CN) <- colnames(cors_GE_CN) <- orgs

pdf("../../figures_GRCh37/heatmap_correlation_rnaseq_CN_GE.pdf")
print(pheatmap::pheatmap(cors_GE_CN, cluster_rows = F, cluster_cols = F))
dev.off()


# joint_counts_CN$

  ggplot(joint_counts_CN_normalised_highlyvar, aes(x=nearestGeneCN.value,
                                                   y=counts.counts.value))+
  geom_point(alpha=0.9)+
  geom_smooth()+
  facet_wrap(.~factor(counts.Var1,
                      levels=c( "119058org", "119148org", "151723org", "151761org", "151773org", "23868org",  "32077org",
                                "50495org", "54059org",  "54276org", '118976org', '119127org', '54327org', '119178org','118947org' )),
             nrow=2, scales = "free")
  
# xx <- joint_counts_CN[!duplicated(joint_counts_CN),] %>%
#   dplyr::select(PDO, DESeq.count, CN.gene_name) %>%
#   filter(!is.na(DESeq.count) & !is.na(CN.gene_name)) %>% distinct()
# sum(duplicated(xx))
# sum(duplicated(joint_counts_CN))
# sum(duplicated(joint_counts_CN[!duplicated(joint_counts_CN),]))
# frequency_genes <- joint_counts_CN %>% filter(!duplicated(joint_counts_CN))%>%
#   dplyr::select(PDO, DESeq.count, CN.gene_name) %>%
#   ungroup() %>% filter(!is.na(DESeq.count)) %>% dplyr::select(CN.gene_name) %>% table() %>% sort(decreasing = T)
# 
# names(frequency_genes[frequency_genes>15])
# NPIPA7
# joint_counts_CN %>% filter(CN.gene_name == 'NPIPA7', PDO=='PDO3')

tacna <- read.table("~/Desktop/Degree of TACNA50/TCGA_degree_of_TACNA50.txt", sep = "\t",
                    header=T, comment.char="#",
                    na.strings=".", stringsAsFactors=FALSE,
                    quote="", fill=FALSE)
colnames(tacna)

ggplot((df_gene_characteristics[,!duplicated(t(df_gene_characteristics))] %>%
               filter(df_average_bottomCN.Gene %in% tacna$SYMBOL[tacna$Degree.of.TACNA > .9])),
       aes(y=r2_normCNnormDESeq, x=df_average_bottomCN.average_comparison_CN_DESeq,
           col=log2(averageDESeq2)))+
  geom_point()+
  geom_line(aes(group=r2_raw_DESeq2.Gene), alpha=0.1)+
  scale_y_continuous(trans = "log2")+
  geom_smooth()+geom_density_2d()
ggsave("../../figures_GRCh37/r2normCNnormDESeq_vs_averagebottomCN_highTACNA.pdf", width = 6, height = 6)


ggplot((df_gene_characteristics[,!duplicated(t(df_gene_characteristics))]),
       aes(y=r2_normCNnormDESeq, x=df_average_bottomCN.average_comparison_CN_DESeq,
           col=log2(averageDESeq2)))+
  geom_point(alpha=0.2)+
  geom_line(aes(group=r2_raw_DESeq2.Gene), alpha=0.1)+
  scale_y_continuous(trans = "log2")+
  geom_smooth()+geom_density_2d()
ggsave("../../figures_GRCh37/r2normCNnormDESeq_vs_averagebottomCN_2.pdf", width = 6, height = 6)


ggplot(joint_counts_CN %>% filter(CN.value > 2, !is.na(PDO)), aes(x=scaled_centered_weighted_CN,
                                                     y=scaled_centered_DESeq))+
  geom_smooth()+geom_density_2d()+facet_wrap(.~PDO)
ggsave("../../figures_GRCh37/scatterplot_normCNweighted_normDESeq_contour.pdf", width = 6, height = 6)


ggplot(joint_counts_CN %>% filter(CN.value > 2, !is.na(PDO),
                                  CN.gene_name %in% tacna$SYMBOL[tacna$Degree.of.TACNA > .9]), aes(x=scaled_centered_weighted_CN,
                                                     y=scaled_centered_DESeq))+
  geom_smooth()+geom_density_2d()+facet_wrap(.~PDO)
ggsave("../../figures_GRCh37/scatterplot_normCNweighted_normDESeq_contour_highTACNA.pdf", width = 6, height = 6)

ggplot(joint_counts_CN %>%filter(CN.gene_name == 'TTK'), aes(x=CN.L1, y=counts.value))+geom_point()

auc <- read.table("../../../survival_analysis/data/20200419-AUC-organoids_PFI.csv", sep = ',', header = T)
auc <- read.table("../../../survival_analysis/data/20200419-AUC-organoids.csv", sep = ',', header = T)

MAD1 <- joint_counts_CN %>% filter(CN.gene_name == 'MAD1L1')
MAD1_and_drugs <- cbind.data.frame(MAD1, auc[match(MAD1$PDO, renaming$PDO[match( auc$id, gsub("org", "", renaming$ID),)]),])[,c('DESeq.count', "auc_ll5.AZD8186",
                                                          "auc_ll5.AZD2014", "auc_ll5.AZD5363",
                                                          "auc_ll5.AZD6738",     "auc_ll5.AZD1775",
                                                          "auc_ll5.Paclitaxel",  "auc_ll5.Oxaliplatin",
                                                          "auc_ll5.Doxorubicin", "auc_ll5.Gemcitabine",
                                                          "auc_ll5.APR.246")]
pairs(MAD1_and_drugs)


deseq_counts_match_with_renormalised <- deseq_counts0[,match(colnames(deseq_counts_renormalised0), colnames(deseq_counts0))]
rownames(deseq_counts_match_with_renormalised) <- rownames(deseq_counts_renormalised0)
deseq_counts_match_with_renormalised <- deseq_counts_match_with_renormalised[match(rownames(deseq_counts_renormalised0), rownames(deseq_counts_match_with_renormalised)),]
dim(deseq_counts_match_with_renormalised)
dim(deseq_counts_renormalised0)
plot(log(unlist(as.vector(deseq_counts_match_with_renormalised))), log(unlist(as.vector(deseq_counts_renormalised0))))
abline(coef=c(0,1), col='blue')


df_pdo3_pdo9 <- data.frame(PDO3=(deseq_counts0[,gsub('-', '.', renaming$sampleNameRNAseq[renaming$PDO == 'PDO3'])]),
                           PDO9=(deseq_counts0[,gsub('-', '.', renaming$sampleNameRNAseq[renaming$PDO == 'PDO9'])]))
rownames(df_pdo3_pdo9) = deseq_counts0$X
df_pdo3_pdo9$off_diag_1 <- log(df_pdo3_pdo9[,1])-log(df_pdo3_pdo9[,2]) > 2
df_pdo3_pdo9$off_diag_2 <- log(df_pdo3_pdo9[,1])-log(df_pdo3_pdo9[,2]) < -2
df_pdo3_pdo9$col = interaction(df_pdo3_pdo9$off_diag_1, df_pdo3_pdo9$off_diag_2)

data(c2BroadSets) ## from GSVAdata`

df_pdo3_pdo9$HRD = ifelse(rownames(df_pdo3_pdo9) %in% t2g_38$ensembl_gene_id[match(geneIds(c2BroadSets['KEGG_HOMOLOGOUS_RECOMBINATION']@.Data[[1]]), t2g_38$entrezgene_id)],
                          yes = t2g_38$external_gene_name[match(rownames(df_pdo3_pdo9), t2g_38$ensembl_gene_id)], no = NA)

remove_na <- function(i){
  .x <- i[!is.na(i)]
  names(.x) <- names(i)[!is.na(i)]
  .x
}

ggplot(df_pdo3_pdo9, aes(x=PDO3, y=PDO9, col=col, label=HRD))+geom_point()+
  scale_x_continuous(trans = 'log2')+scale_y_continuous(trans = 'log2')+
  geom_abline(slope = 1, intercept = 0, lty='dashed')+
  scale_fill_manual(values = c('blue', 'red', 'green'))+
  geom_label_repel()
ggsave("../../figures_GRCh37/PDO3_PDO9_HRD.pdf")


library(goseq)
goResults_upper <- goseq(nullp((remove_na(setNames(object = as.numeric(df_pdo3_pdo9$off_diag_1), nm = rownames(df_pdo3_pdo9)))),
                               "hg19", "ensGene"), "hg19","ensGene", test.cats=c("GO:BP"))
goResults_lower <- goseq(nullp((remove_na(setNames(object = as.numeric(df_pdo3_pdo9$off_diag_2), nm = rownames(df_pdo3_pdo9)))),
                               "hg19", "ensGene"), "hg19","ensGene", test.cats=c("GO:BP"))

pdf("../../figures_GRCh37/PDO3_PDO9_GO.pdf", width = 10, height = 4)
grid.arrange(goResults_upper %>% 
  top_n(10, wt=-over_represented_pvalue) %>% 
  mutate(hitsPerc=numDEInCat*100/numInCat) %>% 
  ggplot(aes(x=hitsPerc, 
             y=term, 
             colour=over_represented_pvalue, 
             size=numDEInCat)) +
  geom_point() +
  expand_limits(x=0) +
  labs(x="Hits (%)", y="GO term", colour="p value", size="Count")+
    ggtitle('Upper diagonal genes (enriched in PDO9)'),
  goResults_lower %>% 
  top_n(10, wt=-over_represented_pvalue) %>% 
  mutate(hitsPerc=numDEInCat*100/numInCat) %>% 
  ggplot(aes(x=hitsPerc, 
             y=term, 
             colour=over_represented_pvalue, 
             size=numDEInCat)) +
  geom_point() +
  expand_limits(x=0) +
  labs(x="Hits (%)", y="GO term", colour="p value", size="Count")+
    ggtitle('Lower diagonal genes (enriched in PDO3)'),
  nrow=1)
dev.off()


genes_pdo3_pdo9_upper <- deseq_counts0[match(remove_na(rownames(df_pdo3_pdo9)[df_pdo3_pdo9$off_diag_1]), deseq_counts0$X),gsub('-', '.', renaming$sampleNameRNAseq[match(remove_na(unique(joint_counts_CN$PDO)), renaming$PDO)])]
genes_pdo3_pdo9_lower <- deseq_counts0[match(remove_na(rownames(df_pdo3_pdo9)[df_pdo3_pdo9$off_diag_2]), deseq_counts0$X),gsub('-', '.', renaming$sampleNameRNAseq[match(remove_na(unique(joint_counts_CN$PDO)), renaming$PDO)])]
colnames(genes_pdo3_pdo9_upper) <- renaming$PDO[match(gsub("[.]", "-", colnames(genes_pdo3_pdo9_upper)), renaming$sampleNameRNAseq)]
colnames(genes_pdo3_pdo9_lower) <- renaming$PDO[match(gsub("[.]", "-", colnames(genes_pdo3_pdo9_lower)), renaming$sampleNameRNAseq)]

pheatmap(as.matrix(log(genes_pdo3_pdo9_upper+0.01)))
pheatmap(as.matrix(log(genes_pdo3_pdo9_lower+0.01)))

# library(topGO)
# sampleGOdata <- new("topGOdata",
#                     description = "Simple session", ontology = "BP",
#                     allGenes = setNames(rep(1, sum(df_pdo3_pdo9$off_diag_1 | df_pdo3_pdo9$off_diag_2, na.rm = T)),
#                                         remove_na(c(rownames(df_pdo3_pdo9)[df_pdo3_pdo9$off_diag_1],
#                                            rownames(df_pdo3_pdo9)[df_pdo3_pdo9$off_diag_2]))),
#                     geneSel = function(i){i %in% remove_na(rownames(df_pdo3_pdo9)[df_pdo3_pdo9$off_diag_1])},
#                     nodeSize = 10,
#                     annot = annFUN.db, affyLib = affyLib)

give_goterm_list_of_genes <- function(genes_in_list_ens, genes_not_in_list_ens){
  require(topGO)
  .xentrez_list = t2g_38[match(genes_in_list_ens, t2g_38$ensembl_gene_id),'entrezgene_id']
  .xentrez_list = .xentrez_list[!is.na(.xentrez_list)]
  .xentrez = t2g_38[match(genes_not_in_list_ens, t2g_38$ensembl_gene_id),'entrezgene_id']
  .xentrez = .xentrez[!is.na(.xentrez)]
  
  .xres_gosimenrich = GOSim::GOenrichment(genesOfInterest = as.character(.xentrez_list), allgenes = as.character(.xentrez))
  .xres_gosimenrich
}

gosim_upper <- give_goterm_list_of_genes(genes_in_list_ens = remove_na(c(rownames(df_pdo3_pdo9)[df_pdo3_pdo9$off_diag_1])),
                                         genes_not_in_list_ens = remove_na(c(rownames(df_pdo3_pdo9)[df_pdo3_pdo9$off_diag_1], rownames(df_pdo3_pdo9)[df_pdo3_pdo9$off_diag_2])))

gosim_lower <- give_goterm_list_of_genes(genes_in_list_ens = remove_na(c(rownames(df_pdo3_pdo9)[df_pdo3_pdo9$off_diag_2])),
                                         genes_not_in_list_ens = remove_na(c(rownames(df_pdo3_pdo9)[df_pdo3_pdo9$off_diag_1], rownames(df_pdo3_pdo9)[df_pdo3_pdo9$off_diag_2])))


gosim_upper$GOTerms$Term[order(gosim_upper$p.values)]
gosim_lower$GOTerms$Term[order(gosim_lower$p.values)]

gosim_upper$GOTerms$go_id[order(gosim_upper$p.values)[1]]
gosim_upper$GOTerms$Term[order(gosim_upper$p.values)[1]]
gosim_upper$GOTerms$Term[order(gosim_upper$p.values)[1]]

df_pdo3_pdo9[,"Genes_GO:0045471_ethanol"] = ifelse(rownames(df_pdo3_pdo9) %in% t2g_38$ensembl_gene_id[match(gosim_upper$genes[[gosim_upper$GOTerms$go_id[order(gosim_upper$p.values)][1]]], t2g_38$entrezgene_id)],
                                           yes = t2g_38$external_gene_name[match(rownames(df_pdo3_pdo9), t2g_38$ensembl_gene_id)], no=NA)

df_pdo3_pdo9[,"Genes_GO:0021700"] = ifelse(rownames(df_pdo3_pdo9) %in% t2g_38$ensembl_gene_id[match(gosim_lower$genes[[1]], t2g_38$entrezgene_id)],
                                           yes = t2g_38$external_gene_name[match(rownames(df_pdo3_pdo9), t2g_38$ensembl_gene_id)], no=NA)
df_pdo3_pdo9[,"Genes_allgo_upper"] = ifelse(rownames(df_pdo3_pdo9) %in% t2g_38$ensembl_gene_id[match(unlist(sapply(1:length(gosim_upper$GOTerms$go_id), function(j) gosim_upper$genes[[j]])), t2g_38$entrezgene_id)],
                                           yes = t2g_38$external_gene_name[match(rownames(df_pdo3_pdo9), t2g_38$ensembl_gene_id)], no=NA)
df_pdo3_pdo9[,"Genes_allgo_lower"] = ifelse(rownames(df_pdo3_pdo9) %in% t2g_38$ensembl_gene_id[match(unlist(sapply(1:length(gosim_lower$GOTerms$go_id), function(j) gosim_lower$genes[[j]])), t2g_38$entrezgene_id)],
                                            yes = t2g_38$external_gene_name[match(rownames(df_pdo3_pdo9), t2g_38$ensembl_gene_id)], no=NA)


ggplot(df_pdo3_pdo9, aes(x=PDO3, y=PDO9, col=col, label=`Genes_GO:0021700`))+geom_point()+
  scale_x_continuous(trans = 'log2')+scale_y_continuous(trans = 'log2')+
  geom_abline(slope = 1, intercept = 0, lty='dashed')+
  scale_fill_manual(values = c('blue', 'red', 'green'))+
  geom_label_repel()

ggplot(df_pdo3_pdo9, aes(x=PDO3, y=PDO9, col=col, label=`Genes_allgo_upper`))+geom_point()+
  scale_x_continuous(trans = 'log2')+scale_y_continuous(trans = 'log2')+
  geom_abline(slope = 1, intercept = 0, lty='dashed')+
  scale_fill_manual(values = c('blue', 'red', 'green'))+
  geom_label_repel()

ggplot(df_pdo3_pdo9, aes(x=PDO3, y=PDO9, col=col, label=`Genes_allgo_lower`))+geom_point()+
  scale_x_continuous(trans = 'log2')+scale_y_continuous(trans = 'log2')+
  geom_abline(slope = 1, intercept = 0, lty='dashed')+
  scale_fill_manual(values = c('blue', 'red', 'green'))+
  geom_label_repel()

ggplot(df_pdo3_pdo9, aes(x=PDO3, y=PDO9, col=col, label=`Genes_GO:0045471_ethanol`))+geom_point()+
  scale_x_continuous(trans = 'log2')+scale_y_continuous(trans = 'log2')+
  geom_abline(slope = 1, intercept = 0, lty='dashed')+
  scale_fill_manual(values = c('blue', 'red', 'green'))+
  geom_label_repel()


# library(KEGGREST)
# KEGGREST::keggList()

library(ReactomePA)
enrichment_reactome_upper <- ReactomePA::enrichPathway(gene = remove_na(t2g_38$entrezgene_id[match(rownames(df_pdo3_pdo9)[df_pdo3_pdo9$off_diag_1], t2g_38$ensembl_gene_id)]),
                          organism = "human")
enrichment_reactome_lower <- ReactomePA::enrichPathway(gene = remove_na(t2g_38$entrezgene_id[match(rownames(df_pdo3_pdo9)[df_pdo3_pdo9$off_diag_2], t2g_38$ensembl_gene_id)]),
                                                 organism = "human")
enrichment_reactome_upper@result
summary(enrichment_reactome_upper)[,1:6]
summary(enrichment_reactome_lower)[,1:6]

##---------
## Correlation of CN status in all genome
# dcast_joint_counts_CN0 <- joint_counts_CN0[,c('CN.gene_name', 'CN.value', 'CN.L1')]
# dcast_joint_counts_CN0 = outer(unique(dcast_joint_counts_CN0$CN.gene_name),
#                                unique(dcast_joint_counts_CN0$CN.gene_name),
#                                Vectorize(function(i,j){
#   .x1 <- dcast_joint_counts_CN0 %>% filter(CN.gene_name == i) %>% dplyr::select(CN.value)
#   .x2 <- dcast_joint_counts_CN0 %>% filter(CN.gene_name == j) %>% dplyr::select(CN.value)
#   return(cor(x = .x2[,1], y = .x1[match(rownames(.x2), rownames(.x1)),1]))
#                                             }))

deseq_counts$chrom = ag$seq_name[match(deseq_counts$Gene, ag$gene_name)]
head(deseq_counts)
table(deseq_counts$chrom)
deseq_counts$PDO= renaming$PDO[match(deseq_counts$DESEq.sample, renaming$ID)]

ggplot(deseq_counts, aes(x=PDO, y=DESeq.count))+geom_boxplot()
ggplot(deseq_counts %>% filter(grepl('PDO', PDO)),
       aes(x=PDO, y=DESeq.count))+geom_boxplot()
ggsave("../../figures_GRCh37/counts_persample_DESeq2.pdf", width=10, height=10)

ggplot(joint_counts_CN, aes(x=PDO, y=counts.value))+geom_boxplot()
ggsave("../../figures_GRCh37/counts_persample_TPM.pdf", width=10, height=10)

ggplot(joint_counts_CN0, aes(x=PDO, y=counts.value))+geom_boxplot()
ggsave("../../figures_GRCh37/counts_persample_TPM_2.pdf", width=10, height=10)
