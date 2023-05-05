rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

### files in  all_CN_states_per_gene/ created using script get_CN_of_genes.R

# library(tsne)
library(umap)
library(ggrepel)
library(dplyr)
library(grid)
library(reshape2)

dendrograminputclr <- readRDS("/Users/morril01/Documents/PhD/other_repos/Vias_Brenton/copy_number_analysis_organoids/robjects/cluster_clades_imput002clr.RDS")

ref_dir <- ""
data_flder <- "/Users/morril01/Documents/PhD/other_repos/Vias_Brenton/RNASeq_and_CN/20191218_ViasM_BJ_orgaBrs/Data/"
gtf.file <- file.path(data_flder,
                      "Homo_sapiens.GRCh37.87.gtf.gz")
sqlite_file <- 'Homo_sapiens.GRCh37.87.sqlite'
sqlite_path <- file.path(data_flder, sqlite_file)
if(!file.exists(sqlite_path)) {
  ## generate the SQLite database file
  ensembldb::ensDbFromGtf(gtf=gtf.file, path = '', outfile=sqlite_file)
}
EnsDb.Hsapiens.v87 <- ensembldb::EnsDb(sqlite_file)

# Genes, used to annotated the TPM matrix to send to Maria
ag <- ensembldb::genes(EnsDb.Hsapiens.v87, filter=list(AnnotationFilter::GeneBiotypeFilter('protein_coding')), return.type="DataFrame") 
ag_subsetchrom <- ag[!(ag$seq_name %in% c("MT", "X", "Y")) & !grepl("GL",ag$seq_name),]

flder <- "/Users/morril01/Documents/PhD/other_repos/Vias_Brenton/RNASeq_and_CN/20191218_ViasM_BJ_orgaBrs/output/all_CN_states_per_gene/"
fles <- list.files(flder, full.names = T)
fles <- fles[grepl('RDS', fles)]
length(fles)
gene_CN <- lapply(fles, readRDS)

basename(fles)[grepl('early', (fles))]

cn <- lapply(gene_CN, `[`, 2)
all_cn <- (do.call('cbind', cn))
colnames(all_cn) <- gsub(".RDS", "", basename(fles))
rownames(all_cn) <- make.names(ag_subsetchrom$gene_name, unique = T)

var_all_cn <- apply(all_cn, 1, var)
all_cn_var <- all_cn[var_all_cn>0,]

all_cn_var[ag_subsetchrom$seq_name == "8",]
all_cn[ag_subsetchrom$seq_name == "8",]


# image(as(all_cn, 'matrix'))
# image(as(all_cn[1:100,], 'matrix'))
# image(as(all_cn_var[1:100,], 'matrix'))


pheatmap::pheatmap(as(all_cn[1:100,], 'matrix'))
most_var <- as(all_cn[1:100,], 'matrix')
hclust_most_var <- hclust(dist(t(most_var)))
cutree_hclust_most_var <- cutree(hclust_most_var, k=2)
plot(hclust_most_var)
table(cutree_hclust_most_var)
remove_samples <- c('IM_311PS', 'IM_312PS', 'JBLAB-4211', 'JBLAB-4263', 'JBLAB-4263PS',
                    'JBLAB-4264', 'JBLAB-4264PS', 'JBLAB-4265') #names(cutree_hclust_most_var)[which(cutree_hclust_most_var == 2)]


##remove 12 samples which seem quite weird
all_cn <- all_cn[,!(colnames(all_cn) %in% remove_samples)]
all_cn_var <- all_cn_var[,!(colnames(all_cn_var) %in% remove_samples)]

min(all_cn) ## there shouldn't be any negative values!

dim(all_cn)

min_genes <- apply(all_cn, 1, min)
max_genes <- apply(all_cn, 1, max)
which_neg <- (min_genes) < 0
table(which_neg)
plot(which_neg, type='h')
plot(min_genes, type='h')
plot(max_genes, type='h')

# head(gene_CN[[1]])
# gene_CN

#------------------------------------------------------------------------#
#------------------------------------------------------------------------#
## UMAP
# umap_res <- umap(t(all_cn_var[1:2000,]))
umap_res <- umap(t(all_cn_var))
dev.off()

table(dendrograminputclr)
dendrograminputclr[dendrograminputclr == 1] <- 'S4-rich'
dendrograminputclr[dendrograminputclr == 2] <- 'S3-rich'
df_umap <- cbind.data.frame(umap_res$layout,
                 sample=colnames(all_cn_var),
                 clade=as.character(dendrograminputclr[match(colnames(all_cn_var), names(dendrograminputclr))]),
                 MYC=unlist(all_cn_var['MYC',]),
                 BRCA1=unlist(all_cn_var['BRCA1',]),
                 BRCA2=unlist(all_cn_var['BRCA2',]),
                 CCNE1=unlist(all_cn_var['CCNE1',]),
                 CDKN2D=unlist(all_cn_var['CDKN2D',])
)
categorise_samples <- function(i){
  if(grepl('JBLAB', i)){
    'JBLAB'
  }else if(grepl('IM', i)){
    'IM'
  }else if(grepl('TCGA', i)){
    'TCGA'
  }else{
    'ICGC'
  }
}
df_umap$dataset <- sapply(rownames(df_umap), categorise_samples)

ggplot(df_umap,
       aes(x=`1`, y=`2`, label=sample, shape=dataset))+
  geom_point( alpha=0.4, aes( color=clade, size=MYC))+
  # geom_label_repel()+
  theme_bw()
ggsave("figures/umap_segmentsMYC.pdf", width = 6.2, height = 5)

saveRDS(df_umap, "robjects/umap_ploidy_genes.RDS")

df_umap[df_umap$`1` > 80,]

ggplot(df_umap,
       aes(x=`1`, y=`2`, label=sample, shape=dataset))+
  geom_point( alpha=0.4, aes( color=clade, size=BRCA1))+
  theme_bw()
ggsave("figures/umap_segmentsBRCA1.pdf", width = 6.2, height = 5)

ggplot(df_umap,
       aes(x=`1`, y=`2`, label=sample, shape=dataset))+
  geom_point( alpha=0.4, aes( color=clade, size=CCNE1))+
  theme_bw()
ggsave("figures/umap_segmentsCCNE1.pdf", width = 6.2, height = 5)

tikzDevice::tikz("figures/umap_segmentsCCNE1.tex", width = 4.2, height = 4)
ggplot(df_umap,
       aes(x=`1`, y=`2`, label=sample, shape=dataset))+
  geom_point( alpha=0.4, aes( color=clade, size=CCNE1))+
  labs(shape='Dataset', size='CCNE1', col='', x='UMAP 1', y='UMAP 2')+
  theme_bw()+theme(legend.position = "bottom",legend.box="vertical", legend.margin=margin())
dev.off()

ggplot(df_umap,
       aes(x=`1`, y=`2`, label=sample, shape=dataset))+
  geom_point( alpha=0.4, aes( color=clade, size=CDKN2D))+
  theme_bw()
ggsave("figures/umap_segmentsCDKN2D.pdf", width = 6.2, height = 5)

ggplot(df_umap,
       aes(x=`1`, y=`2`, label=sample, shape=dataset))+
  geom_point( alpha=0.4, aes( color=clade, size=BRCA2))+
  theme_bw()
ggsave("figures/umap_segmentsBRCA2.pdf", width = 6.2, height = 5)


pca_excluding_samples <- function(features_mat, samples_to_exclude, size_point=10){
  bool_exclude <- (colnames(features_mat) %in% samples_to_exclude)
  .pca_res <- prcomp(t(features_mat[,!bool_exclude]), center = T, scale. = T)
  ggplot(cbind.data.frame(.pca_res$x[,1:2],
                          sample=rownames(.pca_res$x[,1:2]),
                          clade=as.character(dendrograminputclr[match(colnames(features_mat), names(dendrograminputclr))])[!bool_exclude]),
         aes(x=PC1, y=PC2,  label=sample))+
    geom_point(size=size_point, alpha=0.2, aes( color=clade))+geom_label_repel()+theme_bw()
}

pca_res <- prcomp(t(all_cn_var), center = T, scale. = T)

ggplot(cbind.data.frame(pca_res$x[,1:2],
                        sample=rownames(pca_res$x[,1:2]),
                        clade=as.character(dendrograminputclr[match(colnames(all_cn_var), names(dendrograminputclr))])),
       aes(x=PC1, y=PC2,  label=sample))+
  geom_point(size=10, alpha=0.2, aes( color=clade))+geom_label_repel()+theme_bw()
ggsave("figures/pca_segments.pdf", width = 6.2, height = 5)

pca_excluding_samples(features_mat = all_cn_var, samples_to_exclude = c('IM_35', 'IM_29'), size_point = 2)
ggsave("figures/pca_segments_excludedsamples.pdf", width = 6.2, height = 5)


subset_genes_of_interest = c('MYC', 'CCNE1', 'PIK3CA', 'TERT', 'KRAS', 'PTEN', 'RB1', 'AKT1',
                             'AKT2', 'PARP1', 'PARP2', 'ATM', 'ATR', 'WEE1', 'TOP1', 'TUBB1',
                             'AKT3', 'CCND1', 'CCND2', 'CCND3', 'CDKN2A', 'CDKN2B', 'MECOM', 'CDK12')

all_cn_subset_GoI <- all_cn[subset_genes_of_interest,]
ComplexHeatmap::Heatmap(all_cn_subset_GoI)

cor(all_cn_subset_GoI)
remove_cols_with_no_variance <- function(i)  t(i[,-1])[,!(apply(t(i[,-1]), 2, var) == 0)]

dendrograminputclr <- dendrograminputclr[match(colnames(all_cn_var), names(dendrograminputclr))]
colnames(all_cn_subset_GoI) ==  names(dendrograminputclr)
cor_CN_genes_of_interest <- cor(remove_cols_with_no_variance(all_cn_subset_GoI))
cor_CN_genes_of_interest_clade <- lapply(as.character(unique(dendrograminputclr)), function(clade_lab){
  cor(remove_cols_with_no_variance(all_cn_subset_GoI[,which(as.character(dendrograminputclr) == clade_lab)]))
})

ComplexHeatmap::Heatmap(cor_CN_genes_of_interest,
                        cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
                          grid.text(round(cor_CN_genes_of_interest[i, j], 1), x, y)
                        })

library(gridExtra)
ch1 <- ComplexHeatmap::Heatmap(cor_CN_genes_of_interest_clade[[1]], cluster_columns = F, cluster_rows = F,
                        cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
                          grid.text(round(cor_CN_genes_of_interest_clade[[1]][i, j], 1), x, y)
                        }, column_title = as.character(unique(dendrograminputclr))[1])
ch2 <- ComplexHeatmap::Heatmap(cor_CN_genes_of_interest_clade[[2]], cluster_columns = F, cluster_rows = F,
                        cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
                          grid.text(round(cor_CN_genes_of_interest_clade[[2]][i, j], 1), x, y)
                        }, column_title = as.character(unique(dendrograminputclr))[2])
ch3 <- ComplexHeatmap::Heatmap(cor_CN_genes_of_interest_clade[[3]],
                               cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
                                 grid.text(round(cor_CN_genes_of_interest_clade[[2]][i, j], 1), x, y)
                               })

ch1+ch2
plot(as.vector(cor_CN_genes_of_interest_clade[[1]]), as.vector(cor_CN_genes_of_interest_clade[[2]]))
abline(coef = c(0,1), lty='dashed')


## correlation of the most variable genes
cor_CN_genes_most_var <- cor(t(all_cn[order(var_all_cn, decreasing = T)[1:100],]))
ComplexHeatmap::Heatmap(cor_CN_genes_most_var)

ggplot(melt(all_cn_subset_GoI['MYC',]), aes(x=value))+geom_density()+#geom_vline(xintercept = all_cn_subset_GoI['MYC',1])
  geom_segment(aes(x = value, xend=value, y=-Inf, yend=Inf), alpha=0.2)+theme_bw()+scale_x_continuous(trans = "log2")

umap_genes_GoI <- umap(t(all_cn_subset_GoI))
df_umap_selectedgenes <- cbind.data.frame(umap_genes_GoI$layout,
                            sample=colnames(all_cn_var),
                            clade=as.character(dendrograminputclr[match(colnames(all_cn_var), names(dendrograminputclr))]),
                            MYC=unlist(all_cn_var['MYC',]),
                            BRCA1=unlist(all_cn_var['BRCA1',]),
                            BRCA2=unlist(all_cn_var['BRCA2',])
)
ggplot(df_umap_selectedgenes,
       aes(x=`1`, y=`2`, label=sample))+
  geom_point( alpha=0.2, aes( color=clade, size=BRCA1))+
  theme_bw()

umap_genes2 <- umap((all_cn_subset_GoI))
plot(umap_genes$layout)
plot(umap_genes2$layout) ## ??? looks very weird

give_umap_per_chrom <- function(chrom, genes_to_remove=NULL){
  dta <- t(all_cn_var[ag_subsetchrom$seq_name == chrom,])
  if(!is.null(genes_to_remove)){
    dta <- dta[,!(colnames(dta) %in% genes_to_remove)]
  }
  dta <- dta[,!grepl("NA.", colnames(dta))]
  na_idx <- which(colnames(dta) == "NA")
  if(length(na_idx)>0)  dta <- dta[,-na_idx]
  if(!is.null(dim(dta))){
    .umap_res <- umap(dta)
    return(cbind.data.frame(.umap_res$layout,
                                sample=colnames(all_cn_var),
                                clade=as.character(dendrograminputclr[match(colnames(all_cn_var), names(dendrograminputclr))]),
                                MYC=unlist(all_cn_var['MYC',]),
                                BRCA1=unlist(all_cn_var['BRCA1',]),
                                BRCA2=unlist(all_cn_var['BRCA2',]))
    )
  }else{
    NA
  }
}

pdf("figures/umap_segments_per_chrom.pdf", width = 4.2, height = 3)
sapply(gtools::mixedsort(unique(ag_subsetchrom$seq_name)), function(chrom_it){
  .umap <-  give_umap_per_chrom(chrom_it)
  cat('Chromosome ', chrom_it, '\n')
  if(!is.na(.umap)){
    print(ggplot(.umap,
           aes(x=`1`, y=`2`, label=sample))+
      geom_point( alpha=0.2, aes( color=clade, size=MYC))+
      theme_bw()+ggtitle(paste0('Chromosome ', chrom_it)))
  }else{
    print(ggplot()+ggtitle(paste0('Chromosome ', chrom_it)))
  }
})
dev.off()

system("open figures")

##----- multiple regression to see which genes determine the two classes
# cor_CN_genes_of_interest_clade <- lapply(as.character(unique(dendrograminputclr)), function(clade_lab){
#   cor(remove_cols_with_no_variance(

all_cn_var_clades <- all_cn_var[,!(is.na(dendrograminputclr))]
dendrograminputclr_no_na <- dendrograminputclr[!(is.na(dendrograminputclr))]
dim(all_cn_var_clades)
length(dendrograminputclr_no_na)
all_cn_var_clades_df <- (data.frame( dendrograminputclr_no_na, t(all_cn_var_clades)))
all_cn_var_clades_df[1:4,1:4]
all_cn_var_clades_df$dendrograminputclr_no_na <- factor(all_cn_var_clades_df$dendrograminputclr_no_na)

all_cn_var_clades_df_bin <- cbind.data.frame(dendrograminputclr_no_na=all_cn_var_clades_df[,1], apply(all_cn_var_clades_df[,-1], 2, function(i) as.numeric(i<2)))

dim(all_cn_var_clades_df)
dim(all_cn_var_clades_df_bin)

data_for_SVM <- all_cn_var_clades_df
data_for_SVM <- all_cn_var_clades_df_bin

dendrograminputclr_no_na
training_set_idx <- sample(1:nrow(data_for_SVM), size = floor(.7*nrow(data_for_SVM)))
training_set = data_for_SVM[training_set_idx,]
validation_set = data_for_SVM[!(1:nrow(data_for_SVM) %in% training_set_idx),]

library(e1071)
library(caret)

var_genes <- apply(training_set, 2, var)
subset_cols <- c(1, order(var_genes[-1], decreasing = T)[1:1000]) ## selecting the top X most variable genes #1:ncol(training_set)
length(subset_cols)
svmfit = svm(dendrograminputclr_no_na ~ ., data = training_set[,subset_cols], kernel = "linear", cost = 10, scale = FALSE)
print(svmfit)
# pheatmap::pheatmap(svmfit$SV)
# pheatmap::pheatmap(svmfit$coef0)

validation_predicted = predict(svmfit, validation_set[,subset_cols])
confusionMatrix(validation_predicted, validation_set$dendrograminputclr_no_na)

dim(svmfit$coefs)
dim(svmfit$SV)

order_coefs_svm <- order(svmfit$coefs, decreasing = T)
order_coefs_svm
colnames(training_set)

ggplot(all_cn_var_clades_df[,c('dendrograminputclr_no_na', 'MYC')], aes(x=dendrograminputclr_no_na, y = MYC))+
  geom_boxplot()+geom_jitter()+scale_y_continuous(trans = "log2")+theme_bw()

ggplot(all_cn_var_clades_df[,c('dendrograminputclr_no_na', 'CCNE1')], aes(x=dendrograminputclr_no_na, y = CCNE1))+
  geom_boxplot()+geom_jitter()+scale_y_continuous(trans = "log2")+theme_bw()

ggplot(all_cn_var_clades_df[,c('dendrograminputclr_no_na', 'CDKN2D')], aes(x=dendrograminputclr_no_na, y = CDKN2D))+
  geom_boxplot()+geom_jitter()+scale_y_continuous(trans = "log2")+theme_bw()

## ttest between the two groups
idx_s3 <- (all_cn_var_clades_df$dendrograminputclr_no_na == "S3-rich")
idx_s4 <- (all_cn_var_clades_df$dendrograminputclr_no_na == "S4-rich")

ttests_between_groups <- sapply(2:ncol(all_cn_var_clades_df), function(gene_idx){
  t.test(all_cn_var_clades_df[idx_s3,gene_idx], all_cn_var_clades_df[idx_s4,gene_idx])
})

ttests_between_groups[,1]
dim(ttests_between_groups)

## pval
order_pval_ttest <- order(unlist(ttests_between_groups[3,]))
order_pval_ttest

colnames(all_cn_var_clades_df)[-1][order_pval_ttest[1:100]]

plot(unlist(ttests_between_groups[1,]),
     -log10(unlist(ttests_between_groups[3,]))) ## volcano plot

genes_ttest_clades <- cbind.data.frame(tstat=unlist(ttests_between_groups[1,]),
                                       minlogpval=-log10(unlist(ttests_between_groups[3,])),
                                       minlogpvaladj=-log10(p.adjust(unlist(ttests_between_groups[3,]))),
                                       gene=colnames(all_cn_var_clades_df)[-1],
                                       chrom=ag$seq_name[match(colnames(all_cn_var_clades_df)[-1], ag$gene_name)])

ggplot(genes_ttest_clades, aes(x=factor(chrom, level=gtools::mixedsort(unique(chrom))), y=minlogpvaladj))+
  geom_boxplot()+geom_jitter(alpha=0.2)+theme_bw()+
  geom_hline(yintercept = -log10(0.05), lty='dashed', col='blue')
ggsave(filename = "figures/ttest_two_groups_per_chrom.pdf", width = 7)

ggplot(genes_ttest_clades,
       aes(x=tstat, y=minlogpvaladj, col=chrom,
           label=ifelse( gene %in% genes_ttest_clades[order(genes_ttest_clades$minlogpvaladj, decreasing = T)[1:10],'gene'],
                         yes = gene, no = NA)))+
       # aes(x=tstat, y=minlogpvaladj, col=chrom, label=ifelse(tstat<-26, gene, NA)))+
       # aes(x=tstat, y=minlogpvaladj, col=chrom, label=ifelse(gene == "TP53", gene, NA)))+
  geom_point()+theme_bw()+
  geom_label_repel(max.overlaps = 100)
ggsave("figures/ttest_two_groups_volcano.pdf", width = 10)

head(melt(all_cn_var, measure.vars = NULL))

all_cn_var_df <- (melt(as(all_cn_var, 'matrix')))
all_cn_var_df$chrom=ag$seq_name[match(all_cn_var_df$Var1, ag$gene_name)]
# ggplot(all_cn_var_df, aes(x=factor(chrom, level=gtools::mixedsort(unique(chrom))), ## too large
#                                   y=value))+geom_boxplot()+geom_jitter(alpha=0.2)+theme_bw()
all_cn_var_df$clade <- dendrograminputclr[match(all_cn_var_df$Var2, names(dendrograminputclr))]

dim(all_cn_var_df)
all_cn_var_df_summary <- all_cn_var_df %>% dplyr::group_by(Var1) %>% dplyr::summarise(mean_val=mean(value))
all_cn_var_df_summary2 <- all_cn_var_df %>% dplyr::group_by(Var1, clade) %>% dplyr::summarise(mean_val=mean(value))
dim(all_cn_var_df_summary)
all_cn_var_df_summary$chrom=ag$seq_name[match(all_cn_var_df_summary$Var1, ag$gene_name)]
all_cn_var_df_summary2$chrom=ag$seq_name[match(all_cn_var_df_summary2$Var1, ag$gene_name)]

ggplot(all_cn_var_df_summary, aes(x=as.numeric(chrom), group=chrom,
                                  y=mean_val, label=ifelse(mean_val<1.5, as.character(Var1), NA)))+
  geom_boxplot()+geom_jitter(alpha=0.2, width = 0.2)+geom_label_repel(nudge_x = 0)+theme_bw()

all_cn_var_df_summary2[is.na(all_cn_var_df_summary2$chrom),]
ggplot(all_cn_var_df_summary2[!is.na(all_cn_var_df_summary2$chrom),],
       aes(x=clade, group=clade, col=clade,
           y=mean_val, label=ifelse(mean_val<1.5, as.character(Var1), NA)))+
  # geom_boxplot()+
  facet_wrap(.~as.numeric(chrom), scales = "free_x")+
  geom_jitter(alpha=0.2, width = 0.2, aes(group=interaction(chrom, clade)))+
  # geom_label_repel(nudge_x = 0)+
  theme_bw()+labs(x='Clade and chromosome', y='Absolute CN')
ggsave("figures/average_CN_genes_per_chrom.pdf", width = 10, height = 10)
all_cn_var_df_summary2
# umap_genes <- umap((all_cn_var))
# plot(umap_genes$layout)
# plot(umap_genes$layout, xlim(-100, 100))
# plot(umap_genes$layout, xlim = c(-20, 20),  ylim = c(-20, 20), pch=19, cex=0.1) ## not very interesting

#------------------------------------------------------------------------#
#------------------------------------------------------------------------#
## primary and relapse
## read in relapse/primary information
load("../../../other_repos/britroc-1/data/britroc_30kb_signature_data.rds")
patient.meta


df_umap_britroc_rel_arx <- df_umap
df_umap_britroc_rel_arx <- df_umap_britroc_rel_arx[(grepl("*PS$", df_umap_britroc_rel_arx$sample)),]
df_umap_britroc_rel_arx$group = patient.meta$group[match(gsub('PS', '', df_umap_britroc_rel_arx$sample), patient.meta$SAMPLE_ID)]
df_umap_britroc_rel_arx$ploidy = patient.meta$ploidy[match(gsub('PS', '', df_umap_britroc_rel_arx$sample), patient.meta$SAMPLE_ID)]
df_umap_britroc_rel_arx$patient = patient.meta$PATIENT_ID[match(gsub('PS', '', df_umap_britroc_rel_arx$sample), patient.meta$SAMPLE_ID)]


## sensitivity and resistance
clinbrit <- read.csv("../../../other_repos/britroc-1/data/restricted/clinical/britroc_cohort_patient_data.tsv", sep = "\t")
clinbrit2 <- read.csv("../../../other_repos/britroc-1/data/restricted/clinical/britroc_patient_treatment_data.tsv", sep = "\t")
# df_umap_britroc_rel_arx$age <- clinbrit$age[match(as.character(df_umap_britroc_rel_arx$patient), paste0("BRITROC-", clinbrit$britroc_number))]
# df_umap_britroc_rel_arx$plat <- clinbrit$plat[match(as.character(df_umap_britroc_rel_arx$patient), paste0("BRITROC-", clinbrit$britroc_number))]
# df_umap_britroc_rel_arx$os <- clinbrit$os[match(as.character(df_umap_britroc_rel_arx$patient), paste0("BRITROC-", clinbrit$britroc_number))]
# df_umap_britroc_rel_arx$sensitivity <- clinbrit2$pt_sensitivity_at_reg[match(as.character(df_umap_britroc_rel_arx$patient), paste0("BRITROC-", clinbrit2$id))]
df_umap_britroc_rel_arx <- cbind(df_umap_britroc_rel_arx,
                                 clinbrit[match(as.character(df_umap_britroc_rel_arx$patient), paste0("BRITROC-", clinbrit$britroc_number)),])
df_umap_britroc_rel_arx <- cbind(df_umap_britroc_rel_arx,
                                 clinbrit2[match(as.character(df_umap_britroc_rel_arx$patient), paste0("BRITROC-", clinbrit2$idr)),])
ggplot(df_umap_britroc_rel_arx,
       aes(x=`1`, y=`2`, label=sample, shape=dataset))+
  geom_point( alpha=0.4, aes(size=age))+
  theme_bw()
ggplot(df_umap_britroc_rel_arx,
       aes(x=`1`, y=`2`, label=sample, shape=dataset))+
  geom_point( alpha=0.4, aes(size=1, col=sensitivity))+
  theme_bw()
ggplot(df_umap_britroc_rel_arx,
       aes(x=`1`, y=`2`, label=plat, shape=dataset))+
  geom_point( alpha=0.4, aes(size=1, col=sensitivity))+
  theme_bw()
ggplot(df_umap_britroc_rel_arx,
       aes(x=`1`, y=`2`, label=plat, shape=dataset))+
  geom_point( alpha=0.4, aes(size=os, col=sensitivity))+
  theme_bw()

df_umap_britroc_rel_arx <- df_umap_britroc_rel_arx[,!duplicated(colnames(df_umap_britroc_rel_arx))]

df_umap_britroc_rel_arx$int_end
df_umap_britroc_rel_arx$patient = as.character(df_umap_britroc_rel_arx$patient)
df_umap_britroc_rel_arx$patient_only_matched = as.character(df_umap_britroc_rel_arx$patient_only_matched)
pdf("figures/umap_all_colnames_clinical_.pdf")
for(cln in colnames(df_umap_britroc_rel_arx)){
  print(cln)
  if(!all(is.na(df_umap_britroc_rel_arx[,cln]))){
    .x <- ggplot(df_umap_britroc_rel_arx,
           aes(x=`1`, y=`2`, label=sample))+
      geom_point( alpha=0.4, aes(col=get(cln)))+
      theme_bw()+ggtitle(cln)
    if(!(typeof(df_umap_britroc_rel_arx[,cln]) %in% c("numeric", "integer", "double"))){
      if(length(unique(df_umap_britroc_rel_arx[,cln])) > 20){
      .x <- .x+guides(col=FALSE)
      }
    }
    print(.x)
  }
}
dev.off()

# patients_both_samples <- unique(df_umap_britroc_rel_arx$patient)[which(sapply(unique(df_umap_britroc_rel_arx$patient), function(i) all(c('arx', 'rlps') %in% df_umap_britroc_rel_arx$group[df_umap_britroc_rel_arx$patient == i])))]
# patients_both_samples <- patients_both_samples[!is.na(patients_both_samples)]
patients_both_samples <- as.character(unique(patient.meta$PATIENT_ID)[which(sapply(unique(patient.meta$PATIENT_ID), function(i) all(c('arx', 'rlps') %in% patient.meta$group[patient.meta$PATIENT_ID == i])))])

df_umap_britroc_rel_arx$patient_only_matched <- ifelse( (df_umap_britroc_rel_arx$patient %in% patients_both_samples),
                                                       yes = as.character(df_umap_britroc_rel_arx$patient),
                                                       no=NA)
ggplot()+
  geom_point(data = df_umap, aes(x=`1`, y=`2`), alpha=0.02, size=5)+
  geom_point(data = df_umap_britroc_rel_arx, aes(x=`1`, y=`2`, shape=group))+
  theme_bw()+
  geom_line(data = df_umap_britroc_rel_arx[!is.na(df_umap_britroc_rel_arx$patient_only_matched),],
            aes(x=`1`, y=`2`,group=patient, col=gsub('BRITROC-', 'BR', as.character(patient))))+
  # geom_label_repel(data = df_umap_britroc_rel_arx[!is.na(df_umap_britroc_rel_arx$patient_only_matched),],
  #                  aes(x=`1`, y=`2`, col=patient, group=patient, label=patient_only_matched), alpha=0.8,
  #                  max.overlaps = 10, size=2)+
  # theme(legend.position = "bottom", )+
  labs(col='', shape='')+guides(col=guide_legend(ncol=3))
ggsave("figures/umap_britroc_arx_rlps.pdf", width = 6.5, height = 6)

df_umap_britroc_rel_arx$group <- factor(df_umap_britroc_rel_arx$group, levels=c('arx', 'rlps'))
df_umap_britroc_rel_arx <- df_umap_britroc_rel_arx[order(df_umap_britroc_rel_arx$group),] ## sort by group
df_umap_britroc_rel_arx <- df_umap_britroc_rel_arx[with(df_umap_britroc_rel_arx, order(patient, group)), ]

## change_name_for_multisample_patients
# new_patient_name <- lapply(unique(df_umap_britroc_rel_arx$patient), function(patient_it){
new_patient_name <- lapply(patients_both_samples, function(patient_it){

combinations_list= c()
  ct = 1
  for(arx_it in which(patient.meta[patient.meta_subset$PATIENT_ID == patient_it,'group'] == "arx")){
    for(rlps_it in which(patient.meta_subset[patient.meta_subset$PATIENT_ID == patient_it,'group'] == "rlps")){
      sam_arx <- paste0(patient.meta_subset[patient.meta_subset$PATIENT_ID == patient_it,][arx_it,'SAMPLE_ID'])
      sam_rlps <- paste0(patient.meta_subset[patient.meta_subset$PATIENT_ID == patient_it,][rlps_it,'SAMPLE_ID'])
      combinations_list= c(combinations_list, paste0(sam_arx, "/", sam_rlps))
      ct = ct+1
    }
  }
  combinations_list
})
change_name_for_multisample_patients <- do.call('rbind', lapply(do.call('c', new_patient_name), function(pair_samples){
  .x <- df_umap_britroc_rel_arx[match(paste0(strsplit(pair_samples, '/')[[1]], 'PS'), df_umap_britroc_rel_arx$sample),]
  .x$pair <- pair_samples
  .x
}))
change_name_for_multisample_patients
change_name_for_multisample_patients$pair[is.na(change_name_for_multisample_patients$patient_only_matched)] <- NA
# change_name_for_multisample_patients <- change_name_for_multisample_patients[with(change_name_for_multisample_patients, order(pair, group)), ]
change_name_for_multisample_patients <- change_name_for_multisample_patients[order(change_name_for_multisample_patients$`2`),]

ggplot()+
  geom_point(data = df_umap, aes(x=`1`, y=`2`), alpha=0.02, size=5)+
  geom_point(data = df_umap_britroc_rel_arx, aes(x=`1`, y=`2`, shape=group))+
  theme_bw()+
  geom_path(data = df_umap_britroc_rel_arx[!is.na(df_umap_britroc_rel_arx$patient_only_matched),],
            aes(x=`1`, y=`2`,group=patient, colour=group ), arrow = arrow(length=unit(0.15,"cm"), ends="first"))+
  labs(col='', shape='')+guides(col=guide_legend(ncol=3))
ggsave("figures/umap_britroc_arx_rlps_2.pdf", width = 6.5, height = 6)


tikzDevice::tikz("figures/umap_britroc_arx_rlps_2.tex", width = 4.2, height = 3.5)
ggplot()+
  geom_point(data = df_umap, aes(x=`1`, y=`2`), alpha=0.02, size=5)+
  geom_point(data = df_umap_britroc_rel_arx, aes(x=`1`, y=`2`, shape=group))+
  theme_bw()+
  geom_path(data = df_umap_britroc_rel_arx[!is.na(df_umap_britroc_rel_arx$patient_only_matched),],
            aes(x=`1`, y=`2`,group=patient, colour=group ), arrow = arrow(length=unit(0.15,"cm"), ends="first"))+
  labs(col='', shape='', x='UMAP 1', y='UMAP 2')+guides(col=guide_legend(ncol=3))+
  theme(legend.position = "bottom")
dev.off()


ggplot()+
  geom_point(data = df_umap, aes(x=`1`, y=`2`), alpha=0.02, size=5)+
  geom_point(data = change_name_for_multisample_patients, aes(x=`1`, y=`2`, shape=group))+
  theme_bw()+
  geom_path(data = change_name_for_multisample_patients,
            aes(x=`1`, y=`2`,group=pair, colour=group ), arrow = arrow(length=unit(0.15,"cm"), ends="first"))+
  labs(col='', shape='')+guides(col=guide_legend(ncol=3))
ggsave("figures/umap_britroc_arx_rlps_3.pdf", width = 6.5, height = 6)

change_name_for_multisample_patients$pair_specific_patients <- change_name_for_multisample_patients$pair
change_name_for_multisample_patients$interesting_pairs <- as.character(sapply((change_name_for_multisample_patients$pair_specific_patients), function(i){
  .x <- change_name_for_multisample_patients[which(change_name_for_multisample_patients$pair_specific_patients == i),]
  (.x[.x$group == "arx", 2] > 2) & (.x[.x$group == "rlps", 2] < 2)
}))
change_name_for_multisample_patients$interesting_pairsINV <- as.character(sapply((change_name_for_multisample_patients$pair_specific_patients), function(i){
  .x <- change_name_for_multisample_patients[which(change_name_for_multisample_patients$pair_specific_patients == i),]
  (.x[.x$group == "arx", 2] < 0) & (.x[.x$group == "rlps", 2] > 2)
}))
change_name_for_multisample_patients$interesting_pairs = (change_name_for_multisample_patients$interesting_pairs == "TRUE")
change_name_for_multisample_patients$interesting_pairs
# change_name_for_multisample_patients$pair_specific_patients[change_name_for_multisample_patients$patient %in% ]
## for which patients there is a change of ploidy?
ggplot()+
  geom_point(data = df_umap, aes(x=`1`, y=`2`), alpha=0.02, size=5)+
  geom_point(data = change_name_for_multisample_patients[change_name_for_multisample_patients$interesting_pairs,],
             aes(x=`1`, y=`2`, shape=group))+
  theme_bw()+
  geom_path(data = change_name_for_multisample_patients[change_name_for_multisample_patients$interesting_pairs,],
            aes(x=`1`, y=`2`,group=pair, colour=group ), arrow = arrow(length=unit(0.15,"cm"), ends="first"))+
  # labs(col='', shape='', x='UMAP dim 1', y='UMAP dim 2')+guides(col=guide_legend(ncol=3))+facet_wrap(.~pair)
  labs(col='', shape='', x='UMAP dim 1', y='UMAP dim 2')+guides(col=guide_legend(ncol=3))+
  facet_wrap(.~factor(patient, levels=gtools::mixedsort(as.character(unique(patient)))), nrow=2)
ggsave("figures/umap_britroc_arx_rlps_3_facet.pdf", width = 7.5, height = 6)

patients_with_change <- change_name_for_multisample_patients[change_name_for_multisample_patients$interesting_pairs,]
patients_with_change <- patients_with_change[order(patients_with_change$patient),]
table(patients_with_change$patient)
View(patients_with_change)

## patient who has WGD in arx but not rlps
# change_name_for_multisample_patients[(change_name_for_multisample_patients$`1` < 0) & (change_name_for_multisample_patients$`2` > 0),]
change_name_for_multisample_patients[which(change_name_for_multisample_patients$interesting_pairsINV == "TRUE"),]
## BRITROC-67

dim(df_umap_britroc_rel_arx[!is.na(df_umap_britroc_rel_arx$patient_only_matched),])

table(sapply(unique(df_umap_britroc_rel_arx$patient), function(i){
  all(c('arx', 'rlps') %in% df_umap_britroc_rel_arx$group[df_umap_britroc_rel_arx$patient == i])
}))


### Are there differences in the patients who undergo WGD?
patients_who_undergo_WGD <- (unique(change_name_for_multisample_patients[change_name_for_multisample_patients$interesting_pairs,]$patient))
patients_who_dont_undergo_WGD <- unique(df_umap_britroc_rel_arx$patient[!(df_umap_britroc_rel_arx$patient %in% patients_who_undergo_WGD)])
clinbrit_patients_who_undergo_WGD <- clinbrit[match(gsub("BRITROC-", "", patients_who_undergo_WGD), clinbrit$britroc_number),]
clinbrit_patients_who_dont_undergo_WGD <- clinbrit[match(gsub("BRITROC-", "", patients_who_dont_undergo_WGD), clinbrit$britroc_number),]

## need to compare patients who undergo WGD to patients who remain without WGD

boxplot(list(clinbrit_patients_who_undergo_WGD$age,
             clinbrit_patients_who_dont_undergo_WGD$age))

table(sapply(unique(patient.meta$PATIENT_ID), function(i){
  all(c('arx', 'rlps') %in% patient.meta$group[patient.meta$PATIENT_ID == i])
}))

df_umap[which.max(df_umap$`1`),]
df_umap[(df_umap_britroc_rel_arx$`1`)>50,]
df_umap[(df_umap$`1`)>100,]

all_cn_var[,'JBLAB-4264PS']
all_cn_var[,'JBLAB-4264']


df_umap$early = grepl("^early*", df_umap$sample)
ggplot()+
  geom_point(data = df_umap[df_umap$early,], aes(x=`1`, y=`2`, col=early), alpha=0.7, size=3)+
  geom_point(data = df_umap[!df_umap$early,], aes(x=`1`, y=`2`, col=early), alpha=0.2, size=3, col='grey')+
  geom_label_repel(data = df_umap, aes(x=`1`, y=`2`, label=gsub("early_", "", ifelse(df_umap$early, df_umap$sample, NA))), alpha=0.7, size=3)+
  theme_bw()+
  labs(col='', shape='')+guides(col=F)
ggsave("figures/umap_early_samples.pdf", width = 6.5, height = 6)
system("open figures")



pdf("figures/umap_segments_per_chrom_early.pdf", width = 4.2, height = 3)
sapply(gtools::mixedsort(unique(ag_subsetchrom$seq_name)), function(chrom_it){
  .umap <-  give_umap_per_chrom(chrom_it)
  cat('Chromosome ', chrom_it, '\n')
  if(!is.na(.umap)){
    .umap$early = grepl("^early*", .umap$sample)
    ggplot()+
      geom_point(data = .umap[.umap$early,], aes(x=`1`, y=`2`, col=early), alpha=0.7, size=3)+
      geom_point(data = .umap[!.umap$early,], aes(x=`1`, y=`2`, col=early), alpha=0.2, size=3, col='grey')+
      geom_label_repel(data = .umap, aes(x=`1`, y=`2`, label=gsub("early_", "", ifelse(.umap$early, .umap$sample, NA))), alpha=0.7, size=2)+
      theme_bw()+
      labs(col='', shape='')+guides(col=F)+ggtitle(paste0('Chromosome ', chrom_it))
  }else{
    print(ggplot()+ggtitle(paste0('Chromosome ', chrom_it)))
  }
})
dev.off()


.umap4 <-  give_umap_per_chrom(4)
.umap4$early = grepl("^early*", .umap4$sample)
ggplot()+
  geom_point(data = .umap4[.umap4$early,], aes(x=`1`, y=`2`, col=early), alpha=0.7, size=3)+
  geom_point(data = .umap4[!.umap4$early,], aes(x=`1`, y=`2`, col=early), alpha=0.2, size=3, col='grey')+
  geom_label_repel(data = .umap4, aes(x=`1`, y=`2`, label=gsub("early_", "", ifelse(.umap4$early, .umap4$sample, NA))), alpha=0.7, size=2)+
  theme_bw()+
  labs(col='', shape='')+guides(col=F)+ggtitle(paste0('Chromosome ', 4))
ggsave("figures/umap_early_samples_chrom4.pdf", width = 6.5, height = 6)


chroms_vec_all_cn <- ag[match(rownames(all_cn), ag$gene_name),]
all_cn_chrom4 <- all_cn[which(chroms_vec_all_cn$seq_name == "4"),]
all_cn_chrom4
bool_early_col <- grepl('^early_', colnames(all_cn))
ggplot(melt(list(early=as(all_cn_chrom4[,bool_early_col], 'matrix'),
               late=as(all_cn_chrom4[,!bool_early_col], 'matrix'))),
       aes(x=value, col=L1))+
  geom_density()+scale_x_continuous(trans = "log10")+theme_bw()

melt_chrom4_early_late <- melt(list(early=as(all_cn_chrom4[,bool_early_col], 'matrix'),
          late=as(all_cn_chrom4[,!bool_early_col], 'matrix')))
ggplot(melt_chrom4_early_late,
       aes(x=value, col=L1))+
  geom_density()+scale_x_continuous(trans = "log10")+theme_bw()


melt_chrom4_early_late_summary_genes <- melt_chrom4_early_late %>% group_by(Var1, L1) %>% summarize(mean_cn=mean(value))
head(melt_chrom4_early_late_summary_genes)
melt_chrom4_early_late_summary_genes_dcast <- (dcast(melt_chrom4_early_late_summary_genes, Var1~L1, value.var = "mean_cn"))
melt_chrom4_early_late_summary_genes_dcast$logR <- log2(melt_chrom4_early_late_summary_genes_dcast$late)-log2(melt_chrom4_early_late_summary_genes_dcast$early)

min(all_cn)
min(melt_chrom4_early_late_summary_genes_dcast$early)
plot(melt_chrom4_early_late_summary_genes_dcast$early, melt_chrom4_early_late_summary_genes_dcast$late)
abline(coef = c(0,1), lty='dashed')

genes_negative_vals <- ag[match(melt_chrom4_early_late_summary_genes_dcast[melt_chrom4_early_late_summary_genes_dcast$early < 0,'Var1'], ag$gene_name),]

hist(melt_chrom4_early_late_summary_genes_dcast$logR)

head(melt_chrom4_early_late_summary_genes)

length(unique(chroms_vec_all_cn$gene_id))

## even removing these weird genes with negative values, there is the clear difference in the umap
.umap4_v2 <-  give_umap_per_chrom(4, genes_to_remove = genes_negative_vals$gene_name)
dim(.umap4_v2)
dim(.umap4)
.umap4_v2$early = grepl("^early*", .umap4_v2$sample)
ggplot()+
  geom_point(data = .umap4_v2[.umap4_v2$early,], aes(x=`1`, y=`2`, col=early), alpha=0.7, size=3)+
  geom_point(data = .umap4_v2[!.umap4_v2$early,], aes(x=`1`, y=`2`, col=early), alpha=0.2, size=3, col='grey')+
  geom_label_repel(data = .umap4_v2, aes(x=`1`, y=`2`, label=gsub("early_", "", ifelse(.umap4_v2$early, .umap4_v2$sample, NA))), alpha=0.7, size=2)+
  theme_bw()+
  labs(col='', shape='')+guides(col=F)+ggtitle(paste0('Chromosome ', 4))


melt_chrom4_early_late_summary_genes_dcast_order <- order(melt_chrom4_early_late_summary_genes_dcast$logR, na.last = NA)
length(melt_chrom4_early_late_summary_genes_dcast_order)
dim(melt_chrom4_early_late_summary_genes_dcast)
melt_chrom4_early_late_summary_genes_dcast[melt_chrom4_early_late_summary_genes_dcast_order[1:10],]
most_lost_early <- melt_chrom4_early_late_summary_genes_dcast[length(melt_chrom4_early_late_summary_genes_dcast_order):(length(melt_chrom4_early_late_summary_genes_dcast_order)-10),]


ggplot(melt_chrom4_early_late[melt_chrom4_early_late$Var1 %in% most_lost_early$Var1,], aes(x=L1, y=value))+geom_boxplot()+
  geom_jitter()+theme_bw()+facet_wrap(.~Var1)

cn_chrom4 <- (dcast(melt_chrom4_early_late, Var1~Var2, value.var = "value"))
rownames(cn_chrom4) <- cn_chrom4[,1]; cn_chrom4 <- cn_chrom4[,-1]
cn_chrom4 <- cn_chrom4[!(rownames(cn_chrom4) %in% genes_negative_vals$gene_name),]
cn_chrom4 <- cn_chrom4[apply(cn_chrom4, 1, var)>0,]
dim(cn_chrom4)
table(ag$seq_name == "4")
.pca_res_chrom4 <- prcomp(t(cn_chrom4), center = T, scale. = T)
plot(.pca_res_chrom4$x[,1:2], col=c('blue', 'orange')[factor(df_umap$early)])


melt_chrom4_early_late$location <- ag$gene_seq_start[match(melt_chrom4_early_late$Var1, ag$gene_name)]
melt_chrom4_early_late_summary_genes_dcast$location <- ag$gene_seq_start[match(melt_chrom4_early_late_summary_genes_dcast$Var1, ag$gene_name)]

ggplot(melt_chrom4_early_late_summary_genes_dcast, aes(x=location, y=logR,
                                                       col=(logR>0.65),
                                                       label=ifelse(logR>0.65, as.character(Var1), NA)))+geom_point()+
  geom_label_repel()+theme_bw()

all_cn_var_early <- all_cn_var[,grepl("^early_", colnames(all_cn_var))]
all_cn_var_late <- all_cn_var[,!grepl("^early_", colnames(all_cn_var))]
ttests_between_groups_early_late <- sapply(1:nrow(all_cn_var_late), function(gene_idx){
  t.test(all_cn_var_early[gene_idx,], all_cn_var_late[gene_idx,])
})


ttests_between_groups_early_late
order_pval_ttest_early_late <- order(unlist(ttests_between_groups_early_late[3,]))
order_pval_ttest_early_late

rownames(all_cn_var)[order_pval_ttest_early_late[1:100]]

df_ttest_early_late <- data.frame(stat=unlist(ttests_between_groups_early_late[1,]),
           mpval=-log10(unlist(ttests_between_groups_early_late[3,])),
           gene=rownames(all_cn_var))
df_ttest_early_late$chrom <- ag$seq_name[match(df_ttest_early_late$gene, ag$gene_name)]

ggplot(df_ttest_early_late,
       aes(x=stat, y=mpval, label=ifelse(mpval>20, gene, NA)))+geom_point()+theme_bw()+
  geom_label_repel()
ggsave("figures/ttest_early_late_volcano.pdf", width = 6.5, height = 6)

ggplot(df_ttest_early_late[df_ttest_early_late$mpval > 40,],
       aes(x=stat, y=mpval, col=chrom, label=ifelse(mpval>20, gene, NA)))+geom_point()+theme_bw()+
  geom_label_repel()+lims(y=c(40, 65))
ggsave("figures/ttest_early_late_volcano_zoom.pdf", width = 6.5, height = 6)


ggplot(melt_chrom4_early_late[melt_chrom4_early_late$Var1 %in% c('TUBB8'),], aes(x=L1, y=value))+geom_boxplot()+
  geom_jitter()+theme_bw()+facet_wrap(.~Var1)

## volcano plot

ggplot(melt(list(early=all_cn_var_early['TUBB8',],
               late=all_cn_var_late['TUBB8',])),
       aes(x=L1, y=value))+geom_boxplot()+geom_jitter()+theme_bw()+ggtitle('TUBB8')
ggsave("figures/TUBB8_early_late.pdf", width = 3, height = 3)

ggplot(melt(list(early=all_cn_var_early['CDKN1C',],
                 late=all_cn_var_late['CDKN1C',])),
       aes(x=L1, y=value))+geom_boxplot()+geom_jitter()+theme_bw()+ggtitle('CDKN1C')

subset_genes_of_interest = c('MYC', 'CCNE1', 'PIK3CA', 'TERT', 'KRAS', 'PTEN', 'RB1', 'AKT1',
                             'AKT2', 'PARP1', 'PARP2', 'ATM', 'ATR', 'WEE1', 'TOP1', 'TUBB1',
                             'AKT3', 'CCND1', 'CCND2', 'CCND3', 'CDKN2A', 'CDKN2B', 'MECOM', 'CDK12')
ggplot(df_ttest_early_late,
       aes(x=stat, y=mpval, label=ifelse(gene %in% subset_genes_of_interest, gene, NA)))+geom_point()+theme_bw()+
  geom_label_repel(max.overlaps = 40)
ggsave("figures/ttest_early_late_volcano_GoI.pdf", width = 6.5, height = 6)

ggplot(df_ttest_early_late,
       aes(x=stat, y=mpval, col=factor(chrom, levels=gtools::mixedsort(unique(df_ttest_early_late$chrom))),
           label=ifelse(stat<-18, gene, NA)))+geom_point()+theme_bw()+labs(col='')#+ geom_label_repel(max.overlaps = 40)
ggsave("figures/ttest_early_late_volcano_chrom.pdf", width = 6.5, height = 6)

df_ttest_early_late$chrom <- factor(df_ttest_early_late$chrom, levels=gtools::mixedsort(unique(df_ttest_early_late$chrom)))


ggplot(df_ttest_early_late,
       aes(x=stat, y=mpval,
           col=chrom, #shape=(chrom %in% levels(df_ttest_early_late$chrom)[c(T,F)]),
           label=ifelse(stat<15, gene, NA)))+geom_point()+theme_bw()+
  geom_label_repel(max.overlaps = 40)+lims(x=c(-20, -11))
ggsave("figures/ttest_early_late_volcano_annotated.pdf", width = 6.5, height = 6)

                          
ggplot(df_ttest_early_late,
       aes(x=stat, y=mpval, col=chrom, label=ifelse(grepl("^POTE*", gene), gene, NA)))+geom_point()+theme_bw()+
  geom_label_repel()
ggsave("figures/ttest_early_late_volcano_POTE.pdf", width = 6.5, height = 6)


all_cn_melt <- melt(as(all_cn, 'matrix'))
all_cn_melt_sample_summary <- (all_cn_melt) %>% dplyr::group_by(Var2) %>% dplyr::summarise(mean_cn = mean(value), var_cn=var(value))

head(all_cn_melt_sample_summary)
genes_ttest_clades$gene[which.max(genes_ttest_clades$minlogpvaladj)]

genes_ttest_clades[order(genes_ttest_clades$minlogpvaladj, decreasing = T)[1:10],]

all_cn_melt_sample_summary_df <- data.frame(mean_cn=all_cn_melt_sample_summary$mean_cn, var_cn=all_cn_melt_sample_summary$var_cn,
                                            early=factor(grepl('^early*', all_cn_melt_sample_summary$Var2)),
                                            sample=all_cn_melt_sample_summary$Var2)
all_cn_melt_sample_summary_df$POTEG_cn <- unlist(all_cn['POTEG',match(all_cn_melt_sample_summary_df$sample, colnames(all_cn))])
all_cn_melt_sample_summary_df$CHRNB2_cn <- unlist(all_cn['CHRNB2',match(all_cn_melt_sample_summary_df$sample, colnames(all_cn))])
ggplot(all_cn_melt_sample_summary_df, aes(x=mean_cn, y=var_cn, shape=early, col=CHRNB2_cn))+geom_point()+theme_bw()+
  scale_x_continuous(trans = "log2")+scale_y_continuous(trans = "log2")


##-------------------------------------------------------------------------------------##
## Question: is there an increase of LoH regions from primary to relapse?
## in BriTROC-1, how many genes are LoH?

patient.meta_subset <- patient.meta
patient.meta_subset0 <- match(paste0(as.character(patient.meta_subset$SAMPLE_ID), 'PS'), colnames(all_cn))
patient.meta_subset <- patient.meta_subset[-which(is.na(patient.meta_subset0)),]
all_cn_britroc <- all_cn[, match(paste0(as.character(patient.meta_subset$SAMPLE_ID), 'PS'), colnames(all_cn))]
dim(all_cn_britroc)

df_loss_britroc <- melt(list(arx=as(colMeans(all_cn_britroc[,patient.meta_subset$group == "arx"] < 2), 'matrix'),
          rlps=as(colMeans(all_cn_britroc[,patient.meta_subset$group == "rlps"] < 2), 'matrix')))
ggplot(df_loss_britroc,
       aes(x=value, col=L1))+geom_density()+theme_bw()

hist(colMeans(all_cn_britroc[patient.meta_subset$group == "arx",] < 2))
hist(colMeans(all_cn_britroc[patient.meta_subset$group == "rlps",] < 2))

changes_in_losses <- lapply(patients_both_samples, function(patient_it){
  ct = 1
  losses_list= list()
  cat(paste0('patient_it= "', patients_both_samples, '"', ";\n"))
  for(arx_it in which(patient.meta_subset[patient.meta_subset$PATIENT_ID == patient_it,'group'] == "arx")){
    for(rlps_it in which(patient.meta_subset[patient.meta_subset$PATIENT_ID == patient_it,'group'] == "rlps")){
      cat(paste0('arx_it= ', arx_it, '; rlps_it=', rlps_it, ";\n"))
      
      sam_arx <- paste0(patient.meta_subset[patient.meta_subset$PATIENT_ID == patient_it,][arx_it,'SAMPLE_ID'], 'PS')
      sam_rlps <- paste0(patient.meta_subset[patient.meta_subset$PATIENT_ID == patient_it,][rlps_it,'SAMPLE_ID'], 'PS')
      cat(paste0('sam_arx= "', sam_arx, '"; sam_rlps="', sam_rlps, '"\n'))
      losses_list[[paste0(sam_arx, '/', sam_rlps)]]= apply(cbind(all_cn[,as.character(sam_arx)] < 2,
            all_cn[,as.character(sam_rlps)] < 2), 1, paste0, collapse='')
      ct = ct+1
    }
  }
  return(losses_list)
})
changes_in_losses_thresh2 <- lapply(patients_both_samples, function(patient_it){
  ct = 1
  losses_list= list()
  for(arx_it in which(patient.meta_subset[patient.meta_subset$PATIENT_ID == patient_it,'group'] == "arx")){
    for(rlps_it in which(patient.meta_subset[patient.meta_subset$PATIENT_ID == patient_it,'group'] == "rlps")){
      sam_arx <- paste0(patient.meta_subset[patient.meta_subset$PATIENT_ID == patient_it,][arx_it,'SAMPLE_ID'], 'PS')
      sam_rlps <- paste0(patient.meta_subset[patient.meta_subset$PATIENT_ID == patient_it,][rlps_it,'SAMPLE_ID'], 'PS')
      cat(paste0('sam_arx= "', sam_arx, '"; sam_rlps="', sam_rlps, '"\n'))
      losses_list[[paste0(sam_arx, '/', sam_rlps)]]= apply(cbind(all_cn[,as.character(sam_arx)] < 1.2,
                                                                 all_cn[,as.character(sam_rlps)] <1.2), 1, paste0, collapse='')
      ct = ct+1
    }
  }
  return(losses_list)
})
changes_in_losses_thresh3 <- lapply(patients_both_samples, function(patient_it){
  ct = 1
  losses_list= list()
  for(arx_it in which(patient.meta_subset[patient.meta_subset$PATIENT_ID == patient_it,'group'] == "arx")){
    for(rlps_it in which(patient.meta_subset[patient.meta_subset$PATIENT_ID == patient_it,'group'] == "rlps")){
      sam_arx <- paste0(patient.meta_subset[patient.meta_subset$PATIENT_ID == patient_it,][arx_it,'SAMPLE_ID'], 'PS')
      sam_rlps <- paste0(patient.meta_subset[patient.meta_subset$PATIENT_ID == patient_it,][rlps_it,'SAMPLE_ID'], 'PS')
      cat(paste0('sam_arx= "', sam_arx, '"; sam_rlps="', sam_rlps, '"\n'))
      losses_list[[paste0(sam_arx, '/', sam_rlps)]]= apply(cbind((all_cn[,as.character(sam_arx)] < 1.2) & (all_cn[,as.character(sam_arx)] > 0.8),
                                                                 (all_cn[,as.character(sam_rlps)] < 1.2) & (all_cn[,as.character(sam_rlps)] > 0.8) ), 1, paste0, collapse='')
      ct = ct+1
    }
  }
  return(losses_list)
})
changes_in_CN_gene <- lapply(patients_both_samples, function(patient_it){
  ct = 1
  losses_list= list()
  cat(paste0('patient_it= "', patients_both_samples, '"', ";\n"))
  for(arx_it in which(patient.meta_subset[patient.meta_subset$PATIENT_ID == patient_it,'group'] == "arx")){
    for(rlps_it in which(patient.meta_subset[patient.meta_subset$PATIENT_ID == patient_it,'group'] == "rlps")){
      cat(paste0('arx_it= ', arx_it, '; rlps_it=', rlps_it, ";\n"))
      
      sam_arx <- paste0(patient.meta_subset[patient.meta_subset$PATIENT_ID == patient_it,][arx_it,'SAMPLE_ID'], 'PS')
      sam_rlps <- paste0(patient.meta_subset[patient.meta_subset$PATIENT_ID == patient_it,][rlps_it,'SAMPLE_ID'], 'PS')
      cat(paste0('sam_arx= "', sam_arx, '"; sam_rlps="', sam_rlps, '"\n'))
      losses_list[[paste0(sam_arx, '/', sam_rlps)]]= (all_cn[,as.character(sam_rlps)] - all_cn[,as.character(sam_arx)])
      ct = ct+1
    }
  }
  return(losses_list)
})

colnames(all_cn)[(grepl("IM_214", colnames(all_cn)))]

changes_in_CN_gene_log <- lapply(as.character(patients_both_samples), function(patient_it){
  ct = 1
  losses_list= list()
  cat(paste0('patient_it= "', patients_both_samples, '"', "\n"))
  
  for(arx_it in which(patient.meta_subset[patient.meta_subset$PATIENT_ID == patient_it,'group'] == "arx")){
    for(rlps_it in which(patient.meta_subset[patient.meta_subset$PATIENT_ID == patient_it,'group'] == "rlps")){
      cat(paste0('arx_it= ', arx_it, '; rlps_it=', rlps_it, "\n"))
      sam_arx <- paste0(patient.meta_subset[patient.meta_subset$PATIENT_ID == patient_it,][arx_it,'SAMPLE_ID'], 'PS')
      sam_rlps <- paste0(patient.meta_subset[patient.meta_subset$PATIENT_ID == patient_it,][rlps_it,'SAMPLE_ID'], 'PS')
      cat(paste0('sam_arx= "', sam_arx, '"; sam_rlps="', sam_rlps, '"\n'))
      losses_list[[paste0(sam_arx, '/', sam_rlps)]]= (log2(all_cn[,as.character(sam_rlps)]) - log2(all_cn[,as.character(sam_arx)]))
      ct = ct+1
    }
  }
  return(losses_list)
})
CN_arx_gene <- lapply(patients_both_samples, function(patient_it){
  ct = 1
  losses_list= list()
  for(arx_it in which(patient.meta_subset[patient.meta_subset$PATIENT_ID == patient_it,'group'] == "arx")){
    for(rlps_it in which(patient.meta_subset[patient.meta_subset$PATIENT_ID == patient_it,'group'] == "rlps")){
      sam_arx <- paste0(patient.meta_subset[patient.meta_subset$PATIENT_ID == patient_it,][arx_it,'SAMPLE_ID'], 'PS')
      sam_rlps <- paste0(patient.meta_subset[patient.meta_subset$PATIENT_ID == patient_it,][rlps_it,'SAMPLE_ID'], 'PS')
      losses_list[[paste0(sam_arx, '/', sam_rlps)]]= (all_cn[,as.character(sam_arx)])
      ct = ct+1
    }
  }
  return(losses_list)
})

give_stats_gene_CN_changes <- function(subset_patients){
  .res <- lapply(subset_patients, function(patient_it){
    ct = 1
    losses_list= list()
    for(arx_it in which(patient.meta_subset[as.character(patient.meta_subset$PATIENT_ID) == patient_it,'group'] == "arx")){
      for(rlps_it in which(patient.meta_subset[as.character(patient.meta_subset$PATIENT_ID) == patient_it,'group'] == "rlps")){
        sam_arx <- paste0(patient.meta_subset[patient.meta_subset$PATIENT_ID == patient_it,][arx_it,'SAMPLE_ID'], 'PS')
        sam_rlps <- paste0(patient.meta_subset[patient.meta_subset$PATIENT_ID == patient_it,][rlps_it,'SAMPLE_ID'], 'PS')
        cat(paste0('sam_arx= "', sam_arx, '"; sam_rlps="', sam_rlps, '"\n'))
        losses_list[[paste0(sam_arx, '/', sam_rlps)]]= cbind(min1p2=apply(cbind(all_cn[,as.character(sam_arx)] < 1.2,
                                                                                all_cn[,as.character(sam_rlps)] <1.2), 1, paste0, collapse=''),
                                                             min2=apply(cbind(all_cn[,as.character(sam_arx)] < 2,
                                                                              all_cn[,as.character(sam_rlps)] < 2), 1, paste0, collapse=''),
                                                             difference_in_cn=(all_cn[,as.character(sam_rlps)] - all_cn[,as.character(sam_arx)]),
                                                             logR=(log2(all_cn[,as.character(sam_rlps)]) - log2(all_cn[,as.character(sam_arx)])),
                                                             arxCN= (all_cn[,as.character(sam_arx)]),
                                                             gene=rownames(all_cn),
                                                             patient=patient_it,
                                                             comparison=paste0(sam_arx, '/', sam_rlps)
        )
        ct = ct+1
      }
    }
    return(losses_list)
  })
  return(.res)
}

wrapper_give_stats_gene_CN_changes <- function(subset_patients){
  .x <- give_stats_gene_CN_changes(subset_patients)
  .x = do.call('c', .x)
  .x = do.call('rbind.data.frame', .x)
  .x$difference_in_cn <- as.numeric(.x$difference_in_cn )
  .x$logR <- as.numeric(.x$logR )
  .x$arxCN <- as.numeric(.x$arxCN )
  return(.x)
}

stats_gene_CN_changes_WGDsamples <- wrapper_give_stats_gene_CN_changes(c('BRITROC-209', 'BRITROC-216', 'BRITROC-23', 'BRITROC-241',
                                                                         'BRITROC-267', 'BRITROC-274', 'BRITROC-74'))
stats_gene_CN_changes_nonWGDsamples <- wrapper_give_stats_gene_CN_changes(as.character(patients_both_samples)[-match(c('BRITROC-209', 'BRITROC-216', 'BRITROC-23', 'BRITROC-241',
                                                                                                                       'BRITROC-267', 'BRITROC-274', 'BRITROC-74'), as.character(patients_both_samples))])

dim(stats_gene_CN_changes_WGDsamples)

table(sapply(changes_in_losses, length))
changes_in_losses_lists <- do.call('c', changes_in_losses)
changes_in_losses_thresh2_lists <- do.call('c', changes_in_losses_thresh2)
changes_in_losses_thresh3_lists <- do.call('c', changes_in_losses_thresh3)
changes_in_CN_gene_lists <- do.call('c', changes_in_CN_gene)
changes_in_CN_gene_log_lists <- do.call('c', changes_in_CN_gene_log)
CN_arx_gene_lists <- do.call('c', CN_arx_gene)

length(changes_in_losses_lists)
length(changes_in_losses)

changes_in_losses_lists <- do.call('rbind', changes_in_losses_lists)
changes_in_losses_thresh2_lists <- do.call('rbind', changes_in_losses_thresh2_lists)
changes_in_losses_thresh3_lists <- do.call('rbind', changes_in_losses_thresh3_lists)

changes_in_losses_lists_equal <- (apply(changes_in_losses_lists, 2, function(i) (i == "FALSEFALSE") | (i == "TRUETRUE") ))
dim(changes_in_losses_lists_equal)
dim(changes_in_losses_lists)
changes_in_losses_lists_increaseloss <- (apply(changes_in_losses_lists, 2, function(i) (i == "FALSETRUE") ))
dim(changes_in_losses_lists_increaseloss)

# pheatmap::pheatmap(changes_in_losses_lists_equal) ## too large

## we classify the genes in each pair of archival/relapse into
## (a) both with CN < 2
## (b) none with CN < 2
## (c) and (d) one of each
changes_in_losses_lists_table <- apply(changes_in_losses_lists, 1, function(i){
  table(factor(i, levels=c('TRUETRUE', 'FALSETRUE', 'TRUEFALSE', 'FALSEFALSE')))
  })
changes_in_losses_lists_table[,1:5]
pdf("figures/genes_lower2CN.pdf", width = 12.5, height = 6)
print(pheatmap::pheatmap(sweep(changes_in_losses_lists_table, 2, colSums(changes_in_losses_lists_table), '/')))
dev.off()
## the samples are split into two:
## the ones with a lot of CN in general (FALSEFALSE)
## the ones with many losses (although they could also be diploid)
## there is a third, smaller, category, of TRUEFALSE samples, i.e. the ones that undergo WGD (confirmed with umap results)
hist(changes_in_CN_gene_lists[[1]])
changes_in_losses_thresh2_lists_table <- apply(changes_in_losses_thresh2_lists, 1, function(i){
  table(factor(i, levels=c('TRUETRUE', 'FALSETRUE', 'TRUEFALSE', 'FALSEFALSE')))
})
changes_in_losses_thresh3_lists_table <- apply(changes_in_losses_thresh3_lists, 1, function(i){
  table(factor(i, levels=c('TRUETRUE', 'FALSETRUE', 'TRUEFALSE', 'FALSEFALSE')))
})
changes_in_losses_thresh2_lists_table[,1:5]
pdf("figures/genes_lower1p2CN.pdf", width = 12.5, height = 6)
print(pheatmap::pheatmap(sweep(changes_in_losses_thresh2_lists_table, 2, colSums(changes_in_losses_thresh2_lists_table), '/')))
dev.off()

pdf("figures/genes_roughly1.pdf", width = 12.5, height = 6)
print(pheatmap::pheatmap(sweep(changes_in_losses_thresh3_lists_table, 2, colSums(changes_in_losses_thresh3_lists_table), '/')))
dev.off()

## so, what are the genes with a tendency to lose LoH from early to late?
changes_in_CN_gene_lists_mean <- sapply(changes_in_CN_gene_lists, mean)
changes_in_CN_gene_lists_sd <- sapply(changes_in_CN_gene_lists, sd)
changes_in_CN_gene_lists_med <- sapply(changes_in_CN_gene_lists, median)

changes_in_CN_gene_lists_mean_pergene <- apply(do.call('rbind', changes_in_CN_gene_lists), 2, mean, na.rm=T)
changes_in_CN_gene_lists_med_pergene <- apply(do.call('rbind', changes_in_CN_gene_lists), 2, median, na.rm=T)
changes_in_CN_gene_lists_sd_pergene <- apply(do.call('rbind', changes_in_CN_gene_lists), 2, sd, na.rm=T)

changes_in_CN_gene_log_lists_mean_pergene <- apply(do.call('rbind', changes_in_CN_gene_log_lists), 2, mean, na.rm=T)

CN_arx_gene_lists_mean_pergene <- apply(do.call('rbind', CN_arx_gene_lists), 2, mean, na.rm=T)

plot(changes_in_CN_gene_lists_mean, changes_in_CN_gene_lists_sd)
plot(changes_in_CN_gene_lists_mean, changes_in_CN_gene_lists_med)

plot(changes_in_CN_gene_lists_mean_pergene, changes_in_CN_gene_lists_sd_pergene)

changes_per_gene <- cbind.data.frame(mean_changegene=changes_in_CN_gene_lists_mean_pergene,
                 median_changegene=changes_in_CN_gene_lists_med_pergene,
                 sd_changegene=changes_in_CN_gene_lists_sd_pergene,
                 gene=rownames(all_cn),
                 logR_mean_changes=changes_in_CN_gene_log_lists_mean_pergene,
                 arx_CN=CN_arx_gene_lists_mean_pergene)

ggplot(changes_per_gene, aes(x=mean_changegene, y=sd_changegene, label=ifelse(abs(mean_changegene)>.5, gene, NA)))+
  geom_point()+geom_label_repel(max.overlaps = 40)+theme_bw()

ggplot(changes_per_gene, aes(x=median_changegene, y=sd_changegene, label=ifelse(abs(median_changegene)>.15, gene, NA)))+
  geom_point()+geom_label_repel(max.overlaps = 40)+theme_bw()

ggplot(changes_per_gene, aes(x=mean_changegene, y=sd_changegene, label=ifelse(gene %in% subset_genes_of_interest, gene, NA)))+
  geom_point()+geom_label_repel(max.overlaps = 40)+theme_bw()

ggplot(changes_per_gene, aes(x=median_changegene, y=sd_changegene, label=ifelse(gene %in% subset_genes_of_interest, gene, NA)))+
  geom_point()+geom_label_repel(max.overlaps = 40)+theme_bw()

ggplot(changes_per_gene, aes(x=median_changegene, y=sd_changegene, label=ifelse(gene %in% c('CDKN1C'), gene, NA)))+
  geom_point()+geom_label_repel(max.overlaps = 40)+theme_bw()

ggplot(changes_per_gene, aes(x=mean_changegene, y=logR_mean_changes, label=ifelse(abs(mean_changegene)>2, gene, NA)))+
  geom_point()+geom_label_repel(max.overlaps = 40)+theme_bw()

ggplot(changes_per_gene, aes(x=arx_CN, y=mean_changegene, label=ifelse(abs(mean_changegene)>2, gene, NA)))+
  geom_point()+geom_label_repel(max.overlaps = 40)+theme_bw()

ggplot(changes_per_gene, aes(x=arx_CN, y=logR_mean_changes, label=ifelse((logR_mean_changes < -0.1) & (arx_CN)<2, gene, NA)))+
  geom_point()+geom_label_repel(max.overlaps = 40)+theme_bw()

## what about the samples that undergo WGD - do they have some genes that are lost more than others?
stats_gene_CN_changes_WGDsamples_summary = stats_gene_CN_changes_WGDsamples %>% dplyr::group_by(gene) %>% 
  dplyr::summarise(mean_difference_in_cn = mean(difference_in_cn, na.rm = T), mean_logR = mean(logR, na.rm = T),
                   mean_arx_cn=mean(arxCN, na.rm = T))

normalise <- function(i) i/sum(i)

## we have a vert slight increase of genes that are now LoH, but didn't use to be
normalise(table(factor(stats_gene_CN_changes_WGDsamples$min1p2, levels=unique(stats_gene_CN_changes_WGDsamples$min1p2))))
normalise(table(factor(stats_gene_CN_changes_nonWGDsamples$min1p2, levels=unique(stats_gene_CN_changes_nonWGDsamples$min1p2))))
# apply(stats_gene_CN_changes_WGDsamples$min1p2, 2, function(i) table(factor(i, levels=unique(stats_gene_CN_changes_WGDsamples$min1p2))))

head(stats_gene_CN_changes_WGDsamples_summary)
ggplot(stats_gene_CN_changes_WGDsamples, )


stats_gene_CN_changes_nonWGDsamples

length(changes_in_losses_lists)
library(igraph)

table_to_graph <- function(tab, ...){
  adj <- matrix(0, ncol=5, nrow=5)
  colnames(adj) <- rownames(adj) <- c('diploid', 'arx loss', 'arx no loss', 'rlps loss', 'rlps no loss')
  adj['diploid','arx loss'] <- sum(tab[grepl('^TRUE', names(tab))])
  adj['diploid','arx no loss'] <- sum(tab[grepl('^FALSE', names(tab))])
  adj['arx loss', 'rlps loss'] <- sum(tab[grepl('^TRUETRUE', names(tab))])
  adj['arx no loss', 'rlps loss'] <- sum(tab[grepl('^FALSETRUE', names(tab))])
  adj['arx loss', 'rlps no loss'] <- sum(tab[grepl('^TRUEFALSE', names(tab))])
  adj['arx no loss', 'rlps no loss'] <- sum(tab[grepl('^FALSEFALSE', names(tab))])
  adj <- sweep(adj, 1, rowSums(adj), '/')
  adj[is.na(adj)] <- 0
  # for(i in 1:5){for(j in 1:i){adj[i,j] = adj[j,i]}}
  NodeList <- cbind(x=c(1,1,3,2,3), y=c(0,2.5,2.5,6,0))
  graph <- igraph::graph_from_adjacency_matrix(adj, weighted = T, mode = "directed")
  # LO = layout_nicely(graph); LO[,1] <- LO[,1]+3; LO[,2] <- LO[,2]+3
  plot(graph, edge.width=igraph::E(graph)$weight*6, layout=NodeList,
       edge.label = paste0('  ', round(igraph::E(graph)$weight, 2), '        '), ...)
}
table(changes_in_losses_lists[1,])

table_to_graph(table(changes_in_losses_lists[1,]))

##' first cluster samples according to the type of changes that they undergo, and then
##' create this adjacency matrix and graph for each of the sample subgroups

## compositional
tree_changes_in_losses_lists_table <- (hclust(dist(compositions::clr(t(changes_in_losses_lists_table)))))
plot(tree_changes_in_losses_lists_table) ## 3 groups
tree_changes_in_losses_thresh2_lists_table <- (hclust(dist(compositions::clr(t(changes_in_losses_thresh2_lists_table)))))
plot(tree_changes_in_losses_thresh2_lists_table) ## 2 groups
tree_changes_in_losses_thresh3_lists_table <- (hclust(dist(compositions::clr(t(changes_in_losses_thresh3_lists_table)))))
plot(tree_changes_in_losses_thresh3_lists_table) ## 3 groups (including one outlier)
tree_changes_in_losses_lists_table_cutree <- cutree(tree_changes_in_losses_lists_table, k=3)
tree_changes_in_losses_thresh2_lists_table_cutree <- cutree(tree_changes_in_losses_thresh2_lists_table, k=2)
tree_changes_in_losses_thresh3_lists_table_cutree <- cutree(tree_changes_in_losses_thresh3_lists_table, k=3)
tree_changes_in_losses_lists_table$labels == names(tree_changes_in_losses_lists_table_cutree)
tree_changes_in_losses_lists_table$labels = (tree_changes_in_losses_lists_table_cutree)
plot(tree_changes_in_losses_lists_table)

dev.off()

pdf("figures/loss_graph_stratify.pdf", width = 12.5, height = 4)
par(mfrow=c(1,3))
for(k_it in 1:3){
  subset_clustering_loss <- tree_changes_in_losses_lists_table_cutree[(tree_changes_in_losses_lists_table_cutree == k_it)]
  table_to_graph(table(as.vector(changes_in_losses_lists[names(subset_clustering_loss),])),
                 main=paste0('Sample group ', k_it, ' (', length(subset_clustering_loss), ')'),
                 vertex.size=40, vertex.label.color = "black", vertex.color = "white",
                 edge.label.family='mono', vertex.label.family='mono', edge.color='black')
}
dev.off()

pdf("figures/loss_graph_stratify_thresh1p2.pdf", width = 12, height = 7)
par(mfrow=c(1,2))
for(k_it in 1:2){
  subset_clustering_loss <- tree_changes_in_losses_thresh2_lists_table_cutree[(tree_changes_in_losses_thresh2_lists_table_cutree == k_it)]
  table_to_graph(table(as.vector(changes_in_losses_thresh2_lists[names(subset_clustering_loss),])),
                 main=paste0('Sample group ', k_it, ' (', length(subset_clustering_loss), ')'),
                 vertex.size=40, vertex.label.color = "black", vertex.color = "white",
                 edge.label.family='mono', vertex.label.family='mono', edge.color='black')
}
dev.off()

pdf("figures/loss_graph_stratify_roughly1.pdf", width = 12, height = 7)
par(mfrow=c(1,2))
for(k_it in 1:2){
  subset_clustering_loss <- tree_changes_in_losses_thresh3_lists_table_cutree[(tree_changes_in_losses_thresh3_lists_table_cutree == k_it)]
  table_to_graph(table(as.vector(changes_in_losses_thresh3_lists[names(subset_clustering_loss),])),
                 main=paste0('Sample group ', k_it, ' (', length(subset_clustering_loss), ')'),
                 vertex.size=40, vertex.label.color = "black", vertex.color = "white",
                 edge.label.family='mono', vertex.label.family='mono', edge.color='black')
}
dev.off()

table(changes_in_losses_thresh3[[2]][[1]] ==   changes_in_losses_thresh3[[1]][[1]])

##-------------------------------------------------------------------------------------##
## Question: is there genes of consitent increase/decrease from primary to relapse?

##' data (computed above): difference in CN between archival and relapse (relapse - archival)
##' for any archival-relapse pair, including the same sample multiple times if they are multi-sample

changes_in_CN_gene_df <- do.call('rbind', changes_in_CN_gene_lists)

rownames(changes_in_CN_gene_df) <- names(changes_in_CN_gene_lists)
colnames(changes_in_CN_gene_df) <- rownames(all_cn)

colmeans_difference_genes <- colMeans(changes_in_CN_gene_df)
sort(colmeans_difference_genes)

colmeans_difference_genes_df <- data.frame(colmeans_difference_genes=colmeans_difference_genes,
                                           colmedians_difference_genes= apply(changes_in_CN_gene_df, 2, median),
                                           gene=colnames(changes_in_CN_gene_df))

head(colmeans_difference_genes_df)

library(viridis)
library(jcolors)
library(RColorBrewer)
ggplot(colmeans_difference_genes_df, aes(x= colmeans_difference_genes, col=(colmeans_difference_genes),
                                         label=ifelse(abs(colmeans_difference_genes)>0.5, gene, NA)))+
  geom_density()+
  geom_label_repel(aes(y=0), size=2, max.overlaps = Inf, box.padding = .4, nudge_y = 1)+theme_bw()+
  labs(x='Average difference in CN between archival and relapse',
       y='Density')+guides(col=F)+
  # scale_color_viridis() +
  scale_color_jcolors_contin(palette = "pal3")
ggsave("figures/most_extreme_differences_archival_relapse.pdf", width = 5, height = 4)

changes_in_CN_gene_df[,grepl('ERBB2', colnames(changes_in_CN_gene_df))]
boxplot(changes_in_CN_gene_df[,which(colnames(changes_in_CN_gene_df) == 'ERBB2')]) ## one extreme case

colnames(changes_in_CN_gene_df)[grepl('ERBB2', colnames(changes_in_CN_gene_df))]

View(colmeans_difference_genes_df[abs(colmeans_difference_genes_df$colmeans_difference_genes) > 0.5, ] %>% dplyr::arrange(colmeans_difference_genes))

colmeans_difference_genes_df$chrom = ag$seq_name[match(colmeans_difference_genes_df$gene, ag$gene_name)]
table(colmeans_difference_genes_df$chrom)

colmeans_difference_genes_df$chrom <- factor(as.numeric(colmeans_difference_genes_df$chrom ))
ggplot(colmeans_difference_genes_df, aes(x=colmeans_difference_genes, y=colmedians_difference_genes,
                                         col=factor(as.numeric(chrom)),
                                         label=ifelse(abs(colmedians_difference_genes) > 0.18, gene, NA)
       ))+geom_point(shape=1)+geom_label_repel(max.overlaps = Inf, size=2)+theme_bw()+
  scale_color_manual(values = colorRampPalette(brewer.pal(12, "Accent"))(24))+
  labs(col='Chromosome', x='Mean CN difference',  y='Median CN difference')
ggsave("figures/most_extreme_differences_archival_relapse_mean_median.pdf", width = 5, height = 4)

ggplot(colmeans_difference_genes_df, aes(x=colmeans_difference_genes, y=colmedians_difference_genes,
                                         col=chrom,
                                         label=ifelse(abs(colmedians_difference_genes) > 0.18, gene, NA)
))+geom_vline(xintercept = 0, lty='dashed')+geom_hline(yintercept = 0, lty='dashed')+
  geom_point(shape=1)+geom_label_repel(max.overlaps = Inf, size=2)+theme_bw()+
  scale_color_manual(values = colorRampPalette(brewer.pal(12, "Accent"))(24))+
  labs(col='Chromosome', x='Mean CN difference',  y='Median CN difference')+facet_wrap(.~chrom)
ggsave("figures/most_extreme_differences_archival_relapse_mean_median_facet.pdf", width = 8, height = 8)

ggplot(colmeans_difference_genes_df[which(colmeans_difference_genes_df$chrom == 17),], aes(x=colmeans_difference_genes, y=colmedians_difference_genes,
                                         col=chrom,
                                         label=ifelse(abs(colmedians_difference_genes) > 0.15, gene, NA)
))+geom_vline(xintercept = 0, lty='dashed')+geom_hline(yintercept = 0, lty='dashed')+
  geom_point(shape=1)+geom_label_repel(max.overlaps = Inf, size=2)+theme_bw()+
  scale_color_manual(values = colorRampPalette(brewer.pal(12, "Accent"))(24))+
  labs(col='Chromosome', x='Mean CN difference',  y='Median CN difference')+facet_wrap(.~chrom)


# selected_chrom17 <- colmeans_difference_genes_df[colmeans_difference_genes_df$chrom == 17,][order(colmeans_difference_genes_df[colmeans_difference_genes_df$chrom == 17,]$colmedians_difference_genes)[1:10],'gene']
selected_chrom17 <- c(colmeans_difference_genes_df[colmeans_difference_genes_df$chrom == 17,][order(colmeans_difference_genes_df[colmeans_difference_genes_df$chrom == 17,]$colmedians_difference_genes, decreasing = T)[1:10],'gene'])
ggplot(melt(changes_in_CN_gene_df[,which(colnames(changes_in_CN_gene_df) %in% selected_chrom17)]),
       aes(x=Var2, y=value))+geom_boxplot()+geom_jitter()+theme_bw()+
  theme(axis.text.x=element_text(angle = 45, hjust = 1, vjust=1))
  

heatmap(changes_in_CN_gene_df)

##---------------
subset_genes_of_interest = c('MYC', 'CCNE1', 'PIK3CA', 'TERT', 'KRAS', 'PTEN', 'RB1', 'AKT1',
                             'AKT2', 'PARP1', 'PARP2', 'ATM', 'ATR', 'WEE1', 'TOP1', 'TUBB1',
                             'AKT3', 'CCND1', 'CCND2', 'CCND3', 'CDKN2A', 'CDKN2B', 'MECOM', 'CDK12')

head(all_cn[,1:2])

NA_for_zero <- function(i){
  i[is.na(i)] <- 0
  i
}


cors <- cor(t(all_cn[subset_genes_of_interest,]), use = "pairwise.complete.obs")

pheatmap::pheatmap(log(NA_for_zero(cors)+1))

all_cn[c('CDKN2A', 'CDKN2B'),1:10]

plot(log(sort(unlist(all_cn['MYC',]))))
abline(h=log(2))

all_cn['MYC',]

#------------

## add organoids segments

average_lower1p5 = colMeans(all_cn < 1.5)
anyecna = apply(all_cn, 2, function(i) any(i >= 8))
fractionecna = apply(all_cn, 2, function(i) mean(i >= 8, na.rm=T))
fraction2or1vicinity= apply(all_cn, 2, function(i) mean( ((i >= 0.8) & (i <= 1.2)) | ((i >= 1.8) & (i <= 2.2)), na.rm=T))
table(anyecna)

# ??? COPY NUMBER OF -10??
sort(unlist(all_cn))[1:10]
quantiles = apply(all_cn, 2, quantile, seq(0, 1, length.out=20), na.rm=T)
# quantiles <- quantiles[1:10,]
quantilespca <- prcomp(t(quantiles), scale. = T)
quantilestsne <- Rtsne::Rtsne(t(quantiles))
library(ggplot2)
library(ggrepel)
ggplot(cbind.data.frame(quantilestsne$Y),
       aes(x=`1`, y=`2`))+geom_point()+
  theme_bw()

ggplot(cbind.data.frame(quantilespca$x[,1:2]),
       aes(x=PC1, y=PC2))+geom_point()+
  theme_bw()

df_pca <- cbind.data.frame(quantilespca$x[,1:2], averagelow=average_lower1p5,
                           anyecna=anyecna, fractionecna=fractionecna,
                           fraction2or1vicinity=fraction2or1vicinity)
df_tsne <- cbind.data.frame(quantilestsne$Y[,1:2], averagelow=average_lower1p5,
                            anyecna=anyecna, fractionecna=fractionecna,
                            fraction2or1vicinity=fraction2or1vicinity)

brcaStatus = readRDS(file.path('../../../other_repos/cnsigs_Initial_submission/survival_analysis/from_ruben/survival_models/TCGA_OVBRCAonly_Exposures_and_BRCA_Status_plusGene.rds'))
brcaStatus$Sample <- as.character(brcaStatus$Sample)
brcaStatus$Status <- as.character(brcaStatus$Status)
brcaStatus$Gene <- as.character(brcaStatus$Gene)
brcaStatus$Signature <- as.character(brcaStatus$Signature)
df_tsne$brca = brcaStatus$Status[match(rownames(df_tsne), brcaStatus$Sample)]
df_tsne$BRCA2 = brcaStatus$Gene[match(rownames(df_tsne), brcaStatus$Sample)] == 'BRCA2'
df_tsne$Signature = brcaStatus$Signature[match(rownames(df_tsne), brcaStatus$Sample)]
df_tsne$brca_v2 = df_tsne$brca
df_tsne$brca_v2[df_tsne$brca_v2 %in% c('LOH in BRCA1/2') ] = NA
df_tsne$fraction_ecDNA = apply(all_cn, 2, function(i) mean(i>8))
df_tsne$MYC = unlist(all_cn['MYC',])

library(gridExtra)
library(reshape2)
grid.arrange(
  ggplot(df_pca,
         aes(x=PC1, y=PC2, col=averagelow))+geom_point()+
    theme_bw()+lims(x=c(-20, 50), y=c(-10, 10)),
  ggplot(df_pca,
         aes(x=PC1, y=PC2, col=anyecna))+geom_point()+
    theme_bw()+lims(x=c(-20, 50), y=c(-10, 10)),
  ggplot(df_pca,
         aes(x=PC1, y=PC2, col=log(fractionecna+0.001)))+geom_point()+
    theme_bw()+lims(x=c(-20, 50), y=c(-10, 10)),
  ggplot(df_pca,
         aes(x=PC1, y=PC2, col=log(fraction2or1vicinity)))+geom_point()+
    theme_bw()+lims(x=c(-20, 50), y=c(-10, 10)))#+lims(x=c(-5, 10), y=c(-5, 3)))

pdf("figures/tsne_brca_ecdna.pdf", width = 12, height = 2.5)
cowplot::plot_grid(
  # ggplot(df_tsne,
  #      aes(x=`1`, y=`2`, col=averagelow))+geom_point()+
  # theme_bw(),#+lims(x=c(-20, 50), y=c(-10, 10)),
  ggplot(df_tsne,
         aes(x=`1`, y=`2`, col=anyecna))+geom_point()+
    theme_bw()+#theme(legend.position = "bottom")+
    labs(col='pecDNA', x='TSNE1', y='TSNE2'),#+lims(x=c(-20, 50), y=c(-10, 10)),
  # ggplot(df_tsne,
  #        aes(x=`1`, y=`2`, col=log(fractionecna+0.001)))+geom_point()+
  #   theme_bw(),#+lims(x=c(-20, 50), y=c(-10, 10)),
  ggplot(df_tsne,
         aes(x=`1`, y=`2`, col=log(fraction2or1vicinity)))+geom_point()+
    theme_bw()+#theme(legend.position = "bottom")+
    labs(col='Fraction in CN={1,2}', x='TSNE1', y='TSNE2'),#+lims(x=c(-20, 50), y=c(-10, 10)),
  # ggplot(df_tsne,
  #        aes(x=`1`, y=`2`, col=factor(brca)))+geom_point()+
  #   theme_bw(),#+lims(x=c(-20, 50), y=c(-10, 10))
  ggplot()+
    geom_point(data=df_tsne[-which(df_tsne$BRCA2),],
               aes(x=`1`, y=`2`), alpha=0.2, col='grey')+
    geom_point(data=df_tsne[which(df_tsne$BRCA2),],
               aes(x=`1`, y=`2`, col=factor(BRCA2), shape=brca_v2))+
    theme_bw()+labs(col='BRCA2 mut', x='TSNE1', y='TSNE2')  # ggplot(df_tsne,
  #        aes(x=`1`, y=`2`, col=factor(brca_v2)))+geom_point()+
  #   theme_bw()+labs(col='BRCA1/2 status', x='TSNE1', y='TSNE2')#+theme(legend.position = "bottom")+
    #guides(col=guide_legend(nrow=5,byrow=TRUE))
  , ncol=3, rel_widths = c(2.7, 3,3.2)
)
dev.off()

ggplot(df_tsne,
       aes(x=`1`, y=`2`, col=cut(fraction_ecDNA, 10)))+geom_point()+
  theme_bw()+facet_wrap(.~cut(fraction_ecDNA, 10))
  # labs(col='Fraction in CN={1,2}', x='TSNE1', y='TSNE2')

ggplot(df_tsne,
       aes(x=`1`, y=`2`, col=log(MYC)))+geom_point()+
  theme_bw()+facet_wrap(.~cut(fraction_ecDNA, 10))
# labs(col='Fraction in CN={1,2}', x='TSNE1', y='TSNE2')

ggplot(df_tsne,
       aes(x=`1`, y=`2`, col=(MYC)>8))+geom_point()+
  theme_bw()+facet_wrap(.~cut(fraction_ecDNA, 10))

ggplot(df_tsne,
       aes(x=`1`, y=`2`, col=quantilestsne$Y[,1]))+geom_point()+
  theme_bw()+facet_wrap(.~cut(fraction_ecDNA, 10))


pairs(t(quantiles))

ggplot(df_tsne,
       aes(x=`1`, y=`2`, col=quantiles[9,]))+geom_point()+
  theme_bw()+facet_wrap(.~cut(fraction_ecDNA, 10))

ggplot(df_tsne,
       aes(x=`1`, y=`2`, col=Signature))+geom_point()+
  theme_bw()+#theme(legend.position = "bottom")+
  labs(col='Fraction in CN={1,2}', x='TSNE1', y='TSNE2')

ggplot(df_pca,
       aes(x=PC1, y=PC2, col=anyecna))+geom_point()+
  theme_bw()+lims(x=c(-20, 50), y=c(-10, 10))+facet_wrap(.~anyecna)

ggplot(melt(list(df_pca$averagelow[df_pca$anyecna], df_pca$averagelow[!df_pca$anyecna])),
       aes(x=L1, group=L1, y=value))+
  geom_boxplot()+geom_jitter()+theme_bw()


hclustquantiles <- hclust(dist(t(quantiles)))
# hclustquantiles$labels = rep('', length(hclustquantiles$labels))
plot(hclustquantiles)
cutree_quantiles <- cutree(hclustquantiles, 10)

df_pca$cutree_quantiles = cutree_quantiles[match(rownames(df_pca), names(cutree_quantiles))]

ggplot(df_pca,
       aes(x=PC1, y=PC2, col=factor(cutree_quantiles)))+geom_point()+
  theme_bw()+lims(x=c(-20, 50), y=c(-10, 10))

var_cn = apply(all_cn, 2, var)
fraction1vicinity= apply(all_cn, 2, function(i) mean( (i >= 0.8) & (i <= 1.2), na.rm=T))

ggplot(df_tsne,
       aes(x=`1`, y=`2`,
           # col=quantiles[9,]
           # col=log(var_cn)
           col=fraction1vicinity
           ))+geom_point()+
  theme_bw()+#theme(legend.position = "bottom")+
  labs(col='Fraction in CN=1', x='TSNE1', y='TSNE2')

ggplot(melt(quantiles), aes(x=Var1, y=value, group=Var2))+geom_line(alpha=0.2)+
  theme_bw()#+scale_y_continuous(trans = "log2")
ggplot(melt(quantiles), aes(x=Var1, y=value, group=Var2))+geom_line(alpha=0.2, col='blue')+
  theme_bw()+lims(y=c(-5, 10))
