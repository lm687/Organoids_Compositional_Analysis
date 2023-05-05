rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(reshape2)
library(ggplot2)
library(ggdendro)
library(pheatmap)
library(readxl)
library(gridExtra)
library(ggplotify)
library(cowplot)
library(ggh4x) ## nested facets
library(dplyr)
library(EnvStats)
library(GenomicRanges)

source("helper.R")

load_ag_genes <- function(){
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
  # EnsDb.Hsapiens.v87 <- ensembldb::EnsDb(sqlite_path)
  EnsDb.Hsapiens.v87 <- ensembldb::EnsDb(sqlite_file)
  # Genes, used to annotated the TPM matrix to send to Maria
  ag <- ensembldb::genes(EnsDb.Hsapiens.v87, filter=list(AnnotationFilter::GeneBiotypeFilter('protein_coding')), return.type="DataFrame") 
  ag_subsetchrom <- ag[!(ag$seq_name %in% c("MT", "X", "Y")) & !grepl("GL",ag$seq_name),]
  return(ag_subsetchrom)
}

load_chromlens <- function(){
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
  # EnsDb.Hsapiens.v87 <- ensembldb::EnsDb(sqlite_path)
  EnsDb.Hsapiens.v87 <- ensembldb::EnsDb(sqlite_file)
  chromosomes_metadata <- ensembldb::seqinfo(EnsDb.Hsapiens.v87)
  chromlens <- (cbind.data.frame(Chrom=seqnames(chromosomes_metadata),
                                 Length=seqlengths(chromosomes_metadata)))
  return(chromlens)
}

cached_ag_subsetchrom <- T
if(!cached_ag_subsetchrom){
  ag_subsetchrom <- load_ag_genes()
  saveRDS(ag_subsetchrom, "../robjects/ag_subsetchrom.RDS")
}else{
  ag_subsetchrom <- readRDS("../robjects/ag_subsetchrom.RDS")
}
##----------------------------------------------------------------------------
##----------------------------------------------------------------------------

renaming <- readxl::read_excel("../../RNASeq_DE_resistant_sensitive/files/PDOnameProperSample_sWGS_RNAseq.xlsx")

organoid_list = c('118976org', '119148orgb', '23868org')

read_UID_orgs <- function(org){
  org_tab = read.csv(paste0("../UID-", org, ".csv"))
  org_tab
  org_clean = org_tab[org_tab[,1:4][,'num_noisy'] == 0,]
  org_clean
} ## what was done to create the absCN_clean_ files

read_UID_orgs_raw <- function(org){
  org_tab = read.csv(paste0("../UID-", org, ".csv"))
  org_tab
}

absCN0 = lapply(paste0("../data/absCN_clean_", organoid_list, ".RDS"), readRDS)
absCN = lapply(paste0("../data/absCN_clean_", organoid_list, ".RDS"), readRDS)
names(absCN) = organoid_list
names(absCN) = gsub("orgb", "org", names(absCN))
names(absCN) = renaming$PDO[match(names(absCN), renaming$ID)]
saveRDS(absCN, "../robjects/fig2_absCN.RDS")
organoid_list = renaming$PDO[match(gsub("orgb", "org", organoid_list), renaming$ID)]


# image(absCN[[1]])
# 
# plot(hclust(dist(absCN[[org_it]])))

# org_it = '23868org'

## remove the outliers
outliers = list()
outliers$`PDO3` = c(12, 10, 11, 7, 5, 6, 4, 3, 1, 2, 8, 9, 158)
outliers$`PDO6` = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 23, 377, 378, 379, 380, 381)
outliers$`PDO2` = c(1,2,3, 4, 5, 6, 7, 8)

saveRDS(outliers, "../robjects/scDNA_outliers.RDS")

# plot(hclust(dist(absCN[[org_it]][-outliers[[org_it]],])))

for(org_it in organoid_list){
  org_it_id <- renaming$ID[match(org_it, renaming$PDO)]
  if(org_it == 'PDO6'){
    org_it_id <- "119148orgb"
  }
  rownames(absCN[[org_it]]) <- read_UID_orgs(org_it_id)$node_id
  absCN[[org_it]] = absCN[[org_it]][-outliers[[org_it]],]
  
  absCN[[org_it]] = absCN[[org_it]][,rowSums(apply(absCN[[org_it]], 1, is.na)) == 0]
  
  ## select only top variable areas of the genome
  # absCN[[org_it]] = absCN[[org_it]][,order(apply(absCN[[org_it]] , 2, var), decreasing = T)[1:1000]]
  ## doesn't work too well
  
  absCN[[org_it]][absCN[[org_it]] > 14] = 14
  
  mat_breaks <- c(-Inf, 0:6, seq(8, 14, by=2)) - 0.01 ## adding a small number because otherwise the binning is done wrong
  col_list <- c("#2670af", "#2670af", "#8daecf", "#eaeaea", "#ffd2b6", "#f5b788", "#f39d5f", "#f17e37", "#ee1200", "#b60000", "#7a0e04", "#401004", "#000000")
  
  # mat_breaks <- c(-Inf, 0:6, seq(8, 14, by=2)) #seq(min(org_clean, na.rm = T), max(org_clean, na.rm = T), length.out = 10)
  # col_list <- c("#2670af", "#8daecf", "#eaeaea", "#ffd2b6", "#f5b788", "#f39d5f", "#f17e37", "#ee1200", "#b60000", "#7a0e04", "#401004", "#000000")
  annotation_chroms = data.frame(row.names = colnames(absCN[[org_it]]),
                                 chrom=clean_chrom(gsub("\\..*", "", colnames(absCN[[org_it]]))))
  
  ph = pheatmap::pheatmap(absCN[[org_it]]+0.02, show_colnames = FALSE, show_rownames = FALSE,
                     color             = col_list,
                     breaks            = mat_breaks, cluster_cols = FALSE,
                     annotation_col = annotation_chroms, annotation_legend=F, main=org_it)
  
  ph$gtable$grobs[[6]]$children[[1]]$vjust = 0.45
  sexchrom_bool = annotation_chroms$chrom %in% c('X', 'Y')
  annotation_chroms_nosexchrom <- data.frame(row.names=rownames(annotation_chroms)[!sexchrom_bool],
                                             chrom=annotation_chroms$chrom[!sexchrom_bool])
  ph_nosexchrom = pheatmap::pheatmap(absCN[[org_it]][,!sexchrom_bool]+0.02, show_colnames = FALSE, show_rownames = FALSE,
                          color             = col_list,
                          breaks            = mat_breaks, cluster_cols = FALSE,
                          annotation_col = annotation_chroms_nosexchrom, annotation_legend=F, main=org_it)
  ph_nosexchrom$gtable$grobs[[6]]$children[[1]]$vjust = 0.45
  saveRDS(list(mat_breaks=mat_breaks,col_list=col_list), paste0("../robjects/fig2_colours.RDS"))
  saveRDS(absCN[[org_it]], paste0("../robjects/fig2_subclonal_hclust", org_it, "_2.RDS"))
  saveRDS(absCN[[org_it]][,!sexchrom_bool], paste0("../robjects/fig2_subclonal_hclust_nosexchrom", org_it, "_2.RDS"))
  saveRDS(ph, paste0("../robjects/fig2_subclonal_hclust", org_it, ".RDS"))
  saveRDS(ph_nosexchrom, paste0("../robjects/fig2_subclonal_hclust_nosexchrom", org_it, ".RDS"))
  # ph$gtable$grobs[[1]]$gp <- gpar(lwd = 5)
  # ph$gtable$grobs[[2]]$gp <- gpar(col = 'blue')
  
  plot(as.dendrogram(ph$tree_row))

  pdf(paste0("../plots/subclonal_hclust", org_it, ".pdf"))
  print(ph)
  dev.off()
  
}

ggplot(melt(t(absCN[[org_it]])), aes(y=Var1, x=factor(Var2, levels=ph$tree_row$order), fill=value))+
  geom_tile()+
  scale_fill_gradientn(colours=col_list, breaks=mat_breaks) ## alternative plotting

data_orgs <- lapply(organoid_list, function(org_it){
  .x <- readRDS(paste0("../robjects/fig2_subclonal_hclust", org_it, "_2.RDS"))
  ## remove sex chroms
  chrom=clean_chrom(gsub("\\..*", "", colnames(.x)))
  .x[,!(chrom %in% c('X', 'Y'))]
})
names(data_orgs) <- organoid_list

#-------------------------------------------------------------------------#
### general structure of the populations
hclusts_orgs <- lapply(names(absCN), function(nm) hclust(dist((as(absCN[[nm]][,!sexchrom_bool], 'matrix')))))
names(hclusts_orgs) <- names(absCN)

hm <- list()
hm2 <- list()
for(org_it in names(hclusts_orgs)){
  hm[[org_it]] <- pheatmap(absCN[[org_it]][,!sexchrom_bool]+0.02,  show_colnames = FALSE, show_rownames = FALSE,
           color             = col_list,
           breaks            = mat_breaks, cluster_cols = FALSE,
           annotation_col = annotation_chroms_nosexchrom, annotation_legend=F, main=org_it)
  hm2[[org_it]] <- pheatmap(absCN[[org_it]][,!sexchrom_bool],  show_colnames = FALSE, show_rownames = FALSE,
                 color             = col_list,
                 breaks            = mat_breaks, cluster_cols = FALSE,
                 annotation_col = annotation_chroms_nosexchrom, annotation_legend=F, main=org_it)
  saveRDS(hm2[[org_it]], paste0("../robjects/pheatmap_", org_it, ".RDS"))
}

pdf(paste0("../plots/subclonal_hclust_dendrograms.pdf"), width = 5, height = 2.5)
sapply(names(hclusts_orgs), function(org_it){
  print(cowplot::plot_grid(hm[[org_it]][[4]], 
                   hm2[[org_it]][[4]],
                   ggdendrogram(hm2[[org_it]][[1]])+theme(axis.title.x=element_blank(),
                                                axis.text.x=element_blank(),
                                                axis.ticks.x=element_blank())+ coord_flip(),
                   ggdendrogram(hclusts_orgs[[org_it]])+theme(axis.title.x=element_blank(),
                                                              axis.text.x=element_blank(),
                                                              axis.ticks.x=element_blank())+ coord_flip(), ncol = 4 ))
})
dev.off()
hclusts_orgs[[org_it]]$labels

## analysis of clusters

## PDO2
plot(hm2[[org_it]]$tree_row)
kclades_pdo2 <- 4
table(cutree(hm2$PDO2$tree_row, k = kclades_pdo2))
## two major clades, of 42 and 30 cells
hm2lab <- hm2
hm2lab$PDO2$tree_row$labels = cutree(hm2[['PDO2']]$tree_row, k = kclades_pdo2)
hm2lab$PDO2$tree_row$labels[!duplicated(hm2lab$PDO2$tree_row$labels)] = paste0(hm2lab$PDO2$tree_row$labels[!duplicated(hm2lab$PDO2$tree_row$labels)], ' Clade')

saveRDS(dendextend::color_branches(as.dendrogram(hm2lab$PDO2$tree_row),
                                   clusters = cutree((hm2$PDO2$tree_row), k = kclades_pdo2)[labels(as.dendrogram(hm2$PDO2$tree_row))]),
        "../robjects/dendextend_clades_PDO2.RDS")
pdf("../plots/dendrogram_scDNA_PDO2.pdf")
par(mfrow=c(1,1), mar=c(2,0,0,0))
plot(dendextend::color_branches(as.dendrogram(hm2lab$PDO2$tree_row),
                                clusters = cutree((hm2$PDO2$tree_row), k = kclades_pdo2)[labels(as.dendrogram(hm2$PDO2$tree_row))]))
dev.off()
remove_labels <- function(i){
  labels(i) <- NULL
  i
}

pdf("../plots/dendrogram_scDNA_PDO2_v2.pdf", height = 2, width = 2)
par(mfrow=c(1,1), mar=c(0,0,0,0))
plot(remove_labels(dendextend::color_branches(as.dendrogram((hm2lab$PDO2$tree_row)),
  clusters = cutree((hm2$PDO2$tree_row), k = kclades_pdo2)[labels(as.dendrogram(hm2$PDO2$tree_row))])))
text(24,16, 'A')
text(60,16, 'B')
dev.off()

.xxxxx <- dendextend::remove_nodes_nodePar(dendextend::color_branches(as.dendrogram((hm2lab$PDO2$tree_row)),
                                                                      clusters = cutree((hm2$PDO2$tree_row), k = kclades_pdo2)[labels(as.dendrogram(hm2$PDO2$tree_row))],
                                                                      labels=F))
labels(.xxxxx) <- NULL

library(ggdendro)
ggdendro::ggdendrogram(hm2lab$PDO2$tree_row)

dendro_data(hm2lab$PDO2$tree_row, )

# ggplot(melt(list(absCN[[org_it]][,!sexchrom_bool][cutree(hm2[[org_it]]$tree_row, k = 4) == 1,],
#      absCN[[org_it]][,!sexchrom_bool][cutree(hm2[[org_it]]$tree_row, k = 4) == 4,])),
#      aes(x=Var2, y=value, group=interaction(Var2, L1)))+
#   geom_boxplot()

absCN_list_PDO2 <- list(absCN$PDO2[,!sexchrom_bool][cutree(hm2$PDO2$tree_row, k = 4) == 1,],
                        absCN$PDO2[,!sexchrom_bool][cutree(hm2$PDO2$tree_row, k = 4) == 4,])

write.table(absCN_list_PDO2[[1]], "../robjects/table_absCNscDNA_PDO2_cladeA_1.txt", sep = "\t", quote = F, row.names = F)
write.table(absCN_list_PDO2[[2]], "../robjects/table_absCNscDNA_PDO2_cladeB_4.txt", sep = "\t", quote = F, row.names = F)

absCN_df_clades_PDO2 <- cbind.data.frame(clade1=melt(absCN_list_PDO2[[1]]) %>% group_by(Var2) %>% dplyr::summarise(mean=mean(value), sd=sd(value), quant005=quantile(value, probs = 0.05), quant95=quantile(value, probs = 0.95)),
                 clade4=melt(absCN_list_PDO2[[2]]) %>% group_by(Var2) %>% dplyr::summarise(mean=mean(value), sd=sd(value), quant005=quantile(value, probs = 0.05), quant95=quantile(value, probs = 0.95)))
stopifnot(absCN_df_clades_PDO2$clade1.Var2 == absCN_df_clades_PDO2$clade4.Var2)
absCN_df_clades_PDO2$chrom <- clean_chrom(gsub("\\..*", "", (absCN_df_clades_PDO2$clade1.Var2)))
absCN_df_clades_PDO2$confint <- (absCN_df_clades_PDO2$clade1.quant005 > absCN_df_clades_PDO2$clade4.quant95) | absCN_df_clades_PDO2$clade4.quant005 > absCN_df_clades_PDO2$clade1.quant95
ggplot(absCN_df_clades_PDO2,
      aes(x=clade1.mean, y=clade4.mean))+geom_point(aes(col=confint))+
  geom_errorbar(aes(xmin=clade1.mean-clade1.sd, xmax=clade1.mean+clade1.sd,
                    ymin=clade4.mean-clade4.sd, ymax=clade4.mean+clade4.sd))+
  geom_abline(slope = 1, intercept = 0, lty='dashed')+theme_bw()+ggtitle('Clade differences in PDO2')+
labs(x='CN mean across cells in Clade 1', y='CN mean across cells in Clade 2')+
  facet_wrap(.~factor(as.numeric(chrom)))+guides(col=FALSE)
ggsave(paste0("../plots/subclonal_clades_cells_summarised_chromosomes_PDO2.pdf"), width = 5, height = 5)

pvals_clade1_4 <- as.numeric(sapply(1:ncol(absCN_list_PDO2[[1]]), function(idx_col) try(t.test(absCN_list_PDO2[[1]][,idx_col],
                                                                                               absCN_list_PDO2[[2]][,idx_col])$p.value)))
absCN_df_clades_PDO2$clade1.Var2 <- factor(absCN_df_clades_PDO2$clade1.Var2, levels=absCN_df_clades_PDO2$clade1.Var2[order(pvals_clade1_4)])
absCN_df_clades_PDO2$pvalstat1_4 <- p.adjust(pvals_clade1_4) < 0.05

ggplot(absCN_df_clades_PDO2[order(pvals_clade1_4),],
       aes(x=clade1.Var2, y= clade1.mean-clade4.mean, col=factor(ifelse(pvalstat1_4, as.numeric(chrom), NA))))+
  geom_point()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank())+geom_abline(slope = 0, intercept = 0, lty='dashed')+
  theme(legend.position = "bottom")+labs(col='Chromosome')
ggsave(paste0("../plots/subclonal_clades_cells_PDO2_pvals_clade1_4.pdf"), width = 5, height = 5)

ggplot(absCN_df_clades_PDO2[order(pvals_clade1_4),],
       aes(x=clade1.Var2, y= clade1.mean-clade4.mean, col=factor(ifelse(pvalstat1_4, as.numeric(chrom), NA))))+
  geom_point()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank())+geom_abline(slope = 0, intercept = 0, lty='dashed')+
  theme(legend.position = "bottom")+labs(col='Chromosome')+facet_wrap(.~as.numeric(chrom))
ggsave(paste0("../plots/subclonal_clades_cells_PDO2_pvals_clade1_4_facet.pdf"), width = 10, height = 10)

add_names <- function(i,j){
  names(i) <- j
  i
}
ggplot(absCN_df_clades_PDO2,
       aes(x=clade1.Var2, y= clade1.mean-clade4.mean, col=factor(ifelse(pvalstat1_4, as.numeric(chrom), NA))))+
  geom_bar(stat = "identity")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        strip.background = element_blank(), strip.text = element_blank())+geom_abline(slope = 0, intercept = 0, lty='dashed')+
  theme(legend.position = "bottom")+labs(col='Chromosome')+
  facet_grid (.~ as.numeric(chrom), scales = "free_x", space = "free_x")+
  guides(col=F)+labs(y='PDO2\nCN (A) - CN (B)')+
  scale_color_manual(values = add_names(unique(ph$gtable$grobs[[4]]$gp$fill[,1])[1:22], as.character(1:22)))
ggsave(paste0("../plots/subclonal_clades_cells_PDO2_pvals_clade1_4_facet_bar.pdf"), width = 8, height = 1.5)


onerow_dummy <- absCN_df_clades_PDO2[sapply(unique(absCN_df_clades_PDO2$chrom), function(i) which(absCN_df_clades_PDO2$chrom == i)[1]),]
ggplot(absCN_df_clades_PDO2[absCN_df_clades_PDO2$pvalstat1_4,],
       aes(x=clade1.mean, y= clade4.mean))+
  geom_rect(data=onerow_dummy, mapping=aes(xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=2), color=NA, fill="blue", alpha=0.3)+
  geom_rect(data=onerow_dummy, mapping=aes(xmin=-Inf, xmax=Inf, ymin=2, ymax=Inf), color=NA, fill="yellow", alpha=0.3)+
  geom_point()+
  geom_errorbar(aes(xmin=clade1.mean-clade1.sd, xmax=clade1.mean+clade1.sd,
                    ymin=clade4.mean-clade4.sd, ymax=clade4.mean+clade4.sd))+
  geom_abline(slope = 1, intercept = 0, lty='dashed')+theme_bw()+
  geom_abline(slope = 0, intercept = 0, lty='dashed')+
  theme(legend.position = "bottom")+labs(col='Chromosome')+facet_wrap(.~as.numeric(chrom), scales = "free")
ggsave(paste0("../plots/subclonal_clades_cells_summarised_chromosomes_PDO2_DA_regions.pdf"), width = 10, height = 10)

## PDO3
plot(hm2$PDO3$tree_row)
kclades_pdo3 <- 7
table(cutree(hm2[['PDO3']]$tree_row, k = kclades_pdo3))

plot(dendextend::color_branches(as.dendrogram(hm2$PDO3$tree_row),
                                clusters = cutree(as.dendrogram(hm2$PDO3$tree_row), k = kclades_pdo3)[labels(as.dendrogram(hm2$PDO3$tree_row))]))

hm2lab$PDO3$tree_row$labels = cutree(hm2[['PDO3']]$tree_row, k = kclades_pdo3)
hm2lab$PDO3$tree_row$labels[!duplicated(hm2lab$PDO3$tree_row$labels)] = paste0(hm2lab$PDO3$tree_row$labels[!duplicated(hm2lab$PDO3$tree_row$labels)], ' Clade')

saveRDS(dendextend::color_branches(as.dendrogram(hm2lab$PDO3$tree_row),
                                   clusters = cutree((hm2$PDO3$tree_row), k = kclades_pdo3)[labels(as.dendrogram(hm2$PDO3$tree_row))]),
        "../robjects/dendextend_clades_PDO3.RDS")

pdf("../plots/dendrogram_scDNA_PDO3.pdf")
par(mfrow=c(1,1), mar=c(3,0,0,0))
plot(dendextend::color_branches(as.dendrogram(hm2lab$PDO3$tree_row),
                                clusters = cutree((hm2$PDO3$tree_row), k = kclades_pdo3)[labels(as.dendrogram(hm2$PDO3$tree_row))]))
dev.off()

pdf("../plots/dendrogram_scDNA_PDO3_v2.pdf", height = 2, width = 2)
par(mfrow=c(1,1), mar=c(0,0,0,0))
plot(remove_labels(dendextend::color_branches(as.dendrogram(hm2lab$PDO3$tree_row),
                                              clusters = cutree((hm2$PDO3$tree_row), k = kclades_pdo3)[labels(as.dendrogram(hm2$PDO3$tree_row))])))
text(24,18, 'B')
text(78,18, 'A')
text(135,18, 'C')
dev.off()

## three major clades, labelled 3 (40 cells), 4 (52 cells) and 6 (48 cells)
absCN_list_PDO3 <- list(absCN$PDO3[,!sexchrom_bool][cutree(hm2$PDO3$tree_row, k = kclades_pdo3) == 3,],
                        absCN$PDO3[,!sexchrom_bool][cutree(hm2$PDO3$tree_row, k = kclades_pdo3) == 4,],
                        absCN$PDO3[,!sexchrom_bool][cutree(hm2$PDO3$tree_row, k = kclades_pdo3) == 6,])

write.table(absCN_list_PDO3[[1]], "../robjects/table_absCNscDNA_PDO3_cladeA_3.txt", sep = "\t", quote = F, row.names = F)
write.table(absCN_list_PDO3[[2]], "../robjects/table_absCNscDNA_PDO3_cladeB_4.txt", sep = "\t", quote = F, row.names = F)
write.table(absCN_list_PDO3[[3]], "../robjects/table_absCNscDNA_PDO3_cladeC_6.txt", sep = "\t", quote = F, row.names = F)

sapply(absCN_list_PDO3, dim)[1,]
absCN_df_clades_PDO3 <- cbind.data.frame(do.call('cbind', lapply(absCN_list_PDO3, function(i) melt(i) %>% group_by(Var2) %>% dplyr::summarise(mean=mean(value), sd=sd(value), quant005=quantile(value, probs = 0.05), quant95=quantile(value, probs = 0.95)))))
colnames(absCN_df_clades_PDO3) <- paste0(rep(c('clade3.', 'clade4.', 'clade6.'), each=ncol(absCN_df_clades_PDO3)/3), colnames(absCN_df_clades_PDO3))
absCN_df_clades_PDO3$chrom <- clean_chrom(gsub("\\..*", "", (absCN_df_clades_PDO3$clade3.Var2)))

pdf(paste0("../plots/subclonal_clades_cells_summarised_chromosomes_PDO3.pdf"), width = 5, height = 10)
grid.arrange(ggplot(absCN_df_clades_PDO3,
                    aes(x=clade3.mean, y=clade4.mean))+geom_point(aes(col=1))+
               geom_errorbar(aes(xmin=clade3.mean-clade3.sd, xmax=clade3.mean+clade3.sd,
                                 ymin=clade4.mean-clade4.sd, ymax=clade4.mean+clade4.sd))+
               geom_abline(slope = 1, intercept = 0, lty='dashed')+theme_bw()+ggtitle('Clade differences in PDO3')+
               labs(x='CN mean across cells in Clade 3', y='CN mean across cells in Clade 4')+
               facet_wrap(.~factor(as.numeric(chrom)))+guides(col=FALSE),
             ggplot(absCN_df_clades_PDO3,
                    aes(x=clade3.mean, y=clade6.mean))+geom_point(aes(col=1))+
               geom_errorbar(aes(xmin=clade3.mean-clade3.sd, xmax=clade3.mean+clade3.sd,
                                 ymin=clade6.mean-clade6.sd, ymax=clade6.mean+clade6.sd))+
               geom_abline(slope = 1, intercept = 0, lty='dashed')+theme_bw()+ggtitle('Clade differences in PDO3')+
               labs(x='CN mean across cells in Clade 3', y='CN mean across cells in Clade 6')+
               facet_wrap(.~factor(as.numeric(chrom)))+guides(col=FALSE)
)
dev.off()

pvals_clade3_4 <- as.numeric(sapply(1:ncol(absCN_list_PDO3[[1]]), function(idx_col) try(t.test(absCN_list_PDO3[[1]][,idx_col],
                                                                   absCN_list_PDO3[[2]][,idx_col])$p.value)))
absCN_df_clades_PDO3$clade3.Var2 <- factor(absCN_df_clades_PDO3$clade3.Var2, levels=absCN_df_clades_PDO3$clade3.Var2[order(pvals_clade3_4)])
absCN_df_clades_PDO3$pvalstat3_4 <- p.adjust(pvals_clade3_4) < 0.05

ggplot(absCN_df_clades_PDO3[order(pvals_clade3_4),],
       aes(x=clade3.Var2, y= clade3.mean-clade4.mean, col=factor(ifelse(pvalstat3_4, as.numeric(chrom), NA))))+
  geom_point()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank())+geom_abline(slope = 0, intercept = 0, lty='dashed')+
  theme(legend.position = "bottom")+labs(col='Chromosome')
ggsave(paste0("../plots/subclonal_clades_cells_PDO3_pvals_clade3_4.pdf"), width = 5, height = 5)

ggplot(absCN_df_clades_PDO3[order(pvals_clade3_4),],
       aes(x=clade3.Var2, y= clade3.mean-clade4.mean, col=factor(ifelse(pvalstat3_4, as.numeric(chrom), NA))))+
  geom_point()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank())+geom_abline(slope = 0, intercept = 0, lty='dashed')+
  theme(legend.position = "bottom")+labs(col='Chromosome')+facet_wrap(.~as.numeric(chrom))
ggsave(paste0("../plots/subclonal_clades_cells_PDO3_pvals_clade3_4_facet.pdf"), width = 10, height = 10)

chrom_to_num <- function(i){
  i$chrom <- as.numeric(i$chrom)
  i
}
one_to_22_chrom <- absCN_df_clades_PDO3$chrom %in% as.character(1:22)
ggplot(chrom_to_num(absCN_df_clades_PDO3[one_to_22_chrom,]),
       aes(x=clade3.Var2, y= clade3.mean-clade4.mean, col=factor(ifelse(pvalstat3_4, as.numeric(chrom), NA))))+
  geom_bar(stat = "identity")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        strip.background = element_blank(), strip.text = element_blank())+geom_abline(slope = 0, intercept = 0, lty='dashed')+
  theme(legend.position = "bottom")+labs(col='Chromosome')+
  facet_grid (.~ (chrom), scales = "free_x", space = "free_x")+
  guides(col=F)+labs(y='PDO3\nCN (A) - CN (B)')+
  scale_color_manual(values = add_names(unique(ph$gtable$grobs[[4]]$gp$fill[,1])[1:22], as.character(1:22)))
ggsave(paste0("../plots/subclonal_clades_cells_PDO3_pvals_clade3_4_facet_bar.pdf"), width = 8, height = 1.5)


onerow_dummy2 <- absCN_df_clades_PDO3[sapply(unique(absCN_df_clades_PDO3$chrom), function(i) which(absCN_df_clades_PDO3$chrom == i)[1]),]
ggplot(absCN_df_clades_PDO3[absCN_df_clades_PDO3$pvalstat3_4,],
       aes(x=clade3.mean, y= clade4.mean))+
  geom_rect(data=onerow_dummy2, mapping=aes(xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=2), color=NA, fill="blue", alpha=0.3)+
  geom_rect(data=onerow_dummy2, mapping=aes(xmin=-Inf, xmax=Inf, ymin=2, ymax=Inf), color=NA, fill="yellow", alpha=0.3)+
  geom_point()+
  geom_errorbar(aes(xmin=clade3.mean-clade3.sd, xmax=clade3.mean+clade3.sd,
                    ymin=clade4.mean-clade4.sd, ymax=clade4.mean+clade4.sd))+
  geom_abline(slope = 1, intercept = 0, lty='dashed')+theme_bw()+
  geom_abline(slope = 0, intercept = 0, lty='dashed')+
  theme(legend.position = "bottom")+labs(col='Chromosome')+facet_wrap(.~as.numeric(chrom), scales = "free")
ggsave(paste0("../plots/subclonal_clades_cells_summarised_chromosomes_PDO3_DA_regions_3v4.pdf"), width = 10, height = 10)

pvals_clade3_6 <- as.numeric(sapply(1:ncol(absCN_list_PDO3[[1]]), function(idx_col) try(t.test(absCN_list_PDO3[[1]][,idx_col],
                                                                                               absCN_list_PDO3[[3]][,idx_col])$p.value)))
absCN_df_clades_PDO3$clade3.Var2 <- factor(absCN_df_clades_PDO3$clade3.Var2, levels=absCN_df_clades_PDO3$clade3.Var2[order(pvals_clade3_6)])
absCN_df_clades_PDO3$pvalstat3_6 <- p.adjust(pvals_clade3_6) < 0.05

ggplot(absCN_df_clades_PDO3[order(pvals_clade3_6),],
       aes(x=clade3.Var2, y= clade3.mean-clade6.mean, col=factor(ifelse(pvalstat3_6, as.numeric(chrom), NA))))+
  geom_point()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank())+geom_abline(slope = 0, intercept = 0, lty='dashed')+
  theme(legend.position = "bottom")+labs(col='Chromosome')
ggsave(paste0("../plots/subclonal_clades_cells_PDO3_pvals_clade3_6.pdf"), width = 5, height = 5)

ggplot(absCN_df_clades_PDO3[order(pvals_clade3_6),],
       aes(x=clade3.Var2, y= clade3.mean-clade6.mean, col=factor(ifelse(pvalstat3_6, as.numeric(chrom), NA))))+
  geom_point()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank())+geom_abline(slope = 0, intercept = 0, lty='dashed')+
  theme(legend.position = "bottom")+labs(col='Chromosome')+facet_wrap(.~as.numeric(chrom))
ggsave(paste0("../plots/subclonal_clades_cells_PDO3_pvals_clade3_6_facet.pdf"), width = 10, height = 10)


ggplot(chrom_to_num(absCN_df_clades_PDO3[one_to_22_chrom,]),
       aes(x=clade3.Var2, y= clade3.mean-clade6.mean, col=factor(ifelse(pvalstat3_6, as.numeric(chrom), NA))))+
  geom_bar(stat = "identity")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+geom_abline(slope = 0, intercept = 0, lty='dashed')+
  theme(legend.position = "bottom",
        strip.background = element_blank(), strip.text = element_blank())+labs(col='Chromosome')+
  facet_grid (.~ (chrom), scales = "free_x", space = "free_x")+
  guides(col=F)+labs(y='PDO3\nCN (A) - CN (C)')+
  scale_color_manual(values = add_names(unique(ph$gtable$grobs[[4]]$gp$fill[,1])[1:22], as.character(1:22)))
ggsave(paste0("../plots/subclonal_clades_cells_PDO3_pvals_clade3_6_facet_bar.pdf"), width = 8, height = 1.5)


ggplot(absCN_df_clades_PDO3[absCN_df_clades_PDO3$pvalstat3_6,],
       aes(x=clade3.mean, y= clade6.mean))+
  geom_rect(data=onerow_dummy2, mapping=aes(xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=2), color=NA, fill="blue", alpha=0.3)+
  geom_rect(data=onerow_dummy2, mapping=aes(xmin=-Inf, xmax=Inf, ymin=2, ymax=Inf), color=NA, fill="yellow", alpha=0.3)+
  geom_point()+
  geom_errorbar(aes(xmin=clade3.mean-clade3.sd, xmax=clade3.mean+clade3.sd,
                    ymin=clade6.mean-clade6.sd, ymax=clade6.mean+clade6.sd))+
  geom_abline(slope = 1, intercept = 0, lty='dashed')+theme_bw()+
  geom_abline(slope = 0, intercept = 0, lty='dashed')+
  theme(legend.position = "bottom")+labs(col='Chromosome')+facet_wrap(.~as.numeric(chrom), scales = "free")
ggsave(paste0("../plots/subclonal_clades_cells_summarised_chromosomes_PDO3_DA_regions_3v6.pdf"), width = 10, height = 10)


## PDO6
plot(hm2$PDO6$tree_row)
kclades_pdo6 <- 6
table(cutree(hm2[['PDO6']]$tree_row, k = kclades_pdo6))

hm2lab$PDO6$tree_row$labels = cutree(hm2[['PDO6']]$tree_row, k = kclades_pdo6)
# hm2lab$PDO6$tree_row$labels[duplicated(hm2lab$PDO6$tree_row$labels)] = 'hee'
hm2lab$PDO6$tree_row$labels[!duplicated(hm2lab$PDO6$tree_row$labels)] = paste0(hm2lab$PDO6$tree_row$labels[!duplicated(hm2lab$PDO6$tree_row$labels)], ' Clade')

saveRDS(dendextend::color_branches(as.dendrogram(hm2lab$PDO6$tree_row),
                                   clusters = cutree((hm2$PDO6$tree_row), k = kclades_pdo6)[labels(as.dendrogram(hm2$PDO6$tree_row))]),
        "../robjects/dendextend_clades_PDO6.RDS")

pdf("../plots/dendrogram_scDNA_PDO6.pdf")
par(mfrow=c(1,1), mar=c(2,0,0,0))
plot(dendextend::color_branches(as.dendrogram(hm2lab$PDO6$tree_row),
                                clusters = cutree((hm2$PDO6$tree_row), k = kclades_pdo6)[labels(as.dendrogram(hm2$PDO6$tree_row))]))
dev.off()
plot(as.dendrogram(hm2$PDO6$tree_row))

pdf("../plots/dendrogram_scDNA_PDO6_v2.pdf", height = 2, width = 2)
par(mfrow=c(1,1), mar=c(0,0,0,0))
plot(remove_labels(dendextend::color_branches(as.dendrogram(hm2lab$PDO6$tree_row),
                                              clusters = cutree((hm2$PDO6$tree_row), k = kclades_pdo6)[labels(as.dendrogram(hm2$PDO6$tree_row))])))
text(30,25, 'C')
text(160,20, 'A')
text(270,25, 'B')
dev.off()

## group 1 (158 cells), 2 (145 cells) and 3 (49 cells)
absCN_list_PDO6 <- list(absCN$PDO6[,!sexchrom_bool][cutree(hm2$PDO6$tree_row, k = kclades_pdo6) == 1,],
                        absCN$PDO6[,!sexchrom_bool][cutree(hm2$PDO6$tree_row, k = kclades_pdo6) == 2,],
                        absCN$PDO6[,!sexchrom_bool][cutree(hm2$PDO6$tree_row, k = kclades_pdo6) == 3,])

write.table(absCN_list_PDO6[[1]], "../robjects/table_absCNscDNA_PDO6_cladeA_1.txt", sep = "\t", quote = F, row.names = F)
write.table(absCN_list_PDO6[[2]], "../robjects/table_absCNscDNA_PDO6_cladeB_2.txt", sep = "\t", quote = F, row.names = F)
write.table(absCN_list_PDO6[[3]], "../robjects/table_absCNscDNA_PDO6_cladeC_3.txt", sep = "\t", quote = F, row.names = F)


absCN_df_clades_PDO6 <- cbind.data.frame(do.call('cbind', lapply(absCN_list_PDO6, function(i) melt(i) %>% group_by(Var2) %>%
          dplyr::summarise(mean=mean(value), sd=sd(value), quant005=quantile(value, probs = 0.05),
                           quant95=quantile(value, probs = 0.95)))))
colnames(absCN_df_clades_PDO6) <- paste0(rep(c('clade1.', 'clade2.', 'clade3.'), each=ncol(absCN_df_clades_PDO6)/3), colnames(absCN_df_clades_PDO6))
absCN_df_clades_PDO6$chrom <- clean_chrom(gsub("\\..*", "", (absCN_df_clades_PDO6$clade3.Var2)))

pdf(paste0("../plots/subclonal_clades_cells_summarised_chromosomes_PDO6.pdf"), width = 5, height = 10)
grid.arrange(ggplot(absCN_df_clades_PDO6,
                    aes(x=clade1.mean, y=clade2.mean))+geom_point(aes(col=1))+
               geom_errorbar(aes(xmin=clade1.mean-clade1.sd, xmax=clade1.mean+clade1.sd,
                                 ymin=clade2.mean-clade2.sd, ymax=clade2.mean+clade2.sd))+
               geom_abline(slope = 1, intercept = 0, lty='dashed')+theme_bw()+ggtitle('Clade differences in PDO6')+
               labs(x='CN mean across cells in Clade 1', y='CN mean across cells in Clade 2')+
               facet_wrap(.~factor(as.numeric(chrom)))+guides(col=FALSE),
             ggplot(absCN_df_clades_PDO6,
                    aes(x=clade1.mean, y=clade3.mean))+geom_point(aes(col=1))+
               geom_errorbar(aes(xmin=clade1.mean-clade1.sd, xmax=clade1.mean+clade1.sd,
                                 ymin=clade3.mean-clade3.sd, ymax=clade3.mean+clade3.sd))+
               geom_abline(slope = 1, intercept = 0, lty='dashed')+theme_bw()+ggtitle('Clade differences in PDO6')+
               labs(x='CN mean across cells in Clade 1', y='CN mean across cells in Clade 3')+
               facet_wrap(.~factor(as.numeric(chrom)))+guides(col=FALSE)
)
dev.off()

pvals_clade1_2 <- as.numeric(sapply(1:ncol(absCN_list_PDO6[[1]]), function(idx_col) try(t.test(absCN_list_PDO6[[1]][,idx_col],
                                                                                               absCN_list_PDO6[[2]][,idx_col])$p.value)))
absCN_df_clades_PDO6$clade1.Var2 <- factor(absCN_df_clades_PDO6$clade1.Var2, levels=absCN_df_clades_PDO6$clade1.Var2[order(pvals_clade1_2)])
absCN_df_clades_PDO6$pvalstat1_2 <- p.adjust(pvals_clade1_2) < 0.05

ggplot(absCN_df_clades_PDO6[order(pvals_clade1_2),],
       aes(x=clade1.Var2, y= clade1.mean-clade2.mean, col=factor(ifelse(pvalstat1_2, as.numeric(chrom), NA))))+
  geom_point()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank())+geom_abline(slope = 0, intercept = 0, lty='dashed')+
  theme(legend.position = "bottom")+labs(col='Chromosome')
ggsave(paste0("../plots/subclonal_clades_cells_PDO6_pvals_clade1_2.pdf"), width = 5, height = 5)

ggplot(absCN_df_clades_PDO6[order(pvals_clade1_2),],
       aes(x=clade1.Var2, y= clade1.mean-clade2.mean, col=factor(ifelse(pvalstat1_2, as.numeric(chrom), NA))))+
  geom_point()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank())+geom_abline(slope = 0, intercept = 0, lty='dashed')+
  theme(legend.position = "bottom")+labs(col='Chromosome')+facet_wrap(.~as.numeric(chrom))
ggsave(paste0("../plots/subclonal_clades_cells_PDO6_pvals_clade1_2_facet.pdf"), width = 10, height = 10)

one_to_22_chrom_pdo6 <- absCN_df_clades_PDO6$chrom %in% as.character(1:22)
ggplot(chrom_to_num(absCN_df_clades_PDO6[one_to_22_chrom_pdo6,]),
       aes(x=clade1.Var2, y= clade1.mean-clade2.mean, col=factor(ifelse(pvalstat1_2, as.numeric(chrom), NA))))+
  geom_bar(stat = "identity")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+geom_abline(slope = 0, intercept = 0, lty='dashed')+
  theme(legend.position = "bottom",
        strip.background = element_blank(), strip.text = element_blank())+labs(col='Chromosome')+
  facet_grid (.~ (chrom), scales = "free_x", space = "free_x")+
  guides(col=F)+labs(y='PDO6\nCN (A) - CN (B)')+
  scale_color_manual(values = add_names(unique(ph$gtable$grobs[[4]]$gp$fill[,1])[1:22], as.character(1:22)))
ggsave(paste0("../plots/subclonal_clades_cells_PDO6_pvals_clade1_2_facet_bar.pdf"), width = 8, height = 1.5)


onerow_dummy3 <- absCN_df_clades_PDO6[sapply(unique(absCN_df_clades_PDO6$chrom), function(i) which(absCN_df_clades_PDO6$chrom == i)[1]),]
ggplot(absCN_df_clades_PDO6[absCN_df_clades_PDO6$pvalstat1_2,],
       aes(x=clade1.mean, y= clade2.mean))+
  geom_rect(data=onerow_dummy3, mapping=aes(xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=2), color=NA, fill="blue", alpha=0.3)+
  geom_rect(data=onerow_dummy3, mapping=aes(xmin=-Inf, xmax=Inf, ymin=2, ymax=Inf), color=NA, fill="yellow", alpha=0.3)+
  geom_point()+
  geom_errorbar(aes(xmin=clade1.mean-clade1.sd, xmax=clade1.mean+clade1.sd,
                    ymin=clade2.mean-clade2.sd, ymax=clade2.mean+clade2.sd))+
  geom_abline(slope = 1, intercept = 0, lty='dashed')+theme_bw()+
  geom_abline(slope = 0, intercept = 0, lty='dashed')+
  theme(legend.position = "bottom")+labs(col='Chromosome')+facet_wrap(.~as.numeric(chrom), scales = "free")
ggsave(paste0("../plots/subclonal_clades_cells_summarised_chromosomes_PDO6_DA_regions_1v2.pdf"), width = 10, height = 10)


pvals_clade1_3 <- as.numeric(sapply(1:ncol(absCN_list_PDO6[[1]]), function(idx_col) try(t.test(absCN_list_PDO6[[1]][,idx_col],
                                                                                               absCN_list_PDO6[[3]][,idx_col])$p.value)))
absCN_df_clades_PDO6$clade1.Var2 <- factor(absCN_df_clades_PDO6$clade1.Var2, levels=absCN_df_clades_PDO6$clade1.Var2[order(pvals_clade1_3)])
absCN_df_clades_PDO6$pvalstat1_3 <- p.adjust(pvals_clade1_3) < 0.05

ggplot(absCN_df_clades_PDO6[order(pvals_clade1_3),],
       aes(x=clade1.Var2, y= clade1.mean-clade3.mean, col=factor(ifelse(pvalstat1_3, as.numeric(chrom), NA))))+
  geom_point()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank())+geom_abline(slope = 0, intercept = 0, lty='dashed')+
  theme(legend.position = "bottom")+labs(col='Chromosome')
ggsave(paste0("../plots/subclonal_clades_cells_PDO6_pvals_clade1_3.pdf"), width = 5, height = 5)


ggplot(absCN_df_clades_PDO6[order(pvals_clade1_3),],
       aes(x=clade1.Var2, y= clade1.mean-clade3.mean, col=factor(ifelse(pvalstat1_3, as.numeric(chrom), NA))))+
  geom_point()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank())+geom_abline(slope = 0, intercept = 0, lty='dashed')+
  theme(legend.position = "bottom")+labs(col='Chromosome')+facet_wrap(.~as.numeric(chrom))
ggsave(paste0("../plots/subclonal_clades_cells_PDO6_pvals_clade1_3_facet.pdf"), width = 10, height = 10)

ggplot(chrom_to_num(absCN_df_clades_PDO6[one_to_22_chrom_pdo6,]),
       aes(x=clade1.Var2, y= clade1.mean-clade3.mean, col=factor(ifelse(pvalstat1_3, as.numeric(chrom), NA))))+
  geom_bar(stat = "identity")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+geom_abline(slope = 0, intercept = 0, lty='dashed')+
  theme(legend.position = "bottom",
        strip.background = element_blank(), strip.text = element_blank())+labs(col='Chromosome')+
  facet_grid (.~ (chrom), scales = "free_x", space = "free_x")+
  guides(col=F)+labs(y='PDO6\nCN (A) - CN (C)')+
  scale_color_manual(values = add_names(unique(ph$gtable$grobs[[4]]$gp$fill[,1])[1:22], as.character(1:22)))
ggsave(paste0("../plots/subclonal_clades_cells_PDO6_pvals_clade1_3_facet_bar.pdf"), width = 8, height = 1.5)


ggplot(absCN_df_clades_PDO6[absCN_df_clades_PDO6$pvalstat1_3,],
       aes(x=clade1.mean, y= clade3.mean))+
  geom_rect(data=onerow_dummy3, mapping=aes(xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=2), color=NA, fill="blue", alpha=0.3)+
  geom_rect(data=onerow_dummy3, mapping=aes(xmin=-Inf, xmax=Inf, ymin=2, ymax=Inf), color=NA, fill="yellow", alpha=0.3)+
  geom_point()+
  geom_errorbar(aes(xmin=clade1.mean-clade1.sd, xmax=clade1.mean+clade1.sd,
                    ymin=clade3.mean-clade3.sd, ymax=clade3.mean+clade3.sd))+
  geom_abline(slope = 1, intercept = 0, lty='dashed')+theme_bw()+
  geom_abline(slope = 0, intercept = 0, lty='dashed')+
  theme(legend.position = "bottom")+labs(col='Chromosome')+facet_wrap(.~as.numeric(chrom), scales = "free")
ggsave(paste0("../plots/subclonal_clades_cells_summarised_chromosomes_PDO6_DA_regions_1v3.pdf"), width = 10, height = 10)


remove_na <- function(i){i[!is.na(i)]}

par(mfrow=c(1,2))
image(absCN$PDO6[cutree(hm2$PDO6$tree_row, k = kclades_pdo6) == 1,remove_na(match(as.character(absCN_df_clades_PDO6$clade1.Var2[absCN_df_clades_PDO6$pvalstat1_3]), colnames(absCN$PDO6)))])
image(absCN$PDO6[cutree(hm2$PDO6$tree_row, k = kclades_pdo6) == 3,remove_na(match(as.character(absCN_df_clades_PDO6$clade1.Var2[absCN_df_clades_PDO6$pvalstat1_3]), colnames(absCN$PDO6)))])

#-------------------------------------------------------------------------#
## Subclonality in genes of interest

names_org_num <- c("23868org", "118976org", "119148orgb")
names_org_num2 <- c("23868org", "118976org", "119148org")

## this part of the analysis is done with 20kb bins, instead of the 500kb bins from above
## load these binned CN
absCN_bed = lapply(names_org_num,
                   function(i) read.table(paste0("../UID-", i, ".bed")))
names(absCN_bed) <- names_org_num2
absCN_bed_granges = lapply(names_org_num2, function(pdo){
  as(cbind.data.frame(chrom=absCN_bed[[pdo]]$V1, start=absCN_bed[[pdo]]$V2, end=absCN_bed[[pdo]]$V3,
                                    segVal=absCN_bed[[pdo]]$V5, names=absCN_bed[[pdo]]$V4), 'GRanges')
}) ## Granges must be in env
names(absCN_bed_granges) <- names_org_num2
names(absCN_bed_granges) <- names(absCN_bed) <- renaming$PDO[match(names(absCN_bed_granges), renaming$ID)]

## intercept with genes of interest
subset_genes = c('MYC', 'CCNE1', 'PIK3CA', 'TERT', 'KRAS', 'PTEN', 'RB1', 'AKT1',
                 'AKT2', 'PARP1', 'PARP2', 'ATM', 'ATR', 'WEE1', 'TOP1', 'TUBB1', 'ZWINT', 'ERBB2')

self_select <- function(i) i[i]
table(cutree(hm2[['PDO2']]$tree_row, k = kclades_pdo2))
table(cutree(hm2[['PDO3']]$tree_row, k = kclades_pdo3))
table(cutree(hm2[['PDO6']]$tree_row, k = kclades_pdo6))
length(unique(absCN_bed_granges$PDO2[absCN_bed_granges$PDO2$names %in% names(self_select(cutree(hm2[['PDO2']]$tree_row, k = kclades_pdo2) == 1)),]$names))
length(unique(absCN_bed_granges$PDO2[absCN_bed_granges$PDO2$names %in% names(self_select(cutree(hm2[['PDO2']]$tree_row, k = kclades_pdo2) == 4)),]$names))

saveRDS((absCN_bed_granges$PDO2[absCN_bed_granges$PDO2$names %in% names(self_select(cutree(hm2[['PDO2']]$tree_row, k = kclades_pdo2) == 1)),]),
        "../robjects/table_absCNscDNA_PDO2_cladeA_1.RDS")
saveRDS((absCN_bed_granges$PDO2[absCN_bed_granges$PDO2$names %in% names(self_select(cutree(hm2[['PDO2']]$tree_row, k = kclades_pdo2) == 4)),]),
        "../robjects/table_absCNscDNA_PDO2_cladeB_4.RDS")
saveRDS((absCN_bed_granges$PDO3[absCN_bed_granges$PDO3$names %in% names(self_select(cutree(hm2[['PDO3']]$tree_row, k = kclades_pdo3) == 3)),]),
        "../robjects/table_absCNscDNA_PDO3_cladeA_3.RDS")
saveRDS((absCN_bed_granges$PDO3[absCN_bed_granges$PDO3$names %in% names(self_select(cutree(hm2[['PDO3']]$tree_row, k = kclades_pdo3) == 4)),]),
        "../robjects/table_absCNscDNA_PDO3_cladeB_4.RDS")
saveRDS((absCN_bed_granges$PDO3[absCN_bed_granges$PDO3$names %in% names(self_select(cutree(hm2[['PDO3']]$tree_row, k = kclades_pdo3) == 6)),]),
        "../robjects/table_absCNscDNA_PDO3_cladeC_6.RDS")
saveRDS((absCN_bed_granges$PDO6[absCN_bed_granges$PDO6$names %in% names(self_select(cutree(hm2[['PDO6']]$tree_row, k = kclades_pdo6) == 1)),]),
        "../robjects/table_absCNscDNA_PDO6_cladeA_1.RDS")
saveRDS((absCN_bed_granges$PDO6[absCN_bed_granges$PDO6$names %in% names(self_select(cutree(hm2[['PDO6']]$tree_row, k = kclades_pdo6) == 2)),]),
        "../robjects/table_absCNscDNA_PDO6_cladeB_2.RDS")
saveRDS((absCN_bed_granges$PDO6[absCN_bed_granges$PDO6$names %in% names(self_select(cutree(hm2[['PDO6']]$tree_row, k = kclades_pdo6) == 3)),]),
        "../robjects/table_absCNscDNA_PDO6_cladeC_3.RDS")

# PDO2_cladeA_1granges <- readRDS("../robjects/table_absCNscDNA_PDO2_cladeA_1.RDS")
# PDO2_cladeB_4granges <- readRDS("../robjects/table_absCNscDNA_PDO2_cladeB_4.RDS")


dcast(data.frame(PDO2_cladeA_1granges[seqnames(PDO2_cladeA_1granges) == 1]), .~names, value.var = "segVal")


ag_subsetchrom <- ag_subsetchrom[ag_subsetchrom$gene_name %in% subset_genes,]
subset_genes %in% ag_subsetchrom$gene_name
max(table(ag_subsetchrom$gene_name))
ag_subsetchrom <- as(cbind.data.frame(chrom=ag_subsetchrom$seq_name,
                                      start=ag_subsetchrom$gene_seq_start, end=ag_subsetchrom$gene_seq_end,
                                      names=ag_subsetchrom$gene_name), 'GRanges')

source("../../copy_number_analysis_organoids/helper_functions_granges.R")

if(!cached_ag_subsetchrom){
  chromlens <- load_chromlens()
  saveRDS(chromlens, "../robjects/chromlensGRCh37.87.RDS")
}else{
  chromlens <- readRDS("../robjects/chromlensGRCh37.87.RDS")
}

cached_average_CN_per_gene_selected_genes <- T
if(cached_average_CN_per_gene_selected_genes){
  average_CN_per_gene <- readRDS("../robjects/average_CN_per_gene_selected_genes_all_cells.RDS")
}else{
  average_CN_per_gene <- lapply(names(absCN_bed_granges), function(sample_idx){ lapply(sort(unique(absCN_bed_granges[[sample_idx]]$names)), function(cell_idx){
    give_CN_per_gene_v2(segment_arg = absCN_bed_granges[[sample_idx]][absCN_bed_granges[[sample_idx]]$names == cell_idx,],
                        gr_genes = as(ag_subsetchrom, 'GRanges'))
  })
  })
  names(average_CN_per_gene) <- names(absCN_bed_granges)
  saveRDS(average_CN_per_gene, "../robjects/average_CN_per_gene_selected_genes_all_cells.RDS")
  
}
average_CN_per_gene_melt <- lapply(average_CN_per_gene, melt)

min(absCN_bed$PDO2$V4)
max(absCN_bed$PDO2$V4)
min(absCN_bed_granges$PDO2$names)
min(absCN_bed$PDO3$V4)
max(absCN_bed$PDO3$V4)
min(absCN_bed$PDO6$V4)
max(absCN_bed$PDO6$V4)

(read_UID_orgs_raw(renaming$ID[match(org_it, renaming$PDO)]))
min((read_UID_orgs_raw(renaming$ID[match(organoid_list[1], renaming$PDO)]))$node_id)
max((read_UID_orgs_raw(renaming$ID[match(organoid_list[1], renaming$PDO)]))$node_id)

average_CN_per_gene_melt$PDO2$L1 <- sort(unique(absCN_bed_granges$PDO2$names))[average_CN_per_gene_melt$PDO2$L1]
average_CN_per_gene_melt$PDO3$L1 <- sort(unique(absCN_bed_granges$PDO3$names))[average_CN_per_gene_melt$PDO3$L1]
average_CN_per_gene_melt$PDO6$L1 <- sort(unique(absCN_bed_granges$PDO6$names))[average_CN_per_gene_melt$PDO6$L1]

## select only cells which are included in the clonal analysis (filtered)
for(org_it in names(average_CN_per_gene_melt)){
  average_CN_per_gene_melt[[org_it]] <- average_CN_per_gene_melt[[org_it]][average_CN_per_gene_melt[[org_it]]$L1 %in% rownames(absCN[[org_it]]),]
}

## find which cell belongs to which clone

table(cutree(hm2[['PDO2']]$tree_row, k = kclades_pdo2))
average_CN_per_gene_melt_PDO2 <- lapply(c(1, 4), function(idx_clade) average_CN_per_gene_melt$PDO2[average_CN_per_gene_melt$PDO2$L1 %in% rownames(absCN$PDO2[(cutree(hm2[['PDO2']]$tree_row, k = kclades_pdo2)) == idx_clade,]),])
ggplot(melt(average_CN_per_gene_melt_PDO2), aes(x=gene, y=value, group=interaction(gene,L1), col=factor(L1)))+
  geom_boxplot()+scale_y_continuous(trans = "log2")+theme_bw()+labs(y='CN (log2)')
ggsave(paste0("../plots/subclonal_clades_cells_PDO2_genes_of_interest.pdf"), width = 10, height = 2.5)

library(ComplexHeatmap)
ComplexHeatmap::Heatmap(apply(cbind(dcast(average_CN_per_gene_melt$PDO2[average_CN_per_gene_melt$PDO2$L1 %in% rownames(absCN$PDO2[(cutree(hm2[['PDO2']]$tree_row, k = kclades_pdo2)) == 1,]),], gene ~ L1, value.var = "value")[,-1],
      dcast(average_CN_per_gene_melt$PDO2[average_CN_per_gene_melt$PDO2$L1 %in% rownames(absCN$PDO2[(cutree(hm2[['PDO2']]$tree_row, k = kclades_pdo2)) == 4,]),], gene ~ L1, value.var = "value")[,-1]), 2, as.numeric),
      top_annotation = HeatmapAnnotation(df=data.frame(rep(c('cladeA', 'cladeB'), table(cutree(hm2[['PDO2']]$tree_row, k = kclades_pdo2))[c(1,4)]))))

table(cutree(hm2[['PDO3']]$tree_row, k = kclades_pdo3))
average_CN_per_gene_melt_PDO3 <- lapply(c(3, 4, 6), function(idx_clade) average_CN_per_gene_melt$PDO3[average_CN_per_gene_melt$PDO3$L1 %in% rownames(absCN$PDO3[(cutree(hm2[['PDO3']]$tree_row, k = kclades_pdo3)) == idx_clade,]),])
ggplot(melt(average_CN_per_gene_melt_PDO3), aes(x=gene, y=value, group=interaction(gene,L1), col=factor(L1)))+
  geom_boxplot()+scale_y_continuous(trans = "log2")+theme_bw()+labs(y='CN (log2)')
ggsave(paste0("../plots/subclonal_clades_cells_PDO3_genes_of_interest.pdf"), width = 10, height = 2.5)

table(cutree(hm2[['PDO6']]$tree_row, k = kclades_pdo6))
average_CN_per_gene_melt_PDO6 <- lapply(c(1, 2, 3), function(idx_clade) average_CN_per_gene_melt$PDO6[average_CN_per_gene_melt$PDO6$L1 %in% rownames(absCN$PDO6[(cutree(hm2[['PDO6']]$tree_row, k = kclades_pdo6)) == idx_clade,]),])
ggplot(melt(average_CN_per_gene_melt_PDO6), aes(x=gene, y=value, group=interaction(gene,L1), col=factor(L1)))+
  geom_boxplot()+scale_y_continuous(trans = "log2")+theme_bw()+labs(y='CN (log2)')
ggsave(paste0("../plots/subclonal_clades_cells_PDO6_genes_of_interest.pdf"), width = 10, height = 2.5)


length(table(average_CN_per_gene_melt$PDO2$L1))
for(org_it in names(average_CN_per_gene_melt)){
  ggplot(average_CN_per_gene_melt[[org_it]], aes(x=factor(gene, levels=unlist(order_first_by_second_column(average_CN_per_gene_melt[[org_it]] %>% group_by(gene) %>% summarise(median(value))))),
                                            y=value))+geom_boxplot()+
    theme_bw()+geom_abline(slope = 0, intercept = 2, lty='dashed')+
    labs(x='Gene', y='CN')
  ggsave(paste0("../plots/CN_", org_it, "_genes_of_interest.pdf"), width = 10, height = 2.5)
}
for(org_it in names(average_CN_per_gene_melt)){
  ggplot(average_CN_per_gene_melt[[org_it]], aes(x=factor(gene, levels=unlist(order_first_by_second_column(average_CN_per_gene_melt[[org_it]] %>% group_by(gene) %>% summarise(median(value))))),
                                                 y=value))+geom_boxplot()+
    theme_bw()+geom_abline(slope = 0, intercept = log2(2), lty='dashed')+
    labs(x='Gene', y='CN')+scale_y_continuous(trans = "log2")
  ggsave(paste0("../plots/CN_", org_it, "_genes_of_interest_log2.pdf"), width = 10, height = 2.5)
}

replace_factors <- function(i){
  sapply(i, function(j){
    if(j == 'FALSEFALSE'){
      'Diploid'
    }else if(j == 'FALSETRUE'){
      'Amplification'
    }else if(j == 'TRUEFALSE'){
      'Loss'
    }
  })
}
for(org_it in names(average_CN_per_gene_melt)){
  ggplot(average_CN_per_gene_melt[[org_it]], aes(x=factor(gene, levels=unlist(order_first_by_second_column(average_CN_per_gene_melt[[org_it]] %>% group_by(gene) %>% summarise(median(value))))),
                                               fill=replace_factors(paste0((value < 2), (value > 2))) ))+geom_boxplot()+
  geom_bar()+
  # facet_wrap(.~factor(gene, levels=unlist(order_first_by_second_column(average_CN_per_gene_melt[[org_it]] %>% group_by(gene) %>% summarise(median(value))))), nrow=1)+
  theme_bw()+
  labs(x='Gene', y='CN', fill='Status')+#+scale_y_continuous(trans = "log2")
    ggtitle(org_it)+theme(legend.position = "bottom")+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave(paste0("../plots/CN_", org_it, "_genes_of_interest_barplot.pdf"), width =6, height = 2.5)
}

order_first_by_second_column <- function(i) i[order(i[,2]),1]


# sort(unique(absCN_bed_granges$PDO2$names)) ## 0-indexed
# sort(unique(absCN_bed_granges$PDO3$names)) ## 0-indexed
# sort(unique(absCN_bed_granges$PDO6$names)) ## 0-indexed
# 
# 
# dim(absCN$PDO3)
# 
# 
# stopifnot(renaming$ID[match(names(absCN_bed_granges), renaming$PDO)] == names_org_num2)
# 
# keep_cells <- sapply(names_org_num, function(i) as.numeric(rownames(read_UID_orgs(i)))) ## 1-indexed
# sapply(keep_cells, length)
# names(keep_cells) <- names(absCN_bed_granges)
# keep_cells$PDO2 <- keep_cells$PDO2[-outliers$PDO2]
# keep_cells$PDO3 <- keep_cells$PDO3[-outliers$PDO3]
# keep_cells$PDO6 <- keep_cells$PDO6[-outliers$PDO6]
# sapply(keep_cells, length) ## removal of outliers
# keep_cells_0index <- sapply(keep_cells, function(i) i-1)
# 
# sapply(keep_cells, length)
# c(nrow(absCN$PDO2), nrow(absCN$PDO3), nrow(absCN$PDO6))
# c(length(unique(absCN_bed_granges$PDO2$names)),
# length(unique(absCN_bed_granges$PDO3$names)),
# length(unique(absCN_bed_granges$PDO6$names)))
# 
# ## zero-indexed
# # min(average_CN_per_gene_melt$PDO2$L1)
# # min(average_CN_per_gene_melt$PDO3$L1)
# # min(average_CN_per_gene_melt$PDO6$L1)
# average_CN_per_gene_melt$PDO2 <- average_CN_per_gene_melt$PDO2[average_CN_per_gene_melt$PDO2$L1 %in% keep_cells_0index$PDO2,]
# average_CN_per_gene_melt$PDO3 <- average_CN_per_gene_melt$PDO3[average_CN_per_gene_melt$PDO3$L1 %in% keep_cells_0index$PDO3,]
# average_CN_per_gene_melt$PDO6 <- average_CN_per_gene_melt$PDO6[average_CN_per_gene_melt$PDO6$L1 %in% keep_cells_0index$PDO6,]
# 
# c(length(unique(average_CN_per_gene_melt$PDO2$L1)),
# length(unique(average_CN_per_gene_melt$PDO3$L1)),
# length(unique(average_CN_per_gene_melt$PDO6$L1)))
# sapply(keep_cells, length)
# c(nrow(absCN$PDO2), nrow(absCN$PDO3), nrow(absCN$PDO6))

### HERE

#-------------------------------------------------------------------------#


ggplot(melt(list(apply(data_orgs[[1]], 1, var),
     apply(data_orgs[[2]], 1, var),
     apply(data_orgs[[3]], 1, var))), aes(x=value, col=L1, group=L1))+geom_density()


ggplot(melt(data_orgs[[1]][,1:100]), aes(x=Var2, y=value))+geom_violin()+geom_jitter(size=0.2)

data_orgs_centered <- lapply(data_orgs, function(j) sweep(j, 2, apply(j, 2, mean), '-'))
names(data_orgs_centered) <- organoid_list

ggplot(melt(data_orgs_centered[[1]][,1:100]), aes(x=Var2, y=value))+geom_violin()+geom_jitter(size=0.2)

# apply(data_orgs_centered[[1]]

ggplot(melt(data_orgs_centered[[1]]), aes(x=Var2, y=value))+geom_violin()

data_orgs_centered_sd <- lapply(data_orgs_centered, function(j){
  cbind.data.frame(pos=colnames(j),
                 sd_bin=apply(j, 2, sd))})

ggplot(melt(data_orgs_centered_sd),
       aes(x=pos, y=value, group=L1, col=L1))+geom_line()+
  labs(y='Standard deviation')

## compute confidence intervals
data_orgs_centered_CI <- data.frame(do.call('rbind', lapply(data_orgs_centered, function(j) (t(apply(j, 2, quantile, c(0.05, 0.95)))))))
data_orgs_centered_CI$L1 = rep(organoid_list, sapply(data_orgs_centered, ncol))
data_orgs_centered_CI$pos=rownames(data_orgs_centered_CI)

put_limits <- function(x, upperlim=1, lowerlim=-1){
  a <- x
  if(x > upperlim) a <- upperlim
  if(x < lowerlim) a <- lowerlim
  return(a)
}

data_orgs_centered_CI$censoredX5 = sapply(data_orgs_centered_CI$X5., put_limits)
data_orgs_centered_CI$censoredX95 = sapply(data_orgs_centered_CI$X95., put_limits)
data_orgs_centered_CI$chrom = clean_chrom(gsub("\\..*", "", rownames(data_orgs_centered_CI)))
data_orgs_centered_CI$chrom = factor(data_orgs_centered_CI$chrom, levels=c(as.character(1:22)))
data_orgs_centered_CI$pos_num = as.numeric(sapply(rownames(data_orgs_centered_CI), function(i) strsplit(i, '[.]')[[1]][2]))
ggplot(data_orgs_centered_CI)+
  geom_ribbon(aes(x=pos_num, ymin=censoredX5, ymax=censoredX95, group=L1))+
  geom_point(aes(x=pos_num, y=0, col=( (censoredX95-censoredX5)< 0.1 )), size=0.1)+
  facet_nested(L1~chrom, scales = "free", space="free_x")+theme_bw()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  lims(y=c(-1, 1))+theme(legend.position = "bottom")+
  labs(x='90% confidence intervals of the scaled CN values of bins, across organoids',
       y='Value of scaled CN of bin')+labs(col='Clonal bin')
ggsave(paste0("../plots/subclonal_hclust_CI_across_genome_2.pdf"), width = 14)

x = data_orgs_centered_CI %>% filter(L1=='PDO3', chrom==7)
plot(x$censoredX5, type='l')

pnorm(q = 0.95)
cbind.data.frame(org=organoid_list, percentage_above_083sd=sapply(data_orgs_centered_sd, function(i) mean((i$sd_bin > 0.83))))

## plot the absolute CN and the variace, because there might be heteroscedasticity
data_orgs_mean_CN <- lapply(data_orgs, function(j){
  cbind.data.frame(pos=colnames(j),
                   mean=colMeans(j))})

plot((data_orgs_mean_CN[[1]]$mean),
     (data_orgs_centered_sd[[1]]$sd_bin)
     )

df_heteroscedasticity <- data.frame(do.call('rbind', lapply(1:3, function(i) cbind(data_orgs_mean_CN[[i]]$mean,
                              data_orgs_centered_sd[[i]]$sd_bin))))
colnames(df_heteroscedasticity) <- c('Mean', 'Sd')
df_heteroscedasticity$PDO <- rep(organoid_list, sapply(1:3, function(i) length(data_orgs_centered_sd[[i]]$sd_bin)))
head(df_heteroscedasticity)

ggplot(df_heteroscedasticity, aes(x=Mean, y=Sd))+geom_point()+
  geom_smooth(se = FALSE, method = lm)+
  facet_wrap(.~PDO)+theme_bw()
ggsave(paste0("../plots/subclonal_hclust_CI_across_genome_heteroscedasticity.pdf"), width = 5, height = 2.5)

mean_sd_lm <- lapply(organoid_list, function(org) lm(Sd ~ Mean,
                         data = df_heteroscedasticity[df_heteroscedasticity$PDO == org,]))
names(mean_sd_lm) <- organoid_list
sapply(mean_sd_lm, coef)

## for each bin, in each organoid, compute the observed standard deviation, and se
## how it compares to the expected standard deviation, with a Chi squared test
pvals_chisq <- lapply(organoid_list, function(org){
  sapply(1:sum(df_heteroscedasticity$PDO == org), function(pos_it){
    .current <- df_heteroscedasticity[df_heteroscedasticity$PDO == org,][1,]
    .expected_sd <- coef(mean_sd_lm[[org]])[1] + .current$Mean*coef(mean_sd_lm[[org]])[2]
    .expected_var <- .expected_sd**2
    .res_chi <- EnvStats::varTest(data_orgs[[org]][,pos_it], sigma.squared = .expected_var, alternative="greater")
    .res_chi$p.value
  })
})
names(pvals_chisq) <- organoid_list
xtable::xtable(data.frame(fraction_heterogeneous=sapply(pvals_chisq, function(j) mean(j < 0.05))))

data_orgs_centered_CI$pvalstat = NA ## init
for(org_it in organoid_list){
  data_orgs_centered_CI[data_orgs_centered_CI$L1 == org_it,'pvalstat'] = pvals_chisq[[org_it]] < 0.05
}
sum(is.na(data_orgs_centered_CI$pvalstat))
table(data_orgs_centered_CI$pvalstat)

ggplot(data_orgs_centered_CI)+
  geom_ribbon(aes(x=pos_num, ymin=censoredX5, ymax=censoredX95, group=L1))+
  geom_point(aes(x=pos_num, y=0, col=( pvalstat )), size=0.1)+
  facet_nested(L1~chrom, scales = "free", space="free_x")+theme_bw()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  lims(y=c(-1, 1))+theme(legend.position = "bottom")+
  labs(x='90% confidence intervals of the scaled CN values of bins, across organoids',
       y='Value of scaled CN of bin')+labs(col='Subclonal bin')
ggsave(paste0("../plots/subclonal_hclust_CI_across_genome_3_test.pdf"), width = 14)


tiled_ecDNA_table_absCNscDNA_PDO2_cladeA_1 <- get_tiled_ecDNA(table_absCNscDNA_PDO2_cladeA_1_split)
saveRDS(tiled_ecDNA_table_absCNscDNA_PDO2_cladeA_1, "../robjects/tiled_ecDNA_table_absCNscDNA_PDO2_cladeA_1.RDS")

tiled_ecDNA_table_absCNscDNA_PDO2_cladeB_4 <- get_tiled_ecDNA(table_absCNscDNA_PDO2_cladeB_4_split)
saveRDS(tiled_ecDNA_table_absCNscDNA_PDO2_cladeB_4, "../robjects/tiled_ecDNA_table_absCNscDNA_PDO2_cladeB_4.RDS")

tiled_ecDNA_table_absCNscDNA_PDO3_cladeA_3 <- get_tiled_ecDNA(table_absCNscDNA_PDO3_cladeA_3_split)
saveRDS(tiled_ecDNA_table_absCNscDNA_PDO3_cladeA_3, "../robjects/tiled_ecDNA_table_absCNscDNA_PDO3_cladeA_3.RDS")

tiled_ecDNA_table_absCNscDNA_PDO3_cladeB_4 <- get_tiled_ecDNA(table_absCNscDNA_PDO3_cladeB_4_split)
saveRDS(tiled_ecDNA_table_absCNscDNA_PDO3_cladeB_4, "../robjects/tiled_ecDNA_table_absCNscDNA_PDO3_cladeB_4.RDS")

tiled_ecDNA_table_absCNscDNA_PDO3_cladeC_6 <- get_tiled_ecDNA(table_absCNscDNA_PDO3_cladeC_6_split)
saveRDS(tiled_ecDNA_table_absCNscDNA_PDO3_cladeC_6, "../robjects/tiled_ecDNA_table_absCNscDNA_PDO3_cladeC_6.RDS")

tiled_ecDNA_table_absCNscDNA_PDO6_cladeA_1 <- get_tiled_ecDNA(table_absCNscDNA_PDO6_cladeA_1_split)
saveRDS(tiled_ecDNA_table_absCNscDNA_PDO6_cladeA_1, "../robjects/tiled_ecDNA_table_absCNscDNA_PDO6_cladeA_1.RDS")

tiled_ecDNA_table_absCNscDNA_PDO6_cladeB_2 <- get_tiled_ecDNA(table_absCNscDNA_PDO6_cladeB_2_split)
saveRDS(tiled_ecDNA_table_absCNscDNA_PDO6_cladeB_2, "../robjects/tiled_ecDNA_table_absCNscDNA_PDO6_cladeB_2.RDS")

tiled_ecDNA_table_absCNscDNA_PDO6_cladeC_3 <- get_tiled_ecDNA(table_absCNscDNA_PDO6_cladeC_3_split)
saveRDS(tiled_ecDNA_table_absCNscDNA_PDO6_cladeC_3, "../robjects/tiled_ecDNA_table_absCNscDNA_PDO6_cladeC_3.RDS")


##----------------------------------------------------------------------------
##----------------------------------------------------------------------------

##----------------------------------------------------------------------------
##----------------------------------------------------------------------------

### ecDNA
### finding areas of extreme CN values

library(Rtsne)
library(ComplexHeatmap)
library(ggrepel)

##----------------

tiled_ecDNA_table_absCNscDNA_PDO2_cladeA_1 <- readRDS("../robjects/tiled_ecDNA_table_absCNscDNA_PDO2_cladeA_1.RDS")
tiled_ecDNA_table_absCNscDNA_PDO2_cladeB_4 <- readRDS("../robjects/tiled_ecDNA_table_absCNscDNA_PDO2_cladeB_4.RDS")
tiled_ecDNA_table_absCNscDNA_PDO3_cladeA_3 <- readRDS("../robjects/tiled_ecDNA_table_absCNscDNA_PDO3_cladeA_3.RDS")
tiled_ecDNA_table_absCNscDNA_PDO3_cladeB_4 <- readRDS("../robjects/tiled_ecDNA_table_absCNscDNA_PDO3_cladeB_4.RDS")
tiled_ecDNA_table_absCNscDNA_PDO3_cladeC_6 <- readRDS("../robjects/tiled_ecDNA_table_absCNscDNA_PDO3_cladeC_6.RDS")
tiled_ecDNA_table_absCNscDNA_PDO6_cladeA_1 <- readRDS("../robjects/tiled_ecDNA_table_absCNscDNA_PDO6_cladeA_1.RDS")
tiled_ecDNA_table_absCNscDNA_PDO6_cladeB_2 <- readRDS("../robjects/tiled_ecDNA_table_absCNscDNA_PDO6_cladeB_2.RDS")
tiled_ecDNA_table_absCNscDNA_PDO6_cladeC_3 <- readRDS("../robjects/tiled_ecDNA_table_absCNscDNA_PDO6_cladeC_3.RDS")

table_absCNscDNA_PDO2_cladeA_1 <- readRDS("../robjects/table_absCNscDNA_PDO2_cladeA_1.RDS")
table_absCNscDNA_PDO2_cladeB_4 <- readRDS("../robjects/table_absCNscDNA_PDO2_cladeB_4.RDS")
table_absCNscDNA_PDO3_cladeA_3 <- readRDS("../robjects/table_absCNscDNA_PDO3_cladeA_3.RDS")
table_absCNscDNA_PDO3_cladeB_4 <- readRDS("../robjects/table_absCNscDNA_PDO3_cladeB_4.RDS")
table_absCNscDNA_PDO3_cladeC_6 <- readRDS("../robjects/table_absCNscDNA_PDO3_cladeC_6.RDS")
table_absCNscDNA_PDO6_cladeA_1 <- readRDS("../robjects/table_absCNscDNA_PDO6_cladeA_1.RDS")
table_absCNscDNA_PDO6_cladeB_2 <- readRDS("../robjects/table_absCNscDNA_PDO6_cladeB_2.RDS")
table_absCNscDNA_PDO6_cladeC_3 <- readRDS("../robjects/table_absCNscDNA_PDO6_cladeC_3.RDS")

##----------------

split_by_name <- function(i){
  nmes <- unique(i$names)
  .x <- lapply(nmes, function(name) i[i$names == name])
  names(.x) <- nmes
  .x
}

table_absCNscDNA_PDO2_cladeA_1_split <- split_by_name(table_absCNscDNA_PDO2_cladeA_1[table_absCNscDNA_PDO2_cladeA_1$segVal > 8,])
table_absCNscDNA_PDO2_cladeB_4_split <- split_by_name(table_absCNscDNA_PDO2_cladeB_4[table_absCNscDNA_PDO2_cladeB_4$segVal > 8,])

table_absCNscDNA_PDO3_cladeA_3_split <- split_by_name(table_absCNscDNA_PDO3_cladeA_3[table_absCNscDNA_PDO3_cladeA_3$segVal > 8,])
table_absCNscDNA_PDO3_cladeB_4_split <- split_by_name(table_absCNscDNA_PDO3_cladeB_4[table_absCNscDNA_PDO3_cladeB_4$segVal > 8,])
table_absCNscDNA_PDO3_cladeC_6_split <- split_by_name(table_absCNscDNA_PDO3_cladeC_6[table_absCNscDNA_PDO3_cladeC_6$segVal > 8,])

table_absCNscDNA_PDO6_cladeA_1_split <- split_by_name(table_absCNscDNA_PDO6_cladeA_1[table_absCNscDNA_PDO6_cladeA_1$segVal > 8,])
table_absCNscDNA_PDO6_cladeB_2_split <- split_by_name(table_absCNscDNA_PDO6_cladeB_2[table_absCNscDNA_PDO6_cladeB_2$segVal > 8,])
table_absCNscDNA_PDO6_cladeC_3_split <- split_by_name(table_absCNscDNA_PDO6_cladeC_3[table_absCNscDNA_PDO6_cladeC_3$segVal > 8,])

tiled_genome <- GenomicRanges::tileGenome(seqlengths(BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19), tilewidth = 20000)
tiled_genome <- unlist(tiled_genome)
tiled_genome <- tiled_genome[!grepl("_", seqnames(tiled_genome)),]
tiled_genome <- tiled_genome[!(seqnames(tiled_genome) %in% c('chrM', 'chrMT')),]
# seqnames(tiled_genome) <- seqnames(tiled_genome)[names(seqlengths(tiled_genome)) %in% unique(unique(as.character(seqnames(tiled_genome))))]
# tiled_genome <- trim(tiled_genome)
# levels(seqnames(tiled_genome)) <- as.character(unique(seqnames(tiled_genome)))
# levels(seqnames(tiled_genome)) <- NULL
seqlevels(tiled_genome) <- unique( as.character(seqnames(tiled_genome)))
# seqnames(tiled_genome) <- Rle(values = unique(as.character(seqnames(tiled_genome))), lengths = seqlengths(tiled_genome)[names(seqlengths(tiled_genome)) %in% as.character(seqnames(tiled_genome))])

chromlens=data.frame(Chrom=names(seqlengths(tiled_genome))[names(seqlengths(tiled_genome)) %in% unique(unique(as.character(seqnames(tiled_genome))))],
                     Length=seqlengths(tiled_genome)[names(seqlengths(tiled_genome)) %in% unique(unique(as.character(seqnames(tiled_genome))))])

remove_names <- function(i){
  i$names <- NULL
  i
}
# names(tiled_genome) <- 1:length(tiled_genome)

tiled_genome$names <-  1:length(tiled_genome)
table_absCNscDNA_PDO2_cladeA_1_split
get_tiled_ecDNA <- function(tiled_genome_arg){
  tiled_ecDNA_table_absCNscDNA_PDO2_cladeA_1 <- lapply(tiled_genome_arg, function(tab){
    cat('New cell\n')
    give_CN_per_gene_v2(gr_genes = tiled_genome,
                   segment_arg = remove_names(tab),
                   transform_segment_dataframe = F, only_compute_overlaps=T)
  })
  names(tiled_ecDNA_table_absCNscDNA_PDO2_cladeA_1) <- names(tiled_genome_arg)
  return(tiled_ecDNA_table_absCNscDNA_PDO2_cladeA_1)
}

##----------------


##----------------
list_tiles <- list(tiled_ecDNA_table_absCNscDNA_PDO2_cladeA_1,
     tiled_ecDNA_table_absCNscDNA_PDO2_cladeB_4,
     tiled_ecDNA_table_absCNscDNA_PDO3_cladeA_3,
     tiled_ecDNA_table_absCNscDNA_PDO3_cladeB_4,
     tiled_ecDNA_table_absCNscDNA_PDO3_cladeC_6,
     tiled_ecDNA_table_absCNscDNA_PDO6_cladeA_1,
     tiled_ecDNA_table_absCNscDNA_PDO6_cladeB_2,
     tiled_ecDNA_table_absCNscDNA_PDO6_cladeC_3)

annotation_cells=data.frame(clade=rep(c('PDO2 A', 'PDO2 B', 'PDO3 A', 'PDO3 B', 'PDO3 C', 'PDO6 A', 'PDO6 B', 'PDO6 C'), sapply(list_tiles, length)))

tiled_clades_matrix <- lapply(list_tiles, function(j){
  .x <- sapply(j, function(i) i[,2])
  .x[is.na(.x)] <- 0
  .x
})
tiled_clades_matrix_sumNA <- sapply(tiled_clades_matrix, function(i) rowSums(i))
tiled_clades_matrix_sumNA_allzero <- rowSums(tiled_clades_matrix_sumNA) == 0
tiled_clades_matrix_pruned <- lapply(1:length(tiled_clades_matrix), function(i){
  tiled_clades_matrix[[i]][!tiled_clades_matrix_sumNA_allzero,]
})

##----------------

rm(tiled_ecDNA_table_absCNscDNA_PDO2_cladeA_1)
rm(tiled_ecDNA_table_absCNscDNA_PDO2_cladeB_4)
rm(tiled_ecDNA_table_absCNscDNA_PDO3_cladeA_3)
rm(tiled_ecDNA_table_absCNscDNA_PDO3_cladeB_4)
rm(tiled_ecDNA_table_absCNscDNA_PDO3_cladeC_6)
rm(tiled_ecDNA_table_absCNscDNA_PDO6_cladeA_1)
rm(tiled_ecDNA_table_absCNscDNA_PDO6_cladeB_2)
rm(list_tiles)
rm(tiled_clades_matrix)
# rm(tiled_clades_matrix_pruned)
rm(absCN_list_PDO2)
rm(absCN_list_PDO3)
rm(absCN_list_PDO6)
rm(absCN)
rm(absCN0)
rm(data_orgs)
rm(tiled_clades_matrix_sumNA)
rm(tiled_ecDNA_table_absCNscDNA_PDO6_cladeC_3)

##----------------

sapply(tiled_clades_matrix_pruned, dim)
tiled_clades_matrix_pruned_cbind <- do.call('cbind', tiled_clades_matrix_pruned)
# image(tiled_clades_matrix_pruned_cbind) ## TAKES TOO LONG - only regions in which at least in one cell there is a putative ecDNA region

umap_tiled_clades_matrix_pruned <- umap::umap(t(tiled_clades_matrix_pruned_cbind))
ggplot(cbind.data.frame(umap_tiled_clades_matrix_pruned$layout, col=rep(c('PDO2 A', 'PDO2 B', 'PDO3 A', 'PDO3 B', 'PDO3 C', 'PDO6 A', 'PDO6 B', 'PDO6 C'), sapply(tiled_clades_matrix_pruned, ncol))),
       aes(x=`1`, y=`2`, col=col))+
  geom_point()+theme_bw()

tsne_tiled_clades_matrix_pruned <- Rtsne::Rtsne(t(tiled_clades_matrix_pruned_cbind), labels=rep(c('PDO2 A', 'PDO2 B', 'PDO3 A', 'PDO3 B', 'PDO3 C', 'PDO6 A', 'PDO6 B', 'PDO6 C'), sapply(tiled_clades_matrix_pruned, ncol)))
ggplot(cbind.data.frame(tsne_tiled_clades_matrix_pruned$Y, col=rep(c('PDO2 A', 'PDO2 B', 'PDO3 A', 'PDO3 B', 'PDO3 C', 'PDO6 A', 'PDO6 B', 'PDO6 C'), sapply(tiled_clades_matrix_pruned, ncol))),
       aes(x=`1`, y=`2`, col=col))+
  geom_point()+theme_bw()
ggsave("../plots/tsne_ecDNA_scDNA.pdf")

## annotate the bins
length(tiled_clades_matrix_sumNA_allzero) == length(tiled_genome)
names_tiles <- apply(data.frame(tiled_genome)[,1:3], 1, paste0, collapse='_')
rownames(tiled_clades_matrix_pruned_cbind) <- names_tiles[!tiled_clades_matrix_sumNA_allzero]

tiled_clades_matrix_pruned_cbind_topvar <- apply(tiled_clades_matrix_pruned_cbind, 1, var)
tiled_clades_matrix_pruned_cbind_topvar <- tiled_clades_matrix_pruned_cbind[order(tiled_clades_matrix_pruned_cbind_topvar, decreasing = T)[1:50],]
colnames(tiled_clades_matrix_pruned_cbind_topvar) = paste0(annotation_cells$clade, ' - ', colnames(tiled_clades_matrix_pruned_cbind_topvar))
rownames(annotation_cells) = colnames(tiled_clades_matrix_pruned_cbind_topvar)
  
annotation_cells$clade <- gsub(" ", "_", annotation_cells$clade)
col_fun = circlize::colorRamp2(c(0, 5, 10), c("blue", "white", "red"))

annotation_cells_annotation <- ComplexHeatmap::HeatmapAnnotation(df=annotation_cells,
                  col=list(clade=c(PDO2_A='green', PDO2_B='#00ff7f', PDO3_A='#e4717a', PDO3_B='#f06949',
                           PDO3_C='#ff0038', PDO6_A='#4b60c6', PDO6_B='#00008b', PDO6_C='#003cd5')))

pdf("../plots/ecDNA_heatmap_topvar.pdf", width = 9)
ComplexHeatmap::Heatmap(tiled_clades_matrix_pruned_cbind_topvar,
                        top_annotation = annotation_cells_annotation,
                        show_column_names = F)

dev.off()

## random selection of positions (nice, sometimes!)
pdf("../plots/ecDNA_heatmap_randomselection.pdf", width = 9)
ComplexHeatmap::Heatmap(tiled_clades_matrix_pruned_cbind[sample(1:nrow(tiled_clades_matrix_pruned_cbind), 100),],
                        top_annotation = annotation_cells_annotation,
                        show_column_names = F)
dev.off()


# plot(hclust(dist(tiled_clades_matrix_pruned_cbind[order(rowSums(tiled_clades_matrix_pruned_cbind),decreasing = T)[1:100])))

pdf("../plots/ecDNA_heatmap_topsum.pdf", width = 9)
ComplexHeatmap::Heatmap(tiled_clades_matrix_pruned_cbind[order(rowSums(tiled_clades_matrix_pruned_cbind),decreasing = T)[1:200],],
                        top_annotation = annotation_cells_annotation,
                        show_column_names = F)
dev.off()

give_annotation_chroms <- function(regions_vec){
  unique_chroms <- gsub("_.*", "", regions_vec)
  ComplexHeatmap::rowAnnotation(df=unique_chroms,
                                col=list(chrom=give_names(viridisLite::magma(nrow(chromlens)), chromlens$Chrom)[unique(unique_chroms)]))
}
give_names <- function(i,j){
  names(i) <- j
  i
}

## for each organoid, get the areas of highest CN
tiled_clades_matrix_pruned_cbind_orgtop <- sapply(c('PDO2', 'PDO3', 'PDO6'),
      function(org){
        .xx <- tiled_clades_matrix_pruned_cbind[,gsub("_.*", "", annotation_cells$clade) == org]
        order(rowSums(.xx),decreasing = T)[1:200]
      })
tiled_clades_matrix_pruned_cbind_orgtop <- tiled_clades_matrix_pruned_cbind[unique(as.vector(tiled_clades_matrix_pruned_cbind_orgtop)),]
tiled_clades_matrix_pruned_cbind_orgtop_anno <- give_annotation_chroms(rownames(tiled_clades_matrix_pruned_cbind_orgtop))

pdf("../plots/ecDNA_heatmap_toporg.pdf", width = 9)
ComplexHeatmap::Heatmap(tiled_clades_matrix_pruned_cbind_orgtop,
                        top_annotation = annotation_cells_annotation,
                        right_annotation = tiled_clades_matrix_pruned_cbind_orgtop_anno,
                        show_column_names = F, show_row_names = F)
dev.off()

pdf("../plots/ecDNA_heatmap_chrom8.pdf", width = 9)
ComplexHeatmap::Heatmap(tiled_clades_matrix_pruned_cbind[(grepl('chr8_', rownames(tiled_clades_matrix_pruned_cbind))),],
                        top_annotation = annotation_cells_annotation,
                        show_column_names = F, show_row_names = F, cluster_rows = F)
dev.off()

#------------------------------
## annotation of genes in these regions
add_colnames <- function(i){
  colnames(i) <- c('seqnames', 'start', 'end')
  i
}
modify_ag_colnames <- function(i){
  colnames(i) <- gsub("gene_", "", colnames(i))
  colnames(i)[colnames(i) == 'seq_name'] = 'seqname'
  i <- i[,c('seqname', 'seq_start', 'seq_end', 'symbol')]
  colnames(i)[colnames(i) == 'symbol'] = 'names'
  i$segVal = NA
  i
}
tiled_clades_matrix_pruned_cbind <- tiled_clades_matrix_pruned_cbind[!grepl('chrX', rownames(tiled_clades_matrix_pruned_cbind)),]
tiled_clades_matrix_pruned_cbind <- tiled_clades_matrix_pruned_cbind[!grepl('chrY', rownames(tiled_clades_matrix_pruned_cbind)),]

tiled_clades_matrix_pruned_cbind_granges <- as(add_colnames(data.frame(t(sapply(gsub("chr", "", rownames(tiled_clades_matrix_pruned_cbind)), function(i) strsplit(i, '_')[[1]])))), "GRanges")
tiled_clades_matrix_pruned_cbind_granges
tiled_clades_matrix_pruned_cbind_granges$segVal=NA

ag_subsetchrom_granges <- as(modify_ag_colnames(ag_subsetchrom), 'GRanges')
ag_subsetchrom_granges$segVal = ag_subsetchrom_granges$names
tiled_clades_matrix_pruned_cbind_granges_disjoin <- give_CN_per_gene_v2(
  segment_arg = ag_subsetchrom_granges, ## second in disjoin
  gr_genes = tiled_clades_matrix_pruned_cbind_granges, ## first in disjoin
  only_compute_overlaps = T,
  transform_segment_dataframe=F,
  add_columns=T, value_is_CN=F)
length(unique(tiled_clades_matrix_pruned_cbind_granges_disjoin$CN_val))

## add the cell in which the 
annotation_genes_ecDNA <- tiled_clades_matrix_pruned_cbind_granges_disjoin$CN_val[match(1:length(tiled_clades_matrix_pruned_cbind_granges), tiled_clades_matrix_pruned_cbind_granges_disjoin$revmap)]

length(annotation_genes_ecDNA)
dim(tiled_clades_matrix_pruned_cbind)

## most common genes with ecDNA
sort(table(annotation_genes_ecDNA), decreasing = T)[1:30]

## annotate the genes of the clear ecDNA in chrom 8
table(annotation_genes_ecDNA[gsub("_.*", "", rownames(tiled_clades_matrix_pruned_cbind)) == 'chr8'])
## there isn't MYC??

sum(annotation_genes_ecDNA == 'MYC', na.rm = T)
annotation_genes_ecDNA

tiled_clades_matrix_pruned_cbind[gsub("_.*", "", rownames(tiled_clades_matrix_pruned_cbind)) == 'chr8',]

tiled_clades_matrix_pruned_cbind_genes <- tiled_clades_matrix_pruned_cbind
rownames(tiled_clades_matrix_pruned_cbind_genes) <- make.unique(annotation_genes_ecDNA)
# pdf("../plots/ecDNA_heatmap_chrom8_genes.pdf", width = 9)
# ComplexHeatmap::Heatmap(tiled_clades_matrix_pruned_cbind_genes[(grepl('chr8_', rownames(tiled_clades_matrix_pruned_cbind))),],
#                         top_annotation = annotation_cells_annotation,
#                         show_column_names = F, show_row_names = T, cluster_rows = F)
# dev.off()

give_annotation_genes_of_interest <- function(arg_genes, arg_vec_genes_interest){
  arg_genes[!(arg_genes %in% arg_vec_genes_interest)] = 'other'
  ComplexHeatmap::rowAnnotation(df=(arg_genes),
                                col=list(gene=give_names(viridisLite::magma(length((arg_genes))),
                                                          (arg_genes))))
}


pdf("../plots/ecDNA_heatmap_chrom8_genes.pdf", width = 9)
ComplexHeatmap::Heatmap(tiled_clades_matrix_pruned_cbind_genes[(grepl('chr8_', rownames(tiled_clades_matrix_pruned_cbind))),],
                        top_annotation = annotation_cells_annotation,
                        right_annotation = give_annotation_genes_of_interest(annotation_genes_ecDNA[(grepl('chr8_', rownames(tiled_clades_matrix_pruned_cbind)))], c('MYC')),
                        show_column_names = F, show_row_names = F, cluster_rows = F)
dev.off()

idx_chrom_8_PDO3 <- 2870:2970
pdf("../plots/ecDNA_heatmap_chrom8_genes2.pdf", width = 9)
ComplexHeatmap::Heatmap(tiled_clades_matrix_pruned_cbind_genes[idx_chrom_8_PDO3,grepl('PDO3', annotation_cells$clade)],
                        top_annotation = ComplexHeatmap::HeatmapAnnotation(df=annotation_cells[grepl('PDO3', annotation_cells$clade),],
                                                                           col=list(clade=c(PDO3_A='#e4717a', PDO3_B='#f06949',
                                                                                            PDO3_C='#ff0038'))),
                        right_annotation = give_annotation_genes_of_interest(annotation_genes_ecDNA[idx_chrom_8_PDO3], c('MYC')),
                        show_column_names = F, show_row_names = F, cluster_rows = F)
dev.off()

##-----------------------------------------------------------------------------------##
## Analysis of the MYC ecDNA
## in each cell there might be mutiple ecDNA rings

plot(density(rowSums(tiled_clades_matrix_pruned_cbind_genes[idx_chrom_8_PDO3,grepl('PDO3', annotation_cells$clade)] - 2)))

tiled_clades_matrix_pruned_cbind_genes_region_of_interest_chrom8 <- (tiled_clades_matrix_pruned_cbind_genes[idx_chrom_8_PDO3,grepl('PDO3', annotation_cells$clade)])
plot(sort(rowSums(tiled_clades_matrix_pruned_cbind_genes_region_of_interest_chrom8)))
plot(sort(rowMeans(tiled_clades_matrix_pruned_cbind_genes_region_of_interest_chrom8)))

plot(density(rowMeans(tiled_clades_matrix_pruned_cbind_genes_region_of_interest_chrom8)))

# flexmix::flexmix(data.frame(f=rowMeans(tiled_clades_matrix_pruned_cbind_genes_region_of_interest_chrom8)),
#                  data =  )

plot(tiled_clades_matrix_pruned_cbind_genes_region_of_interest_chrom8[,1], ylim=c(0,140))
abline(v=which(annotation_genes_ecDNA[idx_chrom_8_PDO3] == 'MYC'), col='red')
sapply(seq(10, 120, by = 10), function(i) abline(h=i, type='dashed'))
# points(tiled_clades_matrix_pruned_cbind_genes_region_of_interest_chrom8[,1]-29, ylim=c(0,140), col='blue')
# points(tiled_clades_matrix_pruned_cbind_genes_region_of_interest_chrom8[,1]-54, ylim=c(0,140), col='green')
# points(tiled_clades_matrix_pruned_cbind_genes_region_of_interest_chrom8[,1]-67, ylim=c(0,140), col='green')
# points(tiled_clades_matrix_pruned_cbind_genes_region_of_interest_chrom8[,1]-83, ylim=c(0,140), col='green')
# points(tiled_clades_matrix_pruned_cbind_genes_region_of_interest_chrom8[,1]-94, ylim=c(0,140), col='green')
points(tiled_clades_matrix_pruned_cbind_genes_region_of_interest_chrom8[,1]-103, ylim=c(0,140), col='green')
points(tiled_clades_matrix_pruned_cbind_genes_region_of_interest_chrom8[,1]-117, ylim=c(0,140), col='green')

as.numeric(names(uniq_pddo3_1))-29
uniq_pddo3_1
uniq_pddo3_1 <- table(tiled_clades_matrix_pruned_cbind_genes_region_of_interest_chrom8[,1])
ComplexHeatmap::Heatmap(outer(names(uniq_pddo3_1), names(uniq_pddo3_1), function(i,j) as.numeric(i)-as.numeric(j)))
plot(density((tiled_clades_matrix_pruned_cbind_genes_region_of_interest_chrom8[,1])))
sapply(1:ncol(tiled_clades_matrix_pruned_cbind_genes_region_of_interest_chrom8),
       function(i) lines(density((tiled_clades_matrix_pruned_cbind_genes_region_of_interest_chrom8[,i]))))


plot(tiled_clades_matrix_pruned_cbind_genes_region_of_interest_chrom8[,2], ylim=c(0,140))
plot(tiled_clades_matrix_pruned_cbind_genes_region_of_interest_chrom8[,3], ylim=c(0,140))
plot(tiled_clades_matrix_pruned_cbind_genes_region_of_interest_chrom8[,4], ylim=c(0,140))
plot(tiled_clades_matrix_pruned_cbind_genes_region_of_interest_chrom8[,5], ylim=c(0,140))
plot(tiled_clades_matrix_pruned_cbind_genes_region_of_interest_chrom8[,6], ylim=c(0,140))



# segment_arg = ag_subsetchrom_granges
# gr_genes = tiled_clades_matrix_pruned_cbind_granges
# only_compute_overlaps = T
# transform_segment_dataframe=F
# add_columns=T
# value_is_CN=F
##-----------------------------------------------------------------------------------##

## General landscape of ecDNA, with annotations of genes
dim(tiled_clades_matrix_pruned_cbind_genes)
length(annotation_genes_ecDNA)
annotation_cells$clade

## aggregate by clade
tiled_clades_matrix_pruned_cbind_genes_averaged_clade <- sapply(unique(annotation_cells$clade), function(cl) rowMeans(tiled_clades_matrix_pruned_cbind_genes[,annotation_cells$clade == cl]))
tiled_clades_matrix_pruned_cbind_genes_averaged_clade_melt <- melt(tiled_clades_matrix_pruned_cbind_genes_averaged_clade)
tiled_clades_matrix_pruned_cbind_genes_averaged_clade_melt$x = as.numeric(tiled_clades_matrix_pruned_cbind_genes_averaged_clade_melt$Var1)
tiled_clades_matrix_pruned_cbind_genes_averaged_clade_melt$label = tiled_clades_matrix_pruned_cbind_genes_averaged_clade_melt$Var1
top_ecdna_bin_number <- 1000
tiled_clades_matrix_pruned_cbind_genes_averaged_clade_melt$label[!(tiled_clades_matrix_pruned_cbind_genes_averaged_clade_melt$Var1 %in% rownames(tiled_clades_matrix_pruned_cbind_genes_averaged_clade)[order(rowSums(tiled_clades_matrix_pruned_cbind_genes_averaged_clade), decreasing = T)[1:top_ecdna_bin_number]])] <- NA
tiled_clades_matrix_pruned_cbind_genes_averaged_clade_melt$label2 <- tiled_clades_matrix_pruned_cbind_genes_averaged_clade_melt$label
tiled_clades_matrix_pruned_cbind_genes_averaged_clade_melt$label2[grepl('^NA.', tiled_clades_matrix_pruned_cbind_genes_averaged_clade_melt$label2)] <- NA
unique(tiled_clades_matrix_pruned_cbind_genes_averaged_clade_melt$label2)
tiled_clades_matrix_pruned_cbind_genes_averaged_clade_melt$label2[duplicated(gsub("\\..*", "", tiled_clades_matrix_pruned_cbind_genes_averaged_clade_melt$label2))] = NA

## add chrom and position
tiled_clades_matrix_pruned_cbind_genes_averaged_clade_melt

add_pos <- data.frame(t(sapply(rownames(tiled_clades_matrix_pruned_cbind), function(i) strsplit(i, '_')[[1]]))[tiled_clades_matrix_pruned_cbind_genes_averaged_clade_melt$x,])
add_pos$X2 <- as.numeric(add_pos$X2)
add_pos$X3 <- as.numeric(add_pos$X3)
tiled_clades_matrix_pruned_cbind_genes_averaged_clade_melt = cbind.data.frame(tiled_clades_matrix_pruned_cbind_genes_averaged_clade_melt,
                                                                              add_pos)

tiled_clades_matrix_pruned_cbind_genes_averaged_clade_melt$X1 <- factor(tiled_clades_matrix_pruned_cbind_genes_averaged_clade_melt$X1, levels=gtools::mixedsort(unique(tiled_clades_matrix_pruned_cbind_genes_averaged_clade_melt$X1)))
tiled_clades_matrix_pruned_cbind_genes_averaged_clade_melt <- tiled_clades_matrix_pruned_cbind_genes_averaged_clade_melt[!is.na(tiled_clades_matrix_pruned_cbind_genes_averaged_clade_melt$X1),]
ggplot(tiled_clades_matrix_pruned_cbind_genes_averaged_clade_melt, aes(x=X2, y=value, group=Var2,
                                                                       label=label2))+
  geom_ribbon(aes(ymin=0, ymax=value,  col=Var2,
                  fill=gsub("_.*", "", Var2)), alpha=0.2)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  geom_vline(data=tiled_clades_matrix_pruned_cbind_genes_averaged_clade_melt[!is.na(tiled_clades_matrix_pruned_cbind_genes_averaged_clade_melt$label),],
             aes(xintercept=X2), lty='dashed', alpha=0.2)+
  geom_label_repel(data=tiled_clades_matrix_pruned_cbind_genes_averaged_clade_melt[!is.na(tiled_clades_matrix_pruned_cbind_genes_averaged_clade_melt$label2) & (tiled_clades_matrix_pruned_cbind_genes_averaged_clade_melt$Var2 == 'PDO2_A'),],
                   aes(y=200), max.overlaps = Inf, size=2)+
  theme_bw()+labs(fill='PDO', col='Clade')+
  facet_wrap(.~X1, scales = "free_x")+
  theme(legend.position = "bottom")
ggsave("../plots/ecDNA_clade_average_genome_annotated.pdf", height = 20, width = 15)

give_ribbon_plot_clade_chrom_subset <- function(chrom_vec){
  .tiled_clades_matrix_pruned_cbind_genes_averaged_clade_melt <- tiled_clades_matrix_pruned_cbind_genes_averaged_clade_melt[tiled_clades_matrix_pruned_cbind_genes_averaged_clade_melt$X1 %in% chrom_vec,]
  ggplot(.tiled_clades_matrix_pruned_cbind_genes_averaged_clade_melt, aes(x=X2, y=value, group=Var2,
                                                                         label=label2))+
    geom_ribbon(aes(ymin=0, ymax=value,  col=Var2,
                    fill=gsub("_.*", "", Var2)), alpha=0.2)+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())+
    geom_vline(data=.tiled_clades_matrix_pruned_cbind_genes_averaged_clade_melt[!is.na(.tiled_clades_matrix_pruned_cbind_genes_averaged_clade_melt$label),],
               aes(xintercept=X2), lty='dashed', alpha=0.2)+
    geom_label_repel(data=.tiled_clades_matrix_pruned_cbind_genes_averaged_clade_melt[!is.na(.tiled_clades_matrix_pruned_cbind_genes_averaged_clade_melt$label2) & (.tiled_clades_matrix_pruned_cbind_genes_averaged_clade_melt$Var2 == 'PDO2_A'),],
                     aes(y=200), max.overlaps = Inf, size=2)+
    theme_bw()+labs(fill='PDO', col='Clade')+
    facet_wrap(.~X1, scales = "free_x")+
    theme(legend.position = "bottom")
  
}
give_ribbon_plot_clade_chrom_subset('chr19')

pheatmap(cor(apply(dcast(tiled_clades_matrix_pruned_cbind_genes_averaged_clade_melt[,c('value', 'Var2', 'x')],
      Var2~x), 2, as.numeric)))

## checking possibly interesting genes
tiled_clades_matrix_pruned_cbind_genes_averaged_clade_melt[!is.na(tiled_clades_matrix_pruned_cbind_genes_averaged_clade_melt$label2) & (tiled_clades_matrix_pruned_cbind_genes_averaged_clade_melt$Var2 == 'PDO2_A'),]

##-----------------------------------------------------------------------------------##

##----------------

tiled_ecDNA_table_absCNscDNA_PDO2_cladeA_1_all <- sapply(tiled_ecDNA_table_absCNscDNA_PDO2_cladeA_1, function(i) i[,2])
tiled_ecDNA_table_absCNscDNA_PDO2_cladeA_1_all[is.na(tiled_ecDNA_table_absCNscDNA_PDO2_cladeA_1_all)] <- 0
image(tiled_ecDNA_table_absCNscDNA_PDO2_cladeA_1_all)
tiled_ecDNA_table_absCNscDNA_PDO2_cladeA_1_all_bool_na <- apply(tiled_ecDNA_table_absCNscDNA_PDO2_cladeA_1_all, 2, function(i) i ==0)

image(tiled_ecDNA_table_absCNscDNA_PDO2_cladeA_1_all[rowSums(tiled_ecDNA_table_absCNscDNA_PDO2_cladeA_1_all_bool_na)<ncol(tiled_ecDNA_table_absCNscDNA_PDO2_cladeA_1_all_bool_na),])

example_tile[!is.na(example_tile$gene),]

