
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
library(umap)
library(gridExtra)
library(jcolors)

set.seed(1325)

if(load_data){
  load("../objects/umap_image.RData")
}else{
  
  # load("~/Documents/PhD/other_repos/robjects/umap_TCGA_HGSOC.RData")
  
  
  ## Load file output from DESeq2, which has been created with 1_run_DE.sh
  load("../objects/deaObjectFile")
  deObj = `~response`
  
  ## results from DE analysis between resistant and sensitives
  results <- readRDS("../objects/resultsDESeq_TCGA.RDS")
  rownames_short = sapply(rownames(results), function(i) strsplit(i, '[.]')[[1]][1])
  
  ## read GRCh38 conversion
  t2g = readRDS("../../copy_number_analysis_organoids/robjects/t2g2.RDS")
  coding_genes = readRDS("../../copy_number_analysis_organoids/robjects/coding_genes.RDS")
  gene_conversion = t2g[match(rownames_short, t2g$ensembl_gene_id),]
  dim(gene_conversion)
  dim(results)
  
  query_GDC <- readRDS("../objects/query")
  
  results = cbind(gene_conversion, results)
  results = results[order(results$log2FoldChange),]
  
  counts_DESeq_TCGA = DESeq2::counts(deObj, normalized=TRUE)
  # counts_DESeq_TCGA <- counts_DESeq_TCGA[!grepl('__', rownames(counts_DESeq_TCGA)), ]
  counts_DESeq_TCGA_raw = DESeq2::counts(deObj, normalized=FALSE)
  # counts_DESeq_TCGA_raw <- counts_DESeq_TCGA_raw[!grepl('__', rownames(counts_DESeq_TCGA_raw)), ]
  
  
  ## Clustering of TCGA
  
  umap_TCGA <- umap(counts_DESeq_TCGA)
  umap_2 <- umap(t(counts_DESeq_TCGA))
  
  annotation_umap <- cbind.data.frame(Gene=rownames(counts_DESeq_TCGA))
  annotation_umap$replication <- sapply(coding_genes$description, grepl, pattern = 'replication' )[match(rownames(counts_DESeq_TCGA), coding_genes$ensembl_gene_id)]
  # annotation_umap$consensusTMB <- (coding_genes$external_gene_name[match(rownames(counts_DESeq_TCGA), coding_genes$ensembl_gene_id)] %in% unique(unlist(ConsensusTMB_OV)))[match(rownames(counts_DESeq_TCGA), coding_genes$ensembl_gene_id)]
  # plot_umap(umap_TCGA, factor(annotation_umap$replication, levels=c('FALSE', 'TRUE')))
  # plot_umap(umap_TCGA, factor(annotation_umap$consensusTMB))
  
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
  change_names <- function(i){
    names(i) = t2g$external_gene_name[match(gsub("\\..*","",names(i)), t2g$ensembl_gene_id)]
    i
  }
  matched_consensusTME <- ConsensusTME::consensusTMEAnalysis(change_rownames(counts_DESeq_TCGA_raw), cancerType = "OV")
  
  ## clinical ov data
  TCGA_genes <- readRDS("../../../../other_repos/cnsigs_Initial_submission/survival_analysis/from_ruben/survival_models/TCGA_OVBRCAonly_Exposures_and_BRCA_Status_plusGene.rds")
  cluster_fig1 <- readRDS("../../copy_number_analysis_organoids/robjects/dendrograminputclr_tree.RDS")
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
  
  brca <- readRDS("../../../cnsigs_Initial_submission/survival_analysis/from_ruben/survival_models/TCGA_OVBRCAonly_Exposures_and_BRCA_Status_plusGene.rds")
  
  df_umap_2 <- data.frame(umap_2$layout, mean_exprs=log(colMeans(counts_DESeq_TCGA)),
                          matched_clinical,WGD=matched_wgs,matched_queryGDC,
                          matched_exposures,
                          group_clr=matched_group_WGD,
                          genes=matched_group_genes,
                          consensusTME=t(matched_consensusTME))
  colnames(df_umap_2)
  df_umap_2$brca_status = brca$Status[match(df_umap_2$case_submitter_id, brca$Sample)]
  
  
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
  ggplot(df_umap_2,
         aes(x=X1, y=X2, col=brca_status))+geom_point() ## no separation by brca status
  
  
  ggplot(cbind(df_umap_2, gene=counts_DESeq_TCGA[grepl('ENSG00000136997', rownames(counts_DESeq_TCGA)),]),
         aes(x=X1, y=X2, col=gene))+geom_point() ## not separated by MYC
  
  ggplot(cbind(df_umap_2, gene=counts_DESeq_TCGA[grepl('ENSG00000175054', rownames(counts_DESeq_TCGA)),]),
         aes(x=X1, y=X2, col=gene))+geom_point() ## not separated by ATR
  ggplot(cbind(df_umap_2, gene=counts_DESeq_TCGA[grepl('ENSG00000095585', rownames(counts_DESeq_TCGA)),]),
         aes(x=X1, y=X2, col=gene))+geom_point() ## not separated by BLNK
  
  to_ens <- function(i){
    t2g$ensembl_gene_id[match(i, t2g$external_gene_name)]
  }
  ggplot(cbind(df_umap_2, gene=counts_DESeq_TCGA[grepl(to_ens('ESR1'), rownames(counts_DESeq_TCGA)),]),
         aes(x=X1, y=X2, col=gene))+geom_point() ## not separated by ER
  ggplot(cbind(df_umap_2, gene=counts_DESeq_TCGA[grepl(to_ens('PGR'), rownames(counts_DESeq_TCGA)),]),
         aes(x=X1, y=X2, col=gene))+geom_point() ## not separated by PR
  ggplot(cbind(df_umap_2, gene=counts_DESeq_TCGA[grepl(to_ens('RECQL'), rownames(counts_DESeq_TCGA)),]),
         aes(x=X1, y=X2, col=log(gene)))+geom_point() ## not separated by RECQL
  
  pdf("../figures/figures_umap/genes_gradient_umap.pdf", width = 10, height = 3)
  grid.arrange(ggplot(cbind(df_umap_2, gene=counts_DESeq_TCGA[grepl(to_ens('RNVU1-7'), rownames(counts_DESeq_TCGA)),]),
         aes(x=X1, y=X2, col=log(gene)))+geom_point()+labs(col='RNVU1-7')+scale_colour_jcolors_contin("pal3")+theme(legend.position = 'bottom'),
  ggplot(cbind(df_umap_2, gene=counts_DESeq_TCGA[grepl(to_ens('RNVU1-18'), rownames(counts_DESeq_TCGA)),]),
         aes(x=X1, y=X2, col=log(gene)))+geom_point()+labs(col='RNVU1-18')+scale_colour_jcolors_contin("pal3")+theme(legend.position = 'bottom'),
  ggplot(cbind(df_umap_2, gene=counts_DESeq_TCGA[grepl(to_ens('AL355075.4'), rownames(counts_DESeq_TCGA)),]),
         aes(x=X1, y=X2, col=log(gene)))+geom_point()+labs(col='AL355075.4')+scale_colour_jcolors_contin("pal3")+theme(legend.position = 'bottom'),
  ggplot(cbind(df_umap_2, gene=counts_DESeq_TCGA[grepl(to_ens('AL355075.4'), rownames(counts_DESeq_TCGA)),]),
         aes(x=X1, y=X2, col=log(gene)))+geom_point()+labs(col='AL355075.4')+scale_colour_jcolors_contin("pal3")+theme(legend.position = 'bottom'), nrow=1)
  dev.off()
  
  ggplot(cbind(df_umap_2),
         aes(x=X1, y=X2, col=s1))+geom_point() ## not separated by s1
  ggplot(cbind(df_umap_2),
         aes(x=X1, y=X2, col=s3))+geom_point() ## not separated by s3
  ggplot(cbind(df_umap_2),
         aes(x=X1, y=X2, col=s4))+geom_point() ## not separated by s4
  
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
      ggrepel::geom_label_repel(data=dfSigInlier, aes(x=  .data$vx, y = .data$vy, label=labls), max.overlaps = 40)+
      geom_point(data = dfSigOutlier,  aes(x = .data$vx, y = .data$vy),
                 colour = colOut, shape = dfSigOutlier$shape) + 
      ggrepel::geom_label_repel(data=dfSigOutlier, aes(x=  .data$vx, y = .data$vy, label=labls), max.overlaps = 40)+
      theme(legend.position = "none") + ggtitle("Volcano Plot") + 
      labs(x = "Log of Fold Change", y = "-log10 of FDR")+
      ylim(c(0, 1.8))
    return(p)
  }
  
  pdf("../figures/figures_umap/volcano_plot_branches.pdf")
  vlcano(results_DE_umap)
  dev.off()
  
  table(results_DE_umap$padj < 0.00000001)
  xtable::print.xtable(xtable::xtable(print(head(results_DE_umap[order(results_DE_umap$padj, decreasing = F),c('external_gene_name', 'log2FoldChange', 'padj')]), row.names = F)), row.names = F)
  top_umap_DE <- change_rownames((DESeq2::counts(mat, normalized=F)[which(results_DE_umap$padj < 0.000000001),]))
  pheatmap::pheatmap(top_umap_DE,
                     annotation_col = groups_umap, cluster_cols = T, scale = "row",
                     main = 'Genes with adjusted p-value < 0.05', show_colnames = FALSE)
  
  top_umap_DE2 <- change_rownames(head(DESeq2::counts(mat, normalized=T)[order(results_DE_umap$padj),], n=35))
  
  pdf("../figures/figures_umap/top_DE_genes_branches.pdf", height = 5)
  print(pheatmap::pheatmap(log(top_umap_DE2+0.001),
                     annotation_col = groups_umap, cluster_cols = F, scale = "row",
                     main = 'Genes with adjusted p-value < 0.05 (top 35)', show_colnames = FALSE))
  dev.off()
  
  change_rownames2 <- function(i){
    rownames(i) <- gsub("\\..*","",rownames(i))
    i
  }
  
  top_umap_DE3 <- change_rownames(change_rownames2(DESeq2::counts(mat, normalized=T))[match((head(results_DE_umap[results_DE_umap$padj < 0.05,'ensembl_gene_id'], n=100)), gsub("\\..*","",rownames(mat))),])
  top_umap_DE3 <- top_umap_DE3[!is.na(rownames(top_umap_DE3)),]
  top_umap_DE3_cor <- cor(t(top_umap_DE3))
  
  pdf("../figures/figures_umap/top_DE_genes_branches_cor.pdf", height = 11, width=11)
  print(pheatmap::pheatmap(top_umap_DE3_cor))
  dev.off()
  
  
  ### Analysis of bottom group
  rownames(groups_umap)[groups_umap$grouping_umap == 'Group2']
  
  ## for each gene, do the correlation of umap1 dimension, and the expression of the gene
  
  gradient_per_group <- lapply(c('Group1', 'Group2'), function(group_it){
    counts_group2 <- DESeq2::counts(mat, normalized=T)
    counts_group2 <- counts_group2[,rownames(groups_umap)[groups_umap$grouping_umap == group_it]]
    
    counts_group2 <- counts_group2[,match(names(sort(umap_2$layout[,1])), colnames(counts_group2))]
    counts_group2 <- counts_group2[,!is.na(colnames(counts_group2))]
    num_samples_group2 <- ncol(counts_group2)
    cor_gradient2 <- apply(counts_group2, 1, function(i) cor(x = 1:num_samples_group2, y=i,method = "kendall"))
    ## check the top correlated genes
    top_cor_bottombranch <- change_names(head(sort(cor_gradient2, decreasing = T), n=20))
    return(list(counts_group2, cor_gradient2, top_cor_bottombranch))
  })
  names(gradient_per_group) <- c('Group1', 'Group2')
  
  counts_group1 <- gradient_per_group$Group1[[1]]
  cor_gradient1 <- gradient_per_group$Group1[[2]]
  top_cor_upperbranch <- gradient_per_group$Group1[[3]]
  
  counts_group2 <- gradient_per_group$Group2[[1]]
  cor_gradient2 <- gradient_per_group$Group2[[2]]
  top_cor_bottombranch <- gradient_per_group$Group2[[3]]
  
  plot(density(cor_gradient2[!is.na(cor_gradient2)]))
  
  norm_DESeq2_counts <- change_rownames(DESeq2::counts(mat, normalized=T))
  norm_DESeq2_counts_rowmean <- rowMeans(norm_DESeq2_counts)
  
  ggplot(cbind(df_umap_2, gene=counts_DESeq_TCGA[grepl(to_ens('RN7SKP230'), rownames(counts_DESeq_TCGA)),]),
         aes(x=X1, y=X2, col=log(gene)))+geom_point()+labs(col='RN7SKP230')+scale_colour_jcolors_contin("pal3")+theme(legend.position = 'bottom')
  
  ggplot(cbind(df_umap_2, gene=counts_DESeq_TCGA[grepl(to_ens('AL355075.4'), rownames(counts_DESeq_TCGA)),]),
         aes(x=X1, y=X2, col=log(gene)))+geom_point()+labs(col='AL355075.4')+scale_colour_jcolors_contin("pal3")+theme(legend.position = 'bottom')
  
  ## compare the expression of these highly correlated genes to the average expression of all genes
  df_bottombranch <- cbind.data.frame(name=names(norm_DESeq2_counts_rowmean),
                   averagexprs=norm_DESeq2_counts_rowmean,
                   x2=ifelse(test = names(norm_DESeq2_counts_rowmean) %in% names(top_cor_bottombranch)[!is.na(names(top_cor_bottombranch))],
                             yes = norm_DESeq2_counts_rowmean, no = NA),
                   label=ifelse(test = names(norm_DESeq2_counts_rowmean) %in% names(top_cor_bottombranch)[!is.na(names(top_cor_bottombranch))],
                          yes = names(norm_DESeq2_counts_rowmean), no = NA))
  ggplot(df_bottombranch,
         aes(x=averagexprs, label=label))+
    geom_vline(aes(xintercept = x2), data = df_bottombranch)+
    geom_label_repel(aes(y=0, col=label), max.overlaps = 30)+
    scale_x_continuous(trans = "log2")+
    geom_density()
  ## they don't have particularly high exposure  
  
  xtable::xtable(cbind(Correlation=top_cor_bottombranch))
  xtable::xtable(cbind(Correlation=top_cor_upperbranch))
  # save.image("~/Documents/PhD/other_repos/robjects/umap_TCGA_HGSOC.RData")
  
  umap_TPM <- umap(t(sweep(counts_DESeq_TCGA_raw, 2, colSums(counts_DESeq_TCGA_raw), '/')*1e6))
  plot(umap_TPM$layout[,1:2])
  
  df_umap_TPM <- data.frame(umap_TPM$layout, mean_exprs=log(colMeans(counts_DESeq_TCGA)),
                          matched_clinical,WGD=matched_wgs,matched_queryGDC,
                          matched_exposures,
                          group_clr=matched_group_WGD,
                          genes=matched_group_genes,
                          consensusTME=t(matched_consensusTME),
                          brca_status=brca$Status[match(df_umap_2$case_submitter_id, brca$Sample)])
  
  ggplot(cbind(df_umap_TPM, gene=counts_DESeq_TCGA[grepl(to_ens('RN7SKP230'), rownames(counts_DESeq_TCGA)),]),
         aes(x=X1, y=X2, col=log(gene)))+geom_point()+labs(col='RN7SKP230')+scale_colour_jcolors_contin("pal3")+theme(legend.position = 'bottom')
  ggplot(df_umap_TPM,
         aes(x=X1, y=X2, col=brca_status))+geom_point()+
    theme(legend.position = 'bottom')
  ggplot(cbind(df_umap_TPM, gene=counts_DESeq_TCGA[grepl(to_ens('AL355075.4'), rownames(counts_DESeq_TCGA)),]),
         aes(x=X1, y=X2, col=log(gene)))+geom_point(size=1)+labs(col='AL355075.4')+
    scale_colour_jcolors_contin("pal3")+theme(legend.position = 'bottom')
  ggsave("../figures/figures_umap/genes_gradient_umapTPM_AL355075_4.pdf", width = 2.5, height = 3)
  
  pdf("../figures/figures_umap/genes_gradient_umapTPM.pdf", width = 10, height = 3)
  grid.arrange(ggplot(cbind(df_umap_TPM, gene=counts_DESeq_TCGA[grepl(to_ens('RNVU1-7'), rownames(counts_DESeq_TCGA)),]),
                      aes(x=X1, y=X2, col=log(gene)))+geom_point()+labs(col='RNVU1-7')+scale_colour_jcolors_contin("pal3")+theme(legend.position = 'bottom'),
               ggplot(cbind(df_umap_TPM, gene=counts_DESeq_TCGA[grepl(to_ens('RNVU1-18'), rownames(counts_DESeq_TCGA)),]),
                      aes(x=X1, y=X2, col=log(gene)))+geom_point()+labs(col='RNVU1-18')+scale_colour_jcolors_contin("pal3")+theme(legend.position = 'bottom'),
               ggplot(cbind(df_umap_TPM, gene=counts_DESeq_TCGA[grepl(to_ens('AL355075.4'), rownames(counts_DESeq_TCGA)),]),
                      aes(x=X1, y=X2, col=log(gene)))+geom_point()+labs(col='AL355075.4')+scale_colour_jcolors_contin("pal3")+theme(legend.position = 'bottom'),
               ggplot(cbind(df_umap_TPM, gene=counts_DESeq_TCGA[grepl(to_ens('AL355075.4'), rownames(counts_DESeq_TCGA)),]),
                      aes(x=X1, y=X2, col=log(gene)))+geom_point()+labs(col='AL355075.4')+scale_colour_jcolors_contin("pal3")+theme(legend.position = 'bottom'), nrow=1)
  dev.off()
  
  ggplot(cbind(df_umap_2, gene=counts_DESeq_TCGA[grepl(to_ens('AL355075.4'), rownames(counts_DESeq_TCGA)),]),
         aes(x=X1, y=X2, col=log(gene)))+geom_point()+labs(col='AL355075.4')+scale_colour_jcolors_contin("pal3")+theme(legend.position = 'bottom')
  
  ## what is the top correlated transcript, that is not annotated?
  which.max(cor_gradient2)
  
  
  ### Repeated elements?
  counts_DESeq_TCGA_all = DESeq2::counts(deObj, normalized=TRUE)
  repeated_or_problematic <- counts_DESeq_TCGA_all[grepl('__', rownames(counts_DESeq_TCGA_all)),]
  pairs(t(repeated_or_problematic))
  ggplot(melt(repeated_or_problematic['__alignment_not_unique',]),
         aes(x=value))+geom_density()#+scale_x_continuous(trans = "log2")
  
  ##Gene: RNVU1-7 ENSG00000206585
  
  ggplot(cbind.data.frame(alignment_not_unique=counts_DESeq_TCGA_all['__alignment_not_unique',],
                          RNVU1_7=counts_DESeq_TCGA_all[grepl('ENSG00000206585',
                              rownames(counts_DESeq_TCGA_all)),]),
         aes(x=alignment_not_unique, y=RNVU1_7))+
    geom_smooth()+
    geom_point()+scale_x_continuous(trans = "log2")+scale_y_continuous(trans = "log2")+
    theme_bw()
  ggsave("../figures/figures_umap/alignment_not_unique_and_RNVU1_7.pdf", width = 5, height = 5)
  
  save.image("../objects/umap_image.RData")
}

