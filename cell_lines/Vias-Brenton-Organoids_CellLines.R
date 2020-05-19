#<begin_omit>```{r}
rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
Sys.setenv(LANG='en')
#<end_omit>```{r}

#<begin_chunk>```{r, libraries,message=FALSE, cache=TRUE}
library(grid)
library(gridExtra)
library(dendextend)
library(ggrepel)
library(cowplot)
library(compositions)
library(CompSign)
source("../../../CDA_in_Cancer/code/functions/meretricious/pretty_plots/prettySignatures.R")
prevdata = 'Redefinition'
prevdata = 'NatGen'
#<end_chunk>

#<begin_text>
## The background data that we are using is `r prevdata`.
#<end_text>

#<begin_chunk>```{r,fig.height=4.5,messages=FALSE, echo=FALSE,warning=FALSE, cache=TRUE}
#org <- as(read.csv("data/CIOV_cell_lines_signature_exposures.csv", row.names = 1), 'matrix')
org <- as(read.csv("data/CIOV_cell_lines_signature_exposures_p7.csv", row.names = 1), 'matrix')

createBarplot(as(org, 'matrix'), remove_labels = FALSE, verbatim = FALSE, angle_rotation_axis = 45, order_labels = names(sort(org[,1]))) + 
  ggtitle('Exposures for the organoids')
#<end_chunk>

#<begin_chunk>```{r, natgen_data,include=FALSE, cache=TRUE}
natgen <- list()
natgen_metadata <- list()
## prevdata == 'NatGen'
load("../../../CDA_in_Cancer/data/Robj/image_NatGen_rmd.RData")
natgen0 <- as.matrix(sig_data_unorm[,1:7])
natgen_metadata[[1]] <- sig_data_unorm[,8:ncol(sig_data_unorm)]
## Geoff
## Normalisation is not done in such a way that rows add up to 1. Re-normalising
natgen[[1]] <- sweep(natgen0, 1, rowSums(natgen0), '/')
## last exposures from Ruben

id_previous_samples <- 1
natgen_barplt1 <- createBarplot(natgen[[id_previous_samples]], remove_labels = TRUE, verbatim = FALSE, 
                                order_labels = rownames(natgen[[id_previous_samples]])[(order(natgen[[id_previous_samples]][,1]))]) +
  ggtitle('Original')
natgen_barplt2 <- createBarplot(natgen[[id_previous_samples]], remove_labels = TRUE, verbatim = FALSE, 
                                order_labels = rownames(natgen[[id_previous_samples]])[(order(natgen[[id_previous_samples]][,1]))]) +
  ggtitle('Re-normalised')
#grid.arrange(natgen_barplt1, natgen_barplt2)

# natgen_barplt_perstudy <- list()
# for(i in 1:length(unique(natgen_metadata$study))){
#   natgen_barplt_perstudy[[i]] <- createBarplot(natgen[natgen_metadata$study == unique(natgen_metadata$study)[i],],
#                                                remove_labels = TRUE, verbatim = FALSE)+
#     ggtitle(paste0('Re-normalised\n', unique(natgen_metadata$study)[i] ))
# }
# plot_grid(plotlist=natgen_barplt_perstudy)

## }else if(prevdata == 'Redefinition'){
natgen[[2]] <- readRDS("data/Export-matrix_OV_Sigs_on_TCGA-OV_12112019.rds")
natgen_metadata[[2]] <- data.frame(study=rep('Previous', nrow(natgen[[2]])), stringsAsFactors = FALSE)

#<end_chunk>

#<begin_chunk>```{r, include=FALSE, message=FALSE, cache=TRUE}
org_barplot <- createBarplot(org, remove_labels = FALSE, verbatim = FALSE, angle_rotation_axis = 45, order_labels = names(sort(org[,1]))) + 
  ggtitle('Exposures for the organoids')
no1_natgen1 <- createBarplot(natgen[[1]], remove_labels = TRUE, verbatim = FALSE, 
                                              order_labels = rownames(natgen[[1]])[(order(natgen[[1]][,1]))]) +
  ggtitle('Original')
no1_natgen2 <- createBarplot(natgen[[2]], remove_labels = TRUE, verbatim = FALSE, 
                                              order_labels = rownames(natgen[[2]])[(order(natgen[[2]][,1]))]) +
  ggtitle('Original')
grid.arrange(org_barplot, no1_natgen1, no1_natgen2)
#<end_chunk>

#<begin_text>
##' ## PCA
##' ### PCA in compositional data
##' 
##' In the book Analysing compositional data with R they say that PCA should be done on clr-transformed data.
##' Here I am using robust zeroes: for zero exposures, the centered log-ratios are set to zero (as opposed to -Inf).
##' The plot done with (biplot(princomp(acomp(x)))) is the same as plotting princomp(as(clr(x), 'matrix'))
#<end_text>

#<begin_chunk>```{r, clr_funs,include=FALSE,eval=TRUE, cache=TRUE}
clr_vec <- function(x){
  log(x) - mean(log(x))
}
clr_mat <- function(X){
  .res <- t(apply(X, 1, clr_vec))
  stopifnot(dim(.res) == dim(X))
  .res
}
#<end_chunk>

#<begin_chunk>```{r, clr,include=FALSE,eval=TRUE, cache=TRUE}
  ## there were no zeroes
natgen_clr <- list()
for(i in 1:2){
  cat('Zeroes:',sum(natgen[[i]] == 0),'\n')
#  natgen_clr[[i]] <- clr_mat(natgen[[i]])
  natgen_clr[[i]] <- as(compositions::clr(natgen[[i]]), 'matrix')
}

org_clr <- clr_mat(org)
org_clr_robustzeroes <- as(compositions::clr(org), 'matrix')
rownames(org_clr_robustzeroes) <- rownames(org_clr) <- paste0('Cell line ', rownames(org_clr))
#<end_chunk>


#<begin_chunk>```{r, cols,include=FALSE, cache=TRUE}
n <- 60
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unique(unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))))

#<end_chunk>

#<begin_text>
##' #### Projecting
#<end_text>

#<begin_chunk>```{r, princomp,echo=FALSE, cache=TRUE}

createPCA_projectorganoids <- function(input_matrix, annotation, annotation2, return_df=FALSE, labels_active=TRUE){
  bool_organoids <- grepl('Cell line', rownames(input_matrix))
  prcomp_all <- princomp(input_matrix[!bool_organoids,])
  df_prcomp_exposures <- data.frame(prcomp_all$scores[,1:2], study=annotation[!bool_organoids], labels=NA)
  df_prcomp_exposures_org <- data.frame(predict(prcomp_all, (input_matrix[bool_organoids,]))[,1:2], 'Organoid', rownames(input_matrix[bool_organoids,]))
  colnames(df_prcomp_exposures_org) <- c('Comp.1', 'Comp.2', 'study', 'labels')
  df_prcomp_exposures <- rbind(df_prcomp_exposures, df_prcomp_exposures_org)
  df_prcomp_exposures$labels <- gsub('Sample ', '', df_prcomp_exposures$labels) ##here
  
  myColors <- col_vector[1:length(unique(df_prcomp_exposures$study))]
  names(myColors) <- unique(df_prcomp_exposures$study)
  #ggthemr('flat dark')
  #ggthemr_reset()
  # set_swatch(myColors)
  
  var_explained1 <- (prcomp_all$sdev**2)/sum(prcomp_all$sdev**2)
  
  if(return_df){
    return(prcomp_all)
  }else{
    ggplot(df_prcomp_exposures, aes(x=Comp.1, y=Comp.2, col=study))+
      geom_point() +
      geom_label_repel(aes(label=labels))+
      ggtitle("PCA of both datasets with projection")+
      labs(x=paste0('PC1 (', round(var_explained1[1], 2)*100, '%)'),
           y=paste0('PC2 (', round(var_explained1[2], 2)*100, '%)'))+
      theme_dark()+ theme(legend.position = "bottom")
  }
}
#<end_chunk>


#<begin_chunk>```{r, pca_from_scratch,include=TRUE, echo=FALSE, cache=TRUE}
createPCA_fromscratch <- function(input_matrix, annotation, annotation2, return_df=FALSE, labels_active=TRUE){
  prcomp_all_clr <- princomp(input_matrix)
  
  df_prcomp_exposures_clr <- data.frame(prcomp_all_clr$scores[,1:2],
                                        study=annotation,
                                        bool_any_zeroes=annotation2,
                                        labels=rownames(input_matrix))
  df_prcomp_exposures_clr$labels[!grepl('Cell line', df_prcomp_exposures_clr$labels)] <- NA
  df_prcomp_exposures_clr[,'labels'] <- gsub("Cell line ", "", df_prcomp_exposures_clr$labels)
  var_explained2 <- (prcomp_all_clr$sdev**2)/sum(prcomp_all_clr$sdev**2)
  
  if(return_df){
    prcomp_all_clr
  }else{
    if(labels_active){
      ggplot(df_prcomp_exposures_clr, aes(x=Comp.1, y=Comp.2, col=interaction(bool_any_zeroes, study), label=labels))+ geom_point() + geom_label_repel()+
        labs(x=paste0('PC1 (', round(var_explained2[1], 2)*100, '%)'),
             y=paste0('PC2 (', round(var_explained2[2], 2)*100, '%)'))
    }else{
      ggplot(df_prcomp_exposures_clr, aes(x=Comp.1, y=Comp.2, col=interaction(bool_any_zeroes, study)))+ geom_point() +
        labs(x=paste0('PC1 (', round(var_explained2[1], 2)*100, '%)'),
             y=paste0('PC2 (', round(var_explained2[2], 2)*100, '%)'))
    }
  }
}


for(i in 1:2){
  print(createPCA_projectorganoids(input_matrix = rbind(natgen_clr[[i]],org_clr_robustzeroes),
                                                  annotation = c(natgen_metadata[[i]]$study, rep('Cell Line', nrow(org_clr))),
                                                  annotation2 = c(rep(FALSE, dim(natgen_metadata[[i]])[1]),
                                                                  unlist(apply(org, 1, function(i) any(i == 0)))),
                                                  labels_active = TRUE)+ theme_dark()+ theme(legend.position = "bottom") +
    ggtitle(paste0('PCA of both datasets with projection with robust zeroes, datatset=', i)))
  print(createPCA_fromscratch(input_matrix = rbind(natgen_clr[[i]],org_clr_robustzeroes),
                        annotation = c(natgen_metadata[[i]]$study, rep('Cell Line', nrow(org_clr))),
                        annotation2 = c(rep(FALSE, dim(natgen_metadata[[i]])[1]),
                                        unlist(apply(org, 1, function(i) any(i == 0)))),
                        labels_active = TRUE)+ theme_dark()+ theme(legend.position = "bottom") +
    ggtitle(paste0('PCA created from scratch with robust zeroes, dataset=', i)))
}
#<end_chunk>

#<begin_text>
##' Overall, there are two groups (JB126 and 2259). All JB126 are **extremely** similar except that p22 has a non-zero exposure for 25, whereas all others have a zero exposure. Then, for 2259 there are two groups: p3 and p20 (which have a zero exposure of S6 and a non-zero exposure of S5) and p7 and p13 (opposite scenario; non-zero for S6 and zero for S5).
#<end_text>

#<begin_text>
#' ### Loadings for the PCAs
#<end_text>

#<begin_chunk>```{r, loadings, fig.height=4, echo=FALSE, cache=TRUE}
pcas_with_projection <- list()
pcas_from_scratch <- list()
for(i in 1:2){
  pcas_with_projection[[i]] <- createPCA_projectorganoids(input_matrix = rbind(natgen_clr[[i]],org_clr_robustzeroes),
                                   annotation = c(natgen_metadata[[i]]$study, rep('Cell Line', nrow(org_clr))),
                                   annotation2 = c(rep(FALSE, dim(natgen_metadata[[i]])[1]),
                                                   unlist(apply(org, 1, function(i) any(i == 0)))),
                                   labels_active = TRUE, return_df = TRUE)
  pcas_from_scratch[[i]] <- createPCA_fromscratch(input_matrix = rbind(natgen_clr[[i]],org_clr_robustzeroes),
                              annotation = c(natgen_metadata[[i]]$study, rep('Cell Line', nrow(org_clr))),
                              annotation2 = c(rep(FALSE, dim(natgen_metadata[[i]])[1]),
                                              unlist(apply(org, 1, function(i) any(i == 0)))),
                              labels_active = TRUE, return_df = TRUE)
}

par(mfrow=c(1,2))
for(i in 1:2){
  barplot(pcas_with_projection[[i]]$loadings[,1], main='Loadings of the\nfirst principal component')
  barplot(pcas_with_projection[[i]]$loadings[,2], main='Loadings of the\nfirst principal component')
  barplot(pcas_from_scratch[[i]]$loadings[,1], main='Loadings of the\nsecond principal component')
  barplot(pcas_from_scratch[[i]]$loadings[,2], main='Loadings of the\nsecond principal component')
}
#<end_chunk>


#<begin_chunk>```{r, dendrogram_aitchisondistance,echo=FALSE, cache=TRUE}
par(mfrow=c(1,2))
pdf("results/dendrogram.pdf")
names_prev_datasets <- c('NatGen dataset', 'New OV exposures for SNP TCGA')
for(idx in 1){#1:2){
  organoid_metadata <- cbind.data.frame(study=rep('organoids', nrow(org_clr_robustzeroes)), age=NA, age.cat=NA, stringsAsFactors=FALSE)
  rownames(organoid_metadata) <- rownames(org_clr_robustzeroes)
  if(idx==1){
    all_metadata <- rbind(cbind(natgen_metadata[[idx]]$study), cbind(study=organoid_metadata$study))
  }else{
    all_metadata <- rbind(natgen_metadata[[idx]], cbind(study=organoid_metadata$study))
  }
  all_clr <- rbind(natgen_clr[[idx]], org_clr_robustzeroes)
  rownames(all_metadata) <- rownames(all_clr)
  rm_infinite <- apply(all_clr, 1, function(x) any(is.infinite(x)))
  cat(which(rm_infinite), 'removed due to infinite values')
  all_clr_clean <- all_clr[!rm_infinite,]
  
  which(rm_infinite)
  
  dendro_all <- as.dendrogram(hclust(dist(all_clr_clean)))
  levels_study <- levels(factor(all_metadata[labels(dendro_all),'study']))
  levels_study
  which_level_organoids <- which(grepl('organoids', levels_study))
  cols <- rep(NA, length(levels_study))
  cols[which_level_organoids] <- 'blue' #'#88E9A2'
  cols[-which_level_organoids] <- c('#FFA07A', '#FA8072', '#E9967A', '#F08080')
  labels_colors(dendro_all) <- cols[factor(all_metadata[labels(dendro_all),'study'])]
  labels_org_bool <- labels_colors(dendro_all) == 'blue' #'#88E9A2'
  # labels(dendro_all)[labels_org_bool] <- rep('●', sum(labels_org_bool))
  labels(dendro_all)[!labels_org_bool] <- rep('•', sum(!labels_org_bool))
  labels(dendro_all)[!labels_org_bool] <- rep(NA, sum(!labels_org_bool))
  labels(dendro_all) <- gsub('Cell line ', '', labels(dendro_all))
  cex_labels <- rep(1, length(labels_org_bool))
  cex_labels[labels_org_bool] <- 0.9
  dendro_all <- set(dendro_all, "labels_cex", cex_labels)
  plot(dendro_all, cex=0.2, cex.main=1, main=paste0('Dendrogram based on the exposures\n(Aitchison distance)\n', names_prev_datasets[idx]))
}
dev.off()
#<end_chunk>
