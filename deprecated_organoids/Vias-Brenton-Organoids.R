#<begin_text>
##'
##' - **Important!** I believe from Geoff that Ruben has updated exposures for TCGA, BRITROC, etc. samples
##'
##' - How do organoids from patients reflect the generality of OV patients?
##'
##' - Given signatures of CN, can we say something about how this compares to the TCGA cohort?
##' 
##' Some sort of clustering of patients and organoids and see where organoids fall -- based on signatures, and also based on raw CN profiles.
##' 
##' **To do**: for (1) signatures, and (2) CN profiles, do (1) clustering, (2) PCA (i.e. 4 plots).
##' 
##' (Get TCGA, ICGC, BRITROC from Geoff, and compare this to organoids).
##' 
##' Multisample OV: no difference of note.
##' 
##' *note* in latest version of signatures s1 has been removed. Additionally, now the signatures I got are not
##' normalised.
##' 
##' 
#<end_text>

#<begin_omit>```{r}
rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
Sys.setenv(LANG='en')
#<end_omit>```{r}

#<begin_chunk>```{r, libraries,message=FALSE}
library(grid)
library(gridExtra)
library(dendextend)
library(ggrepel)
library(cowplot)
library(compositions)
library(CompSign)
#<end_chunk>

#<begin_text>
##' ## Loading data
##' ### Data from organoids
#<end_text>

#<begin_chunk>```{r,fig.height=2.5,message=FALSE, echo=FALSE,warning=FALSE}
#x <- read.table("data/organoids_signature_exposuresDom.csv", sep = ',', header = T)
first_version <- FALSE
x <- readRDS("data/organoid_exposures_Aug21.rds")
source("../../../CDA_in_Cancer/code/functions/meretricious/pretty_plots/prettySignatures.R")

if(first_version){
  org <- x[,-1]
}else{
  org <- x
}
org <- as(org, 'matrix')
rownames(org) <- paste0('Sample ', 1:nrow(org))

createBarplot(org, remove_labels = TRUE, order_labels = names(sort(org[,1]))) + 
  ggtitle('Exposures for the organoids')


org_nonnormalised <- org
org <- sweep(org, 1, rowSums(org), '/')
createBarplot(org, remove_labels = TRUE, order_labels = names(sort(org[,1]))) + 
  ggtitle('Exposures for the organoids')
#<end_chunk>

#<begin_text>
##' ### Data from Nat Gen paper
#<end_text>

#<begin_chunk>```{r, natgen_data,include=FALSE}
load("../../../CDA_in_Cancer/data/Robj/image_NatGen_rmd.RData")
natgen <- as.matrix(sig_data_unorm[,1:7])
natgen_metadata <- sig_data_unorm[,8:ncol(sig_data_unorm)]
## Normalisation is not done in such a way that rows add up to 1. Re-normalising
natgen_renormalised <- sweep(natgen, 1, rowSums(natgen), '/')

natgen_barplt1 <- createBarplot(natgen, remove_labels = TRUE,
                                order_labels = rownames(natgen)[(order(natgen[,1]))]) +
  ggtitle('Original')
natgen_barplt2 <- createBarplot(natgen_renormalised, remove_labels = TRUE,
                                order_labels = rownames(natgen)[(order(natgen[,1]))]) +
  ggtitle('Re-normalised')
grid.arrange(natgen_barplt1, natgen_barplt2)
#<end_chunk>

#<begin_chunk>```{r, barplots_per_study,echo=FALSE,message=FALSE,warning=FALSE}
natgen_barplt_perstudy <- list()
for(i in 1:length(unique(natgen_metadata$study))){
  natgen_barplt_perstudy[[i]] <- createBarplot(natgen_renormalised[natgen_metadata$study == unique(natgen_metadata$study)[i],],
                                  remove_labels = TRUE)+
    ggtitle(paste0('Re-normalised\n', unique(natgen_metadata$study)[i] ))
}
plot_grid(plotlist=natgen_barplt_perstudy)
#<end_chunk>


#<begin_text>
##' ## PCA
##' ### PCA in compositional data
##' 
##' In the book Analysing compositional data with R they say that PCA should be done on clr-transformed data.
##' Zeroes are an issue (see last section)
#<end_text>

#<begin_chunk>
org_barplot <- createBarplot(org, remove_labels = TRUE, order_labels = names(sort(org[,1]))) + 
  ggtitle('Exposures for the organoids')
no1_natgen <- natgen_barplt1 <- createBarplot(natgen_renormalised, remove_labels = TRUE,
                                order_labels = rownames(natgen)[(order(natgen[,1]))]) +
  ggtitle('Original')
grid.arrange(org_barplot, no1_natgen)
#<end_chunk>


#<begin_text>
##' ### PCA based on the signatures
##' The plot done with (biplot(princomp(acomp(x)))) is the same as plotting princomp(as(clr(x), 'matrix'))
#<end_text>

#<begin_chunk>```{r, clr,include=FALSE,eval=TRUE}
clr_vec <- function(x){
  log(x) - mean(log(x))
}
clr_mat <- function(X){
  .res <- t(apply(X, 1, clr_vec))
  stopifnot(dim(.res) == dim(X))
  .res
}

natgen_clr <- clr_mat(natgen_renormalised)
org_clr <- clr_mat(org)
org_clr_robustzeroes <- as(compositions::clr(org), 'matrix')
rownames(org_clr_robustzeroes) <- rownames(org_clr) <- paste0('Organoid ', rownames(org_clr))
#<end_chunk>


#<begin_chunk>```{r, cols,include=FALSE}
n <- 60
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unique(unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))))

#<end_chunk>

#<begin_text>
##' #### Projecting
#<end_text>

#<begin_chunk>```{r, princomp,echo=FALSE}
prcomp_all <- princomp(as(clr(natgen_renormalised), 'matrix'))
#plot(prcomp_all)

df_prcomp_exposures <- data.frame(prcomp_all$scores[,1:2], study=natgen_metadata$study, labels=NA)
df_prcomp_exposures_org <- data.frame(predict(prcomp_all, (org_clr_robustzeroes))[,1:2], 'Organoid', rownames(org))
colnames(df_prcomp_exposures_org) <- c('Comp.1', 'Comp.2', 'study', 'labels')
df_prcomp_exposures <- rbind(df_prcomp_exposures, df_prcomp_exposures_org)
df_prcomp_exposures$labels <- gsub('Sample ', '', df_prcomp_exposures$labels) ##here

myColors <- col_vector[1:length(unique(df_prcomp_exposures$study))]
names(myColors) <- unique(df_prcomp_exposures$study)
#ggthemr('flat dark')
#ggthemr_reset()
# set_swatch(myColors)

var_explained1 <- (prcomp_all$sdev**2)/sum(prcomp_all$sdev**2)

ggplot(df_prcomp_exposures, aes(x=Comp.1, y=Comp.2, col=study))+
  geom_point() +
  geom_label_repel(aes(label=labels))+
  ggtitle("PCA of both datasets with projection")+
  labs(x=paste0('PC1 (', round(var_explained1[1], 2)*100, '%)'),
       y=paste0('PC2 (', round(var_explained1[2], 2)*100, '%)'))+
  theme_dark()

#  scale_color_discrete(myColors)

# plot(prcomp_all$scores[,1:2], col=factor(natgen_metadata$study), pch=4)
# points(t(prcomp_all$loadings %*% t(org))[,1:2], pch=19, col='red')
## just to point out that they are the same

#<end_chunk>

#<begin_text>
##' #### What is different in these 'underrepresented' samples?
##' Conclusion: it seems as though it's signature 3, the relative abundance of which is never high in organoid samples.
##' I am comparing
##' - the barplots of the exposures
##' - CLR (centered log-ratio) of signature 3 is high in the underrepresented samples
##' - the ratio of the sums of different signatures, e.g. the ratio of 1+3+5 vs 2+4+6+7.
##' - ILR (isometric log-ratio) when splitting the dataset into s3 and all other signatures. It is the log-ratio of the exposure to signature 3 and the geometric mean of all other exposures.
#<end_text>

#<begin_chunk>```{r, underrepresented, fig.height=6, echo=FALSE, eval=TRUE}
selected_underrepresented_right <- natgen_renormalised[which(df_prcomp_exposures$Comp.1 > max(df_prcomp_exposures[df_prcomp_exposures$study == 'Organoid','Comp.1'])),]
selected_underrepresented_left <- natgen_renormalised[which(df_prcomp_exposures$Comp.1 < min(df_prcomp_exposures[df_prcomp_exposures$study == 'Organoid','Comp.1'])),]
grid.arrange(createBarplot(as(selected_underrepresented_right, 'matrix'), remove_labels = TRUE)+ggtitle('Underrepresented (right)'),
             createBarplot(as(selected_underrepresented_left, 'matrix'), remove_labels = TRUE)+ggtitle('Underrepresented (left)'),
             createBarplot(org, remove_labels = TRUE)+ggtitle('Organoids'))
dev.off()


#<end_chunk>

#<begin_text>
#' The *ALR* is the additive log-ratio transformation, where we divide each element by one of the other elements, and then take the log
#' (equivalently, subtract the log of one of the parts to all other logs):
#' $$alr(x)_n= \Bigg(\log\Big(\frac{x_1}{x_n}\Big),\ \log\Big(\frac{x_2}{x_n}\Big), \, ..., \ \log\Big(\frac{x_{n-1}}{x_n}\Big)\Bigg) $$
#'Boxplots for the alr of 
#' (1) right underrepresented samples
#' (2) left underrepresented samples
#' (3) all samples (i.e., non-organoid) including those in (1) and (2)

#<end_text>

#<begin_chunk>```{r, underrepresented_both_boxplot, fig.height=6, echo=FALSE, eval=TRUE}
# par(mfrow=c(2,3))
# for(i in 2:7){
#   boxplot(log(selected_underrepresented_right[,i]/selected_underrepresented_right[,1]),
#           log(selected_underrepresented_left[,i]/selected_underrepresented_left[,1]),
#           log(natgen_renormalised[,i]/natgen_renormalised[,1]), names = c('underrepresented right',
#                                                                           'underrepresented left',
#                                                                           'all Natgen'),
#           main=paste0('ALR signature ', i, ' (wrt signature 1)'))
# }

melt_alr <- melt(lapply(2:7, function(i) list(underrepresented_right=log(selected_underrepresented_right[,i]/selected_underrepresented_right[,1]),
                                  underrepresented_left=log(selected_underrepresented_left[,i]/selected_underrepresented_left[,1]),
                                  all_Natgen=log(natgen_renormalised[,i]/natgen_renormalised[,1]),
                                  organoids=tomatrix(alr(org, ivar = 1))[,i-1])))
melt_alr$L2 <- factor(melt_alr$L2,  levels=c('underrepresented_right',
                     'underrepresented_left',
                     'all_Natgen', 'organoids'))
melt_alr$L1 <- melt_alr$L1+1
ggplot(melt_alr, aes(x=L2, y=value))+geom_boxplot()+facet_wrap(.~L1)+theme_light()+
  theme(axis.text.x = element_text(angle = 30, hjust = 1))+ ggtitle('ALR with respect to signature 1')

#<end_chunk>

#<begin_text>
#' This is with the *CLR* (centered log-ratio transformation):
#' #' $$clr(x)_n= \Bigg(\log\Big(\frac{x_1}{\textrm{gm}(x)}\Big),\ \log\Big(\frac{x_2}{\textrm{gm}(x)}\Big), \, ..., \ \log\Big(\frac{x_n}{\textrm{gm}(x)}\Big)\Bigg) $$
#' where $\textrm{gm}(x)$ is the geometric mean.
#' Same order:
#' (1) right underrepresented samples
#' (2) left underrepresented samples
#' (3) all samples (i.e., non-organoid) including those in (1) and (2)
#<end_text>

#<begin_chunk>```{r, underrepresented2, fig.height=6, echo=FALSE, eval=TRUE}
## is it signature 3?
# par(mfrow=c(2,4))
# for(i in 1:7){
#   boxplot(clr_mat(selected_underrepresented_right)[,i],
#           clr_mat(selected_underrepresented_left)[,i],
#           org_clr_robustzeroes, main=paste0('CLR of signature ', i))
# }
melt_clr <- melt(lapply(1:7, function(i) list(underrepresented_right=clr_mat(selected_underrepresented_right)[,i],
                                              underrepresented_left=clr_mat(selected_underrepresented_left)[,i],
                                              all_Natgen=natgen_clr[,i],
                                              organoids=org_clr_robustzeroes[,i])))
melt_clr[melt_clr$L1 %in% c(4, 6, 7),'grouping'] <- 'group 1'
melt_clr[melt_clr$L1 %in% c(1, 3, 5),'grouping'] <- 'group 2'
melt_clr$L2 <- factor(melt_clr$L2,  levels=c('underrepresented_right',
                                             'underrepresented_left',
                                             'all_Natgen', 'organoids'))
ggplot(melt_clr, aes(x=L2, y=value, col=grouping))+geom_boxplot()+facet_wrap(.~L1)+theme_light()+
  theme(axis.text.x = element_text(angle = 30, hjust = 1))+ ggtitle('CLR for each of the signatures')
ggplot(melt_clr, aes(x=L2, y=value, col=grouping))+geom_boxplot()+facet_wrap(.~L1, nrow=1)+theme_light()+
  theme(axis.text.x = element_text(angle = -30, hjust = 0))+ ggtitle('CLR for each of the signatures')+labs(x="")#+ guides(col=FALSE)
ggsave("../../First-Year-Report/FYR2/figures/clr_organoids.pdf", height = 2.7, width = 8)
#<end_chunk>

plot(prcomp(ilr(org))$x[,1:2])

#<begin_text>
#' The colouring is based on some apparent patterns.
#' 
#' 
#' ### Loadings
#' Looking at the loadings. In particular, looking for components in the first and second PC
#<end_text>

#<begin_chunk>```{r, loadings, fig.height=8, echo=FALSE}
par(mfrow=c(2,1))
barplot(prcomp_all$loadings[,1], main='Loadings of the\nfirst principal component')
barplot(prcomp_all$loadings[,2], main='Loadings of the\nsecond principal component')
#<end_chunk>

#<begin_text>
#' Signatures 3 and 6 seem to be quite important for the underrepresented groups
#' 
#' ### Amalgamation
#' Here I am amalgamating, i.e. adding together exposures
#'
#' Boxplots for the amalgamated ratios: e.g. the first one is the ratio between 1 and 2+3+...+7. Note the ratio for 3 (second plot).
#' 
#' The boxplots are (1) underrepresented samples (those further to the right, i.e. first principal component larger than the rightmost organoid sample) (subset of (2)), (2) all Nature Genetics samples, (3) organoid samples.
#' 
#' 
#<end_text>

#<begin_chunk>```{r, underrepresented_amalgamation, fig.height=5, echo=FALSE, eval=TRUE}
## amalgamating: everything against s3
## actually I should use the one from the package
# amalgamate <- function(subset, which_to_amalgamate){
#   keep <- subset[,!(1:ncol(subset) %in% which_to_amalgamate)]
#   keep_colnames <- colnames(subset)[which(!(1:ncol(subset) %in% which_to_amalgamate))]
#   notkeep_colnames <- colnames(subset)[which((1:ncol(subset) %in% which_to_amalgamate))]
#   .tmp <- cbind( keep, rowSums(subset[,which_to_amalgamate]))
#   colnames(.tmp) <- c(keep_colnames,
#                       paste0(notkeep_colnames, collapse='+')) 
#   .tmp
# }

ratio_amalgamation <- function(subset, which_to_amalgamate){
  .tmp_amalgamation <- amalgamate(subset, which_to_amalgamate)
  .tmp_amalgamation[,1]/.tmp_amalgamation[,2]
}

boxplot_amalgamation <- function(partition){
  # boxplot(log(ratio_amalgamation(selected_underrepresented_right, partition)),
  #         log(ratio_amalgamation(selected_underrepresented_left, partition)),
  #         log(ratio_amalgamation(natgen_renormalised, partition)),
  #         log(ratio_amalgamation(org, partition)),
  #         names = c("underrepresented right", "underrepresented left", "all natgen", "organoids"),
  #         main=paste0(paste0((1:7)[!((1:7) %in% partition)], collapse="+"), " | all others"))
  
  .title <- paste0('[ ', paste0((1:7)[!((1:7) %in% partition)], collapse="+"), " | all others ]")
  .tmp_ratios <- melt(list(underrepresented_right=log(ratio_amalgamation(selected_underrepresented_right, partition)),
                                                underrepresented_left=log(ratio_amalgamation(selected_underrepresented_left, partition)),
                                                all_Natgen=log(ratio_amalgamation(natgen_renormalised, partition)),
                                                organoids=log(ratio_amalgamation(org, partition))))
  .tmp_ratios$L1 <- factor(.tmp_ratios$L1,  levels=c('underrepresented_right',
                                               'underrepresented_left',
                                               'all_Natgen', 'organoids'))
  ggplot(.tmp_ratios, aes(x=L1, y=value))+geom_boxplot()+theme_light()+
    theme(axis.text.x = element_text(angle = 30, hjust = 1))+ ggtitle(paste0('Ratio for partitioning ', .title))
}

#<end_chunk>

#<begin_chunk>```{r, warning=FALSE, echo=FALSE, eval=TRUE}
par(mfrow=c(1,1))
boxplot_amalgamation(c(1:7)[-c(5, 3)])
#<end_chunk>


#<begin_chunk>```{r, warning=FALSE, echo=FALSE, eval=FALSE}
## building an ilr where the first split is between 3 and everything else
ilr_3 <- ilr(x = rbind(natgen_renormalised, org),
    V = gsi.buildilrBase(W = t(matrix(c(1, 1, -1, 1, 1, 1, 1,
                                      -1, 1 ,0, 1, 1, 1, 1,
                                      0, -1, 0, -1, 1, 1, 1,
                                      0, 0, 0, 0,  -1, 1, 1,
                                      0, 0, 0, 0,  0, -1, 1,
                                      0, 0, 0, 0,  0, 0, -1), ncol = 7, byrow = TRUE))))
## this can't be correct; there are only NaN in last column

ilr_split <- melt(split(ilr_3[,1], f = as.factor(c(natgen_metadata$study, rep('Organoid', dim(org)[1])))))
ggplot(ilr_split, aes(x=value, col=L1))+geom_density()+ggtitle('ILR of s3 vs all others')

#<end_chunk>

#<begin_chunk>```{r, ilr_2, warning=FALSE, echo=FALSE, eval=FALSE}
## here; ongoing
plot(prcomp(as(ilr(x = rbind(natgen_renormalised, org)), 'matrix'))$x[,1:2])

melt(split(prcomp(as(ilr(x = rbind(natgen_renormalised, org)), 'matrix'))$x[,1:2],
           f = as.factor(c(natgen_metadata$study, rep('Organoid', dim(org)[1])))))
ggplot(melt(prcomp(as(ilr(x = rbind(natgen_renormalised, org)), 'matrix'))$x[,1:2]),
       aes(x=value, col=L1))+geom_density()+ggtitle('ILR of s3 vs all others')

#<end_chunk>
#<begin_text>
##' #### Creating a new PCA
##' The results are basically the same
#<end_text>

#<begin_chunk>```{r, princomp_clr,include=TRUE,eval=TRUE,echo=FALSE}

## input matrix is the already clr-transformed matrix
createPCA_fromscratch <- function(input_matrix, annotation, annotation2, return_df=FALSE, labels_active=TRUE){
  prcomp_all_clr <- princomp(input_matrix)

  df_prcomp_exposures_clr <- data.frame(prcomp_all_clr$scores[,1:2],
                                        study=annotation,
                                        bool_any_zeroes=annotation2,
                                        labels=rownames(input_matrix))
    df_prcomp_exposures_clr$labels[!grepl('Sample', df_prcomp_exposures_clr$labels)] <- NA
    df_prcomp_exposures_clr[,'labels'] <- gsub("Organoid Sample", "", df_prcomp_exposures_clr$labels)
  var_explained2 <- (prcomp_all_clr$sdev**2)/sum(prcomp_all_clr$sdev**2)
  
  if(return_df){
    df_prcomp_exposures_clr
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

createPCA_fromscratch(input_matrix = rbind(natgen_clr,org_clr_robustzeroes),
                      annotation = c(natgen_metadata$study, rep('Organoid', nrow(org_clr))),
                      annotation2 = c(rep(FALSE, dim(natgen_metadata)[1]), unlist(apply(org, 1, function(i) any(i == 0)))))+
  ggtitle('PCA created from scratch with robust zeroes')

#<end_chunk>

#<begin_chunk>```{r, raw_pca}
## raw pca, with no compositional aspect
raw_pca_df <- createPCA_fromscratch(input_matrix = rbind(natgen_renormalised,org),
                                    annotation = c(natgen_metadata$study, rep('Organoid', nrow(org_clr))),
                                    annotation2 = c(rep(FALSE, dim(natgen_metadata)[1]), unlist(apply(org, 1, function(i) any(i == 0)))),
                                    return_df = TRUE)
## see what the vertices are
extreme_cases <- rbind(rbind(natgen,org)[which.min(raw_pca_df$Comp.2),],
       rbind(natgen,org)[which.max(raw_pca_df$Comp.2),],
       rbind(natgen,org)[which.min(raw_pca_df$Comp.1),])
rownames(extreme_cases) <- paste0('Extreme sample ', 1:3)
createBarplot(extreme_cases, remove_labels = TRUE)
bpextreme1 <- createBarplot(t(as.matrix(extreme_cases[1,])), remove_labels = TRUE)
bpextreme2 <- createBarplot(t(as.matrix(extreme_cases[2,])), remove_labels = TRUE)
bpextreme3 <- createBarplot(t(as.matrix(extreme_cases[3,])), remove_labels = TRUE)
pca_raw <- createPCA_fromscratch(input_matrix = rbind(natgen_renormalised,org),
                      annotation = c(natgen_metadata$study, rep('Organoid', nrow(org_clr))),
                      annotation2 = c(rep(FALSE, dim(natgen_metadata)[1]), unlist(apply(org, 1, function(i) any(i == 0)))), labels_active = FALSE)+
  ggtitle('PCA created from scratch without any transformation')

blankPanel<-grid.rect(gp=gpar(col="white"))
grid.arrange(blankPanel, bpextreme2, bpextreme3, pca_raw, blankPanel, bpextreme1, ncol=2)
#<end_chunk>

#<begin_text>
#' Pseudocounts analysis (not shown)
#'The lower the pseudocount, the more the organoid samples are left out of the PCA of all other samples.
#' That is because there is quite a good representation of signature landscapes in the NatGen samples, so
#' that anything with more non-zero exposures is likely to be represented
#<end_text>

#<begin_chunk>```{r, pseucocount_PCA,include=FALSE,eval=TRUE,echo=FALSE,fig.height=8}
plt_pc1 <- createPCA_fromscratch(input_matrix = rbind(natgen_clr,clr_mat(addPseudoCounts(org, 1e-7))),
                                 annotation = c(natgen_metadata$study, rep('Organoid', nrow(org_clr))),
                                 annotation2 = c(rep(FALSE, dim(natgen_metadata)[1]), unlist(apply(org, 1, function(i) any(i == 0)))))+
  ggtitle('PCA created from scratch with pseudocount of 1e-7')

plt_pc2 <- createPCA_fromscratch(input_matrix = rbind(natgen_clr,clr_mat(addPseudoCounts(org, 1e-8))),
                                 annotation = c(natgen_metadata$study, rep('Organoid', nrow(org_clr))),
                                 annotation2 = c(rep(FALSE, dim(natgen_metadata)[1]), unlist(apply(org, 1, function(i) any(i == 0)))))+
  ggtitle('PCA created from scratch with pseudocount of 1e-8')

plt_pc3 <- createPCA_fromscratch(input_matrix = rbind(natgen_clr,clr_mat(addPseudoCounts(org, 1e-10))),
                                 annotation = c(natgen_metadata$study, rep('Organoid', nrow(org_clr))),
                                 annotation2 = c(rep(FALSE, dim(natgen_metadata)[1]), unlist(apply(org, 1, function(i) any(i == 0)))))+
  ggtitle('PCA created from scratch with pseudocount of 1e-10')
  
grid.arrange(plt_pc1, plt_pc2, plt_pc3)
#<end_chunk>


#<begin_chunk>```{r, subcomposition_pca, eval=FALSE, echo=FALSE}
apply(org, 2, function(i) sum(i == 0))

subcomp_minus <- list()
rm_subcomp <- list()
subcomp_minus[['4']] <- clr_mat(close_data(rbind(natgen_renormalised,org)[,-4]))
subcomp_minus[['1']] <- clr_mat(close_data(rbind(natgen_renormalised,org)[,-1]))

for(i in names(subcomp_minus)){
  rm_subcomp[[i]] <- apply(apply(subcomp_minus[[i]], 1, is.infinite), 2, any)
  subcomp_minus[[i]] <- subcomp_minus[[i]][!rm_subcomp[[i]],]
}

createPCA_fromscratch(input_matrix = subcomp_minus[['4']],
                      annotation = c(natgen_metadata$study, rep('Organoid', nrow(org_clr)))[!rm_subcomp[['4']]],
                      annotation2 = c(rep(FALSE, dim(natgen_metadata)[1]), unlist(apply(org, 1, function(i) any(i == 0))))[!rm_subcomp[['4']]])+
  ggtitle('PCA created from scratch with a subcomposition (no 4)')

createPCA_fromscratch(input_matrix = subcomp_minus[['1']],
                      annotation = c(natgen_metadata$study, rep('Organoid', nrow(org_clr)))[!rm_subcomp[['1']]],
                      annotation2 = c(rep(FALSE, dim(natgen_metadata)[1]), unlist(apply(org, 1, function(i) any(i == 0))))[!rm_subcomp[['1']]])+
  ggtitle('PCA created from scratch with a subcomposition (no 1)')

#<end_chunk>
#<begin_text>
## Conclusion: the PCA from this subcomposition is very similar to the PCA of all signatures made with robust zeroes.
#<end_text>

#<begin_text>
##' Exposures from Nat Gen don't have zeroes; for organoids many do.
#<end_text>

#<begin_text>
##' ### PCA based on the CN profiles
##' To do (need to get data from Dilrini)
#<end_text>

#<begin_text>
##' ## Dendrograms
##' ### Dendrogram based on the signatures
#<end_text>

#<begin_chunk>```{r, dendrogram_aitchisondistance,echo=FALSE}
organoid_metadata <- cbind.data.frame(study=rep('organoids', nrow(org_clr_robustzeroes)), age=NA, age.cat=NA)
rownames(organoid_metadata) <- rownames(org_clr_robustzeroes)
all_metadata <- rbind(natgen_metadata, organoid_metadata)
all_clr <- rbind(natgen_clr, org_clr_robustzeroes)
rm_infinite <- apply(all_clr, 1, function(x) any(is.infinite(x)))
cat(which(rm_infinite), 'removed due to infinite values')
all_clr_clean <- all_clr[!rm_infinite,]

which(rm_infinite)
#all_clr_clean <- all_clr_clean[c(1:2, grep( 'Organoid', rownames(all_clr_clean))),]

dendro_all <- as.dendrogram(hclust(dist(all_clr_clean)))
levels_study <- levels(factor(all_metadata[labels(dendro_all),'study']))
levels_study
which_level_organoids <- which(grepl('organoids', levels_study))
cols <- rep(NA, length(levels_study))
cols[which_level_organoids] <- '#88E9A2'
cols[-which_level_organoids] <- c('#FFA07A', '#FA8072', '#E9967A', '#F08080')
labels_colors(dendro_all) <- cols[factor(all_metadata[labels(dendro_all),'study'])]
labels_org_bool <- labels_colors(dendro_all) == '#88E9A2'
labels(dendro_all)[labels_org_bool] <- rep('●', sum(labels_org_bool))
labels(dendro_all)[!labels_org_bool] <- rep('•', sum(!labels_org_bool))
labels(dendro_all)[!labels_org_bool] <- rep(NA, sum(!labels_org_bool))
cex_labels <- rep(1, length(labels_org_bool))
cex_labels[labels_org_bool] <- 2
dendro_all <- set(dendro_all, "labels_cex", cex_labels)
plot(dendro_all, cex=0.4, cex.lab=4, main='Dendrogram based on the exposures\n(Aitchison distance)')

table(labels(dendro_all))
#<end_chunk>

#<begin_text>
##' ### Dendrogram based on the CN profiles
#<end_text>


#<begin_text>
##' ## Comments
##' ### Handling of zeroes and robust zeroes
##' The exposures for organoids have some zero values, which is not the case for
##' the ones from Nature Genetics. I could either use the clr in its natural sense so that
##' more than half the exposures are Inf/NaN, or use the robust zeroes approach from the
##' package compositions.
##' Here I have used robust zeroes.
##' 
##' Below are the samples we would use if we didn't (because they contain some zero).
#<end_text>

#<begin_chunk>```{r, echo=FALSE}
all_clr_nonrobust <- rbind(natgen_clr, org_clr)
rm_infinite_nonrobust <- apply(all_clr_nonrobust, 1, function(x) any(is.infinite(x)))
cat(paste0(names(which(rm_infinite_nonrobust)), collapse = ', '), 'removed due to infinite values')

apply(org, 1, function(x) sum(x == 0)/nrow(org))
table(unlist(apply(org, 1, function(x) which(x == 0))))
#<end_chunk>



#<begin_chunk>```{r, echo=FALSE, eval=FALSE}

# image(cov(sweep(natgen_renormalised[,-7], 1, natgen_renormalised[,7], '/')))
#<end_chunk>

#<begin_text>
##' More on zeroes: talk with Dom
##' 
##' Missingicompositions: check this section of the compositions package
##' ILR is with imputation
##' 
##' Patterns of missing signatures: 100, 101, 110, … etc
##' 
##' - Ilr on pseudo count means that you do a clustering on missing data
##' - Logistic SVD with the dichotomised matrix is the way to go!
##' - Signal is picked up when you dichotomise it, not with normal irl, which totally misses it
##' 
##' There is a laplace prior on zeroes for the Lasso
##' 
##' To do:
##' - Get ratings of signatures from Geoff
##' - See what sort of dats patterns there are
##' It’s problematic not to have any zero when you are doing the dichotomised approach
#<end_text>

#<begin_text>
dichot <- apply(org, 2, function(i) i == 0)
dichot <- apply(dichot, 2, as.numeric)
rownames(dichot) <- rownames(org); colnames(dichot) <- colnames(org)
pheatmap(dichot, cluster_rows = FALSE, cluster_cols = FALSE)
#<end_text>

#<begin_text>
##' ## Alternative transformations
##' Alternative transformations might deal better with zeroes.
##' ### Tsagris, Preston, and Wood (2011)
##' Alternative transformation of the data following section 3.2 from [this paper](https://www.tandfonline.com/doi/full/10.1080/01621459.2014.990563#aHR0cHM6Ly93d3cudGFuZGZvbmxpbmUuY29tL2RvaS9wZGYvMTAuMTA4MC8wMTYyMTQ1OS4yMDE0Ljk5MDU2Mz9uZWVkQWNjZXNzPXRydWVAQEAw)
#<end_text>


#<begin_chunk>```{r, alt_transformation,echo=FALSE, eval=TRUE, include=FALSE}
tsagris_transformation <- function(u, alpha){
  ## u is a compositional vector
  ## alpha is a scalar
  p <- length(u)
  1/alpha*( (u**alpha) / (1/p)*(u**alpha) - rep(1,p))
}

tsagris_transformation_matrix <- function(U, alpha){
  t(apply(U, 1, tsagris_transformation, alpha))
}

# Umat <- MCMCpack::rdirichlet(10, rep(1,4))
# tsagris_transformation_matrix(Umat, alpha=1.03)
# tsagris_transformation(Umat[1,], alpha=1.03)

createPCA_fromscratch(input_matrix = tsagris_transformation_matrix(rbind(natgen,org), alpha = .8),
                      annotation = c(natgen_metadata$study, rep('Organoid', nrow(org_clr))),
                      annotation2 = c(rep(FALSE, dim(natgen_metadata)[1]), unlist(apply(org, 1, function(i) any(i == 0)))))+
  ggtitle('PCA created from scratch with tsagris transformation, alpha=.8')

createPCA_fromscratch(input_matrix = tsagris_transformation_matrix(rbind(natgen,org), alpha = 1.03),
                                 annotation = c(natgen_metadata$study, rep('Organoid', nrow(org_clr))),
                                 annotation2 = c(rep(FALSE, dim(natgen_metadata)[1]), unlist(apply(org, 1, function(i) any(i == 0)))))+
  ggtitle('PCA created from scratch with tsagris transformation, alpha=1.03')

createPCA_fromscratch(input_matrix = tsagris_transformation_matrix(rbind(natgen,org), alpha = 1.2),
                      annotation = c(natgen_metadata$study, rep('Organoid', nrow(org_clr))),
                      annotation2 = c(rep(FALSE, dim(natgen_metadata)[1]), unlist(apply(org, 1, function(i) any(i == 0)))))+
  ggtitle('PCA created from scratch with tsagris transformation, alpha=1.2')

# image(cov(sweep(natgen_renormalised[,-7], 1, natgen_renormalised[,7], '/')))
#<end_chunk>

#<begin_text>
#' ## CLR pairs
#<end_text>

#<begin_chunk>```{r, clr_pairs,echo=FALSE, eval=TRUE, include=TRUE}
clr_rbind <- list()
clr_rbind[[1]] <- rbind(natgen_clr,clr(org))
clr_rbind[[2]] <- rbind(natgen_clr,clr_mat(addPseudoCounts(org, 1e-4)))
clr_rbind[[3]] <- rbind(natgen_clr,clr_mat(addPseudoCounts(org, 1e-7)))

for(i in 1:3){
  pairs(clr_rbind[[i]], col=c("#4a1d96", "#22ae6b")[factor(c(rep(1, nrow(natgen_clr)),
                              c(rep(2, nrow(org)))))],
        pch=c(15, 2)[as.factor(c(rep(FALSE, dim(natgen_metadata)[1]), unlist(apply(org, 1, function(i) any(i == 0)))))],
        main='clr pairs')
}
#<end_chunk>


