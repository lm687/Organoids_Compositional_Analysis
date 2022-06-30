### Files related to WGD clades
###            ../check_WGD.R
###            TCGA_HGSOC_dendrogram.R
###            two_cluster_exposures_modelling.R

rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(TMB)
library(ggplot2)
library(cowplot)
library(tikzDevice)
library(ggrepel)
library(reshape2)
library(gridExtra)
library(dplyr)
library(readxl)
library(dendextend)
library(ggdendro)

source("../helper_functions.R")
source("../../../../GlobalDA/code/2_inference_TMB/helper_TMB.R")
source("../../../mcneish/code/models/helper_TMB.R")
source("../../../britroc-1/code/models/helper/functions.R")

cluster_fig1 <- readRDS("../../copy_number_analysis_organoids/robjects/dendrograminputclr_tree.RDS")
exposures_fig1 <- readRDS("../../copy_number_analysis_organoids/robjects/allnatgen_UpdatedExposures.RDS")

dim(exposures_fig1)

.cutree <- cutree(cluster_fig1, k = 2)
table(.cutree)

rowSums(exposures_fig1)

.cutree <- .cutree[match(rownames(exposures_fig1), names(.cutree))]
length(.cutree)
dim(exposures_fig1)

cluster_fig1

#------------------------------------------------------------------#
#---------------------------- plotting ----------------------------#
#------------------------------------------------------------------#

dendrograminputclr <- readRDS("../../copy_number_analysis_organoids/robjects/dendrograminputclr.RDS")
heatmapinputclr <- readRDS("../../copy_number_analysis_organoids/robjects/heatmapinputclr_with_ticks.RDS")

e1 <- dendrograminputclr#+geom_label_repel(label.size = NA)
e2 <- heatmapinputclr+  theme(axis.title.x=element_text(), axis.title.y=element_text(angle=90))+
  labs(x='Public tumour datasets (TCGA, PCAWG and BriTROC) and organoids', y='Copy number signature activity')+
  theme(axis.text.y = element_text(size = 8), axis.title.y = element_text(vjust=3))
# e <- plot_grid(e1, e2, nrow = 2, labels = c('e'), rel_widths = c(3,5))
e2_v2 <- e2#+geom_label_repel(fill='white', col='black',nudge_x = 0, nudge_y = -5,
           #                  aes(y=0,label=ifelse(grepl('PDO', Var2) & Var1 == 's1',
           #                                       as.character(Var2), NA)))+
  # scale_y_continuous(expand=c(0,0))+lims(y=c(-1.6, 1))
e1_v2 <- e1; e1_v2$layers[[2]] = NULL; e1_v2$scales$scales[[2]] <- NULL
e1_v2$theme$plot.margin[3] = unit(0.0, "cm")
e1_v2$theme$plot.margin[2] = unit(-25, "points")
e1_v2$theme$plot.margin
# e1_v2 <- e1_v2+scale_y_continuous(expand=c(0,0))
e1_v2

e2_v2$theme$plot.margin[1] = unit(0, "cm")
# e_v2 <- plot_grid(plot_grid(plot.new(), e1_v2, rel_widths=c(0.07, 5), rel_heights = c(3,5)),
e_v2 <- plot_grid(plot_grid(e1_v2),
                  e2_v2, nrow = 2)
e_v2

##-------------------------------------------------------------------------------------------##
## Adding annotation from cohorts
annotation_heatmap <- (give_annotation_from_names(as.character(e2_v2$data$Var2)))
annotation_heatmap <- (reshape2::melt(cbind.data.frame(idx=1:length(annotation_heatmap), anno=annotation_heatmap)))
annotation_heatmap <- cbind.data.frame(annotation_heatmap, t(apply(sapply(annotation_heatmap$anno, function(i) (i == unique(annotation_heatmap$anno))), 2, as.numeric)))
colnames(annotation_heatmap)[-c(1:3)] <- unique(annotation_heatmap$anno)
colnames(annotation_heatmap)[3] <- 'idx'
head(annotation_heatmap)
annotation_heatmap <- reshape2::melt(annotation_heatmap[,-c(1:2)], id.vars='idx')
annotation_heatmap$variable <- as.character(annotation_heatmap$variable)
annotation_heatmap <- annotation_heatmap[annotation_heatmap$value == 1,]
annotation_heatmap$variable[annotation_heatmap$variable == 'BriTROC-1'] <- 'BriTROC'
annotation_cohort <- ggplot(annotation_heatmap,
       aes(x=idx, y=variable), col='black')+geom_tile()+theme_bw()+guides(fill='none')+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        # axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text=element_text(size=8))+
  coord_cartesian(expand = FALSE)

pdf( '../../../../CDA_in_Cancer/text/thesis/figures/general_compositional/partial_fig1_2_with_annopdf.pdf', height = 5, width = 6)
plot_grid(plot_grid(e1_v2),
          e2_v2, annotation_cohort, nrow = 3, rel_heights = c(3,3,1))
dev.off()
tikz( '../../../../CDA_in_Cancer/text/thesis/figures/general_compositional/partial_fig1_2_with_anno.tex', height = 5, width = 6)
plot_grid(plot_grid(e1_v2),
          e2_v2, annotation_cohort, nrow = 3, rel_heights = c(3,3,1))
dev.off()
##-------------------------------------------------------------------------------------------##

# system("open ../../../../CDA_in_Cancer/text/thesis/figures/")
tikz( '../../../../CDA_in_Cancer/text/thesis/figures/general_compositional/partial_fig1_2.tex', height = 4, width = 6)
e_v2
dev.off()

pdf( '../../../../CDA_in_Cancer/text/thesis/figures/general_compositional/partial_fig1_2_image.pdf',
     height = 4, width = 6)
e_v2
dev.off()


#------------------------------------------------------------------#
#-------------------------- model (zeros) -------------------------#
#------------------------------------------------------------------#

# bernoulli for zeros
TMB::compile("../../../britroc-1/code/models/tmb_RE/tmb_correlated_multinom_2_allFE_b.cpp", "-std=gnu++17")
dyn.load(dynlib("../../../britroc-1/code/models/tmb_RE/tmb_correlated_multinom_2_allFE_b"))
# partial ILR for nonzero exposures
TMB::compile("../../../britroc-1/code/models/tmb_RE/tmb_MVN_partial_ILR_FEb.cpp", "-std=gnu++17")
dyn.load(dynlib("../../../britroc-1/code/models/tmb_RE/tmb_MVN_partial_ILR_FEb"))
# partial ILR for nonzero exposures, with correlations
TMB::compile("../../../britroc-1/code/models/tmb_RE/tmb_MVN_partial_ILR_FEe.cpp", "-std=gnu++17")
dyn.load(dynlib("../../../britroc-1/code/models/tmb_RE/tmb_MVN_partial_ILR_FEe"))


exposures_zeros = apply(exposures_fig1 == 0, 2, as.numeric)

X <- cbind(1, as.numeric(.cutree == 2))

## Zeros
TMB_data_zeros = list(Y = exposures_zeros,
                      num_sigs = ncol(exposures_zeros),
                      x = X,
                      z=diag(nrow(X)))

fraction_zeros = sapply(list(clade_Hs3=exposures_zeros[TMB_data_zeros$x[,2] == 0,],
                             clade_Hs4=exposures_zeros[TMB_data_zeros$x[,2] == 1,]), colMeans)
fraction_zeros

bar_zeros <- ggplot(melt(fraction_zeros), aes(x=Var1, y=value, fill=Var2))+
  geom_bar(stat = "identity", position=position_dodge(width=0.5), alpha=0.7)+
  theme_bw()+labs(x='Signatures', y='Fraction of zero exposures', col='Clade')

pdf("../figures/LN_clusters_exposures_WGD/clades_barplot_zeros.pdf",
    width = 3, height = 1.8, onefile = F)
bar_zeros 
dev.off()

tikzDevice::tikz("../figures/LN_clusters_exposures_WGD/clades_barplot_zeros.tex",
    width = 3, height = 1.8, onefile = F)
bar_zeros 
dev.off()

# results_zeros = give_res_report(TMB_data_zeros, d=ncol(exposures_zeros))
results_zeros = wrapper_run_TMB_use_nlminb(model = 'bernoulliFE', object = TMB_data_zeros)
results_zeros

# wald_generalised(v = select_slope_2(results_zeros$report$par.fixed,v=F), sigma = results_zeros$report$cov.fixed[c(F,T), c(F,T)])
wald_generalised(v = select_slope_2(results_zeros$par.fixed,v=F), sigma = results_zeros$cov.fixed[c(F,T), c(F,T)])

pdf("../figures/LN_clusters_exposures_WGD/clades_estimate_beta_slope_zeros.pdf",
    width = 3, height = 1.8, onefile = F)
plot(plot_betas(results_zeros, names_cats = paste0('s', 1:7), xlab = '', return_plot=T))
dev.off()

colMeans(TMB_data_zeros$Y)

xtable::xtable(do.call('rbind', lapply(split_rows(df = TMB_data_zeros$Y,
                                                  f = TMB_data_zeros$x[,2]), colMeans)))

fitted_logR <- TMB_data_zeros$x %*% matrix(results_zeros$par.fixed, nrow=2)
fitted_vals = apply(fitted_logR, 2, function(i) exp(i)/(1+exp(i)) ) ## transform to probabilities
colnames(fitted_vals) = colnames(exposures_zeros)
rownames(fitted_vals) = rownames(exposures_zeros)

ggplot(melt(fitted_vals[order(TMB_data_zeros$x[,2]),]))+
  geom_tile( aes( x = Var2, y = Var1, fill = value )  )+
  scale_fill_gradient(low="blue", high="red")+theme_bw()+
  geom_point(data = melt(exposures_zeros[order(TMB_data_zeros$x[,2]),])  %>% filter(value == 1), aes(x = Var2, y = Var1),
             col='white',  shape=18, size=.3)+
  # ggtitle('Observed zeros and\nfitted probabilities')+
  labs(x='Signatures' , y='Samples', fill='Probability')+
  theme(legend.position = "bottom")+
  scale_x_discrete(expand = c(0, 0))+  scale_y_discrete(expand = c(0, 0))
ggsave("../figures/LN_clusters_exposures_WGD/clades_fitted_zeros.pdf", width = 3, height = 3)

fitted_probs = t(fitted_vals[c(which(TMB_data_zeros$x[,2] == 0)[1],
                               which(TMB_data_zeros$x[,2] == 1)[1]),]); colnames(fitted_probs) = c('clade_Hs3', 'clade_Hs4')
## fitted probability of zero vs actual fraction of zeros
comparison_probs = (dcast(melt(list(fraction_zeros=fraction_zeros,
                                    fitted_probs=fitted_probs), measure.vars=c('L1')), Var1+Var2~L1))

cbind(fraction_zeros=fraction_zeros, fitted_probs=fitted_probs)

ggplot(comparison_probs, aes(x=fraction_zeros, y=fitted_probs, col=Var1, shape=Var2))+geom_point()+
  geom_abline(coef=c(0,1), lty='dashed')+
  #ggtitle('Observed fraction of zero and fitted probability')+
  theme_bw()+labs(y='Fitted probabilities',
                  x='Fraction of zeros')+theme(legend.title =element_blank())
ggsave("../figures/LN_clusters_exposures_WGD/clades_fitted_zeros_observed_fractions.pdf", height=3, width=3.8)

ggplot(comparison_probs, aes(x=fraction_zeros, y=fitted_probs, col=Var1, shape=Var2))+
  geom_point(size=3)+
  geom_abline(coef=c(0,1), lty='dashed')+#ggtitle('Observed fraction of zero and fitted probability')+
  theme_bw()+facet_wrap(.~Var2)+labs(y='Fitted probabilities',
                                     x='Fraction of zeros')+theme(legend.title =element_blank())
ggsave("../figures/LN_clusters_exposures_WGD/clades_fitted_zeros_observed_fractions_v2.pdf", height=3, width=6)

## comparison per patient
comparison_exposures_clades <- cbind(fitted=melt(fitted_vals[order(TMB_data_zeros$x[,2]),]),
      obs=melt(exposures_zeros[order(TMB_data_zeros$x[,2]),]),
      cohort=rep(give_annotation_from_names(rownames(exposures_fig1[order(TMB_data_zeros$x[,2]),])), 7))
ggplot(comparison_exposures_clades, aes(x=fitted.value, y=obs.value))+
  geom_point(alpha=0.2)+theme_bw()+facet_wrap(.~cohort)


split_exposures_per_cohort <- split_rows(exposures_fig1, give_annotation_from_names(rownames(exposures_fig1)))
names(split_exposures_per_cohort) <- unique(give_annotation_from_names(rownames(exposures_fig1)))
lapply(split_exposures_per_cohort, min)
## only in TCGA and organoids are there zeros

#------------------------------------------------------------------#
#--------------------------- model (ILR) --------------------------#
#------------------------------------------------------------------#

partialILR = give_partial_irl(exposures_fig1)
as(compositions::ilr(exposures_fig1), 'matrix')
partialILR

irl_robust = as(compositions::ilr(exposures_fig1), 'matrix')

plot(unlist(irl_robust), unlist(partialILR)) ## remember that zeros are treated differently in the inference!

.keep = (rowSums(partialILR == 0) < (ncol(exposures_fig1) - 2) )
partialILR = partialILR[.keep,]  ## if a sample has too many zeros (i.e. all but 1 or 2) remove
colnames(partialILR) = paste0('ILR', 1:ncol(partialILR))
dim(partialILR)

TMB_data_partialILR = list(Y = partialILR,
                           d = ncol(partialILR),
                           n = sum(.keep),
                           x = cbind(1, as.numeric(.cutree[match(rownames(partialILR), names(.cutree))] == 1)),
                           z = diag(sum(.keep)))

# results_ILR = give_res_report(TMB_data_partialILR, d = (ncol(partialILR)), model='tmb_MVN_partial_ILR_FEb')
# results_ILR = wrapper_run_TMB_use_nlminb(TMB_data_partialILR, model='tmb_MVN_partial_ILR_FEb') ## no cor
results_ILR = wrapper_run_TMB_use_nlminb(TMB_data_partialILR, model='tmb_MVN_partial_ILR_FEe') ## cor

results_ILR
results_ILR$pdHess

#----
pdf("../figures/LN_clusters_exposures_WGD/clades_estimate_beta_slope_ILR.pdf",
    width = 4, height = 2.5, onefile = F)
plot(plot_betas(results_ILR, names_cats = paste0('ILR', 1:6), xlab = '', return_plot = T))
dev.off()

## s3 decreases, s1 possibly decreases
##  s4 increases, s2 and s6 increase
## s5, s7: unclear

exp(python_like_select_rownames(summary(results_ILR), 'logsd')[,1])

colMeans(exposures_fig1 == 0)
colMeans(partialILR == 0)

ggplot(melt(list(Hs3=exposures_fig1[rownames(TMB_data_partialILR$Y),][TMB_data_partialILR$x[,2] == 0,],
                 Hs4=exposures_fig1[rownames(TMB_data_partialILR$Y),][TMB_data_partialILR$x[,2] == 1,])),
       aes(x=value, col=L1))+geom_density()+facet_wrap(.~Var2, nrow=2, scales = "free_x")+
  scale_color_brewer(palette = "Set2")+theme_bw()+theme(legend.position = "bottom")
ggsave("../figures/LN_clusters_exposures_WGD/clades_densities_change.pdf", width = 5, height = 3.4)

wald_generalised(v = select_slope_2(results_ILR$par.fixed,v=F), sigma = results_ILR$cov.fixed[c(F,T), c(F,T)])
wald_TMB_wrapper(results_ILR)

## Put the results of the two models together
ggplot(cbind.data.frame(label=colnames(exposures_zeros),
                        col=(colnames(exposures_zeros) == 's1'), ## s1 is used as baseline, somewhat
                        partialILR=rbind(0, python_like_select_rownames(summary(results_ILR), 'beta')[c(F,T),])[match(colnames(exposures_zeros), colnames(exposures_fig1)),],
                        zeros=python_like_select_rownames(summary(results_zeros), 'beta')[c(F,T),]),
       aes(x=partialILR.Estimate, y=-zeros.Estimate, label=label))+geom_point()+
  geom_hline(yintercept = 0, col='blue', lty='dashed', size=0.3)+
  geom_vline(xintercept = 0, col='blue', lty='dashed', size=0.3)+
  geom_errorbar(aes(xmin=`partialILR.Estimate`-`partialILR.Std. Error`,xmax=`partialILR.Estimate`+`partialILR.Std. Error`), width=.1)+
  geom_errorbar(aes(ymin=-`zeros.Estimate`-`zeros.Std. Error`, ymax=-`zeros.Estimate`+`zeros.Std. Error`), width=.1)+
  geom_label_repel(aes(col=col))+facet_wrap(.~label, nrow=1)+guides(col='none')+
  labs(x='Partial ILR estimate', y='Bernoulli estimate')+
  # ggtitle('Summary of both models')+
  theme_bw()
ggsave("../figures/LN_clusters_exposures_WGD/clades_summary_both_sets_coefs.pdf", width = 8, height = 1.8)

##---------------------------------------------------------------------------------------------------------##
## Imputation

TMB::compile("../../../../other_repos/britroc-1/code/models/tmb_RE/tmb_MVN_with_mean_2.cpp",  "-std=gnu++17")
dyn.load(dynlib("../../../../other_repos/britroc-1/code/models/tmb_RE/tmb_MVN_with_mean_2"))
mod_model_name = "mvnbetacor"

TMB::compile("../../../../other_repos/britroc-1/code/models/tmb_RE/tmb_MVN_with_mean_2_scaledsdpergroup.cpp",  "-std=gnu++17")
dyn.load(dynlib("../../../../other_repos/britroc-1/code/models/tmb_RE/tmb_MVN_with_mean_2_scaledsdpergroup"))

TMB_data_imput1 = list(Y = as(compositions::ilr(impute(exposures_fig1, inputation_value = 1e-2)), 'matrix'),
                       d = 6,
                       n = nrow(exposures_fig1),
                       x = cbind(1, 1-(.cutree-1)), ## so that 1=WGD, 0=non-WGD
                       z= give_z_matrix_from_labels(1:nrow(exposures_fig1)))

imput1 <- wrapper_run_TMB_use_nlminb(model = 'mvnbetacor', object = TMB_data_imput1,
                                             use_nlminb = T)
plot_betas(imput1, names_cats = paste0('ILR', 1:6))

imput1overdisp <- wrapper_run_TMB_use_nlminb(model = 'tmb_MVN_with_mean_2_scaledsdpergroup', object = TMB_data_imput1,
                                     use_nlminb = T)
pdf("../figures/LN_clusters_exposures_WGD/clades_estimate_beta_slope_imputationILR.pdf",
    width = 4, height = 2.5, onefile = F)
plot(plot_betas(imput1overdisp, names_cats = paste0('ILR', 1:6), xlab = '', return_plot = T))
dev.off()


## extremely similar results in all cases
pdf("../figures/LN_clusters_exposures_WGD/clades_estimate_beta_slope_comparison.pdf",
    width = 4, height = 5.5, onefile = F)
grid.arrange(plot_betas(results_ILR, title='partial ILR', names_cats = paste0('ILR', 1:6)),
             plot_betas(imput1, title = 'imput 1e-2', names_cats = paste0('ILR', 1:6)),
             plot_betas(imput1overdisp, title = 'imput 1e-2, overdisp', names_cats = paste0('ILR', 1:6)), nrow=3)
dev.off()

##---------------------------------------------------------------------------------------------------------##
## Pending: how to go back to fitted values from partial ILR results? I guess we have to put NA in many places
## (but I don't think there is much a thing as a fitted value... unless I put a zero everywhere else)

# inverse_partialILR <- function(i, d){
#   apply(TMB_data_partialILR$Y, 1, function(j){
#     .which_zero = as.vector(which(j==0))
#     .basis_row = give_partial_ilr_basis(.which_zero, d)
#     
#   })
# }

fitted_logR_partialILR <- TMB_data_partialILR$x %*% matrix(results_ILR$par.fixed[grepl('beta', names(results_ILR$par.fixed))], nrow=2)
fitted_vals_partialILR = apply(cbind(fitted_logR_partialILR, 0), 2, function(i) exp(i)/(1+exp(i)) ) ## transform to probabilities
colnames(fitted_vals_partialILR) = colnames(exposures_fig1) #paste0('ILR', 1:ncol(partialILR))
colnames(exposures_fig1)
rownames(fitted_vals_partialILR) = rownames(TMB_data_partialILR$Y)
fitted_vals_partialILR

.sftmax_partialILR <- cbind(partialILR, 0)
colnames(.sftmax_partialILR) <- colnames(exposures_fig1)
partialILR_melt <- (melt(apply(.sftmax_partialILR, 2, function(i) exp(i)/(1+exp(i)) ) ))
fitted_vals_partialILR_melt <- (melt(apply(fitted_vals_partialILR, 2, function(i) exp(i)/(1+exp(i)) ) ))

ggplot(cbind.data.frame(partialILR_melt,
                        values_observed=fitted_vals_partialILR_melt[match(paste0(partialILR_melt$Var1, partialILR_melt$Var2),
      paste0(fitted_vals_partialILR_melt$Var1, fitted_vals_partialILR_melt$Var2)),3]), aes(x=value, y=values_observed))+
  geom_point(alpha=0.2)+facet_wrap(.~Var2)+theme_bw()
ggsave("../figures/LN_clusters_exposures_WGD/clades_fitted_partialILR.pdf", width = 3, height = 4)

partialILR_melt[partialILR_melt$Var2 == 's7',]


#------------------------------------------------------------------#
#------------- model on independent WGD classification ------------#
#------------------------------------------------------------------#

## Using the two external groups (WGD and not) as covariates
wgd_hasse <- read.table("../../../../CDA_in_Cancer/data/TCGA_WGD/Haase2019_TCGA.giScores.wgd.txt", header = T)
wgd_icgc <- read_xlsx("../../../../CDA_in_Cancer/data/pcawg_supplementary_tables/41586_2019_1907_MOESM4_ESM.xlsx")
wgd_icgc_withWGD <- sapply(unique(wgd_icgc$samplename), function(i) any(wgd_icgc[wgd_icgc$samplename == i,'type'] == 'WGD'))
wgd_icgc_withWGD <- cbind.data.frame(sample=names(wgd_icgc_withWGD), WGD=wgd_icgc_withWGD)

unique(wgd_icgc_withWGD$WGD)
unique(wgd_hasse$wgd)

## using both ICGC and TCGA data. 633 samples in total
wgd <- rbind.data.frame(cbind(wgd_hasse$patient, wgd_hasse$wgd),
      cbind(wgd_icgc_withWGD$sample, wgd_icgc_withWGD$WGD))
wgd <- wgd[match(rownames(partialILR), wgd$V1),]
wgd <- wgd[!is.na(wgd$V1),]

length(wgd$V1)
length(unique(wgd$V1))

exposures_WGD = exposures_fig1[remove_na(match(wgd$V1, rownames(exposures_fig1))),]
exposures_zeros_WGD = apply(exposures_WGD == 0, 2, as.numeric)
rownames(exposures_zeros_WGD) <- rownames(exposures_WGD)
partialILR_WGD <- partialILR[remove_na(match(wgd$V1, rownames(partialILR))),]
dim(exposures_zeros_WGD)
dim(partialILR_WGD)
exposures_zeros_WGD[apply(exposures_zeros_WGD, 1, function(i) any(is.na(i))),]

stopifnot(dim(partialILR_WGD)[1] == dim(exposures_zeros_WGD)[1])
stopifnot(all(rownames(partialILR_WGD) == wgd$V1))
# exposures_fig1_subset <- exposures_fig1[match(wgd$V1, rownames(exposures_fig1)),]

TMB_data_partialILR_WGD = list(Y = partialILR_WGD,
                           d = ncol(partialILR_WGD),
                           n = nrow(partialILR_WGD),
                           x = cbind(1, as.numeric(factor(wgd$V2))-1),
                           z=diag(nrow(partialILR_WGD)))
# TMB_data_partialILR_WGD_probs = list(Y = exposures_fig1_subset,
#                                d = ncol(partialILR),
#                                n = nrow(partialILR),
#                                x = cbind(1, as.numeric(factor(wgd$V2))-1),
#                                z=diag(nrow(exposures_fig1_subset))) ## dummy, for the number of individuals

# results_ILR = give_res_report(TMB_data_partialILR_WGD, d = (ncol(partialILR)), model='tmb_MVN_partial_ILR_FEb')
results_ILR_WGDnocor = wrapper_run_TMB_use_nlminb(TMB_data_partialILR_WGD, model='tmb_MVN_partial_ILR_FEb')
results_ILR_WGD = wrapper_run_TMB_use_nlminb(TMB_data_partialILR_WGD, model='tmb_MVN_partial_ILR_FEe')

results_ILR_WGD
results_ILR_WGD$pdHess

plot_betas(results_ILR_WGD)
pdf("../figures/LN_clusters_exposures_WGD/WGD_betas_partialILR.pdf", width = 3, height = 2, onefile = F)
plot(plot_betas(results_ILR_WGD, names_cats = paste0('ILR', 1:6), xlab = '', return_plot = T))
dev.off()

#-----------------
## Zeros
TMB_data_zeros_WGD = list(Y = exposures_zeros_WGD,
                          num_sigs = ncol(exposures_zeros_WGD),
                          x = cbind(1, as.numeric(factor(wgd$V2))-1),
                          z=diag(nrow(ncol(exposures_zeros_WGD))))
fraction_zeros_WGD = sapply(list(nonWGD=exposures_zeros_WGD[TMB_data_zeros_WGD$x[,2] == 0,],
                                 WGD=exposures_zeros_WGD[TMB_data_zeros_WGD$x[,2] == 1,]), colMeans)
fraction_zeros_WGD

exposures_zeros_WGD[TMB_data_zeros_WGD$x[,2] == 0,]

bar_zeros_WGD <- ggplot(melt(fraction_zeros_WGD), aes(x=Var1, y=value, fill=Var2))+
  geom_bar(stat = "identity", position=position_dodge(width=0.5), alpha=0.7)+
  theme_bw()+labs(x='Signatures', y='Fraction of zero exposures', col='Clade')
bar_zeros_WGD

colSums(fraction_zeros_WGD)
rowSums(fraction_zeros_WGD)

pdf("../figures/LN_clusters_exposures_WGD/WGD_barplot_zeros.pdf",
    width = 3, height = 1.8, onefile = F)
bar_zeros_WGD 
dev.off()

results_zeros_WGD = wrapper_run_TMB_use_nlminb(model = 'bernoulliFE', object = TMB_data_zeros_WGD)
results_zeros_WGD

# TMB::compile("../../../britroc-1/code/models/tmb_RE/mvn_beta_no_cor.cpp", "-std=gnu++17")
# dyn.load(dynlib("../../../britroc-1/code/models/tmb_RE/mvn_beta_no_cor"))
# 
# res_ALR_impute <- wrapper_run_TMB_use_nlminb(model = 'mvn_beta_no_cor', object = give_ALR(imputate_TMB(TMB_data_partialILR_WGD_probs, 1e-2)),
#                            use_nlminb = T)
# plot_betas(res_ALR_impute, names=(paste0('ALR', 1:6)))
# 
# wald_TMB_wrapper(results_ILR)
# wald_TMB_wrapper(res_ALR_impute)
# 
# ggplot(data.frame(beta_slope=python_like_select_name(res_ALR_impute$par.fixed, 'beta')[c(F,T)],
#                   total_exposure=colSums(exposures_fig1[,-7]),
#                   sig=paste0('s', 1:6)), aes(x=beta_slope, y=total_exposure, label=sig))+
#   geom_point()+geom_label()+theme_bw()+scale_y_continuous(trans = "log2")

results_zeros_ILR_WGD <- cbind.data.frame(label=colnames(exposures_zeros_WGD),
                                          col=(colnames(exposures_zeros) == 's1'), ## s1 is used as baseline, somewhat
                                          partialILR=rbind(0, python_like_select_rownames(summary(results_ILR_WGD), 'beta')[c(F,T),])[match(colnames(exposures_zeros_WGD), colnames(exposures_fig1)),],
                                      zeros=python_like_select_rownames(summary(results_zeros_WGD), 'beta')[c(F,T),])
ggplot(results_zeros_ILR_WGD,
       aes(x=partialILR.Estimate, y=-zeros.Estimate, label=label, col=col))+geom_point()+
  geom_hline(yintercept = 0, col='blue', lty='dashed', size=0.3)+
  geom_vline(xintercept = 0, col='blue', lty='dashed', size=0.3)+
  geom_errorbar(aes(xmin=`partialILR.Estimate`-`partialILR.Std. Error`,xmax=`partialILR.Estimate`+`partialILR.Std. Error`),
                width=.1)+
  geom_errorbar(aes(ymin=-`zeros.Estimate`-`zeros.Std. Error`, ymax=-`zeros.Estimate`+`zeros.Std. Error`), width=.1)+
  geom_label_repel()+facet_wrap(.~label, nrow=1, scales = "free")+
  # ggtitle('Summary of both models')+
  theme_bw()+
  labs(x='Partial ILR estimate', y='Bernoulli estimate')+guides(col=F)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave("../figures/LN_clusters_exposures_WGD/WGD_summary_both_sets_coefs.pdf", width = 10, height = 2)

ggplot(results_zeros_ILR_WGD,
       aes(x=partialILR.Estimate, y=-zeros.Estimate, label=label, col=col))+geom_point()+
  geom_hline(yintercept = 0, col='blue', lty='dashed', size=0.3)+
  geom_vline(xintercept = 0, col='blue', lty='dashed', size=0.3)+
  geom_errorbar(aes(xmin=`partialILR.Estimate`-`partialILR.Std. Error`,xmax=`partialILR.Estimate`+`partialILR.Std. Error`),
                width=.1)+
  geom_errorbar(aes(ymin=-`zeros.Estimate`-`zeros.Std. Error`, ymax=-`zeros.Estimate`+`zeros.Std. Error`), width=.1)+
  geom_label_repel()+facet_wrap(.~label, nrow=1)+
  theme_bw()+
  labs(x='Partial ILR estimate', y='Bernoulli estimate')+guides(col=F)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave("../figures/LN_clusters_exposures_WGD/WGD_summary_both_sets_coefs_sharedaxes.pdf", width = 10, height = 2)

ggplot(melt(list(nonWGD=exposures_WGD[TMB_data_partialILR_WGD$x[,2] == 0,],
                 WGD=exposures_WGD[TMB_data_partialILR_WGD$x[,2] == 1,])),
       aes(x=value, col=L1))+geom_density()+facet_wrap(.~Var2, nrow=1, scales = "free_x")+
  scale_color_brewer(palette = "Set2")+theme_bw()+theme(legend.position = "bottom")+
  labs(col='', x='Exposure')
ggsave("../figures/LN_clusters_exposures_WGD/WGD_densities_change.pdf", width = 10, height = 2.3)

## compare fit to observed
fitted_logR_WGD <- TMB_data_zeros_WGD$x %*% matrix(results_zeros_WGD$par.fixed, nrow=2)
fitted_vals_WGD = apply(fitted_logR_WGD, 2, function(i) exp(i)/(1+exp(i)) ) ## transform to probabilities
colnames(fitted_vals_WGD) = colnames(exposures_zeros_WGD)
rownames(fitted_vals_WGD) = rownames(exposures_WGD)

min(fitted_vals_WGD)
max(fitted_vals_WGD)
fitted_probs_WGD = t(fitted_vals[c(which(TMB_data_zeros_WGD$x[,2] == 0)[1],
                                   which(TMB_data_zeros_WGD$x[,2] == 1)[1]),]); colnames(fitted_probs_WGD) = c('nonWGD', 'WGD')
round(fitted_probs_WGD, digits = 3)
## fitted probability of zero vs actual fraction of zeros
comparison_probs_WGD = (dcast(melt(list(fraction_zeros=fraction_zeros_WGD,
                                    fitted_probs=fitted_probs_WGD), measure.vars=c('L1')), Var1+Var2~L1))
ggplot(comparison_probs_WGD, aes(x=fraction_zeros, y=fitted_probs, col=Var1, shape=Var2))+geom_point()+
  geom_abline(coef=c(0,1), lty='dashed')+
  #ggtitle('Observed fraction of zero and fitted probability')+
  theme_bw()+labs(y='Fitted probabilities',
                  x='Fraction of zeros')+theme(legend.title =element_blank())
ggsave("../figures/LN_clusters_exposures_WGD/WGD_fitted_zeros_observed_fractions.pdf", height=3, width=3.8)

max(melt_fitted_zeros_WGD$value)
min(melt_fitted_zeros_WGD$value)

melt_fitted_zeros_WGD <- reshape2::melt(fitted_vals_WGD)
## sort by group
reshuffled_WGD_labels <- unlist(lapply(sapply(c(0,1), function(i) which(TMB_data_zeros_WGD$x[,2] == i)), function(j) j[sample(x = 1:length(j), size = length(j), replace = F)]))
# reshuffled_WGD_labels <- order(TMB_data_zeros_WGD$x[,2])
melt_fitted_zeros_WGD$Var1 <- factor(melt_fitted_zeros_WGD$Var1, levels=rownames(exposures_WGD)[reshuffled_WGD_labels])
fitted_zeros_WGD <- ggplot(melt_fitted_zeros_WGD)+
  geom_tile( aes( x = Var2, y = Var1, fill = value )  )+
  scale_fill_gradient(low="blue", high="red")+
  # geom_point(data = reshape2::melt(exposures_zeros_WGD)  %>% filter(value == 1), aes(x = Var2, y = Var1),
  #            col='white',  shape=18, size=.3)+
  geom_tile(data = reshape2::melt(exposures_zeros_WGD)  %>% filter(value == 1), aes(x = Var2, y = Var1),
             col='white', width=.3)+
  # ggtitle('Observed zeros and fitted probabilities')+
  labs(x='Signatures' , y='Samples', fill='Probability')+
  theme(legend.position = "bottom", axis.title.y=element_blank(),
        axis.text.y=element_blank(), axis.ticks.y=element_blank())+labs(col='Probability of zero')
fitted_zeros_WGD
ggsave("../figures/LN_clusters_exposures_WGD/WGD_fitted_zeros.pdf", width = 2.5, height = 3)

tikzDevice::tikz("../figures/LN_clusters_exposures_WGD/WGD_fitted_zeros.tex")
fitted_zeros_WGD
dev.off()

#-------------
TMB_data_imput1_WGD = list(Y = as(compositions::ilr(impute(exposures_WGD, inputation_value = 1e-2)), 'matrix'),
                       d = 6,
                       n = nrow(exposures_WGD),
                       x =  cbind(1, as.numeric(factor(wgd$V2))-1), ## so that 1=WGD, 0=non-WGD
                       z= give_z_matrix_from_labels(1:nrow(exposures_WGD)))

imput1_WGD <- wrapper_run_TMB_use_nlminb(model = 'mvnbetacor', object = TMB_data_imput1_WGD,
                                     use_nlminb = T)
plot_betas(imput1_WGD, names_cats = paste0('ILR', 1:6))

imput1overdisp_WGD <- wrapper_run_TMB_use_nlminb(model = 'tmb_MVN_with_mean_2_scaledsdpergroup',
                                                 object = TMB_data_imput1_WGD, use_nlminb = T)
plot_betas(imput1overdisp_WGD)

wald_TMB_wrapper(results_ILR_WGD)

table(exposures_zeros_WGD[,4], exposures_zeros_WGD[,6])
fisher.test(table(exposures_zeros_WGD[,4], exposures_zeros_WGD[,6]))

ggplot(cbind.data.frame(exposures_zeros_WGD[,c(4,6)], group=TMB_data_zeros_WGD$x[,2]), aes(x=s4, y=s6))+
  geom_point()+facet_wrap(.~group)
ggplot(cbind.data.frame(exposures_WGD[,c(4,6)], group=TMB_data_zeros_WGD$x[,2],
                        cohort=give_annotation_from_names(rownames(exposures_WGD))),
       aes(x=s4, y=s6, col=cohort))+
  geom_point()+facet_wrap(.~group+cohort)+theme_bw()

pdf("../figures/LN_clusters_exposures_WGD/WGD_estimate_beta_slope_comparison.pdf",
    width = 4, height = 6.5, onefile = F)
grid.arrange(plot_betas(results_ILR_WGD, title='partial ILR cor', names_cats = paste0('ILR', 1:6)),
             plot_betas(results_ILR_WGDnocor, title='partial ILR no cor', names_cats = paste0('ILR', 1:6)),
             plot_betas(imput1_WGD, title = 'imput 1e-2', names_cats = paste0('ILR', 1:6)),
             plot_betas(imput1overdisp_WGD, title = 'imput 1e-2, overdisp', names_cats = paste0('ILR', 1:6)), nrow=4)
dev.off()

plot_lambdas(imput1overdisp_WGD)

#-------------

#-------------

give_amalgamation(as(exposures_WGD[TMB_data_zeros_WGD$x[,2] == 1,], 'matrix'), c(list('s3'), list('s4'), list('s1', 's2', 's5', 's6', 's7')))
give_ternary(give_amalgamation(as(exposures_WGD[TMB_data_zeros_WGD$x[,2] == 1,], 'matrix'), c(list('s3'), list('s4'), list(c('s1', 's2', 's5', 's6', 's7')))))

## 
give_ternary(as(exposures_WGD[TMB_data_zeros_WGD$x[,2] == 1,], 'matrix')[,c(1, 3, 4)])
give_ternary(as(exposures_WGD[TMB_data_zeros_WGD$x[,2] == 0,], 'matrix')[,c(1, 3, 4)])

pdf("../figures/LN_clusters_exposures_WGD/WGD_simplex_groupWGD.pdf", width = 3, height = 3)
give_ternary_v2(as(exposures_WGD[TMB_data_zeros_WGD$x[,2] == 1,], 'matrix')[,c(1, 3, 4)], legend_off=T)
dev.off()
pdf("../figures/LN_clusters_exposures_WGD/WGD_simplex_groupnonWGD.pdf", width = 3, height = 3)
give_ternary_v2(as(exposures_WGD[TMB_data_zeros_WGD$x[,2] == 0,], 'matrix')[,c(1, 3, 4)], legend_off=T)
dev.off()
#-------------

#------------------------------------------------------------------#
#------------------- with new britroc exposures -------------------#
#------------------------------------------------------------------#


new_britroc_exposures <- load_britroc_exposures()

exposures_fig1

dendroimputclr_onlynewbritroc = give_dendrogram_generalised(as(compositions::clr(impute(new_britroc_exposures, 1e-2)), 'matrix'),
                                                            modify_labels=F, keep_only_PDO = F)

dendroimputclr_all_oldbritroc = give_dendrogram_generalised(as(compositions::clr(impute(exposures_fig1, 1e-2)), 'matrix'),
                                                            modify_labels=F, keep_only_PDO = F)

match(rownames(new_britroc_exposures), rownames(exposures_fig1_new_britroc))

rownames(exposures_fig1_new_britroc)[grepl('JBLAB', rownames(exposures_fig1_new_britroc))]
rownames(exposures_fig1_new_britroc)[grepl('IM_', rownames(exposures_fig1_new_britroc))]

## missing from natgen
rownames(new_britroc_exposures)[is.na(match(rownames(new_britroc_exposures), rownames(exposures_fig1)))]
## present in natgen (very few!)
rownames(new_britroc_exposures)[!is.na(match(rownames(new_britroc_exposures), rownames(exposures_fig1)))]
rownames(new_britroc_exposures)[grepl('IM_', rownames(new_britroc_exposures))]

exposures_fig1_new_britroc <- exposures_fig1
exposures_fig1_new_britroc <- exposures_fig1_new_britroc[!grepl('IM_', rownames(exposures_fig1_new_britroc)),]
exposures_fig1_new_britroc <- rbind(exposures_fig1_new_britroc, new_britroc_exposures)

## 265 samples
give_heatmap_and_dendro(dendrogram_arg = dendroimputclr_onlynewbritroc, exposures_arg = new_britroc_exposures,
                        extra_expand = 0, extra_expand_v2 = 0)
## a fraction of high s4 is found together with the high s3 groups, and the outgroup is the high s3 and low s1


## 710 samples
give_heatmap_and_dendro(dendrogram_arg = dendroimputclr_all_oldbritroc, exposures_arg = exposures_fig1,
                        extra_expand = 0, extra_expand_v2 = 0)
## using the old britroc exposures, you get the 1/3 to 2/3 split in high s3 and high s4


## 924 samples
dendroimputclr_all_newbritroc <- give_dendrogram_generalised(as(compositions::clr(impute(exposures_fig1_new_britroc, 1e-2)), 'matrix'),
                            modify_labels=F, keep_only_PDO = F)
give_heatmap_and_dendro(dendrogram_arg = dendroimputclr_all_newbritroc, exposures_arg = exposures_fig1_new_britroc,
                        extra_expand = 0, extra_expand_v2 = 0)
## updated dendrogram: s3 and s4 are split. there are two clades, each with the 1/3 and 2/3 proportions

wrapper_give_heatmap(exposures_fig1_new_britroc, 1e-2)
wrapper_give_heatmap(exposures_fig1_new_britroc, 1e-3)

wrapper_give_heatmap(exposures_fig1[grepl('TCGA', rownames(exposures_fig1)),], 1e-2)
wrapper_give_heatmap(exposures_fig1[grepl('TCGA', rownames(exposures_fig1)),], 1e-2)

annotation_cohort_exposures_fig1 <- give_annotation_from_names(rownames(exposures_fig1))

table(annotation_cohort_exposures_fig1)
table(give_annotation_from_names(rownames(exposures_fig1_new_britroc)))

wrapper_give_heatmap(exposures_fig1[annotation_cohort_exposures_fig1 == 'ICGC',], 1e-2)

give_annotation_from_names(rownames(exposures_fig1))

#------------------------------------------------------------------#
#------------------- Fitting of all exposures ---------------------#
#------------------------------------------------------------------#

TMB_data_imputall = list(Y = as(compositions::ilr(impute(exposures_fig1, inputation_value = 1e-2)), 'matrix'),
                       d = 6,
                       n = nrow(exposures_fig1),
                       x = as.matrix(rep(1, length(.cutree))),
                       z= give_z_matrix_from_labels(1:nrow(exposures_fig1))) ## dummy z

imputall <- wrapper_run_TMB_use_nlminb(model = 'mvnbetacorONEGROUP', object = TMB_data_imputall,
                                     use_nlminb = T, num_covariates = 1)

plot_betas(imputall)


pheatmap::pheatmap(give_rowcolnames(L_to_cov_wrapper_TMB_v2(TMB_obj = imputall, name_cov_par = "cov_RE", d = 6),
                                 paste0('ILR', 1:6)))

cov_all <- L_to_cov_wrapper_TMB_v2(TMB_obj = imputall, name_cov_par = "cov_RE", d = 6)
diag(cov_all) <- exp(python_like_select_name(imputall$par.fixed, 'logs_sd_RE'))

simulated_all <- apply(compositions::ilrInv(mvtnorm::rmvnorm(2000, mean = python_like_select_name(imputall$par.fixed, 'beta'), sigma = cov_all)), 2, as.numeric)
pairs(simulated_all)
pairs(impute(exposures_fig1, inputation_value = 1e-2))

ggplot(cbind.data.frame(rbind(simulated_all[,c(5,7)], impute(exposures_fig1, inputation_value = 1e-2)[,c(5,7)]),
                 col=c(rep('s', nrow(simulated_all)), rep('o', nrow(exposures_fig1)))),
       aes(x=`5`, y=`7`, col=col))+geom_point()+theme_bw()

ggplot()+
  geom_point(data = cbind.data.frame(simulated_all[,c(5,7)]),
             aes(x=`5`, y=`7`), color='blue', alpha=0.06, size=5)+
  geom_point(data = cbind.data.frame(impute(exposures_fig1, inputation_value = 1e-2)[,c(5,7)]),
             aes(x=s5, y=s7), color='red', pch=2, size=1)+
theme_bw()

par(mfrow=c(1,2))
give_ternary(simulated_all[,1:3], main='simulated')
give_ternary(normalise_rw(impute(exposures_fig1, inputation_value = 1e-2)[,1:3]), main='observed')

#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################

## ------------------------------------##
# ## Modelling with logistic normal
# library(TMB)
# 
# source("../../../GlobalDA/code/2_inference_TMB/helper_TMB.R")
# source("../../../other_repos/mcneish/code/models/functions.R")
# # bernoulli for zeros
# TMB::compile("../../../other_repos/mcneish/code/models/tmb_correlated_multinom_2_allFE_b.cpp", "-std=gnu++17")
# dyn.load(dynlib("../../../other_repos/mcneish/code/models/tmb_correlated_multinom_2_allFE_b"))
# # partial ILR for nonzero exposures
# TMB::compile("../../../other_repos/mcneish/code/models/tmb_MVN_partial_ILR_FEb.cpp", "-std=gnu++17")
# dyn.load(dynlib("../../../other_repos/mcneish/code/models/tmb_MVN_partial_ILR_FEb"))
# 
# exposures = exposures_fig1#all_natgen$UpdatedExposures
# ## select only TCGA, as in BriTROC and organoids and ICGC we don't have these zeros
# exposures <- exposures[grepl('TCGA', rownames(exposures)),]
# exposures_zeros = apply(exposures == 0, 2, as.numeric)
# rownames(exposures_zeros) <- rownames(exposures)
# colMeans(exposures_zeros)
# grouping <- split_s3s4tree_df$L1[match(rownames(exposures), split_s3s4tree_df$value)]
# 
# ## Zeros
# ## the levels (High_s3 and High_s4) indicate that the coefficients are for high-s4 wrt high-s3
# TMB_data_zeros = list(Y = exposures_zeros,
#                       num_sigs = ncol(exposures_zeros),
#                       x = cbind(1, as.numeric(factor(grouping, levels=c('High_s3', 'High_s4')))-1))
# 
# xtable::xtable(sapply(list(early=exposures_zeros[TMB_data_zeros$x[,2] == 0,],
#                            late=exposures_zeros[TMB_data_zeros$x[,2] == 1,]), colMeans))
# 
# 
# ## There seems to be a huge change in s3, and perhaps s4 too
# fraction_zeros = sapply(list(Group1=exposures_zeros[TMB_data_zeros$x[,2] == 0,],
#                              Group2=exposures_zeros[TMB_data_zeros$x[,2] == 1,]), colMeans)
# fraction_zeros
# 
# 
# # give_res_report = function(TMB_data, d, model='tmb_correlated_multinom_2_allFE_b'){
# #   TMB_params = list(beta = (matrix(runif(d*2, min = -4, max = 4),
# #                                    nrow = 2, byrow=TRUE))
# #   )
# #   
# #   if(model == 'tmb_MVN_partial_ILR_FEb'){
# #     TMB_params = c(TMB_params, list(logsd = rep(0, d)))
# #   }
# #   
# #   obj <- MakeADFun(data = TMB_data, parameters = TMB_params, DLL=model)
# #   obj$hessian <- TRUE
# #   opt <- do.call("optim", obj)
# #   opt
# #   opt$hessian ## <-- FD hessian from optim
# #   report = sdreport(obj)
# #   
# #   return(list(TMB_params=TMB_params, report=report))
# #   
# # }
# 
# results_zeros = give_res_report(TMB_data_zeros, d=ncol(exposures_zeros))
# results_zeros
# 
# wald_generalised(v = select_slope_2(results_zeros$report$par.fixed,v=F), sigma = results_zeros$report$cov.fixed[c(F,T), c(F,T)])
# 
# pdf("figures/estimate_beta_zeros_bernoulli.pdf", width = 7, height = 3)
# grid.arrange(ggplot(cbind.data.frame(sig=colnames(exposures_zeros), beta=summary(results_zeros$report)[c(T,F),]), aes(x=sig, y=`beta.Estimate`))+
#                geom_point()+
#                geom_errorbar(aes(ymin=`beta.Estimate`-`beta.Std. Error`, ymax=`beta.Estimate`+`beta.Std. Error`), width=.1)+
#                ggtitle('Intercepts'),
#              ggplot(cbind.data.frame(sig=colnames(exposures_zeros), beta=summary(results_zeros$report)[c(F,T),]), aes(x=sig, y=`beta.Estimate`))+
#                geom_point()+
#                geom_errorbar(aes(ymin=`beta.Estimate`-`beta.Std. Error`, ymax=`beta.Estimate`+`beta.Std. Error`), width=.1)+
#                ggtitle('Slopes'), nrow=1
# )
# dev.off()
# 
# 
# # split_df_rw <- function(x, f){
# #   f_levels <- unique(f)
# #   lapply(f_levels, function(f_it) x[which(f == f_it),] )
# # }
# # fitted_probs = sapply(split_df_rw(fitted_vals, TMB_data_zeros$x[,2]), function(i)  colMeans(apply(i, 2, function(j) j ==0)))
# # colnames(fitted_probs) = colnames(fraction_zeros)
# # ## fitted probability of zero vs actual fraction of zeros
# # comparison_probs = (dcast(melt(list(fraction_zeros=fraction_zeros,
# #                                     fitted_probs=fitted_probs), measure.vars=c('L1')), Var1+Var2~L1))
# # 
# # rbind(fraction_zeros, 1-fitted_probs)
# # 
# # ggplot(comparison_probs, aes(x=fraction_zeros, y=fitted_probs, col=Var1, shape=Var2))+geom_point()+
# #   geom_abline(coef=c(0,1), lty='dashed')+ggtitle('Observed fraction of zero and fitted probability')
# # ggsave("../figures/fitted_zeros_observed_fractions.pdf", height=5, width=6)

##---------------------------------------------------------------------------------##

##---------------------------------------------------------------------------------##

exposures_TCGA <- exposures[grepl('TCGA', rownames(exposures)),]
pairs(exposures_TCGA/exposures_TCGA[,1])


pairs(cbind.data.frame(Rest=rowSums(exposures_TCGA[,-c(3,4)]),
                       s3=exposures_TCGA[,3],
                       s4=exposures_TCGA[,4]), col=factor(TMB_data_zeros$x[,2]))

amalgamation = cbind(rowSums(exposures_TCGA[,-c(3,4)]),
                     exposures_TCGA[,3],
                     exposures_TCGA[,4])
colnames(amalgamation) = c('Rest', 's3', 's4')
stopifnot(dim(exposures_TCGA)[1] == dim(TMB_data_zeros$x)[1])
amalgamation1 = amalgamation[TMB_data_zeros$x[,2] == 0,]
amalgamation2 = amalgamation[TMB_data_zeros$x[,2] == 1,]

# png("../figures/simplex_early.png", width = 5, height = 5, units = 'in', res = 300)
pdf("figures/simplex_group1.pdf", width = 5, height = 5)
par(mar=c(0,0,0,0))
TernaryPlot(atip = colnames(amalgamation)[1], btip = colnames(amalgamation)[2], ctip = colnames(amalgamation)[3],
            grid.lines = 0, grid.col = NULL)
dens <- TernaryDensity(amalgamation1, resolution = 10L)

cls_legend = rbind(viridisLite::viridis(48L, alpha = 0.6),
                   seq(from = 0, to = 47, by=1))
legend(x=-0.4,y=1.08,
       fill = cls_legend[1,][c(T,F,F,F,F)],
       legend = round(as.numeric(cls_legend[2,][c(T,F,F,F,F)])/sum(dens['z',]), 2), ncol=5,
       y.intersp=0.8,x.intersp=0.5,text.width=0.1, cex=0.9, bty = "n")
ColourTernary(dens)
TernaryPoints(amalgamation1, col = 'red', pch = '.', cex=5)
TernaryDensityContour(amalgamation1, resolution = 30L)
dev.off()

# png("../figures/simplex_late.png", width = 5, height = 5, units = 'in', res = 300)
pdf("figures/simplex_group2.pdf", width = 5, height = 5)
par(mar=c(0,0,0,0))
TernaryPlot(atip = colnames(amalgamation)[1], btip = colnames(amalgamation)[2], ctip = colnames(amalgamation)[3], grid.lines = 0, grid.col = NULL)
dens <- TernaryDensity(amalgamation2, resolution = 10L)

cls_legend = rbind(viridisLite::viridis(48L, alpha = 0.6),
                   seq(from = 0, to = 47, by=1))
legend(x=-0.4,y=1.08,
       fill = cls_legend[1,][c(T,F,F,F,F)],
       legend = round(as.numeric(cls_legend[2,][c(T,F,F,F,F)])/sum(dens['z',]), 2), ncol=5,
       y.intersp=0.8,x.intersp=0.5,text.width=0.1, cex=0.9, bty = "n")
ColourTernary(dens)
TernaryPoints(amalgamation2, col = 'red', pch = '.', cex=5)
TernaryDensityContour(amalgamation2, resolution = 30L)
dev.off()

#---------------------------------------
# ## done above
# ## do the same for the WGD classes
# wgd1 = read.table("../../../../CDA_in_Cancer/data/TCGA_WGD/Haase2019_TCGA.giScores.wgd.txt", header = T)
# wgd1 = wgd1[!duplicated(wgd1),]
# wgd1 = wgd1[wgd1$cancer_type == "OV",]
# 
# grouping_WGD <- wgd1$wgd[match(rownames(exposures_fig1), wgd1$patient)]
# exposures_WGD <- exposures_fig1[!is.na(grouping_WGD),]
# grouping_WGD <- grouping_WGD[!is.na(grouping_WGD)]
# exposures_zeros_WGD = apply(exposures_WGD == 0, 2, as.numeric)
# 
# ## Zeros
# TMB_data_zeros_WGD = list(Y = exposures_zeros_WGD,
#                           num_sigs = ncol(exposures_zeros_WGD),
#                           x = cbind(1, 1-as.numeric(grouping_WGD)),
#                           z=diag(length(grouping_WGD)))
# 
# xtable::xtable(sapply(list(WGD=exposures_zeros_WGD[TMB_data_zeros_WGD$x[,2] == 0,],
#                            nonWGD=exposures_zeros_WGD[TMB_data_zeros_WGD$x[,2] == 1,]), colMeans))
# 
# 
# ## There seems to be a huge change in s3, and perhaps s4 too
# # fraction_zeros_WGD = sapply(list(nonWGD=exposures_zeros_WGD[TMB_data_zeros_WGD$x[,2] == 0,],
# #                                  WGD=exposures_zeros_WGD[TMB_data_zeros_WGD$x[,2] == 1,]), colMeans)
# # fraction_zeros_WGD
# 
# 
# # results_zeros_WGD = give_res_report(TMB_data_zeros_WGD, d=ncol(exposures_zeros_WGD))
# results_zeros_WGD = wrapper_run_TMB_use_nlminb(model = 'bernoulliFE', object = TMB_data_zeros_WGD)
# 
# results_zeros_WGD
# 
# wald_generalised(v = select_slope_2(results_zeros$par.fixed,v=F), sigma = results_zeros$cov.fixed[c(F,T), c(F,T)])
# 
# plot_betas(results_zeros_WGD)
# pdf("figures/estimate_beta_zeros_bernoulli_WGD.pdf", width = 7, height = 3)
# plot_betas(results_zeros_WGD)
# # grid.arrange(ggplot(cbind.data.frame(sig=colnames(exposures_zeros_WGD), beta=summary(results_zeros_WGD$report)[c(T,F),]), aes(x=sig, y=`beta.Estimate`))+
# #                geom_point()+
# #                geom_errorbar(aes(ymin=`beta.Estimate`-`beta.Std. Error`, ymax=`beta.Estimate`+`beta.Std. Error`), width=.1)+
# #                ggtitle('Intercepts'),
# #              ggplot(cbind.data.frame(sig=colnames(exposures_zeros_WGD), beta=summary(results_zeros_WGD$report)[c(F,T),]), aes(x=sig, y=`beta.Estimate`))+
# #                geom_point()+
# #                geom_errorbar(aes(ymin=`beta.Estimate`-`beta.Std. Error`, ymax=`beta.Estimate`+`beta.Std. Error`), width=.1)+
# #                ggtitle('Slopes'), nrow=1
# # )
# dev.off()
# 
# 
# ## fit
# fitted_logR_WGD <- TMB_data_zeros_WGD$x %*% matrix(results_zeros_WGD$par.fixed, nrow=2)
# fitted_vals_WGD = apply(fitted_logR_WGD, 2, function(i) exp(i)/(1+exp(i)) ) ## transform to probabilities
# colnames(fitted_vals_WGD) = colnames(exposures_zeros_WGD)
# rownames(fitted_vals_WGD) = rownames(exposures_zeros_WGD)
# 
# ggplot(melt(fitted_vals_WGD[order(TMB_data_zeros_WGD$x[,2]),]))+
#   geom_raster( aes( x = Var2, y = Var1, fill = value )  )+
#   scale_fill_gradient(low="blue", high="red")+
#   geom_point(data = melt(exposures_zeros[order(TMB_data_zeros_WGD$x[,2]),])  %>% filter(value == 1), aes(x = Var2, y = Var1),
#              col='white',  shape=18, size=.3)+ggtitle('Observed zeros and\nfitted probabilities')+labs(x='Signatures' , y='Samples')+
#   theme(legend.position = "bottom")
# ggsave("../figures/LN_clusters_exposures_WGD/fitted_zeros_WGD.pdf", width = 3, height = 4)
# 
# 
# ##---------------------------------------------------------------------------------##
# ## Partial ILR
# 
# partialILR_WGD = give_partial_irl(exposures_WGD)
# irl_robust_WGD = as(compositions::ilr(exposures_WGD), 'matrix')
# 
# plot(unlist(irl_robust_WGD), unlist(partialILR_WGD)) ## remember that zeros are treated differently in the inference!
# 
# .keep = (rowSums(partialILR_WGD == 0) < (ncol(exposures) - 2) )
# partialILR_WGD = partialILR_WGD[.keep,]  ## if a sample has too many zeros (i.e. all but 1 or 2) remove
# colnames(partialILR_WGD) = paste0('ILR', 1:ncol(partialILR_WGD))
# 
# TMB_data_partialILR_WGD = list(Y = partialILR_WGD,
#                                d = ncol(partialILR_WGD),
#                                n = sum(.keep),
#                                x = cbind(1,1-as.numeric(grouping_WGD[.keep])-1))
# 
# results_ILR_WGD = give_res_report(TMB_data_partialILR_WGD, d = (ncol(partialILR_WGD)), model='tmb_MVN_partial_ILR_FEb')
# results_ILR_WGD
# 
# pdf("figures/estimate_beta_ILR_WGD.pdf", width = 7, height = 3)
# grid.arrange(ggplot(cbind.data.frame(sig=colnames(partialILR_WGD),
#                                      beta=python_like_select_rownames(summary(results_ILR_WGD$report), 'beta')[c(T,F),]),
#                     aes(x=sig, y=`beta.Estimate`))+
#                geom_point()+
#                geom_errorbar(aes(ymin=`beta.Estimate`-`beta.Std. Error`, ymax=`beta.Estimate`+`beta.Std. Error`), width=.1)+
#                ggtitle('Intercepts'),
#              ggplot(cbind.data.frame(sig=colnames(partialILR_WGD),
#                                      beta=python_like_select_rownames(summary(results_ILR_WGD$report), 'beta')[c(F,T),]),
#                     aes(x=sig, y=`beta.Estimate`))+
#                geom_point()+
#                geom_errorbar(aes(ymin=`beta.Estimate`-`beta.Std. Error`, ymax=`beta.Estimate`+`beta.Std. Error`), width=.1)+
#                ggtitle('Slopes'), nrow=1)
# dev.off()
# 
# ggplot(melt(list(group1=exposures[TMB_data_partialILR_WGD$x[,2] == 0,],
#                  group2=exposures[TMB_data_partialILR_WGD$x[,2] == 1,])),
#        aes(x=value, col=L1))+geom_density()+facet_wrap(.~Var2, nrow=1, scales = "free_x")+
#   scale_color_brewer(palette = "Set2")
# ggsave("figures/densities_change_WGD.pdf", width = 10, height = 2.4)
# 
# wald_generalised(v = select_slope_2(results_ILR$report$par.fixed,v=F), sigma = results_ILR$report$cov.fixed[c(F,T), c(F,T)])
# 
# ## Put the results of the two models together
# results_zeros_ILR_WGD <- cbind.data.frame(label=colnames(exposures_zeros_WGD),
#                                           partialILR=rbind(NA, python_like_select_rownames(summary(results_ILR_WGD$report), 'beta')[c(F,T),])[match(colnames(exposures_zeros_WGD), colnames(exposures_WGD)),],
#                                           zeros=python_like_select_rownames(summary(results_zeros_WGD$report), 'beta')[c(F,T),])
# ggplot(results_zeros_ILR_WGD,
#        aes(x=partialILR.Estimate, y=-zeros.Estimate, label=label))+geom_point()+
#   geom_hline(yintercept = 0, col='blue', lty='dashed', size=0.3)+
#   geom_vline(xintercept = 0, col='blue', lty='dashed', size=0.3)+
#   geom_errorbar(aes(xmin=`partialILR.Estimate`-`partialILR.Std. Error`,xmax=`partialILR.Estimate`+`partialILR.Std. Error`),
#                 width=.1)+
#   geom_errorbar(aes(ymin=-`zeros.Estimate`-`zeros.Std. Error`, ymax=-`zeros.Estimate`+`zeros.Std. Error`), width=.1)+
#   geom_label_repel()+facet_wrap(.~label, nrow=1, scales = "free")+
#   ggtitle('Summary of both models')
# ggsave("figures/summary_both_sets_coefs_WGD.pdf", width = 10, height = 2.4)

##---------------------------------------------------------------------------------##

amalgamation_WGD = cbind(rowSums(exposures_WGD[,-c(3,4)]),
                         exposures_WGD[,3],
                         exposures_WGD[,4])
colnames(amalgamation_WGD) = c('Rest', 's3', 's4')
amalgamation1_WGD = amalgamation_WGD[TMB_data_zeros_WGD$x[,2] == 0,]
amalgamation2_WGD = amalgamation_WGD[TMB_data_zeros_WGD$x[,2] == 1,]

# png("../figures/simplex_early.png", width = 5, height = 5, units = 'in', res = 300)
pdf("figures/simplex_group1_WGD.pdf", width = 5, height = 5)
par(mar=c(0,0,0,0))
TernaryPlot(atip = colnames(amalgamation_WGD)[1], btip = colnames(amalgamation_WGD)[2], ctip = colnames(amalgamation_WGD)[3],
            grid.lines = 0, grid.col = NULL)
dens <- TernaryDensity(amalgamation1_WGD, resolution = 10L)

cls_legend = rbind(viridisLite::viridis(48L, alpha = 0.6),
                   seq(from = 0, to = 47, by=1))
legend(x=-0.4,y=1.08,
       fill = cls_legend[1,][c(T,F,F,F,F)],
       legend = round(as.numeric(cls_legend[2,][c(T,F,F,F,F)])/sum(dens['z',]), 2), ncol=5,
       y.intersp=0.8,x.intersp=0.5,text.width=0.1, cex=0.9, bty = "n")
ColourTernary(dens)
TernaryPoints(amalgamation1_WGD, col = 'red', pch = '.', cex=5)
TernaryDensityContour(amalgamation1_WGD, resolution = 30L)
dev.off()

# png("../figures/simplex_late.png", width = 5, height = 5, units = 'in', res = 300)
pdf("figures/simplex_group2_WGD.pdf", width = 5, height = 5)
par(mar=c(0,0,0,0))
TernaryPlot(atip = colnames(amalgamation_WGD)[1], btip = colnames(amalgamation_WGD)[2], ctip = colnames(amalgamation_WGD)[3], grid.lines = 0, grid.col = NULL)
dens <- TernaryDensity(amalgamation2_WGD, resolution = 10L)

cls_legend = rbind(viridisLite::viridis(48L, alpha = 0.6),
                   seq(from = 0, to = 47, by=1))
legend(x=-0.4,y=1.08,
       fill = cls_legend[1,][c(T,F,F,F,F)],
       legend = round(as.numeric(cls_legend[2,][c(T,F,F,F,F)])/sum(dens['z',]), 2), ncol=5,
       y.intersp=0.8,x.intersp=0.5,text.width=0.1, cex=0.9, bty = "n")
ColourTernary(dens)
TernaryPoints(amalgamation2_WGD, col = 'red', pch = '.', cex=5)
TernaryDensityContour(amalgamation2_WGD, resolution = 30L)
dev.off()

table(give_annotation_from_names(rownames(exposures_fig1)))

save.image("../image_code_WGD_modelling.RData")

