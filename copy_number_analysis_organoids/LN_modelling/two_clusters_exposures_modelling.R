rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(TMB)
library(ggplot2)
library(cowplot)
library(tikzDevice)
library(ggrepel)

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
  scale_y_continuous(expand=c(0,0))+lims(y=c(-1.6, 1))
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

# system("open ../../../../CDA_in_Cancer/text/thesis/figures/")
tikz( '../../../../CDA_in_Cancer/text/thesis/figures/general_compositional/partial_fig1_2.tex', height = 4, width = 6)
e_v2
dev.off()

#----------------------------
source("../../../mcneish/code/models/functions.R")
source("../../../mcneish/code/models/helper_TMB.R")
# bernoulli for zeros
TMB::compile("../../../mcneish/code/models/tmb_correlated_multinom_2_allFE_b.cpp", "-std=gnu++17")
dyn.load(dynlib("../../../mcneish/code/models/tmb_correlated_multinom_2_allFE_b"))
# partial ILR for nonzero exposures
TMB::compile("../../../mcneish/code/models/tmb_MVN_partial_ILR_FEb.cpp", "-std=gnu++17")
dyn.load(dynlib("../../../mcneish/code/models/tmb_MVN_partial_ILR_FEb"))

exposures_zeros = apply(exposures_fig1 == 0, 2, as.numeric)

X <- cbind(1, as.numeric(.cutree == 2))

## Zeros
TMB_data_zeros = list(Y = exposures_zeros,
                      num_sigs = ncol(exposures_zeros),
                      x = X)

fraction_zeros = sapply(list(cluster1=exposures_zeros[TMB_data_zeros$x[,2] == 0,],
                             cluster2=exposures_zeros[TMB_data_zeros$x[,2] == 1,]), colMeans)
fraction_zeros
ggplot(melt(fraction_zeros), aes(x=Var1, y=value, fill=Var2))+geom_bar(stat = "identity",
                                                                       position=position_dodge(width=0.5))

give_res_report = function(TMB_data, d, model='tmb_correlated_multinom_2_allFE_b'){
  TMB_params = list(beta = (matrix(runif(d*2, min = -4, max = 4),
                                   nrow = 2, byrow=TRUE))
  )
  
  if(model == 'tmb_MVN_partial_ILR_FEb'){
    TMB_params = c(TMB_params, list(logsd = rep(0, d)))
  }
  
  obj <- MakeADFun(data = TMB_data, parameters = TMB_params, DLL=model)
  obj$hessian <- TRUE
  opt <- do.call("optim", obj)
  opt
  opt$hessian ## <-- FD hessian from optim
  report = sdreport(obj)
  
  return(list(TMB_params=TMB_params, report=report))
  
}

results_zeros = give_res_report(TMB_data_zeros, d=ncol(exposures_zeros))
results_zeros

wald_generalised(v = select_slope_2(results_zeros$report$par.fixed,v=F), sigma = results_zeros$report$cov.fixed[c(F,T), c(F,T)])

pdf("../figures/LN_clusters_exposures_WGD/estimate_beta_slope_zeros.pdf", width = 6, height = 3)
grid.arrange(ggplot(cbind.data.frame(sig=colnames(exposures_zeros), beta=summary(results_zeros$report)[c(T,F),]), aes(x=sig, y=`beta.Estimate`))+
  geom_point()+
  geom_errorbar(aes(ymin=`beta.Estimate`-`beta.Std. Error`, ymax=`beta.Estimate`+`beta.Std. Error`), width=.1)+
  ggtitle('Intercepts')+theme_bw(),
  ggplot(cbind.data.frame(sig=colnames(exposures_zeros), beta=summary(results_zeros$report)[c(F,T),]), aes(x=sig, y=`beta.Estimate`))+
  geom_point()+
  geom_errorbar(aes(ymin=`beta.Estimate`-`beta.Std. Error`, ymax=`beta.Estimate`+`beta.Std. Error`), width=.1)+
  ggtitle('Slopes')+theme_bw(), ncol=2)
dev.off()


fitted_logR <- TMB_data_zeros$x %*% matrix(results_zeros$report$par.fixed, nrow=2)
fitted_vals = apply(fitted_logR, 2, function(i) exp(i)/(1+exp(i)) ) ## transform to probabilities
colnames(fitted_vals) = colnames(exposures_zeros)
rownames(fitted_vals) = rownames(exposures_zeros)

ggplot(melt(fitted_vals))+
  geom_raster( aes( x = Var2, y = Var1, fill = value )  )+
  scale_fill_gradient(low="blue", high="red")+
  geom_point(data = melt(exposures_zeros)  %>% filter(value == 1), aes(x = Var2, y = Var1),
             col='white',  shape=18, size=.3)+ggtitle('Observed zeros and\nfitted probabilities')+labs(x='Signatures' , y='Samples')+
  theme(legend.position = "bottom")
ggsave("../figures/LN_clusters_exposures_WGD/fitted_zeros.pdf", width = 3, height = 4)

fitted_probs = t(fitted_vals[c(1, sum(TMB_data_zeros$x[,2] == 0)+1),]); colnames(fitted_probs) = c('cluster1', 'cluster2')
## fitted probability of zero vs actual fraction of zeros
comparison_probs = (dcast(melt(list(fraction_zeros=fraction_zeros,
                                    fitted_probs=fitted_probs), measure.vars=c('L1')), Var1+Var2~L1))

rbind(fraction_zeros, 1-fitted_probs)

ggplot(comparison_probs, aes(x=fraction_zeros, y=fitted_probs, col=Var1, shape=Var2))+geom_point()+
  geom_abline(coef=c(0,1), lty='dashed')+ggtitle('Observed fraction of zero and fitted probability')+
  theme_bw()
ggsave("../figures/LN_clusters_exposures_WGD/fitted_zeros_observed_fractions.pdf", height=5, width=6)

##-----------------

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
                           x = cbind(1, as.numeric(.cutree[match(rownames(partialILR), names(.cutree))] == 1)))

results_ILR = give_res_report(TMB_data_partialILR, d = (ncol(partialILR)), model='tmb_MVN_partial_ILR_FEb')
results_ILR
results_ILR$report$pdHess

#----
pdf("../figures/LN_clusters_exposures_WGD/estimate_beta_slope_ILR.pdf", width = 6, height = 3)
grid.arrange(ggplot(cbind.data.frame(sig=colnames(partialILR), beta=python_like_select_rownames(summary(results_ILR$report), 'beta')[c(T,F),]),
       aes(x=sig, y=`beta.Estimate`))+
  geom_point()+
  geom_errorbar(aes(ymin=`beta.Estimate`-`beta.Std. Error`, ymax=`beta.Estimate`+`beta.Std. Error`), width=.1)+
  ggtitle('Intercepts')+theme_bw(),
  ggplot(cbind.data.frame(sig=colnames(partialILR), beta=python_like_select_rownames(summary(results_ILR$report), 'beta')[c(F,T),]),
       aes(x=sig, y=`beta.Estimate`))+
  geom_point()+
  geom_errorbar(aes(ymin=`beta.Estimate`-`beta.Std. Error`, ymax=`beta.Estimate`+`beta.Std. Error`), width=.1)+
  ggtitle('Slopes')+theme_bw(), ncol=2)
dev.off()

exp(python_like_select_rownames(summary(results_ILR$report), 'logsd')[,1])

colMeans(exposures_fig1 == 0)
colMeans(partialILR == 0)

ggplot(melt(list(cluster1=exposures_fig1[TMB_data_partialILR$x[,2] == 0,],
                 cluster2=exposures_fig1[TMB_data_partialILR$x[,2] == 1,])),
       aes(x=value, col=L1))+geom_density()+facet_wrap(.~Var2, nrow=1, scales = "free_x")+
  scale_color_brewer(palette = "Set2")
ggsave("../figures/LN_clusters_exposures_WGD/densities_change.pdf", width = 10, height = 2.4)

wald_generalised(v = select_slope_2(results_ILR$report$par.fixed,v=F), sigma = results_ILR$report$cov.fixed[c(F,T), c(F,T)])

## Put the results of the two models together
ggplot(cbind.data.frame(label=colnames(exposures_zeros),
                        partialILR=rbind(NA, python_like_select_rownames(summary(results_ILR$report), 'beta')[c(F,T),])[match(colnames(exposures_zeros), colnames(exposures_fig1)),],
                        zeros=python_like_select_rownames(summary(results_zeros$report), 'beta')[c(F,T),]),
       aes(x=partialILR.Estimate, y=-zeros.Estimate, label=label))+geom_point()+
  geom_hline(yintercept = 0, col='blue', lty='dashed', size=0.3)+
  geom_vline(xintercept = 0, col='blue', lty='dashed', size=0.3)+
  geom_errorbar(aes(xmin=`partialILR.Estimate`-`partialILR.Std. Error`,xmax=`partialILR.Estimate`+`partialILR.Std. Error`), width=.1)+
  geom_errorbar(aes(ymin=-`zeros.Estimate`-`zeros.Std. Error`, ymax=-`zeros.Estimate`+`zeros.Std. Error`), width=.1)+
  geom_label()+facet_wrap(.~label, nrow=1)+
  ggtitle('Summary of both models')+theme_bw()
ggsave("../figures/LN_clusters_exposures_WGD/summary_both_sets_coefs.pdf", width = 8, height = 2.4)


##----
## Pending: how to go back to fitted values from partial ILR results? I guess we have to put NA in many places
## (but I don't think there is much a thing as a fitted value... unless I put a zero everywhere else)

fitted_logR_partialILR <- TMB_data_partialILR$x %*% matrix(results_ILR$report$par.fixed[grepl('beta', names(results_ILR$report$par.fixed))], nrow=2)
fitted_vals_partialILR = apply(cbind(fitted_logR_partialILR, 0), 2, function(i) exp(i)/(1+exp(i)) ) ## transform to probabilities
colnames(fitted_vals_partialILR) = colnames(exposures_fig1) #paste0('ILR', 1:ncol(partialILR))
colnames(exposures_fig1)
rownames(fitted_vals_partialILR) = rownames(TMB_data_partialILR$Y)
fitted_vals_partialILR


fitted_vals_partialILR

.sftmax_partialILR <- cbind(partialILR, 0)
colnames(.sftmax_partialILR) <- colnames(exposures_fig1)
partialILR_melt <- (melt(apply(.sftmax_partialILR, 2, function(i) exp(i)/(1+exp(i)) ) ))
fitted_vals_partialILR_melt <- (melt(apply(fitted_vals_partialILR, 2, function(i) exp(i)/(1+exp(i)) ) ))

ggplot(cbind.data.frame(partialILR_melt,
                        values_observed=fitted_vals_partialILR_melt[match(paste0(partialILR_melt$Var1, partialILR_melt$Var2),
      paste0(fitted_vals_partialILR_melt$Var1, fitted_vals_partialILR_melt$Var2)),3]), aes(x=value, y=values_observed))+
  geom_point(alpha=0.2)+facet_wrap(.~Var2)+theme_bw()
ggsave("../figures/LN_clusters_exposures_WGD/fitted_partialILR.pdf", width = 3, height = 4)

partialILR_melt[partialILR_melt$Var2 == 's7',]



