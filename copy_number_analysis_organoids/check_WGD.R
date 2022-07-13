rm(list = ls())
### Files related to WGD clades
###            ../check_WGD.R
###            TCGA_HGSOC_dendrogram.R
###            two_cluster_exposures_modelling.R

setwd("~/Documents/PhD/other_repos/Vias_Brenton/figures/")

## for fitting the model, check two_clusters_exposures_modelling.R

require(reshape2)
require(ggplot2)
require(ggrepel)
require(cowplot)
require(gridExtra)
require(latex2exp)
library(cowplot)
library(ggpubr)
library(survminer)
library(survival)
library(CompSign)
library(caret)
library(readxl)
library(ggdendro)

source("../copy_number_analysis_organoids/helper_functions.R")
# source("../../../../GlobalDA/code/2_inference_TMB/helper_TMB.R")
# source("../../../mcneish/code/models/helper_TMB.R")
source("../../../other_repos/britroc-1/code/models/helper/functions.R")


exposures <- readRDS("../copy_number_analysis_organoids/robjects/exposures.RDS")
dendrograminputclr <- readRDS("../copy_number_analysis_organoids/robjects/dendrograminputclr.RDS")
heatmapinputclr <- readRDS("../copy_number_analysis_organoids/robjects/heatmapinputclr.RDS")

heatmapinputclr

wgd_hasse <- read.table("../../../CDA_in_Cancer/data/TCGA_WGD/Haase2019_TCGA.giScores.wgd.txt", header = T)
wgd_icgc <- read_xlsx("../../../CDA_in_Cancer/data/pcawg_supplementary_tables/41586_2019_1907_MOESM4_ESM.xlsx")
wgd_icgc_withWGD <- sapply(unique(wgd_icgc$samplename), function(i) any(wgd_icgc[wgd_icgc$samplename == i,'type'] == 'WGD'))
wgd_icgc_withWGD <- cbind.data.frame(sample=names(wgd_icgc_withWGD), WGD=wgd_icgc_withWGD)

heatmap_data <- heatmapinputclr$data
heatmap_data = heatmap_data[heatmap_data$Var2 %in% wgd_hasse$patient,]
head(heatmap_data)

heatmap_data$WGD = wgd_hasse$wgd[match(heatmap_data$Var2, wgd_hasse$patient)]

heatmap_data_ICGC = heatmapinputclr$data
heatmap_data_ICGC$WGD = wgd_icgc_withWGD$WGD[match(heatmap_data_ICGC$Var2, wgd_icgc_withWGD$sample)]
heatmap_data_ICGC = heatmap_data_ICGC[!is.na(heatmap_data_ICGC$WGD),]

ggplot(heatmap_data, aes(x=Var2, fill=Var1, y=value))+geom_bar(stat = "identity")+
  scale_fill_brewer(palette="Dark2")+facet_wrap(.~ifelse(WGD, 'non-WGD', 'WGD'), nrow=2)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.title = element_blank())+
  # ggtitle('WGD')+
  guides(fill='none')
ggsave("../copy_number_analysis_organoids/figures/WGD_exposures.pdf", height = 2.5, width = 2.5)

## adding the information about cohorts
heatmap_data

heatmap_data$WGD <- ifelse(heatmap_data$WGD, 'WGD', 'non-WGD')

tikzDevice::tikz("../copy_number_analysis_organoids/figures/WGD_exposures_tikz.tex", height = 3, width = 2.5)

ggplot(heatmap_data, aes(x=Var2, fill=Var1, y=value))+geom_bar(stat = "identity")+
  scale_fill_brewer(palette="Dark2")+facet_wrap(.~WGD, nrow=2)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.title = element_blank())+
  ggtitle('TCGA exposures')+guides(fill='none')+labs(y='Exposure')
dev.off()

system("open ../copy_number_analysis_organoids/figures/")

ggplot(heatmap_data_ICGC, aes(x=Var2, fill=Var1, y=value))+geom_bar(stat = "identity")+
  scale_fill_brewer(palette="Dark2")+facet_wrap(.~ifelse(WGD, 'non-WGD', 'WGD'), nrow=2)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.title = element_blank())+
  # ggtitle('WGD (ICGC samples)')+
  guides(fill='none')
ggsave("../copy_number_analysis_organoids/figures/WGD_exposures_ICGC.pdf", height = 2.5, width = 2.5)

heatmap_data_ICGC$WGD <- ifelse(heatmap_data_ICGC$WGD, 'WGD', 'non-WGD')

tikzDevice::tikz("../copy_number_analysis_organoids/figures/WGD_exposures_ICGC_tikz.tex", height = 3, width = 2.5)
ggplot(heatmap_data_ICGC, aes(x=Var2, fill=Var1, y=value))+geom_bar(stat = "identity")+
  scale_fill_brewer(palette="Dark2")+facet_wrap(.~WGD, nrow=2)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.title = element_blank(),
        panel.background = element_blank())+
  ggtitle('ICGC exposures')+guides(fill='none')+labs(y='Exposure')
dev.off()



## SVM

library(e1071)
give_svm_res <- function(df_exposures){
  exposures <- dcast(df_exposures, formula =  Var2 ~ Var1, value.var = 'value')
  rownames(exposures) <- exposures[,1]; exposures <- exposures[,-1]
  clr <- compositions::clr(impute(exposures, 1e-2))
  sample_WGD <- data.frame(df_exposures) %>% dplyr::filter(Var1 == 's1') %>% dplyr::select(WGD) %>% unlist()
  training_set_idx <- sample(1:nrow(clr), size = floor(.7*nrow(clr)))
  training_set = clr[training_set_idx,]
  validation_set = clr[!(1:nrow(clr) %in% training_set_idx),]
  training_dat = data.frame(training_set, y = as.factor(sample_WGD[training_set_idx]))
  validation_dat = data.frame(validation_set, y = as.factor(sample_WGD[!(1:nrow(clr) %in% training_set_idx)]))
  svmfit = svm(y ~ ., data = training_dat, kernel = "linear", cost = 10, scale = FALSE)
  print(svmfit)
  
  validation_predicted = predict(svmfit, validation_dat)
  
  print(confusionMatrix(validation_predicted, validation_dat$y))
  return(list(svmfit=svmfit, confusionMatrix=confusionMatrix(validation_predicted, validation_dat$y),
              training_dat=training_dat))
}

length(unique(heatmap_data_ICGC$Var2))
length(unique(heatmap_data$Var2))
table(heatmap_data_ICGC[!duplicated(heatmap_data_ICGC$Var2),'WGD'])
table(heatmap_data[!duplicated(heatmap_data$Var2),'WGD'])

svmfit_rep <- replicate(give_svm_res(df_exposures = heatmap_data), n = 20, simplify = F) ## TCGA
svmfit_ICGC_rep <- replicate(give_svm_res(df_exposures = heatmap_data_ICGC), n=20, simplify = F) ## ICGC
svmfit_both_rep <- replicate(give_svm_res(df_exposures = rbind(heatmap_data, heatmap_data_ICGC)), n=20, simplify = F) ## both

saveRDS(list(svmfit_rep=svmfit_rep, svmfit_ICGC_rep=svmfit_ICGC_rep, svmfit_both_rep=svmfit_both_rep),
        file = "../copy_number_analysis_organoids/robjects/SVM_WGD.RDS")

mean(sapply(svmfit_rep, function(j) j$confusionMatrix$byClass['Sensitivity']))
mean(sapply(svmfit_rep, function(j) j$confusionMatrix$byClass['Specificity']))
mean(sapply(svmfit_rep, function(j) j$confusionMatrix$byClass['Balanced Accuracy']))
mean(sapply(svmfit_ICGC_rep, function(j) j$confusionMatrix$byClass['Balanced Accuracy']))


svmfit <- svmfit_rep[[1]]$svmfit
training_dat <- svmfit_rep[[1]]$training_dat
plot(x = svmfit, data = training_dat, formula = s1 ~ s4, fill = TRUE, grid = 50, slice = list())
pdf("../copy_number_analysis_organoids/figures/WGD_exposures_svm1.pdf", height = 4, width = 4)
plot(x = svmfit, data = training_dat, formula = s4 ~ s3, fill = TRUE, grid = 50, slice = list())
dev.off()
pdf("../copy_number_analysis_organoids/figures/WGD_exposures_svm2.pdf", height = 4, width = 4)
plot(x = svmfit, data = training_dat, formula = s3 ~ s5, fill = TRUE, grid = 50, slice = list())
dev.off()
plot(x = svmfit, data = training_dat, formula = s5 ~ s6, fill = TRUE, grid = 50, slice = list())

pheatmap::pheatmap(svmfit$SV)
pheatmap::pheatmap(svmfit$coefs, cluster_cols = F, cluster_rows = F)

add_rownames <- function(i){
  rownames(i) <- paste0('Dim', 1:nrow(i))
  i
}

myplotSVM <- function (x, data, formula = NULL, fill = TRUE, grid = 50, slice = list(), 
                       symbolPalette = palette(), svSymbol = "x", dataSymbol = "o", main='SVM Classification Plot',
                       ...) {
  if (x$type < 3) {
    if (is.null(formula) && ncol(data) == 3) {
      formula <- formula(delete.response(terms(x)))
      formula[2:3] <- formula[[2]][2:3]
    }
    if (is.null(formula)) 
      stop("missing formula.")
    if (fill) {
      sub <- model.frame(formula, data)
      xr <- seq(min(sub[, 2]), max(sub[, 2]), length.out = grid)
      yr <- seq(min(sub[, 1]), max(sub[, 1]), length.out = grid)
      l <- length(slice)
      if (l < ncol(data) - 3) {
        slnames <- names(slice)
        slice <- c(slice, rep(list(0), ncol(data) - 3 - 
                                l))
        names <- labels(delete.response(terms(x)))
        names(slice) <- c(slnames, names[!names %in% 
                                           c(colnames(sub), slnames)])
      }
      for (i in names(which(vapply(data, is.factor, NA)))) if (!is.factor(slice[[i]])) {
        levs <- levels(data[[i]])
        lev <- if (is.character(slice[[i]])) 
          slice[[i]]
        else levs[1]
        fac <- factor(lev, levels = levs)
        if (is.na(fac)) 
          stop(paste("Level", dQuote(lev), "could not be found in factor", 
                     sQuote(i)))
        slice[[i]] <- fac
      }
      lis <- c(list(yr), list(xr), slice)
      names(lis)[1:2] <- colnames(sub)
      new <- expand.grid(lis)[, labels(terms(x))]
      preds <- predict(x, new)
      filled.contour(xr, yr, matrix(as.numeric(preds), 
                                    nrow = length(xr), byrow = TRUE), plot.axes = {
                                      axis(1)
                                      axis(2)
                                      colind <- as.numeric(model.response(model.frame(x, 
                                                                                      data)))
                                      dat1 <- data[-x$index, ]
                                      dat2 <- data[x$index, ]
                                      coltmp1 <- symbolPalette[colind[-x$index]]
                                      coltmp2 <- symbolPalette[colind[x$index]]
                                      points(formula, data = dat1, pch = dataSymbol, 
                                             col = coltmp1)
                                      points(formula, data = dat2, pch = svSymbol, 
                                             col = coltmp2)
                                    }, levels = 1:(length(levels(preds)) + 1), key.axes = axis(4, 
                                                                                               1:(length(levels(preds))) + 0.5, labels = levels(preds), 
                                                                                               las = 3), plot.title = title(main = main, 
                                                                                                                            xlab = names(lis)[2], ylab = names(lis)[1]), 
                     ...)
    }
    else {
      plot(formula, data = data, type = "n", ...)
      colind <- as.numeric(model.response(model.frame(x, 
                                                      data)))
      dat1 <- data[-x$index, ]
      dat2 <- data[x$index, ]
      coltmp1 <- symbolPalette[colind[-x$index]]
      coltmp2 <- symbolPalette[colind[x$index]]
      points(formula, data = dat1, pch = dataSymbol, col = coltmp1)
      points(formula, data = dat2, pch = svSymbol, col = coltmp2)
      invisible()
    }
  }
}

previous_mar <- par()$mar

tikzDevice::tikz("../copy_number_analysis_organoids/figures/WGD_exposures_SVM.tex", height = 2.0, width = 3.1)
par(mar=c(4,4,1,2))
myplotSVM(x = svmfit, data = training_dat, formula = s4 ~ s3, fill = TRUE, grid = 50, slice = list(),
     col=c('#e7feff', '#ffdead'), symbolPalette=c('black', 'red'), main='')
dev.off()

tikzDevice::tikz("../copy_number_analysis_organoids/figures/WGD_exposures_SVM_ICGC.tex", height = 2.0, width = 3.1)
par(mar=c(4,4,1,2))
myplotSVM(x = svmfit_ICGC_rep[[1]]$svmfit, data = svmfit_ICGC_rep[[1]]$training_dat, formula = s4 ~ s3, fill = TRUE, grid = 50, slice = list(),
          col=c('#e7feff', '#ffdead'), symbolPalette=c('black', 'red'), main='')
dev.off()

pheatmap::pheatmap(svmfit$SV, annotation_row = add_rownames(svmfit$coefs))

# ## TRYING TO FIND ALL THE EXPOSURES OF ALL BRITROC SAMPLES (?????)
# ## predict status of britroc samples
# CN_features<-readRDS("../../britroc-cnsignatures-bfb69cd72c50/data/feat_sig_mat.rds")#extractCopynumberFeatures(hq_CN,cores=num_cores)
# CN_components<-readRDS("../../britroc-cnsignatures-bfb69cd72c50/data/component_parameters.rds")#extractCopynumberFeatures(hq_CN,cores=num_cores)
# readRDS("../../britroc-cnsignatures-bfb69cd72c50/manuscript_Rmarkdown/data/")#extractCopynumberFeatures(hq_CN,cores=num_cores)
# source("../../britroc-cnsignatures-bfb69cd72c50/main_functions.R")
# source("../../britroc-cnsignatures-bfb69cd72c50/helper_functions.R")
# 
# britroc_sample_component_matrix<-generateSampleByComponentMatrix(CN_features,CN_components,cores=1,subcores=num_cores)
# 
# 
# quantifySignatures<-function(sample_by_component,component_by_signature=NULL){
#   if(is.null(component_by_signature))
#   {
#     component_by_signature<-readRDS(paste(this_path,"data/feat_sig_mat.rds",sep="/"))
#   }
#   signature_by_sample<-YAPSA::LCD(t(sample_by_component),
#                                   YAPSA:::normalize_df_per_dim(component_by_signature,2))
#   signature_by_sample<-normaliseMatrix(signature_by_sample)
#   signature_by_sample
# }
# 
# 
# ### WHERE ARE THE BRITROC SAMPLES IN THIS DENDROGRAM???
# sig_data = readRDS("../copy_number_analysis_organoids/data/sig_data_unorm.RDS")
# sig_data = readRDS("../../britroc-cnsignatures-bfb69cd72c50/manuscript_Rmarkdown/data/sig")
# View(rownames(sig_data))
# dim(exposures_all)
# sum(grepl('^IM_', rownames(exposures_all)))
# sum(grepl('^TCGA-', rownames(exposures_all)))
# sum(grepl('^PDO', rownames(exposures_all)))

exposures_all <- dcast(heatmapinputclr$data, formula =  Var2 ~ Var1, value.var = 'value')
rownames(exposures_all) <- exposures_all[,1]; exposures_all <- exposures_all[,-1]
clr_all <- compositions::clr(impute(exposures_all, 1e-2))

predict_im_britroc <- predict(svmfit_rep[[1]]$svmfit, data.frame(clr_all[grepl('^IM_', rownames(clr_all)),]))
table(predict_im_britroc) ## subset of britroc samples from 2018 natgen


## get exposures of new britroc samples
# give_exposures <- function(){
#   load("../../britroc-1/data/britroc_30kb_signature_data.rds")
#   exposures = t(sig_quants)
#   exposures
# }

exposures_britroc <- load_britroc_exposures("../../../other_repos/britroc-1/data/britroc_30kb_signature_data.rds")
patient.meta <- load_britroc_meta("../../../other_repos/britroc-1/data/britroc_30kb_signature_data.rds")
patient.meta$SAMPLE_ID <- as.character(patient.meta$SAMPLE_ID)
patient.meta$PATIENT_ID <- as.character(patient.meta$PATIENT_ID)
clr_britroc<- compositions::clr(impute(exposures_britroc, 1e-2))
clr_britroc_predict1 <- predict(svmfit_rep[[1]]$svmfit, data.frame(clr_britroc))
clr_britroc_predict2 <- predict(svmfit_ICGC_rep[[1]]$svmfit, data.frame(clr_britroc))

prediction_WGD_britroc <- cbind.data.frame(svm_TCGA = clr_britroc_predict1, svm_ICGC = clr_britroc_predict2)
table(prediction_WGD_britroc)
saveRDS(prediction_WGD_britroc, "../copy_number_analysis_organoids/robjects/prediction_WGD_britroc.RDS")

## load early
late_samples <- read.csv("../../mcneish/data/sig_data_by_samples_lena.csv")
late_samples <- late_samples[late_samples$stage == "late_stage_cohort",]
late_samples <- late_samples[,paste0('s', 1:7)]

early_samples <- read.csv("../../mcneish/data/sig_data_by_samples_lena.csv")
early_samples <- early_samples[early_samples$stage == "early_stage_cohort",]
early_samples <- early_samples[,paste0('s', 1:7)]

clr_early_predict1 <- predict(svmfit_rep[[1]]$svmfit, data.frame(compositions::clr(impute(early_samples, 1e-2))))
clr_early_predict2 <- predict(svmfit_ICGC_rep[[1]]$svmfit, data.frame(compositions::clr(impute(early_samples, 1e-2))))
clr_early_predictboth <- predict(svmfit_both_rep[[1]]$svmfit, data.frame(compositions::clr(impute(early_samples, 1e-2))))

clr_late_predict1 <- predict(svmfit_rep[[1]]$svmfit, data.frame(compositions::clr(impute(late_samples, 1e-2))))
clr_late_predict2 <- predict(svmfit_ICGC_rep[[1]]$svmfit, data.frame(compositions::clr(impute(late_samples, 1e-2))))
clr_late_predictboth <- predict(svmfit_both_rep[[1]]$svmfit, data.frame(compositions::clr(impute(late_samples, 1e-2))))

table(tcga=clr_early_predict1, icgc=clr_early_predict2)

table(tcga=clr_late_predict1, icgc=clr_late_predict2)

table(clr_early_predictboth)
table(clr_late_predictboth)

extra_expand <- 0.2
extra_expand_v2 <- 0.2
dendro_all_britroc <- give_dendrogram_from_imputation(impute_VALUE = 1e-2,
                                                      exposures = rbind(exposures_all[!grepl('^IM_', rownames(exposures_all)),],
                                                                        exposures_britroc), return_grob = T)
plot(dendro_all_britroc)

table(clr_britroc_predict1, clr_britroc_predict2)

ggplot(melt(list(TCGA=sapply(svmfit_rep, function(i) i$confusionMatrix$overall['Accuracy']),
     ICGC=sapply(svmfit_ICGC_rep, function(i) i$confusionMatrix$overall['Accuracy']))), aes(x=L1, y=value))+
  geom_boxplot()+geom_jitter()+theme_bw()

dim(svmfit$coefs)
dim(svmfit$SV)

plot(exposures$s3,
     exposures$s4)

ilr <- compositions::ilr(impute(exposures, 1e-2))
library(mclust)
fit_mvn <- mclust::mvn(modelName = "Ellipsoidal", data = ilr)

sigma_est <- fit_mvn$parameters$variance$Sigma

paletteLength <- 50
myColor <- colorRampPalette(c("yellow", "white", "blue"))(paletteLength)
myBreaks <- c(seq(min(sigma_est), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(sigma_est)/paletteLength, max(sigma_est), length.out=floor(paletteLength/2)))
rownames(sigma_est) <- colnames(sigma_est) <- paste0('ILR', 1:6)

pheatmap::pheatmap(sigma_est, cluster_rows = F, cluster_cols = F,
                   color=myColor, breaks=myBreaks)

simulated_ilr <- mvtnorm::rmvnorm(1000, mean = fit_mvn$parameters$mean, sigma = sigma_est)
simulated_exp <- compositions::clrInv(compositions::ilr2clr(simulated_ilr))
colnames(simulated_exp) <- colnames(exposures)

prcomp <- prcomp(rbind(exposures, simulated_exp))
ggplot(cbind.data.frame(prcomp$x[,c(1,2)], source= c(rep('real', nrow(exposures)),
                                              rep('sim', nrow(simulated_exp)))),
       aes(x=PC1, y=PC2, col=source, group=source))+geom_point(alpha=0.1)+geom_density_2d()


ggplot(cbind.data.frame(exposures), aes(x=s3, y=s4, col=s7))+geom_point()


## ------------------------------------------------------------------------------------- ##
## With imputation
library('robCompositions')
library('zCompositions')


all_exposures <- rbind(exposures_all[!grepl('^IM_', rownames(exposures_all)),],
                       exposures_britroc)
dim(all_exposures)
View(all_exposures)
robimpute <- robCompositions::imputeBDLs(x = all_exposures, dl=matrix(rep(0.05, 7), nrow=1),
                                         maxit=10,eps=0.1,R=10,method="subPLS") ## dosesn't work


multLNimpute <- zCompositions::multLN(X = all_exposures, dl = matrix(rep(0.05, 7), nrow=1), label = 0)
multReplimpute <- zCompositions::multRepl(X = all_exposures, dl = matrix(rep(0.05, 7), nrow=1), label = 0)
rownames(multLNimpute) <- rownames(multReplimpute) <- rownames(all_exposures)

min(multLNimpute)
min(all_exposures)
plot(unlist(multLNimpute), unlist(all_exposures))
plot(unlist(multReplimpute), unlist(all_exposures))
plot(unlist(multReplimpute), unlist(multLNimpute))
multReplimpute == multLNimpute

dendro_all_britroc_multLNimpute <- give_dendrogram_from_imputation(exposures = as(multLNimpute, 'matrix'),
                                                                   impute_VALUE = 0, return_grob = T, plot=T)

dendro_all_britroc_multReplimpute <- give_dendrogram_from_imputation(exposures = as(multReplimpute, 'matrix'),
                                                                   impute_VALUE = 0, return_grob = T, plot=T)

cowplot::plot_grid(dendro_all_britroc,
                   dendro_all_britroc_multLNimpute,
                   dendro_all_britroc_multReplimpute)

####### Barplotso of exposures, sorted
exposures_allv2 <- exposures_all
rownames(exposures_allv2) <- exposures_allv2$Var2
exposures_allv2 <- exposures_allv2[,-1]

createBarplot(matrix_exposures = as.matrix(exposures_allv2),
              order_labels = as.character(rownames(exposures_allv2))[order(exposures_allv2$s1)]) ## continuum

createBarplot(matrix_exposures = as.matrix(exposures_allv2),
              order_labels = as.character(rownames(exposures_allv2))[order(exposures_allv2$s2)],
              levels_signatures = paste0('s', 1:7)[c(c(1, 3:7, 2))]) ## sparser; white interesting

createBarplot(matrix_exposures = as.matrix(exposures_allv2),
              order_labels = as.character(rownames(exposures_allv2))[order(exposures_allv2$s3)],
              levels_signatures = paste0('s', 1:7)[c(c(1:2, 4:7, 3))]) ## sparser similar to s2

createBarplot(matrix_exposures = as.matrix(exposures_allv2),
              order_labels = as.character(rownames(exposures_allv2))[order(exposures_allv2$s4)],
              levels_signatures = paste0('s', 1:7)[c(c(1:3, 5:7, 4))]) ## sparser, similar to s2/s3

createBarplot(matrix_exposures = as.matrix(exposures_allv2),
              order_labels = as.character(rownames(exposures_allv2))[order(exposures_allv2$s5)],
              levels_signatures = paste0('s', 1:7)[c(c(1:4, 6:7, 5))]) ## almost everywhere

createBarplot(matrix_exposures = as.matrix(exposures_allv2),
              order_labels = as.character(rownames(exposures_allv2))[order(exposures_allv2$s6)],
              levels_signatures = paste0('s', 1:7)[c(c(1:5, 7, 6))]) ## the groups that don't have it are the ones that "look weird"

createBarplot(matrix_exposures = as.matrix(exposures_allv2),
              order_labels = as.character(rownames(exposures_allv2))[order(exposures_allv2$s7)],
              levels_signatures = paste0('s', 1:7)[c(c(1:7))]) ## continuum. weird

##--------------------------------------------------------

# > table(clr_early_predictboth)
# clr_early_predictboth
# FALSE  TRUE 
# 17     3 
# > table(clr_late_predictboth)
# clr_late_predictboth
# FALSE  TRUE 
# 32    40 

clr_britroc_arx <- compositions::clr(impute(exposures_britroc[as.character(patient.meta$SAMPLE_ID[patient.meta$group == "arx"]),], 1e-2))
clr_britroc_rlps <- compositions::clr(impute(exposures_britroc[as.character(patient.meta$SAMPLE_ID[patient.meta$group == "rlps"]),], 1e-2))
clr_britroc_arx_pred <- predict(svmfit_both_rep[[1]]$svmfit, data.frame(clr_britroc_arx))
clr_britroc_rlps_pred <- predict(svmfit_both_rep[[1]]$svmfit, data.frame(clr_britroc_rlps))
table(clr_britroc_arx_pred)
table(clr_britroc_arx_pred)/length(clr_britroc_arx_pred)

diploid_archivals <- names(clr_britroc_arx_pred)[as.character(clr_britroc_arx_pred) == "FALSE"]
length(diploid_archivals) ## 65 samples
## which of the diploid archivals have relapse?
diploid_archivals_patients <- as.character(patient.meta$PATIENT_ID[patient.meta$SAMPLE_ID %in% diploid_archivals])
length(diploid_archivals_patients) ## 65 patients with 
diploid_archivals_patients_with_rlps <- unique(as.character(patient.meta$PATIENT_ID[c(patient.meta$PATIENT_ID %in% diploid_archivals_patients) & (patient.meta$group == "rlps")]))
length((diploid_archivals_patients_with_rlps)) ## 31 patients with diploid archival have some relapse

diploid_archivals_patients_with_rlps_pred_arx <- sapply(diploid_archivals_patients_with_rlps, function(i) table(clr_britroc_arx_pred[patient.meta$SAMPLE_ID[(patient.meta$PATIENT_ID == i) & patient.meta$group == "arx"]]))
table(apply(diploid_archivals_patients_with_rlps_pred_arx, 2, function(i) sum(i>0)))
diploid_archivals_patients_with_rlps_pred_rlps <- sapply(diploid_archivals_patients_with_rlps, function(i) table(clr_britroc_rlps_pred[patient.meta$SAMPLE_ID[(patient.meta$PATIENT_ID == i) & patient.meta$group == "rlps"]]))

table(apply(diploid_archivals_patients_with_rlps_pred_rlps, 2, function(i) sum(i>0))) ## 6 samples with different ploidy in their relapses
different_ploidy_relapses <- unique(names(which(apply(diploid_archivals_patients_with_rlps_pred_rlps, 2, function(i) sum(i>0)) > 1)))

table(apply(diploid_archivals_patients_with_rlps_pred_rlps[,!(colnames(diploid_archivals_patients_with_rlps_pred_rlps) %in% different_ploidy_relapses)],
      2, function(i) which(i>0)))

## the patients who undergo WGD, according to the SVM
unique(names(which(apply(diploid_archivals_patients_with_rlps_pred_rlps[,!(colnames(diploid_archivals_patients_with_rlps_pred_rlps) %in% different_ploidy_relapses)],
      2, function(i) which(i>0)) == 2)))

# "BRITROC-242" "BRITROC-65"  "BRITROC-94" : these undergo WGD according to the SVM

#' 100% diploid -----> 0.15 WGD in early (3/20)
#'              |
#'              |-----> 0.85 diploid in early (17/20)
#'                                  |
#'                                  |----> 0.48 WGD in late archival (61/126)
#'                                  |
#'                                  |----> 0.52 diploid in late archival (65/126)
#'                                          |
#'                                          |-----> 31 patients with some relapse. No archival samples have contradictory ploidy
#'                                                    |
#'                                                    |----> 0.0645 (2 patients, 6 samples) with different ploidy in their relapses, ("BRITROC-274" "BRITROC-74")
#'                                                    |
#'                                                    |----> 0.838 (26 patients) remain diploid
#'                                                    |
#'                                                    |----> 0.0967 (3 patients) undergo WGD

#' "BRITROC-274": relapses of different ploidy (archival all the same). change according to UMAP
#' "BRITROC-74": relapses of different ploidy (archival all the same). change according to UMAP
#' "BRITROC-242": archival diploid, relapse WGD, according to SVM
#' "BRITROC-65": archival diploid, relapse WGD, according to SVM
#' "BRITROC-94": archival diploid, relapse WGD, according to SVM
#' 
#' BRITROC-23: according to UMAP
#' BRITROC-209: according to UMAP
#' BRITROC-216: : according to UMAP
#' BRITROC-241: according to UMAP
#' BRITROC-267: according to UMAP

give_samples_from_patient <- function(patient, group){
  patient.meta$SAMPLE_ID[(patient.meta$PATIENT_ID == patient) & (patient.meta$group == group)]
}

patient.meta[patient.meta$PATIENT_ID == "BRITROC-242",] ## no clear change in ploidy

c(unique(as.character(clr_britroc_arx_pred[give_samples_from_patient("BRITROC-242", "arx")])),
     unique(as.character(clr_britroc_rlps_pred[give_samples_from_patient("BRITROC-242", "rlps")])))

lapply(c("BRITROC-242", "BRITROC-65", "BRITROC-94"), function(pat){
  c(unique(as.character(clr_britroc_arx_pred[give_samples_from_patient(pat, "arx")])),
    unique(as.character(clr_britroc_rlps_pred[give_samples_from_patient(pat, "rlps")])))
}) ## samples that SVM finds tetraploid in relapse, but UMAP doesn't

lapply(c("BRITROC-23", "BRITROC-209", "BRITROC-216", "BRITROC-241", "BRITROC-267"), function(pat){
  c(unique(as.character(clr_britroc_arx_pred[give_samples_from_patient(pat, "arx")])),
    unique(as.character(clr_britroc_rlps_pred[give_samples_from_patient(pat, "rlps")])))
}) ## samples that SVM finds diploid in relapse, but UMAP detects a change

