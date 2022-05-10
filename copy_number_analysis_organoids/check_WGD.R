rm(list = ls())

setwd("~/Documents/PhD/other_repos/Vias_Brenton/figures/")

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

source("../copy_number_analysis_organoids/helper_functions.R")

exposures <- readRDS("../copy_number_analysis_organoids/robjects/exposures.RDS")
dendrograminputclr <- readRDS("../copy_number_analysis_organoids/robjects/dendrograminputclr.RDS")
heatmapinputclr <- readRDS("../copy_number_analysis_organoids/robjects/heatmapinputclr.RDS")

heatmapinputclr

wgd_hasse <- read.table("../../../CDA_in_Cancer/data/TCGA_WGD/Haase2019_TCGA.giScores.wgd.txt", header = T)


heatmap_data <- heatmapinputclr$data
heatmap_data = heatmap_data[heatmap_data$Var2 %in% wgd_hasse$patient,]
head(heatmap_data)

heatmap_data$WGD = wgd_hasse$wgd[match(heatmap_data$Var2, wgd_hasse$patient)]

ggplot(heatmap_data, aes(x=Var2, fill=Var1, y=value))+geom_bar(stat = "identity")+
  scale_fill_brewer(palette="Dark2")+facet_wrap(.~WGD, nrow=2)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.title = element_blank())+
  ggtitle('WGD')
ggsave("../copy_number_analysis_organoids/figures/WGD_exposures.pdf", height = 4, width = 7)


## SVM
exposures <- dcast(heatmap_data, formula =  Var2 ~ Var1, value.var = 'value')
rownames(exposures) <- exposures[,1]; exposures <- exposures[,-1]
clr <- compositions::clr(impute(exposures, 1e-2))

library(e1071)
sample_WGD <- data.frame(heatmap_data) %>% dplyr::filter(Var1 == 's1') %>% dplyr::select(WGD) %>% unlist()
training_set_idx <- sample(1:nrow(clr), size = floor(.7*nrow(clr)))
training_set = clr[training_set_idx,]
validation_set = clr[!(1:nrow(clr) %in% training_set_idx),]
training_dat = data.frame(x=training_set, y = as.factor(sample_WGD[training_set_idx]))
validation_dat = data.frame(x=validation_set, y = as.factor(sample_WGD[!(1:nrow(clr) %in% training_set_idx)]))
svmfit = svm(y ~ ., data = training_dat, kernel = "linear", cost = 10, scale = FALSE)
print(svmfit)

validation_predicted = predict(svmfit, validation_dat)

confusionMatrix(validation_predicted, validation_dat$y)

plot(x = svmfit, data = training_dat, formula = x.s1 ~ x.s4, fill = TRUE, grid = 50, slice = list())
pdf("../copy_number_analysis_organoids/figures/WGD_exposures_svm1.pdf", height = 4, width = 4)
plot(x = svmfit, data = training_dat, formula = x.s4 ~ x.s3, fill = TRUE, grid = 50, slice = list())
dev.off()
pdf("../copy_number_analysis_organoids/figures/WGD_exposures_svm2.pdf", height = 4, width = 4)
plot(x = svmfit, data = training_dat, formula = x.s3 ~ x.s5, fill = TRUE, grid = 50, slice = list())
dev.off()
plot(x = svmfit, data = training_dat, formula = x.s5 ~ x.s6, fill = TRUE, grid = 50, slice = list())

pheatmap::pheatmap(svmfit$SV)
pheatmap::pheatmap(svmfit$coefs, cluster_cols = F, cluster_rows = F)

add_rownames <- function(i){
  rownames(i) <- paste0('Dim', 1:nrow(i))
  i
}

pheatmap::pheatmap(svmfit$SV, annotation_row = add_rownames(svmfit$coefs))


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


