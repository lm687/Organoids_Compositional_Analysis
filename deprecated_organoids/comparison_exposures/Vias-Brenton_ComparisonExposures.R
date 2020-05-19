#<begin_omit>```{r}
rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
Sys.setenv(LANG='en')
#<end_omit>```{r}

#<begin_chunk>```{r, libraries,message=FALSE}
library(ggrepel)
#library(cowplot)
library(compositions)
#library(CompSign)
source("../../../../CDA_in_Cancer/code/functions/meretricious/pretty_plots/prettySignatures.R")
#<end_chunk>

#<begin_text>
## Comparison of the old exposures (from Nat Gen) and the new ones defined by Ruben, where the ASCAT segmentation has been changed, and normal segments have been removed.
#<end_text>

#<begin_chunk>```{r, natgen_data,include=TRUE}

load("../../../../CDA_in_Cancer/data/Robj/image_NatGen_rmd.RData")
natgen <- as.matrix(sig_data_unorm[,1:7])
## Normalisation is not done in such a way that rows add up to 1. Re-normalising
natgen_renormalised <- sweep(natgen, 1, rowSums(natgen), '/')

new_exposures <- readRDS("../data/3_tcga_exposuresYAPSA.rds")
new_exposures <- readRDS("../data/Export-matrix_OV_Sigs_on_TCGA-OV_12112019.rds")
rownames(new_exposures) <- sapply(rownames(new_exposures), function(u) paste0(strsplit(u, '-')[[1]][1:3], collapse = '-'))
#<end_chunk>


#<begin_chunk>```{r, match,include=TRUE}
natgen_renormalised <- natgen_renormalised[grep('TCGA', rownames(natgen_renormalised)),]
rownames(natgen_renormalised) <- sapply(rownames(natgen_renormalised), function(i) paste0(strsplit(i, '-')[[1]][1:3], collapse = '-'))
mtch <- match(rownames(new_exposures), rownames(natgen_renormalised))
melt_comparison <- (melt(list(new_exposures, natgen_renormalised[mtch,])))
melt_comparison <- melt_comparison[!is.na(melt_comparison$Var1),]
melt_comparison[which(duplicated(melt_comparison[,c(1,2,4)])),]
melt_comparison <- dcast(melt_comparison, Var1+Var2~L1, value.var = 'value')
colnames(melt_comparison)[3:4] <- c('NewExposures', 'OldExposures')
#<end_chunk>

#<begin_chunk>```{r, scatterplot,include=TRUE}
ggplot(melt_comparison, aes(x=NewExposures, y=OldExposures))+geom_point()+facet_wrap(.~Var2, scales = "free", nrow=2)+
  geom_abline(intercept = 0, slope = 1)
#<end_chunk>


