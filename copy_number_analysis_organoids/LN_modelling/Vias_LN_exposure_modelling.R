## Infence using stan_fit_LNM for data in the simplex
## Simple LN, no fixed or random effects

rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

set.seed(1234)

library(rstan)
library(uuid)
library(ggplot2)
library(optparse)
library(reshape2)
library(pheatmap)
library(bayesplot)
source("../../../../../GlobalDA/code/2_inference/helper/helper_DA_stan.R")
source("../../../../../GlobalDA/code/3_analysis/helper/helper_analyse_posteriors.R")
source("../../../../../GlobalDA/code/2_inference_TMB/helper_TMB.R")
Nits = 2500

model_file_name = "../../../../../ProjectOvarianMultisampleTree/code/files_analysis/modelling_LN/stan_fit_simple_LN.stan"
rstan_options(auto_write = TRUE)
stanc(model_file_name)

#-------------------------------------------------------------------------------------------#
give_short_names = function(i){
  if(grepl('TCGA', i)){
    paste0(strsplit(i, split = '-')[[1]][1:3], collapse='-')
  }else{
    i
  }
}

org <- as(readRDS("../data/organoid_exposures_Aug21.rds"), 'matrix')
rownames(org) <- paste0('Sample ', 1:nrow(org))
names_orgs = readxl::read_xlsx("data/NewOrganoidNaming.xlsx")

natgen = natgen_metadata = list()
sig_data = readRDS("../data/sig_data_unorm.RDS")
sig_data = cbind(sweep(sig_data[,1:7], 1, rowSums(sig_data[,1:7]), '/'),
                 sig_data[,8:ncol(sig_data)])
natgen[[1]] <- as.matrix(sig_data[,1:7])
natgen_metadata[[1]] <- sig_data[,8:10]
natgen[[2]] <- readRDS("../data/Export-matrix_OV_Sigs_on_TCGA-OV_12112019.rds")
natgen_metadata[[2]] <- cbind.data.frame(study=rep('Updated TCGA', nrow(natgen[[2]])), age=NA, age.cat=NA, stringsAsFactors = FALSE, row.names=rownames(natgen[[2]]))
names(natgen_metadata) = names(natgen) = c('ExposuresNatGen', 'UpdatedExposures')

#------------ Only keep TCGA samples which are of good enough quality------------#
summary_ascat = read.table("../data/summary.ascatTCGA.penalty70.txt", header = TRUE, stringsAsFactors = FALSE)
good_tcga = summary_ascat$name[summary_ascat$dCIN]
good_tcga = good_tcga[!is.na(good_tcga)]
bool_tcga= lapply(natgen, function(i) grepl('TCGA', rownames(i)))
rm_na = function(df) !apply(df, 1, function(rw) all(is.na(rw)))

## modify the dataframes
for(version in 1:2){
  rm_bad_samples = !(sapply(rownames(natgen[[version]])[bool_tcga[[version]]], give_short_names) %in% good_tcga)
  natgen[[version]][bool_tcga[[version]],][ rm_bad_samples,] <- NA
  natgen[[version]] = natgen[[version]][rm_na(natgen[[version]]),]
  natgen_metadata[[version]][bool_tcga[[version]],][ rm_bad_samples,] <- NA
  natgen_metadata[[version]] = natgen_metadata[[version]][rm_na(natgen_metadata[[version]]),]
}
rownames(natgen[[1]]) = as.character(sapply(rownames(natgen[[1]]), give_short_names))

# add the non-TCGA samples to natgen2
bool_tcga = lapply(natgen, function(i) grepl('TCGA', rownames(i))) ## re-compute
natgen[[2]] = rbind(natgen[[2]], natgen[[1]][!bool_tcga[[1]],])

exposures1 = rbind(org, natgen[[1]])
exposures2 = rbind(org, natgen[[2]])

exposures1_no_zero = exposures1
exposures1_no_zero[which(exposures1 == 0)] = 1e-4
exposures1_no_zero = normalise_rw(exposures1_no_zero)

exposures2_no_zero = exposures2
exposures2_no_zero[which(exposures2 == 0)] = 1e-4
exposures2_no_zero = normalise_rw(exposures2_no_zero)

exposures = exposures1_no_zero

#-------------------------------------------------------------------------------------------#

stan_data = list(n=nrow(exposures1_no_zero),
                 d = ncol(exposures1_no_zero),
                 W = exposures1_no_zero)

params = c('mu', 'Sigma')

fit_stan <- stan(file = model_file_name, data = stan_data,
                 iter = Nits, chains = 4, cores = 4, thin = 1, pars = params,
                 control = list(stepsize=3, adapt_delta=.95, max_treedepth=15))

# pairs(fit_stan)

## Good convergence
max(bayesplot::rhat(fit_stan))
posterior2 = as.matrix(fit_stan)

# saveRDS(posterior2, file = paste0("../../../out/robj/inference/simple_", sample_name, "_posteriors.RDS"))

names_betas = colnames(posterior2)[grep("beta", colnames(posterior2))]
p <- bayesplot::mcmc_trace(posterior2,  pars = names_betas,
                           facet_args = list(nrow = 1, labeller = label_parsed))+ theme(text = element_text(size=10))
p + facet_text(size = 6)

## looking at the correlation between inferred parameters (in this case, just betas)
## there is a clear correlations between the intersect and the corresponding slope, for each of the (sub)features index by 1, ..., d-1
pdf(paste0("../figures/LM_modelling/exposures1_pairs_plot.pdf"))
pairs(fit_stan)
dev.off()
# ## beta coefficients for intersect
# pairs(fit_stan, pars = names(fit_stan)[grep('beta\\[1,', names(fit_stan))], text.panel = "Coefficients for the intercept")
# ## beta coefficients for slope
# pairs(fit_stan, pars = names(fit_stan)[grep('beta\\[2,', names(fit_stan))], text.panel = "Coefficients for the slope")

## Reminder that beta[1,x] corresponds to the intercept, and beta[2,x] to the slopes
## plot all but LP
pdf(paste0("../figures/LM_modelling/exposures1_pars.pdf"))
bayesplot::mcmc_areas(posterior2, pars = colnames(posterior2)[-ncol(posterior2)])+ggtitle('Slope')
dev.off()

pdf(paste0("../figures/LM_modelling/exposures1_parcoord.pdf"), width = 12)
bayesplot::mcmc_parcoord(posterior2)+ggtitle('Slope')
dev.off()

pdf(paste0("../figures/LM_modelling/exposures1_traceplot.pdf"), width = 18, height = 8)
traceplot(fit_stan, pars = c("mu", "Sigma"), inc_warmup = TRUE, nrow = 2)
dev.off()
# rstanarm::posterior_vs_prior(fit_stan)
# 
## let's simulate under the model with the inferred parameters and compare it to the actual data
size_sim = 400
stopifnot(size_sim %% 2 == 0)
idx_posteriors = sample(1:nrow(posterior2), size = size_sim, replace = FALSE)

simulate_from_model = function(idx_posterior){
  posterior_row = posterior2[idx_posterior,]
  mu = python_like_select_name(posterior_row, 'mu')
  multiv = t(apply(t(matrix(mu)), 1, mvtnorm::rmvnorm, n = 1,
                   sigma =  matrix(python_like_select_name(posterior_row, 'Sigma'),
                                   ncol = ncol(exposures)-1)))
  return(softmax(c(multiv, 0)))
}

sim_results = t(sapply(idx_posteriors, simulate_from_model))
colnames(sim_results) = colnames(exposures)

pdf(paste0("../figures/LM_modelling/exposures1_comparison_with_simulation.pdf"))
pairs(rbind(sim_results, exposures), col=c('#080a4d', 'red')[rep(c(1,2), c(nrow(sim_results), nrow(exposures)))], pch=19, cex=0.1)
dev.off()

## examine betas
hist(posterior2[,'beta[2,1]'])
hist(posterior2[,'beta[2,2]'])
hist(posterior2[,'beta[2,3]'])
hist(posterior2[,'beta[2,4]'])

colMeans(python_like_select_colnames(posterior2, 'beta\\[2'))

## credible intervals
