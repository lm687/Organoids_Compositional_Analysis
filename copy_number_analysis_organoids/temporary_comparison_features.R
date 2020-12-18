
setwd("/Users/morril01/Documents/PhD/other_repos/b_tape/Vias_Brenton/copy_number_analysis_organoids/")
rm(list = ls())

require(dplyr)
require(gridExtra)
source("../../../../other_repos/BriTROC-cnsignatures-bfb69cd72c50/main_functions.R")
source("helper_functions.R")

pcawg_CN_features = readRDS("data/pcawg_CN_features.rds")
organoids_absolute_copynumber = readRDS("data/organoid_absolute_CN.rds")
organoids_CN_features = extractCopynumberFeatures(organoids_absolute_copynumber)
tcga_CN_features = readRDS("data/tcga_CN_features.rds")
BriTROC_CN_features = readRDS("data/BriTROC_CN_features.rds")

#------------------------------------------------------------------------------------------------#

create_distrib_df = function(feature, name_value_col){
  df = list()
  df[['pcawg']] = pcawg_CN_features[[feature]] %>% group_by(ID) %>% mutate(ct1_num = as.numeric(get(name_value_col))) %>% summarize(mean(ct1_num))
  df[['org']] = organoids_CN_features[[feature]] %>% group_by(ID) %>% mutate(ct1_num = as.numeric(get(name_value_col))) %>% summarize(mean(ct1_num))
  df[['tcga']] = tcga_CN_features[[feature]] %>% group_by(ID) %>% mutate(ct1_num = as.numeric(get(name_value_col))) %>% summarize(mean(ct1_num))
  df[['BriTROC']] = BriTROC_CN_features[[feature]] %>% group_by(ID) %>% mutate(ct1_num = as.numeric(get(name_value_col))) %>% summarize(mean(ct1_num))
  return(df)
}

remove_infty = function(i){
  i = i[!is.infinite(i)]
  return(i)
}

#------------------------------------------------------------------------------------------------#

distrib_segsize = create_distrib_df('segsize', 'value')
distrib_bp10MB = create_distrib_df('bp10MB', 'value')
distrib_osCN = create_distrib_df('osCN', 'value')
distrib_bpchrarm = create_distrib_df('bpchrarm', 'ct1')
distrib_changepoint = create_distrib_df('changepoint', 'value')
distrib_copynumber = create_distrib_df('copynumber', 'value')

#------------------------------------------------------------------------------------------------#

## 1/6 breakpoints per 10MB
give_joint_histogram(list(distrib_bp10MB[['org']]$`mean(ct1_num)`,
                          distrib_bp10MB[['pcawg']]$`mean(ct1_num)`,
                          distrib_bp10MB[['tcga']]$`mean(ct1_num)`,
                          distrib_bp10MB[['BriTROC']]$`mean(ct1_num)`), no_colour=FALSE)

give_joint_histogram(list(distrib_bp10MB[['org']]$`mean(ct1_num)` %>% log,
                          distrib_bp10MB[['pcawg']]$`mean(ct1_num)` %>% log,
                          distrib_bp10MB[['tcga']]$`mean(ct1_num)` %>% log,
                          distrib_bp10MB[['BriTROC']]$`mean(ct1_num)` %>% log), no_colour=FALSE)

t.test(distrib_bp10MB[['org']]$`mean(ct1_num)` %>% log,
       c(distrib_bp10MB[['pcawg']]$`mean(ct1_num)` %>% log,
         distrib_bp10MB[['tcga']]$`mean(ct1_num)` %>% log,
         distrib_bp10MB[['BriTROC']]$`mean(ct1_num)` %>% log))

#------------------------------------------------------------------------------------------------#
## 2/6 segment size
give_joint_histogram(list(distrib_segsize[['org']]$`mean(ct1_num)`,
                          distrib_segsize[['pcawg']]$`mean(ct1_num)`,
                          distrib_segsize[['tcga']]$`mean(ct1_num)`,
                          distrib_segsize[['BriTROC']]$`mean(ct1_num)`), no_colour=FALSE)

give_joint_histogram(list(distrib_segsize[['org']]$`mean(ct1_num)` %>% log,
                          distrib_segsize[['pcawg']]$`mean(ct1_num)` %>% log,
                          distrib_segsize[['tcga']]$`mean(ct1_num)` %>% log,
                          distrib_segsize[['BriTROC']]$`mean(ct1_num)` %>% log), no_colour=FALSE)

t.test(distrib_segsize[['org']]$`mean(ct1_num)` %>% log,
       c(distrib_segsize[['pcawg']]$`mean(ct1_num)` %>% log,
         distrib_segsize[['tcga']]$`mean(ct1_num)` %>% log,
         distrib_segsize[['BriTROC']]$`mean(ct1_num)` %>% log))

#------------------------------------------------------------------------------------------------#
## 3/6 oscillating CN
grid.arrange(give_joint_histogram(list(distrib_osCN[['org']]$`mean(ct1_num)`,
                          distrib_osCN[['pcawg']]$`mean(ct1_num)`,
                          distrib_osCN[['tcga']]$`mean(ct1_num)`,
                          distrib_osCN[['BriTROC']]$`mean(ct1_num)`), no_colour=FALSE)+ggtitle('Raw'),
give_joint_histogram(list(distrib_osCN[['org']]$`mean(ct1_num)` %>% log,
                          distrib_osCN[['pcawg']]$`mean(ct1_num)` %>% log,
                          distrib_osCN[['tcga']]$`mean(ct1_num)` %>% log,
                          distrib_osCN[['BriTROC']]$`mean(ct1_num)` %>% log), no_colour=FALSE)+ggtitle('Log transform'), ncol=2)

t.test(distrib_osCN[['org']]$`mean(ct1_num)` %>% log,
       remove_infty(c(distrib_osCN[['pcawg']]$`mean(ct1_num)` %>% log,
         distrib_osCN[['tcga']]$`mean(ct1_num)` %>% log,
         distrib_osCN[['BriTROC']]$`mean(ct1_num)` %>% log)))

#------------------------------------------------------------------------------------------------#
## 4/6 num of breakpoints per chromosome arm
give_joint_histogram(list(distrib_bpchrarm[['org']]$`mean(ct1_num)`,
                          distrib_bpchrarm[['pcawg']]$`mean(ct1_num)`,
                          distrib_bpchrarm[['tcga']]$`mean(ct1_num)`,
                          distrib_bpchrarm[['BriTROC']]$`mean(ct1_num)`), no_colour=FALSE)

give_joint_histogram(list(distrib_bpchrarm[['org']]$`mean(ct1_num)` %>% log,
                          distrib_bpchrarm[['pcawg']]$`mean(ct1_num)` %>% log,
                          distrib_bpchrarm[['tcga']]$`mean(ct1_num)` %>% log,
                          distrib_bpchrarm[['BriTROC']]$`mean(ct1_num)` %>% log), no_colour=FALSE)

## Same distribution of num of breakpoints per chromosome arm
t.test(distrib_bpchrarm[['org']]$`mean(ct1_num)` %>% log,
       c(distrib_bpchrarm[['pcawg']]$`mean(ct1_num)` %>% log,
         distrib_bpchrarm[['tcga']]$`mean(ct1_num)` %>% log,
         distrib_bpchrarm[['BriTROC']]$`mean(ct1_num)` %>% log))


#------------------------------------------------------------------------------------------------#
## 5/6 num of changepoints
give_joint_histogram(list(distrib_changepoint[['org']]$`mean(ct1_num)`,
                          distrib_changepoint[['pcawg']]$`mean(ct1_num)`,
                          distrib_changepoint[['tcga']]$`mean(ct1_num)`,
                          distrib_changepoint[['BriTROC']]$`mean(ct1_num)`), no_colour=FALSE)

give_joint_histogram(list(distrib_changepoint[['org']]$`mean(ct1_num)` %>% log,
                          distrib_changepoint[['pcawg']]$`mean(ct1_num)` %>% log,
                          distrib_changepoint[['tcga']]$`mean(ct1_num)` %>% log,
                          distrib_changepoint[['BriTROC']]$`mean(ct1_num)` %>% log), no_colour=FALSE)

## Same distribution of num of breakpoints per chromosome arm
t.test(distrib_changepoint[['org']]$`mean(ct1_num)` %>% log,
       c(distrib_changepoint[['pcawg']]$`mean(ct1_num)` %>% log,
         distrib_changepoint[['tcga']]$`mean(ct1_num)` %>% log,
         distrib_changepoint[['BriTROC']]$`mean(ct1_num)` %>% log))

#------------------------------------------------------------------------------------------------#
## 6/6 copy number
give_joint_histogram(list(distrib_copynumber[['org']]$`mean(ct1_num)`,
                          distrib_copynumber[['pcawg']]$`mean(ct1_num)`,
                          distrib_copynumber[['tcga']]$`mean(ct1_num)`,
                          distrib_copynumber[['BriTROC']]$`mean(ct1_num)`), no_colour=FALSE)

give_joint_histogram(list(distrib_copynumber[['org']]$`mean(ct1_num)` %>% log,
                          distrib_copynumber[['pcawg']]$`mean(ct1_num)` %>% log,
                          distrib_copynumber[['tcga']]$`mean(ct1_num)` %>% log,
                          distrib_copynumber[['BriTROC']]$`mean(ct1_num)` %>% log), no_colour=FALSE)

## Same distribution of num of breakpoints per chromosome arm
t.test(distrib_copynumber[['org']]$`mean(ct1_num)` %>% log,
       c(distrib_copynumber[['pcawg']]$`mean(ct1_num)` %>% log,
         distrib_copynumber[['tcga']]$`mean(ct1_num)` %>% log,
         distrib_copynumber[['BriTROC']]$`mean(ct1_num)` %>% log))

