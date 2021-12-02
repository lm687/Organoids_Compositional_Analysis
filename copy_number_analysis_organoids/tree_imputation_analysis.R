rm(list = ls())
setwd("~/Documents/PhD/other_repos/Vias_Brenton/copy_number_analysis_organoids/")
source("helper_functions.R")
library(CompSign)
library(ggplot2)
library(ggrepel)
library(ggdendro)
library(reshape2)
library(gridExtra)
library(TreeDist)
library(jcolors)
library(pheatmap)
library(cowplot)

## read data
all_natgen0 <- readRDS("robjects/allnatgen_UpdatedExposures.RDS")
all_natgen <- list()
which_natgen = 'UpdatedExposures'
all_natgen[[which_natgen]] <- all_natgen0

## parameters for plotting
extra_expand <- 0.2
extra_expand_v2 <- 0.2

###########################
## original
maximput <- 10^-1
lenghtout <- 30
## with more values
# maximput <- 10^-0
# lenghtout <- 30
###########################
maximput_char <- gsub("[.]", "", as.character(maximput))

imputation_values_0 <- seq(-3, log10(maximput), length.out = lenghtout)
imputation_values <- sapply(imputation_values_0, function(i) 10^{i})
imputation_values_large <- sapply(seq(-3, log10(maximput), length.out = 100), function(i) 10^{i})
imputation_values

dendrogams_imputation <- lapply(imputation_values, give_distance_from_imputation)
dendrogams_imputation_large <- lapply(imputation_values_large, give_distance_from_imputation)

## metrics; mean distance between organoids vs mean distance from organoids to primary
## same for median
## same for max?
## some other metric of clustering?

comparison_orgs_TCGA <- lapply(dendrogams_imputation, function(i){
  .dists <- as(i, 'matrix')
  dists_within_orgs <- .dists[grepl('PDO', rownames(.dists)),grepl('PDO', colnames(.dists))]
  diag(dists_within_orgs) <- NA
  dists_with_primary <- .dists[!grepl('PDO', rownames(.dists)),grepl('PDO', colnames(.dists))]
  diag(dists_with_primary) <- NA
  
  first_split <- cutree(hclust(i), k = 2)
  
  ## find the size of the first clade without any orgs
  split_without_orgs <- first_split; .k <- 2
  while(length(table(split_without_orgs[grepl('PDO', names(split_without_orgs))])) == .k){
    .k <- .k+1
    split_without_orgs <- cutree(hclust(i), k = .k)
  }
  ## find the clade without organoids
  first_empty_clade_size <- table(split_without_orgs)[names(table(split_without_orgs))[!(names(table(split_without_orgs)) %in% names(table(split_without_orgs[grepl('PDO', names(split_without_orgs))])))]]
  
  
  cbind.data.frame(mean=c(mean(dists_within_orgs, na.rm = T),
    mean(dists_with_primary, na.rm = T)),
    median=c(median(dists_within_orgs, na.rm = T),
      median(dists_with_primary, na.rm = T)),
    max=c(max(dists_within_orgs, na.rm = T),
      max(dists_with_primary, na.rm = T)),
    min=c(min(dists_within_orgs, na.rm = T),
      min(dists_with_primary, na.rm = T)),
    largest_clade_size=max(table(first_split)),
    num_orgs_in_largest_clade=max(table(first_split[grepl('PDO', names(first_split))])),
    size_first_empty_clade=first_empty_clade_size)
})

summary_orgs_in_dendro <- cbind.data.frame(ratio_mean_dist=sapply(comparison_orgs_TCGA, function(i) i['mean'][1,]/i['mean'][2,]),
                                ratio_median_dist=sapply(comparison_orgs_TCGA, function(i) i['median'][1,]/i['median'][2,]),
                                ratio_max_dist=sapply(comparison_orgs_TCGA, function(i) i['max'][1,]/i['max'][2,]),
                                ratio_min_dist=sapply(comparison_orgs_TCGA, function(i) i['min'][1,]/i['min'][2,]),
                                largest_clade_size=sapply(comparison_orgs_TCGA, function(i) i['largest_clade_size'][1,]),
                                num_orgs_in_largest_clade=sapply(comparison_orgs_TCGA, function(i) i['num_orgs_in_largest_clade'][1,]),
                                first_empty_clade_size=sapply(comparison_orgs_TCGA, function(i) i['size_first_empty_clade'][1,]))
                                
summary_orgs_in_dendro$imputation_values = imputation_values
summary_orgs_in_dendro_melt <- melt(summary_orgs_in_dendro, id.vars='imputation_values')
summary_orgs_in_dendro_melt$facet <- paste0((summary_orgs_in_dendro_melt$variable %in% c('largest_clade_size')), 
  (summary_orgs_in_dendro_melt$variable %in% c('num_orgs_in_largest_clade')),
  (summary_orgs_in_dendro_melt$variable %in% c('first_empty_clade_size')))

ggplot(summary_orgs_in_dendro_melt, aes(x=imputation_values, y = value,
                                   col=variable, group=variable))+
  geom_line()+geom_point()+scale_x_continuous(trans = "log10")+theme_bw()+
  # facet_wrap(.~facet, scales = "free_y")+
  facet_wrap(.~variable, scales = "free_y", nrow=1)+
  theme(legend.position = "bottom")
ggsave(paste0("figures/metrics_dendrogram_imputation_maximput_char", maximput_char, ".pdf"), height = 3, width = 10)

## now, compute the similarity of trees, as imputation grows
trees <- lapply(dendrogams_imputation, function(i){
  ape::as.phylo.hclust(hclust(i))
})
trees_2 <- lapply(dendrogams_imputation, function(i){
  (hclust(i))
})
trees_3 <- lapply(dendrogams_imputation_large, function(i){ ## simpler version
  (hclust(i))
})

NyeSimilarity(trees[[1]], trees[[2]])

nyesimilarities <- matrix(NA, length(trees), length(trees))
for(i in 1:length(trees)){
  for(j in 1:i){
    nyesimilarities[i,j] = nyesimilarities[j,i] <- NyeSimilarity(trees[[i]], trees[[j]])
  }
}

##' compute the fraction of samples which are found in the same clade from one
##' imputation value to the next (i.e. we need to compare trees, like with the Nye similarities)
fraction_samples_in_same_major_clade <- matrix(NA, length(trees_3), length(trees_3))
for(i in 1:length(trees_3)){
  for(j in 1:i){
    .cutreei <- cutree(trees_3[[i]], k = 2)
    .cutreej <- cutree(trees_3[[j]], k = 2)
    fraction_samples_in_same_major_clade[i,j] = max(c(mean(.cutreei == .cutreej), mean(!.cutreei == .cutreej)))
  }
}

relative_size_clusters <- sapply(trees_3, function(i){
  .cutreei <- cutree(i, k = 2)
  max(table(.cutreei))/length(.cutreei)
})

plot(trees_3[[which.max(relative_size_clusters)]])

## smart way: not including the outgroup that we sometimes get
fraction_samples_in_same_major_clade_smart <- matrix(NA, length(trees_3), length(trees_3))
for(i in 1:length(trees_3)){
  for(j in 1:i){
    cat(i, '/', length(trees_3), '\t', j, '/', i, '\n')
    .cutreei <- cutree(trees_3[[i]], k = 2)
    .cutreej <- cutree(trees_3[[j]], k = 2)
    
    if((min(table(.cutreei))/length(.cutreei) < 0.1) | (min(table(.cutreej))/length(.cutreej) < 0.1)){
      .tab <- table(.cutreej)
      if(min(table(.cutreei))/length(.cutreei) < 0.1){
        samples_remove <- names(which(.cutreei == names(which.min(.tab))))
      }
      if(min(table(.cutreej))/length(.cutreej) < 0.1){
        samples_remove <- names(which(.cutreej == names(which.min(.tab))))
      }
      .i <- as(dendrogams_imputation_large[[i]],  'matrix')
      .j <- as(dendrogams_imputation_large[[j]],  'matrix')

      new_tree_i <- hclust(dist(.i[!(rownames(.i) %in% samples_remove),!(colnames(.i) %in% samples_remove)]))
      new_tree_j <- hclust(dist(.j[!(rownames(.j) %in% samples_remove),!(colnames(.j) %in% samples_remove)]))
      .cutreei <- cutree(new_tree_i, k = 2)
      .cutreej <- cutree(new_tree_j, k = 2)
    }
    
    fraction_samples_in_same_major_clade_smart[i,j] = max(c(mean(.cutreei == .cutreej), mean(!.cutreei == .cutreej)))
  }
}

colnames(nyesimilarities) <- rownames(nyesimilarities) <- imputation_values
colnames(fraction_samples_in_same_major_clade) <- rownames(fraction_samples_in_same_major_clade) <- imputation_values_large
colnames(fraction_samples_in_same_major_clade_smart) <- rownames(fraction_samples_in_same_major_clade_smart) <- imputation_values_large

image(nyesimilarities)

ggplot(melt(nyesimilarities), aes(x=Var1, y=Var2, fill=value))+geom_tile()+
  scale_x_continuous(trans = "log10")+scale_y_continuous(trans = "log10")+
  theme_bw()+labs(x='Imputation value', y='Imputation value')+
  guides(fill=guide_legend(title = 'NyeSimilarity'))+
  scale_fill_jcolors_contin("pal3", reverse = TRUE, bias = 2.25) 
ggsave(paste0("figures/nyesimilarities_dendrogram_imputation_maximput_char", maximput_char, ".pdf"),
       height = 4.5, width = 5.5)

fraction_samples_in_same_major_clade
ggplot(melt(fraction_samples_in_same_major_clade), aes(x=Var1, y=Var2, fill=value))+geom_tile()+
  scale_x_continuous(trans = "log10")+scale_y_continuous(trans = "log10")+
  theme_bw()+labs(x='Imputation value', y='Imputation value')+
  # guides(fill=guide_legend(title = 'NyeSimilarity'))+
  scale_fill_jcolors_contin("pal3", reverse = TRUE, bias = 2.25) 

ggplot(melt(fraction_samples_in_same_major_clade_smart), aes(x=Var1, y=Var2, fill=value))+geom_tile()+
  scale_x_continuous(trans = "log10")+scale_y_continuous(trans = "log10")+
  theme_bw()+labs(x='Imputation value', y='Imputation value')+
  # guides(fill=guide_legend(title = 'NyeSimilarity'))+
  scale_fill_jcolors_contin("pal3", reverse = TRUE, bias = 2.25) 

ggplot(data.frame(fraction_in_same_clade=rev(fraction_samples_in_same_major_clade[nrow(fraction_samples_in_same_major_clade),]),
                  imput=imputation_values_large,
                  size_largest=relative_size_clusters),
       aes(x=imput, y=fraction_in_same_clade))+geom_line()+
  geom_point(aes(size=size_largest, col=size_largest))+
  scale_x_continuous(trans = "log10")+theme_bw()



ggplot(data.frame(fraction_in_same_clade=sapply(1:(nrow(fraction_samples_in_same_major_clade)-1),
                                                function(i) fraction_samples_in_same_major_clade[i+1,i]),
                  imput=imputation_values_large[-1],
                  size_largest_change=relative_size_clusters[-length(relative_size_clusters)]-relative_size_clusters[-1]),
       aes(x=imput, y=fraction_in_same_clade))+geom_line()+
  # geom_point(aes(size=size_largest_change, col=size_largest))+
  scale_x_continuous(trans = "log10")+theme_bw()

ggplot(data.frame(fraction_in_same_clade=sapply(1:(nrow(fraction_samples_in_same_major_clade_smart)-1),
                                                function(i) fraction_samples_in_same_major_clade_smart[i+1,i]),
                  imput=imputation_values_large[-1],
                  size_largest_change=relative_size_clusters[-length(relative_size_clusters)]-relative_size_clusters[-1]),
       aes(x=imput, y=fraction_in_same_clade))+geom_line()+
  # geom_point(aes(size=size_largest_change, col=size_largest))+
  scale_x_continuous(trans = "log10")+theme_bw()

ggplot(data.frame(fraction_in_same_clade=rev(fraction_samples_in_same_major_clade[nrow(fraction_samples_in_same_major_clade),]),
                  size_largest=relative_size_clusters),
       aes(x=size_largest, y=fraction_in_same_clade))+
  geom_point()
ggplot(data.frame(fraction_in_same_clade=rev(fraction_samples_in_same_major_clade_smart[nrow(fraction_samples_in_same_major_clade_smart),]),
                  size_largest=relative_size_clusters),
       aes(x=size_largest, y=fraction_in_same_clade))+
  geom_point()+theme_bw()

k_large <- 30
## find the big clades that don't have any organoids, for each imputation value
pdf(paste0("figures/exposures_from_empty_clades_per_impoutation_value_maximput_char", maximput_char, ".pdf"))
samples_in_clades_without_orgs <- lapply(1:length(trees_2), function(i_idx){
  i <- trees_2[[i_idx]]
  .cutree_large <- cutree(i, k = k_large)
  clades_without_orgs <- (1:k_large)[!((1:k_large) %in% .cutree_large[grepl('PDO', names(.cutree_large))])]
  samples_in_clades <- all_natgen[[which_natgen]][names(.cutree_large[.cutree_large %in% clades_without_orgs]),]
  print(createBarplot(samples_in_clades)+
    ggtitle(paste0('Imputation:', imputation_values[i_idx], '\nNumber of clades with no orgs: ', 
                  length(clades_without_orgs))))
  return(samples_in_clades)
})
dev.off()

samples_in_clades_without_orgs_names <- sapply(samples_in_clades_without_orgs, rownames)
samples_in_clades_without_orgs_names

names(table(unlist(samples_in_clades_without_orgs_names)))
plot(sort(table(unlist(samples_in_clades_without_orgs_names))))

binary_sample_in_clade_without_orgs <- sapply(samples_in_clades_without_orgs_names, function(i) rownames(all_natgen[[which_natgen]]) %in% i)

dev.off()
pdf(paste0("figures/samples_in_empty_clades_bool_maximput_char", maximput_char, '.pdf'), height = 4.5, width = 5.5)
print(pheatmap(apply(binary_sample_in_clade_without_orgs, 2, as.numeric), cluster_cols=F))
dev.off()

##' The conclusion here is that, for the majority of samples, and unless the imputation
##' value is around ~1/3 (~ 0.004)
imputation_values[length(imputation_values)/3]
##' the samples that appear in areas which don't include any organoids are the same, i.e.
##' those that have a zero everywhere in the heatmap above (the central blue fringe)

hist(rowSums(binary_sample_in_clade_without_orgs))

exposures_generally_not_in_organoid_clades <- all_natgen[[which_natgen]][which(rowSums(binary_sample_in_clade_without_orgs) > 20),]
exposures_generally_not_in_organoid_clades <- exposures_generally_not_in_organoid_clades[!grepl('PDO', rownames(exposures_generally_not_in_organoid_clades)),]

createBarplot(exposures_generally_not_in_organoid_clades, remove_labels = T)
ggsave(paste0("figures/constantly_underrepresented_samples_maximput_char", maximput_char, '.pdf'), height = 4.5, width = 5.5)
## and now need to see if that's what I wrote in the paper. Looks like exposures with high s3 and s7

dendrogram_underrepresented <- give_dendrogram_from_imputation(plot = T, exposures = exposures_generally_not_in_organoid_clades,
                                                               impute_VALUE = 0.01, return_grob = T,
                                                               expand_vec = c(0.05, 0, 0.05, 0))

pdf(paste0("figures/constantly_underrepresented_samples_tree_maximput_char", maximput_char, '.pdf'),
    height = 6.5, width = 10.5)
cowplot::plot_grid(dendrogram_underrepresented)
dev.off()

