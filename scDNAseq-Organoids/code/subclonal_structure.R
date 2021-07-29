rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(reshape2)
library(ggplot2)
library(ggdendro)
library(pheatmap)
library(readxl)
library(gridExtra)
library(ggplotify)
library(cowplot)
library(ggh4x) ## nested facets
library(dplyr)
library(EnvStats)

source("helper.R")


renaming <- readxl::read_excel("../../RNASeq_DE_resistant_sensitive/files/PDOnameProperSample_sWGS_RNAseq.xlsx")

organoid_list = c('118976org', '119148orgb', '23868org')

absCN = lapply(paste0("../data/absCN_clean_", organoid_list, ".RDS"), readRDS)
names(absCN) = organoid_list
names(absCN) = gsub("orgb", "org", names(absCN))
names(absCN) = renaming$PDO[match(names(absCN), renaming$ID)]
saveRDS(absCN, "../robjects/fig2_absCN.RDS")
organoid_list = renaming$PDO[match(gsub("orgb", "org", organoid_list), renaming$ID)]


# image(absCN[[1]])
# 
# plot(hclust(dist(absCN[[org_it]])))

# org_it = '23868org'

## remove the outliers
outliers = list()
outliers$`PDO3` = c(12, 10, 11, 7, 5, 6, 4, 3, 1, 2, 8, 9, 158)
outliers$`PDO6` = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 23, 377, 378, 379, 380, 381)
outliers$`PDO2` = c(1,2,3, 4, 5, 6, 7, 8)

saveRDS(outliers, "../robjects/scDNA_outliers.RDS")

# plot(hclust(dist(absCN[[org_it]][-outliers[[org_it]],])))

for(org_it in organoid_list){
  absCN[[org_it]] = absCN[[org_it]][-outliers[[org_it]],]
  
  absCN[[org_it]] = absCN[[org_it]][,rowSums(apply(absCN[[org_it]], 1, is.na)) == 0]
  
  ## select only top variable areas of the genome
  # absCN[[org_it]] = absCN[[org_it]][,order(apply(absCN[[org_it]] , 2, var), decreasing = T)[1:1000]]
  ## doesn't work too well
  
  absCN[[org_it]][absCN[[org_it]] > 14] = 14
  
  mat_breaks <- c(-Inf, 0:6, seq(8, 14, by=2)) - 0.01 ## adding a small number because otherwise the binning is done wrong
  col_list <- c("#2670af", "#2670af", "#8daecf", "#eaeaea", "#ffd2b6", "#f5b788", "#f39d5f", "#f17e37", "#ee1200", "#b60000", "#7a0e04", "#401004", "#000000")
  
  # mat_breaks <- c(-Inf, 0:6, seq(8, 14, by=2)) #seq(min(org_clean, na.rm = T), max(org_clean, na.rm = T), length.out = 10)
  # col_list <- c("#2670af", "#8daecf", "#eaeaea", "#ffd2b6", "#f5b788", "#f39d5f", "#f17e37", "#ee1200", "#b60000", "#7a0e04", "#401004", "#000000")
  annotation_chroms = data.frame(row.names = colnames(absCN[[org_it]]),
                                 chrom=clean_chrom(gsub("\\..*", "", colnames(absCN[[org_it]]))))
  
  ph = pheatmap::pheatmap(absCN[[org_it]]+0.02, show_colnames = FALSE, show_rownames = FALSE,
                     color             = col_list,
                     breaks            = mat_breaks, cluster_cols = FALSE,
                     annotation_col = annotation_chroms, annotation_legend=F, main=org_it)
  
  ph$gtable$grobs[[6]]$children[[1]]$vjust = 0.45
  sexchrom_bool = annotation_chroms$chrom %in% c('X', 'Y')
  annotation_chroms_nosexchrom <- data.frame(row.names=rownames(annotation_chroms)[!sexchrom_bool],
                                             chrom=annotation_chroms$chrom[!sexchrom_bool])
  ph_nosexchrom = pheatmap::pheatmap(absCN[[org_it]][,!sexchrom_bool]+0.02, show_colnames = FALSE, show_rownames = FALSE,
                          color             = col_list,
                          breaks            = mat_breaks, cluster_cols = FALSE,
                          annotation_col = annotation_chroms_nosexchrom, annotation_legend=F, main=org_it)
  ph_nosexchrom$gtable$grobs[[6]]$children[[1]]$vjust = 0.45
  saveRDS(list(mat_breaks=mat_breaks,col_list=col_list), paste0("../robjects/fig2_colours.RDS"))
  saveRDS(absCN[[org_it]], paste0("../robjects/fig2_subclonal_hclust", org_it, "_2.RDS"))
  saveRDS(ph, paste0("../robjects/fig2_subclonal_hclust", org_it, ".RDS"))
  saveRDS(ph_nosexchrom, paste0("../robjects/fig2_subclonal_hclust_nosexchrom", org_it, ".RDS"))
  # ph$gtable$grobs[[1]]$gp <- gpar(lwd = 5)
  # ph$gtable$grobs[[2]]$gp <- gpar(col = 'blue')
  
  plot(as.dendrogram(ph$tree_row))
  # quartz()
  # print(grid.arrange(ggdendrogram(as.dendrogram(ph$tree_row), no.margin = TRUE)+
  #                theme(axis.line=element_blank(),axis.text.x=element_blank(),
  #                      axis.text.y=element_blank(),axis.ticks=element_blank(),
  #                      # axis.title.x=element_blank(),
  #                      axis.title.y=element_blank(),legend.position="none")+scale_x_continuous(expand = c(0,0)),
  #              ggplot(melt(t(absCN[[org_it]])), aes(y=Var1, x=factor(Var2, levels=ph$tree_row$order), fill=value))+geom_tile()+
  #   # scale_colour_steps(breaks=mat_breaks,value=col_list)
  #     scale_fill_gradientn(colours=col_list, breaks=mat_breaks)+
  #   theme(axis.line=element_blank(),axis.text.x=element_blank(),
  #         axis.text.y=element_blank(),axis.ticks=element_blank(),
  #         axis.title.x=element_blank(),
  #         axis.title.y=element_blank(),legend.position="none",
  #         panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
  #         panel.grid.minor=element_blank(),plot.background=element_blank()), nrow=2, top=org_it))
  pdf(paste0("../plots/subclonal_hclust", org_it, ".pdf"))
  print(ph)
  dev.off()
  # dendro <- ggdendrogram(as.dendrogram(ph$tree_row), no.margin = TRUE)+
  #                  theme(axis.line=element_blank(),axis.text.x=element_blank(),
  #                        axis.text.y=element_blank(),axis.ticks=element_blank(),
  #                        # axis.title.x=element_blank(),
  #                        axis.title.y=element_blank(),legend.position="none")+scale_x_continuous(expand = c(0,0))+
  #   coord_flip()
  # cowplot::plot_grid(ggplotify::as.grob(ph), dendro)
  # dev.off()
  
}

# quartz()
ggplot(melt(t(absCN[[org_it]])), aes(y=Var1, x=factor(Var2, levels=ph$tree_row$order), fill=value))+
  geom_tile()+
  # scale_colour_steps(breaks=mat_breaks,value=col_list)
  scale_fill_gradientn(colours=col_list, breaks=mat_breaks)
# 
# 

data_orgs <- lapply(organoid_list, function(org_it){
  .x <- readRDS(paste0("../robjects/fig2_subclonal_hclust", org_it, "_2.RDS"))
  ## remove sex chroms
  chrom=clean_chrom(gsub("\\..*", "", colnames(.x)))
  .x[,!(chrom %in% c('X', 'Y'))]
})
names(data_orgs) <- organoid_list

ggplot(melt(list(apply(data_orgs[[1]], 1, var),
     apply(data_orgs[[2]], 1, var),
     apply(data_orgs[[3]], 1, var))), aes(x=value, col=L1, group=L1))+geom_density()


ggplot(melt(data_orgs[[1]][,1:100]), aes(x=Var2, y=value))+geom_violin()+geom_jitter(size=0.2)

data_orgs_centered <- lapply(data_orgs, function(j) sweep(j, 2, apply(j, 2, mean), '-'))
names(data_orgs_centered) <- organoid_list

ggplot(melt(data_orgs_centered[[1]][,1:100]), aes(x=Var2, y=value))+geom_violin()+geom_jitter(size=0.2)

# apply(data_orgs_centered[[1]]

ggplot(melt(data_orgs_centered[[1]]), aes(x=Var2, y=value))+geom_violin()

data_orgs_centered_sd <- lapply(data_orgs_centered, function(j){
  cbind.data.frame(pos=colnames(j),
                 sd_bin=apply(j, 2, sd))})

ggplot(melt(data_orgs_centered_sd),
       aes(x=pos, y=value, group=L1, col=L1))+geom_line()+
  labs(y='Standard deviation')

## compute confidence intervals
data_orgs_centered_CI <- data.frame(do.call('rbind', lapply(data_orgs_centered, function(j) (t(apply(j, 2, quantile, c(0.05, 0.95)))))))
data_orgs_centered_CI$L1 = rep(organoid_list, sapply(data_orgs_centered, ncol))
data_orgs_centered_CI$pos=rownames(data_orgs_centered_CI)

put_limits <- function(x, upperlim=1, lowerlim=-1){
  a <- x
  if(x > upperlim) a <- upperlim
  if(x < lowerlim) a <- lowerlim
  return(a)
}

data_orgs_centered_CI$censoredX5 = sapply(data_orgs_centered_CI$X5., put_limits)
data_orgs_centered_CI$censoredX95 = sapply(data_orgs_centered_CI$X95., put_limits)
data_orgs_centered_CI$chrom = clean_chrom(gsub("\\..*", "", rownames(data_orgs_centered_CI)))
data_orgs_centered_CI$chrom = factor(data_orgs_centered_CI$chrom, levels=c(as.character(1:22)))
data_orgs_centered_CI$pos_num = as.numeric(sapply(rownames(data_orgs_centered_CI), function(i) strsplit(i, '[.]')[[1]][2]))
ggplot(data_orgs_centered_CI)+
  geom_ribbon(aes(x=pos_num, ymin=censoredX5, ymax=censoredX95, group=L1))+
  geom_point(aes(x=pos_num, y=0, col=( (censoredX95-censoredX5)< 0.1 )), size=0.1)+
  facet_nested(L1~chrom, scales = "free", space="free_x")+theme_bw()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  lims(y=c(-1, 1))+theme(legend.position = "bottom")+
  labs(x='90% confidence intervals of the scaled CN values of bins, across organoids',
       y='Value of scaled CN of bin')+labs(col='Clonal bin')
ggsave(paste0("../plots/subclonal_hclust_CI_across_genome_2.pdf"), width = 14)

x = data_orgs_centered_CI %>% filter(L1=='PDO3', chrom==7)
plot(x$censoredX5, type='l')

pnorm(q = 0.95)
cbind.data.frame(org=organoid_list, percentage_above_083sd=sapply(data_orgs_centered_sd, function(i) mean((i$sd_bin > 0.83))))

## plot the absolute CN and the variace, because there might be heteroscedasticity
data_orgs_mean_CN <- lapply(data_orgs, function(j){
  cbind.data.frame(pos=colnames(j),
                   mean=colMeans(j))})

plot((data_orgs_mean_CN[[1]]$mean),
     (data_orgs_centered_sd[[1]]$sd_bin)
     )

df_heteroscedasticity <- data.frame(do.call('rbind', lapply(1:3, function(i) cbind(data_orgs_mean_CN[[i]]$mean,
                              data_orgs_centered_sd[[i]]$sd_bin))))
colnames(df_heteroscedasticity) <- c('Mean', 'Sd')
df_heteroscedasticity$PDO <- rep(organoid_list, sapply(1:3, function(i) length(data_orgs_centered_sd[[i]]$sd_bin)))
head(df_heteroscedasticity)

ggplot(df_heteroscedasticity, aes(x=Mean, y=Sd))+geom_point()+
  geom_smooth(se = FALSE, method = lm)+
  facet_wrap(.~PDO)+theme_bw()
ggsave(paste0("../plots/subclonal_hclust_CI_across_genome_heteroscedasticity.pdf"), width = 5, height = 2.5)

mean_sd_lm <- lapply(organoid_list, function(org) lm(Sd ~ Mean,
                         data = df_heteroscedasticity[df_heteroscedasticity$PDO == org,]))
names(mean_sd_lm) <- organoid_list
sapply(mean_sd_lm, coef)

## for each bin, in each organoid, compute the observed standard deviation, and se
## how it compares to the expected standard deviation, with a Chi squared test
pvals_chisq <- lapply(organoid_list, function(org){
  sapply(1:sum(df_heteroscedasticity$PDO == org), function(pos_it){
    .current <- df_heteroscedasticity[df_heteroscedasticity$PDO == org,][1,]
    .expected_sd <- coef(mean_sd_lm[[org]])[1] + .current$Mean*coef(mean_sd_lm[[org]])[2]
    .expected_var <- .expected_sd**2
    .res_chi <- EnvStats::varTest(data_orgs[[org]][,pos_it], sigma.squared = .expected_var, alternative="greater")
    .res_chi$p.value
  })
})
names(pvals_chisq) <- organoid_list
xtable::xtable(data.frame(fraction_heterogeneous=sapply(pvals_chisq, function(j) mean(j < 0.05))))
