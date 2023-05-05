getwd()

present1 <- natgen$ExposuresNatGen > 0.05
present2 <- natgen$UpdatedExposures > 0.05
mean_cooc_1 <- outer(1:7, 1:7, Vectorize(function(i,j) {
  mean(apply(present1[,c(i,j)], 1, function(k) k[1] == k[2] ))
}))

mean_cooc_2 <- outer(1:7, 1:7, Vectorize(function(i,j) {
  mean(apply(present2[,c(i,j)], 1, function(k) k[1] == k[2] ))
}))

colnames(mean_cooc_2) <- rownames(mean_cooc_2) <- colnames(mean_cooc_1) <- rownames(mean_cooc_1) <- paste0('s', 1:7)

ph1 <- pheatmap(mean_cooc_1, main='Old exposures')
ph2 <- pheatmap(mean_cooc_2, main = 'Updated exposures')

pdf("/Users/morril01/Documents/PhD/other_repos/Vias_Brenton/copy_number_analysis_organoids/figures/coocurrence_signatures.pdf", height = 3.3)
cowplot::plot_grid(plotlist = list(ph1[[4]],ph2[[4]]))
dev.off()
