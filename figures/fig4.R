rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

set.seed(234)

library(ggplotify)
renaming <- readxl::read_excel("../RNASeq_DE_resistant_sensitive/files/PDOnameProperSample_sWGS_RNAseq.xlsx")
ascites = readRDS("../copy_number_analysis_organoids/robjects/fig4_ascites.RDS")
AUC_all_df = readRDS("../survival_analysis/robjects/AUC_all_df.RDS")
PDS_PDO = readRDS("../survival_analysis/robjects/fig4_PDS_PDO.RDS")

b <- ggplot(melt(ascites, id.vars=c('sample', 'bool_ascites', 'sample_paired')), aes(x=bool_ascites, y=value, fill=variable))+geom_bar(stat = "identity")+
  facet_wrap(.~sample, scales = "free_x", nrow=2)+
  scale_fill_brewer(palette="Dark2", name = NULL)+theme(legend.position="bottom")+labs(x=NULL, y=NULL)+
  guides(fill=guide_legend(nrow=1))


PDS_PDO$sample <- renaming$PDO[match(as.character(PDS_PDO$sample), gsub("org", "", renaming$ID))]
c <- ggplot(data = PDS_PDO, aes(auc_ll5.sph, auc_ll5.org, colour = sample))+
  geom_point() +
  labs(x="Patient-derived spheroids", y="Patient-derived organoids")+
  theme_bw() +
  geom_smooth(method=lm, se=FALSE)+facet_wrap(.~sample, ncol=1)+theme(legend.position = "bottom")

d <- pheatmap(data.frame(AUC_all_df[,-(c(1:7, 21:22))]), annotation_row = AUC_all_df %>% select(PFI),
              annotation_colors = list(Resistant='red', Sensitive='purple'))

a <- cowplot::ggdraw()+draw_image("~/Desktop/fig4a.png", scale = 1)

pdf("fig4.pdf", height = 9, width = 9)
plot_grid(plot_grid(a, b, labels = c('a', 'b')),
          plot_grid(c, ggplotify::as.grob(d), labels = c('c', 'd')),
          label_size = 12, ncol=1, scale = 0.9, rel_heights = c(1.5,2))
dev.off()




