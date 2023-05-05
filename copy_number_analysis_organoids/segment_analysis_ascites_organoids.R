rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(ggrepel)
source("../../../other_repos/britroc-cnsignatures-bfb69cd72c50/main_functions.R")
source("../../../other_repos/britroc-cnsignatures-bfb69cd72c50/helper_functions.R")


## analysis of segments of ascites and organoids
# ascites_segments <- readRDS("data/20220511BH_ascites_absoluteCN_bestfit.rds")
ascites_segments1 <- readRDS("data/20220629BH_ascites_absoluteCN_bestfit.rds")
ascites_segments2 <- readRDS("data/20220718BH_absoluteCN_bestfit.rds")
colnames(ascites_segments1@assayData$copynumber)
colnames(ascites_segments2@assayData$copynumber)
organoids_segments <- readRDS("data/organoid_absolute_CN.rds")

## get copy number in bins
ascites_segments_segtable <- getSegTable(ascites_segments1)
# View(ascites_segments_segtable)
organoids_segments_segtable <- getSegTable(organoids_segments)

ascites_segments_assay <- cbind(ascites_segments1@assayData$copynumber,
                                ascites_segments2@assayData$copynumber)
organoids_segments_assay <- organoids_segments@assayData$copynumber

all(rownames(ascites_segments_assay) == rownames(organoids_segments_assay))

ascc <- readxl::read_excel("data/matching_ascites_samples_Lena.xlsx")
ascc$sWGS_barcode <- gsub("-", "_", ascc$sWGS_barcode)
colnames(ascites_segments_assay)[colnames(ascites_segments_assay) == "16421.D705tp_D503tp"] = "16421.D705_D503"
colnames(ascites_segments_assay)[colnames(ascites_segments_assay) == "16421.D705tp_D501tp"] = "16421.D705_D501"

paste0(ascc$sWGS_SLX, '.', ascc$sWGS_barcode)[is.na(match(colnames(ascites_segments_assay), paste0(ascc$sWGS_SLX, '.', ascc$sWGS_barcode)))]

colnames(ascites_segments_assay)[!(colnames(ascites_segments_assay) %in% paste0(ascc$sWGS_SLX, '.', ascc$sWGS_barcode))]
paste0(ascc$sWGS_SLX, '.', ascc$sWGS_barcode)[!(paste0(ascc$sWGS_SLX, '.', ascc$sWGS_barcode) %in% colnames(ascites_segments_assay) )]

colnames(ascites_segments_assay) <- paste0(ascc$Derived_organoid[match(colnames(ascites_segments_assay), paste0(ascc$sWGS_SLX, '.', ascc$sWGS_barcode))],
       '_', colnames(ascites_segments_assay))

colnames(organoids_segments_assay) <- paste0(ascc$Derived_organoid[match(gsub("org", "", colnames(organoids_segments_assay)), ascc$CInumber)], '_', colnames(organoids_segments_assay))

all_segments_assay <- cbind(ascites_segments_assay, organoids_segments_assay)
all_segments_assay <- all_segments_assay[!apply(all_segments_assay, 1, function(i) all(is.na(i))),]

pca_segments <- prcomp(t(all_segments_assay))
umap_segments <- umap::umap(t(all_segments_assay))
pca_segments

ggplot(cbind.data.frame(name=rownames(pca_segments$x), pca_segments$x[,1:2],
       type=ifelse(grepl("org", rownames(pca_segments$x)),'organoid', 'ascites'),
       col=gsub("\\_.*","",rownames(pca_segments$x))),
       aes(x=PC1, y=PC2, label=name,col=col, shape=type))+
  geom_point()+
  geom_label_repel()+theme_bw()
ggsave("figures/ascites_organoids_matching_PCA.pdf", width = 8, height = 8)


ggplot(cbind.data.frame(name=rownames(pca_segments$x), umap_segments$layout,
                        type=ifelse(grepl("org", rownames(pca_segments$x)),'organoid', 'ascites'),
       col=gsub("\\_.*","",rownames(pca_segments$x))),
       aes(x=`1`, y=`2`, label=name, col=col, shape=type))+
  geom_point()+
  geom_label_repel()+theme_bw()
ggsave("figures/ascites_organoids_matching_umap.pdf", width = 8, height = 8)

dim(ascites_segments)
dim(organoids_segments)
dim(ascites_segments_assay)
dim(organoids_segments_assay)

organoids_segments_assay <- organoids_segments_assay[,sapply(paste0('PDO', 1:18, '_'), function(i){ 
  grep(i,(colnames(organoids_segments_assay)))
  })]
ascites_segments_assay <- ascites_segments_assay[,sapply(paste0('PDO', 1:18, '_'), function(i){ 
  grep(i,(colnames(organoids_segments_assay)))
})]

cbind(colnames(organoids_segments_assay), colnames(ascites_segments_assay))
  
j=1
plts_scatter <- lapply(1:18, function(j) ggplot(data.frame(org=organoids_segments_assay[,j], asc=ascites_segments_assay[,j]),
       aes(y=org, x=asc))+
  geom_density_2d_filled()+theme_bw()+#+scale_y_continuous(trans = "log")+scale_x_continuous(trans = "log")
  geom_smooth(method = "lm")+labs(x='CN in ascites', y='CN in organoid')+ggtitle(paste0('PDO', j))+guides(fill='none'))

pdf("figures/scatter_ascites_organoids.pdf")
do.call('grid.arrange', plts_scatter)
dev.off()
