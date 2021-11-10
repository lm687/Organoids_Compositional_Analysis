# This is a github repo for the analyses of the paper "High grade serous ovarian cancer organoids as models of chromosomal instability"
Lena Morrill, 2021

## Summary of paper
There are four main figures. All the figures are created using the scripts inside the folder `figures`.

## Summary of folders and files
I only include the files which were finally used for the paper

- **copy_number_analysis_organoids** Analysis of copy number exposures, and of segments
	- `copy_number_analysis_organoids.Rmd`: R markdown file
- **RNASeq_DE_resistant_sensitive** Analysis of DE for sensitive vs resistant in TCGA samples
	- analysis_scripts
		- `3_analysis_organoids_PCA_subset_samples.R`: PCA for RNA-Seq data
		- `DGE_Genepathways_RnaSeqPipAucWo3PBias_MVorganoids_CMS20210518.Rmd`: pathway enrichment analysis
- **RNASeq_and_CN** Comparison of copy number and gene expression
	- 20191218_ViasM_BJ_orgaBrs
		- Scripts
			- `analyse_joint_counts_CN.R`
- **scDNAseq-Organoids** Scripts to create single cell DNA plots
	- code
		- `plotting_nosexchrom.R`: plotting scDNA data without sex chromosomes
		- `subclonal_structure.R`: clustering of scDNA profiles
- **survival_analysis** Scripts to create plots of the survival of organoids, and miscelaneous plots
	- code
		- `organoid_survival2.rmd` R markdown
		- `response_drugs.R`
		- `Survival.Rmd`
- **figures** folder to create the final version of the figures; explained above
	- `fig1.R`
	- `fig2.R` 
	- `fig3.R`
	- `fig4.R`

## Order in which to run the files
I would run the files in the order they appear above.

## Other folders for related analyses, but not part of the paper
- **not_used** additional related analyses that haven't made it to the paper
	- **cell_lines** Analysis of cell lines (deprecated)


## Input files
- **copy_number_analysis_organoids**
	- `copy_number_analysis_organoids.Rmd`
		- `data/organoid_exposures.rds`
		- `data/NewOrganoidNaming.csv`
		- `data/sig_data_unorm.RDS`
		- `data/Export-matrix_OV_Sigs_on_TCGA-OV_12112019.rds`
		- `data/summary.ascatTCGA.penalty70.txt`
		- `data/pcawg_CN_features.rds`
		- `data/tcga_CN_features.rds`
		- `../../cnsignatures/manuscript_Rmarkdown/data/BriTROC_absolute_copynumber.rds`
		- `data/6_TCGA_Signatures_on_BRITROC/0_BRITROC_absolute_CN.rds`
		- `data/organoid_absolute_CN.rds`
		- `data/BriTROC_CN_features.rds`
		- `data/CN_Calls_ABSOLUTE_PCAWG/OV-AU.segments.raw.rds`
		- `data/CN_Calls_ABSOLUTE_PCAWG/OV-US.segments.raw.rds`
		- `../RNASeq_DE_resistant_sensitive/files/20191218_ViasM_BJ_orgaBrs_tpm.csv`
		- `../RNASeq_DE_resistant_sensitive/files/PDOnameProperSample_sWGS_RNAseq.xlsx`
		- `data/asorg_PDO.csv` signature exposures for ascites and organoids
		- `data/ascites_exposures_20210125.rds`
- **RNASeq_DE_resistant_sensitive**
	- analysis_scripts
		- `3_analysis_organoids_PCA_subset_samples.R`
		- `DGE_Genepathways_RnaSeqPipAucWo3PBias_MVorganoids_CMS20210518.Rmd`
- **RNASeq_and_CN** 
	- 20191218_ViasM_BJ_orgaBrs
		- Scripts
			- `analyse_joint_counts_CN.R`
- **scDNAseq-Organoids**
	- code
		- `plotting_nosexchrom.R`
		- `subclonal_structure.R`
- **survival_analysis**
	- code
		- `organoid_survival2.rmd` 
		- `response_drugs.R`
		- `Survival.Rmd`
- **figures**
	- `fig1.R`
		- `../copy_number_analysis_organoids/robjects/exposures.RDS`
		- `../copy_number_analysis_organoids/robjects/dendrograminputclr.RDS`
		- `../copy_number_analysis_organoids/robjects/heatmapinputclr_with_ticks.RDS`
		- `../copy_number_analysis_organoids/robjects/heatmapinputclr.RDS`
		- `../copy_number_analysis_organoids/robjects/rank_nsegments.RDS`
		- `../copy_number_analysis_organoids/robjects/rank_ploidy.RDS`
		- `../copy_number_analysis_organoids/robjects/rank_nsegments_df.RDS`
		- `../copy_number_analysis_organoids/robjects/rank_ploidy_df.RDS`
		- `../survival_analysis/robjects/km_as_one.RDS`
		- `../survival_analysis/data/OrganoidSurvival.csv`
		- `../survival_analysis/data/OrganoidCulturesSurvival_091117OR.xlsx`
		- `../survival_analysis/data/OrganoidCulturesSurvival_091117OR_17082021.xlsx`
	- `fig2.R` 
		- `../scDNAseq-Organoids/robjects/fig2_absCN.RDS`
		- `../scDNAseq-Organoids/robjects/fig2_subclonal_hclust_nosexchromPDO2.RDS`
		- `../scDNAseq-Organoids/robjects/fig2_subclonal_hclust_nosexchromPDO3.RDS`
		- `../scDNAseq-Organoids/robjects/fig2_subclonal_hclust_nosexchromPDO6.RDS`
		- `../scDNAseq-Organoids/robjects/fig2_subclonal_hclust_nosexchromPDO2_2.RDS`
		- `../scDNAseq-Organoids/robjects/fig2_subclonal_hclust_nosexchromPDO3_2.RDS`
		- `../scDNAseq-Organoids/robjects/fig2_subclonal_hclust_nosexchromPDO6_2.RDS`
		- `../scDNAseq-Organoids/robjects/fig2_colours.RDS`
		- `../copy_number_analysis_organoids/data/absolute_profiles/PDO2_annotated_2.pdf`
		- `../copy_number_analysis_organoids/data/absolute_profiles/PDO3_annotated_2.pdf`
		- `../copy_number_analysis_organoids/data/absolute_profiles/PDO6_annotated_2.pdf`
		- `fig2_chrom_annotation_crop_nosexchrom.png`
	- `fig3.R`
		- `../RNASeq_and_CN/20191218_ViasM_BJ_orgaBrs/output/output_GRCh37/fig3_df_gene_characteristics.RDS`
		- `../RNASeq_and_CN/20191218_ViasM_BJ_orgaBrs/output/output_GRCh37/fig3_df_average_bottomCN.RDS`
		- `../RNASeq_and_CN/20191218_ViasM_BJ_orgaBrs/output/output_GRCh37/joint_counts_CN_subset.RDS`
		- `../RNASeq_DE_resistant_sensitive/objects/fig4_pca_with_gsva_annotation_NC.RDS`
		- `../RNASeq_DE_resistant_sensitive/objects/fig4_pca_with_gsva_annotation_NC_prcomp.RDS`
		- `../RNASeq_DE_resistant_sensitive/objects/fig4_df_colmeans_deseqcounts_correlation_tcga_org.RDS`
		- `../RNASeq_DE_resistant_sensitive/objects/fig3_ssgsea_repair.RDS`
	- `fig4.R`
		- `../RNASeq_DE_resistant_sensitive/files/PDOnameProperSample_sWGS_RNAseq.xlsx`
		- `../RNASeq_and_CN/20191218_ViasM_BJ_orgaBrs/Input/samples.csv`
		- `../copy_number_analysis_organoids/data/Book1_ascites.xlsx`
		- `../copy_number_analysis_organoids/robjects/fig4_ascites.RDS`
		- `../survival_analysis/robjects/AUC_all_df.RDS`
		- `../survival_analysis/robjects/fig4_PDS_PDO.RDS`
		- `../copy_number_analysis_organoids/robjects/exposures.RDS`
		- `../RNASeq_DE_resistant_sensitive/objects/fig4_fgseaResTidy.RDS`


<!-- ## Symlinks
- copy_number_analysis_organoids/BriTROC-cnsignatures-symlink: symlink to repository [CNsignatures](https://bitbucket.org/britroc/cnsignatures/src/master/) -->

## Dependencies
### Specific packages
- You have to download the bitbucket repo [CNsignatures](https://bitbucket.org/britroc/cnsignatures/src/master/) and save it in the same folder where this repo is, i.e.
	 
	                 |--- this folder
	Mother folder ---| 
	                 |--- CNsignatures

- [CompSign package](https://github.com/lm687/CompSign)
- rnaseqRpkg: internal CRUK RNA-Seq pipeline package. Send me (`lm687 at cam.ac.uk`) an email about it

### General R packages
```
library(AnnotationHub)
library(Biobase)
library(CNTools)
library(CompSign)
library(DESeq2)
library(EnvStats)
library(GSVA)
library(GSVAdata)
library(GenomicRanges)
library(MASS)
library(QDNAseq)
library(RColorBrewer)
library(ReactomePA)
library(biomaRt)
library(compositions)
library(cowplot)
library(dendextend)
library(dplyr)
library(fgsea)
library(ggdendro)
library(ggh4x)
library(ggplot2)
library(ggplotify)
library(ggrepel)
library(ggthemr)
library(grid)
library(gridExtra)
library(jcolors)
library(latex2exp)
library(lsa)
library(parallel)
library(pheatmap)
library(readxl)
library(reshape2)
library(tidyverse)
library(viridis)
require(GSVA)
require(GSVAdata)
require(biomaRt)
require(dplyr)
require(ggplot2)
require(ggrepel)
require(gridExtra)
require(jcolors)
require(pheatmap)
require(reshape2)
```
