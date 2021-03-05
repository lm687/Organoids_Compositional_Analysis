##

rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

require(ggplot2)
require(ggrepel)
require(dplyr)
require(reshape2)
require(biomaRt)

##------------------------------------------------------------------------------------------------------------##
joint_counts_CN = readRDS("../output/joint_counts_CN_TPM_20210301210511.RDS")

head(joint_counts_CN)

ggplot(joint_counts_CN %>% filter(CN.L1 == '119025org' ), aes(x=CN.value, y=counts.value))+geom_point()


ggplot(joint_counts_CN, aes(x=CN.value, y=counts.value))+geom_point()+facet_wrap(.~CN.L1, scales = "free")+theme_bw()
ggsave("~/Desktop/Vias/figures/joint_counts_CN_TPM_all.pdf", width=10, height=10)

ggplot(joint_counts_CN, aes(x=CN.value, y=counts.value))+geom_point()+facet_wrap(.~CN.L1, scales = "free")+theme_bw()+
  scale_x_continuous(trans = "log2")+scale_y_continuous(trans = "log2")
ggsave("~/Desktop/Vias/figures/joint_counts_CN_TPM_all_loglog.pdf", width=10, height=8)

topvariable = read.table("~/Desktop//Vias/top_variable.txt", comment.char = "#")

ggplot(joint_counts_CN %>% filter(CN.gene_name %in% topvariable$V1), aes(x=CN.value, y=counts.value, col=CN.gene_name))+
  geom_point()

pdf("~/Desktop/Vias/figures/joint_counts_CN_TPM_topvar.pdf")
for(i in topvariable$V1){
  print(ggplot(joint_counts_CN %>% filter(CN.gene_name == i), aes(x=CN.value, y=counts.value, label=CN.L1))+
          geom_point()+geom_smooth()+ggtitle(i)+geom_label_repel()+theme_bw())
}
dev.off()
##------------------------------------------------------------------------------------------------------------##

##------------------------------------------------------------------------------------------------------------##
## now add DESeq counts
deseq_counts = read.table("/Users/morril01/Documents/PhD/other_repos/b_tape/Vias_Brenton/RNASeq_DE_resistant_sensitive/files/counts_norm.csv",
                          sep=',', header = T)

deseq_counts = (melt(deseq_counts))
colnames(deseq_counts) = c('Gene', 'DESEq.sample', 'DESeq.count')

# mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
# ensembl <- useEnsembl(biomart = "genes")
# ensembl <- useDataset(dataset = "hsapiens_gene_ensembl", mart = ensembl)
# 
# t2g <- getBM(
#   attributes = c('ensembl_gene_id', 'external_gene_name'),
#   values = deseq_counts$Gene,
#   filter = 'ensembl_gene_id',
#   mart = ensembl, useCache = FALSE)
# saveRDS(t2g, file = "~/Desktop/t2g.RDS")
t2g = readRDS("~/Desktop/t2g.RDS")

genename_ensmid_matched = t2g[match(deseq_counts$Gene, t2g$ensembl_gene_id),'external_gene_name']
deseq_counts$Gene = genename_ensmid_matched

renaming = readxl::read_xlsx("/Users/morril01/Documents/PhD/other_repos/b_tape/Vias_Brenton/RNASeq_DE_resistant_sensitive/files/PDOnameProperSample_sWGS_RNAseq.xlsx")
match(paste0(joint_counts_CN$CN.L1, joint_counts_CN$CN.gene_name))

deseq_counts$DESEq.sample = renaming$ID[match(gsub("[.]", "-", deseq_counts$DESEq.sample), renaming$sampleNameRNAseq)]

joint_counts_CN = cbind.data.frame(joint_counts_CN, deseq_counts[match(paste0(joint_counts_CN$CN.L1, joint_counts_CN$CN.gene_name), paste0(deseq_counts$DESEq.sample, deseq_counts$Gene)),])
joint_counts_CN[4,]

saveRDS(joint_counts_CN, file = "../output/joint_counts_CN_TPM_DESeq2.RDS")

ggplot(joint_counts_CN, aes(x=CN.value, y=DESeq.count))+geom_point()+facet_wrap(.~CN.L1, scales = "free")+theme_bw()
ggsave("~/Desktop/Vias/figures/joint_counts_CN_DESeq_all.pdf", width=10, height=10)

ggplot(joint_counts_CN, aes(x=CN.value, y=DESeq.count))+geom_point()+facet_wrap(.~CN.L1, scales = "free")+theme_bw()+
  scale_x_continuous(trans = "log2")+scale_y_continuous(trans = "log2")
ggsave("~/Desktop/Vias/figures/joint_counts_CN_DESeq_all_loglog.pdf", width=10, height=8)

pdf("~/Desktop/Vias/figures/joint_counts_CN_DESeq_topvar.pdf")
for(i in topvariable$V1){
  print(ggplot(joint_counts_CN %>% filter(CN.gene_name == i), aes(x=CN.value, y=DESeq.count, label=CN.L1))+
          geom_point()+geom_smooth()+ggtitle(i)+geom_label_repel()+theme_bw())
}
dev.off()


plot(joint_counts_CN$counts.value, joint_counts_CN$DESeq.count)

##------------------------------------------------------------------------------------------------------------##

ggplot(joint_counts_CN %>% filter(CN.gene_name == "CCNE1"), aes(x=CN.value, y=counts.value, label=CN.L1))+
  geom_point()+geom_smooth()+ggtitle("CCNE1")+geom_label_repel()+theme_bw()


ggplot(joint_counts_CN %>% filter(CN.gene_name == "CDK12"), aes(x=CN.value, y=counts.value, label=CN.L1))+
  geom_point()+geom_smooth()+ggtitle("CDK12")+geom_label_repel()+theme_bw()


##------------------------------------------------------------------------------------------------------------##
corDfAll = read.csv("../GexVsCnGwDeseq2/corStatTable.csv")
joint_counts_CN_TPM_subset = readRDS("../output/joint_counts_CN_TPM_subset.RDS")

## only genes included in Stephane's analysis
ggplot(joint_counts_CN %>% filter(CN.gene_name %in% joint_counts_CN_TPM_subset$nearestGeneCN.gene),
       aes(x=CN.value, y=DESeq.count))+geom_point()+facet_wrap(.~CN.L1, scales = "free")+theme_bw()+
  scale_x_continuous(trans = "log2")+scale_y_continuous(trans = "log2")+
  geom_abline(slope = 1, intercept = 0, col='blue', lty='dashed')

ggplot(joint_counts_CN %>% filter(CN.gene_name %in% joint_counts_CN_TPM_subset$nearestGeneCN.gene),
       aes(x=CN.value, y=DESeq.count))+geom_point()+theme_bw()+
  scale_x_continuous(trans = "log2")+scale_y_continuous(trans = "log2")+
  geom_abline(slope = 1, intercept = 0, col='blue', lty='dashed')+geom_smooth()


