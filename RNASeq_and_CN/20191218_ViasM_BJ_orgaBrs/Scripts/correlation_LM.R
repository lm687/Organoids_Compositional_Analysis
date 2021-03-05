rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

absCnDf <- read.table( "../Input/organoidAbsolute_SegTable.txt", sep="\t", header=TRUE)

## get absolute copy number
absCnDf

gr_genes = GenomicRanges::makeGRangesFromDataFrame(df = cbind.data.frame(start=ag$gene_seq_start, end=ag$gene_seq_end, seqnames=ag$seq_name))
gr_genes
gr_genes$name_gene = ag$gene_name

head(gr_genes)
