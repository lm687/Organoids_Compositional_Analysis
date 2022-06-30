## excerpts copied from gexVsCnGw_LM.R

rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
R.version

.libPaths( "/mnt/scratcha/fmlab/morril01/software/miniconda3/envs/vias_cor/lib/R/library/")
.libPaths( )

library(GenomicRanges, lib.loc = "/mnt/scratcha/fmlab/morril01/software/miniconda3/envs/vias_cor/lib/R/library/")
library(ensembldb, lib.loc = "/mnt/scratcha/fmlab/morril01/software/miniconda3/envs/vias_cor/lib/R/library/")
library(parallel, lib.loc = "/mnt/scratcha/fmlab/morril01/software/miniconda3/envs/vias_cor/lib/R/library/")
source("../../../copy_number_analysis_organoids/helper_functions_granges.R")

gtf.file <- file.path("../Data/", "Homo_sapiens.GRCh37.87.gtf.gz")
sqlite_file <- 'Homo_sapiens.GRCh37.87.sqlite'
sqlite_path <- file.path("../Data/", sqlite_file)

if(!file.exists(sqlite_path)) {
  ## generate the SQLite database file
  ensembldb::ensDbFromGtf(gtf=gtf.file, path = ref_dir, outfile=sqlite_file)
}
EnsDb.Hsapiens.v87 <- ensembldb::EnsDb(sqlite_path)

# Genes, used to annotated the TPM matrix to send to Maria
ag <- ensembldb::genes(EnsDb.Hsapiens.v87, filter=list(AnnotationFilter::GeneBiotypeFilter('protein_coding')), return.type="DataFrame") 
ag

chromosomes_metadata <- ensembldb::seqinfo(EnsDb.Hsapiens.v87)
chromlens <- (cbind.data.frame(Chrom=seqnames(chromosomes_metadata),
                               Length=seqlengths(chromosomes_metadata)))
head(chromlens)

ag_subsetchrom <- ag[!(ag$seq_name %in% c("MT", "X", "Y")) & !grepl("GL",ag$seq_name),]
unique(ag_subsetchrom$seq_name)
gr_genes = as(paste0(ag_subsetchrom$seq_name, ':', ag_subsetchrom$gene_seq_start, '-', ag_subsetchrom$gene_seq_end), "GRanges")
gr_genes


levels(seqnames(gr_genes))
seqnames(gr_genes) <- droplevels(seqnames(gr_genes))
levels(seqnames(gr_genes))
unique(seqnames(gr_genes))

## get segments from the unified version created in copy_number_analysis_organoids.Rmd
segments_Britroc <- readRDS("../../../copy_number_analysis_organoids/data/clean_segtables/segtables_BriTROC_absolute_copynumber.RDS")
segments_ICGCAU <- readRDS("../../../copy_number_analysis_organoids/data/clean_segtables/segtables_ICGC_absolute_copynumber_AU.RDS")
segments_ICGCUS <- readRDS("../../../copy_number_analysis_organoids/data/clean_segtables/segtables_ICGC_absolute_copynumber_US.RDS")
segments_TCGA <- readRDS("../../../copy_number_analysis_organoids/data/clean_segtables/segtables_TCGA_absolute_copynumber.RDS")
segments_early <- readRDS("../../../copy_number_analysis_organoids/data/early_segTable.rds")
names_early <- paste0('early_', names(segments_early))
segments_early <- lapply(segments_early, function(i){
  i$start = as.numeric(i$start)
  i$end = as.numeric(i$end)
  i$segVal = as.numeric(i$segVal)
  i
})
names(segments_early) <- names_early
segments_Britroc_Phil <- read.table("../../../copy_number_analysis_organoids/data/britroc_30kb_ds_absCopyNumber_segmentTable.tsv",
                                    stringsAsFactors = F, header = T)
britroc1_phil_samples <- unique(segments_Britroc_Phil$sample)
segments_Britroc_Phil <- lapply(britroc1_phil_samples, function(sam) segments_Britroc_Phil[segments_Britroc_Phil$sample == sam,-ncol(segments_Britroc_Phil)])
names(segments_Britroc_Phil) <- britroc1_phil_samples
head(segments_TCGA[[1]])
head(segments_Britroc_Phil[[1]])
head(segments_early[[1]]$segVal)
head(segments_Britroc_Phil[[1]]$segVal)

table(names(segments_Britroc_Phil) %in% names(segments_Britroc))
table(names(segments_Britroc) %in% names(segments_Britroc_Phil))

names(segments_Britroc_Phil) <- paste0(names(segments_Britroc_Phil), 'PS')

sapply(list(segments_Britroc, segments_Britroc_Phil, segments_ICGCAU, segments_ICGCUS, segments_TCGA), length)

all_segs <- c(segments_Britroc, segments_Britroc_Phil, segments_ICGCAU, segments_ICGCUS, segments_TCGA, segments_early)
length(all_segs)


names_run <- list.files("../output/output_GRCh37/all_CN_states_per_gene/")
names_not_run <- names(all_segs)[!(names(all_segs) %in% gsub(".RDS", "", names_run))]

head(all_segs$IM_100)
head(all_segs$`TCGA-04-1331`)

head(all_segs[[sample_name_it]])

change_column_names <- function(i){
  colnames(i)[colnames(i) == "startpos"] <- "start"
  colnames(i)[colnames(i) == "endpos"] <- "end"
  i
}

names_not_run

CN_averages = lapply(names_not_run, function(sample_name_it){
  CN_averages <- give_CN_per_gene(segment_arg = change_column_names(all_segs[[sample_name_it]]))
  averaged_CN_df = cbind.data.frame(gene_name=(ag_subsetchrom$symbol),
                                    CN=CN_averages$CN_bin_averaged)
  saveRDS(averaged_CN_df, paste0("../output/output_GRCh37/all_CN_states_per_gene/", sample_name_it, ".RDS"))
})


hist(averaged_CN_df$CN)

