rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
Sys.setenv(LANG='en')

library(MetaGxOvarian)
library(TCGAbiolinks)

cached = TRUE

### Clinical data ####
## this file had to be modifed: replace "'--" by "--"
# biospecimen_data <- read.table("biospecimen.cases_selection.2020-05-04/sample.tsv",  sep = "\t", fill = TRUE,
#                                stringsAsFactors = FALSE, header = TRUE)


#### Read in counts data ####
flder = "../data/"
subfolders = list.files(flder)

if(cached){
  counts = readRDS("../objects/counts")
  query = readRDS("../objects/query")
}

if(!cached){
counts = lapply(subfolders, function(sf){
  fle_counts = list.files(paste0(flder, sf))
  fle_counts = fle_counts[grepl("htseq.counts", fle_counts)]
  print(fle_counts)
  fle = paste0(flder, sf, "/", fle_counts)
  if(grepl("[.]gz", fle_counts)){
    system(paste0("gunzip -d ", fle))
    fle_unzipped = gsub(".gz", "", fle)
  }else{
    fle_unzipped = fle
  }
  .x = read.table(fle_unzipped, stringsAsFactors = FALSE)
  system(paste0("gzip  ", fle_unzipped))
  .x
})
}

head(counts[[1]])

# match(subfolders, biospecimen_data$sample_id)
# match(subfolders, biospecimen_data$case_id)
# biospecimen_data$sample_type

if(!cached){
  query <- GDCquery(project = c("TCGA-OV"),
                                   data.category = "Transcriptome Profiling")
  
  saveRDS(query, "../objects/query")
  saveRDS(counts, "../objects/counts")
}

cases_sorted = query$results[[1]]['cases.submitter_id'][match(subfolders, query$results[[1]][,'id']),]
sample_type_sorted = query$results[[1]]['sample_type'][match(subfolders, query$results[[1]][,'id']),]
table(sample_type_sorted) ## there are no matched normals?
stopifnot(length(sample_type_sorted) == length(subfolders))

files_df = query$results[[1]][match(subfolders, unlist(query$results[[1]]['id'])),c('id', 'sample.submitter_id')]
length(counts)
dim(files_df)[1]
saveRDS(object = files_df, file = "../objects/list_files.RDS")

