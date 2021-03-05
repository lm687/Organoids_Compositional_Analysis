
rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
Sys.setenv(LANG='en')

library(MetaGxOvarian)
library(TCGAbiolinks)

cached = TRUE

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

if(!cached){
  query <- GDCquery(project = c("TCGA-OV"),
                    data.category = "Transcriptome Profiling")
  
  saveRDS(query, "../objects/query")
  saveRDS(counts, "../objects/counts")
}

head(counts[[1]])
counts_bind = do.call('cbind', lapply(counts, function(i) i[,2]))
colnames(counts_bind) = subfolders
rownames(counts_bind) = counts[[1]][,1]


cases_sorted = query$results[[1]]['cases.submitter_id'][match(subfolders, query$results[[1]][,'id']),]
sample_type_sorted = query$results[[1]]['sample_type'][match(subfolders, query$results[[1]][,'id']),]
table(sample_type_sorted) ## there are no matched normals?
stopifnot(length(sample_type_sorted) == length(subfolders))

files_df = query$results[[1]][match(subfolders, unlist(query$results[[1]]['id'])),c('id', 'sample.submitter_id')]
length(counts)
dim(files_df)[1]
#saveRDS(object = files_df, file = "../objects/list_files.RDS")

## Save files for DE
## sample file
## raw counts file

## read treatment group
treatment = read.table("files/TCGA_OV_ClinData_forMariaLena.csv", sep = ",", header = TRUE, stringsAsFactors = FALSE)

# complete response VS progressive disease 
treatment = treatment[(treatment$patient.follow_ups.follow_up.primary_therapy_outcome_success) %in% c("complete remission/response", "progressive disease"),]
table(treatment$patient.follow_ups.follow_up.primary_therapy_outcome_success)
treatment$patient.follow_ups.follow_up.primary_therapy_outcome_success[treatment$patient.follow_ups.follow_up.primary_therapy_outcome_success == "complete remission/response"] = "complete_remission_or_response"
treatment$patient.follow_ups.follow_up.primary_therapy_outcome_success[treatment$patient.follow_ups.follow_up.primary_therapy_outcome_success == "progressive disease"] = "progressive_disease"
counts_bind = counts_bind[,match(treatment$id, colnames(counts_bind))]
stopifnot(treatment$id == colnames(counts_bind))

samples_df = data.frame(sample=treatment$id,
           response=treatment$patient.follow_ups.follow_up.primary_therapy_outcome_success)
write.table(samples_df,
file = "files/sample_sheet.csv", quote = FALSE, sep = ",", row.names = FALSE)
#sapply(counts, function(i) i[,1])
write.table(counts_bind, file = "files/table_raw_counts.csv", quote = FALSE, sep = ",", col.names = TRUE)


