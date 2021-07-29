## for DE between resistant and sensitive (TCGA)
Rscript --vanilla other_scripts/runDEAnalysis.R --sampleSheet files/sample_sheet.csv --model ~response --counts_raw files/table_raw_counts.csv --deaObjectFile ../objects/deaObjectFile

## for DE between the two clusters in the data (organoids)
Rscript --vanilla other_scripts/runDEAnalysis.R --sampleSheet files/sample_sheet_organoids_cluster.csv --model ~cluster --counts_raw ../../RNASeq_and_CN/20191218_ViasM_BJ_orgaBrs/RnaSeqPip/counts/counts_raw_subsetno3pbias.csv --deaObjectFile ../objects/deaObjectFile_organoids_cluster

