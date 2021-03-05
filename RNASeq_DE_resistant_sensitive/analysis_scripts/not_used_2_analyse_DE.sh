

Rscript --vanilla other_scripts/reportDEAnalysis.R --model response --numerator "complete_remission_or_response" --denominator "progressive_disease" --numeratorName "complete remission/response" --denominatorName "progressive disease" --threshold 0.05 --deaObject ../objects/deaObjectFile --differential_all ../objects/differential_all --differential_sig ../objects/differential_sig --summary_counts ../objects/summary_counts --sampleSheet files/sample_sheet.csv  --refflat files/refFlat.txt --geneNames files/Homo_sapiens.GRCh37.87.gff3 


