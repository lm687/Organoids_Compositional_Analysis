library('GEOquery')
gse <- getGEO("GSE165897", GSEMatrix = TRUE)
show(gse)

gse$`GSE165897-GPL16791_series_matrix.txt.gz`$geo_accession
filePaths = getGEOSuppFiles("GSE165897")
filePaths



http://genomicsclass.github.io/book/pages/GEOquery.html