library(gsheet)
aa <- gsheet2tbl('https://docs.google.com/spreadsheets/d/1naFASG3iOCZ2MOdLj9hb4fokTeV-OBRCNOBnU8qx3-w/edit#gid=2117194230', sheetid = "sWGS (EGA) FINAL")

sort(table(aa$...9), dec=F)
all(table(aa$...9) == 1) ## check that there are no duplicated SLX numbers
sort(table(aa$...11))
