remove.packages("googlesheets")

googlesheets::gs_url("https://docs.google.com/spreadsheets/d/1ulHgNIMsoFECt8cim2va5pTXu3LaCJlqaiSfvTCPAxs/edit#gid=0")

install.packages('gsheet')
library(gsheet)
summary_file = gsheet2tbl('https://docs.google.com/spreadsheets/d/1ulHgNIMsoFECt8cim2va5pTXu3LaCJlqaiSfvTCPAxs/edit#gid=0"')

apply(summary_file, 2, function)
\Checkmark
xtable::xtable(summary_file)
