# call with Rscript
# See also https://readxl.tidyverse.org/articles/articles/readxl-workflows.html
require(dplyr)
require(readxl)
require(readr)

readxl::read_xlsx("NewOrganoidNaming.xlsx") %>%
  write_csv("names_orgs.csv")
