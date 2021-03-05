## Differential expression of *sensitive* vs *resistant* OV samples in TCGA
### Lena Morrill, May 2020


Run code from the folder `analysis_scripts`

#### Dependencies
- Must have package `rnaseqRpkg` from the CRUK core facilities installed

#### Files to download

- refFlat file downloaded from http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/
- Human-readable names of genes from ftp://ftp.ensembl.org/pub/grch37/current/gff3/homo_sapiens/

#### Running code

1. To get the data from TCGA in the right shape:

  ```0_preparing_files.R```

2. For DE pipeline from counts data:

  ```sh 1_run_DE.sh```

3. To analyse the DE results:
``` sh 2_analyse_DE.sh ```
