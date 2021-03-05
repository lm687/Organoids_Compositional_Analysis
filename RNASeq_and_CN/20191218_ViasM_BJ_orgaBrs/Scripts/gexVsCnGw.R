# Context:
 #-------------
# Maria needs to check the link between copy number of gene expression for a manuscript on organoids ready for submission (little more information provided)

# Aim:
 #-------------
# Check correlation between absolute copy number and gene expression in 18 samples.
# For a set of genes Maria is interested in,
# or across the genome depending on Maria's decision

# Note:
 #-------------
# the code comprises chunks for featureCounts and or kallisto text output that
# we do  not need, but want to keep until the project is finished. 
# 20200705: adapting code to DESeq2 matrix (replacing the kallisto TPM matrix)

# Data sets
 #-------------

# Sample sheet
# samples.xlsx: sample sheet with the two sets of sample IDs used: one for CN, the other for gene expression
# samples.csv: text version of samples.xlsx

# CN:
# organoidAbsolute_SegTable.txt: absolute copy number for 18 samples (columns) and 30 kb bins along autosomes 
# organoidAbsolute_geneList.txt: table in the long format with genes of interest their coordinates, absolute copy number ('segVal') and sample ('sample')

# gene expression
# featureCounts, now replaced with kallisto 
# kallisto output in Kallisto2
# 2020705:  DESeq2 normalised counts

#countType <- "kallisto"
countType <- "deseq2"

# Tool
 #-------------

# cluster # /home/bioinformatics/software/R/R-3.6.1/bin/R
# osx # R version 3.6.1

# Functions
 #-------------

# Packages
 #-------------
library(dplyr)
library(GenomicRanges)
library(ggplot2)
library(ggrepel)

# Variables
 #-------------
projName <- "20191218_ViasM_BJ_orgaBrs"
projDir <- "~/MyProjectsSvn/SvnRepoForStdRnaSeq/Brenton/ViasM/20191218_ViasM_BJ_orgaBrs/trunk"
#projDir <- "/mnt/scratchb/bioinformatics/baller01/20191218_ViasM_BJ_orgaBrs"
setwd(projDir)

# Assess correlation between gene expression and copy number.
 #-------------

# Read sample sheet in
 #-------------
# samples.csv has three columns:
# organoid_name;
# JBLAB-number; # eg JBLAB-19902
# OV04_number
tmpFn <- "Data/samples.csv"
splSht <- read.table(tmpFn, sep=";", header=TRUE)
splSht <- splSht %>% mutate(SampleName=gsub("-","",JBLAB.number)) 

# Annotation
 #-------------

setwd(projDir)

# GRCh38
if(FALSE)
{
ref_dir <- 'KallistoIdx/homo_sapiens'
gtf.file <- file.path(ref_dir, "Homo_sapiens.GRCh38.96.gtf.gz")
sqlite_file <- 'Homo_sapiens.GRCh38.96.sqlite'
sqlite_path <- file.path(ref_dir, sqlite_file)

if(!file.exists(sqlite_path)) {
    ## generate the SQLite database file
    ensembldb::ensDbFromGtf(gtf=gtf.file, path=ref_dir, outfile=sqlite_file)
}
EnsDb.Hsapiens.v96 <- ensembldb::EnsDb(sqlite_path)

# Genes, used to annotate the TPM matrix to send to Maria
ag <- ensembldb::genes(EnsDb.Hsapiens.v96, filter=list(AnnotationFilter::GeneBiotypeFilter('protein_coding')), return.type="DataFrame") 
ag
}

# GRCh37
# Not from kallisto
ref_dir <- 'Data'
gtf.file <- file.path(ref_dir, "Homo_sapiens.GRCh37.87.gtf.gz")
sqlite_file <- 'Homo_sapiens.GRCh37.87.sqlite'
sqlite_path <- file.path(ref_dir, sqlite_file)

if(!file.exists(sqlite_path)) {
    ## generate the SQLite database file
    ensembldb::ensDbFromGtf(gtf=gtf.file, path=ref_dir, outfile=sqlite_file)
}
EnsDb.Hsapiens <- ensembldb::EnsDb(sqlite_path)

# Genes, used to annotate the TPM matrix to send to Maria
ag <- ensembldb::genes(EnsDb.Hsapiens, filter=list(AnnotationFilter::GeneBiotypeFilter('protein_coding')), return.type="DataFrame") 
ag


# Read CNV data in:
 #-------------
# 1: the whole table
tmpFn <- "Data/organoidAbsolute_SegTable.txt"
absCnDf <- read.table(tmpFn, sep="\t", header=TRUE)
print(head(absCnDf[,-1]))
#hist(absCnDf[,-1], n=50)
#absCnNbNa <- apply(absCnDf[,-1], 1, function(x){sum(is.na(x))})
	
absCnVar <- apply(absCnDf[,-1], 1, var)
summary(absCnVar)
#hist(absCnVar, n=50)
absCnVarOrd <- order(absCnVar)

absCnMean <- apply(absCnDf[,-1], 1, mean)
summary(absCnMean)

absCnMedian <- apply(absCnDf[,-1], 1, median)
summary(absCnMedian)

#Find most variable segment
absCnVarQ90 <- quantile(absCnVar, probs=0.90)
absCnVarQ90rd <- round(quantile(absCnVar, probs=0.90),1)

table(absCnVar >= absCnVarQ90, absCnVar >= absCnVarQ90rd)

segToKeep <- absCnVar >= absCnVarQ90

absCnDfOrig <- absCnDf

absCnDf <- absCnDfOrig[segToKeep,]

# Make GRanges from first column
 #-------------
# Get gene list
# Get overlap

# derive data frame to feed makeGRangesFromDataFrame()
tmpl <- strsplit(as.character(absCnDf[,1]),"[:-]")
cnSeg <- do.call(rbind, lapply(tmpl,as.numeric))
colnames(cnSeg) <- c("chr", "start", "end")
# make GRanges, keeping the CN values.
cnSegRg <- makeGRangesFromDataFrame(cbind(cnSeg,absCnDf[,-1]), keep.extra.columns=TRUE)

cnSegRgOrig <- cnSegRg
cnSegRg <- reduce(cnSegRgOrig, ignore.strand=T)

tmpFn <- sprintf("%s.csv", "cnSegRg")
##write.csv(data.frame(cnSegRg), file=tmpFn)

# find nearest gene:
# first get GRanges for genes:
agRg <-  ag %>% data.frame() %>% dplyr::rename(chrom = seq_name,
		      start = gene_seq_start,
		      end = gene_seq_end) %>%
	makeGRangesFromDataFrame(keep.extra.columns=TRUE)

# find nearest:
x <- nearest(cnSegRg, agRg, select=c("all"), ignore.strand=TRUE)

genesToGet <- mcols(agRg[subjectHits(x),])$gene_id

# split the coord into chr, start and end to filter on:
absCnDf2 <- absCnDf %>% mutate(coord = X) %>%
	tidyr::separate(X, c("chrom", "start", "end"), sep = "[:-]")
absCnDf2$chromNum <- as.numeric(absCnDf2$chrom)
absCnDf2$startNum <- as.numeric(absCnDf2$start)
absCnDf2$endNum <- as.numeric(absCnDf2$end)

# split the coord into chr, start and end to filter on:
absCnDf3 <- absCnDfOrig %>% mutate(coord = X) %>%
	tidyr::separate(X, c("chrom", "start", "end"), sep = "[:-]")
absCnDf3$chromNum <- as.numeric(absCnDf3$chrom)
absCnDf3$startNum <- as.numeric(absCnDf3$start)
absCnDf3$endNum <- as.numeric(absCnDf3$end)

# Read kallisto TMP matrix in
# or the DESeq2 normalised count
 #-------------------------

if(countType == "kallisto")
{
	tmpFn <- sprintf("%s/TpmMat/%s_tpm.csv", projDir, projName)
	tpmMat <- read.table(tmpFn, header=TRUE, sep=",")
} else if(countType == "deseq2")
{
	tmpFn <- sprintf("%s/RnaSeqPip/counts/counts_norm.csv", projDir)
	tpmMat <- read.table(tmpFn, header=TRUE, sep=",")
	# one row per gene (ENSG) and one column per sample (JBLAB ID)
	# same as for kallisto above, but without the gene annotation
	# so also get the TPM file to fetch annotation
	tmpFn <- sprintf("%s/TpmMat/%s_tpm.csv", projDir, projName)
	tpmKal <- read.table(tmpFn, header=TRUE, sep=",")
	kalCol <- c("gene_id", "gene_name","entrezid","gene_biotype","gene_seq_start","gene_seq_end","seq_name","seq_strand","symbol")
	tpmMat <- tpmMat %>%
		dplyr::rename(gene_id = X) %>%	
		dplyr::left_join(tpmKal[,kalCol], by="gene_id")
	# replace column names
	tmpInd <- which(grepl("JBLAB", colnames(tpmMat)))
	splNames <- gsub("\\.", "-", colnames(tpmMat)[tmpInd])
	tmpInd2 <- match(splNames, splSht$JBLAB.number)
	splNamesNew <- splSht[tmpInd2, "organoid_name"]
	colnames(tpmMat)[tmpInd] <- splNamesNew
} else
{
	stop("Unknown countType '%s'", countType)
}

# correlation gene expression and CNV.
 #-------------------------

setwd(projDir)

# format the TPM matrix
 #-------------

# geneDf uses 'sample' eg 119025org, aka 'organoid_name' in splSht

genesToKeep <- genesToGet
tpmMatOrig <- tpmMat 

# Filter genes of interest
if(countType == "kallisto")
{
	tpmMat <- tpmMatOrig %>% dplyr::filter(gene_id %in% genesToKeep) %>%
	mutate(gene_stg = sprintf("%s::%s", gene_id, gene_name)) %>%
	select(gene_stg, splSht$SampleName)
} else {
	tpmMat <- tpmMatOrig %>% dplyr::filter(gene_id %in% genesToKeep) %>%
	mutate(gene_stg = sprintf("%s::%s", gene_id, gene_name)) %>%
	#select(gene_name, grep("org$", colnames(tpmMatOrig)))
	select(gene_stg, grep("org$", colnames(tpmMatOrig)))
}

# Write the TPM matrix for genes of interest to file:
tmpFn <- sprintf("%s/TpmMat/%s_InGeneList_tpm.csv", projDir, projName)
##write.csv(tpmMat, file=tmpFn, row.names = FALSE)

# convert TPM matrix to long format to match that keeping CN sent by Maria:
##tpmMatLong <- tpmMat %>% tidyr::pivot_longer(-gene_name, names_to = "SampleName", values_to = "tpm")
tpmMatLong <- tpmMat %>% tidyr::pivot_longer(-gene_stg, names_to = "SampleName", values_to = "count")

# add cn to gene to make geneDf

# for the first gene in region
# Read gene list in. Not any more.
# Now parse CN table.
 #-------------

#setwd("GexVsCnGw") TPM
setwd("GexVsCnGwDeseq2") # DESeq2

# have list to keep outcome
# one entry per region 
myls <- vector("list", length = length(cnSegRg)) 
# to keep a data frame with one row per gene and columns keeping gene
# coordinates, and correlation estimate and p-value.

names(myls) <- data.frame(cnSegRg) %>%
	mutate(regionName = sprintf("%s-%s-%s-%s",
				seqnames, start, end,
				width/1000
				)) %>%
	pull(regionName)

# For each region:
for (i in 1:length(cnSegRg)) # 
{
	tmpBit <- sprintf("reg_%s", names(myls)[i])

	# find nearest genes:
	#-------------
	x <- nearest(cnSegRg[i,], agRg, select=c("all"), ignore.strand=TRUE)
	print(agRg[subjectHits(x),])

	genesToGet <- mcols(agRg[subjectHits(x),]) %>%
		data.frame() %>%
		mutate(gene_stg = sprintf("%s::%s", gene_id, gene_name)) %>%
		pull(gene_stg)
	print(genesToGet)

	# have list for that region
	# with one entry per gene

	myls2 <- vector("list", length = length(genesToGet)) 
	names(myls2) <- genesToGet

	# get coordinates of that region
	#-------------
	# to filter segments on.
	chr <- as.numeric(as.vector(seqnames(cnSegRg[i,]))) # see also runValue()
	# have flanking bits for plot:
	l1 <- start(cnSegRg[i,]) - 1000000
	l2 <- end(cnSegRg[i,]) + 1000000

	# recall: absCnDf2 is same as absCnDf but with coord split in chr, start and end
	# to filter on
	# recall: absCnDf is the df with one segement per row with coord and cn with one
	# column per sample, for the regions of interest, not the whole table.
	# now subset that for the current region:
	absCnDf2Sub <- absCnDf2 %>% dplyr::filter( chromNum == chr & l1 <= startNum & endNum <= l2)
	print(absCnDf2Sub[1:2,])

	absCnDf3Sub <- absCnDf3 %>% dplyr::filter( chromNum == chr & l1 <= startNum & endNum <= l2)
	print(absCnDf3Sub[1:2,])

	# convert to long format to have a segment-sample pair per row and
	# corresponding cn:
	#cnSegLong <- absCnDf2Sub %>%
	cnSegLong <- absCnDf3Sub %>%
		dplyr::select(-chrom, -start, -end, -chromNum, -startNum, -endNum) %>%
		tidyr::pivot_longer(-coord, names_to = "SampleName", values_to = "cn")
	print(cnSegLong[1:2,])

	# split segment coord into chr, start and end
	df <- cnSegLong %>% tidyr::separate(coord, c("chrom", "start", "end"), sep = "[:-]") %>%
		mutate(chrom=as.numeric(chrom),
			start=as.numeric(start),
			end=as.numeric(end))
	# and compute mid point to use in plot of CN along region:
	df <- df %>% mutate(x = start + (end - start + 1))
	print(df[1:2,])
	
	# plot CN along region	
	p <- ggplot(df, aes(x=x, y=cn, group = SampleName, colour=SampleName)) +
		geom_line() +
		geom_point()
	p <- p + geom_vline(xintercept = start(cnSegRg[i,]))
	p <- p + geom_vline(xintercept = end(cnSegRg[i,]))
	tmpStart <- start(agRg[subjectHits(x),])
	tmpEnd <- end(agRg[subjectHits(x),])
	segDf <- data.frame(x1 = tmpStart, x2 = tmpEnd, y1 = -0.5, y2=-0.5, col="black")
	segDf$SampleName <- "gene"
	segDf$gene_name <- mcols(agRg[subjectHits(x),])$gene_name
	segDf$mid <- segDf$x1 + (segDf$x2 - segDf$x1)/2
	p <- p + geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2),
			      size = 3,
			      #arrow = arrow(length = unit(0.1, "inches")), #
			      #would need strand
			      data = segDf)

	p <- p + geom_label_repel(aes(x=mid, label=gene_name),
				  y=-0.5,
				  data = segDf,
				show.legend = FALSE)
	p <- p + xlab(sprintf("chr %s", chr))
	p <- p + ylab("absolute copy number")
	tmpTitle <- sprintf("%s:%s-%s (width: %d kb)",
				chr,
				start(cnSegRg[i,]),
				end(cnSegRg[i,]),
				width(cnSegRg[i,])/1000
				)
	p <- p + ggtitle(tmpTitle)
	p <- p + theme(plot.title = element_text(hjust = 0.5))
	#print(p)

	# write to file:
	tmpFn <- sprintf("%s_cn.png", tmpBit)
	print(sprintf("tmpFn: %s", tmpFn))
	ggsave(tmpFn, plot = p, width = 2*5, height = 1*5, units = c("in"))

	# for each gene
	for (iGene in 1:length(agRg[subjectHits(x),]))
	{
		# have GRange for that gene to find which segment it overlaps:
		tmpGeneGr <- agRg[subjectHits(x),][iGene,]
		geneStg <- sprintf("%s::%s",
				   mcols(tmpGeneGr)$gene_id,
				   mcols(tmpGeneGr)$gene_name)
		
		# find overlapping segments:
		y <- findOverlaps(tmpGeneGr, cnSegRgOrig)
		nbSegOver <- length(y)
		# The region of interest may not overlap any gene
		# and the gene in the vicinity may not overlap any of the
		# segments available
		# Then use the segments in the region of interest.
		if(length(y) == 0)
		{
			y <- findOverlaps(cnSegRg[i,], cnSegRgOrig)
			nbSegOver <- -length(y)
			print(subjectHits(y))
		}
		# retrieve cn for these overlapped segments:
		cnSegOver <- cnSegRgOrig[subjectHits(y),]
		# compute the average cn across these segments:
		aveCn <- apply(mcols(cnSegOver),2,mean)
		# have average for each sample in data frame to join with tpm
		aveCnDf <- data.frame(aveCn) %>%
			tibble::rownames_to_column("sample") %>%
			mutate(sample = gsub("^X", "", sample)) %>%
			dplyr::rename(segVal = aveCn)
		aveCnDf$Gene <- sprintf("%s::%s", tmpGeneGr$gene_id, tmpGeneGr$gene_name)
		# rename to use code written for Maria's gene of interest:
		geneDf <- aveCnDf
		print(head(geneDf))

		# merge with CN matrix, using {gene_name, SampleName} IDs.
		#-------------
		# tpm uses JBLAB IDs, eg JBLAB19902
		# cn uses another ID, eg X119025org
		# so first merge the cn table with the sample sheet
		geneDfOrig <- geneDf
		geneDfOrig$sample <- as.character(geneDfOrig$sample)

		geneDf <- geneDfOrig %>% mutate(organoid_name = sample) %>%
			left_join(splSht[, c("organoid_name", "SampleName", "OV04_number")], by="organoid_name")

		# now merge the cn and tpm matrices,
		# in long format, so with gene-sample pair
		geneDf2 <- geneDf %>%
		       	select(Gene, segVal, SampleName, OV04_number, organoid_name) %>%
			mutate(gene_stg = Gene) %>%
			mutate(SampleNameOrig = SampleName) %>%
			mutate(SampleName = organoid_name) %>%
			left_join(tpmMatLong, by=c("gene_stg", "SampleName")) %>%
			dplyr::select(-gene_stg, -organoid_name)
		print(head(geneDf2))

		#countColumn <- "tpm"
		countColumn <- "count"
		geneDf2.nbNa <- geneDf2 %>%
			group_by(Gene) %>%
			select(one_of(countColumn)) %>%
			summarise(nbNa = sum(is.na(.)))

		genesIn <- geneDf2.nbNa %>% filter(nbNa == 0) %>%
			pull(Gene)
	
		# skip if no data for that gene.
		if(length(genesIn)==0)
		{
			ctDf <- data.frame("GeneName"=mcols(tmpGeneGr)$gene_name,
				"EnsemblId"=mcols(tmpGeneGr)$gene_id, 
				"esti"=NA,
				"pval"=NA,
				"nbSegOver"=0,
				"plot"=NA)
			# add empty count and cn for sample
			countAndCn <- data.frame(t(matrix(rep(NA, 2*nrow(geneDf2)))))
			colnames(countAndCn) <- c(paste("segVal", geneDf2$SampleName, sep="_"),
							paste(countColumn, geneDf2$SampleName, sep="_"))
			ctDf <- cbind(ctDf, countAndCn)
			myls2[[geneStg]] <- ctDf
			next
		}

		# plot single gene
		library(ggplot2)
		tmpdf <- geneDf2 %>% filter(Gene == geneStg)
		# de-factor OV04_number 
		tmpdf$OV04_number <- as.character(tmpdf$OV04_number)
		tmpdf$expVal <- tmpdf[,countColumn]

		# correlation
		# linear regression - skip 
		# lm(geneDf2$tpm ~ geneDf2$segVal)

		print(geneDf2)
		ct <- cor.test(geneDf2$segVal, geneDf2[,countColumn], alternative = "greater",
			 method = "kendall", exact = FALSE)
		##ct$estimate
		##ct$p.value
		##ct$method

		# plot
		p <- ggplot(tmpdf, aes(x=segVal, y=expVal, label=SampleName))
		p <- p + geom_smooth(method = lm, se = TRUE)
		p <- p + geom_point(aes(colour=OV04_number)) # , size=0.5
		p <- p + geom_text_repel(aes(colour=OV04_number), show.legend = FALSE) # , size = 1
		p <- p + xlab("absolute copy number")
		p <- p + ylab(countColumn)
		# title: gene name and coord, cor estimate and p-value
		tmpCoord <- sprintf("%s:%s-%s",
				    as.numeric(as.vector(seqnames(tmpGeneGr))),
				    start(tmpGeneGr),
				    end(tmpGeneGr)
				)
		tmpTile <- sprintf("%s\n%s\ntau: %.2f, p-val: %.2g, nb seg: %d",
				  geneStg,
				  tmpCoord,
				  round(ct$estimate,2),
				  ct$p.value,
				nbSegOver)
		p <- p + ggtitle(tmpTile)
		p <- p + theme(plot.title = element_text(hjust = 0.5))
		# file
		geneStg2 <- sprintf("%s_%s",
				   mcols(tmpGeneGr)$gene_id,
				   mcols(tmpGeneGr)$gene_name)
		
		tmpFn <- sprintf("%s_%s.png", tmpBit, geneStg2)
		print(sprintf("tmpFn: %s", tmpFn))
		ggsave(tmpFn, plot = p, width = 2*5, height = 1*5, units = c("in"))

		# data.frame
		ctDf <- data.frame("GeneName"=mcols(tmpGeneGr)$gene_name,
				   "EnsemblId"=mcols(tmpGeneGr)$gene_id, 
				   "esti"=ct$estimate,
				"pval"=ct$p.value,
				"nbSegOver"=nbSegOver,
				"plot"=tmpFn)
		# add count and cn:
		head(tmpdf)
		tmpdfWide <- tmpdf %>%
			select(SampleName, segVal, all_of(countColumn)) %>%
			tidyr::pivot_wider(names_from = "SampleName",
					values_from = c(segVal, all_of(countColumn))
			)
		head(tmpdfWide)
		head(ctDf)
		ctDf <- cbind(ctDf, tmpdfWide)
		head(ctDf)
		# keep that in list:	
		myls2[[geneStg]] <- ctDf
	}
	print(myls2)
	lapply(myls2, dim)
	corDf <- do.call(rbind, myls2) # assumes tmpdfWide will always be the same
	corDf$regionName <- names(myls)[i]
	myls[[i]] <- corDf

} # end loop for region

#print(myls)

corDfAll <- do.call(rbind, myls)
corDfAll <- corDfAll %>% mutate(esti = round(esti,2),
				pval = sprintf("%.2g", pval))
tmpFn <- sprintf("%s.csv", "corStatTable")
write.csv(corDfAll, file=tmpFn)
