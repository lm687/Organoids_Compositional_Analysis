## this file is copied from
## /home/bioinformatics/software/pipelines/rnaseq/rnaseq-3.0.current/R/runDEAnalysis.R 

suppressMessages(library(optparse))
suppressMessages(library(rnaseqRpkg))

parseCommandLine <- function() {
  parser <- OptionParser(description="Perform a differential binding analysis using DESeq2")
  parser <- add_option(parser,"--sampleSheet",
                       help="sample sheet")
  parser <- add_option(parser,"--model",
                       help="model for DE analysis (in text form)")
  parser <- add_option(parser,"--counts_raw",
                       help="table of raw counts")
  parser <- add_option(parser,"--deaObjectFile",
                       help="file to write DEA object to")
  opts = parse_args(parser,positional_arguments=0)
  return(opts$options)
}



loadCountTable <- function(fn) {
  data <- read.csv(fn,row.names=1,header=TRUE,check.names=FALSE)
  return(data)
}

loadRefFlat <- function(fn) {
  data <- read.table(fn,sep="\t",quote="",header=FALSE,
                     colClasses=c("character","NULL","character","character",
                                  "numeric","numeric","NULL","NULL","NULL",
                                  "NULL","NULL"))
  colnames(data) <- c('Gene','Chrom','Strand','Left','Right')
  return(data)
}

loadGeneNames <- function(fn) {
  data <- read.table(fn,sep="\t",quote="",header=FALSE)
  colnames(data) <- c('Gene','Symbol','Name')
  return(data)
}

annotateGenes <- function(mat,refflat,geneNames) {
  refflat_ordered <- refflat[match(rownames(mat),refflat$Gene),]
  geneNames_ordered <- geneNames[match(rownames(mat),geneNames$Gene),]
  mat$chrom <- refflat_ordered$Chrom
  mat$strand <- refflat_ordered$Strand
  mat$left <- refflat_ordered$Left
  mat$right <- refflat_ordered$Right
  mat$symbol <- geneNames_ordered$Symbol
  mat$name <- geneNames_ordered$Name
  return(mat)
}

makeVolcanoPlot <- function(results,target_dir,prefix) {
  fn <- file.path(target_dir,sprintf("%s_VolcanoPlot.pdf",prefix))
  pdf(fn)
  print(plotVolcanogg(results))
  dev.off()
}

makeMAplot <- function(results,target_dir,prefix) {
  fn <- file.path(target_dir,sprintf("%s_MAplot.pdf",prefix))
  pdf(fn)
  print(plotMAgg(results))
  dev.off()
}
#makeMAplot <- function(results,target_dir,prefix,cutoffs) {
#  fn <- file.path(target_dir,sprintf("%s_MAplot.pdf",prefix))
#  pdf(fn)
#  xdata <- log2(results$baseMean)
#  ydata <- results$log2FoldChange
#  sig <- results$padj
#  smoothScatter(xdata,ydata,main=sprintf("%s MA Plot",target_dir),xlab="log2baseMean",
#                ylab="log2FC")
#  abline(h=0,col="grey",lwd=2)
#  cols <- rainbow(n=length(cutoffs),start=2/6,end=0)
#  for(i in 1:length(cutoffs)) {
#    cutoff <- cutoffs[i]
#    points(xdata[sig<=cutoff],ydata[sig<=cutoff],col=cols[i],pch=19,cex=0.30)
#  }
#  dev.off()
#}

makeHeatmap <- function(deMat,report,count,target_dir,prefix,all_samples,target,num,denom) {
  if (all_samples) {
    label <- "allSamples"
  } else {
    label <- "testSamples"
    wanted <- colData(deMat)[,target] %in% c(num,denom)
    deMat <- deMat[,wanted]
  }
  fn <- file.path(target_dir,sprintf("%s_top200_heatmap_%s.pdf",
                                     prefix,label))
  
  pdf(fn)
  print(plotTopNSigHeatmap(report,deMat,target,topN=count))
  dev.off()
}

makePCA <- function(deMat,model,report,target_dir,prefix,target,numerator,denominator,samples) {
  fn <- file.path(target_dir,sprintf("%s_PCA.pdf",prefix))
  pdf(fn)
  print(plotPCAbyContrast(report,target,deMat,samples))
  dev.off()
}

getMean <- function(data,factor,values,sample) {
  cols <- sample[,factor] %in% values
  data <- data[,cols,drop=FALSE]
  sums <- apply(data,c(1),sum)
  return(sums/sum(cols))
}

getMeans <- function(mat,factor,numerator,denominator,samples) {
  data <- counts(mat,normalized=TRUE)
  numData <- getMean(data,factor,numerator,samples)
  denData <- getMean(data,factor,denominator,samples)
  allData <- getMean(data,factor,c(numerator,denominator),samples)
  df <- data.frame(allData,numData,denData)
  colnames(df) <- c("Both_Mean",sprintf("%s_Mean",numerator),sprintf("%s_Mean",denominator))
  return(df)
}

addMeans <- function(mat,res,factor,numerator,denominator,samples) {
  means <- getMeans(mat,factor,numerator,denominator,samples)
  means$FC = 2^res$log2FoldChange
  res <- cbind(res,means)
  res <- res[c(1,7,8,9,2,10,3,4,5,6)]
  return(res)
}

writeReport <- function(deMat,model,target,denominator,numerator,refflat,geneNames,destDir,subDir,fnPrefix,samples) {
  target_dir <- file.path(destDir,subDir)
  fnBase <- sprintf("%s_DEA_%s-vs-%s",fnPrefix,numerator,denominator)
  if (!dir.exists(target_dir)) {
    dir.create(target_dir)
  }
  report <- results(deMat,c(target,numerator,denominator),alpha=0.05,format='DataFrame')
  report_df <- as.data.frame(addMeans(deMat,report,target,numerator,denominator,samples))
  for (column in (1:6)) {
    report_df[,column] <- round(report_df[,column],3)
  }
  report_df <- annotateGenes(report_df,refflat,geneNames)
  report_df$location <- with(report_df,paste(chrom,":",left,"-",right,sep=""))
  fn_all <- file.path(destDir,subDir,sprintf("%s_all.csv",fnBase))
  write.csv(report_df,file=fn_all)
  fn_filt <- file.path(destDir,subDir,sprintf("%s_filtered.csv",fnBase))
  report_filt <- report_df[!is.na(report_df$padj)&report_df$padj<0.05,]
  write.csv(report_filt,file=fn_filt)
  makeMAplot(report,target_dir,fnBase)
  makeVolcanoPlot(report,target_dir,fnBase)
  if (nrow(report_filt) > 0) {
    makeHeatmap(deMat,report_filt,200,target_dir,fnBase,all_samples=TRUE,target,numerator,denominator)
    makeHeatmap(deMat,report_filt,200,target_dir,fnBase,all_samples=FALSE,target,numerator,denominator)
    makePCA(deMat,model,report,target_dir,fnBase,target,numerator,denominator,samples)
  }
  numeratorCount <- sum(samples[[target]] == numerator)
  denominatorCount <- sum(samples[[target]] == denominator)
  upcount <- sum(report_filt$log2FoldChange > 0)
  downcount <- sum(report_filt$log2FoldChange < 0)
  return(c(upcount,downcount,numeratorCount,denominatorCount))
}

runDEanalysis <- function(model,samples,counts,refflat,geneNames,outputDir) {
  #  print(sprintf("outputDir: %s",outputDir))
  fsamples <- sampleSheet2Factors(samples)
  print(fsamples)
  head(counts)
  print(model)
  mat <- suppressMessages(DESeqDataSetFromMatrix(countData=counts,colData=fsamples,design=model))
  mat <- estimateSizeFactors(mat)
  mat <- suppressMessages(estimateDispersions(mat,fitType="local"))
  tryCatch( mat <- nbinomWaldTest(mat),
            error = function(e) {message("failed in nbinom")})
  terms <- rev(labels(terms(model)))
  fnPrefix <- paste(terms,sep="-",collapse="-")
  target <- terms[1]
  values <- levels(fsamples[,target])
  nval <- length(values)
  summaries <- list()
  for (i in 1:(nval-1)) {
    for (j in (i+1):nval) {
      subdir <- sprintf("%s_DEA_%s_vs_%s",fnPrefix,values[j],values[i])
      #      print(fnPrefix)
      #      print(values[j])
      #      print(values[i])
      #      print(sprintf("subDir: %s",subdir))
      upAndDown <- writeReport(mat,model,target,values[i],values[j],refflat,geneNames,outputDir,subdir,fnPrefix,samples)
      res <- c('Model'=as.character(model),'Numerator'=values[j],'Denominator'=values[i],'Up-regulated'=upAndDown[1],'Down-regulated'=upAndDown[2],'NumCount'=upAndDown[3],'DenomCount'=upAndDown[4])
      summaries <- c(summaries,list(res))
    }
  }
  return(list(mat,summaries))
}


opts <- parseCommandLine()
samplesFN <- opts$sampleSheet
modelStr <- opts$model
countsFN <- opts$counts_raw
destination <- opts$deaObjectFile

deaObj <- runDEAnalysis(sampleSheetFN = samplesFN, modelString = modelStr, countsFN=countsFN)
assign(modelStr,deaObj)
save(list=modelStr,file=destination)
