# copied from /home/bioinformatics/software/pipelines/rnaseq/rnaseq-3.0.current/R/utilities.R
# Various utility functions for plotting, loading data and so on.
suppressMessages(library(gplots))
suppressMessages(library(DESeq2))
suppressMessages(library(ggplot2))
suppressMessages(library(ggfortify))
suppressMessages(library(genefilter))
suppressMessages(library(formula.tools))
suppressMessages(library(hexbin))
suppressMessages(library(RColorBrewer))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(dplyr))

###############################################################################
#
# Constants we need here and there
#

# columns of sample sheet that should not be interpreted as factors
NON_FACTOR_COLUMNS = c("SampleName","Barcode","SLXId","FileName")

###############################################################################
#
# Miscellaneous Text Functions
# 

truncateStr <- function(strs,n) {
  locs <- !is.na(strs) & nchar(strs) > n
  strs[locs] = paste(substr(strs[locs],1,n),"...",sep="")
  return(strs)
}

toPercent <- function(n) {
  return(sprintf("%.1f%%",n*100))
}

toShorter <- function(n) {
  return(sprintf("%.2f",n))
}

friendlyNum <- function(n) {
  return(gsub("b$","",format(structure(n, class="object_size"), units="auto")))
}

###############################################################################
#
# Functions for loading counts, sample sheets, etc.
#

loadCounts <- function(fn) {
  data <- read.csv(fn,header=TRUE,row.names=1,check.names=FALSE)
  return(data)
}

loadCountStats <- function(fn) {
  data <- read.csv(fn,header=TRUE,check.names=FALSE,stringsAsFactors=FALSE)
  data$Sample <- as.character(data$Sample)
  data <- data[c(14,1:13,15)] # put sample name first
  data$Sample = as.character(data$Sample)
  return(data)
}

loadCoverageData <- function(fn) {
  df <- read.csv(fn,colClasses=c(NA,NA,"character"))
  return(df)
}

loadDEobjects <- function(fnstr) {
  fns <- strsplit(fnstr,",",fixed=TRUE)[[1]]
  data <- list()
  for (fn in fns) {
    names <- load(fn)
    for (name in names) {
      data[[name]] <- get(name)
    }
  }
  return(data)
}

loadSampleSheet <- function(fn) {
  data <- read.csv(fn,header=TRUE,stringsAsFactors=FALSE)
  data$SampleName = as.character(data$SampleName)
  return(data)
}

sampleSheet2Factors <- function(sheet) {
  for (name in colnames(sheet)) {
    if (!(name %in% NON_FACTOR_COLUMNS)) {
      sheet[,name] <- factor(sheet[,name],levels=unique(sheet[,name]))
    }
  }
  return(sheet)
}

loadModels <- function(modelfile) {
  data <- readLines(con=modelfile)
  data <- data[!grepl("^\\s*$",data)]
  models <- lapply(data,as.formula)
  return(models)
}

###############################################################################
#
# Filtering functions
#

getTopNByPadj <- function(res,n) {
  data <- as.data.frame(res)
  data$sample <- rownames(data)
  data <- arrange(data,padj)
  names <- dplyr::slice(data,1:n)$sample
  return(names)
}

getTopNByVar <- function(counts,n) {
  rv <- rowVars(counts)
  sel <- order(rv,decreasing=TRUE)[seq_len(min(n,length(rv)))]
  return(sel)
}

###############################################################################
#
# Plotting functions
#

plotPCAgg <- function(factor,data,samples,topn=500) {
  main <- sprintf("Colour: %s",factor)
  sel <- getTopNByVar(data,topn)
  pca <- prcomp(t(data[sel,]),scale.=TRUE)
  summ <- summary(pca)$importance
  xlab <- tryCatch(sprintf("PC1: %.1f%% variance",summ[2,1]*100),
                   error = function(err) { sprintf("PC1: unknown variance")})
  ylab <- tryCatch(sprintf("PC2: %.1f%% variance",summ[2,2]*100),
                   error = function(err) { sprintf("PC2: unknown variance")})
  samples[,factor] <- as.character(samples[,factor])
  return(autoplot(pca,
                  data=samples,
                  colour=factor,
                  main=main))
}

plotPCAbyContrast <- function(res,factor,data,samples,small=FALSE) {
  main <- sprintf("PCA Plot: by %s",factor)
  res <- res[!is.na(res$padj)&res$padj<0.05,]
  if (nrow(res) < 1) {
    return(NULL)
  }
  keepers <- rownames(res)
  rdata <- data[keepers,]
  vstdata <- tryCatch(varianceStabilizingTransformation(rdata,
                                                        blind=FALSE,
                                                        fitType="local"),
                      error=function(err) { 
                        return(rlog(rdata,blind=FALSE,fitType="local"))
                      })
  pca <- prcomp(t(assay(vstdata)),scale.=TRUE)
  summ <- summary(pca)$importance
  xlab <- tryCatch(sprintf("PC1: %.1f%% variance",summ[2,1]*100),
                   error = function(err) { sprintf("PC1: unknown variance")})
  ylab <- tryCatch(sprintf("PC2: %.1f%% variance",summ[2,2]*100),
                   error = function(err) { sprintf("PC2: unknown variance")})
  samples[,factor] <- as.character(samples[,factor])
  p = tryCatch({
    if (small) {
      p <- autoplot(pca,
                    data=samples,
                    colour=factor)
    } else {
      p <- autoplot(pca,
                    data=samples,
                    colour=factor,
                    main=main)
    }
  }, error=function(err) { return(NULL) })
  return(p)
}

plotMAhex <- function(res,small=FALSE,bins=75) {
  df <- data.frame(mean=res$baseMean,
                   lfc=res$log2FoldChange,
                   isDE=ifelse(is.na(res$padj),FALSE,res$padj<0.05))
  colNonSig <- "gray32"
  colSig <- "red3"
  colOut <- "blue"
  dotShape <- 20
  downShape <- 25
  upShape <- 24
  countThreshold <- 150
  
  ylim = c(-1,1) * quantile(abs(df$lfc[is.finite(df$lfc)]), probs=0.995) * 1.1
  df$py <- pmax(ylim[1], pmin(ylim[2], df$lfc))
  df$outlier <- df$lfc < ylim[1] | df$lfc > ylim[2]
  
  dfNonSig <- df[!df$isDE,]
  dfSig <- df[df$isDE,]
  dfSigInlier <- dfSig[!dfSig$outlier,]
  dfSigOutlier <- dfSig[dfSig$outlier,]
  dfSigOutlier$pch <- ifelse(dfSigOutlier$lfc<0, downShape, upShape)
  
  hexNonSig <- hexbin(log10(dfNonSig$mean),dfNonSig$py,xbins=bins)
  dataNonSig <- data.frame(mean=10^(hexNonSig@xcm),py=hexNonSig@ycm)
  
  if (nrow(dfSigInlier) > countThreshold) {
    hexSig <- hexbin(log10(dfSigInlier$mean),dfSigInlier$py,xbins=bins)
    dataSig <- data.frame(mean=10^hexSig@xcm,py=hexSig@ycm)
  } else {
    dataSig <- dfSigInlier
  }
  
  if (small) {
    p <- ggplot(dataNonSig,aes(x=mean,y=py)) +
      geom_point(colour=colNonSig,pch=dotShape) +
      geom_point(data=dataSig,aes(x=mean,y=py),colour=colSig,shape=dotShape) +
      geom_point(data=dfSigOutlier,aes(x=mean,y=py),colour=colOut,shape=dfSigOutlier$pch) +
      scale_x_log10() +
      theme(legend.position="none") +
      labs(x="Mean of Normalised Counts",y="Log of Fold Change")
  } else {
    p <- ggplot(dataNonSig,aes(x=mean,y=py)) +
      geom_point(colour=colNonSig,pch=dotShape) +
      geom_point(data=dataSig,aes(x=mean,y=py),colour=colSig,shape=dotShape) +
      geom_point(data=dfSigOutlier,aes(x=mean,y=py),colour=colOut,shape=dfSigOutlier$pch) +
      scale_x_log10() +
      theme(legend.position="none") +
      labs(x="Mean of Normalised Counts",y="Log of Fold Change")
  }
  return(p)
}


plotMAgg <- function(res) {
  df <- data.frame(mean=res$baseMean,
                   lfc=res$log2FoldChange,
                   isDE=ifelse(is.na(res$padj),FALSE,res$padj<0.05))
  colNonSig <- "gray32"
  colSig <- "red3"
  colOut <- "blue"
  ylim = c(-1,1) * quantile(abs(df$lfc[is.finite(df$lfc)]), probs=0.995) * 1.1
  df$outlier <- df$lfc < ylim[1] | df$lfc > ylim[2]
  py <- pmax(ylim[1], pmin(ylim[2], df$lfc))
  pch <- ifelse(df$lfc<ylim[1], 25, ifelse(df$lfc>ylim[2], 24, 20))
  col <- ifelse(df$outlier,colOut,ifelse(df$isDE, colSig, colNonSig))
  df$py <- py
  p <- ggplot(df,aes(x=mean,y=py)) +
    geom_point(colour=col,shape=pch) +
    scale_x_log10() +
    labs(x="Mean of Normalised Counts",y="Log of Fold Change")
  return(p)
}

plotClusteringHeatMap <- function(counts,samples) {
  mat <- as.matrix(dist(t(counts)))
  rownames(mat) <- samples$SampleName
  colnames(mat) <- samples$SampleName
  heatmap.2(mat,
            margins=c(14,14),
            trace="none",
            key.xlab="Distance")
}

plotClusteringHeatmap1 <- function(counts,samples,targets,small=FALSE) {
  d <- as.matrix(dist(t(counts)))
  colours <- mkColourMaps(samples,targets)
  df <- data.frame(samples[,targets])
  colnames(df) <- targets
  ha <- HeatmapAnnotation(df=df,col=colours,show_legend=FALSE,show_annotation_name=TRUE)
  if (small) {
    hm <- Heatmap(d,name="Distance",
                  top_annotation=ha,
                  show_row_names=FALSE,show_column_names=FALSE)
  } else {
    hm <- Heatmap(d,name="Distance",
                  top_annotation=ha,show_row_names=FALSE,show_column_names=TRUE,
                  column_title="Sample Clustering",
                  column_title_side="top",
                  column_names_gp = gpar(fontsize=10))
  }
  hm <- hm + rowAnnotation(df=df,col=colours)
  return(hm)
}


plotClusteringHeatmap2 <- function(res,mat,target,samples,small=FALSE) {
  keepers <- rownames(res)[!is.na(res$padj)&res$padj < 0.05]
  kmat <- assay(mat)[keepers,]
  d <- as.matrix(dist(t(kmat)))
  #  categories <- unique(samples[,target])
  #  ccount <- length(categories)
  #  cols <- brewer.pal(ccount,"RdBu")
  #  names(cols) <- categories
  #  colours <- list(type=cols)
  #  names(colours) <- c(target)
  colours <- mkColourMaps(samples,c(target))
  df <- data.frame(samples[,target])
  colnames(df) <- c(target)
  ha <- HeatmapAnnotation(df=df,col=colours,show_legend=FALSE,show_annotation_name=TRUE)
  if (small) {
    hm <- Heatmap(d,name="Distance",
                  top_annotation=ha,
                  show_row_names=FALSE,show_column_names=FALSE)
  } else {
    hm <- Heatmap(d,name="Distance",
                  top_annotation=ha,show_row_names=FALSE,show_column_names=TRUE,
                  column_title="Sample Clustering",
                  column_title_side="top",
                  column_names_gp = gpar(fontsize=10))
  }
  hm <- hm + rowAnnotation(df=df,col=colours)
  return(hm)
}

plotTopNHeatMap <- function(counts,samples,colour_factors,topN=200) {
  sel <- getTopNByVar(counts,topN)
  data <- as.matrix(counts[sel,])
  
  # set up row, column annotations
  df <- data.frame(row.names=samples$SampleName)
  for (f in colour_factors) {
    df[,f] = samples[,f]
  }
  colours <- mkColourMaps(samples,colour_factors)
  ha <- HeatmapAnnotation(df=df,col=colours,show_annotation_name=TRUE)
  title <- sprintf("Heat Map: Top %d Most Variable Genes",topN)
  p <- Heatmap(data,name="log(count)",top_annotation=ha,show_row_names=FALSE,
               show_column_names=TRUE,column_title=title,
               column_title_side="top",
               column_names_gp=gpar(fontsize=10))
  return(p)
}

mkColourMaps <- function(data,targets,palette="RdBu") {
  colours <- list()
  for (target in targets) {
    categories <- unique(data[,target])
    ccount <- length(categories)
    if (ccount == 1) {
      cols <- c("#0000FF")
    } else if (ccount == 2) {
      cols <- c("#FF0000","#0000FF")
    } else if (ccount > 11) {
      cmap <- colorRampPalette(brewer.pal(11,palette))
      cols <- cmap(ccount)
    } else {
      cols <- brewer.pal(ccount,palette)
    }
    names(cols) <- categories
    colours <- c(colours,list(cols))
  }
  names(colours) <- targets
  return(colours)
}

plotTopNSigHeatmap <- function(res,mat,target,topN=200) {
  names <- getTopNByPadj(res,topN)
  keepers <- mat[names]
  data <- assay(keepers)
  data <- as.matrix(data)
  data[data==0] = 1 # so log2 doesn't give -inf
  data <- log2(data)
  samples <- colData(mat)
  ha_df <- data.frame(samples[,target])
  colnames(ha_df) <- c(target)
  colours <- mkColourMaps(samples,c(target))
  ha <- HeatmapAnnotation(df=ha_df,col=colours,
                          show_legend=TRUE,show_annotation_name=TRUE,
                          annotation_name_gp = gpar(fontsize=10))
  p <- Heatmap(data,
               show_row_names=FALSE,
               show_column_names=TRUE,
               name="log(Read Count)",
               top_annotation=ha,
               column_title=sprintf("Top %s Most DE Genes",topN),
               column_title_side="top",
               column_names_gp = gpar(fontsize=10))
  return(p)
}

plotVolcano <- function(res,small=FALSE,bins=75) {
  df <- data.frame(lfc=res$log2FoldChange,
                   pv=-log10(res$padj),
                   isDE=ifelse(is.na(res$padj),FALSE,res$padj<0.05))
  colNonSig <- "gray32"
  colSig <- "red3"
  colOut <- "blue"
  dotShape <- 20
  upShape <- 24
  outShape <- 19
  
  xlim <- c(-1,1) * quantile(abs(df$lfc[is.finite(df$lfc)]), probs=0.995) * 1.1
  ylim <- c(0,quantile(df$pv[is.finite(df$pv)],probs=0.995)*1.1)
  df$vy <- pmin(df$pv,ylim[2])
  df$vx <- pmax(xlim[1],pmin(xlim[2],df$lfc))
  df$outY <- df$pv > ylim[2]
  df$outX <- df$lfc < xlim[1] | df$lfc > xlim[2]
  df$outXY <- df$outX & df$outY
  
  dfNonSig <- df[!df$isDE,]
  dfNonSig$vy[is.na(dfNonSig$vy)] = 0
  # need to catch errors: in case all y values are 0, hexbin fails
  dataNonSig <- tryCatch({
    hexNonSig <- hexbin(dfNonSig$vx,dfNonSig$vy,xbins=bins)
    data.frame(vx=hexNonSig@xcm,vy=hexNonSig@ycm)
  }, error = function(e) { return(dfNonSig) })
  
  dfSigInlier <- df[df$isDE & !df$outX & !df$outY,]
  if (length(df$SigInlier$vx) > 0) {
    hexSig <- hexbin(dfSigInlier$vx,dfSigInlier$vy,xbins=bins)
    dataSig <- data.frame(vx=hexSig@xcm,vy=hexSig@ycm)
  } else {
    dataSig <- dfSigInlier
  }
  
  dfSigOutlier <- df[df$isDE & (df$outX | df$outY),]
  dfSigOutlier$shape <- ifelse(dfSigOutlier$outXY,outShape,ifelse(dfSigOutlier$outY,upShape,dotShape))
  
  if (small) {
    p <- ggplot() +
      geom_point(data=dataNonSig,aes(x=vx,y=vy),colour=colNonSig,pch=dotShape) +
      geom_point(data=dataSig,aes(x=vx,y=vy),colour=colSig,shape=dotShape) +
      geom_point(data=dfSigOutlier,aes(x=vx,y=vy),colour=colOut,shape=dfSigOutlier$shape) +
      theme(legend.position="none") +
      labs(x="Log of Fold Change",y="-log10 of FDR")
  } else {
    p <- ggplot() +
      geom_point(data=dataNonSig,aes(x=vx,y=vy),colour=colNonSig,pch=dotShape) +
      geom_point(data=dataSig,aes(x=vx,y=vy),colour=colSig,shape=dotShape) +
      geom_point(data=dfSigOutlier,aes(x=vx,y=vy),colour=colOut,shape=dfSigOutlier$shape) +
      theme(legend.position="none") +
      ggtitle("Volcano Plot") +
      labs(x="Log of Fold Change",y="-log10 of FDR")
  }
  
  return(p)
}

plotVolcanogg <- function(res) {
  df <- data.frame(lfc=res$log2FoldChange,
                   pv=-log10(res$padj),
                   isDE=ifelse(is.na(res$padj),FALSE,res$padj<0.05))
  colNonSig <- "gray32"
  colSig <- "red3"
  colOut <- "blue"
  dotShape <- 20
  upShape <- 24
  outShape <- 19
  
  xlim <- c(-1,1) * quantile(abs(df$lfc[is.finite(df$lfc)]), probs=0.995) * 1.1
  ylim <- c(0,quantile(df$pv[is.finite(df$pv)],probs=0.995)*1.1)
  df$vy <- pmin(df$pv,ylim[2])
  df$vx <- pmax(xlim[1],pmin(xlim[2],df$lfc))
  df$outY <- df$pv > ylim[2]
  df$outX <- df$lfc < xlim[1] | df$lfc > xlim[2]
  df$outXY <- df$outX & df$outY
  
  dfNonSig <- df[!df$isDE,]
  dfNonSig$vy[is.na(dfNonSig$vy)] = 0
  
  dfSigInlier <- df[df$isDE & !df$outX & !df$outY,]
  
  dfSigOutlier <- df[df$isDE & (df$outX | df$outY),]
  dfSigOutlier$shape <- ifelse(dfSigOutlier$outXY,outShape,ifelse(dfSigOutlier$outY,upShape,dotShape))
  
  p <- ggplot() +
    geom_point(data=dfNonSig,aes(x=vx,y=vy),colour=colNonSig,pch=dotShape) +
    geom_point(data=dfSigInlier,aes(x=vx,y=vy),colour=colSig,shape=dotShape) +
    geom_point(data=dfSigOutlier,aes(x=vx,y=vy),colour=colOut,shape=dfSigOutlier$shape) +
    theme(legend.position="none") +
    ggtitle("Volcano Plot") +
    labs(x="Log of Fold Change",y="-log10 of FDR")
  
  return(p)
}

###############################################################################
#
# Prepare tables for "top N" lists
#
dumpTopN <- function(model,numerator,denominator,keepers) {
  tag <- paste(rev(labels(terms(model))),sep="-",collapse="-")
  fn <- sprintf("%s/%s_DEA_%s_vs_%s/%s_DEA_%s-vs-%s_filtered.csv",
                opts$deaDir,
                tag,numerator,denominator,tag,numerator,denominator)
  df <- read.csv(fn,header=TRUE,row.names=1)
  real_keepers = keepers[keepers %in% rownames(df)]
  ndf <- df[real_keepers,c(3,4,5,9,15,16)]
  names <- colnames(ndf)
  colnames(ndf) <- c(names[1],names[2],"logFC","padj","symbol","name")
  ndf$name <- truncateStr(ndf$name,35)
  return(kable(ndf,caption="Top DE Sites",format="pandoc"))
}
