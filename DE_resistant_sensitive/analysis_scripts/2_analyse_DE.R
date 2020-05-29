
rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(rnaseqRpkg)
library(biomaRt)

## Load file output from DESeq2
load("../objects/deaObjectFile")
deObj = `~response`
results <- results(deObj, c("response", "complete_remission_or_response", "progressive_disease"), 
                  alpha = 0.05, format = "DataFrame")
rownames_short = sapply(rownames(results), function(i) strsplit(i, '[.]')[[1]][1])

## Re-name
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
gene_conversion <- getBM(attributes = c("external_gene_name", "ensembl_gene_id"),
                 filters = "ensembl_gene_id", values = rownames_short,
                 mart = mart)
gene_conversion = gene_conversion[match(rownames_short, gene_conversion$ensembl_gene_id),]
dim(gene_conversion)
dim(results)


results = cbind(gene_conversion, results)
results = results[order(results$log2FoldChange),]

highDA = results[1:20,]
lowDA = highDA[nrow(highDA):(nrow(highDA)-20),]

ggplot(cbind.data.frame(log2fc= results$log2FoldChange, pval=-log10(results$padj)))+
  geom_point(aes(x=log2fc, y=pval))

vlcano = function (res) {
  df <- data.frame(labls = res$external_gene_name, lfc = res$log2FoldChange, pv = -log10(res$padj), 
                   isDE = ifelse(is.na(res$padj), FALSE, res$padj < 0.05))
  colNonSig <- "gray32"
  colSig <- "red3"
  colOut <- "blue"
  dotShape <- 20
  upShape <- 24
  outShape <- 19
  xlim <- c(-1, 1) * quantile(abs(df$lfc[is.finite(df$lfc)]), 
                              probs = 0.995) * 1.1
  ylim <- c(0, quantile(df$pv[is.finite(df$pv)], probs = 0.995) * 
              1.1)
  df$vy <- pmin(df$pv, ylim[2])
  df$vx <- pmax(xlim[1], pmin(xlim[2], df$lfc))
  df$outY <- df$pv > ylim[2]
  df$outX <- df$lfc < xlim[1] | df$lfc > xlim[2]
  df$outXY <- df$outX & df$outY
  dfNonSig <- df[!df$isDE, ]
  dfNonSig$vy[is.na(dfNonSig$vy)] = 0
  dfSigInlier <- df[df$isDE & !df$outX & !df$outY, ]
  dfSigOutlier <- df[df$isDE & (df$outX | df$outY), ]
  dfSigOutlier$shape <- ifelse(dfSigOutlier$outXY, outShape, 
                               ifelse(dfSigOutlier$outY, upShape, dotShape))
    p <- ggplot() +
      geom_point(data = dfNonSig, aes(x = .data$vx,   y = .data$vy),
                 colour = colNonSig, pch = dotShape) + 
      geom_point(data = dfSigInlier, aes(x = .data$vx, y = .data$vy), 
               colour = colSig, shape = dotShape) +
      ggrepel::geom_label_repel(data=dfSigInlier, aes(x=  .data$vx, y = .data$vy, label=labls))+
      geom_point(data = dfSigOutlier,  aes(x = .data$vx, y = .data$vy),
               colour = colOut, shape = dfSigOutlier$shape) + 
      ggrepel::geom_label_repel(data=dfSigOutlier, aes(x=  .data$vx, y = .data$vy, label=labls))+
      theme(legend.position = "none") + ggtitle("Volcano Plot") + 
      labs(x = "Log of Fold Change", y = "-log10 of FDR")+
      ylim(c(0, 1.8))
  return(p)
}
vlcano(results)
ggsave("../objects/volcano_plot.pdf", width = 20, height = 20)

