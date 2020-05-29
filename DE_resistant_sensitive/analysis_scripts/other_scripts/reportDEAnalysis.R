## copied from  /home/bioinformatics/software/pipelines/rnaseq/rnaseq-3.0.current/R/reportDEAnalysis.R 

suppressMessages(library(optparse))
suppressMessages(library(rnaseqRpkg))
suppressMessages(library(readr))

parseCommandLine <- function() {
  parser <- OptionParser(description="Report on a contrast, given a DE Object.")
  parser <- add_option(parser,"--model",help="model for DE analysis")
  parser <- add_option(parser,"--numerator",help="numerator for DE report")
  parser <- add_option(parser,"--denominator",help="denominator for DE report")
  parser <- add_option(parser,"--numeratorName",help="numerator for DE report")
  parser <- add_option(parser,"--denominatorName",help="denominator for DE report")
  parser <- add_option(parser,"--threshold",help="p-value significance cut-off")
  parser <- add_option(parser,"--refflat",help="Refflat-formatted gene list")
  parser <- add_option(parser,"--geneNames",help="Human-readable names of genes")
  parser <- add_option(parser,"--deaObject",help="DEA object to work from")
  parser <- add_option(parser,"--differential_all",help="Output file, all")
  parser <- add_option(parser,"--differential_sig",help="Output file, filtered for significance.")
  parser <- add_option(parser,"--summary_counts",help="Summary of differential counts")
  parser <- add_option(parser,"--sampleSheet",help="Sample Sheet")
  opts = parse_args(parser,positional_arguments=0)
  return(opts$options)
}

opts <- parseCommandLine()
threshold <- as.numeric(opts$threshold)
fnPrefix <- sprintf("%s_DEA_%s-vs-%s",opts$model,opts$numerator,opts$denominator)
targetDir <- dirname(opts$differential_all)
results <- reportDEAnalysis(opts$model,opts$numerator,opts$denominator,opts$refflat,opts$geneNames,opts$deaObject,threshold,opts$sampleSheet)

# use write.csv to get column labels (Don't tell Hadley.)
write.csv(results$report,file = opts$differential_all)
write.csv(results$filtered,file = opts$differential_sig)
reportDEAplots(results$samples,results$deaObj,results$report,results$filtered,targetDir,fnPrefix,opts$model,opts$numerator,opts$denominator,200)

annotations <- data.frame("Model" = opts$model, "Numerator" = opts$numerator, "Denominator" = opts$denominator)
summary <- cbind(annotations,results$stats)
print(summary)
write_csv(summary,opts$summary_counts)