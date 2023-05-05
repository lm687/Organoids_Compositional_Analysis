give_CN_per_gene <- function(segment_arg, gr_genes=gr_genes){
  ## in this case GR_bins is gr_genes
  
  GR_bins <- gr_genes
  gr_CN = as(data.frame(segment_arg), "GRanges")
  GR_bins <- trim(GR_bins)
  gr_CN <- trim(gr_CN)
  
  values(gr_CN) = segment_arg[,'segVal']
  
  ## do per sample
  gr_CN_org = gr_CN
  seqlevels(gr_CN) ## necessary for mcolAsRleList
  seqlengths(gr_CN) = as.numeric(as.vector(chromlens$Length[match(seqlevels(gr_CN), gsub("chr", "", chromlens$Chrom))])) ## necessary for mcolAsRleList
  sub("_*", "", seqlevels(GR_bins)) %in% seqlevels(gr_CN) ##20220526
  seqlevels(GR_bins) = seqlevels(gr_CN) ## necessary for mcolAsRleList
  seqlengths(GR_bins) = seqlengths(gr_CN) ## necessary for mcolAsRleList
  
  # GR_bins$X = 
  GR_bins$segVal = 0
  # gr_CN_org$X <- NULL
  full_GR <- c(GR_bins, gr_CN_org)
  disjoint_gr <- GenomicRanges::disjoin(full_GR, with.revmap=TRUE, ignore.strand=TRUE)
  disjoint_gr <- trim(disjoint_gr)
  
  disjoint_gr_revmap_first = sapply(disjoint_gr$revmap, function(i) i[1])
  
  ## keep only those in which there is an overlap with a gene, i.e. revmap has to include the idxs 1:length(gr_genes)
  disjoint_gr = disjoint_gr[disjoint_gr_revmap_first %in% 1:length(GR_bins),]
  disjoint_gr_revmap_first = sapply(disjoint_gr$revmap, function(i) i[1])
  disjoint_gr_revmap_second = sapply(disjoint_gr$revmap, function(i) i[2])
  table(is.na(disjoint_gr_revmap_second))
  
  ## remove sections which only contain non-genes
  disjoint_gr[is.na(disjoint_gr_revmap_second),]
  
  ## get the idx of the copy number segments
  idx_gr = disjoint_gr_revmap_second-length(GR_bins)
  is.na(idx_gr)
  
  idx_gr[idx_gr <= 0 ] = NA
  
  disjoint_gr$CN_val = gr_CN_org$X[idx_gr]
  disjoint_gr$CN_val[is.na(disjoint_gr$CN_val)] = 2
  
  disjoint_gr$revmap = disjoint_gr_revmap_first
  seqlevels(disjoint_gr) = as.character(unique(seqnames(disjoint_gr)))
  seqlengths(disjoint_gr) = as.numeric(as.vector(chromlens$Length[match(seqlevels(disjoint_gr), gsub("chr", "", chromlens$Chrom))]))
  seqlevels(GR_bins) = seqlevels(disjoint_gr)
  seqlengths(GR_bins) = seqlengths(disjoint_gr)
  
  CN_RleList = mcolAsRleList(disjoint_gr, "CN_val") ## that takes forever even when length(disjoint_gr) is a VERY low value
  averageCN = binnedAverage(bins = GR_bins, numvar = CN_RleList, varname = "CN_bin_averaged")
  
  return(averageCN)
}

give_CN_per_gene_v2 <- function(segment_arg, gr_genes=gr_genes, transform_segment_dataframe=T, only_compute_overlaps=F, name_val_col='segVal', add_columns=F, value_is_CN=T){
  ## in this case GR_bins is gr_genes
  
  GR_bins <- gr_genes
  if(transform_segment_dataframe){
    gr_CN = as(data.frame(segment_arg), "GRanges")
  }else{
    gr_CN <- segment_arg
  }
  GR_bins <- trim(GR_bins)
  gr_CN <- trim(gr_CN)
  
  values(gr_CN) = segment_arg[,name_val_col]
  
  ## do per sample
  gr_CN_org = gr_CN
  seqlevels(gr_CN) ## necessary for mcolAsRleList
  seqlengths(gr_CN) = as.numeric(as.vector(chromlens$Length[match(seqlevels(gr_CN), gsub("chr", "", chromlens$Chrom))])) ## necessary for mcolAsRleList
  seqlevels(GR_bins) = seqlevels(gr_CN) ## necessary for mcolAsRleList
  seqlengths(GR_bins) = seqlengths(gr_CN) ## necessary for mcolAsRleList
  
  gr_CN_orgraw <- gr_CN_org
  gr_CN_orgraw$X <- NULL ## 20220526
  gr_CN_orgraw@elementMetadata[,name_val_col] <- NULL ## 20220526
  full_GR <- c(GR_bins, gr_CN_orgraw)
  
  if(add_columns){
    ## otherwise there is no match between the two, because they don't have the same metadata
    gr_CN_orgraw$names <- NA
    gr_CN_orgraw$segVal <- NA
  }
  disjoint_gr <- GenomicRanges::disjoin(full_GR, with.revmap=TRUE, ignore.strand=TRUE)
  print(length(disjoint_gr))
  # disjoint_gr <- trim(disjoint_gr) ## commented out 20220527
  
  disjoint_gr_revmap_first = sapply(disjoint_gr$revmap, function(i) i[1])
  # disjoint_gr_revmap_length = sapply(disjoint_gr$revmap, length)
  # disjoint_gr[disjoint_gr_revmap_length > 1]
  # disjoint_gr[disjoint_gr_revmap_length > 1]
  
  ## keep only those in which there is an overlap with a gene, i.e. revmap has to include the idxs 1:length(gr_genes)
  disjoint_gr = disjoint_gr[disjoint_gr_revmap_first %in% 1:length(GR_bins),]
  
  ## it might be that there are multiple entries in revmap but they all belong to genes, without any overlap.
  ## remove these cases
  only_genes_revmap <- sapply(disjoint_gr$revmap, function(i) all(i %in% 1:length(GR_bins)))
  disjoint_gr <- disjoint_gr[!only_genes_revmap]
  
  ## otherwise it means that there isn't a single overlap
  stopifnot(max(unlist(disjoint_gr$revmap)) > length(GR_bins))
  
  ## if there are more than two overlaps, split in in different rows
  which_overlaps <- sapply(disjoint_gr$revmap, length) > 2
  overlaps <- disjoint_gr[which_overlaps,]
  disjoint_gr <- disjoint_gr[!which_overlaps]
  cat('Number of entries with overlaps: ', length(disjoint_gr), '\n')
  
  ## in some cases there might be more than 2 revmaps, if there are overlaps in either of the two input granges
  
  ## split
  ## this takes too long
  if(length(overlaps) > 300){
    warning('This is going to take a long time\n')
  }
  if(length(overlaps)>0){
    additional_split_overlaps <- do.call('c', sapply(1:length(overlaps), function(idx){
      i = overlaps[idx,]
      ## get which indices are from genes, and which from annotated regions
      which_genes = which(i$revmap[[1]] <= length(GR_bins))
      which_not_genes = which(i$revmap[[1]] > length(GR_bins))
      do.call('c', sapply(which_genes, function(j){
        .x <- i
        sapply(which_not_genes, function(idx_not_genes){
          .x$revmap[[1]] <- c(.x$revmap[[1]][j], .x$revmap[[1]][which_not_genes])
          .x
        })
      }))
    }))
    cat('Concatenaning non-overlapped segments and split overlap segments\n')
    disjoint_gr <- c(disjoint_gr, additional_split_overlaps)
  }
  
  disjoint_gr_revmap_first = sapply(disjoint_gr$revmap, function(i) i[1])
  
  disjoint_gr_revmap_second = sapply(disjoint_gr$revmap, function(i) i[2])
  lengths_revmap <- sapply(disjoint_gr$revmap, length)
  
  table(is.na(disjoint_gr_revmap_second))
  
  ## remove sections which only contain non-genes
  # disjoint_gr[is.na(disjoint_gr_revmap_second),] ##???
  
  ## get the idx of the copy number segments (or the segments of interest)
  idx_gr = disjoint_gr_revmap_second-length(GR_bins)
  is.na(idx_gr)
  
  idx_gr[idx_gr <= 0 ] = NA
  
  disjoint_gr$CN_val = gr_CN_org@elementMetadata[,name_val_col][idx_gr] ### difference
  if(value_is_CN){
    replacement_val = 2
  }else{
    replacement_val= ''
  }
  disjoint_gr$CN_val[is.na(disjoint_gr$CN_val)] = replacement_val
  
  previous_revmap <- disjoint_gr$revmap
  disjoint_gr$revmap = disjoint_gr_revmap_first
  seqlevels(disjoint_gr) = as.character(unique(seqnames(disjoint_gr)))
  seqlengths(disjoint_gr) = as.numeric(as.vector(chromlens$Length[match(seqlevels(disjoint_gr), gsub("chr", "", chromlens$Chrom))]))
  seqlevels(GR_bins) = seqlevels(disjoint_gr)
  seqlengths(GR_bins) = seqlengths(disjoint_gr)
  
  # CN_RleList = mcolAsRleList(disjoint_gr, "CN_val") ## that takes forever even when length(disjoint_gr) is a VERY low value
  # averageCN = binnedAverage(bins = GR_bins, numvar = CN_RleList, varname = "CN_bin_averaged")
  # averageCN = binnedAverage(bins = GR_bins, numvar = disjoint_gr$revmap, varname = "CN_bin_averaged")
  
  if(value_is_CN){
    if(only_compute_overlaps){
      ## more optimal
      ## relevant revmap for the first part of indices (i.e. the ones relating to bins, e.g. genes)
      ## we want to map each of these relevant bins
      stopifnot(length(disjoint_gr_revmap_first) == length(disjoint_gr_revmap_second))
      relevant_revmap <- unique(disjoint_gr$revmap[!is.na(disjoint_gr_revmap_second)])
      length(relevant_revmap)
      ## now aggregating all segments that belong to the same gene
      averageCN2 <- data.frame(gene=rep(NA, length(unique(disjoint_gr$revmap))),
                               averageCN=rep(NA, length(unique(disjoint_gr$revmap))))
      dim(averageCN2)
      length(length(unique(disjoint_gr$revmap)))
      length(disjoint_gr_revmap_second)
      length(disjoint_gr_revmap_first)
      length(disjoint_gr$revmap)
      length(unique(disjoint_gr$revmap))
      # stopifnot(dim(averageCN2)[1] == length(disjoint_gr_revmap_second))
      GR_bins$names[unique(disjoint_gr$revmap) %in% relevant_revmap]
      averageCN2[unique(disjoint_gr$revmap) %in% relevant_revmap, 'gene'] = GR_bins$names[unique(disjoint_gr$revmap) %in% relevant_revmap]
      averageCN2[unique(disjoint_gr$revmap) %in% relevant_revmap, 'averageCN'] = sapply(unique(disjoint_gr$revmap)[unique(disjoint_gr$revmap) %in% relevant_revmap], function(i){
        cat(i, '\n')
        sum(width(disjoint_gr[disjoint_gr$revmap == i,]) * disjoint_gr[disjoint_gr$revmap == i,]$CN_val)/
          (sum(width(disjoint_gr[disjoint_gr$revmap == i,])))})
    }else{
      averageCN2 <- cbind.data.frame(gene=GR_bins$names[unique(disjoint_gr$revmap)],
                                     averageCN=sapply(unique(disjoint_gr$revmap), function(i){
                                       sum(width(disjoint_gr[disjoint_gr$revmap == i,]) * disjoint_gr[disjoint_gr$revmap == i,]$CN_val)/
                                         (sum(width(disjoint_gr[disjoint_gr$revmap == i,])))}))
    }
  }else{
    ## we do not want to compute any average
    disjoint_gr$second_revmap = disjoint_gr_revmap_second
    averageCN2 = disjoint_gr
    
  }
  
  return(averageCN2)
}


