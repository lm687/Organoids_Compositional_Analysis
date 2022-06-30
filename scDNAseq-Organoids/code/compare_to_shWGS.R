rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(GenomicRanges)
library(Biobase)
library(QDNAseq)
library(ggplot2)
library(cowplot)
library(parallel)
source("../../../cnsignatures/main_functions.R")
source("../../../cnsignatures/helper_functions.R")
source("../../copy_number_analysis_organoids/helper_functions.R")
source("helper.R")
library(dplyr)
library(reshape2)
# library(ggplot2)
# library(ggdendro)
# library(pheatmap)
# library(readxl)
# library(gridExtra)
# library(ggplotify)
# library(cowplot)

renaming <- readxl::read_excel("../../RNASeq_DE_resistant_sensitive/files/PDOnameProperSample_sWGS_RNAseq.xlsx")

organoid_list = c('PDO2', 'PDO3', 'PDO6')
organoid_list_large = renaming$PDO; organoid_list_large = organoid_list_large[grepl('PDO', organoid_list_large)]

absCN = list()
absCN_granges = list() ## mean over all organoids
absCN_granges_perorg = list() ## for each organoid independently
for(org_it in organoid_list){
  absCN[[org_it]] = readRDS(paste0("../robjects/fig2_subclonal_hclust", org_it, "_2.RDS"))
  coords = sapply(colnames(absCN[[org_it]]), function(i){
    .x = strsplit(i, "[.]")[[1]][-1]
    if(length(.x) == 6){
      as.numeric(c(paste0(.x[1:3], collapse = ""), paste0(.x[4:6], collapse = "")))
    }else if(length(.x) == 4){
      as.numeric(c(paste0(.x[1], collapse = ""), paste0(.x[2:4], collapse = "")))
    }
    })
  absCN_granges[[org_it]] = as(cbind.data.frame(chrom=clean_chrom(gsub("\\..*", "", colnames(absCN[[org_it]]))),
                                             start=coords[1,], end=coords[2,],
                                             CN_value=colMeans(absCN[[org_it]])), "GRanges")
  absCN_granges_perorg[[org_it]] = as(cbind.data.frame(chrom=rep(as.character(seqnames(absCN_granges[[org_it]])), each=nrow(absCN[[org_it]])),
                                                       start=rep(as.character(start(absCN_granges[[org_it]])), each=nrow(absCN[[org_it]])),
                                                       end=rep(as.character(end(absCN_granges[[org_it]])), each=nrow(absCN[[org_it]])),
                                                       CN_value=as.vector(absCN[[org_it]]),
                                                       cell=rep(1:nrow(absCN[[org_it]]), ncol(absCN[[org_it]]))), "GRanges")
  
}

organoid_absolute_CN = readRDS("../../copy_number_analysis_organoids/data/organoid_absolute_CN.rds")
segtables_organoids_absolute_copynumber = lapply(renaming$ID[match(organoid_list_large, renaming$PDO)],
                                                 function(samplename) getSegTable(organoid_absolute_CN[,samplename]))
names(segtables_organoids_absolute_copynumber) = organoid_list_large
names(absCN_granges_perorg) = organoid_list

segtables_organoids_absolute_copynumber = sapply(segtables_organoids_absolute_copynumber, as, "GRanges")

segtables_organoids_absolute_copynumber$PDO2
absCN_granges$PDO2

get_df_comparison = function(pdo){
  disjoin =   GenomicRanges::disjoin(c(segtables_organoids_absolute_copynumber[[pdo]],
                                       absCN_granges[[pdo]]), with.revmap=TRUE, ignore.strand=TRUE)
  
  ## keep only segments with both types of data
  disjoin = disjoin[sapply(disjoin$revmap, length) == 2,]
  revmap_idx1 = sapply(disjoin$revmap, `[`, 1)
  revmap_idx2 = sapply(disjoin$revmap, `[`, 2) - max(revmap_idx1)
  
  CN_sWGS_scDNA = cbind.data.frame(CN_shWGS=as.numeric(segtables_organoids_absolute_copynumber[[pdo]]$segVal[revmap_idx1]),
                   CN_mean_scDNA=absCN_granges[[pdo]]$CN_value[revmap_idx2])
  
  return(CN_sWGS_scDNA)
}

df_comparison = lapply(organoid_list, get_df_comparison); names(df_comparison) = organoid_list
a <- ggplot(df_comparison$PDO2, aes(x=CN_shWGS, y=CN_mean_scDNA))+
  geom_abline(slope = 1, intercept = 0, lty='dashed', alpha=0.8)+geom_point()+ggtitle("PDO2")
b <- ggplot(df_comparison$PDO3, aes(x=CN_shWGS, y=CN_mean_scDNA))+
  geom_abline(slope = 1, intercept = 0, lty='dashed', alpha=0.8)+geom_point()+ggtitle("PDO3")
c <- ggplot(df_comparison$PDO6, aes(x=CN_shWGS, y=CN_mean_scDNA))+
  geom_abline(slope = 1, intercept = 0, lty='dashed', alpha=0.8)+geom_point()+ggtitle("PDO6")

pdf("../plots/comparison_sWGS_scDNA_log2.pdf", width = 11, height = 4)
cowplot::plot_grid(a+scale_x_continuous(trans = "log2")+scale_y_continuous(trans = "log2"),
                   b+scale_x_continuous(trans = "log2")+scale_y_continuous(trans = "log2"),
                   c+scale_x_continuous(trans = "log2")+scale_y_continuous(trans = "log2"), nrow=1)
dev.off()

pdf("../plots/comparison_sWGS_scDNA.pdf", width = 11, height = 4)
cowplot::plot_grid(a,
                   b,
                   c, nrow=1)
dev.off()

#----------------------------------------------------------------------------------------------------#

## now with bed files instead of csv
renaming[match(organoid_list, renaming$PDO),'ID'] %>% unlist
absCN_bed = lapply(c("23868org", "118976org", "119148orgb"),
                   function(i) read.table(paste0("../UID-", i, ".bed")))
names(absCN_bed) = organoid_list

# pdo = "PDO2"
get_df_comparison_2 = function(pdo){
  
  .x_grouped_df = as(cbind.data.frame(chrom=absCN_bed[[pdo]]$V1, start=absCN_bed[[pdo]]$V2, end=absCN_bed[[pdo]]$V3,
                                      segVal=absCN_bed[[pdo]]$V5, names=absCN_bed[[pdo]]$V4), 'GRanges')
  .dis = GenomicRanges::disjoin(.x_grouped_df, with.revmap=T)
  .dis$segVal = sapply(.dis$revmap, function(k) mean(.x_grouped_df$segVal[k[[1]]]))
  
  disjoin =   GenomicRanges::disjoin(c(segtables_organoids_absolute_copynumber[[pdo]],
                                       .dis), with.revmap=TRUE, ignore.strand=TRUE)
  disjoin = disjoin[sapply(disjoin$revmap, length) == 2,]
  revmap_idx1 = sapply(disjoin$revmap, `[`, 1)
  revmap_idx2 = sapply(disjoin$revmap, `[`, 2) - max(revmap_idx1)
  
  CN_sWGS_scDNA = cbind.data.frame(CN_shWGS=as.numeric(segtables_organoids_absolute_copynumber[[pdo]]$segVal[revmap_idx1]),
                                   CN_mean_scDNA=.dis$segVal[revmap_idx2])

  return(list(CN_sWGS_scDNA, ggplot(CN_sWGS_scDNA, aes(x=CN_shWGS, y=CN_mean_scDNA))+
    geom_abline(slope = 1, intercept = 0, lty='dashed', alpha=0.8)+geom_point()))
}

df_comparison_from_bed = lapply(organoid_list, get_df_comparison_2); names(df_comparison_from_bed) = organoid_list
a_frombed <- ggplot(df_comparison_from_bed$PDO2[[1]], aes(x=CN_shWGS, y=CN_mean_scDNA))+
  geom_abline(slope = 1, intercept = 0, lty='dashed', alpha=0.8)+geom_point()+ggtitle("PDO2")
b_frombed <- ggplot(df_comparison_from_bed$PDO3[[1]], aes(x=CN_shWGS, y=CN_mean_scDNA))+
  geom_abline(slope = 1, intercept = 0, lty='dashed', alpha=0.8)+geom_point()+ggtitle("PDO3")
c_frombed <- ggplot(df_comparison_from_bed$PDO6[[1]], aes(x=CN_shWGS, y=CN_mean_scDNA))+
  geom_abline(slope = 1, intercept = 0, lty='dashed', alpha=0.8)+geom_point()+ggtitle("PDO6")

pdf("../plots/comparison_sWGS_scDNA_frombed.pdf", width = 11, height = 4)
cowplot::plot_grid(a_frombed,#+scale_x_continuous(trans = "log2")+scale_y_continuous(trans = "log2"),
                   b_frombed,#+scale_x_continuous(trans = "log2")+scale_y_continuous(trans = "log2"),
                   c_frombed,#+scale_x_continuous(trans = "log2")+scale_y_continuous(trans = "log2"),
                   nrow=1)
dev.off()

#---------------------------------------------------------------------------------------#
## in regions of interest, plot the single cell and the sgWGS together, for scDNA in 20kb


absCN_bed.granges = lapply(organoid_list, function(org){
  .x = as(cbind.data.frame(chromosome=absCN_bed[[org]]$V1, start=absCN_bed[[org]]$V2,
                    end=absCN_bed[[org]]$V3),
   "GRanges")
  .x$cell = absCN_bed[[org]]$V4
  .x$CN_value = absCN_bed[[org]]$V5
  .x$event_confidence = absCN_bed[[org]]$V6
  return(.x)
})
names(absCN_bed.granges) = organoid_list
absCN_bed.granges

chrlen = readRDS("../../copy_number_analysis_organoids/data/chrlen.RDS")

# list_params = list(segtab <- segtables_organoids_absolute_copynumber$PDO2,
#                    chrom_AOI=7
#                    start_AOI=1
#                    end_AOI=chrlen %>% filter(V1 == "chr7") %>% dplyr::select(V2) %>% as.integer()
#                    subclonal_line1=2
#                    subclonal_line2=3
#                    title='Chrom 7 in PDO2')

list_params = list(list(segtab=segtables_organoids_absolute_copynumber$PDO2,
                        PDO='PDO2',
                            chrom_AOI=7, start_AOI=1,
                            end_AOI=chrlen %>% filter(V1 == "chr7") %>% dplyr::select(V2) %>% as.integer(),
                            subclonal_line1=2, subclonal_line2=3, title='Chrom 7 in PDO2'),
list(segtab=segtables_organoids_absolute_copynumber$PDO2,
     PDO='PDO2',
     chrom_AOI=10, start_AOI=1,
                            end_AOI=135534747/4,
                            subclonal_line1=2, subclonal_line2=3, title='Chrom 10 (start) in PDO2'),
list(segtab=segtables_organoids_absolute_copynumber$PDO2,
     PDO='PDO2',
     chrom_AOI=10,
                            start_AOI=chrlen %>% filter(V1 == "chr10") %>% dplyr::select(V2) %>% as.integer() / 4,
                            end_AOI=chrlen %>% filter(V1 == "chr10") %>% dplyr::select(V2) %>% as.integer(),
                            subclonal_line1=2, subclonal_line2=3, title='Chrom 10 (start) in PDO2'),
list(segtab=segtables_organoids_absolute_copynumber$PDO3,
     PDO='PDO3',
     chrom_AOI=5,
                            start_AOI=51200001,
                            end_AOI=66560000,
                            subclonal_line1=2, subclonal_line2=1, title='Chrom 5 (start) in PDO3'),
list(segtab=segtables_organoids_absolute_copynumber$PDO3,
     PDO='PDO3',
     chrom_AOI=5,
                            start_AOI=chrlen %>% filter(V1 == "chr5") %>% dplyr::select(V2) %>% as.integer() * 2/3,
                            end_AOI=chrlen %>% filter(V1 == "chr5") %>% dplyr::select(V2) %>% as.integer(),
                            subclonal_line1=2, subclonal_line2=3, title='Chrom 5 (end) in PDO3'),
list(segtab=segtables_organoids_absolute_copynumber$PDO6,
     PDO='PDO6',
     chrom_AOI=4,
                            start_AOI=71680001,
                            end_AOI=81920000,
                            subclonal_line1=2, subclonal_line2=3, title='Chrom 4 (end) in PDO6'),
list(segtab=segtables_organoids_absolute_copynumber$PDO6,
     PDO='PDO6',
     chrom_AOI=7,
     start_AOI=15360001,
     end_AOI=20480000,
     subclonal_line1=1, subclonal_line2=2, title='Chrom 7 (start) in PDO6'),
list(segtab=segtables_organoids_absolute_copynumber$PDO3,
     PDO='PDO3',
     chrom_AOI=5,
     start_AOI=51200001,
     end_AOI=56320000,
     subclonal_line1=1, subclonal_line2=2, title='Chrom 5 (start) in PDO3'),
list(segtab=segtables_organoids_absolute_copynumber$PDO6,
     PDO='PDO6',
     chrom_AOI=10,
     start_AOI=10240001,
     end_AOI=15360000,
     subclonal_line1=3, subclonal_line2=15, title='Chrom 10 (start) in PDO6'),
list(segtab=segtables_organoids_absolute_copynumber$PDO6,
     PDO='PDO6',
     chrom_AOI=2,
     start_AOI=138240001,
     end_AOI=143360000,
     subclonal_line1=2, subclonal_line2=3, title='Chrom 2 in PDO6')
)

list_params[length(list_params)][[1]]$start_AOI/chrlen[2,'V2']
list_params[length(list_params)][[1]]$end/chrlen[2,'V2']


# give_plot_region = function(region_list_params, abs_sc_GR, return_plot){
#   segtabGR <- as(data.frame(region_list_params$segtab),"GRanges")
#   region_subclonal = as(cbind.data.frame(chromosome=region_list_params$chrom_AOI, start=region_list_params$start_AOI,
#                                          end=region_list_params$end_AOI),
#                         "GRanges")
#   
#   region_subclonal
#   
#   disjoin_region_subclonal_sWGS <- GenomicRanges::disjoin(c(region_subclonal, segtabGR),
#                                                           with.revmap=TRUE, ignore.strand=TRUE)
#   
#   names_cells <- sort(unique(abs_sc_GR[[region_list_params$PDO]]$cell))
#   disjoin_region_subclonal_scDNA <- lapply(names_cells, function(idx_cell){
#     GenomicRanges::disjoin(c(
#     region_subclonal,
#     abs_sc_GR[[region_list_params$PDO]][abs_sc_GR[[region_list_params$PDO]]$cell == idx_cell,]),
#     with.revmap=TRUE, ignore.strand=TRUE)
#   })
#   names(disjoin_region_subclonal_scDNA) = names_cells
#   
#   disjoin_region_subclonal_sWGS = disjoin_region_subclonal_sWGS[sapply(disjoin_region_subclonal_sWGS$revmap, function(i) 1 %in% i),]
#   disjoin_region_subclonal_scDNA = lapply(disjoin_region_subclonal_scDNA, function(i) i[sapply(i$revmap, function(i) 1 %in% i),])
#   
#   # disjoin_region_subclonal_sWGS
#   # disjoin_region_subclonal_scDNA
#   # 
#   
#   CN_sWGS = sapply(disjoin_region_subclonal_sWGS$revmap, function(i){
#     if(!is.na(i[2])){
#       as(data.frame(region_list_params$segtab),"GRanges")$segVal[i[2]]
#     }else{
#       2
#     }
#   })
#   disjoin_region_subclonal_sWGS$CNval = as.numeric(CN_sWGS)
#   
#   disjoin_region_subclonal_scDNA= lapply(disjoin_region_subclonal_scDNA, function(j){
#     j[width(j) > 1]
#   })
#   len_firstGR = length(region_subclonal)
#   CN_scDNA = mclapply(names(disjoin_region_subclonal_scDNA), function(cell_name){
#     j = disjoin_region_subclonal_scDNA[[cell_name]]
#     sapply(j$revmap, function(i){
#     # if(!is.na(i[2])){
#     if(length(i) > 1){
#       as.numeric(abs_sc_GR[[region_list_params$PDO]][abs_sc_GR[[region_list_params$PDO]]$cell == cell_name,][i[2]-len_firstGR]$CN_value)
#     }else{
#       2
#     }
#     })
#   }) ### this can be sped up!!!
#   event_conf = mclapply(names(disjoin_region_subclonal_scDNA), function(cell_name){
#     j = disjoin_region_subclonal_scDNA[[cell_name]]
#     sapply(j$revmap, function(i){
#       # if(!is.na(i[2])){
#       if(length(i) > 1){
#         as.numeric(abs_sc_GR[[region_list_params$PDO]][abs_sc_GR[[region_list_params$PDO]]$cell == cell_name,][i[2]-len_firstGR]$event_confidence)
#       }else{
#         2
#       }
#     })
#   })
#   
#   for(i in 1:length(disjoin_region_subclonal_scDNA)){
#     disjoin_region_subclonal_scDNA[[i]]$CNval = as.numeric(CN_scDNA[[i]])
#     disjoin_region_subclonal_scDNA[[i]]$event_confidence = event_conf[[i]]
#   }
#   
#   plt1 <- ggplot(cbind.data.frame(start=start(disjoin_region_subclonal_sWGS),
#                           end=end(disjoin_region_subclonal_sWGS),
#                           CNval=disjoin_region_subclonal_sWGS$CNval))+
#     geom_rect(aes(xmin=start, xmax=end, ymin=0, ymax=CNval), alpha=0.2)+
#     geom_abline(slope = 0, intercept = region_list_params$subclonal_line1, color='blue', lty='dashed')+
#     geom_abline(slope = 0, intercept = region_list_params$subclonal_line2, color='blue', lty='dashed')+
#     ggtitle(region_list_params$title)
#   
#   
#   disjoin_region_subclonal_scDNA_df <- lapply(disjoin_region_subclonal_scDNA, function(i){
#     cbind.data.frame(start=start(i),
#                    end=end(i),
#                    CNval=i$CNval,
#                    event_confidence=i$event_confidence)
#   })
#   for(i in 1:length(disjoin_region_subclonal_scDNA_df)){
#     disjoin_region_subclonal_scDNA_df[[i]] = cbind(disjoin_region_subclonal_scDNA_df[[i]], names_cells[i])
#   }
#   disjoin_region_subclonal_scDNA_df = do.call('rbind', disjoin_region_subclonal_scDNA_df)
#   
#   plt2 <- ggplot(disjoin_region_subclonal_scDNA_df)+
#     geom_rect(aes(xmin=start, xmax=end, ymin=0, ymax=CNval,group=`names_cells[i]`, col=`names_cells[i]`),fill=NA, alpha=0.1, size=0.1)+
#     geom_abline(slope = 0, intercept = region_list_params$subclonal_line1, color='blue', lty='dashed')+
#     geom_abline(slope = 0, intercept = region_list_params$subclonal_line2, color='blue', lty='dashed')+
#     geom_rect(data=cbind.data.frame(start=start(disjoin_region_subclonal_sWGS),
#                                     end=end(disjoin_region_subclonal_sWGS),
#                                     CNval=disjoin_region_subclonal_sWGS$CNval),
#               aes(xmin=start, xmax=end, ymin=0, ymax=CNval), col='red', fill=NA, alpha=0.01)+
#     ggtitle(region_list_params$title)+
#     scale_y_continuous(trans = "log2")
#   if(return_plot){
#     return(list(plt1, plt2))
#   }else{
#     return(list(scDNA=disjoin_region_subclonal_scDNA_df, sWGS=disjoin_region_subclonal_sWGS))
#   }
# }

give_plot_region_v2 = function(region_list_params, abs_sc_GR, return_plot, include_event_confidence=T, include_sc=T){
  
  #------------------------------------------------------------------------------#
  
  segtabGR <- as(data.frame(region_list_params$segtab),"GRanges")
  region_subclonal = as(cbind.data.frame(chromosome=region_list_params$chrom_AOI, start=region_list_params$start_AOI,
                                         end=region_list_params$end_AOI),
                        "GRanges")
  len_firstGR <- length(region_subclonal)
  
  #------------------------------------------------------------------------------#
  
  disjoin_region_subclonal_sWGS <- GenomicRanges::disjoin(c(region_subclonal, segtabGR),
                                                          with.revmap=TRUE, ignore.strand=TRUE)
  disjoin_region_subclonal_sWGS = disjoin_region_subclonal_sWGS[sapply(disjoin_region_subclonal_sWGS$revmap, function(i) 1 %in% i),]
  CN_sWGS = sapply(disjoin_region_subclonal_sWGS$revmap, function(i){
    if(!is.na(i[2])){
      as.numeric(as.character(segtabGR[i[2]-len_firstGR]$segVal))
    }else{
      2
    }
  })
  disjoin_region_subclonal_sWGS$CNval = CN_sWGS
  
  disjoin_region_subclonal_scDNA <- 
    GenomicRanges::disjoin(c(
      region_subclonal,
      abs_sc_GR[[region_list_params$PDO]]),
      with.revmap=TRUE, ignore.strand=TRUE)

  disjoin_region_subclonal_scDNA = disjoin_region_subclonal_scDNA[sapply(disjoin_region_subclonal_scDNA$revmap, function(i) 1 %in% i),]
  
  #------------------------------------------------------------------------------#
  
  if(include_sc){
    CN_scDNA = lapply(disjoin_region_subclonal_scDNA$revmap, function(i){
        # if(!is.na(i[2])){
        if(length(i) > 1){
          as.numeric(abs_sc_GR[[region_list_params$PDO]][i[-1]-len_firstGR]$CN_value)
        }else{
          2
        }
      })
    if(include_event_confidence){
      event_conf = lapply(disjoin_region_subclonal_scDNA$revmap, function(i){
        if(length(i) > 1){
          as.numeric(abs_sc_GR[[region_list_params$PDO]][i[-1]-len_firstGR]$event_confidence)
        }else{
          NA
        }
      })
    }
    cell = lapply(disjoin_region_subclonal_scDNA$revmap, function(i){
      if(length(i) > 1){
        as.numeric(abs_sc_GR[[region_list_params$PDO]][i[-1]-len_firstGR]$cell)
      }else{
        NA
      }
    })
  }
  
  #------------------------------------------------------------------------------#
  
  plt1 <- ggplot(cbind.data.frame(start=start(disjoin_region_subclonal_sWGS),
                                  end=end(disjoin_region_subclonal_sWGS),
                                  CNval=disjoin_region_subclonal_sWGS$CNval))+
    geom_rect(aes(xmin=start, xmax=end, ymin=0, ymax=CNval), alpha=0.2)+
    geom_abline(slope = 0, intercept = region_list_params$subclonal_line1, color='blue', lty='dashed')+
    geom_abline(slope = 0, intercept = region_list_params$subclonal_line2, color='blue', lty='dashed')+
    ggtitle(region_list_params$title)
  
  #------------------------------------------------------------------------------#
  
  if(include_sc){
    disjoin_region_subclonal_scDNA_df = do.call('rbind', lapply(1:length(disjoin_region_subclonal_scDNA), function(idx_row){
      as.data.frame(rep(disjoin_region_subclonal_scDNA[idx_row], length(cell[[idx_row]])))
    }))
    disjoin_region_subclonal_scDNA_df$cell = (unlist(cell))
    disjoin_region_subclonal_scDNA_df$CNval = (unlist(CN_scDNA))
    if(include_event_confidence){
      disjoin_region_subclonal_scDNA_df$event_confidence = (unlist(event_conf))
    }else{
      disjoin_region_subclonal_scDNA_df$event_confidence = NA
    }
    
    #------------------------------------------------------------------------------#
    
    plt2 <- ggplot(disjoin_region_subclonal_scDNA_df)+
      geom_rect(aes(xmin=start, xmax=end, ymin=0, ymax=CNval,group=cell, col=cell),fill=NA, alpha=0.1, size=0.1)+
      geom_abline(slope = 0, intercept = region_list_params$subclonal_line1, color='blue', lty='dashed')+
      geom_abline(slope = 0, intercept = region_list_params$subclonal_line2, color='blue', lty='dashed')+
      geom_rect(data=cbind.data.frame(start=start(disjoin_region_subclonal_sWGS),
                                      end=end(disjoin_region_subclonal_sWGS),
                                      CNval=disjoin_region_subclonal_sWGS$CNval),
                aes(xmin=start, xmax=end, ymin=0, ymax=CNval), col='red', fill=NA, alpha=0.01)+
      ggtitle(region_list_params$title)+
      scale_y_continuous(trans = "log2")
  }else{
    plt2 <- ggplot()
    disjoin_region_subclonal_scDNA_df = NA
  }
  
  #------------------------------------------------------------------------------#
  
  if(return_plot){
    return(list(plt1, plt2))
  }else{
    return(list(scDNA=disjoin_region_subclonal_scDNA_df, sWGS=disjoin_region_subclonal_sWGS))
  }
}

# for(k in 1:length(absCN_granges)){
#   absCN_granges[[k]]$CN_value = absCN_granges[[k]]$segVal
# }

# region_list_params=list_params[[2]]; abs_sc_GR=absCN_bed.granges; return_plot=T
# all_plots_scDNA_sWGS_regions = mclapply(list_params, give_plot_region, abs_sc_GR = absCN_bed.granges)
all_plots_scDNA_sWGS_regions_v2 = lapply(list_params[2], give_plot_region_v2, abs_sc_GR = absCN_bed.granges, return_plot=T, include_event_confidence=T)
# all_plots_scDNA_sWGS_regions_largebins = mclapply(list_params, give_plot_region, abs_sc_GR=absCN_granges_perorg)
all_plots_scDNA_sWGS_regions_largebins_v2 = mclapply(list_params[1], give_plot_region_v2, abs_sc_GR=absCN_granges_perorg, return_plot=T, include_event_confidence=F)

# pdf("../plots/comparison_sWGS_scDNA_specificregions.pdf", width = 10, height = 4)
# for(k in 1:length(all_plots_scDNA_sWGS_regions)){
#   print(plot_grid(all_plots_scDNA_sWGS_regions[[k]][[1]], all_plots_scDNA_sWGS_regions[[k]][[2]], nrow=1, rel_widths=c(0.3, 0.7)))
# }
# dev.off()

pdf("../plots/comparison_sWGS_scDNA_specificregions_v2.pdf", width = 10, height = 4)
for(k in 1:length(all_plots_scDNA_sWGS_regions_v2)){
  print(plot_grid(all_plots_scDNA_sWGS_regions_v2[[k]][[1]], all_plots_scDNA_sWGS_regions_v2[[k]][[2]], nrow=1, rel_widths=c(0.3, 0.7)))
}
dev.off()


# pdf("../plots/comparison_sWGS_scDNA_specificregions_largebins.pdf", width = 10, height = 4)
# for(k in 1:length(all_plots_scDNA_sWGS_regions_largebins)){
#   print(plot_grid(all_plots_scDNA_sWGS_regions_largebins[[k]][[1]], all_plots_scDNA_sWGS_regions_largebins[[k]][[2]], nrow=1, rel_widths=c(0.3, 0.7)))
# }
# dev.off()
# 
# pdf("../plots/comparison_sWGS_scDNA_specificregions_largebins_normalscale.pdf", width = 10, height = 4)
# for(k in 1:length(all_plots_scDNA_sWGS_regions_largebins)){
#   print(plot_grid(all_plots_scDNA_sWGS_regions_largebins[[k]][[1]]+scale_y_continuous(trans = "identity"),
#                   all_plots_scDNA_sWGS_regions_largebins[[k]][[2]]+scale_y_continuous(trans = "identity"),
#                   nrow=1, rel_widths=c(0.3, 0.7)))
# }
# dev.off()
pdf("../plots/comparison_sWGS_scDNA_specificregions_largebins_normalscale_v2.pdf", width = 10, height = 4)
for(k in 1:length(all_plots_scDNA_sWGS_regions_largebins_v2)){
  print(plot_grid(all_plots_scDNA_sWGS_regions_largebins_v2[[k]][[1]]+scale_y_continuous(trans = "identity"),
                  all_plots_scDNA_sWGS_regions_largebins_v2[[k]][[2]]+scale_y_continuous(trans = "identity"),
                  nrow=1, rel_widths=c(0.3, 0.7)))
}
dev.off()


## looking at specific ones in detail
all_plots_scDNA_sWGS_regions_ex1 = lapply(list_params[2], give_plot_region_v2, abs_sc_GR=absCN_bed.granges, return_plot = F)
all_plots_scDNA_sWGS_regions_ex1b = lapply(list_params[1], give_plot_region_v2, abs_sc_GR=absCN_bed.granges, return_plot = F)
all_plots_scDNA_sWGS_regions_ex2 = lapply(list_params[5], give_plot_region_v2, abs_sc_GR=absCN_bed.granges, return_plot = F)
all_plots_scDNA_sWGS_regions_ex3 = lapply(list_params[6], give_plot_region_v2, abs_sc_GR=absCN_bed.granges, return_plot = F)
all_plots_scDNA_sWGS_regions_ex3b = lapply(list_params[7], give_plot_region_v2, abs_sc_GR=absCN_bed.granges, return_plot = F)
all_plots_scDNA_sWGS_regions_ex2b = lapply(list_params[8], give_plot_region_v2, abs_sc_GR=absCN_bed.granges, return_plot = F)
all_plots_scDNA_sWGS_regions_ex3c = lapply(list_params[9], give_plot_region_v2, abs_sc_GR=absCN_bed.granges, return_plot = F)
all_plots_scDNA_sWGS_regions_ex3d = lapply(list_params[10], give_plot_region_v2, abs_sc_GR=absCN_bed.granges, return_plot = F)
all_plots_scDNA_sWGS_regions_ex3d_largebins = lapply(list_params[10], give_plot_region_v2, abs_sc_GR=absCN_granges_perorg, return_plot = F, include_event_confidence = F)


plot(density(all_plots_scDNA_sWGS_regions_ex1[[1]][[1]]$event_confidence))

todf = function(disjoin_region_subclonal_sWGS){
  cbind.data.frame(start=start(disjoin_region_subclonal_sWGS),
                 end=end(disjoin_region_subclonal_sWGS),
                 CNval=disjoin_region_subclonal_sWGS$CNval)
}

disjoin_region_subclonal_sWGS_ex1 = todf(all_plots_scDNA_sWGS_regions_ex1[[1]][[2]])

saveRDS(list(region1=list(all_plots_scDNA_sWGS_regions_ex1[[1]],list_params[2]),
             region1b=list(all_plots_scDNA_sWGS_regions_ex1b[[1]], list_params[1]),
             region2=list(all_plots_scDNA_sWGS_regions_ex2[[1]], list_params[5]),
             region3=list(all_plots_scDNA_sWGS_regions_ex3[[1]],list_params[6]),
             region3b=list(all_plots_scDNA_sWGS_regions_ex3b[[1]],list_params[7]),
             region2b=list(all_plots_scDNA_sWGS_regions_ex2b[[1]],list_params[8]),
             region3c=list(all_plots_scDNA_sWGS_regions_ex3c[[1]],list_params[9]),
             region3d=list(all_plots_scDNA_sWGS_regions_ex3d[[1]],list_params[10])),
file = "../robjects/fig2_regions_subclonalb.RDS")

ggplot(all_plots_scDNA_sWGS_regions_ex1[[1]][[1]] %>% filter(event_confidence > 10)# %>%
         # filter(`names_cells[i]` == 10)
                                                             )+
  geom_rect(aes(xmin=start, xmax=end, ymin=0, ymax=CNval,group=`names_cells[i]`), col='black',
            fill=NA, alpha=0.1, size=0.1)+
  # geom_step(aes(x=start, y=CNval,group=`names_cells[i]`, col=`names_cells[i]`),
  #           alpha=1)+
  # geom_abline(slope = 0, intercept = region_list_params$subclonal_line1, color='blue', lty='dashed')+
  # geom_abline(slope = 0, intercept = region_list_params$subclonal_line2, color='blue', lty='dashed')+
  geom_rect(data=disjoin_region_subclonal_sWGS_ex1,
            aes(xmin=start, xmax=end, ymin=0, ymax=CNval), col='red', fill=NA, alpha=0.01)+
  ggtitle("Chrom 10 (start) in PDO2")


ggplot()+
  geom_rect(data = all_plots_scDNA_sWGS_regions_ex3d_largebins[[1]][[1]],
            aes(xmin=start, xmax=end, ymin=0, ymax=CNval,group=cell), col='blue',
            fill=NA, alpha=0.1, size=0.1)+
  geom_rect(data=data.frame(all_plots_scDNA_sWGS_regions_ex3d_largebins[[1]][[2]]),
            aes(xmin=start, xmax=end, ymin=0, ymax=CNval), col='red', fill=NA, alpha=0.01)+
  ggtitle("")

ggplot()+
  geom_rect(data = all_plots_scDNA_sWGS_regions_ex3d[[1]][[1]],
            aes(xmin=start, xmax=end, ymin=0, ymax=CNval,group=cell), col='black',
            fill=NA, alpha=0.1, size=0.1)+
  geom_rect(data = all_plots_scDNA_sWGS_regions_ex3d_largebins[[1]][[1]],
            aes(xmin=start, xmax=end, ymin=0, ymax=CNval,group=cell), col='blue',
            fill=NA, alpha=0.1, size=0.1)+
  geom_rect(data=data.frame(all_plots_scDNA_sWGS_regions_ex3d_largebins[[1]][[2]]),
            aes(xmin=start, xmax=end, ymin=0, ymax=CNval), col='red', fill=NA, alpha=0.01)+
  ggtitle("")

ggplot()+
  geom_rect(data = all_plots_scDNA_sWGS_regions_ex3d[[1]][[1]] %>%
              filter(event_confidence > 10),
            aes(xmin=start, xmax=end, ymin=0, ymax=CNval,group=cell), col='black',
            fill=NA, alpha=0.1, size=0.1)
  
.xx <- all_plots_scDNA_sWGS_regions_ex3d[[1]][[1]] %>%
  filter(event_confidence > 10)
plot(0, xlim=c(min(.xx$start), max(.xx$end)), ylim=c(min(.xx$CNval), max(.xx$CNval)))
for(i in 1:nrow(.xx)){
  segments(x0 = .xx$start[i], x1 = .xx$end[i], y0 = .xx$CNval[i], y1 = .xx$CNval[i], col=alpha(colour = "black", alpha = 0.01))
}

.xx2 <- .xx %>% group_by(cell) %>% dplyr::summarise(sum(CNval*width)/sum(width))
plot(density(.xx2$`sum(CNval * width)/sum(width)`))

## a single cell in chrom 2 of PDO6
ggplot(absCN_bed$PDO6 %>% filter(V1 == 2, V4 %in% 100:102), aes(x=V2, y=V5, col=V4))+
  geom_step()+lims(y=c(-1,10))+
  geom_segment(aes(x = specific_regions[[2]][[2]][[1]]$start_AOI, xend= specific_regions[[2]][[2]][[1]]$end_AOI,
                   y=-0.5, yend=-0.5))

# .xx <- all_plots_scDNA_sWGS_regions_ex3d[[1]][[1]] %>%
#   filter(event_confidence > 10)
# plot(0, xlim=c(min(.xx$start), max(.xx$end)), ylim=c(min(.xx$CNval), max(.xx$CNval)))
# for(i in 1:nrow(.xx)){
#   segments(x0 = .xx$start[i], x1 = .xx$end[i], y0 = .xx$CNval[i], y1 = .xx$CNval[i], col=alpha(colour = "black", alpha = 0.01))
# }


#---------------------------------------------------------------------------------------#
## Look at the genes in the regions of interest

specific_regions <- readRDS("../../scDNAseq-Organoids/robjects/fig2_regions_subclonalb.RDS")
specific_regions[[2]] ## pdo2
specific_regions[[6]] ## pdo3
specific_regions[[8]] ## pdo6

gtf.file <- file.path("../../RNASeq_and_CN/20191218_ViasM_BJ_orgaBrs/Data/Homo_sapiens.GRCh37.87.gtf.gz")
sqlite_file <- '../../RNASeq_and_CN/20191218_ViasM_BJ_orgaBrs/Data/Homo_sapiens.GRCh37.87.sqlite'
sqlite_path <- file.path(sqlite_file)

if(!file.exists(sqlite_path)) {
  ## generate the SQLite database file
  ensembldb::ensDbFromGtf(gtf=gtf.file, path=ref_dir, outfile=sqlite_file)
}
EnsDb.Hsapiens.v87 <- ensembldb::EnsDb(sqlite_file)

# Genes, used to annotated the TPM matrix to send to Maria
ag <- ensembldb::genes(EnsDb.Hsapiens.v87, filter=list(AnnotationFilter::GeneBiotypeFilter('protein_coding')), return.type="DataFrame") 
ag

specific_regions[[2]]
specific_regions[[2]]

genes = ag[ag$seq_name == 4,]
df_chrom = data.frame(apply(as(segtables_organoids_absolute_copynumber$PDO6[seqnames(segtables_organoids_absolute_copynumber$PDO6) == 4], "data.frame"),
                 2, as.numeric))
df_chrom = rbind(c(df_chrom[1,1], -Inf, df_chrom[1,2], NA, 2), df_chrom, c(df_chrom[nrow(df_chrom),1], df_chrom[nrow(df_chrom), 3], -Inf, NA, NA))

ggplot(df_chrom)+geom_step(aes(x=start,  y=segVal))+geom_abline(slope = 0, intercept = 2, col='red', lty='dashed')+
  geom_segment(data = data.frame(genes), aes(x=gene_seq_start, xend=gene_seq_end, y=0, yend=0.5))+
  geom_segment(data = data.frame(genes) %>% filter(gene_name %in% c('SLC34A2', 'IL2', 'PPARGC1A')),
               aes(x=gene_seq_start, xend=gene_seq_end, y=0, yend=0.5), col='blue')

list_AOI <- list(specific_regions[[2]],
     specific_regions[[6]],
     specific_regions[[8]])
genes_in_AOI = lapply(list_AOI,
 function(i) data.frame(ag) %>%
   filter (seq_name == i[[2]][[1]]$chrom_AOI, gene_seq_start > i[[2]][[1]]$start_AOI,
           gene_seq_end < i[[2]][[1]]$end_AOI)  )

sapply(genes_in_AOI, nrow)

plots_genes <- lapply(1:3, function(i) ggplot(list_AOI[[i]][[1]][[1]])+geom_step(aes(x=start,  y=CNval))+
         geom_abline(slope = 0, intercept = 2, col='red', lty='dashed')+
         geom_segment(data = data.frame(genes_in_AOI[[i]]), aes(x=gene_seq_start, xend=gene_seq_end, y=0, yend=0.5)))#+
         # geom_segment(data = data.frame(genes_in_AOI[[i]]) %>% filter(gene_name %in% c('SLC34A2', 'IL2', 'PPARGC1A')),
         #       aes(x=gene_seq_start, xend=gene_seq_end, y=0, yend=0.5), col='blue'))

plots_genes2 <- lapply(1:3, function(i) ggplot()+
                        geom_abline(slope = 0, intercept = 2, col='red', lty='dashed')+
                        geom_segment(data = data.frame(genes_in_AOI[[i]]), aes(x=gene_seq_start, xend=gene_seq_end, y=0, yend=0.5)))
                        # geom_segment(data = data.frame(genes_in_AOI[[i]]) %>% filter(gene_name %in% c('SLC34A2', 'IL2', 'PPARGC1A')),
                        #              aes(x=gene_seq_start, xend=gene_seq_end, y=0, yend=0.5), col='blue'))

plots_genes[[1]]
plots_genes[[2]]
plots_genes[[3]]
plots_genes2[[1]]
plots_genes2[[2]]
plots_genes2[[3]]

pdo2specific <- ggplot(specific_regions[[2]][[1]][[1]] %>% filter(event_confidence > 5, width >1) %>%
                         group_by(cell) %>% 
                         mutate(mean_CN = mean(CNval)))+
  geom_rect(aes(xmin=start, xmax=end, ymin=0, ymax=CNval,group=cell),
            fill='black', alpha=0.008, size=0.001)+
  geom_step(data=as(specific_regions[[2]][[1]][[2]], "data.frame"),
            aes(x=start, y=CNval), col='red', fill=NA)+
  ggtitle(specific_regions[[2]][[2]][[1]]$title)+
  lims(y=c(0,6))
pdo2specific

add_extra_rows_step_percell <- function(df2){
  do.call('rbind', lapply(unique(df2$cell), function(i){
    df = df2[df2$cell == i,]
    df = df[order(df$start),]
    rbind.data.frame(c(df[1,1], -Inf, df[1,2], NA, NA, NA, df[1,7]),
                     df,
                     c(df[nrow(df),1], df[nrow(df),3], Inf, NA, NA, NA, df[nrow(df),7]))
  }))
}

pdo6specific <- ggplot(add_extra_rows_step_percell(specific_regions[[8]][[1]][[1]] %>%
                         filter(event_confidence > 10, width >1)))+
  geom_step(aes(x=start, y=CNval,group=cell))+
  geom_step(data=add_extra_rows_step(as(specific_regions[[8]][[1]][[2]], "data.frame")),
            aes(x=start, y=CNval), col='red', fill=NA)+
  geom_rect(data = all_plots_scDNA_sWGS_regions_ex3d_largebins[[1]][[1]],
            aes(xmin=start, xmax=end, ymin=0, ymax=CNval,group=cell), col='blue',
            fill=NA, alpha=0.1, size=0.1)+
  ggtitle(specific_regions[[8]][[2]][[1]]$title)#+lims(y=c(0,4))
pdo6specific

# give_subset_in_interval = function(gr_obj, start_int, end_int, chrom_int){
#   .region_subclonal = as(cbind.data.frame(chromosome=chrom_int, start=start_int,
#                                          end=end_int), "GRanges")
#   
#   .disj = disjoin(c(.region_subclonal, gr_obj),  with.revmap=TRUE, ignore.strand=TRUE)
#   .firstrevmap = sapply(.disj$revmap, `[`, 1)
#   .disj_subset = .disj[.firstrevmap == 1,] ## in interval
#   .disj_subset2 <- gr_obj[unlist(.disj_subset$revmap)[-1]-1,]
#   return(.disj_subset2)
# }

# subset_in_interval6 <- give_subset_in_interval(gr_obj = absCN_bed.granges$PDO6,
#                                                start_int = list_params[[6]]$start_AOI,
#                                                end_int = list_params[[6]]$end_AOI,
#                                                chrom_int=list_params[[6]]$chrom_AOI)
# subset_in_interval6_df <- data.frame(subset_in_interval6)
# subset_in_interval6_df$start = sapply(subset_in_interval6_df$start, function(i) max(i, list_params[[6]]$start_AOI))
# subset_in_interval6_df$end = sapply(subset_in_interval6_df$end, function(i) min(i, list_params[[6]]$end_AOI))
# 
# ggplot(subset_in_interval6_df %>% filter(cell < 50, event_confidence > 10))+
#   geom_rect(aes(xmin=start, xmax=end, ymin=0, ymax=CN_value,group=cell), col='black',
#             fill=NA, alpha=0.01, size=0.1)+
#   ggtitle("Chrom 4 in PDO6")+
#   lims(y=c(0,13))


ZWINT_PD06 <- list(segtab=segtables_organoids_absolute_copynumber$PDO6,
     PDO='PDO6',
     chrom_AOI=as.numeric(ag[ag$symbol == 'ZWINT','seq_name']),
     start_AOI=ag[ag$symbol == 'ZWINT','gene_seq_start'],
     end_AOI=ag[ag$symbol == 'ZWINT','gene_seq_end'],
     subclonal_line1=NA, subclonal_line2=NA, title='ZWINT in PDO6')
ZWINT_PD010 <- list(segtab=segtables_organoids_absolute_copynumber$PDO10,
                   PDO='PDO10',
                   chrom_AOI=as.numeric(ag[ag$symbol == 'ZWINT','seq_name']),
                   start_AOI=ag[ag$symbol == 'ZWINT','gene_seq_start'],
                   end_AOI=ag[ag$symbol == 'ZWINT','gene_seq_end'],
                   subclonal_line1=NA, subclonal_line2=NA, title='ZWINT in PD10')
MIEN1_PD014 <- list(segtab=segtables_organoids_absolute_copynumber$PDO14,
                    PDO='PDO14',
                    chrom_AOI=as.numeric(ag[ag$symbol == 'MIEN1','seq_name']),
                    start_AOI=ag[ag$symbol == 'MIEN1','gene_seq_start'],
                    end_AOI=ag[ag$symbol == 'MIEN1','gene_seq_end'],
                    subclonal_line1=NA, subclonal_line2=NA, title='MIEN1 in PD14')

give_params_gene <- function(gene, pdo){
  list(segtab=segtables_organoids_absolute_copynumber[[pdo]],
                      PDO=pdo,
                      chrom_AOI=as.numeric(ag[ag$symbol == gene,'seq_name']),
                      start_AOI=ag[ag$symbol == gene,'gene_seq_start'],
                      end_AOI=ag[ag$symbol == gene,'gene_seq_end'],
                      subclonal_line1=NA, subclonal_line2=NA, title=paste0(gene, ' in ', pdo))
}

# MYC_PD06 <- list(segtab=segtables_organoids_absolute_copynumber$PDO14,
#                     PDO='PDO6',
#                     chrom_AOI=as.numeric(ag[ag$symbol == 'ZWINT','seq_name']),
#                     start_AOI=ag[ag$symbol == 'MIEN1','gene_seq_start'],
#                     end_AOI=ag[ag$symbol == 'MIEN1','gene_seq_end'],
#                     subclonal_line1=NA, subclonal_line2=NA, title='MIEN1 in PD14')

# region_list_params = ZWINT_PD06
# abs_sc_GR=absCN_granges_perorg
# return_plot = F
# include_event_confidence = FALSE
all_plots_scDNA_sWGS_regions_ZWINT_PDO6 = give_plot_region_v2(region_list_params = ZWINT_PD06,
                                                         abs_sc_GR=absCN_granges_perorg,
                                                         return_plot = T, include_event_confidence = FALSE)
all_plots_scDNA_sWGS_regions_ZWINT_PDO6[[1]]
all_plots_scDNA_sWGS_regions_ZWINT_PDO6[[2]]

all_plots_scDNA_sWGS_regions_ZWINT_PDO10 = give_plot_region_v2(region_list_params = ZWINT_PD010,
                                                         abs_sc_GR=absCN_granges_perorg,
                                                         return_plot = T, include_event_confidence = FALSE, include_sc = FALSE)
all_plots_scDNA_sWGS_regions_ZWINT_PDO10[[1]]

all_plots_scDNA_sWGS_regions_MIEN1_PD014 = give_plot_region_v2(region_list_params = MIEN1_PD014,
                                                               abs_sc_GR=absCN_granges_perorg,
                                                               return_plot = T, include_event_confidence = FALSE, include_sc = FALSE)
all_plots_scDNA_sWGS_regions_MIEN1_PD014[[1]]

give_plot_region_v2(region_list_params = give_params_gene(gene = 'ZWINT', pdo = 'PDO2'),
                    abs_sc_GR=absCN_granges_perorg,
                    return_plot = T, include_event_confidence = FALSE)

give_plot_region_v2(region_list_params = give_params_gene(gene = 'ZWINT', pdo = 'PDO5'),
                    abs_sc_GR=absCN_granges_perorg,
                    return_plot = T, include_event_confidence = FALSE, include_sc = F)[[1]]

give_plot_region_v2(region_list_params = give_params_gene(gene = 'ZWINT', pdo = 'PDO6'),
                    abs_sc_GR=absCN_granges_perorg,
                    return_plot = T, include_event_confidence = FALSE)

ZWINT_PDO6 <- give_plot_region_v2(region_list_params = give_params_gene(gene = 'ZWINT', pdo = 'PDO6'),
                    abs_sc_GR=absCN_granges_perorg,
                    return_plot = F, include_event_confidence = FALSE)

ZWINT_neighbour_PDO6 <- give_plot_region_v2(region_list_params = list(segtab=segtables_organoids_absolute_copynumber[["PDO6"]],
                                                                     PDO="PDO6",
                                                                     chrom_AOI=10,
                                                                     start_AOI=58110001,
                                                                     end_AOI=58140000,
                                                                     subclonal_line1=NA, subclonal_line2=NA, title=paste0('Neighbouring 30kb region of ZWINT in PDO6')),
                                  abs_sc_GR=absCN_granges_perorg,
                                  return_plot = F, include_event_confidence = FALSE)
ZWINT_neighbour_PDO6[[2]]
segtables_organoids_absolute_copynumber$PDO6[seqnames(segtables_organoids_absolute_copynumber$PDO6) == 10,][10]

ag[ag$gene_name == 'ZWINT',]
x[grepl('ZWINT', x$Gene),]
ZWINT_PDO6[[2]]
# x[grepl('Gene', x$Gene),]
head(x[grepl('^10:', rownames(x)),][-(1:1555),])

max(as.numeric(segtables_organoids_absolute_copynumber$PDO6$segVal))
# [174]       10 58110001-58710000      * |    37.5732330396688


##############################

## plotting Chrom 8 in PDO3
chrom8_PDO3_sc = absCN_granges_perorg$PDO3[seqnames(absCN_granges_perorg$PDO3) == 8,]
chrom8_PDO3_sWGD <- segtables_organoids_absolute_copynumber$PDO3[seqnames(segtables_organoids_absolute_copynumber$PDO3) == 8,]
chrom8_PDO3_sWGD_df <- data.frame(chrom8_PDO3_sWGD)
chrom8_PDO3_sWGD_df$segVal = as.numeric(chrom8_PDO3_sWGD_df$segVal)
chrom8_PDO3_sWGD_df$start = as.numeric(chrom8_PDO3_sWGD_df$start)
chrom8_PDO3_sWGD_df$end = as.numeric(chrom8_PDO3_sWGD_df$end)
colnames(chrom8_PDO3_sWGD_df)[colnames(chrom8_PDO3_sWGD_df) == 'segVal'] = 'CN_value'

chrom8_PDO3_sc_df <- data.frame(chrom8_PDO3_sc)
chrom8_PDO3_sc_df$CN_value = as.numeric(chrom8_PDO3_sc_df$CN_value)
chrom8_PDO3_sc_df$cell = as.numeric(chrom8_PDO3_sc_df$cell)
chrom8_PDO3_sc_df$start = as.numeric(chrom8_PDO3_sc_df$start)
chrom8_PDO3_sc_df$end = as.numeric(chrom8_PDO3_sc_df$end)

length(chrom8_PDO3_sc)
length(chrom8_PDO3_sWGD)
chrom8_PDO3_sc

chrom8_PDO3_sWGD
plt1 <- ggplot(chrom8_PDO3_sWGD_df)+
  # geom_rect(aes(xmin=start, xmax=end, ymin=0, ymax=CNval,group=`names_cells[i]`), col='black',
  #           fill=NA, alpha=0.1, size=0.1)+
  geom_step(aes(x=start, y=CN_value),
            alpha=1)+
  theme_bw()+
  geom_hline(yintercept = 2, lty='dashed', col='blue')+
  labs(x='Genomic position of chromosome 8',y='Absolute CN value in PDO3')
plt1$theme$plot.margin[2] = unit(.8, "cm")
plt1
ggsave("../plots/chromothripsis_PDO3_onlyWGD.pdf", width = 5, height = 3)


ggplot(chrom8_PDO3_sc_df)+
  geom_step(aes(x=start, y=CN_value, group=cell,col=cell),
            alpha=1)+
  labs(x='Genomic position')

ggplot(rbind.data.frame(cbind.data.frame(chrom8_PDO3_sc_df[,c('start', 'end', 'CN_value', 'cell')], id='Single cell'),
                        cbind.data.frame(chrom8_PDO3_sWGD_df[,c('start', 'end', 'CN_value')], cell=-1, id='sWGS')))+
  geom_step(aes(x=start, y=CN_value, group=cell,col=id),
            alpha=1)+#facet_wrap(.~id, nrow=2)+
  labs(x='Genomic position of chromosome 8', y='Absolute CN value in PDO3')+
  labs(col='')+
  theme_bw()
  # guides(col=F)
ggsave("../plots/chromothripsis_PDO3.pdf", width = 4.5, height = 3)

## chromosome 8 for PDO9 (twin sample of PDO3)
absCN_granges_perorg$PDO9[seqnames(absCN_granges_perorg$PDO9) == 8,]
