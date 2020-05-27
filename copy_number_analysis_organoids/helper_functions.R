#' Create the barplot given a matrix of exposures
createBarplot <- function(matrix_exposures, angle_rotation_axis = 0, order_labels=NULL,
                          remove_labels=FALSE, levels_signatures=NULL, includeMelt=NULL, Melt=NULL, verbose=TRUE){
  #' error due to it not being a matrix:  
  #'    No id variables; using all as measure variables
  #'    Rerun with Debug
  #'    Error in `[.data.frame`(.mat, , "Var1") : undefined columns selected \\
  #' (use tomatrix())
  
  
  if(is.null(colnames(matrix_exposures))) stop('columns must have names')
  require(reshape2)
  require(ggplot2)
  library(RColorBrewer)
  if(verbose){
    cat(paste0('Creating plot... it might take some time if the data are large. Number of samples: ', nrow(matrix_exposures), '\n'))
  }
  
  if( (!is.null(order_labels) & typeof(order_labels) == "logical")){if(!order_labels){cat('WARNING: Order labels is either a vector with desired order or NULL, not bool')}}
  
  if(!is.null(levels_signatures)){
    library(ggthemes)
    ggthemes_data$economist
  }else{
    levels_signatures <- colnames(matrix_exposures)
  }
  n <- 60
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unique(unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))))
  
  myColors <- col_vector[1:length(levels_signatures)]
  names(myColors) <- levels_signatures
  if(is.null(order_labels)) order_labels = rownames(matrix_exposures)
  if(!is.null(includeMelt)){
    cat('For whatever reason sometimes the melt does not work. Here it is passed as argument.')
    .mat <- Melt
  }else{
    .mat <- melt(matrix_exposures)
  }
  .mat[,'Var1'] <- factor(.mat[,'Var1'], levels=order_labels)
  .mat[,'Var2'] <- factor(.mat[,'Var2'], levels=levels_signatures)
  ###rownames(.mat) <- rownames(matrix_exposures) ### new
  
  if(!remove_labels){
    if(!is.null(levels_signatures)){
      ggplot(.mat, aes(x=Var1, y=value, fill=factor(Var2, levels=levels_signatures)))+
        geom_bar(stat = 'identity')+
        theme(axis.text.x = element_text(angle = angle_rotation_axis, hjust = 1))+
        #theme(axis.title.x=element_blank(),
        #      axis.text.x=element_blank(),
        #      axis.ticks.x=element_blank())+
        guides(fill=guide_legend(title="Signature"))+
        scale_fill_manual(name = "grp",values = myColors)
    }else{
      ggplot(.mat, aes(x=Var1, y=value, fill=factor(Var2, levels=levels_signatures)))+
        geom_bar(stat = 'identity')+
        theme(axis.text.x = element_text(angle = angle_rotation_axis, hjust = 1))+
        #theme(axis.title.x=element_blank(),
        #      axis.text.x=element_blank(),
        #      axis.ticks.x=element_blank())+
        guides(fill=guide_legend(title="Signature"))
    }
  }else{
    if(!is.null(levels_signatures)){
      ggplot(.mat, aes(x=Var1, y=value, fill=factor(Var2, levels=levels_signatures)))+
        geom_bar(stat = 'identity')+
        theme(axis.text.x = element_text(angle = angle_rotation_axis, hjust = 1))+
        theme(axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank())+
        guides(fill=guide_legend(title="Signature"))+
        scale_fill_manual(name = "grp",values = myColors)
    }else{
      ggplot(.mat, aes(x=Var1, y=value, fill=factor(Var2, levels=levels_signatures)))+
        geom_bar(stat = 'identity')+
        theme(axis.text.x = element_text(angle = angle_rotation_axis, hjust = 1))+
        theme(axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank())+
        guides(fill=guide_legend(title="Signature"))
    }
  }
}


give_joint_histogram = function(some_list, bins=30, no_colour=FALSE){
  require(reshape2)
  require(ggplot2)
  .x = reshape2::melt(some_list)
  if(no_colour){
    ggplot(.x, aes(x=value, group=L1), alpha=0.2)+geom_histogram(bins = bins)+facet_wrap(.~L1, ncol=1, scales='free_y')
  }else{
    ggplot(.x, aes(x=value, fill=L1, group=L1), alpha=0.2)+geom_histogram(bins = bins)+facet_wrap(.~L1, ncol=1, scales='free_y')
  }
}


extractCopynumberFeatures = function(CN_data, cores = 1)
{
  #get chromosome lengths
  chrlen<-read.table(paste(this_path,"data/hg19.chrom.sizes.txt",sep="/"),sep="\t",stringsAsFactors = F)[1:24,]
  
  #get centromere locations
  gaps<-read.table(paste(this_path,"data/gap_hg19.txt",sep="/"),sep="\t",header=F,stringsAsFactors = F)
  centromeres<-gaps[gaps[,8]=="centromere",]
  
  if(cores > 1) {
    require(foreach)
    doMC::registerDoMC(cores)
    
    temp_list = foreach::foreach(i=1:6) %dopar% {
      if(i == 1){
        list(segsize = getSegsize(CN_data) )
      } else if (i == 2) {
        list(bp10MB = getBPnum(CN_data,chrlen) )
      } else if (i == 3) {
        list(osCN = getOscilation(CN_data,chrlen) )
      } else if (i == 4) {
        list(bpchrarm = getCentromereDistCounts(CN_data,centromeres,chrlen) )
      } else if (i == 5) {
        list(changepoint = getChangepointCN(CN_data) )
      } else {
        list(copynumber = getCN(CN_data) )
      }
      
    }
    unlist( temp_list, recursive = FALSE )
  } else {  
    
    segsize<-getSegsize(CN_data)
    bp10MB<-getBPnum(CN_data,chrlen)
    osCN<-getOscilation(CN_data,chrlen)
    bpchrarm<-getCentromereDistCounts(CN_data,centromeres,chrlen)
    changepoint<-getChangepointCN(CN_data)
    copynumber<-getCN(CN_data)
    
    list(segsize=segsize,bp10MB=bp10MB,osCN=osCN,bpchrarm=bpchrarm,changepoint=changepoint,copynumber=copynumber)
  }
  
}


getSampNames<-function(abs_profiles)
{
  if(class(abs_profiles)=="QDNAseqCopyNumbers")
  {
    samps<-colnames(abs_profiles)
  }
  else
  {
    samps<-names(abs_profiles)
  }
  samps
}

getBPnum2<-function(abs_profiles,chrlen){
  out<-c()
  cln = colnames(abs_profiles[[1]])
  abs_profiles = do.call('rbind', abs_profiles)
  colnames(abs_profiles) = cln
  # samps<-getSampNames(abs_profiles)
  # for(i in samps)
  # {
  segTab = abs_profiles#[[i]]
  chrs<-unique(segTab[,'chromosome'])
  allBPnum<-c()
  for(c in chrs)
  {
    currseg<-segTab[segTab[,'chromosome']==c,]
    intervals<-seq(1,chrlen[chrlen[,1]== c,2]+10000000,10000000)
    if(is.null(nrow(currseg))){
      res <- hist(as.numeric(currseg['end'][-1]),breaks=intervals,plot=FALSE)$counts
    }else{
      res <- hist(as.numeric(currseg[,'end'][-nrow(currseg)]),breaks=intervals,plot=FALSE)$counts
    }
    allBPnum<-c(allBPnum,res)
  }
  out<-rbind(out,cbind(ID=rep(i,length(allBPnum)),value=allBPnum))
  # out
  # }
  # rownames(out)<-NULL
  data.frame(out,stringsAsFactors = F)
}

wrapper_get_bp = function(df){
  for(i in 1:length(df)){
    df[[i]][,'chromosome'] = paste0('chr', df[[i]][,'chromosome'])
  }
  names(df) = paste0('df', 1:length(df))
  getBPnum2(df, chrlen = chrlen)
}