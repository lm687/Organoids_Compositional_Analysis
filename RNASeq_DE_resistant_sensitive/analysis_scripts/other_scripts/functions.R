give_linebreaks = function(strng, max_line){
  splt = strsplit(strng, " ")[[1]]
  if(length(splt)>1){
    which_good = which(sapply(1:length(splt), function(i) nchar(paste0(splt[1:i], collapse=' ')) ) < max_line)
    which_good = which_good[length(which_good)]
    if(length(which_good)==0){
      return(paste0(paste0(splt[1], collapse=" "), "\n",
                    give_linebreaks(paste0(splt[2:length(splt)], collapse=" "), max_line)))
    }else{
      if(which_good == length(splt)){
        return(paste0(splt[1:which_good], collapse=" "))
      }else{
        ## continue
        return(paste0(paste0(splt[1:which_good], collapse=" "), "\n",
                      give_linebreaks(paste0(splt[(which_good+1):length(splt)], collapse=" "), max_line)))
      }
    }
  }else{
    return(strng)
  }
}