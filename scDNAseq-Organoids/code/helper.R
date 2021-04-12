clean_chrom = function(i){
  sapply(i, function(j){
    if(j %in% c('X', 'Y')){
      j
    }else{
      ## autosomal
      substr(j, 2, 1000)
    }
  })
}