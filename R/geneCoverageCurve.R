geneCoverageView <- function(count, isRowGene=TRUE){
  count = as.matrix(count)
  count[count!=0] = 1
  if(isRowGene){
    tmp = rowSums(count)
    coverage = sapply(0:ncol(count), function(x) length(which(tmp>x)))
    res = data.frame(Cutoff=0:ncol(count), Number=coverage)
  }
  else{
    tmp = colSums(count)
    coverage = sapply(0:nrow(count), function(x) length(which(tmp>x)))
    res = data.frame(Cutoff=0:nrow(count), Number=coverage)
  }
  p = ggplot(res,aes(Cutoff, Number))
  p = p + geom_line()
  p = p + theme_bw(12)
  p = p + labs(x="Number of samples", y="Number of expressed genes")
  p
  return(p)
}
