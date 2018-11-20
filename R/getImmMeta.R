getImmMeta <- function(mat, filter = TRUE){
  immGenes = list(CD3 = c('CD3D','CD3E','CD3G'),
                  CYT = c('GZMA','GZMB','PRF1'),
                  CD8 = c('CD8A','CD8B'),
                  CD4 = c('CD4'),
                  T0 = c('IL2','IL3','IL4','IL5','IL6','IL10','CSF2','IFNG','LTA','TNF'),
                  T1 = c('IFNG','TNF','LTA','IL3','CSF2','IL2'),
                  T2 = c('IL4','IL5','IL6','IL10','IL13'))
  CDres <- lapply(immGenes, function(x){
    x = intersect(x, rownames(mat))
    colMeans(mat[x, , drop = FALSE])
  })
  CDres = as.data.frame(CDres)
  CDres$CD8_CD4 = CDres$CD8 - CDres$CD4
  #' Select cells with high CD3 expression(T cell marker), non-zero expression
  #' of CYT, and unequal expression of CD8 and CD4
  if(filter){
    cd3Idx1 =  CDres$CD3 > mean(CDres$CD3)
    CDres = CDres[cd3Idx1, ]
    cytIdx1 = CDres$CYT != 0
    CDres = CDres[cytIdx1, ]
    cd8mcd4Idx1 = CDres$CD8_CD4 != 0
    CDres = CDres[cd8mcd4Idx1, ]
  }

  return(CDres)
}
