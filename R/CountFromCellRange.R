CountFromCellRange <- function(pipestance_path = "PBMC4k/"){
  GeneBCMatrix <- load_cellranger_matrix(pipestance_path)
  expr = exprs(GeneBCMatrix) # expression matrix
  return(expr)
}

