#' AssignCellType
#'
#' Assign Cell Type for single cells
#'
#' @docType methods
#' @name AssignCellType
#' @rdname AssignCellType
#' @aliases assigncelltype
#'
#' @param markers A data frame which includes differential expressed genes in each cell cluster.
#' @param org "hsa" or "mmu".
#' @param celltype.markers A list, which customize marker of all different cell type.
#'
#' @return A Vector of assigned cell type names.
#'
#' @author Wubing Zhang
#'
#' @examples
#'
#' @export
#'
AssignCellType <- function(markers, org = "hsa", celltype.markers = NA){
  if(is.na(celltype.markers)){
    data("Merged_CellMarker")
    celltype.markers = Merged_CellMarker
  }

  if(org == "mmu"){
    data("HumanGene_MouseGene_Mapping")
    geneMap = HumanGene_MouseGene_Mapping
    celltype.markers <- lapply(celltype.markers, function(x){
      geneMap[geneMap$HGNC.symbol%in%x, "MGI.symbol"]
    })
  }
  CT = sapply(as.integer(levels(markers$cluster)), function(x){
    idx = markers$cluster==x
    tmpFC = markers$avg_logFC[idx]
    names(tmpFC) = markers$gene[idx]
    score_cluster = sapply(celltype.markers, function(y){
      score = sum(tmpFC[y], na.rm = TRUE) / log2(length(y)+1)
      return(score)
    })
    if(max(score_cluster, na.rm = TRUE)>0.6)
      CT = names(score_cluster)[which.max(score_cluster)]
    else
      CT = "Others"
  })
  return(CT)
}
