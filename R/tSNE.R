#' tSNE
#'
#' tSNE plot
#'
#' @docType methods
#' @name tSNE
#' @rdname tSNE
#' @aliases tsne
#'
#' @param SeuratObj A Seurat object.
#' @param dims.use Which dimensions to use as input features.
#' @param org "hsa" or "mmu".
#' @param outdir The output directory.
#'
#' @return A list contrains two objects, including \code{SeuratObj} and \code{markers}.
#'
#' @author Wubing Zhang
#'
#' @examples
#'
#' @import Seurat
#' @import dplyr
#'
#' @export

tSNE <- function(SeuratObj, dims.use = 1:10, org = "hsa", outdir = ".")
{
  date = format(Sys.time(), "%b.%d.%Y")
  #=========tSNE===========
  message("t-SNE analysis ...")
  SeuratObj <- FindClusters(object = SeuratObj, reduction.type = "pca", dims.use = dims.use,
						    resolution = 0.6, print.output = 0, save.SNN = TRUE)
  SeuratObj <- RunTSNE(object = SeuratObj, dims.use = dims.use, do.fast = TRUE)
  png(file.path(outdir, paste0(SeuratObj@project.name, "_origIdent_tSNE_", date, ".png")),
	  res=300, width=6, height=4, units = "in")
  TSNEPlot(object = SeuratObj, do.label = TRUE, pt.size = 0.5, group.by = "res.0.6")
  dev.off()

  #=========identify marker===========
  message("Find marker genes ...")
  pbmc.markers <- NULL
  pbmc.markers <- FindAllMarkers(object = SeuratObj, only.pos = TRUE, min.pct = 0.1)
  saveRDS(pbmc.markers, file.path(outdir, paste0(SeuratObj@project.name, "_DiffMarkers_", date, ".rds")))
  saveRDS(SeuratObj, file.path(outdir, paste0(SeuratObj@project.name, "_SeuratObj_", date, ".rds")))

  idx = c(sapply(levels(pbmc.markers$cluster), function(x) which(pbmc.markers$cluster==x)[1:50]))
  pbmc.markers <- pbmc.markers[idx, ]
  current.cluster.ids = as.integer(levels(pbmc.markers$cluster))
  new.cluster.ids = AssignCellType(pbmc.markers, org = org)
  SeuratObj@meta.data$assign.ident = SeuratObj@ident[rownames(SeuratObj@meta.data)]
  SeuratObj@meta.data$assign.ident <- plyr::mapvalues(x = SeuratObj@meta.data$assign.ident,
											 from = current.cluster.ids, to = new.cluster.ids)

  png(file.path(outdir, paste0(SeuratObj@project.name, "_assignIdent_tSNE_", date, ".png")),res=300,
                                                         width=6, height=4, units = "in")
  TSNEPlot(object = SeuratObj, do.label = TRUE, pt.size = 0.5, group.by = "assign.ident")
  dev.off()
  return(list(SeuratObj=SeuratObj, markers = pbmc.markers))
}

