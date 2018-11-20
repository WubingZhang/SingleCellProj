#' Create Seurat object
#'
#' Run Seurat pipeline and reture Seurat object
#'
#' @docType methods
#' @name SeuratPipeline
#' @rdname SeuratPipeline
#' @aliases seuratpipeline
#'
#' @param tpmMat Matrix like data object, or a file path of data.
#' @param proj Project name (string).
#' @param min.c Include genes with detected expression in at least this many cells.
#' Will subset the raw.data matrix as well. To reintroduce excluded genes, create a
#' new object with a lower cutoff.
#' @param min.g Include cells where at least this many genes are detected.
#' @param normalization.method Method for cell normalization. Default is no normalization.
#' In this case, run NormalizeData later in the workflow. As a shortcut, you can specify
#' a normalization method (i.e. LogNormalize) here directly.
#' @param do.scale In object@scale.data, perform row-scaling (gene-based z-score). FALSE
#' by default. In this case, run ScaleData later in the workflow. As a shortcut, you can
#' specify do.scale = TRUE (and do.center = TRUE) here.
#' @param do.center In object@scale.data, perform row-centering (gene-based centering).
#' @param vars.to.regress Variables to be regressed.
#' @param outdir The output directory.
#' @param ... Other useful paramters in method `CreateSeuratObject`.
#'
#' @return A Seurat object.
#'
#' @author Wubing Zhang
#'
#' @examples
#'
#' @import Seurat
#'
#' @export
#'
SeuratPipeline <- function(tpmMat, proj, min.c = 3, min.g = 500, max.g = 20000,
                           normalization.method = NULL, do.scale = TRUE, do.center = FALSE,
                           vars.to.regress = c("nUMI"), outdir = ".", ...)
{
  date = format(Sys.time(), "%b.%d.%Y")
  #=========QC========
  message("Check gene and cell coverage ...")
  nGene = apply(tpmMat, 2, function(x) length(x[x>0]))
  nCell = apply(tpmMat, 1, function(x) length(x[x>0]))
  pdf(file.path(outdir, paste0(proj,"_QC_Coverage_", date, ".pdf")),width=8,height=4.5)
  par(mfrow=c(1,2))
  plot(1:ncol(tpmMat),sort(nGene),pch=16,col="blue",ylab="Number of Genes Expressed",xlab="Cells",main="Cell Filter")
  abline(h=min.g,lwd=2,lty=2);text(ncol(tpmMat)/2,min.g+max(nGene)*0.05,paste0("n = ",min.g))
  plot(1:nrow(tpmMat),sort(nCell),pch=16,col="blue",ylab="Number of Cells Expressed",xlab="Genes",main="Gene Filter")
  abline(h=min.c,lwd=2,lty=2);text(nrow(tpmMat)/2,min.c+max(nCell)*0.05,paste0("n = ",min.c))
  dev.off()

  SeuratObj <- CreateSeuratObject(tpmMat, project = proj, min.cells = min.c, min.genes = min.g, ...)
  mito.genes <- grep("^MT-", rownames(SeuratObj@data), value = TRUE, ignore.case = TRUE)
  ercc.genes <- grep("^ERCC", rownames(SeuratObj@data), value = TRUE, ignore.case = TRUE)
  percent.mito <- colSums(SeuratObj@data[mito.genes, ])/colSums(SeuratObj@data)
  percent.ercc <- colSums(SeuratObj@data[ercc.genes, ])/colSums(SeuratObj@data)
  SeuratObj <- AddMetaData(SeuratObj, percent.mito, "percent.mito")
  SeuratObj <- AddMetaData(SeuratObj, percent.ercc, "percent.ercc")
  p1<-VlnPlot(SeuratObj, c("percent.mito"), nCol = 1)
  ggsave(file.path(outdir, paste0(proj,"_QC_Spikein_", date, ".pdf")), p1, width = 6, height = 4.5)

  #=========Filter========
  message("Filter cells and find variable genes ...")
  SeuratObj <- SubsetData(SeuratObj, subset.name = "percent.mito", accept.high = 0.05)
  SeuratObj <- SubsetData(SeuratObj, subset.name = "percent.ercc", accept.high = 0.05)
  SeuratObj <- FilterCells(object = SeuratObj, subset.names = c("nGene"),
                           low.thresholds = min.g, high.thresholds = max.g)

  SeuratObj <- NormalizeData(object = SeuratObj, normalization.method = normalization.method, scale.factor = 10000)
  SeuratObj <- FindVariableGenes(object = SeuratObj, mean.function = ExpMean, dispersion.function = LogVMR,
                                 x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, do.plot = FALSE)
  if(do.scale) SeuratObj <- ScaleData(object = SeuratObj, vars.to.regress = vars.to.regress, do.center = do.center)

  #=========PCA===========
  message("PCA analysis ...")
  SeuratObj <- RunPCA(object = SeuratObj, pc.genes = SeuratObj@var.genes,
                      do.print = FALSE, rev.pca = TRUE)
  # VizPCA(object = SeuratObj, pcs.use = 1:2)
  # PCAPlot(object = SeuratObj, dim.1 = 1, dim.2 = 2)
  # PCHeatmap(object = SeuratObj, pc.use = 1:12, cells.use = 500, do.balanced = TRUE,
  #           label.columns = FALSE, use.full = FALSE)
  # SeuratObj <- JackStraw(object = SeuratObj, num.replicate = 100, display.progress = FALSE)
  # JackStrawPlot(object = SeuratObj, PCs = 1:12)
  p3 = PCElbowPlot(object = SeuratObj)
  ggsave(file.path(outdir, paste0(proj,"_PCElbowPlot_", date, ".pdf")), p3, width = 5, height = 4)

  return(SeuratObj)
}

