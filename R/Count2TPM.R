#' Transform raw count to TPM
#'
#' @docType methods
#' @name Count2TPM
#' @rdname Count2TPM
#'
#' @param countMat A matrix-like object.
#' @param idType Ensembl, Symbol, or Entrez.
#' @param org One of "hsa", "mmu", "bta", "cfa", "ptr", "rno" and "ssc".
#'
#' @return A TPM matrix.
#'
#' @examples
#'
#' @import biomaRt
#' @export

Count2TPM <- function(countMat, idType = "Ensembl", org="hsa")
{
  datasets = paste0(c("hsapiens", "mmusculus", "btaurus", "cfamiliaris",
                      "ptroglodytes", "rnorvegicus", "sscrofa"), "_gene_ensembl")
  type = c("ensembl_gene_id", "entrezgene", "hgnc_symbol", "start_position", "end_position")
  if(org=="mmu") type[3] = "mgi_symbol"
  # listEnsemblArchives()
  # listMarts()
  # listAttributes()
  ds = datasets[grepl(org, datasets)]
  mart <- useMart(host = "www.ensembl.org", biomart = 'ENSEMBL_MART_ENSEMBL', dataset = ds)
  ensembl = getBM(attributes=type, mart = mart)
  ensembl$Length <- abs(ensembl$end_position - ensembl$start_position)
  if(toupper(idType) == "ENSEMBL"){
    tmp = ensembl[match(rownames(countMat),ensembl$ensembl_gene_id), 3]
    countMat = countMat[!duplicated(tmp), ]
    rownames(countMat) = tmp[!duplicated(tmp)]
    len <- ensembl[match(rownames(countMat),ensembl[,3]), "Length"]
  }
  else if(toupper(idType) == "SYMBOL")
    len <- ensembl[match(rownames(countMat), ensembl[,3]), "Length"]
  else if(toupper(idType) == "ENTREZ")
    len <- ensembl[match(rownames(countMat), ensembl[,2]), "Length"]
  else
    stop("Please input right type of gene name, such as Ensembl or gene Symbol ...")

  na_idx = which(is.na(len))
  if(length(na_idx)>0){
    warning(paste0("Omit ", length(na_idx), " genes of which length is not available !"))
    countMat = countMat[!is.na(len),]
    len = len[!is.na(len)]
  }
  tmp <- countMat / len
  TPM <- 1e6 * t(t(tmp) / colSums(tmp))
  TPM = TPM[!duplicated(rownames(TPM)),]
  return(TPM)
}
