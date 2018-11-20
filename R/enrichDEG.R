enrichDEG <- function(res, prefix = NULL, width = 6, height=3.5){
  require(MAGeCKFlute)
  universe = TransGeneID(rownames(res))
  # #=======GSEA============
  # genelist = -log2(res$adj.P.Val)
  # tmp = res$logFC
  # tmp[tmp>0] = 1; tmp[tmp<0] = -1;
  # genelist = genelist * tmp
  # genelist = genelist[!is.na(universe)]
  # names(genelist) = universe[!is.na(universe)]
  # GSEA.res = enrich.GSE(genelist)
  # p1 = EnrichedGSEView(GSEA.res, color = "#e41a1c")
  ##===ORT enrichment======
  universe = universe[!is.na(universe)]
  upgene = rownames(res)[res$logFC>0.1 & res$adj.P.Val<0.05]
  upgene = TransGeneID(upgene)
  upgene = upgene[!is.na(upgene)]
  if(!is.null(upgene)){
    enrichUp = enrich.ORT(upgene, universe)
    p.up = EnrichedView(enrichUp@result, color = "#e41a1c", termNum = 10, charLength = 30)
  }else{
    enrichUp = NULL
    p.up = NULL
  }

  downgene = rownames(res)[res$logFC< -0.1 & res$adj.P.Val<0.05]
  downgene = TransGeneID(downgene)
  downgene = downgene[!is.na(downgene)]
  if(!is.null(downgene)){
    enrichDown = enrich.ORT(downgene, universe)
    p.down = EnrichedView(enrichDown@result, color = "#377eb8", termNum = 10, charLength = 30)
  }else{
    enrichDown = NULL
    p.down = NULL
  }

  res = list(p.up = p.up, p.down = p.down, enrichUp = enrichUp, enrichDown = enrichDown)
  if(!is.null(prefix)){
    if(!is.null(p.up))
      ggsave(plot=p.up, filename=paste0(prefix, "_ORTup.png"),
             units = "in", dpi=600, width = width, height = height)
    if(!is.null(p.down))
      ggsave(plot=p.down, filename=paste0(prefix, "_ORTdown.png"),
             units = "in", dpi=600, width = width, height = height)
  }
  return(res)
}
