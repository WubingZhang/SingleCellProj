#Plot square
viewCyto <- function(gg, title=NULL, x_cutoff = c(-1,1), y_cutoff=c(1.5,1.5),
                     filename = NULL, ...){
  requireNamespace("ggExtra", quietly=TRUE) || stop("need ggExtra package")
  gg$color = "Others"
  idx1 = gg$CD8_CD4 < x_cutoff[1]
  idx2 = gg$CD8_CD4 > x_cutoff[2]
  idx3 = gg$CYT < y_cutoff[1]
  idx4 = gg$CYT > y_cutoff[2]
  gg$color[idx1 & idx3] = "nonCytoCD4"
  gg$color[idx1 & idx4] = "CytoCD4"
  gg$color[idx2 & idx3] = "nonCytoCD8"
  gg$color[idx2 & idx4] = "CytoCD8"

  mycolor=c("Others"="aliceblue",  "CytoCD4"="#ff7f00", "nonCytoCD8"="#005CB7",
             "CytoCD8"="#984ea3", "nonCytoCD4"="#4daf4a" )
  p=ggplot(gg,aes(x=CD8_CD4, y=CYT, color=color))
  p=p+geom_point(shape=".", alpha=1/1, size = 1)
  p=p+scale_color_manual(values=mycolor)
  p=p+geom_jitter(size = 1)
  p=p+labs(title=title)
  p=p+geom_vline(xintercept = x_cutoff[x_cutoff>min(gg$CD8_CD4) & x_cutoff<max(gg$CD8_CD4)], linetype = "dotted")
  p=p+geom_hline(yintercept = y_cutoff[y_cutoff>min(gg$CYT) & y_cutoff<max(gg$CYT)], linetype = "dotted")
  p=p+annotate("text",color="red",x=min(gg$CD8_CD4), y=max(gg$CYT),hjust = 0,
               label=paste0("CytoCD4: ", length(which(gg$color=="CytoCD4"))))
  p=p+annotate("text",color="red",x=min(gg$CD8_CD4), y=min(gg$CYT), hjust = 0,
               label=paste0("nonCytoCD4: ", length(which(gg$color=="nonCytoCD4"))))
  p=p+annotate("text",color="red",x=max(gg$CD8_CD4), y=max(gg$CYT),hjust = 1,
               label=paste0("CytoCD8: ", length(which(gg$color=="CytoCD8"))))
  p=p+annotate("text",color="red",x=max(gg$CD8_CD4), y=min(gg$CYT), hjust = 1,
               label=paste0("nonCytoCD8: ", length(which(gg$color=="nonCytoCD8"))))
  p = p + theme(legend.key = element_rect(fill = "transparent"))
  p = p + theme(text = element_text(colour="black",size = 14, family = "Helvetica"),
                plot.title = element_text(hjust = 0.5, size=18),
                axis.text = element_text(colour="gray10"))
  p = p + theme(axis.line = element_line(size=0.5, colour = "black"),
                panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.border = element_blank(), panel.background = element_blank())
  p = p + theme(legend.position = "none")
  p=suppressWarnings(ggExtra::ggMarginal(p, type="histogram",bins=200,color="#0062A7"))
  p$data = gg

  if(!is.null(filename)){
    ggsave(plot=p, filename=filename, units = "in", dpi=600, ...)
  }
  return(p)
}
