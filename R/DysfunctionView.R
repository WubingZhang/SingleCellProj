#' Cell surface expression of TIM-3 and PD1 can be used to partition CD8+
#' T cells into three different groups: double negative, double positive,
#' and single positive.
#'
DysfunctionView <- function(expr, xcutoff, ycutoff, fill="CYT", title=NULL, filename = NULL, ...){
  gg = as.data.frame(expr)
  idx = gg$HAVCR2==0 | gg$PDCD1==0
  gg = gg[!idx,]
  # xcutoff = c(median(gg$HAVCR2)-scale*sd(gg$HAVCR2),
  #             median(gg$HAVCR2)+scale*sd(gg$HAVCR2))
  # ycutoff = c(median(gg$PDCD1)-scale*sd(gg$PDCD1),
  #             median(gg$PDCD1)+scale*sd(gg$PDCD1))
  gg$Group = "Others"
  gg$Group[gg$HAVCR2>xcutoff[2]&gg$PDCD1>ycutoff[2]] = "DP"
  gg$Group[gg$HAVCR2<xcutoff[1]&gg$PDCD1<ycutoff[1]] = "DN"
  gg$Group[gg$HAVCR2<xcutoff[1]&gg$PDCD1>ycutoff[2]] = "SP.PD1"
  gg$Group[gg$HAVCR2>xcutoff[2]&gg$PDCD1<ycutoff[1]] = "SP.TIM3"

  mycolour=c("Others"="aliceblue",  "DP"="#ff7f00", "DN"="#4daf4a",
             "SP.PD1"="#984ea3", "SP.TIM3"="SlateBlue" )
  p=ggplot(gg, aes(x=HAVCR2, y=PDCD1))
  p=p+geom_point(aes(color=Group, fill=gg[,fill]), alpha=1, shape=21, size = 1)
  p=p+scale_color_manual(values = mycolour, guide=FALSE)
  p=p+scale_fill_gradient(low = "gray90", high = "blue")
  p=p+geom_hline(yintercept = ycutoff, size=0.5, linetype="dashed")
  p=p+geom_vline(xintercept = xcutoff, size=0.5, linetype="dashed")
  p=p+annotate("text",color="#ff7f00",x=max(gg$HAVCR2), y=max(gg$PDCD1),hjust = 1, vjust=1,
               label=paste0("DP: ", length(which(gg$Group=="DP"))))
  p=p+annotate("text",color="#4daf4a",x=min(gg$HAVCR2), y=min(gg$PDCD1),hjust = 0, vjust=0,
               label=paste0("DN: ", length(which(gg$Group=="DN"))))
  p=p+annotate("text",color="#984ea3",x=min(gg$HAVCR2), y=max(gg$PDCD1),hjust = 0, vjust=1,
               label=paste0("SP.PD1: ", length(which(gg$Group=="SP.PD1"))))
  p=p+annotate("text",color="SlateBlue",x=max(gg$HAVCR2), y=min(gg$PDCD1),hjust = 1, vjust=0,
               label=paste0("SP.TIM3: ", length(which(gg$Group=="SP.TIM3"))))
  p=p+labs(x="TIM3", y="PD1", fill=fill, title=title)

  if(!is.null(filename)){
    ggsave(plot=p, filename=filename, units = "in", dpi=600, ...)
  }
  return(p)
}
