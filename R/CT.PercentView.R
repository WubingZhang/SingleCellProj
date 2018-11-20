CT.PercentView <- function(gg, ylab = "Count", pos = c("fill", "stack")){
  p = ggplot(gg, aes(Sample, weight=Count, fill=CellType))
  p = p + geom_bar(position = pos[1])
  p = p + labs(x = NULL, y=ylab)
  p = p + theme(plot.title = element_text(size=12))
  p = p + theme(text = element_text(colour="black",size = 14, family = "Helvetica"),
                plot.title = element_text(hjust = 0.5, size=18),
                axis.text = element_text(colour="gray10"))
  p = p + theme(axis.line = element_line(size=0.5, colour = "black"),
                axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
                panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.border = element_blank(), panel.background = element_blank())
  p
}
