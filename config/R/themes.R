# custom.levels<-names(pretty.level.names)
# gg.labs.name<-"Age group"
# gg.title.groups<-"age groups"
mag.taxonomy.barplot.theme <- theme(axis.text.x = element_text(size=13),
                                      axis.text.y = element_text(size=13), # size of y axis ticks
                                      axis.title = element_text(size = 15), # size of axis names
                                      plot.caption = element_text(size=13), # size of plot caption
                                      legend.title = element_text(size = 15), # size of legend title
                                    legend.position = "right",
                                    panel.grid.minor = element_blank(),
                                    panel.grid.major = element_blank()
)

#' ## Plot MAG genome sizes. 
mag.boxplot.theme <- theme(axis.text.x = element_text(angle=45,size=14,hjust=1),# rotate 
      # the x-axis labels by 45 degrees and shift to the right
      axis.title = element_text(size = 14), # size of axis names
      plot.title = element_text(size = 10), 
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank()
)

diffabund.plot.theme <-theme(
  # plot.margin=unit(c(1,1,1,2), 'cm'),
  axis.title = element_text(size = 10),
  axis.title.y = element_text(size = 12),
  axis.text.y = ggtext::element_markdown(size=10),
  axis.text.x = element_text(size=10),
  strip.text.x = ggtext::element_markdown(size=10),
  plot.title = element_text(size = 13),
  legend.text = element_text(size = 10),
  legend.title = element_text(size = 12),
  legend.position = "right",
  panel.grid.minor = element_blank(),
  panel.grid.major = element_blank())
