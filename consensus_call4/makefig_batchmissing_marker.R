#!/usr/bin/Rscript
#Crosslink Copyright (C) 2016 NIAB EMR see included NOTICE file for details
#compare LOD-distance plot from recomb versus global optimisation

library(ggplot2)
library(gridExtra)
library(gtable)
library(png)

theme = theme(
  #panel.background = element_rect(fill="white"),
  #axis.ticks = element_line(colour=NA),
  #panel.grid = element_line(colour="grey"),
  axis.text.y = element_text(colour="black"),
  axis.text.x = element_text(colour="black",angle=45,hjust=1),
  text = element_text(size=8, family="Arial"),
  title = element_text(size=12, family="Arial")
)

inpfile = "figs/batch_missing_marker.csv"
outfile = "figs/batch_missing_marker.png"

dat = read.table(inpfile,col.names=c("sample","batch","missing","total"))

dat$pc_missing = dat$missing / dat$total

plt = ggplot(dat, aes(batch, pc_missing)) +
      geom_boxplot() +
      ylab("Prop. Missing Data") +
      xlab("Batch") +
      ggtitle("Proportion of Missing Data per Marker by Batch") + 
      theme

#      scale_y_log10() +   

ggsave(file=outfile,plot=plt,dpi=600)
