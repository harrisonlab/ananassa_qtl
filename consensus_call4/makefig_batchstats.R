#!/usr/bin/Rscript
#Crosslink Copyright (C) 2016 NIAB EMR see included NOTICE file for details
#compare LOD-distance plot from recomb versus global optimisation

library(ggplot2)
library(gridExtra)
library(gtable)
library(png)

#args = commandArgs(trailingOnly = TRUE)
#inpfile = args[1]  
#outfile = args[2]

inpfile = "figs/batch_stats.csv"
outfile = "figs/batch_stats.png"

theme = theme(
  #panel.background = element_rect(fill="white"),
  #axis.ticks = element_line(colour=NA),
  #panel.grid = element_line(colour="grey"),
  axis.text.y = element_text(colour="black"),
  axis.text.x = element_text(colour="black",angle=45,hjust=1),
  text = element_text(size=8, family="Arial"),
  title = element_text(size=12, family="Arial")
)

dat = read.table(inpfile,col.names=c("batch","class","freq"))

dat$class <- factor(dat$class,
    levels = c("PolyHighResolution","NoMinorHom","MonoHighResolution","OTV","CallRateBelowThreshold","Other"))
    
plt = ggplot(dat, aes(batch, fill=class, weight=freq)) +
      geom_bar() +
      ylab("Markers") +
      xlab("Batch") +
      ggtitle("Marker Classifications by Batch") + 
      scale_fill_discrete(guide=guide_legend(reverse=TRUE)) + 
      theme

ggsave(file=outfile,plot=plt,dpi=600)
