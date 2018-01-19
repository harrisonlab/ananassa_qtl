#!/usr/bin/Rscript

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

inpfile = "figs/popn_stats.csv"
outfile = "figs/popn_stats.png"
outfile2 = "figs/popn_stats_usable.png"

dat = read.table(inpfile,col.names=c("popn","class","freq"))

dat$class <- factor(dat$class,
    levels = c("Maternal","Paternal","Shared","ImpossibleGenotype","MissingData","SegDistorted","MissingParental","Monomorphic"))

#, order = -as.numeric(class)
plt = ggplot(dat, aes(popn, fill=class, weight=freq)) +
      geom_bar() +
      ylab("Markers") +
      xlab("Population") +
      ggtitle("Marker Classification by Population") +
      scale_fill_discrete(guide=guide_legend(reverse=TRUE)) + 
      theme

ggsave(file=outfile,plot=plt,dpi=600)

dat = subset(dat, class=='Maternal'|class=='Paternal'|class=='Shared')
dat$class <- factor(dat$class, levels = c("Maternal","Paternal","Shared"))

plt = ggplot(dat, aes(popn, fill=class, weight=freq)) +
      geom_bar() +
      ylab("Markers") +
      xlab("Population") +
      ggtitle("Usable Markers by Population") +
      scale_fill_discrete(guide=guide_legend(reverse=TRUE)) + 
      theme

ggsave(file=outfile2,plot=plt,dpi=600)
