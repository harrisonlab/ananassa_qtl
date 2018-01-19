#!/usr/bin/Rscript

#
# synteny plot of consensus map versus all biparental maps
# run from ~/octoploid_mapping/consensus_map4/consensus_map
#

library(ggplot2)

#construct list of linkage group names
hg <- 1:7
sg <- c('A','B','C','D')
#lglist <- character()
#for (x in hg) for (y in sg) lglist[length(lglist)+1] <- paste0(x,y)

#blank starting dataframe
df <- data.frame(marker=character(),lg=character(),hg=character(),sg=character(),conpos=numeric(),mapname=character(),mappos=numeric(),parent=character())

for (x in hg)
{
    for (y in sg)
    {
        lg <- paste0(x,y)
        #read consensus map positions
        dfcon <- read.csv(paste0(lg,'.map'),stringsAsFactor=F,header=T)
        colnames(dfcon) <- c("marker","conpos")
        
        #read individual map positions (with fragments merged together)
        dfmap <- read.csv(paste0(lg,'.merge'),stringsAsFactor=F,header=T)
        colnames(dfmap)[1] <- 'marker'
        mlist <- colnames(dfmap)
        mlist <- mlist[2:length(mlist)]
        
        for (mp in mlist)
        {
            #get the parental origin
            parent <- 'maternal'
            if(substr(mp,1,1) == 'p') parent <- 'paternal'
            
            #population name, eg RGxHA
            name <- substr(mp,2,6)
            
            #select markers with calls in this map
            wh <- !is.na(dfmap[[mp]])
            dfsub <- data.frame(marker=dfmap$marker[wh])
            dfsub$mappos <- dfmap[[mp]][wh]
            
            #normalise map positions
            dfsub$mappos <- dfsub$mappos / max(dfsub$mappos)
            
            #find the consensus map positions of these markers for this lg only
            df2 <- merge(dfcon,dfsub)
            df2$lg <- lg
            df2$mapname <- name
            df2$parent <- parent 
            df2$sg <- y
            df2$hg <- x
            
            df <- rbind(df,df2)
        }
    }
}

gg <- ggplot(df,aes(x=conpos,y=mappos,colour=parent)) +
      ylab("Biparental Position (normalised cM)") +
      xlab("Consensus Position (normalised cM)") +
      theme_bw() +
      geom_point(size=0.2) +
      facet_grid(mapname ~ lg,scales="fixed") + 
      scale_x_continuous(breaks=c(0.0,0.5)) +
      theme(legend.position = "right",
            panel.margin = unit(0, "lines"),
            axis.text=element_text(size=8))

ggsave("comparative_plot_conmap4.png",plot=gg,units="in",dpi=200,width=14.0,height=8.0)

gg <- ggplot(df,aes(x=conpos,colour=mapname,fill=mapname)) +
      ylab("Marker Count") +
      xlab("Consensus Position (normalised cM)") +
      theme_bw() +
      geom_histogram(size=0.2,bins=50) +
      facet_grid(sg ~ hg,scales="fixed") + 
      scale_x_continuous(breaks=c(0.0,0.5)) +
      theme(legend.position = "right",
            panel.margin = unit(0, "lines"),
            axis.text=element_text(size=8))

ggsave("marker_density_conmap4.png",plot=gg,units="in",dpi=200,width=14.0,height=8.0)

mplist <- c("RGxHA","FLxCH","EMxFE","CAxDO","CAxCF")

for ( mp in mplist )
{
    gg <- ggplot(df[df$mapname == mp,],aes(x=conpos,y=mappos,colour=parent)) +
          ylab(paste0(mp," Position (normalised cM)")) +
          xlab("Consensus Position (normalised cM)") +
          theme_bw() +
          geom_point(size=0.2) +
          facet_grid(sg ~ hg,scales="fixed") + 
          scale_x_continuous(breaks=c(0.0,0.5)) +
          theme(legend.position = "right",
                panel.margin = unit(0, "lines"),
                axis.text=element_text(size=8))

    ggsave(paste0(mp,"_vs_conmap4.png"),plot=gg,units="in",dpi=200,width=14.0,height=8.0)
}
