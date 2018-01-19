#!/usr/bin/Rscript

#
# merge the linkage group fragments for all maps where it is fragmented
# merge by finding the offsets between fragments which lead to the minimum
# residual SSE when a linear regression model is fitted to the merged version of the map
# with respect to the current consensus map positions of the markers
#

args = commandArgs(trailingOnly=TRUE)
inpfile = args[1]
mapfile = args[2]
outfile = args[3]

#calculate sum of squared error of a linear regression model fitted to the merged
#LG fragments versus their consensus map positions, assuming the given fragment offset positions
calc_sse <- function(offsets,df2)
{
    #new column to receive the merged positions of the LG
    df2$tmp = NA
    i=1
    #offsets[1] = 0.0
    
    for( mp in fraglist )
    {
        #assuming each marker is only present in at most one of the fragments
        wh = !is.na(df2[[mp]])
        df2$tmp[wh] = df2[[mp]][wh] + offsets[i]
        i = i + 1
    }

    #plot(df$comb,df$tmp)

    model = lm(tmp ~ mean, df2)
    df2$pred = predict(model,newdata=df2)
    sse = sum((df2$tmp-df2$pred)^2,na.rm=TRUE)
    return(sse)
}

#original, individual map positions
df = read.csv(inpfile,row.names = 1)

#current consensus map positions
con = read.csv(mapfile,row.names = 1)

maplist = colnames(df)
nmaps = length(maplist)

out <- subset(con,select=c('mean'))

i <- 1
while(i <= nmaps)
{
    #get basename of map i
    basename <- substr(maplist[i],1,10)
    
    #find all other maps which need to be merged with it
    fraglist = character()
    while(i <= nmaps)
    {
        if(substr(maplist[i],1,10) != basename) break
        fraglist[length(fraglist)+1] <- maplist[[i]]
        i <- i + 1
    }
    
    if(length(fraglist) == 1)
    {
        #no merging required, retain map unmodified
        out[[ fraglist[[1]] ]] = df[[ fraglist[[1]] ]]
        next
    }
    
    df2 = subset(df,select=fraglist)
    df2$mean = con$mean

    #find the minimum sse arrangement of the lg fragments
    #with respect to their consensus map positions
    offsets = numeric(length(fraglist))
    res <- optim(offsets,calc_sse,gr=NULL,df2)
    
    #create merged version of the map
    offsets = res[[1]] - min(res[[1]])
    
    #merged map has suffix of 'M' instead of 1, 2 etc
    name <- paste0(basename,'M')
    
    out[[ name ]] = NA
    
    j=1
    for( mp in fraglist )
    {
        wh = !is.na(df[[mp]])
        out[[name]][wh] = df[[mp]][wh] + offsets[j]
        j = j + 1
    }
}

out <- subset(out, select=-mean )

write.csv(out,outfile,quote=FALSE)
