#!/usr/bin/Rscript

#
# following Phil's iterative method
# only fit the model to the real data
# use the predictions to calculate the mean positions
# extended to included different weights per map
# apply different weights to different maps
#

#require(splines)

args = commandArgs(trailingOnly=TRUE)
inpfile = args[1]
niters = args[3]
weight_file = args[4]
outfile = args[5]

df = read.csv(inpfile,row.names = 1)

cols = colnames(df)
ncols = length(cols)

maplist = cols
nmaps = length(maplist)

#create copy of map data
df2 = subset(df,select=maplist)

#create map weights
weights = read.csv(weight_file,header=F,row.names=1,col.names=c('map','weight'))
weights$count = 0
all_weights = numeric(length=nmaps)
for (i in 1:nmaps)
{
    mp = substr(maplist[i],1,6)
    weights[[mp,"count"]] = weights[[mp,"count"]] + 1
}

for (i in 1:nmaps)
{
    mp = substr(maplist[i],1,6)
    all_weights[i] = weights[[mp,"weight"]]
    
    #divide map weight by number of fragments to compensate for making multiple predictions
    all_weights[i] = all_weights[i] / weights[[mp,"count"]]
}

for (it in 1:niters)
{
    #==calc mean normalised marker positions across all maps==
    #normalise map data
    for (mp in maplist)
    {
        low = min(df2[[mp]],na.rm=T)
        high = max(df2[[mp]],na.rm=T)
        df2[[mp]] = (df2[[mp]] - low) / (high - low)
    }

    #take mean of normalised positions, ignoring missing values
    #df$mean = rowMeans(df2,na.rm=T)
    df$mean = apply(df2, 1, function(row) weighted.mean(row, all_weights,na.rm=TRUE))

    #==fit lines / splines and predict map positions from mean positions for all markers
    sse = 0

    #for each map
    for (mp in maplist)
    {
        df$x = df$mean
        df$y = df[[mp]]
        
        #fit model
        model = lm(y ~ x, df)
        
        #predict values
        df$yp = predict(model,newdata=df)
        
        #find sse of predicted versus non-missing values
        sse = sse + sum((df$y-df$yp)^2,na.rm=TRUE)

        #copy across unnormalised original map data again
        df2[[mp]] = df[[mp]]
        
        #replace only missing values with the predictions
        wh = is.na(df[[mp]])
        df2[[mp]][wh] = df$yp[wh]
    }
}

df_out = subset(df,select=c('mean'))
write.csv(df_out,outfile,quote=FALSE)
