#!/usr/bin/Rscript

#
# following Phil's iterative method
# only fit the model to the real data
# use the predictions to calculate the mean positions
# extended to included different weights per map
# apply different weights to different maps and to extrapolated markers
#

require(splines)

#source("~/rjv_mnt/cluster/git_repos/axiom_strawberry/consensus_map4/fit_linear_real_weighted_test.R")

use_masking = "TRUE"

inpfile = "3B.csv"
plot_map = "mRExRE_3B.1"

weight_file = "../popnsizes"
niters = 500

maskfile = FALSE
if(use_masking == "TRUE")
{
    maskfile = paste(inpfile,".mask",sep="")
    mask = read.csv(maskfile,row.names = 1)
}

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
    
    if(maskfile == FALSE)
    {
        #divide map weight by number of fragments to compensate for making multiple predictions
        all_weights[i] = all_weights[i] / weights[[mp,"count"]]
    }
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
        
        #set masked predictions to NA
        if(maskfile != FALSE)
        {
            wh = mask[[mp]]
            df2[[mp]][wh] = NA
        }
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
        
        #linear only
        model = lm(y ~ x, df)
        df$yp = predict(model,newdata=df)
        
        #add line to plot showing predicted values
        df_tmp = data.frame(x=seq(0,1.0001,0.01))
        df_tmp$y = predict(model, newdata = df_tmp)
        
        #find sse of predicted versus non-missing values
        sse = sse + sum((df$y-df$yp)^2,na.rm=TRUE)

        #plot map data versus mean
        if(mp == plot_map)
        {
            plot(df$x,df$y)
            title(main=mp)
            lines(df_tmp$x,df_tmp$y)
        }

        #copy across unnormalised original map data again
        df2[[mp]] = df[[mp]]
        
        #replace missing values with the predictions
        wh = is.na(df[[mp]])
        df2[[mp]][wh] = df$yp[wh]
        
    }
    
    print(sse)
}

#df_out = subset(df,select=c('mean'))
#write.csv(df_out,outfile,quote=FALSE)
