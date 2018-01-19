
df <- data.frame(comb=1:100,map1=NA,map2=NA,map3=NA)
df$map1[1:30]   <- 5*1:30
df$map2[31:60]  <- 2*1:30
df$map3[61:100] <- 3*1:40
model1 <- lm(map1 ~ comb,df)
model2 <- lm(map2 ~ comb,df)
model3 <- lm(map3 ~ comb,df)
plot(df$comb,df$map1)
points(df$comb,df$map2)
points(df$comb,df$map3)

frags <- c("map1","map2","map3")
n_frags = 3

calc_sse <- function(offsets)
{
    df$tmp = NA
    i=1
    offsets[1] = 0.0
    for( mp in frags )
    {
        wh = !is.na(df[[mp]])
        df$tmp[wh] = df[[mp]][wh] + offsets[i]
        i = i + 1
    }

    #plot(df$comb,df$tmp)

    model = lm(tmp ~ comb, df)
    df$pred = predict(model,newdata=df)
    sse = sum((df$tmp-df$pred)^2,na.rm=TRUE)
    return(sse)
}

plot_merged <- function(offsets)
{
    df$tmp = NA
    i=1
    offsets[1] = 0.0
    for( mp in frags )
    {
        wh = !is.na(df[[mp]])
        df$tmp[wh] = df[[mp]][wh] + offsets[i]
        i = i + 1
    }

    plot(df$comb,df$tmp)

}

#calc_sse(c(0,150,200))

res <- optim(c(0,0,0),calc_sse)

plot_merged(res[[1]])
