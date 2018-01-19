#!/usr/bin/Rscript

#
# interactively select marker blocks where the vesca positions are linear with consensus map
# then use correlation coeff to decide whether to invert and mean consensus map position to order
# also allows manual selection of regions to invert in one subgenome only
# to account for genuine differences between vesca and individual ananassa homeologs
#

logfile <- "autofix.log"
inpfile <- "corrected_chromosomes_2017-04-27.csv"
outfile <- "autofixed_vesca.csv"

#catch error where kendall's tau cannot be calculated
#return dummy, nonsignificant values
my_corr_test <- function(x,y)
{
    res <- tryCatch(
        {
            r <- cor.test(x,y,method="kendall")
            p <- as.double(r[[3]])
            t <- as.double(r[[4]])
            c(p,t) #pvalue and tau value
        },
        error = function(x){ return ( c(1,0) ) }, #values to return in case kendall corr cannot be calculated
        warning = function(x){ return ( c(1,0) ) },
        finally = {})

    return ( res )
}


cat("start\n",file=logfile,append=FALSE)

#load dataframe output from correct_chromosome.R
df <- read.csv(inpfile,sep=",",header=T)
df$X <- NULL

#filter out everything where vesca chromosome does not agree with homeolog group
df <- df[df$hg==df$chr,]

#define linear regions, determine if they are inverted
#find centre points, determine best ordering
#for ( cchr in 1:28 )
for ( vchr in 1:7 )
{
    cat("vesca chromosome",vchr,"\n",file=logfile,append=T)

    xx <- NULL #view range
    bb <- c()  #breakpoint locations

    while(TRUE)
    {
        #restrict to one vesca chromosome
        wh <- df$chr==vchr

        plot(df$bp[wh],df$cm[wh],
             xlim=xx,
             col=df$sg[wh],
             xlab=paste0("vesca",vchr),
             ylab=paste0("consensus map hg",vchr))
        abline(v=bb)
        legend("topleft",legend=c(1,2,3,4),col=c("black","red","green","blue"),pch=1)

        cat("z - zoom, b - set breakpoint(s), c - cancel breakpoint, o - invert octoploid, a - autocorrect, x - next chrm\n")
        cmd <- scan(what=character(),n=1)

        if(cmd == "z") {
            #zoom
            loc <- locator(n=2)
            if(length(loc$x) != 2) {
                #reset to default view
                xx <- NULL
            } else {
                xx <- loc$x
            }
        } else if(cmd == "x") {
            #next chr
            break
        } else if(cmd == "b") {
            #set new breakpoint(c)
            loc <- locator(n=99)
            if(length(loc$x) < 1) next #abort setting a new breakpoint

            bb <- c(bb,loc$x)
            #cat(bb,"\n")
        } else if(cmd == "c") {
            #cancel last breakpoint
            bb <- head(bb,-1)
        } else if(cmd == "o") {
            #invert *octoploid* position within the range
            if(length(bb) != 2) next
            cat("which subgenome (1-4)\n")
            sub <- scan(what=integer(),n=1)
            x1 <- min(bb)
            x2 <- max(bb)

            cat(x1,x2,"\n")

            wh <- df$chr==vchr&df$bp>=x1&df$bp<=x2&df$sg==sub

            omin = min(df$cm[wh])
            omax = max(df$cm[wh])
            df$cm[wh] <- omin + omax - df$cm[wh]

        } else if(cmd == "a") {
            #autocorrect
            if(length(bb) < 1) next

            bb <- sort(bb)

            minbp <- min(df$bp[df$chr==vchr])
            maxbp <- max(df$bp[df$chr==vchr])

            #start and end of the interval (vesca bps)
            bp1 <- c(minbp , bb)
            bp2 <- c(bb , maxbp + 1)

            ww <- c() #which markers are in the interval
            mm <- data.frame(frag=integer(),pos=double(),size=double()) #mean consensus map position of the interval
            for ( i in 1:length(bp1) )
            {
                #log info about the interval
                cat("frag",i,bp1[[i]],"-",bp2[[i]],file=logfile,append=T)

                #select markers in this interval
                wh <- df$chr==vchr&df$bp>=bp1[[i]]&df$bp<bp2[[i]]
                ww[[length(ww)+1]] <- wh

                #calculate kendall's tau (correlation) and pvalue
                res <- my_corr_test(df$bp[wh],df$cm[wh])

                #test for significant correlation and negative Tau
                if(res[[1]] <= 0.01 && res[[2]] < 0.0)
                {
                    #invert marker positions within this interval
                    df$bp[wh] <- bp1[[i]] + bp2[[i]] - df$bp[wh]
                    cat(" invert\n",file=logfile,append=T)
                } else {
                    cat(" noinvert\n",file=logfile,append=T)
                }


                #log the marker ids
                #write(as.character(df$mid[wh]),file=logfile,append=T,ncolumns=9999)
                cat(as.character(df$mid[wh]),"\n",file=logfile,append=T)

                #find mean consensus map position within this interval
                mm[nrow(mm)+1,] <- c(i,mean(df$cm[wh]),bp2[[i]]-bp1[[i]])
            }

            #sort frags by mean position
            mm <- mm[order(mm$pos),]
            base <- minbp
            cat("new ordering:",file=logfile,append=T)
            for ( i in 1:nrow(mm) )
            {
                #which fragment
                ff <- mm$frag[i]

                cat(" ",ff,file=logfile,append=T)

                #select markers in this interval
                wh <- ww[[ff]]

                df$bp[wh] <- df$bp[wh] - bp1[[ff]] + base

                base <- base + mm$size[i]
            }
            cat("\n",file=logfile,append=T)

            bb <- c()
        } else {
            cat("unknown command\n")
        }
    }
}

write("end",file=logfile,append=TRUE)

write.csv(df,file=outfile,quote=FALSE)

new_map <- data.frame(mid=df$mid,chr=df$chr,bp=df$bp)
new_map <- new_map[order(new_map$chr,new_map$bp),]

write.table(new_map,file="vesca_fixed.csv",quote=FALSE,row.names=F,col.names=F,sep=',')
