#!/usr/bin/Rscript

#
# interactively select regions where the vesca chromosome disagrees with the consensus map LG
# move the affected marker block to the end of the putative correct vesca chromosome
#

logfile <- "correct_chromosome.log"
outfile <- "corrected_chromosomes.csv"
write("start",file=logfile,append=FALSE)

#load all data: consensus map and vesca v2.0  SNP coords
con <- read.csv("consensus_map_2017-03-20.csv",sep=",",header=F)
colnames(con) <- c("mid","lg","cm")

ves <- read.csv("v2_accurate_snp_posn_no_dups.csv",sep=",",header=F)
colnames(ves) <- c("mid","chr","bp")

#create numerical LG identifiers
con$hg <- as.integer(substr(con$lg,1,1))
con$sg <- sapply(substr(con$lg,2,2),utf8ToInt)-64
con$chrno <- con$sg + (con$hg-1)*4

#create cumulative consensus map position
con$concum <- con$cm * 0.99 + con$chrno

#create cumulative vesca position
#ves$vescum <- ves$bp + (ves$chr-1) * 4.5e7

df <- merge(con,ves)

#find translocations to wrong LG
for ( vchr in c(1,2,3,4,5,6,7) )
{
    #restrict to one vesca chromosome at a time
    df1 <- df[df$chr==vchr,]

    #default view range
    xx <- NULL
    vv <- c()

    while(TRUE)
    {
        plot(df1$bp,df1$concum,xlim=xx,col=df1$hg%%4+1)
        abline(v=vv)

        cat("MARKER TRANSLOCATIONS TO WRONG LG\n")
        cat("z - zoom, m - mark range, x - next chrm\n")
        cmd <- scan(what=character(),n=1)

        if(cmd == "z") {
            #zoom into a selected region
            loc <- locator(n=2)
            if(length(loc) != 2) {
                #reset to default view
                xx <- NULL
            } else {
                xx <- loc$x
            }
        } else if(cmd == "m") {
            #select a block of markers
            loc <- locator(n=2)
            if(length(loc) != 2) {
                #start again
                vv <- c()
            } else {
                #add new range
                vv[[length(vv)+1]] <- min(loc$x)
                vv[[length(vv)+1]] <- max(loc$x)
            }
        } else if(cmd == "x") {
            #done marking the present chromosome
            break
        }
    }

    #move selected range(s) to the most likely vesca chromosome, tacked on the end for now
    while(TRUE)
    {
        if(length(vv) == 0) break

        x1 <- vv[[1]]
        x2 <- vv[[2]]
        vv <- vv[-c(1,2)]

        #grab markers contained in the interval
        chunk <- df1[df1$bp>=x1&df1$bp<=x2,]

        #find most likely true chromosome from consensus map
        truechr <- names(sort(table(chunk$hg),decreasing=T))[1]

        #grab only those mapping to the true chromosome
        chunk <- df1[df1$bp>=x1 & df1$bp<=x2 & df1$hg==truechr,]

        #record to logfile
        write(as.character(chunk$mid),file=logfile,append=T,ncolumns=9999) #marker ids
        write(as.numeric(c(x1,x2,vchr,truechr)),file=logfile,append=T,ncolumns=9999) #startbp,endbp,orig chr,true chr

        #find the start and range of bp positions
        start <- min(chunk$bp)
        end <- max(chunk$bp)
        rng <- end - start

        #find end of current true chr
        curr_max <- max(df$bp[df$chr==truechr])

        #put selected markers at the end of the truechr
        wh <- df$bp>=x1 & df$bp<=x2 & df$chr==vchr
        df$chr[wh] <- truechr
        df$bp[wh] <- df$bp[wh] - start + curr_max + 50000
    }
}

write("end",file=logfile,append=TRUE)

write.csv(df,file=outfile,quote=FALSE)
