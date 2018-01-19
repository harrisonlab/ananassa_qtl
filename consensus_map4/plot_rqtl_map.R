#input file format:
#markers,,
#m1,1A,0.0
#m2,1A,0.1... etc

library(qtl)

#no need for any genotype data in the file! only marker name, lg and map position
data <- read.cross("csvr", ".", "rqtl_map.csv")

png(filename="output_map.png",width=1000,height=800)
plotMap(data)
dev.off()
