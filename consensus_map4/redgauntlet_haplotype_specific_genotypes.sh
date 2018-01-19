#!/bin/bash

#
# extract redgauntlet genotype specific calls
#

#run from /home/vicker/octoploid_mapping/consensus_map4/popn_RGxHA/redgauntlet_haplotypes

#==========bash
#affy call codes, with missing data
cat ../mhr_phr_nmh_merged.tsv | awk '{print $1,$1,$2}' | tail -n +2 > rg_calls.tsv

cat ../mhr_phr_nmh_merged.tsv | awk '{print $1,$3}' | tail -n +2 > ha_calls.tsv
cat ./ha_calls.tsv | probe_to_snp.py > ha_conv.tsv

#phased locfiles
cat ../map/bestrenamedgrps2/*.loc > bestrenamedgrp2_all.loc


#=============R
options(width=120)
df <- read.csv("bestrenamedgrp2_all.loc",sep="",header=F)
df <-data.frame(snpid=df$V1,mtype=df$V2,phase=df$V3)

df2 <- read.csv("rg_snpid_clusters_probeid_rgcall.tsv",sep="",header=F)
colnames(df2) <- c("snpid","clustering","probeid","rgcall")

all <- merge(df,df2,all=T)

df3 <- read.csv("ha_conv.tsv",sep="",header=F)
colnames(df3) <- c("snpid","hacall")

all <- merge(all,df3,all=T)

all$rghap1 = NA
all$rghap2 = NA

#rows where RG is het and HA is AA and matphase is 0
wh = !is.na(all$mtype)&all$mtype=="<lmxll>"&all$hacall==0&all$phase=="{0-}"
all$rghap1[wh]="A"
all$rghap2[wh]="B"

#rows where RG is het and HA is BB and matphase is 0
wh = !is.na(all$mtype)&all$mtype=="<lmxll>"&all$hacall==2&all$phase=="{0-}"
all$rghap1[wh]="B"
all$rghap2[wh]="A"

#rows where RG is het and HA is AA and matphase is 1
wh = !is.na(all$mtype)&all$mtype=="<lmxll>"&all$hacall==0&all$phase=="{1-}"
all$rghap1[wh]="B"
all$rghap2[wh]="A"

#rows where RG is het and HA is BB and matphase is 1
wh = !is.na(all$mtype)&all$mtype=="<lmxll>"&all$hacall==2&all$phase=="{1-}"
all$rghap1[wh]="A"
all$rghap2[wh]="B"

#nnxnp and RG is 0
wh = !is.na(all$mtype) & all$mtype=="<nnxnp>" & all$rgcall==0
all$rghap1[wh]="A"
all$rghap2[wh]="A"

#nnxnp and RG is 2
wh = !is.na(all$mtype) & all$mtype=="<nnxnp>" & all$rgcall==2
all$rghap1[wh]="B"
all$rghap2[wh]="B"

#mhr and RG is 0
wh = all$clustering=="MHR" & all$rgcall==0
all$rghap1[wh]="A"
all$rghap2[wh]="A"

#mhr and RG is 2
wh = all$clustering=="MHR" & all$rgcall==2
all$rghap1[wh]="B"
all$rghap2[wh]="B"

#hkxhk and matphase is "0"
wh = !is.na(all$mtype) & all$mtype=="<hkxhk>" & (all$phase=="{00}" | all$phase=="{01}")
all$rghap1[wh]="A"
all$rghap2[wh]="B"

#hkxhk and matphase is "1"
wh = !is.na(all$mtype) & all$mtype=="<hkxhk>" & (all$phase=="{10}" | all$phase=="{11}")
all$rghap1[wh]="B"
all$rghap2[wh]="A"

#get consensus map positions
conmap <- read.csv("../../consensus_map3.csv",sep=",",header=F)
colnames(conmap) <- c("snpid","lg_con","pos_con")
all <- merge(all,conmap,all=T)

#get rgxha map positions
rgxhamap <- read.csv("../map/bestrenamedgrps2_map.csv",sep=",",header=F)
colnames(rgxhamap) <- c("snpid","lg_rgxha","pos_rgxha")
all <- merge(all,rgxhamap,all=T)

#now contains duplicated entries for those markers which map to more than one location
#on consensus map
write.csv(all,file="maindataframe.csv",quote=F,row.names=F)
