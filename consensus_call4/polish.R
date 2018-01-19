#!/usr/bin/Rscript

#import SNPolisher library, do not comment this line out
library("SNPolisher")
library("methods") #otherwise Rscript complains that function "is" is not found

ps2snp="IStraw90.r1.ps2snp_map.ps"
call_threshold=as.numeric(Sys.getenv("POLISH_CALL_THRESHOLD"))
analysis_dir=Sys.getenv("ISTRAW90_DIR")
varz=3
clustermin=5

#=============re-call genotypes using a different threshold, may convert some to NOCALL
Ps_CallAdjust(callFile="genotype_step2/AxiomGT1.calls.txt",
              confidenceFile="genotype_step2/AxiomGT1.confidences.txt",
              threshold=call_threshold,
              outputFile="genotype_step2/recalled.txt")

#========calculate SNP quality metrics
dir.create("metrics",showWarnings = FALSE)
Ps_Metrics(posteriorFile="genotype_step2/AxiomGT1.snp-posteriors.txt",
           callFile="genotype_step2/recalled.txt",
           output.metricsFile="metrics/metrics.txt")

#========classify SNPs
Ps_Classification(metricsFile="metrics/metrics.txt",
                  ps2snpFile=paste0(analysis_dir,"/",ps2snp),
                  SpeciesType="Polyploid",
                  output.dir="classification",
                  output.converted=TRUE)

#========perform additional classification to try to catch plots with more than 3 clusters
Ps_Classification_Supplemental(performanceFile="classification/Ps.performance.txt",
                               summaryFile="genotype_step2/AxiomGT1.summary.txt",
                               callFile="genotype_step2/recalled.txt",
                               posteriorFile="genotype_step2/AxiomGT1.snp-posteriors.txt",
                               output.dir="supplemental_classification",
                               AA.varX.Z.cut=varz, AA.varY.Z.cut=varz,
                               AB.varX.Z.cut=varz, AB.varY.Z.cut=varz,
                               BB.varX.Z.cut=varz, BB.varY.Z.cut=varz,
                               clustermin=clustermin)
