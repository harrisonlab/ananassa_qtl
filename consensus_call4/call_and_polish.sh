#!/bin/bash

#
# call genotypes from CEL files using apt-tools, SNP polisher and custom python scripts
# get to the point where genotypes are in joinmap compatible format
# and split up into the individual populations
#

set -eu

input_list=$1

analysis_dir=${ISTRAW90_DIR}
genoqc_file="IStraw90.r1.apt-geno-qc.AxiomQC1.xml"
models="IStraw90.r1.generic_prior.txt"  #generic priors
#MODELS="IStraw90.r1.AxiomGT1.models"   #SNP specific priors
min_dqc=${MIN_DQC}
min_call_rate=${MIN_CALL_RATE}

step1_file="IStraw90_96orMore_Step1.r1.apt-probeset-genotype.AxiomGT1.xml"
step2_file="IStraw90_96orMore_Step2.r1.apt-probeset-genotype.AxiomGT1.xml"


#calculate DQC values
#0.82 was the recommended threshold value (see axiom_genotyping_solution_analysis_guide.pdf)
apt-geno-qc \
    --analysis-files-path="${analysis_dir}" \
    --xml-file="${analysis_dir}/${genoqc_file}" \
    --cel-files="${input_list}" \
    --out-file="apt_geno_qc.out" \
    --log-file="apt_geno_qc.log"
    
#create new list of cel files with dish QC >= (eg)0.82
cat apt_geno_qc.out\
    | grep -v cel_files\
    | awk -F '\t' "\$18 >= ${min_dqc} {print \$1}"\
    > passed_dqc
echo cel_files > passed_dqc_cel_files
cat ${input_list}\
    | fgrep -f passed_dqc\
    >> passed_dqc_cel_files

#generate preliminary genotype calls
apt-probeset-genotype \
    --analysis-files-path="${analysis_dir}" \
    --xml-file="${analysis_dir}/${step1_file}" \
    --cel-files="passed_dqc_cel_files" \
    --out-dir="genotype_step1" \
    --log-file="genotype_step1/step1_qclog.log" \
    --read-models-brlmmp="${models}"

#create new list of cel files with call rate >= (eg)96
cat genotype_step1/AxiomGT1.report.txt \
    | grep -v -e '^#' -e cel_files \
    | awk -F '\t' "\$3 >= ${min_call_rate} {print \$1}"\
    > genotype_step1/passed_qc
echo cel_files > passed_callqc_cel_files
cat passed_dqc_cel_files\
    | fgrep -f genotype_step1/passed_qc\
    >> passed_callqc_cel_files

#perform final genotype calling on the filtered set of CEL files
apt-probeset-genotype \
    --analysis-files-path="${analysis_dir}" \
    --xml-file="${analysis_dir}/${step2_file}" \
    --cel-files="passed_callqc_cel_files" \
    --out-dir="genotype_step2" \
    --log-file="genotype_step2/step2_qclog.log" \
    --read-models-brlmmp="${models}"\
    --summaries \
    --write-models \
    --cc-chp-output

#run SNP polisher
polish.R

#extract PHR and NMH markers, only the best probeset per marker
calls="genotype_step2/recalled.txt"
class="classification/Ps.performance.txt"
supp="supplemental_classification/Ps.performance.uniquevar.txt"
    
mkdir -p markers
rm -f markers/*

MARKER_CLASS=(PolyHighResolution NoMinorHom CallRateBelowThreshold MonoHighResolution OTV Other)

for mclass in ${MARKER_CLASS[*]} ; do
    echo ${mclass}
    
    #extract best-per-marker probeset ids
    cat ${class}\
        | grep ${mclass}\
        | awk '$17==1 {print $1}'\
        > markers/tmp${mclass}.ps
        
    #retain only those which also passed supplemental classfication
    cat ${supp}\
        | grep ${mclass}\
        | fgrep -w -f markers/tmp${mclass}.ps\
        | awk '{print $1}'\
        > markers/${mclass}.ps
    rm markers/tmp${mclass}.ps
        
    #extract the genotype calls for these markers
    head -n 1  ${calls}                                    >  markers/${mclass}.tsv
    tail -n +2 ${calls} | fgrep -w -f markers/${mclass}.ps >> markers/${mclass}.tsv
done

#combine PHR and NMH markers, modify marker name to reflect classification
#transpose marker file to one sample per row rather than per column
head -n 1  markers/PolyHighResolution.tsv                     >  markers/phr_and_nmh.tsv
tail -n +2 markers/PolyHighResolution.tsv | sed 's/^AX/PHR/g' >> markers/phr_and_nmh.tsv
tail -n +2 markers/NoMinorHom.tsv         | sed 's/^AX/NMH/g' >> markers/phr_and_nmh.tsv
transpose_tsv.py markers/phr_and_nmh.tsv > markers/phr_and_nmh_trans.tsv

#combine MHR, PHR and NMH markers, do not modify marker name
head -n 1  markers/PolyHighResolution.tsv                      >  markers/mhr_phr_nmh.tsv
tail -n +2 markers/PolyHighResolution.tsv  | sed 's/^AX/PHR/g' >> markers/mhr_phr_nmh.tsv
tail -n +2 markers/NoMinorHom.tsv          | sed 's/^AX/NMH/g' >> markers/mhr_phr_nmh.tsv
tail -n +2 markers/MonoHighResolution.tsv  | sed 's/^AX/MHR/g' >> markers/mhr_phr_nmh.tsv
transpose_tsv.py markers/mhr_phr_nmh.tsv > markers/mhr_phr_nmh_trans.tsv
