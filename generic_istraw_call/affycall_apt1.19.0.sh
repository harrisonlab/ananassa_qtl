#!/bin/bash

#
# call genotypes using APT 1.19.0 and SNPolisher
# this script only deals with a single batch of samples at once
# you must create the list of input cel files first
#
# usage: affycall_apt1.19.0.sh <input_celfile_list> \
#                              istraw35|istraw90 \
#                              <min_dish_qc> <min_call_rate> <call_threshold> \
#                              <output_dir>
#

set -eu

input_list=$1
datatype=$2
min_dqc=$3
min_call_rate=$4
call_threshold=$5
outdir=$6

#scripts directory
export PATH=/home/vicker/git_repos/axiom_strawberry/generic_istraw_call:${PATH}
export PATH=/home/vicker/programs/apt-1.19.0-x86_64-intel-linux/bin:${PATH}

#choose datatype
if [ "${datatype}" == "istraw35" ] ; then
    analysis_dir="/home/groups/harrisonlab/raw_data/raw_celfiles/metadata/Axiom_IStraw35"
    ps2snp="Axiom_IStraw35.r1.ps2snp_map.ps"
    genoqc_file="Axiom_IStraw35.r1.apt-geno-qc.AxiomQC1.xml"
    models="Axiom_IStraw35.r1.generic_prior.txt"  #generic priors
    #models="Axiom_IStraw35.r1.AxiomGT1.models"   #SNP specific priors
    step1_file="Axiom_IStraw35_96orMore_Step1.r1.apt-genotype-axiom.AxiomGT1.apt2.xml"
    step2_file="Axiom_IStraw35_96orMore_Step2.r1.apt-genotype-axiom.AxiomGT1.apt2.xml"

elif [ "${datatype}" == "istraw90" ] ; then
    analysis_dir="/home/groups/harrisonlab/project_files/affymetrix_strawberry/istraw90"
    ps2snp="IStraw90.r1.ps2snp_map.ps"
    genoqc_file="IStraw90.r1.apt-geno-qc.AxiomQC1.xml"
    models="IStraw90.r1.generic_prior.txt"  #generic priors
    #models="IStraw90.r1.AxiomGT1.models"   #SNP specific priors
    step1_file="IStraw90_96orMore_Step1.r1.apt-probeset-genotype.AxiomGT1.xml"
    step2_file="IStraw90_96orMore_Step2.r1.apt-probeset-genotype.AxiomGT1.xml"

else
    echo unknown datatype ${datatype}
    exit 1
fi

mkdir -p ${outdir}
cd ${outdir}
cp ../${input_list} .

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
apt-genotype-axiom \
    --analysis-files-path="${analysis_dir}" \
    --arg-file="${analysis_dir}/${step1_file}" \
    --cel-files="./passed_dqc_cel_files" \
    --out-dir="genotype_step1" \
    --log-file="genotype_step1/step1_qclog.log" \
    --dual-channel-normalization=true

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
apt-genotype-axiom \
    --analysis-files-path="${analysis_dir}" \
    --arg-file="${analysis_dir}/${step2_file}" \
    --cel-files="./passed_callqc_cel_files" \
    --out-dir="genotype_step2" \
    --log-file="genotype_step2/step2_qclog.log" \
    --dual-channel-normalization=true \
    --summaries \
    --write-models
    #--cc-chp-output


#run SNP polisher
polish.R ${ps2snp} ${analysis_dir} ${call_threshold}

#extract PHR and NMH markers, only the best probeset per marker
calls="genotype_step2/recalled.txt"
class="classification/Ps.performance.txt"
supp="supplemental_classification/Ps.performance.uniquevar.txt"

mkdir -p markers
rm -f markers/*

MARKER_CLASS=(PolyHighResolution NoMinorHom CallRateBelowThreshold MonoHighResolution OTV Other)

for mclass in ${MARKER_CLASS[*]}
do
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

#combine MHR, PHR and NMH markers
head -n 1  markers/PolyHighResolution.tsv                      >  markers/mhr_phr_nmh.tsv
tail -n +2 markers/PolyHighResolution.tsv  | sed 's/^AX/PHR/g' >> markers/mhr_phr_nmh.tsv
tail -n +2 markers/NoMinorHom.tsv          | sed 's/^AX/NMH/g' >> markers/mhr_phr_nmh.tsv
tail -n +2 markers/MonoHighResolution.tsv  | sed 's/^AX/MHR/g' >> markers/mhr_phr_nmh.tsv

#transpose marker file to one sample per row rather than per column
transpose_tsv.py markers/mhr_phr_nmh.tsv > markers/mhr_phr_nmh_trans.tsv
