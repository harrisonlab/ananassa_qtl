#!/bin/bash

#
# filter haplocode csv files to retian only popseq progeny
# run from ~/octoploid_mapping/consensus_map4/popn_RGxHA/map

LIST=popseq_progeny_filter_list

for x in haplotypes/??_haplocodes_cleaned.csv
do
    echo ${x}
    #~/git_repos/rjvbio/transpose_csv.py ${x} | fgrep -f ${LIST} > ${x}_trans
    ~/git_repos/rjvbio/transpose_csv.py ${x}_trans > ${x/cleaned/popseq}
done
