#!/bin/bash

#
# combine all individual maps into a consensus map
# this version deals with fragmented LGs by assigning each marker to only one fragment
# and masking out the predictions from the other fragments
# but this seems to cause convergence to suboptimal solutions
#

#run from /home/vicker/octoploid_mapping/consensus_map4

set -eu

export PATH=/home/vicker/git_repos/axiom_strawberry/consensus_map4:${PATH}
export PATH=/home/vicker/rjv_mnt/cluster/git_repos/axiom_strawberry/consensus_map4:${PATH}

inpdirs=bestrenamedgrps
outdir=consensus_map
min_markers=10
min_map_length=4.0
niters=400

mkdir -p ${outdir}

if false ; then ##################################

#combine individual maps into a single matrix
create_map_matrix.py   ${inpdirs}   ${outdir}   ${min_markers}   ${min_map_length}

#create first version of consensus map where all LG fragments produce predictions
#for all missing markers
for csvfile in ${outdir}/*.csv
do
    outfile=${csvfile/csv/out}
    fit_linear_real_weighted.R ${csvfile} FALSE ${niters} popnsizes ${outfile} &
done

sleep 4
wait

#assign missing markers to only one LG fragment so that only the nearest fragment
#is predicting its position
assign_to_fragment.py ${outdir}

#create final version of consensus map where each map only produces one prediction
#for each missing marker
for csvfile in ${outdir}/*.csv
do
    outfile=${csvfile/csv/out2}
    fit_linear_real_weighted.R ${csvfile} TRUE ${niters} popnsizes ${outfile} &
done

sleep 4
wait

fi ###################################################

#combine all consensus lgs into a single csv map file
conmap2csv.py consensus_map.csv ${outdir}/*.csv
