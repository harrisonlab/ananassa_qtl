#!/bin/bash

#
# combine all individual maps into a consensus map
# this version deals with fragmented LGs by explicitly merging
# them based on their initial consensus map positions
#

#run from /home/vicker/octoploid_mapping/consensus_map4

set -eu

export PATH=/home/vicker/git_repos/axiom_strawberry/consensus_map4:${PATH}
export PATH=/home/vicker/rjv_mnt/cluster/git_repos/axiom_strawberry/consensus_map4:${PATH}

inpdirs=bestrenamedgrps2
outdir=consensus_map
min_markers=10
min_map_length=4.0
niters=300

nmerges=10

NODES="-N2,5,6,7,8,9"
MAXJOBS="-L40"
MEM="-M2"

mkdir -p ${outdir}


#combine individual maps into a single matrix
create_map_matrix.py   ${inpdirs}   ${outdir}   ${min_markers}   ${min_map_length}



#create first version of consensus map where all LG fragments produce predictions
#for all missing markers
rm -f joblist
for csvfile in ${outdir}/*.csv
do
    mapfile=${csvfile/csv/map}
    JOBNAME="${csvfile//\//_}${RANDOM}"
    grid_run -J${JOBNAME} ${NODES} ${MAXJOBS} ${MEM} \
        "fit_linear_real_weighted.R ${csvfile} FALSE ${niters} popnsizes ${mapfile}" \
        >> joblist
done

#sleep 4
#wait
grid_wait -Ljoblist


for i in $(seq 1 ${nmerges})
do
    rm -f joblist
    for csvfile in ${outdir}/*.csv
    do
        echo ${csvfile}
        mapfile=${csvfile/csv/map}
        mergefile=${csvfile/csv/merge}

        #merge all fragmented LGs based on their current consensus map positions
        JOBNAME="merge${RANDOM}"
        grid_run -J${JOBNAME} ${NODES} ${MAXJOBS} ${MEM} \
            "consensus_merge_frags.R ${csvfile} ${mapfile} ${mergefile}" \
            >> joblist

        #create refined version of consensus map
        JOBNAME2="fit${RANDOM}"
        grid_run -J${JOBNAME2} -W${JOBNAME} ${NODES} ${MAXJOBS} ${MEM} \
            "fit_linear_real_weighted.R ${mergefile} FALSE ${niters} popnsizes ${mapfile}" \
            >> joblist
    done

    grid_wait -Ljoblist
    #sleep 4
    #wait
done

#combine all consensus lgs into a single csv map file
conmap2csv.py consensus_map4.csv ${outdir}/*.map
