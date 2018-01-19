#!/bin/bash

set -eu

inpfile=$1
outfile=$2
SEED=$3

export CL_MAP_RANDOMISE=0
export CL_MAP_CYCLES=50
export CL_GA_ITERS=2000
export CL_GA_OPTIMISEMETH=2
export CL_GA_SKIPORDER1=1

export CL_GA_USEMST=0
export CL_GA_MINLOD=7.0
export CL_GA_MSTNONHK=0
export CL_GA_PROBHOP=0.3333
export CL_GA_MAXHOP=1.0
export CL_GA_PROBMOVE=0.3333
export CL_GA_MAXMOVESEG=1.0
export CL_GA_MAXMOVEDIST=1.0
export CL_GA_PROBINV=0.5
export CL_GA_MAXSEG=1.0
export CL_GIBBS_SAMPLES=500
export CL_GIBBS_BURNIN=20
export CL_GIBBS_PERIOD=1
export CL_GIBBS_PROBSEQUEN=0.0
export CL_GIBBS_PROBUNIDIR=1.0
export CL_GIBBS_MINPROB1=0.1
export CL_GIBBS_MINPROB2=0.0
export CL_GIBBS_TWOPT1=0.1
export CL_GIBBS_TWOPT2=0.0

crosslink_map\
    --inp=${inpfile} --out=${outfile} --map=${outfile/loc/map}\
    --randomise_order=${CL_MAP_RANDOMISE} --seed=${SEED}\
    --ga_gibbs_cycles=${CL_MAP_CYCLES}\
    --ga_iters=${CL_GA_ITERS} --ga_optimise_meth=${CL_GA_OPTIMISEMETH} --ga_skip_order1=${CL_GA_SKIPORDER1}\
    --ga_use_mst=${CL_GA_USEMST} --ga_minlod=${CL_GA_MINLOD} --ga_mst_nonhk=${CL_GA_MSTNONHK}\
    --ga_prob_hop=${CL_GA_PROBHOP} --ga_max_hop=${CL_GA_MAXHOP}\
    --ga_prob_move=${CL_GA_PROBMOVE} --ga_max_mvseg=${CL_GA_MAXMOVESEG} --ga_max_mvdist=${CL_GA_MAXMOVEDIST}\
    --ga_prob_inv=${CL_GA_PROBINV} --ga_max_seg=${CL_GA_MAXSEG}\
    --gibbs_samples=${CL_GIBBS_SAMPLES} --gibbs_burnin=${CL_GIBBS_BURNIN} --gibbs_period=${CL_GIBBS_PERIOD}\
    --gibbs_prob_sequential=${CL_GIBBS_PROBSEQUEN} --gibbs_prob_unidir=${CL_GIBBS_PROBUNIDIR}\
    --gibbs_min_prob_1=${CL_GIBBS_MINPROB1} --gibbs_min_prob_2=${CL_GIBBS_MINPROB2}\
    --gibbs_twopt_1=${CL_GIBBS_TWOPT1} --gibbs_twopt_2=${CL_GIBBS_TWOPT2}

