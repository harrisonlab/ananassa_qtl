#filter out marker typing errors and other bad markers, rename to standard Affx SNP marker names
function filter_type_errors
{
    export CL_GROUP_MINLOD=7.0
    export CL_GROUP_MATPATLOD=10.0
    export CL_MATPAT_WEIGHTS="01P03"
    export CL_GROUP_REDUNLOD=20.0
    export CL_GROUP_KNN=0
    export CL_GROUP_LOGFILE="fixtypes.log"
    
    cl_fixtypes.sh initial.loc /dev/null /dev/null
    cat fixtypes.log | grep 'type corrected' | cut -d' ' -f4 > type_error_pids
    cat ./initial.loc | fgrep -v -f type_error_pids > filtered.loc
    convert_probe2snp.py filtered.loc ${ISTRAW90_DIR}/${PS2SNPFILE} renamed.loc
    cat ./renamed.loc | fgrep -v -f bad_marker_ids > filtered2.loc
}

function knn_impute
{
    export CL_GROUP_MINLOD=7.0
    export CL_GROUP_KNN=3
    
    cl_knnimpute.sh   filtered2.loc   imputed.loc   /dev/null
}

function grouping1
{
    export CL_GROUP_MINLOD=4.5
    locfile=imputed.loc
    outdir=imputedgrps
    
    rm -rf ${outdir}
    mkdir -p ${outdir}
    cl_group.sh         ${locfile}  ${outdir} ${CL_GROUP_MINLOD}
    cl_phase.sh         ${outdir}  ${outdir}
    cl_mappos.sh        ${outdir}  ${outdir}
    cl_map2csv.sh       ${outdir}  ${outdir}_map.csv
}

function remove_redundancy
{
    export inpdir=imputedgrps
    export outdir=uniqgrps
    export CL_GROUP_MINLOD=7.0
    export CL_GROUP_REDUNLOD=7.0
    
    rm -rf ${outdir}
    mkdir -p ${outdir}
    for locfile in ${inpdir}/*.loc
    do
        redfile=${locfile/\.loc/.redun}
        outfile=${locfile/${inpdir}/${outdir}}
        cl_findredun.sh   ${locfile}   ${redfile}   /dev/null
        cl_extract.sh     ${locfile}   ${redfile}   ${outfile}
    done
}

function order_maps
{
    inpdir=uniqgrps
    outdir=orderedgrps
    
    export CL_MAP_RANDOMISE=1
    export CL_MAP_CYCLES=30
    export CL_GA_ITERS=200000
    export CL_GA_OPTIMISEMETH=0
    export CL_GA_SKIPORDER1=0
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
    export CL_GIBBS_SAMPLES=300
    export CL_GIBBS_BURNIN=20
    export CL_GIBBS_PERIOD=1
    export CL_GIBBS_PROBSEQUEN=0.0
    export CL_GIBBS_PROBUNIDIR=1.0
    export CL_GIBBS_MINPROB1=0.1
    export CL_GIBBS_MINPROB2=0.0
    export CL_GIBBS_TWOPT1=0.1
    export CL_GIBBS_TWOPT2=0.0
    
    rm -rf ${outdir}
    mkdir -p ${outdir}
    for locfile in ${inpdir}/*.loc
    do
        outfile=${locfile/${inpdir}/${outdir}}
        crosslink_map\
            --inp=${locfile} --out=${outfile}\
            --ga_gibbs_cycles=${CL_MAP_CYCLES}\
            --ga_iters=${CL_GA_ITERS} --ga_optimise_meth=${CL_GA_OPTIMISEMETH} --ga_skip_order1=${CL_GA_SKIPORDER1}\
            --ga_use_mst=${CL_GA_USEMST} --ga_minlod=${CL_GA_MINLOD} --ga_mst_nonhk=${CL_GA_MSTNONHK}\
            --ga_prob_hop=${CL_GA_PROBHOP} --ga_max_hop=${CL_GA_MAXHOP}\
            --ga_prob_move=${CL_GA_PROBMOVE} --ga_max_mvseg=${CL_GA_MAXMOVESEG} --ga_max_mvdist=${CL_GA_MAXMOVEDIST}\
            --ga_prob_inv=${CL_GA_PROBINV} --ga_max_seg=${CL_GA_MAXSEG}\
            --gibbs_samples=${CL_GIBBS_SAMPLES} --gibbs_burnin=${CL_GIBBS_BURNIN} --gibbs_period=${CL_GIBBS_PERIOD}\
            --gibbs_prob_sequential=${CL_GIBBS_PROBSEQUEN} --gibbs_prob_unidir=${CL_GIBBS_PROBUNIDIR}\
            --gibbs_min_prob_1=${CL_GIBBS_MINPROB1} --gibbs_min_prob_2=${CL_GIBBS_MINPROB2}\
            --gibbs_twopt_1=${CL_GIBBS_TWOPT1} --gibbs_twopt_2=${CL_GIBBS_TWOPT2}&
    done
    
    sleep 4
    wait

    cl_mappos.sh        ${outdir}  ${outdir}
    cl_map2csv.sh       ${outdir}  ${outdir}_map.csv
}

function reinsert_redundant
{
    export CL_GIBBS_SAMPLES=300
    export CL_GIBBS_BURNIN=20
    export CL_GIBBS_PERIOD=1
    export CL_GIBBS_PROBSEQUEN=0.0
    export CL_GIBBS_PROBUNIDIR=1.0
    export CL_GIBBS_MINPROB1=0.1
    export CL_GIBBS_MINPROB2=0.0
    export CL_GIBBS_TWOPT1=0.1
    export CL_GIBBS_TWOPT2=0.0
    inpdir=orderedgrps
    reddir=imputedgrps
    outdir=redungrps
    
    rm -rf ${outdir}
    mkdir -p ${outdir}
    cat ${reddir}/*.redun > ${outdir}/all_redun
    cl_reinsert_loc.sh   ${inpdir}   imputed.loc   ${outdir}/all_redun   ${outdir}   /dev/null
    cl_reinsert_map.sh   ${inpdir}   ${outdir}/all_redun   ${outdir}
    cl_map2csv.sh        ${outdir}   ${outdir}_map.csv
}

function rename_lgs
{
    inpdir=redungrps
    outdir=renamedgrps
    refmap=/home/vicker/octoploid_mapping/our6plates_plus_RGxHAros/rgxha_map2/snpids.csv
    
    rm -rf ${outdir}
    mkdir -p ${outdir}
    cl_adjustlgs3.py    ${inpdir} ${refmap} ${outdir}
    cl_map2csv.sh       ${outdir}  ${outdir}_map.csv
}
