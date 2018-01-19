#functions for 020_build_maps.sh

function fix_type_errors_inner
{
    cl_fixtypes.sh initial.loc fixtype.loc /dev/null
    convert_probe2snp.py fixtype.loc ${ISTRAW90_DIR}/${PS2SNPFILE} renamed.loc
}

#fix marker typing errors, then rename to standard Affx SNP marker names
function fix_type_errors
{
    source ${CONFDIR}/defaults
    export CL_GROUP_LOGFILE="fixtypes.log"
    
    for popn in $*
    do
        popndir=popn_${popn}
        mkdir -p ${popndir}/map
        cp ${popndir}/final.loc ${popndir}/map/initial.loc
        cd ${popndir}/map
        fix_type_errors_inner &
        cd ../..
    done
    
    echo waiting for fix_type_errors_inner
    sleep 4
    wait
}

function knn_impute
{
    source ${CONFDIR}/defaults
    
    for popn in $*
    do
        popndir=popn_${popn}
        cd ${popndir}/map
        cl_knnimpute.sh   renamed.loc   imputed.loc   /dev/null &
        cd ../..
    done
    
    echo waiting for cl_knnimpute.sh
    sleep 4
    wait
}

function grouping_inner
{
    popn=$1
    
    cl_group.sh         ${locfile}  ${outdir}    ${CL_GROUP_MINLOD}
    cl_phase.sh         ${outdir}  ${outdir}
    cl_mappos.sh        ${outdir}  ${outdir}
    cl_map2csv.sh       ${outdir}  ${outdir}_map.csv
}

function grouping
{
    source ${CONFDIR}/defaults

    export locfile=$1
    export outdir=$2
    shift 2

    rm -rf ${outdir}
    mkdir -p ${outdir} 

    for popn in $*
    do
        export CL_GROUP_MINLOD=$(grep -e "${popn}" ${CONFDIR}/group_minlod | cut -d' ' -f2)
        popndir=popn_${popn}
        cd ${popndir}/map
        grouping_inner ${popn} &
        cd ../..
    done
    
    echo waiting for grouping_inner
    sleep 4
    wait
}

#work out which provisional LGs are fusions of one or more true groups
#by comparison to a reference map
function find_fused_groups
{
    inpdir=$1
    shift
    
    export refmapfile=$(readlink -f ${CONFDIR}/refmap.csv) #get full path to reference map file
    
    for popn in $*
    do
        export MINCOUNT=$(grep -e "${popn}" ${CONFDIR}/fusedgroup | cut -d' ' -f2)
        export MINPROP=$(grep -e "${popn}" ${CONFDIR}/fusedgroup | cut -d' ' -f3)
        popndir=popn_${popn}
        cd ${popndir}/map/${inpdir}
        find_fused_groups.py ${refmapfile} ${MINCOUNT} ${MINPROP} *.loc > fused_groups &
        cd ../../..
    done
    
    echo waiting for find_fused_groups.py
    sleep 4
    wait
}

function filter_bad_markers_inner
{
    rm -f lod_stats
    for locfile in *.loc
    do
        rflodfile=${locfile/\.loc/.rflod}
        crosslink_rflod --inp=${locfile} --out=${rflodfile} --min_lod=${RFLOD_BAD}
        compare_grouping.py ${rflodfile} ${refmapfile} >> lod_stats
    done
    
    find_bad_markers.py lod_stats fused_groups ${popn} ${K_RATIO} ${K_CONST} > bad_markers
    fgrep -e ${popn} ${extrabadfile} | cut -d' ' -f 2 >> bad_markers
    cut -d' ' -f1 bad_markers > bad_marker_ids
    fgrep -f bad_marker_ids -v ../imputed.loc > ../filtered.loc
}

#find which markers appear to be falsely holding together homeologous lgs
function filter_bad_markers
{
    export refmapfile=$(readlink -f ${CONFDIR}/refmap.csv) #get full path to reference map file
    export extrabadfile=$(readlink -f ${CONFDIR}/extra_bad_ids) #get full path to extra bad ids file
    
    for popn in $*
    do
        export popn
        export RFLOD_BAD=$(grep -e "${popn}" ${CONFDIR}/rflod_bad | cut -d' ' -f2)
        export K_RATIO=$(grep -e "${popn}" ${CONFDIR}/badmarker_criteria | cut -d' ' -f2)
        export K_CONST=$(grep -e "${popn}" ${CONFDIR}/badmarker_criteria | cut -d' ' -f3)
        popndir=popn_${popn}
        cd ${popndir}/map/imputedgrps
        filter_bad_markers_inner ${popn} &
        cd ../../..
    done
    
    echo waiting for filter_bad_markers_inner
    sleep 4
    wait
}

function merge_split_groups_inner1
{
    rm -rf tmpmergegrps
    mkdir -p tmpmergegrps

    cat ${inpdir}/split_groups | \
    while read line
    do
        lg=$(echo ${line} | cut -d' ' -f1)
        frags=$(echo ${line} | cut -d' ' -f2 | sed 's/,/ /g')

        #merge all frags from the same lg together
        for frag in ${frags}
        do
            #prefix marker names with fragment name
            cat ${inpdir}/${frag}.loc | sed "s/^/${frag}-/g" >> tmpmergegrps/${lg}.loc
        done
    done
}

function merge_split_groups_inner2
{
    cl_phase.sh tmpmergegrps tmpmergegrps

    cat ${inpdir}/split_groups | \
    while read line
    do
        lg=$(echo ${line} | cut -d' ' -f1)
        frags=$(echo ${line} | cut -d' ' -f2 | sed 's/,/ /g')

        #calculate 2-pt rf and lod
        crosslink_rflod --inp=tmpmergegrps/${lg}.loc\
                        --out=tmpmergegrps/${lg}.rflod\
                        --min_lod=${RFLOD_MIN}

        #create merged LGs
        merge_split_groups.py tmpmergegrps/${lg}.rflod ${inpdir} ${outdir} ${lg} "${frags}"
    done

    cl_phase.sh    ${outdir}  ${outdir}
    cl_mappos.sh   ${outdir}  ${outdir}
    cl_map2csv.sh  ${outdir}  ${outdir}_map.csv
}

#work out which true LGs are split into two or more pieces wrt a reference map
function merge_split_groups
{
    export inpdir=$1
    export outdir=$2
    shift 2
    
    export refmapfile=$(readlink -f ${CONFDIR}/refmap.csv) #get full path to reference map file

    #find which fragments belong to which true LG
    for popn in $*
    do
        popndir=popn_${popn}
        cd ${popndir}/map/${inpdir}
        find_split_groups.py ${refmapfile} *.loc > split_groups &
        cd ../../..
    done
    
    echo waiting for find_split_groups.py
    sleep 4
    wait

    #merge fragments which belong to the same true LG and have significant LOD
    for popn in $*
    do
        export RFLOD_MIN=$(grep -e "${popn}" ${CONFDIR}/mergegroups | cut -d' ' -f2)
        popndir=popn_${popn}
        cd ${popndir}/map
        rm -rf ${outdir}
        mkdir -p ${outdir}
        merge_split_groups_inner1
        merge_split_groups_inner2 &
        cd ../..
    done
    
    echo waiting for merge_split_groups_inner2
    sleep 4
    wait
}

function refine_order1_inner
{
    export CL_PARALLEL_JOBS=24
    export CL_IGNORE_PREVIOUS=1
    export CL_MAP_RANDOMISE=0
    export CL_GA_SKIPORDER1=0
    export CL_MAP_CYCLES=500
    export CL_GA_USEMST=30
    export CL_GA_ITERS=50000
    cl_refine_order.sh   ${inpdir}   ${outdir}   20   24   /dev/null #refine map ordering using trial and error
}

#refine order using recombination event counting
function refine_order1
{
    export inpdir=$1
    export outdir=$2
    shift 2
    
    source ${CONFDIR}/defaults
    
    #find which fragments belong to which true LG
    for popn in $*
    do
        popndir=popn_${popn}
        cd ${popndir}/map
        refine_order1_inner &
        cd ../..
    done
    
    echo waiting for refine_order1_inner
    sleep 4
    wait

    for popn in $*
    do
        popndir=popn_${popn}
        cd ${popndir}/map
        cl_mappos.sh   ${outdir}  ${outdir}
        cl_map2csv.sh  ${outdir}  ${outdir}_map.csv
        cd ../..
    done
}

function refine_order2_inner
{
    export CL_PARALLEL_JOBS=30
    export CL_GA_OPTIMISEMETH=2
    export CL_IGNORE_PREVIOUS=1
    export CL_MAP_RANDOMISE=0
    export CL_GA_SKIPORDER1=1
    export CL_MAP_CYCLES=50
    export CL_GA_USEMST=30
    export CL_GA_ITERS=10000
    cl_refine_order.sh   ${inpdir}   ${outdir}   5   30   /dev/null
}

#refine order using "global" optimisation method
function refine_order2
{
    export inpdir=$1
    export outdir=$2
    shift 2
    
    source ${CONFDIR}/defaults
    
    #find which fragments belong to which true LG
    for popn in $*
    do
        popndir=popn_${popn}
        cd ${popndir}/map
        refine_order2_inner &
        cd ../..
    done
    
    echo waiting for refine_order2_inner
    sleep 4
    wait

    for popn in $*
    do
        popndir=popn_${popn}
        cd ${popndir}/map
        cl_mappos.sh   ${outdir}  ${outdir}
        cl_map2csv.sh  ${outdir}  ${outdir}_map.csv
        cd ../..
    done
}

function remove_redundancy_inner
{
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

function remove_redundancy
{
    export inpdir=$1
    export outdir=$2
    shift 2
    
    source ${CONFDIR}/defaults
    
    #find which fragments belong to which true LG
    for popn in $*
    do
        export CL_GROUP_MINLOD=$(grep -e "${popn}" ${CONFDIR}/group_minlod | cut -d' ' -f2)
        export CL_GROUP_REDUNLOD=$(grep -e "${popn}" ${CONFDIR}/group_minlod | cut -d' ' -f2)
        popndir=popn_${popn}
        cd ${popndir}/map
        remove_redundancy_inner &
        cd ../..
    done
    
    echo waiting for remove_redundancy_inner
    sleep 4
    wait
}
