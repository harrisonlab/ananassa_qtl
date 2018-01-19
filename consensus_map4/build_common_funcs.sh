function optimise_order
{
    inpdir=uniqgrps

    REPS=20
    MAX_JOBS=40
    NODES=2,4,5,6,7,8,9,10
    MEM=2

    rm -f optimorder_joblist
    rm -rf optim_*_grps

    for TYPE in mstrec ranrec mstglo ranglo ranglo2
    do
        for SEED in $(seq 1 ${REPS})
        do
            outdir=optim_${TYPE}_${SEED}_grps
            rm -rf ${outdir}
            mkdir -p ${outdir}

            for locfile in ${inpdir}/*.loc
            do
                lg=$(basename ${locfile} .loc)
                grid_run -J${TYPE}${SEED}_${lg} -L${MAX_JOBS} -N${NODES} -M${MEM} \
                    "optimise_order.sh ${locfile} ${TYPE} ${SEED} ${outdir} ${lg}" \
                    >> optimorder_joblist
            done
        done
    done
}

function optimise_order2
{
    inpdir=bestrenamedgrps

    REPS=20
    MAX_JOBS=40
    NODES=2,4,5,6,7,8,9,10
    MEM=3

    rm -f optimorder2_joblist
    rm -rf optim2_*_grps

    for TYPE in mstrec ranrec mstglo ranglo ranglo2
    do
        for SEED in $(seq 1 ${REPS})
        do
            outdir=optim2_${TYPE}_${SEED}_grps
            rm -rf ${outdir}
            mkdir -p ${outdir}

            for locfile in ${inpdir}/*.loc
            do
                lg=$(basename ${locfile} .loc)
                grid_run -J${TYPE}_${SEED}_${lg} -L${MAX_JOBS} -N${NODES} -M${MEM} \
                    "optimise_order2.sh ${locfile} ${TYPE} ${SEED} ${outdir} ${lg}" \
                    >> optimorder2_joblist
            done
        done
    done
}

function global_order
{
    inpdir=bestgrps2
    outdir=globalgrps

    rm -rf ${outdir}
    mkdir -p ${outdir}

    for inpfile in ${inpdir}/*.loc
    do
        lg=$(basename ${inpfile} .loc)
        outfile=${outdir}/${lg}.loc
        #grid_run -Jglobal_${lg} -L${MAX_JOBS} -N${NODES} -M${MEM}
        optimise_order_global.sh ${inpfile} ${outfile} 1000 &
        # >> optimorder2_joblist
    done

    sleep 4
    wait

    cl_map2csv.sh       ${outdir}  ${outdir}_map.csv
}

function pick_best_order
{
    outdir=bestgrps

    rm -rf ${outdir}
    mkdir -p ${outdir}

    pick_best_order.py ${outdir} optim_*_grps
    cl_map2csv.sh      ${outdir}  ${outdir}_map.csv
}

function pick_best_order2
{
    outdir=bestgrps2

    rm -rf ${outdir}
    mkdir -p ${outdir}

    pick_best_order.py ${outdir} optim2_*_grps
    cl_map2csv.sh      ${outdir}  ${outdir}_map.csv
}

function reinsert_redundant2
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

    inpdir=bestgrps
    reddir=imputedgrps
    outdir=bestredungrps

    rm -rf ${outdir}
    mkdir -p ${outdir}
    cat ${reddir}/*.redun > ${outdir}/all_redun
    cl_reinsert_loc.sh   ${inpdir}   imputed.loc   ${outdir}/all_redun   ${outdir}   /dev/null
    cl_reinsert_map.sh   ${inpdir}   ${outdir}/all_redun   ${outdir}
    cl_map2csv.sh        ${outdir}   ${outdir}_map.csv
}

function rename_lgs2
{
    inpdir=bestredungrps
    outdir=bestrenamedgrps
    refmap=/home/vicker/octoploid_mapping/our6plates_plus_RGxHAros/rgxha_map2/snpids.csv

    rm -rf ${outdir}
    mkdir -p ${outdir}
    cl_adjustlgs3.py    ${inpdir} ${refmap} ${outdir}
    cl_map2csv.sh       ${outdir}  ${outdir}_map.csv
}

function rename_lgs3
{
    inpdir=bestgrps2
    outdir=bestrenamedgrps2
    refmap=../../../our6plates_plus_RGxHAros/rgxha_map2/snpids.csv

    rm -rf ${outdir}
    mkdir -p ${outdir}
    cl_adjustlgs3.py    ${inpdir} ${refmap} ${outdir}
    cl_map2csv.sh       ${outdir}  ${outdir}_map.csv
}

function sort_columns
{
    inpdir=bestrenamedgrps2
    outdir=popnordergrps
    currorderfile=progeny_names
    neworderfile=conf/popn_order

    rm -rf ${outdir}
    mkdir -p ${outdir}

    for locfile in ${inpdir}/*.loc
    do
        lg=$(basename ${locfile} .loc)
        mapfile=${locfile/loc/map}
        sort_columns.py ${locfile} ${currorderfile} ${neworderfile} ${outdir}/${lg}.loc
        cp ${mapfile} ${outdir}
    done
}

function extract_mapbased_haplotypes
{
    inpdir=popnordergrps
    outdir=haplotypes
    mapfile=bestrenamedgrps2_map.csv

    #rm -rf ${outdir}
    mkdir -p ${outdir}

    for locfile in ${inpdir}/*.loc
    do
        lg=$(basename ${locfile} .loc)
        mapbased_haplotypes2.py \
            ${locfile} \
            ${mapfile} \
            parental_genotypes.tsv \
            conf/popn_order \
            ${outdir}/${lg}
        colour_haplotypes_array2.py \
            ${outdir}/${lg}_haplotypes.csv \
            ${outdir}/${lg}_haplocodes.csv \
            focal_snps \
            ${outdir}/${lg}_coloured.xlsx
    done
}

function extract_mapbased_haplotypes2ped
{
    origdir=origdata_ped                 #original genotype calls in ped format
    locdir=popnordergrps                 #phased loc file from crosslink
    outdir=haplotypes_ped                #output "phased" ped file and a plink map file
    mapfile=bestrenamedgrps2_map.csv     #map order

    #rm -rf ${outdir}
    mkdir -p ${outdir}

    haplotypes2ped.py 

}

function convert_to_joinmap
{
    outdir=joinmap
    inpdir=popnordergrps

    mkdir -p ${outdir}
    cl_map2joinmap.sh   ${inpdir}    ${outdir}/redun.map
    cl_loc2joinmap.sh   ${inpdir}    ${outdir}/redun.loc

    cat ./bestgrps/*.map | grep -v '^group' | cut -f1 | sort > uniq_loci
    cat ./bestgrps2/*.map | grep -v '^group' | cut -f1 | sort > all_loci
    n_redun=$(cat all_loci|wc --lines)
    n_uniq=$(cat uniq_loci|wc --lines)
    cat all_loci | fgrep -v -f uniq_loci > redun_loci

    cat ${outdir}/redun.map | fgrep -v -f redun_loci > ${outdir}/uniq.map
    cat ${outdir}/redun.loc | fgrep -v -f redun_loci \
                            | sed 's/redun/uniq/g' \
                            | sed "s/${n_redun}/${n_uniq}/g" \
                            > ${outdir}/uniq.loc

    echo 'individual names:' >> ${outdir}/redun.loc
    cat conf/popn_order      >> ${outdir}/redun.loc

    echo 'individual names:' >> ${outdir}/uniq.loc
    cat conf/popn_order      >> ${outdir}/uniq.loc
}

function get_names
{
    #affy genotype codes of parents
    get_parental_genotypes.py ../mhr_phr_nmh_merged.tsv > parental_genotypes.tsv

    #progeny names in current ordering
    get_progeny_names.py ../mhr_phr_nmh_merged.tsv \
        | cut -d_ -f4 \
        | sed 's/\.CEL//g' \
        > progeny_names
}

function make_uniq_popnorder_locs
{
    #unique population order loci in csv with progeny names in header
    echo -n 'marker,type,phase,'                                 > uniq_popnorder_genotypes.csv
    cat conf/popn_order | tr '\n' ',' | sed 's/,$/\n/g'         >> uniq_popnorder_genotypes.csv
    cat popnordergrps/*.loc | fgrep -f ./uniq_loci | tr ' ' ',' >> uniq_popnorder_genotypes.csv

    #all loci in csv with progeny names
    echo -n 'marker,type,phase,'                                 > all_popnorder_genotypes.csv
    cat conf/popn_order | tr '\n' ',' | sed 's/,$/\n/g'         >> all_popnorder_genotypes.csv
    cat popnordergrps/*.loc | tr ' ' ','                        >> all_popnorder_genotypes.csv
}

function make_affycode_popnorder_locs
{
    for fname in haplotypes/*_haplotypes.csv
    do
        extract_ordered_affycodes.py ${fname} > ${fname/_haplotypes/_affycodes}
    done

    head -n 1 haplotypes/1A_affycodes.csv > all_affycodes.csv
    tail -q -n +2 haplotypes/*_affycodes.csv >> all_affycodes.csv

    head -n 1 all_affycodes.csv > uniq_affycodes.csv
    cat all_affycodes.csv | fgrep -f ./uniq_loci >> uniq_affycodes.csv
}
