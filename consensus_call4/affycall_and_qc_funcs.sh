#bash functions for 000_affycall_and_qc.sh

#run the affymetrix calling pipeline (APT + SNPolisher)
function run_affy_pipeline
{
    for popndir in popn_*
    do
        #call genotypes
        cd ${popndir}
        call_and_polish.sh input_file_list
        cd ..
    done
}

function classify_samples
{
    for popndir in popn_*
    do
        popn=${popndir/popn_/}
        params=$(cat ${CONFDIR}/similarity | grep ${popn} | cut -d ' ' -f 2)
        
        cd ${popndir}
    
        mat=$(cat matpat_samples | cut -d' ' -f1)
        pat=$(cat matpat_samples | cut -d' ' -f2)
        
        #plot all samples versus parents, classify as obvious non progeny, candidate progeny or parental (inc replicates)
        samples_vs_parents.py\
            markers/mhr_phr_nmh_trans.tsv\
            "${mat}"\
            "${pat}"\
            "${params}"\
            ../figs/${popn}_progeny_vs_parents.png\
            ${popn}\
            > sample_classif &
            
        cd ..
    done
    
    echo waiting for samples_vs_parents.py
    wait
}

function order_populations
{
    #exclude rogues, order samples so that parental samples are at the top
    #process CAxDO to retain only the F1 hybrid and progeny, exclude Camerosa and Dover
    for popndir in popn_*
    do
        popn=${popndir/popn_/}
        
        cd ${popndir}
        
        if [ "${popn}" == "CAxDO" ]
        then
            #process CAxDO differently
            order_populations_caxdo.py\
                sample_classif\
                markers/mhr_phr_nmh_trans.tsv\
                mhr_phr_nmh_trans.tsv\
                parental_samples &
        else
            #process F1 outcross population
            order_populations.py\
                sample_classif\
                markers/mhr_phr_nmh_trans.tsv\
                mhr_phr_nmh_trans.tsv\
                parental_samples &
        fi
        
        cd ..
    done

    echo waiting for order_populations
    wait
}

function merge_duplicate_samples
{
    for popndir in popn_*
    do
        popn=${popndir/popn_/}

        #merge duplicate samples
        duplicate_check_popn.py\
            ${popndir}/mhr_phr_nmh_trans.tsv\
            ${popndir}/parental_samples\
            ${popndir} \
            figs/${popn}_dup_plot.png \
            > ${popndir}/out 2> ${popndir}/err &
    done

    echo waiting for duplicate_check_popn.py
    wait
}

function recode_markers
{
    for popndir in popn_*
    do
        popn=${popndir/popn_/}
        
        #transpose back to one marker per row
        transpose_tsv.py ${popndir}/mhr_phr_nmh_trans_merged.tsv\
            > ${popndir}/mhr_phr_nmh_merged.tsv

        #recode genotype calls wrt parental calls
        recode_markers.py ${popndir}/mhr_phr_nmh_merged.tsv ${popn}\
            > ${popndir}/initial.loc &
    done
    
    echo waiting for recode_markers.py
    wait
}

function qc_filter_markers
{
    for popndir in popn_*
    do
        popn=${popndir/popn_/}
        
        rm -f ${popndir}/qc_stats.csv
        
        #filter out impossible calls
        if [ "${popn}" == "CAxDO" ]
        then
            grep -v -e ' xx' -e '<lmxll>' -e '<nnxnp>' ${popndir}/initial.loc > ${popndir}/noxx.loc
        else
            grep -v ' xx' ${popndir}/initial.loc > ${popndir}/noxx.loc
        fi
        
        echo ${popn} ImpossibleGenotype $(( $(cat ${popndir}/initial.loc | wc --lines) - $(cat ${popndir}/noxx.loc | wc --lines) ))\
            >> ${popndir}/qc_stats.csv

        #missing-data filter by marker
        filter_markers.py 0.0 ${MARKER_PC_MISSING} ${popndir}/noxx.loc > ${popndir}/missfil.loc
        echo ${popn} MissingData $(( $(cat ${popndir}/noxx.loc | wc --lines) - $(cat ${popndir}/missfil.loc | wc --lines) ))\
            >> ${popndir}/qc_stats.csv

        #segregation filter by marker
        filter_markers.py ${MARKER_MIN_PVALUE} 1.0 ${popndir}/missfil.loc > ${popndir}/final.loc
        echo ${popn} SegDistorted $(( $(cat ${popndir}/missfil.loc | wc --lines) - $(cat ${popndir}/final.loc | wc --lines) ))\
            >> ${popndir}/qc_stats.csv
        
        echo ${popn} Maternal $(cat ${popndir}/final.loc | grep -c '<lmxll>' | cat) >> ${popndir}/qc_stats.csv
        echo ${popn} Paternal $(cat ${popndir}/final.loc | grep -c '<nnxnp>' | cat) >> ${popndir}/qc_stats.csv
        echo ${popn} Shared   $(cat ${popndir}/final.loc | grep -c '<hkxhk>' | cat) >> ${popndir}/qc_stats.csv
    done
}

function collect_batch_stats
{
    mkdir -p figs
    rm -f figs/batch_stats.csv
    MARKER_CLASS=(PolyHighResolution NoMinorHom CallRateBelowThreshold MonoHighResolution OTV Other)
    
    for popndir in popn_*
    do
        popn=${popndir/popn_/}

        for mclass in ${MARKER_CLASS[*]}
        do
            echo ${popn} ${mclass} $(tail -n +2 ${popndir}/markers/${mclass}.tsv | wc --lines)\
                >> figs/batch_stats.csv
        done
        
        call_rate_per_row.py ${popndir}/markers/mhr_phr_nmh_trans.tsv ${popn}\
            > ${popndir}/markers/batch_missing_sample.csv &
            
        call_rate_per_row.py ${popndir}/markers/mhr_phr_nmh.tsv ${popn}\
            > ${popndir}/markers/batch_missing_marker.csv &
    done
    
    echo waiting for call_rate_per_row.py
    wait

    cat popn_*/markers/batch_missing_sample.csv > figs/batch_missing_sample.csv
    cat popn_*/markers/batch_missing_marker.csv > figs/batch_missing_marker.csv
    
    list2table.py figs/batch_stats.csv > figs/batch_stats_table.csv
}

function collect_popn_stats
{
    rm -f figs/popn_stats.csv
    
    for popndir in popn_*
    do
        popn=${popndir/popn_/}
        
        #stats output by recode_markers.py
        cat ${popndir}/mhr_phr_nmh_merged.tsv.out\
            | grep -v Informative\
            >> figs/popn_stats.csv
            
        #stats output by affycall_and_qc_funcs.sh:qc_filter_markers
        cat ${popndir}/qc_stats.csv\
            >> figs/popn_stats.csv
    done
    
    list2table.py figs/popn_stats.csv > figs/popn_stats_table.csv
}
