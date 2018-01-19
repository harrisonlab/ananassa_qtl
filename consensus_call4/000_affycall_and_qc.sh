#!/bin/bash

#
# do initial genotype calling and quality control
#

set -eu

#add scripts folder to path
source conf/includes

#source script functions
export SCRIPTDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

#source ${SCRIPTDIR}/bash_utils.sh
source ${SCRIPTDIR}/affycall_and_qc_funcs.sh

mkdir -p figs

create_sample_db.py
output_sample_lists.py
run_affy_pipeline
classify_samples
order_populations
merge_duplicate_samples
recode_markers
update_sample_db.py
qc_filter_markers
collect_batch_stats
collect_popn_stats

makefig_popnstats.R
makefig_batchstats.R
makefig_batchmissing_sample.R
makefig_batchmissing_marker.R

