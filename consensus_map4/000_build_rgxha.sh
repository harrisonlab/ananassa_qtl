#!/bin/bash

#
# build rgxha map
#

set -eu

#add scripts folder to path
source conf/includes

#source script functions
export SCRIPTDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

source ${SCRIPTDIR}/build_common_funcs.sh
source ${SCRIPTDIR}/build_rgxha_funcs.sh

filter_type_errors
knn_impute
grouping1
remove_redundancy
order_maps1
reinsert_redundant
rename_lgs
optimise_order
pick_best_order
reinsert_redundant2
rename_lgs2
optimise_order2
pick_best_order2
rename_lgs3
get_names
make_rgxha_order.py progeny_names > conf/popn_order
sort_columns

convert_to_joinmap
extract_mapbased_haplotypes

make_uniq_popnorder_locs
make_affycode_popnorder_locs
