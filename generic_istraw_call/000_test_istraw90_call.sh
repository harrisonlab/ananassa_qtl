#!/bin/bash

#
# test the new single batch script using the old istraw90 pipeline
# working directory: /home/vicker/octoploid_mapping/test_istraw90
#

export PATH=/home/vicker/git_repos/rjv:${PATH} #grid_run
export PATH=/home/vicker/git_repos/axiom_strawberry/generic_istraw_call:${PATH}

#create list of celfiles test files
echo cel_files > input_files
mysql -h mongo -u vicker -p$(cat /home/vicker/passwords/mysql_mongo_vicker) -D strawberry_samples <<XXX | tail -n +2 | tr '\t' '/' >> input_files
select path,file from sample where batch like 'istraw90_batch1_plate1'
XXX

#call genotypes using APT and SNPolisher
grid_run -Jaffycall -M4 affycall_apt1.16.1.sh input_files istraw90 0.82 96 0.01 affy_calls
