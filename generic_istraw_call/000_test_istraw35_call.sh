#!/bin/bash

#
# test the new single batch script using the a subset of istraw35 samples pipeline
# working directory: /home/vicker/octoploid_mapping/test_istraw35
#

export PATH=/home/vicker/git_repos/rjv:${PATH} #grid_run
export PATH=/home/vicker/git_repos/axiom_strawberry/generic_istraw_call:${PATH}

#create list of celfiles test files
#echo cel_files > input_files
#mysql -h mongo -u vicker -p$(cat /home/vicker/passwords/mysql_mongo_vicker) -D strawberry_samples <<XXX | tail -n +2 | tr '\t' '/' >> input_files
#select path,file from sample where batch like 'istraw35_%' limit 96
#XXX

#call genotypes using APT and SNPolisher
grid_run -Jaffycall -M4 affycall_apt1.19.0.sh input_files istraw35 0.82 96 0.01 affy_calls
