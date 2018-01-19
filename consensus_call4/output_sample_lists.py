#!/usr/bin/python

#
# output lists of sample filenames, one list per mapping population
#

popn_list =\
[
    'RGxHA',
    'EMxFE',
    'FLxCH',
    'BSxEL',
    'HOxKO',
    'P1xP2',
    'DAxMO',
    'CAxDO',
    'CAxCF',
]

from pymongo import MongoClient
import os

import sys
sys.path.append('/home/vicker/git_repos/rjv')
from rjvmongo import *

cl = MongoClient('192.168.1.100', 27017)
db = cl.strawberry_genetics
db.authenticate('vicker', mongo_password())

for popn in popn_list:
    popn_dir = 'popn_'+popn
    list_file = popn_dir + '/input_file_list'
    
    if not os.path.exists(popn_dir): os.makedirs(popn_dir)
    
    f = open(list_file,'wb')
    f.write('cel_files\n')
    for doc in db.samps.find({'popn':popn}):
        fname = doc['dir'] + '/' + doc['celfile']
        f.write(fname + '\n')
    f.close()

    parent_file = popn_dir + '/matpat_samples'
    
    f = open(parent_file,'wb')
    mat = db.samps.find_one({'popn':popn,'type':'MAT'})
    pat = db.samps.find_one({'popn':popn,'type':'PAT'})
    f.write(mat['celfile'] + ' ' + pat['celfile'] +'\n')
    f.close()
