#!/usr/bin/python

#
# add new information to sample db to reflect rogue and duplicate information
# obtained from calling the samples
#

from pymongo import MongoClient
from os import path
from glob import glob

cl = MongoClient('192.168.1.100', 27017)
db = cl.strawberry_genetics
password = open('/home/vicker/passwords/mongo_vicker').read().strip()
db.authenticate('vicker', password)

db.create_collection('variety')

for fname in glob('popn_*'):
    popn = path.basename(fname).split('_')[-1]
    
    f = open(fname + '/sample_classif')
    for line in f:
        tok = line.strip().split()
        celfile = tok[0]
        _type = tok[1].upper()
        
        db.samps.find_one_and_update({'celfile':celfile},{'$set':{'auxtype':_type}})
    f.close()

    f = open(fname + '/merged_samples')
    for line in f:
        tok = line.strip().split()
        merged_name = tok[0].replace('_merged','')
        samples = tok[1:]
        
        db.variety.insert({'representative':merged_name,'samples':samples})
        
        for celfile in samples: db.samps.find_one_and_update({'celfile':celfile},{'$set':{'dup':merged_name}})
    f.close()

cl.close()
