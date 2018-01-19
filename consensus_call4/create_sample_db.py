#!/usr/bin/python

#
# load list of sample names into mongodb
# classify by mapping population membership where appropriate
#

glob_list = \
[
    '/home/vicker/octoploid_mapping/caxcf/celfiles/*.CEL',
    '/home/vicker/octoploid_mapping/xxxxxxxxx/daxmo/symlinks/*.CEL',
    '/home/vicker/octoploid_mapping/xxxxxxxxx/P150761_p1/cels/*.CEL',
    '/home/vicker/octoploid_mapping/crag/*.CEL',
    '/home/groups/harrisonlab/raw_data/raw_celfiles/istraw90_6plates_RGxHA_EMxFE_FLxCH_BSxEL_etc/symlinks/*.CEL',
    '/home/groups/harrisonlab/raw_data/raw_celfiles/istraw90_2plates_P150761Q151138_20160722/symlinks/*.CEL',
    '/home/groups/harrisonlab/raw_data/raw_celfiles/istraw90_4plates_rosbreed/symlinks/known/*.CEL',
]

popn_parents = \
{
    'Redgauntlet':['RGxHA','MAT'],
    'Hapil':['RGxHA','PAT'],
    'Flamenco':['FLxCH','MAT'],
    'Chandler':['FLxCH','PAT'],
    'Emily':['EMxFE','MAT'],
    'Fenella':['EMxFE','PAT'],
    'BSP14':['BSxEL','MAT'],
    'Elvira':['BSxEL','PAT'],
}

popn_codes =\
{
    'RGxHA':None,
    'EMxFE':None,
    'FLxCH':None,
    'BSxEL':None,
}

from glob import glob
from pymongo import MongoClient
from os import path

cl = MongoClient('192.168.1.100', 27017)
db = cl.strawberry_genetics
password = open('/home/vicker/passwords/mongo_vicker').read().strip()
db.authenticate('vicker', password)

for gl in glob_list:
    for fname in glob(gl):
        _dir = path.dirname(fname)
        _base = path.basename(fname)
        tok = _base.split('_')
        batch = tok[0]
        plate = tok[1]
        well = tok[2]
        sample = '_'.join(tok[3:])
        assert sample.endswith('.CEL')
        sample = sample[:-4]

        popn = None

        if batch == 'inra':
            if sample in ['318','319','320']:
                popn = ['CAxCF','MAT']
            elif sample in ['125','218','323','324','322']:
                popn = ['CAxCF','PAT']
            else:
                popn = ['CAxCF','PRO']

        elif batch == 'crag':
            if sample in ['P_Camarosa_2H_2','P_Camarosa_2H_1','1F_Camarosa']:
                popn = ['CAxDO','MAT']
            elif sample in ['1F_E_Dover']:
                popn = ['CAxDO','PAT']
            else:
                popn = ['CAxDO','PRO']

        elif batch == 'xxxx':
            if '44R40_rep' in sample:
                popn = ['P1xP2','MAT']
            elif '56T325_rep' in sample:
                popn = ['P1xP2','PAT']
            else:
                popn = ['P1xP2','PRO']

        elif batch == 'DAxMO':
            if sample == 'Darselect':
                popn = ['DAxMO','MAT']
            elif sample == 'Monterey':
                popn = ['DAxMO','PAT']
            else:
                popn = ['DAxMO','PRO']

        elif sample in popn_parents:        #eg Redgauntlet
            popn = popn_parents[sample]

        elif sample[:5] in popn_codes:    #eg RGxHAxxxxx
            popn = [sample[:5],'PRO']

        elif batch == 'ROS':
            if sample == 'Holiday':
                popn = ['HOxKO','MAT']
            elif sample == 'Korona':
                popn = ['HOxKO','PAT']
            elif sample.startswith('H-'):
                popn = ['HOxKO','PRO']

        doc = {
                "dir":_dir,
                "celfile":_base,
                "batch":batch,
                "plate":plate,
                "well":well,
                "sample":sample,
              }

        if popn != None:
            doc['popn'] = popn[0]
            doc['type'] = popn[1]

        print _base,_dir
        print batch,plate,well,sample,popn
        print

        db.samps.insert_one(doc)

#db.samps
