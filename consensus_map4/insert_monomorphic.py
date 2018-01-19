#!/usr/bin/python

#
# insert monomorphic markers
# process single lg at a time
#

import sys

htypein = sys.argv[1]  #lacking monomorphic eg inpdir/1A_haplotype.csv
hcodein = sys.argv[2]
parfile = sys.argv[3]  #original parental genotype calls
mapfile = sys.argv[4]
htypeout = sys.argv[5] #with monomorphic inserted
hcodeout = sys.argv[6]

lg = htypein.split('/')[-1].split('_')[0]

#load parental genotype calls
parcall = {}
f = open(parfile)
for line in f:
    tok = line.strip().split('\t')
    uid = tok[0]
    mat = tok[1]
    pat = tok[2]
    parcall[uid] = [mat,pat]
f.close()

#load map positions
posn = {}
f = open(mapfile)
for line in f:
    uid,lg,pos = line.strip().split(',')
    posn[uid] = [lg,float(pos)]
f.close()

#load EMxFE markers
marker = {}
f = open(htypein)
header = f.readline().strip().split(',')

nprogeny = (len(header) - 8 ) /2
for line in f:
    tok = line.strip().split(',')
    uid,pos,mtype,phase = tok[0:4]
    
    if not uid in posn: continue
    
    #replace map position
    pos = posn[uid][1]
    
    call = tok[4:]
    
    marker[uid] = [pos,mtype,phase,call]
f.close()

f = open(hcodein)
header2 = f.readline().strip().split(',')
assert header == header2
for line in f:
    tok = line.strip().split(',')
    uid,pos,mtype,phase = tok[0:4]
    
    if not uid in posn: continue

    code = tok[4:]
    
    marker[uid].append(code)
f.close()

#print len(marker)

#insert all other markers from the same lg
for uid in posn:
    #wrong lg
    if posn[uid][0] != lg: continue
    
    #marker already included
    if uid in marker: continue
    
    if not uid in parcall:
        #print 'missing parental calls'
        continue
    
    matpat = parcall[uid]
    
    if '-1' in matpat:
        #print 'missing call'
        continue

    if '1' in matpat:
        #print 'het call'
        continue
        
    if uid not in posn:
        #print 'not in map'
        continue
    
    mat = matpat[0]
    pat = matpat[1]
    
    if mat == '0': mat = 'A'
    else:          mat = 'B'
    
    if pat == '0': pat = 'A'
    else:          pat = 'B'

    #pos,mtype,phase,call,code
    pos = posn[uid][1]
    
    if mat == pat:
        mtype = '<aaxaa>'
    else:
        mtype = '<aaxbb>'
    
    phase = '{--}'
    
    call = [ mat,mat,pat,pat ]
    code = [ mat,mat,pat,pat ]
    
    for i in xrange(nprogeny):
        call.append(mat)
        call.append(pat)
        code.append('-')
        code.append('-')
        
    marker[uid] = [pos,mtype,phase,call,code]
    #print uid

#sort by map order and output to new files

data = [ [uid] + marker[uid] for uid in marker ]
data.sort(key=lambda x:x[1])

for row in data: row[1] = str(int(round(row[1])))

fout = open(htypeout,'wb')
fout.write(','.join(header)+'\n')
for row in data: fout.write(','.join(row[:4]) + ',' + ','.join(row[4]) + '\n')
fout.close()

fout = open(hcodeout,'wb')
fout.write(','.join(header)+'\n')
for row in data: fout.write(','.join(row[:4]) + ',' + ','.join(row[5]) + '\n')
fout.close()
