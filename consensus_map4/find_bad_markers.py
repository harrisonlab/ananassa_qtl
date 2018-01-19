#!/usr/bin/python

'''
plot lod information
'''

import sys
import matplotlib.pyplot as plt

#figure size, inches
width = 5.0
height = 5.0
dpi = 600

#plot config
linewidth = 0.2
alpha = 0.5

inpname = sys.argv[1]
fusedgrps = sys.argv[2]
popn = sys.argv[3]
k_ratio = float(sys.argv[4])
k_const = float(sys.argv[5])
outname = 'lod_stats.png'

#load info about which groups are fused in the current grouping wrt the reference map
fused = {}
f = open(fusedgrps) #newlg reflg reflg [reflg...]
for line in f:
    tok = line.strip().split()[1:]
    for i,x in enumerate(tok):
        for y in tok[i+1:]:
            if x < y:#sort alphabetically
                fused[x+':'+y] = None
            else:
                fused[y+':'+x] = None
f.close()

f = open(inpname)
xx = []
yy = []
cc = []
fig = plt.figure(figsize=(width, height)) 

for line in f:
    uid,lg1,sum1,n1,lg2,sum2,n2 = line.strip().split()
    sum1 = float(sum1)
    sum2 = float(sum2)
    n1 = int(n1)
    n2 = int(n2)

    if lg1 > lg2: lg1,lg2 = lg2,lg1 #sort alphabetically
    lgpair = lg1+':'+lg2
    
    if sum1 < sum2: sum1,sum2 = sum2,sum1
    xx.append(sum1)
    yy.append(sum2)
    
    cc.append('#0000ff') #not classified as 'bad'
    
    if sum2 > k_const and sum2/sum1 > k_ratio:
        plt.text(sum1,sum2,'%s-%s'%(lg1,lg2))
        if lgpair in fused:
            cc[-1]= '#ff0000' #classified as 'bad' using current thresholds
            print uid,lg1,lg2

f.close()

minx = min(xx)
maxx = max(xx)
miny = min(yy)
maxy = max(yy)

plt.scatter(xx,yy,c=cc,alpha=alpha,linewidth=linewidth)
plt.plot([minx, maxx], [k_const, k_const], color='k', linestyle='-', linewidth=2)
plt.plot([minx, maxx], [minx*k_ratio, maxx*k_ratio], color='k', linestyle='-', linewidth=2)
plt.xlabel('sum of LODs to first LG')
plt.ylabel('sum of LODs to second LG')
plt.title(popn)
    
plt.savefig(outname, dpi=dpi) #, bbox_inches='tight'
