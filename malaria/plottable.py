#!/usr/bin/env python
# grep : table | tr -s ': _' ' ' | awk '{print $3, $4, $6}' |
# us this instead
# grep : table | tr -s ':_' ' ' | awk '{ print $4,$5 }' | python ../lml_alpha.py
import sys
import matplotlib
from matplotlib.pyplot import cm 
#matplotlib.use('PDF')
import matplotlib.pyplot as plt
#from matplotlib.backends.backend_pdf import PdfFile
import numpy as np
from matplotlib import gridspec

def column(matrix,i):
        return [float(row[i]) for row in matrix]




def read(file):
    with open(file) as f:
        content = f.readlines()
    content = [float(x.strip()) for x in content] 
    return content
    
data = sys.stdin.readlines()
alldata = []
oldj1=0
c=[]
for i in data:
    j = [float(x) for x in i.strip().split()]
    if oldj1 == j[0]:
        c.append([j[1],j[2]])
    else:
        alldata.append(c)
        c = []
        c.append([j[1],j[2]])
        oldj1 = j[0]
alldata.append(c)

alldata.pop(0)
ii = 0.1
for i in  alldata:
    i.sort()
    print ii, i
    ii += 0.1
    
color=iter(cm.rainbow(np.linspace(0,1,20)))
c=next(color)
plt.figure()

ii=0.1
gs1 = gridspec.GridSpec(1,1)
count = 0
for i in alldata:
    x,y = zip(*i)
    ax = plt.subplot(gs1[count])
    count += 1
    ax.scatter(x,y,color=c)
    ax.plot(x,y,color=c, label=ii)
    ii += 0.1
    c=next(color)
#plt.ylim(-50,0)    
#plt.legend()
plt.show()
