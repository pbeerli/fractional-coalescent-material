#!/usr/bin/env python
from collections import Counter

import sys
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfFile


#read the data
def read_data():
	"""
	this function reads any data
	"""
	data =[]
	for line in sys.stdin:
		if len(line)<=1:
			break
		newline = line.rstrip('\r\n')
		v = newline.split(' ')
		data.append(v)
	return data





#change data into right format
#plot data
def hist(x):
	xheader="X"
	yheader="Freq."
	plt.figure()
	plt.hist(x)
	plt.xlabel(xheader)
	plt.ylabel(yheader)
	plt.savefig('hist.pdf', format='pdf')


def column(matrix,i):
	return [float(row[i]) for row in matrix]

if sys.argv[1]=='-c':
	cont=True
else:
	cont=False
	min = float(sys.argv[1])
	max = float(sys.argv[2])
	if len(sys.argv)==4:
		freq=True
	else:
		freq=False
#this is the main program
data = read_data()
try:
	l=(len(data[0]))
except:
	l=0

if l==0:
	xdata = data
else:
	xdata = sorted(column(data,0))

grouped = {}
for elem in data:
    key = float(elem[0])
    if key in grouped:
        grouped[key].append(float(elem[1]))
    else:
        grouped[key]=[]
        grouped[key].append(float(elem[1]))            
print grouped

import numpy as np


# seven subplots sharing both x/y axes

f, ax = plt.subplots(7, sharex=True, sharey=True)
for v,i in zip(sorted(grouped.keys()),ax):
        n, bins, patches = i.hist(grouped[v],[0.35,0.45,0.55,0.65,0.75,0.85,0.95,1.05],facecolor='black',)
        print n
        print bins
        i.plot(bins,color='k')
        i.set_yticks([]) 
f.subplots_adjust(hspace=0)
f.text(0.5, 0.04, r'Estimated $\hat{\alpha}$', ha='center')
f.text(0.05, 0.5, r'$\alpha$', ha='center')
f.text(0.1, 0.83, r'0.4', ha='center')
f.text(0.1, 0.715, r'0.5', ha='center')
f.text(0.1, 0.6, r'0.6', ha='center')
f.text(0.1, 0.485, r'0.7', ha='center')
f.text(0.1, 0.371, r'0.8', ha='center')
f.text(0.1, 0.257, r'0.9', ha='center')
f.text(0.1, 0.145, r'1.0', ha='center')

plt.xlim([0.35,1.05])
plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
plt.savefig('hist.pdf', format='pdf')

