#!/usr/bin/env python
# use
# grep ^"      1     " real10/outfile* | tr -s '_:' ' '
# in mittag-leffler dir real10
# to get the plot
import sys
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfFile


#read the data
def read_data(filename):
	data =[]
        f = open(filename,'r')
	for line in f:
		if len(line)<=1:
			break
		newline = line.rstrip('\r\n')
		v = newline.split()
		data.append(v)
        f.close()
	return data





#change data into right format
#plot data
def plot(x, y, thetas):
	xheader=r"$\alpha$"
	yheader=r"$\ln mL$"
	xmin = min(x)
	ymin = min(y)
	xmax = max(x)
	ymax = max(y)
        fig, ax = plt.subplots()
        # Twin the x-axis twice to make independent y-axes.
        axes = [ax, ax.twinx()] #, ax.twinx()]
        # Make some space on the right side for the extra y-axis.
        fig.subplots_adjust(right=0.75)
        # Move the last y-axis spine over to the right by 20% of the width of the axes
        #axes[-1].spines['right'].set_position(('axes', 1.2))

        # To make the border of the right-most axis visible, we need to turn the frame
        # on. This hides the other plots, however, so we need to turn its fill off.
        axes[-1].set_frame_on(True)
        axes[-1].patch.set_visible(False)

        # And finally we get to plot things...
        #colors = ('Green', 'Red', 'Blue')
        axes[0].plot(x,y, marker='o', linestyle='-', color='k')
        axes[0].plot(x[-5],y[-5], marker='s', ms=8, linestyle=' ', color='k')
        axes[0].set_ylabel(r'$\ln$ marginal likelihood', color='k')
        axes[0].tick_params(axis=y, colors='k')
        axes[0].set_ylim(-10.,0.2)
        axes[1].plot(x[-11:],thetas[-11:], marker='o',linestyle=':', color='k')
        axes[1].plot(x[-5],thetas[-5], marker='s', ms=8, linestyle=' ', color='k')
        axes[1].set_ylabel(r'Mutation-scaled population size $\Theta$', color='k')
        axes[1].tick_params(axis=y, colors='k')
        axes[0].set_xlabel(xheader)
        axes[0].set_xlim(0.49,1.01)

        
	#plt.figure()
	#plt.plot(x,y,'b')
	#plt.plot(x,y,'b',marker='o',ms=5)
	#plt.axis([0.5,xmax,-10,ymax+0.2])
	#plt.xlabel(xheader)
	#plt.ylabel(yheader)
	plt.savefig('plot.pdf', format='pdf')


def column(matrix,i, f):
	return [f(row[i]) for row in matrix]


#this is the main program
data = read_data('table')
xdata = column(data,0,float)
ydata = column(data,2,float)
data = read_data('thetatable')
thetas = column(data,1,float)
plot(xdata,ydata,thetas)


