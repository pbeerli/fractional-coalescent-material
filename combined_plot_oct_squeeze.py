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
def plot(x, yinfluenza, yplasm, yhumpback,  thetasi, thetasp, thetash):
	xheader=r"$\alpha$"
	yheader=r"$\ln mL$"
	#xmin = min(x)
	#ymin = min(y)
	#xmax = max(x)
	#ymax = max(y)

        font = {'family': 'sans-serif',
        'color':  'black',
        'weight': 'bold',
        'size': 16,
        }

        
        fig, axall = plt.subplots(1,3,figsize=(15,3.8))
        fig.subplots_adjust(wspace=0, hspace=0)
        #fig = plt.figure()
        ax1 = axall[0] #fig.add_axes([0.1, 0.1, 0.25, 0.5],
              #             ylim=(-10,0.2), xlim=(0.4,1))
        ax2 = axall[1] #fig.add_axes([0.366, 0.1, 0.25, 0.5],
                       #    ylim=(-10,0.2),yticklabels=[], xlim=(0.4,1))
        ax3 = axall[2] #fig.add_axes([0.632, 0.1, 0.25, 0.5],
                       #    ylim=(-10,0.2), yticklabels=[], xlim=(0.4,1))
        ax1.tick_params(axis='both', which='major', labelsize=15)
        ax2.tick_params(axis='both', which='major', labelsize=15)
        ax3.tick_params(axis='both', which='major', labelsize=15)
        ax1.set_ylim(-10,0.2)
        ax2.set_ylim(-10,0.2)
        ax3.set_ylim(-10,0.2)
        ax1.set_xticks([0.4, 0.5,0.6,0.7,0.8,0.9,1.0])
        ax1.set_xticklabels(['','0.5','','0.7','','0.9',''])
        ax2.set_xticks([0.4, 0.5,0.6,0.7,0.8,0.9,1.0])
        ax2.set_xticklabels(['','0.5','','0.7','','0.9',''])
        ax3.set_xticks([0.4, 0.5,0.6,0.7,0.8,0.9,1.0])
        ax3.set_xticklabels(['','0.5','','0.7','','0.9',''])
        # Twin the x-axis twice to make independent y-axes.
        #ax = axall[0]
        #axes = axall[0]
        axes1 = [ax1, ax1.twinx()] #, ax.twinx()]
        axes2 = [ax2, ax2.twinx()] #, ax.twinx()]
        axes3 = [ax3, ax3.twinx()] #, ax.twinx()]
        axes1[1].set_yticks=[]
        axes2[1].set_yticks=[]
        axes3[1].set_yticks=['x']
        axes1[1].set_yticklabels=[]
        axes2[1].set_yticklabels=[]
        axes3[1].set_yticklabels=[]
        plt.setp(axes1[1].get_yticklabels(), visible=False)
        plt.setp(axes2[0].get_yticklabels(), visible=False)
        plt.setp(axes2[1].get_yticklabels(), visible=False)
        plt.setp(axes3[0].get_yticklabels(), visible=False)
        ax1.text(0.42, -0.5, r'A',fontdict=font)
        ax2.text(0.4, -0.5, r'B',fontdict=font)
        ax3.text(0.4, -0.5, r'C',fontdict=font)
        # Make some space on the right side for the extra y-axis.
        #fig.subplots_adjust(right=0.75)
        # Move the last y-axis spine over to the right by 20% of the width of the axes
        #axes[-1].spines['right'].set_position(('axes', 1.2))

        axes3[1].set_ylabel(r'$\Theta$')
        # To make the border of the right-most axis visible, we need to turn the frame
        # on. This hides the other plots, however, so we need to turn its fill off.
        #axes1[1].set_frame_on(False)
        #axes1[1].patch.set_visible(False)
        #axes2[1].set_frame_on(False)
        #axes2[1].patch.set_visible(False)
        axes3[1].set_frame_on(True)
        axes3[1].patch.set_visible(False)

        # And finally we get to plot things...
        #colors = ('Green', 'Red', 'Blue')
        axes1[0].plot(x,yplasm, marker='o', ms=10,linestyle='-', color='k')
        #axes1[0].plot(x[-3],yplasm[-3], marker='s', ms=8, linestyle=' ', color='k')
        #
        axes2[0].plot(x,yinfluenza, marker='o', ms=10, linestyle='-', color='k')
        #axes2[1].plot(x[-5],yinfluenza[-5], marker='s', ms=8, linestyle=' ', color='k')
        #
        axes3[0].plot(x,yhumpback, marker='o', ms=10, linestyle='-', color='k')
        #axes3[1].plot(x[-2],yhumpback[-2], marker='s', ms=8, linestyle=' ', color='k')
        #
        axes1[0].set_ylabel(r'$\ln$ marginal likelihood', color='k',fontsize=15)
        axes1[0].tick_params(axis=yplasm, colors='k',ms=10)
        #axes1[0].set_ylim(-10.,0.2)
        axes1[0].set_xlim(0.4,1.0)

        axes1[1].plot(x[-11:],thetasp[-11:], marker='o',ms=10,markerfacecolor="white",
                      markeredgecolor='k', linestyle=':', color='k')
        axes1[1].plot(x[-7],thetasp[-7], marker='s' ,markerfacecolor="white",
                      markeredgecolor='k', ms=10, linestyle=' ', color='k')
        #
        axes2[1].plot(x[-11:],thetasi[-11:], marker='o',markerfacecolor="white",
                      markeredgecolor='k', ms=10,linestyle=':', color='k')
        axes2[1].plot(x[-7],thetasi[-7], marker='s' ,markerfacecolor="white",
                      markeredgecolor='k',ms=10, linestyle=' ', color='k')
        #
        axes3[1].plot(x[-11:],thetash[-11:], marker='o',ms=10,markerfacecolor="white",
                      markeredgecolor='k',linestyle=':', color='k')
        axes3[1].plot(x[-2],thetash[-2], marker='s' ,markerfacecolor="white",
                      markeredgecolor='k',ms=10, linestyle=' ', color='k')
        axes3[1].set_xticklabels(['','0.5','','0.7','','0.9',''])
        axes3[1].tick_params(axis='both', which='major', labelsize=15)
        axes3[1].set_ylabel(r'Mutation-scaled population size $\Theta$', color='k',fontsize=15)
        #axes[1].tick_params(axis=yplasm, colors='k')
        axes2[0].set_xlabel(xheader,fontsize=15)
        #axes[0].set_xlim(0.49,1.01)
        axes1[0].set_yticks=[]
        axes2[0].set_yticks=[]
        axes3[0].set_yticks=[]
        axes1[0].set_yticklabels=[]
        axes2[0].set_yticklabels=[]
        axes3[0].set_yticklabels=[]

        plt.tight_layout()        
	#plt.figure()
	#plt.plot(x,y,'b')
	#plt.plot(x,y,'b',marker='o',ms=5)
	#plt.axis([0.5,xmax,-10,ymax+0.2])
	plt.xlabel(xheader)
	#plt.ylabel(yheader)
	plt.savefig('plot2.pdf', format='pdf')


def column(matrix,i, f):
	return [f(row[i]) for row in matrix]


#this is the main program
datap = read_data('plasmodium_dec2017brerun/bf-tabler')
datai = read_data('influenza_rerun_oct27four/bf-tabler')
datah = read_data('humpback_oct272018_rerun/bf-tabler')
xdata = column(datap,0,float)
ydatap = column(datap,1,float)
ydatai = column(datai,1,float)
ydatah = column(datah,1,float)
datap = read_data('plasmodium_dec2017brerun/theta-tabler')
thetasp = column(datap,1,float)
datai = read_data('influenza_rerun_oct27four/theta-tabler')
thetasi = column(datai,1,float)
datah = read_data('humpback_oct272018_rerun/theta-tabler')
thetash = column(datah,1,float)
plot(xdata,ydatap,ydatai,ydatah,thetasp,thetasi,thetash)


