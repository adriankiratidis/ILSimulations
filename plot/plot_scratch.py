#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
from pylab import *
from matplotlib import rc
import math

#Here we make use of the python plotting library matplotlib.
#Documentation can be found here http://matplotlib.org/

#Set the rc parameters to fix the font and allow tex to be used.
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

#Read in the data from "data.txt".  This will read in the whole file 
#and if necessary arrange the data into a rank 2 array.
#For example, if data.txt contained
#1   4.3   0.2
#2   8.6   0.1
#3   7.4   0.4
#4   8.9   0.4
#5   9.9   0.3
#then data[2,1] = 7.4 that is the 2+1 row and the 1+1 coloumn.
#(Default python array indicies start at 0 unlike fortran which starts at 1.) 
data_28 = np.loadtxt("../bin/test-n_plus_028.txt")
data_20 = np.loadtxt("../bin/test-n_plus_020.txt")
data_10 = np.loadtxt("../bin/test-n_plus_010.txt")
data_05 = np.loadtxt("../bin/test-n_plus_005.txt")

x_28 = data_28[:,0]
y_28 = data_28[:,1]

x_20 = data_20[:,0]
y_20 = data_20[:,1]

x_10 = data_10[:,0]
y_10 = data_10[:,1]

x_05 = data_05[:,0]
y_05 = data_05[:,1]

#err = data[:,2]

y_28[0:5] = -0.1
y_20[0:5] = -0.1
y_10[0:5] = -0.1
y_05[0:5] = -0.1

y_28[95:101] = -0.1
y_20[95:101] = -0.1
y_10[95:101] = -0.1
y_05[95:101] = -0.1

#A description of the available plotting characters and colours 
#to be used in place of 'rx' can be found here...
#http://matplotlib.org/api/pyplot_api.html#matplotlib.pyplot.plot
#To plot data with errorbars use 
#plt.errorbar(x,y,err,fmt='bo',label="Some Label")
#If you want to plot data but not show a key in the legend use
#label='_nolegend_'
plt.plot(x_28,y_28,'rx',label=r'$n_{b}\sigma^{3} = 0.28$')
plt.plot(x_20,y_20,'bo',label=r'$n_{b}\sigma^{3} = 0.20$')
plt.plot(x_10,y_10,'gs',label=r'$n_{b}\sigma^{3} = 0.10$')
plt.plot(x_05,y_05,'yx',label=r'$n_{b}\sigma^{3} = 0.05$')

#Set the axis labels.  Labelpad option controls the spacing between actual axis and axis label.  The r option tells python to interpret as a raw string literal.
plt.xlabel(r"$z/\sigma$",labelpad=10)
plt.ylabel(r"$n\sigma^{3}$",labelpad=5)

#Set the axis limits if you want to specify a region of interest.
#The default is auto zoom.
plt.ylim(0.0,0.4)
#plt.xlim(0.304,0.311)

#Change the position and label of the axis ticks.
#xtickloc = [ -1.5,-1.0,-0.5,0.0,0.5,1.0,1.5 ]
#xtickval = [ '','John','Paul','Adam', 'Bill', 'Dave']
#plt.xticks(xtickloc,xtickval)

#ytickloc = [ -1.5,-1.0,-0.5,0.0,0.5,1.0,1.5 ]
#ytickval = [ '','John','Paul','Adam', 'Bill', 'Dave']
#plt.yticks(ytickloc,ytickval)

#Uncomment if legend required.
plt.legend(loc='best',ncol=1, numpoints=1, frameon=False)

#Uncomment if title is required
#plt.title(r"Some Title")

savefig("testingplot.pdf",bbox_inches='tight')

#Open a window and show the plot
plt.show()
