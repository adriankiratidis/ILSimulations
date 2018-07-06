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
data_contactthm = np.loadtxt("../testing-normal-pressure-left-wall.txt")
data_negativederiv = np.loadtxt("../testing-negative_deriv_of_potential.txt")

x_contactthm = data_contactthm[:,0]
y_contactthm = data_contactthm[:,1]

x_negativederiv = data_negativederiv[:,0]
y_negativederiv = data_negativederiv[:,1]

#err = data[:,2]

#A description of the available plotting characters and colours 
#to be used in place of 'rx' can be found here...
#http://matplotlib.org/api/pyplot_api.html#matplotlib.pyplot.plot
#To plot data with errorbars use 
#plt.errorbar(x,y,err,fmt='bo',label="Some Label")
#If you want to plot data but not show a key in the legend use
#label='_nolegend_'
#plt.plot(x_plus,y_plus,'rx',label=r'$n_{+}$')

#y_negativederiv = y_negativederiv - y_negativederiv[-1]
#y_contactthm = y_contactthm - y_contactthm[-1]

y_contactthm = y_contactthm - y_contactthm[-1]

#y_contactthm = y_contactthm * (y_negativederiv[0]/y_contactthm[0])

plt.plot(x_contactthm[2:],y_contactthm[2:],'bo',label='Pressure from contact theorem')
plt.plot(x_negativederiv[2:], y_negativederiv[2:],'gs',label='Pressure from derivative of potential')

#Set the axis labels.  Labelpad option controls the spacing between actual axis and axis label.  The r option tells python to interpret as a raw string literal.
plt.xlabel(r"$\sigma$",labelpad=10)
plt.ylabel(r"${P}$",labelpad=5)

#Set the axis limits if you want to specify a region of interest.
#The default is auto zoom.
#plt.ylim(-0.02,0.06)
#plt.xlim(10,30)

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

savefig("Pressure_comparison_Neutral_Dimer.pdf",bbox_inches='tight')

#Open a window and show the plot
plt.show()
