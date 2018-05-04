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
#data_plus = np.loadtxt("../bin/test-n_plus.txt")
data_neutral = np.loadtxt("../bin/test-n_neutral_separation10.txt")
data_neutral20 = np.loadtxt("../bin/test-n_neutral_separation20.txt")
data_neutral30 = np.loadtxt("../bin/test-n_neutral_separation30.txt")
data_neutral40 = np.loadtxt("../bin/test-n_neutral_separation40.txt")
data_neutral50 = np.loadtxt("../bin/test-n_neutral_separation50.txt")
data_neutral60 = np.loadtxt("../bin/test-n_neutral_separation60.txt")
#data_minus = np.loadtxt("../bin/test-n_minus.txt")

#x_plus = data_plus[:,0]
#y_plus = data_plus[:,1]

x_neutral = data_neutral[:,0]
y_neutral = data_neutral[:,1]

x_neutral20 = data_neutral20[:,0]
y_neutral20 = data_neutral20[:,1]

x_neutral30 = data_neutral30[:,0]
y_neutral30 = data_neutral30[:,1]

x_neutral40 = data_neutral40[:,0]
y_neutral40 = data_neutral40[:,1]

x_neutral50 = data_neutral50[:,0]
y_neutral50 = data_neutral50[:,1]

x_neutral60 = data_neutral60[:,0]
y_neutral60 = data_neutral60[:,1]

#x_minus = data_minus[:,0]
#y_minus = data_minus[:,1]

#err = data[:,2]

#A description of the available plotting characters and colours 
#to be used in place of 'rx' can be found here...
#http://matplotlib.org/api/pyplot_api.html#matplotlib.pyplot.plot
#To plot data with errorbars use 
#plt.errorbar(x,y,err,fmt='bo',label="Some Label")
#If you want to plot data but not show a key in the legend use
#label='_nolegend_'
#plt.plot(x_plus,y_plus,'rx',label=r'$n_{+}$')
plt.plot(x_neutral,y_neutral,'bo',label=r'$n^{b}_{0}\sigma=0.35; d = 10\sigma$')
plt.plot(x_neutral20 + 10,y_neutral20,'rx',label=r'$n^{b}_{0}\sigma=0.35; d = 20\sigma$')
plt.plot(x_neutral30 + 30,y_neutral30,'g^',label=r'$n^{b}_{0}\sigma=0.35; d = 30\sigma$')
plt.plot(x_neutral40 + 60,y_neutral40,'y*',label=r'$n^{b}_{0}\sigma=0.35; d = 40\sigma$')
plt.plot(x_neutral50 + 100,y_neutral50,'bs',label=r'$n^{b}_{0}\sigma=0.35; d = 50\sigma$')
plt.plot(x_neutral60 + 150,y_neutral60,'m>',label=r'$n^{b}_{0}\sigma=0.35; d = 60\sigma$')
#plt.plot(x_minus,y_minus,'gs',label=r'$n_{b}\sigma = 0.04$')

#Set the axis labels.  Labelpad option controls the spacing between actual axis and axis label.  The r option tells python to interpret as a raw string literal.
plt.xlabel(r"$h/\sigma$",labelpad=10)
plt.ylabel(r"$n\sigma^{3}$",labelpad=5)

#Set the axis limits if you want to specify a region of interest.
#The default is auto zoom.
plt.ylim(0.348,0.358)
plt.xlim(0,350)

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

savefig("density_profiles_at_varying_separation_nb35.pdf",bbox_inches='tight')

#Open a window and show the plot
plt.show()
