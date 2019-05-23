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

separation = 10

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
#data_0 = np.loadtxt("../testing-n_" + charge + "_separation" + str(separation) + ".00000.txt")

data_plus_1 = np.loadtxt("../testing-n_plus_separation" + str(separation) + ".00000charge1.txt")
# data_plus_2 = np.loadtxt("../testing-n_plus_separation" + str(separation) + ".00000charge2.txt")
# data_plus_3 = np.loadtxt("../testing-n_plus_separation" + str(separation) + ".00000charge3.txt")
# data_plus_4 = np.loadtxt("../testing-n_plus_separation" + str(separation) + ".00000charge4.txt")
# data_plus_5 = np.loadtxt("../testing-n_plus_separation" + str(separation) + ".00000charge5.txt")
# data_plus_6 = np.loadtxt("../testing-n_plus_separation" + str(separation) + ".00000charge6.txt")
# data_plus_7 = np.loadtxt("../testing-n_plus_separation" + str(separation) + ".00000charge7.txt")
# data_plus_8 = np.loadtxt("../testing-n_plus_separation" + str(separation) + ".00000charge8.txt")


data_minus_1 = np.loadtxt("../testing-n_minus_separation" + str(separation) + ".00000charge1.txt")
# data_minus_2 = np.loadtxt("../testing-n_minus_separation" + str(separation) + ".00000charge2.txt")
# data_minus_3 = np.loadtxt("../testing-n_minus_separation" + str(separation) + ".00000charge3.txt")
# data_minus_4 = np.loadtxt("../testing-n_minus_separation" + str(separation) + ".00000charge4.txt")
# data_minus_5 = np.loadtxt("../testing-n_minus_separation" + str(separation) + ".00000charge5.txt")
# data_minus_6 = np.loadtxt("../testing-n_minus_separation" + str(separation) + ".00000charge6.txt")
# data_minus_7 = np.loadtxt("../testing-n_minus_separation" + str(separation) + ".00000charge7.txt")
# data_minus_8 = np.loadtxt("../testing-n_minus_separation" + str(separation) + ".00000charge8.txt")


x_p_1 = data_plus_1[:,0]
y_p_1 = data_plus_1[:,1]

# x_p_2 = data_plus_2[:,0]
# y_p_2 = data_plus_2[:,1]

# x_p_3 = data_plus_3[:,0]
# y_p_3 = data_plus_3[:,1]

# x_p_4 = data_plus_4[:,0]
# y_p_4 = data_plus_4[:,1]

# x_p_5 = data_plus_5[:,0]
# y_p_5 = data_plus_5[:,1]

# x_p_6 = data_plus_6[:,0]
# y_p_6 = data_plus_6[:,1]

# x_p_7 = data_plus_7[:,0]
# y_p_7 = data_plus_7[:,1]

# x_p_8 = data_plus_8[:,0]
# y_p_8 = data_plus_8[:,1]




x_m_1 = data_minus_1[:,0]
y_m_1 = data_minus_1[:,1]

# x_m_2 = data_minus_2[:,0]
# y_m_2 = data_minus_2[:,1]

# x_m_3 = data_minus_3[:,0]
# y_m_3 = data_minus_3[:,1]

# x_m_4 = data_minus_4[:,0]
# y_m_4 = data_minus_4[:,1]

# x_m_5 = data_minus_5[:,0]
# y_m_5 = data_minus_5[:,1]

# x_m_6 = data_minus_6[:,0]
# y_m_6 = data_minus_6[:,1]

# x_m_7 = data_minus_7[:,0]
# y_m_7 = data_minus_7[:,1]

# x_m_8 = data_minus_8[:,0]
# y_m_8 = data_minus_8[:,1]





#A description of the available plotting characters and colours 
#to be used in place of 'rx' can be found here...
#http://matplotlib.org/api/pyplot_api.html#matplotlib.pyplot.plot
#To plot data with errorbars use 
#plt.errorbar(x,y,err,fmt='bo',label="Some Label")
#If you want to plot data but not show a key in the legend use
#label='_nolegend_'
plt.plot(x_p_1+0.02,y_p_1,'rx',label='+ve 1st')
# plt.plot(x_p_2,y_p_2,'gx',label='+ve 2nd')
# plt.plot(x_p_3,y_p_3,'bx',label='+ve 3rd')
# plt.plot(x_p_4,y_p_4,'yx',label='+ve 4th')
# plt.plot(x_p_5,y_p_5,'kx',label='+ve 5th')
# plt.plot(x_p_6,y_p_6,'rx',label='+ve 6rd')
# plt.plot(x_p_7,y_p_7,'gx',label='+ve 7th')
# plt.plot(x_p_8,y_p_8,'bx',label='+ve 8th')


plt.plot(x_m_1,y_m_1,'ro',label='-ve 1st')
# plt.plot(x_m_2,y_m_2,'go',label='-ve 2nd')
# plt.plot(x_m_3,y_m_3,'bo',label='-ve 3rd')
# plt.plot(x_m_4,y_m_4,'yo',label='-ve 4th')
# plt.plot(x_m_5,y_m_5,'ko',label='-ve 5th')
# plt.plot(x_m_6,y_m_6,'ro',label='-ve 6rd')
# plt.plot(x_m_7,y_m_7,'go',label='-ve 7th')
# plt.plot(x_m_8,y_m_8,'bo',label='-ve 8th')


#plt.plot(x_0,y_0,'rx',label='converged profile')
#plt.plot(x_1+0.02,y_1,'bx',label=r'$n_{+}$')
#plt.plot(x_2+0.04,y_2,'go',label=r'$n_{0}$')
#plt.plot(x_3,y_3,'rs',label=r'$n_{-}$')
#plt.plot(x_4,y_4,'ch',label=r'$n_{s}$')
#plt.plot(x_5,y_5,'m*',label='iteration 5')
#plt.plot(x_6,y_6,'y<',label='iteration 6')
# plt.plot(x_7,y_7,'k>',label='iteration 7')
# plt.plot(x_8,y_8,'w^',label='iteration 8')
# plt.plot(x_9,y_9,'bD',label='iteration 9')

#Set the axis labels.  Labelpad option controls the spacing between actual axis and axis label.  The r option tells python to interpret as a raw string literal.
plt.xlabel(r"$z/\sigma$",labelpad=10)
plt.ylabel(r"$n\sigma^{3}$",labelpad=5)

#Set the axis limits if you want to specify a region of interest.
#The default is auto zoom.
#plt.ylim(-0.05,1.1)
#plt.xlim(0.5,4.5)

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

savefig("Dimer-iterations.pdf",bbox_inches='tight')

#Open a window and show the plot
plt.show()
