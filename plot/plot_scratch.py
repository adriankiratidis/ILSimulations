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

separation = 212

charge = "minus"
chargenum = 1
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

# ".00000charge" + str(chargenum) + 
# data_1 = np.loadtxt("../bin/test-n_" + charge + "_separation" + str(separation) + "iteration31.txt")
# data_2 = np.loadtxt("../bin/test-n_" + charge + "_separation" + str(separation) + "iteration32.txt")
# data_3 = np.loadtxt("../bin/test-n_" + charge + "_separation" + str(separation) + "iteration33.txt")
# data_4 = np.loadtxt("../bin/test-n_" + charge + "_separation" + str(separation) + "iteration34.txt")
# data_5 = np.loadtxt("../bin/test-n_" + charge + "_separation" + str(separation) + "iteration35.txt")
# data_6 = np.loadtxt("../bin/test-n_" + charge + "_separation" + str(separation) + "iteration36.txt")
# data_7 = np.loadtxt("../bin/test-n_" + charge + "_separation" + str(separation) + "iteration37.txt")
# data_8 = np.loadtxt("../bin/test-n_" + charge + "_separation" + str(separation) + "iteration38.txt")
# data_9 = np.loadtxt("../bin/test-n_" + charge + "_separation" + str(separation) + "iteration39.txt")
# data_10 = np.loadtxt("../bin/test-n_" + charge + "_separation" + str(separation) + "iteration40.txt")
# data_11 = np.loadtxt("../bin/test-n_" + charge + "_separation" + str(separation) + "iteration41.txt")
# data_12 = np.loadtxt("../bin/test-n_" + charge + "_separation" + str(separation) + "iteration42.txt")

data_1 = np.loadtxt("../run_results/compare_epsilonLJ/charge/35.51-35.51/testing3-n_" + charge + "_separation" + str(separation) + ".00000charge" + str(chargenum) + "iteration604.txt")
data_2 = np.loadtxt("../run_results/compare_epsilonLJ/charge/35.51-35.51/testing3-n_" + charge + "_separation" + str(separation) + ".00000charge" + str(chargenum) + "iteration605.txt")
data_3 = np.loadtxt("../run_results/compare_epsilonLJ/charge/35.51-35.51/testing3-n_" + charge + "_separation" + str(separation) + ".00000charge" + str(chargenum) + "iteration606.txt")
data_4 = np.loadtxt("../run_results/compare_epsilonLJ/charge/35.51-35.51/testing3-n_" + charge + "_separation" + str(separation) + ".00000charge" + str(chargenum) + "iteration607.txt")
data_5 = np.loadtxt("../run_results/compare_epsilonLJ/charge/35.51-35.51/testing3-n_" + charge + "_separation" + str(separation) + ".00000charge" + str(chargenum) + "iteration608.txt")
data_6 = np.loadtxt("../run_results/compare_epsilonLJ/charge/35.51-35.51/testing3-n_" + charge + "_separation" + str(separation) + ".00000charge" + str(chargenum) + "iteration609.txt")
data_7 = np.loadtxt("../run_results/compare_epsilonLJ/charge/35.51-35.51/testing3-n_" + charge + "_separation" + str(separation) + ".00000charge" + str(chargenum) + "iteration610.txt")
# data_8 = np.loadtxt("../testing-n_" + charge + "_separation" + str(separation) + ".00000charge" + str(chargenum) + "iteration8.txt")
# data_9 = np.loadtxt("../testing-n_" + charge + "_separation" + str(separation) + ".00000charge" + str(chargenum) + "iteration9.txt")
# data_10 = np.loadtxt("../testing-n_" + charge + "_separation" + str(separation) + ".00000charge" + str(chargenum) + "iteration10.txt")
# data_11 = np.loadtxt("../testing-n_" + charge + "_separation" + str(separation) + ".00000charge" + str(chargenum) + "iteration11.txt")
# data_12 = np.loadtxt("../testing-n_" + charge + "_separation" + str(separation) + ".00000charge" + str(chargenum) + "iteration12.txt")



# data_4 = np.loadtxt("../testing-n_" + charge + "_separation" + str(separation) + "000charge1iteration84.txt")
# data_5 = np.loadtxt("../testing-n_" + charge + "_separation" + str(separation) + "000charge1iteration85.txt")
# data_6 = np.loadtxt("../testing-n_" + charge + "_separation" + str(separation) + "000charge1iteration86.txt")
# data_7 = np.loadtxt("../testing-n_" + charge + "_separation" + str(separation) + "000charge1iteration87.txt")
# data_8 = np.loadtxt("../testing-n_" + charge + "_separation" + str(separation) + "000charge1iteration88.txt")
# data_9 = np.loadtxt("../testing-n_" + charge + "_separation" + str(separation) + "000charge1iteration89.txt")
# data_10 = np.loadtxt("../testing-n_" + charge + "_separation" + str(separation) + ".00000charge" + str(chargenum) + "iteration10.txt")
# data_11 = np.loadtxt("../testing-n_" + charge + "_separation" + str(separation) + ".00000charge" + str(chargenum) + "iteration11.txt")
# data_12 = np.loadtxt("../testing-n_" + charge + "_separation" + str(separation) + ".00000charge" + str(chargenum) + "iteration12.txt")
# data_13 = np.loadtxt("../testing-n_" + charge + "_separation" + str(separation) + ".00000charge" + str(chargenum) + "iteration13.txt")
# data_14 = np.loadtxt("../testing-n_" + charge + "_separation" + str(separation) + ".00000charge" + str(chargenum) + "iteration14.txt")
# data_15 = np.loadtxt("../testing-n_" + charge + "_separation" + str(separation) + ".00000charge" + str(chargenum) + "iteration15.txt")

#x_0 = data_0[:,0]
#y_0 = data_0[:,1]

x_1 = data_1[:,0]
y_1 = data_1[:,1]

x_2 = data_2[:,0]
y_2 = data_2[:,1]

x_3 = data_3[:,0]
y_3 = data_3[:,1]

x_4 = data_4[:,0]
y_4 = data_4[:,1]

x_5 = data_5[:,0]
y_5 = data_5[:,1]


x_6 = data_6[:,0]
y_6 = data_6[:,1]

x_7 = data_7[:,0]
y_7 = data_7[:,1]

# x_8 = data_8[:,0]
# y_8 = data_8[:,1]

# x_9 = data_9[:,0]
# y_9 = data_9[:,1]

# x_10 = data_10[:,0]
# y_10 = data_10[:,1]


# x_11 = data_11[:,0]
# y_11 = data_11[:,1]

# x_12 = data_12[:,0]
# y_12 = data_12[:,1]

# x_13 = data_13[:,0]
# y_13 = data_13[:,1]

# x_14 = data_14[:,0]
# y_14 = data_14[:,1]

# x_15 = data_15[:,0]
# y_15 = data_15[:,1]



# x_6 = data_6[:,0]
# y_6 = data_6[:,1]

# x_7 = data_7[:,0]
# y_7 = data_7[:,1]

# x_8 = data_8[:,0]
# y_8 = data_8[:,1]

# x_9 = data_9[:,0]
# y_9 = data_9[:,1]


#A description of the available plotting characters and colours 
#to be used in place of 'rx' can be found here...
#http://matplotlib.org/api/pyplot_api.html#matplotlib.pyplot.plot
#To plot data with errorbars use 
#plt.errorbar(x,y,err,fmt='bo',label="Some Label")
#If you want to plot data but not show a key in the legend use
#label='_nolegend_'
#plt.plot(x_0,y_0,'rx',label='converged profile')
plt.plot(x_1,y_1,'bx',label='iteration1')
plt.plot(x_2,y_2,'go',label='iteration2')
plt.plot(x_3,y_3,'rs',label='iteration3')
plt.plot(x_4,y_4,'ch',label='iteration4')
plt.plot(x_5+0.2,y_5,'m*',label='iteration5')
plt.plot(x_6+0.25,y_6,'kx',label='iteration6')
plt.plot(x_7,y_7,'gx',label='iteration7')
# plt.plot(x_8,y_8,'ro',label='iteration8')
# plt.plot(x_9,y_9,'ch',label='iteration9')
#plt.plot(x_10,y_10,'m*',label='iteration10')
#plt.plot(x_11,y_11,'bx',label='iteration11')
#plt.plot(x_12,y_12,'go',label='iteration12')
# plt.plot(x_13,y_13,'rs',label='iteration13')
# plt.plot(x_14,y_14,'ch',label='iteration14')
# plt.plot(x_15,y_15,'m*',label='iteration15')


#Set the axis labels.  Labelpad option controls the spacing between actual axis and axis label.  The r option tells python to interpret as a raw string literal.
plt.xlabel(r"$z/\sigma$",labelpad=10)
plt.ylabel(r"$n\sigma^{3}$",labelpad=5)

#Set the axis limits if you want to specify a region of interest.
#The default is auto zoom.
#plt.ylim(-0.02,0.05)
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
