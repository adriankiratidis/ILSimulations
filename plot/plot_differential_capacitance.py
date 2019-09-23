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

data = np.loadtxt("/home/adriank/postdoc/2018/branch/ILSimulations/testing14-a2-electric_potential_and_charge.txt")

separation = data[:, 0]
sigma_left = data[:, 1]
sigma_right = data[:, 2]
potential = data[:, 3]
differential_capacitance = list(data[:, 3])
#print(range(len(potential))[1:len(potential) - 1]) 
#print(sigma_left)
#print(potential)
#print("")
for i in range(len(potential))[1:len(potential) - 1]:
    #print("calculating")
    #print(i)
    differential_capacitance[i] = (sigma_left[i+1] - sigma_left[i-1])/(potential[i+1] - potential[i-1])
    #print((sigma_left[i+1] - sigma_left[i-1])/(potential[i+1] - potential[i-1]))
    
print(potential[1:len(potential) - 1])
print("")
print(differential_capacitance[1:len(differential_capacitance) - 1])
plt.plot(potential[1:len(potential) - 1], differential_capacitance[1:len(differential_capacitance) - 1],'bo', label='DC')
#Set the axis labels.  Labelpad option controls the spacing between actual axis and axis label.  The r option tells python to interpret as a raw string literal.
plt.xlabel(r"U (V)",labelpad=10)
plt.ylabel(r"DC",labelpad=5)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

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

#savefig("DC.pdf",bbox_inches='tight')
#savefig("Potential_35.51-"+ELJ_wall[i]+"-"+charges[j]+"_C4MIM+_TFSI-_model1_density2_plus_halfplus.pdf",bbox_inches='tight')
#savefig("Potential_35.51-"+ELJ_wall[0]+"-"+charges[j]+IL[i]+"_density1_a2.pdf",bbox_inches='tight')

#Open a window and show the plot
plt.show()
