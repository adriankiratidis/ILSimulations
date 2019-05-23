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

#data_contactthm = np.loadtxt("../run_results/testing-normal-pressure-left-wall.txt5-6")
#data_negativederiv = np.loadtxt("../run_results/testing-negative_deriv_of_potential.txt5-6")

#data_contactthm2 = np.loadtxt("../run_results/compare_epsilonLJ/nocharge/20-100/testing-normal-pressure-left-wall.txt")
#data_negativederiv2 = np.loadtxt("../run_results/compare_epsilonLJ/nocharge/20-100/testing-negative_deriv_of_potential.txt")

#data_contactthm2 = np.loadtxt("../run_results/compare_epsilonLJ/charge4-12/35.51-88.78/testing6-potential-per-unit-area.txt")


ELJ_wall=["0", "17.76", "53.27", "71.02", "88.78"]
ELJ_wall_index=["", "2", "4", "5", "6"]

for i in range(len(ELJ_wall)):
        
    data_contactthm1 = np.loadtxt("../run_results/compare_epsilonLJ/long_range_attempt/charge5/35.51-"+ELJ_wall[i]+"/testing"+ELJ_wall_index[i]+"-potential-per-unit-area.txt")
    data_contactthm2 = np.loadtxt("../run_results/compare_epsilonLJ/long_range_attempt/charge6/35.51-"+ELJ_wall[i]+"/testing"+ELJ_wall_index[i]+"-potential-per-unit-area.txt")
    data_contactthm3 = np.loadtxt("../run_results/compare_epsilonLJ/long_range_attempt/charge7/35.51-"+ELJ_wall[i]+"/testing"+ELJ_wall_index[i]+"-potential-per-unit-area.txt")
    
    #x_contactthm = np.concatenate((data_contactthm[:,0], data_contactthm2[:,0]), axis=0)
    #y_contactthm = np.concatenate((data_contactthm[:,1], data_contactthm2[:,1]), axis=0)
    #y_contactthm = y_contactthm - y_contactthm[-1]
    #x_negativederiv = np.concatenate((data_negativederiv[:,0], data_negativederiv2[:,0]), axis=0)
    #y_negativederiv = np.concatenate((data_negativederiv[:,1] + y_contactthm[len(data_negativederiv[:,1]) - 1], data_negativederiv2[:,1]), axis=0)
    
    x_contactthm = data_contactthm1[:,0]
    y_contactthm = data_contactthm1[:,1]*(10**20)*2*pi
    y_contactthm = y_contactthm - y_contactthm[14]

    x_contactthm2 = data_contactthm2[:,0]
    y_contactthm2 = data_contactthm2[:,1]*(10**20)*2*pi
    y_contactthm2 = y_contactthm2 - y_contactthm2[36]

    x_contactthm3 = data_contactthm3[:,0]
    y_contactthm3 = data_contactthm3[:,1]*(10**20)*2*pi
    y_contactthm3 = y_contactthm3 - y_contactthm3[142]

    
    #x_negativederiv = data_negativederiv2[:,0]
    #y_negativederiv = data_negativederiv2[:,1]


    # for ij in range(len(y_contactthm)):
    # #!if(abs(y_contactthm[i]) <= 0 ):
    # #    print abs(y_contactthm[i])
    #     y_contactthm[ij] = math.log(abs(y_contactthm[ij]) if abs(y_contactthm[ij]) > 0 else 1)


            
    # for ik in range(len(y_contactthm2)):
    # #!if(abs(y_contactthm[i]) <= 0 ):
    # #    print abs(y_contactthm[i])
    #     y_contactthm2[ik] = math.log(abs(y_contactthm2[ik]) if abs(y_contactthm2[ik]) > 0 else 1)
  

    #y_contactthm = y_contactthm - y_contactthm[-1]
    
    plt.plot(x_contactthm[0:14],y_contactthm[0:14],'bo', label=r"\textrm{Calculated every }$5\sigma$")
    plt.plot(x_contactthm2[0:36],y_contactthm2[0:36],'rx', label=r"\textrm{Calculated every }$2\sigma$")
    plt.plot(x_contactthm3[0:142],y_contactthm3[0:142],'go', label=r"\textrm{Calculated every }$0.5\sigma$")
    
    #y_contactthm = y_contactthm * (y_negativederiv[0]/y_contactthm[0])
    #y_negativederiv = y_contactthm - y_contactthm[-1]
    

        
    
    #y_contactthm = math.log(y_contactthm)
    
    #plt.plot(x_contactthm[2:],y_contactthm[2:],'bo')
        
    #Set the axis labels.  Labelpad option controls the spacing between actual axis and axis label.  The r option tells python to interpret as a raw string literal.
    plt.xlabel(r"$z/\sigma$",labelpad=10)
    plt.ylabel(r"${2\pi\Delta{E_{s}}}(N/m)$",labelpad=5)
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
    
    #savefig("Potential_35.51-35.51-charge_a2_TEST_C4MIMTFSI_model2.pdf",bbox_inches='tight')
    #savefig("Potential_35.51-88.78-charge_a2_C4MIMTFSI_model2.pdf",bbox_inches='tight')
    savefig("Potential_35.51-"+ELJ_wall[i]+"_a2_C4MIMTFSI_model2_LONGRANGE.pdf",bbox_inches='tight')
    
    #Open a window and show the plot
    plt.show()
