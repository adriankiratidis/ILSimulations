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


charges=["CentreCentrePotential/minuswalls"]
#charges=["nocharge4-12"]
ELJ_wall=["100"]
ELJ_wall_index=[""]
particle_wall_interaction="0"
density = "0.005"

for j in range(len(charges)):
    for i in range(len(ELJ_wall)):
        data_contactthm2_1 = np.loadtxt("../run_results/compare_epsilonLJ/CentreCentrePotential/minuswalls/35.51-0_Heptamer_SingleSphere-0.005-10--0.003125-hs1.0CHARGE/testing2-a2-potential-per-unit-area.txt")
        data_contactthm2_2 = np.loadtxt("../run_results/compare_epsilonLJ/CentreCentrePotential/minuswalls/35.51-35.51_Heptamer_SingleSphere-0.005-10--0.003125-hs1.0CHARGE/testing2-a2-potential-per-unit-area.txt")
        data_contactthm2_3 = np.loadtxt("../run_results/compare_epsilonLJ/CentreCentrePotential/minuswalls/35.51-53.27_Heptamer_SingleSphere-0.005-10--0.003125-hs1.0CHARGE/testing2-a2-potential-per-unit-area.txt")
        data_contactthm2_4 = np.loadtxt("../run_results/compare_epsilonLJ/CentreCentrePotential/minuswalls/35.51-88.78_Heptamer_SingleSphere-0.005-10--0.003125-hs1.0CHARGE/testing2-a2-potential-per-unit-area.txt")


        
        x_contactthm_1 = data_contactthm2_1[:,0]
        y_contactthm_1 = data_contactthm2_1[:,1]*(10**20)*2*pi
        y_contactthm_1 = y_contactthm_1 - y_contactthm_1[-1]
        
        x_contactthm_2 = data_contactthm2_2[:,0]
        y_contactthm_2 = data_contactthm2_2[:,1]*(10**20)*2*pi
        y_contactthm_2 = y_contactthm_2 - y_contactthm_2[-1]
        
        x_contactthm_3 = data_contactthm2_3[:,0]
        y_contactthm_3 = data_contactthm2_3[:,1]*(10**20)*2*pi
        y_contactthm_3 = y_contactthm_3 - y_contactthm_3[-1]
        
        x_contactthm_4 = data_contactthm2_4[:,0]
        y_contactthm_4 = data_contactthm2_4[:,1]*(10**20)*2*pi
        y_contactthm_4 = y_contactthm_4 - y_contactthm_4[-1]
        
        
        plt.plot(x_contactthm_1[2:len(y_contactthm_1)],(y_contactthm_1[2:len(y_contactthm_1)]),'bo', label='Potential C = 10, epsilon-pw = 0')
        plt.plot(x_contactthm_2[2:len(y_contactthm_2)],(y_contactthm_2[2:len(y_contactthm_2)]),'rx', label='Potential C = 10, epsilon-pw = 35.51')
        plt.plot(x_contactthm_3[2:len(y_contactthm_3)],(y_contactthm_3[2:len(y_contactthm_3)]),'gs', label='Potential C = 10, epsilon-pw = 53.27')
        plt.plot(x_contactthm_4[2:len(y_contactthm_4)],(y_contactthm_4[2:len(y_contactthm_4)]),'k^', label='Potential C = 10, epsilon-pw = 88.78')

                
        #Set the axis labels.  Labelpad option controls the spacing between actual axis and axis label.  The r option tells python to interpret as a raw string literal.
        plt.xlabel(r"$h/\sigma$",labelpad=10)
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

        #savefig("Potential-Heptamer_SingleSphere-density"+ density + "35.51-"+particle_wall_interaction+".pdf",bbox_inches='tight')
        
        # savefig("Potential-Heptamer_SingleSphere-density"+ density + "35.51-"+particle_wall_interaction+"at_maxima2.pdf",bbox_inches='tight')
        savefig("Potential-Heptamer_SingleSphere-density"+ density + "35.51-fixedC_different_epsilon_negativewalls.pdf",bbox_inches='tight')
        
        #Open a window and show the plot
        plt.show()
