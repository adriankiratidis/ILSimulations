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


charges=["CentreCentrePotential"]
#charges=["nocharge4-12"]
ELJ_wall=["100"]
ELJ_wall_index=[""]

for j in range(len(charges)):
    for i in range(len(ELJ_wall)):
        #data_contactthm2 = np.loadtxt("../run_results/compare_epsilonLJ/"+charges[j]+"/100-"+ELJ_wall[i]+"_C4MIM_BF4--TESTING-NEW-HS/testing"+ELJ_wall_index[i]+"-a2-potential-per-unit-area.txt")
        #data_contactthm2 = np.loadtxt("../run_results/compare_epsilonLJ/"+charges[j]+"/35.51-53.27_C4MIM_BF4--0.005-0-0.00-LongTEST10/testing"+ELJ_wall_index[i]+"-a2-potential-per-unit-area.txt")
        data_contactthm2 = np.loadtxt("../run_results/compare_epsilonLJ/"+charges[j]+"/35.51-53.27_C4MIM_BF4--0.005-0-0.00-LongLong/testing2"+ELJ_wall_index[i]+"-a2-potential-per-unit-area.txt")


        

# for j in range(len(charges)):
#     for i in range(len(IL)):
#         data_contactthm2 = np.loadtxt("../run_results/compare_epsilonLJ/"+charges[j]+"4-12/35.51-"+ELJ_wall[0]+IL[i]+"/testing"+IL_index[i]+"-a2-potential-per-unit-area.txt")
        
        x_contactthm = data_contactthm2[:,0]
        y_contactthm = data_contactthm2[:,1]*(10**20)*2*pi
        y_contactthm = y_contactthm - y_contactthm[-1]

        
        #A description of the available plotting characters and colours 
        #to be used in place of 'rx' can be found here...
        #http://matplotlib.org/api/pyplot_api.html#matplotlib.pyplot.plot
        #To plot data with errorbars use 
        #plt.errorbar(x,y,err,fmt='bo',label="Some Label")
        #If you want to plot data but not show a key in the legend use
        #label='_nolegend_'
        #plt.plot(x_plus,y_plus,'rx',label=r'$n_{+}$')
        

        #y_contactthm = y_contactthm - y_contactthm[-1]
        
        
        
        #y_contactthm = y_contactthm * (y_negativederiv[0]/y_contactthm[0])
        #y_negativederiv = y_contactthm - y_contactthm[-1]
        
        #for ij in range(len(y_contactthm)):
            #!if(abs(y_contactthm[i]) <= 0 ):
            #    print abs(y_contactthm[i])
        #    y_contactthm[ij] = math.log(abs(y_contactthm[ij]) if abs(y_contactthm[ij]) > 0 else 1)
        
        #y_contactthm = math.log(y_contactthm)
        
        #plt.plot(x_contactthm[2:],y_contactthm[2:],'bo')
        plt.plot(x_contactthm[2:],y_contactthm[2:],'bo')

        
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
        
        
        #savefig("Potential_35.51-"+ELJ_wall[i]+"-"+charges[j]+"_C4MIM+_TFSI-_model1_density2_plus_halfplus.pdf",bbox_inches='tight')
        #savefig("Potential_35.51-"+ELJ_wall[0]+"-"+charges[j]+IL[i]+"_density1_a2.pdf",bbox_inches='tight')
        
        #Open a window and show the plot
        plt.show()
