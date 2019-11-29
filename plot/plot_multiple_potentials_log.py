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

for j in range(len(charges)):
    for i in range(len(ELJ_wall)):
        data_contactthm2_0 = np.loadtxt("../run_results/compare_epsilonLJ/"+charges[j]+"/35.51-" + particle_wall_interaction + "_PositiveMinusSpheres-0.005-0--0.003125-hs1.0CHARGE/testing"+ELJ_wall_index[i]+"-a2-potential-per-unit-area.txt")
        data_contactthm2_1 = np.loadtxt("../run_results/compare_epsilonLJ/"+charges[j]+"/35.51-" + particle_wall_interaction + "_PositiveMinusSpheres-0.0005-0--0.003125-hs1.0CHARGE/testing"+ELJ_wall_index[i]+"-a2-potential-per-unit-area.txt")
        data_contactthm2_2 = np.loadtxt("../run_results/compare_epsilonLJ/"+charges[j]+"/35.51-" + particle_wall_interaction + "_PositiveMinusSpheres-0.00005-0--0.003125-hs1.0CHARGE/testing"+ELJ_wall_index[i]+"-a2-potential-per-unit-area.txt")
        data_contactthm2_3 = np.loadtxt("../run_results/compare_epsilonLJ/"+charges[j]+"/35.51-" + particle_wall_interaction + "_PositiveMinusSpheres-0.000005-0--0.003125-hs1.0CHARGE/testing"+ELJ_wall_index[i]+"-a2-potential-per-unit-area.txt")



        data_contactthm2_7 = np.loadtxt("../run_results/compare_epsilonLJ/"+charges[j]+"/35.51-" + particle_wall_interaction + "_PositiveMinusSpheres-0.005-20--0.003125-hs1.0CHARGE/testing2"+ELJ_wall_index[i]+"-a2-potential-per-unit-area.txt")
        data_contactthm2_8 = np.loadtxt("../run_results/compare_epsilonLJ/"+charges[j]+"/35.51-" + particle_wall_interaction + "_PositiveMinusSpheres-0.0005-100--0.003125-hs1.0CHARGE/testing6"+ELJ_wall_index[i]+"-a2-potential-per-unit-area.txt")
        data_contactthm2_9 = np.loadtxt("../run_results/compare_epsilonLJ/"+charges[j]+"/35.51-" + particle_wall_interaction + "_PositiveMinusSpheres-0.00005-100--0.003125-hs1.0CHARGE/testing6"+ELJ_wall_index[i]+"-a2-potential-per-unit-area.txt")
        data_contactthm2_10 = np.loadtxt("../run_results/compare_epsilonLJ/"+charges[j]+"/35.51-" + particle_wall_interaction + "_PositiveMinusSpheres-0.000005-100--0.003125-hs1.0CHARGE/testing2"+ELJ_wall_index[i]+"-a2-potential-per-unit-area.txt")

        
        x_contactthm_0 = data_contactthm2_0[:,0]
        y_contactthm_0 = data_contactthm2_0[:,1]*(10**20)*2*pi
        y_contactthm_0 = y_contactthm_0 - y_contactthm_0[-1]
        
        x_contactthm_1 = data_contactthm2_1[:,0]
        y_contactthm_1 = data_contactthm2_1[:,1]*(10**20)*2*pi
        y_contactthm_1 = y_contactthm_1 - y_contactthm_1[-1]
        
        x_contactthm_2 = data_contactthm2_2[:,0]
        y_contactthm_2 = data_contactthm2_2[:,1]*(10**20)*2*pi
        y_contactthm_2 = y_contactthm_2 - y_contactthm_2[-1]
        
        x_contactthm_3 = data_contactthm2_3[:,0]
        y_contactthm_3 = data_contactthm2_3[:,1]*(10**20)*2*pi
        y_contactthm_3 = y_contactthm_3 - y_contactthm_3[-1]
        
        # x_contactthm_4 = data_contactthm2_4[:,0]
        # y_contactthm_4 = data_contactthm2_4[:,1]*(10**20)*2*pi
        # y_contactthm_4 = y_contactthm_4 - y_contactthm_4[-1]
        
        # x_contactthm_5 = data_contactthm2_5[:,0]
        # y_contactthm_5 = data_contactthm2_5[:,1]*(10**20)*2*pi
        # y_contactthm_5 = y_contactthm_5 - y_contactthm_5[-1]
        
        # x_contactthm_6 = data_contactthm2_6[:,0]
        # y_contactthm_6 = data_contactthm2_6[:,1]*(10**20)*2*pi
        # y_contactthm_6 = y_contactthm_6 - y_contactthm_6[-1]
        
        x_contactthm_7 = data_contactthm2_7[:,0]
        y_contactthm_7 = data_contactthm2_7[:,1]*(10**20)*2*pi
        y_contactthm_7 = y_contactthm_7 - y_contactthm_7[-1]

        x_contactthm_8 = data_contactthm2_8[:,0]
        y_contactthm_8 = data_contactthm2_8[:,1]*(10**20)*2*pi
        y_contactthm_8 = y_contactthm_8 - y_contactthm_8[-1]

        x_contactthm_9 = data_contactthm2_9[:,0]
        y_contactthm_9 = data_contactthm2_9[:,1]*(10**20)*2*pi
        y_contactthm_9 = y_contactthm_9 - y_contactthm_9[-1]

        x_contactthm_10 = data_contactthm2_10[:,0]
        y_contactthm_10 = data_contactthm2_10[:,1]*(10**20)*2*pi
        y_contactthm_10 = y_contactthm_10 - y_contactthm_10[-1]

        # x_contactthm_11 = data_contactthm2_11[:,0]
        # y_contactthm_11 = data_contactthm2_11[:,1]*(10**20)*2*pi
        # y_contactthm_11 = y_contactthm_11 - y_contactthm_11[-1]

        # x_contactthm_12 = data_contactthm2_12[:,0]
        # y_contactthm_12 = data_contactthm2_12[:,1]*(10**20)*2*pi
        # y_contactthm_12 = y_contactthm_12 - y_contactthm_12[-1]

        
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

        #print type(y_contactthm_0)
        
        #print type
        
        #plt.plot(x_contactthm_0[2:len(y_contactthm_0)-2],(np.log(np.array(np.absolute(y_contactthm_0[2:len(y_contactthm_0)-2]), dtype=float32))),'bo', label='Potential C = 0, n=0.005')
        #plt.plot(x_contactthm_1[2:len(y_contactthm_1)-2],(np.log(np.array(np.absolute(y_contactthm_1[2:len(y_contactthm_1)-2]), dtype=float32))),'rx', label='Potential C = 0, n=0.0005')
        #plt.plot(x_contactthm_2[2:len(y_contactthm_0)-2],(np.log(np.array(np.absolute(y_contactthm_2[2:len(y_contactthm_2)-2]), dtype=float32))),'gs', label='Potential C = 0, n=0.00005')
        #plt.plot(x_contactthm_3[2:len(y_contactthm_0)-2],(np.log(np.array(np.absolute(y_contactthm_3[2:len(y_contactthm_3)-2]), dtype=float32))),'k^', label='Potential C = 0, n=0.000005')
        
        # plt.plot(x_contactthm_4[2:len(y_contactthm_0)-2],(np.log(np.array(np.absolute(y_contactthm_4[2:len(y_contactthm_4)-2]), dtype=float32))),'c<', label='Potential C = 5')
        # plt.plot(x_contactthm_5[2:len(y_contactthm_0)-2],(np.log(np.array(np.absolute(y_contactthm_5[2:len(y_contactthm_5)-2]), dtype=float32))),'m>', label='Potential C = 10')
        # plt.plot(x_contactthm_6[2:len(y_contactthm_0)-2],(np.log(np.array(np.absolute(y_contactthm_6[2:len(y_contactthm_6)-2]), dtype=float32))),'yH', label='Potential C = 20')
        
        plt.plot(x_contactthm_7[2:len(y_contactthm_0)-2],(np.log(np.array(np.absolute(y_contactthm_7[2:len(y_contactthm_7)-2]), dtype=float32))),'g+', label='Potential C = 20, n=0.005')
        plt.plot(x_contactthm_8[2:len(y_contactthm_8)-2],(np.log(np.array(np.absolute(y_contactthm_8[2:len(y_contactthm_8)-2]), dtype=float32))),'c<', label='Potential C = 100, n=0.0005')
        plt.plot(x_contactthm_9[2:len(y_contactthm_9)-2],(np.log(np.array(np.absolute(y_contactthm_9[2:len(y_contactthm_9)-2]), dtype=float32))),'yH', label='Potential C = 100, n=0.00005')
        plt.plot(x_contactthm_10[2:len(y_contactthm_10)-2],(np.log(np.array(np.absolute(y_contactthm_10[2:len(y_contactthm_10)-2]), dtype=float32))),'ro', label='Potential C = 100, n=0.000005')

        # plt.plot(x_contactthm_11[2:],y_contactthm_11[2:],'mH', label='Potential C = 140')
        # plt.plot(x_contactthm_12[2:],y_contactthm_12[2:],'y-', label='Potential C = 160')
                
        #Set the axis labels.  Labelpad option controls the spacing between actual axis and axis label.  The r option tells python to interpret as a raw string literal.
        plt.xlabel(r"$h/\sigma$",labelpad=10)
        plt.ylabel(r"${2\pi\log|\Delta{E_{s}}}|(N/m)$",labelpad=5)
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        
        #Set the axis limits if you want to specify a region of interest.
        #The default is auto zoom.
        plt.ylim(-28.0,-2)
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
        
        #savefig("Log_Potential_C4MIM_BF4_Negative_Wall_charge_C_Comparison-pw" + particle_wall_interaction + ".pdf",bbox_inches='tight')
        #savefig("Potential_35.51-"+ELJ_wall[i]+"-"+charges[j]+"_C4MIM+_TFSI-_model1_density2_plus_halfplus.pdf",bbox_inches='tight')
        #savefig("Potential_35.51-"+ELJ_wall[0]+"-"+charges[j]+IL[i]+"_density1_a2.pdf",bbox_inches='tight')
        savefig("Log_Potential_PositiveMinusSpheres_35.51-" + particle_wall_interaction + "Cnonzero_NegativeWallCharge-0.003125.pdf",bbox_inches='tight')
        
        #Open a window and show the plot
        plt.show()
