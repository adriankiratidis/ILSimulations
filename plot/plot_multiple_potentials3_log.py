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
        data_contactthm2_1 = np.loadtxt("../run_results/compare_epsilonLJ/CentreCentrePotential/minuswalls/35.51-" + particle_wall_interaction + "_C4MIM+_TFSI-_model1-0.005-0-r4-0.0-hs1.0CHARGE-LONGRANGE30/testing3-a2-potential-per-unit-area.txt")
        #data_contactthm2_2 = np.loadtxt("../run_results/compare_epsilonLJ/CentreCentrePotential/minuswalls/35.51-" + particle_wall_interaction + "_C4MIM+_TFSI-_model1-0.005-0-r4--0.003125-hs1.0CHARGE-LONGRANGE30/testing-a2-potential-per-unit-area.txt")

        data_contactthm2_12 = np.loadtxt("../run_results/compare_epsilonLJ/CentreCentrePotential/minuswalls/35.51-" + particle_wall_interaction + "_C4MIM+_TFSI-_model1-0.005-60-r4-0.0-hs1.0CHARGE-LONGRANGE30/testing4-a2-potential-per-unit-area.txt")
        #data_contactthm2_22 = np.loadtxt("../run_results/compare_epsilonLJ/CentreCentrePotential/minuswalls/35.51-" + particle_wall_interaction + "_C4MIM+_TFSI-_model1-0.005-160-r4--0.003125-hs1.0CHARGE-LONGRANGE30/testing2-a2-potential-per-unit-area.txt")





        
        #data_contactthm2_3 = np.loadtxt("../run_results/compare_epsilonLJ/CentreCentrePotential/minuswalls/35.51-" + particle_wall_interaction + "_C4MIM+_TFSI-_model1-0.005-0--0.003125-hs1.0CHARGE/testing2-a2-potential-per-unit-area.txt")
        #data_contactthm2_4 = np.loadtxt("../run_results/compare_epsilonLJ/CentreCentrePotential/minuswalls/35.51-" + particle_wall_interaction + "_C4MIM+_TFSI-_model1-0.005-0.9524--0.003125/testing4-a2-potential-per-unit-area.txt")
        #data_contactthm2_5 = np.loadtxt("../run_results/compare_epsilonLJ/CentreCentrePotential/minuswalls/35.51-" + particle_wall_interaction + "_C4MIM+_TFSI-_model1-0.005-2.3810--0.003125/testing5-a2-potential-per-unit-area.txt")
        #data_contactthm2_6 = np.loadtxt("../run_results/compare_epsilonLJ/CentreCentrePotential/minuswalls/35.51-" + particle_wall_interaction + "_C4MIM+_TFSI-_model1-0.005-4.7619--0.003125/testing6-a2-potential-per-unit-area.txt")
        #data_contactthm2_7 = np.loadtxt("../run_results/compare_epsilonLJ/CentreCentrePotential/minuswalls/35.51-" + particle_wall_interaction + "_C4MIM+_TFSI-_model1-0.005-9.5238--0.003125/testing7-a2-potential-per-unit-area.txt")
        #data_contactthm2_8 = np.loadtxt("../run_results/compare_epsilonLJ/CentreCentrePotential/minuswalls/35.51-" + particle_wall_interaction + "_C4MIM+_TFSI-_model1-0.005-19.0476--0.003125/testing8-a2-potential-per-unit-area.txt")
        #data_contactthm2_9 = np.loadtxt("../run_results/compare_epsilonLJ/CentreCentrePotential/minuswalls/35.51-" + particle_wall_interaction + "_C4MIM+_TFSI-_model1-0.005-28.5714--0.003125/testing9-a2-potential-per-unit-area.txt")
        #data_contactthm2_10 = np.loadtxt("../run_results/compare_epsilonLJ/CentreCentrePotential/minuswalls/35.51-" + particle_wall_interaction + "_C4MIM+_TFSI-_model1-0.005-38.0952--0.003125/testing10-a2-potential-per-unit-area.txt")
        #data_contactthm2_11 = np.loadtxt("../run_results/compare_epsilonLJ/CentreCentrePotential/minuswalls/35.51-" + particle_wall_interaction + "_C4MIM+_TFSI-_model1-0.005-57.1429--0.003125/testing11-a2-potential-per-unit-area.txt")
        #data_contactthm2_12 = np.loadtxt("../run_results/compare_epsilonLJ/CentreCentrePotential/minuswalls/35.51-" + particle_wall_interaction + "_C4MIM+_TFSI-_model1-0.005-66.6667--0.003125/testing12-a2-potential-per-unit-area.txt")
        #data_contactthm2_13 = np.loadtxt("../run_results/compare_epsilonLJ/CentreCentrePotential/minuswalls/35.51-" + particle_wall_interaction + "_C4MIM+_TFSI-_model1-0.005-76.1908--0.003125/testing13-a2-potential-per-unit-area.txt")
 

        
# for j in range(len(charges)):76.1908
#     for i in range(len(IL)):
#         data_contactthm2 = np.loadtxt("../run_results/compare_epsilonLJ/"+charges[j]+"4-12/35.51-"+ELJ_wall[0]+IL[i]+"/testing"+IL_index[i]+"-a2-potential-per-unit-area.txt")
        
        x_contactthm_1 = data_contactthm2_1[:,0]
        y_contactthm_1 = data_contactthm2_1[:,1]*(10**20)*2*pi
        y_contactthm_1 = y_contactthm_1 - y_contactthm_1[-1]
        
        #x_contactthm_2 = data_contactthm2_2[:,0]
        #y_contactthm_2 = data_contactthm2_2[:,1]*(10**20)*2*pi
        #y_contactthm_2 = y_contactthm_2 - y_contactthm_2[-1]

        x_contactthm_12 = data_contactthm2_12[:,0]
        y_contactthm_12 = data_contactthm2_12[:,1]*(10**20)*2*pi
        y_contactthm_12 = y_contactthm_12 - y_contactthm_12[-1]
        
        #x_contactthm_22 = data_contactthm2_22[:,0]
        #y_contactthm_22 = data_contactthm2_22[:,1]*(10**20)*2*pi
        #y_contactthm_22 = y_contactthm_22 - y_contactthm_22[-1]

        
        #x_contactthm_3 = data_contactthm2_3[:,0]
        #y_contactthm_3 = data_contactthm2_3[:,1]*(10**20)*2*pi
        #y_contactthm_3 = y_contactthm_3 - y_contactthm_3[-1]
        
        #x_contactthm_4 = data_contactthm2_4[:,0]
        #y_contactthm_4 = data_contactthm2_4[:,1]*(10**20)*2*pi
        #y_contactthm_4 = y_contactthm_4 - y_contactthm_4[-1]
        
        #x_contactthm_5 = data_contactthm2_5[:,0]
        #y_contactthm_5 = data_contactthm2_5[:,1]*(10**20)*2*pi
        #y_contactthm_5 = y_contactthm_5 - y_contactthm_5[-1]
        
        #x_contactthm_6 = data_contactthm2_6[:,0]
        #y_contactthm_6 = data_contactthm2_6[:,1]*(10**20)*2*pi
        #y_contactthm_6 = y_contactthm_6 - y_contactthm_6[-1]
        
        #x_contactthm_7 = data_contactthm2_7[:,0]
        #y_contactthm_7 = data_contactthm2_7[:,1]*(10**20)*2*pi
        #y_contactthm_7 = y_contactthm_7 - y_contactthm_7[-1]

        #x_contactthm_8 = data_contactthm2_8[:,0]
        #y_contactthm_8 = data_contactthm2_8[:,1]*(10**20)*2*pi
        #y_contactthm_8 = y_contactthm_8 - y_contactthm_8[-1]

        #x_contactthm_9 = data_contactthm2_9[:,0]
        #y_contactthm_9 = data_contactthm2_9[:,1]*(10**20)*2*pi
        #y_contactthm_9 = y_contactthm_9 - y_contactthm_9[-1]

        #x_contactthm_10 = data_contactthm2_10[:,0]
        #y_contactthm_10 = data_contactthm2_10[:,1]*(10**20)*2*pi
        #y_contactthm_10 = y_contactthm_10 - y_contactthm_10[-1]

        #x_contactthm_11 = data_contactthm2_11[:,0]
        #y_contactthm_11 = data_contactthm2_11[:,1]*(10**20)*2*pi
        #y_contactthm_11 = y_contactthm_11 - y_contactthm_11[-1]

        #x_contactthm_12 = data_contactthm2_12[:,0]
        #y_contactthm_12 = data_contactthm2_12[:,1]*(10**20)*2*pi
        #y_contactthm_12 = y_contactthm_12 - y_contactthm_12[-1]

        #x_contactthm_13 = data_contactthm2_13[:,0]
        #y_contactthm_13 = data_contactthm2_13[:,1]*(10**20)*2*pi
        #y_contactthm_13 = y_contactthm_13 - y_contactthm_13[-1]

        
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
        
        
        #Print every nth one
        n = 18
        
        #print(y_contactthm_1[2:len(y_contactthm_1)].info())
        
        x_plot = [x_contactthm_1[x] for x in range(0, len(x_contactthm_1), n)]
        y_plot = [np.log(abs(y_contactthm_1[x])) for x in range(0, len(y_contactthm_1), n)]
        plt.plot(x_plot,y_plot,'bo', markersize=10, label='Potential C = 0')

        #x_plot = [x_contactthm_2[x] for x in range(0, len(x_contactthm_2), n)]
        #y_plot = [np.log(abs(y_contactthm_2[x])) for x in range(0, len(y_contactthm_2), n)]
        #plt.plot(x_plot,y_plot,'rx', markersize=10, label='Potential C = 0')
        
        #print(x_contactthm_12)

        x_plot = [x_contactthm_12[x] for x in range(0, len(x_contactthm_12), n)]
        y_plot = [np.log(abs(y_contactthm_12[x])) for x in range(0, len(y_contactthm_12), n)]
        plt.plot(x_plot,y_plot,'ks', markersize=10, label='Potential C = 60')

        #x_plot = [x_contactthm_22[x] for x in range(0, len(x_contactthm_22), n)]
        #y_plot = [np.log(abs(y_contactthm_22[x])) for x in range(0, len(y_contactthm_22), n)]
        #plt.plot(x_plot,y_plot,'g*', markersize=10, label='Potential C = 160')
        

        #plt.plot(x_contactthm_3[2:len(y_contactthm_3)],(y_contactthm_3[2:len(y_contactthm_3)]),'gs', label='Potential C = 0, n=0.000504, (Same total density as C=200)')
        #plt.plot(x_contactthm_4[2:len(y_contactthm_4)],(y_contactthm_4[2:len(y_contactthm_4)]),'k^', label='Potential C = 40')
        #plt.plot(x_contactthm_5[2:len(y_contactthm_5)],(y_contactthm_5[2:len(y_contactthm_5)]),'c<', label='Potential C = 60')
        #plt.plot(x_contactthm_6[2:len(y_contactthm_6)],(y_contactthm_6[2:len(y_contactthm_6)]),'m>', label='Potential C = 70')
        #plt.plot(x_contactthm_7[2:len(y_contactthm_7)],(y_contactthm_7[2:len(y_contactthm_7)]),'yH', label='Potential C = 77.02')
        #plt.plot(x_contactthm_8[2:len(y_contactthm_8)],(y_contactthm_8[2:len(y_contactthm_8)]),'ko', label='Potential C = 120')
        #plt.plot(x_contactthm_9[2:len(y_contactthm_9)],(y_contactthm_9[2:len(y_contactthm_9)]),'bx', label='Potential C = 150')
        #plt.plot(x_contactthm_10[2:len(y_contactthm_10)],(y_contactthm_10[2:len(y_contactthm_10)]),'ro', label='Potential C = 175')
        #plt.plot(x_contactthm_11[2:len(y_contactthm_11)],(y_contactthm_11[2:len(y_contactthm_11)]),'g^', label='Potential C = 200')
        #plt.plot(x_contactthm_12[2:len(y_contactthm_12)],(y_contactthm_12[2:len(y_contactthm_12)]),'y<', label='Potential C = 500')
        #plt.plot(x_contactthm_13[2:len(y_contactthm_13)],(y_contactthm_13[2:len(y_contactthm_13)]),'y<', label='Potential C = 500')

        
                
        #Set the axis labels.  Labelpad option controls the spacing between actual axis and axis label.  The r option tells python to interpret as a raw string literal.
        plt.xlabel(r"$h/\sigma$",labelpad=10)
        plt.ylabel(r"$\log|{2\pi\Delta{E_{s}}}|(N/m)$",labelpad=5)
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

        savefig("LongLogNeutralWalls-C60.pdf",bbox_inches='tight')
        #savefig("CompareLogPotentials-C4MIM+_TFSI-_model1-C0-C160_compare-negativewalls.pdf",bbox_inches='tight')



        
        #savefig("Potential-C4MIM+_TFSI-_model1density"+ density + "35.51-"+particle_wall_interaction+".pdf",bbox_inches='tight')
        
        # savefig("Potential-C4MIM+_TFSI-_model1density"+ density + "35.51-"+particle_wall_interaction+"at_maxima2.pdf",bbox_inches='tight')
        #savefig("Potential-C4MIM+_TFSI-_model1density"+ density + "35.51-"+particle_wall_interaction+"at_maxima_second_iteration_scheme.pdf",bbox_inches='tight')
        
        #Open a window and show the plot
        plt.show()
