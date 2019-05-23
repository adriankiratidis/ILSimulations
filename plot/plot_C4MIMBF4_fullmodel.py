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
charge = 1

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


charges=["unequalcharge"]
ELJ_wall=["0", "17.76", "35.51", "53.27", "71.02", "88.78"]
ELJ_wall_index=["", "2", "3", "4", "5", "6"]

for j in range(len(charges)):
    for i in range(len(ELJ_wall)):
        data_1 = np.loadtxt("../run_results/compare_epsilonLJ/"+charges[j]+"4-12/35.51-"+ELJ_wall[i]+"_C4MIM+_TFSI-_model1_plus_halfplus/testing"+ELJ_wall_index[i]+"-a2-n_plus_separation" +  str(separation) + ".00000charge" + str(charge) + ".txt")
        data_2 = np.loadtxt("../run_results/compare_epsilonLJ/"+charges[j]+"4-12/35.51-"+ELJ_wall[i]+"_C4MIM+_TFSI-_model1_plus_halfplus/testing"+ELJ_wall_index[i]+"-a2-n_neutral_separation" +  str(separation) + ".00000charge" + str(charge) + ".txt")
        data_3 = np.loadtxt("../run_results/compare_epsilonLJ/"+charges[j]+"4-12/35.51-"+ELJ_wall[i]+"_C4MIM+_TFSI-_model1_plus_halfplus/testing"+ELJ_wall_index[i]+"-a2-n_minus_separation" + str(separation) + ".00000charge" + str(charge) + ".txt")
        

# charges=["charge"]
# ELJ_wall=["53.27"]
# IL=["_C2MIM+_TFSI-_model2", "_C4MIM+_TFSI-_model2", "_C6MIM+_TFSI-_model2", "_C8MIM+_TFSI-_model2", "_C10MIM+_TFSI-_model2", "_C2MIM_BF4-", "_C4MIMBF4", "_C6MIM_BF4-", "_C8MIM_BF4-", "_C10MIM_BF4-", "_C2MIM+_TFSI-_model1", "_C4MIM+_TFSI-_model1", "_C6MIM+_TFSI-_model1", "_C8MIM+_TFSI-_model1", "_C10MIM+_TFSI-_model1"]
# IL_index=["", "4", "2", "3", "4", "", "4", "2", "3", "4", "", "5", "2", "3", "4"]


# for j in range(len(charges)):
#     for i in range(len(IL)):
#         data_1 = np.loadtxt("../run_results/compare_epsilonLJ/"+charges[j]+"4-12/35.51-"+ELJ_wall[0]+IL[i]+"/testing"+IL_index[i]+"-a2-n_plus_separation" +  str(separation) + ".00000charge" + str(charge) + ".txt")
#         data_2 = np.loadtxt("../run_results/compare_epsilonLJ/"+charges[j]+"4-12/35.51-"+ELJ_wall[0]+IL[i]+"/testing"+IL_index[i]+"-a2-n_neutral_separation" +  str(separation) + ".00000charge" + str(charge) + ".txt")
#         data_3 = np.loadtxt("../run_results/compare_epsilonLJ/"+charges[j]+"4-12/35.51-"+ELJ_wall[0]+IL[i]+"/testing"+IL_index[i]+"-a2-n_minus_separation" + str(separation) + ".00000charge" + str(charge) + ".txt")


        
        
        x_1 = data_1[:,0]
        y_1 = data_1[:,1]
        
        x_2 = data_2[:,0]
        y_2 = data_2[:,1]
        
        x_3 = data_3[:,0]
        y_3 = data_3[:,1]
        
        # x_4 = data_4[:,0]
        # y_4 = data_4[:,1]
        

        
        #A description of the available plotting characters and colours 
        #to be used in place of 'rx' can be found here...
        #http://matplotlib.org/api/pyplot_api.html#matplotlib.pyplot.plot
        #To plot data with errorbars use 
        #plt.errorbar(x,y,err,fmt='bo',label="Some Label")
        #If you want to plot data but not show a key in the legend use
        #label='_nolegend_'
        #plt.plot(x_0,y_0,'rx',label='converged profile')
        plt.plot(x_1+0.02,y_1,'bx',label=r'$n_{+}$')
        plt.plot(x_2+0.04,y_2,'go',label=r'$n_{0}$')
        plt.plot(x_3,y_3,'rs',label=r'$n_{-}$')
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
        
        
        #savefig("C4MIM_BF4-35.51-"+ELJ_wall[i]+"-"+charges[j]+"_a2_density2.pdf",bbox_inches='tight')
        #savefig("Density_35.51-"+ELJ_wall[0]+"-"+charges[j]+"_C4MIM+_TFSI-_model1_density2_plus_halfplus.pdf",bbox_inches='tight')
        
        savefig("Density_35.51-"+ELJ_wall[i]+"-"+charges[j]+"_C4MIM+_TFSI-_model1_density2_plus_halfplus.pdf",bbox_inches='tight')
        #savefig("Density_35.51-"+ELJ_wall[0]+"-"+charges[j]+IL[i]+"SEPARATION5.pdf",bbox_inches='tight')
        
        #Open a window and show the plot
        plt.show()
