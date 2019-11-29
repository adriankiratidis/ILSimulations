#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
from pylab import *
from matplotlib import rc
import math
import matplotlib

#Here we make use of the python plotting library matplotlib.
#Documentation can be found here http://matplotlib.org/

#Set the rc parameters to fix the font and allow tex to be used.
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

matplotlib.rcParams.update({'font.size': 17})
#separation = 4
#charge = 1

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

separation_array = [3.02 + 0.02*x for x in range(50*12)]
separation_array = ['%.2f' % a for a in separation_array]
sum_1_density=[]
sum_2_density=[]
x=[]
for sep in separation_array:

    separation=sep
    r4r8="r4"
    wall_charge="-0.003125"
    n_positive_beads = 5.0
    n_neutral_beads = 19.0
    n_negative_beads = 1.0

    #
    #Hexamer_SingleSphere
    #C4MIM_BF4-
    #C4MIM+_TFSI-_model1

    #../run_results/compare_epsilonLJ/nocharge4-12/35.51-88.78_C4MIM+_TFSI-_model1-3
    data_1 = np.loadtxt("../run_results/compare_epsilonLJ/CentreCentrePotential/minuswalls/35.51-0_C4MIM+_TFSI-_model1-0.005-0-"+r4r8+"-"+wall_charge+"-hs1.0CHARGE-LONGRANGE30/testing-a2-n_plus_separation"+separation+"000charge1.txt")
    data_2 = np.loadtxt("../run_results/compare_epsilonLJ/CentreCentrePotential/minuswalls/35.51-0_C4MIM+_TFSI-_model1-0.005-0-"+r4r8+"-"+wall_charge+"-hs1.0CHARGE-LONGRANGE30/testing-a2-n_neutral_separation"+separation+"000charge1.txt")
    data_3 = np.loadtxt("../run_results/compare_epsilonLJ/CentreCentrePotential/minuswalls/35.51-0_C4MIM+_TFSI-_model1-0.005-0-"+r4r8+"-"+wall_charge+"-hs1.0CHARGE-LONGRANGE30/testing-a2-n_minus_separation"+separation+"000charge1.txt")


    #data_1_end = np.loadtxt("../run_results/compare_epsilonLJ/CentreCentrePotential/minuswalls/35.51-53.27_C4MIM+_TFSI-_model1-0.005-20-"+r4r8+"-"+wall_charge+"-hs1.0CHARGE/testing2-a2-n_plus_separation"+separation+"000charge1.txt")
    #data_2_end = np.loadtxt("../run_results/compare_epsilonLJ/CentreCentrePotential/minuswalls/35.51-53.27_C4MIM+_TFSI-_model1-0.005-20-"+r4r8+"-"+wall_charge+"-hs1.0CHARGE/testing2-a2-n_neutral_separation"+separation+"000charge1.txt")
    #data_3_end = np.loadtxt("../run_results/compare_epsilonLJ/CentreCentrePotential/minuswalls/35.51-53.27_C4MIM+_TFSI-_model1-0.005-20-"+r4r8+"-"+wall_charge+"-hs1.0CHARGE/testing2-a2-n_minus_separation"+separation+"000charge1.txt")


    data_1_end2 = np.loadtxt("../run_results/compare_epsilonLJ/CentreCentrePotential/minuswalls/35.51-0_C4MIM+_TFSI-_model1-0.005-20-"+r4r8+"-"+wall_charge+"-hs1.0CHARGE-LONGRANGE30/testing5-a2-n_plus_separation"+separation+"000charge1.txt")
    data_2_end2 = np.loadtxt("../run_results/compare_epsilonLJ/CentreCentrePotential/minuswalls/35.51-0_C4MIM+_TFSI-_model1-0.005-20-"+r4r8+"-"+wall_charge+"-hs1.0CHARGE-LONGRANGE30/testing5-a2-n_neutral_separation"+separation+"000charge1.txt")
    data_3_end2 = np.loadtxt("../run_results/compare_epsilonLJ/CentreCentrePotential/minuswalls/35.51-0_C4MIM+_TFSI-_model1-0.005-20-"+r4r8+"-"+wall_charge+"-hs1.0CHARGE-LONGRANGE30/testing5-a2-n_minus_separation"+separation+"000charge1.txt")


    x_1 = data_1[:,0]
    y_1t = data_1[:,1]
    y_1 = [e/n_positive_beads for e in y_1t]

    x_2 = data_2[:,0]
    y_2t = data_2[:,1]
    y_2 = [e/n_neutral_beads for e in y_2t]

    x_3 = data_3[:,0]
    y_3t = data_3[:,1]
    y_3 = [e/n_negative_beads for e in y_3t]

    sum_1_density.append((max(y_1) + max(y_2) + max(y_3))/3.0)
    
    # x_1_end = data_1_end[:,0]
    # y_1_end = data_1_end[:,1]
    
    # x_2_end = data_2_end[:,0]
    # y_2_end = data_2_end[:,1]

    # x_3_end = data_3_end[:,0]
    # y_3_end = data_3_end[:,1]


    x_1_end2 = data_1_end2[:,0]
    y_1_end2t = data_1_end2[:,1]
    y_1_end2 = [e/n_positive_beads for e in y_1_end2t]

    x_2_end2 = data_2_end2[:,0]
    y_2_end2t = data_2_end2[:,1]
    y_2_end2 = [e/n_neutral_beads for e in y_2_end2t]

    x_3_end2 = data_3_end2[:,0]
    y_3_end2t = data_3_end2[:,1]
    y_3_end2 = [e/n_negative_beads for e in y_3_end2t]

    sum_2_density.append((max(y_1_end2) + max(y_2_end2) + max(y_3_end2))/3.0)

    x.append(separation)
    
    #A description of the available plotting characters and colours 
    #to be used in place of 'rx' can be found here...
    #http://matplotlib.org/api/pyplot_api.html#matplotlib.pyplot.plot
    #To plot data with errorbars use 
    #plt.errorbar(x,y,err,fmt='bo',label="Some Label")
    #If you want to plot data but not show a key in the legend use
    #label='_nolegend_'
    #plt.plot(x_0,y_0,'rx',label='converged profile')

n=20

x_plot = [x[a] for a in range(0, len(x), n)]
y_1_plot = [sum_1_density[a] for a in range(0, len(sum_1_density), n)]
y_2_plot = [sum_2_density[a] for a in range(0, len(sum_2_density), n)]

plt.plot(x_plot, y_1_plot, 'bx', markersize = 9,label=r'C=0')
plt.plot(x_plot, y_2_plot, 'rs', markersize = 9,label=r'C=20')


#plt.plot(x_1_plot,y_1_plot,'bx', markersize = 9,label=r'$n_{+}$')
#plt.plot(x_2_plot,y_2_plot,'go', markersize = 9,label=r'$n_{0}$')
#plt.plot(x_3_plot,y_3_plot,'rs', markersize = 9,label=r'$n_{-}$')


#Set the axis labels.  Labelpad option controls the spacing between actual axis and axis label.  The r option tells python to interpret as a raw string literal.
plt.xlabel(r"$z/\sigma$",labelpad=10, fontsize=20)
plt.ylabel(r"$n\sigma^{3}$",labelpad=5, fontsize=20)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
#Set the axis limits if you want to specify a region of interest.
#The default is auto zoom.
#plt.ylim(-0.05,1.1)
#plt.xlim(0,14)

#plt.axvline(x=8,color='black')

#Change the position and label of the axis ticks.
#xtickloc = [ 0,3,6,9,12,15,18, 21, 24, 27 ]
#xtickval = [ '3','6','9','12', '15', '18', '21', '24', '27', '30']
#plt.xticks(xtickloc,xtickval)

#ytickloc = [ -1.5,-1.0,-0.5,0.0,0.5,1.0,1.5 ]
#ytickval = [ '','John','Paul','Adam', 'Bill', 'Dave']
#plt.yticks(ytickloc,ytickval)

#Uncomment if legend required.
plt.legend(loc='best',ncol=1, numpoints=1, frameon=False)

#Uncomment if title is required
#plt.title(r"Some Title")

savefig("CompareMaximumDensities-C60-negative-walls.pdf", bbox_inches='tight')

#Open a window and show the plot
plt.show()
