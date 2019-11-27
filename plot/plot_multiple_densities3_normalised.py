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

separation="9.98"
r4r8="r8"
wall_charge="0.0"
n_positive_beads = 1.0
n_neutral_beads = 5.0
n_negative_beads = 1.0

#
#Hexamer_SingleSphere
#C4MIM_BF4-
#C4MIM+_TFSI-_model1

#../run_results/compare_epsilonLJ/nocharge4-12/35.51-88.78_Hexamer_SingleSphere-3
data_1 = np.loadtxt("../run_results/compare_epsilonLJ/CentreCentrePotential/minuswalls/35.51-0_Hexamer_SingleSphere-0.005-0-"+r4r8+"-"+wall_charge+"-hs1.0CHARGE/testing-a2-n_plus_separation"+separation+"000charge1.txt")
data_2 = np.loadtxt("../run_results/compare_epsilonLJ/CentreCentrePotential/minuswalls/35.51-0_Hexamer_SingleSphere-0.005-0-"+r4r8+"-"+wall_charge+"-hs1.0CHARGE/testing-a2-n_neutral_separation"+separation+"000charge1.txt")
data_3 = np.loadtxt("../run_results/compare_epsilonLJ/CentreCentrePotential/minuswalls/35.51-0_Hexamer_SingleSphere-0.005-0-"+r4r8+"-"+wall_charge+"-hs1.0CHARGE/testing-a2-n_minus_separation"+separation+"000charge1.txt")


#data_1_end = np.loadtxt("../run_results/compare_epsilonLJ/CentreCentrePotential/minuswalls/35.51-53.27_Hexamer_SingleSphere-0.005-20-"+r4r8+"-"+wall_charge+"-hs1.0CHARGE/testing2-a2-n_plus_separation"+separation+"0000charge1.txt")
#data_2_end = np.loadtxt("../run_results/compare_epsilonLJ/CentreCentrePotential/minuswalls/35.51-53.27_Hexamer_SingleSphere-0.005-20-"+r4r8+"-"+wall_charge+"-hs1.0CHARGE/testing2-a2-n_neutral_separation"+separation+"0000charge1.txt")
#data_3_end = np.loadtxt("../run_results/compare_epsilonLJ/CentreCentrePotential/minuswalls/35.51-53.27_Hexamer_SingleSphere-0.005-20-"+r4r8+"-"+wall_charge+"-hs1.0CHARGE/testing2-a2-n_minus_separation"+separation+"0000charge1.txt")


data_1_end2 = np.loadtxt("../run_results/compare_epsilonLJ/CentreCentrePotential/minuswalls/35.51-0_Hexamer_SingleSphere-0.005-60-"+r4r8+"-"+wall_charge+"-hs1.0CHARGE/testing4-a2-n_plus_separation"+separation+"000charge1.txt")
data_2_end2 = np.loadtxt("../run_results/compare_epsilonLJ/CentreCentrePotential/minuswalls/35.51-0_Hexamer_SingleSphere-0.005-60-"+r4r8+"-"+wall_charge+"-hs1.0CHARGE/testing4-a2-n_neutral_separation"+separation+"000charge1.txt")
data_3_end2 = np.loadtxt("../run_results/compare_epsilonLJ/CentreCentrePotential/minuswalls/35.51-0_Hexamer_SingleSphere-0.005-60-"+r4r8+"-"+wall_charge+"-hs1.0CHARGE/testing4-a2-n_minus_separation"+separation+"000charge1.txt")


x_1 = data_1[:,0]
y_1t = data_1[:,1]
y_1 = [e/n_positive_beads for e in y_1t]

x_2 = data_2[:,0]
y_2t = data_2[:,1]
y_2 = [e/n_neutral_beads for e in y_2t]

x_3 = data_3[:,0]
y_3t = data_3[:,1]
y_3 = [e/n_negative_beads for e in y_3t]

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



#A description of the available plotting characters and colours 
#to be used in place of 'rx' can be found here...
#http://matplotlib.org/api/pyplot_api.html#matplotlib.pyplot.plot
#To plot data with errorbars use 
#plt.errorbar(x,y,err,fmt='bo',label="Some Label")
#If you want to plot data but not show a key in the legend use
#label='_nolegend_'
#plt.plot(x_0,y_0,'rx',label='converged profile')

n=10
x_1_plot = [x_1[x] for x in range(0, len(x_1), n)]
x_2_plot = [x_2[x] for x in range(0, len(x_2), n)]
x_3_plot = [x_3[x] for x in range(0, len(x_3), n)]

y_1_plot = [y_1[x] for x in range(0, len(y_1), n)]
y_2_plot = [y_2[x] for x in range(0, len(y_2), n)]
y_3_plot = [y_3[x] for x in range(0, len(y_3), n)]

plt.plot(x_1_plot,y_1_plot,'bx', markersize = 9,label=r'$n_{+}$')
plt.plot(x_2_plot,y_2_plot,'go', markersize = 9,label=r'$n_{0}$')
plt.plot(x_3_plot,y_3_plot,'rs', markersize = 9,label=r'$n_{-}$')

#plt.plot(x_1_end,y_1_end,'c<',label=r'$n_{+}$ (C=20)')
#plt.plot(x_2_end,y_2_end,'k>',label=r'$n_{0}$ (C=20)')
#plt.plot(x_3_end,y_3_end,'y*',label=r'$n_{-}$ (C=20)')


x_1_end2_plot = [x_1_end2[x] for x in range(0, len(x_1_end2), n)]
x_2_end2_plot = [x_2_end2[x] for x in range(0, len(x_2_end2), n)]
x_3_end2_plot = [x_3_end2[x] for x in range(0, len(x_3_end2), n)]

y_1_end2_plot = [y_1_end2[x] for x in range(0, len(y_1_end2), n)]
y_2_end2_plot = [y_2_end2[x] for x in range(0, len(y_2_end2), n)]
y_3_end2_plot = [y_3_end2[x] for x in range(0, len(y_3_end2), n)]

plt.plot(x_1_end2_plot,y_1_end2_plot,'y^',markersize = 9, label=r'$n_{+}$ (C=60)')
plt.plot(x_2_end2_plot,y_2_end2_plot,'k*',markersize = 9, label=r'$n_{0}$ (C=60)')
plt.plot(x_3_end2_plot,y_3_end2_plot,'m<',markersize = 9, label=r'$n_{-}$ (C=60)')


#Set the axis labels.  Labelpad option controls the spacing between actual axis and axis label.  The r option tells python to interpret as a raw string literal.
plt.xlabel(r"$z/\sigma$",labelpad=10, fontsize=20)
plt.ylabel(r"$n\sigma^{3}$",labelpad=5, fontsize=20)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.tick_params(direction="in", top=True, right=True)
#Set the axis limits if you want to specify a region of interest.
#The default is auto zoom.
#plt.ylim(-0.05,1.1)
#plt.ylim(0,0.038)
#plt.ylim(0,0.014)

#plt.axvline(x=8,color='black')

#Change the position and label of the axis ticks.
#xtickloc = [ -1.5,-1.0,-0.5,0.0,0.5,1.0,1.5 ]
#xtickval = [ '','John','Paul','Adam', 'Bill', 'Dave']
#plt.xticks(xtickloc,xtickval)

#ytickloc = [ -1.5,-1.0,-0.5,0.0,0.5,1.0,1.5 ]
#ytickval = [ '','John','Paul','Adam', 'Bill', 'Dave']
#plt.yticks(ytickloc,ytickval)

#Uncomment if legend required.
plt.legend(loc='best', ncol=2, handletextpad=0.1, labelspacing=0.4, columnspacing=0.5, numpoints=1, frameon=False)
#plt.legend(loc='center', bbox_to_anchor=(0.5, 0.4), ncol=2, handletextpad=0.1, labelspacing=0.4, columnspacing=0.5, numpoints=1, frameon=False)

#Uncomment if title is required
#plt.title(r"Some Title")

#savefig("CompareDifferentDensities.pdf", bbox_inches='tight)'
savefig("NormalisedDensities-Hexamer_SingleSphere-C60-"+r4r8+"-wallcharge"+wall_charge+"sep"+separation+".pdf", bbox_inches='tight')

#savefig("NormalisedDensities-Hexamer_SingleSphere--density0.005-epsilon35.51-0-sep"+separation+"_"+r4r8+"-"+ wall_charge+".pdf", bbox_inches='tight')
#savefig("Densities-Hexamer_SingleSphere--density0.005-epsilon35.51-0-sep8.0neutral_walls.pdf", bbox_inches='tight')
#savefig("Densities-Hexamer_SingleSphere--density0.005-epsilon35.51-0-sep8.0neutral_walls.pdf", bbox_inches='tight')

#savefig("Densities_Hexamer_SingleSphere_Negative_Wall_Charge-C20-sep3.5-pw88.78.pdf", bbox_inches='tight')
#savefig("Single_density_plot_35.51-0-nocharge_C10MIM+_TFSI-_model2_density1_a2.pdf", bbox_inches='tight')
#savefig("C4MIM_TFSI_model2-35.51-0-charge_a2_TEST.pdf",bbox_inches='tight')
#savefig("C2MIM+_TFSI-_35.51-0_model2_TEST_plus_minus.pdf",bbox_inches='tight')

#Open a window and show the plot
plt.show()
