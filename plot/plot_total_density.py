#!/usr/bin/python
import sys
import numpy as np
import matplotlib.pyplot as plt
from pylab import *
from matplotlib import rc
import math
import matplotlib
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

#Here we make use of the python plotting library matplotlib.
#Documentation can be found here http://matplotlib.org/

#Set the rc parameters to fix the font and allow tex to be used.
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

matplotlib.rcParams.update({'font.size': 17})

#Read in the data from "data.txt".  This will read in the whole file 
#and if necessary arrange the data into a rank 2 array.
#For example, if data.txt contained
#1   4.3   0.2
#2   8.6   0.1
#3   7.4   0.4
#4   8.9   0.4
#5   9.9   0.3
#then data[2,1] = 7.4 that is the 2+1 row and the 1+1 coloumn.

fig, ax = plt.subplots()

def find_total_integrated_density(y_plus, y_neutral, y_minus):
    y_accumulate = 0
    for i in range(len(y_plus) - 1):
        y_accumulate = y_accumulate + (((y_plus[i] + y_plus[i + 1]) / 2.0) * (x_plus[i+1] - x_plus[i]) )
        y_accumulate = y_accumulate + (((y_neutral[i] + y_neutral[i + 1]) / 2.0) * (x_neutral[i+1] - x_neutral[i]) )
        y_accumulate = y_accumulate + (((y_minus[i] + y_minus[i + 1]) / 2.0) * (x_minus[i+1] - x_minus[i]) )

    return y_accumulate

#Negative wall charge
possibleCs = ['0','10','20','40','60']
associated_indicies = ['','2','3','4','5']

#Positive wall charge
#possibleCs = ['0','10','20','30','40']
#associated_indicies = ['','3','5','7','9']

#Neutral wall charge
possibleCs = ['0','20','40','60','80','100','120','140','160','180','200']
associated_indicies = ['','2','3','4','5','6','7','8','9','10','']

#possibleCs = ['80']
#associated_indicies = ['5']

#possibleCs = ['0']
#associated_indicies = ['2']


separation = "8.0"
for iC_index, iC in enumerate(possibleCs):
    data_plus = np.loadtxt("../run_results/compare_epsilonLJ/CentreCentrePotential/minuswalls/35.51-0_Heptamer_SingleSphere-0.005-" + iC + "-0.0-hs1.0CHARGE/testing"+associated_indicies[iC_index]+"-a2-n_plus_separation"+separation+"0000charge1.txt")
    data_neutral = np.loadtxt("../run_results/compare_epsilonLJ/CentreCentrePotential/minuswalls/35.51-0_Heptamer_SingleSphere-0.005-" + iC + "-0.0-hs1.0CHARGE/testing"+associated_indicies[iC_index]+"-a2-n_neutral_separation"+separation+"0000charge1.txt")
    data_minus = np.loadtxt("../run_results/compare_epsilonLJ/CentreCentrePotential/minuswalls/35.51-0_Heptamer_SingleSphere-0.005-" + iC + "-0.0-hs1.0CHARGE/testing"+associated_indicies[iC_index]+"-a2-n_minus_separation"+separation+"0000charge1.txt")

    

    x_plus = data_plus[:,0]
    y_plus = data_plus[:,1]

    x_neutral = data_neutral[:,0]
    y_neutral = data_neutral[:,1]

    x_minus = data_minus[:,0]
    y_minus = data_minus[:,1]

    print("")
    print("iC = " + iC)
    print(find_total_integrated_density(y_plus, y_neutral, y_minus))
    print("")


    plt.plot(x_plus,y_plus,'c<',label=r'$n_{+}$')
    plt.plot(x_neutral,y_neutral,'k>',label=r'$n_{0}$')
    plt.plot(x_minus,y_minus,'y*',label=r'$n_{-}$')

plt.show()
    
#plt.axvline(x=0.5, color='k', linestyle='--')
#plt.hlines(0.0, 0.0, 5.0, color='b', linestyle='-')
#plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

#A description of the available plotting characters and colours 
#to be used in place of 'rx' can be found here...
#http://matplotlib.org/api/pyplot_api.html#matplotlib.pyplot.plot
#To plot data with errorbars use 
#plt.errorbar(x,y,err,fmt='bo',label="Some Label")
#If you want to plot data but not show a key in the legend use
#label='_nolegend_'
#plt.plot(x_0,y_0,'rx',label='converged profile')
#plt.plot(x_1+0.02,y_1,'bx',label=r'$n_{+}$')
#plt.plot(x_2+0.04,y_2,'go',label=r'$n_{0}$')
#plt.plot(x_3,y_3,'rs',label=r'$n_{-}$')
#plt.plot(x_4,y_4,'ch',label=r'$n_{s}$')
#plt.plot(x_5,y_5,'m*',label='iteration 5')
#plt.plot(x_6,y_6,'y<',label='iteration 6')
# plt.plot(x_7,y_7,'k>',label='iteration 7')
# plt.plot(x_8,y_8,'w^',label='iteration 8')
# plt.plot(x_9,y_9,'bD',label='iteration 9')

#Set the axis labels.  Labelpad option controls the spacing between actual axis and axis label.  The r option tells python to interpret as a raw string literal.
#plt.xlabel(r'$\epsilon_{LJ}^{pw}/\rho\sigma^{3}\epsilon_{LJ}^{pp}$',labelpad=10, fontsize=20)
#plt.ylabel(r"$D_{\sigma_{T}=0}/\sigma$",labelpad=5, fontsize=20)
#plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

# axins = zoomed_inset_axes(ax, 1.5, loc=1) # zoom-factor: 2.5, location: upper-left
# axins.plot(x_plus[0:len(y_plus)/2], accumulate_to_plot[0:len(y_plus)/2], 'rx')
# x1, x2, y1, y2 = 1, 3, -0.0005, 0.0005 # specify the limits
# axins.set_xlim(x1, x2) # apply the x-limits
# axins.set_ylim(y1, y2) # apply the y-limits
# plt.yticks(visible=False)
# plt.xticks(visible=False)
# mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="0.5")



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
#plt.legend(loc='best',ncol=1, numpoints=1, frameon=False)

#Uncomment if title is required
#plt.title(r"Some Title")

#savefig("Distance_vs_epsilon.pdf", bbox_inches='tight')
#savefig("Single_density_plot_35.51-53.27-nocharge_C10MIM+_TFSI-_model2_density1_a2.pdf", bbox_inches='tight')
#savefig("C4MIM_TFSI_model2-35.51-35.51-charge_a2_TEST.pdf",bbox_inches='tight')
#savefig("C2MIM+_TFSI-_35.51-35.51_model2_TEST_plus_minus.pdf",bbox_inches='tight')

#Open a window and show the plot
#plt.show()
