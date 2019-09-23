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

separation = 
#35.51-35.51_C4MIM+_TFSI-_model1_density2
data_plus = np.loadtxt("../run_results/compare_epsilonLJ/charge4-12/35.51-0_Heptamer_SingleSphere-0.005-45-0.003125-hs1.0CHARGE/testing10-a2-n_plus_separation"+separation+"0000charge1.txt")
data_neutral = np.loadtxt("../run_results/compare_epsilonLJ/charge4-12/35.51-0_Heptamer_SingleSphere-0.005-45-0.003125-hs1.0CHARGE/testing10-a2-n_neutral_separation"+separation+"0000charge1.txt")
data_minus = np.loadtxt("../run_results/compare_epsilonLJ/charge4-12/35.51-0_Heptamer_SingleSphere-0.005-45-0.003125-hs1.0CHARGE/testing10-a2-n_minus_separation"+separation+"0000charge1.txt")

data_plus2 = np.loadtxt("../run_results/compare_epsilonLJ/charge4-12/35.51-0_Heptamer_SingleSphere-0.005-45-0.003125-hs1.0CHARGE/testing10-a2-n_plus_separation"+separation+"0000charge1.txt")
data_neutral2 = np.loadtxt("../run_results/compare_epsilonLJ/charge4-12/35.51-0_Heptamer_SingleSphere-0.005-45-0.003125-hs1.0CHARGE/testing10-a2-n_neutral_separation"+separation+"0000charge1.txt")
data_minus2 = np.loadtxt("../run_results/compare_epsilonLJ/charge4-12/35.51-0_Heptamer_SingleSphere-0.005-45-0.003125-hs1.0CHARGE/testing10-a2-n_minus_separation"+separation+"0000charge1.txt")

data_plus3 = np.loadtxt("../run_results/compare_epsilonLJ/charge4-12/35.51-0_Heptamer_SingleSphere-0.005-45-0.003125-hs1.0CHARGE/testing10-a2-n_plus_separation"+separation+"0000charge1.txt")
data_neutral3 = np.loadtxt("../run_results/compare_epsilonLJ/charge4-12/35.51-0_Heptamer_SingleSphere-0.005-45-0.003125-hs1.0CHARGE/testing10-a2-n_neutral_separation"+separation+"0000charge1.txt")
data_minus3 = np.loadtxt("../run_results/compare_epsilonLJ/charge4-12/35.51-0_Heptamer_SingleSphere-0.005-45-0.003125-hs1.0CHARGE/testing10-a2-n_minus_separation"+separation+"0000charge1.txt")

data_plus4 = np.loadtxt("../run_results/compare_epsilonLJ/charge4-12/35.51-0_Heptamer_SingleSphere-0.005-45-0.003125-hs1.0CHARGE/testing10-a2-n_plus_separation"+separation+"0000charge1.txt")
data_neutral4 = np.loadtxt("../run_results/compare_epsilonLJ/charge4-12/35.51-0_Heptamer_SingleSphere-0.005-45-0.003125-hs1.0CHARGE/testing10-a2-n_neutral_separation"+separation+"0000charge1.txt")
data_minus4 = np.loadtxt("../run_results/compare_epsilonLJ/charge4-12/35.51-0_Heptamer_SingleSphere-0.005-45-0.003125-hs1.0CHARGE/testing10-a2-n_minus_separation"+separation+"0000charge1.txt")

data_plus5 = np.loadtxt("../run_results/compare_epsilonLJ/charge4-12/35.51-0_Heptamer_SingleSphere-0.005-45-0.003125-hs1.0CHARGE/testing10-a2-n_plus_separation"+separation+"0000charge1.txt")
data_neutral5 = np.loadtxt("../run_results/compare_epsilonLJ/charge4-12/35.51-0_Heptamer_SingleSphere-0.005-45-0.003125-hs1.0CHARGE/testing10-a2-n_neutral_separation"+separation+"0000charge1.txt")
data_minus5 = np.loadtxt("../run_results/compare_epsilonLJ/charge4-12/35.51-0_Heptamer_SingleSphere-0.005-45-0.003125-hs1.0CHARGE/testing10-a2-n_minus_separation"+separation+"0000charge1.txt")

data_plus6 = np.loadtxt("../run_results/compare_epsilonLJ/charge4-12/35.51-0_Heptamer_SingleSphere-0.005-45-0.003125-hs1.0CHARGE/testing10-a2-n_plus_separation"+separation+"0000charge1.txt")
data_neutral6 = np.loadtxt("../run_results/compare_epsilonLJ/charge4-12/35.51-0_Heptamer_SingleSphere-0.005-45-0.003125-hs1.0CHARGE/testing10-a2-n_neutral_separation"+separation+"0000charge1.txt")
data_minus6 = np.loadtxt("../run_results/compare_epsilonLJ/charge4-12/35.51-0_Heptamer_SingleSphere-0.005-45-0.003125-hs1.0CHARGE/testing10-a2-n_minus_separation"+separation+"0000charge1.txt")

data_plus7 = np.loadtxt("../run_results/compare_epsilonLJ/charge4-12/35.51-0_Heptamer_SingleSphere-0.005-45-0.003125-hs1.0CHARGE/testing10-a2-n_plus_separation"+separation+"0000charge1.txt")
data_neutral7 = np.loadtxt("../run_results/compare_epsilonLJ/charge4-12/35.51-0_Heptamer_SingleSphere-0.005-45-0.003125-hs1.0CHARGE/testing10-a2-n_neutral_separation"+separation+"0000charge1.txt")
data_minus7 = np.loadtxt("../run_results/compare_epsilonLJ/charge4-12/35.51-0_Heptamer_SingleSphere-0.005-45-0.003125-hs1.0CHARGE/testing10-a2-n_minus_separation"+separation+"0000charge1.txt")

data_plus8 = np.loadtxt("../run_results/compare_epsilonLJ/charge4-12/35.51-0_Heptamer_SingleSphere-0.005-45-0.003125-hs1.0CHARGE/testing10-a2-n_plus_separation"+separation+"0000charge1.txt")
data_neutral8 = np.loadtxt("../run_results/compare_epsilonLJ/charge4-12/35.51-0_Heptamer_SingleSphere-0.005-45-0.003125-hs1.0CHARGE/testing10-a2-n_neutral_separation"+separation+"0000charge1.txt")
data_minus8 = np.loadtxt("../run_results/compare_epsilonLJ/charge4-12/35.51-0_Heptamer_SingleSphere-0.005-45-0.003125-hs1.0CHARGE/testing10-a2-n_minus_separation"+separation+"0000charge1.txt")

data_plus9 = np.loadtxt("../run_results/compare_epsilonLJ/charge4-12/35.51-0_Heptamer_SingleSphere-0.005-45-0.003125-hs1.0CHARGE/testing10-a2-n_plus_separation"+separation+"0000charge1.txt")
data_neutral9 = np.loadtxt("../run_results/compare_epsilonLJ/charge4-12/35.51-0_Heptamer_SingleSphere-0.005-45-0.003125-hs1.0CHARGE/testing10-a2-n_neutral_separation"+separation+"0000charge1.txt")
data_minus9 = np.loadtxt("../run_results/compare_epsilonLJ/charge4-12/35.51-0_Heptamer_SingleSphere-0.005-45-0.003125-hs1.0CHARGE/testing10-a2-n_minus_separation"+separation+"0000charge1.txt")

data_plus10 = np.loadtxt("../run_results/compare_epsilonLJ/charge4-12/35.51-0_Heptamer_SingleSphere-0.005-45-0.003125-hs1.0CHARGE/testing10-a2-n_plus_separation"+separation+"0000charge1.txt")
data_neutral10 = np.loadtxt("../run_results/compare_epsilonLJ/charge4-12/35.51-0_Heptamer_SingleSphere-0.005-45-0.003125-hs1.0CHARGE/testing10-a2-n_neutral_separation"+separation+"0000charge1.txt")
data_minus10 = np.loadtxt("../run_results/compare_epsilonLJ/charge4-12/35.51-0_Heptamer_SingleSphere-0.005-45-0.003125-hs1.0CHARGE/testing10-a2-n_minus_separation"+separation+"0000charge1.txt")






x_plus = data_plus[:,0]
y_plus = data_plus[:,1]
y_minus = data_minus[:,1]
accumulate_to_plot = []



#x_minus = data_minus[:,0]
#y_minus = data_minus[:,1]

y_accumulate = 0.0

for i in range(len(y_plus) - 1):
    y_accumulate = y_accumulate + (((y_plus[i] + y_plus[i + 1]) / 2.0) * (x_plus[i + 1] - x_plus[i]))





plt.plot(1.0, y_accumulate, 'gs', label=r"$[\mathrm{C}_{4}\mathrm{MIM}^+][\mathrm{BF}_4^-]$")




plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))


#Set the axis labels.  Labelpad option controls the spacing between actual axis and axis label.  The r option tells python to interpret as a raw string literal.
plt.xlabel(r"$z/\sigma$",labelpad=10, fontsize=20)
plt.ylabel(r"$$",labelpad=5, fontsize=20)


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

#savefig("Single_density_plot_35.51-53.27-nocharge_C10MIM+_TFSI-_model2_density1_a2.pdf", bbox_inches='tight')
#savefig("C4MIM_TFSI_model2-35.51-35.51-charge_a2_TEST.pdf",bbox_inches='tight')
#savefig("C4MIM_TFSI1_running_charge_density_density2.pdf",bbox_inches='tight')

#Open a window and show the plot
plt.show()
