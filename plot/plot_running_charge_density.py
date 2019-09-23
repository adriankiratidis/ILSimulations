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
#e = 1.6021766E-19
e = 1
sigma = 2.4
wall_charge = (1.0/320.0) * e
bead_charge_plus = 0.2 * e
bead_charge_minus_BF4 =  0.2 * e
bead_charge_minus_TFSI1 =  e
bead_charge_minus_TFSI2 =  0.25 * e
n_discretisation_z = 50.0


#35.51-35.51_C4MIM+_TFSI-_model1_density2
data_plus_BF4 = np.loadtxt("../run_results/compare_epsilonLJ/charge4-12/35.51-88.78_C4MIM_BF4-_density2/testing6-a2-n_plus_separation10.00000charge1.txt")
data_minus_BF4 = np.loadtxt("../run_results/compare_epsilonLJ/charge4-12/35.51-88.78_C4MIM_BF4-_density2/testing6-a2-n_minus_separation10.00000charge1.txt")

data_plus_TFSI1 = np.loadtxt("../run_results/compare_epsilonLJ/charge4-12/35.51-88.78_C4MIM+_TFSI-_model1_density2/testing6-a2-n_plus_separation10.00000charge1.txt")
data_minus_TFSI1 = np.loadtxt("../run_results/compare_epsilonLJ/charge4-12/35.51-88.78_C4MIM+_TFSI-_model1_density2/testing6-a2-n_minus_separation10.00000charge1.txt")

data_plus_TFSI2 = np.loadtxt("../run_results/compare_epsilonLJ/charge4-12/35.51-88.78_C4MIM+_TFSI-_model2_density2/testing6-a2-n_plus_separation10.00000charge1.txt")
data_minus_TFSI2 = np.loadtxt("../run_results/compare_epsilonLJ/charge4-12/35.51-88.78_C4MIM+_TFSI-_model2_density2/testing6-a2-n_minus_separation10.00000charge1.txt")

x_plus_BF4 = data_plus_BF4[:,0]
y_plus_BF4 = data_plus_BF4[:,1]
y_minus_BF4 = data_minus_BF4[:,1]
accumulate_to_plot_BF4 = []

x_plus_TFSI1 = data_plus_TFSI1[:,0]
y_plus_TFSI1 = data_plus_TFSI1[:,1]
y_minus_TFSI1 = data_minus_TFSI1[:,1]
accumulate_to_plot_TFSI1 = []

x_plus_TFSI2 = data_plus_TFSI2[:,0]
y_plus_TFSI2 = data_plus_TFSI2[:,1]
y_minus_TFSI2 = data_minus_TFSI2[:,1]
accumulate_to_plot_TFSI2 = []

#x_minus = data_minus[:,0]
#y_minus = data_minus[:,1]

y_accumulate_BF4 = wall_charge
found_negative_number_BF4 = False
first_negative_number_occurs_at_BF4 = 0

y_accumulate_TFSI1 = wall_charge
found_negative_number_TFSI1 = False
first_negative_number_occurs_at_TFSI1 = 0

y_accumulate_TFSI2 = wall_charge
found_negative_number_TFSI2 = False
first_negative_number_occurs_at_TFSI2 = 0


for i in range(len(y_plus_BF4) - 1):
    y_accumulate_BF4 = y_accumulate_BF4 + (((y_plus_BF4[i] + y_plus_BF4[i + 1]) / 2.0) * (sigma / 50.0) * bead_charge_plus / sigma**3) - \
                   (((y_minus_BF4[i] + y_minus_BF4[i + 1]) / 2.0) * (sigma / 50.0) * bead_charge_minus_BF4 / sigma**3)

    if y_accumulate_BF4 < 0:
        if found_negative_number_BF4 == False:
            first_negative_number_occurs_at_BF4 = x_plus_BF4[i]
            found_negative_number_BF4 = True
    
    accumulate_to_plot_BF4.append(y_accumulate_BF4)


print "first negative number_BF4 = ", first_negative_number_occurs_at_BF4


# for i in range(len(y_plus_TFSI1) - 1):
#     y_accumulate_TFSI1 = y_accumulate_TFSI1 + (((y_plus_TFSI1[i] + y_plus_TFSI1[i + 1]) / 2.0) * (sigma / 50.0) * bead_charge_plus / sigma**3) - \
#                    (((y_minus_TFSI1[i] + y_minus_TFSI1[i + 1]) / 2.0) * (sigma / 50.0) * bead_charge_minus_TFSI1 / sigma**3)

#     if y_accumulate_TFSI1 < 0:
#         if found_negative_number_TFSI1 == False:
#             first_negative_number_occurs_at_TFSI1 = x_plus_TFSI1[i]
#             found_negative_number_TFSI1 = True
    
#     accumulate_to_plot_TFSI1.append(y_accumulate_TFSI1)

# print "first negative number_TFSI1 = ", first_negative_number_occurs_at_TFSI1


# for i in range(len(y_plus_TFSI2) - 1):
#     y_accumulate_TFSI2 = y_accumulate_TFSI2 + (((y_plus_TFSI2[i] + y_plus_TFSI2[i + 1]) / 2.0) * (sigma / 50.0) * bead_charge_plus / sigma**3) - \
#                    (((y_minus_TFSI2[i] + y_minus_TFSI2[i + 1]) / 2.0) * (sigma / 50.0) * bead_charge_minus_TFSI2 / sigma**3)

#     if y_accumulate_TFSI2 < 0:
#         if found_negative_number_TFSI2 == False:
#             first_negative_number_occurs_at_TFSI2 = x_plus_TFSI2[i]
#             found_negative_number_TFSI2 = True
    
#     accumulate_to_plot_TFSI2.append(y_accumulate_TFSI2)

# print "first negative number_TFSI2 = ", first_negative_number_occurs_at_TFSI2

#accumulate_to_plot_BF4 = [x - accumulate_to_plot_BF4[len(y_plus_BF4)/2] for x in accumulate_to_plot_BF4]

#plt.plot(x_plus_TFSI2[0:len(y_plus_TFSI2)/2], accumulate_to_plot_TFSI2[0:len(y_plus_TFSI2)/2], 'bo', label=r"$[\mathrm{C}_{4}\mathrm{MIM}^+][\mathrm{TFSI}^-]_2$")
#plt.plot(x_plus_TFSI1[0:len(y_plus_TFSI1)/2], accumulate_to_plot_TFSI1[0:len(y_plus_TFSI1)/2], 'rx', label=r"$[\mathrm{C}_{4}\mathrm{MIM}^+][\mathrm{TFSI}^-]_1$")
plt.plot(x_plus_BF4[0:len(y_plus_BF4)/2], accumulate_to_plot_BF4[0:len(y_plus_BF4)/2], 'gs', label=r"$[\mathrm{C}_{4}\mathrm{MIM}^+][\mathrm{BF}_4^-]$")




plt.axvline(x=0.5, color='k', linestyle='--')
plt.hlines(0.0, 0.0, 5.0, color='k', linestyle='-')
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

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
plt.xlabel(r"$z/\sigma$",labelpad=10, fontsize=20)
plt.ylabel(r"$\sigma_{T}\,\,(e / \AA^2)$",labelpad=5, fontsize=20)


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
plt.legend(loc='best',ncol=1, numpoints=1, frameon=False)

#Uncomment if title is required
#plt.title(r"Some Title")

#savefig("Single_density_plot_35.51-53.27-nocharge_C10MIM+_TFSI-_model2_density1_a2.pdf", bbox_inches='tight')
#savefig("C4MIM_TFSI_model2-35.51-35.51-charge_a2_TEST.pdf",bbox_inches='tight')
#savefig("C4MIM_TFSI1_running_charge_density_density2.pdf",bbox_inches='tight')

#Open a window and show the plot
plt.show()
