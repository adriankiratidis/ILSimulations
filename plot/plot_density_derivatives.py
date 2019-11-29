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

C_max = "45"
separation = "8.5"
epsilon_pw = "0"
density = "0.005"
particle = "Heptamer_SingleSphere"

data_plus_zero = np.loadtxt("/home/adriank/postdoc/2018/branch/ILSimulations/run_results/compare_epsilonLJ/CentreCentrePotential/minuswalls/35.51-" + epsilon_pw + "_" + particle + "-" + density + "-0-0.003125-hs1.0CHARGE/testing-a2-n_plus_separation"+ separation +"0000charge1.txt")
data_neutral_zero = np.loadtxt("/home/adriank/postdoc/2018/branch/ILSimulations/run_results/compare_epsilonLJ/CentreCentrePotential/minuswalls/35.51-" + epsilon_pw + "_" + particle + "-" + density + "-0-0.003125-hs1.0CHARGE/testing-a2-n_neutral_separation" + separation + "0000charge1.txt")
data_minus_zero = np.loadtxt("/home/adriank/postdoc/2018/branch/ILSimulations/run_results/compare_epsilonLJ/CentreCentrePotential/minuswalls/35.51-" + epsilon_pw + "_" + particle + "-" + density + "-0-0.003125-hs1.0CHARGE/testing-a2-n_minus_separation" + separation + "0000charge1.txt")


data_plus_separation = np.loadtxt("/home/adriank/postdoc/2018/branch/ILSimulations/run_results/compare_epsilonLJ/CentreCentrePotential/minuswalls/35.51-" + epsilon_pw + "_" + particle + "-" + density + "-" + C_max + "-0.003125-hs1.0CHARGE/testing10-a2-n_plus_separation" + separation + "0000charge1.txt")
data_neutral_separation = np.loadtxt("/home/adriank/postdoc/2018/branch/ILSimulations/run_results/compare_epsilonLJ/CentreCentrePotential/minuswalls/35.51-" + epsilon_pw + "_" + particle + "-" + density + "-" + C_max + "-0.003125-hs1.0CHARGE/testing10-a2-n_neutral_separation" + separation + "0000charge1.txt")
data_minus_separation = np.loadtxt("/home/adriank/postdoc/2018/branch/ILSimulations/run_results/compare_epsilonLJ/CentreCentrePotential/minuswalls/35.51-" + epsilon_pw + "_" + particle + "-" + density + "-" + C_max + "-0.003125-hs1.0CHARGE/testing10-a2-n_minus_separation" + separation + "0000charge1.txt")


x_plus_zero = data_plus_zero[:,0]
y_plus_zero = data_plus_zero[:,1]
y_plus_zero_derivative = list(data_plus_zero[:,1])

x_neutral_zero = data_neutral_zero[:,0]
y_neutral_zero = data_neutral_zero[:,1]
y_neutral_zero_derivative = list(data_neutral_zero[:,1])

x_minus_zero = data_minus_zero[:,0]
y_minus_zero = data_minus_zero[:,1]
y_minus_zero_derivative = list(data_minus_zero[:,1])

###########################################################################

x_plus_separation = data_plus_separation[:,0]
y_plus_separation = data_plus_separation[:,1]
y_plus_separation_derivative = list(data_plus_separation[:,1])

x_neutral_separation = data_neutral_separation[:,0]
y_neutral_separation = data_neutral_separation[:,1]
y_neutral_separation_derivative = list(data_neutral_separation[:,1])

x_minus_separation = data_minus_separation[:,0]
y_minus_separation = data_minus_separation[:,1]
y_minus_separation_derivative = list(data_minus_separation[:,1])


#for i in range(len(potential))[1:len(potential) - 1]:
#    differential_capacitance[i] = (sigma_left[i+1] - sigma_left[i-1])/(potential[i+1] - potential[i-1])

for i in range(len(x_plus_zero))[1:len(x_plus_zero) - 1]:
    y_plus_zero_derivative[i] = (y_plus_zero[i+1] - y_plus_zero[i-1])/(x_plus_zero[i+1] - x_plus_zero[i-1])
    y_neutral_zero_derivative[i] = (y_neutral_zero[i+1] - y_neutral_zero[i-1])/(x_neutral_zero[i+1] - x_neutral_zero[i-1])
    y_minus_zero_derivative[i] = (y_minus_zero[i+1] - y_minus_zero[i-1])/(x_minus_zero[i+1] - x_minus_zero[i-1])
    
    y_plus_separation_derivative[i] = (y_plus_separation[i+1] - y_plus_separation[i-1])/(x_plus_separation[i+1] - x_plus_separation[i-1])
    y_neutral_separation_derivative[i] = (y_neutral_separation[i+1] - y_neutral_separation[i-1])/(x_neutral_separation[i+1] - x_neutral_separation[i-1])
    y_minus_separation_derivative[i] = (y_minus_separation[i+1] - y_minus_separation[i-1])/(x_minus_separation[i+1] - x_minus_separation[i-1])

    




    
#plt.plot(potential[1:len(potential) - 1], differential_capacitance[1:len(differential_capacitance) - 1],'bo', label='DC')   

plt.plot(x_plus_zero, y_plus_zero_derivative,'bo', label='+; C=0')
plt.plot(x_neutral_zero, y_neutral_zero_derivative,'rx', label='0; C=0')
plt.plot(x_minus_zero, y_minus_zero_derivative,'g^', label='-; C=0')

plt.plot(x_plus_separation, y_plus_separation_derivative,'c^', label='+; C="max"')
plt.plot(x_neutral_separation, y_neutral_separation_derivative,'m<', label='0; C="max"')
plt.plot(x_minus_separation, y_minus_separation_derivative,'k>', label='-; C="max"')





#Set the axis labels.  Labelpad option controls the spacing between actual axis and axis label.  The r option tells python to interpret as a raw string literal.
plt.xlabel(r"$z/\sigma$",labelpad=10, fontsize=20)
plt.ylabel(r"$\frac{d(n\sigma^{3})}{dz}$",labelpad=5, fontsize=20)
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

savefig("DerivativeOfDensity-" + particle + "-epsilon_pw" + epsilon_pw + "-" + "density" + density + "-separation" + separation + ".pdf",bbox_inches='tight')
#savefig("DC.pdf",bbox_inches='tight')
#savefig("Potential_35.51-"+ELJ_wall[i]+"-"+charges[j]+"_C4MIM+_TFSI-_model1_density2_plus_halfplus.pdf",bbox_inches='tight')
#savefig("Potential_35.51-"+ELJ_wall[0]+"-"+charges[j]+IL[i]+"_density1_a2.pdf",bbox_inches='tight')

#Open a window and show the plot
plt.show()
