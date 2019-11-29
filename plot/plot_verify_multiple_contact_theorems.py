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
ELJ_wall=["0"]
ELJ_wall_index=[""]
particle_wall_epsilon="88.78"

for j in range(len(charges)):
   for i in range(len(ELJ_wall)):
      data_contactthm2_0 = np.loadtxt("../run_results/compare_epsilonLJ/"+charges[j]+"/35.51-" + particle_wall_epsilon + "_C4MIM_BF4--0.005-0--0.003125/testing"+ELJ_wall_index[i]+"-a2-normal-pressure-left-wall.txt")
      data_negativederiv2_0 = np.loadtxt("../run_results/compare_epsilonLJ/"+charges[j]+"/35.51-" + particle_wall_epsilon + "_C4MIM_BF4--0.005-0--0.003125/testing"+ELJ_wall_index[i]+"-a2-negative_deriv_of_potential.txt")

      data_contactthm2_1 = np.loadtxt("../run_results/compare_epsilonLJ/"+charges[j]+"/35.51-" + particle_wall_epsilon + "_C4MIM_BF4--0.005-0.2381--0.003125/testing2"+ELJ_wall_index[i]+"-a2-normal-pressure-left-wall.txt")
      data_negativederiv2_1 = np.loadtxt("../run_results/compare_epsilonLJ/"+charges[j]+"/35.51-" + particle_wall_epsilon + "_C4MIM_BF4--0.005-0.2381--0.003125/testing2"+ELJ_wall_index[i]+"-a2-negative_deriv_of_potential.txt")

      data_contactthm2_2 = np.loadtxt("../run_results/compare_epsilonLJ/"+charges[j]+"/35.51-" + particle_wall_epsilon + "_C4MIM_BF4--0.005-0.4762--0.003125/testing3"+ELJ_wall_index[i]+"-a2-normal-pressure-left-wall.txt")
      data_negativederiv2_2 = np.loadtxt("../run_results/compare_epsilonLJ/"+charges[j]+"/35.51-" + particle_wall_epsilon + "_C4MIM_BF4--0.005-0.4762--0.003125/testing3"+ELJ_wall_index[i]+"-a2-negative_deriv_of_potential.txt")

      data_contactthm2_3 = np.loadtxt("../run_results/compare_epsilonLJ/"+charges[j]+"/35.51-" + particle_wall_epsilon + "_C4MIM_BF4--0.005-0.9524--0.003125/testing4"+ELJ_wall_index[i]+"-a2-normal-pressure-left-wall.txt")
      data_negativederiv2_3 = np.loadtxt("../run_results/compare_epsilonLJ/"+charges[j]+"/35.51-" + particle_wall_epsilon + "_C4MIM_BF4--0.005-0.9524--0.003125/testing4"+ELJ_wall_index[i]+"-a2-negative_deriv_of_potential.txt")

      data_contactthm2_4 = np.loadtxt("../run_results/compare_epsilonLJ/"+charges[j]+"/35.51-" + particle_wall_epsilon + "_C4MIM_BF4--0.005-2.3810--0.003125/testing5"+ELJ_wall_index[i]+"-a2-normal-pressure-left-wall.txt")
      data_negativederiv2_4 = np.loadtxt("../run_results/compare_epsilonLJ/"+charges[j]+"/35.51-" + particle_wall_epsilon + "_C4MIM_BF4--0.005-2.3810--0.003125/testing5"+ELJ_wall_index[i]+"-a2-negative_deriv_of_potential.txt")

      data_contactthm2_5 = np.loadtxt("../run_results/compare_epsilonLJ/"+charges[j]+"/35.51-" + particle_wall_epsilon + "_C4MIM_BF4--0.005-4.7619--0.003125/testing6"+ELJ_wall_index[i]+"-a2-normal-pressure-left-wall.txt")
      data_negativederiv2_5 = np.loadtxt("../run_results/compare_epsilonLJ/"+charges[j]+"/35.51-" + particle_wall_epsilon + "_C4MIM_BF4--0.005-4.7619--0.003125/testing6"+ELJ_wall_index[i]+"-a2-negative_deriv_of_potential.txt")

      data_contactthm2_6 = np.loadtxt("../run_results/compare_epsilonLJ/"+charges[j]+"/35.51-" + particle_wall_epsilon + "_C4MIM_BF4--0.005-9.5238--0.003125/testing7"+ELJ_wall_index[i]+"-a2-normal-pressure-left-wall.txt")
      data_negativederiv2_6 = np.loadtxt("../run_results/compare_epsilonLJ/"+charges[j]+"/35.51-" + particle_wall_epsilon + "_C4MIM_BF4--0.005-9.5238--0.003125/testing7"+ELJ_wall_index[i]+"-a2-negative_deriv_of_potential.txt")

      # data_contactthm2_7 = np.loadtxt("../run_results/compare_epsilonLJ/"+charges[j]+"/35.51-" + particle_wall_epsilon + "_C4MIM_BF4--0.005-19.0476--0.003125/testing8"+ELJ_wall_index[i]+"-a2-normal-pressure-left-wall.txt")
      # data_negativederiv2_7 = np.loadtxt("../run_results/compare_epsilonLJ/"+charges[j]+"/35.51-" + particle_wall_epsilon + "_C4MIM_BF4--0.005-19.0476--0.003125/testing8"+ELJ_wall_index[i]+"-a2-negative_deriv_of_potential.txt")

      
   
      x_contactthm_0 = data_contactthm2_0[:,0]
      y_contactthm_0 = data_contactthm2_0[:,1]*(10**30)
      y_contactthm_0 = y_contactthm_0 - y_contactthm_0[-1]
      x_negativederiv_0 = data_negativederiv2_0[:,0]
      y_negativederiv_0 = data_negativederiv2_0[:,1]*(10**30)
      
      x_contactthm_1 = data_contactthm2_1[:,0]
      y_contactthm_1 = data_contactthm2_1[:,1]*(10**30)
      y_contactthm_1 = y_contactthm_1 - y_contactthm_1[-1]
      x_negativederiv_1 = data_negativederiv2_1[:,0]
      y_negativederiv_1 = data_negativederiv2_1[:,1]*(10**30)
      
      x_contactthm_2 = data_contactthm2_2[:,0]
      y_contactthm_2 = data_contactthm2_2[:,1]*(10**30)
      y_contactthm_2 = y_contactthm_2 - y_contactthm_2[-1]
      x_negativederiv_2 = data_negativederiv2_2[:,0]
      y_negativederiv_2 = data_negativederiv2_2[:,1]*(10**30)
      
      x_contactthm_3 = data_contactthm2_3[:,0]
      y_contactthm_3 = data_contactthm2_3[:,1]*(10**30)
      y_contactthm_3 = y_contactthm_3 - y_contactthm_3[-1]
      x_negativederiv_3 = data_negativederiv2_3[:,0]
      y_negativederiv_3 = data_negativederiv2_3[:,1]*(10**30)
      
      x_contactthm_4 = data_contactthm2_4[:,0]
      y_contactthm_4 = data_contactthm2_4[:,1]*(10**30)
      y_contactthm_4 = y_contactthm_4 - y_contactthm_4[-1]
      x_negativederiv_4 = data_negativederiv2_4[:,0]
      y_negativederiv_4 = data_negativederiv2_4[:,1]*(10**30)
      
      x_contactthm_5 = data_contactthm2_5[:,0]
      y_contactthm_5 = data_contactthm2_5[:,1]*(10**30)
      y_contactthm_5 = y_contactthm_5 - y_contactthm_5[-1]
      x_negativederiv_5 = data_negativederiv2_5[:,0]
      y_negativederiv_5 = data_negativederiv2_5[:,1]*(10**30)
      
      x_contactthm_6 = data_contactthm2_6[:,0]
      y_contactthm_6 = data_contactthm2_6[:,1]*(10**30)
      y_contactthm_6 = y_contactthm_6 - y_contactthm_6[-1]
      x_negativederiv_6 = data_negativederiv2_6[:,0]
      y_negativederiv_6 = data_negativederiv2_6[:,1]*(10**30)

      # x_contactthm_7 = data_contactthm2_7[:,0]
      # y_contactthm_7 = data_contactthm2_7[:,1]*(10**30)
      # y_contactthm_7 = y_contactthm_7 - y_contactthm_7[-1]
      # x_negativederiv_7 = data_negativederiv2_7[:,0]
      # y_negativederiv_7 = data_negativederiv2_7[:,1]*(10**30)
        

      
      #A description of the available plotting characters and colours 
      #to be used in place of 'rx' can be found here...
      #http://matplotlib.org/api/pyplot_api.html#matplotlib.pyplot.plot
      #To plot data with errorbars use 
      #plt.errorbar(x,y,err,fmt='bo',label="Some Label")
      #If you want to plot data but not show a key in the legend use
      #label='_nolegend_'
      #plt.plot(x_plus,y_plus,'rx',label=r'$n_{+}$')
       
      #y_negativederiv = y_negativederiv - y_negativederiv[-1]
      #y_contactthm = y_contactthm - y_contactthm[-1]
       
       
       
      #y_contactthm = y_contactthm * (y_negativederiv[0]/y_contactthm[0])
      #y_negativederiv = y_contactthm - y_contactthm[-1]
       
      #for i in range(len(y_contactthm)):
      #    #!if(abs(y_contactthm[i]) <= 0 ):
      #    #    print abs(y_contactthm[i])
      #    y_contactthm[i] = math.log(abs(y_contactthm[i]) if abs(y_contactthm[i]) > 0 else 1)
       
      #y_contactthm = math.log(y_contactthm)

              
      #plt.plot(x_contactthm_zero[2:],y_contactthm_zero[2:],'kx',label='Pressure from contact theorem 0')
      #plt.plot(x_negativederiv_zero[2:], y_negativederiv_zero[2:],'y^',label='Pressure from derivative of potential 0')

       
      plt.plot(x_contactthm_0[2:],y_contactthm_0[2:],'bo',label='Pressure from contact theorem C = 0')
      plt.plot(x_negativederiv_0[2:], y_negativederiv_0[2:],'gs',label='Pressure from derivative of potential C = 0')

      plt.plot(x_contactthm_1[2:],y_contactthm_1[2:],'rx',label='Pressure from contact theorem C = 0.5')
      plt.plot(x_negativederiv_1[2:], y_negativederiv_1[2:],'k^',label='Pressure from derivative of potential C = 0.5')

      plt.plot(x_contactthm_2[2:],y_contactthm_2[2:],'y<',label='Pressure from contact theorem C = 1')
      plt.plot(x_negativederiv_2[2:], y_negativederiv_2[2:],'c>',label='Pressure from derivative of potential C = 1')

      plt.plot(x_contactthm_3[2:],y_contactthm_3[2:],'mH',label='Pressure from contact theorem C = 2')
      plt.plot(x_negativederiv_3[2:], y_negativederiv_3[2:],'k*',label='Pressure from derivative of potential C = 2')

      plt.plot(x_contactthm_4[2:],y_contactthm_4[2:],'bo',label='Pressure from contact theorem C = 5')
      plt.plot(x_negativederiv_4[2:], y_negativederiv_4[2:],'gs',label='Pressure from derivative of potential C = 5')

      plt.plot(x_contactthm_5[2:],y_contactthm_5[2:],'r<',label='Pressure from contact theorem C = 10')
      plt.plot(x_negativederiv_5[2:], y_negativederiv_5[2:],'k>',label='Pressure from derivative of potential C = 10')

      plt.plot(x_contactthm_6[2:],y_contactthm_6[2:],'y^',label='Pressure from contact theorem C = 20')
      plt.plot(x_negativederiv_6[2:], y_negativederiv_6[2:],'kH',label='Pressure from derivative of potential C = 20')

      # plt.plot(x_contactthm_7[2:],y_contactthm_7[2:],'rx',label='Pressure from contact theorem C = 40')
      # plt.plot(x_negativederiv_7[2:], y_negativederiv_7[2:],'bs',label='Pressure from derivative of potential C = 40')

      
      #Set the axis labels.  Labelpad option controls the spacing between actual axis and axis label.  The r option tells python to interpret as a raw string literal.
      plt.xlabel(r"$h/\sigma$",labelpad=10, fontsize=20)
      plt.ylabel(r"${P_{int}}(N/m^2)$",labelpad=5, fontsize=20)
       
      plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0), fontsize=20)
    
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
      plt.legend(loc='best',ncol=1, numpoints=1, frameon=False, fontsize='xx-small')

      #Uncomment if title is required
      #plt.title(r"Some Title")
       
      #savefig("Contact_theorem_for_paper_35.51-53.17-charge4-12_C10MIM+_TFSI-_model2.pdf",bbox_inches='tight')
      #savefig("Contact_theorem_35.51-"+ELJ_wall[0]+"-"+charges[j]+IL[i]+"_density1_a2.pdf",bbox_inches='tight')
      #savefig("Contact_theorem_35.51-"+ELJ_wall[i]+"-"+charges[j]+"_C4MIM+_TFSI-_model1_density2_plus_halfplus.pdf",bbox_inches='tight')
      #savefig("Contact_theorem_for_paper_35.51-"+ELJ_wall[i]+"-"+charges[j]+"_C10MIM+_TFSI-_model2.pdf",bbox_inches='tight')
      #savefig("Contact_theorem_for_paper_35.51-53.17-nocharge4-12_C10MIM+_TFSI-_model2.pdf", bbox_inches='tight')
      #savefig("Contact_theorem_C4MIM_BF4_Negative_Wall_charge_C_Comparison.pdf",bbox_inches='tight')

      #Open a window and show the plot
      plt.show()
