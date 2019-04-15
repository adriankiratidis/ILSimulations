#!/usr/bin/python
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.collections import PolyCollection
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib import colors as mcolors
import numpy as np
import math
from pylab import *

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
#35.51-53.27_C4MIM+_TFSI-_model2_density2
#C4MIM_BF4-_density2
charges = ["charge4-12"]
#list_of_particles=["_C4MIM+_TFSI-_model2"]
#list_of_particles=["_C4MIM+_TFSI-_model2_density2"]
#list_of_particles=["_C4MIM_BF4-_density2", "_C4MIMBF4"]
list_of_particles=["_C4MIM_BF4-_density2"]

for i in range(len(charges)):
    for j in range(len(list_of_particles)):
        #print charges[i], list_of_particles[j]
        
        data_contactthm1 = np.loadtxt("../run_results/compare_epsilonLJ/"+charges[i]+"/35.51-0"+list_of_particles[j]+"/testing-a2-potential-per-unit-area.txt")
        data_contactthm2 = np.loadtxt("../run_results/compare_epsilonLJ/"+charges[i]+"/35.51-17.76"+list_of_particles[j]+"/testing2-a2-potential-per-unit-area.txt")
        data_contactthm3 = np.loadtxt("../run_results/compare_epsilonLJ/"+charges[i]+"/35.51-35.51"+list_of_particles[j]+"/testing3-a2-potential-per-unit-area.txt")
        data_contactthm4 = np.loadtxt("../run_results/compare_epsilonLJ/"+charges[i]+"/35.51-53.27"+list_of_particles[j]+"/testing4-a2-potential-per-unit-area.txt")
        data_contactthm5 = np.loadtxt("../run_results/compare_epsilonLJ/"+charges[i]+"/35.51-71.02"+list_of_particles[j]+"/testing5-a2-potential-per-unit-area.txt")
        data_contactthm6 = np.loadtxt("../run_results/compare_epsilonLJ/"+charges[i]+"/35.51-88.78"+list_of_particles[j]+"/testing6-a2-potential-per-unit-area.txt")

# for i in range(len(charges)):
#     for j in range(len(list_of_particles)):

#         data_contactthm1 = np.loadtxt("../run_results/compare_epsilonLJ/"+charges[i]+"4-12/35.51-0_"+list_of_particles[j]+"/testing-a2-potential-per-unit-area.txt")
#         data_contactthm2 = np.loadtxt("../run_results/compare_epsilonLJ/"+charges[i]+"4-12/35.51-17.76_"+list_of_particles[j]+"/testing2-a2-potential-per-unit-area.txt")
#         data_contactthm3 = np.loadtxt("../run_results/compare_epsilonLJ/"+charges[i]+"4-12/35.51-35.51_"+list_of_particles[j]+"/testing3-a2-potential-per-unit-area.txt")
#         data_contactthm4 = np.loadtxt("../run_results/compare_epsilonLJ/"+charges[i]+"4-12/35.51-53.27_"+list_of_particles[j]+"/testing4-a2-potential-per-unit-area.txt")
#         data_contactthm5 = np.loadtxt("../run_results/compare_epsilonLJ/"+charges[i]+"4-12/35.51-71.02_"+list_of_particles[j]+"/testing5-a2-potential-per-unit-area.txt")
#         data_contactthm6 = np.loadtxt("../run_results/compare_epsilonLJ/"+charges[i]+"4-12/35.51-88.78_"+list_of_particles[j]+"/testing6-a2-potential-per-unit-area.txt")


        
        x_contactthm = list(data_contactthm1[:,0])
        y_contactthm = data_contactthm1[:,1]*(10**20)*2*pi*(10**2)
        y_contactthm = list(y_contactthm - y_contactthm[-1])
        y_contactthm[0] = 0.0
        
        x_contactthm2 = list(data_contactthm2[:,0])
        y_contactthm2 = data_contactthm2[:,1]*(10**20)*2*pi*(10**2)
        y_contactthm2 = list(y_contactthm2 - y_contactthm2[-1])
        y_contactthm2[0] = 0.0
        
        x_contactthm3 = list(data_contactthm3[:,0])
        y_contactthm3 = data_contactthm3[:,1]*(10**20)*2*pi*(10**2)
        y_contactthm3 = list(y_contactthm3 - y_contactthm3[-1])
        y_contactthm3[0] = 0.0
        
        x_contactthm4 = list(data_contactthm4[:,0])
        y_contactthm4 = data_contactthm4[:,1]*(10**20)*2*pi*(10**2)
        y_contactthm4 = list(y_contactthm4 - y_contactthm4[-1])
        y_contactthm4[0] = 0.0
        
        x_contactthm5 = list(data_contactthm5[:,0])
        y_contactthm5 = data_contactthm5[:,1]*(10**20)*2*pi*(10**2)
        y_contactthm5 = list(y_contactthm5 - y_contactthm5[-1])
        y_contactthm5[0] = 0.0
        
        x_contactthm6 = list(data_contactthm6[:,0])
        y_contactthm6 = data_contactthm6[:,1]*(10**20)*2*pi*(10**2)
        y_contactthm6 = list(y_contactthm6 - y_contactthm6[-1])
        y_contactthm6[0] = 0.0
        
        
        #y_contactthm2 = y_contactthm.copy()
        #for i in range(len(y_contactthm2)):
        #    y_contactthm2[i] = 0.0 
        
        
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        
        
        #def cc(arg):
        #    return arg#mcolors.to_rgba(arg, alpha=0.6)
        
        #xs = np.arange(0, 10, 0.4)
        #verts = []
        #zs = [0.0]
        #verts.append(list(zip(x_contactthm, y_contactthm)))
        
        #for z in zs:
        #    ys = np.random.rand(len(xs))
        #    ys[0], ys[-1] = 0, 0
        #    verts.append(list(zip(xs, ys)))
        
        
        
        #print type(verts), len(verts)#, type(verts[1])
        #print verts
        #verts=
        
        xs = np.arange(0, 10, 0.4)
        verts = []
        verts2 = []
        verts3 = []
        verts4 = []
        verts5 = []
        verts6 = []
        
        zs = [0.0]
        zs2 = [0.5]
        zs3 = [1.0]
        zs4 = [1.5]
        zs5 = [2.0]
        zs6 = [2.5]
        

        #for z in zs:
        #    ys = np.random.rand(len(xs))
        #    ys[0], ys[-1] = 0, 0
        #    verts.append(list(zip(xs, ys)))
        verts.append(list(zip(x_contactthm, y_contactthm)))
        verts2.append(list(zip(x_contactthm2, y_contactthm2)))
        verts3.append(list(zip(x_contactthm3, y_contactthm3)))
        verts4.append(list(zip(x_contactthm4, y_contactthm4)))
        verts5.append(list(zip(x_contactthm5, y_contactthm5)))
        verts6.append(list(zip(x_contactthm6, y_contactthm6)))
        
        #ax.plot(x_contactthm, y_contactthm, zs=0, zdir='y', color='b')
        #ax.plot(x_contactthm, y_contactthm2, zs=0, zdir='y', color='b')
        #ax.fill_between(x_contactthm, y_contactthm, y_contactthm2, facecolor='blue', alpha=0.5)
        
        #poly = PolyCollection(verts, facecolors=['r', 'g', 'y', 'c', 'm', 'b'])
        poly = PolyCollection(verts, facecolors=['r'])
        poly2 = PolyCollection(verts2, facecolors=['r'])
        poly3 = PolyCollection(verts3, facecolors=['r'])
        poly4 = PolyCollection(verts4, facecolors=['r'])
        poly5 = PolyCollection(verts5, facecolors=['r'])
        poly6 = PolyCollection(verts6, facecolors=['r'])
        
        
        poly.set_alpha(1.0)
        poly2.set_alpha(0.75)
        poly3.set_alpha(0.55)
        poly4.set_alpha(0.4)
        poly5.set_alpha(0.25)
        poly6.set_alpha(0.15)
        
        ax.add_collection3d(poly, zs=zs, zdir='y')
        ax.add_collection3d(poly2, zs=zs2, zdir='y')
        ax.add_collection3d(poly3, zs=zs3, zdir='y')
        ax.add_collection3d(poly4, zs=zs4, zdir='y')
        ax.add_collection3d(poly5, zs=zs5, zdir='y')
        ax.add_collection3d(poly6, zs=zs6, zdir='y')
        
        
        ax.set_xlabel(r'$h/\sigma$', labelpad=1.5, fontsize=17)
        ax.set_xlim3d(4, 12)
        ax.tick_params(axis='x', which='major', pad=-2)
        
        ax.set_ylabel(r'$\epsilon_{LJ}^{pw}/\rho\sigma^{3}\epsilon_{LJ}^{pp}$', labelpad=4, fontsize=17)
        ax.set_ylim3d(-0.2, 2.7)
        ax.tick_params(axis='y', which='major', pad=-2)
        
        ax.set_zlabel(r'${2\pi\Delta{E_{s}}}(N/m)\left(\times 10^{-2}\right)$',labelpad=8, fontsize=17)
        #ax.set_zlim3d(-0.0003, 0.0011)
        ax.set_zlim3d(min(min(y_contactthm), min(y_contactthm2), min(y_contactthm3), min(y_contactthm4), min(y_contactthm5), min(y_contactthm6)), max(max(y_contactthm), max(y_contactthm2), max(y_contactthm3), max(y_contactthm4),max(y_contactthm5), max(y_contactthm6)))
        ax.tick_params(axis='z', which='major', pad=4)
        #ax.autoscale_view(tight=None)
        
        #plt.ticklabel_format(style='sci', axis='z', scilimits=(0,0))
        
        savefig("paper_revision_plots/WATERFALL_PLOT_POTENTIAL_"+charges[i]+"_"+list_of_particles[j]+".pdf",bbox_inches='tight')
        
        plt.show()
        
