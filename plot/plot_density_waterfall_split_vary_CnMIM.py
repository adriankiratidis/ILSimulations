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

charges = ["charge_minuswalls"] 

#list_of_particles=["MIM+_TFSI-_model2"] #nocharge and charge 4-a2 on C4
#list_of_particles=["MIM+_TFSI-_model1"] #no a2 on C4 for nocharge, for charge 5-a2
#list_of_particles=["MIM_BF4-"] #a2 on C4
#list_of_particles=["MIM_BF4-_minus_minus"] #a2 on C4

#list_of_particles=["MIM_BF4-_minus_minus","MIM+_TFSI-_model1_minus_minus"]
#list_of_particles=["MIM+_TFSI-_model2_density2_minus_minus"]

#list_of_particles=["MIM+_TFSI-_model1_density2_minus_minus", "MIM+_TFSI-_model2_minus_minus"]
#list_of_particles=["MIM+_TFSI-_model2_density2"]
list_of_particles=["MIM_BF4-_density2_minus_minus"]
#list_of_particles=["MIM+_TFSI-_model1_density2_minus_minus"]

for i in range(len(charges)):
    for j in range(len(list_of_particles)):

        data1_plus = np.loadtxt("../run_results/compare_epsilonLJ/"+charges[i]+"4-12/35.51-53.27_C2"+list_of_particles[j]+"/testing-a2-n_plus_separation10.00000charge1.txt")
        data1_minus = np.loadtxt("../run_results/compare_epsilonLJ/"+charges[i]+"4-12/35.51-53.27_C2"+list_of_particles[j]+"/testing-a2-n_minus_separation10.00000charge1.txt")

        data2_plus = np.loadtxt("../run_results/compare_epsilonLJ/"+charges[i]+"4-12/35.51-53.27_C4"+list_of_particles[j]+"/testing2-a2-n_plus_separation10.00000charge1.txt")
        data2_minus = np.loadtxt("../run_results/compare_epsilonLJ/"+charges[i]+"4-12/35.51-53.27_C4"+list_of_particles[j]+"/testing2-a2-n_minus_separation10.00000charge1.txt")
        
        data3_plus = np.loadtxt("../run_results/compare_epsilonLJ/"+charges[i]+"4-12/35.51-53.27_C6"+list_of_particles[j]+"/testing3-a2-n_plus_separation10.00000charge1.txt")
        data3_minus = np.loadtxt("../run_results/compare_epsilonLJ/"+charges[i]+"4-12/35.51-53.27_C6"+list_of_particles[j]+"/testing3-a2-n_minus_separation10.00000charge1.txt")
        
        data4_plus = np.loadtxt("../run_results/compare_epsilonLJ/"+charges[i]+"4-12/35.51-53.27_C8"+list_of_particles[j]+"/testing4-a2-n_plus_separation10.00000charge1.txt")
        data4_minus = np.loadtxt("../run_results/compare_epsilonLJ/"+charges[i]+"4-12/35.51-53.27_C8"+list_of_particles[j]+"/testing4-a2-n_minus_separation10.00000charge1.txt")
        
        data5_plus = np.loadtxt("../run_results/compare_epsilonLJ/"+charges[i]+"4-12/35.51-53.27_C10"+list_of_particles[j]+"/testing5-a2-n_plus_separation10.00000charge1.txt")
        data5_minus = np.loadtxt("../run_results/compare_epsilonLJ/"+charges[i]+"4-12/35.51-53.27_C10"+list_of_particles[j]+"/testing5-a2-n_minus_separation10.00000charge1.txt")

        
        x1_plus = list(data1_plus[:,0])
        y1_plus = list(data1_plus[:,1])
        #y1_plus[0] = 0.0
        
        x1_minus = list(data1_minus[:,0])
        y1_minus = data1_minus[:,1]+0.3
        #y1_plus[0] = 0.0
        
        x2_plus = list(data2_plus[:,0])
        y2_plus = data2_plus[:,1]
        #y2_plus[0] = 0.0
        
        x2_minus = list(data2_minus[:,0])
        y2_minus = data2_minus[:,1]+0.3
        #y2_minus[0] = 0.0
        
        x3_plus = list(data3_plus[:,0])
        y3_plus = data3_plus[:,1]
        #y3_plus[0] = 0.0
        
        x3_minus = list(data3_minus[:,0])
        y3_minus = data3_minus[:,1]+0.3
        #y3_minus[0] = 0.0
        
        x4_plus = list(data4_plus[:,0])
        y4_plus = data4_plus[:,1]
        #y4_plus[0] = 0.0
        
        x4_minus = list(data4_minus[:,0])
        y4_minus = data4_minus[:,1]+0.3
        #y4_minus[0] = 0.0
        
        x5_plus = list(data5_plus[:,0])
        y5_plus = data5_plus[:,1]
        #y5_plus[0] = 0.0
        
        x5_minus = list(data5_minus[:,0])
        y5_minus = data5_minus[:,1]+0.3
        #y5_minus[0] = 0.0
        
        
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        
        zeros=[0.0]*len(y1_plus)
        half=[0.5]*len(y2_plus)
        ones=[1.0]*len(y3_plus)
        onehalf=[1.5]*len(y4_plus)
        twos=[2.0]*len(y5_plus)

        
        xs = np.arange(0, 10, 0.4)
        
        #ax.scatter(x1_plus, zeros, y1_plus, c='k')
        #ax.scatter(x1_minus, zeros, y1_minus, c='k')
        
        #ax.scatter(x2_plus, half, y2_plus, c='k')
        #ax.scatter(x2_minus, half, y2_minus, c='k')
        
        #ax.scatter(x3_plus, ones, y3_plus, c='k')
        #ax.scatter(x3_minus, ones, y3_minus, c='k')
        
        #ax.scatter(x4_plus, onehalf, y4_plus, c='k')
        #ax.scatter(x4_minus, onehalf, y4_minus, c='k')
        
        #ax.scatter(x5_plus, twos, y5_plus, c='k')
        #ax.scatter(x5_minus, twos, y5_minus, c='k')
        
        #ax.scatter(x6_plus, twohalf, y6_plus, c='k')
        #ax.scatter(x6_minus, twohalf, y6_minus, c='k')
        
        verts = []
        verts2 = []
        verts3 = []
        verts4 = []
        verts5 = []
        
        pverts = []
        pverts2 = []
        pverts3 = []
        pverts4 = []
        pverts5 = []
        
        zs = [0.0]
        zs2 = [0.5]
        zs3 = [1.0]
        zs4 = [1.5]
        zs5 = [2.0]
        
        # y1_minus[len(y1_minus)/2] = 0.0
        # y2_minus[len(y1_minus)/2] = 0.0
        # y3_minus[len(y1_minus)/2] = 0.0
        # y4_minus[len(y1_minus)/2] = 0.0
        # y5_minus[len(y1_minus)/2] = 0.0
        # y6_minus[len(y1_minus)/2] = 0.0
        # verts.append(list(zip(x1_minus[len(x1_minus)/2:len(x1_minus)], y1_minus[len(x1_minus)/2:len(x1_minus)])))
        # verts2.append(list(zip(x2_minus[len(x1_minus)/2:len(x1_minus)], y2_minus[len(x1_minus)/2:len(x1_minus)])))
        # verts3.append(list(zip(x3_minus[len(x1_minus)/2:len(x1_minus)], y3_minus[len(x1_minus)/2:len(x1_minus)])))
        # verts4.append(list(zip(x4_minus[len(x1_minus)/2:len(x1_minus)], y4_minus[len(x1_minus)/2:len(x1_minus)])))
        # verts5.append(list(zip(x5_minus[len(x1_minus)/2:len(x1_minus)], y5_minus[len(x1_minus)/2:len(x1_minus)])))
        # verts6.append(list(zip(x6_minus[len(x1_minus)/2:len(x1_minus)], y6_minus[len(x1_minus)/2:len(x1_minus)])))
        
        verts.append(list(zip(x1_minus, y1_minus)))
        verts2.append(list(zip(x2_minus, y2_minus)))
        verts3.append(list(zip(x3_minus, y3_minus)))
        verts4.append(list(zip(x4_minus, y4_minus)))
        verts5.append(list(zip(x5_minus, y5_minus)))
        
        
        poly = PolyCollection(verts, facecolors=['g'])
        poly2 = PolyCollection(verts2, facecolors=['g'])
        poly3 = PolyCollection(verts3, facecolors=['g'])
        poly4 = PolyCollection(verts4, facecolors=['g'])
        poly5 = PolyCollection(verts5, facecolors=['g'])
        
        
        poly.set_alpha(0.25)
        poly2.set_alpha(0.4)
        poly3.set_alpha(0.55)
        poly4.set_alpha(0.7)
        poly5.set_alpha(0.85)
        
        poly.set_alpha(1.0)
        poly2.set_alpha(0.85)
        poly3.set_alpha(0.7)
        poly4.set_alpha(0.55)
        poly5.set_alpha(0.4)
        
        poly.set_alpha(1.0)
        poly2.set_alpha(0.65)
        poly3.set_alpha(0.45)
        poly4.set_alpha(0.3)
        poly5.set_alpha(0.2)

        
        pverts.append(list(zip(x1_plus, y1_plus)))
        pverts2.append(list(zip(x2_plus, y2_plus)))
        pverts3.append(list(zip(x3_plus, y3_plus)))
        pverts4.append(list(zip(x4_plus, y4_plus)))
        pverts5.append(list(zip(x5_plus, y5_plus)))
        
        
        ppoly = PolyCollection(pverts, facecolors=['b'])
        ppoly2 = PolyCollection(pverts2, facecolors=['b'])
        ppoly3 = PolyCollection(pverts3, facecolors=['b'])
        ppoly4 = PolyCollection(pverts4, facecolors=['b'])
        ppoly5 = PolyCollection(pverts5, facecolors=['b'])
        
        
        # ppoly.set_alpha(0.25)
        # ppoly2.set_alpha(0.4)
        # ppoly3.set_alpha(0.55)
        # ppoly4.set_alpha(0.7)
        # ppoly5.set_alpha(0.85)
        
        # ppoly.set_alpha(1.0)
        # ppoly2.set_alpha(0.85)
        # ppoly3.set_alpha(0.7)
        # ppoly4.set_alpha(0.55)
        # ppoly5.set_alpha(0.4)

        ppoly.set_alpha(1.0)
        ppoly2.set_alpha(0.65)
        ppoly3.set_alpha(0.45)
        ppoly4.set_alpha(0.3)
        ppoly5.set_alpha(0.2)
        
        
        
        ax.add_collection3d(ppoly, zs=zs, zdir='y')
        ax.add_collection3d(ppoly2, zs=zs2, zdir='y')
        ax.add_collection3d(ppoly3, zs=zs3, zdir='y')
        ax.add_collection3d(ppoly4, zs=zs4, zdir='y')
        ax.add_collection3d(ppoly5, zs=zs5, zdir='y')
        
        ax.add_collection3d(poly, zs=zs, zdir='y')
        ax.add_collection3d(poly2, zs=zs2, zdir='y')
        ax.add_collection3d(poly3, zs=zs3, zdir='y')
        ax.add_collection3d(poly4, zs=zs4, zdir='y')
        ax.add_collection3d(poly5, zs=zs5, zdir='y')

        
        labels = [item.get_text() for item in ax.get_zticklabels()] + ['u', 'u']
        #print labels
        labels[0] = '0.0'
        labels[1] = '0.1'
        labels[2] = '0.2'
        labels[3] = '0.0'
        labels[4] = '0.1'
        #labels[5] = '0.05'
        #labels[6] = '0.02'
        #labels[7] = '0.05'

        
        ax.set_zticklabels(labels, fontsize=11)
        
        labelsx = [item.get_text() for item in ax.get_xticklabels()]
        labelsx[0] = '0'
        labelsx[1] = '2'
        labelsx[2] = '4'
        labelsx[3] = '6'
        labelsx[4] = '8'
        labelsx[5] = '10'
        
        ax.set_xticklabels(labelsx, fontsize=11)
        
        labelsy = [item.get_text() for item in ax.get_yticklabels()]
        labelsy[0] = r'$C_{2}$'
        labelsy[1] = r'$C_{4}$'
        labelsy[2] = r'$C_{6}$'
        labelsy[3] = r'$C_{8}$'
        labelsy[4] = r'$C_{10}$'
    
        ax.set_yticklabels(labelsy, fontsize=11)
        
        
        
        
        ax.set_xlabel(r'$z/\sigma$', labelpad=1.5, fontsize=18)
        ax.set_xlim3d(0, 10)
        ax.tick_params(axis='x', which='major', pad=-2)
        
        ax.set_ylabel('cation tail', labelpad=4, fontsize=15.5)
        ax.set_ylim3d(-0.2, 2.2)
        ax.tick_params(axis='y', which='major', pad=-2)
        
        ax.set_zlabel(r'$n\sigma^3$',labelpad=8, fontsize=18)
        #ax.set_zlim3d(-0.0003, 0.0011)
        ax.set_zlim3d(-0.001, 0.41)
        
        ax.tick_params(axis='z',which='major', pad=4)
        
        #plt.ticklabel_format(style='sci', axis='z', scilimits=(0,0))
        
        savefig("DENSITY_PLOT_VARY_CnMIM_"+charges[i]+"_"+list_of_particles[j]+".pdf",bbox_inches='tight')
        
        plt.show()
        
