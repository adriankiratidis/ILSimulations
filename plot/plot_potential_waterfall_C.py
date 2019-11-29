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

#rcParams['ytick.major.pad']='80'

#35.51-53.27_C4MIM+_TFSI-_model2_density2
#C4MIM_BF4-_density2
charges = ["CentreCentrePotential/minuswalls"]

#list_of_particles=["MIM+_TFSI-_model2_density2_minus_minus"]
#list_of_particles=["MIM+_TFSI-_model2_density2"]
#list_of_particles=["MIM_BF4-"]
#list_of_particles=["MIM+_TFSI-_model1"]

list_of_particles=["Hexamer_SingleSphere"]
separation = "8.0"
r4r8="r8"
wall_charge="0.0"

for i in range(len(charges)):
    for j in range(len(list_of_particles)):

        data_contactthm_p1 = np.loadtxt("../run_results/compare_epsilonLJ/"+charges[i]+"/35.51-0_"+list_of_particles[j]+"-0.005-0-r8-"+wall_charge+"-hs1.0CHARGE/testing-a2-potential-per-unit-area.txt")
        data_contactthm_p2 = np.loadtxt("../run_results/compare_epsilonLJ/"+charges[i]+"/35.51-0_"+list_of_particles[j]+"-0.005-40-r8-"+wall_charge+"-hs1.0CHARGE/testing3-a2-potential-per-unit-area.txt")
        data_contactthm_p3 = np.loadtxt("../run_results/compare_epsilonLJ/"+charges[i]+"/35.51-0_"+list_of_particles[j]+"-0.005-80-r8-"+wall_charge+"-hs1.0CHARGE/testing5-a2-potential-per-unit-area.txt")
        data_contactthm_p4 = np.loadtxt("../run_results/compare_epsilonLJ/"+charges[i]+"/35.51-0_"+list_of_particles[j]+"-0.005-120-r8-"+wall_charge+"-hs1.0CHARGE/testing7-a2-potential-per-unit-area.txt")
        data_contactthm_p5 = np.loadtxt("../run_results/compare_epsilonLJ/"+charges[i]+"/35.51-0_"+list_of_particles[j]+"-0.005-180-r8-"+wall_charge+"-hs1.0CHARGE/testing10-a2-potential-per-unit-area.txt")


        data_contactthm_n1 = np.loadtxt("../run_results/compare_epsilonLJ/"+charges[i]+"/35.51-0_"+list_of_particles[j]+"-0.005-0-r4-"+wall_charge+"-hs1.0CHARGE/testing-a2-potential-per-unit-area.txt")
        data_contactthm_n2 = np.loadtxt("../run_results/compare_epsilonLJ/"+charges[i]+"/35.51-0_"+list_of_particles[j]+"-0.005-20-r4-"+wall_charge+"-hs1.0CHARGE/testing2-a2-potential-per-unit-area.txt")
        data_contactthm_n3 = np.loadtxt("../run_results/compare_epsilonLJ/"+charges[i]+"/35.51-0_"+list_of_particles[j]+"-0.005-40-r4-"+wall_charge+"-hs1.0CHARGE/testing3-a2-potential-per-unit-area.txt")
        data_contactthm_n4 = np.loadtxt("../run_results/compare_epsilonLJ/"+charges[i]+"/35.51-0_"+list_of_particles[j]+"-0.005-60-r4-"+wall_charge+"-hs1.0CHARGE/testing4-a2-potential-per-unit-area.txt")




        
        data_contactthm_m1 = np.loadtxt("../run_results/compare_epsilonLJ/"+charges[i]+"/35.51-0_C4MIM+_TFSI-_model1-0.005-0-"+r4r8+"-"+wall_charge+"-hs1.0CHARGE/testing-a2-potential-per-unit-area.txt")
        data_contactthm_m2 = np.loadtxt("../run_results/compare_epsilonLJ/"+charges[i]+"/35.51-0_C4MIM+_TFSI-_model1-0.005-40-"+r4r8+"-"+wall_charge+"-hs1.0CHARGE/testing3-a2-potential-per-unit-area.txt")
        data_contactthm_m3 = np.loadtxt("../run_results/compare_epsilonLJ/"+charges[i]+"/35.51-0_C4MIM+_TFSI-_model1-0.005-80-"+r4r8+"-"+wall_charge+"-hs1.0CHARGE/testing5-a2-potential-per-unit-area.txt")
        #data_contactthm_m4 = np.loadtxt("../run_results/compare_epsilonLJ/"+charges[i]+"/35.51-0_C4MIM+_TFSI-_model1-0.005-60-"+r4r8+"-"+wall_charge+"-hs1.0CHARGE/testing4-a2-potential-per-unit-area.txt")

        
        
        x_contactthm_p1 = list(data_contactthm_p1[:,0])
        y_contactthm_p1 = data_contactthm_p1[:,1]*(10**20)*2*pi*(10**2)
        y_contactthm_p1 = list(y_contactthm_p1 - y_contactthm_p1[-1])
        y_contactthm_p1[0] = 0.0
        
        x_contactthm_p2 = list(data_contactthm_p2[:,0])
        y_contactthm_p2 = data_contactthm_p2[:,1]*(10**20)*2*pi*(10**2)
        y_contactthm_p2 = list(y_contactthm_p2 - y_contactthm_p2[-1])
        y_contactthm_p2[0] = 0.0
        
        x_contactthm_p3 = list(data_contactthm_p3[:,0])
        y_contactthm_p3 = data_contactthm_p3[:,1]*(10**20)*2*pi*(10**2)
        y_contactthm_p3 = list(y_contactthm_p3 - y_contactthm_p3[-1])
        y_contactthm_p3[0] = 0.0

        x_contactthm_p4 = list(data_contactthm_p4[:,0])
        y_contactthm_p4 = data_contactthm_p4[:,1]*(10**20)*2*pi*(10**2)
        y_contactthm_p4 = list(y_contactthm_p4 - y_contactthm_p4[-1])
        y_contactthm_p4[0] = 0.0

        x_contactthm_p5 = list(data_contactthm_p5[:,0])
        y_contactthm_p5 = data_contactthm_p5[:,1]*(10**20)*2*pi*(10**2)
        y_contactthm_p5 = list(y_contactthm_p5 - y_contactthm_p5[-1])
        y_contactthm_p5[0] = 0.0



        

        x_contactthm_n1 = list(data_contactthm_n1[:,0])
        y_contactthm_n1 = data_contactthm_n1[:,1]*(10**20)*2*pi*(10**2)
        y_contactthm_n1 = list(y_contactthm_n1 - y_contactthm_n1[-1])
        y_contactthm_n1[0] = 0.0
        
        x_contactthm_n2 = list(data_contactthm_n2[:,0])
        y_contactthm_n2 = data_contactthm_n2[:,1]*(10**20)*2*pi*(10**2)
        y_contactthm_n2 = list(y_contactthm_n2 - y_contactthm_n2[-1])
        y_contactthm_n2[0] = 0.0
        
        x_contactthm_n3 = list(data_contactthm_n3[:,0])
        y_contactthm_n3 = data_contactthm_n3[:,1]*(10**20)*2*pi*(10**2)
        y_contactthm_n3 = list(y_contactthm_n3 - y_contactthm_n3[-1])
        y_contactthm_n3[0] = 0.0

        x_contactthm_n4 = list(data_contactthm_n4[:,0])
        y_contactthm_n4 = data_contactthm_n4[:,1]*(10**20)*2*pi*(10**2)
        y_contactthm_n4 = list(y_contactthm_n4 - y_contactthm_n4[-1])
        y_contactthm_n4[0] = 0.0


        x_contactthm_m1 = list(data_contactthm_m1[:,0])
        y_contactthm_m1 = data_contactthm_m1[:,1]*(10**20)*2*pi*(10**2)
        y_contactthm_m1 = list(y_contactthm_m1 - y_contactthm_m1[-1])
        y_contactthm_m1[0] = 0.0
        
        x_contactthm_m2 = list(data_contactthm_m2[:,0])
        y_contactthm_m2 = data_contactthm_m2[:,1]*(10**20)*2*pi*(10**2)
        y_contactthm_m2 = list(y_contactthm_m2 - y_contactthm_m2[-1])
        y_contactthm_m2[0] = 0.0
        
        x_contactthm_m3 = list(data_contactthm_m3[:,0])
        y_contactthm_m3 = data_contactthm_m3[:,1]*(10**20)*2*pi*(10**2)
        y_contactthm_m3 = list(y_contactthm_m3 - y_contactthm_m3[-1])
        y_contactthm_m3[0] = 0.0

        
        fig = plt.figure()
        ax = fig.gca(projection='3d')
    
        xs = np.arange(0, 10, 0.4)
        #Beginning of p
        verts_p1 = []
        verts_p2 = []
        verts_p3 = []
        verts_p4 = []
        verts_p5 = []

        zs_p1 = [2.0]
        zs_p2 = [4.0]
        zs_p3 = [6.0]
        zs_p4 = [8.0]
        zs_p5 = [10.0]
        
        verts_p1.append(list(zip(x_contactthm_p1, y_contactthm_p1)))
        verts_p2.append(list(zip(x_contactthm_p2, y_contactthm_p2)))
        verts_p3.append(list(zip(x_contactthm_p3, y_contactthm_p3)))
        verts_p4.append(list(zip(x_contactthm_p4, y_contactthm_p4)))
        verts_p5.append(list(zip(x_contactthm_p5, y_contactthm_p5)))
        
        poly_p1 = PolyCollection(verts_p1, facecolors=['b'])
        poly_p2 = PolyCollection(verts_p2, facecolors=['g'])
        poly_p3 = PolyCollection(verts_p3, facecolors=['r'])
        poly_p4 = PolyCollection(verts_p4, facecolors=['c'])
        poly_p5 = PolyCollection(verts_p5, facecolors=['k'])
        
        poly_p1.set_alpha(1.0)
        poly_p2.set_alpha(0.8)
        poly_p3.set_alpha(0.6)
        poly_p4.set_alpha(0.4)
        poly_p5.set_alpha(0.3)

        ax.add_collection3d(poly_p1, zs=zs_p1, zdir='y')
        ax.add_collection3d(poly_p2, zs=zs_p2, zdir='y')
        ax.add_collection3d(poly_p3, zs=zs_p3, zdir='y')
        ax.add_collection3d(poly_p4, zs=zs_p4, zdir='y')
        ax.add_collection3d(poly_p5, zs=zs_p5, zdir='y')

        #Beginning of n
        verts_n1 = []
        verts_n2 = []
        verts_n3 = []
        verts_n4 = []

        zs_n1 = [22.0]
        zs_n2 = [24.0]
        zs_n3 = [26.0]
        zs_n4 = [28.0]
        
        verts_n1.append(list(zip(x_contactthm_n1, y_contactthm_n1)))
        verts_n2.append(list(zip(x_contactthm_n2, y_contactthm_n2)))
        verts_n3.append(list(zip(x_contactthm_n3, y_contactthm_n3)))
        verts_n4.append(list(zip(x_contactthm_n4, y_contactthm_n4)))
        
        poly_n1 = PolyCollection(verts_n1, facecolors=['b'])
        poly_n2 = PolyCollection(verts_n2, facecolors=['g'])
        poly_n3 = PolyCollection(verts_n3, facecolors=['r'])
        poly_n4 = PolyCollection(verts_n4, facecolors=['c'])
        
        poly_n1.set_alpha(1.0)
        poly_n2.set_alpha(0.75)
        poly_n3.set_alpha(0.5)
        poly_n4.set_alpha(0.2)

        ax.add_collection3d(poly_n1, zs=zs_n1, zdir='y')
        ax.add_collection3d(poly_n2, zs=zs_n2, zdir='y')
        ax.add_collection3d(poly_n3, zs=zs_n3, zdir='y')
        ax.add_collection3d(poly_n4, zs=zs_n4, zdir='y')


        #Beginning of m
        verts_m1 = []
        verts_m2 = []
        verts_m3 = []

        zs_m1 = [25.0]
        zs_m2 = [27.0]
        zs_m3 = [29.0]
        
        verts_m1.append(list(zip(x_contactthm_m1, y_contactthm_m1)))
        verts_m2.append(list(zip(x_contactthm_m2, y_contactthm_m2)))
        verts_m3.append(list(zip(x_contactthm_m3, y_contactthm_m3)))
        
        poly_m1 = PolyCollection(verts_m1, facecolors=['b'])
        poly_m2 = PolyCollection(verts_m2, facecolors=['g'])
        poly_m3 = PolyCollection(verts_m3, facecolors=['r'])
        
        poly_m1.set_alpha(1.0)
        poly_m2.set_alpha(0.75)
        poly_m3.set_alpha(0.5)

        #ax.add_collection3d(poly_m1, zs=zs_m1, zdir='y')
        #ax.add_collection3d(poly_m2, zs=zs_m2, zdir='y')
        #ax.add_collection3d(poly_m3, zs=zs_m3, zdir='y')




        

        
        
        # labels = [item.get_text() for item in ax.get_zticklabels()]
        # #print labels
        # labels[0] = '0.0'
        # labels[1] = '0.1'
        # labels[2] = '0.2'
        # labels[3] = '0.0'
        # labels[4] = '0.1'
        # labels[5] = '0.2'
        # #labels[5] = '0.15'
        
        # ax.set_zticklabels(labels, fontsize=11)
        
        labelsx = [item.get_text() for item in ax.get_xticklabels()] + ['u', 'u', 'u']
        labelsx[0] = '4'
        labelsx[1] = '5'
        labelsx[2] = '6'
        labelsx[3] = '7'
        labelsx[4] = '8'
        labelsx[5] = '9'
        labelsx[6] = '10'
        #labelsx[7] = '11'
        #labelsx[8] = '12'
        
        ax.set_xticklabels(labelsx, fontsize=11)
        
        labelsy = [item.get_text() for item in ax.get_yticklabels()] + ['u', 'u', 'u']
        labelsy[0] = r'$C_{r_{8}}=0$'
        labelsy[1] = r'$C_{r_{8}}=40$'
        labelsy[2] = r'$C_{r_{8}}=80$'
        labelsy[3] = r'$C_{r_{8}}=120$'
        labelsy[4] = r'$C_{r_{8}}=180$'

        #for tick in ax.get_xaxis().get_major_ticks():
        #    tick.set_pad(80)
        #    tick.label1 = tick._get_text1()
    
        #ax.set_yticklabels(labelsy, fontsize=11)
  
        ax.set_xlabel(r'$h/\sigma$', labelpad=1.5, fontsize=17)
        ax.set_xlim3d(3, 10)
        ax.tick_params(axis='x', which='major', pad=-2)
        
        ax.set_ylabel('C value', labelpad=4, fontsize=15.5)
        ax.set_ylim3d(-0.2, 30.2)
        ax.tick_params(axis='y', which='minor', pad=-2)
        
        ax.set_zlabel(r'${2\pi\Delta{E_{s}}}(N/m)\left(\times 10^{-2}\right)$',labelpad=8, fontsize=17)
        ax.set_zlim3d(min(min(y_contactthm_p1), min(y_contactthm_p2), min(y_contactthm_p3), min(y_contactthm_p4), min(y_contactthm_p5)), \
                      max(max(y_contactthm_p1), max(y_contactthm_p2), max(y_contactthm_p3), max(y_contactthm_p4), max(y_contactthm_p5)))
        ax.tick_params(axis='z', which='major', pad=4)
        ax.view_init(22, -28 )
        plt.minorticks_on()
        plt.grid(b=True, which='minor', axis='y', color='#666666', linestyle='-')

        #savefig("Test3D.pdf",bbox_inches='tight')
        savefig("TESTING3D.pdf")
        
        plt.show()
        
