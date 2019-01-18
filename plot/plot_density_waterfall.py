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

data1_plus = np.loadtxt("../run_results/compare_epsilonLJ/charge4-12/35.51-0_C4MIM+_TFSI-_model1_density2/testing-a2-n_plus_separation10.00000charge1.txt")
data1_minus = np.loadtxt("../run_results/compare_epsilonLJ/charge4-12/35.51-0_C4MIM+_TFSI-_model1_density2/testing-a2-n_minus_separation10.00000charge1.txt")

data2_plus = np.loadtxt("../run_results/compare_epsilonLJ/charge4-12/35.51-17.76_C4MIM+_TFSI-_model1_density2/testing2-a2-n_plus_separation10.00000charge1.txt")
data2_minus = np.loadtxt("../run_results/compare_epsilonLJ/charge4-12/35.51-17.76_C4MIM+_TFSI-_model1_density2/testing2-a2-n_minus_separation10.00000charge1.txt")

data3_plus = np.loadtxt("../run_results/compare_epsilonLJ/charge4-12/35.51-35.51_C4MIM+_TFSI-_model1_density2/testing3-a2-n_plus_separation10.00000charge1.txt")
data3_minus = np.loadtxt("../run_results/compare_epsilonLJ/charge4-12/35.51-35.51_C4MIM+_TFSI-_model1_density2/testing3-a2-n_minus_separation10.00000charge1.txt")

data4_plus = np.loadtxt("../run_results/compare_epsilonLJ/charge4-12/35.51-53.27_C4MIM+_TFSI-_model1_density2/testing4-a2-n_plus_separation10.00000charge1.txt")
data4_minus = np.loadtxt("../run_results/compare_epsilonLJ/charge4-12/35.51-53.27_C4MIM+_TFSI-_model1_density2/testing4-a2-n_minus_separation10.00000charge1.txt")

data5_plus = np.loadtxt("../run_results/compare_epsilonLJ/charge4-12/35.51-71.02_C4MIM+_TFSI-_model1_density2/testing5-a2-n_plus_separation10.00000charge1.txt")
data5_minus = np.loadtxt("../run_results/compare_epsilonLJ/charge4-12/35.51-71.02_C4MIM+_TFSI-_model1_density2/testing5-a2-n_minus_separation10.00000charge1.txt")

data6_plus = np.loadtxt("../run_results/compare_epsilonLJ/charge4-12/35.51-88.78_C4MIM+_TFSI-_model1_density2/testing6-a2-n_plus_separation10.00000charge1.txt")
data6_minus = np.loadtxt("../run_results/compare_epsilonLJ/charge4-12/35.51-88.78_C4MIM+_TFSI-_model1_density2/testing6-a2-n_minus_separation10.00000charge1.txt")

x1_plus = list(data1_plus[:,0])
y1_plus = list(data1_plus[:,1])
#y1_plus[0] = 0.0

x1_minus = list(data1_minus[:,0])
y1_minus = data1_minus[:,1]
#y1_plus[0] = 0.0

x2_plus = list(data2_plus[:,0])
y2_plus = data2_plus[:,1]
#y2_plus[0] = 0.0

x2_minus = list(data2_minus[:,0])
y2_minus = data2_minus[:,1]
#y2_minus[0] = 0.0

x3_plus = list(data3_plus[:,0])
y3_plus = data3_plus[:,1]
#y3_plus[0] = 0.0

x3_minus = list(data3_minus[:,0])
y3_minus = data3_minus[:,1]
#y3_minus[0] = 0.0

x4_plus = list(data4_plus[:,0])
y4_plus = data4_plus[:,1]
#y4_plus[0] = 0.0

x4_minus = list(data4_minus[:,0])
y4_minus = data4_minus[:,1]
#y4_minus[0] = 0.0

x5_plus = list(data5_plus[:,0])
y5_plus = data5_plus[:,1]
#y5_plus[0] = 0.0

x5_minus = list(data5_minus[:,0])
y5_minus = data5_minus[:,1]
#y5_minus[0] = 0.0

x6_plus = list(data6_plus[:,0])
y6_plus = data6_plus[:,1]
#y6_plus[0] = 0.0

x6_minus = list(data6_minus[:,0])
y6_minus = data6_minus[:,1]
#y6_minus[0] = 0.0


#y_contactthm2 = y_contactthm.copy()
#for i in range(len(y_contactthm2)):
#    y_contactthm2[i] = 0.0 
#print y1_plus

fig = plt.figure()
ax = fig.gca(projection='3d')

zeros=[0.0]*len(y1_plus)
half=[0.5]*len(y2_plus)
ones=[1.0]*len(y3_plus)
onehalf=[1.5]*len(y4_plus)
twos=[2.0]*len(y5_plus)
twohalf=[2.5]*len(y6_plus)



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

ax.scatter(x1_plus, zeros, y1_plus, c='k')
#ax.scatter(x1_minus, zeros, y1_minus, c='k')

ax.scatter(x2_plus, half, y2_plus, c='k')
#ax.scatter(x2_minus, half, y2_minus, c='k')

ax.scatter(x3_plus, ones, y3_plus, c='k')
#ax.scatter(x3_minus, ones, y3_minus, c='k')

ax.scatter(x4_plus, onehalf, y4_plus, c='k')
#ax.scatter(x4_minus, onehalf, y4_minus, c='k')

ax.scatter(x5_plus, twos, y5_plus, c='k')
#ax.scatter(x5_minus, twos, y5_minus, c='k')

ax.scatter(x6_plus, twohalf, y6_plus, c='k')
#ax.scatter(x6_minus, twohalf, y6_minus, c='k')

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


verts.append(list(zip(x1_minus, y1_minus)))
verts2.append(list(zip(x2_minus, y2_minus)))
verts3.append(list(zip(x3_minus, y3_minus)))
verts4.append(list(zip(x4_minus, y4_minus)))
verts5.append(list(zip(x5_minus, y5_minus)))
verts6.append(list(zip(x6_minus, y6_minus)))

poly = PolyCollection(verts, facecolors=['k'])
poly2 = PolyCollection(verts2, facecolors=['k'])
poly3 = PolyCollection(verts3, facecolors=['k'])
poly4 = PolyCollection(verts4, facecolors=['k'])
poly5 = PolyCollection(verts5, facecolors=['k'])
poly6 = PolyCollection(verts6, facecolors=['k'])


poly.set_alpha(0.3)
poly2.set_alpha(0.4)
poly3.set_alpha(0.5)
poly4.set_alpha(0.6)
poly5.set_alpha(0.7)
poly6.set_alpha(0.8)

ax.add_collection3d(poly, zs=zs, zdir='y')
ax.add_collection3d(poly2, zs=zs2, zdir='y')
ax.add_collection3d(poly3, zs=zs3, zdir='y')
ax.add_collection3d(poly4, zs=zs4, zdir='y')
ax.add_collection3d(poly5, zs=zs5, zdir='y')
ax.add_collection3d(poly6, zs=zs6, zdir='y')


ax.set_xlabel(r'$z/\sigma$', labelpad=1.5, fontsize=13)
ax.set_xlim3d(0, 10)
ax.tick_params(axis='x', which='major', pad=-2)

ax.set_ylabel(r'$\epsilon_{LJ}^{pw}/\rho\sigma^{3}\epsilon_{LJ}^{pp}$', labelpad=4, fontsize=13)
#ax.set_ylim3d(-0.2, 2.7)
ax.tick_params(axis='y', which='major', pad=-2)

ax.set_zlabel(r'$n\sigma^3$',labelpad=8, fontsize=13)
#ax.set_zlim3d(-0.0003, 0.0011)
ax.set_zlim3d(-0.001, 0.08)
ax.tick_params(axis='z', which='major', pad=4)

#plt.ticklabel_format(style='sci', axis='z', scilimits=(0,0))

savefig("DENSITY_TESTING.pdf",bbox_inches='tight')

plt.show()

