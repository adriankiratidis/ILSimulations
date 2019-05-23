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

separation = 10
charge = 1

data_1 = np.loadtxt("../bin/test-n_plus_separation" +  str(separation) + ".txt")
data_2 = np.loadtxt("../bin/test-n_minus_separation" + str(separation) + ".txt")

x_1 = data_1[:,0]
y_1 = data_1[:,1]

x_2 = data_2[:,0]
y_2 = data_2[:,1]

#A description of the available plotting characters and colours 
#to be used in place of 'bx' etc can be found here...
#http://matplotlib.org/api/pyplot_api.html#matplotlib.pyplot.plot
plt.plot(x_1+0.02,y_1,'bx',label=r'$n_{+}$')
plt.plot(x_2+0.04,y_2,'go',label=r'$n_{0}$')

#Set the axis labels.  Labelpad option controls the spacing between actual axis and axis label.  The r option tells python to interpret as a raw string literal.
plt.xlabel(r"$z/\sigma$",labelpad=10)
plt.ylabel(r"$n\sigma^{3}$",labelpad=5)

#Uncomment if legend required.
plt.legend(loc='best',ncol=1, numpoints=1, frameon=False)

#Uncomment if title is required
#plt.title(r"Some Title")

savefig("PositiveNegativeSpheres.pdf",bbox_inches='tight')

#Open a window and show the plot
plt.show()
