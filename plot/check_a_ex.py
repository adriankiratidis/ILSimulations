#!/usr/bin/python

import numpy as np
import math
import matplotlib.pyplot as plt


Y_special = 9.0
Psi = 24.0
W = -1.0
Z = -3.0
X = -8.0

# x = nbar * (hs_diameter**3)
x = np.arange(0.05, 0.4, 0.01)
eta_bar = (4.0 * math.atan(1.0) / 6.0 ) * x

GetAExDerivIntegrand = ( (-1.0 / 8.0) * ((np.log(1.0 - eta_bar)/(eta_bar**2)) + (1.0/eta_bar)) ) + \
( ((-Y_special)/8.0) * ((1.0/eta_bar) + (1.0/(1.0 - eta_bar)) + (1.0/((1.0 - eta_bar)**2)) + \
(1.0/((1.0 - eta_bar)**3)) + (np.log(1.0 - eta_bar)/eta_bar**2)) ) + \
(((Psi + (2.0 * Z) - X)/8.0) * (1.0 / ((1.0 - eta_bar)**3))) + \
((W/8.0) * ((1.0/((1.0 - eta_bar)**3)) - (2.0/((1.0 - eta_bar)**2)) + \
(1.0/(1.0 - eta_bar)))) + \
((Psi/8.0) * ((1.0/((1.0 - eta_bar)**2)) - (1.0/((1.0 - eta_bar)**3))))

GetAEx = (( (1.0 - eta_bar) / eta_bar ) * np.log(1.0 - eta_bar) ) + 1.0 + \
( Y_special * ( (np.log(1.0 - eta_bar) * ( ((1.0 - eta_bar)/(eta_bar)) + 1.0 )) - (1.0/(1.0 - eta_bar)) -\
(1.0/(2.0 * ((1.0 - eta_bar)**2))) + (5.0/2.0) ) ) + \
( (Psi + (2.0 * Z) - X) * (((1.0)/(2.0 * ((1.0 - eta_bar)**2))) - 0.5) ) + \
(W * ((1.0/(2.0 * ((1.0 - eta_bar)**2))) - (2.0/(1.0 - eta_bar)) - np.log(1.0 - eta_bar) + 1.5 )) +\
(Psi * ((1.0 / (1.0 - eta_bar)) - (1.0 / (2.0 * ((1.0 - eta_bar)**2))) - 0.5) )






# red dashes, blue squares and green triangles
plt.plot(x, GetAExDerivIntegrand, 'b>')
plt.plot(x, GetAEx, 'gs')
plt.show()
