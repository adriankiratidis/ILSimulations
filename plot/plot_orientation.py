from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Make data
u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)
x = 1.2 * np.outer(np.cos(u), np.sin(v)) + 12
y = 1.2 * np.outer(np.sin(u), np.sin(v)) + 12
z = 1.2 * np.outer(np.ones(np.size(u)), np.cos(v))

# Plot the surface with charge
ax.plot_surface(x, y, z+25*0.048, color='g')
ax.plot_surface(x+0.197539*0.048, y+0.357552*0.048, z+75*0.048, color='g')
ax.plot_surface(x+9.94453*0.048, y+2.342506*0.048, z+124*0.048, color='b')
ax.plot_surface(x+46.01593*0.048, y+26.1004*0.048, z+149*0.048, color='b')

# Plot the surface without charge
# ax.plot_surface(x, y, z+25*0.048, color='g')
# ax.plot_surface(x+2.088*0.048, y+2.084434*0.048, z+75*0.048, color='g')
# ax.plot_surface(x+16.3*0.048, y+2.08738*0.048, z+124*0.048, color='b')
# ax.plot_surface(x+43.97*0.048, y+26.1004*0.048, z+89*0.048, color='b')

# Plot the surface without charge fixed at peak
# ax.plot_surface(x, y, z+89*0.048, color='g')
# ax.plot_surface(x+19.96935*0.048, y+45.80879*0.048, z+89*0.048, color='g')
# ax.plot_surface(x+7.142617*0.048, y+94.12238806*0.048, z+88*0.048, color='b')
# ax.plot_surface(x+39.02221*0.048, y+132.5309*0.048, z+88*0.048, color='b')



ax.set_xlabel('X')
ax.set_xlim3d(0, 24)

ax.set_ylabel('Y')
ax.set_ylim3d(0, 24)

ax.set_zlabel('Z')
ax.set_zlim3d(0, 24)

plt.show()
