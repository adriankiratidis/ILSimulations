import matplotlib.pyplot as plt
import numpy as np

t = np.arange(0.16, 1.5, 0.0001)
s = -1.6471*(1.0/t)**6 + 0.6471*np.exp(14.3863*(1-t))
#s = 0.6471*np.exp(14.3863*(1-t))
plt.plot(t, s)

plt.xlabel('r')
plt.ylabel('Potential')
#plt.title('About as simple as it gets, folks')
#plt.grid(True)
#plt.savefig("test.png")
plt.show()
