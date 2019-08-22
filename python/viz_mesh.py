import matplotlib.pyplot as plt
import numpy as np

p = np.loadtxt("Points.txt")
t = np.loadtxt("Facets.txt")
plt.triplot(p[:,0],p[:,1],t-1)
plt.show()
