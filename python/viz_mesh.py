import matplotlib.pyplot as plt
import numpy as np

p = np.loadtxt("POINTS0010.TXT")
t = np.loadtxt("FACETS0010.TXT")
plt.triplot(p[:,0],p[:,1],t-1)
plt.show()
