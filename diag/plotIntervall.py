import numpy as np 
import matplotlib.pyplot as plt 

x, y , z= np.genfromtxt('intervallHalb.txt', unpack = True)

N = np.linspace(0, len(x), len(x))

plt.plot(N, x, label = "x")
plt.plot(N, y, label = "y")
plt.plot(N, z, label = "z")

plt.legend()
plt.savefig('plotIntervall.pdf')