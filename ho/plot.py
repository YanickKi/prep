import numpy as np 
import matplotlib.pyplot as plt 

t, f = np.genfromtxt('data.txt', unpack = True)

plt.plot(t,f)
plt.savefig('plot.pdf')