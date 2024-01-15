import matplotlib.pyplot as plt
import numpy as np
import sys
import matplotlib.animation as animation

eV = 1.602176634e-19


np.set_printoptions(threshold=sys.maxsize)

N, T, Tc = np.genfromtxt('InputData.txt')

N = int(N)
T = int(T)

a = np.loadtxt('groundstates2d.txt')
a = a.reshape((T,N,N))
temperatures, magnetization = np.genfromtxt('magnetization2d.txt', unpack = True)
temperatures, energies = np.genfromtxt('energies2d.txt', unpack = True)
temperatures, specificHeat = np.genfromtxt('specificHeat2d.txt', unpack = True)
temperatures, susceptibility = np.genfromtxt('susceptibility2d.txt', unpack = True)


print(len(a))
print(len(a[0]))
print(len(a[0][0]))


def init():
    global fig, ax, im
    fig = plt.figure()
    ax = plt.axes()
    im = ax.imshow(a[0], cmap = 'RdBu', vmin = -1, vmax = +1)
    im.set_data(a[0])
    return

def update(t):
 im.set_data(a[t])
 ax.set_title(temperatures[t]/Tc)
 print(t)
 return [im]

init()

anim = animation.FuncAnimation(
                               fig, 
                               update,
                               frames = len(a), 
                               repeat = False,
                               )

anim.save('test_anim.mp4', fps=2 , extra_args=['-vcodec', 'libx264'])




#for t in range(T):
#    plt.title(temperatures[t]/Tc)
#    im = ax.imshow(a[t], animated=True)
#    if t == 0:
#        ax.imshow(a[0])  # show an initial one first
#    ims.append([im])
#
#ani = animation.ArtistAnimation(fig, ims, interval=50, blit=True,
#                                repeat_delay=1000)

figMagn = plt.figure()

plt.plot(temperatures/Tc, magnetization, '.')
plt.xlabel(r'$T/T_c$')
plt.ylabel(r'magnetization')
plt.savefig('MagnvsTemp.pdf')

figEnergies = plt.figure()

plt.plot(temperatures/Tc, energies/eV, '.')
plt.xlabel(r'$T/T_c$')
plt.ylabel(r'energy per spin/ev')
plt.savefig('EnergyvsTemp.pdf')

figSpecificHeat = plt.figure()

plt.plot(temperatures/Tc, specificHeat, '.')
plt.xlabel(r'$T/T_c$')
plt.ylabel(r'$C$')
plt.savefig('specificHeatvsTemp.pdf')

figSusceptibility = plt.figure()

plt.plot(temperatures/Tc, susceptibility, '.')
plt.xlabel(r'$T/T_c$')
plt.ylabel(r'$\chi$')
plt.savefig('SusceptibilityvsTemp.pdf')