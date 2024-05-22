import matplotlib.pyplot as plt
import numpy as np


# Read the DOS data
E_DOS_4 = np.loadtxt('4.dat',skiprows=1)
E_DOS_5 = np.loadtxt('5.dat',skiprows=1)
E_DOS_6 = np.loadtxt('6.dat',skiprows=1)
E_DOS_7 = np.loadtxt('7.dat',skiprows=1)
E_DOS_8 = np.loadtxt('8.dat',skiprows=1)

fig,ax = plt.subplots(figsize=(3,3.5))

ax.plot(E_DOS_4[:,0],E_DOS_4[:,1],linewidth=1.0,label='dxy')
ax.plot(E_DOS_5[:,0],E_DOS_5[:,1],linewidth=1.0,label='dyz')
ax.plot(E_DOS_6[:,0],E_DOS_6[:,1],linewidth=1.0,label='dz2')
ax.plot(E_DOS_7[:,0],E_DOS_7[:,1],linewidth=1.0,label='dxz')
ax.plot(E_DOS_8[:,0],E_DOS_8[:,1],linewidth=1.0,label='dx2-y2')

ax.set_xlabel('Energy (eV)')
ax.set_ylabel('DOS')
plt.legend()

plt.savefig('DOS.png',bbox_inches='tight',dpi=300)













