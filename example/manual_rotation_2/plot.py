import matplotlib.pyplot as plt
import numpy as np


files = ['6.dat','8.dat']
files_dn = ['6-dn.dat','8-dn.dat']
color = ['tab:blue','tab:orange']
linestyle = ['-','--']
legend = ['$d_{z^2}$','$d_{x^2-y^2}$']
fig,ax = plt.subplots(figsize=(1.8,3.5))

for i in range(2):
    E_DOS = np.loadtxt(files[i])
    ax.plot(E_DOS[:,0],E_DOS[:,1],color=color[i], linestyle=linestyle[i], linewidth=1.0)

for i in range(2):
    E_DOS = np.loadtxt(files_dn[i])
    ax.plot(E_DOS[:,0],-E_DOS[:,1],color=color[i],label=legend[i], linestyle=linestyle[i], linewidth=1.0)

ax.set_xlabel('Energy (eV)')
ax.set_ylabel('DOS')
plt.legend()

plt.savefig('DOS_eg.png',bbox_inches='tight',dpi=300)

files = ['4.dat','5.dat','7.dat']
files_dn = ['4-dn.dat','5-dn.dat','7-dn.dat']
color = ['tab:green','tab:red','tab:purple']
linestyle = ['-','--',':']
legend = ['$d_{xy}$','$d_{xz}$','$d_{yz}$']
fig,ax = plt.subplots(figsize=(1.8,3.5))

for i in range(3):
    E_DOS = np.loadtxt(files[i])
    ax.plot(E_DOS[:,0],E_DOS[:,1],color=color[i], linestyle=linestyle[i], linewidth=1.0)

for i in range(3):
    E_DOS = np.loadtxt(files_dn[i])
    ax.plot(E_DOS[:,0],-E_DOS[:,1],color=color[i],label=legend[i], linestyle=linestyle[i], linewidth=1.0)

ax.set_xlabel('Energy (eV)')
ax.set_ylabel('DOS')
plt.legend()

plt.savefig('DOS_t2g.png',bbox_inches='tight',dpi=300)
