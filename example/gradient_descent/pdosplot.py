from pymatgen.io.vasp import Vasprun
from pymatgen.electronic_structure.plotter import DosPlotter
from pymatgen.io.vasp import Procar
from pymatgen.electronic_structure.plotter import Spin
import numpy as np
import matplotlib.pyplot as plt
import sys

filein = (sys.argv[1])
ion = int(sys.argv[2])
dos_min = float(sys.argv[3])
dos_max = float(sys.argv[4])
samples = int(sys.argv[5])
sigma = float(sys.argv[6])
#0.001

p = Procar(filein)


def energy(kpoint,band, spin):
    if spin == 0:       
        return p.eigenvalues[Spin.up][kpoint,band]
    elif spin == 1:
        return p.eigenvalues[Spin.down][kpoint,band]

#p.data in form {spin: np.array accessed with (k-point index, band index, ion index, orbital index)}
dataup = p.data[Spin.up]
datadown = p.data[Spin.down]
data = [dataup, datadown]


nbands = p.nbands
nkpoints = p.nkpoints
nions = p.nions
nspins = p.nspins
weights = p.weights


occupancies = [p.occupancies[Spin.up], p.occupancies[Spin.down]]

#find fermi energy
E_fermi_min = -1000000000
E_fermi_max = 1000000000
for ikpt in range(nkpoints):
    for ibnd in range(nbands):
        for ispn in range(nspins):
            if occupancies[ispn][ikpt,ibnd] > 0.5:
                if energy(ikpt, ibnd,ispn) > E_fermi_min:
                    E_fermi_min = energy(ikpt, ibnd,ispn)
            if occupancies[ispn][ikpt,ibnd] < 0.5:
                if energy(ikpt, ibnd,ispn) < E_fermi_max:
                    E_fermi_max = energy(ikpt, ibnd,ispn)
E_fermi = 0.5*(E_fermi_max+E_fermi_min)




fig, ax = plt.subplots(figsize = (4,3.5))

#dictionary so both spins of the same orbital have the same color in plot
orb_color = {4:'b', 5:'g', 6:'r', 7:'y', 8:'m'}

#find dos at "samples" points, equally spaced from dos_min to dos_max
#repeat for each orbital
def get_dos(E, orbital_num, spin):
    dos = 0
    for ikpt in range(nkpoints):
        for ibnd in range(nbands):
            a = E-energy(ikpt,ibnd,spin)
            #gaussian smearing
            dos += np.exp(-0.5*(a/sigma)**2)*weights[ikpt]*(data[spin][ikpt,ibnd,ion,orbital_num])
    return dos

for orbital_num in range(4,9):
    E_DOS = []
    for E in np.linspace(dos_min+E_fermi,dos_max+E_fermi,samples):
        dos = get_dos(E,orbital_num, 0)
        #  print(E,dos)
        E_DOS.append([E,dos])
    E_DOS = np.array(E_DOS)
    ax.plot(E_DOS[:,0]-E_fermi,E_DOS[:,1],color = orb_color[orbital_num], linewidth=0.5)
#in case of spin splitting
if nspins == 2:
    for orbital_num in range(4,9):
        E_DOS = []
        for E in np.linspace(dos_min+E_fermi,dos_max+E_fermi,samples):
            dos = get_dos(E,orbital_num, 1)
            #  print(E,dos)
            E_DOS.append([E,-dos])
        E_DOS = np.array(E_DOS)
        ax.plot(E_DOS[:,0]-E_fermi,E_DOS[:,1],color = orb_color[orbital_num], linewidth=0.5)

ax.set_xlabel('Energy (eV)')
ax.set_ylabel('DOS')


plt.axvline(x = 0, color = 'k', linestyle = 'dashed', linewidth=0.3)
plt.legend(["dxy","dyz","dz2","dxz","dx2-y2"])
plt.savefig(filein+'_DOS.png',bbox_inches='tight',dpi=300)
