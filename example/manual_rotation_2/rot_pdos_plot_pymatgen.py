import numpy as np
from rotsph import rotsph
import sys
from pymatgen.electronic_structure.core import Spin
from pymatgen.io import vasp

#---------------------------------------
# read data
#---------------------------------------
# read the input
atom_idx = int(sys.argv[1])
orbital = int(sys.argv[2])

if orbital < 4:
    orb_idx = 1
elif orbital < 9: 
    orb_idx = 2
else:
    orb_idx = 3

z_angle_degree = float(sys.argv[3])
y_angle_degree = float(sys.argv[4])
x_angle_degree = float(sys.argv[5])

dos_min = float(sys.argv[6])
dos_max = float(sys.argv[7])
step = int(sys.argv[8])
sigma = float(sys.argv[9])
spin = float(sys.argv[10])

if spin == 1:
    spin = Spin.up
elif spin == -1:
    spin = Spin.down
else:
    print('Error: spin should be either 1 or -1')
    exit()

print('Reading data...')
#kpt, kpt_weight, bnd_energy, bnd_occ, proj_data = read_PROCAR12('PROCAR',report=False)
procar = vasp.Procar('PROCAR')
eigenval = vasp.outputs.Eigenval('EIGENVAL')

proj_data = procar.phase_factors[spin] 
kpt = np.array(eigenval.kpoints)
kpt_weight = np.array(eigenval.kpoints_weights)
bnd_energy = eigenval.eigenvalues[spin]
bnd_occ = bnd_energy[:,:,1]
bnd_energy = bnd_energy[:,:,0]

assert atom_idx < proj_data.shape[2]
#  assert orbital in ['p','d','f']
assert orbital in range(1,16) 

#---------------------------------------
# generate transformation matrix
#---------------------------------------
#  print('Generating transformation matrix...')
z_angle = z_angle_degree/180.*np.pi
rot_z =np.array([[np.cos(z_angle),-np.sin(z_angle),0],\
        [np.sin(z_angle),np.cos(z_angle),0],[0,0,1]])
y_angle = y_angle_degree/180.*np.pi
rot_y=np.array([[np.cos(y_angle),0,np.sin(y_angle)],[0,1,0],\
        [-np.sin(y_angle),0,np.cos(y_angle)]])
x_angle = x_angle_degree/180.*np.pi
rot_x=np.array([[1,0,0],[0,np.cos(x_angle),-np.sin(x_angle)],\
        [0,np.sin(x_angle),np.cos(x_angle)]])

# construct the rotation matrix of coordinate system
rot_mat = np.dot(np.dot(rot_z,rot_y),rot_x)

# get the transformation matrix for d orbitals
T_mat = rotsph.get_R_mat(orb_idx,rot_mat)

#---------------------------------------
# rotating PDOS
#---------------------------------------
print('Rotating PDOS...')
# initialize array for rotated weight
weight_rotated = np.zeros((len(kpt),len(bnd_energy[0])),dtype=complex)

# extract the weight of d orbitals
weight_rotated_atom = proj_data[:,:,atom_idx,orb_idx**2:orb_idx**2+2*orb_idx+1]
# extract T_mat colum for dz
T_mat_dz = T_mat[orbital-orb_idx**2,:]

# perform the rotation
weight_rotated = np.einsum('ijl,l',weight_rotated_atom,T_mat_dz)
# square the weight
weight_rotated = weight_rotated * np.conj(weight_rotated)


#---------------------------------------
# construct the DOS
#---------------------------------------
def get_dos(E, bnd_energy, kpt_weight, proj_data, sigma):
    dos = 0
    for ikpt in range(kpt.shape[0]):
        for ibnd in range(bnd_energy.shape[1]):
            a = E-bnd_energy[ikpt,ibnd]
            dos += np.exp(-0.5*(a/sigma)**2)*kpt_weight[ikpt]*np.real(proj_data[ikpt,ibnd])
    return dos


#dos = get_dos(1, bnd_energy, kpt_weight, weight_rotated, 0.1)
#print(1,dos)
E_DOS = []
for E in np.linspace(dos_min,dos_max,step):
    dos = get_dos(E, bnd_energy, kpt_weight, weight_rotated, sigma)
    #  print(E,dos)
    E_DOS.append([E,dos])

E_DOS = np.array(E_DOS)


#---------------------------------------
# plot data
#---------------------------------------
print('Plotting data...')
import matplotlib.pyplot as plt

fig,ax = plt.subplots(figsize=(1.8,3.5))

ax.plot(E_DOS[:,0],E_DOS[:,1],color='black',linewidth=1.0)
ax.set_xlabel('Energy (eV)')
ax.set_ylabel('DOS')

plt.savefig('DOS.png',bbox_inches='tight',dpi=300)

#---------------------------------------
# write data
#---------------------------------------
print('Writing data...')
with open('pdos.dat','w') as file:
    #  file.write('plot "pdos.dat" w l\n')
    for E in range(E_DOS.shape[0]):
        file.write(f"{E_DOS[E,0]:+10.4f} {E_DOS[E,1]:+10.4f} \n")


#---------------------------------------
# write integrated data
#---------------------------------------
#  for ibnd in range(len(bnd_energy[ikpt])):
#      print(f"{ibnd:10d} {np.sum(weight_rotated[:,ibnd].real):+10.4f}")
#
#  print(f"{z_angle/np.pi*180:10.4f} {x_angle/np.pi*180:10.4f} {np.max(weight_rotated[0,:].real):+10.4f}")




