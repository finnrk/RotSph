import numpy as np
from rotsph import rotsph
import sys

def read_PROCAR12(file_path,report=False):
    with open(file_path, 'r') as file:
        data = file.readlines()
    header = data[1].split()
    num_kpt = int(header[3])
    num_bnd = int(header[7])
    num_ion = int(header[11])
    if report:
        print(num_kpt, num_bnd, num_ion)

    bnd_line = 2*num_ion + 7
    kpt_line = num_bnd*(bnd_line)+1 

    # initialize arrays
    kpt = [] # [kpt,3]
    kpt_weight = [] # [kpt]
    bnd_energy = [] # [kpt,bnd]
    bnd_occ = []   # [kpt,bnd]
    proj_data = [] # [kpt,bnd,ion,orbital]


    start_idx = 2
    for ikpt in range(num_kpt):
        if report:
            print('Reading kpt:',ikpt+1,'/',num_kpt)
        start_idx += ikpt*kpt_line
        kpt.append(np.array(data[start_idx+1].split()[3:6],dtype=float))
        kpt_weight.append(float(data[start_idx+1].split()[8]))
        start_idx += 2
        bnd_energy_tmp = []
        bnd_occ_tmp = []
        proj_data_bnd = []
        for ibnd in range(num_bnd):
            if report:
                print('Reading band:',ibnd+1,'/',num_bnd)
            start_idx += ibnd*bnd_line
            bnd_energy_tmp.append(np.array(data[start_idx+1].split()[4],dtype=float))
            bnd_occ_tmp.append(np.array(data[start_idx+1].split()[7],dtype=float))

            # move to IORBIT=12 phase line
            start_idx += 6+num_ion

            proj_data_ion = []
            for iion in range(num_ion):
                if report:
                    print('Reading ion:',iion+1,'/',num_ion)
                proj_data_tmp = np.array(
                    data[start_idx+iion].split()[1:-1],dtype=float)
                # check if we have even number of data (due to complex)
                if len(proj_data_tmp)%2==0:
                    max_num_orb = int(len(proj_data_tmp)/2)
                else:
                    print('Error: odd number of data')
                    exit()
                proj_data_tmp = proj_data_tmp.reshape((max_num_orb,2))
                # cast to complex
                proj_real = proj_data_tmp[:,0]
                proj_imag = proj_data_tmp[:,1]
                proj_data_tmp = proj_real + 1j*proj_imag

                # put into main array
                proj_data_ion.append(proj_data_tmp)

            # move back
            start_idx -= 6+num_ion
            start_idx -= ibnd*bnd_line

            # put into band array
            proj_data_bnd.append(proj_data_ion)

        # put into band array
        proj_data.append(proj_data_bnd)

        # move back
        start_idx -= ikpt*kpt_line

        # put into main array
        bnd_energy.append(bnd_energy_tmp)
        bnd_occ.append(bnd_occ_tmp)

    return np.array(kpt,dtype=float), \
           np.array(kpt_weight,dtype=float), \
           np.array(bnd_energy,dtype=float), \
           np.array(bnd_occ,dtype=float), \
           np.array(proj_data,dtype=complex)

#---------------------------------------
# automatically finding the best rotation
#---------------------------------------
def build_n_mat(proj_data, bnd_occ, kpt_weight,l, atom_idx):
    n_mat = np.zeros([2*l+1,2*l+1],dtype=complex)
    for m in np.arange(l**2,l**2+(2*l+1)):
        for m_ in np.arange(l**2,l**2+(2*l+1)):
            for ikpt in range(bnd_occ.shape[0]):
                for ibnd in range(bnd_occ.shape[1]):
                    n_mat[m-l**2,m_-l**2] += bnd_occ[ikpt,ibnd]*proj_data[ikpt,ibnd,atom_idx,m]*np.conj(proj_data[ikpt,ibnd,atom_idx,m_])*kpt_weight[ikpt]
                    #  n_mat[m-l**2,m_-l**2] += proj_data[ikpt,ibnd,atom_idx,m]*np.conj(proj_data[ikpt,ibnd,atom_idx,m_])*kpt_weight[ikpt]

    return n_mat

#---------------------------------------
# read data
#---------------------------------------
print('Reading data...')
kpt, kpt_weight, bnd_energy, bnd_occ, proj_data = read_PROCAR12('PROCAR',report=False)


# read the input
atom_idx = int(sys.argv[1])
assert atom_idx < proj_data.shape[2]
orbital = int(sys.argv[2])
#  assert orbital in ['p','d','f']
assert orbital in range(1,16) 
if orbital < 4:
    orb_idx = 1
elif orbital < 9: 
    orb_idx = 2
else:
    orb_idx = 3

# build the n matrix, see Eq. A1 in 10.1103/PhysRevMaterials.5.104402
n_mat = build_n_mat(proj_data, bnd_occ, kpt_weight, 2, 0)
# diagonalize n_mat to get the eigenvec which is the T mat.
# note that the eigenvec is column vector, so we need to transpose it to get the row vector
T_mat = np.linalg.eig(n_mat)[1].T
print('T_mat:',T_mat)
print('Eigenvalues:',np.linalg.eig(n_mat)[0])

# due to the diagonalization, we lose the orbtial order, 
# and there's no way to restore the original order.

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
# construct the k-path
#---------------------------------------
print('Constructing k-path...')
from ase.io import read

atoms = read('POSCAR',format='vasp')
rec_lat = atoms.cell.reciprocal()

kpath = [0]
for ikpt in np.arange(1,len(kpt)):
    k_point = np.dot(kpt[ikpt],rec_lat)
    k_point_prev = np.dot(kpt[ikpt-1],rec_lat)
    k_diff = k_point - k_point_prev
    # absolute value of the difference
    k_diff_abs = np.linalg.norm(k_diff)
    kpath.append(kpath[ikpt-1]+k_diff_abs)

#---------------------------------------
# plot data
#---------------------------------------
print('Plotting data...')
import matplotlib.pyplot as plt

fig,ax = plt.subplots(figsize=(1.8,3.5))

for ibnd in range(len(bnd_energy[0])):
    ax.plot(kpath,bnd_energy[:,ibnd],'k-')
    ax.scatter(kpath,bnd_energy[:,ibnd],s=weight_rotated[:,ibnd].real*50,c="red")
    for high_symmetry_point in kpath[::30]:
        ax.axvline(x=high_symmetry_point,color='gray',linestyle='--')
    ax.set_xticks(kpath[::30])
    ax.set_xticklabels([])
    ax.set_ylabel('Energy [eV]')
    ax.set_xlim([0,kpath[-1]])
    ax.set_ylim([-2,12])

plt.savefig('band_diag.png',bbox_inches='tight',dpi=300)
