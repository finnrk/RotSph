import numpy as np
import rotsph
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


if len(sys.argv) > 3:
    z_angle_degree = float(sys.argv[3])
    y_angle_degree = float(sys.argv[4])
    x_angle_degree = float(sys.argv[5])
else:
    z_angle_degree, y_angle_degree, x_angle_degree = 0,0,0

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

plt.savefig(f'band{orbital}.png',bbox_inches='tight',dpi=300)

#---------------------------------------
# write data
#---------------------------------------
# print('Writing data...')
# with open('pband.dat','w') as file:
#     file.write('plot "pband.dat" u 1:2:3 w points lt 1 pt 10 ps variable\n')
#     for ikpt in range(len(kpt)):
#         file.write("\n")
#         for ibnd in range(len(bnd_energy[ikpt])):
#             file.write(f"{kpath[ikpt]:+10.4f} {bnd_energy[ikpt][ibnd]:+10.4f} {weight_rotated[ikpt,ibnd].real:+10.4f} \n")


#---------------------------------------
# write integrated data
#---------------------------------------
#  for ibnd in range(len(bnd_energy[ikpt])):
#      print(f"{ibnd:10d} {np.sum(weight_rotated[:,ibnd].real):+10.4f}")
#  
#  print(f"{z_angle/np.pi*180:10.4f} {x_angle/np.pi*180:10.4f} {np.max(weight_rotated[0,:].real):+10.4f}")




