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

z_angle_degree = float(sys.argv[3])
y_angle_degree = float(sys.argv[4])
x_angle_degree = float(sys.argv[5])

dos_min = float(sys.argv[6])
dos_max = float(sys.argv[7])
step = int(sys.argv[8])
sigma = float(sys.argv[9])

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


#  dos = get_dos(1, bnd_energy, kpt_weight, weight_rotated, 0.1)
#  print(1,dos)
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
    file.write('plot "pdos.dat" w l\n')
    for E in range(E_DOS.shape[0]):
        file.write(f"{E_DOS[E,0]:+10.4f} {E_DOS[E,1]:+10.4f} \n")


#---------------------------------------
# write integrated data
#---------------------------------------
#  for ibnd in range(len(bnd_energy[ikpt])):
#      print(f"{ibnd:10d} {np.sum(weight_rotated[:,ibnd].real):+10.4f}")
#
#  print(f"{z_angle/np.pi*180:10.4f} {x_angle/np.pi*180:10.4f} {np.max(weight_rotated[0,:].real):+10.4f}")




