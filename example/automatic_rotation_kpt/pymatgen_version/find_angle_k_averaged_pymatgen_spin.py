import numpy as np
from rotsph import rotsph
import sys
from pymatgen.electronic_structure.core import Spin
from pymatgen.io import vasp


#---------------------------------------
# automatically finding the best rotation
#---------------------------------------
def get_max(proj_data, atom_idx=0, orbital='d', alpha=0, beta=0, gamma=0, half_range=10, dtheta=1):

    # sensible range
    if half_range >= 180:
        half_range = 180
        # y angle only needs to be -90 to 90
        half_range_y = 90
    elif half_range >= 90:
        half_range_y = 90
    else:
        half_range_y = half_range

    # convert to radian
    z1_1 = (alpha-half_range)/180.*np.pi; z1_2=(alpha+half_range)/180.*np.pi
    z2_1 = (beta -half_range_y)/180.*np.pi; z2_2=(beta +half_range_y)/180.*np.pi
    z3_1 = (gamma-half_range)/180.*np.pi; z3_2=(gamma+half_range)/180.*np.pi

    # initialize array
    data=[]
    data_tmp=[]

    for z_angle in np.linspace(z1_1,z1_2,2*int(half_range/dtheta)+1):
        for y_angle in np.linspace(z2_1,z2_2,2*int(half_range_y/dtheta)+1):
            for x_angle in np.linspace(z3_1,z3_2,2*int(half_range/dtheta)+1):
                #---------------------------------------
                # generate transformation matrix
                #---------------------------------------
                #  print('Generating transformation matrix...')
                #  z_angle = 0.0
                #  z_angle = -np.pi/4.
                rot_z =np.array([[np.cos(z_angle),-np.sin(z_angle),0],\
                        [np.sin(z_angle),np.cos(z_angle),0],[0,0,1]])
                #  y_angle = 0
                rot_y=np.array([[np.cos(y_angle),0,np.sin(y_angle)],[0,1,0],\
                        [-np.sin(y_angle),0,np.cos(y_angle)]])
                #  x_angle = 0.0
                #  x_angle = np.arctan(2**0.5/1)
                rot_x=np.array([[1,0,0],[0,np.cos(x_angle),-np.sin(x_angle)],\
                        [0,np.sin(x_angle),np.cos(x_angle)]])

                # construct the rotation matrix of coordinate system
                rot_mat = np.dot(np.dot(rot_z,rot_y),rot_x)


                #---------------------------------------
                # rotating PDOS
                #---------------------------------------
                #  print('Rotating PDOS...')
                # initialize array for rotated weight
                weight_rotated = np.zeros((len(kpt),len(bnd_energy[0])),dtype=complex)

                if orbital == 'p':
                    assert proj_data.shape[3] >= 4
                    # get the transformation matrix for d orbitals
                    T_mat = rotsph.get_R_mat(1,rot_mat)
                    # extract the weight of d orbitals (4:10) of the first atom (0)
                    weight_rotated_atom = proj_data[:,:,atom_idx,1:4]
                    # extract T_mat colum for dz (2)
                    # we just need to fit one orbtial at gamma to get the optimum rotation
                    #  T_mat = T_mat[1,:]
                elif orbital == 'd':
                    assert proj_data.shape[3] >= 9
                    # get the transformation matrix for d orbitals
                    T_mat = rotsph.get_R_mat(2,rot_mat)
                    # extract the weight of d orbitals (4:10) of the first atom (0)
                    weight_rotated_atom = proj_data[:,:,atom_idx,4:9]
                    # extract T_mat colum for dz (2)
                    # we just need to fit one orbtial at gamma to get the optimum rotation
                    #  T_mat = T_mat[2,:]
                elif orbital == 'f':
                    assert proj_data.shape[3] >= 16
                    # get the transformation matrix for d orbitals
                    T_mat = rotsph.get_R_mat(3,rot_mat)
                    # extract the weight of d orbitals (4:10) of the first atom (0)
                    weight_rotated_atom = proj_data[:,:,atom_idx,9:16]
                    # extract T_mat colum for dz (2)
                    # we just need to fit one orbtial at gamma to get the optimum rotation
                    #  T_mat = T_mat[3,:]

                # perform the rotation
                weight_rotated = np.einsum('ijl,kl',weight_rotated_atom,T_mat)
                # square the weight
                weight_rotated = weight_rotated * np.conj(weight_rotated)

                # find max proj value of each orbital across all bands and kpts.
                max_val = np.max(weight_rotated[:,:,:].real,axis=(0,1))
                # and sum up, use this as metric
                max_val = np.sum(max_val)

                #  number_of_significant_bands = np.sum(weight_rotated[0,:].real>0.9*max_val)
                #  metric = np.max(weight_rotated[0,:].real)/\
                #           number_of_significant_bands
                #  metric = np.max(weight_rotated[0,:].real)*\
                #           (len(weight_rotated[0,:])-number_of_significant_bands)

                data.append([z_angle/np.pi*180,\
                             y_angle/np.pi*180,\
                             x_angle/np.pi*180,\
                             max_val])
                             #  np.max(weight_rotated[0,:].real)])
                             #  np.max(weight_rotated[0,:].real),\
                             #  np.argmax(weight_rotated[0,:].real),\
                             #  number_of_significant_bands,\
                             #  metric])


    data = np.array(data,dtype=float)

    return data[np.argmax(data, axis=0)[3]]


#---------------------------------------
# read data
#---------------------------------------
#print('Reading data...')
#kpt, kpt_weight, bnd_energy, bnd_occ, proj_data = read_PROCAR12('PROCAR',report=False)

# read the input
atom_idx = int(sys.argv[1])
orbital  = sys.argv[2]

if orbital < 4:
    orb_idx = 1
elif orbital < 9: 
    orb_idx = 2
else:
    orb_idx = 3

spin = float(sys.argv[3])

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
assert orbital in ['p','d','f']
#assert orbital in range(1,16) 

#---------------------------------------

alpha    = 1
beta     = 1
gamma    = 1

dtheta     = [60,  10, 2, 0.25, 0.05, 0.005]
half_range = [180, 30, 5, 1,    0.1, 0.025]

print('Optimizing angle...')

for iter in range(4):
    print(f"iter: {iter}; accuracy: {dtheta[iter]} degrees")
    result = get_max(proj_data, atom_idx, orbital, alpha,beta,gamma, half_range[iter], dtheta[iter])
    alpha = result[0]
    beta = result[1]
    gamma = result[2]
    #  print(result[3])


# pz can be set to be along px or py as they are identical... p orbtials are weird.
#  print(result[:])
print(f'''
   Euler_alpha:  {result[0]:5.2f} degrees
   Euler_beta :  {result[1]:5.2f} degrees 
   Euler_gamma:  {result[2]:5.2f} degrees''')

# pz along actual pz
#  print(get_max(proj_data, 3, 'p', 45,0,54.73, 1, 1))

#---------------------------------------
# generate transformation matrix
#---------------------------------------
#  print('Generating transformation matrix...')
z_angle = alpha/180.*np.pi
rot_z =np.array([[np.cos(z_angle),-np.sin(z_angle),0],\
        [np.sin(z_angle),np.cos(z_angle),0],[0,0,1]])
y_angle = beta/180.*np.pi
rot_y=np.array([[np.cos(y_angle),0,np.sin(y_angle)],[0,1,0],\
        [-np.sin(y_angle),0,np.cos(y_angle)]])
x_angle = gamma/180.*np.pi
rot_x=np.array([[1,0,0],[0,np.cos(x_angle),-np.sin(x_angle)],\
        [0,np.sin(x_angle),np.cos(x_angle)]])

# construct the rotation matrix of coordinate system
rot_mat = np.dot(np.dot(rot_z,rot_y),rot_x)
T_mat = rotsph.get_R_mat(2,rot_mat)
#  print(T_mat)
