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
        for y_angle in np.linspace(z2_1,z2_2,2*int(half_range/dtheta)+1):
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

                # find max proj value of each orbital at gamma point (0) across all bands
                max_val = np.max(weight_rotated[0,:,:].real,axis=0)
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
print('Reading data...')
kpt, kpt_weight, bnd_energy, bnd_occ, proj_data = read_PROCAR12('PROCAR',report=False)

#---------------------------------------
atom_idx = int(sys.argv[1])
assert atom_idx < proj_data.shape[2]
orbital  = sys.argv[2]
assert orbital in ['p','d','f']

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
   Eular_alpha:  {result[0]:5.2f} degrees
   Eular_beta :  {result[1]:5.2f} degrees 
   Eular_gamma:  {result[2]:5.2f} degrees''')

# pz along actual pz
#  print(get_max(proj_data, 3, 'p', 45,0,54.73, 1, 1))

