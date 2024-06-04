from pymatgen.core import Lattice, Structure, Molecule
from scipy.spatial.transform import Rotation as R
import numpy as np
import sys
import warnings
with warnings.catch_warnings():
    warnings.filterwarnings("ignore", message="divide by zero encountered in divide")


atom_idx = int(sys.argv[1])

structure = Structure.from_file("POSCAR")

# find octhedral axes
for r in np.linspace(1,5,100):
    neighbors = structure.get_neighbors(structure[atom_idx],r=r)
    if len(neighbors) == 6:
        break
    elif len(neighbors) > 6:
        print("Too many neighbors.")
        print("Not octahedral symmetry.")
        exit()
		
# measure distance
dist = np.array([structure[0].distance(neighbor) for neighbor in neighbors])

# use the longest and shortest distance to determine z and x
longest = dist.max()
shortest = dist.min()
half = (longest + shortest)/2
if np.sum(dist > half)<3:
    xy_plane = np.where(dist < half)[0]
    z_axis = np.where(dist > half)[0]
elif np.sum(dist < half)<3:
    z_axis = np.where(dist < half)[0]
    xy_plane = np.where(dist > half)[0]
else:
    print("Cannot determine z and x axes.")
    print("Not octahedral symmetry (NOT even approximately).")
    xy_plane = np.where(dist < half)[0]
    z_axis = np.where(dist > half)[0]

z_atom = neighbors[z_axis[0]]
z = z_atom.coords - structure[0].coords
x_atom = neighbors[xy_plane[0]]
x = x_atom.coords - structure[0].coords
#  print(z); print(x)

#  print(structure[0].coords); print(z_atom.coords); print(x_atom.coords)

# calculate b
y = np.cross(z,x)
# check if axes are orthogonal
if not np.isclose(np.dot(z,x),0.0,atol=1e-12):
    # force a to be orthogonal to z?
	if True:
		x=-np.cross(z,y)
	else:
		raise ValueError("a and c not orthogonal.")
		

# normalize axes
x = x/np.linalg.norm(x)
y = y/np.linalg.norm(y)
z = z/np.linalg.norm(z)

# construct the final frame
frame_target = np.array([x,y,z])

# construct the rotation matrix
rot_mat = np.dot(np.linalg.inv(np.eye(3)), frame_target)
#  print(rot_mat)

# fix reflection
if np.linalg.det(rot_mat) < 0:
    rot_mat = -rot_mat
#  print(rot_mat)

# get Euler angles
# Filterout User warning about Gimbo lock.
with warnings.catch_warnings():
    warnings.simplefilter("ignore", category=UserWarning)
    Euler_angles = R.from_matrix(rot_mat.T).as_euler('ZYX', degrees=True)
    print(Euler_angles)

		
