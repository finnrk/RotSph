import numpy as np
from scipy.spatial.transform import Rotation as R
import sys
import warnings
with warnings.catch_warnings():
    warnings.filterwarnings("ignore", message="divide by zero encountered in divide")

'''
Usage:
python construct_euler.py c0 c1 c2 a0 a1 a2 [force_orth]
'''

# reading x and y axis from command line
c=np.zeros(3)
a=np.zeros(3)
c[0] = float(sys.argv[1]);c[1] = float(sys.argv[2]);c[2] = float(sys.argv[3])
a[0] = float(sys.argv[4]);a[1] = float(sys.argv[5]);a[2] = float(sys.argv[6])

# do we want to force a to be orthogonal to c?
if len(sys.argv) > 7:
	force_orth = True
else:
	force_orth = False

# calculate b
b = np.cross(c,a)

# check if axes are orthogonal
if not np.isclose(np.dot(c,a),0.0,atol=1e-12):
    # force a to be orthogonal to c?
	if force_orth:
		a=-np.cross(c,b)
	else:
		raise ValueError("a and c not orthogonal.")

# normalize axes
a = a/np.linalg.norm(a)
b = b/np.linalg.norm(b)
c = c/np.linalg.norm(c)

# construct the final frame
frame_target = np.array([a,b,c])

# construct the rotation matrix
rot_mat = np.dot(np.linalg.inv(np.eye(3)), frame_target)

# fix reflection
if np.linalg.det(rot_mat) < 0:
    rot_mat = -rot_mat

print(np.dot(rot_mat, np.eye(3)))

# get Euler angles
# Filterout User warning about Gimbo lock.
with warnings.catch_warnings():
    warnings.simplefilter("ignore", category=UserWarning)
    Euler_angles = R.from_matrix(rot_mat.T).as_euler('ZYX', degrees=True)
    print(Euler_angles)
