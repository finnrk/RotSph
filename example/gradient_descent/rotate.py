from pymatgen.io.vasp import Procar
from pymatgen.electronic_structure.plotter import Spin
import numpy as np
import rotsph
from multiprocessing import Pool
import sys





#assign default parameters         
ion = 0
filein = 'PROCAR'
fileout = None
# read parameters from command line
if len(sys.argv) > 1:
    filein = (sys.argv[1])
    if len(sys.argv) > 2:
        ion = int(sys.argv[2])
        if len(sys.argv) > 3:
            fileout = sys.argv[3]        


#read projections from PROCAR file
procar = Procar(filein)
nbands = procar.nbands
nkpoints = procar.nkpoints
nions = procar.nions
nspins = procar.nspins

if nspins == 1:
    projs2 = [procar.data[Spin.up]]
    phases = [procar.phase_factors[Spin.up]]
else:
    projs2 = [procar.data[Spin.up], procar.data[Spin.down]]
    phases = [procar.phase_factors[Spin.up], procar.phase_factors[Spin.down]]


#normalise phases
with np.nditer(phases[0], op_flags=['readwrite']) as it:
    for x in it:
        if np.abs(x) != 0 and not np.isnan(x):
            x[...] = x/np.abs(x)
if nspins == 2:
    with np.nditer(phases[1], op_flags=['readwrite']) as it:
        for x in it:
            if np.abs(x) != 0 and not np.isnan(x):
                x[...] = x/np.abs(x)

#calculate complex projections as product of magnitude and phase
if nspins == 1:
    projs = [np.sqrt(projs2[0])*phases[0]]
else:
    projs = [np.sqrt(projs2[0])*phases[0],np.sqrt(projs2[1])*phases[1]]


#--------------------
# random sampling
#--------------------

print("Starting random sampling")

# must create 1-paramater version of metric for use in threading
def metric_1_param(q):
   return rotsph.metric(projs, ion, q)


num_samples = 5000
#generate rotations using quaternions randomly distributed over S^3
unnormalisedqs = np.random.normal(0,1,size=(4,num_samples))
normalisedqs = [np.array(unnormalisedqs[:, i])/np.linalg.norm(unnormalisedqs[:, i]) for i in range(num_samples)]


#evaluating the metric corresponding to each quaternion using multithreading
with Pool(20) as pool:
    metric_list = pool.map(metric_1_param, normalisedqs)
minmet = 1000000
#find rotation corresponding to the smallest metric
for i in range(num_samples):
    if metric_list[i] < minmet:
        minmet = metric_list[i]
        qfin = normalisedqs[i]



#convert from quaternion to euler angles
[a,b,g] = rotsph.quaternion_to_euler(*qfin)


#-----------------
#gradient descent
#-----------------

a,b,g = rotsph.gradient_descent(projs, ion, a,b,g)

#------------------------
#write output PROCAR
#------------------------

#rotation matrices for projections

if fileout is not None:
    rotsph.write_procar(procar, a,b,g, fileout)



#reverse rotation from S'->S to S->S'
[w,x,y,z] = rotsph.euler_to_quaternion(a,b,g)
a,b,g = rotsph.quaternion_to_euler(-w,x,y,z)


#note we swap a and g to convert from xyz to zyx convention  (so that these angles can be copy and pasted into pbandplot.py)
print("\n\nEuler angles: \nalpha, beta, gamma: ", g*180/np.pi,b*180/np.pi,a*180/np.pi)

print("\n\nRotation matrix: ")
for i in np.transpose(rotsph.euler_to_matrix(a,b,g).T):
    print(' '.join(map(str, i)))
print("\n")

