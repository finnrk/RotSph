from pymatgen.io.vasp import Procar
from pymatgen.electronic_structure.plotter import Spin
import numpy as np
import time
import rotsph
from multiprocessing import Pool
import sys





"""
 The function which quantifies the suitabilty of an orientation. We aim to minimise this.
 Parameter: 
 q : np.array of size (4,1)
     A quaternion, q, which parametrises the rotation.

 Returns: 
 float 
    A value >= 1 where 1 is the metric of a completely diagonal system.
 """
def metric(q):
    out = 0
    #find rotation matrices for spherical harmonics / projectors
    rot_mat = rotsph.quaternion_to_matrix(*q)
    pmat = rotsph.get_R_mat(1,rot_mat.T)
    dmat = rotsph.get_R_mat(2,rot_mat.T)

    total_projection = 0

    #calculate metric - to be minimised
    for s in range(nspins):
        for k in np.linspace(0, nkpoints-1,min(7,nkpoints)):
            k = int(k)
            for b in range(nbands):
                #rotate projectors
                pin = projs[s][k,b,ion,1:4]
                pout = np.square(np.abs(np.matmul(pmat, pin)))
                ptot = sum(pout)

                din = projs[s][k,b,ion,4:]
                dout = np.square(np.abs(np.matmul(dmat, din)))
                dtot = sum(dout)
                total_projection += ptot + dtot

            #increase output for off-diagonal terms
            #define the specific function used for the metric (can be any concave function)
                metfun = lambda x : np.sqrt(x)
                if dtot >0:
                    dout /= dtot
                    for i in dout:
                        out += metfun(i)*dtot
                if ptot >0:
                    pout /= ptot
                    for i in pout:
                        out += metfun(i)*ptot

    #ideal metric would be 1.
    return out/(total_projection)

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
p = Procar(filein)
nbands = p.nbands
nkpoints = p.nkpoints
nions = p.nions
nspins = p.nspins



if nspins == 1:
    projs2 = [p.data[Spin.up]]
    phases = [p.phase_factors[Spin.up]]
else:
    projs2 = [p.data[Spin.up], p.data[Spin.down]]
    phases = [p.phase_factors[Spin.up], p.phase_factors[Spin.down]]






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

#calculate complex rojections as product of magnitude and phase
if nspins == 1:
    projs = [np.sqrt(projs2[0])*phases[0]]
else:
    projs = [np.sqrt(projs2[0])*phases[0],np.sqrt(projs2[1])*phases[1]]



time_sample_start = time.time()

#generate rotations using quaternions randomly distributed over S^3
numpoints = 5000
print("Starting random sampling")
unnormalisedqs = np.random.normal(0,1,size=(4,numpoints))
normalisedqs = [np.array(unnormalisedqs[:, i])/np.linalg.norm(unnormalisedqs[:, i]) for i in range(numpoints)]


#evaluating the metric corresponding to each quaternion using multithreading
threadstarttime = time.time()
if __name__ == '__main__':
    with Pool(20) as pool:
        metric_list = pool.map(metric, normalisedqs)
threadendtime = time.time()


minmet = 1000000
#find rotation corresponding to the smallest metric
for i in range(numpoints):
    if metric_list[i] < minmet:
        minmet = metric_list[i]
        qfin = normalisedqs[i]

time_sample_end = time.time()
time_sample = time_sample_end-time_sample_start


#convert from quaternion to euler angles
[a,b,g] = rotsph.quaternion_to_euler(*qfin)


#-----------------
#gradient descent
#-----------------
print("Starting gradient descent")
grada, gradb, gradg = 200,200,200
eta = 0.1
da = 0.00001
db = 0.00001
dg = 0.00001
a1, b1, g1 = a+1,b+1,g+1
time_sample_start = time.time()

count = 0
while np.abs(a-a1) > 1e-4 or  np.abs(b-b1) > 1e-4 or np.abs(g-g1) > 1e-4:
    count += 1
    a1,b1,g1 = a,b,g
    grada = (metric(rotsph.euler_to_quaternion(a+da,b,g))-metric(rotsph.euler_to_quaternion(a-da,b,g)))/(2*da)
    a = a - grada*eta
    gradb = (metric(rotsph.euler_to_quaternion(a,b+db,g))-metric(rotsph.euler_to_quaternion(a,b-db,g)))/(2*db)
    b = b - gradb*eta
    gradg = (metric(rotsph.euler_to_quaternion(a,b,g+dg))-metric(rotsph.euler_to_quaternion(a,b,g-dg)))/(2*dg)
    g = g - gradg*eta
    eta = eta/(1+0.01*count)
time_sample_end = time.time()


#------------------------
#write output PROCAR
#------------------------

#rotation matrices for projections
rot_mat = rotsph.euler_to_matrix(a,b,g)
pmat = rotsph.get_R_mat(1,rot_mat.T)
dmat = rotsph.get_R_mat(2,rot_mat.T)




#write output PROCAR
if fileout is not None:
    with open(fileout, "w") as w:
        w.write("PROCAR lm decomposed + phase; opt. projector for interval (eV)   -15.0000000    15.0000000\n")
        for spin in range(nspins):
            w.write(f'# of k-points:  {nkpoints}         # of bands:   {nbands}         # of ions:    {nions}\n')
            for kpointwrite in range(nkpoints):
                w.write(f'\n k-point     {kpointwrite+1} :  {p.kpoints[kpointwrite,0]} {p.kpoints[kpointwrite,1]} {p.kpoints[kpointwrite,2]}      weight = {p.weights[kpointwrite]}\n\n')
                for bandwrite in range(nbands):
                    tot = np.zeros(10)
                    w.write(f'band      {bandwrite+1} # energy  {p.eigenvalues[Spin.up if spin == 0 else 1][kpointwrite,bandwrite]} # occ.   {np.abs(p.occupancies[Spin.up if spin == 0 else 1][kpointwrite,bandwrite])}\n')
                    w.write("\nion      s     py     pz     px    dxy    dyz    dz2    dxz  x2-y2    tot\n")
                    #store calculated projections so they do not need to be recalculated
                    sout_list = []
                    pout_list = []
                    dout_list = []
                    total_list = []
                    #write square of projections
                    for ion in range(nions):
                        s = projs[spin][kpointwrite,bandwrite,ion,0]
                        pin = projs[spin][kpointwrite,bandwrite,ion,1:4]
                        pout = np.matmul(pmat, pin)
                        din = projs[spin][kpointwrite,bandwrite,ion,4:]
                        dout = np.matmul(dmat, din)
                        sout_list.append(s)
                        pout_list.append(pout)
                        dout_list.append(dout)
                        s2 = np.abs(s)**2
                        dout2 = np.round(np.abs(dout)**2,3)
                        pout2 = np.round(np.abs(pout)**2,3)
                        #total projection onto all orbitals for a given ion
                        total = s2 + np.sum(pout2) + np.sum(dout2)
                        total_list.append(total)
                        tot[0] += s2
                        tot[-1] += total
                        poutstring = ''
                        for i in range(3):
                            poutstring = poutstring + "{:.3f}".format(pout2[i]) + "  "
                            tot[1+i] += pout2[i]
                        doutstring = ''
                        for i in range(5):
                            doutstring = doutstring + "{:.3f}".format(dout2[i]) + "  "
                            tot[4+i] += dout2[i]
                        out = ''.join(("   ", " "*(1-int((ion+1)/10)), str(ion+1), "  ","{:.3f}".format(s2), "  ", poutstring, doutstring, "{:.3f}".format(total), "\n"))
                        w.write(out)
                        if ion == nions-1:
                            out = 'tot    '
                            for i in range(10):
                                out += "{:.3f}".format(tot[i]) + "  "
                            w.write((out+ "\n"))

                    #write complex phases - normalised to the projection squared
                    #I can't find what they are normalised to in the PROCAR file outputted by VASP
                    w.write("ion          s             py             pz             px             dxy            dyz            dz2            dxz          dx2-y2  \n")
                    for ion in range(nions):
                        #use stored projections instead of re-calculating
                        s = sout_list[ion]
                        pout = pout_list[ion]
                        dout = dout_list[ion]
                        total = total_list[ion]
                        #fixes formatting issue for s = 0 - 0j
                        if s.imag == 0 and s.real == 0:
                            s = 0
                        sstring = " "*int(2+np.floor(0.5*np.sign(s.real))) + "{:.3f}".format(s.real) + " "*int(2+np.floor(0.5*np.sign(s.imag))) + "{:.3f}".format(s.imag)
                        poutstring = ''
                        for i in range(3):
                            poutstring = poutstring + " "*int(3+np.floor(0.5*np.sign(pout[i].real))) + "{:.3f}".format(pout[i].real) + " "*int(2+np.floor(0.5*np.sign(pout[i].imag))) + "{:.3f}".format(pout[i].imag)
                        doutstring = ''
                        for i in range(5):
                            doutstring = doutstring + " "*int(3+np.floor(0.5*np.sign(dout[i].real))) + "{:.3f}".format(dout[i].real) + " "*int(2+np.floor(0.5*np.sign(dout[i].imag))) + "{:.3f}".format(dout[i].imag)
                        out = ''.join(("   ", " "*(1-int((ion+1)/10)), str(ion+1), sstring, poutstring, doutstring, "   ", "{:.3f}".format(total), "\n"))
                        w.write(out)

                        #can use totals calculated in the previous part
                        if ion == nions-1:
                            out = 'charge '
                            for i in range(10):
                                out += "{:.3f}".format(tot[i]) + "          "
                            w.write((out+ "\n\n"))

    time_write_end = time.time()



timetakengrad = time_sample_end-time_sample_start
    


# print("\n\nEuler angles in radians: \nalpha, beta, gamma: ", a,b,g)
# print("\n\nEuler angles in degrees: \nz,y,x: ", g*180/np.pi,b*180/np.pi,a*180/np.pi)



# print("\nQuaternion: ", rotsph.euler_to_quaternion(a,b,g))

#reverse rotation from S'->S to S->S'
[w,x,y,z] = rotsph.euler_to_quaternion(a,b,g)
a,b,g = rotsph.quaternion_to_euler(-w,x,y,z)


#note we swap a and g to convert from xyz to zyx convention  (so that these angles can be copy and pasted into pbandplot.py)
print("\n\nEuler angles: \nalpha, beta, gamma: ", g*180/np.pi,b*180/np.pi,a*180/np.pi)

print("\n\nRotation matrix: ")
for i in np.transpose(rotsph.euler_to_matrix(a,b,g).T):
    print(' '.join(map(str, i)))
print("\n")

# inputmat = np.array([[ 0.46436496, -0.48844994,  0.7387705 ],
#   [ 0.76637311, -0.19646982, -0.61161415],
#   [ 0.443889,    0.85018602,  0.28310119]])

# inputmat = np.array([[   0.6923367,  0.5989676, -0.4023776],
#   [-0.5374349,  0.8001396,  0.2663462],
#    [0.4814910,  0.0318505,  0.8758721]])

# print("\n\nRotation matrix: ")
# for i in np.transpose(inputmat.T):
#     print(' '.join(map(str, i)))
# print("\n")


# for i in np.matmul(inputmat, rotsph.euler_to_matrix(a,b,g)):
#     print(' '.join(map(str, i)))
# print("\n")

# for i in np.matmul(rotsph.euler_to_matrix(a,b,g),inputmat):
#     print(' '.join(map(str, i)))
# print("\n")

# for i in np.matmul(inputmat, np.identity(3)):
#     print(' '.join(map(str, i)))
# print("\n")

# for i in np.matmul(np.identity(3),inputmat):
#     print(' '.join(map(str, i)))
# print("\n")


# print("Time taken on mesh search:", time_sample)
# print("Time taken on gradient descent:", timetakengrad)
# if write_to_file:
#     print("Time taken on writing to file:", time_write_end-time_write_start)
