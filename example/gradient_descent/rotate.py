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
write_to_file = 0
filein = 'PROCAR'
fileout = 'PROCAR_OUT'
# read parameters from command line
if len(sys.argv) > 1:
    ion = int(sys.argv[1])
    if len(sys.argv) > 2:
        write_to_file = int(sys.argv[2])
        if len(sys.argv) > 3:
            filein = sys.argv[3]
            if len(sys.argv) > 4:
                fileout = sys.argv[4]        


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
        if np.abs(x) != 0:
            x[...] = x/np.abs(x)
if nspins == 2:
    with np.nditer(phases[1], op_flags=['readwrite']) as it:
        for x in it:
            if np.abs(x) != 0:
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
    with Pool(20) as p:
        metric_list = p.map(metric, normalisedqs)
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

#gradient descent
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



#write output PROCAR


#rotation matrices for projections
rot_mat = rotsph.euler_to_matrix(a,b,g)
pmat = rotsph.get_R_mat(1,rot_mat.T)
dmat = rotsph.get_R_mat(2,rot_mat.T)



bandwrite = 0
kpointwrite = 0
ionwrite = 0
spinwrite = -1

#clear output file
if write_to_file:
    #whether we're writing squares of projections or the complex part
    writesquares = True 
    time_write_start = time.time()
    open(fileout, 'w').close()
    with open(filein) as r:
        with open(fileout, "a") as w:

            for line in list(r):
                if len(line.split()) > 0:
                    if line.split()[0] == "ion":
                        w.write(line)
                    elif line.split()[0] == "k-point":
                        kpointwrite = int(line.split()[1])-1
                        w.write(line)
                    elif line.split()[0] == "band":
                        bandwrite = int(line.split()[1])-1
                        w.write(line)
                        #reset totals
                        tot = np.zeros(10)
                    elif line.split()[0].isnumeric():
                        ionwrite = int(line.split()[0])-1
                        s = projs[spinwrite][kpointwrite,bandwrite,ionwrite,0]
                        pin = projs[spinwrite][kpointwrite,bandwrite,ionwrite,1:4]
                        pout = np.matmul(pmat, pin)
                        din = projs[spinwrite][kpointwrite,bandwrite,ionwrite,4:]
                        dout = np.matmul(dmat, din)

                        s2 = np.abs(s)**2
                        dout2 = np.round(np.abs(dout)**2,3)
                        pout2 = np.round(np.abs(pout)**2,3)
                        total = s2 + np.sum(pout2) + np.sum(dout2)
                        #write square of projections
                        if writesquares:
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

                            out = ''.join(("   ", " "*(1-int((ionwrite+1)/10)), str(ionwrite+1), "  ","{:.3f}".format(s2), "  ", poutstring, doutstring, "{:.3f}".format(total), "\n"))
                            w.write(out)
                            if ionwrite == nions-1:
                                out = 'tot    '
                                for i in range(10):
                                    out += "{:.3f}".format(tot[i]) + "  "
                                w.write((out+ "\n"))

                            #write complex phases - normalised to the projection squared
                            #I can't find what they are normalised to in the files outputted by VASP

                        else:
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


                            out = ''.join(("   ", " "*(1-int((ionwrite+1)/10)), str(ionwrite+1), sstring, poutstring, doutstring, "   ", "{:.3f}".format(total), "\n"))
                            w.write(out)

                            #can use totals calculated in the previous part
                            if ionwrite == nions-1:
                                out = 'charge '
                                for i in range(10):
                                    out += "{:.3f}".format(tot[i]) + "          "
                                w.write((out+ "\n"))

                        if ionwrite == nions-1:
                            ionwrite = 0
                            writesquares = not writesquares
                    elif line.split()[0] == "#":
                        w.write(line)
                        spinwrite += 1
                    elif line.split()[0] == "PROCAR":
                        w.write(line)
                else: 
                    w.write(line)

    time_write_end = time.time()



timetakengrad = time_sample_end-time_sample_start
    


print("\n\nEuler angles: \nalpha, beta, gamma: ", a,b,g)
print("\nQuaternion: ", rotsph.euler_to_quaternion(a,b,g))
print("\nRotation matrix: ")
for i in np.transpose(rotsph.euler_to_matrix(a,b,g)):
    print(' '.join(map(str, i)))
print("\n")


print("Time taken on mesh search:", time_sample)
print("Time taken on gradient descent:", timetakengrad)
if write_to_file:
    print("Time taken on writing to file:", time_write_end-time_write_start)
