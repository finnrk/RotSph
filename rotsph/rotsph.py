import numpy as np
from pymatgen.io.vasp import Procar
from pymatgen.electronic_structure.plotter import Spin

# REF: J. Phys. Chem. 1996, 100, 15, 6342–6347

#% build small u, v, w
def delta(i,j):
    '''
    kronecker delta function.
    '''
    assert type(i) is int;
    assert type(j) is int;
    if i == j:
        return 1
    else:
        return 0

def u_mat(l,m,m_):
    '''
    Calculate u matrix.
    '''
    assert type(l) is int;
    assert type(m) is int;
    assert type(m_) is int;
    #
    if np.abs(m_) < l:
        return ((l+m)*(l-m)/((l+m_)*(l-m_)))**(0.5)
    elif np.abs(m_) == l:
        return ((l+m)*(l-m)/(2*l*(2*l-1)))**(0.5)
    else:
        print('Error in generating u_matrix');
        exit()

def v_mat(l,m,m_):
    '''
    Calculate v matrix.
    '''
    assert type(l) is int;
    assert type(m) is int;
    assert type(m_) is int;
    #
    if np.abs(m_) < l:
        return 0.5*( \
                   ((1+delta(m,0))*(l+np.abs(m)-1)*(l+np.abs(m))) \
                   /((l+m_)*(l-m_)) \
                   )**0.5 \
                  *(1-2*delta(m,0))
    elif np.abs(m_) == l:
        return 0.5*(\
                   ((1+delta(m,0))*(l+np.abs(m)-1)*(l+np.abs(m))) \
                   /(2*l*(2*l-1))\
                   )**0.5 \
                  *(1-2*delta(m,0))
    else:
        print('Error in generating v_matrix');
        exit()

def w_mat(l,m,m_):
    '''
    Calculate w matrix.
    '''
    assert type(l) is int;
    assert type(m) is int;
    assert type(m_) is int;
    #
    if np.abs(m_) < l:
        return -0.5*(((l-np.abs(m)-1)*(l-np.abs(m)))/((l+m_)*(l-m_)))**0.5 \
                   *(1-delta(m,0))
    elif np.abs(m_) == l:
        return -0.5*(((l-np.abs(m)-1)*(l-np.abs(m)))/(2*l*(2*l-1)))**0.5 \
                   *(1-delta(m,0))
    else:
        print('Error in generating v_matrix');
        exit()

#% build big U V W
def p_mat(i,l,mu,m_, rot):
    '''
    Calculate p matrix element.
    '''
    assert type(i) is int;
    assert type(l) is int;
    assert type(mu) is int;
    assert type(m_) is int;
    
    if np.abs(m_) < l:
        return rot[i][0]*R_mat(l-1,mu,m_, rot)
    elif m_==l:
        if (np.abs(mu)>l-1 or np.abs(m_-1)>l-1 or np.abs(-m_+1)>l-1):
            # In REF, mu and m_ are allowed to be larger than l.
            # However, Eq. (2.5b) in REF, if m is larger than l, P_l^m(t) = 0
            # hence we are setting everything to zero if abs(mu) or 
            # abs(m_-1) and abs(m_-1)) are larger than l.
            return 0
        else:
            # P functions defined in Table 1 for m'=l is not correct.
            # The correct form can be found in Eq. 7.4a.
            return rot[i][1]*R_mat(l-1,mu,m_-1, rot) \
                  -rot[i][-1]*R_mat(l-1,mu,-m_+1, rot)
    elif m_==-l:
        if (np.abs(mu)>l-1 or np.abs(m_+1)>l-1 or np.abs(-m_-1)>l-1):
            # In REF, mu and m_ are allowed to be larger than l.
            # However, Eq. (2.5b) in REF, if m is larger than l, P_l^m(t) = 0
            # hence we are setting everything to zero if abs(mu) or 
            # abs(m_-1) and abs(m_-1)) are larger than l.
            return 0
        else:
            # P functions defined in Table 1 for m'=-l is not correct.
            # The correct form can be found in Eq. 7.4b.
            return rot[i][1]*R_mat(l-1,mu,m_+1, rot) \
                  +rot[i][-1]*R_mat(l-1,mu,-m_-1, rot)

def bU_mat(l,m,m_, rot):
    '''
    Calculate big U matrix.
    '''
    assert type(l) is int;
    assert type(m) is int;
    assert type(m_) is int;
    
    if m==0:
        return p_mat(0,l,0,m_, rot)
    elif m>0:
        return p_mat(0,l,m,m_, rot)
    elif m<0:
        return p_mat(0,l,m,m_, rot)

def bV_mat(l,m,m_, rot):
    '''
    Calculate big V matrix.
    '''
    assert type(l) is int;
    assert type(m) is int;
    assert type(m_) is int;
    
    if m==0:
        return p_mat(1,l,1,m_, rot) + p_mat(-1,l,-1,m_, rot)
    elif m>0:
        return p_mat(1,l,m-1,m_, rot)*(1+delta(m,1))**0.5 \
              -p_mat(-1,l,-m+1,m_, rot)*(1-delta(m,1))
    elif m<0:
        return p_mat(1,l,m+1,m_, rot)*(1-delta(m,-1)) \
              +p_mat(-1,l,-m-1,m_, rot)*(1+delta(m,-1))**0.5

def bW_mat(l,m,m_, rot):
    '''
    Calculate big W matrix.
    '''
    if m==0:
        return 0
    elif m>0:
        return p_mat(1,l,m+1,m_, rot) + p_mat(-1,l,-m-1,m_, rot)
    elif m<0:
        return p_mat(1,l,m-1,m_, rot) - p_mat(-1,l,-m+1,m_, rot)

#% recursive function
def R_mat(l,m,m_, rot):
    '''
    Recursively calculate transformation matrix.
    '''
    if l==0:
        return 1
    elif l==1:
        if (np.abs(m)>l or np.abs(m_)>l):
            return 0
        else:
            return rot[m][m_]
    else:
        return u_mat(l,m,m_)*bU_mat(l,m,m_, rot) \
              +v_mat(l,m,m_)*bV_mat(l,m,m_, rot) \
              +w_mat(l,m,m_)*bW_mat(l,m,m_, rot)

#%
def reorder_rot(rot_mat):
    '''
    Normalize and reorder the input rotation matrix:
    e.g., [[1,0,0],[0,1,0],[0,0,1]] -> [[0,0,1],[0,1,0],[1,0,0]]

    Here the reordered rot_matrix is returned as a dict to use custome indeces.
    '''
    assert np.allclose(np.dot(rot_mat.T, rot_mat), np.eye(3)) == True, \
        'ERROR: rot_mat is not unitary.'
    
    # first we renormalize rot_mat
    rot_mat = rot_mat/np.linalg.norm(rot_mat,axis=1)
    
    # Note: In the REF, the rotation matrix is defined as
    #                                (r_00, r_01, r_02)
    #        (x', y', z') = (x, y, z)|r_10, r_11, r_12|
    #                                (r_20, r_21, r_22)
    #       However, our input rotation mat is defined (in accordance with the
    #       Euler angle) so that:
    #        (a1')   (r'_00, r'_01, r'_02)(a1)
    #        (a2') = |r'_10, r'_11, r'_12||a2|
    #        (a3')   (r'_20, r'_21, r'_22)(a3)
    #       I.e., np.dot(rot_mat, lat_mat.T) = new_lat_mat.T
    #       Or:
    #        (a1', b1', c1')   (r'_00, r'_01, r'_02)(a1, b2, c3)
    #        (a2', b2', c2') = |r'_10, r'_11, r'_12||a2, b2, c3| 
    #        (a3', b3', c3')   (r'_20, r'_21, r'_22)(a3, b2, c3) 
    # 
    #       Because of this, a transpose is needed to translate
    #       r' to r. 
    rot_mat = rot_mat.T
    
    # NOTE: Also, in REF, the axes are re-ordered: xyz -> yzx
    rot = {
        -1: {
            -1: rot_mat[1,1],
            0 : rot_mat[1,2],
            1 : rot_mat[1,0]
        },
        0: {
            -1: rot_mat[2,1],
            0 : rot_mat[2,2],
            1 : rot_mat[2,0]
        },
        1: {
            -1: rot_mat[0,1],
            0 : rot_mat[0,2],
            1 : rot_mat[0,0]
        }
    }
    return rot


def get_R_mat(l, rot_mat):
    '''
    orbital order:
          m:  -3         |-2   | -1   | 0   | 1    | 2         | 3
    - l = 0:             |     |      | s   |      |           |
    - l = 1:             |     | p_y  | p_z | p_x  |           |
    - l = 2:             |d_xy | d_yz | d_z2| d_xz | d_x2-y2   |
    - l = 3:  f_y(3x2-y2)|f_xyz| f_yz2| f_z3| f_xz2| f_z(x2-y2)| f_x(x2-3y2)

    Order of real spherical harmonics taken from: 
    https://en.wikipedia.org/wiki/Table_of_spherical_harmonics

    The rotated orbtial (\hat{Y_m}) can be constructed from the original orbital 
    (Y_m) as follows:

    \hat{Y_m} = \sum_{m'} R_{mm'} Y_m' 

    Here the R is calculated for each index l m m_.
    '''
    rot = reorder_rot(rot_mat)
    # l = 2
    R = np.empty([2*l+1,2*l+1])
    for m in range(2*l+1):
        for m_ in range(2*l+1):
            R[m,m_] = np.real(R_mat(l,m-l,m_-l, rot))
    return R

#%% test
# rot_mat=np.matrix([[0,1,0],[1,0,0],[0,0,1]])
# rot_mat=np.matrix([[1,1,0],[1,-1,0],[0,0,1]])
# rot_mat=np.matrix([[0,0,1],[0,1,0],[1,0,0]])
# rot_mat=np.matrix([[1,0,0],[0,1,0],[0,0,1]])
# rot_mat=np.matrix([[-1,0,1],[0,1,0],[1,0,1]])

# rot_mat=np.matrix([[1,0,1],[0,1,0],[-1,0,1]])
# R1 = get_R_mat(3,rot_mat)
# print(R1)

# rot_mat=np.matrix([[-1,0,1],[0,1,0],[1,0,1]])
# R2 = get_R_mat(3,rot_mat)
# print(R2)

def quaternion_to_euler(w, x, y, z):
    """
    Convert a normalised quaternion into euler angles (roll, pitch, yaw) in Body 3-2-1 sequence.

    Parameters:
    w, x, y, z : floats 
        The components of a normalised quaternion which characterises a rotation.

    Returns:
    tuple (roll, pitch, yaw) of floats
    roll
        Rotation angle around the x-axis in radians.
    pitch
        Rotation angle around the y-axis in radians.
    yaw 
        Rotation angle around the z-axis in radians.

    """
    t0 = +2.0 * (w * x + y * z)
    t1 = +1.0 - 2.0 * (x**2 + y**2)
    roll_x = np.arctan2(t0, t1)

    t2 = +2.0 * (w * y - z * x)
    t2 = +1.0 if t2 > +1.0 else t2
    t2 = -1.0 if t2 < -1.0 else t2
    pitch_y = np.arcsin(t2)

    t3 = +2.0 * (w * z + x * y)
    t4 = +1.0 - 2.0 * (y**2 + z**2)
    yaw_z = np.arctan2(t3, t4)

    return roll_x, pitch_y, yaw_z  # in radians

def euler_to_quaternion(roll, pitch, yaw):
    """
    Convert Euler angles (roll, pitch, yaw) in Body 3-2-1 sequence to a normalised quaternion.
    
    Parameters:
    roll : float
        Rotation angle around the x-axis in radians.
    pitch : float
        Rotation angle around the y-axis in radians.
    yaw : float
        Rotation angle around the z-axis in radians.
        
    Returns:
    tuple
        A tuple representing the quaternion (w, x, y, z).
    """
    cy = np.cos(yaw * 0.5)
    sy = np.sin(yaw * 0.5)
    cp = np.cos(pitch * 0.5)
    sp = np.sin(pitch * 0.5)
    cr = np.cos(roll * 0.5)
    sr = np.sin(roll * 0.5)

    w = cr * cp * cy + sr * sp * sy
    x = sr * cp * cy - cr * sp * sy
    y = cr * sp * cy + sr * cp * sy
    z = cr * cp * sy - sr * sp * cy

    return (w, x, y, z)


def quaternion_to_matrix(w,x,y,z):
    """
    Converts a normalised quaternion to a rotation matrix

    Parameters:
    w, x, y, z : floats 
        The components of a normalised quaternion which characterises a rotation.

    Returns:
    numpy 3x3 array
        orthogonal rotation matrix
    
    """
    return np.array([[1 - 2*(y**2 + z**2),  2.0*(x*y+w*z),   2.0*(x*z-w*y)],   [2.0*(x*y-w*z), 1 - 2*(x**2 + z**2),  2.0*(w*x+y*z)],  [2.0*(w*y+x*z),    2.0*(y*z-w*x),     1 - 2*(x**2 + y**2)]]).T


def euler_to_matrix(a,b,g):

    """
    Convert Euler angles (roll, pitch, yaw) in Body 3-2-1 sequence to a rotation matrix.
    
    Parameters:
    a : float
        Rotation angle around the x-axis in radians.
    b : float
        Rotation angle around the y-axis in radians.
    g : float
        Rotation angle around the z-axis in radians.

    Returns:
    numpy 3x3 array
        orthogonal rotation matrix  
        """

    return quaternion_to_matrix(*euler_to_quaternion(a,b,g))



def metric(projs, ion, q):
    """
    Quantifies the suitabilty of an orientation. Smaller metric corresponds to a "better" orientation.
    Parameter: 
    q : np.array of size (4,1)
        A quaternion, q, which parametrises the rotation.

    Returns: 
    float 
        A value >= 1 where 1 is the metric of a completely diagonal system.
    """
    out = 0
    #find rotation matrices for spherical harmonics / projectors
    rot_mat = quaternion_to_matrix(*q)
    pmat = get_R_mat(1,rot_mat.T)
    dmat = get_R_mat(2,rot_mat.T)

    total_projection = 0

    nspins = len(projs)
    nkpoints, nbands, nion, norbitals =  np.shape(projs[0])
    #calculate metric - to be minimised
    for s in range(nspins):
        for k in np.linspace(0, nkpoints-1,min(23,nkpoints)):
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
    #print(out/total_projection)
    return out/(total_projection)


def gradient_descent(projs, ion, a,b,g):
    """
    Performs gradient descent to find the rotation which locally minimises metric.
    Parameters: 
    projs : list of length nspins (1 or 2) of np.arrays of size [nkpoints, nbands, nions, norbitals]
        Contains the complex projections of the wavefunction onto each atomic orbital.
        Accessed as projs[spin][kpoint, band, ion, orbital]
    ion : int
        the ion onto which the metric should be minimised
    Euler angles in Body 3-2-1 sequence:
        a : float
            Rotation angle around the x-axis in radians.
        b : float
            Rotation angle around the y-axis in radians.
        g : float
            Rotation angle around the z-axis in radians.
    Returns: 
    tuple
        the euler angles a, b, g which minimise the metric 
    """
    print("Starting gradient descent")
    grada, gradb, gradg = 200,200,200
    eta = 0.1
    da = 0.00001
    db = 0.00001
    dg = 0.00001
    a1, b1, g1 = a+1,b+1,g+1

    count = 0
    while np.abs(a-a1) > 1e-3 or  np.abs(b-b1) > 1e-3 or np.abs(g-g1) > 1e-3:
        count += 1
        a1,b1,g1 = a,b,g
        grada = (metric(projs, ion, euler_to_quaternion(a+da,b,g))-metric(projs, ion, euler_to_quaternion(a-da,b,g)))/(2*da)
        a = a - grada*eta
        gradb = (metric(projs, ion, euler_to_quaternion(a,b+db,g))-metric(projs, ion, euler_to_quaternion(a,b-db,g)))/(2*db)
        b = b - gradb*eta
        gradg = (metric(projs, ion, euler_to_quaternion(a,b,g+dg))-metric(projs, ion, euler_to_quaternion(a,b,g-dg)))/(2*dg)
        g = g - gradg*eta
        eta = eta/(1+0.01*count)
    return a,b,g


def write_procar(procar,a,b,g, fileout):
    """
    Writes a rotated Procar object to a PROCAR file at a specified location
    Parameters: 
    procar : pymatgen.io.vasp.Procar
        The input projections, read from an input PROCAR file
    Euler angles in Body 3-2-1 sequence:
        a : float
            Rotation angle around the x-axis in radians.
        b : float
            Rotation angle around the y-axis in radians.
        g : float
            Rotation angle around the z-axis in radians.
    fileout : string
        location to write file to
    """

    if procar.nspins == 1:
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
    if procar.nspins == 2:
        with np.nditer(phases[1], op_flags=['readwrite']) as it:
            for x in it:
                if np.abs(x) != 0 and not np.isnan(x):
                    x[...] = x/np.abs(x)

    #calculate complex projections as product of magnitude and phase
    if procar.nspins == 1:
        projs = [np.sqrt(projs2[0])*phases[0]]
    else:
        projs = [np.sqrt(projs2[0])*phases[0],np.sqrt(projs2[1])*phases[1]]


    rot_mat = euler_to_matrix(a,b,g)
    pmat = get_R_mat(1,rot_mat.T)
    dmat = get_R_mat(2,rot_mat.T)
    with open(fileout, "w") as w:
        #check if interval is always correct
        w.write("PROCAR lm decomposed + phase; opt. projector for intervall (eV)   -15.0000000    15.0000000\n")
        for spin in range(procar.nspins):
            w.write(f'# of k-points:  {procar.nkpoints}         # of bands:   {procar.nbands}         # of ions:    {procar.nions}\n')
            for kpointwrite in range(procar.nkpoints):
                w.write(f'\n k-point     {kpointwrite+1} :  {procar.kpoints[kpointwrite,0]} {procar.kpoints[kpointwrite,1]} {procar.kpoints[kpointwrite,2]}      weight = {procar.weights[kpointwrite]}\n\n')
                for bandwrite in range(procar.nbands):
                    tot = np.zeros(10)
                    w.write(f'band      {bandwrite+1} # energy  {procar.eigenvalues[Spin.up if spin == 0 else Spin.down][kpointwrite,bandwrite]} # occ.   {np.abs(procar.occupancies[Spin.up if spin == 0 else Spin.down][kpointwrite,bandwrite])}\n')
                    w.write("\nion      s     py     pz     px    dxy    dyz    dz2    dxz  x2-y2    tot\n")
                    #store calculated projections so they do not need to be recalculated
                    sout_list = []
                    pout_list = []
                    dout_list = []
                    total_list = []
                    for ion in range(procar.nions):
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
                        #write square of projections
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
                        if ion == procar.nions-1:
                            out = 'tot    '
                            for i in range(10):
                                out += "{:.3f}".format(tot[i]) + "  "
                            w.write((out+ "\n"))

                            #write complex phases - normalised to the projection squared
                            #I can't find what they are normalised to in the files outputted by VASP
                    w.write("ion          s             py             pz             px             dxy            dyz            dz2            dxz          dx2-y2  \n")
                    for ion in range(procar.nions):
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
                        if ion == procar.nions-1:
                            out = 'charge '
                            for i in range(10):
                                out += "{:.3f}".format(tot[i]) + "          "
                            w.write((out+ "\n\n"))