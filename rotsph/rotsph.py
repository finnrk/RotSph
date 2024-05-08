import numpy as np

#% build small u, v, w
def delta(i,j):
    '''
    kronecker delta function.
    '''
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
        return 0.5*(((1+delta(m,0))*(l+np.abs(m)-1)*(l+np.abs(m)))/((l+m_)*(l-m_)))**0.5*(1-2*delta(m,0))
    elif np.abs(m_) == l:
        return 0.5*(((1+delta(m,0))*(l+np.abs(m)-1)*(l+np.abs(m)))/(2*l*(2*l-1)))**0.5*(1-2*delta(m,0))
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
        return -0.5*(((l-np.abs(m)-1)*(l-np.abs(m)))/((l+m_)*(l-m_)))**0.5*(1-delta(m,0))
    elif np.abs(m_) == l:
        return -0.5*(((l-np.abs(m)-1)*(l-np.abs(m)))/(2*l*(2*l-1)))**0.5*(1-delta(m,0))
    else:
        print('Error in generating v_matrix');
        exit()

#% build big U V W
def p_mat(i,l,mu,m_, rot):
    '''
    Calculate p matrix element.
    '''
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
            # P functions defeined Table 1 for m'=l is not correct.
            # The correct form can be found in Eq. 7.4a.
            return rot[i][1]*R_mat(l-1,mu,m_-1, rot)-rot[i][-1]*R_mat(l-1,mu,-m_+1, rot)
    elif m_==-l:
        if (np.abs(mu)>l-1 or np.abs(m_+1)>l-1 or np.abs(-m_-1)>l-1):
            # In REF, mu and m_ are allowed to be larger than l.
            # However, Eq. (2.5b) in REF, if m is larger than l, P_l^m(t) = 0
            # hence we are setting everything to zero if abs(mu) or 
            # abs(m_-1) and abs(m_-1)) are larger than l.
            return 0
        else:
            # P functions defeined Table 1 for m'=l is not correct.
            # The correct form can be found in Eq. 7.4b.
            return rot[i][1]*R_mat(l-1,mu,m_+1, rot)+rot[i][-1]*R_mat(l-1,mu,-m_-1, rot)

def bU_mat(l,m,m_, rot):
    '''
    Calculate big U matrix.
    '''
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
    if m==0:
        return p_mat(1,l,1,m_, rot) + p_mat(-1,l,-1,m_, rot)
    elif m>0:
        return p_mat(1,l,m-1,m_, rot)*(1+delta(m,1))**0.5 - \
               p_mat(-1,l,-m+1,m_, rot)*(1-delta(m,1))
    elif m<0:
        return p_mat(1,l,m+1,m_, rot)*(1-delta(m,-1)) + \
               p_mat(-1,l,-m-1,m_, rot)*(1+delta(m,-1))**0.5

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
        return u_mat(l,m,m_)*bU_mat(l,m,m_, rot) + \
               v_mat(l,m,m_)*bV_mat(l,m,m_, rot) + \
               w_mat(l,m,m_)*bW_mat(l,m,m_, rot)

#%
def reorder_rot(rot_mat):
    '''
    Normalize and reorder the input rotation matrix:
    e.g., [[1,0,0],[0,1,0],[0,0,1]] -> [[0,0,1],[0,1,0],[1,0,0]]
    '''
    # first we renormalize rot_mat
    rot_mat = rot_mat/np.linalg.norm(rot_mat,axis=1)
    # Note: rotation matrix is defined as
    #                                (r_00, r_01, r_02)
    #        (x', y', z') = (x, y, z)|r_10, r_11, r_12|
    #                                (r_20, r_21, r_22)
    #       However, our input is so that we can write:
    #        (a1', b1', c1')   (r'_00, r'_01, r'_02)(a1, a2, a3)
    #        (a2', b2', c2') = |r'_10, r'_11, r'_12||b1, b2, b3| 
    #        (a3', b3', c3')   (r'_20, r'_21, r'_22)(c1, c2, c3) 
    #       or, in other words:
    #        (a1')   (r'_00, r'_01, r'_02)(a1)
    #        (a2') = |r'_10, r'_11, r'_12||a2|
    #        (a3')   (r'_20, r'_21, r'_22)(a3)
    #       hence, a transpose is needed here: r=r'.T
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
    - l = 0:  s
    - l = 1:  p_y p_z p_x
    - l = 2:  d_xy d_yz d_z2 d_xz d_x2-y2
              ...
    Order taken from: https://en.wikipedia.org/wiki/Table_of_spherical_harmonics#Real_spherical_harmonics

    The rotated orbtial (\hat{Y_m}) can be constructed from the original orbital (Y_m) as follows:
    \hat{Y_m} = \sum_{m'} R_{mm'} Y_m' 
    '''
    rot = reorder_rot(rot_mat)
    # l = 2
    R = np.empty([2*l+1,2*l+1])
    for i in range(2*l+1):
        for j in range(2*l+1):
            # print("main",l,i-l,j-l)
            R[i,j] = np.real(R_mat(l,i-l,j-l, rot))
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
