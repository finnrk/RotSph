import numpy as np

#% build small u, v, w
def delta(i,j):
    if i == j:
        return 1
    else:
        return 0

def u_mat(l,m,m_):
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
    # print("p")
    if np.abs(m_) < l:
        return rot[i][0]*R_mat(l-1,mu,m_, rot)
    elif m_==l:
        if (np.abs(mu)>l-1 or np.abs(m_-1)>l-1 or np.abs(-m_+1)>l-1):
            return 0
        else:
            return rot[i][1]*R_mat(l-1,mu,m_-1, rot)-rot[i][-1]*R_mat(l-1,mu,-m_+1, rot)
    elif m_==-l:
        if (np.abs(mu)>l-1 or np.abs(m_+1)>l-1 or np.abs(-m_-1)>l-1):
            return 0
        else:
            return rot[i][1]*R_mat(l-1,mu,m_+1, rot)+rot[i][-1]*R_mat(l-1,mu,-m_-1, rot)

def bU_mat(l,m,m_, rot):
    # print("U")
    if m==0:
        return p_mat(0,l,0,m_, rot)
    elif m>0:
        return p_mat(0,l,m,m_, rot)
    elif m<0:
        return p_mat(0,l,m,m_, rot)

def bV_mat(l,m,m_, rot):
    # print("V")
    if m==0:
        return p_mat(1,l,1,m_, rot) + p_mat(-1,l,-1,m_, rot)
    elif m>0:
        return p_mat(1,l,m-1,m_, rot)*(1+delta(m,1))**0.5 - \
               p_mat(-1,l,-m+1,m_, rot)*(1-delta(m,1))
    elif m<0:
        return p_mat(1,l,m+1,m_, rot)*(1-delta(m,-1)) + \
               p_mat(-1,l,-m-1,m_, rot)*(1+delta(m,-1))**0.5

def bW_mat(l,m,m_, rot):
    # print("W")
    if m==0:
        return 0
    elif m>0:
        return p_mat(1,l,m+1,m_, rot) + p_mat(-1,l,-m-1,m_, rot)
    elif m<0:
        return p_mat(1,l,m-1,m_, rot) - p_mat(-1,l,-m+1,m_, rot)

#% recursive function
def R_mat(l,m,m_, rot):
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
    #e.g., rot_mat=np.matrix([[0,0,1],[0,1,0],[1,0,0]])
    # first we renormalize rot_mat
    rot_mat = rot_mat/np.linalg.norm(rot_mat,axis=1)
    # NOTE: a weird convention is used in the original paper
    # NOTE: here we transpose traditional rotation matrix.
    rot_mat = rot_mat.T
    # NOTE: Also, need to some axes: xyz -> yzx
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
    \hat{Y_m} = \sum_{m} R_{mm'} Y_m' 
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
