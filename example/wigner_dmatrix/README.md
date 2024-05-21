## Copmaring to Wigner's D-matrix

Wigner's D-matrix rotates a set of complex (pure) spherical harmonics to a
desired local frame:

$$
\tilde{Y}\_{l}^{m'} = \sum_{m} D_{m'm} Y\_{l}^{m} 
$$

Since we know the transformation matrix $U$ to go from [the complex spherical 
harmonics to the real ones](https://en.wikipedia.org/wiki/Table_of_spherical_harmonics#Real_spherical_harmonics).
we can write:

$$
\begin{aligned}
\tilde{Y}\_{l}^{n'} &= \sum_{n} D^{n'}\_{n} Y_{l}^{n}\\
&= \sum_{n} D^{n'}\_{n} \left( \sum_{b }\mathbb{I}^{n}\_{b} Y_{l}^{b} \right)\\
&= \sum_{n} D^{n'}\_{n} \left[ \sum_{b} \left( \sum_{c}  U_{c}^{n}U_{b}^{c\dagger} \right) Y_{l}^{b} \right]\\
&= \sum_{n} D^{n'}\_{n}  \sum_{c}  U_{c}^{n} \sum_{b} U_{b}^{c\dagger} Y_{l}^{b} \\
\tilde{Y}\_{l}^{n'} &= \sum_{n} D^{n'}\_{n}  \sum_{c}  U_{c}^{n} Y_{l,c} \\
\sum_{a} \tilde{Y}\_{l}^{a} \mathbb{I}^{n'}\_{a} &= \sum_{n} D^{n'}\_{n}  \sum_{c}  U_{c}^{n} Y_{l,c} \\
\sum_{a} \tilde{Y}\_{l}^{a} \left( \sum_{d}  U_{a}^{d\dagger}U_{d}^{n'} \right)  &= \sum_{n} D^{n'}\_{n}  \sum_{c}  U_{c}^{n} Y_{l,c} \\
\sum_{d} \tilde{Y}\_{l,d} U_{d}^{n'}  &= \sum_{n} D^{n'}\_{n}  \sum_{c}  U_{c}^{n} Y_{l,c} \\
\sum_{n'} U^{d\dagger}\_{n'}\sum_{d} \tilde{Y}\_{l,d} U_{d}^{n'}  &= \sum_{n'} U^{d\dagger}\_{n'}\sum_{n} D^{n'}\_{n}  \sum_{c}  U_{c}^{n} Y_{l,c} \\\sum_{n'} \tilde{Y}_{l,d} &= \sum_{n'} U^{d\dagger}_{n'}\sum_{n} D^{n'}_{n}  \sum_{c}  U_{c}^{n} Y_{l,c} 
\end{aligned}
$$

Hence, the transformation matrix is:

$$
\mathbf{T} = \mathbf{U^{\dagger}DU}.
$$

This has been implemented in [pyscf](https://pyscf.org/_modules/pyscf/symm/Dmatrix.html).
We can use it to verify that our implementation here is correct.

> [!NOTE]
> pyscf adopts a [z-y-z convention of the Euler angles](https://en.wikipedia.org/wiki/Euler_angles#Chained_rotations_equivalence).
> This means in pyscf, alpha rotates by z axis, beta rotates by the new y axis 
> and gamma rotates by the final z axis. However, we adopt the 
> [Tait–Bryan notation](https://en.wikipedia.org/wiki/Euler_angles#Tait–Bryan_angles),
> specifically, [the z-y-x conention (alpha-z, beta-y, gamma-x)](https://en.wikipedia.org/wiki/Rotation_matrix#General_3D_rotations).


With pyscf, we test the Euler angles of (-45, 54.74, 90):
```python
from pyscf.symm import Dmatrix

np.set_printoptions(suppress=True, precision=4)

print(Dmatrix.Dmatrix(3,-45/180*np.pi,54.74/180*np.pi,90/180*np.pi).T)
```
we get output of:
```python
array([[-0.3402,  0.6667, -0.2635, -0.4304,  0.2635,  0.    , -0.3402],
       [ 0.4082,  0.    ,  0.527 , -0.    ,  0.527 ,  0.3335, -0.4082],
       [ 0.2635, -0.0001, -0.6124,  0.3332,  0.6124,  0.    ,  0.2635],
       [-0.3043, -0.7454, -0.2356, -0.385 ,  0.2356,  0.    , -0.3043],
       [-0.4565,  0.    ,  0.1178, -0.    ,  0.1178,  0.7453,  0.4565],
       [ 0.4714, -0.0001,  0.0001, -0.7454, -0.0001, -0.    ,  0.4714],
       [ 0.3535, -0.    , -0.4565,  0.    , -0.4565,  0.5773, -0.3535]])
```

With rotsph, the equivalent Euler angles is (45,0,54.74):
```python
z_angle = 45/180*np.pi
rot_z =np.array([[np.cos(z_angle),-np.sin(z_angle),0],\
        [np.sin(z_angle),np.cos(z_angle),0],[0,0,1]])
y_angle = 0/180*np.pi
rot_y=np.array([[np.cos(y_angle),0,np.sin(y_angle)],[0,1,0],\
        [-np.sin(y_angle),0,np.cos(y_angle)]])
x_angle = 54.74/180*np.pi
rot_x=np.array([[1,0,0],[0,np.cos(x_angle),-np.sin(x_angle)],\
        [0,np.sin(x_angle),np.cos(x_angle)]])

# construct the rotation matrix of coordinate system
rot_mat = np.dot(np.dot(rot_z,rot_y),rot_x)

# get the transformation matrix for d orbitals
T_mat = rotsph.get_R_mat(2,rot_mat)

np.set_printoptions(suppress=True, precision=4)
print(T_mat)
```
and get:
```
array([[-0.3402,  0.6667, -0.2635, -0.4304,  0.2635,  0.    , -0.3402],
       [ 0.4082, -0.    ,  0.527 ,  0.    ,  0.527 ,  0.3335, -0.4082],
       [ 0.2635, -0.0001, -0.6124,  0.3332,  0.6124, -0.    ,  0.2635],
       [-0.3043, -0.7454, -0.2356, -0.385 ,  0.2356, -0.    , -0.3043],
       [-0.4565, -0.    ,  0.1178,  0.    ,  0.1178,  0.7453,  0.4565],
       [ 0.4714, -0.0001,  0.0001, -0.7454, -0.0001,  0.    ,  0.4714],
       [ 0.3535, -0.    , -0.4565,  0.    , -0.4565,  0.5773, -0.3535]])
```
