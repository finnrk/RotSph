## Optimized Wannier projection [BaTiO₃]

As I've demonstrated in automatic rotation, optimally rotated atomic orbitals
can help finding the best orbtial features of a wavefunction. Apart of helping
analyzing the band structure, this can also help choosing the best initial
projections used for Wannierization (using Wannier90).

Again, we are going to be using BaTiO₃ as an example. Here, we are going to use 
the _diag_ structure where c-axis is aligned along the [1,-1,1] direction with 
Eular angles of [45, 0, 54.73] degrees.

> [!NOTE]
> Here, we adopt the [Tait–Bryan notation](https://en.wikipedia.org/wiki/Euler_angles#Tait–Bryan_angles),
> specifically, [the z-y-x conention (alpha-z, beta-y, gamma-x)](https://en.wikipedia.org/wiki/Rotation_matrix#General_3D_rotations).


### Using un-rotated atomic orbtials

If we plot the projected bands calculaed using the default `LORBIT=11` method,
we get the following:

<p align="center">
<img src="https://github.com/Chengcheng-Xiao/RotSph/blob/master/example/optimized_wannier_projection/images/band_diag.png?raw=true" width="50%" height="50%">
</p>

Now, if we want to wannierze the three bands between 5.8 eV and 8 eV, from these
projected band structures, it is very hard to decide which orbtial to be used to
construct the initial guess.

The best choice we have here is to use $d\_{z^2}$, $d_{xz}$ and $d_{x2-y2}$
orbitals of Ti atom.

The spread of the Wannier functions calculated using these orbtials as initial 
guess are (projection-only Wannier functions):
```
Initial State
 WF centre and spread    1  (  1.857135,  0.975222,  2.884264 )     2.30076325
 WF centre and spread    2  (  1.811067,  1.005568,  2.831399 )     2.07963346
 WF centre and spread    3  (  1.822488,  1.001925,  2.876755 )     2.11022391
 Sum of centres and spreads (  5.490690,  2.982715,  8.592418 )     6.49062062

     0     0.649E+01     0.0000000000        6.4906206224       0.00  <-- CONV
       O_D=      0.1084170 O_OD=      0.1777149 O_TOT=      6.4906206 <-- SPRD
```


### Using optimized atomic orbtials

Since we know the optimized local frame can be constructed using Eular angles of 
[45, 0, 54.73] degrees. The rotated z and x axis can be found using the
following script (`get_axis.py`):
```python
import numpy as np
import sys

z_angle_degree = float(sys.argv[1])
y_angle_degree = float(sys.argv[2])
x_angle_degree = float(sys.argv[3])

z_angle = z_angle_degree/180.*np.pi
rot_z =np.array([[np.cos(z_angle),-np.sin(z_angle),0],\
        [np.sin(z_angle),np.cos(z_angle),0],[0,0,1]])
y_angle = y_angle_degree/180.*np.pi
rot_y=np.array([[np.cos(y_angle),0,np.sin(y_angle)],[0,1,0],\
        [-np.sin(y_angle),0,np.cos(y_angle)]])
x_angle = x_angle_degree/180.*np.pi
rot_x=np.array([[1,0,0],[0,np.cos(x_angle),-np.sin(x_angle)],\
        [0,np.sin(x_angle),np.cos(x_angle)]])

# construct the rotation matrix of coordinate system
rot_mat = np.dot(np.dot(rot_z,rot_y),rot_x)

# get the axis of the new coordinate system
Cart = np.array([[1,0,0],[0,1,0],[0,0,1]])

# new axis needs to be transposed
new_axis = np.dot(rot_mat,Cart).T

print(new_axis)
```
by `python get_axis.py 45 0 54.75`. 

If we plot the projected bands calculaed using the optimized atomic orbitals,
we get the following:

<p align="center">
<img src="https://github.com/Chengcheng-Xiao/RotSph/blob/master/example/optimized_wannier_projection/images/band_rotsph.png?raw=true" width="50%" height="50%">
</p>

Again, we want to wannierze the three bands between 5.8 eV and 8 eV and from 
these projected band structures, we see very clear that these bands are compoesd
from $d_{xy}$, $d_{yz}$ and $d_{xz}$ orbitals of Ti atom.

The spread of the Wannier functions calculated using these (rotated) orbtials as 
initial guess are (projection-only Wannier functions):
```
 Initial State
  WF centre and spread    1  (  1.784574,  1.039964,  2.818440 )     1.91953326
  WF centre and spread    2  (  1.853310,  0.971228,  2.887160 )     2.15364331
  WF centre and spread    3  (  1.853301,  0.971237,  2.887142 )     2.15353038
  Sum of centres and spreads (  5.491185,  2.982429,  8.592742 )     6.22670695

      0     0.623E+01     0.0000000000        6.2267069529       0.00  <-- CONV
        O_D=      0.0088859 O_OD=      0.0133322 O_TOT=      6.2267070 <-- SPRD
```

And we see that the spread is better than the one we get from before. 

I feel that this is not the typical use case as the difference is not prominent
and the band structrues and the imag/real ration of the Wannier functions are
pretty similar. I suspect a greater difference would appear for systems in which
a smooth Wannier gauge is harder to find and only when we use a good initial
projection will the calculation make sense. 

### Steps to reproduce
0. go to `diag` dir.
1. go to `wannier_xy_yz_xz` dir and run VASP.
2. run `wannier90.x wannier90`.
3. rinse and repeat for `wannier_z2_xz_x2` dir.
