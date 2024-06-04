## Manual mode [CrI₃]
Here, we try to reproduce the results reported in https://doi.org/10.1016/j.commatsci.2020.109820.

1. Run VASP.
2. Calculate local a axis (Cr[1]-I[5]) and c axis (Cr[1]-I[2]):
```
c = Cr[1]-I[5]  1.23734 -1.87958 1.5594
a = Cr[1]=I[2] -2.24445 -0.13639 1.57024
```

<p align="center">
<img src="https://github.com/Chengcheng-Xiao/RotSph/blob/master/example/manual_rotation_2/images/structure.png?raw=true" width="25%" height="25%">
</p>


3. Calculate Euler angle:
```
python construct_euler.py 1.23734 -1.87958 1.5594  -2.24445 -0.13639 1.57024 T
```
4. run `export PYTHONPATH="${PYTHONPATH}:/path/to/repo"`
5. Calculate pDOS:
```
for i in  4 5 6 7 8
do
python rot_pdos_plot_pymatgen.py 0 $i -176.04162059  -35.31274179  -45.73258177 -9 2 100 0.1 -1 && imgcat DOS.png && mv pdos.dat $i-dn.dat
done
```
6. plot everything:
```
python plot.py
```

With this, we get the $e_g$ DOS:

<p align="center">
<img src="https://github.com/Chengcheng-Xiao/RotSph/blob/master/example/manual_rotation_2/images/DOS_eg.png?raw=true" width="25%" height="25%">
</p>

and $t_{2g}$ DOS:

<p align="center">
<img src="https://github.com/Chengcheng-Xiao/RotSph/blob/master/example/manual_rotation_2/images/DOS_t2g.png?raw=true" width="25%" height="25%">
</p>

However, we see that these plots are still slightly different from Fig. 4 in 
10.1016/j.commatsci.2020.109820. This is because here the z-axis is not aligned
along the long axis of the octahedron

### Using the Euler angle finder for octahedrons
Alternatively, one can use the `eg_t2g_Euler_angle.py`:
```
python eg_t2g_Euler_angle.py 0
                             ^ atom index
```
and it will spit out the proper Euler angles, in this case
```
[-117.12133235  -35.10745263   45.58423345]
```
which properly alignes the z axis along the long axis of the Cr-I octahedron,
and gives us:

<p align="center">
<img src="https://github.com/Chengcheng-Xiao/RotSph/blob/master/example/manual_rotation_2/images/DOS_eg_1.png?raw=true" width="25%" height="25%">
</p>

and $t_{2g}$ DOS:

<p align="center">
<img src="https://github.com/Chengcheng-Xiao/RotSph/blob/master/example/manual_rotation_2/images/DOS_t2g_1.png?raw=true" width="25%" height="25%">
</p>

This time, the plots looks exactly like the ones in Fig.4 of 
10.1016/j.commatsci.2020.109820.
