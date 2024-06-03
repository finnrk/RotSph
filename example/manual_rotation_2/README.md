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

