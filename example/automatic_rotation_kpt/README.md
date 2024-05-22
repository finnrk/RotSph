## Automatic rotation [TiO₂]

For TiO₂, we use a sructure where the bonds of the corner Ti atoms are aligned
along the x, y and z axis (Cartesian). With this, we can plot the $e_g$ and
$t_{2g}$ orbtials splitting of the corner atoms as:

<p align="center">
<img src="https://github.com/Chengcheng-Xiao/RotSph/blob/master/example/automatic_rotation_kpt/images/dos_corner.png?raw=true" width="50%" height="50%">
</p>

which is almost identical to Fig. 2 of [this paper](https://journals.aps.org/prb/pdf/10.1103/PhysRevB.95.195159).

However, if we want to plot the same for the centere Ti atom, rotation is
needed. If we use the old `find_angles.py` which only finds the best angle for
states at $\Gamma$ point, we get weird plots which doesn't resemble that of the
$e_g$ and $t_{2g}$ orbitals:

<p align="center">
<img src="https://github.com/Chengcheng-Xiao/RotSph/blob/master/example/automatic_rotation_kpt/images/dos_wrong.png?raw=true" width="50%" height="50%">
</p>

This is because $d_{z^2}$ (and other orbitals) are indeed not oriented along the 
bonds:

<p align="center">
<img src="https://github.com/Chengcheng-Xiao/RotSph/blob/master/example/automatic_rotation_kpt/images/wfc.png?raw=true" width="50%" height="50%">
</p>

To help analyze similar systems, I've updated the `find_angles.py` to
`find_angles_kpt.py` so that it is now trying to find the best angle across all 
k-points. With this, we get the following plots:

<p align="center">
<img src="https://github.com/Chengcheng-Xiao/RotSph/blob/master/example/automatic_rotation_kpt/images/dos_rotsph.png?raw=true" width="50%" height="50%">
</p>

#### Steps to reproduce
0. go to `TiO2` dir and run a single-point calculation.
5. run `export PYTHONPATH="${PYTHONPATH}:/path/to/repo"`
6. run `python find_angles.py 1 d`. (for the second atom and for d orbitals)
7. run `python rot_pdos_plot.py 1 4 -143.75 -90 -81.25 -3 13 100 0.1` (for the second atom and for dxy orbital)
8. run `python find_angles_kpt.py 1 d`. (for the second atom and for d orbitals)
9. run `python rot_pdos_plot.py 1 4 -98 91 -143 -3 13 100 0.1` (for the second atom and for dxy orbital)
