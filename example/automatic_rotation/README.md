## Automatic rotation [BaTiO₃]

Just like in the manual rotation, here we have two BaTiO₃ structures, one with 
c-axis aligned along the z-direction ("_orth_") and the other with the c-axis 
aligned along the [1,-1,1] direction ("_diag_") with Eular angles of 
[45, 0, 54.73] degrees.

Whilst we know the actual rotation matrix used to obtain the _diag_ structure,
here, I'll show an way that can automatically search for the best Eular angles.

### Method

The method used here is based on the observation that the best projector 
(atomic orbitals) will yield the highest projection coeffs (squared).

In other workds, this means we want to find the Eular angles that gives the
highest total projection coeffs (squared) of an $l$-shell for a set of
wavefunctions at $\Gamma$ point (in theory this can also be extended to be
k-dependent, but for the moment we shall stick to only use $\Gamma$ point).

The reason I only focused on a single $l$-shell is the fact that the
orthogonality is retained via the radial part of the wavefunctions instead of
the angular part. And using a different rotation frame for different $l$-shells
should also be valid.

Whilst in theory we could find the best Eular angles by performing an 
analytical analysis on the angle dependence of the projection coeffs, due to 
the normalization factor difference (e.g. https://www.vasp.at/forum/viewtopic.php?t=18618)
I'll stick with an exhaust search of all the possible angles.

Some improvement might also be done with symmetry of the real-spherical
harmonics and the local environment, but for the moment, I'll ignore these.

### Testing 

With the _orth_ structure, we have the following projected bands:

<p align="center">
<img src="https://github.com/Chengcheng-Xiao/RotSph/blob/master/example/automatic_rotation/images/band_orth.png?raw=true" width="50%" height="50%">
</p>

With the _diag_ structure, we have

<p align="center">
<img src="https://github.com/Chengcheng-Xiao/RotSph/blob/master/example/automatic_rotation/images/band_diag.png?raw=true" width="50%" height="50%">
</p>

Now, using `find_angles.py` file (`python find_angles.py 0 d`), we get the following Eular angles:
```
   Eular_alpha:  135.50 degrees
   Eular_beta :  -54.75 degrees
   Eular_gamma:  -180.50 degrees
```

And with these angles, we can re-calculate the projected band structures using
`rot_pband_plot.py` (`python rot_pband_plot.py 0 4 135.5 -54.75 -180.50`). And we get:

<p align="center">
<img src="https://github.com/Chengcheng-Xiao/RotSph/blob/master/example/automatic_rotation/images/band_rotsph.png?raw=true" width="50%" height="50%">
</p>

Comparing our rotated projected band structrue to the one we got for the _orth_
structure, we realize that $d_{xz}$ and $d_{yz}$ swapped places. But this is fine as
they are symmetry-related.

#### Steps to reproduce
0. go to `diag` dir.
1. run a single-point calculation in `scf` folder to get `CHGCAR`.
2. copy `CHGCAR` to `band` folder, change `LORBIT=11` to `LORBIT=12` and run VASP again.
4. copy `find_angles.py` and `rot_pband_plot.py` to `band` folder.
5. run `export PYTHONPATH="${PYTHONPATH}:/path/to/repo"`
6. run `python find_angles.py 0 d`. (for the 0-th atom and for d orbitals)
7. run `python rot_pband_plot.py 0 4 135.5 -54.75 -180.50` (for the 0-th atom and for dxy orbital)


