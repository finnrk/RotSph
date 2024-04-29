## BaTiO₃

### Problem

Here I demonstrate the effect of canonical porjection using using two BaTiO₃
(BTO) sturctures, one with c-axis aligned along the z-direction ("orth") and 
one with the c-axis aligned along the [1,1,1] direction ("diag").

Using canonical projection, for the orth structure, we obtain the following
projected band structure:

<img src="https://github.com/Chengcheng-Xiao/RotSph/blob/master/example/BTO/images/band_orth.png?raw=true" width="50%" height="50%">

And if we plot the wavefunction of #25 band at $\Gamma$ point, we see that it
indeed is a $d_{z^2}$ orbtial oriented along the z-direction:

<img src="https://github.com/Chengcheng-Xiao/RotSph/blob/master/example/BTO/images/wfc_25_gamma.png?raw=true" width="50%" height="50%">

Now if we run the same calculation but with the diag structure, we see
completely different projected band structure:

<img src="https://github.com/Chengcheng-Xiao/RotSph/blob/master/example/BTO/images/band_diag.png?raw=true" width="50%" height="50%">

But if we plot the wavefunction of #25 band at $\Gamma$ point, the $d_{z^2}$
feature is still there.

#### Steps to reproduce
0. go to `orth` dir.
1. run a single-point calculation in `scf` folder to get `CHGCAR`.
2. copy `CHGCAR` to `band` folder and run VASP again.
3. plot band structure using [pyband](https://github.com/QijingZheng/pyband) 
with `pyband -y -2 12 -z 0 --occ '1' --spd '6'`.
4. Rinse and repeat for `diag` dir.

### Solution

We can use the transformation matrix generated with `rotsph` code to obtain the
optimally rotated projected bandstructure. For this to work, we need to run VASP
with `IORBIT=12` (so that the actual projection coeffs along with its phase 
factors are printed out).

<img src="https://github.com/Chengcheng-Xiao/RotSph/blob/master/example/BTO/images/band_rotsph.png?raw=true" width="50%" height="50%">

#### Steps to reproduce
0. go to `orth` dir.
1. run a single-point calculation in `scf` folder to get `CHGCAR`.
2. copy `CHGCAR` to `band` folder and run VASP again.
3. copy `rot_pband.py` to `band` folder.
4. run `export PYTHONPATH="${PYTHONPATH}:/path/to/repo"`     
5. run `python rot_pband.py`.

*Please take a look into the `rot_pband.py` file and make sure you understand
what is going on.*



