## BaTiO₃

### Problem

Here I demonstrate the effect of canonical porjection using using two BaTiO₃
(BTO) sturctures, one with c-axis aligned along the z-direction ("_orth_") and
one with the c-axis aligned along the [1,1,1] direction ("_diag_").

Using canonical projection, for the orth structure, we obtain the following
projected band structure:

<p align="center">
<img src="https://github.com/Chengcheng-Xiao/RotSph/blob/master/example/BTO/images/band_orth.png?raw=true" width="50%" height="50%">
</p>

And if we plot the wavefunction of #25 band at $\Gamma$ point, we see that it
indeed is a $d_{z^2}$ orbital oriented along the z-direction:

<p align="center">
<img src="https://github.com/Chengcheng-Xiao/RotSph/blob/master/example/BTO/images/wfc_25_gamma.png?raw=true" width="50%" height="50%">
</p>

Now if we run the same calculation but with the _diag_ structure, we get a
completely different projected band structure:

<p align="center">
<img src="https://github.com/Chengcheng-Xiao/RotSph/blob/master/example/BTO/images/band_diag.png?raw=true" width="50%" height="50%">
</p>

But if we plot the wavefunction of #25 band at $\Gamma$ point, the $d_{z^2}$
feature is still there. 

The difference in these two projected band structures is caused by the use of 
atomic orbital aligned along the canonical Cartesian axes that are not rotated 
along with the crystal axes.

#### Steps to reproduce
0. go to `orth` dir.
1. run a single-point calculation in `scf` folder to get `CHGCAR`.
2. copy `CHGCAR` to `band` folder and run VASP again.
3. plot band structure using [pyband](https://github.com/QijingZheng/pyband)
with `pyband -y -2 12 -z 0 --occ '1' --spd '6'`.
4. Rinse and repeat for `diag` dir.

### Solution

We can use the transformation matrix generated with `RotSph` code to obtain the
optimally rotated projected bandstructure. For this to work, we need to run VASP
with `LORBIT=12` so that the projection coeffs along with its phase
factors are printed out (not only the square of the projection coeffs).

Once we have the projection coeffs (e.g. to atomic orbital $j$) for band $i$ at
$k$-point ($c_{ik,j}=\braket{\phi_{ik}|p_{j}}$), we can calculate the
rotated projection coeffs ($\tilde{c}\_{ik,j'}=\braket{\phi_{ik}|p_{j'}}$)
using the transformation matrix. Specifically:

$$
\begin{aligned}
\tilde{c}\_{ik,j'} &= \braket{\phi_{ik}|(\hat{R}|p)\_{j'}} \\
&= \int \phi_{ik}(\mathbf{r}) \sum_{j} R_{j'j} p_{j}(\mathbf{r}) d\mathbf{r} \\
&= \sum_{j} R_{j'j} \braket{\phi_{ik}|p_{j}} \\
&= \sum_{j} R_{j'j} c_{ik,j},
\end{aligned}
$$

and for the projected band structure, we need its square $|\tilde{c}_{ik,j'}|^2$.

The resulting projected band sturctures calculated using the _diag_ structure now
looks exactly like the one calculated using the _orth_ structure:

<p align="center">
<img src="https://github.com/Chengcheng-Xiao/RotSph/blob/master/example/BTO/images/band_rotsph.png?raw=true" width="50%" height="50%">
</p>

#### Steps to reproduce
0. go to `diag` dir.
1. run a single-point calculation in `scf` folder to get `CHGCAR`.
2. copy `CHGCAR` to `band` folder, change `LORBIT=11` to `LORBIT=12` and run VASP again.
4. copy `rot_pband.py` to `band` folder.
5. run `export PYTHONPATH="${PYTHONPATH}:/path/to/repo"`
6. run `python rot_pband.py`.

*`rot_pband.py` reads in the projection coeffs* $c\_{ik,j}$ and calculates the 
rotated $\tilde{c}\_{ik,j'}$ *and then plot the projected band structures.*
