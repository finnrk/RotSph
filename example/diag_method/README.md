## Diagonalization based rotation [BaTiO₃]

Just like in the manual rotation, here we have two BaTiO₃ structures, one with 
c-axis aligned along the z-direction ("_orth_") and the other with the c-axis 
aligned along the [1,-1,1] direction ("_diag_") with Eular angles of 
[45, 0, 54.73] degrees.

Unlike the rotation mathod used in this repo, here, an alternative method based 
on diagonalizing the projected occupation matrix proposed by [Ruchika Mahajan
et. al.](https://doi.org/10.1103/PhysRevMaterials.5.104402) is used to get the 
projected bands.

### Method

The general projected occupation matrix is defined as:

$$
u^{IJ\sigma}\_{mm'} = \sum_{v,\mathbf{k}} f_{v,\mathbf{k}}^{\sigma} \braket{\psi_{v,\mathbf{k}}^{\sigma}|\phi_{m'}^J} \braket{\phi_{m}^I|\psi_{v,\mathbf{k}}^{\sigma}}
$$

where $\sigma$ represents the spin, $I$ and $J$ are the atom indices, $m$ and 
$m'$ are the magnetic quantum numbers, $v$ is the band index, $\mathbf{k}$ is 
the k-point, $\psi$ represent the Bloch states calculated from DFT and $\phi$ 
is the atomic orbital.

The diagonalization method used here is based on the diagonalizing the diagonal
part (i.e., $I=J$) of the general projected occupation matrix: 
$u^{II\sigma}_{mm'} $

The transformation matrix is then constructed from the eigenvectors.

Comparing to the rotation based method, this method has several down-sides when
it comes to chemical analysis.

1. The orbitals are no longer labled due to the process of diagonalization. One
   can only find the corresponding orbitals by plotting out the rotated atomic
   orbitals.
3. The diagonalization based method utilizes the complete degrees of freedom 
   available. For example, for $d$ orbitals, the degrees of freedom is 5. 
   Because of this, some orbtials might be oriented differently from the
   canonical notation. For example, $d_{xz}$ and $d_{yz}$ might be rotated by
   any degrees along $z$ axis while the rest of the $d$ orbitals are still
   in their canonical position.

### Testing 

With the _orth_ structure, we have the following projected bands:

<p align="center">
<img src="https://github.com/Chengcheng-Xiao/RotSph/blob/master/example/diag_method/images/band_orth.png?raw=true" width="50%" height="50%">
</p>

With the _diag_ structure, we have

<p align="center">
<img src="https://github.com/Chengcheng-Xiao/RotSph/blob/master/example/diag_method/images/band_diag.png?raw=true" width="50%" height="50%">
</p>

Using diagonalization based method, we can re-calculate the projected band 
structures using `find_angles_diag.py` (`python rot_pband_plot.py 0 4`):

<p align="center">
<img src="https://github.com/Chengcheng-Xiao/RotSph/blob/master/example/diag_method/images/band_rotsph.png?raw=true" width="50%" height="50%">
</p>

We see that while the projected bands resemble that of the _orth_ structure, the
order is completely changed.
<!-- Comparing our rotated projected band structrue to the one we got for the _orth_ -->
<!-- structure, we realize that $d_{xz}$ and $d_{yz}$ swapped places. But this is fine as -->
<!-- they are symmetry-related. -->

#### Steps to reproduce
0. go to `diag` dir.
1. run a single-point calculation in `scf` folder to get `CHGCAR`.
2. copy `CHGCAR` to `band` folder, change `LORBIT=11` to `LORBIT=12` and run VASP again.
4. copy `find_angles_diag.py` to `band` folder.
5. run `export PYTHONPATH="${PYTHONPATH}:/path/to/repo"`
6. run `python find_angles.py 0 4`. (for the 0-th atom and for $d_{xy}$(m=-2) orbitals)
