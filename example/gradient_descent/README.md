## Automatic rotation [BaTiO₃]

Similarly to in automatic\_rotation, we have two BaTio3 structures. One has the Ti-O pairs aligned with the x, y, z axes ("_aligned_"), 
while the other has has its lattice vectors rotated by a randomly generated matrix ("_rand_").

> [!NOTE]
> Here, we adopt the [Tait–Bryan notation](https://en.wikipedia.org/wiki/Euler_angles#Tait–Bryan_angles),
> specifically, [the z-y-x convention (alpha-z, beta-y, gamma-x)](https://en.wikipedia.org/wiki/Rotation_matrix#General_3D_rotations).

Whilst we know the actual rotation matrix used to obtain the _aligned_ structure from the _rand_ structure,
here, I'll show a way that can automatically search for the best Euler angles.

### Method

We define a metric to quantify the degree to which a band is projected onto multiple orbitals.
The metric is smallest when each band is only projected onto a single orbital and is largest when each band is projected equally onto all orbitals. We aim to minimise this metric.

We construct this metric by taking a function $f(x)$, and summing together $f(x)$ for the square projections $x$ of each k-point, band and $l$ onto each orbital.
We desire a function $f$ such that $f(a+b) < f(a) + f(b)$, (satified by any concave function). I used $f(x) = \sqrt{x}$.

We then evaluate the metric at randomly sampled orientations to find one in 
the viscinity of the optimal orientation, before performing a gradient descent to find the orientation to the desired accuracy.


### Testing 

With the _aligned_ structure, we have the following projected bands in the order $d_{xy}$, $d_{yz}$, $d_{z^2}$, $d_{xz}$, $d_{x^2-y^2}$:

<p align="center">
<img src="https://github.com/finnrk/RotSph/blob/master/example/gradient_descent/images/bands_aligned.png?raw=true" width="80%" height="80%">
</p>

With the _rand_ structure, we have:

<p align="center">
<img src="https://github.com/finnrk/RotSph/blob/master/example/gradient_descent/images/bands_unaligned.png?raw=true" width="80%" height="80%">
</p>

Now, using `rotate.py` file (`python rotate.py PROCAR_UNALIGN 1`), we get the following Euler angles:
```
alpha, beta, gamma:  -56.11458181966653 15.448190216500603 155.32619084727617
```

And with these angles, we can re-calculate the projected band structures using
`rot_pband_plot.py` (`python pbandplot.py 1 4 -56.11458181966653 15.448190216500603 155.32619084727617`). And we get:

<p align="center">
<img src="https://github.com/finnrk/RotSph/blob/master/example/gradient_descent/images/bands_realigned.png?raw=true" width="80%" height="80%">
</p>

We see that the splitting between the $e_g$ and $t_{2g}$ orbitals is recovered.
<!-- Comparing our rotated projected band structrue to the one we got for the _rand_ -->
<!-- structure, we realize that $d_{xz}$ and $d_{yz}$ swapped places, But this is fine as -->
<!-- they are symmetry-related. Similarly, x2y2 and z2 swapped. These are related, since they are the two orbitals which make up the eg orbitals-->

#### Steps to reproduce
1. go to `batio3/rand` dir.
2. run a single-point calculation in `scf` folder to get `CHGCAR`.
3. copy `CHGCAR` to `band` folder and run VASP again.
4. copy `rotate.py` and `pbandplot.py` to `band` folder.
5. run `export PYTHONPATH="${PYTHONPATH}:/path/to/repo"`
6. run `python rotate.py PROCAR 1`. (for the 1st atom)
7. run `python pbandplot.py 1 4 -56.11458181966653 15.448190216500603 155.32619084727617` (for the 1st atom and for dxy orbital)

## Creating a rotated PROCAR file [Hexaquacopper]

Similarly to in the previous section, we have two Hexaquacopper structures. One has the elongated [Jahn-Teller](https://en.wikipedia.org/wiki/Jahn%E2%80%93Teller_effect)
ligand along the z axis, with the other two ligands aligned with the x and y axes ("_aligned_"), 
while the other has had its ions rotated by a randomly generated matrix ("_rand_").

To have `rotate.py` create a `PROCAR` file for the rotated frame, include the name of the output file as an extra parameter in the command line.

### Testing
With the _aligned_ structure, we have the following projected density of states:

<p align="center">
<img src="https://github.com/finnrk/RotSph/blob/master/example/gradient_descent/images/dos_aligned.png?raw=true" width="80%" height="80%">
</p>

With the _rand_ structure, we have:

<p align="center">
<img src="https://github.com/finnrk/RotSph/blob/master/example/gradient_descent/images/dos_unaligned.png?raw=true" width="80%" height="80%">
</p>

Now, using `rotate.py` file (`python rotate.py PROCAR -1 PROCAR_OUT`), we generate a new `PROCAR` file.

With this new file, we can replot the projected density of states using
`pdosplot.py` (`python pdosplot.py PROCAR_OUT -1 -1.5 0.5 600 0.001`). And we get:

<p align="center">
<img src="https://github.com/finnrk/RotSph/blob/master/example/gradient_descent/images/dos_realigned.png?raw=true" width="80%" height="80%">
</p>

#### Steps to reproduce
1. go to `molecule/rand` dir.
2. run VASP to get `PROCAR`
3. run `export PYTHONPATH="${PYTHONPATH}:/path/to/repo"`
4. run `python rotate.py PROCAR -1 PROCAR_OUT` (for the last atom)
5. run `python pdosplot.py PROCAR_OUT -1 -1.5 0.5 600 0.001`

The command line arguments of `pdosplot.py` are:
```
Last atom | dos E min | dos E max | dos points | smearing Sigma|
-1        | -1.5      | 0.5       | 600        | 0.001         |
```
Note that for this particular example we require a very small smearing width so that the dz2 and dx2-y2 orbitals can be distinguished at the fermi energy.

