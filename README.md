# Rotation Matrix of Real Spherical Harmonics
This repo contains a script to calculate the transformation matrix for real 
spherical harmonics ($Y_{lm}$) under different basis frame using the 
[real-space rotation matrix](https://en.wikipedia.org/wiki/Rotation_matrix) of 
the basis set.

I.e., the rotated real spherical harmonics ($\tilde Y_{lm}$) can be expressed 
as a linear combination of the original real spherical harmonics ($Y_{lm}$),

$$
\tilde Y_{lm'} = \sum_{m} R_{m'm} Y_{lm}.
$$

Where $R$ is the transformation matrix.

The rotation matrix ($r$) that rotates the Cartesian axes ($e$) is defined zs:

$$
\tilde e_{i} = \sum_j r_{ij} e_j,
$$

where $\tilde e$ is the Cartesian axes after rotation.

For technical details, see: [J. Phys. Chem. 1996, 100, 15, 6342–634](https://pubs.acs.org/doi/10.1021/jp953350u)

## Requirement
- numpy

## Usage

```python
from rotsph import rotsph
import numpy as np

# basis rotation matrix
r=np.matrix([[1,0,1],[0,1,0],[-1,0,1]])

# for f orbtials, we use "3".
R = rotsph.get_R_mat(3,r)

print(R)
```

## Result Checker
See `plot_rotated_orbitals.nb` (Mathematica Notebook).

## Use Case
In normal DFT codes, orbital projections are usually done with real atomic 
orbitals that are oriented on the canonical Cartesian axes. However, this may 
not be compatible with the crystal symmetries and may render difficulties in 
bonding analysis.

The transfomration matrix generated here can be used to construct optimally 
oriented local orbitals and can be applied to the canonical projection 
coefficients to get the projection coefficients on the rotated frame. 
I.e., get projection coefficients obtained by projecting wavefunctions on to 
atomic orbitals on the rotated frame.

This feature has been implemented (see example/manual_rotation) and an routine
that performes automatic search of the best Eular angles is also implemented
(see example/automatic_rotation).

## License
RotSph is released under the MIT license.

