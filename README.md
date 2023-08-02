# Rotation Matrix of Real Spherical Harmonics
This repo contains a script to calculate the rotation matrix to real spherical harmonics ($Y_{lm}$) under the [real-space rotation matrix](https://en.wikipedia.org/wiki/Rotation_matrix)($R$).

I.e. the rotated real spherical harmonics ($\tilde Y_{lm}$) can be expressed as a linear combination of the original real spherical harmonics ($Y_{lm}$),

$$
\tilde Y_{lm} = \sum_{m} Y_{lm} R_{mm'}.
$$

The real-space rotation matrix ($r$) that rotates the Cartesian axes ($e$) is defined through:

$$
\tilde e_{i} = \sum_j r_{ij} e_j,
$$

where $\tilde e$ is the Cartesian axes after rotation.

For details, see: [J. Phys. Chem. 1996, 100, 15, 6342–634](https://pubs.acs.org/doi/10.1021/jp953350u)

## Requirement
- numpy

## Usage

```python
from rotsph import rotsph
import numpy as np

r=np.matrix([[1,0,1],[0,1,0],[-1,0,1]])
R = rotsph.get_R_mat(3,r)

print(R)
```

## Result Checker
See `plot_rotated_orbitals.nb`.

## Use Case
In normal DFT codes, orbital projections are usually done with real atomic orbitals that are oriented on the canonical Cartesian axes. 
However, this may not be compatible with the crystal symmetries and may render difficulties in bonding analysis.

This rotation matrix can be used to construct optimally oriented local orbitals and can be applied to the canonical projection coefficients to get the projection coefficients on the rotated frame. I.e., get projection coefficients obtained by projecting wavefunctions on to atomic orbitals on the rotated frame.

Example: coming soon.

## Disclaimer
The correctness of this code has not been thoroughly tested.
Make sure you KNOW WHAT YOU ARE DOING!
