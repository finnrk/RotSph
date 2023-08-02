# Rotation Matrix of Real Spherical Harmonics
This repo contains a script to calculate the rotation materix to real spherical harmonics ($Y_{lm}$) under the [real-space rotation matrix](https://en.wikipedia.org/wiki/Rotation_matrix)($R$).

I.e. the rotated real spherical harmonics is,
$$
\tilde{Y}_{lm'} = \sum_{m} Y_{lm} R_{mm'}
$$

The real-space rotation matrix ($r$) that rotates the Cartesian axes ($\bm{e}$) is defined through:

$$
\bm{\tilde{e}}_i =\sum_j r_ij \bm{e}_j.
$$

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

## Disclaimer
The correctness of this code has not been thoroughly tested.
Make sure you KNOW WHAT YOU ARE DOING!
