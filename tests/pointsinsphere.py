
#!/usr/bin/env python

import lattice
import numpy as np
from pymatgen.io.ase import AseAtomsAdaptor

lat = lattice.Lattice([[9, 0, 2.3], [0, -12, 7.2], [-4, 1, 8]])
# lat =  lattice.Lattice([[9, 0, 0],[0, -12, 0],[0, 0, 8]])

coords = [1, 1, 1], [3, 4, -5.6], [1, 1, 1.2], [-3.4, 2.2, -9.8]
cart_coords = lat.get_cartesian_coords(coords)
print(cart_coords)

result = lat.get_points_in_sphere(coords, [9, -3.1, 0.8], r=9.3)
dist = [dist for (fcoords, dist, i, img) in result]
fcoords = [fcoords for (fcoords, dist, i, img) in result]
print()
print("OUT:")
print(dist)
print(fcoords)
