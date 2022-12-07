# Shift gridpoints file

import numpy as np
import ase, ase.io
import argparse
from pathlib import Path

# Have not been tested:
def coords_within_cell_boundaries(coords, atoms, cart_coord=False):
    """Checks whether the coordinates are inside a unitcell."""

    cell = atoms.get_cell()
    a_cell, b_cell, c_cell = cell[0,:], cell[1,:], cell[2,:]

    if cart_coord:
        frac_coords = np.linalg.solve(cell.T, coords.T).T
    else:
        frac_coords = coords
    
    is_in_cell =  ((0.0 <= frac_coords[:,0] <= 1.0).all() and (0.0 <= frac_coords[2] <= 1.0).all() and (0.0 <= frac_coords[2] <= 1.0).all())

    return is_in_cell


parser = argparse.ArgumentParser()
parser.add_argument("coords", action="store", type=str, help="Text file containing coordinates of the grid.")
parser.add_argument("atoms", action="store", type=str, help="File containing the unitcell.")

args = parser.parse_args()

coords_file = args.coords
atoms_file = args.atoms

# Read coordinates and atoms object
coords = np.loadtxt(coords_file)
atoms = ase.io.read(atoms_file)

# shift coordinates
a,b,c = atoms.get_cell()[0,:], atoms.get_cell()[1,:], atoms.get_cell()[2,:]
scale = 0.00001
coords_shifted = coords + scale * (a + b + c)

np.savetxt(coords_file + "_shifted.txt", coords_shifted)

