#!/usr/bin/python3 -> check this thing

# THIS IS UNFINISHED: TODO: Properly structure the argument parser. Add read section for the input file to get atoms obj. Then add section to compute spacing
import numpy as np
import ase
import ase.io
import argparse
from pathlib import Path

parser = argparse.ArgumentParser()

parser.add_argument("file", help="File from which to read atoms object and optionally the grid.")
parser.add_argument("-g", type=int, nargs="*", default=None, help="The grid discretization (i.e. number of points of the grid in each dimension) of the unit cell in atoms object from file, for which the spacing should be computed.")


args = parser.parse_args()

stem = Path()

# To compute spacing, starting point:
        a,b,c = np.linalg.norm(lgps.get_cell()[0,:]), np.linalg.norm(lgps.get_cell()[1,:]), np.linalg.norm(lgps.get_cell()[2,:])
        nx,ny,nz = ma.ceil(a/spacing_x), ma.ceil(b/spacing_y), ma.ceil(c/spacing_z)
        true_spacing = (a/nx, b/ny, c/nz)
        print('Desired spacing: ', (spacing_x, spacing_y, spacing_z),' True spacing: ', true_spacing)
