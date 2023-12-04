#!/usr/bin/env python

import numpy as np
import ase, ase.io
from ase.io.cube import read_cube, write_cube
from ase.io.cif import read_cif
from ase.io.xsf import read_xsf
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("grid1", action="store", type=str, help="File containing first grid to be compared.")
parser.add_argument("grid2", action="store", type=str, help="File containing second grid to be compared.")
#parser.add_argument("-f", "--format", action="store", nargs="+", options=["cube", "cif", "xsf"], default=None, type=str, help="File containing second grid to be compared.")


args = parser.parse_args()

if args.grid1[-3:] == "cif":
     atoms1, grid1 = read_cif(args.grid1, read_data=True)
elif args.grid1[-3:] == "xsf":
     atoms1, grid1 = read_xsf(args.grid1, read_data=True)
else:
     with open(args.grid1, 'r') as g1:
         grid1_content = read_cube(g1, read_data=True)
         atoms1 = grid1_content["atoms"]
         grid1 = grid1_content["data"]


if args.grid2[-3:] == "cif":
     atoms2, grid2 = read_cif(args.grid2, read_data=True)
elif args.grid2[-3:] == "xsf":
     atoms2, grid2 = read_xsf(args.grid2, read_data=True)
else:
     with open(args.grid2, 'r') as g2:
         grid2_content = read_cube(g2, read_data=True)
         atoms2 = grid2_content["atoms"]
         grid2 = grid2_content["data"]
    

diff_grid = grid1 - grid2
abs_diff_grid = np.abs(diff_grid)

with open("diff_grid.cube", 'w') as cube1:
    write_cube(cube1, atoms1, data=diff_grid)

with open("abs_diff_grid.cube", 'w') as cube2:
    write_cube(cube2, atoms1, data=abs_diff_grid) 


