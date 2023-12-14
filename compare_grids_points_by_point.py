#!/usr/bin/env python


import ase
from ase.io.cube import read_cube
import numpy as np
import argparse
import unitcellsampling.sample as sample

parser = argparse.ArgumentParser()

parser.add_argument("grid1", type=str, action='store', help="The 1st grid (cube-file) to be compared. Need to have same atoms object and grid dimensions.")
parser.add_argument("grid2", type=str, action='store', help="The 1st grid (cube-file) to be compared. Need to have same atoms object and grid dimensions.")
parser.add_argument("--vdw", type=float, action='store', help="The vdw cutoff to apply to the respective grids when doing comparison, i.e. exclude points within the cutoff.")
#parser.add_argument("--midvox", action='store_ture', help="Whether to use midvox ")
args = parser.parse_args()


with open(args.grid1, 'r') as f1:
    cont1 = read_cube(f1)
    atoms1 = cont1["atoms"]
    grid1 = cont1["data"]

with open(args.grid2, 'r') as f2:
    cont2 = read_cube(f2)
    atoms2 = cont2["atoms"]
    grid2 = cont2["data"]


assert (grid2.shape == grid1.shape), "Grids have different shapes!"
assert atoms1 == atoms2, "Atoms are not the same!!!"

sampler1 = sample.UnitCellSampler(atoms1)
sampler2 = sample.UnitCellSampler(atoms2)

gridpoints1, incl1 = sampler1.generate_grid_vectors(grid1.shape, vdw_scale=args.vdw, midvox=False)
gridpoints2, incl2 = sampler2.generate_grid_vectors(grid2.shape, vdw_scale=args.vdw, midvox=False)

assert (incl1 == incl2).all(), "incl1 and incl2 are different!!!"

grid1_greater_than_grid2 = np.array(grid1[incl1] > grid2[incl2], dtype=bool)
if grid1_greater_than_grid2.all():
    print("Grid 1 is greater than grid 2 at all points.")
    print("I.e. grid 1 >= grid 2 everywhere")

grid2_greater_than_grid1 = np.array(grid1[incl1] < grid2[incl2], dtype=bool)

grids_equal = np.array(grid1[incl1] == grid2[incl2], dtype=bool)

if grid2_greater_than_grid1.all():
    print("Grid 2 is greater than grid 1 at all points.")
    print("I.e. grid 1 =< grid 2 everywhere")

if grid1_greater_than_grid2.all() and grid2_greater_than_grid1.all():
    print("Grids are equal everywhere!")

if (not grid1_greater_than_grid2.all()) and (not grid2_greater_than_grid1.all()):
    print("No. of points where grid1 < grid2: ", np.sum(grid2_greater_than_grid1) - np.sum(grids_equal))
    print("No. of points where grid1 > grid2: ", np.sum(grid1_greater_than_grid2) - np.sum(grids_equal))
    print("No. of points where grid1 = grid2: ", np.sum(grids_equal))

    print("Total no. of points: ", np.sum(incl1))
