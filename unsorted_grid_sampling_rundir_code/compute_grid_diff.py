import numpy as np
import ase, ase.io
from ase.io.cube import read_cube
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
     atoms1, grid1 = read_cube(args.grid1, read_data=True)


if args.grid2[-3:] == "cif":
     atoms2, grid2 = read_cif(args.grid2, read_data=True)
elif args.grid2[-3:] == "xsf":
     atoms2, grid2 = read_xsf(args.grid2, read_data=True)
else:
     atoms2, grid2 = read_cube(args.grid2, read_data=True)



