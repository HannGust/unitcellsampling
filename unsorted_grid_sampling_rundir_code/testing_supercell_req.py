import ase 
import ase.io
from ase.io import read
from ase.io.cube import read_cube
from pathlib import Path
from preparatory_fcns import compute_req_supercells

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("file", type=str, action='store', help='File containing atoms object of unitcell (must be readable with ase.io.read).')
parser.add_argument("cutoff", type=float, action='store', help='Cutoff to be tested in determining the supercell size.')

args = parser.parse_args()

file = Path(args.file)
cutoff = args.cutoff

atoms = read(file)


assert isinstance(atoms, ase.Atoms), "atoms is not ase.Atoms object. Something went wrong after reading file."

supercell_size = compute_req_supercells(atoms, cutoff=cutoff)

print("Determined supercell size: ", supercell_size)

