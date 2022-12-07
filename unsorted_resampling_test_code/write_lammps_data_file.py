import ase
import ase.io
from ase.io.cube import read_cube
from ase.io import read
from ase.io.lammpsdata import write_lammps_data, read_lammps_data

import sys
from pathlib import Path

import argparse

from read_atom_types import read_lammpsdata_atom_info

from preparatory_fcns import unitcell_to_supercell_frac_coords, \
                             remove_nonframework_cations_fancy, \
                             compute_req_supercells
#
parser  = argparse.ArgumentParser(description="Can take a cif, cube or lammpsdata file, read an atoms object, add a sampling atom and write out a new lammpsdatafile.")
parser.add_argument('-a','--atom', type=str, action='store', default=None, help='Atom to sample with. Deafult: None')
parser.add_argument('-c','--coord', type=float, action='store', nargs=3, default=(0.0, 0.0, 0.0), help='Coordinates to put sample atom at. Deafult: (0.0, 0.0, 0.0)')
parser.add_argument('--lammpsdata', type=str, action='store', default=None, help='Lammps data file to read from.')
parser.add_argument('-d', '--del_sample_atoms', action='store_true', help='If specified, removes all atoms of the same kind as the sample atom from the input structure before adding the sample atom.')
parser.add_argument('--cif', type=str, action='store', default=None, help='Cif file to read from.')
parser.add_argument('--cube', type=str, action='store', default=None, help='Cube file to read from.')

args = parser.parse_args()
#

if args.cif:
    with open(args.cif, 'r') as f:
        atoms = read(f)
elif args.cube:
    with open(args.cube, 'r') as f:
        atoms = read_cube(f)
elif args.lammpsdata:
    with open(args.lammpsdata, 'r') as f:
        atoms = read_lammps_data(f)
        () = read_lammpsdata_atom_info(f)
else:
    raise Exception("Need an input file to read from.")


if args.a:
    sample_atom = ase.Atom(args.a, list(args.c))
else:
    sample_atom = None


if len(sys.argv) == 2:
    inp = sys.argv[1]
    add_atom = False
elif len(sys.argv) > 2:
    inp = sys.argv[1]
    add_atom = True
    print('Adding sample Na atom to structure at [0,0,0].')
else:
    raise Exception("Need an input.")


with open(inp, 'r') as fp:
    atoms = read(inp)
    assert isinstance(atoms,ase.Atoms), "atoms is not ase.Atoms object!!!"



if add_atom:
    atoms = atoms + sample_atom
    out = Path(inp).stem + "_with_Na.lmpdata"
else:
    out = Path(inp).stem + ".lmpdata"



with open(out, 'w') as fp:
    write_lammps_data(out, atoms)
