#!/usr/bin/env python

from unitcellsampling import symmetry
import ase.io
import ase
from ase.spacegroup import get_spacegroup

import pymatgen
from pymatgen.core import Structure, Lattice
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

import numpy as np
import argparse
import os

def print_line(leng=79):
    print(leng*"-")

### Properties to fetch and print: ###
#label_to_prop = {"Formula (as is)":lambda x : x.get_chemical_formula(),
                 #"Number of electrons":lambda x : sum(x.get_atomic_numbers()),
                 #"Electron parity":lambda x : "EVEN" * int(sum(x.get_atomic_numbers()) % 2 == 0) + "ODD" * int(sum(x.get_atomic_numbers()) % 2 == 1),
                 #"Cell parameters":lambda x : dict(zip(("a","b","c","alpha", "beta", "gamma"), x.cell.cellpar())),
                 #"Cell volume":lambda x : x.get_volume(),
                 #"Spacegroup symbol":lambda x : ase.spacegroup.get_spacegroup(x).symbol,
                 #"Spacegroup number":lambda x : ase.spacegroup.get_spacegroup(x).no, 
                 #"Spacegroup order":lambda x : ase.spacegroup.get_spacegroup(x).nsymop
                 #}


parser = argparse.ArgumentParser()
parser.add_argument("-f", "--files", type=str, action='store', nargs='+', default=None, help="Files containing the structures to be analysed.")
parser.add_argument("-d", "--directory", type=str, action='store', default=None, help="Directory containing the structures to be searched among.")
#parser.add_argument("-o", "--output-directory", type=str, action='store', default=None, help="Directory to put cif files in.")

args = parser.parse_args()


assert args.files or args.directory and not (args.files and args.directory), "Error: One and only one of --files and --directory must be given!"


if args.files:
    files = args.files
elif args.directory:
    dr = args.directory
    assert os.path.exists(dr) and os.path.isdir(dr), "Error: {} was not found or is not a directory!".format(dr)
    files = [os.path.join(dr,i) for i in os.listdir(dr)]
else:
   raise Exception("Either --files or --directory arguments must be set!!! This should not have happened!")

with open("num_Li_table.txt", "a") as f:
    f.write("   ".join(("Structure","#Li Prim. Cell.","#Li Conv. Cell","Primitive == Conventional?"))+"\n")
    f.write("-"*len("   ".join(("Structure","#Li Prim. Cell.","#Li Conv. Cell","Primitive == Conventional?")))+"\n")

files.sort()


for f in files:
#    """do the thing"""
   structure = str(f)
   init_unitcell = ase.io.read(f)
   
   # Compute:
   is_conv = symmetry.is_conventional_cell(init_unitcell)
   sg_analyzer = SpacegroupAnalyzer(AseAtomsAdaptor.get_structure(init_unitcell))
   conv_cell = sg_analyzer.get_conventional_standard_structure()
   conv_unitcell = AseAtomsAdaptor.get_atoms(conv_cell)
   #
   init_num_Li = np.sum(init_unitcell.symbols == "Li")
   conv_num_Li = np.sum(conv_unitcell.symbols == "Li")
   #
   props = {}
   props["Is conventional cell"] = is_conv
   props["Init chemical formula"] = init_unitcell.get_chemical_formula()
   props["Init num Li"] = init_num_Li
   props["Init Cell parameters"] = dict(zip(("a","b","c","alpha", "beta", "gamma"), init_unitcell.cell.cellpar()))

   props["Conv. chemical formula"] = conv_unitcell.get_chemical_formula()
   props["Conv num Li"] = conv_num_Li
   props["Conv. Cell parameters"] = dict(zip(("a","b","c","alpha", "beta", "gamma"), conv_unitcell.cell.cellpar()))


   
   print_line()
   print("Structure: ", structure)
   print_line()


   for k in props.keys():
       print(k, " : ", props[k])
   print_line()

   with open("num_Li_table.txt", "a") as f:
       f.write("    ".join((str(structure), str(init_num_Li), str(conv_num_Li), str(is_conv)))+"\n")




