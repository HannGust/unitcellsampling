#!/usr/bin/env python

"""A script to be able to print properties from structure files, 
e.g. atomic, structural and spacegroup-symmetric properties. This is
to be able to more quickly and easily choose or select among sturctures."""
# OBS: WORK IN PROGRESS! NOT FINISHED!
# MERELY AN IDEA 

import ase
import ase.io
import ase.spacegroup
import os
import argparse

def print_line(leng=79):
    print(leng*"-")

### Properties to fetch and print: ###
label_to_prop = {"Formula (as is)":lambda x : x.get_chemical_formula(),
                 "Number of electrons":lambda x : sum(x.get_atomic_numbers()),
                 "Electron parity":lambda x : "EVEN" * int(sum(x.get_atomic_numbers()) % 2 == 0) + "ODD" * int(sum(x.get_atomic_numbers()) % 2 == 1),
                 "Cell parameters":lambda x : dict(zip(("a","b","c","alpha", "beta", "gamma"), x.cell.cellpar())),
                 #"Cell volume":lambda x : x.get_volume(),
                 #"Spacegroup symbol":lambda x : ase.spacegroup.get_spacegroup(x).symbol,
                 #"Spacegroup number":lambda x : ase.spacegroup.get_spacegroup(x).no, 
                 #"Spacegroup order":lambda x : ase.spacegroup.get_spacegroup(x).nsymop
                 }



parser = argparse.ArgumentParser()
parser.add_argument("-f", "--files", type=str, action='store', nargs='+', default=None, help="Files containing the structures to be analysed.")
parser.add_argument("-d", "--directory", type=str, action='store', default=None, help="Directory containing the structures to be searched among.")

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

for f in files:
#    """do the thing"""
   structure = str(f)
   atoms = ase.io.read(f)
   props = {}
   props = {key: label_to_prop[key](atoms) for key in list(label_to_prop.keys())}
   print_line()
   print("Structure: ", structure)
   print_line()
   for k in props.keys():
       print(k, " : ", props[k])
   print_line()
   
   
 

