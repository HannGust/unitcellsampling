#!/usr/bin/env python3

# File to generate a single point calculation from which to start grid sampling

import sys
import ase
#form ase.spacegroup import Spacegroup
from ase.io import read
from ase.calculators import cp2k
from ase.calculators.cp2k import CP2K
import os
from pathlib import Path

print(sys.argv)

file = sys.argv[1]

atoms = read(file)

print(atoms.get_chemical_formula())
print(atoms.get_global_number_of_atoms())
#print(atoms.get_charges())
print(atoms.get_initial_charges())
print(sum(atoms.get_initial_charges()))
print(atoms.info["spacegroup"].symbol)

# Preliminary approximate computation of charges
charge_dict = {'O':-2, 'Si':4, 'Al':3, 'Na':1, 'Li':1}
comp_charge = sum([charge_dict[a] for a in atoms.get_chemical_symbols()])
print("Approximately estimated charge (no Li+ nor Na+): ", comp_charge)
print("Approximately estimated charge with q=+1 sampling ion: ", comp_charge+1)

