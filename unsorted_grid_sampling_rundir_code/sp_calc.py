#!/usr/bin/env python3

# File to generate a single point calculation from which to start grid sampling
# Adds an atom/ion of the type we want to sample with (typically Li+, Na+)

import sys
import ase
from ase.spacegroup import get_spacegroup
from ase.io import read
from ase.calculators import cp2k
from ase.calculators.cp2k import CP2K
import os
from pathlib import Path
import argparse

parser = argparse.ArgumentParser(description='Run inital single point calculation with DFT (PBE/PBE-GTH/DZVP-MOLOPT-SR-GTH), for grid sampling.\n Optionally adds a sampling atom at (0,0,0).')
#parser.add_argument('file', metavar='filename', type=str, action='store', help="Name of the cif-file containing the structure, without sample atom/ion.")
#parser.add_argument('name', metavar='calcname', type=str, action='store', help="Desired name for the calculation (applied to generated output files, directories, etc.).")
#parser.add_argument('method', metavar='method', type=str, action='store', choices=['pbe', 'ff', 'ff_sholl'], help="Method to calculate the energy during grid sampling.")
parser.add_argument('--noatom', action='store_false', help="Don't add a sample atom before calculating; only calculate the structure as-is.")
parser.add_argument('--atom', type=str, action='store', help="For specifying the sample atom/ion, to be placed at (0,0,0). (Not implemented yet)", default="Na")
parser.add_argument('--nocalc', action='store_true', help="Does not perform calculation, but does everything before that, i.e. reads input and estimates charge.")
parser.add_argument('--charge', type=int, action='store', nargs=1, help="For specifying the total charge.", default=0)
parser.add_argument('--aq', action='store_true', help="Sets the total charge of the system to the estimated one. Overides \"--charge\"")


args = parser.parse_args()

file = args.file


print("Arguments: ",sys.argv)


atoms = read(file)

sg = get_spacegroup(atoms,symprec=1.0e-05)
#atoms.info["spacegroup"]=sg

print("Chemical formula: ", atoms.get_chemical_formula())
print("Total num. atoms: ", atoms.get_global_number_of_atoms())
print("Init. Spacegroup: ", atoms.info["spacegroup"].symbol)
print("Found Spacegroup: ", sg.symbol)



if args.noatom:
    print('')
    print('Adding sampling Na atom to structure.')
    sample_atom = ase.Atoms('Na', positions=[(0.0, 0.0, 0.0)])
    atoms.extend(sample_atom)
    print(atoms.get_chemical_formula())
    print(atoms.get_global_number_of_atoms())
    print(atoms.info["spacegroup"].symbol)



proj_name = Path(file).stem

print('')
print('Attempting to run a single point calculation now, of the structure in the inputfile (cif).')

## Here lets attempt a single point calculation using ase
CP2K.command = "env OMP_NUM_TREADS=4 cp2k_shell" # The cp2k command

# unclear if we need this - might be useful in more complex cases of course, but for now I'll sart by defining calcs via ase
#inp_temp = """
#&GLOBAL
#   PROJECT {}
#   RUN_TYPE ENERGY_FORCE
#   PRINT_LEVEL MEDIUM
#&END GLOBAL
#&FORCE_EVAL
#   &SUBSYS
#       
#       &TOPOLOGY
#           COORD_FILE_FORMAT CIF
#           COORD_FILE_NAME {}
#       &END TOPOLOGY
#&END FORCE_EVAL
#"""

print("Project name: ", proj_name)
inp = ""
#"""&GLOBAL
#   PROJECT {}
#   RUN_TYPE ENERGY_FORCE
#   PRINT_LEVEL MEDIUM
#&END GLOBAL
#""".format(proj_name)

#

# set path to single point calc directory, and change to there

p = Path('.')
calcdir = p / 'single_point_calc'

if not (os.path.exists(calcdir) and os.path.isdir(calcdir)):
    os.mkdir(calcdir)

os.chdir(calcdir)

# Preliminary approximate computation of charges
charge_dict = {'O':-2, 'Si':4, 'Al':3, 'Na':1, 'Li':1}
comp_charge = sum([charge_dict[a] for a in atoms.get_chemical_symbols()])
print("Approximately estimated charge: ", comp_charge)



# Set the total charge
if args.aq:
    charge = comp_charge
else:
    charge = args.charge

print("Total charge set to ", charge)

# construct calculator, and compute energies
calc = CP2K(inp=inp, label=proj_name, basis_set="DZVP-MOLOPT-SR-GTH", charge=charge, pseudo_potential='GTH-PBE',
                max_scf=600, uks=True, xc='PBE', print_level='LOW')

atoms.set_calculator(calc)

if args.nocalc:
    print("--nocalc specified ; no calculation performed.")
else:
    energies = atoms.get_potential_energy()


