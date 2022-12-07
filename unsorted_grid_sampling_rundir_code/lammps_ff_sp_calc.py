#!/usr/bin/env python3

# File to generate a single point calculation from which to start grid sampling
# Adds an atom/ion of the type we want to sample with (typically Li+, Na+)

import sys
import ase
from ase.spacegroup import get_spacegroup
from ase.io import read
from ase.calculators import cp2k
from ase.calculators.cp2k import CP2K
from ase.calculators.lammpsrun import LAMMPS
from ase.calculators.lammpslib import LAMMPSlib

import os
from pathlib import Path
import argparse

os.environ['OMP_NUM_THREADS'] = '4'
#os.environ['OMP_NU]


parser = argparse.ArgumentParser(description='Run inital single point calculation with a Force-Field (), for grid sampling.\n Optionally adds a sampling atom at (0,0,0).')
parser.add_argument('file', metavar='filename', type=str, action='store', help="Name of the cif-file containing the structure, without sample atom/ion.")
parser.add_argument('--noatom', action='store_false', help="Don't add a sample atom before calculating; only calculate the structure as-is.")
parser.add_argument('--atom', type=str, action='store', help="For specifying the sample atom/ion, to be placed at (0,0,0). (Not implemented yet)", default="Na")
parser.add_argument('--nocalc', action='store_true', help="Does not perform calculation, but does everything before that, i.e. reads input and estimates charge.")
parser.add_argument('--proj', type=str, action='store', help='Specifies an addition to the project name, in addition to using the filename.')
#parser.add_argument('--charge', type=int, action='store', nargs=1, help="For specifying the total charge.", default=0)
#parser.add_argument('--aq', action='store_true', help="Sets the total charge of the system to the estimated one. Overides \"--charge\"")


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

if args.proj:
    proj_name = proj_name + "_" + args.proj



print('')
print('Attempting to run a single point calculation now, of the structure in the inputfile (cif).')

## Here lets attempt a single point calculation using ase

print("Project name: ", proj_name)
proceed_var = input("Do you want to proceed? (y/n) (def: n): ")

if not (proceed_var == "y" or proceed_var == "Y"):
    print("Aborting.")
    sys.exit("Aborted by user.")



# set path to single point calc directory, and change to there

p = Path('.')
calcdir = p / 'lammps_ff_single_point_calc'

if not (os.path.exists(calcdir) and os.path.isdir(calcdir)):
    os.mkdir(calcdir)

os.chdir(calcdir)

# Preliminary approximate computation of charges
charge_dict = {'O':-2, 'Si':4, 'Al':3, 'Na':1, 'Li':1} # Naive charges (oxidation numbers) /H
charge_dict_from_ff = {'O':-0.9364, 'Si':1.8708, 'Al':1.7906, 'Na':0.9094, 'Li':0.9094} # Essentially those atomic charges from ff by Sholl et al 2021


comp_charge = sum([charge_dict[a] for a in atoms.get_chemical_symbols()])
comp_charge_ff = sum([charge_dict_from_ff[a] for a in atoms.get_chemical_symbols()])
print("Approximately estimated (naive) charge: ", comp_charge)
print("Approximately estimated ff-charge: ", comp_charge_ff)

# Setting the initial charges here; hopefully lammps obtains the initial charges from ase ...
atoms.set_initial_charges(charges=[charge_dict_from_ff[a] for a in atoms.get_chemical_symbols()])
print('Initial charges, array: ', atoms.get_initial_charges())


# Set the total charge for FF not sure we need it
#if args.aq:
#    charge = comp_charge
#else:
#    charge = args.charge
#
#print("Total charge set to ", charge)

# construct calculator, and compute energies

log_file = proj_name+".log"

# Set the atom type numbers for lammps
lammps_atom_types = {'Si':1, 'Al':2, 'O':3, 'B':4, 'Na':5, 'Li':6} # Note B is dummy for now - it is originally Si-O-Al while type 3 is Si-O-Si / H

# lammps header to control units etc.
lammps_header = ['units real', 
                 'atom_style full'] 


# Add setting of atomic charges in amendments 
lammps_amendments = ['set type 1 charge 1.8708', 
                 'set type 2 charge 1.7906', 
                 'set type 3 charge -0.9364', 
                 'set type 4 charge -1.1427', 
                 'set type 5 charge 0.9094'] 

# lammps commmands defining the force field in [SOURCE Sholl et al 2021] for the framework and cations
lammps_cmds = ['pair_style hybrid/overlay buck 11.0 coul/long 11.0', 
                  'pair_modify tail yes', 
                  'kspace_style ewald 1.0e-6', 
                  'pair_coeff * * coul/long', 
                  'pair_coeff 3 5 buck 110905.00 0.25094 1821.66']  # Define the lattice-cation interaction /H 
# First line i modified from 'pair_style hybrid/overlay lj/cut 11.0 buck 11.0 coul/long 11.0', 


# This is just for testing with a very simple lammps lennard jones ff
#lammps_cmds = ['pair_style lj/cut 8.0', 'pair_coeff * * 10.0 2.0'] ## Just define something simple to see if it works as intended .. /H  
#

calc = LAMMPSlib(lmpcmds=lammps_cmds, atom_types=lammps_atom_types, lammps_header=lammps_header, amendments=lammps_amendments, log_file=log_file)

atoms.set_calculator(calc)

if args.nocalc:
    print("--nocalc specified ; no calculation performed.")
else:
    energies = atoms.get_potential_energy()


