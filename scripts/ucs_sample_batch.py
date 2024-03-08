#!/usr/bin/env python

"""CLI script to UCS sampling, intended for batch sampling of grids.
Called by the ucs_batch_run.py to run individual batches constructed by that script.
Can otherwise be used to sample a set of coordinates in a structure read from file,
as is, without any preprocessing. """

from pathlib import Path
from unitcellsampling import sample

import ase.io
import ase


# energy calcs
from ase.calculators.lj import LennardJones as LJ
from ase.calculators.cp2k import CP2K
from ase.calculators.lammpsrun import LAMMPS
from ase.calculators.lammpslib import LAMMPSlib


# Should change this, to autogenerate the dict as well
#from unitcellsampling.autocreated_methods import *
from unitcellsampling import autocreated_methods # This is better 

# Methods for 1075 of the structures in the rest of the dataset:
from unitcellsampling import autocreated_methods_structures_mcloud_rest_WO_duplabels_nolog_nosubdircalc 

# predefined by Ben
from unitcellsampling.energy_calculator import cp2k_dft

# New cp2k dft methods and wrappers
from unitcellsampling import cp2k_calculators
from unitcellsampling.cp2k_calculators import pass_cp2k_charge
from unitcellsampling.cp2k_calculators import cp2k_calculator_from_input
from unitcellsampling.cp2k_calculators import cp2k2ucs

# from unitcellsampling.cp2k_calculators import pass_cp2k_input

# File handling
from ase.io.cube import write_cube
from ase.io.xsf import write_xsf # Added to test if its better /H 

from unitcellsampling.decorators import subdir_calc
import os
import argparse
import numpy as np
import datetime



# Created: 2024-03-07
# Last update: 2024-03-07 

### Some global logic control variables ###
save_energies_as_txt = True # Adds saving the enrgies as .txt files with np.savetxt, useful for debug or testing perhaps
###

# Define and import methods like in the main cli script

#method_list = ['pbe', 'lammps_lj', 'lammps_lj_coul', 'ff_boulfelfel', 'ff_boulfelfel_buck', 'ff_garcia_sanches', 'ase_lj', '54189', '73679']
method_list = []

# Adding CP2K dft methods 
# This is the simple predefined one, should move it...
dft_methods = {'cp2k_dft':cp2k_dft}

# These are from the cp2k_calculators module
cp2k_dft_methods = {}
cp2k_calc_dict = cp2k_calculators.__dict__

for key in cp2k_calc_dict.keys():
    if key[0:4] == "cp2k":
        cp2k_dft_methods[key] = cp2k_calc_dict[key]


automethods = {}

# Adding UFF + REPEAT force field methods from the 84 intitial structures
autocreated_methods_dict = autocreated_methods.__dict__

for key in autocreated_methods_dict.keys():
    if key[0] == "m" and key[-4:]=="auto":
        automethods["A"+key.strip("m_auto")] = autocreated_methods_dict[key]


# Adding UFF + REPEAT force field methods for 1075 structures in the larger dataset
mcloud_rest_full_dict = autocreated_methods_structures_mcloud_rest_WO_duplabels_nolog_nosubdircalc.__dict__

#mcloud_rest_method_dict = {}
for key in mcloud_rest_full_dict.keys():
    if key[0] == "m" and key[-4:]=="auto":
        automethods["A"+key.strip("m_auto")] = mcloud_rest_full_dict[key]

# Extending the methods list:
method_list.extend(automethods.keys())

method_list.extend(dft_methods.keys())

method_list.extend(cp2k_dft_methods.keys())

# TODO: Remove preprocessing arguments: --nosym, --vdw, --mic-cutoff, --ra, --sg, --guc, --midvox, --conv
# Additionally, remove the following arguments: --wfn, --space, --grid, --name
### Definition and parsing of arguments
parser = argparse.ArgumentParser(description='Energy sampling of a (periodic) solid system with an added atom/ion on the grid.') 
parser.add_argument('file', metavar='filename', type=str, action='store', help="Name of the structure file containing the unit cell structure to be sampled. Formats: cif, cube, xsf, ...")
#parser.add_argument('-m', '--method', type=str, action='store', default='ff_boulfelfel', choices=['pbe', 'ff', 'ff_boulfelfel', 'ff_boulfelfel_buck'], help="Method to calculate the energy during grid sampling.") # should method be optional or not? No, probably not optional? Or maybe optinal with default?
parser.add_argument('method', type=str, action='store', metavar='METHOD', choices=method_list, help="Method to calculate the energy during grid sampling. Example: AXXXXX for UFF force field with REPEAT charges for structure with number XXXXX (implemented for some but not all), cp2k_dft for simple UKS dft with CP2K. For specific cp2k-dft methods, see:" + str(cp2k_dft_methods.keys()))
parser.add_argument('--coords-from-file', type=str, action='store', help="File containing the cartesian coordinates of the grid points that should be sampled. Mandatory argument. Accepted formats are .txt or .npy, both as written by numpy's .savetxt() and save() methods.")
#parser.add_argument('-n','--name', metavar='jobname', type=str, action='store', default=None, help="Desired name for the calculation (applied to generated output files, directories, etc.).")
#parser.add_argument('-w', '--wfn', type=str, action='store', default=None, help="Specify the initial wfn-file for a DFT calculation.")
parser.add_argument('-a', '--atom', type=str, action='store', default='Li', help="Specify the atom/ion used for sampling.")
#parser.add_argument('-g', '--grid', type=int, action='store', default=[10], nargs='+', help="Specify the number of grid points in each dimension (or cubic grid) (mutually exclusive with \"--space\").")
#parser.add_argument('-s', '--space', type=float, action='store', default=None, nargs='+', help="Specify the spacing between the grid points in each dimension (mutually exclusive with \"--grid\").")
parser.add_argument('--cp2k_q', '--cp2k-total-charge', type=str, action='store', default=None, help="Specify the total charge of the structure that is to be sampled. Observe that this is for the full, final structure that is actually sampled, not the input structure. This charge is for instance passed to CP2K DFT calculators to set the number of electrons. Options: int - to set the charge explicitly, \"auto\" for automatic determination based of cp2k_aq, None for no determination or setting of charge (e.g. if present in cp2k input template already).")
#parser.add_argument('--tc', '--total-charge', type=float, action='store', default=None, help="Specify the total charge of the structure that is to be sampled. Observe that this is for the final structure that is actually sampled, not the input structure. This charge is for instance passed to CP2K dft calculators to set the number of electrons.")
#NOT NEEDED: parser.add_argument('--cp2k_aq', '--cp2k-atom-charge', type=int, action='store', default=1, help="Specify the charge of the sampling atom/ion, for the purpose of determining the CP2K total charge. This charge is only used in the automatic determination of the charge that are passed to CP2K DFT calculators. Thus, if cp2k total charge is given as an integer, or if cp2k is not used for sampling, this input is redundant.")

parser.add_argument('--cp2k_template', type=str, action='store', help="Specify the cp2k input file template. Only important for methods using a cp2k input template.")

parser.add_argument('--cp2k_wfn_mode', type=str, action='store', default=None, choices=[None, 'off', 'on', 'seq'], help="Specify the cp2k wfn guess mode, i.e. how the SCF_GUESS is set. By default, it is off. \"on\" or \"seq\" uses a simple sequential strategy, in which each subsequent calculation uses the resulting wfn file from the previous calculation. Only relevant for methods using DFT cp2k methods.")

parser.add_argument('--cp2k_shell_reset_freq', type=int, action='store', default=500, help="Specify the frequency with which the CP2K-shell is re-instantiated, i.e. reset. If this frequency is N, the cp2k shell is killed and restarted every N:th calculation. This is too avoid the logical unit error arising from too many io-units assigned in the cp2k fortran source code. Too disable this, either set it to a higher number than the total number of calculations, or set it to a value <=0. In the latter case the program will automatically set it to a high value so that no reset is performed. Default: 500")

# CP2K command argument - REDUNDANT FOR NOW
#parser.add_argument('--cp2k_cmd', '--cp2k-command', type=str, action='store', default=default_cp2k_cmd, help="Specify the CP2K-command that is used in the ASE-CP2K calculator interface to start and run the CP2K program. Default: ")


args = parser.parse_args()

### End of parser definition

# Print raw input
def print_raw_input(args):
    input_dict = args.__dict__
    print("--------------------------------------------------")
    print("-                   Raw input                    -")
    print("--------------------------------------------------")
    for key,item in input_dict.items():
        print(key,":",item)

    print("--------------------------------------------------")
    print("-                End of raw input                -")
    print("--------------------------------------------------")


def read_coordinates(coord_file_path:str):
    if coord_file_path[-4:] == ".txt":
        coords = np.loadtxt(coord_file_path, dtype=np.float64)
    elif coord_file_path[-4:] == ".npy":
        coords = np.load(coord_file_path)
    else:
        raise ValueError("ERROR: Coordinate file \""+str(coord_file_path)+"\" has unsupported extension. Must be \".txt\" or \".npy\".")
    return coords


### Script start ###
# Print header and date
print("UCS batch sampler")
print("-----------------")
print("Batch sampling run initiated")
print("Date and time: ", datetime.datetime.now().ctime())
print()

# The print raw parsed input
print_raw_input(args)
#

# Settings to control script functions
xsf_output = False
#

## Set inputfile, jobname, method, wfn_file (if method is pbe) and atom variables from input arguments
structure_file = Path(args.file)
method = args.method
coord_file = str(args.coords_from_file)

# don't know if I should keep this?
#if args.name:
#    jobname = args.name

# if method == 'pbe' and not args.wfn:
#     raise Exception()
# if method == 'pbe' and args.wfn:
#     wfn_file = args.wfn
#     restart_file_name = wfn_file # Quite unneccessary / H 

atom = args.atom

## Check existence of inputfile, extract the name of the inputfile and directories
## Construct filenames - the indata file (.cif) is given by structure_file (the relative or absolute path to it)
## Should take the indata-name and extract the non-path name without the .cif ending, which can be done with pathlib.Path 
if not structure_file.exists():
    raise FileNotFoundError('The given inputfile \'' + str(structure_file) + '\' was not found.')

input_basename = structure_file.name
while len(Path(input_basename).suffixes) > 0:
    input_basename = Path(input_basename).stem # input_basename is here meant to be just the name of the input file without parent directories and extensions: /foo/bar.cif -> input_basename = bar
print("Input basename: ", input_basename)


if not Path(coord_file).exists():
    raise FileNotFoundError('The given coordinate file \'' + str(structure_file) + '\' was not found.')


# Set the input filenames and directory, and read input
indir = structure_file.parent
infile =  structure_file.name 

# READ FILES - structure from cif and coordinates from file:
structure = ase.io.read(Path(indir, infile))

coords = read_coordinates(coord_file)

# Check shape:
if not (coords.shape == (len(coords),3)):
    print("Coordinates was of shape ", str(coords.shape),", and not (-1,3): reshaping to (-1,3)...")
    coords = coords.reshape(-1,3)

### End of parsing of arguments

# Current implementation: Check CP2K.command first, second ASE_CP2K_COMMAND
# note that CP2K.command seems to be initiated as None.
if CP2K.command is not None:
    raise Exception('CP2K.command class variable was not None. Please use \
            environment variable $ASE_CP2K_COMMAND to set CP2K command!')

if method in cp2k_dft_methods.keys():
    try:
        ase_cp2k_command = os.environ['ASE_CP2K_COMMAND']
    except KeyError:
        raise Exception('Environment variable ASE_CP2K_COMMAND was not set!\
                Please set this to the desired CP2K command!')
    assert isinstance(ase_cp2k_command, str)
    print("CP2K command used is the ASE_CP2K_COMMAND environment "\
            "variable.")
    print("$ASE_CP2K_COMMAND = ", ase_cp2k_command)


# Here, fetch caluclation directory from environment variable
print()
print("Fetching calculation directory from environment (UCS_CALCULATION_DIR)")
try:
    calc_dir = os.environ["UCS_CALCULATION_DIR"]
except KeyError:
    raise KeyError("ERROR: Could not fetch calculation directory from UCS_CALCULATION_DIR environment variable. The vaiable is likely unset! Please set it to specify the calculation directory.")

print("Calculation directory retrieved from UCS_CALCULATION_DIR enviroment variable")
print("Calculation directory set to: ", calc_dir)
print("Full path: ", os.path.realpath(calc_dir))
print()

#######################################################################

# Make the sampler object
sampler = sample.UnitCellSampler(structure) # DONE
#sampler.n_frac = (nx, ny, nz) # Do not need these for batch sampling I believe
#sampler.n_supercell = num_cells 

print("Setting up energy calculation method...")
# Here, do method selection as in main script (but call sampler after logic section)
# Generated Force-field methods (mostly UFF+REPEAT)
if method in automethods.keys():
    calc_method = automethods[method]

elif method in dft_methods.keys():
    calc_method = dft_methods[method]

# CP2K methods!!! Needs testing!!! 
elif method in cp2k_dft_methods.keys() and method != "cp2k_calculator_from_input":
    cp2k_calc = cp2k_dft_methods[method]

    # TODO: This has been updated in cp2k_calculators.py. Testing needs to be done
    if "requires_charge" in cp2k_calc.__dict__.keys() and cp2k_calc.requires_charge:
        cp2k_calc = pass_cp2k_charge(cp2k_calc, args.cp2k_q)

    calc_method = cp2k_calc

# CP2K from input
elif method in cp2k_dft_methods.keys() and method == "cp2k_calculator_from_input":
    # First: Check the cp2k_template is set, and that it is a file.
    if args.cp2k_template is None:
        raise Exception("Method " + str(method) + " requires a cp2k template \
              input file to be specified through the --cp2k_template argument.")

    cp2k_template_file  = Path(args.cp2k_template) 

    if not cp2k_template_file.exists():
        raise Exception("CP2K template input file " + str(cp2k_template_file) 
                       + " was not found!")

    # Then parse the file into a str passable to the calculator:
    with open(cp2k_template_file, 'r') as f_inp:
        parsed_cp2k_input = "".join(f_inp.readlines()) 

    # Now create the calculator
    cp2k_kwargs = dict(print_level="MEDIUM",
                       basis_set="DZVP-MOLOPT-SR-GTH",
                       pseudo_potential="auto")

    if args.cp2k_q is not None:
        print("Passing total charge = ", int(args.cp2k_q), " to CP2K calculator")
        cp2k_kwargs["charge"] = int(args.cp2k_q)
    else:
        print("--cp2k_q is None - charge is not passed, but that present"
              + " in the inputfile template is used if present.")


    cp2k_calc_from_inp = cp2k_calculator_from_input(parsed_cp2k_input, 
                                                    **cp2k_kwargs)
    # Check charge (sanity check):
    print("Sanity check, cp2k calculator charge: ",
          cp2k_calc_from_inp.__dict__["parameters"]["charge"])

    # Manage wfn file usage:
    if args.cp2k_wfn_mode in ['on', 'seq']:
        wfn_restart = True
    else:
        wfn_restart = False


    # Manage and control the cp2k shell reset frequency setting
    assert isinstance(args.cp2k_shell_reset_freq, int), "ERROR: Shell reset frequency must be integer!!!"
    if args.cp2k_shell_reset_freq <= 0:
        print("CP2K-shell frequency is <= 0. Setting it to 1 + total number of grid points, to disable it.")
        args.cp2k_shell_reset_freq = int(coords.shape[0]) + 1 # Set it to more than the total number of grid points
        print("CP2K-shell reset frequency was set to: ", args.cp2k_shell_reset_freq)
        print("This should be 1 + total number of grid points!")

    print("Shell reset frequency:", args.cp2k_shell_reset_freq)

    # Turn the constructed calculator into a UCS compatible one
    # i.e. singature atoms -> float (energy)
    cp2k_calc = cp2k2ucs(cp2k_calc_from_inp, 
                    base_calc_dir=calc_dir,
                    base_label="cp2k",
                    update_mode="label",
                    restart=wfn_restart,
                    shell_reset_freq=args.cp2k_shell_reset_freq)

    # Print the wfn mode setting:
    if wfn_restart:
        print("wfn_restart = ", wfn_restart,
              ": Each calculation uses wfn from previous calculation",
              " as scf_guess - this is an experimental setting!!!")
    else: 
        print("wfn_restart = ", wfn_restart,
              ": SCF_GUESS is not managed by the sampler",
              " - wfn files are not reused.")

    calc_method = cp2k_calc
    #energies = sampler.calculate_energies(grid_points=supercell_cart_grid,
    #        method=cp2k_calc, atom=atom, exploit_symmetry=use_sym)


else:
    print("No default method defined yet.")
    raise Exception('No method.')


# Here do the sampling after method selection:
# - NOTE: Since preprocing of coordinates and structure is already done,
# and coordinates are simply obtained from file, we do not incorporate symmetry
# or normalization in this case.
print("Starting grid sampling of batch...")
energies = sampler.calculate_energies(grid_points=coords,
            method=calc_method,
            atom=atom,
            exploit_symmetry=False,
            normalize=False)


# OK here we ned to handle output of the sampled points:
# We should perhaps print a list of them, but could also use the indices to print a grid
##

# Write .npy file of energies:
energies_npy_fname ="energies.npy"
np.save(energies_npy_fname, energies)
print("Energies saved to file: ", energies_npy_fname)

if save_energies_as_txt:
    energies_txt_fname ="energies.txt"
    np.savetxt(energies_txt_fname,energies)
    print("Energies saved to txt-file: ", energies_txt_fname)


# NOTE: For now, unless we know the entire grid size (e.g. it is passed to this script)
    # we cannot put the energies into a cube-file. 
    # One possibility is to use the indices in the indices file to produce a cube file.
    # This is left as a TODO!

#### Below here should be the cube writing functonality, if ever implemented ####
#fill_value = np.nan_to_num(np.inf)
#energies_full = np.full(unitcell_included.shape, fill_value=fill_value, dtype=np.float64)
#print("unitcell_included", unitcell_included, np.sum(unitcell_included))
#energies_full[unitcell_included] = energies 
##

# cube_filename = ".".join((calc_name, 'cube'))
# xsf_filename = ".".join((calc_name, 'xsf'))
# print()

# Disabled here: Printing of energy grid values in the output file.
#print(*np.array2string(energies.reshape((len(energies),1)), formatter={"float_kind":lambda x: "%.7f" % x }))
#for row in energies_full:
#    print(row)


# with open(cube_filename, 'w') as fp:
#     write_cube(fp, unitcell, data=energies_full.reshape(
#         unitcell_ucs.included_grid_vectors.shape))

# print("Grid written to cube-file: ", cube_filename)

# if xsf_output:
#     with open(xsf_filename, 'w') as fp:
#         write_xsf(fp, unitcell, data=energies_full.reshape(
#             unitcell_ucs.included_grid_vectors.shape))

#     print("Grid written to xsf-file: ", xsf_filename)
