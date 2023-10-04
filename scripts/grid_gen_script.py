#!/usr/bin/env python

# This is based on the unitcellsampling/examples/dft_grid_gen/grid_gen .. script named grid_gen.py by Benjamin Bolbrinker. Modified by Hannes Gustafsson.

from pathlib import Path
from unitcellsampling import sample
from unitcellsampling import symmetry
import ase.io
import ase
from ase.spacegroup import get_spacegroup

import pymatgen
from pymatgen.core import Structure, Lattice
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

# energy calcs
from ase.calculators.lj import LennardJones as LJ
from ase.calculators.cp2k import CP2K
from ase.calculators.lammpsrun import LAMMPS
from ase.calculators.lammpslib import LAMMPSlib

from unitcellsampling.lammps_calc_from_inp import parser_lammps_mel_inp
from unitcellsampling.lammps_calc_from_inp import lammps_method_from_data

from unitcellsampling.special_methods import struct_73679_ff

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
# from unitcellsampling.cp2k_calculators import pass_cp2k_input

# File handling
from ase.io.cube import write_cube
from ase.io.xsf import write_xsf # Added to test if its better /H 

from unitcellsampling.decorators import subdir_calc
import os
import argparse
import numpy as np
import math as ma

from unitcellsampling.preparatory_fcns import unitcell_to_supercell_frac_coords, \
                                              remove_nonframework_cations_fancy, \
                                              compute_req_supercells

# Defaults

# Last update: 2023-10-04 

### TODO: Create input argument handling, to make usage much smoother. Some particular changes:
     # -  TODO: make the  program either recognize the filename from the input file automatically, and name the "project" accordingly
     #          or make the desired project name an input so that the variables all are set in conjuction /H
     # -  TODO: Make the subdirectory for the gridsampling also be project name-derived, and check so that it cannot overide anything by default /H

# TODO: Actually: Move the setting of a project name and setting of environment variables to an exterior runscript - MAYBE, but arguably could be set here too.

# TODO: Test this
# TODO: Move to cp2k-module
# Function determining the charge automatically based on charge of sampling
# ion and how many was removed
def determine_total_cp2k_charge(neutral_unitcell, unitcell_to_sample, 
        sample_ion_symb, sample_ion_charge=1):
    """Determines the total charge of the unitcell that is sampled based on 
    the number and charge of the removed sampling ions.
    Assumes the difference between the neutral unitcell or supercell and the
    cell used in the sampling is only in the number of sample ions."""
    # Type assertions
    assert isinstance(sample_ion_symb, str), "sample_ion_symb must be chemical symbol of smapling ion, of type str."
    assert isinstance(neutral_unitcell, (ase.Atoms, list)), "neutral_unitcell must be Atoms or list object."
    assert isinstance(unitcell_to_sample, (ase.Atoms, list)), "unitcell_to_sample must be Atoms or list object."
    #

    # Make them into list of str-type chemical symbols
    if isinstance(neutral_unitcell, ase.Atoms):
        neutral_atoms_symb = neutral_unitcell.get_chemical_symbols()
    elif all([isinstance(i, str) for i in neutral_unitcell]):
        neutral_atoms_symb = neutral_unitcell

    if isinstance(unitcell_to_sample, ase.Atoms):
        sample_atoms_symb = unitcell_to_sample.get_chemical_symbols()
    elif all([isinstance(i, str) for i in unitcell_to_sample]):
        sample_atoms_symb = unitcell_to_sample
    #
    # Compute the number of sampling ions in neutral and sample cell
    # Add one for the sample ion that will be placed in the sampler
    n_init_sample_ions = len([ion for ion in neutral_atoms_symb if ion == sample_ion_symb])
    n_sampling_ions = len([ion for ion in sample_atoms_symb if ion == sample_ion_symb]) + 1

    # Compute charge:
    # Total charge q = -q_ion * ions_removed = -sample_ion_charge * (n_init_sample_ions - n_sampling_ions) 
    # = sample_ion_charge * (n_sampling_ions - n_init_sample_ions)
    tot_q = (n_sampling_ions - n_init_sample_ions) * sample_ion_charge

    # Print info
    print("Automatically computing CP2K charge...")
    print("No. of sampling ions in intial structure: ", str(n_init_sample_ions))
    print("No of sampling ions in sampling structure: ", str(n_sampling_ions))
    print("Sampling ion charge: ", sample_ion_charge)
    print("Total CP2K charge determined to be: ", str(tot_q))
    #
    return int(tot_q)


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

### Definition and parsing of arguments
parser = argparse.ArgumentParser(description='Energy sampling of a (periodic) solid system with an added atom/ion on the grid.') 
parser.add_argument('file', metavar='filename', type=str, action='store', help="Name of the structure file containing the unit cell structure to be sampled. Formats: cif, cube, xsf, ...")
#parser.add_argument('-m', '--method', type=str, action='store', default='ff_boulfelfel', choices=['pbe', 'ff', 'ff_boulfelfel', 'ff_boulfelfel_buck'], help="Method to calculate the energy during grid sampling.") # should method be optional or not? No, probably not optional? Or maybe optinal with default?
parser.add_argument('method', type=str, action='store', metavar='METHOD', choices=method_list, help="Method to calculate the energy during grid sampling. Example: AXXXXX for UFF force field with REPEAT charges for structure with number XXXXX (implemented for some but not all), cp2k_dft for simple UKS dft with CP2K. For specific cp2k-dft methods, see:" + str(cp2k_dft_methods.keys()))
parser.add_argument('-n','--name', metavar='jobname', type=str, action='store', default=None, help="Desired name for the calculation (applied to generated output files, directories, etc.).")
parser.add_argument('-w', '--wfn', type=str, action='store', default=None, help="Specify the initial wfn-file for a DFT calculation.")
parser.add_argument('-a', '--atom', type=str, action='store', default='Li', help="Specify the atom/ion used for sampling.")
parser.add_argument('-g', '--grid', type=int, action='store', default=[10], nargs='+', help="Specify the number of grid points in each dimension (or cubic grid) (mutually exclusive with \"--space\").")
parser.add_argument('-s', '--space', type=float, action='store', default=None, nargs='+', help="Specify the spacing between the grid points in each dimension (mutually exclusive with \"--grid\").")
parser.add_argument('--vdw', type=float, action='store', default=0.0, help="Specify the fraction of the van der Waals radius that should be excluded from the sampled volume around each atom in the host framework structure.")
parser.add_argument('--cp2k_q', '--cp2k-total-charge', type=int, action='store', default=None, help="Specify the total charge of the structure that is to be sampled. Observe that this is for the full, final structure that is actually sampled, not the input structure. This charge is for instance passed to CP2K DFT calculators to set the number of electrons. If None, it will be automatically determined.")
#parser.add_argument('--tc', '--total-charge', type=float, action='store', default=None, help="Specify the total charge of the structure that is to be sampled. Observe that this is for the final structure that is actually sampled, not the input structure. This charge is for instance passed to CP2K dft calculators to set the number of electrons.")
parser.add_argument('--cp2k_aq', '--cp2k-atom-charge', type=int, action='store', default=1, help="Specify the charge of the sampling atom/ion, for the purpose of determining the CP2K total charge. This charge is only used in the automatic determination of the charge that are passed to CP2K DFT calculators. Thus, if cp2k total charge is given as an integer, or if cp2k is not used for sampling, this input is redundant.")

# CP2K command argument - REDUNDANT FOR NOW
#parser.add_argument('--cp2k_cmd', '--cp2k-command', type=str, action='store', default=default_cp2k_cmd, help="Specify the CP2K-command that is used in the ASE-CP2K calculator interface to start and run the CP2K program. Default: ")
# Minimum-image-convention cutoff argument
parser.add_argument('--mic-cutoff', type=float, action='store', default=0.0, help="Specify the cut-off used to construct a supercell obeying the minimum image convention.")

parser.add_argument('--nosym', action='store_false', help="Turns off usage of spacegroup symmetry. Default is to apply symmetry to save the number of required calculations.")
parser.add_argument('--ra', action='store_true', help="Specify whether to remove all atoms of the type that is used for sampling from the structure, before doing the sampling.")
parser.add_argument('--sg', type=int, default=None, action='store', help="Manually specify the spacegroup to use for symmetry. Default is None, in which case spacegroup will be automatically determined from the structure.")
parser.add_argument('--guc', '--gemmiunitcell', action='store_true', help="If given, will pass unit cell information to gemmi. This is currently only for testing, to see if symmetry will be better handled in certain cases with primitive unitcells. TEST.")
parser.add_argument('--conv', '--conventional-cell', action='store_true', help="If given, will check whether unit cell is conventional cell, and if not, determine the conventional cell based on the input structure and then sample the found conventional cell instead.")
parser.add_argument('--midvox', action='store_true', help="Specifies that sampling should be done in the center of voxels rather than in the corners, as is default. Corresponds to shifting the coordinates of the grid points with 0.5 * (1/na, 1/nb, 1/nc), where ni is the number of grid points in direction of unit cell vector i.")

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

# Before anything, print raw parsed input
print_raw_input(args)
#

# Settings to control script functions
xsf_output = False
#

## Set inputfile, jobname, method, wfn_file (if method is pbe) and atom variables from input arguments
input_file = Path(args.file)
method = args.method

use_sym = args.nosym


if args.name:
    jobname = args.name

if method == 'pbe' and not args.wfn:
    raise Exception()
if method == 'pbe' and args.wfn:
    wfn_file = args.wfn
    restart_file_name = wfn_file # Quite unneccessary / H 

atom = args.atom

## Check existence of inputfile, extract the name of the inputfile and directories
## Construct filenames - the indata file (.cif) is given by input_file (the relative or absolute path to it)
## Should take the indata-name and extract the non-path name without the .cif ending, which can be done with pathlib.Path 
if not input_file.exists():
    raise FileNotFoundError('The given inputfile \'' + str(input_file) + '\' was not found.')

input_basename = input_file.name
while len(Path(input_basename).suffixes) > 0:
    input_basename = Path(input_basename).stem # input_basename is here meant to be just the name of the input file without parent directories and extensions: /foo/bar.cif -> input_basename = bar
print("Input basename: ", input_basename)

# Set the input filenames and directory, and read input
indir = input_file.parent
infile =  input_file.name 

# TODO: Change this reading to the full given path for infile??
input_unitcell = ase.io.read(Path(indir, infile))
init_unitcell = input_unitcell.copy()

# Change the unitcell to the conventional cell if symmetry 
# is used or if conv argument is specified
# This is a bit convoluted as we may be able to use ase
# directly I have realized now...
if use_sym or args.conv:
    if not symmetry.is_conventional_cell(init_unitcell):
        sg_analyzer = SpacegroupAnalyzer(AseAtomsAdaptor.get_structure(init_unitcell)) 
        conv_cell = sg_analyzer.get_conventional_standard_structure()
        init_unitcell = AseAtomsAdaptor.get_atoms(conv_cell)

#atoms = ase.io.read(input_file) # Like so? 

if args.ra:
    print("Removing sampling atoms", str(atom),"from structure.")
    atoms_temp = remove_nonframework_cations_fancy(init_unitcell, ase.Atom(atom))
    unitcell = atoms_temp
else:
    unitcell = init_unitcell

if use_sym:
    if args.sg:
        spacegroup = args.sg
    else:
        spacegroup = get_spacegroup(unitcell) 

## Set the number of points in the grid in each dimension (or equivalently, the mesh size)
# TODO: If symmetr is used, check/determine nearest shape that is compatible with spacegroup
if not args.space:
    if len(args.grid) == 1:
        nx,ny,nz = args.grid * 3
    elif len(args.grid) == 2:
        raise Exception("--grid must be either 1 or 3 int inputs.")
    else:
        nx,ny,nz = args.grid[0:3]
    
    a,b,c = np.linalg.norm(unitcell.get_cell()[0,:]),\
            np.linalg.norm(unitcell.get_cell()[1,:]),\
            np.linalg.norm(unitcell.get_cell()[2,:])
    true_spacing = (a/nx, b/ny, c/nz)

    # There should be a symmetry adaption here probably ... 

    print('True spacing: ', true_spacing)
    print('(nx, ny, nz) =', (nx, ny, nz))

else:
    # Small section to compute nx, ny, nz i.e. the # of points in each direction in case the --space option is used
    # This however is not needed - sampler.generate_grid_vectors has this functionality already, so we can just pass arguments to this
    # Although, that option does parheps something odd - it might ue interpolation to improve the spacing... 
    if args.space:
        if len(args.space) == 1:
            spacing_x, spacing_y, spacing_z = args.space * 3
        elif len(args.space) == 2:
            raise Exception("--space must be either 1 or 3 float inputs.")
        else:
            spacing_x, spacing_y, spacing_z =  args.space[0:3]

        a,b,c = np.linalg.norm(unitcell.get_cell()[0,:]),\
                np.linalg.norm(unitcell.get_cell()[1,:]),\
                np.linalg.norm(unitcell.get_cell()[2,:])
        nx,ny,nz = ma.ceil(a/spacing_x), ma.ceil(b/spacing_y), ma.ceil(c/spacing_z)
        
        if use_sym:
            nx,ny,nz = symmetry.find_spacegroup_compatible_gridshape((nx,ny,nz), spacegroup, search_denser=True)
            print("Determined grid spacing compatible with spacegroup.")

        true_spacing = (a/nx, b/ny, c/nz)
        print('Desired spacing: ', (spacing_x, spacing_y, spacing_z),' True spacing: ', true_spacing)
        print('(nx, ny, nz) =', (nx, ny, nz))

print("vdW cutoff factor: ", args.vdw)

### End of parsing of arguments

# Now handle naming of the calculation files: Chosen format: inputfile_method_gridsize_jobname
# Could addversion if spacing is given, to indicate spacing?
#calc_name = '_'.join((input_basename, method, 'x'.join((str(nx), str(ny), str(nz))), jobname)) # For grid size
if args.name:
    calc_name = '_'.join((input_basename, method, 'x'.join(tuple(map(str, np.around(true_spacing,4)))), jobname))
else:
    calc_name = '_'.join((input_basename, method, 'x'.join(tuple(map(str, np.around(true_spacing,4))))))



### Set the directories and related/relevant environment vars
#calc_dir = './UCS_CALCS_s0.4/' + calc_name
#calc_dir = './ucs_out_mel_data/' + calc_name
work_dir = '.'
os.environ['UCS_WORK_DIR'] = work_dir


# CP2K Things:
# Examples:
#CP2K.command = "env OMP_NUM_THREADS=32 srun cp2k_shell.psmp" Does it work if srun runs within here?
#CP2K.command = "env OMP_NUM_THREADS=4 cp2k_shell"   ## Should be the same / H
# NOTES: CP2K() in ase can be passed a cp2k_cmd option. If it isn't, it checks in order the
# 1. the CP2K.command class variable, 2. the $ASE_CP2K_COMMAND 3. cp2k_shell
# So I could either set the CP2K.command class variable or I can set the ASE_CP2K_COMMAND environment variable.
# For now, externally setting the ASE_CP2K_COMMAND environment vatiable seems most reasonable.
# However, for full control I should check that the CP2K.command variable is not set already.
#CP2K.command = args.cp2k_cmd # Though input: Check default as well

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
    print("CP2K command used is the ASE_CP2K_COMMAND environment\
            variable.")
    print("$ASE_CP2K_COMMAND = ", ase_cp2k_command)


#if not os.path.isdir(Path(calc_dir)):
#    os.mkdir(calc_dir)                ## Might be better later to export these in terminal bash script /H
#os.environ['UCS_CALCULATION_DIR'] = calc_dir  #'./LTA_lammps_ff_grid_test' #./LTA4A_reduced_grid_gen'

# We leave this to the enviroment setting for now
#os.environ['OMP_NUM_THREADS'] = '4'


# CP2K stuff, for DFT 
# Path to previously relaxed single point calculation with an atom placed at
# the first sampled position in this case (0, 0, 0)
#restart_file_name = "./single_point_calc/LTA4A_reduced-RESTART.wfn"  ## This is the old line /H 


#####################################################################
# Create supercell from unitcell depending on cutoff !!!
####################################################################
# This codeblock is taken from custom_lammps_grid_gen.py

# Metainfo: unitcell = unitcell atoms object whose sampling atoms have been removed if desired
#           sampler = UnitcellSampler object for the unitcell

# First construct supercell
#cutoff = 12.5 # Force cutoff in Å
cutoff =  args.mic_cutoff # Force cutoff in Å, which is used to determine supercell compliant with MIC
print("Force cutoff used to determine supercell size: ", cutoff)
num_cells = compute_req_supercells(unitcell, cutoff)


print("num_cells in supercell: ", num_cells)

supercell_from_init_unitcell = ase.build.make_supercell(
            init_unitcell,
            np.diag(num_cells), wrap=True
            )

supercell_from_unitcell_wo_ions = ase.build.make_supercell(
            unitcell,
            np.diag(num_cells), wrap=True
            )

# Here goes charge determination for cp2k:
# TODO: Think about which method to use - supercell seems better now
# Also, test it, see that it works...
if not args.cp2k_q or args.cp2k_q == "auto":
    # Determine charge per unitcell:
    cp2k_charge_per_unitcell = determine_total_cp2k_charge(init_unitcell, 
            unitcell, args.atom, args.cp2k_aq)

    # Then multiply with number of unit cells in supercell,
    # and remove one for each
    cp2k_charge_for_supercell_alt = ma.prod(num_cells) * cp2k_charge_per_unitcell \
                              - (ma.prod(num_cells) - 1) * args.cp2k_aq

    # Copmuting from supercells now
    print()
    print("Computing CP2K total charge based on supercells")
    cp2k_charge_for_supercell = determine_total_cp2k_charge(
            supercell_from_init_unitcell, supercell_from_unitcell_wo_ions, 
            args.atom, args.cp2k_aq)

    assert cp2k_charge_for_supercell_alt == cp2k_charge_for_supercell, "Total CP2K charge not equal with different methods!!! Check charge determination..."

    # Set args namespace
    args.cp2k_q = cp2k_charge_for_supercell
    


print("supercell from uc params: ", supercell_from_unitcell_wo_ions.get_cell_lengths_and_angles())
print("supercell from uc cell: ", supercell_from_unitcell_wo_ions.get_cell())

print("\n")
print("Spacegroup, unitcell: ", ase.spacegroup.get_spacegroup(unitcell, 1.0e-6))
print("Spacegroup, supercell: ", ase.spacegroup.get_spacegroup(supercell_from_unitcell_wo_ions, 1.0e-6))

# "Middle" unit cell (we don't need this though, since pbc)
#uc_indices = [int((i//2)-1) if i % 2 == 0 else int((i-1)//2) for i in num_cells]
uc_indices = (0, 0, 0)


# Now generate grid for unitcell:
unitcell_ucs = sample.UnitCellSampler(unitcell) # DONE

if (args.midvox and use_sym):
    print("WARNING: Both symmetry and midvoxel sampling enabled. Their compatibility has NOT been properly tested!!!")

unitcell_grid, unitcell_included = unitcell_ucs.generate_grid_vectors((nx, ny, nz), vdw_scale=args.vdw, midvox=args.midvox) # DONE

unitcell_grid = unitcell_grid.reshape(-1,3) # DONE
unitcell_included = unitcell_included.reshape(-1) # DONE

# Convert to fractional coordinates, and convert to supercell grid

print("Shape check, grid: ", np.array(unitcell_grid).shape, np.array(unitcell_grid).T.shape) # DONE
unitcell_frac_grid = np.linalg.solve(np.array(unitcell.get_cell()).T, np.array(unitcell_grid).T).T # DONE

print("unitcell frac grid shape: ", unitcell_frac_grid.shape) #DONE
supercell_frac_grid = unitcell_to_supercell_frac_coords(unitcell_frac_grid[unitcell_included], num_cells, unitcell_ind=(0,0,0)) #DONE

# Convert to cartesian supercell grid
supercell_cart_grid = supercell_frac_grid @ np.array(supercell_from_unitcell_wo_ions.get_cell()) # DONE
print(type(supercell_cart_grid)) # DONE


#######################################################################

## Old version of grid_gen_script.py:
#sampler = sample.UnitCellSampler(lgps)
#sampler.generate_grid_vectors(n_frac=(nx, ny, nz), vdw_scale=args.vdw)
#sampler.spacegroup = args.sg # Set spacegroup for the sampler

# New, supercell sampler
sampler = sample.UnitCellSampler(supercell_from_unitcell_wo_ions) # DONE
sampler.n_frac = (nx, ny, nz)
sampler.n_supercell = num_cells 

if args.sg:
    sampler.spacegroup = args.sg # Set spacegroup for the sampler

if args.guc:
    print("Passing unitcell information (input structure) to UCS and gemmi.")
    sampler.gemmi_unitcell = unitcell

print("Spacegroup input: ", args.sg)
print("Symmetry: ", use_sym)


#if method == 'pbe':
    # For DFT:
#    energies = sampler.calculate_energies(grid_points=supercell_cart_grid,
#        method=dft_pbe, atom=atom, exploit_symmetry=use_sym)

#elif method == '73679':
#    energies = sampler.calculate_energies(grid_points=supercell_cart_grid,
#            method=struct_73679_ff, atom=atom, exploit_symmetry=use_sym)

if method in automethods.keys():
    energies = sampler.calculate_energies(grid_points=supercell_cart_grid,
            method=automethods[method], atom=atom, exploit_symmetry=use_sym)

elif method in dft_methods.keys():
    energies = sampler.calculate_energies(grid_points=supercell_cart_grid,
            method=dft_methods[method], atom=atom, exploit_symmetry=use_sym)

# CP2K methods!!! Needs testing!!! 
elif method in cp2k_dft_methods.keys():
    cp2k_calc = cp2k_dft_methods[method]

    # TODO: This has been updated in cp2k_calculators.py. Testing needs to be done
    if "requires_charge" in cp2k_calc.__dict__.keys() and cp2k_calc.requires_charge:
        cp2k_calc = pass_cp2k_charge(cp2k_calc, args.cp2k_q)

    energies = sampler.calculate_energies(grid_points=supercell_cart_grid,
            method=cp2k_calc, atom=atom, exploit_symmetry=use_sym)

else:
    print("No default method defined yet.")
    raise Exception('No method.')

print("Included grid vectors shape: ", unitcell_ucs.included_grid_vectors.shape)

##
fill_value = np.nan_to_num(np.inf)
energies_full = np.full(unitcell_included.shape, fill_value=fill_value, dtype=np.float64)
print("unitcell_included", unitcell_included, np.sum(unitcell_included))
energies_full[unitcell_included] = energies 
##

cube_filename = ".".join((calc_name, 'cube'))
xsf_filename = ".".join((calc_name, 'xsf'))
print()
#print(*np.array2string(energies.reshape((len(energies),1)), formatter={"float_kind":lambda x: "%.7f" % x }))
for row in energies_full:
    print(row)


with open(cube_filename, 'w') as fp:
    write_cube(fp, unitcell, data=energies_full.reshape(
        unitcell_ucs.included_grid_vectors.shape))

if xsf_output:
    with open(xsf_filename, 'w') as fp:
        write_xsf(fp, unitcell, data=energies_full.reshape(
            unitcell_ucs.included_grid_vectors.shape))

