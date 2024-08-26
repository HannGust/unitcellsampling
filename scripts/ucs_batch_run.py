#!/usr/bin/env python

"""Runscript for unitcellsampling that can batch a grid calculation into 
smaller batches of points, and setup and run independent parallel calculations
of the batches. Very useful in a cluster environment."""


# What does this script need to do?
# It needs to do the preprocessing of the structure,
# then find the grid points to-be sampled, then batch them,
# then create a hiearchical directory structure for the
# batch runs, and put the relevant files in their place there.

# Lastly, it should submit the calculations as parallel jobs.

# Update how the calculation directory is determined/set? - > environment variable/cli argument?

# E.g. UCS_CALCULATION_DIR , UCS_WORK_DIR  ?

# ??? We have an out dir. Maybe create a work dir?

# "GLOBAL" File hierarchy
# WORK_DIR/UCS_CALC/
#                  /BATCH_1 ...
#                  /BATCH_N
#         /UCS_OUT/batch_1.out ...
#         
#
# CALC_DIR/BATCH_1
#         /BATCH_2
#         ...
#         /BATCH_N/cp2k_1.inp /cp2k_1.out ...         

# OR :
# "LOCAL" File hierarchy
# WORK_DIR/BATCH_1/ .out, .cube, .coords, .indices, .cif , cp2k_template.inp, jobsubmit.sh, py-script.py calc/
# NOTE: We can have WORK_DIR/structure.cif , cp2k_template.inp and even py-script.py or jobscript? these will be the same



# I have chosen to go for a "LOCAL" directory structure

# Also, in the first edition, first goal, aim to use the same grid_gen_script as before, and just pass relevant arguments and structure as is.
# Then, you only need to implement a --points-from-file or something akin to that, and some way of 
# maybe having better control over the working directory and the calculation directory.

# I probably have to incorporate the UCS_CALCULATION_DIR and UCS_WORK_DIR environment variables again.

# pseudo code:

# 1. parse arguments for grid_gen_script and batching configs
#    e.g. batch size, jobsubmit.sh-template.

# 2. Prepare arguments to pass to each individual batch calculation.

# 3. Use symmetry and volume exclusion to produce list of points to sample

# 4. Decompose the list of points to sample into batches with desired size

# 5. Prepare a jobscript based on template.

# 6. Create directory structure:
#       - Create main job directory (based on project name or so)
#       - For each batch, create a subdirectory that acts working directory
#       - Move all necessary files and info to the respective batch directory

# 7. Run the batches separately. Decend into each batch directory one by one
#    and (optionally) run the jobs by calling the batch-jobsubmit.sh script.


# information flow: 
# Structure, method, points to sample, 



# Reusing code from grid_gen_script here:

from pathlib import Path
from unitcellsampling import sample
from unitcellsampling import symmetry
import ase.io
import ase
from ase.spacegroup import get_spacegroup
import gemmi

import pymatgen
from pymatgen.core import Structure, Lattice
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

# Should change this, to autogenerate the dict as well
#from unitcellsampling.autocreated_methods import *
from unitcellsampling import autocreated_methods # This is better 

# Methods for 1075 of the structures in the rest of the dataset:
from unitcellsampling import autocreated_methods_structures_mcloud_rest_WO_duplabels_nolog_nosubdircalc 

from unitcellsampling.energy_calculator import cp2k_dft, Li13Si4_test

# New cp2k dft methods and wrappers
from unitcellsampling import cp2k_calculators
from unitcellsampling.cp2k_calculators import pass_cp2k_charge
from unitcellsampling.cp2k_calculators import cp2k_calculator_from_input
from unitcellsampling.cp2k_calculators import cp2k2ucs

# File handling
from ase.io.cube import write_cube
from ase.io.xsf import write_xsf # Added to test if its better /H 

import os
import sys
import shutil
import subprocess
import argparse
import numpy as np
import math as ma

import datetime

from unitcellsampling.preparatory_fcns import unitcell_to_supercell_frac_coords, \
                                              remove_nonframework_cations_fancy, \
                                              compute_req_supercells



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
    return float(tot_q)


#########################################################################
########       Hardcoded options for the batching process:       ########
#########################################################################

#### Logic switches

batch_coord_txt_format = True # Switch to write coordinate/indices files in the .txt format if True
batch_coord_npy_format = False # Switch to write coordinate/indices files in the .txt format if True

assert batch_coord_txt_format or batch_coord_txt_format, "Error: You must switch on at least one format option (either txt or npy) for the coordinate/indices files."

#### Generic filenames for the batch run
batch_structure = "structure.cif"
batch_jobscript_file = "jobsubmit.sh"
batch_indices_file = "indices"
batch_coord_file = "coords"
batch_frac_coord_file = "frac_coords"
batch_cp2k_template = "cp2k_template.inp"
batch_sampler_pyscript = "ucs_sample_batch.py"
batch_sample_outfile = "ucs_sample_batch.out"

#### Set path to pyscript
batcher_script_path = os.path.realpath(sys.argv[0])
batcher_script_dir = os.path.dirname(batcher_script_path)
batcher_pyscript_path = os.path.join(batcher_script_dir, batch_sampler_pyscript)

#########################################################################


#########################################################################
######## HERE BEGINS PREPROCESSING AS DONE IN grid_gen_script.py ########
#########################################################################
### Method section is redundant- actually no, just needed to be passed to batches

#method_list = ['pbe', 'lammps_lj', 'lammps_lj_coul', 'ff_boulfelfel', 'ff_boulfelfel_buck', 'ff_garcia_sanches', 'ase_lj', '54189', '73679']
method_list = ['Li13Si4_test']

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
parser.add_argument('--rc', type=float, action='store', default=0.0, help="Specify a cutoff radius in Å that should be excluded from the sampled volume around each atom in the host framework structure.")
parser.add_argument('--vdw', type=float, action='store', default=0.0, help="Specify a fraction of the van der Waals radius that should be excluded from the sampled volume around each atom in the host framework structure.")
parser.add_argument('--cp2k_q', '--cp2k-total-charge', type=str, action='store', default=None, help="Specify the total charge of the structure that is to be sampled. Observe that this is for the full, final structure that is actually sampled, not the input structure. This charge is for instance passed to CP2K DFT calculators to set the number of electrons. Options: int - to set the charge explicitly, \"auto\" for automatic determination based of cp2k_aq, None for no determination or setting of charge (e.g. if present in cp2k input template already).")
#parser.add_argument('--tc', '--total-charge', type=float, action='store', default=None, help="Specify the total charge of the structure that is to be sampled. Observe that this is for the final structure that is actually sampled, not the input structure. This charge is for instance passed to CP2K dft calculators to set the number of electrons.")
parser.add_argument('--cp2k_aq', '--cp2k-atom-charge', type=int, action='store', default=1, help="Specify the charge of the sampling atom/ion, for the purpose of determining the CP2K total charge. This charge is only used in the automatic determination of the charge that are passed to CP2K DFT calculators. Thus, if cp2k total charge is given as an integer, or if cp2k is not used for sampling, this input is redundant.")

parser.add_argument('--cp2k_template', type=str, action='store', help="Specify the cp2k input file template. Only important for methods using a cp2k input template.")

parser.add_argument('--cp2k_wfn_mode', type=str, action='store', default=None, choices=[None, 'off', 'on', 'seq'], help="Specify the cp2k wfn guess mode, i.e. how the SCF_GUESS is set. By default, it is off. \"on\" or \"seq\" uses a simple sequential strategy, in which each subsequent calculation uses the resulting wfn file from the previous calculation. Only relevant for methods using DFT cp2k methods.")

parser.add_argument('--cp2k_basis_set', type=str, action='store', default=None, help="Specify the cp2k basis set, which will be applied to all atoms. By default, it is None, and in this case the sampler will not tamper with the basis set option for any atoms - it has to be manually specified in the template, if applicable. Only relevant for methods using DFT cp2k methods. Can be e.g. DVZP-MOLOPT-SR-GTH, DZV-MOLOPT-GTH, DZV-MOLOPT-SR-GTH.")

parser.add_argument('--cp2k_pseudo_pot', type=str, action='store', default=None, help="Specify the cp2k pseudo potential, which will be applied to all atoms. By default, it is None, and in this case the sampler will not tamper with the pseudo potential  option for any atoms - it has to be manually specified in the template, if applicable. Only relevant for methods using DFT cp2k methods. Can be e.g. GTH-PBE, or auto, which enables automatic selection within the ase cp2k calculator.")

parser.add_argument('--cp2k_print_level', type=str, action='store', default="MEDIUM", choices=[None, "SILENT", "LOW", "MEDIUM", "HIGH", "DEBUG"], help="Specify the global print level in the cp2k output. If None, it has to be specified in the template, if applicable, or the CP2K default will apply. Only relevant for methods using DFT cp2k methods. Default: MEDIUM")

parser.add_argument('--cp2k_shell_reset_freq', type=int, action='store', default=500, help="Specify the frequency with which the CP2K-shell is re-instantiated, i.e. reset. If this frequency is N, the cp2k shell is killed and restarted every N:th calculation. This is too avoid the logical unit error arising from too many io-units assigned in the cp2k fortran source code. Too disable this, either set it to a higher number than the total number of calculations, or set it to a value <=0. In the latter case the program will automatically set it to a high value so that no reset is performed. Default: 500")

# CP2K command argument - REDUNDANT FOR NOW
#parser.add_argument('--cp2k_cmd', '--cp2k-command', type=str, action='store', default=default_cp2k_cmd, help="Specify the CP2K-command that is used in the ASE-CP2K calculator interface to start and run the CP2K program. Default: ")
# Minimum-image-convention cutoff argument
parser.add_argument('--mic-cutoff', type=float, action='store', default=0.0, help="Specify the cut-off used to construct a supercell obeying the minimum image convention.")
parser.add_argument('--ra', action='store_true', help="Specify whether to remove all atoms of the type that is used for sampling from the structure, before doing the sampling.")
parser.add_argument('--conv', '--conventional-cell', action='store_true', help="If given, will check whether unit cell is conventional cell, and if not, determine the conventional cell based on the input structure and then sample the found conventional cell instead.")

parser.add_argument('--nosym', action='store_false', help="Turns off usage of spacegroup symmetry. Default is to apply symmetry to save the number of required calculations.")
parser.add_argument('--sg', type=int, default=None, action='store', help="Manually specify the spacegroup to use for symmetry, by specifying the corresponding index in the gemmi table (as obtained from gemmi.spacegroup_table_itb()). Note that it is the python index, e.g. starting from 0. Default is None, in which case spacegroup will be automatically determined from the structure. Note that it is assumed that the structure already matches this spacegroup, as this bypasses the matching procedure.")
#parser.add_argument('--guc', '--gemmiunitcell', action='store_true', help="If given, will pass unit cell information to gemmi. This is currently only for testing, to see if symmetry will be better handled in certain cases with primitive unitcells. TEST.")
parser.add_argument('--sym-spglib-std', type=str, action='store', default="allow", choices=["off", "allow", "on"], help="Specifies the setting of the spglib standardization, in the processing of structure to match with spacegroup. Can be on, in which it is enforced, or off, in which case it is not applied, and allow, in which case it is applied only if initial spacegroup matching fails.")
parser.add_argument('--sym-change-basis', type=str, action='store', default="allow", choices=["off", "allow", "on"], help="Specifies how to apply the change of basis in symmetry structure standardization, after spglib standardization. Can be on, in which it is enforced, or off, in which case it is not applied, and allow, in which case it is applied only if initial spacegroup matching without it fails. NOTE: Can only be applied after spglib standardization, and thus this setting must be at most that of sym-spglib-std, e.g. off if it is off, and off or allow if it is allow.")
# TODO: add this
#parser.add_argument('--sym-standardize', action='store_true', help="If given, will perform symmetry standardization of structure and grid to match the spacegroup, even if symmetry is not applied.")

parser.add_argument('--midvox', action='store_true', help="Specifies that sampling should be done in the center of voxels rather than in the corners, as is default. Corresponds to shifting the coordinates of the grid points with 0.5 * (1/na, 1/nb, 1/nc), where ni is the number of grid points in direction of unit cell vector i.")


# Add Batcher arguments
batch_args_group = parser.add_argument_group(title="UCS batch arguments", description="Defines and controls the behaviour of the batching procedure, in a UCS batch run.")
batch_args_group.add_argument("--batch-size", type=int, action="store", default=500, help="Sets the batch size, i.e. how many grid point calculations assigned to each batch, at most.")
batch_args_group.add_argument("--jobscript-cmd", type=str, action="store", default=None, choices=[None, "bash", "sbatch"], help="The command used to launch the jobscripts of the batches, e.g. \"bash\" or \"sbatch\".")
batch_args_group.add_argument("--jobscript-template", type=str, action="store", default=None, help="The run script template that is used to generate runscripts to start the calculations. This should be a script that is runnable with the option given in by --jobscript-cmd. It should contain any extra setup needed to run the sampling. The actual line to run the sampler should not be present, as it will be determined and appended to this template.")
batch_args_group.add_argument("--no-run", action="store_true", help="If specified, will not automatically start the jobs after dividing into batches and setting up the calculations.")



# Parse arguments
args = parser.parse_args()

### End of parser definition

# Log raw input definition
def log_raw_input(logfile, args):
    input_dict = args.__dict__
    logfile.write("--------------------------------------------------\n")
    logfile.write("-                   Raw input                    -\n")
    logfile.write("--------------------------------------------------\n")
    for key,item in input_dict.items():
        logfile.write(key+" : "+str(item)+"\n")

    logfile.write("--------------------------------------------------\n")
    logfile.write("-                End of raw input                -\n")
    logfile.write("--------------------------------------------------\n")



# Check batch argument logic
if not args.no_run and not args.jobscript_cmd:
    raise ValueError("ERROR: Need to supply a jobscript command unless --no-run is specified.")


if (not args.jobscript_template or not os.path.isfile(args.jobscript_template)) and not args.no_run:
    raise FileNotFoundError("ERROR: Jobscript template was not given or not found. A valid template file must be given!")

# Then update jobscript_template to contain FULL path
full_jobscript_path = os.path.realpath(args.jobscript_template)
args.jobscript_template = full_jobscript_path

# Also update cp2k_template if it exists to a FULL path
if args.cp2k_template:
    full_cp2k_tmpl_path = os.path.realpath(args.cp2k_template)
    args.cp2k_template = full_cp2k_tmpl_path

if args.batch_size < 1:
    raise ValueError("ERROR: The batch size must be larger than 0.")


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
#print("Input basename: ", input_basename)

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
    #print("Removing sampling atoms", str(atom),"from structure.")
    atoms_temp = remove_nonframework_cations_fancy(init_unitcell, ase.Atom(atom))
    unitcell = atoms_temp
else:
    unitcell = init_unitcell

if use_sym:
    if args.sg:
        # TODO: Change so that this grabs the spacegroup entry from the gemmi-table. -> DONE BUT TEST? AND UPDATE THE CLI DOCS
        # Grab the desired entry from the gemmi table. Note indexation! It is the pythonic index 
        spacegroup = list(gemmi.spacegroup_table_itb())[args.sg]
        #spacegroup = args.sg # OLD! DEPRECATE!
    else:
        # TODO: Change here so that a gemmi.SpaceGroup is obtained -> DONE, but test.
        # TODO: Make CLI arguments to get more control over symmetry adaption here
        _OLD_spacegroup = get_spacegroup(unitcell) # OLD! DEPRECATE. KEEP ONLY FOR PRINTING INFO.

        spgrp_matching_unitcell, spacegroup, basis_change_tuple = symmetry.prepare_matching_structure_and_spacegroup(atoms=unitcell,
                                                                                      clean=True,
                                                                                      wrap=False,
                                                                                      to_primitive=False,
                                                                                      no_idealize=True,
                                                                                      spglib_standardize=args.sym_spglib_std,
                                                                                      change_basis=args.sym_change_basis)
        unitcell = spgrp_matching_unitcell
else:
    # If symmetry should not be used, set spacegroup to None
    spacegroup = None
## Set the number of points in the grid in each dimension (or equivalently, the mesh size)
## and if symmetry is used, make sure grid and spacegroup are compatible
# TODO: If symmetry is used, check/determine nearest shape that is compatible with spacegroup -> DONE, see todo below
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
    
    # Set init spacing from input
    spacing_x, spacing_y, spacing_z = a/nx, b/ny, c/nz

    # There should be a symmetry adaption here probably ... 
    # TODO: Here, check the compatibility -> DONE, but test.
    if use_sym:
        sym_grid_shape = symmetry.find_spacegroup_compatible_gridshape((nx,ny,nz),
                                                      spacegroup,
                                                      a,
                                                      b,
                                                      c,
                                                      search_extent=10)
        nx, ny, nz = sym_grid_shape

    true_spacing = (a/nx, b/ny, c/nz)

    #print('True spacing: ', true_spacing)
    #print('(nx, ny, nz) =', (nx, ny, nz))

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
        
        # TODO: Fix this, i.e. make sure the new version of the grid shape compatibility function works -> DONE, but test!
        if use_sym:
            sym_grid_shape = symmetry.find_spacegroup_compatible_gridshape((nx,ny,nz),
                                                                     spacegroup,
                                                                     a,
                                                                     b,
                                                                     c,
                                                                     search_extent=10)
            nx, ny, nz = sym_grid_shape
            #print("Determined grid spacing compatible with spacegroup.")

        true_spacing = (a/nx, b/ny, c/nz)
        #print('Desired spacing: ', (spacing_x, spacing_y, spacing_z),' True spacing: ', true_spacing)
        #print('(nx, ny, nz) =', (nx, ny, nz))

#print("vdW cutoff factor: ", args.vdw)

### End of parsing of arguments

# Now handle naming of the calculation files: Chosen format: inputfile_method_gridsize_jobname
# Could addversion if spacing is given, to indicate spacing?
#calc_name = '_'.join((input_basename, method, 'x'.join((str(nx), str(ny), str(nz))), jobname)) # For grid size
if args.name:
    calc_name = '_'.join((input_basename, method, 'x'.join(tuple(map(str, np.around(true_spacing,4)))), jobname))
else:
    calc_name = '_'.join((input_basename, method, 'x'.join(tuple(map(str, np.around(true_spacing,4))))))


##### Make the overarching batch calculation directory, and batch log ####

# Note that paths of files read after this point 
# need to be converted to real absolute paths before changing directory!

batch_calc_dir = calc_name + "_batch_run"
os.mkdir(batch_calc_dir)
os.chdir(batch_calc_dir)
batch_wd = os.getcwd()
batch_log = open(batch_wd+"/ucs_batch.log", 'w')

# Write to log
batch_log.write("UCS batch calculation log file\n")
batch_log.write("------------------------------\n")
batch_log.write("".join(["Date and time: ", str(datetime.datetime.now().ctime()), "\n"]))
batch_log.write("".join(["Invoked with: ", " ".join(sys.argv[:]), "\n"]))
batch_log.write("".join(["Batch directory: ", str(batch_wd),"\n"]))
batch_log.write("\n\n")
batch_log.write("UCS batch settings:\n")
batch_log.write("-------------------\n")
batch_log.write("".join(["Batch size: ", str(args.batch_size), "\n"]))
batch_log.write("".join(["Jobscript template: ", str(args.jobscript_template), "\n"]))
batch_log.write("".join(["Jobscript command: ", str(args.jobscript_cmd), "\n"]))
batch_log.write("".join(["--no-run flag: ", str(args.no_run),"\n"]))
batch_log.write("-------------------\n\n")


batch_log.write("Printing output from UCS preprocessing (similar to grid_gen_script.py)\n")
batch_log.write("======================================================================\n")
log_raw_input(batch_log, args)
batch_log.write("\n\n\n")
if args.ra:
    batch_log.write("Removed sampling atoms " + str(atom) + " from structure.\n")


if use_sym:
    batch_log.write('Determined grid shape compatible with spacegroup.\n')
batch_log.write('Initial spacing: ' + str((spacing_x, spacing_y, spacing_z)) + ' True spacing: ' + str(true_spacing) + "\n")
batch_log.write('Input grid shape: (nx, ny, nz) = ' + str((nx, ny, nz)) + ' True grid shape: (nx, ny, nz) = ' + str((nx, ny, nz)) + "\n")


batch_log.write("Radial cutoff [Å]: " + str(args.rc) + "\n")
batch_log.write("vdW cutoff factor: " + str(args.vdw) + "\n")

if use_sym:
    batch_log.write("Symmetry is used.\n")
if use_sym or args.conv:
    batch_log.write("Determined conventional cell for sampling.\n")
if args.midvox:
    batch_log.write("Midvox sampling enabled.\n")



### Set the directories and related/relevant environment vars
#work_dir = '.'
#os.environ['UCS_WORK_DIR'] = work_dir


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
batch_log.write("Force cutoff used to determine supercell size: "+ str(cutoff)+"\n")
num_cells = compute_req_supercells(unitcell, cutoff)


batch_log.write("num_cells in supercell: " + str(num_cells) + "\n")

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
if args.cp2k_q == "auto":
    # Determine charge per unitcell:
    cp2k_charge_per_unitcell = determine_total_cp2k_charge(init_unitcell, 
            unitcell, args.atom, args.cp2k_aq)

    # Then multiply with number of unit cells in supercell,
    # and remove one for each
    cp2k_charge_for_supercell_alt = ma.prod(num_cells) * cp2k_charge_per_unitcell \
                              - (ma.prod(num_cells) - 1) * args.cp2k_aq

    # Computing from supercells now
    batch_log.write("\n")
    batch_log.write("Computing CP2K total charge based on supercells\n")
    cp2k_charge_for_supercell = determine_total_cp2k_charge(
            supercell_from_init_unitcell, supercell_from_unitcell_wo_ions, 
            args.atom, args.cp2k_aq)

    assert cp2k_charge_for_supercell_alt == cp2k_charge_for_supercell, "Total CP2K charge not equal with different methods!!! Check charge determination..."

    # Set args namespace
    args.cp2k_q = cp2k_charge_for_supercell
    


batch_log.write("supercell from uc params: " + str(supercell_from_unitcell_wo_ions.get_cell_lengths_and_angles()) + "\n")
batch_log.write("supercell from uc cell: " + str(supercell_from_unitcell_wo_ions.get_cell()) + "\n")

batch_log.write("\n")
# NOTE: Printing of the new spacegroup information has been updated (in prev. update of batch analzer).
# NOTE: Added debug printing/comparison of old spacegroup determination/information prior to update, and fixed logic, i.e. only when using symmetry
### THIS IS PRINTING OF OLD SYMMETRY INFORMATION: ###
if use_sym:
    batch_log.write("\n")
    batch_log.write("DEBUG: Printing old, deprecated symmetry information below.\n")
    batch_log.write("DEBUG (OLD SYMMETRY INFO): OLD SPACEGROUP, UNITCELL, number and name: " + str(_OLD_spacegroup.no) + ", " + str(_OLD_spacegroup.symbol) + " \n") # OLD, DEPRECATE!
    
    _OLD_GEMMI_spacegroup = gemmi.find_spacegroup_by_number(_OLD_spacegroup) # OLD! Deprecate this
    try:
        _OLD_GEMMI_idx = list(gemmi.spacegroup_table_itb()).index(_OLD_GEMMI_spacegroup)
    except:
        _OLD_GEMMI_idx = "(NOT FOUND)"
    
    SYM_DEBUG_CMPR_old_new_spgrps = symmetry.compare_gemmi_spacegroups(_OLD_GEMMI_spacegroup, spacegroup)
    SYM_DEBUG_CMPR_old_new_ops = symmetry.are_symops_equal(symmetry.gemmi_spacegroup_to_symops_list(_OLD_GEMMI_spacegroup),
                                                           symmetry.gemmi_spacegroup_to_symops_list(spacegroup))
    
    batch_log.write("DEBUG (OLD SYMMETRY INFO): OLD/NEW SPACEGROUP COMPARISONS: groups equal = " + str(SYM_DEBUG_CMPR_old_new_spgrps) + " ;;  operations equal = "  + str(SYM_DEBUG_CMPR_old_new_ops)+ "\n") # OLD, DEPRECATE!
    batch_log.write("DEBUG (OLD SYMMETRY INFO): OLD SPACEGROUP INFO, gemmi: table index " + str(_OLD_GEMMI_idx) + " ;; "  + str(_OLD_GEMMI_spacegroup)+ "\n") # OLD, DEPRECATE!
    
    batch_log.write("\n")
### 
#batch_log.write("DEBUG (OLD SYMMETRY INFO): Spacegroup, unitcell: " + str(ase.spacegroup.get_spacegroup(unitcell, 1.0e-6)) + "\n") # OLD, DEPRECATE!
#batch_log.write("DEBUG (OLD SYMMETRY INFO): Spacegroup, supercell: " + str(ase.spacegroup.get_spacegroup(supercell_from_unitcell_wo_ions, 1.0e-6)) + "\n") # OLD, DEPRECATE
if use_sym:
    if spacegroup is not None:
        spacegroup_gemmi_table_idx = list(gemmi.spacegroup_table_itb()).index(spacegroup)
    else:
        spacegroup_gemmi_table_idx = "Not found."
        spacegroup = "NONE"
    batch_log.write("Spacegroup, unitcell:  table index " + str(spacegroup_gemmi_table_idx) + " ; " + str(spacegroup) + "\n")
    
    try:
        tmp_supercell_spacegroup = symmetry.find_matching_gemmi_spacegroup_from_spglib_dataset(
                                      symmetry.spglib.get_symmetry_dataset(
                                      symmetry.atoms_to_spglib_cell_tuple(supercell_from_unitcell_wo_ions)
                                      )
                                      )
        assert len(tmp_supercell_spacegroup) == 1 or (len(tmp_supercell_spacegroup) >= 1 
                                                      and all([gsg.number == 68 for gsg in tmp_supercell_spacegroup]))
        
        tmp_supercell_spacegroup = tmp_supercell_spacegroup[0]
        tmp_supercell_spacegroup_idx = list(gemmi.spacegroup_table_itb()).index(tmp_supercell_spacegroup)
    except:
        tmp_supercell_spacegroup = "Not found."
        tmp_supercell_spacegroup_idx = "NONE"
    
    batch_log.write("Spacegroup, supercell: table index " + str(tmp_supercell_spacegroup_idx) + " ; " + str(tmp_supercell_spacegroup) + "\n")
else:
    batch_log.write("Symmetry is not used. Spacegroup is not determined." + "\n")

# "Middle" unit cell (we don't need this though, since pbc)
#uc_indices = [int((i//2)-1) if i % 2 == 0 else int((i-1)//2) for i in num_cells]
uc_indices = (0, 0, 0)


# Now generate grid for unitcell:
unitcell_ucs = sample.UnitCellSampler(unitcell) # DONE

if (args.midvox and use_sym):
    batch_log.write("WARNING: Both symmetry and midvoxel sampling enabled. Their compatibility has NOT been properly tested!!!\n")

unitcell_grid, unitcell_included = unitcell_ucs.generate_grid_vectors((nx, ny, nz), 
                                                                      cutoff_radii=args.rc, 
                                                                      vdw_scale=args.vdw, 
                                                                      midvox=args.midvox) # DONE

unitcell_grid = unitcell_grid.reshape(-1,3) # DONE
unitcell_included = unitcell_included.reshape(-1) # DONE

# Convert to fractional coordinates, and convert to supercell grid

batch_log.write("Shape check, grid: " + str(np.array(unitcell_grid).shape) + " " + str(np.array(unitcell_grid).T.shape)+"\n") # DONE
unitcell_frac_grid = np.linalg.solve(np.array(unitcell.get_cell()).T, np.array(unitcell_grid).T).T # DONE

batch_log.write("unitcell frac grid shape: " + str(unitcell_frac_grid.shape) + "\n") #DONE
# Full grid:
supercell_frac_grid_full = unitcell_to_supercell_frac_coords(unitcell_frac_grid, num_cells, unitcell_ind=(0,0,0)) #DONE

# Then filtered:
supercell_frac_grid = unitcell_to_supercell_frac_coords(unitcell_frac_grid[unitcell_included], num_cells, unitcell_ind=(0,0,0)) #DONE

# Convert to cartesian supercell grid
# Full grid:
supercell_cart_grid_full = supercell_frac_grid_full @ np.array(supercell_from_unitcell_wo_ions.get_cell())

# Then filtered
supercell_cart_grid = supercell_frac_grid @ np.array(supercell_from_unitcell_wo_ions.get_cell()) # DONE



assert (supercell_cart_grid_full[unitcell_included] == supercell_cart_grid).all(), "Cartesian supercell grids are wrong somehow!"
assert (supercell_frac_grid_full[unitcell_included] == supercell_frac_grid).all(), "Fractional supercell grids are wrong somehow!"

#######################################################################


batch_log.write("Spacegroup input: " + str(args.sg) + "\n")
batch_log.write("Symmetry: " + str(use_sym) + "\n")


# CP2K from input - This is kept for the exception hanlding!
if method in cp2k_dft_methods.keys() and method == "cp2k_calculator_from_input":
    # First: Check the cp2k_template is set, and that it is a file.
    if args.cp2k_template is None:
        raise Exception("Method " + str(method) + " requires a cp2k template \
              input file to be specified through the --cp2k_template argument.")

    cp2k_template_file  = Path(args.cp2k_template) 

    if not cp2k_template_file.exists():
        raise Exception("CP2K template input file " + str(cp2k_template_file) 
                       + " was not found!")


    # TODO: This has been updated in cp2k_calculators.py. Testing needs to be done
    if args.cp2k_q is not None:
        batch_log.write("Passing total charge = "+ str(int(args.cp2k_q))+ " to CP2K calculator\n")
    else:
        batch_log.write("--cp2k_q is None - charge is not passed, but that present"
              + " in the inputfile template is used if present.\n")


    # Manage wfn file usage:
    if args.cp2k_wfn_mode in ['on', 'seq']:
        wfn_restart = True
        batch_log.write("cp2k_wfn_mode = " + str(args.cp2k_wfn_mode) + ": Reusing consecutive wfn-files as initial guess.\n")
    else:
        wfn_restart = False
        batch_log.write("cp2k_wfn_mode = " + str(args.cp2k_wfn_mode) + ": Disabled. Wfn files are not used.\n")


    if args.cp2k_basis_set is not None:
        batch_log.write("CP2K basis set specified: " + str(args.cp2k_basis_set) + "\n")
    else:
        batch_log.write("No CP2K basis set provided (= None): Basis set taken to be specified in input template or implicitly in the method.\n")
    
    if args.cp2k_pseudo_pot is not None:
        batch_log.write("CP2K pseudo potential specified: " + str(args.cp2k_pseudo_pot) + "\n")
    else:
        batch_log.write("No CP2K pseudo potential provided (= None): Pseudo potential taken to be specified in input template or implicitly in the method.\n")

    if args.cp2k_print_level is not None:
        batch_log.write("CP2K print level: " + str(args.cp2k_print_level) + "\n")
    else:
        batch_log.write("CP2K print level not set in sampler.\n")

    # Manage and control the cp2k shell reset frequency setting
    assert isinstance(args.cp2k_shell_reset_freq, int), "ERROR: Shell reset frequency must be integer!!!"
    if args.cp2k_shell_reset_freq <= 0:
        batch_log.write("CP2K-shell frequency is <= 0. Setting it to 1 + total number of grid points, to disable it.\n")
        args.cp2k_shell_reset_freq = int(nx + ny + nz) + 1 # Set it to more than the total number of grid points
        batch_log.write("CP2K-shell reset frequency was set to: " + str(args.cp2k_shell_reset_freq)+"\n")
        batch_log.write("This should be 1 + total number of grid points!\n")

    batch_log.write("Shell reset frequency:"+str(args.cp2k_shell_reset_freq)+"\n")

batch_log.write("===================== END OF PREPROCESSING OUTPUT =============================\n\n\n")


################################################################################
########## DETERMINING GRID POINT COORDINATES AND INDICES, AND BATCHING ########
################################################################################

# Notes: In the preprocessing, we get:
# unitcell_grid             =  The coordinates of the full grid in cartesian coordinates in a single unit cell, shape (-1, 3)
# unitcell_included         =  The mask that give the included grid points, that are not excluded by vdw criterion, shape (-1)
# unitcell_frac_grid        =  Like unitcell_grid but in fractional coordinates of the unitcell 
# supercell_frac_grid_full  =  The full grid points but in fractional coordinates of the determined supercell 
# supercell_cart_grid_full  =  Like supercell_frac_grid_full but in cartesian coordinates
# supercell_frac_grid       =  The INCLUDED grid points but in fractional coordinates of the determined supercell 
# supercell_cart_grid       =  Like supercell_frac_grid but in cartesian coordinates

# We only want the coordinates to be sampled
# We could write out all of them too?
# Here however, it is better to just collect the ones that are sampled, that is the intent

# Initialize variables:
cart_grid_coord_list = []
frac_grid_coord_list = []
indices_list = []

# code adapted from sample.py:

def get_cartesian_coords(fractional_coords, cell_vectors) -> np.ndarray:
    return np.dot(fractional_coords, cell_vectors)


def get_fractional_coords(cartesian_coords, cell_vectors) -> np.ndarray:
    return np.dot(cartesian_coords, np.linalg.inv(cell_vectors))


# intialise counters
n_excluded_total = 0
n_included_total = 0
n_exploited_symmetry = 0
n_calculations_total = 0
total_points = np.size(unitcell_included)

# Counters for the specific spherical cutoff types
n_radial_cutoff_included = np.count_nonzero(unitcell_ucs.cutoff_included)
n_radial_cutoff_excluded = np.count_nonzero(np.logical_not(unitcell_ucs.cutoff_included))
n_vdw_included = np.count_nonzero(unitcell_ucs.vdw_included)
n_vdw_excluded = np.count_nonzero(np.logical_not(unitcell_ucs.vdw_included))

batch_log.write("Determining actual number of required calculations, and corresponding grid points...\n")

# TODO: Here we need to update the setting of the spacegroup!!! -> DONE, but test!
if use_sym:
    bool_grid = gemmi.Int8Grid(nx,ny,nz)
    #bool_grid.spacegroup = gemmi.find_spacegroup_by_number(spacegroup) OLD! Deprecate this
    bool_grid.spacegroup = spacegroup # New! This is gemmi.SpaceGroup instance

    # NOTE: Added for constructing/printing of symmetry cube
    symID_grid = gemmi.FloatGrid(nx,ny,nz)
    symID_grid.spacegroup = spacegroup


    for idx, grid_point in enumerate(supercell_cart_grid_full): # originally looped over unitcell_grid
        if not unitcell_included[idx]:
            n_excluded_total += 1
            continue
        n_included_total += 1

        # old grid point index getter:
        #grid_point_index = get_fractional_coords(
        #    grid_point, unitcell.cell[:])
        #grid_point_index[0] = grid_point_index[0]*nx
        #grid_point_index[1] = grid_point_index[1]*ny
        #grid_point_index[2] = grid_point_index[2]*nz
        #grid_point_index = np.array(
        #    np.around(grid_point_index), dtype=int)

        # New using numpy
        grid_point_index_new = np.hstack(np.unravel_index(idx, (nx,ny,nz)))

        # Sanity check
        #assert (grid_point_index == grid_point_index_new).all(), "New and old grid point indices does not match!!!"

        # check if sampled before:
        if bool_grid.get_value(*grid_point_index_new):
            n_exploited_symmetry += 1
            continue

        # Store index and cartesian grid point coordinate, since we sample it
        indices_list.append(grid_point_index_new)
        cart_grid_coord_list.append(grid_point)
        frac_grid_coord_list.append(supercell_frac_grid_full[idx,:])

        assert (grid_point == get_cartesian_coords(supercell_frac_grid_full[idx,:], supercell_from_unitcell_wo_ions.cell[:])).all(), "Error: Cartesian and fractional coord mismatch!?"

        #atoms = self.atoms + ase.Atom(atom, grid_point)

        bool_grid.set_value(*grid_point_index_new, 1)
        bool_grid.symmetrize_max()
        n_calculations_total += 1

        # NOTE: Added for constructing/printing of symmetry cube
        symID_grid.set_value(*grid_point_index_new, n_calculations_total)

    # NOTE: Moved line below, symmetrization of symmetry cube, to outside loop, as it is not needed to be inside and it is faster
    symID_grid.symmetrize_max()
    indices_list = np.array(indices_list, dtype=int)
    cart_grid_coord_list = np.array(cart_grid_coord_list, dtype=np.float64)
    frac_grid_coord_list = np.array(frac_grid_coord_list, dtype=np.float64)

    # NOTE: Printing symmetry cube here:
    # print symmetry information
    sym_name = calc_name+'_symInfo'
    cube_filename = ".".join((sym_name, 'cube'))
    #with open(cube_filename, 'w') as fp:
    #    write_cube(fp, unitcell, data=np.array(symID_grid))
    data = np.array(symID_grid, dtype='int')
    symmetry_cube = open(cube_filename,'w')
    symmetry_cube.write('symmetry information file, units: Angstrom?\n')
    symmetry_cube.write('--------------------------------\n')
    symmetry_cube.write(''.join(["{:d}".format(1),' ','0.000000 0.000000 0.000000\n']))
    symmetry_cube.write(''.join(["{:d}".format(nx),' ',"{:1.6f}".format(true_spacing[0]),' ',"{:1.6f}".format(0.),' ',"{:1.6f}".format(0.),'\n']))
    symmetry_cube.write(''.join(["{:d}".format(ny),' ',"{:1.6f}".format(0.),' ',"{:1.6f}".format(true_spacing[1]),' ',"{:1.6f}".format(0.),'\n']))
    symmetry_cube.write(''.join(["{:d}".format(nz),' ',"{:1.6f}".format(0.),' ',"{:1.6f}".format(0.),' ',"{:1.6f}".format(true_spacing[2]),'\n']))
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                symmetry_cube.write(''.join(["{:d}".format(data[i,j,k]),'\n']))
    symmetry_cube.close()

else:
    # Here only the previous vdw exclusion is performed
    indices_list = np.argwhere(unitcell_included.reshape(nx,ny,nz))
    indices_list_2 = np.hstack(np.unravel_index(np.argwhere(unitcell_included), (nx,ny,nz))) # Just to see if I understand this correctly
    assert (indices_list == indices_list_2).all(), "Error: Indices list with diff methods not the same!!!"
    cart_grid_coord_list = supercell_cart_grid
    frac_grid_coord_list = supercell_frac_grid

    # Compute counts
    n_excluded_total = total_points - np.count_nonzero(unitcell_included)
    n_included_total = np.count_nonzero(unitcell_included)
    n_exploited_symmetry = 0
    n_calculations_total = n_included_total

# Info print to batch-log:
batch_log.write("Total grid points: "+str(total_points)+"\n")
batch_log.write("Points included by total spherical cutoffs: "+str(n_included_total)+"\n")
batch_log.write("Points excluded by total spherical cutoffs: "+str(n_excluded_total)+"\n")
batch_log.write("Points excluded by radial cutoff: "+str(n_radial_cutoff_excluded)+"\n")
batch_log.write("Points excluded by vdW cutoff: "+str(n_vdw_excluded)+"\n")
batch_log.write("Points excluded by symmetry: "+str(n_exploited_symmetry)+"\n")
batch_log.write("Actual number of calculations: "+str(n_calculations_total)+"\n")
batch_log.write("\n")
#batch_log.write("Grid point coordinates of points\n")
# 

### NOW DO THE BATCHING ###
batch_log.write("Determining the number of batches...\n")
# Determine number of batches 
n_batches = int(n_calculations_total) // args.batch_size

if not (n_batches * args.batch_size == n_calculations_total):
    n_batches += 1

assert n_batches * args.batch_size >= int(n_calculations_total), "ERROR: Number of batches times batch size is somehow less than total number of required calculations!!!"

batch_log.write("Batch size: "+str(args.batch_size)+"\n")
batch_log.write("Number of calculations: "+str(n_calculations_total)+"\n")
batch_log.write("Number of batches determined: "+str(n_batches)+"\n")
batch_log.write("\n")


#################################################
########     Prepare the runscript!     #########
#################################################

# Run script should be adapted to run the calculation eiher locally or in whatever environment
# What the batcher is supplying is the template lines for the calculation dir, and 
# setting of the passed parameters that are 

# We need to pass the relevant arguments to the script running the sampling
### First, store arguments:
full_arguments_list = sys.argv[1:]
# Need to pass: file - determined, method - pass, name - None?, filename with coords - determined, cp2k_q - calculated, cp2k_template - determined
#"structure.cif"
#"coords.txt"
#"cp2k_template.inp"
# method - from method here? so: args.method
# other cp2k arguments?!?
## Needs to be constructed
batch_log.write("Determining and constructing arguments to pass to batch grid sampler\n")

# Determine format for the passed coordinate file
if batch_coord_npy_format:
    coord_ext = ".npy"
else:
    coord_ext = ".txt"


args_to_pass = ["../structure.cif",
                args.method,
                "--coords-from-file", batch_coord_file+coord_ext,
                "-a", args.atom]

# Then determine and pass required cp2k arguments:  
if args.cp2k_template is not None:
    args_to_pass.extend(["--cp2k_template", "../cp2k_template.inp"])
    batch_log.write("Passing cp2k_template.\n")
    
if args.cp2k_q is not None:
    args_to_pass.extend(["--cp2k_q", str(int(args.cp2k_q))])
    batch_log.write("Passing cp2k_q.\n")

# TODO: Remove passing of cp2k_aq - it isn't needed!
#if args.cp2k_aq is not None and ("--cp2k_aq" in full_arguments_list or "--cp2k-atom-charge" in full_arguments_list):
#    args_to_pass.extend(["--cp2k_aq", str(args.cp2k_aq)])
#    batch_log.write("Passing cp2k_aq (cp2k-atom-charge).\n")

if args.cp2k_wfn_mode is not None:
    args_to_pass.extend(["--cp2k_wfn_mode", str(args.cp2k_wfn_mode)])
    batch_log.write("Passing cp2k_wfn_mode.\n")

# TODO; Add passing of cp2k_basis_set, pseudo pot, and print level
if args.cp2k_basis_set is not None:
    args_to_pass.extend(["--cp2k_basis_set", str(args.cp2k_basis_set)])
    batch_log.write("Passing cp2k_basis_set.\n")

if args.cp2k_pseudo_pot is not None:
    args_to_pass.extend(["--cp2k_pseudo_pot", str(args.cp2k_pseudo_pot)])
    batch_log.write("Passing cp2k_pseudo_pot.\n")

if args.cp2k_print_level is not None and "--cp2k_print_level" in full_arguments_list:
    args_to_pass.extend(["--cp2k_print_level", str(args.cp2k_print_level)])
    batch_log.write("Passing cp2k_print_level.\n")

# TODO: add passing of cp2k shell reset frequency here
if args.cp2k_shell_reset_freq is not None and "--cp2k_shell_reset_freq" in full_arguments_list:
    args_to_pass.extend(["--cp2k_shell_reset_freq", str(args.cp2k_shell_reset_freq)])
    batch_log.write("Passing cp2k_shell_reset_freq.\n")



formatted_args_to_pass = " ".join(args_to_pass)
batch_log.write("\n")
batch_log.write("Full argument list that is passed: "+str(args_to_pass)+"\n")
batch_log.write("Formatted argument string: "+formatted_args_to_pass+"\n")



batch_log.write("Making jobscript from template... \n")
# Read the template
with open(args.jobscript_template, 'r') as jsf:
    jobscript_tmp = "".join(jsf.readlines())


# Add necessary lines
jobscript_tmp += "\n\n"
jobscript_tmp += "# Automatically added lines by the ucs batcher:\n"
#jobscript_tmp += "# Here for testing and debugging:\n"   # Debug CHANGE HERE REMOVE THESE
#jobscript_tmp += "echo Hello from BATCH {batchnr}\n"     # Debug
#jobscript_tmp += "exit 0 \n"                             # Debug
jobscript_tmp += "export UCS_CALCULATION_DIR=./calc\n"
#jobscript_tmp += "export UCS_WORK_DIR=.\n"
jobscript_tmp += "python -u {pyscript} {sampler_args} > {outfile}\n".format(pyscript="../"+batch_sampler_pyscript,
                                                             sampler_args=formatted_args_to_pass,
                                                             outfile=batch_sample_outfile)


jobscript_tmp += "\n# End of automatically added lines \n"


#############################################################################
######## Create directory hierarchy for batch runs, move/write files ########
#############################################################################

# First place files in the batch run main directory
# E.g.:
# cp2k template, cif file, python script symlink


# Write cif containing processed structure to be sampled:
batch_log.write("Writing processed structure to batch directory...\n")
ase.io.write(batch_structure, supercell_from_unitcell_wo_ions)

# Create symlink to the sampling script that is called to compute each batch
batch_log.write("Placing link to sampler in batch directory...\n")
os.symlink(batcher_pyscript_path, os.path.join(batch_wd, batch_sampler_pyscript))

# Copy the cp2k_template if necessary
if args.cp2k_template is not None:
    batch_log.write("Copying cp2k template to batch directory...\n")
    shutil.copy(args.cp2k_template, os.path.join(batch_wd, batch_cp2k_template))
#
batch_log.write("Done transfering files to main batch directory.\n\n")

batch_log.write("Generating batch subdirectory hierarchy and preparing individual batch calculations...\n")
# Now iteratively make each batch directory
# and place relevant files there
for batchnr in range(n_batches):
    batch_subdir_name = batch_wd + "/BATCH" + str(batchnr)
    batch_log.write("Creating "+str(batch_subdir_name)+" ...\n")
    os.mkdir(batch_subdir_name)
    os.chdir(batch_subdir_name)

    if batchnr == (n_batches - 1):
        # If last batch!
        batch_slice = slice(batchnr*args.batch_size, None)
    else:
        batch_slice = slice(batchnr*args.batch_size, (batchnr+1)*args.batch_size)

    batch_indices = np.copy(indices_list[batch_slice])
    batch_cart_grid_coords = np.copy(cart_grid_coord_list[batch_slice])
        
    
    # NOTE: Numpy can save with save(), .tofile() and savetxt(). I believe documentation suggests save or savetxt should be used?
    # Save yields .npy files ; savetxt .txt ; tofile depends
    if batch_coord_txt_format:
        np.savetxt(batch_indices_file+".txt", batch_indices)
        np.savetxt(batch_coord_file+".txt", batch_cart_grid_coords)
    
    if batch_coord_npy_format:
        np.save(batch_indices_file+".npy", batch_indices)
        np.save(batch_coord_file+".npy", batch_cart_grid_coords)

    batch_log.write("Coordinates written.\n")
    
    # write jobscript to the folder
    with open(batch_jobscript_file, "w") as js:
        js.write(jobscript_tmp) # What it should be
        #js.write(jobscript_tmp.format(batchnr=batchnr)) # Test/Debug CHANGE HERE
    batch_log.write("Jobscript written.\n")
    batch_log.write("BATCH"+str(batchnr)+" done.\n\n")


    # Change back to create next batch dir
    os.chdir(batch_wd)

print("Batch run setup completed.")
batch_log.write("\nBatch run setup completed.\n")

##########################################################
######## Lastly run the calculations if desired:  ########
##########################################################

if args.no_run: # Intuitive
    # Don't run ...
    print("Argument --no-run specified: The batch sampling jobs are not started.")
    batch_log.write("Argument --no-run specified: The batch sampling jobs are not started.\n")
else:
    print("Starting the batch grid sampling jobs...\n")
    batch_log.write("\nStarting the batch grid sampling jobs...\n")
    batch_log.write("The jobscript command given: " + args.jobscript_cmd + "\n")
    run_cmd = [str(args.jobscript_cmd), "./" + batch_jobscript_file]
    batch_log.write("Starting each batch sampling through the following command from the respective subdirectory: \"" + " ".join(run_cmd) + "\"\n")
    
    #run ... 
    for batchnr in range(n_batches):
        batch_subdir_name = batch_wd+"/BATCH"+str(batchnr)
        os.chdir(batch_subdir_name)
        #batch_log.write("Simulating running batch: " + str(batchnr) + " from directory: " + str(os.getcwd()) + "\nShould be: " + str(batch_subdir_name)+"\n")
        batch_log.write("Running batch " + str(batchnr) + " from: " + str(os.getcwd())+ "\n")
        #pwd_proc = subprocess.run(["pwd"], check=True)
        #print(pwd_proc.stdout)
        subprocess.run(run_cmd) # Here it is run!
        batch_log.write("BATCH "+str(batchnr)+" submitted.\n")

        os.chdir(batch_wd)
    
    print("Batch run submission completed.")
    batch_log.write("\n\nBatch run submission completed.\n")


# Finally close the logfile
batch_log.close()
