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
from unitcellsampling.autocreated_methods import *

# Methods for 1075 of the structures in the rest of the dataset:
from unitcellsampling import autocreated_methods_structures_mcloud_rest_WO_duplabels_nolog_nosubdircalc 

# predefined by Ben
from unitcellsampling.energy_calculator import cp2k_dft



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

# Last update: 

### TODO: Create input argument handling, to make usage much smoother. Some particular changes:
     # -  TODO: make the  program either recognize the filename from the input file automatically, and name the "project" accordingly
     #          or make the desired project name an input so that the variables all are set in conjuction /H
     # -  TODO: Make the subdirectory for the gridsampling also be project name-derived, and check so that it cannot overide anything by default /H


# TODO: Actually: Move the setting of a project name and setting of environment variables to an exterior runscript - MAYBE, but arguably could be set here too.
# TODO: Move the definition the the method arguments to this list here, for easier editing
method_list = ['pbe', 'lammps_lj', 'lammps_lj_coul', 'ff_boulfelfel', 'ff_boulfelfel_buck', 'ff_garcia_sanches', 'ase_lj', '54189', '73679']

dft_methods = {'cp2k_dft':cp2k_dft}

automethods = {'A54189':m54189_auto,
  'A54209':m54209_auto,
  'A54297':m54297_auto,
  'A54449':m54449_auto,
  'A54683':m54683_auto,
  'A54837':m54837_auto,
  'A54865':m54865_auto,
  'A54879':m54879_auto,
  'A54884':m54884_auto,
  'A55184':m55184_auto,
  'A55319':m55319_auto,
  'A55983':m55983_auto,
  'A56568':m56568_auto,
  'A56627':m56627_auto,
  'A57347':m57347_auto,
  'A57382':m57382_auto,
  'A57448':m57448_auto,
  'A57608':m57608_auto,
  'A57644':m57644_auto,
  'A57761':m57761_auto,
  'A59631':m59631_auto,
  'A59632':m59632_auto,
  'A59715':m59715_auto,
  'A59948':m59948_auto,
  'A60450':m60450_auto,
  'A61111':m61111_auto,
  'A61237':m61237_auto,
  'A61329':m61329_auto,
  'A63257':m63257_auto,
  'A63279':m63279_auto,
  'A63600':m63600_auto,
  'A63663':m63663_auto,
  'A63922':m63922_auto,
  'A64383':m64383_auto,
  'A64396':m64396_auto,
  'A64754':m64754_auto,
  'A64867':m64867_auto,
  'A65350':m65350_auto,
  'A65898':m65898_auto,
  'A66685':m66685_auto,
  'A66876':m66876_auto,
  'A66985':m66985_auto,
  'A67415':m67415_auto,
  'A67418':m67418_auto,
  'A67779':m67779_auto,
  'A67785':m67785_auto,
  'A67806':m67806_auto,
  'A67999':m67999_auto,
  'A68042':m68042_auto,
  'A68103':m68103_auto,
  'A68120':m68120_auto,
  'A68242':m68242_auto,
  'A68674':m68674_auto,
  'A68894':m68894_auto,
  'A68960':m68960_auto,
  'A69116':m69116_auto,
  'A69642':m69642_auto,
  'A70003':m70003_auto,
  'A70442':m70442_auto,
  'A70681':m70681_auto,
  'A71368':m71368_auto,
  'A71681':m71681_auto,
  'A71925':m71925_auto,
  'A71990':m71990_auto,
  'A72018':m72018_auto,
  'A72372':m72372_auto,
  'A72631':m72631_auto,
  'A72652':m72652_auto,
  'A73298':m73298_auto,
  'A73679':m73679_auto,
  'A73766':m73766_auto,
  'A74395':m74395_auto,
  'A74865':m74865_auto,
  'A75073':m75073_auto,
  'A75138':m75138_auto,
  'A75630':m75630_auto,
  'A75961':m75961_auto,
  'A76276':m76276_auto,
  'A76595':m76595_auto,
  'A77893':m77893_auto,
  'A77899':m77899_auto,
  'A78086':m78086_auto,
  'A78099':m78099_auto,
  'A78355':m78355_auto}


# Adding methods for 1075 structures in the larger dataset
mcloud_rest_full_dict = autocreated_methods_structures_mcloud_rest_WO_duplabels_nolog_nosubdircalc.__dict__

#mcloud_rest_method_dict = {}
for key in mcloud_rest_full_dict.keys():
    if key[0] == "m" and key[-4:]=="auto":
        automethods["A"+key.strip("m_auto")] = mcloud_rest_full_dict[key]

# Extending the methods list:
method_list.extend(automethods.keys())

method_list.extend(dft_methods.keys())

### Definition and parsing of arguments
parser = argparse.ArgumentParser(description='Energy sampling of a (periodic) solid system with an added atom/ion on the grid.') 
parser.add_argument('file', metavar='filename', type=str, action='store', help="Name of the structure file containing the unit cell structure to be sampled. Formats: cif, cube, xsf, ...")
#parser.add_argument('-m', '--method', type=str, action='store', default='ff_boulfelfel', choices=['pbe', 'ff', 'ff_boulfelfel', 'ff_boulfelfel_buck'], help="Method to calculate the energy during grid sampling.") # should method be optional or not? No, probably not optional? Or maybe optinal with default?
parser.add_argument('method', type=str, action='store', choices=method_list, help="Method to calculate the energy during grid sampling.")
parser.add_argument('-n','--name', metavar='jobname', type=str, action='store', default=None, help="Desired name for the calculation (applied to generated output files, directories, etc.).")
parser.add_argument('-w', '--wfn', type=str, action='store', default=None, help="Specify the initial wfn-file for a DFT calculation.")
parser.add_argument('-a', '--atom', type=str, action='store', default='Li', help="Specify the atom/ion used for sampling.")
parser.add_argument('-g', '--grid', type=int, action='store', default=[10], nargs='+', help="Specify the number of grid points in each dimension (or cubic grid) (mutually exclusive with \"--space\").")
parser.add_argument('-s', '--space', type=float, action='store', default=None, nargs='+', help="Specify the spacing between the grid points in each dimension (mutually exclusive with \"--grid\").")
parser.add_argument('--vdw', type=float, action='store', default=0.0, help="Specify the fraction of the van der Waals radius that should be excluded from the sampled volume around each atom in the host frameworkstructure.")
parser.add_argument('--nosym', action='store_false', help="Turns off usage of spacegroup symmetry. Default is to apply symmetry to save the number of required calculations.")
parser.add_argument('--ra', action='store_true', help="Specify whether to remove all atoms of the type that is used for sampling from the structure, before doing the sampling.")
parser.add_argument('--sg', type=int, default=None, action='store', help="Manually specify the spacegroup to use for symmetry. Default is None, in which case spacegroup will be automatically determined from the structure.")
parser.add_argument('--guc', '--gemmiunitcell', action='store_true', help="If given, will pass unit cell information to gemmi. This is currently only for testing, to see if symmetry will be better handled in certain cases with primitive unitcells. TEST.")
parser.add_argument('--conv', '--conventional-cell', action='store_true', help="If given, will check whether unit cell is conventional cell, and if not, determine the conventional cell based on the input structure and then sample the found conventional cell instead.")
parser.add_argument('--midvox', action='store_true', help="Specifies that sampling should be done in the center of voxels rather than in the corners, as is default. Corresponds to shifting the coordinates of the grid points with 0.5 * (1/na, 1/nb, 1/nc), where ni is the number of grid points in direction of unit cell vector i.")

args = parser.parse_args()

### End of parser definition

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
lgps = ase.io.read(Path(indir, infile))

# Change the unitcell to the conventional cell if symmetry 
# is used or if conv argument is specified
# This is a bit convoluted as we may be able to use ase
# directly I have realized now...
if use_sym or args.conv:
    if not symmetry.is_conventional_cell(lgps):
        sg_analyzer = SpacegroupAnalyzer(AseAtomsAdaptor.get_structure(lgps)) 
        conv_cell = sg_analyzer.get_conventional_standard_structure()
        lgps = AseAtomsAdaptor.get_atoms(conv_cell)

#atoms = ase.io.read(input_file) # Like so? 

if args.ra:
    print("Removing sampling atoms", str(atom),"from structure.")
    atoms_temp = remove_nonframework_cations_fancy(lgps, ase.Atom(atom))
    lgps = atoms_temp

if use_sym:
    if args.sg:
        spacegroup = args.sg
    else:
        spacegroup = get_spacegroup(lgps) 

## Set the number of points in the grid in each dimension (or equivalently, the mesh size)
# TODO: If symmetr is used, check/determine nearest shape that is compatible with spacegroup
if not args.space:
    if len(args.grid) == 1:
        nx,ny,nz = args.grid * 3
    elif len(args.grid) == 2:
        nx,ny= args.grid
        nz = nx
    else:
        nx,ny,nz = args.grid[0:3]
    
    a,b,c = np.linalg.norm(lgps.get_cell()[0,:]),\
            np.linalg.norm(lgps.get_cell()[1,:]),\
            np.linalg.norm(lgps.get_cell()[2,:])
    true_spacing = (a/nx, b/ny, c/nz)
    print('True spacing: ', true_spacing)

else:
    # Small section to compute nx, ny, nz i.e. the # of points in each direction in case the --space option is used
    # This however is not needed - sampler.generate_grid_vectors has this functionality already, so we can just pass arguments to this
    # Although, that option does parheps something odd - it might ue interpolation to improve the spacing... 
    if args.space:
        if len(args.space) == 1:
            spacing_x, spacing_y, spacing_z = args.space * 3
        elif len(args.space) == 2:
            spacing_x, spacing_y = args.space
            spacing_z = spacing_x
        else:
            spacing_x, spacing_y, spacing_z =  args.space[0:3]

        a,b,c = np.linalg.norm(lgps.get_cell()[0,:]), np.linalg.norm(lgps.get_cell()[1,:]), np.linalg.norm(lgps.get_cell()[2,:])
        nx,ny,nz = ma.ceil(a/spacing_x), ma.ceil(b/spacing_y), ma.ceil(c/spacing_z)
        
        if use_sym:
            nx,ny,nz = symmetry.find_spacegroup_compatible_gridshape((nx,ny,nz), spacegroup, search_denser=True)
            print("Determined grid spacing compatible with spacegroup.")

        true_spacing = (a/nx, b/ny, c/nz)
        print('Desired spacing: ', (spacing_x, spacing_y, spacing_z),' True spacing: ', true_spacing)

print("vdW cutoff factor: ", args.vdw)

### End of parsing of arguments

# Now handle naming of the calculation files: Chosen format: inputfile_method_gridsize_jobname
# Could addversion if spacing is given, to indicate spacing?
#calc_name = '_'.join((input_basename, method, 'x'.join((str(nx), str(ny), str(nz))), jobname)) # For grid size
if args.name:
    calc_name = '_'.join((input_basename, method, 'x'.join(tuple(map(str, np.around(true_spacing,4)))), jobname))
else:
    calc_name = '_'.join((input_basename, method, 'x'.join(tuple(map(str, np.around(true_spacing,4))))))



### Set the directories and relted/relevant environment vars
#calc_dir = './UCS_CALCS_s0.4/' + calc_name
#calc_dir = './ucs_out_mel_data/' + calc_name
work_dir = '.'

#CP2K.command = "env OMP_NUM_THREADS=32 srun cp2k_shell.psmp"
CP2K.command = "env OMP_NUM_THREADS=4 cp2k_shell"   ## Should be the same / H
#if not os.path.isdir(Path(calc_dir)):
#    os.mkdir(calc_dir)                ## Might be better later to export these in terminal bash script /H
#os.environ['UCS_CALCULATION_DIR'] = calc_dir  #'./LTA_lammps_ff_grid_test' #./LTA4A_reduced_grid_gen'
os.environ['UCS_WORK_DIR'] = work_dir
os.environ['OMP_NUM_THREADS'] = '4'


# CP2K stuff, for DFT 
# Path to previously relaxed single point calculation with an atom placed at
# the first sampled position in this case (0, 0, 0)
#restart_file_name = "./single_point_calc/LTA4A_reduced-RESTART.wfn"  ## This is the old line /H 

# TODO: Implement method to read lammps parameters from file - decorator maybe
##########  Define methods for computing energy here, like so /H  ########## 

## CP2K DFT (PBE functional)
@subdir_calc
def dft_pbe(atoms: ase.Atoms):
    global restart_file_name

    # Use the same parameter here as in `single_point_calculation/` ## Need to check how to structure this, the input /H
    ## Remember to alter the charge to the correct one
    inp = '''
&FORCE_EVAL
  &DFT
    CHARGE -11
    LSD
    RESTART_FILE_NAME {}
    &SCF
       EPS_SCF 5.00000000E-005
       SCF_GUESS RESTART
       &OT T
           PRECONDITIONER FULL_SINGLE_INVERSE
           MINIMIZER CG
           LINESEARCH 3PNT
       &END OT
    &END SCF
  &END DFT
&END FORCE_EVAL
    '''.format(restart_file_name)

    # Change the restart file to the previous calculation ## As a note, this makes parallellization of calculations impossible (unless using initial calculation) /H
    p = Path('..').glob('**/*')
    dirs = [x for x in p if x.is_dir()]
    restart_file_name = str(
        Path(max(dirs, key=os.path.getmtime), 'LTA4A_reduced-RESTART.wfn'))

    calc = CP2K(inp=inp, pseudo_potential='GTH-PBE',
                max_scf=600, xc='PBE', print_level='LOW')
    atoms.set_calculator(calc)
    return atoms.get_potential_energy()


## CP2K Lennard-Jones
@subdir_calc
def cp2k_lj(atoms: ase.Atoms):
    # TODO: Finish this energy evaluation method (the input is half done at this point...)
    inp = '''
&FORCE_EVAL
  &MM
    &FORCEFIELD
      &NONBONDED
        &LENNARD-JONES
          ATOMS Na Si
          EPSILON
          RCUT
          SIGMA
        &END LENNARD-JONES
        &LENNARD-JONES
          ATOMS Na Si
          EPSILON
          RCUT
          SIGMA
        &END LENNARD-JONES
        &LENNARD-JONES
          ATOMS Na Si
          EPSILON
          RCUT
          SIGMA
        &END LENNARD-JONES

       &OT T
           PRECONDITIONER FULL_SINGLE_INVERSE
           MINIMIZER CG
           LINESEARCH 3PNT
       &END OT
    &END SCF
  &END DFT
&END FORCE_EVAL
    '''.format(restart_file_name)

    # Change the restart file to the previous calculation ## As a note, this makes parallellization of calculations impossible (unless using initial calculation) /H
    p = Path('..').glob('**/*')
    dirs = [x for x in p if x.is_dir()]
    restart_file_name = str(
        Path(max(dirs, key=os.path.getmtime), 'LTA4A_reduced-RESTART.wfn'))

    calc = CP2K(inp=inp, pseudo_potential='GTH-PBE',
                max_scf=600, xc='PBE', print_level='LOW')
    atoms.set_calculator(calc)
    return atoms.get_potential_energy()




# Lammps / ff stuff
## LAMMPS very simple lj force field /H
@subdir_calc
def lammps_lj(atoms: ase.Atoms):

    #lammps_header = ['units metal']
    #lammps_header = ['units real']


    log_file = "job.log"
    
    lammps_lj_sigma = 3.0/(2.0**(1.0/6.0))
    lammps_lj_cutoff = 3.0 * lammps_lj_sigma
    # Debug 
    #print('lammps_lj_sigma: ', lammps_lj_sigma)


    lammps_cmds = ['pair_style lj/cut '+str(lammps_lj_cutoff), 'pair_coeff * *  1.0 '+str(lammps_lj_sigma)+' '+str(lammps_lj_cutoff),
                   'pair_modify shift yes'] ## eps = 1.0 eV sigma = 3.0/(2.0**(1.0/6.0)) Å (for metal units, i.e. default) .. /H  
    
    calc = LAMMPSlib(lmpcmds=lammps_cmds, log_file=log_file)
    atoms.set_calculator(calc)
    return atoms.get_potential_energy()


## LAMMPS simple force-field lj + electrostatics with Ewald summation /H
@subdir_calc
def lammps_lj_coul(atoms: ase.Atoms):
    # TODO: Add electrostatics and Ewald summation to energy calculation, see to it that lj is the same as the simple lj
    # add definition specification of atomic charges DONE
    lammps_header = ['units metal',
                     'atom_style full']

    #lammps_header = ['units real',
    #                 'atom_style full']

    lammps_atom_types = {'Si':1, 'Al':2, 'O':3, 'B':4, 'Na':5, 'Li':6} # Note B is just a dummy for now - it is originally Si-O-Al while type 3 is Si-O-Si / H

    # Simple charges for testing etc.
    """lammps_amendments = ['set type 1 charge 1.9583', 
                         'set type 2 charge 1.0000', 
                         'set type 3 charge -1.0000', 
                         'set type 4 charge 0.0000', 
                         'set type 5 charge 1.0000'] """

    # Redox numbers as charges (if yhou want a small correction to get charge neutrality, put Si to 3.9583)
    lammps_amendments = ['set type 1 charge 4.0000', 
                         'set type 2 charge 3.0000', 
                         'set type 3 charge -2.0000', 
                         'set type 4 charge 0.0000', 
                         'set type 5 charge 1.0000'] 


    # charges from the boulfelfel FF
    """#lammps_amendments = ['set type 1 charge 1.8708', 
                         'set type 2 charge 1.7906', 
                         'set type 3 charge -0.9364', 
                         'set type 4 charge -1.1427', 
                         'set type 5 charge 0.9094'] """
    

    log_file = "job.log"
    
    lammps_lj_sigma = 3.0/(2.0**(1.0/6.0))
    lammps_lj_cutoff = 3.0 * lammps_lj_sigma
    lammps_lj_eps = 1.0

    # Debug 
    print('lammps_lj_sigma: ', lammps_lj_sigma)

    lammps_cmds = ['pair_style lj/cut/coul/long '+str(lammps_lj_cutoff)+' 11.0',
                   'kspace_style ewald 1.0e-6', 
                   'pair_coeff * * 0.0 1.0 1.0', 
                   'pair_coeff 3 5 '+str(lammps_lj_eps)+' '+str(lammps_lj_sigma)+' '+str(lammps_lj_cutoff),
                   'pair_modify shift yes'] ## eps = 1.0 eV sigma = 3.0/(2.0**(1.0/6.0)) Å (for metal units, i.e. default) .. /H  
    
    calc = LAMMPSlib(lmpcmds=lammps_cmds,
                     atom_types=lammps_atom_types,
                     lammps_header=lammps_header,
                     amendments=lammps_amendments,
                     log_file=log_file) 
                     
    atoms.set_calculator(calc)
    return atoms.get_potential_energy()

## LAMMPS force-field: lj + electrostatics with Ewald summation, from Garcia-Sanches et al 2009 
@subdir_calc
def lammps_ff_garcia_sanches(atoms: ase.Atoms):
    # TODO: Enter the correct charges DONE, and LJ-parameters from the article
    lammps_header = ['units metal',
                     'atom_style full']

    #lammps_header = ['units real',
    #                 'atom_style full']

    lammps_atom_types = {'Si':1, 'Al':2, 'O':3, 'B':4, 'Na':5} # Note O 4 is just a dummy for now - it is originally Si-O-Al while type 3 is Si-O-Si / H

    ## Charges in the force-field
    # charges for only Si, i.e. no Al, all O are Si-O-Si:
    #lammps_amendments = ['set type 1 charge 0.78598', 
    #                     'set type 2 charge 0.48598', 
    #                     'set type 3 charge -0.39299', 
    #                     'set type 4 charge 0.00000', 
    #                     'set type 5 charge 0.38340'] # 1 Si q=0.785 98, 2 Al q=0.485 98, 3 O-Si q=-0.392 99, 4 B (O-Al) q=-0.413 84, 5 Na q=0.383 40 

    # Charges for Si:Al = 1:1, i.e. all O atoms are Si-O-Al:
    lammps_amendments = ['set type 1 charge 0.78598', 
                         'set type 2 charge 0.48598', 
                         'set type 3 charge -0.41384', 
                         'set type 4 charge 0.00000', 
                         'set type 5 charge 0.38340'] # 1 Si q=0.785 98, 2 Al q=0.485 98, 3 O-Si q=-0.392 99, 4 B (O-Al) q=-0.413 84, 5 Na q=0.383 40 
 

    log_file = "job.log"
    
    # Parameters from the article Garcia-Sanches (2009): eps/kb=23.000 [K], sigma = 3.400 [Å]
    lammps_lj_sigma = 3.400 # [eV]
    lammps_lj_cutoff = 12.0 # [Å] (mentions cut-off of 12 Å in article; a bit unclear if it refers to all LJ-potentials however)
    lammps_lj_eps = 23.000 * ase.units.kB  # [K] * [eV/K] = [eV]

    # Debug 
    #print('lammps_lj_sigma: ', lammps_lj_sigma)

    lammps_cmds = ['pair_style lj/cut/coul/long '+str(lammps_lj_cutoff)+' 11.0',
                   'kspace_style ewald 1.0e-6', 
                   'pair_coeff * * 0.0 1.0 1.0', 
                   'pair_coeff 3 5 '+str(lammps_lj_eps)+' '+str(lammps_lj_sigma)+' '+str(lammps_lj_cutoff),
                   'pair_modify shift yes'] 
    
    calc = LAMMPSlib(lmpcmds=lammps_cmds,
                     atom_types=lammps_atom_types,
                     lammps_header=lammps_header,
                     amendments=lammps_amendments,
                     log_file=log_file) 
                     
    atoms.set_calculator(calc)
    return atoms.get_potential_energy()





## LAMMPS force field by Boulfelfel et al (2021) J. Phys. Chem. C /H
@subdir_calc
def lammps_ff_boulfelfel(atoms: ase.Atoms):

    log_file = "job.log"
    # Set the atom type numbers for lammps
    lammps_atom_types = {'Si':1, 'Al':2, 'O':3, 'B':4, 'Na':5, 'Li':6} # Note B is just a dummy for now - it is originally Si-O-Al while type 3 is Si-O-Si / H

    # lammps header to control units etc. TODO: CHECK THAT UNITS ARE ACTUALLY CORRECT (I.e. compatible between lammps and ase --> THEY ARE. ASE converts the units based on the unit specification in the lammps input and produces its own output, i.e. in eV and Å units)
    #lammps_header = ['atom_style full'] # no unit specification (i.e. default units) DOES NOT SEEM TO WORK - RESULTS IN ERROR IN THE ASE MODULE

    lammps_header = ['units metal', 
                     'atom_style full'] # For metal units /H 

    #lammps_header = ['units real', 
    #                 'atom_style full'] # Real units (based on compatibility wih ASE?) /H


    # Add setting of atomic charges in amendments 
    lammps_amendments = ['set type 1 charge 1.8708', 
                         'set type 2 charge 1.7906', 
                         'set type 3 charge -0.9364', 
                         'set type 4 charge -1.1427', 
                         'set type 5 charge 0.9094'] 

    # lammps commmands defining the force field in [SOURCE Boulfelfel et al 2021] for the framework and cations (real units)
    """lammps_cmds = ['pair_style hybrid/overlay buck 11.0 coul/long 11.0', 
                   'pair_modify tail yes', 
                   'kspace_style ewald 1.0e-6', 
                   'pair_coeff * * coul/long', 
                   'pair_coeff 3 5 buck 110905.00 0.25094 1821.66']  # For "units real". Define the lattice-cation interaction TODO: Check units, recompute these values in correct units if wrong /H """

    lammps_cmds = ['pair_style hybrid/overlay buck 11.0 coul/long 11.0', 
                   'pair_modify tail yes', 
                   'kspace_style ewald 1.0e-6', 
                   'pair_coeff * * coul/long', 
                   'pair_coeff 3 5 buck 4809.30 0.25094 78.9947']  # For "units metal", default when using ase. Define the lattice-cation interaction TODO: Check units, recompute these values in correct     units if wrong /H """

    # Note: Real units are: Energy = Kcal/mole, Distance = Ångström, Mass = gram/mole, charge = (electron) elementary charge e
    # Default units (when using the ase calculator) are instead of type "metal": Energy = eV, Distance = Å, Mass = gram/mole, charge = e (elementary charge)
    # Kcal/mole -> eV: 1 kcal/mole = (1/Na)*1000*4.184 J = (1/Na)*1000*4.184*(1/e) 1 eV = e* 1 V = |e| J => X J = X/|e| eV
    # Thus:   110905.00 kcal/mole     = 110905.00 * (1/(6.0221408 * 10**23))*1000*4.184*(1/(1.60217663*10**(-19))) eV = 4809.295960999744 eV ≃ 4809.30 eV
    #         0.25094 Å               = 0.25094 Å
    #         1821.66 kcal/mole * Å^6 = 1821.66 * (1/(6.0221408 * 10**23))*1000*4.184*(1/(1.60217663*10**(-19))) eV * Å^6 = 78.99465380564261 eV Å^6 ≃ 78.9947 eV Å^6




    calc = LAMMPSlib(lmpcmds=lammps_cmds, 
                     atom_types=lammps_atom_types, 
                     lammps_header=lammps_header, 
                     amendments=lammps_amendments, 
                     log_file=log_file)

    atoms.set_calculator(calc)
    return atoms.get_potential_energy()


## LAMMPS Boulfelfel et al forcefield except only the Buckingham part (i.e. without electrostatics) /H
@subdir_calc
def lammps_ff_boulfelfel_buck(atoms: ase.Atoms):
 
    log_file = "job.log"
    # Set the atom type numbers for lammps
    lammps_atom_types = {'Si':1, 'Al':2, 'O':3, 'B':4, 'Na':5, 'Li':6} # Note B is just a dummy for now - it is originally Si-O-Al while type 3 is Si-O-Si / H

    # lammps header to control units etc. TODO: CHECK THAT UNITS ARE ACTUALLY CORRECT (I.e. compatible between lammps and ase)
    #lammps_header = ['units metal', 

    lammps_header = ['units real',
                     'atom_style full'] # For real units, compatible with certain force-field parameters, see below /H 

    #lammps_header = ['atom_style full'] # Default units (based on compatibility wih ASE?) /H

    # lammps commmands defining the force field in [SOURCE Boulfelfel et al 2021] for the framework and cations
    lammps_cmds = ['pair_style buck 11.0', 
                   'pair_modify tail yes', 
                   'pair_coeff * * 0.0 1.0 0.0', 
                   'pair_coeff 3 5 110905.00 0.25094 1821.66']  # For "units real". Define the lattice-cation interaction TODO: Check units, recompute these values in correct units if wrong /H """

    """lammps_cmds = ['pair_style buck 11.0', 
                   'pair_modify tail yes',
                   'pair_coeff * * 0.0 1.0 0.0', 
                   'pair_coeff 3 5 4809.30 0.25094 78.9947']  # For "units metal", default when using ase. Define the lattice-cation interaction TODO: Check units, recompute these values in correct units if wrong /H """

    calc = LAMMPSlib(lmpcmds=lammps_cmds, 
                     atom_types=lammps_atom_types, 
                     lammps_header=lammps_header,
                     log_file=log_file)   
    #lammps_header=lammps_header, 

    atoms.set_calculator(calc)
    return atoms.get_potential_energy()

### Do I need subdir_calc here? / H
@subdir_calc
def ase_lj(atoms: ase.Atoms):
    sigma = 3.0/(2.0**(1.0/6.0))
    calc = LJ(sigma=sigma, epsilon=1.0, rc=3*sigma, smooth=False)
    atoms.set_calculator(calc)
    return atoms.get_potential_energy()

###

# Manually defined ff for one of the structures, but with another cut-off just for testing
@subdir_calc
def structure_54189_ff_manually_coded_ALTERED_CUTOFF(atoms: ase.Atoms):
    
    lmpcmds = ["pair_style lj/cut/coul/long 1.0",
               "pair_modify     tail yes mix arithmetic",
               "special_bonds   lj/coul 0.0 0.0 1.0",
               "dielectric      1.0",
               "kspace_style ewald 1.0e-5",
               "pair_coeff 1 1 0.025000 2.183593",
               "pair_coeff 2 2 0.124000 2.461553",
               "pair_coeff 3 3 0.379000 3.813047"] # "box tilt        large", went here before 
    
    atom_types = {"Li":1, "Zn":2, "Ge":3}
    atom_type_masses = {"Li":6.9410, "Zn":65.380, "Ge":72.640}
    log_file = "54189_job.log"

    # TODO: Note that "boundary p p p" and "box tilt large" are automatically added by ase in the inputfile. Can likely remove them here.
    
    lammps_header = ["units real",
                     "atom_style full",
                     "boundary p p p",
                     "box tilt large"]  # Added this line here instead

    amendments = ["set type 1 charge 0.72",
                  "set type 2 charge 1.5",
                  "set type 3 charge -2.95"]
    #post_changebox_cmds = []

    calc = LAMMPSlib(lmpcmds=lmpcmds,
                     atom_types=atom_types,
                     atom_type_masses=atom_type_masses,
                     log_file=log_file,
                     lammps_header=lammps_header,
                     amendments=amendments)
    
    atoms.set_calculator(calc)
    return atoms.get_potential_energy()


##### Lammps method read from input > for now in a separate script
#    lmp_file = 
#    lmp_dict = parser_lammps_mel_inp(lmp_file)

#    lmp_inp_calc = lammps_method_from_data(**lmp_dict)

#####




#####################################################################
# Create supercell from unitcell depending on cutoff !!!
####################################################################
# This codeblock is taken from custom_lammps_grid_gen.py

# Metainfo: lgps = unitcell atoms object whose sampling atoms have been removed if desired
#           sampler = UnitcellSampler object for the unitcell

# First construct supercell
#cutoff = 12.5 # Force cutoff in Å
cutoff = 0.0 # Force cutoff in Å
print("Force cutoff used to determine supercell size: ", cutoff)
num_cells = compute_req_supercells(lgps, cutoff)

print("num_cells in supercell: ", num_cells)

supercell_from_unitcell_wo_ions = ase.build.make_supercell(
            lgps,
            np.diag(num_cells), wrap=True
            )

print("supercell from uc params: ", supercell_from_unitcell_wo_ions.get_cell_lengths_and_angles())
print("supercell from uc cell: ", supercell_from_unitcell_wo_ions.get_cell())

print("\n")
print("Spacegroup, unitcell: ", ase.spacegroup.get_spacegroup(lgps, 1.0e-6))
print("Spacegroup, supercell: ", ase.spacegroup.get_spacegroup(supercell_from_unitcell_wo_ions, 1.0e-6))

# "Middle" unit cell (we don't need this though, since pbc)
#uc_indices = [int((i//2)-1) if i % 2 == 0 else int((i-1)//2) for i in num_cells]
uc_indices = (0, 0, 0)


# Now generate grid for unitcell:
unitcell_ucs = sample.UnitCellSampler(lgps) # DONE

if (args.midvox and use_sym):
    print("WARNING: Both symmetry and midvoxel sampling enabled. Their compatibility has NOT been properly tested!!!")

unitcell_grid, unitcell_included = unitcell_ucs.generate_grid_vectors((nx, ny, nz), vdw_scale=args.vdw, midvox=args.midvox) # DONE

unitcell_grid = unitcell_grid.reshape(-1,3) # DONE
unitcell_included = unitcell_included.reshape(-1) # DONE

# Convert to fractional coordinates, and convert to supercell grid

print("Shape check, grid: ", np.array(unitcell_grid).shape, np.array(unitcell_grid).T.shape) # DONE
unitcell_frac_grid = np.linalg.solve(np.array(lgps.get_cell()).T, np.array(unitcell_grid).T).T # DONE

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
    sampler.gemmi_unitcell = lgps

print("Spacegroup input: ", args.sg)
print("Symmetry: ", use_sym)


if method == 'pbe':
    # For DFT:
    energies = sampler.calculate_energies(grid_points=supercell_cart_grid,
        method=dft_pbe, atom=atom, exploit_symmetry=use_sym)

elif method == 'lammps_lj':
    # For ForceField TEST, simple lj:
    energies = sampler.calculate_energies(grid_points=supercell_cart_grid,
        method=lammps_lj, atom=atom, exploit_symmetry=use_sym)

elif method == 'lammps_lj_coul':
    # For ForceField test, simple lj + electrostatics with Ewald:
    energies = sampler.calculate_energies(grid_points=supercell_cart_grid,
        method=lammps_lj_coul, atom=atom, exploit_symmetry=use_sym)

elif method == 'ff_boulfelfel':
    # For ForceField (Boulfelfel et al, 2021):
    energies = sampler.calculate_energies(grid_points=supercell_cart_grid,
        method=lammps_ff_boulfelfel, atom=atom, exploit_symmetry=use_sym)

elif method == 'ff_boulfelfel_buck':
    # For ForceField (Boulfelfel et al, 2021), but only the Buckingham part:
    energies = sampler.calculate_energies(grid_points=supercell_cart_grid,
        method=lammps_ff_boulfelfel_buck, atom=atom, exploit_symmetry=use_sym)

elif method == 'ff_garcia_sanches':
    # For ForceField (Garcia-Sanches et al, 2009):
    energies = sampler.calculate_energies(grid_points=supercell_cart_grid,
        method=lammps_ff_garcia_sanches, atom=atom, exploit_symmetry=use_sym)

elif method == 'ase_lj':
    energies = sampler.calculate_energies(grid_points=supercell_cart_grid,
        method=ase_lj, atom=atom, exploit_symmetry=use_sym)

elif method == '54189':
    energies = sampler.calculate_energies(grid_points=supercell_cart_grid,
            method=structure_54189_ff_manually_coded_ALTERED_CUTOFF, atom=atom, exploit_symmetry=use_sym)

elif method == '73679':
    energies = sampler.calculate_energies(grid_points=supercell_cart_grid,
            method=struct_73679_ff, atom=atom, exploit_symmetry=use_sym)

elif method in automethods.keys():
    energies = sampler.calculate_energies(grid_points=supercell_cart_grid,
            method=automethods[method], atom=atom, exploit_symmetry=use_sym)

elif method in dft_methods.keys():
    energies = sampler.calculate_energies(grid_points=supercell_cart_grid,
            method=dft_methods[method], atom=atom, exploit_symmetry=use_sym)

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
    write_cube(fp,  lgps, data=energies_full.reshape(
        unitcell_ucs.included_grid_vectors.shape))

if xsf_output:
    with open(xsf_filename, 'w') as fp:
        write_xsf(fp,  lgps, data=energies_full.reshape(
            unitcell_ucs.included_grid_vectors.shape))

