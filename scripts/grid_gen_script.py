#!/usr/bin/env python

# This is based on the unitcellsampling/examples/dft_grid_gen/grid_gen .. script named grid_gen.py by Benjamin Bolbrinker. Modified by Hannes Gustafsson.

from pathlib import Path
from unitcellsampling import sample
import ase.io
import ase

# energy calcs
from ase.calculators.lj import LennardJones as LJ
from ase.calculators.cp2k import CP2K
from ase.calculators.lammpsrun import LAMMPS
from ase.calculators.lammpslib import LAMMPSlib

from unitcellsampling.lammps_calc_from_inp import parser_lammps_mel_inp
from unitcellsampling.lammps_calc_from_inp import lammps_method_from_data

from unitcellsampling.special_methods import struct_73679_ff
# predefined by Ben
#import energy_calculators

# File handling
from ase.io.cube import write_cube
from ase.io.xsf import write_xsf # Added to test if its better /H 

from unitcellsampling.decorators import subdir_calc
import os
import argparse
import numpy as np
import math as ma

from unitcellsampling.preparatory_fcns import unitcell_to_supercell_frac_coords
from unitcellsampling.preparatory_fcns import remove_nonframework_cations_fancy
# Defaults

# Last update: 

### TODO: Create input argument handling, to make usage much smoother. Some particular changes:
     # -  TODO: make the  program either recognize the filename from the input file automatically, and name the "project" accordingly
     #          or make the desired project name an input so that the variables all are set in conjuction /H
     # -  TODO: Make the subdirectory for the gridsampling also be project name-derived, and check so that it cannot overide anything by default /H


# TODO: Actually: Move the setting of a project name and setting of environment variables to an exterior runscript - MAYBE, but arguably could be set here too.
# TODO: Move the definition the the method arguments to this list here, for easier editing
method_list = ['pbe', 'lammps_lj', 'lammps_lj_coul', 'ff_boulfelfel', 'ff_boulfelfel_buck', 'ff_garcia_sanches', 'ase_lj', '54189', '73679']
### Definition and parsing of arguments
parser = argparse.ArgumentParser(description='Energy sampling of a (periodic) solid system with an added cation on the grid. Methods: (PBE/PBE-GTH/DZVP-MOLOPT-SR-GTH), Forcefields (Scholl et al)') 
parser.add_argument('file', metavar='filename', type=str, action='store', help="Name of the cif-file containing the structure, without sample atom/ion.")
#parser.add_argument('-m', '--method', type=str, action='store', default='ff_boulfelfel', choices=['pbe', 'ff', 'ff_boulfelfel', 'ff_boulfelfel_buck'], help="Method to calculate the energy during grid sampling.") # should method be optional or not? No, probably not optional? Or maybe optinal with default?
parser.add_argument('method', type=str, action='store', choices=method_list, help="Method to calculate the energy during grid sampling.")
parser.add_argument('-n','--name', metavar='jobname', type=str, action='store', default=None, help="Desired name for the calculation (applied to generated output files, directories, etc.).")
parser.add_argument('-w', '--wfn', type=str, action='store', default=None, help="Specify the initial wfn-file for a DFT calculation.")
parser.add_argument('-a', '--atom', type=str, action='store', default='Na', help="Specify the atom ((cat)ion) used for sampling.")
parser.add_argument('-g', '--grid', type=int, action='store', default=[10], nargs='+', help="Specify the number of grid points in each dimension (or cubic grid) (mutually exclusive with \"--space\").")
parser.add_argument('-s', '--space', type=float, action='store', default=None, nargs='+', help="Specify the spacing between the grid points in each dimension (mutually exclusive with \"--grid\").")
parser.add_argument('--vdw', type=float, action='store', default=1.0, help="Specify the fraction of the van der Waals radius that should be excluded from the sampled volume around each atom in thei host structure.")
parser.add_argument('--nosym', action='store_false', help="Turns off usage of spacegroup symmetry. Default is to apply symmetry to save the number of required calculations.")
parser.add_argument('--ra', action='store_true', help="Specify whether to remove all atoms of the type that is used for sampling from the structure, before doing the sampling.")
parser.add_argument('--sg', type=int, default=None, action='store', help="Manually specify the spacegroup to use for symmetry. Default is None, in which case spacegroup will be automatically determined from the stucture.")

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
#atoms = ase.io.read(input_file) # Like so? 
if args.ra:
    print("Removing sampling atoms", str(atom),"from structure.")
    atoms_temp = remove_nonframework_cations_fancy(lgps, ase.Atom(atom))
    lgps = atoms_temp


## Set the number of points in the grid in each dimension (or equivalently, the mesh size)
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
        true_spacing = (a/nx, b/ny, c/nz)
        print('Desired spacing: ', (spacing_x, spacing_y, spacing_z),' True spacing: ', true_spacing)


### End of parsing of arguments

# Now handle naming of the calculation files: Chosen format: inputfile_method_gridsize_jobname
# Could addversion if spacing is given, to indicate spacing?
#calc_name = '_'.join((input_basename, method, 'x'.join((str(nx), str(ny), str(nz))), jobname)) # For grid size
if args.name:
    calc_name = '_'.join((input_basename, method, 'x'.join(tuple(map(str, np.around(true_spacing,4)))), jobname))
else:
    calc_name = '_'.join((input_basename, method, 'x'.join(tuple(map(str, np.around(true_spacing,4))))))



### Set the directories and relted/relevant environment vars
#calc_dir = './UCS_CALC/' + calc_name
calc_dir = './ucs_out_mel_data/' + calc_name
work_dir = '.'

#CP2K.command = "env OMP_NUM_THREADS=32 srun cp2k_shell.psmp"
CP2K.command = "env OMP_NUM_THREADS=4 cp2k_shell"   ## Should be the same / H
if not os.path.isdir(Path(calc_dir)):
    os.mkdir(calc_dir)                ## Might be better later to export these in terminal bash script /H
os.environ['UCS_CALCULATION_DIR'] = calc_dir  #'./LTA_lammps_ff_grid_test' #./LTA4A_reduced_grid_gen'
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

#lgps = ase.io.read(Path(indir, infile))

sampler = sample.UnitCellSampler(lgps)
sampler.generate_grid_vectors(n_frac=(nx, ny, nz), vdw_scale=args.vdw)
sampler.spacegroup = args.sg # Set spacegroup for the sampler

print("Spacegroup input: ", args.sg)
print("Symmetry: ", use_sym)
if method == 'pbe':
    # For DFT:
    energies = sampler.calculate_energies(
        method=dft_pbe, atom=atom, exploit_symmetry=use_sym)

elif method == 'lammps_lj':
    # For ForceField TEST, simple lj:
    energies = sampler.calculate_energies(
        method=lammps_lj, atom=atom, exploit_symmetry=use_sym)

elif method == 'lammps_lj_coul':
    # For ForceField test, simple lj + electrostatics with Ewald:
    energies = sampler.calculate_energies(
        method=lammps_lj_coul, atom=atom, exploit_symmetry=use_sym)

elif method == 'ff_boulfelfel':
    # For ForceField (Boulfelfel et al, 2021):
    energies = sampler.calculate_energies(
        method=lammps_ff_boulfelfel, atom=atom, exploit_symmetry=use_sym)

elif method == 'ff_boulfelfel_buck':
    # For ForceField (Boulfelfel et al, 2021), but only the Buckingham part:
    energies = sampler.calculate_energies(
        method=lammps_ff_boulfelfel_buck, atom=atom, exploit_symmetry=use_sym)

elif method == 'ff_garcia_sanches':
    # For ForceField (Garcia-Sanches et al, 2009):
    energies = sampler.calculate_energies(
        method=lammps_ff_garcia_sanches, atom=atom, exploit_symmetry=use_sym)

elif method == 'ase_lj':
    energies = sampler.calculate_energies(
        method=ase_lj, atom=atom, exploit_symmetry=use_sym)

elif method == '54189':
    energies = sampler.calculate_energies(
            method=structure_54189_ff_manually_coded_ALTERED_CUTOFF, atom=atom, exploit_symmetry=use_sym)

elif method == '73679':
    energies = sampler.calculate_energies(
            method=struct_73679_ff, atom=atom, exploit_symmetry=use_sym)

else:
    print("No default method defined yet.")
    raise Exception('No method.')

print("Included grid vectors shape: ", sampler.included_grid_vectors.shape)
cube_filename = ".".join((calc_name, 'cube'))
xsf_filename = ".".join((calc_name, 'xsf'))
print()
#print(*np.array2string(energies.reshape((len(energies),1)), formatter={"float_kind":lambda x: "%.7f" % x }))
for row in energies:
    print(row)


with open(cube_filename, 'w') as fp:
    write_cube(fp,  lgps, data=energies.reshape(
        sampler.included_grid_vectors.shape))

if xsf_output:
    with open(xsf_filename, 'w') as fp:
        write_xsf(fp,  lgps, data=energies.reshape(
            sampler.included_grid_vectors.shape))

