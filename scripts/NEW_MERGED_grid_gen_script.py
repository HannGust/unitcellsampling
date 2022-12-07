#!/usr/bin/env python

# This is based on the unitcellsampling/examples/dft_grid_gen/grid_gen .. script named grid_gen.py by Benjamin Bolbrinker. Modified by Hannes Gustafsson.
#!/usr/bin/env python

# This is based on the unitcellsampling/examples/dft_grid_gen/grid_gen .. script named grid_gen.py by Benjamin Bolbrinker. Modified by Hannes Gustafsson.


from genericpath import isdir
from pathlib import Path
from unitcellsampling import sample
import ase.io
import ase
import ase.build

# energy calculators
from ase.calculators.lammpsrun import LAMMPS
from ase.calculators.lammpslib import LAMMPSlib

from lammps_calc_from_inp import parser_lammps_mel_inp
from lammps_calc_from_inp import lammps_method_from_data

from special_methods import struct_54209_ff, \
                            struct_54297_ff, \
                            struct_54449_ff, \
                            struct_54683_ff, \
                            struct_54837_ff, \
                            struct_54865_ff, \
                            struct_54879_ff, \
                            struct_54884_ff, \
                            struct_55184_ff, \
                            struct_55319_ff


from special_methods import struct_56568_ff, \
                            struct_56627_ff, \
                            struct_57382_ff


# File handling
from ase.io.cube import write_cube
from ase.io.xsf import write_xsf # Added to test if its better /H 
from ase.io.lammpsdata import read_lammps_data, write_lammps_data

from read_atom_types import read_lammpsdata_atom_info

from decorators import subdir_calc
import os
import argparse
import numpy as np
import math as ma

from preparatory_fcns import unitcell_to_supercell_frac_coords, \
                             remove_nonframework_cations_fancy, \
                             compute_req_supercells



from pathlib import Path
from unitcellsampling import sample
import ase.io
import ase

# energy calcs
from ase.calculators.lj import LennardJones as LJ
from ase.calculators.cp2k import CP2K
from ase.calculators.lammpsrun import LAMMPS
from ase.calculators.lammpslib import LAMMPSlib

from lammps_calc_from_inp import parser_lammps_mel_inp
from lammps_calc_from_inp import lammps_method_from_data

from special_methods import struct_73679_ff
# predefined by Ben
#import energy_calculators

# File handling
from ase.io.cube import write_cube
from ase.io.xsf import write_xsf # Added to test if its better /H 

from decorators import subdir_calc
import os
import argparse
import numpy as np
import math as ma

from preparatory_fcns import unitcell_to_supercell_frac_coords
from preparatory_fcns import remove_nonframework_cations_fancy
# Defaults

# Last update: 

### TODO: Create input argument handling, to make usage much smoother. Some particular changes:
     # -  TODO: make the  program either recognize the filename from the input file automatically, and name the "project" accordingly
     #          or make the desired project name an input so that the variables all are set in conjuction /H
     # -  TODO: Make the subdirectory for the gridsampling also be project name-derived, and check so that it cannot overide anything by default /H


# TODO: Actually: Move the setting of a project name and setting of environment variables to an exterior runscript - MAYBE, but arguably could be set here too.
# TODO: Move the definition the the method arguments to this list here, for easier editing
method_list = ['pbe', 'lammps_lj', 'lammps_lj_coul', 'ff_boulfelfel', 
               'ff_boulfelfel_buck', 'ff_garcia_sanches', 'ase_lj']
method_list.extend(['54189',
                    '54209',
                    '54297',
                    '54449',
                    '54683',
                    '54837',
                    '54865',
                    '54879',
                    '54884',
                    '55184',
                    '55319',
                    '56568',
                    '56627',
                    '57382',
                    '73679'])
### Definition and parsing of arguments
parser = argparse.ArgumentParser(description='Energy sampling of a (periodic) solid system with an added cation on the grid.') 
parser.add_argument('file', metavar='filename', type=str, action='store', help="Name of the cif-file containing the structure, without sample atom/ion.")
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

#

# Parser
parser.add_argument('unitcell', type=str, action='store', help="Path to the cif-file containing the unit cell of the system (not supercell!).")
parser.add_argument('infile', metavar='filename', type=str, action='store', help="Path to the LAMMPS inputfile on which to base calculation. (NOT IMPLEMENTED)")
parser.add_argument('-d','--datafile', type=str, action='store', default=None, help="Path to the LAMMPS-datafile containing atomic information, force-field parameters, and atomic coordinates (for supercell).")
parser.add_argument('-n','--name', metavar='jobname', type=str, action='store', default=None, help="Desired name for the calculation (applied to generated output files, directories, etc.).")
parser.add_argument('-a', type=str, action='store', help="Atom to sample with. Default: Li", default="Li")
parser.add_argument('-g', '--grid', type=int, action='store', default=[10], nargs='+', help="Specify the number of grid points in each dimension (or cubic grid) (mutually exclusive with \"--space\").")
parser.add_argument('-s', '--space', type=float, action='store', default=None, nargs='+', help="Specify the spacing between the grid points in each dimension (mutually exclusive with \"--grid\").")
parser.add_argument('-m', '--method', type=str, action='store', default=None, choices=method_list, help="Method to calculate the energy during grid sampling. OBS: Currently specific to structures with the same number.")


args = parser.parse_args()
#


use_sym = args.nosym
if use_sym:
    print("Symmetry on.")
    params_name_addon = ""
else:
    print("Symmetry off.")
    params_name_addon = "nosym_"

params_name_addon += "vdw"+str(args.vdw)

# Handle the unitcell input file and name
unitcell_file = Path(args.unitcell)
if not unitcell_file.exists():
    raise FileNotFoundError('The inputfile \'' + str(unitcell_file) + '\' was not found. Need file containing unitcell.')

input_basename = unitcell_file.name

while len(Path(input_basename).suffixes) > 0:
    input_basename = Path(input_basename).stem # input_basename is here meant to be just the name of the input file without parent directories and extensions: /foo/bar.cif -> input_basename = bar
print("Input basename: ", input_basename)





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





### Here bagins other file


# Preparatory things...
# Remove ions of the type that is used to sample


sample_ion = ase.Atom(args.a)
unitcell = ase.io.read(args.unitcell)

unitcell_wo_ions = remove_nonframework_cations_fancy(unitcell, sample_ion)

#if args.datafile:
#    ref_supercell_wo_ions = remove_nonframework_cations_fancy(ref_supercell, sample_ion)
#    print("ref supercell params: ", ref_supercell_wo_ions.get_cell_lengths_and_angles())
#    print("ref supercell cell: ", ref_supercell_wo_ions.get_cell())
    
#tmp_Li_atom = sample_ion
#tmp_Li_atom.position = list(np.array([0.522, 0.548, 0.466]) @ ref_supercell_atoms_wo_ions.get_cell())
#
#
#ref_supercell_atoms_with_only_one_Li = ref_supercell_atoms_wo_ions + tmp_Li_atom
#with open("lammps_data_ref_supercell_only_one_Li_0.522_0.548_0.466.54189", "w") as f:
#    write_lammps_data(f, ref_supercell_atoms_with_only_one_Li, units="real", atom_style="full")

#wrapped_ref_atoms = ref_supercell_atoms_wo_ions.copy()
#wrapped_ref_atoms.wrap() # -> This doesn't seem to do anything ???
#atoms = remove_nonframework_cations(supercell_atoms, sample_ion)

# Debug - filenames
print("Filenames:")
print("infile: ", args.infile)
print("datafile: ", args.datafile)
print("unitcell (cif): ", args.unitcell)

print("\n\n")

# Debug - various atoms objects
print("Atoms objects (structures):")
#if args.datafile:
#    print("supercell_atoms: ", ref_supercell)
print("sample_ion: ", sample_ion)
print("unitcell: ", unitcell)

#if args.datafile:
#    print("supercell w/o sample ion: ", ref_supercell_wo_ions)
print("unitcell w/o sample ion: ", unitcell_wo_ions)

print("\n\n")
print("Parameters: ")
print("vdw: ", args.vdw)
if args.space:
    print("Grid spacing: ", args.space)
else:
    print("Grid shape/size: ", args.grid)

#############################################################
# Now we should have clean supercell and ordinary unitcell  #
# Clean supercell should be passed to UCS together with     #
# grid which are fractional supercell-coordinates for the   #
# desired grid in one sub-unitcell within the supercell.    #
#############################################################

# we need to compute the number of unitcells in each
# spatial direction of the supercell

# If C = [[cx, 0, 0], [0, cy, 0], [0, 0, cz]], where cx, cy, cz is the number of unitcells in
# every direction, and A is the cell matrix i.e.
# its rows is the cell vectors, then
# B = C @ A => C = B @ A.inv
#############################################################

# Since we have two ways of constructing the supercell, I will
# choose the approach of reading the supercell object from
# file i.e. as input if possible.

#if args.datafile:
#    C = ref_supercell.get_cell() @ np.linalg.inv(unitcell.get_cell())
#    num_cells = np.rint(np.diagonal(np.copy(C)))
#    num_cells = tuple([int(i) for i in num_cells])
#
#    P_mat = np.diag(num_cells)
#    print(ref_supercell.get_masses())
#
#    supercell_from_unitcell_wo_ions = ase.build.make_supercell(unitcell_wo_ions, P_mat, wrap=True)


# Naive:
#supercell_from_unitcell_wo_ions = unitcell_wo_ions.repeat(c)
# Making supercell from unitcell:

cutoff = 12.5 # Force cutoff in Å
print("Force cutoff used to determine supercell size: ", cutoff)
num_cells = compute_req_supercells(unitcell_wo_ions, cutoff)
print("num_cells in supercell: ", num_cells, (num_cells,)*3)

supercell_from_unitcell_wo_ions = ase.build.make_supercell(
            unitcell_wo_ions, 
            np.diag((num_cells, num_cells, num_cells)), wrap=True
            )

print("supercell from uc params: ", supercell_from_unitcell_wo_ions.get_cell_lengths_and_angles())
print("supercell from uc cell: ", supercell_from_unitcell_wo_ions.get_cell())

print("\n")
print("Spacegroup, unitcell: ", ase.spacegroup.get_spacegroup(unitcell_wo_ions, 1.0e-6))
print("Spacegroup, supercell: ", ase.spacegroup.get_spacegroup(supercell_from_unitcell_wo_ions, 1.0e-6))

# "Middle" unit cell (we don't need this though, since pbc)
#uc_indices = [int((i//2)-1) if i % 2 == 0 else int((i-1)//2) for i in num_cells]
uc_indices = (0, 0, 0)
print(uc_indices)


# Debug: Write out supercells to cubefiles for comparison
if False:
    #with open("ref_sc.cube", 'w') as file:
    #    write_cube(file, ref_supercell)
    
    with open("sc_from_uc_"+params_name_addon+".cube", 'w') as file:
        write_cube(file, supercell_from_unitcell_wo_ions)


# Grid
## Set the number of points in the grid in each dimension (or equivalently, the mesh size)
if not args.space:
    if len(args.grid) == 1:
        nx,ny,nz = args.grid * 3
    elif len(args.grid) == 2:
        nx,ny= args.grid
        nz = nx
    else:
        nx,ny,nz = args.grid[0:3]
    
    a,b,c = np.linalg.norm(unitcell.get_cell()[0,:]),\
            np.linalg.norm(unitcell.get_cell()[1,:]),\
            np.linalg.norm(unitcell.get_cell()[2,:])
    true_spacing = (a/nx, b/ny, c/nz)
    print('True spacing: ', true_spacing)
else:
    # Small section to compute nx, ny, nz i.e. the # of points in each direction in case the --space option is used
    # This however is not needed - sampler.generate_grid_vectors has this functionality already, so we can just pass arguments to this
    # Although, that option does parheps something odd - it might ue interpolation to improve the spacing... 

    if len(args.space) == 1:
        spacing_x, spacing_y, spacing_z = args.space * 3
    elif len(args.space) == 2:
        spacing_x, spacing_y = args.space
        spacing_z = spacing_x
    else:
        spacing_x, spacing_y, spacing_z =  args.space[0:3]

    a,b,c = np.linalg.norm(unitcell.get_cell()[0,:]), np.linalg.norm(unitcell.get_cell()[1,:]), \
            np.linalg.norm(unitcell.get_cell()[2,:])

    nx,ny,nz = ma.ceil(a/spacing_x), ma.ceil(b/spacing_y), ma.ceil(c/spacing_z)
    true_spacing = (a/nx, b/ny, c/nz)
    print('Desired spacing: ', (spacing_x, spacing_y, spacing_z),' True spacing: ', true_spacing)


# TODO: Change output naming to general, input dependent, as in other script

# Naming of the calculation
if args.name and args.datafile:
    #calc_name = '_'.join((input_basename, input_datafile_name, 'x'.join(tuple(map(str, np.around(true_spacing,4)))), args.name))
    calc_name = '_'.join((input_basename, 'x'.join(tuple(map(str, np.around(true_spacing,4)))), params_name_addon, args.name))
elif args.datafile:
    #calc_name = '_'.join((input_basename, input_datafile_name, 'x'.join(tuple(map(str, np.around(true_spacing,4))))))
    calc_name = '_'.join((input_basename, 'x'.join(tuple(map(str, np.around(true_spacing,4))))), params_name_addon)
elif args.name:
    calc_name = '_'.join((input_basename, 'x'.join(tuple(map(str, np.around(true_spacing,4)))), params_name_addon, args.name))
else:
    calc_name = '_'.join((input_basename, 'x'.join(tuple(map(str, np.around(true_spacing,4))))), params_name_addon)

if args.method:
    method_str = "mtd" + args.method
    calc_name = '_'.join((calc_name, args.method))

calc_dir = "./ucs_out_mel_data/" + calc_name
if not (Path(calc_dir).exists()):
    os.mkdir(Path(calc_dir))
elif not Path(calc_dir).is_dir():
    raise Exception("Someting already has the name " + str(calc_dir) + ", but it's not a directory.")
elif Path(calc_dir).is_dir():
    print("WARNING: The calculation directory already exists. Proceeds. Be catious.")

work_dir = "."
os.environ["UCS_CALCULATION_DIR"] = calc_dir
os.environ["UCS_WORK_DIR"] = work_dir

# Now generate grid for unitcell:
unitcell_ucs = sample.UnitCellSampler(unitcell_wo_ions)
unitcell_grid, unitcell_included = unitcell_ucs.generate_grid_vectors((nx, ny, nz), vdw_scale=args.vdw)

unitcell_grid = unitcell_grid.reshape(-1,3)
unitcell_included = unitcell_included.reshape(-1)

# Convert to fractional coordinates, and convert to supercell grid
# (rather to unitcell embedded centrally in supercell)
print("Shape check, grid: ", np.array(unitcell_grid).shape, np.array(unitcell_grid).T.shape)
unitcell_frac_grid = np.linalg.solve(np.array(unitcell.get_cell()).T, np.array(unitcell_grid).T).T
print("unitcell frac grid shape: ", unitcell_frac_grid.shape)
supercell_frac_grid = unitcell_to_supercell_frac_coords(unitcell_frac_grid[unitcell_included], (num_cells,)*3, unitcell_ind=(0,0,0))

# Convert to cartesian supercell grid
supercell_cart_grid = supercell_frac_grid @ np.array(supercell_from_unitcell_wo_ions.get_cell())
print(type(supercell_cart_grid))

##### Temporary: Just to have the grids available for plotting
#np.save("supercell_cart_grid_54189.npy", supercell_cart_grid)
#np.save("unitcell_cart_grid_54189.npy", unitcell_grid)
#####

#exit() # Comment/decomment if you want to do or not to do actual calculation
# Finally make sampler for the supercell together with grid
supercell_ucs = sample.UnitCellSampler(supercell_from_unitcell_wo_ions)
supercell_ucs.n_frac = (nx, ny, nz)



### Set the method to be used # 10 struct first
if args.method:
    if args.method == '54209':
        energy_method = struct_54209_ff

    elif args.method == '54297':
        energy_method = struct_54297_ff

    elif args.method == '54449':
        energy_method = struct_54449_ff

    elif args.method == '54683':
        energy_method = struct_54683_ff

    elif args.method == '54837':
        energy_method = struct_54837_ff

    elif args.method == '54865':
        energy_method = struct_54865_ff

    elif args.method == '54879':
        energy_method = struct_54879_ff

    elif args.method == '54884':
        energy_method = struct_54884_ff

    elif args.method == '55184':
        energy_method = struct_55184_ff

    elif args.method == '55319':
        energy_method = struct_55319_ff

    elif args.method == '56568': # Three more struct below here
        energy_method = struct_56568_ff

    elif args.method == '56627':
        energy_method = struct_56627_ff

    elif args.method == '57382':
        energy_method = struct_57382_ff



else: 
    energy_method = structure_54189_ff_manually_coded
    print("Default energy method selected: ", str(energy_method))

print("Energy method used: ", str(energy_method))

# Compute energy grids below. Currently testing different things.

# Compute energies using reference supercell
included_energies_manually_coded_ff = supercell_ucs.calculate_energies(grid_points=supercell_cart_grid,
                                 method=energy_method,
                                 atom=args.a,
                                 exploit_symmetry=use_sym)

fill_value = np.nan_to_num(np.inf)
#energies_ref_sc = np.full(unitcell_included.shape, fill_value=fill_value, dtype=np.float64)
energies = np.full(unitcell_included.shape, fill_value=fill_value, dtype=np.float64)
print("unitcell_included", unitcell_included, np.sum(unitcell_included))
#energies_ref_sc[unitcell_included] = included_energies_manually_coded_ref_sc
energies[unitcell_included] = included_energies_manually_coded_ff

cube_out = '.'.join((calc_name, 'cube'))
xsf_out = '.'.join((calc_name, 'xsf')) # If ever xsf is wanted... 

cube_dir = Path("10_structures_cube_out")

if not (cube_dir.exists() and cube_dir.is_dir()):
    print("Cube directory doesn't exist, so I create it...")
    os.mkdir(cube_dir)

cube_out_path = Path(cube_dir, cube_out) # To store cubes in separate dir
#cube_ref_sc = "TEST_54189_ref_sc_ucs.cube"
#cube_sc_from_uc = "TEST_54189_sc_from_uc_ucs.cube"

with open(cube_out_path, "w") as cf:
    write_cube(cf, unitcell, data=energies.reshape(nx, ny, nz))


print("DONE.")
# TODO: import functions to define lammps potential,
# and use them to define a method to feed to the sampler
# Make sure to feed the sampler

##### Lammps method read from input > for now in a separate script
#lmp_file =
#lmp_dict = parser_lammps_mel_inp(lmp_file)
#lmp_inp_calc = lammps_method_from_data(**lmp_dict)

#####

