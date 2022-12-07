
import numpy as np
import math as ma
import argparse
import ase, ase.io
import os
from pathlib import Path
from unitcellsampling import sample
from preparatory_fcns import remove_nonframework_cations_fancy


#grid_coord_dir = os.getenv("COORD_GRID_GEN_WD")
#
#if (not grid_coord_dir) or (not Path(grid_coord_dir).is_dir()):
#    print("This script utilized environ")
#    grid_coord_dir = "GRID_COORD_WD"

# Parser
parser = argparse.ArgumentParser(description='Generate grid coordinates (cartesian) for sampling of unit cells, with unitcellsampling.')
parser.add_argument('unitcell', type=str, action='store', help="Path to the cif-file containing the unit cell of the system (not supercell!).")
parser.add_argument('-n','--name', metavar='jobname', type=str, action='store', default=None, help="Desired name for the calculation (applied to generated output files, directories, etc.).")
parser.add_argument('-a', type=str, action='store', help="Atom to sample with. Default: Li", default="Li")
parser.add_argument('-g', '--grid', type=int, action='store', default=[10], nargs='+', help="Specify the number of grid points in each dimension (or cubic grid) (mutually exclusive with \"--space\").")
parser.add_argument('-s', '--space', type=float, action='store', default=None, nargs='+', help="Specify the spacing between the grid points in each dimension (mutually exclusive with \"--grid\").")
parser.add_argument('--vdw', type=float, action='store', default=0.01, help="Specify the fraction of the van der Waals radius that should be excluded from the sampled volume around each atom in the host structure.")
parser.add_argument('--shift', type=float, action='store', default=None, help="Specify the fraction of the fraction of the unit cell vectors by which the gridpoints should be shifted. \
    This is intended to solve some problems with lammps complaining about points not being within the boundaries of the given cell. Default is no shift.")
#parser.add_argument('--autoshift', type=float, action='store', default=None, help="Automatically determines the shift that needs to be applied to the coordinates for all of them to be . \


args = parser.parse_args()

sample_ion = ase.Atom(args.a)

# Handle the unitcell input file and name
unitcell_file = Path(args.unitcell)
if not unitcell_file.exists():
    raise FileNotFoundError('The inputfile \'' + str(unitcell_file) + '\' was not found. Need file containing unitcell.')

input_basename = unitcell_file.name

while len(Path(input_basename).suffixes) > 0:
    input_basename = Path(input_basename).stem # input_basename is here meant to be just the name of the input file without parent directories and extensions: /foo/bar.cif -> input_basename = bar
print("Input basename: ", input_basename)

if args.name:
    name = "_".join((input_basename, args.name))
else:
    name = input_basename

unitcell = ase.io.read(args.unitcell)

unitcell_wo_ions = remove_nonframework_cations_fancy(unitcell, sample_ion)

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


unitcell_ucs = sample.UnitCellSampler(unitcell_wo_ions)
unitcell_grid, unitcell_included = unitcell_ucs.generate_grid_vectors((nx, ny, nz), vdw_scale=args.vdw)

unitcell_grid = unitcell_grid.reshape(-1,3)
unitcell_included = unitcell_included.reshape(-1)

unitcell_included_bool = np.array(unitcell_included, dtype=bool)
unitcell_grid_included = unitcell_grid[unitcell_included_bool]

# Shift the grid is desired
if args.shift:
    unitcell_grid_included = unitcell_grid_included + args.shift * (np.sum(unitcell.get_cell(), 1))

# TODO: Add some check to see whether the coordinates are within the cell
#       boundaries or not.


# filenames with paths
if not Path(input_basename).is_dir():
    os.mkdir(Path(input_basename))

coords_txt_name = name+"_coords.txt"
input_data_name = name+"_coord_input_data.txt"
path_to_input_data_file = input_basename+"/"+input_data_name
path_to_coord_txt_file = input_basename+"/"+coords_txt_name

with open(path_to_input_data_file, 'w') as f:
    f.write("Structure: "+input_basename+"\n")
    f.write("Unitcell file path: "+args.unitcell+"\n")

    if args.space:
        f.write("Spacing: "+str((spacing_x, spacing_y, spacing_z))+"\n")
    
    f.write("Grid shape: "+str((nx, ny, nz))+"\n")
    f.write("Num points: "+str(nx*ny*nz)+"\n")
    f.write("VdW parameter: "+str(args.vdw)+"\n")
    
    f.write("Is grid shifted: "+str(bool(args.shift))+"\n")
    if args.shift:
        f.write("Grid shift fraction: "+str(args.shift))

np.savetxt(path_to_coord_txt_file, unitcell_grid_included)

