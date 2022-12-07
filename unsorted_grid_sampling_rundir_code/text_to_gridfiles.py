#!/usr/bin/env/ python

# Text to cube/cif/xsf
import numpy as np
import math as ma
from pathlib import Path
import ase, ase.io
from ase.io.cube import write_cube
import ase.build
import argparse
import itertools as it


def grid_text_to_cube(energyfile:str, 
                      coordfile:str, 
                      unitcell, 
                      grid_shape, 
                      name="txt2cube.cube", 
                      fill_value=np.nan_to_num(np.inf), 
                      return_grids=False, 
                      write_file=True,
                      test_coords=True,
                      shift=None):
    '''
    Takes an energy grid with coordinates from textfiles
    and a unitcell and writes it to a cube file.
    Assumes that the energies are in kcal/mol,
    and converts them into eV.
    '''
    grid_values = np.loadtxt(energyfile)
    coords_computed = np.loadtxt(coordfile)
    
    # Shift grid
    grid_values = grid_values - np.min(grid_values)
    
    # Then convert grid to eV from kcal/mol
    grid_values = grid_values * ase.units.kcal / ase.units.mol
    
    grid = np.empty((np.product(grid_shape)), dtype=np.float64)
    f_coords = np.array(list(
        it.product(
            np.linspace(0, 1, grid_shape[0], endpoint=False), 
            np.linspace(0, 1, grid_shape[1], endpoint=False), 
            np.linspace(0, 1, grid_shape[2], endpoint=False)
            )
            )) 

    coords = f_coords @ unitcell.get_cell()
    
    grid.fill(fill_value)

    # compute indices from coordinates
    indices = np.rint(
        np.linalg.solve(unitcell.get_cell().T , coords_computed.T).T 
        * np.array(grid_shape)
        ).astype(dtype=int)

    raveled_indices  = np.ravel_multi_index((indices[:,0], indices[:,1], indices[:,2]), grid_shape)


    #print((np.array(coords)[raveled_indices])[:4])
    #print(np.array(coords_computed)[:4])
    if test_coords:
        if not shift:
            shift = 0.0
        a,b,c = unitcell.get_cell()[0,:], unitcell.get_cell()[1,:], unitcell.get_cell()[2,:],
        assert (np.abs(np.array(coords)[raveled_indices] - (np.array(coords_computed) - shift * (a + b + c))) < 10**(-5)).all(), "Coords are different."

    grid = grid.flatten()
    grid[raveled_indices] = grid_values

    if write_file:
        with open(name+".cube", 'w') as f:
            write_cube(f, unitcell, data=grid.reshape(grid_shape))

    if return_grids:
        return grid, coords, indices, raveled_indices



def grid_text_to_cif():
    pass

def grid_text_to_xsf():
    pass


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('energy_grid', type=str, action='store', help='Textfile containing energies to read and convert into cube/cif/xsf file.')
    parser.add_argument('coord_file', type=str, action='store', help='Textfile containing coordinates to read and convert into cube/cif/xsf file.')
    parser.add_argument('atoms', type=str, action='store', help='File containing atoms object of the unit cell.')
    parser.add_argument('name', type=str, action='store', help='Name of the output file (not including ".cube").')
    parser.add_argument('-g', '--grid_size', type=int, nargs='+', action='store', help='The number of gridpoints in each direction.')
    parser.add_argument('-s', '--spacing', type=float, nargs='+', action='store',
                        help='The spacing between adjacent gridpoints in each\
                             direction. Computes the grid shape from this in \
                             the same')
    parser.add_argument('--noct', '--nocoordtest', action='store_false',
                        help='Whether to make a test to see that correct\
                             indices have been computed. Defualt True.\
                             May not work for shifted coordinates.')
    parser.add_argument('--shift', type=float, action='store', default=None,
                        help='The shift that has been applied to the\
                             coordinates, in terms of the fraction of\
                             a+b+c that has been added. Only relevant\
                             for coordinate testing. Default is no shift.')
    
    
    parser.add_argument('-f', '--format', type=str, action='store', choices=['cube', 'cif', 'xsf'], default='cube', help='Format of the output file')
    
    args = parser.parse_args()
    
    energyfile = Path(args.energy_grid)
    coordfile = Path(args.coord_file)
    unitcell = ase.io.read(args.atoms)

    # Grid
    ## Set the number of points in the grid in each dimension (or equivalently, the mesh size)
    if not args.spacing:
        if len(args.grid_size) == 1:
            nx,ny,nz = args.grid_size * 3
        else:
            nx,ny,nz = args.grid_size[0:3]
        
        a,b,c = np.linalg.norm(unitcell.get_cell()[0,:]),\
                np.linalg.norm(unitcell.get_cell()[1,:]),\
                np.linalg.norm(unitcell.get_cell()[2,:])
        true_spacing = (a/nx, b/ny, c/nz)
        grid_shape = (nx, ny, nz)
        print('True spacing: ', true_spacing)
    else:
        # Small section to compute nx, ny, nz i.e. the # of points in each direction in case the --space option is used
        # This however is not needed - sampler.generate_grid_vectors has this functionality already, so we can just pass arguments to this
        # Although, that option does parheps something odd - it might ue interpolation to improve the spacing... 
    
        if len(args.spacing) == 1:
            spacing_x, spacing_y, spacing_z = args.spacing * 3
        else:
            spacing_x, spacing_y, spacing_z =  args.spacing[0:3]
    
        a,b,c = np.linalg.norm(unitcell.get_cell()[0,:]), np.linalg.norm(unitcell.get_cell()[1,:]), \
                np.linalg.norm(unitcell.get_cell()[2,:])
    
        nx,ny,nz = ma.ceil(a/spacing_x), ma.ceil(b/spacing_y), ma.ceil(c/spacing_z)
        true_spacing = (a/nx, b/ny, c/nz)
        grid_shape = (nx, ny, nz)
        print('True spacing: ', true_spacing)



    if args.format == 'cube':
        grid_text_to_cube(energyfile,
                          coordfile,
                          unitcell,
                          grid_shape=grid_shape,
                          name=args.name,
                          test_coords=args.noct,
                          shift=args.shift)
    elif args.format == 'cif':
        grid_text_to_cif(energyfile, coordfile, unitcell, grid_shape=grid_shape, name=args.name)
    elif args.format == 'xsf':
        grid_text_to_xsf(energyfile, coordfile, unitcell, grid_shape=grid_shape, name=args.name)
    else:
        argparse.ArgumentError('Invalid format: ' + str(args.format))