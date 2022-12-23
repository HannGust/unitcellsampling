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
from unitcellsampling import sample

# TODO: Possibly make this function able to read input data files 
#       for the grid shape/vdw/spacing parameters. Then the grid indices
#       could be directly reconstructed from that using UCS?

#       Also, make it an option to generate the cooridnates the same way
#       as if it generates them with UCS, from vdw/spacing parameters?


def indices_from_ucs(unitcell, grid_shape, vdw):
    """
    Retrieves the input parameters from the ucs coordinate generator
    and uses the mask and to obtain the sampled indices.
    """
    ucsampler = sample.UnitCellSampler(unitcell)
    dummy_grid, grid_mask = ucsampler.generate_grid_vectors(grid_shape, vdw_scale=vdw)
    
    # Make indices here
    all_indices = np.rint(list(it.product(
                  np.linspace(0, grid_shape[0], grid_shape[0], endpoint=False),
                  np.linspace(0, grid_shape[1], grid_shape[1], endpoint=False), 
                  np.linspace(0, grid_shape[2], grid_shape[2], endpoint=False)))).reshape(grid_shape + (3,))
    
    all_indices_flat = all_indices.reshape(-1,3)
    flat_mask = grid_mask.flatten()
    flat_mask = grid_mask.reshape(-1,1)

    # Debug
    print("all_indices shape:", all_indices.shape)
    print("all_indices initial elements ([0,0,0:4,:]): ", all_indices[0,0,0:4,:])  
    print("all_indices_flat initial 4 rows: ", all_indices_flat[0:4,:])


    indices_ucs = all_indices[grid_mask]
    indices_ucs_from_flat = all_indices_flat[flat_mask]
    assert (indices_ucs.reshape(-1,3) == indices_ucs_from_flat).all(), "Indices from full array and from flattened array and mask are not the same?!?!"
    
    print("Passed the test. Indices are the same.")
    #raveled_indices_ucs = np.ravel_multi_index((indices_ucs[:,0], indices_ucs[:,1], indices_ucs[:,2]), grid_shape)
    return None # raveled_indices, indices_ucs



def indices_from_coord(unitcell, coords_computed, grid_shape):
    """
    Determines the indices of the points that where sampled 
    from the grid size, unitcell and cartesian coordinates of the
    sampled points.
    """
    indices = np.rint(
        np.linalg.solve(unitcell.get_cell().T , coords_computed.T).T 
        * np.array(grid_shape)
        ).astype(dtype=int)

    raveled_indices  = np.ravel_multi_index((indices[:,0], indices[:,1], indices[:,2]), grid_shape)

    return raveled_indices

# TODO later: Note that shift and midvox is not implemented here - thats because the vdw check
# is preformed after coordiantes has been chosen when generated thru ucs

def grid_text_to_cube_ucs(energyfile:str,
                      unitcell, 
                      grid_shape, 
                      name="txt2cube.cube", 
                      fill_value=np.nan_to_num(np.inf), 
                      return_grids=False, 
                      write_file=True,
                      vdw=0.0
                      ):
    """
    Takes an energy grid with coordinates from textfiles
    and a unitcell and writes it to a cube file.
    Assumes that the energies are in kcal/mol,
    and converts them into eV.
    Uses UCS with given input parameters to determine
    the sampled grid points, i.e. which coordinates the energies
    correspond to.
    """
    
    # Read in energy values
    grid_values = np.loadtxt(energyfile)
    
    # Shift energy values so minimum is 0
    grid_values = grid_values - np.min(grid_values)
    
    # Then convert energies to eV from kcal/mol
    grid_values = grid_values * ase.units.kcal / ase.units.mol
    
    # Initialize grid
    #grid = np.empty((np.product(grid_shape)), dtype=np.float64)
    grid = np.empty(grid_shape, dtype=np.float64)
    grid.fill(fill_value)
    
    # Let UCS generate masking array from parameters
    if vdw is None:
        vdw = 0.0 # Need this parameter
    print("vwd= ",vdw)
    
    ucsampler = sample.UnitCellSampler(unitcell)
    dummy_grid, grid_mask = ucsampler.generate_grid_vectors(grid_shape, vdw_scale=vdw)
    
    grid_mask = np.array(grid_mask, dtype=bool)
    # debug
    print("printing grid, masked grid and energy values shapes.")
    print(grid.shape, grid[grid_mask].shape, grid_values.shape)
    assert len(grid_mask.reshape(-1)) == (len(grid_values)), "grid_mask not the same length as grid_values"
    
    grid[grid_mask] = grid_values

    #grid_mask_flat = np.array(grid_mask.reshape(-1), dtype=bool)
    #grid[grid_mask_flat] =  grid_values

    
    if write_file:
        with open(name+".cube", 'w') as f:
            write_cube(f, unitcell, data=grid.reshape(grid_shape))
        print("Writing energy grid to cubefile: ", name+".cube")

    if return_grids:
        return grid
    
    return None


##########

def grid_text_to_cube_from_coords(energyfile:str, 
                      coordfile:str, 
                      unitcell, 
                      grid_shape, 
                      name="txt2cube.cube", 
                      fill_value=np.nan_to_num(np.inf), 
                      return_grids=False, 
                      write_file=True,
                      test_coords=True,
                      shift=None,
                      midvox=False,
                      **kwargs):
    '''
    Takes an energy grid with coordinates from textfiles
    and a unitcell and writes it to a cube file.
    Assumes that the energies are in kcal/mol,
    and converts them into eV.
    Uses cartesian coordinates read from a file to 
    determine the indices of the energy values.
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
    
    # compute indices from coordinates => need to take into account midvoxel sampling
    # Take into account midvox and shift samplig

    if midvox:
        shift = 0.5 * (1/np.array((grid_shape)))

    if (shift is not None) or (shift != 0):
        
        if isinstance(shift, float):
            shift = np.array((shift,) * 3)
        
        elif len(shift) >= 3:
            shift = shift[0:3]
        
        else:
            raise Exception("Shift needs to be float or 1d array with len at least 3")
        
        print("shift: ", shift)
        coords_computed = coords_computed - shift @ unitcell.get_cell()

    indices = np.rint(
        np.linalg.solve(unitcell.get_cell().T , coords_computed.T).T 
        * np.array(grid_shape)
        ).astype(dtype=int)
    
    print(coords_computed.shape)
    print(coords_computed[0:3])
    print(indices.shape)
    print(indices[0:3])
    raveled_indices  = np.ravel_multi_index((indices[:,0], indices[:,1], indices[:,2]), grid_shape)

    
    if test_coords:
        assert (np.abs(np.array(coords)[raveled_indices] - np.array(coords_computed)) < 10**(-5)).all(), "Coords are different at 1.0e-5 absolute error tolerance."

    
    grid = grid.flatten()
    grid[raveled_indices] = grid_values
    
    if write_file:
        with open(name+".cube", 'w') as f:
            write_cube(f, unitcell, data=grid.reshape(grid_shape))
        print("Writing energy grid to cubefile: ", name+".cube") 

    if return_grids:
        return grid, coords, indices, raveled_indices




############

def grid_text_to_cube(energyfile:str, 
                      coordfile:str, 
                      unitcell, 
                      grid_shape, 
                      name="txt2cube.cube", 
                      fill_value=np.nan_to_num(np.inf), 
                      return_grids=False, 
                      write_file=True,
                      test_coords=True,
                      shift=None,
                      use_ucs_indices=False,
                      vdw=0.0,
                      midvox=False,
                      **kwargs):
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
    

    #if points_from_ucs:
    # Let UCS generate masking array from parameters
    if vdw is None:
        vdw = 0.0 # Need this parameter
    print("vwd= ",vdw)
    
    ucsampler = sample.UnitCellSampler(unitcell)
    dummy_grid, grid_mask = ucsampler.generate_grid_vectors(grid_shape, vdw_scale=vdw)
    
    #ucs_grid = np.copy(grid).reshape(grid_shape)
    #ucs_grid[grid_mask] = grid_values
    
    # Make indices here
    all_indices = np.rint(list(it.product(
                  np.linspace(0, grid_shape[0], grid_shape[0], endpoint=False),
                  np.linspace(0, grid_shape[1], grid_shape[1], endpoint=False), 
                  np.linspace(0, grid_shape[2], grid_shape[2], endpoint=False)))
                  ).astype(int).reshape(grid_shape + (3,))
    
    indices_ucs = all_indices[grid_mask]
    print(indices_ucs.shape)
    raveled_indices_ucs = np.ravel_multi_index((indices_ucs[:,0], indices_ucs[:,1], indices_ucs[:,2]), grid_shape)
    
    #else:
    # compute indices from coordinates => need to take into account midvoxel sampling
    indices_alt = np.rint(
        np.linalg.solve(unitcell.get_cell().T , coords_computed.T).T 
        * np.array(grid_shape)
        ).astype(dtype=int)
    
    print(grid_shape)
    print("indices_alt shape: ", indices_alt.shape)
    raveled_indices_alt  = np.ravel_multi_index((indices_alt[:,0], indices_alt[:,1], indices_alt[:,2]), grid_shape)

    assert (raveled_indices_ucs == raveled_indices_alt).all(), "UCS indices and alternative indices do not agree!!!"

    #print((np.array(coords)[raveled_indices])[:4])
    #print(np.array(coords_computed)[:4])
    if test_coords:
        if not midvox:
            if not shift:
                shift = 0.0
            a,b,c = unitcell.get_cell()[0,:], unitcell.get_cell()[1,:], unitcell.get_cell()[2,:]
            assert (np.abs(np.array(coords)[raveled_indices_alt] - (np.array(coords_computed) - shift * (a + b + c))) < 10**(-5)).all(), "Coords are different."
        else:
            a_mv,b_mv,c_mv = (1.0/grid_shape[0]) * unitcell.get_cell()[0,:], (1.0/grid_shape[1]) * unitcell.get_cell()[1,:], (1.0/grid_shape[2]) * unitcell.get_cell()[2,:]
            assert (np.abs(np.array(coords)[raveled_indices_alt] - (np.array(coords_computed) - (a_mv + b_mv + c_mv))) < 10**(-5)).all(), "Coords are different."


    if use_ucs_indices:
        raveled_indices = raveled_indices_ucs
    else:
        raveled_indices = raveled_indices_alt    
    
    grid = grid.flatten()
    grid[raveled_indices] = grid_values

    if write_file:
        with open(name+".cube", 'w') as f:
            write_cube(f, unitcell, data=grid.reshape(grid_shape))
        print("Writing energy grid to cubefile: ", name+".cube") 

    if return_grids:
        return grid, coords, all_indices, raveled_indices



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
    parser.add_argument('--vdw', type=float, action='store', default=0.0, help='The van der Waals radius scaling factor used in exclusion of points too close to atoms. Used when UCS is utulized for determining the indices of the energies.')
    
    parser.add_argument('--midvox', action='store_true', help='Specifies that the coords corresponds to the midpoints of voxels - this overrides the shift option.')
    parser.add_argument('-f', '--format', type=str, action='store', choices=['cube', 'ucs_cube', 'cif', 'xsf'], default='cube', help='Format of the output file')
    
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
        print("cube format, using coordinates from file, chosen.")
        grid_text_to_cube_from_coords(energyfile, 
                              coordfile, 
                              unitcell, 
                              grid_shape=grid_shape, 
                              name=args.name, 
                              write_file=True,
                              test_coords=args.noct,
                              shift=args.shift,
                              midvox=args.midvox)
    elif args.format == 'ucs_cube':
        print("cube format, using UCS is selected. Note that coordfile, shift and midvox are redundant.")
        grid_text_to_cube_ucs(energyfile, 
                              unitcell,
                              grid_shape=grid_shape, 
                              name=args.name,  
                              vdw=args.vdw)
    elif args.format == 'cif':
        grid_text_to_cif(energyfile, coordfile, unitcell, grid_shape=grid_shape, name=args.name)
    elif args.format == 'xsf':
        grid_text_to_xsf(energyfile, coordfile, unitcell, grid_shape=grid_shape, name=args.name)
    else:
        argparse.ArgumentError('Invalid format: ' + str(args.format))
