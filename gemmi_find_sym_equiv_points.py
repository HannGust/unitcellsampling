

import gemmi
import ase
import ase.io
import numpy as np
import itertools as it

from ase.spacegroup import get_spacegroup
import argparse 


def find_symmetrically_equivalent_indices(grid_shape, atoms, spacegroup_no=None):
    """ Finds the indices of symmetrically equivalent points, i.e. the
        equivalence classes of points due to symmetry, in erms of their
        indices. Returns a list of indices for each equivalence class."""

    tot_bol_grid = gemmi.Int8Grid()

    if not spacegroup_no:
        spgroup = gemmi.find_spacegroup_by_number(get_spacegroup(atoms))
    else: 
        spgroup = gemmi.find_spacegroup_by_number(spacegroup_no)

    tot_bol_grid.spacegroup = spgroup

    tot_bol_grid.set_size(*grid_shape)


    eq_class_list = [] 

    for indx in it.product(range(grid_shape[0]), range(grid_shape[1]), range(grid_shape[2])):
        if tot_bol_grid.get_value(*indx) != 0:
            continue
        else:
            # Inititalize a new temporary integer grid with all 0 entries            
            tmp_bol_grid = gemmi.Int8Grid()
            tmp_bol_grid.spacegroup = spgroup
            tmp_bol_grid.set_size(*grid_shape)
            
            # Set the value of the current point to 1, then symmetrize
            tmp_bol_grid.set_value(*indx, 1)
            tmp_bol_grid.symmetrize_max() 

            # Do the same for the total grid (to keep track of which points have been considered)
            tot_bol_grid.set_value(*indx, 1)
            tot_bol_grid.symmetrize_max()

            # Find the list (array) of indices for nonzero entries in the temporary grid, append it to the list of equivalence classes
            eq_class_list.append(np.argwhere(tmp_bol_grid.array))

    return eq_class_list


def evaluate_equivalence_classes(grid, eq_class_list):
    eq_class_values_list = []

    # Make sure grid is an np.array
    if isinstance(grid, gemmi.FloatGrid):
        grid_as_array = np.array(grid.array)
    else:
        grid_as_array = np.array(grid)

    # For all lists of indices in the given equivalence class list, return those values from the grid array by fancy indexing, then add this array to the result list
    for lst in eq_class_list:
        eq_class_values_list.append(grid_as_array[lst[:,0], lst[:,1], lst[:,2]])

    return eq_class_values_list




def main():
    parser = argparse.ArgumentParser(description="Finds and returns the equivalent indices of a grid with a certain shape and spacegroup.")
    parser.add_argument("atoms", action='store', help="File containing the unitcell, from which the spacegroup symmetry is determined.")
    parser.add_argument("--gs", "--grid_shape", action='store', type=int, nargs='+', default=None, help="One or three integers specifying the number of points in each direction of the grid, i.e. the grid shape.")
    parser.add_argumemt("-s", "--spacing", action='store', type=float, nargs='+', default=None, help="One or three floats to specifythe grid spacing in Å. This is used to determine the grid shape, in the same manner as in the grid sampler.")
    parser.add_argument("--gsfc", "--gs-from-cube", action='store', default=None, help="Specifies which cube file to read the grid shape from. If neither this, -s or --gs are given, it will try to do it from the atoms-file input as if it was a cube-file.")

    args = parser.parse_args()

    atoms = ase.io.read(args.atoms)

    if (not args.gs and not args.s):
        # If cube file is given, read grid shape from that
        # Else try the atoms file.
        if args.gsfc:    
            with open(args.gsfc, 'r') as f:
                data = ase.io.cube.read_cube(f, read_data=True)
                grid = np.array(data['data'])
                print(grid)
                
        else:
            with open(args.atoms, 'r') as f:
                data = ase.io.cube.read_cube(f, read_data=True)
                grid = np.array(data['data'])
                print(grid)

        # Set the grid shape
        grid_shape = grid.shape

    elif args.gs and args.s:
        raise Exception("Grid shape and grid spacing are mutually exclusive arguments.")
    else:
        # Nothing in particular

    # Set grid shape to the specified grid shape if given
    if args.gs:
        grid_shape = args.gs

    # Determine grid shape
    if args.s:
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



 



    #grid = gemmi.FloatGrid(data.astype(np.float32))
    #grid.spacegroup = gemmi.get_spacegroup_from_number(get_spacegroup(atoms))

    equiv_indices = find_symmetrically_equivalent_indices(args.grid_shape, atoms) 
    
    print("Symmetry equivalence classes, and points in them: ")
    for eq in equiv_indices:
        print()
        print(np.array2string(eq))

    eq_class_values = evaluate_equivalence_classes(grid, equiv_indices)
    print()
    print("Values fo the points in each equivalence class: ")
    print()
    print(eq_class_values)

if __name__ == '__main__':
    main()
