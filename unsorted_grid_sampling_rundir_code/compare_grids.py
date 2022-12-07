

import ase
from ase.io import read
from ase.io.cube import read_cube, write_cube
import numpy as np
import math as ma
import argparse
import itertools as it
from pathlib import Path

# Reads atoms object and data from a cubefile written in a wierd way by ase so that ase can't read the data.
def custom_read_cube(filename: str) -> tuple:
    filename = Path(filename)

    atoms = read(filename)
    with open(filename,'r') as f:
        contents = f.readlines()

    #print("Debug/test statements:") 
    #print("Total number of lines: ", len(contents))
    #print("Number of lines to skip (until grid data): ", 2 + 4 + len(atoms))
    #print("Number of lines of grid data: ", len(contents) - (2 + 4 + len(atoms)))
    #print("Line 3: ", contents[2])
    #print("Line 4: ", contents[3])

    nx = int(contents[3].split()[0])
    ny = int(contents[4].split()[0])
    nz = int(contents[5].split()[0])
    print("Grid size: ", (nx, ny, nz))
    grid = np.array(list(map(float,contents[2 + 4 + len(atoms):]))).reshape((nx,ny,nz))

    print(grid.shape)

    return (atoms, grid)

def compare_txt_lists(file1, file2, out=False):
    """ 
    Reads shape (n,) or (n, 3) arrays from textfiles
    and computes the difference, absolute difference,
    maximum absolute difference and root mean square deviation.
    """
    grid1 = np.loadtxt(file1)
    grid2 = np.loadtxt(file2)
    
    assert grid1.shape == grid2.shape, "Grid shapes must be the same"

    if len(grid1.squeeze().shape) == 1:
        print("grid shape is (n,)")
        diff = grid1 - grid2
        abs_diff = np.abs(diff)
        max_abs_diff = np.max(abs_diff)
        rmsd = np.sqrt(np.mean(np.square(diff)))
    
    elif len(grid1.squeeze().shape) == 2 and grid1.squeeze().shape[1] == 3:
        print("grid shape is (n,...)")
        diff = grid1 - grid2
        abs_diff = np.linalg.norm(diff, axis=1)
        max_abs_diff = np.max(abs_diff)
        rmsd = np.sqrt(np.mean(np.square(abs_diff)))

    return (diff, abs_diff, max_abs_diff, rmsd)




if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("files", type=str, nargs='+', help="Give the list of files from which to read the grids that should be compared. Observe: The grids needs to contain common coordinates!")
    parser.add_argument("-o", "--out", action="store_true", help="If given, writes the difference and absolute difference (grid_1 - grid_2) at the common coordinates into cube files.")

    args = parser.parse_args()
    
    grids = []
    grid_coords = []
    nx_lst = []
    ny_lst = []
    nz_lst = []
        
    for i,infile in enumerate(args.files):
        print(infile)
        (atoms, grid) = custom_read_cube(infile)
        grids.append(grid)
        print("Grid shape: ", grid.shape)
        nx,ny,nz = np.array(grid).shape
        nx_lst.append(nx)
        ny_lst.append(ny)
        nz_lst.append(nz)
        grid_coords.append((list(it.product(np.linspace(0, 1, nx, endpoint=False), np.linspace(0, 1, ny, endpoint=False), np.linspace(0, 1, nz, endpoint=False)))))
        print('First three grid coords: ', grid_coords[i][0:3])
        #print(len(grid_coords[i]))
        # for now, assume its grids within the same unit cell
        
        
    grid_sets = list(map(set,grid_coords))
    
    grid_intersect = grid_sets[0].intersection(*grid_sets[1:])
    grid_intersect = np.array(list(grid_intersect))
            
        # Begin easy - take two grids
    gx = ma.gcd(nx_lst[0],nx_lst[1])
    gy = ma.gcd(ny_lst[0],ny_lst[1])
    gz = ma.gcd(nz_lst[0],nz_lst[1])
    
    """kx1 = np.linspace(0, shape_list[0][0], 1, endpoint=False)
    kx2 = np.linspace(0, shape_list[1][0], 1, endpoint=False)
    
    ky1 = np.linspace(0, shape_list[0][1], 1, endpoint=False)
    ky2 = np.linspace(0, shape_list[1][1], 1, endpoint=False)
    
    kz1 = np.linspace(0, shape_list[0][2], 1, endpoint=False)
    kz2 = np.linspace(0, shape_list[1][2], 1, endpoint=False)"""
    
    # grids for comparison
    print(round((nx_lst[0]/gx)),
    round((nx_lst[1]/gx)),
    round((ny_lst[0]/gy)),
    round((ny_lst[1]/gy)),
    round((nz_lst[0]/gz)),
    round((nz_lst[1]/gz)))
    
    
    
    grids_within_precision = 1.0e-1
    grids_within_rel_precision = 1.0e-06
    coord_num_decimals = 8
    
    
    
    shift = grids[0][0,0,0] - grids[1][0,0,0]
    
    grid_1 = grids[0][0:nx_lst[0]:round(nx_lst[0]/gx),0:ny_lst[0]:round(ny_lst[0]/gy),0:nz_lst[0]:round(nz_lst[0]/gz)]
    grid_2 = grids[1][0:nx_lst[1]:round(nx_lst[1]/gx),0:ny_lst[1]:round(ny_lst[1]/gy),0:nz_lst[1]:round(nz_lst[1]/gz)]
    grid_1 = grid_1 - shift
    
    coords_1 = np.round(np.array(grid_coords[0]).reshape((nx_lst[0], ny_lst[0], nz_lst[0], 3))[
        0:nx_lst[0]:round(nx_lst[0]/gx),
        0:ny_lst[0]:round(ny_lst[0]/gy),
        0:nz_lst[0]:round(nz_lst[0]/gz)], decimals=coord_num_decimals)
    
    coords_2 = np.round(np.array(grid_coords[1]).reshape((nx_lst[1], ny_lst[1], nz_lst[1], 3))[
        0:nx_lst[1]:round(nx_lst[1]/gx),
        0:ny_lst[1]:round(ny_lst[1]/gy),
        0:nz_lst[1]:round(nz_lst[1]/gz)], decimals=coord_num_decimals)
    
    
    print('Grid shapes: ', grid_1.shape, grid_2.shape)
    print('Coord shapes: ', coords_1.shape, coords_2.shape)
    print('')
    print('')
    
    ### Comparing grids up to some precision
    
    #  Absolute error:
    if ((np.abs(grid_1-grid_2) <= grids_within_precision).all()): # checks for absolute error below threshold
        print('The grids are identical within elementwise ABSOLUTE ERROR of: ', grids_within_precision)
    else:
        print('The grids are different! In the following sense:  MAX(ELEMENTWISE_ABSOLUTE_ERROR) > ', grids_within_precision)
    print('')
    print('MAX(ELEMENTWISE_ABSOLUTE_ERROR) = ', np.max(np.abs(grid_1-grid_2)))
    print('')
    
    print('')
    # Relative error:
    if ((np.abs((grid_1-grid_2)/grid_1) <= grids_within_rel_precision).all()): # checks for relative error below threshold
        print('The grids are identical within elementwise RELATIVE ERROR of: ', grids_within_rel_precision)
    else:
        print('The grids are different! In the following sense:  MAX(ELEMENTWISE_RELATIVE_ERROR) > ', grids_within_rel_precision)
    
    print('')
    print('MAX(ELEMENTWISE_RELATIVE_ERROR) = ', np.max(np.abs((grid_1-grid_2)/grid_1)))
    
    print('')
    print('')
    
    # Comparing coordinates, i.e. that the correct points (the same points in the unit cell) are compared
    if (coords_1 == coords_2).all():
        print('The coordinates are identical. (Precision: rounded to  ', coord_num_decimals,'  decimal places)')
    else:
        print('The coordinates are different!!?!')
        for i,c in enumerate(coords_1.flatten()):
            if coords_1.flatten()[i] != coords_2.flatten()[i]:
                print("Coords differ at: ", i, coords_1.flatten()[i], coords_2.flatten()[i])



    # write difference between grids to a cube file    
    if args.out:
        diff_out_name = "_".join(("difference",) + tuple(Path(file).stem for file in args.files) + (".cube",))
        abs_diff_out_name = "_".join(("abs", diff_out_name))

        difference_grid = grid_1 - grid_2
        abs_difference_grid = np.abs(grid_1 - grid_2)

        with open(diff_out_name, "w") as f1:
            write_cube(f1, atoms, difference_grid)

        with open(abs_diff_out_name, "w") as f2:
            write_cube(f2, atoms, abs_difference_grid)



    
    # Print statements for debugging and testing purposes gathered below
    #print(grid_1.flatten())
    #print(grid_2.flatten())
    #print(coords_1.flatten()[:])
    #print(coords_2.flatten()[:])
    #print(coords_1,coords_2)
    #print(coords_1.shape,coords_2.shape)
    
    #print("Printing the 2 lists of coodinates: ")
    #print(list(map(list,grid_coords[0])))
    #print(list(map(list,grid_coords[1])))
    #print("Grid intersection: ", list(map(list,grid_intersect)))
    
    
