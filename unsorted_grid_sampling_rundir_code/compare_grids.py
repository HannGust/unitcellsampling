

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

def compute_diff_grid(grid1, grid2):
    """Calculates the difference between two grids."""
    grid1 = np.array(grid1)
    grid2 = np.array(grid2)
    assert grid1.shape == grid2.shape, "Grids must have the same shape!"

    return np.subtract(grid1, grid2)

def compute_abs_diff_grid(grid1, grid2):
    return np.abs(compute_diff_grid(grid1, grid2))

def compute_grid_rmsd(grid1, grid2):
    return np.sqrt(np.mean(np.square(compute_diff_grid(grid1, grid2))))

def compute_grid_rel_error(grid1, grid2):
    """Computes number of NaNs, Inf and otherwise computes relative difference between grids."""
    rel_err_grid = np.empty(grid1.shape)
    #print("Is None 1:", rel_err_grid is None)
    rel_err_grid.fill(np.nan_to_num(np.inf))
    #print("Is None 2:", rel_err_grid is None) 
    diff_grid = compute_diff_grid(grid1, grid2)
    #print("Is None 3:", diff_grid is None)
    nonzero_grid1 = grid1[grid1 != 0.0] 

    #print("Is None 4:", nonzero_grid1 is None)
    rel_err_grid[diff_grid == 0.0] = 0.0
    rel_err_grid[grid1 != 0.0] = np.abs(diff_grid[grid1 != 0.0]/grid1[grid1 != 0.0])
    #print("Is None 5:", rel_err_grid is None)
    assert rel_err_grid.shape == grid1.shape, "Relative error grid has different shape from input grid! Something went wrong!"

    return rel_err_grid

def compare_grids(grid1, grid2):
    assert len(grid1.shape) == 1 or len(grid1.shape) == 3, "Grid 1 must be either of shape (n,) or (n,m,l)."
    assert len(grid2.shape) == 1 or len(grid2.shape) == 3, "Grid 2 must be either of shape (n,) or (n,m,l)."   
    assert grid1.shape == grid2.shape, "Grids must have the same shape."

    diff_grid = compute_diff_grid(grid1, grid2)
    abs_diff_grid = compute_abs_diff_grid(grid1, grid2)
    rel_err_grid = compute_grid_rel_error(grid1, grid2) 
    rmsd = compute_grid_rmsd(grid1, grid2)

    return (diff_grid, abs_diff_grid, rel_err_grid, rmsd)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("files", type=str, nargs='+', help="Give the list of files from which to read the grids that should be compared. Observe: The grids needs to contain common coordinates!")
    parser.add_argument("-o", "--allout", action="store_true", help="If given, writes the difference and absolute difference (grid_1 - grid_2) at the common coordinates into cube files.")
    parser.add_argument("--cc", "--customcube", action="store_true", help="Whetherto use a custom cube file reader, instead of ase's cube reader which is the default.")
    parser.add_argument("--do", "--diff-out", action="store_true", help="If given, writes the difference at the common coordinates to a cube file. ")    
    parser.add_argument("--ao", "--abs-out", action="store_true", help="If given, writes the absolute difference at the common coordinates to a cube file. ")    
    parser.add_argument("--ro", "--rel-out", action="store_true", help="If given, writes the realative difference at the common coordinates to a cube file. ") 
    parser.add_argument("-s", "--simple", action="store_true", help="If given, does a very straightforward comparison of the two grids. Requires that they have the same shape. ") 

    args = parser.parse_args()
    if not args.simple:
        grids = []
        grid_coords = []
        nx_lst = []
        ny_lst = []
        nz_lst = []
            
        for i,infile in enumerate(args.files):
            print(infile)
            if args.cc:
                (atoms, grid) = custom_read_cube(infile)
            else:
                with open(infile, "r") as cf:
                    cube_cont = read_cube(cf)
                    atoms = cube_cont["atoms"]
                    grid = cube_cont["data"]

            grids.append(grid)
            print("Grid shape: ", grid.shape)
            nx,ny,nz = np.array(grid).shape
            nx_lst.append(nx)
            ny_lst.append(ny)
            nz_lst.append(nz)
            grid_coords.append((list(it.product(np.linspace(0, 1, nx, endpoint=False), np.linspace(0, 1, ny, endpoint=False), np.linspace(0, 1, nz, endpoint=False)))))
            #print('First three grid coords: ', grid_coords[i][0:3])
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
        print("shift (at [0,0,0] coodrinate):",shift) 
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
        print('AVERAGE(ELEMENTWISE_ABSOLUTE_ERROR) = ', np.mean(np.abs(grid_1-grid_2)))
        print('')
        # Relative error:
        if ((np.abs((grid_1-grid_2)[grid_1 != 0]/grid_1[grid_1 != 0]) <= grids_within_rel_precision).all()): # checks for relative error below threshold
            print('The grids are identical within elementwise RELATIVE ERROR of: ', grids_within_rel_precision)
        else:
            print('The grids are different! In the following sense:  MAX(ELEMENTWISE_RELATIVE_ERROR) > ', grids_within_rel_precision)
        
        print('')
        print('MAX(ELEMENTWISE_RELATIVE_ERROR) = ', np.max(np.abs((grid_1-grid_2)[grid_1 != 0]/grid_1[grid_1 != 0])))
        print('')
        print('AVERAGE(ELEMENTWISE_RELATIVE_ERROR) = ', np.mean(np.abs((grid_1-grid_2)[grid_1 != 0]/grid_1[grid_1 != 0])))
        
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
    
        diff_out_name = "_".join(("difference",) + tuple(Path(file).stem for file in args.files) + (".cube",))
        abs_diff_out_name = "_".join(("abs", diff_out_name))
        rel_diff_out_name = "_".join(("rel", diff_out_name))
     
        difference_grid = grid_1 - grid_2
        abs_difference_grid = np.abs(grid_1 - grid_2)
     
        rel_difference_grid = np.empty(grid_1.shape)
        rel_difference_grid.fill(np.nan_to_num(np.inf))
        rel_difference_grid[grid_1 != 0] = np.abs((grid_1-grid_2)[grid_1 != 0]/grid_1[grid_1 != 0])
        rel_difference_grid[(grid_1-grid_2) == 0] = 0.0 
        rel_diff_reference_1 = compute_grid_rel_error(grid_1, grid_2) 
     
        #rel_diff_reference_2 = compute_grid_rel_error(grids[0], grids[1])
     
    
        if args.allout or args.do:
            with open(diff_out_name, "w") as f1:
                write_cube(f1, atoms, difference_grid)
    
        if args.allout or args.ao:
            with open(abs_diff_out_name, "w") as f2:
                write_cube(f2, atoms, abs_difference_grid)
     
        if args.allout or args.ro:
            assert np.all(rel_diff_reference_1 == rel_difference_grid), "Rel diff not equal to reference 1!!!"
            #assert np.all(rel_diff_reference_2 == rel_difference_grid), "Rel diff not equal to reference 2!!!"
            #assert np.all(rel_diff_reference_1 == rel_diff_reference_2), "Reference 1 not equal to reference 2!!!"
     
            with open(rel_diff_out_name, "w") as f3:
                write_cube(f3, atoms, rel_difference_grid)
     
        
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
        
    else:
        grids = []

        for i,infile in enumerate(args.files):
            print(infile)
            if args.cc:
                (atoms, grid) = custom_read_cube(infile)
            else:
                with open(infile, "r") as cf:
                    cube_cont = read_cube(cf)
                    atoms = cube_cont["atoms"]
                    grid = cube_cont["data"]

            grids.append(grid)
            print("Grid shape: ", grid.shape)
        
        grid1 = grids[0] 
        grid2 = grids[1]

        # Apply energy filtering if desired:
        e_thresh = 10.0**3
        energy_grid_mask = grid1 < e_thresh
        print("Applying energy mask of", e_thresh, "based on grid 1.")
        print("Total points:", grid1.size)
        print("Num points above/below thresh:", np.sum(energy_grid_mask), energy_grid_mask.size - np.sum(energy_grid_mask))
        #
        print("Sanity check, grid maxima after filtering:", np.max(grid1[energy_grid_mask]), np.max(grid2[energy_grid_mask]))
        print("")
         
        diff_grid_1, abs_diff_grid_1, rel_diff_grid_1, rmsd = compare_grids(grid1, grid2)
        #diff_grid, abs_diff_grid, rel_diff_grid, rmsd = compare_grids(grid1[energy_grid_mask], grid2[energy_grid_mask])
        diff_grid, abs_diff_grid, rel_diff_grid = diff_grid_1[energy_grid_mask], abs_diff_grid_1[energy_grid_mask], rel_diff_grid_1[energy_grid_mask]

        masked_rmsd = compute_grid_rmsd(grid1[energy_grid_mask], grid2[energy_grid_mask])
        print()
        print("Grid difference stats: ")
        print("RMSD: ", rmsd)
        print("RMSD, (energy masked): ", masked_rmsd)
        print()
        print("Mean diff: ", np.mean(diff_grid))
        print("Max diff: ", np.max(diff_grid))
        print("Min diff: ", np.min(diff_grid))
        print()
        print("Mean abs diff: ", np.mean(abs_diff_grid))
        print("Max abs diff: ", np.max(abs_diff_grid))
        print("Min abs diff: ", np.min(abs_diff_grid))
        print()
        print("Mean rel diff: ", np.mean(rel_diff_grid))
        print("Max rel diff: ", np.max(rel_diff_grid))
        print("Min rel diff: ", np.min(rel_diff_grid))
 

        diff_out_name = "_".join(("difference",) + tuple(Path(file).stem for file in args.files) + (".cube",))
        abs_diff_out_name = "_".join(("abs", diff_out_name))
        rel_diff_out_name = "_".join(("rel", diff_out_name))
 

        if args.allout or args.do:
             with open(diff_out_name, "w") as f1:
                 write_cube(f1, atoms, diff_grid)
  
        if args.allout or args.ao:
             with open(abs_diff_out_name, "w") as f2:
                 write_cube(f2, atoms, abs_diff_grid)
     
        if args.allout or args.ro:
             #assert np.all(rel_diff_reference_1 == rel_difference_grid), "Rel diff not equal to reference 1!!!"
             #assert np.all(rel_diff_reference_2 == rel_difference_grid), "Rel diff not equal to reference 2!!!"
             #assert np.all(rel_diff_reference_1 == rel_diff_reference_2), "Reference 1 not equal to reference 2!!!"
     
             with open(rel_diff_out_name, "w") as f3:
                 write_cube(f3, atoms, rel_diff_grid)
       
