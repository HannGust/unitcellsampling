
import numpy as np
import ase.io
from ase.io import read
from ase.io.cube import read_cube
from ase.io.cube import write_cube
#from ase.io.xsf import write_xsf
import argparse
from pathlib import Path

### Error functions

def calc_abs_error(grid1, grid2):
    abs_error = np.copy(np.abs(np.array(grid1) - np.array(grid2)))
    return abs_error

def calc_rel_error(test_grid, ref_grid, divide_by_ref=True):
    if divide_by_ref:
       rel_error = np.zeros(ref_grid.shape)
       rel_error[np.logical_not(np.isclose(ref_grid, 0.0))] = np.copy(np.abs(np.divide(
           np.array(ref_grid[np.logical_not(np.isclose(ref_grid, 0.0))]) - np.array(test_grid[np.logical_not(np.isclose(ref_grid, 0.0))]), 
           np.array(ref_grid[np.logical_not(np.isclose(ref_grid, 0.0))]))))

       rel_error[np.logical_and(np.isclose(ref_grid, 0.0), np.isclose(test_grid, 0.0))] = 0.0

       rel_error[np.logical_and(np.isclose(ref_grid, 0.0), np.logical_not(np.isclose(test_grid, 0.0)))] = np.nan
    else:
       rel_error = np.zeros(test_grid.shape)
       rel_error[np.logical_not(np.isclose(test_grid, 0.0))] = np.copy(np.abs(np.divide(
           np.array(test_grid) - np.array(ref_grid), 
           np.array(test_grid))))

       rel_error[np.logical_and(np.isclose(ref_grid, 0.0), np.isclose(test_grid, 0.0))] = 0.0

       rel_error[np.logical_and(np.isclose(test_grid, 0.0), np.logical_not(np.isclose(ref_grid, 0.0)))] = np.nan

    return rel_error

def calc_RMSD(grid1, grid2):
    RMSD = np.sqrt(np.mean(np.square(np.array(grid1) - np.array(grid2))))
    return RMSD

###

def main():
    arg_parser = argparse.ArgumentParser()
    
    arg_parser.add_argument('grids', type=str, action='store', nargs='+', help='Specify grids to compare/compute difference between.')
    arg_parser.add_argument('-co', '--cubeout', action='store_true', help='Specifies that error outputs should be written to a cube file.')
    arg_parser.add_argument('-xo', '--xsfout', action='store_true', help='Specifies that error outputs should be written to a xsf file.')
    arg_parser.add_argument('-nb', '--numbins', type=int, action='store', help='Specifies the number of energy intervals ("bins"), in which to separately compute errors. Done both in linear and log scale.')
    arg_parser.add_argument('-sh', '--shift', action='store', default=0.0, help='Specifies the shift in energy to be applied to grid 2. This is to alleviate the effect of normalization of the grids, i.e. to put them on the same energy scale. Can be float or \'auto\', in which case it is determined by comparing a pair of equivalent points, and letting the shift be the difference between those. Default: 0.0.')
    arg_parser.add_argument('-emax', action='store', type=float, default=None, help='Whether to apply a maximum energy threshold, that is, only consider the points below this energy.')
    arg_parser.add_argument('-emin', action='store', type=float, default=None, help='Whether to apply a minimum energy threshold, that is, only consider the points above this energy.')
    
    #arg_parser.add_argument('--factor')
    
    args = arg_parser.parse_args()
    
    cube_out = args.cubeout
    xsf_out = args.xsfout
    
    if len(args.grids) < 2:
        raise Exception("Need at least 2 arguments.")
    
    filenames = args.grids
    
    #"for infile in filenames:
    #   atoms_list.append(read(infile))
    #   grids.append(read(infile, data=True))
    with open(Path(filenames[0]), 'r') as fp:
        contents1 = read_cube(fp)
        atoms1 = contents1['atoms']
        grid1 = np.array(contents1['data'])
    
    with open(Path(filenames[1]), 'r') as fp:
        contents2 = read_cube(fp)
        atoms2 = contents2['atoms']
        grid2 = np.array(contents2['data'])
    
    assert atoms1 == atoms2, "Atoms objects should be the same."
    
    factor = 2
    print("factor is: ", factor , "NOT IMPLEMENTED YET!")
    
    if args.shift == 'auto':
        # TODO: perhaps can add argument to pick indices of the point?
        shift1_2 = grid1[0,0,0] - grid2[0,0,0]
    else:
        assert isinstance(args.shift, (float,int)), "ERROR: shift must be \'auto\' or otherwise float (or int) type!"
        shift1_2 = float(args.shift)
    
    grid2 = grid2 + shift1_2 # shift to make energies comparable if different normalization
    # General useful information
    print("Energy shift [eV] (added to grid 2): ", shift1_2)
    
    # Debug information
    #print(type(grid1),type(grid2))
    #print(type(grid1[0]),type(grid2[0]))
    #print("g1 num points < 5000: ", np.sum((grid1.flatten() < 5000.0)))
    #print("g2 num points < 5000: ", np.sum((grid2.flatten() < 5000.0)))
    #print("g1 num points > 400 and  < 5000: ", np.sum((400 < grid1.flatten()) * (grid1.flatten() < 5000.0)))
    #print("g2 num points > 400 and < 5000: ", np.sum((400 < grid2.flatten()) * (grid2.flatten() < 5000.0)))
    
    
    # Check the maximum number that are still very high
    #print("Grid1, max entry < 1e300: ", np.max(grid1[grid1 < 1e+300]))
    #print("Grid2, max entry < 1e300: ", np.max(grid2[grid2 < 1e+300]))
    #
    
    
    ### Filter out nans / infs and very large floats
    
    #assert np.all(np.isfinite(grid1)) and np.all(np.isfinite(grid2)), "Grids should be finite."
    
    inf_or_nan_indices = np.logical_or(np.logical_or(np.isinf(grid1), np.isnan(grid1)), np.logical_or(np.isinf(grid2), np.isnan(grid2)))
    
    not_inf_or_nan_indices = np.logical_not(inf_or_nan_indices)
    
    print(not_inf_or_nan_indices.shape, not_inf_or_nan_indices.dtype)
    
    not_overflow_indices = np.less(grid1, 1.0e+300) * np.less(grid2, 1.0e+300)
    
    inf_float = np.nan_to_num(np.inf)
    
    not_close_to_inf_indices = np.logical_not(np.logical_or(np.isclose(grid1,inf_float), np.isclose(grid2,inf_float)))
    
    nonzero_indices = np.logical_and((grid1 != 0.0), (grid2 != 0.0))
    
    no_problem_indices = np.logical_and(not_close_to_inf_indices, np.logical_and(not_inf_or_nan_indices, not_overflow_indices))
    
    no_problem_and_nonzero_indices = np.logical_and(no_problem_indices, nonzero_indices)
   
    #if args.emax is not None:
    # TO BE IMPLEMENTED



    ###
    
    # Is the indexing actally working?
    grid1_max = np.max(grid1[no_problem_indices])
    grid2_max = np.max(grid2[no_problem_indices])
    print("grid 1 filtered max: ", np.max(grid1[no_problem_indices]))
    print("grid 2 filtered max: ", np.max(grid2[no_problem_indices]))
    
    #
    
    abs_error = calc_abs_error(grid1[no_problem_indices], grid2[no_problem_indices])
    rel_error = calc_rel_error(grid1[no_problem_and_nonzero_indices], grid2[no_problem_and_nonzero_indices])
    tot_RMSD = calc_RMSD(grid1[no_problem_indices], grid2[no_problem_indices])
    
    # Check maximum errors ...
    print("Max pointswise abs error: ", np.max(abs_error))
    print("Max pointswise rel error: ", np.max(rel_error))
    #
    
    # Print genreal information - Could print to output file? Or redirect output in bash script
    print("Grids: ", filenames[0], "  -  ", filenames[1])
    print()
    print("Sum Total Abs. Error: ", np.sum(abs_error))
    print("Sum Total Rel. Error: ", np.sum(rel_error)) # Not much meaning
    print("Mean Abs. Error: ", np.mean(abs_error))
    print("Mean Rel. Error: ", np.mean(rel_error)) 
    print("Total RMSD: ", tot_RMSD)
    
    
    ### Print grids/files with errors 
    abs_error_cube_fname = "abs_error_" + Path(filenames[0]).stem + "_VS_" + Path(filenames[1]).stem + ".cube"
    rel_error_cube_fname = "rel_error_" + Path(filenames[0]).stem + "_VS_" + Path(filenames[1]).stem + ".cube"
    abs_error_xsf_fname = "abs_error_" + Path(filenames[0]).stem + "_VS_" + Path(filenames[1]).stem + ".xsf"
    rel_error_xsf_fname = "rel_error_" + Path(filenames[0]).stem + "_VS_" + Path(filenames[1]).stem + ".xsf"
    
    
    if cube_out:
        with open(abs_error_cube_fname, 'w') as fp:
            print("Writing asbolute error to cube-file: ", fp)
            write_cube(fp, atoms=atoms1, data=abs_error)
    
        with open(rel_error_cube_fname, 'w') as fp:
            print("Writing relative error to cube-file: ", fp)
            write_cube(fp, atoms=atoms1, data=rel_error)
    if xsf_out:
        with open(abs_error_xsf_fname, 'w') as fp:
            print("Writing absolute error to xsf-file: ", fp)
            write_xsf(fp, atoms=atoms1, data=abs_error)
    
        with open(rel_error_xsf_fname, 'w') as fp:
            print("Writing relative error to xsf-file: ", fp)
            write_xsf(fp, atoms=atoms1, data=rel_error)
    
    
    #####
    # Do error analysis in separate energy intervals/ranges
    # Find the total energy range, then split into equal size intervals
    # Could also do equal log-scale intervals
    num_nrg_intervals = args.numbins
    
    intervals = np.hstack((np.linspace(np.min([np.min(grid1), np.min(grid2)]), 
                                       np.max([grid1_max, grid2_max]), 
                                       num=num_nrg_intervals)[0:-1].reshape(-1,1), 
                           np.linspace(np.min([np.min(grid1), np.min(grid2)]), 
                                       np.max([grid1_max, grid2_max]), 
                                       num=num_nrg_intervals)[1:].reshape(-1,1)))
    
    
    # Choose lower threshhold, say log10 = 0 => log10(1)
    log_intervals_log = np.hstack((np.linspace(0, 
                                       np.log10(np.max([grid1_max, grid2_max])), 
                                       num=num_nrg_intervals)[0:-1].reshape(-1,1), 
                                   np.linspace(0, 
                                       np.log10(np.max([grid1_max, grid2_max])), 
                                       num=num_nrg_intervals)[1:].reshape(-1,1)))
    
    log_intervals = np.vstack(([0,1], np.power(10, log_intervals_log)))
    
    #print(log_intervals_log)
    #print(log_intervals)
    
    ### Compute the number of points in each interval
    # Linear intervals
    num_points_in_intervals_g1 = [np.count_nonzero(np.logical_and(np.greater(grid1, interval[0]), np.less_equal(grid1, interval[1]))) for interval in intervals]
    num_points_in_intervals_g2 = [np.count_nonzero(np.logical_and(np.greater(grid2, interval[0]), np.less_equal(grid2, interval[1]))) for interval in intervals]
    
    # Log intervals
    num_points_in_log_intervals_g1 = [np.count_nonzero(np.logical_and(np.greater(grid1, interval[0]), np.less_equal(grid1, interval[1]))) for interval in log_intervals]
    num_points_in_log_intervals_g2 = [np.count_nonzero(np.logical_and(np.greater(grid2, interval[0]), np.less_equal(grid2, interval[1]))) for interval in log_intervals]
    ###
    
    # Compute the errors in each energy interval separately
    # Linear intervals
    intervals_RMSD_g1 = [calc_RMSD(grid1[np.logical_and(np.greater(grid1, interval[0]), np.less_equal(grid1, interval[1]))],
        grid2[np.logical_and(np.greater(grid1, interval[0]), np.less_equal(grid1, interval[1]))]) for i,interval in enumerate(intervals) if (num_points_in_intervals_g1[i] != 0)]
    intervals_RMSD_g2 = [calc_RMSD(grid1[np.logical_and(np.greater(grid2, interval[0]), np.less_equal(grid2, interval[1]))],
        grid2[np.logical_and(np.greater(grid2, interval[0]), np.less_equal(grid2, interval[1]))]) for i,interval in enumerate(intervals) if (num_points_in_intervals_g2[i] != 0)]
    
    intervals_mean_abs_err_g1 = [np.mean(calc_abs_error(grid1[np.logical_and(np.greater(grid1, interval[0]), np.less_equal(grid1, interval[1]))],
        grid2[np.logical_and(np.greater(grid1, interval[0]), np.less_equal(grid1, interval[1]))])) for i,interval in enumerate(intervals) if (num_points_in_intervals_g1[i] != 0)]
    
    intervals_mean_abs_err_g2 = [np.mean(calc_abs_error(grid1[np.logical_and(np.greater(grid2, interval[0]), np.less_equal(grid2, interval[1]))],
        grid2[np.logical_and(np.greater(grid2, interval[0]), np.less_equal(grid2, interval[1]))])) for i,interval in enumerate(intervals) if (num_points_in_intervals_g2[i] != 0)]
    
    intervals_mean_rel_err_g1 = [np.mean(calc_rel_error(grid1[np.logical_and(np.greater(grid1, interval[0]), np.less_equal(grid1, interval[1]))],
        grid2[np.logical_and(np.greater(grid1, interval[0]), np.less_equal(grid1, interval[1]))])) for i,interval in enumerate(intervals) if (num_points_in_intervals_g1[i] != 0)]
    
    intervals_mean_rel_err_g2 = [np.mean(calc_rel_error(grid1[np.logical_and(np.greater(grid2, interval[0]), np.less_equal(grid2, interval[1]))],
        grid2[np.logical_and(np.greater(grid2, interval[0]), np.less_equal(grid2, interval[1]))])) for i,interval in enumerate(intervals) if (num_points_in_intervals_g2[i] != 0)]
    
    
    #log_intervals_RMSD_g1 = [calc_RMSD(grid1[np.logical_and(np.greater(grid1, interval[0]), np.less_equal(grid1, interval[1]))],
    #    grid2[np.logical_and(np.greater(grid1, interval[0]), np.less_equal(grid1, interval[1]))]) for interval in log_intervals]
    #log_intervals_RMSD_g2 = [calc_RMSD(grid1[np.logical_and(np.greater(grid2, interval[0]), np.less_equal(grid2, interval[1]))],
    #    grid2[np.logical_and(np.greater(grid2, interval[0]), np.less_equal(grid2, interval[1]))]) for interval in log_intervals]
    
    
    
    print()
    print("Linear intervals:")
    print("Energy intervals: ", intervals)
    print("Num gridpoints in intervals: ")
    print(num_points_in_intervals_g1)
    print(num_points_in_intervals_g2)
    print()
    print("Log intervals:")
    print("Num gridpoints in intervals: ")
    print(num_points_in_log_intervals_g1)
    print(num_points_in_log_intervals_g2)

###

if __name__ == "__main__":
    main()
