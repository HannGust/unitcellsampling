#!/usr/bin/env python

import numpy as np
import ase
import ase.units
from ase.io.cube import read_cube, write_cube
import argparse


def grid_energy_mask(ref_grid, cutoff):
    """Return a boolean mask that yields the grid points with energy below
    the cutoff energy."""
    return ref_grid <= cutoff


def energy_filter_grid(grid, ref_grid=None, cutoff=1.0, fill_value=None,
                       prec=np.float64):
    if ref_grid is None:
        ref_grid = grid
    assert grid.shape == refgrid.shape
    
    if fill_value is None:
        fill_value = np.nan_to_num(np.inf).as_type(prec)

    
    energy_mask = np.logical_not(grid_energy_mask(ref_grid, cutoff))
    
    filtered_grid = np.copy(grid)
    filtered_grid[energy_mask] = fill_value

    return filtered_grid


def define_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('grid', type=str, action='store',
            help='Specify the grid to filter/mask.')
    parser.add_argument('ref_grid', type=str, action='store',
            help='Specify the reference grid used for constructing the '
                  'filter/mask.')
    parser.add_argument('cutoff', type=float, action='store',
            help='The cutoff energy used for filtering. Every grid point in '
            'the reference grid with an energy above this value will be remove '
            'and replace by a fill value.')
    
    parser.add_argument('-f', type=float, action='store', default=None,
            help='The fill value to replace the filtered points with.')
    parser.add_argument('-o', type=str, action='store', default=None,
            help='Name of output file. Default no output file.')
    return parser

def main():
    p = define_parser()
    p.parse_args()

if __name__ == '__main__':
    main()
