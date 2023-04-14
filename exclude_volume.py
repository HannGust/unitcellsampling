""" Program that contains the volume exclusion routine
    and can apply it to existing grids to generate new grids. """

### OBSERVE: NOT FINISHED YET!!! ###

import ase
import numpy as np
import gemmi
import argparse

import unicellsampling.sample
from unicellsampling.sample import UnitCellSampler as UCS

# TODO: Fix the parser, add relevant arguments
parser = argparse.ArgumentParser(description="")
parser.add_argument('grid', action='store', help='The grid from which to read in grid data and atoms object. Can be cube-file.')
parser.add_argument('-o', '--out', action='store', default='vol_excl_grid_out.cube', help='The name of the output file name, to name the cube file for the resulting grid to which volume exclusion has been applied.')

def add_volume_exclusion_arguments(arg_parser:argparse.ArgumentParser):
    arg_parser.add_argument("--type")
    arg_parser.add_argument("-r")
    arg_parser.add_argument("")

    return arg_parser

def exclude_volume(coords:np.ndarray, grid:np.ndarray, excl_type='vdW', rad=None):
    if not isinstance(coords, np.ndarray):
        coords = np.array(coords)

    if not isinstance(grid, np.ndarray):
        grid = np.array(grid)

    assert coords.shape[:-1] == grid.shape[:], "Coordinate array must have compatible shape with energy grid array!"
    # alternative: assert coords.shape == grid.shape + (,3), "Coordinate array must have compatible shape with energy grid array!"


