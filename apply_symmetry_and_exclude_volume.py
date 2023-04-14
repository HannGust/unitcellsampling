""" Program that contains the symmetry and volume exclusion routines
    and can apply it to existing grids to generate new grids. """

### OBSERVE: NOT FINISHED YET!!! ###

import ase
import numpy as np
import gemmi
import argparse

import unicellsampling.sample
from unicellsampling.sample import UnitCellSampler as UCS


parser = argparse.ArgumentParser(description="")
parser.add_argument('grid', action='store', help='The grid from which to read in grid data and atoms object. Can be cube-file.')
parser.add_argument('-o', '--out', action='store', default='sym_vdw_grid_out.cube', help='The name of the output file name, to name the cube fiel fo the resulting grid.')
parser.add_argument('-s', '--sym', action='store_true', help='Specifies that symmetry should be applied to the grid (symmetrization).')


#TODO: Add the identical symmetry routine as in apply_symmetry.py here!
#      Add the identical volume exclusion routine from exclude_volume.py here!  
#      Import whatever functions used in these! 
#      Modify the argument parser to easily be able to enter and modify settings as desired. 
#      Use maybe functions in the other routines to setup the parser and add the relevant arguments?
#      Make unit tests for the rountines separately and test to test them jointly.

