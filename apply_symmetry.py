""" Program that contains the symmetry routine and can apply it to
    existing grids to generate new grids. """

### OBSERVE: NOT FINISHED YET!!! ###

import ase
import numpy as np
import gemmi
import argparse

import unicellsampling.sample
from unicellsampling.sample import UnitCellSampler as UCS


# TODO: Finish and structure the parser - enable automatic and manual setting of spacegroup
parser = argparse.ArgumentParser(description="")
parser.add_argument('grid', action='store', help='The grid from which to read in grid data and atoms object. Can be cube-file.')
parser.add_argument('-o', '--out', action='store', default='sym_grid_out.cube', help='The name of the output file name, to name the cube file for the resulting grid.')

# TODO:Think about whether this is a good approach?
def add_symmetry_arguments(parser:argparse.ArgumentParser):
    """ Add the symmetry settings to an argument parser. """
    parser.add_argument('-s', '--sym', action='store_true', help='Specifies that symmetry should be applied to the grid (symmetrization).')
    return parser


# TODO: Fix this piece of code (the raw symmetry code from the UCS sampler)
if exploit_symmetry:
    bool_grid = gemmi.Int8Grid(
        self.n_frac[0], self.n_frac[1], self.n_frac[2])
    energies_grid = gemmi.FloatGrid(
        self.n_frac[0], self.n_frac[1], self.n_frac[2])

    if self.spacegroup:
        bool_grid.spacegroup = self.spacegroup
        energies_grid.spacegroup = self.spacegroup
    else:
        bool_grid.spacegroup = gemmi.find_spacegroup_by_number(
            get_spacegroup(self.atoms))
        energies_grid.spacegroup = gemmi.find_spacegroup_by_number(
            get_spacegroup(self.atoms))

    # For debugging symmetry / H 
    gemmi_spgrp_from_structure = gemmi.find_spacegroup_by_number(
            get_spacegroup(self.atoms))
    print("UCS: Gemmi spacegroup from sturcture: ", gemmi_spgrp_from_structure)
    print("UCS: Gemmi spacegroup internal: ", self.spacegroup)
    #i

    if self.n_supercell is None:
        self.n_supercell = (1, 1, 1)

    # gemmi unitcell TESTING! WARNING: EXPERIMENTAL CODE!!!
    if self.gemmi_unitcell:
          gemmi_uc = gemmi.UnitCell(*self.atoms.cell.cellpar())
          bool_grid.set_unit_cell(gemmi_uc)
          energies_grid.set_unit_cell(gemmi_uc)
    #

grid_points = np.array(grid_points)

# Logging before main loop
self._log_calculate_energies_before(
    grid_points, included_grid_points, exploit_symmetry)

n_exploited_symmetry = 0
energies = np.empty(grid_points.shape[0], dtype=np.float64)
for idx, grid_point in enumerate(grid_points):
    if not included_grid_points[idx]:
        energies[idx] = np.inf
        #print(np.array2string(get_fractional_coords(
        #      grid_point, self.atoms.cell[:]),
        #    formatter={'float_kind': lambda x: "%.7f" % x}),
        #    "Excluded point")
        continue

    if exploit_symmetry:
