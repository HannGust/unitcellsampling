""" Program that contains the symmetry routine and can apply it to
    existing grids to generate new grids. """

### OBSERVE: NOT FINISHED YET!!! ###

import ase
import numpy as np
import gemmi
import argparse

import unitcellsampling.sample
from unitcellsampling.sample import UnitCellSampler as UCS
from unitcellsampling.sample import get_fractional_coords


# TODO: Finish and structure the parser - enable automatic and manual setting of spacegroup
def init_parser():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('grid', action='store', help='The grid from which to read in grid data and atoms object. Can be cube-file.')
    parser.add_argument('-o', '--out', action='store', default='sym_grid_out.cube', help='The name of the output file name, to name the cube file for the resulting grid.')
    return parser

# TODO:Think about whether this is a good approach?
def add_symmetry_arguments(parser:argparse.ArgumentParser):
    """ Add the symmetry settings to an argument parser. """
    parser.add_argument('-s', '--sym', action='store_true', help='Specifies that symmetry should be applied to the grid (symmetrization).')
    parser.add_argument('-n', '--nfrac', action='store', nargs='+', type=int, default=None, help='Specifies the grid size, if this is not to be decided from the grid itself.')
    parser.add_argument('--sg', '--spacegroup', action='store', default=None, help='Specifies the spacegroup that should be used for symmetrization, if it is not to be decided from the input structure.')
    return parser

###  Scrap code to be reused:

#    assert (atoms is not None) or (spacegroup is not None), "One of atoms or spacegroup must be given!"

# Setting of spacegroup according to how it is done in unitcellsampler
# get_spacegroup is ase's function for this purpose
#        bool_grid.spacegroup = gemmi.find_spacegroup_by_number(
#            get_spacegroup(self.atoms))
#        energies_grid.spacegroup = gemmi.find_spacegroup_by_number(
#            get_spacegroup(self.atoms))

###

### Scrap code probably not to be reused
#    # For debugging symmetry / H 
#    gemmi_spgrp_from_structure = gemmi.find_spacegroup_by_number(
#            get_spacegroup(self.atoms))
#    print("UCS: Gemmi spacegroup from sturcture: ", gemmi_spgrp_from_structure)
#    print("UCS: Gemmi spacegroup internal: ", self.spacegroup)
#    #i
## Logging before main loop
#self._log_calculate_energies_before(
#    grid_points, included_grid_points, exploit_symmetry)v
#            atoms = self.atoms + ase.Atom(atom, grid_point)
#            if callable(method):
#                energies[idx] = method(atoms)
#            elif method == 'dft':
#                energies[idx] = energy_calculator.cp2k_dft(atoms)
#            elif method == 'dftb':
#                raise NotImplementedError
#            elif method == 'xtb':
#                raise NotImplementedError
#                energies[idx] = energy_calculator.cp2k_xtb(atoms)
#            elif method == 'lammps_simple':
#                energies[idx] = energy_calculator.lammps_forcefield_simple(
#                    atoms)
#            elif method == 'lammps':
#                energies[idx] = energy_calculator.lammps_forcefield(atoms)
#            elif method == 'qe':
#                raise NotImplementedError
#            elif method == 'random':
#                energies[idx] = np.random.rand()
#            else:
#                raise ValueError('Invalid value', method, 'for method')
#            # Disabled printing here / H
#            #print(np.array2string(get_fractional_coords(
#            #    grid_point, self.atoms.cell[:]),
#            #    formatter={'float_kind': lambda x: "%.7f" % x}),
#            #    'Calculated energy to',
#            #    energies[idx])
#      if not included_grid_points[idx]:
#            energies[idx] = np.inf
#            #print(np.array2string(get_fractional_coords(
#            #      grid_point, self.atoms.cell[:]),
#            #    formatter={'float_kind': lambda x: "%.7f" % x}),
#            #    "Excluded point")
#            continue
#

###

# This should be in principle the symmetry routine from UCS
def UCS_apply_spgrp_symmetry_to_grid(grid, grid_points, atoms, spacegroup=None, n_frac=None):
    """Applies symmetry to a grid within a unitcell, in principle
    identically as if it was sampled with UnitCellSampler with
    utilization of spacegroup symmetry. Returns the symmetrized
    grid. """

    grid = np.array(grid)
    grid_points = np.array(grid_points)

    grid_energies = grid.flatten()
    grid_points = grid_points.reshape(-1,3)
    
    if n_frac is None:
        n_frac = grid.shape
    elif isinstance(n_frac, int):
        n_frac = (n_frac,) * 3

    assert isinstance(n_frac, (tuple, list)), "n_frac must be tuple or list."
    assert len(n_frac) == 3, "n_frac must have length 3."
    assert all(isinstance(i, int) for i in n_frac)

    bool_grid = gemmi.Int8Grid(n_frac[0], n_frac[1], n_frac[2])
    energies_grid = gemmi.FloatGrid(n_frac[0], n_frac[1], n_frac[2])

    if spacegroup:
        bool_grid.spacegroup = spacegroup
        energies_grid.spacegroup = spacegroup
    else:
        bool_grid.spacegroup = gemmi.find_spacegroup_by_number(
            get_spacegroup(atoms))
        energies_grid.spacegroup = gemmi.find_spacegroup_by_number(
            get_spacegroup(atoms))

    
    n_exploited_symmetry = 0
    energies = np.empty(grid_points.shape[0], dtype=np.float64)

    for idx, grid_point in enumerate(grid_points): 
        grid_point_index = get_fractional_coords(
               grid_point, atoms.cell[:])
        grid_point_index[0] = grid_point_index[0] * n_frac[0]
        grid_point_index[1] = grid_point_index[1] * n_frac[1]
        grid_point_index[2] = grid_point_index[2] * n_frac[2]
        grid_point_index = np.array(
               np.around(grid_point_index), dtype=np.int)

        if bool_grid.get_value(*grid_point_index):
            energies[idx] = energies_grid.get_value(*grid_point_index)
            # Disabled printing here / H
            #print(np.array2string(get_fractional_coords(
            #   grid_point, self.atoms.cell[:]),
            #   formatter={'float_kind': lambda x: "%.7f" % x}),
            #   "Exploited_symmetry:",
            #   energies[idx])
            n_exploited_symmetry += 1
            continue

        energies[idx] = grid_energies[idx]

        bool_grid.set_value(*grid_point_index, 1)
        bool_grid.symmetrize_max()

        energies_grid.set_value(*grid_point_index, energies[idx])
        if energies[idx] <= 0:
            energies_grid.symmetrize_min()
        else:
            energies_grid.symmetrize_max()

        # Logging after main loop
        #self._log_calculate_energies_after(
        #    grid_points, included_grid_points, n_exploited_symmetry)

        # Normalize
        #energies = energies - energies.min()
        #np.nan_to_num(energies, copy=False)

        return energies.reshape(n_frac)

def main():
    #parser = argparse.ArgumentParser(description="")
    #parser.add_argument('grid', action='store', help='The grid from which to read in grid data and atoms object. Can be cube-file.')
    #parser.add_argument('-o', '--out', action='store', default='sym_grid_out.cube', help='The name of the output file name, to name the cube file for the resulting grid.')
    parser = init_parser()
    
    #add_symmetry_arguments(parser)

    args = parser.parse_args()

    with open(args.grid, 'r') as f:
        cube_content = read_cube(f, read_data=True)
        atoms = cube_content['atoms']
        grid = cube_content['data']

    grid = np.array(grid)
    n_frac = grid.shape

    # generate grid coordinates here:
    ucs = UCS(atoms)
    grid_points, inlcuded = ucs.generate_grid_vectors(n_frac=n_frac, vdw_scale=0.0, midvox=args.midvox)
    #

    symmetrized_grid = UCS_apply_spgrp_symmetry_to_grid(grid, grid_points, atoms, spacegroup=None, n_frac=None)

    if args.out:
        with open(args.out, 'w') as f:
            write_cube(f, data=symmetrized_grid)

    return

if __name__ == "__main__":
    main()


#####################################################################
#       This is the raw symmetry code from the UCS sampler          #
#####################################################################
#
#if exploit_symmetry:
#    bool_grid = gemmi.Int8Grid(
#        self.n_frac[0], self.n_frac[1], self.n_frac[2])
#    energies_grid = gemmi.FloatGrid(
#        self.n_frac[0], self.n_frac[1], self.n_frac[2])
#
#    if self.spacegroup:
#        bool_grid.spacegroup = self.spacegroup
#        energies_grid.spacegroup = self.spacegroup
#    else:
#        bool_grid.spacegroup = gemmi.find_spacegroup_by_number(
#            get_spacegroup(self.atoms))
#        energies_grid.spacegroup = gemmi.find_spacegroup_by_number(
#            get_spacegroup(self.atoms))
#
#    # For debugging symmetry / H 
#    gemmi_spgrp_from_structure = gemmi.find_spacegroup_by_number(
#            get_spacegroup(self.atoms))
#    print("UCS: Gemmi spacegroup from sturcture: ", gemmi_spgrp_from_structure)
#    print("UCS: Gemmi spacegroup internal: ", self.spacegroup)
#    #i
#
#    if self.n_supercell is None:
#        self.n_supercell = (1, 1, 1)
#
#    # gemmi unitcell TESTING! WARNING: EXPERIMENTAL CODE!!!
#    if self.gemmi_unitcell:
#          gemmi_uc = gemmi.UnitCell(*self.atoms.cell.cellpar())
#          bool_grid.set_unit_cell(gemmi_uc)
#          energies_grid.set_unit_cell(gemmi_uc)
#    #
#
#grid_points = np.array(grid_points)
#
## Logging before main loop
#self._log_calculate_energies_before(
#    grid_points, included_grid_points, exploit_symmetry)
#
#n_exploited_symmetry = 0
#energies = np.empty(grid_points.shape[0], dtype=np.float64)
#for idx, grid_point in enumerate(grid_points):
#    if not included_grid_points[idx]:
#        energies[idx] = np.inf
#        #print(np.array2string(get_fractional_coords(
#        #      grid_point, self.atoms.cell[:]),
#        #    formatter={'float_kind': lambda x: "%.7f" % x}),
#        #    "Excluded point")
#        continue
#
#    if exploit_symmetry:
#                grid_point_index = get_fractional_coords(
#                    grid_point, self.atoms.cell[:])
#                grid_point_index[0] = grid_point_index[0]*self.n_frac[0]*self.n_supercell[0]
#                grid_point_index[1] = grid_point_index[1]*self.n_frac[1]*self.n_supercell[1]
#                grid_point_index[2] = grid_point_index[2]*self.n_frac[2]*self.n_supercell[2]
#                grid_point_index = np.array(
#                    np.around(grid_point_index), dtype=np.int)
#
#                if bool_grid.get_value(*grid_point_index):
#                    energies[idx] = energies_grid.get_value(*grid_point_index)
#                    # Disabled printing here / H
#                    #print(np.array2string(get_fractional_coords(
#                     #   grid_point, self.atoms.cell[:]),
#                     #   formatter={'float_kind': lambda x: "%.7f" % x}),
#                     #   "Exploited_symmetry:",
#                     #   energies[idx])
#                    n_exploited_symmetry += 1
#                    continue
#
#            atoms = self.atoms + ase.Atom(atom, grid_point)
#            if callable(method):
#                energies[idx] = method(atoms)
#            elif method == 'dft':
#                energies[idx] = energy_calculator.cp2k_dft(atoms)
#            elif method == 'dftb':
#                raise NotImplementedError
#            elif method == 'xtb':
#                raise NotImplementedError
#                energies[idx] = energy_calculator.cp2k_xtb(atoms)
#            elif method == 'lammps_simple':
#                energies[idx] = energy_calculator.lammps_forcefield_simple(
#                    atoms)
#            elif method == 'lammps':
#                energies[idx] = energy_calculator.lammps_forcefield(atoms)
#            elif method == 'qe':
#                raise NotImplementedError
#            elif method == 'random':
#                energies[idx] = np.random.rand()
#            else:
#                raise ValueError('Invalid value', method, 'for method')
#            # Disabled printing here / H
#            #print(np.array2string(get_fractional_coords(
#            #    grid_point, self.atoms.cell[:]),
#            #    formatter={'float_kind': lambda x: "%.7f" % x}),
#            #    'Calculated energy to',
#            #    energies[idx])
#
#            if exploit_symmetry:
#                bool_grid.set_value(*grid_point_index, 1)
#                bool_grid.symmetrize_max()
#                energies_grid.set_value(*grid_point_index, energies[idx])
#                if energies[idx] <= 0:
#                    energies_grid.symmetrize_min()
#                else:
#                    energies_grid.symmetrize_max()
#
#        # Logging after main loop
#        self._log_calculate_energies_after(
#            grid_points, included_grid_points, n_exploited_symmetry)
#
#        # Normalize
#        energies = energies - energies.min()
#        np.nan_to_num(energies, copy=False)
#
#        return energies.reshape(included_grid_points.shape)
#
####################################################################################
