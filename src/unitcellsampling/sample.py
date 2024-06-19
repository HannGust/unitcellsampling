import ase
import numpy as np
import gemmi
from ase.spacegroup import get_spacegroup
import datetime
import time

from unitcellsampling import energy_calculator
from unitcellsampling.volume_exclusion import RadialExcluder, ScaledVdWExcluder


def get_cartesian_coords(fractional_coords, cell_vectors) -> np.ndarray:
    return np.dot(fractional_coords, cell_vectors)


def get_fractional_coords(cartesian_coords, cell_vectors) -> np.ndarray:
    return np.dot(cartesian_coords, np.linalg.inv(cell_vectors))


class UnitCellSampler:

    def __init__(self, atoms: ase.Atoms):
        # basics
        self.atoms = atoms
        self.grid_vectors = None
        self.n_frac = None

        # volume exclusion related
        self.included_grid_vectors = None
        self.cutoff_included = None
        self.vdw_included = None
        self.radial_cutoff_excluder = None
        self.scaled_vdw_excluder = None
        
        # symmetry related
        self.spacegroup = None # / H

        self.gemmi_unitcell = None # IN EXPERIMENTAL TESTING STAGE!!! /H
        # To maybe be introduced:
        self.supercell = None # / H
        self.unitcell = None # / H
        self.n_supercell = None # / H
        self.n_sim_cell = None # / H

    def set_spacegroup(self, sp_number: int):
        self.spacegroup = gemmi.find_spacegroup_by_number(sp_number)


    def calculate_energies(self, grid_points=None,
                           fractional=False,
                           method='random',
                           atom='Li',
                           exploit_symmetry=True,
                           normalize=True):
        """Calculate energy for `atom` on each grid point within a unitcell

        Parameters
        ----------
        grid_points: (m, 3) array_like
            grid points to place `atom` on. If not
            specified you have to call `generate_grid_vectors` first

        fractional: bool
            specifies whether `grid_points` are specified via fractional
            coordinates.

        atom: str
            specifies the atom to be placed on `grid_points`

        method: str or func
            specifies how energies should be computed.
            If function:
                Function should have the following signature:
                func(atoms: ase.Atoms) -> np.float64
            If string:
                'dft': CP2K quickstep calculation
                'dftb': DFTB method. not implemented yet
                'xtb': xTB method. not implemented yet
                'qe': Quantum Espresso: not implemented yet
                'random': fill grid points with random values.

        exploit_symmetry: bool
            use spacegroup symmetry to save computational effort

        Returns
        -------
        np.array
            Real space coordinates of grid points in unit cell.
        """
        self.start = time.time()
        if grid_points is not None:
            # if provided with external grid points
            included_grid_points = np.full(grid_points.shape[0], True)
            
            # These are set for proper logging
            included_radial_cutoff = None
            included_vdw = None

            print("UnitCellSampler: External grid points given - No volume exclusion nor information about volume exclusion is accessed or applied by the sampler.")

        else:
            # Otherwise: Use internal grid points
            included_grid_points = self.included_grid_vectors.flatten()
            grid_points = self.grid_vectors.reshape(
                np.product(self.grid_vectors.shape[:-1]),
                self.grid_vectors.shape[-1])
            
            # These are set for proper logging 
            included_radial_cutoff = self.cutoff_included
            included_vdw = self.vdw_included
            
            print("UnitCellSampler: Internal grid points are used - Accessing information about volume exclusion masks from sampler.")

        if fractional:
            raise Warning('Fractional input not tested. Check results!')
            grid_points = get_cartesian_coords(grid_points, self.atom.cell[:])

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
        self._log_calculate_energies_before(grid_points,
                                            included_grid_points,
                                            exploit_symmetry,
                                            included_radial_cutoff,
                                            included_vdw)

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
                grid_point_index = get_fractional_coords(
                    grid_point, self.atoms.cell[:])
                grid_point_index[0] = grid_point_index[0]*self.n_frac[0]*self.n_supercell[0]
                grid_point_index[1] = grid_point_index[1]*self.n_frac[1]*self.n_supercell[1]
                grid_point_index[2] = grid_point_index[2]*self.n_frac[2]*self.n_supercell[2]
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

            atoms = self.atoms + ase.Atom(atom, grid_point)
            if callable(method):
                energies[idx] = method(atoms)
            elif method == 'dft':
                energies[idx] = energy_calculator.cp2k_dft(atoms)
            elif method == 'dftb':
                raise NotImplementedError
            elif method == 'xtb':
                raise NotImplementedError
                energies[idx] = energy_calculator.cp2k_xtb(atoms)
            elif method == 'lammps_simple':
                energies[idx] = energy_calculator.lammps_forcefield_simple(
                    atoms)
            elif method == 'lammps':
                energies[idx] = energy_calculator.lammps_forcefield(atoms)
            elif method == 'qe':
                raise NotImplementedError
            elif method == 'random':
                energies[idx] = np.random.rand()
            else:
                raise ValueError('Invalid value', method, 'for method')
            # Disabled printing here / H
            #print(np.array2string(get_fractional_coords(
            #    grid_point, self.atoms.cell[:]),
            #    formatter={'float_kind': lambda x: "%.7f" % x}),
            #    'Calculated energy to',
            #    energies[idx])

            if exploit_symmetry:
                bool_grid.set_value(*grid_point_index, 1)
                bool_grid.symmetrize_max()
                energies_grid.set_value(*grid_point_index, energies[idx])
                if energies[idx] <= 0:
                    energies_grid.symmetrize_min()
                else:
                    energies_grid.symmetrize_max()

        # Logging after main loop
        self._log_calculate_energies_after(grid_points,
                                included_grid_points,
                                n_exploited_symmetry,
                                included_radial_cutoff=included_radial_cutoff,
                                included_vdw=included_vdw
                                )

        # Normalize
        if exploit_symmetry:
            print("Symmetry used: Casting to np.float32 before grid normalization.")
            energies = energies.astype(np.float32)
        
        if normalize:
            min_grid_energy = energies.min()
            print("Normalizing grid by shifting minimum energy to 0.0, i.e. subtracting minimum energy.")
            print("Minimum grid energy: ", min_grid_energy)
            print("Added shift: ", -min_grid_energy)
            energies = energies - min_grid_energy
        else:
            print("normalize =", normalize)
            print("Final normalization of grid disabled: minimum energy is NOT shifted to 0.0.")
        np.nan_to_num(energies, copy=False)

        return energies.reshape(included_grid_points.shape)


    def generate_grid_vectors(self, n_frac=(10, 10, 10), abs=None, 
                              cutoff_radii=0.0, 
                              vdw_scale=None, 
                              midvox=False, 
                              charge_scale=None):
        """Get gridpoint vectors in the unitcell for the use of the
        TuTraSt methodology.

        Parameters
        ----------
        atoms
            Atoms object containing the atoms and unitcell.

        n_frac
            3D vector or integer specifying the number of grid points
            for the a, b and c direction in the unit cell respectively.

        abs
            3D vector or float specifying the spacing (in Ångström)
            between grid points for the a, b and c direction in the unit
            cell respectively.

        cutoff_radii
            Single radius or dictionary mapping from atomic symbols to radii
            in Ångström, for use in spherical cutoffs. That is, exclusion of
            points within a radius of the framework atoms. If a single float
            is given, this will be applied to all atoms in the framework.
            If a dictionary is given, it must contain keys corresponding to
            all atom types in the framework, and/or a \"default\" entry. 
            The default entry, if present, will be used for all atom types
            not explicitly specified.
            

        vdw_scale
            Scaling factor or dictionary mapping from atomic symbols to 
            scaling factors, to use for scaling van der Waals radii. This
            is used for a spherical cutoff based on scaled van der Waals
            radii, e.g. exclusion of points within a radius of the framework
            atoms. If a single float is given, this will be applied as a
            van der Waals scaling factor for all atoms in the framework.
            If a dictionary is given, it must contain keys corresponding to
            all atom types in the framework, and/or a \"default\" entry. 
            The default entry, if present, will be used for all atom types
            not explicitly specified.

        midvox
            Setting that specifies that grid points should be the 
            midpoint of voxels, rather than the corners. Corresponds to
            shifting the grid point coordinates with the vector 
            0.5 * (a/nx, b/ny, c/nz), i.e. half of the corresponding 
            grid spacing in each direction.

        charge_scale
            If set vdw_radius will be weighted by the atomic charge via
            this factor (not implemented).

        Returns
        -------
        np.array
            Real space coordinates of grid points in unit cell.
        """

        if charge_scale:
            raise NotImplementedError

        assert abs or n_frac and not abs and n_frac

        if n_frac:
            if isinstance(n_frac, int):
                n_frac = (n_frac,)*3
            assert len(n_frac) == 3
            assert all(type(x) == int for x in n_frac)

        if abs:
            if isinstance(abs, (int, float)):
                abs = (abs,)*3
            n_frac = tuple(
                [int(np.linalg.norm(vec)/abs[idx])
                 for idx, vec in enumerate(self.atoms.get_cell()[:])]
            )

        self.n_frac = n_frac

        a_mesh = np.linspace(
            0, self.atoms.cell[:][0], n_frac[0], endpoint=False)
        b_mesh = np.linspace(
            0, self.atoms.cell[:][1], n_frac[1], endpoint=False)
        c_mesh = np.linspace(
            0, self.atoms.cell[:][2], n_frac[2], endpoint=False)

        # New midvox functionality - TODO: Test it properly.
        if midvox:
            a_mesh += 0.5 * (self.atoms.cell[:][0]/n_frac[0])
            b_mesh += 0.5 * (self.atoms.cell[:][1]/n_frac[1])
            c_mesh += 0.5 * (self.atoms.cell[:][2]/n_frac[2])


        cart_coords = np.empty((len(a_mesh)*len(b_mesh)*len(c_mesh), 3))
        idx = 0
        for a in a_mesh:
            for b in b_mesh:
                for c in c_mesh:
                    vec = a + b + c
                    cart_coords[idx] = vec
                    idx += 1
        
        #### Start block of new vdw code
        #### Here also new cutoff radii implemented
        
        if cutoff_radii is not None:
            # NOTE: Excluder return within_cutfoff, outside_cutoff masks
            cutoff_excluder = RadialExcluder(radii=cutoff_radii)
            cutoff_excluded, cutoff_included = cutoff_excluder.construct_cutoff_filters(self.atoms, 
                                                                                        grid_coords=cart_coords, 
                                                                                        frac_input=False, 
                                                                                        periodic=True)
        else:
            # If not set, set excluder to None and include all points in this mask:
            cutoff_excluder = None
            cutoff_included = np.full(cart_coords.shape[0], fill_value=True, dtype=bool)
            

        ### Here's the vdw exclusion:
        if vdw_scale is not None:
            vdw_excluder = ScaledVdWExcluder(vdw_scaling_map=vdw_scale)
            vdw_excluded, vdw_included = vdw_excluder.construct_cutoff_filters(self.atoms, 
                                                                           grid_coords=cart_coords, 
                                                                           frac_input=False, 
                                                                           periodic=True)
        else:
            # If not set, set excluder to None and include all points in this mask:
            vdw_excluder = None
            vdw_included = np.full(cart_coords.shape[0], fill_value=True, dtype=bool)


        # Combine the masks like so:
        included = np.logical_and(cutoff_included, vdw_included)


        ##### End block of new vdw code

        assert included.shape[0] == cart_coords.shape[0]
        included = included.reshape(n_frac)
        cart_coords = cart_coords.reshape(n_frac + (cart_coords.shape[1],))

        self.grid_vectors = cart_coords
        self.included_grid_vectors = included

        # Adding these filters as well too the sampler 
        self.cutoff_included = cutoff_included.reshape(n_frac)
        self.vdw_included = vdw_included.reshape(n_frac)

        # NOTE: Adding Excluders too
        self.radial_cutoff_excluder = cutoff_excluder
        self.scaled_vdw_excluder = vdw_excluder

        return (cart_coords, included)


    def _log_calculate_energies_before(self, grid_points,
                                    included_grid_points,
                                    exploit_symmetry,
                                    included_radial_cutoff,
                                    included_vdw):
        print('Start calculation of energy grid...')
        print(str(datetime.datetime.now()))
        print()
        print('Calculating energies for the following framework structure (excluding the additional atom):')
        print()
        print('Unit cell:')
        print(np.array2string(self.atoms.get_cell()[:],
                              formatter={'float_kind': lambda x: "%.7f" % x}))
        print()
        print('Scaled coordinates:')
        chemical_symbols = self.atoms.get_chemical_symbols()
        scaled_positions = self.atoms.get_scaled_positions()
        for (symbol, position) in zip(chemical_symbols, scaled_positions):
            print(symbol, np.array2string(position,
                                          formatter={'float_kind': lambda x: "%.7f" % x}))
        print()            


        if exploit_symmetry:
            #print("Using spacegroup", get_spacegroup(self.atoms))
            print("Using spacegroup", self.spacegroup) # / H
            print()
        
        if self.radial_cutoff_excluder is not None and included_radial_cutoff is not None:
            print("Radial cutoff exclusion used:")
            self.radial_cutoff_excluder.print_settings()
            print()

        if self.scaled_vdw_excluder is not None and included_vdw is not None:
            print("Scaled van der Waals exclusion used:")
            self.scaled_vdw_excluder.print_settings()
            print()

        print("Gridpoint mesh to calculate:", self.n_frac)
        print("Total number of grid points:", grid_points.shape[0])
        print("Number of neglected grid points from combined spherical exclusion:",
              np.size(included_grid_points)
              - np.count_nonzero(included_grid_points))
        if included_radial_cutoff is not None:
            print("Number of neglected grid points from cutoff radii:",
                  np.size(included_radial_cutoff) - np.count_nonzero(included_radial_cutoff))         
        if included_vdw is not None:
            print("Number of neglected grid points from scaled van der Waals radii:",
                  np.size(included_vdw) - np.count_nonzero(included_vdw))
        print()
        print("Grid points to calculate:")
        for grid_point in grid_points:
            print(np.array2string(get_fractional_coords(
                grid_point, self.atoms.cell[:]),
                formatter={'float_kind': lambda x: "%.7f" % x}))
        
        print()
        print("Sampling started...")


    def _log_calculate_energies_after(self, grid_points,
                                      included_grid_points,
                                      n_exploited_symmetry,
                                      included_radial_cutoff=None,
                                      included_vdw=None):
        print("Sampling completed.")

        print("=============================================")
        print("==          GRID SAMPLING SUMMARY          ==")
        print("=============================================")
        if self.n_frac:
            print("Total grid shape:", self.n_frac)
        print("Number of total grid points:", grid_points.shape[0])
        actual_num_calcs = (grid_points.shape[0] - n_exploited_symmetry -
              (np.size(included_grid_points)
               - np.count_nonzero(included_grid_points)))
        print("Total number of actual calculations:", actual_num_calcs)
        print("Number of neglected grid points from combined spherical exclusion:",
              np.size(included_grid_points)
              - np.count_nonzero(included_grid_points))
        
        if included_radial_cutoff is not None:
            print("Number of neglected grid points from cutoff radii:",
                  np.size(included_radial_cutoff) - np.count_nonzero(included_radial_cutoff))
            
        if included_vdw is not None:
            print("Number of neglected grid points from scaled van der Waals radii:",
                  np.size(included_vdw) - np.count_nonzero(included_vdw))
            

        print("Number of saved calculations due to symmetry:",
              n_exploited_symmetry)

        total_time = time.time()-self.start

        print("Total calculation time was {:.7f} seconds".format(total_time))
        print("Average time/grid point calculation: {:.7f} s".format(
                                         total_time/actual_num_calcs))
        print("Extrapolated est. time for full grid: {:.7f} s".format(
                   grid_points.shape[0] * total_time/actual_num_calcs))
        print("Estimated time saved: {:.7f} s".format(
            (grid_points.shape[0] * total_time/actual_num_calcs) - total_time))
        print("=============================================")
        print()
