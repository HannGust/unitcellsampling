import ase
import numpy as np
import gemmi
from ase.spacegroup import get_spacegroup
import datetime
import time

from mendeleev import element
from sklearn.neighbors import KDTree
from unitcellsampling import energy_calculator


def get_cartesian_coords(fractional_coords, cell_vectors) -> np.ndarray:
    return np.dot(fractional_coords, cell_vectors)


def get_fractional_coords(cartesian_coords, cell_vectors) -> np.ndarray:
    return np.dot(cartesian_coords, np.linalg.inv(cell_vectors))


class UnitCellSampler:

    def __init__(self, atoms: ase.Atoms):
        self.atoms = atoms
        self.grid_vectors = None
        self.included_grid_vectors = None
        self.n_frac = None
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
                           exploit_symmetry=True):
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
            included_grid_points = np.full(grid_points.shape[0], True)
        else:
            included_grid_points = self.included_grid_vectors.flatten()
            grid_points = self.grid_vectors.reshape(
                np.product(self.grid_vectors.shape[:-1]),
                self.grid_vectors.shape[-1])

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
                # TODO: Fix this indexing routine, so that it works with supercells
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
        self._log_calculate_energies_after(
            grid_points, included_grid_points, n_exploited_symmetry)

        # Normalize
        energies = energies - energies.min()
        np.nan_to_num(energies, copy=False)

        return energies.reshape(included_grid_points.shape)

    def generate_grid_vectors(self, n_frac=(10, 10, 10), abs=None,
                              vdw_scale=0.75, charge_scale=None):
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
            3D vector or integer specifying the spacing (in Ångström)
            between grid points for the a, b and c direction in the unit
            cell respectively.

        vdw_scale
            Scaling factor Van der Waals radius

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

        vdw_radius = {}

        for atom in self.atoms:
            if atom.number not in vdw_radius:
                vdw_radius.update({
                    atom.number: element(
                        atom.symbol).vdw_radius * vdw_scale / 100.
                })  # [pm] / 100 = [Å];

        data_points = np.empty((len(a_mesh)*len(b_mesh)*len(c_mesh), 3))
        idx = 0
        for a in a_mesh:
            for b in b_mesh:
                for c in c_mesh:
                    vec = a + b + c
                    data_points[idx] = vec
                    idx += 1

        mesh_index = KDTree(data_points, metric='euclidean')
        ind = mesh_index.query_radius(
            self.atoms.positions,
            [vdw_radius[atom.number] for atom in self.atoms])

        ind = np.unique(np.hstack(ind))

        included = np.empty(data_points.shape[0], dtype=np.bool)

        idx1 = 0
        for idx2, point in enumerate(data_points):
            if len(ind) == idx1:
                included[idx2] = True
            elif ind[idx1] > idx2:
                included[idx2] = True
            elif ind[idx1] == idx2:
                included[idx2] = False
                idx1 += 1
            else:
                raise Exception('This should not have happened. ' +
                                'Something must have gone terribly wrong! ' +
                                'Check indexing.')

        assert included.shape[0] == data_points.shape[0]
        included = included.reshape(n_frac)
        data_points = data_points.reshape(n_frac + (data_points.shape[1],))

        self.grid_vectors = data_points
        self.included_grid_vectors = included

        return (data_points, included)

    def _log_calculate_energies_before(self, grid_points, included_grid_points, exploit_symmetry):
        print('Start calculation of energy grid...')
        print(str(datetime.datetime.now()))
        print()
        print('Calculating energies for the following structure (excluding the additional atom):')
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

        print('Gridpoint mesh to calculate', self.n_frac)
        print("Grid points to calculate:")

        for grid_point in grid_points:
            print(np.array2string(get_fractional_coords(
                grid_point, self.atoms.cell[:]),
                formatter={'float_kind': lambda x: "%.7f" % x}))

        print("Number of total grid points:", grid_points.shape[0])
        print("Number of neglected grid points due to vdw radii:",
              np.size(included_grid_points)
              - np.count_nonzero(included_grid_points))

    def _log_calculate_energies_after(self, grid_points,
                                      included_grid_points,
                                      n_exploited_symmetry):
        print("Number of saved calculations due to symmetry:",
              n_exploited_symmetry)

        print("Total number of actual calculations:",
              grid_points.shape[0] - n_exploited_symmetry -
              (np.size(included_grid_points)
               - np.count_nonzero(included_grid_points)))

        print("Total calculation time was {:.7f} seconds".format(
            time.time()-self.start))
