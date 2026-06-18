import numpy as np
import ase
import os

import multiprocessing
import psutil

import scipy.constants as const

# from pymatgen import Structure
# import pymatgen.analysis.ewald # Depreciated warning!
import gridmc


class Walker:
    def __init__(self, coords, limits):
        self.coords = np.array(coords)
        self.limits = np.array(limits)

    def step(self, vec: np.array):
        assert (np.abs(vec) < self.limits).all()
        new_coords = self.coords + vec
        for idx, new_coord in enumerate(new_coords):
            if new_coord < 0:
                new_coords[idx] = self.limits[idx] + new_coord
            elif new_coord >= self.limits[idx]:
                new_coords[idx] = new_coord - self.limits[idx]
        assert (new_coords >= 0).all()
        assert (new_coords < self.limits).all()
        return new_coords

    def set_coords(self, coords):
        self.coords = coords


class MonteCarlo:
    def __init__(self, data: np.array, unitcell: np.array, n_atoms: int, perturbation_method='sstep', T=300, seed=None, ewald=False, charge=1.0):
        self.data = data
        self.n_atoms = n_atoms
        self.unitcell = unitcell
        self.charge = charge

        # Seed from OS
        if not seed:
            new_seed = int.from_bytes(os.urandom(3), byteorder='little')
            np.random.seed(new_seed)
        else:
            np.random.seed(seed)

        # Initialize atoms as walkers at random positions
        self.atom_walkers = [Walker(
            coords=np.random.randint(low=(0, 0, 0), high=data.shape),
            limits=data.shape)
            for _ in range(n_atoms)]

        self.visited = np.zeros(data.shape)
        for atom in self.atom_walkers:
            self.visited[atom.coords[0]][atom.coords[1]][atom.coords[2]] = 1

        self.perturbation_method = perturbation_method
        self.T = T
        self.ewald = ewald

        assert len(self.atom_walkers) == n_atoms

    def perturbate(self):
        if self.perturbation_method == 'sstep':
            positions = []
            for walker in self.atom_walkers:
                step = np.random.randint(low=-1, high=2, size=(3,))
                positions.append(walker.step(step))
            return positions
        else:
            raise NotImplementedError()

    def _mic(self, r1, r2, cell):
        # TODO Include r_c
        cell_lengths = [np.linalg.norm(c) for c in cell]
        vector = []
        for comp1, comp2, cl in zip(r1, r2, cell_lengths):
            dist = np.abs(comp1 - comp2)
            while not dist <= cl/2.:
                dist -= cl
            vector.append(dist)
        return vector

    def calculate_energy(self, coords):
        energy = 0

        for idx_1, atom_1 in enumerate(coords):
            energy += self.data[atom_1[0]][atom_1[1]][atom_1[2]]
            if not self.ewald:
                for idx_2, atom_2 in enumerate(coords):
                    if idx_1 > idx_2:
                        r_1 = np.matmul(self.unitcell, atom_1 /
                                        np.array(self.data.shape))
                        r_2 = np.matmul(self.unitcell, atom_2 /
                                        np.array(self.data.shape))
                        coulomb = self.coulomb(np.linalg.norm(
                            self._mic(r_1, r_2, self.unitcell)))
                        energy += coulomb
        if self.ewald:
            energy += gridmc.ewald_summation(
                self.unitcell,
                [(np.array(coord) / np.array(self.data.shape)).tolist()
                 for coord in coords],
                self.charge
            )

        return energy

    def coulomb(self, r):
        z1 = 1
        z2 = 1
        q1 = z1  # ASE units
        q2 = z2  # ASE units

        k = ase.units.m * ase.units.J/(4*np.pi*const.epsilon_0*ase.units.C**2)
        coulomb = k * q1 * q2 / r
        return coulomb

    def boltzmann(self, dE):
        if np.isclose(self.T, 0):
            return 0
        return np.exp(dE / (ase.units.kB * self.T))

    def run(self, n_steps=1000, log_skip=10):
        """
        Parameters
        ----------
        n_steps
            Specify number of steps to perform
        """
        energy = self.calculate_energy(
            [walker.coords for walker in self.atom_walkers])

        if log_skip:
            log_counter = 0
            energy_log = [energy]
            pos_log = [[walker.coords / walker.limits.astype(float)
                        for walker in self.atom_walkers]]

        for step in range(n_steps):
            pos_new = self.perturbate()

            energy_new = self.calculate_energy(pos_new)
            dE = energy - energy_new

            if dE >= 0 or self.boltzmann(dE) > np.random.rand():
                energy = energy_new
                for walker, p_new in zip(self.atom_walkers, pos_new):
                    walker.set_coords(p_new)
            for walker, p_new in zip(self.atom_walkers, pos_new):
                self.visited[walker.coords[0]
                             ][walker.coords[1]][walker.coords[2]] += 1
            if log_skip:
                log_counter += 1
                if log_counter >= log_skip:
                    energy_log.append(energy)
                    pos_log.append([walker.coords / walker.limits.astype(float)
                                    for walker in self.atom_walkers])
                    log_counter = 0
        return energy_log, pos_log, self.visited


def _run_mc(kwargs):
    mc = MonteCarlo(kwargs['data'],
                    kwargs['unit_cell'],
                    kwargs['n_atoms'],
                    kwargs['perturbation_method'],
                    kwargs['T'],
                    kwargs['seed'],
                    kwargs['ewald'],
                    kwargs['charge'])

    energies, positions, visited = mc.run(
        kwargs['mc_steps'], kwargs['log_skip'])
    return {
        'E': energies,
        'r': positions,
        'visited': visited
    }


def loading(data: np.array, unit_cell: np.array, n_atoms: int, n_walkers=1, mc_steps=10000, T=300, log_skip=10, perturbation_method='sstep', seed=None, ewald=False, charge=1.0):
    """Calcuate loading effects and provide new grid for certain loading.

    Parameters
    ----------
    data
        Specify energy isosurface shape(i, j, k) for a single atom
        List or np.array shape(i, j, k)
        Supply in eV
    unit_cell
        Supply in Angstrom
    n_atoms
        Specify number of atoms in unitcell
    Returns
    -------
    np.array
        shape(i, j, k) modified energy grid.
    """
    assert perturbation_method in ['sstep']
    assert n_atoms > 0
    assert n_walkers > 0
    assert mc_steps > 0
    assert T >= 0
    assert log_skip >= 0
    assert np.array(data).shape
    assert np.array(unit_cell).shape == (3, 3)

    kwargs = {
        'data': data,
        'unit_cell': unit_cell,
        'n_atoms': n_atoms,
        'n_walkers': n_walkers,
        'mc_steps': mc_steps,
        'T': T,
        'log_skip': log_skip,
        'perturbation_method': perturbation_method,
        'seed': seed,
        'ewald': ewald,
        'charge': charge
    }

    walker_positions = []
    walker_energies = []
    walker_visited = []

    n_cpus = psutil.cpu_count(logical=True)
    n_processes = n_walkers if n_cpus > n_walkers else n_cpus
    with multiprocessing.Pool(processes=n_processes) as pool:
        walker_dicts = [dic for dic in pool.imap_unordered(_run_mc,
                                                           [kwargs for _ in range(n_walkers)])]

    for i in range(n_walkers):
        walker_positions.append(walker_dicts[i]['r'])
        walker_energies.append(walker_dicts[i]['E'])
        walker_visited.append(walker_dicts[i]['visited'])

    return np.array(walker_energies), np.array(walker_positions), np.array(walker_visited)
