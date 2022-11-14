# import ase.units
# import ase.io.cube
# from pathlib import Path
# import unittest
# import numpy as np
# import functools
# import time
# import gridmc

# from unitcellsampling import isosurface


# class TestCalculate(unittest.TestCase):
#     """Tests for `unitcellsampling` package."""

#     def setUp(self):
#         """Set up test fixtures, if any."""
#         pass
#         # shutil.rmtree(os.getenv('UCS_CALCULATION_DIR'))

#     def tearDown(self):
#         """Tear down test fixtures, if any."""
#         pass

#     def test_walker(self):
#         walker = isosurface.Walker(limits=(79, 52, 85), coords=(0, 0, 0))
#         assert (walker.step([-1, -1, -1]) == np.array([78, 51, 84])).all()
#         assert (walker.step([0, 0, 0]) == np.array([0, 0, 0])).all()
#         assert (walker.step([1, 1, 1]) == np.array([1, 1, 1])).all()
#         assert (walker.step([0, 1, 1]) == np.array([0, 1, 1])).all()
#         assert (walker.step([-1, 1, 1]) == np.array([78, 1, 1])).all()
#         assert (walker.step([-1, -1, 1]) == np.array([78, 51, 1])).all()

#     def test_loading(self):
#         tdc_path = Path('tests/structures/tdc.cube')
#         with open(tdc_path, 'r') as fp:
#             tdc = ase.io.cube.read_cube(fp)

#         loading = 2  # Two atoms
#         start = time.time()
#         walker_energies, walker_scaled_pos, walker_hist = isosurface.loading(
#             tdc['data']*ase.units.kJ / ase.units.mol,
#             tdc['atoms'].get_cell()[:],
#             n_atoms=loading,
#             n_walkers=1,
#             T=1000,
#             mc_steps=1000,
#             log_skip=20,
#             ewald=True)
#         end = time.time()
#         print("Time", end - start)
#         print(walker_hist[0].sum())
#         hist = functools.reduce(lambda h1, h2: h1+h2, walker_hist)

#         with open('tdc_hist.cube', 'w') as fp:
#             ase.io.cube.write_cube(fp, tdc['atoms'], data=hist)
#         atoms = tdc['atoms']

#         # Only observe single walker from now on:
#         with open('tdc_hist_single_walker.cube', 'w') as fp:
#             ase.io.cube.write_cube(fp, tdc['atoms'], data=walker_hist[0])

#         energies = walker_energies[0]
#         scaled_pos = walker_scaled_pos[0]
#         for _ in range(loading):
#             atoms.append(ase.Atom('Li'))

#         traj = ase.io.Trajectory('li.traj', 'w')
#         # traj.write(atoms)
#         for E, pos in zip(energies, scaled_pos):
#             for idx_atom in range(loading):
#                 atoms[-idx_atom -
#                       1].position = np.matmul(atoms.get_cell(), pos[idx_atom])
#             traj.write(atoms, energy=E.sum())

#         assert tdc['data'].any()
#         assert tdc_path.is_file()

#     def test_compare_mmc(self):
#         tdc_path = Path('tests/structures/tdc.cube')
#         with open(tdc_path, 'r') as fp:
#             tdc = ase.io.cube.read_cube(fp)

#         T = 1000
#         loading = 2  # Two atoms
#         mmc_steps = 1000
#         n_walkers = 16
#         point_charge = 1

#         start = time.time()
#         walker_energies, walker_scaled_pos, walker_hist = isosurface.loading(
#             tdc['data']*ase.units.kJ / ase.units.mol,
#             tdc['atoms'].get_cell()[:],
#             n_atoms=loading,
#             n_walkers=n_walkers,
#             T=T,
#             mc_steps=mmc_steps,
#             charge=point_charge,
#             # log_skip=20,
#             ewald=True)
#         end = time.time()
#         print("Time (python)", end - start)

#         # fn correct_energy_grid_metropolis_monte_carlo(
#         # single_particle_grid: Vec<f64>,
#         # grid_shape: Vec<usize>,
#         # unitcell: Vec<Vec<f64>>,
#         # temperature: f64,
#         # point_charge: f64,
#         # n_walkers: u8,
#         # n_runs: u32,
#         # mmc_steps: u64,
#         # ) -> PyResult<Vec<f64>> {
#         start = time.time()

#     # fn correct_energy_grid_metropolis_monte_carlo(
#     #     single_particle_grid: Vec<f64>,
#     #     grid_shape: Vec<usize>,
#     #     unitcell: Vec<Vec<f64>>,
#     #     temperature: f64,
#     #     point_charge: f64,
#     #     n_walkers: u8,
#     #     n_runs: u32,
#     #     mmc_steps: u64,
#     #     equilibrate: u64,
#     #     verbosity_level: u8,
#     #     initial_grid_positions: Vec<Vec<i32>>,
#     # ) -> PyResult<Vec<f64>> {
#         grid = gridmc.correct_energy_grid_metropolis_monte_carlo(
#             (tdc['data']*ase.units.kJ / ase.units.mol).flatten(),
#             tdc['data'].shape,
#             tdc['atoms'].get_cell()[:],
#             T,
#             point_charge,
#             loading,
#             n_walkers,
#             mmc_steps,
#             0,
#             1,
#             [[]]
#         )
#         end = time.time()
#         print("Time (rust)", end - start)

#         def hist2E(hist, dE=0.02):
#             prop = hist/hist.sum()
#             prop = np.array([[[c if c > 1e-9 else 1e-9 for c in b]
#                               for b in a] for a in prop])
#             E_mc = -np.log(prop)*ase.units.kB*T
#             E_mc = E_mc - E_mc.min()
#             return E_mc
#         l1 = functools.reduce(lambda h1, h2: h1+h2, walker_hist)

#         with open('python.cube', 'w') as fp:
#             ase.io.cube.write_cube(fp, tdc['atoms'], data=hist2E(l1))
#         with open('rust.cube', 'w') as fp:
#             ase.io.cube.write_cube(fp, tdc['atoms'], data=np.array(
#                 grid).reshape(tdc['data'].shape))

#         # assert False
