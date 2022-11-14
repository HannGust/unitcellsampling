# import unittest

# import time
# import numpy as np


# import pymatgen.analysis.ewald
# import pymatgen.io.cif as cif
# import unitcellsampling.ewald_sum
# import gridmc
# import ase
# import ase.io

# from pymatgen import Structure


# class TestUnitcellsampling(unittest.TestCase):
#     """Tests for `unitcellsampling` package."""

#     def setUp(self):
#         """Set up test fixtures, if any."""
#         pass

#     def tearDown(self):
#         """Tear down test fixtures, if any."""
#         pass

#     # def test_cubic_convert(self):
#     #     L = 10.0
#     #     wire = ase.Atoms('Au',
#     #                      positions=[[0, L / 2, L / 2],
#     #                                 [L/8, L/3, L/2], [L/4, L/7, L/4]],
#     #                      cell=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
#     #                      pbc=[1, 1, 1])
#     #     ase.gui.view(wire)

#     def test_ewald(self):
#         """Test something."""
#         # # Test 1
#         coords = [[0, 0, 0], [0.0, 0.0, 0.1]]
#         cell = [[10., 0, 0], [0, 10., 0], [0, 0, 10.]]
#         struct = Structure(cell, ["Li", "Li"], coords)
#         struct.add_oxidation_state_by_element({'Li': 1.})
#         # self_energy_pymatgen = unitcellsampling.ewald_sum.EwaldSummation(
#         #     struct, real_space_cut=14.5, eta=2).point_energy
#         # real_pymatgen = unitcellsampling.ewald_sum.EwaldSummation(
#         #     struct, real_space_cut=20.0, eta=4, acc_factor=12).real_space_energy
#         # rec_pymatgen = unitcellsampling.ewald_sum.EwaldSummation(
#         #     struct, recip_space_cut=2., eta=4, acc_factor=12).reciprocal_space_energy
#         total_pymatgen = unitcellsampling.ewald_sum.EwaldSummation(
#             struct).total_energy

#         # assert np.isclose(gridmc.ewald_summation_self(2, 1., 2.), self_energy_pymatgen)
#         # assert np.isclose(gridmc.ewald_summation_real(cell, coords, 1, 20., 2.), real_pymatgen)
#         # assert np.isclose(gridmc.ewald_summation_reciprocal(cell, coords, 2., 1., struct.volume, 4), rec_pymatgen)

#         assert np.isclose(gridmc.ewald_summation(
#             cell, coords, 1), total_pymatgen)

#         # # Test 2
#         coords = [[0, 0, 0], [0.0, 0.0, 0.2]]
#         cell = [[10., 0, 0], [0, 10., 0], [0, 0, 10.]]
#         struct = Structure(cell, ["Li"]*len(coords), coords)
#         struct.add_oxidation_state_by_element({'Li': 1.})
#         # self_energy_pymatgen = unitcellsampling.ewald_sum.EwaldSummation(
#         #     struct, real_space_cut=14.5, eta=2).point_energy
#         # real_pymatgen = unitcellsampling.ewald_sum.EwaldSummation(
#         #     struct, real_space_cut=20.0, eta=4, acc_factor=12).real_space_energy
#         # rec_pymatgen = unitcellsampling.ewald_sum.EwaldSummation(
#         #     struct, recip_space_cut=2., eta=4, acc_factor=12).reciprocal_space_energy
#         total_pymatgen = unitcellsampling.ewald_sum.EwaldSummation(
#             struct).total_energy

#         # assert np.isclose(gridmc.ewald_summation_self(
#         #     len(coords), 1., 2.), self_energy_pymatgen)
#         # assert np.isclose(gridmc.ewald_summation_real(
#         #     cell, coords, 1, 20., 2.), real_pymatgen)
#         # assert np.isclose(gridmc.ewald_summation_reciprocal(
#         #     cell, coords, 2., 1., struct.volume, 4), rec_pymatgen)

#         assert np.isclose(gridmc.ewald_summation(
#             cell, coords, 1), total_pymatgen)

#         # Test 3
#         coords = [[0, 0, 2], [0.0, 0.0, 0.1], [
#             0.5, 0.6, 0.15], [0.9, 0.3, -0.85]]
#         cell = [[-10., 0, 9], [3.4, 10., 9], [-10, 10, -0.8]]
#         struct = Structure(cell, ["Li"]*len(coords), coords)
#         struct.add_oxidation_state_by_element({'Li': 1.})
#         # self_energy_pymatgen = unitcellsampling.ewald_sum.EwaldSummation(
#         #     struct, real_space_cut=14.5, eta=2).point_energy
#         # real_pymatgen = unitcellsampling.ewald_sum.EwaldSummation(
#         #     struct, real_space_cut=20.0, eta=4, acc_factor=12).real_space_energy
#         # rec_pymatgen = unitcellsampling.ewald_sum.EwaldSummation(
#         #     struct, recip_space_cut=.6, eta=4, acc_factor=12).reciprocal_space_energy

#         ewald_pymatgen = pymatgen.analysis.ewald.EwaldSummation(
#             struct).total_energy
#         ewald_rust = gridmc.ewald_summation(cell, coords, 1)
#         # print("RUST COMP", ewald_pymatgen, ewald_rust)
#         # TOTAL EWALD

#         # assert(np.isclose(ewald_pymatgen, ewald_rust))

#         # start = time.time()
#         # for _ in range(1000):
#         #     coords = [[0, 0, 0], [0.0, 0.0, 0.5]]
#         #     struct = Structure([[10., 0, 0], [0, 10., 0], [0, 0, 10.]], [
#         #                        "Li", "Li"], coords)
#         #     struct.add_oxidation_state_by_element({'Li': 1.})
#         #     energy_pymatgen = pymatgen.analysis.ewald.EwaldSummation(
#         #         struct).total_energy
#         # end = time.time()
#         # print("Time, Ewald (pymatgen)", end - start)

#         # start = time.time()
#         # for _ in range(1000):
#         #     coords = [[0, 0, 0], [0.0, 0.0, 0.5]]
#         #     cell = [[10., 0, 0], [0, 10., 0], [0, 0, 10.]]
#         #     gridmc.ewald_summation(cell, coords)

#         # end = time.time()
#         # print("Time, Ewald (rust)", end - start)

#         coords = [[0.000000, 0.000000, 0.000000],
#                   [0.000000, 0.000000, 1.089000],
#                   [1.026719, 0.000000, -0.363000],
#                   [-0.513360, -0.889165, -0.363000],
#                   [-0.513360, 0.889165, -0.363000]]

#         struct = Structure([[11., 0, -2.], [0, 10., 0],
#                             [0, 0, 7.4]], ["Li"]*5, coords)
#         struct.add_oxidation_state_by_element({'Li': 1.})

#         energy_pymatgen = pymatgen.analysis.ewald.EwaldSummation(
#             struct).total_energy
#         energy_unitcellsampling = unitcellsampling.ewald_sum.EwaldSummation(
#             struct).total_energy

#         assert(np.isclose(energy_pymatgen, energy_unitcellsampling))
