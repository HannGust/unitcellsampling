import unittest
from ase.io.cube import read_cube, write_cube
import numpy as np

from unitcellsampling.supercell import create_supercell_grid


class TestUnitcellsampling(unittest.TestCase):
    """Tests for `unitcellsampling` package."""

    def test_supercell_grid(self):
        with open("tests/structures/tdc.cube") as fp:
            cube = read_cube(fp)

        super_cell_grid = create_supercell_grid(
            cube['atoms'].cell, cube['data'], dim=(2, 2, 2))
        assert((np.array(super_cell_grid.shape)
               == np.array(cube['data'].shape)*2).all())
