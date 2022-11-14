#!/usr/bin/env python

"""Tests for `unitcellsampling` package."""

import ase.io
import ase.io.xsf
from ase.spacegroup import get_spacegroup
import unitcellsampling.cli
import os
import unittest

from pathlib import Path
from dotenv import load_dotenv
load_dotenv()


class TestCalculate(unittest.TestCase):
    """Tests for `unitcellsampling` package."""

    def setUp(self):
        """Set up test fixtures, if any."""
        pass

    def tearDown(self):
        """Tear down test fixtures, if any."""
        pass

    def test_exploit_grid_symmetry(self):
        inpath = Path('tests/structures/lgps_dft.cube')
        with open(inpath, 'r') as fp:
            lgps = ase.io.cube.read_cube(fp)

        unitcellsampling.cli.exploit_grid_symmetry(
            lgps['data'], get_spacegroup(lgps['atoms']))

    def test_cli(self):
        inpath = Path('tests/structures/lgps.cif')
        outpath = Path('tests/out/lgps.cube')
        os.system('python unitcellsampling/cli.py --abs 0.5 --out_file '
                  + str(outpath) + ' ' + str(inpath))
