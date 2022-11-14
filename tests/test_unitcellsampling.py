#!/usr/bin/env python

"""Tests for `unitcellsampling` package."""


import unittest
from pathlib import Path

import ase
import ase.io
import ase.io.cube
import ase.visualize
import time

import numpy as np

import unitcellsampling.sample as ucs


class TestUnitcellsampling(unittest.TestCase):
    """Tests for `unitcellsampling` package."""

    def setUp(self):
        """Set up test fixtures, if any."""
        pass

    def tearDown(self):
        """Tear down test fixtures, if any."""
        pass

    def test_sample(self):
        """Test something."""
        indir = Path(Path(__file__).resolve().parent, 'structures')
        outdir = Path(Path(__file__).resolve().parent, 'out')

        au_cube = ase.Atoms(symbols='Au',
                            scaled_positions=[[0.5, 0.5, 1.0]],
                            cell=[
                                [30.0, 0.0, 0.0],
                                [0.0, 30.0, 0.0],
                                [0.0, 0.0, 30.0]
                            ],
                            pbc=True)
        start = time.time()
        end = time.time()
        print("Elapsed time", end-start)
        sampler = ucs.UnitCellSampler(au_cube)
        sample = sampler.generate_grid_vectors(n_frac=20)

        assert (np.isclose(
            sample[0], sampler.generate_grid_vectors(n_frac=20)[0])).all()
        end = time.time()
        print("Elapsed time", end-start)

        sample, included = sample
        sample = sample.reshape(np.product(sample.shape[:-1]), 3)
        included = included.flatten()
        sample = [vec for idx, vec in enumerate(sample) if included[idx]]

        end = time.time()
        print("Elapsed time", end-start)
        [au_cube.append(ase.Atom('H', s)) for s in sample]
        lgps = ase.io.read(str(Path(indir, 'lgps.cif')))
        ase.io.write(str(Path(outdir, 'lgps_test_no_sample.xyz')), lgps)
        with open(Path(outdir, 'test.cube'), 'w') as fp:
            ase.io.cube.write_cube(fp, lgps, data=None)
        sampler = ucs.UnitCellSampler(lgps)
        sample = sampler.generate_grid_vectors(abs=0.8)

        sample, included = sample
        sample = sample.reshape(np.product(sample.shape[:-1]), 3)
        included = included.flatten()
        sample_incl = [vec for idx, vec in enumerate(sample) if included[idx]]
        sample_del = [vec for idx, vec in enumerate(
            sample) if not included[idx]]

        [lgps.append(ase.Atom('H', s)) for s in sample_incl]
        ase.io.write(Path(outdir, 'lgps_test.xyz'), lgps)

        end = time.time()
        print("Elapsed time", end-start)
        lgps = ase.io.read(str(Path(indir, 'lgps.cif')))
        [lgps.append(ase.Atom('H', s)) for s in sample_del]
        ase.io.write(Path(outdir, 'lgps_test_del.xyz'), lgps)
        end = time.time()
        print("Elapsed time", end-start)
