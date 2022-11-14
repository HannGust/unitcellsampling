#!/usr/bin/env python

"""Tests for `unitcellsampling` package."""

import os
import unittest

import ase.io
import numpy as np

import unitcellsampling
import unitcellsampling.sample

import shutil
from pathlib import Path
from dotenv import load_dotenv
load_dotenv()


class TestCalculate(unittest.TestCase):
    """Tests for `unitcellsampling` package."""

    def setUp(self):
        """Set up test fixtures, if any."""
        shutil.rmtree(os.getenv('UCS_CALCULATION_DIR'), ignore_errors=True)

    def tearDown(self):
        """Tear down test fixtures, if any."""
        shutil.rmtree(os.getenv('UCS_CALCULATION_DIR'), ignore_errors=True)

    def test_energy(self):
        """Test potential energy calculation"""
        indir = Path(Path(__file__).resolve().parent, 'structures')
        outdir = Path(Path(__file__).resolve().parent, 'out')

        os.environ['UCS_CALCULATION_DIR'] = str(Path(outdir,
                                                     'ucs_calculation').resolve())
        os.environ['UCS_WORK_DIR'] = str(outdir)

        lgps = ase.io.read(Path(indir, 'lgps.cif'))

        sampler = unitcellsampling.sample.UnitCellSampler(lgps)
        sample = np.array(sampler.generate_grid_vectors(n_frac=3)[0][0][0])

        energies = sampler.calculate_energies(
            grid_points=sample,
            method=unitcellsampling.energy_calculator.lammps_forcefield_simple, atom="Li")

        assert energies.shape == (3,)

        print(energies.shape)
        print(energies)

    def test_lammps_lib(self):
        indir = Path(Path(__file__).resolve().parent, 'structures')
        outdir = Path(Path(__file__).resolve().parent, 'out')

        os.environ['UCS_CALCULATION_DIR'] = str(Path(outdir,
                                                     'ucs_calculation').resolve())
        os.environ['UCS_WORK_DIR'] = str(outdir)
        lgps = ase.io.read(Path(indir, 'lgps.cif'))
        sampler = unitcellsampling.sample.UnitCellSampler(lgps)
        sample = np.array(sampler.generate_grid_vectors(n_frac=4)[0][0][0])
        # energies = sampler.calculate_energies(
        #     grid_points=sample,
        #     method=ucs.energy_calculator.lammps_forcefield, atom="Li")
        # print(energies)

    def test_lammps_simple(self):
        indir = Path(Path(__file__).resolve().parent, 'structures')
        lgps = ase.io.read(Path(indir, 'lgps.cif'))
        sampler = unitcellsampling.sample.UnitCellSampler(lgps)
        sample = np.array(sampler.generate_grid_vectors(n_frac=4)[0][0][0])
        energies = sampler.calculate_energies(
            grid_points=sample,
            method=unitcellsampling.energy_calculator.lammps_forcefield_simple, atom="Li")

        assert energies.shape == (4,)

    def test_lammps_simple_inbuilt(self):
        indir = Path(Path(__file__).resolve().parent, 'structures')
        lgps = ase.io.read(Path(indir, 'lgps.cif'))
        sampler = unitcellsampling.sample.UnitCellSampler(lgps)
        sample = np.array(sampler.generate_grid_vectors(n_frac=4)[0][0][0])
        energies = sampler.calculate_energies(
            grid_points=sample,
            method='lammps_simple', atom="Li")
        energies_costum = sampler.calculate_energies(
            grid_points=sample,
            method=unitcellsampling.energy_calculator.lammps_forcefield_simple,
            atom="Li")
        assert np.isclose(energies, energies_costum).all()
