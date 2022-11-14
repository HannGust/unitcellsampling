"""CLI for unitcellsampling."""

import argparse
import json
import os
import psutil
import re
import subprocess
from pathlib import Path

import ase
import ase.io
import ase.io.cube
import ase.units
from ase.calculators.cp2k import CP2K
from ase.io.trajectory import Trajectory, TrajectoryWriter
from ase.optimize import BFGS
from ase.spacegroup import Spacegroup, get_spacegroup
from ase.visualize import view
import numpy as np

import gridmc
#from ppinterpf import ppinterpf # Cannot find this file, is not python package. Must be own program by Benjamin Bolbrinker // Hannes
import unitcellsampling.sample as ucs
from unitcellsampling.get_partial_ion_charges import parse_ddec_file
from unitcellsampling.get_partial_ion_charges import parse_repeat_file


__author__ = """Benjamin Bolbrinker"""
__email__ = 'benjamin.bolbrinker@kemi.uu.se'
__version__ = '0.1.0'

"""Console script for unitcellsampling."""


def create_traj_from(filename, traj_path='traj.traj', lim=(10000, 20000), view_traj=False):
    traj = TrajectoryWriter(traj_path)
    coords = []
    energies = []
    with open(filename, 'r') as fp:
        while "Printing trajectory:" not in fp.readline():
            pass
        fp.readline()
        fp.readline()
        regex_coord = r',\s+(\d+)\s+(\d+)\s+(\d+)'
        regex_energy = r'\d+,\s([-+]?(?:[0-9]*[.])?[0-9]+(?:[eE][-+]?\d+)?),'
        for line in fp:
            if "End of trajectory" in line:
                break
            coords.append(np.array(
                re.findall(regex_coord, line), dtype='int').tolist())
            energies.append(float(re.match(regex_energy, line).group(1)))

    coords = np.array(coords)

    # TODO: Make exact read from file
    ####
    cell = [max(coords[:, :, 0].flatten()) + 1,
            max(coords[:, :, 1].flatten()) + 1,
            max(coords[:, :, 2].flatten()) + 1]
    ###

    for e, coord in zip(energies[lim[0]:lim[1]], coords[lim[0]:lim[1]]):
        traj.write(
            ase.Atoms('Li' * len(coords[0]), positions=coord, cell=cell), energy=e)
    if view_traj:
        view(Trajectory(traj_path))
    return traj


def view_trajectory():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_file',
                        help='File containing the unitcell to be sampled.')
    parser.add_argument('--start', type=int, default=None,
                        help='Lower limit.')
    parser.add_argument('--end', type=int, default=None,
                        help='Upper limit.')
    args = parser.parse_args()
    create_traj_from(args.input_file, Path(args.input_file).suffix +
                     '.traj', view_traj=True, lim=(args.start, args.end))


def down_sample():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_file',
                        help='File containing the unitcell to be sampled.')
    parser.add_argument('--spacing', type=float, default=0.2,
                        help='Desired spacing.')
    args = parser.parse_args()

    in_path = Path(args.input_file)
    if in_path.suffix == '.xsf':
        atoms = ase.io.read(in_path, read_data=False)
        data = ase.io.read(in_path, read_data=True)
    elif in_path.suffix == '.cube':
        with open(in_path, 'r') as fp:
            cube_data = ase.io.cube.read_cube(fp)
            atoms = cube_data['atoms']
            data = cube_data['data']
    else:
        raise Exception('Unknown input file type.')

    n_frac_a = int(np.linalg.norm(atoms.get_cell()[:][0]) / args.spacing)
    n_frac_b = int(np.linalg.norm(atoms.get_cell()[:][1]) / args.spacing)
    n_frac_c = int(np.linalg.norm(atoms.get_cell()[:][2]) / args.spacing)

#    data_i = ppinterpf.interp_periodic(data, nfrac=(    # Commenting out, as it won't work; program ppinterpf is not accesible
#        n_frac_a, n_frac_b, n_frac_c), method='linear')
    data_i = None # Filler value
    print('WARNING: This function has been altered and does not function properly as the original.')
    print('Expect results to be useless. This has been done due to some functionalities not being available, ')
    print('upon which the original version (by Ben Bolbrinker) relied on. // Hannes Gustafsson')

    out_path = in_path.stem + '_down_sampled_' + str(args.spacing) + '.cube'
    with open(out_path, 'w') as fp:
        ase.io.cube.write_cube(fp, atoms, data_i)


def exploit_grid_symmetry(grid, spacegroup, average='mean'):
    sp = Spacegroup(spacegroup)

    for idx_a, a in enumerate(grid):
        for idx_b, b in enumerate(a):
            for idx_c, c in enumerate(b):
                site = [idx_a/grid.shape[0],
                        idx_b / grid.shape[1], idx_c/grid.shape[2]]
                equivalent_sites = sp.equivalent_sites([site])[0]
                equivalent_sites_idx = np.rint(np.array([
                    equivalent_sites[:, 0]*grid.shape[0],
                    equivalent_sites[:, 1]*grid.shape[1],
                    equivalent_sites[:, 2]*grid.shape[2]
                ]).T).astype(int)

                grid_mean = 0
                for site_idx in equivalent_sites_idx:
                    grid_mean += grid[site_idx[0], site_idx[1], site_idx[2]]
                grid_mean /= equivalent_sites_idx.shape[0]

                for site_idx in equivalent_sites_idx:
                    grid[site_idx[0], site_idx[1], site_idx[2]] = grid_mean
    return grid


def average_grid_values_by_symmetry():
    parser = argparse.ArgumentParser(
        description='Average grid by applying crystal symmetry.')
    parser.add_argument(
        'input_file', help='File containing the single particle energy grid.')
    parser.add_argument('--temperature', type=float,
                        help='Temperature which was used in MMC in Kelvin (default: 1000).', default=1000)

    args = parser.parse_args()
    in_path = Path(args.input_file)
    if in_path.suffix == '.xsf':
        atoms = ase.io.read(in_path, read_data=False)
        data = ase.io.read(in_path, read_data=True)
    elif in_path.suffix == '.cube':
        with open(in_path, 'r') as fp:
            cube_data = ase.io.cube.read_cube(fp)
            atoms = cube_data['atoms']
            data = cube_data['data']
    else:
        raise Exception('Unknown input file type.')
    data = np.exp(-data / (ase.units.kB * args.temperature))

    atoms_spacegroup = get_spacegroup(atoms, symprec=1e-05)
    print(atoms_spacegroup)
    data_new = exploit_grid_symmetry(data, atoms_spacegroup)
    out_path = in_path.stem + '_averaged.cube'

    data_new = -ase.units.kB*args.temperature * np.log(data_new)
    with open(out_path, 'w') as fp:
        ase.io.cube.write_cube(fp, atoms, data_new)


def create_sparse_grid_from(grid, spacing=0.2):
    pass


def run_mmc():
    parser = argparse.ArgumentParser(
        description='Correct energy grid via metropolis monte carlo.')
    parser.add_argument(
        'input_file', help='File containing the single particle energy grid.')
    parser.add_argument('--unit', help='Energy unit in input_file.',
                        type=str, choices=['ev', 'ry', 'kjmol'], default='ev')
    parser.add_argument('--temperature', type=float,
                        help='Temperature to be used in MMC in Kelvin (default: 1000).', default=1000)
    parser.add_argument('--point_charge', type=float,
                        help='Point charge of ions (default: 1.0).', default=1.0)
    parser.add_argument(
        '--n_runs', type=int, help='Number of independent runs. For maximal efficiency a multiple of the usable threads is recommended (default: 64).', default=64)
    parser.add_argument('--mmc_steps', type=int,
                        help='Number Markov permutations in MMC loop (default: 50000).', default=50000)
    parser.add_argument('--equilibrate', type=int,
                        help='Number Markov permutations in MMC loop before keeping track of visited grid points (default: 5000).', default=5000)
    parser.add_argument('--start_from_file', action='store_true',
                        help='Start MMC from Li coordinates in file.')
    parser.add_argument('--n_walkers', type=int,
                        help='Loading (default: 1).', default=1)
    parser.add_argument('--verbose', action='store_true',
                        help='Adjust print level.')
    parser.add_argument('--continue_calculation', action='store_true',
                        help='Continue or restart previous calculation.')
    parser.add_argument('--read_walker_coords_from', type=str, default=None,
                        help='Start MMC calculation from specified walker coordinates. Overwrites n_walkers. If omitted place walkers randomly.')
    parser.add_argument('--read_walker_coords_from_minimum', type=str, default=None,
                        help='Start MMC calculation from walker coordinates with lowest energy. Overwrites n_walkers. If omitted place walkers randomly.')

    args = parser.parse_args()
    in_path = Path(args.input_file)
    if in_path.suffix == '.xsf':
        atoms = ase.io.read(in_path, read_data=False)
        data = ase.io.read(in_path, read_data=True)
    elif in_path.suffix == '.cube':
        with open(in_path, 'r') as fp:
            cube_data = ase.io.cube.read_cube(fp)
            atoms = cube_data['atoms']
            data = cube_data['data']
    else:
        raise Exception('Unknown input file type.')

    if args.unit == 'ev':
        data = data*ase.units.eV
    elif args.unit == 'ry':
        data = data*ase.units.Ry
    elif args.unit == 'kjmol':
        data = data*ase.units.kJ/ase.units.mol
    else:
        raise ValueError

    final_coords_from_last = []
    n_walkers = args.n_walkers

    # TODO: Add search for global minimum
    if args.read_walker_coords_from:
        with open(args.read_walker_coords_from, 'r') as fp:
            while "Printing trajectory:" not in fp.readline():
                pass

            last_line = None
            # TODO: Delete ':' in rust log
            for line in fp:
                if "End of trajectory" not in line:
                    last_line = line
                else:
                    break

            print(last_line)
            regex = r',\s+(\d+)\s+(\d+)\s+(\d+)'
            final_coords_from_last = np.array(
                re.findall(regex, last_line), dtype='int').tolist()
            n_walkers = len(final_coords_from_last)

    # TODO: Add search for global minimum
    if args.read_walker_coords_from_minimum:
        with open(args.read_walker_coords_from_minimum, 'r') as fp:
            while "Printing trajectory:" not in fp.readline():
                pass
            fp.readline()
            fp.readline()
            minimum_e = None
            last_line = None
            # TODO: Delete ':' in rust log
            for line in fp:
                if "End of trajectory" not in line:
                    regex_f_exp = r'\d+,\s?((?:[-+]?\d*\.?\d+)(?:[eE](?:[-+]?\d+))?),'
                    energy = float(re.match(regex_f_exp, line).group(1))
                    if minimum_e is None:
                        minimum_e = energy
                        last_line = line
                    elif minimum_e >= energy:
                        minimum_e = energy
                        last_line = line
                else:
                    break
            print(last_line)
            regex = r',\s+(\d+)\s+(\d+)\s+(\d+)'
            final_coords_from_last = np.array(
                re.findall(regex, last_line), dtype='int').tolist()
            n_walkers = len(final_coords_from_last)

    data_corrected = np.array(
        gridmc.correct_energy_grid_metropolis_monte_carlo(
            data.flatten(),
            data.shape, atoms.get_cell()[:],
            temperature=args.temperature,
            point_charge=args.point_charge,
            n_walkers=n_walkers,
            n_runs=args.n_runs,
            mmc_steps=args.mmc_steps,
            equilibrate=args.equilibrate,
            verbose=args.verbose,
            continue_calculation=args.continue_calculation,
            initial_grid_positions=final_coords_from_last))

    out_path = in_path.stem + '_corrected_' + \
        str(n_walkers) + '_temperature_' + \
        str(int(args.temperature)) + '.cube'
    with open(out_path, 'w') as fp:
        ase.io.cube.write_cube(fp, atoms, data_corrected.reshape(data.shape))


def main():
    """Console script for unitcellsampling."""
    parser = argparse.ArgumentParser()
    parser.add_argument('in_file',
                        help='File containing the unitcell to be sampled.')
    parser.add_argument('--abs', default=2, type=float,
                        help='Specify spacing of gridpoints in Ångström.')
    parser.add_argument('--nfrac', type=int,
                        help='Number of gridpoints along each axis.')
    parser.add_argument('--method', type=str, default='lammps',
                        help='Method for energy calculation on each gridpoint.')
    parser.add_argument('--atomtype', type=str, default='lammps',
                        help='Atom which is placed on each gridpoint.')
    parser.add_argument('--out_file', type=str, default='uc_sampling.cube',
                        help='String specifying the output file (Gaussian cube format).')
    args = parser.parse_args()

    if args.nfrac and args.abs:
        raise ValueError('Specify one on --abs or --nfrac')

    atoms = ase.io.read(args.in_file)
    sampler = ucs.UnitCellSampler(atoms)
    if args.abs:
        sample = sampler.generate_grid_vectors(abs=args.abs)[0]
        energies = sampler.calculate_energies(
            method=args.method)  # From here the SQLite error originates

    energies = energies.reshape(sample.shape[:-1])
    with open(args.out_file, 'w') as fp:
        ase.io.cube.write_cube(fp, atoms, data=energies)


def get_partial_ion_charges():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_file',
                        help='File containing the unitcell to be sampled.')
    parser.add_argument('--fill_ions', type=str, default=None,
                        help='Fill unit cell with ions of specified type. All existing ions of this type will be removed.')
    parser.add_argument('--optimize', type=bool, default=False,
                        help='Optimize structure (default: False).')
    parser.add_argument('--core_electrons', type=str, default='{}',
                        help='JSON string specifying core electrons.')
    parser.add_argument('--method', type=str, default='repeat',
                        choices=['ddec', 'repeat'], help='Method to use for ')
    args = parser.parse_args()

    core_electrons = json.loads(args.core_electrons)

    in_path = Path(args.input_file)
    if in_path.suffix == '.xsf':
        atoms = ase.io.read(in_path, read_data=False)
        data = ase.io.read(in_path, read_data=True)
    elif in_path.suffix == '.cube':
        with open(in_path, 'r') as fp:
            cube_data = ase.io.cube.read_cube(fp)
            atoms = cube_data['atoms']
            data = cube_data['data']
    elif in_path.suffix == '.cif':
        atoms = ase.io.read(in_path)
        assert(not args.fill_ions)
    else:
        raise Exception('Unknown input file type.')

    if args.fill_ions:
        del atoms[[
            atom.index for atom in atoms
            if atom.symbol == args.fill_ions
        ]]
        atoms = fill_with_ions(atoms, data, type=args.fill_ions)

    inp = '''
    &FORCE_EVAL
    &DFT
    {}
    &SCF
    EPS_SCF 5.00000000E-005
    SCF_GUESS RESTART
    &OT T
        PRECONDITIONER FULL_SINGLE_INVERSE
        MINIMIZER CG
        LINESEARCH 3PNT
    &END OT
    &END SCF
    &PRINT
        &E_DENSITY_CUBE
            FILENAME =valence_density.cube
            STRIDE 1
        &END E_DENSITY_CUBE
    &END PRINT
    &END DFT
    &PROPERTIES
    &RESP
        USE_REPEAT_METHOD
        &SPHERE_SAMPLING
            AUTO_VDW_RADII_TABLE UFF
        &END
    &PRINT
        &RESP_CHARGES_TO_FILE
        &END
    &END PRINT
    &END RESP
    &END PROPERTIES
    &END FORCE_EVAL
    '''
    inp = inp.format('') if len(atoms) / 2 else inp.format('LSD')

    # CP2K.command = "env OMP_NUM_THREADS=32 srun cp2k_shell.psmp"
    # CP2K.command = "mpirun cp2k_shell.popt"
    if 'ASE_CP2K_COMMAND' not in os.environ:
        raise Exception(
            'ASE_CP2K_COMMAND not defined. Set environment variable!')

    calc = CP2K(inp=inp, pseudo_potential='GTH-PBE',
                max_scf=600, xc='PBE', print_level='LOW')
    atoms.calc = calc

    if args.optimize:
        dyn = BFGS(atoms, trajectory='traj.traj')
        dyn.run(fmax=0.05)
    else:
        atoms.get_potential_energy()

    if args.method == 'ddec':
        job_control_str = '''<net charge>
0.0 
</net charge>

<periodicity along A, B, and C vectors>
.true.
.true.
.true.
</periodicity along A, B, and C vectors>

<atomic densities directory complete path>
{}
</atomic densities directory complete path>

<number of core electrons>
{}</number of core electrons>
'''
        core_electrons_str = '\n'.join([str(k) + ' ' + str(v)
                                        for k, v in core_electrons.items()])
        core_electrons_str = core_electrons_str if len(
            core_electrons_str) == 1 else core_electrons_str + '\n'
        if 'DDEC_ATOMIC_DENSITY_DIR' not in os.environ:
            raise Exception(
                'DDEC_ATOMIC_DENSITY_DIR not defined. Set environment variable!')
        job_control_str = job_control_str.format(
            os.environ['DDEC_ATOMIC_DENSITY_DIR'], core_electrons_str)

        with open('job_control.txt', 'w') as fp:
            fp.write(job_control_str)
        subprocess.run(['env', 'OMP_NUM_THREADS=' + str(psutil.cpu_count()),
                        'Chargemol_09_26_2017_linux_parallel'])
        ion_charges = np.array([float(charge)
                                for (symbol, charge) in zip(*parse_ddec_file())
                                if symbol == 'Li'])
        print('Standart deviation:', ion_charges.std())
        print('Mean ion charge:', ion_charges.mean())
        return ion_charges.mean()

    elif args.method == 'repeat':
        # [charge for (symbol, charge) in zip(parse_repeat_file())
        #  if symbol == args.fill_ions]
        ion_charges = np.array([float(charge)
                                for (symbol, charge) in zip(*parse_repeat_file())
                                if symbol == 'Li'])
        print('Standart deviation:', ion_charges.std())
        print('Mean ion charge:', ion_charges.mean())
        return ion_charges.mean()


if __name__ == "__main__":
    run_mmc()
