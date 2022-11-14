#!/home/x_benbo/.pyenv/shims/python -u

from pathlib import Path
import unitcellsampling.sample
import ase.io

from ase.calculators.cp2k import CP2K
from ase.io.cube import write_cube

from unitcellsampling.decorators import subdir_calc
import os

#CP2K.command = "env OMP_NUM_THREADS=32 srun cp2k_shell.psmp"
CP2K.command = "env OMP_NUM_THREADS=4 cp2k"

# Path to previously relaxed single point calculation with an atom placed at
# the first sampled position in this case (0, 0, 0)
restart_file_name = "./single_point_calculation/cp2k-RESTART.wfn"


@subdir_calc
def dft_pbe(atoms: ase.Atoms):
    global restart_file_name

    # Use the same parameter here as in `single_point_calculation/`
    inp = '''
&FORCE_EVAL
  &DFT
    CHARGE -19
    LSD
    RESTART_FILE_NAME {}
    &SCF
       EPS_SCF 5.00000000E-005
       SCF_GUESS RESTART
       &OT T
           PRECONDITIONER FULL_SINGLE_INVERSE
           MINIMIZER CG
           LINESEARCH 3PNT
       &END OT
    &END SCF
  &END DFT
&END FORCE_EVAL
    '''.format(restart_file_name)

    # Change the restart file to the previous calculation
    p = Path('..').glob('**/*')
    dirs = [x for x in p if x.is_dir()]
    restart_file_name = str(
        Path(max(dirs, key=os.path.getmtime), 'cp2k-RESTART.wfn'))

    calc = CP2K(inp=inp, pseudo_potential='GTH-PBE',
                max_scf=600, xc='PBE', print_level='LOW')
    atoms.set_calculator(calc)
    return atoms.get_potential_energy()


indir = '.'
lgps = ase.io.read(Path(indir, 'lgps.xsf'))

sampler = unitcellsampling.sample.UnitCellSampler(lgps)
sampler.generate_grid_vectors(n_frac=(16, 16, 24), vdw_scale=1)

energies = sampler.calculate_energies(
    method=dft_pbe, atom="Li", exploit_symmetry=True)

with open('lgps.cube', 'w') as fp:
    write_cube(fp,  lgps, data=energies.reshape(
        sampler.included_grid_vectors.shape))
