import ase

import os
from pathlib import Path


# TODO: Add argument for getting information from 'master' calculation
def subdir_calc(calculate_energy):
    """
    Decorator for performing calculation in folder specified
    by environment variable UCS_CALCULATION_DIR.
    """
    def calculate_in_subdir(atoms: ase.Atoms):
        cwd = os.getcwd()
        calc_dir = os.getenv("UCS_CALCULATION_DIR")
        sub_dir_name = '_'.join(
            ['{:.3f}'.format(c)
             for c in atoms.get_scaled_positions()[-1]])
        sub_dir = Path(calc_dir, sub_dir_name)

        Path(calc_dir).mkdir(exist_ok=True, parents=True)
        sub_dir.mkdir(exist_ok=False)
        os.chdir(sub_dir)

        ase.io.write('atoms.cif', atoms)
        potential_energy = calculate_energy(atoms)

        os.chdir(str(cwd))
        return potential_energy

    return calculate_in_subdir
