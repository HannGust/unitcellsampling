import ase
from ase.calculators.lammpsrun import LAMMPS
from ase.calculators.lammpslib import LAMMPSlib
from ase.calculators.cp2k import CP2K

from unitcellsampling.decorators import subdir_calc


@subdir_calc
def cp2k_dft(atoms: ase.Atoms):
    calc = CP2K(inp='''&FORCE_EVAL
    &DFT
        LSD
    &END DFT
&END FORCE_EVAL''')
    atoms.set_calculator(calc)
    return atoms.get_potential_energy()


def lammps_forcefield_simple(atoms: ase.Atoms):
    # Here no decorator necessary because directory for calculation needed
    atoms.calc = LAMMPS()
    return atoms.get_potential_energy()


@subdir_calc
def lammps_forcefield(atoms: ase.Atoms):
    cmds = ["pair_style eam",
            "pair_coeff * * Cu_u3.eam"]

    atoms.calc = LAMMPSlib(lmpcmds=cmds,
                           log_file='log.log')

    return atoms.get_potential_energy()
