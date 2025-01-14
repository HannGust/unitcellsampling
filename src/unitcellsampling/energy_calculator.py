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


# Test method for Li13Si4 structure - only for testing the symmetry quickly
# Charges and parameters are nonsense, but inspired from UFF forcefield
# For Si, I have used that of Ge. It is similar, but not exactlty, equal to those
# listed in the UFF publication.
# Modified the charges slightly to have charge neutrality.
def Li13Si4_test(atoms: ase.Atoms):

    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
	       "pair_modify tail yes mix arithmetic",
	       "special_bonds lj/coul 0.0 0.0 1.0",
	       "dielectric 1.0",
	       "kspace_style ewald 1.0e-5",
	       "pair_coeff 1 1 0.025000 2.183593",
	       "pair_coeff 2 2 0.379000 3.813047"]
    
    atom_types = {'Li': 1, 'Si': 2}

    atom_type_masses = {'Li': 6.941, 'Si':28.0855}

    log_file = None

    lammps_header = ["units real",
		     "atom_style full",
		     "boundary p p p"]

    amendments = ["set type 1 charge 0.70",
		  "set type 2 charge -2.275"]

    post_changebox_cmds = None

    calc = LAMMPSlib(lmpcmds=lmpcmds,
                     atom_types=atom_types,
                     atom_type_masses=atom_type_masses,
                     log_file=log_file,
                     lammps_header=lammps_header,
                     amendments=amendments,
                     post_changebox_cmds=post_changebox_cmds)
    
    atoms.set_calculator(calc)

    return atoms.get_potential_energy()


# The code below defines a classical lammps force-field calculator,
# which only contains interactions for the Li-ions. In this specific
# case, this constitute electrostatic (Ewald) interactions with 
# charge = +0.5e, and a UFF Li-Li Lennard-Jones potential.
# Can be generalized to any sampling ion, if it it useful ...

global_atom_types = {}
n=100
for i in range(1,n):
    global_atom_types[ase.atom.Atom(i).symbol] = i
   

assert global_atom_types['Li'] == 3 and global_atom_types['H'] == 1 and global_atom_types['U'] == 92, "Something wrong with atom types dict!"

#glob_atom_type_masses = {label:0.00 for label in global_atom_types.keys()}
glob_atom_type_masses = {label:ase.atom.Atom(label).mass for label in global_atom_types.keys()}

Li_charge = 0.50
global_amendment = []
for i in range(1,n):
    if i == 3:
        global_amendment.append("set type {} charge {}".format(i, Li_charge))
    else:
        global_amendment.append("set type {} charge {}".format(i, 0.00))


glob_pair_coeffs = []
for lbl in global_atom_types.keys():
    number = global_atom_types[lbl]

    if number == 3:
        glob_pair_coeffs.append("pair_coeff 3 3 0.044000 2.571134")
    else:
        glob_pair_coeffs.append("pair_coeff {} {} 0.0000 0.0000".format(number,number))


def only_Li_lammps_ewald_uff(atoms: ase.Atoms):
    """Classical force-field defined for all atoms but only containing nonzero
    force parameters (charge and a simple LJ) for Li. Consists of Ewald
    electrostatics as well as a Lennard-Jones potential based on UFF."""

    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
           "pair_modify tail yes mix arithmetic",
           "dielectric 1.0",
           "special_bonds lj/coul 0.0 0.0 1.0",
           "kspace_style ewald 1.0e-5"]
           #"pair_coeff 3 3 0.044000 2.571134"]

    lmpcmds.extend(glob_pair_coeffs)

    atom_types = global_atom_types

    atom_type_masses = glob_atom_type_masses 

    log_file = "only_Li.log"

    lammps_header = ["units real",
             "atom_style full",
             "boundary p p p"]

    amendments = global_amendment

    post_changebox_cmds = None

    calc = LAMMPSlib(lmpcmds=lmpcmds,
                     atom_types=atom_types,
                     atom_type_masses=atom_type_masses,
                     log_file=log_file,
                     lammps_header=lammps_header,
                     amendments=amendments,
                     post_changebox_cmds=post_changebox_cmds)

    atoms.set_calculator(calc)

    return atoms.get_potential_energy()


