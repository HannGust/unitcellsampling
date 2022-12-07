
####################################################################################
#                                                                                  #
# This file contains automatically generated LAMMPSlib-energy calculators for UCS. #
#                                                                                  #
####################################################################################
import ase
from ase.calculators.lammpsrun import LAMMPS
from ase.calculators.lammpslib import LAMMPSlib
from decorators import subdir_calc

# Manually defined ff for structure 54209
@subdir_calc
def struct_54209_ff(atoms: ase.Atoms):
    
    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
               "pair_modify     tail yes mix arithmetic",
               "special_bonds   lj/coul 0.0 0.0 1.0",
               "dielectric      1.0",
               "kspace_style ewald 1.0e-5",
               "pair_coeff 1 1 0.025000 2.183593",
               "pair_coeff 2 2 0.059000 2.819694",
               "pair_coeff 3 3 0.060000 3.118146"] 
    
    atom_types = {"Li":1, "Nb":2, "O":3}

    atom_type_masses = {"Li":6.9410, "Nb":92.906380, "O":15.99940}

    log_file = "54209_job.log"
    
    lammps_header = ["units real",
                     "atom_style full",
                     "boundary p p p",
                     "box tilt large"]

    amendments = ["set type 1 charge 0.92",
                  "set type 2 charge 2.59",
                  "set type 3 charge -1.17"]

    #post_changebox_cmds = []

    calc = LAMMPSlib(lmpcmds=lmpcmds,
                     atom_types=atom_types,
                     atom_type_masses=atom_type_masses,
                     log_file=log_file,
                     lammps_header=lammps_header,
                     amendments=amendments)
    
    atoms.set_calculator(calc)

    return atoms.get_potential_energy()
# Manually defined ff for structure 54297
@subdir_calc
def struct_54297_ff(atoms: ase.Atoms):
    
    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
               "pair_modify     tail yes mix arithmetic",
               "special_bonds   lj/coul 0.0 0.0 1.0",
               "dielectric      1.0",
               "kspace_style ewald 1.0e-5",
               "pair_coeff 1 1 0.025000 2.183593",
               "pair_coeff 2 2 0.036000 2.804549",
               "pair_coeff 3 3 0.060000 3.118146"] 
    
    atom_types = {"Li":1, "Ag":2, "O":3}

    atom_type_masses = {"Li":6.9410, "Ag":107.868200000, "O":15.999400000}

    log_file = "54297_job.log"
    
    lammps_header = ["units real",
                     "atom_style full",
                     "boundary p p p",
                     "box tilt large"]

    amendments = ["set type 1 charge 0.14",
                  "set type 2 charge -0.05",
                  "set type 3 charge -0.09"] # ? charges look a bit wierd

    #post_changebox_cmds = []

    calc = LAMMPSlib(lmpcmds=lmpcmds,
                     atom_types=atom_types,
                     atom_type_masses=atom_type_masses,
                     log_file=log_file,
                     lammps_header=lammps_header,
                     amendments=amendments)
    
    atoms.set_calculator(calc)

    return atoms.get_potential_energy()
# Manually defined ff for structure 54449
@subdir_calc
def struct_54449_ff(atoms: ase.Atoms):
    
    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
               "pair_modify     tail yes mix arithmetic",
               "special_bonds   lj/coul 0.0 0.0 1.0",
               "dielectric      1.0",
               "kspace_style ewald 1.0e-5",
               "pair_coeff 1 1 0.025000 2.183593",
               "pair_coeff 2 2 0.505000 4.008153",
               "pair_coeff 3 3 0.073000 2.530152"] 
    
    atom_types = {"Li":1, "Al":2, "Ir":3}

    atom_type_masses = {"Li":6.9410, "Al":26.981538600, "Ir":192.217000000}

    log_file = "54449_job.log"
    
    lammps_header = ["units real",
                     "atom_style full",
                     "boundary p p p",
                     "box tilt large"]

    amendments = ["set type 1 charge 1.25",
                  "set type 2 charge -5.29",
                  "set type 3 charge 1.39"]

    #post_changebox_cmds = []

    calc = LAMMPSlib(lmpcmds=lmpcmds,
                     atom_types=atom_types,
                     atom_type_masses=atom_type_masses,
                     log_file=log_file,
                     lammps_header=lammps_header,
                     amendments=amendments)
    
    atoms.set_calculator(calc)

    return atoms.get_potential_energy()
# Manually defined ff for structure 54865
@subdir_calc
def struct_54865_ff(atoms: ase.Atoms):
    
    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
               "pair_modify     tail yes mix arithmetic",
               "special_bonds   lj/coul 0.0 0.0 1.0",
               "dielectric      1.0",
               "kspace_style ewald 1.0e-5",
               "pair_coeff 1 1 0.025000 2.183593",
               "pair_coeff 2 2 0.048000 2.582715",
               "pair_coeff 3 3 0.060000 3.118146",
               "pair_coeff 4 4 0.060000 3.118146"] 
    
    atom_types = {"Li":1, "Pd":2, "O":3, "O":4}

    atom_type_masses = {"Li":6.9410, "Pd":106.420000000, "O":15.999400000}

    log_file = "54865_job.log"
    
    lammps_header = ["units real",
                     "atom_style full",
                     "boundary p p p",
                     "box tilt large"]

    amendments = ["set type 1 charge 0.16",
                  "set type 2 charge -0.53",
                  "set type 3 charge 0.11",
                  "set type 4 charge 0.11"] # ? very wierd charges

    #post_changebox_cmds = []

    calc = LAMMPSlib(lmpcmds=lmpcmds,
                     atom_types=atom_types,
                     atom_type_masses=atom_type_masses,
                     log_file=log_file,
                     lammps_header=lammps_header,
                     amendments=amendments)
    
    atoms.set_calculator(calc)

    return atoms.get_potential_energy()
# Manually defined ff for structure 55184
@subdir_calc
def struct_55184_ff(atoms: ase.Atoms):
    
    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
               "pair_modify     tail yes mix arithmetic",
               "special_bonds   lj/coul 0.0 0.0 1.0",
               "dielectric      1.0",
               "kspace_style ewald 1.0e-5",
               "pair_coeff 1 1 0.025000 2.183593",
               "pair_coeff 2 2 0.019000 2.935511",
               "pair_coeff 3 3 0.227000 3.516377",
               "pair_coeff 4 4 0.680000 3.872737"] 
    
    atom_types = {"Li":1, "Sc":2, "Cl":3, "Tl":4}

    atom_type_masses = {"Li":6.9410, "Sc":44.955912000, "Cl":35.453000000, "Tl":204.383300000}

    log_file = "55184_job.log"
    
    lammps_header = ["units real",
                     "atom_style full",
                     "boundary p p p",
                     "box tilt large"]

    amendments = ["set type 1 charge 0.9",
                  "set type 2 charge 1.29",
                  "set type 3 charge -0.49",
                  "set type 4 charge 0.38"]

    #post_changebox_cmds = []

    calc = LAMMPSlib(lmpcmds=lmpcmds,
                     atom_types=atom_types,
                     atom_type_masses=atom_type_masses,
                     log_file=log_file,
                     lammps_header=lammps_header,
                     amendments=amendments)
    
    atoms.set_calculator(calc)

    return atoms.get_potential_energy()
