import ase, ase.io
import numpy as np
from decorators import subdir_calc

# energy calculators
from ase.calculators.lammpsrun import LAMMPS
from ase.calculators.lammpslib import LAMMPSlib

############# THIS IS THE TEMPLATE
# Manually defined ff for one of the structures
@subdir_calc
def structure_54189_ff_manually_coded(atoms: ase.Atoms):
    
    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
               "pair_modify     tail yes mix arithmetic",
               "special_bonds   lj/coul 0.0 0.0 1.0",
               "dielectric      1.0",
               "kspace_style ewald 1.0e-5",
               "pair_coeff 1 1 0.025000 2.183593",
               "pair_coeff 2 2 0.124000 2.461553",
               "pair_coeff 3 3 0.379000 3.813047"] # "box tilt        large", went here before 
    
    atom_types = {"Li":1, "Zn":2, "Ge":3}
    atom_type_masses = {"Li":6.9410, "Zn":65.380, "Ge":72.640}
    log_file = "54189_job.log"

    # TODO: Note that "boundary p p p" and "box tilt large" are automatically added by ase in the inputfile. Can likely remove them here.
    
    lammps_header = ["units real",
                     "atom_style full",
                     "boundary p p p",
                     "box tilt large"]  # Added this line here instead

    amendments = ["set type 1 charge 0.72",
                  "set type 2 charge 1.5",
                  "set type 3 charge -2.95"]
    #post_changebox_cmds = []

    calc = LAMMPSlib(lmpcmds=lmpcmds,
                     atom_types=atom_types,
                     atom_type_masses=atom_type_masses,
                     log_file=log_file,
                     lammps_header=lammps_header,
                     amendments=amendments)
    
    atoms.set_calculator(calc)
    return atoms.get_potential_energy()
############################

##################################
#       Actual methods here      #
##################################

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

#


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

#



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

#



# Manually defined ff for structure 54683
@subdir_calc
def struct_54683_ff(atoms: ase.Atoms):
    
    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
               "pair_modify     tail yes mix arithmetic",
               "special_bonds   lj/coul 0.0 0.0 1.0",
               "dielectric      1.0",
               "kspace_style ewald 1.0e-5",
               "pair_coeff 1 1 0.025000 2.183593",
               "pair_coeff 2 2 0.080000 2.453535",
               "pair_coeff 3 3 0.060000 3.118146"] 
    
    atom_types = {"Li":1, "Pt":2, "O":3}
    atom_type_masses = {"Li":6.9410, "Pt":195.084000000, "O":15.999400000}
    log_file = "54683_job.log"
    
    lammps_header = ["units real",
                     "atom_style full",
                     "boundary p p p",
                     "box tilt large"]

    amendments = ["set type 1 charge 0.72",
                  "set type 2 charge 1.47",
                  "set type 3 charge -0.48"]

    #post_changebox_cmds = []

    calc = LAMMPSlib(lmpcmds=lmpcmds,
                     atom_types=atom_types,
                     atom_type_masses=atom_type_masses,
                     log_file=log_file,
                     lammps_header=lammps_header,
                     amendments=amendments)
    
    atoms.set_calculator(calc)
    return atoms.get_potential_energy()

#



# Manually defined ff for structure 54837
@subdir_calc
def struct_54837_ff(atoms: ase.Atoms):
    
    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
               "pair_modify     tail yes mix arithmetic",
               "special_bonds   lj/coul 0.0 0.0 1.0",
               "dielectric      1.0",
               "kspace_style ewald 1.0e-5",
               "pair_coeff 1 1 0.025000 2.183593",
               "pair_coeff 2 2 0.072000 2.798313",
               "pair_coeff 3 3 0.050000 2.996983"] 
    
    atom_types = {"Li":1, "Hf":2, "F":3}
    atom_type_masses = {"Li":6.9410, "Hf":178.490000000, "F":18.998403200}
    log_file = "54837_job.log"
    
    lammps_header = ["units real",
                     "atom_style full",
                     "boundary p p p",
                     "box tilt large"]

    amendments = ["set type 1 charge 0.92",
                  "set type 2 charge 2.28",
                  "set type 3 charge -0.69"]

    #post_changebox_cmds = []

    calc = LAMMPSlib(lmpcmds=lmpcmds,
                     atom_types=atom_types,
                     atom_type_masses=atom_type_masses,
                     log_file=log_file,
                     lammps_header=lammps_header,
                     amendments=amendments)
    
    atoms.set_calculator(calc)
    return atoms.get_potential_energy()

#



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

#



# Manually defined ff for structure 54879
@subdir_calc
def struct_54879_ff(atoms: ase.Atoms):
    
    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
               "pair_modify     tail yes mix arithmetic",
               "special_bonds   lj/coul 0.0 0.0 1.0",
               "dielectric      1.0",
               "kspace_style ewald 1.0e-5",
               "pair_coeff 1 1 0.025000 2.183593",
               "pair_coeff 2 2 0.072000 2.980056",
               "pair_coeff 3 3 0.518000 3.893227"] 
    
    atom_types = {"Li":1, "Y":2, "Bi":3}
    atom_type_masses = {"Li":6.9410, "Y":88.905850000, "Bi":208.980400000}
    log_file = "54879_job.log"
    
    lammps_header = ["units real",
                     "atom_style full",
                     "boundary p p p",
                     "box tilt large"]

    amendments = ["set type 1 charge 0.98",
                  "set type 2 charge 2.32",
                  "set type 3 charge -2.63"]

    #post_changebox_cmds = []

    calc = LAMMPSlib(lmpcmds=lmpcmds,
                     atom_types=atom_types,
                     atom_type_masses=atom_type_masses,
                     log_file=log_file,
                     lammps_header=lammps_header,
                     amendments=amendments)
    
    atoms.set_calculator(calc)
    return atoms.get_potential_energy()

#



# Manually defined ff for structure 54884
@subdir_calc
def struct_54884_ff(atoms: ase.Atoms):
    
    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
               "pair_modify     tail yes mix arithmetic",
               "special_bonds   lj/coul 0.0 0.0 1.0",
               "dielectric      1.0",
               "kspace_style ewald 1.0e-5",
               "pair_coeff 1 1 0.040000 3.665157",
               "pair_coeff 2 2 0.025000 2.183593",
               "pair_coeff 3 3 0.019000 2.935511",
               "pair_coeff 4 4 0.227000 3.516377"] 
    
    atom_types = {"Rb":1, "Li":2, "Sc":3, "Cl":4}
    atom_type_masses = {"Rb":85.467800000 , "Li":6.9410, "Sc":44.955912000, "Cl":35.453000000}
    log_file = "54884_job.log"
    
    lammps_header = ["units real",
                     "atom_style full",
                     "boundary p p p",
                     "box tilt large"]

    amendments = ["set type 1 charge 0.78",
                  "set type 2 charge 0.95",
                  "set type 3 charge 1.34",
                  "set type 4 charge -0.64"]

    #post_changebox_cmds = []

    calc = LAMMPSlib(lmpcmds=lmpcmds,
                     atom_types=atom_types,
                     atom_type_masses=atom_type_masses,
                     log_file=log_file,
                     lammps_header=lammps_header,
                     amendments=amendments)
    
    atoms.set_calculator(calc)
    return atoms.get_potential_energy()

#



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

#



# Manually defined ff for structure 55319
@subdir_calc
def struct_55319_ff(atoms: ase.Atoms):
    
    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
               "pair_modify     tail yes mix arithmetic",
               "special_bonds   lj/coul 0.0 0.0 1.0",
               "dielectric      1.0",
               "kspace_style ewald 1.0e-5",
               "pair_coeff 1 1 0.025000 2.183593",
               "pair_coeff 2 2 0.015000 2.524807",
               "pair_coeff 3 3 0.069000 3.260689"] 
    
    atom_types = {"Li":1, "Ni":2, "N":3}
    atom_type_masses = {"Li":6.9410, "Ni":58.693400000, "N":14.006700000}
    log_file = "55319_job.log"
    
    lammps_header = ["units real",
                     "atom_style full",
                     "boundary p p p",
                     "box tilt large"]

    amendments = ["set type 1 charge 0.59",
                  "set type 2 charge 0.03",
                  "set type 3 charge -0.63"]

    #post_changebox_cmds = []

    calc = LAMMPSlib(lmpcmds=lmpcmds,
                     atom_types=atom_types,
                     atom_type_masses=atom_type_masses,
                     log_file=log_file,
                     lammps_header=lammps_header,
                     amendments=amendments)
    
    atoms.set_calculator(calc)
    return atoms.get_potential_energy()

#

#######
####### ABOVE IS "10-STRUCT" SUBSET OF METHODS (EXCEPT THE VERY FIRST ONE, WHICH IS THE 54189 STRUCTURE METHOD)
#######


# Manually defined ff for structure 56568
@subdir_calc
def struct_56568_ff(atoms: ase.Atoms):
    
    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
               "pair_modify     tail yes mix arithmetic",
               "special_bonds   lj/coul 0.0 0.0 1.0",
               "dielectric      1.0",
               "kspace_style ewald 1.0e-5",
               "pair_coeff 1 1 0.025000 2.183593",
               "pair_coeff 2 2 0.017000 3.137745",
               "pair_coeff 3 3 0.379000 3.813047"] 
    
    atom_types = {"Li":1, "La":2, "Ge":3}
    atom_type_masses = {"Li":6.9410, "La":138.905470000, "Ge":72.640000000}
    log_file = "56568_job.log"
    
    lammps_header = ["units real",
                     "atom_style full",
                     "boundary p p p",
                     "box tilt large"]

    amendments = ["set type 1 charge 0.91",
                  "set type 2 charge 1.01",
                  "set type 3 charge -2.83"]

    #post_changebox_cmds = []

    calc = LAMMPSlib(lmpcmds=lmpcmds,
                     atom_types=atom_types,
                     atom_type_masses=atom_type_masses,
                     log_file=log_file,
                     lammps_header=lammps_header,
                     amendments=amendments)
    
    atoms.set_calculator(calc)
    return atoms.get_potential_energy()

#


# Manually defined ff for structure 56627
@subdir_calc
def struct_56627_ff(atoms: ase.Atoms):
    
    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
               "pair_modify     tail yes mix arithmetic",
               "special_bonds   lj/coul 0.0 0.0 1.0",
               "dielectric      1.0",
               "kspace_style ewald 1.0e-5",
               "pair_coeff 1 1 0.025000 2.183593",
               "pair_coeff 2 2 0.014000 2.558661",
               "pair_coeff 3 3 0.449000 3.937772"] 
    
    atom_types = {"Li":1, "Co":2, "Sb":3}
    atom_type_masses = {"Li":6.9410, "Co":58.933195000, "Sb":121.760000000}
    log_file = "56627_job.log"
    
    lammps_header = ["units real",
                     "atom_style full",
                     "boundary p p p",
                     "box tilt large"]

    amendments = ["set type 1 charge 0.79",
                  "set type 2 charge 0.37",
                  "set type 3 charge -1.94"]

    #post_changebox_cmds = []

    calc = LAMMPSlib(lmpcmds=lmpcmds,
                     atom_types=atom_types,
                     atom_type_masses=atom_type_masses,
                     log_file=log_file,
                     lammps_header=lammps_header,
                     amendments=amendments)
    
    atoms.set_calculator(calc)
    return atoms.get_potential_energy()

#


# Manually defined ff for structure 57382
@subdir_calc
def struct_57382_ff(atoms: ase.Atoms):
    
    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
               "pair_modify     tail yes mix arithmetic",
               "special_bonds   lj/coul 0.0 0.0 1.0",
               "dielectric      1.0",
               "kspace_style ewald 1.0e-5",
               "pair_coeff 1 1 0.040000 3.665157",
               "pair_coeff 2 2 0.025000 2.183593",
               "pair_coeff 3 3 0.291000 3.746229"] 
    
    atom_types = {"Rb":1, "Li":2, "Se":3}
    atom_type_masses = {"Rb":85.467800000, "Li":6.9410, "Se":78.960000000}
    log_file = "57382_job.log"
    
    lammps_header = ["units real",
                     "atom_style full",
                     "boundary p p p",
                     "box tilt large"]

    amendments = ["set type 1 charge 0.61",
                  "set type 2 charge 0.63",
                  "set type 3 charge -1.25"]

    #post_changebox_cmds = []

    calc = LAMMPSlib(lmpcmds=lmpcmds,
                     atom_types=atom_types,
                     atom_type_masses=atom_type_masses,
                     log_file=log_file,
                     lammps_header=lammps_header,
                     amendments=amendments)
    
    atoms.set_calculator(calc)
    return atoms.get_potential_energy()

#







# Manually defined ff for structure 73679
@subdir_calc
def struct_73679_ff(atoms: ase.Atoms):
    
    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
               "pair_modify     tail yes mix arithmetic",
               "special_bonds   lj/coul 0.0 0.0 1.0",
               "dielectric      1.0",
               "kspace_style ewald 1.0e-5",
               "pair_coeff 1 1 0.025000 2.183593",
               "pair_coeff 2 2 0.505000 4.008153",
               "pair_coeff 3 3 0.053000 2.609442"] 
    
    atom_types = {"Li":1, "Al":2, "Rh":3}
    atom_type_masses = {"Li":6.9410, "Al":26.981538600, "Rh":102.905500000}
    log_file = "73679_job.log"
    
    lammps_header = ["units real",
                     "atom_style full",
                     "boundary p p p",
                     "box tilt large"]

    amendments = ["set type 1 charge 1.23",
                  "set type 2 charge -4.11",
                  "set type 3 charge 0.82"]

    #post_changebox_cmds = []

    calc = LAMMPSlib(lmpcmds=lmpcmds,
                     atom_types=atom_types,
                     atom_type_masses=atom_type_masses,
                     log_file=log_file,
                     lammps_header=lammps_header,
                     amendments=amendments)
    
    atoms.set_calculator(calc)
    return atoms.get_potential_energy()

#
