
####################################################################################
#                                                                                  #
# This file contains automatically generated LAMMPSlib-energy calculators for UCS. #
#                                                                                  #
####################################################################################
import ase
from unitcellsampling.decorators import subdir_calc

# energy calculators
from ase.calculators.lammpsrun import LAMMPS
from ase.calculators.lammpslib import LAMMPSlib

#####   Automatically generated lammps method m54189_auto   #####
@subdir_calc
def m54189_auto(atoms: ase.Atoms):

    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
	       "pair_modify tail yes mix arithmetic",
	       "special_bonds lj/coul 0.0 0.0 1.0",
	       "dielectric 1.0",
	       "kspace_style ewald 1.0e-5",
	       "pair_coeff 1 1 0.025000 2.183593",
	       "pair_coeff 2 2 0.124000 2.461553",
	       "pair_coeff 3 3 0.379000 3.813047"]
    
    atom_types = {'Li': 1, 'Zn': 2, 'Ge': 3}

    atom_type_masses = {'Li': 6.941, 'Zn': 65.38, 'Ge': 72.64}

    log_file = "m54189_auto.log"

    lammps_header = ["units real",
		     "atom_style full",
		     "boundary p p p"]

    amendments = ["set type 1 charge 0.72",
		  "set type 2 charge 1.5",
		  "set type 3 charge -2.95"]

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
######   Automatically generated lammps method m54209_auto   #####
@subdir_calc
def m54209_auto(atoms: ase.Atoms):

    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
	       "special_bonds lj/coul 0.0 0.0 1.0",
	       "dielectric 1.0",
	       "pair_modify tail yes mix arithmetic",
	       "kspace_style ewald 1.0e-5",
	       "pair_coeff 1 1 0.025000 2.183593",
	       "pair_coeff 2 2 0.059000 2.819694",
	       "pair_coeff 3 3 0.060000 3.118146"]
    
    atom_types = {'Li': 1, 'Nb': 2, 'O': 3}

    atom_type_masses = {'Li': 6.941, 'Nb': 92.90638, 'O': 15.9994}

    log_file = "m54209_auto.log"

    lammps_header = ["units real",
		     "atom_style full",
		     "boundary p p p"]

    amendments = ["set type 1 charge 0.92",
		  "set type 2 charge 2.59",
		  "set type 3 charge -1.17"]

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
######   Automatically generated lammps method m54297_auto   #####
@subdir_calc
def m54297_auto(atoms: ase.Atoms):

    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
	       "dielectric 1.0",
	       "pair_modify tail yes mix arithmetic",
	       "special_bonds lj/coul 0.0 0.0 1.0",
	       "kspace_style ewald 1.0e-5",
	       "pair_coeff 1 1 0.025000 2.183593",
	       "pair_coeff 2 2 0.036000 2.804549",
	       "pair_coeff 3 3 0.060000 3.118146"]
    
    atom_types = {'Li': 1, 'Ag': 2, 'O': 3}

    atom_type_masses = {'Li': 6.941, 'Ag': 107.8682, 'O': 15.9994}

    log_file = "m54297_auto.log"

    lammps_header = ["units real",
		     "atom_style full",
		     "boundary p p p"]

    amendments = ["set type 1 charge 0.14",
		  "set type 2 charge -0.05",
		  "set type 3 charge -0.09"]

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
######   Automatically generated lammps method m54449_auto   #####
@subdir_calc
def m54449_auto(atoms: ase.Atoms):

    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
	       "dielectric 1.0",
	       "special_bonds lj/coul 0.0 0.0 1.0",
	       "pair_modify tail yes mix arithmetic",
	       "kspace_style ewald 1.0e-5",
	       "pair_coeff 1 1 0.025000 2.183593",
	       "pair_coeff 2 2 0.505000 4.008153",
	       "pair_coeff 3 3 0.073000 2.530152"]
    
    atom_types = {'Li': 1, 'Al': 2, 'Ir': 3}

    atom_type_masses = {'Li': 6.941, 'Al': 26.9815386, 'Ir': 192.217}

    log_file = "m54449_auto.log"

    lammps_header = ["units real",
		     "atom_style full",
		     "boundary p p p"]

    amendments = ["set type 1 charge 1.25",
		  "set type 2 charge -5.29",
		  "set type 3 charge 1.39"]

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
######   Automatically generated lammps method m54683_auto   #####
@subdir_calc
def m54683_auto(atoms: ase.Atoms):

    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
	       "dielectric 1.0",
	       "pair_modify tail yes mix arithmetic",
	       "special_bonds lj/coul 0.0 0.0 1.0",
	       "kspace_style ewald 1.0e-5",
	       "pair_coeff 1 1 0.025000 2.183593",
	       "pair_coeff 2 2 0.080000 2.453535",
	       "pair_coeff 3 3 0.060000 3.118146"]
    
    atom_types = {'Li': 1, 'Pt': 2, 'O': 3}

    atom_type_masses = {'Li': 6.941, 'Pt': 195.084, 'O': 15.9994}

    log_file = "m54683_auto.log"

    lammps_header = ["units real",
		     "atom_style full",
		     "boundary p p p"]

    amendments = ["set type 1 charge 0.72",
		  "set type 2 charge 1.47",
		  "set type 3 charge -0.48"]

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
######   Automatically generated lammps method m54837_auto   #####
@subdir_calc
def m54837_auto(atoms: ase.Atoms):

    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
	       "dielectric 1.0",
	       "pair_modify tail yes mix arithmetic",
	       "special_bonds lj/coul 0.0 0.0 1.0",
	       "kspace_style ewald 1.0e-5",
	       "pair_coeff 1 1 0.025000 2.183593",
	       "pair_coeff 2 2 0.072000 2.798313",
	       "pair_coeff 3 3 0.050000 2.996983"]
    
    atom_types = {'Li': 1, 'Hf': 2, 'F': 3}

    atom_type_masses = {'Li': 6.941, 'Hf': 178.49, 'F': 18.9984032}

    log_file = "m54837_auto.log"

    lammps_header = ["units real",
		     "atom_style full",
		     "boundary p p p"]

    amendments = ["set type 1 charge 0.92",
		  "set type 2 charge 2.28",
		  "set type 3 charge -0.69"]

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
######   Automatically generated lammps method m54865_auto   #####
@subdir_calc
def m54865_auto(atoms: ase.Atoms):

    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
	       "pair_modify tail yes mix arithmetic",
	       "special_bonds lj/coul 0.0 0.0 1.0",
	       "dielectric 1.0",
	       "kspace_style ewald 1.0e-5",
	       "pair_coeff 1 1 0.025000 2.183593",
	       "pair_coeff 2 2 0.048000 2.582715",
	       "pair_coeff 3 3 0.060000 3.118146",
	       "pair_coeff 4 4 0.060000 3.118146"]
    
    atom_types = {'Li': 1, 'Pd': 2, 'O': 4}

    atom_type_masses = {'Li': 6.941, 'Pd': 106.42, 'O': 15.9994}

    log_file = "m54865_auto.log"

    lammps_header = ["units real",
		     "atom_style full",
		     "boundary p p p"]

    amendments = ["set type 1 charge 0.16",
		  "set type 2 charge -0.53",
		  "set type 3 charge 0.11",
		  "set type 4 charge 0.11"]

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
######   Automatically generated lammps method m54879_auto   #####
@subdir_calc
def m54879_auto(atoms: ase.Atoms):

    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
	       "dielectric 1.0",
	       "special_bonds lj/coul 0.0 0.0 1.0",
	       "pair_modify tail yes mix arithmetic",
	       "kspace_style ewald 1.0e-5",
	       "pair_coeff 1 1 0.025000 2.183593",
	       "pair_coeff 2 2 0.072000 2.980056",
	       "pair_coeff 3 3 0.518000 3.893227"]
    
    atom_types = {'Li': 1, 'Y': 2, 'Bi': 3}

    atom_type_masses = {'Li': 6.941, 'Y': 88.90585, 'Bi': 208.9804}

    log_file = "m54879_auto.log"

    lammps_header = ["units real",
		     "atom_style full",
		     "boundary p p p"]

    amendments = ["set type 1 charge 0.98",
		  "set type 2 charge 2.32",
		  "set type 3 charge -2.63"]

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
######   Automatically generated lammps method m54884_auto   #####
@subdir_calc
def m54884_auto(atoms: ase.Atoms):

    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
	       "dielectric 1.0",
	       "special_bonds lj/coul 0.0 0.0 1.0",
	       "pair_modify tail yes mix arithmetic",
	       "kspace_style ewald 1.0e-5",
	       "pair_coeff 1 1 0.040000 3.665157",
	       "pair_coeff 2 2 0.025000 2.183593",
	       "pair_coeff 3 3 0.019000 2.935511",
	       "pair_coeff 4 4 0.227000 3.516377"]
    
    atom_types = {'Rb': 1, 'Li': 2, 'Sc': 3, 'Cl': 4}

    atom_type_masses = {'Rb': 85.4678, 'Li': 6.941, 'Sc': 44.955912, 'Cl': 35.453}

    log_file = "m54884_auto.log"

    lammps_header = ["units real",
		     "atom_style full",
		     "boundary p p p"]

    amendments = ["set type 1 charge 0.78",
		  "set type 2 charge 0.95",
		  "set type 3 charge 1.34",
		  "set type 4 charge -0.64"]

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
######   Automatically generated lammps method m55184_auto   #####
@subdir_calc
def m55184_auto(atoms: ase.Atoms):

    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
	       "pair_modify tail yes mix arithmetic",
	       "dielectric 1.0",
	       "special_bonds lj/coul 0.0 0.0 1.0",
	       "kspace_style ewald 1.0e-5",
	       "pair_coeff 1 1 0.025000 2.183593",
	       "pair_coeff 2 2 0.019000 2.935511",
	       "pair_coeff 3 3 0.227000 3.516377",
	       "pair_coeff 4 4 0.680000 3.872737"]
    
    atom_types = {'Li': 1, 'Sc': 2, 'Cl': 3, 'Tl': 4}

    atom_type_masses = {'Li': 6.941, 'Sc': 44.955912, 'Cl': 35.453, 'Tl': 204.3833}

    log_file = "m55184_auto.log"

    lammps_header = ["units real",
		     "atom_style full",
		     "boundary p p p"]

    amendments = ["set type 1 charge 0.9",
		  "set type 2 charge 1.29",
		  "set type 3 charge -0.49",
		  "set type 4 charge 0.38"]

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
######   Automatically generated lammps method m55319_auto   #####
@subdir_calc
def m55319_auto(atoms: ase.Atoms):

    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
	       "pair_modify tail yes mix arithmetic",
	       "dielectric 1.0",
	       "special_bonds lj/coul 0.0 0.0 1.0",
	       "kspace_style ewald 1.0e-5",
	       "pair_coeff 1 1 0.025000 2.183593",
	       "pair_coeff 2 2 0.015000 2.524807",
	       "pair_coeff 3 3 0.069000 3.260689"]
    
    atom_types = {'Li': 1, 'Ni': 2, 'N': 3}

    atom_type_masses = {'Li': 6.941, 'Ni': 58.6934, 'N': 14.0067}

    log_file = "m55319_auto.log"

    lammps_header = ["units real",
		     "atom_style full",
		     "boundary p p p"]

    amendments = ["set type 1 charge 0.59",
		  "set type 2 charge -0.63",
		  "set type 2 charge 0.03",
		  "set type 3 charge -0.63"]

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
######   Automatically generated lammps method m55983_auto   #####
@subdir_calc
def m55983_auto(atoms: ase.Atoms):

    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
	       "dielectric 1.0",
	       "pair_modify tail yes mix arithmetic",
	       "special_bonds lj/coul 0.0 0.0 1.0",
	       "kspace_style ewald 1.0e-5",
	       "pair_coeff 1 1 0.025000 2.183593",
	       "pair_coeff 2 2 0.039000 2.933729",
	       "pair_coeff 3 3 0.274000 3.594776",
	       "pair_coeff 4 4 0.039000 2.933729",
	       "pair_coeff 5 5 0.274000 3.594776"]
    
    atom_types = {'Li': 1, 'Au': 4, 'S': 5}

    atom_type_masses = {'Li': 6.941, 'Au': 196.966569, 'S': 32.065}

    log_file = "m55983_auto.log"

    lammps_header = ["units real",
		     "atom_style full",
		     "boundary p p p"]

    amendments = ["set type 1 charge 0.74",
		  "set type 2 charge 0.28",
		  "set type 3 charge -1.02",
		  "set type 4 charge 0.28",
		  "set type 5 charge -1.02"]

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
######   Automatically generated lammps method m56568_auto   #####
@subdir_calc
def m56568_auto(atoms: ase.Atoms):

    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
	       "pair_modify tail yes mix arithmetic",
	       "dielectric 1.0",
	       "special_bonds lj/coul 0.0 0.0 1.0",
	       "kspace_style ewald 1.0e-5",
	       "pair_coeff 1 1 0.025000 2.183593",
	       "pair_coeff 2 2 0.017000 3.137745",
	       "pair_coeff 3 3 0.379000 3.813047"]
    
    atom_types = {'Li': 1, 'La': 2, 'Ge': 3}

    atom_type_masses = {'Li': 6.941, 'La': 138.90547, 'Ge': 72.64}

    log_file = "m56568_auto.log"

    lammps_header = ["units real",
		     "atom_style full",
		     "boundary p p p"]

    amendments = ["set type 1 charge 0.91",
		  "set type 2 charge 1.01",
		  "set type 3 charge -2.83"]

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
######   Automatically generated lammps method m56627_auto   #####
@subdir_calc
def m56627_auto(atoms: ase.Atoms):

    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
	       "special_bonds lj/coul 0.0 0.0 1.0",
	       "pair_modify tail yes mix arithmetic",
	       "dielectric 1.0",
	       "kspace_style ewald 1.0e-5",
	       "pair_coeff 1 1 0.025000 2.183593",
	       "pair_coeff 2 2 0.014000 2.558661",
	       "pair_coeff 3 3 0.449000 3.937772"]
    
    atom_types = {'Li': 1, 'Co': 2, 'Sb': 3}

    atom_type_masses = {'Li': 6.941, 'Co': 58.933195, 'Sb': 121.76}

    log_file = "m56627_auto.log"

    lammps_header = ["units real",
		     "atom_style full",
		     "boundary p p p"]

    amendments = ["set type 1 charge 0.79",
		  "set type 1 charge 0.79",
		  "set type 2 charge 0.37",
		  "set type 2 charge 0.37",
		  "set type 3 charge -1.94",
		  "set type 3 charge -1.94"]

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
######   Automatically generated lammps method m57347_auto   #####
@subdir_calc
def m57347_auto(atoms: ase.Atoms):

    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
	       "pair_modify tail yes mix arithmetic",
	       "special_bonds lj/coul 0.0 0.0 1.0",
	       "dielectric 1.0",
	       "kspace_style ewald 1.0e-5",
	       "pair_coeff 1 1 0.025000 2.183593",
	       "pair_coeff 2 2 0.505000 4.008153",
	       "pair_coeff 3 3 0.044000 2.571134",
	       "pair_coeff 4 4 0.044000 2.571134"]
    
    atom_types = {'Li': 1, 'Al': 2, 'H': 4}

    atom_type_masses = {'Li': 6.941, 'Al': 26.9815386, 'H': 1.00794}

    log_file = "m57347_auto.log"

    lammps_header = ["units real",
		     "atom_style full",
		     "boundary p p p"]

    amendments = ["set type 1 charge 0.96",
		  "set type 2 charge 1.34",
		  "set type 3 charge -0.58",
		  "set type 4 charge -0.58"]

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
######   Automatically generated lammps method m57382_auto   #####
@subdir_calc
def m57382_auto(atoms: ase.Atoms):

    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
	       "special_bonds lj/coul 0.0 0.0 1.0",
	       "pair_modify tail yes mix arithmetic",
	       "dielectric 1.0",
	       "kspace_style ewald 1.0e-5",
	       "pair_coeff 1 1 0.040000 3.665157",
	       "pair_coeff 2 2 0.025000 2.183593",
	       "pair_coeff 3 3 0.291000 3.746229"]
    
    atom_types = {'Rb': 1, 'Li': 2, 'Se': 3}

    atom_type_masses = {'Rb': 85.4678, 'Li': 6.941, 'Se': 78.96}

    log_file = "m57382_auto.log"

    lammps_header = ["units real",
		     "atom_style full",
		     "boundary p p p"]

    amendments = ["set type 1 charge 0.61",
		  "set type 1 charge 0.61",
		  "set type 2 charge 0.63",
		  "set type 2 charge 0.63",
		  "set type 3 charge -1.25",
		  "set type 3 charge -1.25"]

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
######   Automatically generated lammps method m57448_auto   #####
@subdir_calc
def m57448_auto(atoms: ase.Atoms):

    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
	       "pair_modify tail yes mix arithmetic",
	       "dielectric 1.0",
	       "special_bonds lj/coul 0.0 0.0 1.0",
	       "kspace_style ewald 1.0e-5",
	       "pair_coeff 1 1 0.025000 2.183593",
	       "pair_coeff 2 2 0.069000 2.783168",
	       "pair_coeff 3 3 0.069000 3.260689"]
    
    atom_types = {'Li': 1, 'Zr': 2, 'N': 3}

    atom_type_masses = {'Li': 6.941, 'Zr': 91.224, 'N': 14.0067}

    log_file = "m57448_auto.log"

    lammps_header = ["units real",
		     "atom_style full",
		     "boundary p p p"]

    amendments = ["set type 1 charge 1.0",
		  "set type 2 charge 2.8",
		  "set type 3 charge -2.39"]

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
######   Automatically generated lammps method m57608_auto   #####
@subdir_calc
def m57608_auto(atoms: ase.Atoms):

    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
	       "dielectric 1.0",
	       "pair_modify tail yes mix arithmetic",
	       "special_bonds lj/coul 0.0 0.0 1.0",
	       "kspace_style ewald 1.0e-5",
	       "pair_coeff 1 1 0.025000 2.183593",
	       "pair_coeff 2 2 0.072000 2.798313",
	       "pair_coeff 3 3 0.228000 2.537280",
	       "pair_coeff 4 4 0.050000 2.996983"]
    
    atom_types = {'Li': 1, 'Hf': 2, 'Cd': 3, 'F': 4}

    atom_type_masses = {'Li': 6.941, 'Hf': 178.49, 'Cd': 112.411, 'F': 18.9984032}

    log_file = "m57608_auto.log"

    lammps_header = ["units real",
		     "atom_style full",
		     "boundary p p p"]

    amendments = ["set type 1 charge 1.03",
		  "set type 2 charge 2.67",
		  "set type 3 charge 1.7",
		  "set type 4 charge -0.8"]

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
######   Automatically generated lammps method m57644_auto   #####
@subdir_calc
def m57644_auto(atoms: ase.Atoms):

    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
	       "dielectric 1.0",
	       "pair_modify tail yes mix arithmetic",
	       "special_bonds lj/coul 0.0 0.0 1.0",
	       "kspace_style ewald 1.0e-5",
	       "pair_coeff 1 1 0.025000 2.183593",
	       "pair_coeff 2 2 0.085000 2.445517",
	       "pair_coeff 3 3 0.309000 3.768502"]
    
    atom_types = {'Li': 1, 'Be': 2, 'As': 3}

    atom_type_masses = {'Li': 6.941, 'Be': 9.012182, 'As': 74.9216}

    log_file = "m57644_auto.log"

    lammps_header = ["units real",
		     "atom_style full",
		     "boundary p p p"]

    amendments = ["set type 1 charge 0.84",
		  "set type 2 charge 0.71",
		  "set type 3 charge -1.55"]

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
######   Automatically generated lammps method m57761_auto   #####
@subdir_calc
def m57761_auto(atoms: ase.Atoms):

    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
	       "pair_modify tail yes mix arithmetic",
	       "dielectric 1.0",
	       "special_bonds lj/coul 0.0 0.0 1.0",
	       "kspace_style ewald 1.0e-5",
	       "pair_coeff 1 1 0.025000 2.183593",
	       "pair_coeff 2 2 0.505000 4.008153",
	       "pair_coeff 3 3 0.044000 2.571134",
	       "pair_coeff 4 4 0.044000 2.571134"]
    
    atom_types = {'Li': 1, 'Al': 2, 'H': 4}

    atom_type_masses = {'Li': 6.941, 'Al': 26.9815386, 'H': 1.00794}

    log_file = "m57761_auto.log"

    lammps_header = ["units real",
		     "atom_style full",
		     "boundary p p p"]

    amendments = ["set type 1 charge 0.82",
		  "set type 2 charge 1.47",
		  "set type 3 charge -0.57",
		  "set type 4 charge -0.57"]

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
######   Automatically generated lammps method m59631_auto   #####
@subdir_calc
def m59631_auto(atoms: ase.Atoms):

    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
	       "special_bonds lj/coul 0.0 0.0 1.0",
	       "dielectric 1.0",
	       "pair_modify tail yes mix arithmetic",
	       "kspace_style ewald 1.0e-5",
	       "pair_coeff 1 1 0.025000 2.183593",
	       "pair_coeff 2 2 0.005000 3.113691",
	       "pair_coeff 3 3 0.449000 3.937772"]
    
    atom_types = {'Li': 1, 'Cu': 2, 'Sb': 3}

    atom_type_masses = {'Li': 6.941, 'Cu': 63.546, 'Sb': 121.76}

    log_file = "m59631_auto.log"

    lammps_header = ["units real",
		     "atom_style full",
		     "boundary p p p"]

    amendments = ["set type 1 charge 0.7",
		  "set type 2 charge 0.94",
		  "set type 3 charge -2.35"]

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
######   Automatically generated lammps method m59632_auto   #####
@subdir_calc
def m59632_auto(atoms: ase.Atoms):

    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
	       "dielectric 1.0",
	       "pair_modify tail yes mix arithmetic",
	       "special_bonds lj/coul 0.0 0.0 1.0",
	       "kspace_style ewald 1.0e-5",
	       "pair_coeff 1 1 0.025000 2.183593",
	       "pair_coeff 2 2 0.599000 3.976081",
	       "pair_coeff 3 3 0.291000 3.746229"]
    
    atom_types = {'Li': 1, 'In': 2, 'Se': 3}

    atom_type_masses = {'Li': 6.941, 'In': 114.818, 'Se': 78.96}

    log_file = "m59632_auto.log"

    lammps_header = ["units real",
		     "atom_style full",
		     "boundary p p p"]

    amendments = ["set type 1 charge 0.69",
		  "set type 2 charge 0.74",
		  "set type 3 charge -0.72"]

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
######   Automatically generated lammps method m59715_auto   #####
@subdir_calc
def m59715_auto(atoms: ase.Atoms):

    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
	       "pair_modify tail yes mix arithmetic",
	       "dielectric 1.0",
	       "special_bonds lj/coul 0.0 0.0 1.0",
	       "kspace_style ewald 1.0e-5",
	       "pair_coeff 1 1 0.025000 2.183593",
	       "pair_coeff 2 2 0.599000 3.976081",
	       "pair_coeff 3 3 0.567000 3.912827"]
    
    atom_types = {'Li': 1, 'In': 2, 'Sn': 3}

    atom_type_masses = {'Li': 6.941, 'In': 114.818, 'Sn': 118.71}

    log_file = "m59715_auto.log"

    lammps_header = ["units real",
		     "atom_style full",
		     "boundary p p p"]

    amendments = ["set type 1 charge 0.9",
		  "set type 2 charge -1.41",
		  "set type 3 charge 0.52"]

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
######   Automatically generated lammps method m59948_auto   #####
@subdir_calc
def m59948_auto(atoms: ase.Atoms):

    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
	       "dielectric 1.0",
	       "special_bonds lj/coul 0.0 0.0 1.0",
	       "pair_modify tail yes mix arithmetic",
	       "kspace_style ewald 1.0e-5",
	       "pair_coeff 1 1 0.025000 2.183593",
	       "pair_coeff 2 2 0.044000 2.571134",
	       "pair_coeff 3 3 0.053000 2.609442"]
    
    atom_types = {'Li': 1, 'H': 2, 'Rh': 3}

    atom_type_masses = {'Li': 6.941, 'H': 1.00794, 'Rh': 102.9055}

    log_file = "m59948_auto.log"

    lammps_header = ["units real",
		     "atom_style full",
		     "boundary p p p"]

    amendments = ["set type 1 charge 0.45",
		  "set type 2 charge -0.21",
		  "set type 3 charge -0.53"]

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
######   Automatically generated lammps method m60450_auto   #####
@subdir_calc
def m60450_auto(atoms: ase.Atoms):

    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
	       "pair_modify tail yes mix arithmetic",
	       "dielectric 1.0",
	       "special_bonds lj/coul 0.0 0.0 1.0",
	       "kspace_style ewald 1.0e-5",
	       "pair_coeff 1 1 0.025000 2.183593",
	       "pair_coeff 2 2 0.085000 2.445517",
	       "pair_coeff 3 3 0.305000 3.694557"]
    
    atom_types = {'Li': 1, 'Be': 2, 'P': 3}

    atom_type_masses = {'Li': 6.941, 'Be': 9.012182, 'P': 30.973762}

    log_file = "m60450_auto.log"

    lammps_header = ["units real",
		     "atom_style full",
		     "boundary p p p"]

    amendments = ["set type 1 charge 0.41",
		  "set type 2 charge 0.55",
		  "set type 3 charge -0.96"]

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
######   Automatically generated lammps method m61111_auto   #####
@subdir_calc
def m61111_auto(atoms: ase.Atoms):

    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
	       "dielectric 1.0",
	       "special_bonds lj/coul 0.0 0.0 1.0",
	       "pair_modify tail yes mix arithmetic",
	       "kspace_style ewald 1.0e-5",
	       "pair_coeff 1 1 0.030000 2.657551",
	       "pair_coeff 2 2 0.025000 2.183593",
	       "pair_coeff 3 3 0.105000 3.430851",
	       "pair_coeff 4 4 0.060000 3.118146"]
    
    atom_types = {'Na': 1, 'Li': 2, 'C': 3, 'O': 4}

    atom_type_masses = {'Na': 22.98976928, 'Li': 6.941, 'C': 12.0107, 'O': 15.9994}

    log_file = "m61111_auto.log"

    lammps_header = ["units real",
		     "atom_style full",
		     "boundary p p p"]

    amendments = ["set type 1 charge 0.81",
		  "set type 2 charge 0.73",
		  "set type 3 charge 1.3",
		  "set type 4 charge -0.95"]

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
######   Automatically generated lammps method m61237_auto   #####
@subdir_calc
def m61237_auto(atoms: ase.Atoms):

    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
	       "pair_modify tail yes mix arithmetic",
	       "dielectric 1.0",
	       "special_bonds lj/coul 0.0 0.0 1.0",
	       "kspace_style ewald 1.0e-5",
	       "pair_coeff 1 1 0.025000 2.183593",
	       "pair_coeff 2 2 0.124000 2.461553",
	       "pair_coeff 3 3 0.305000 3.694557",
	       "pair_coeff 4 4 0.274000 3.594776"]
    
    atom_types = {'Li': 1, 'Zn': 2, 'P': 3, 'S': 4}

    atom_type_masses = {'Li': 6.941, 'Zn': 65.38, 'P': 30.973762, 'S': 32.065}

    log_file = "m61237_auto.log"

    lammps_header = ["units real",
		     "atom_style full",
		     "boundary p p p"]

    amendments = ["set type 1 charge 0.65",
		  "set type 2 charge 0.77",
		  "set type 3 charge 0.36",
		  "set type 4 charge -0.45"]

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
######   Automatically generated lammps method m61329_auto   #####
@subdir_calc
def m61329_auto(atoms: ase.Atoms):

    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
	       "dielectric 1.0",
	       "special_bonds lj/coul 0.0 0.0 1.0",
	       "pair_modify tail yes mix arithmetic",
	       "kspace_style ewald 1.0e-5",
	       "pair_coeff 1 1 0.045000 4.024190",
	       "pair_coeff 2 2 0.025000 2.183593",
	       "pair_coeff 3 3 0.044000 2.571134",
	       "pair_coeff 4 4 0.060000 3.118146"]
    
    atom_types = {'Cs': 1, 'Li': 2, 'H': 3, 'O': 4}

    atom_type_masses = {'Cs': 132.9054519, 'Li': 6.941, 'H': 1.00794, 'O': 15.9994}

    log_file = "m61329_auto.log"

    lammps_header = ["units real",
		     "atom_style full",
		     "boundary p p p"]

    amendments = ["set type 1 charge 0.71",
		  "set type 2 charge 0.64",
		  "set type 3 charge 0.21",
		  "set type 4 charge -0.87"]

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
######   Automatically generated lammps method m63257_auto   #####
@subdir_calc
def m63257_auto(atoms: ase.Atoms):

    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
	       "special_bonds lj/coul 0.0 0.0 1.0",
	       "pair_modify tail yes mix arithmetic",
	       "dielectric 1.0",
	       "kspace_style ewald 1.0e-5",
	       "pair_coeff 1 1 0.025000 2.183593",
	       "pair_coeff 2 2 0.238000 3.028165",
	       "pair_coeff 3 3 0.505000 4.008153",
	       "pair_coeff 4 4 0.069000 3.260689"]
    
    atom_types = {'Li': 1, 'Ca': 2, 'Al': 3, 'N': 4}

    atom_type_masses = {'Li': 6.941, 'Ca': 40.078, 'Al': 26.9815386, 'N': 14.0067}

    log_file = "m63257_auto.log"

    lammps_header = ["units real",
		     "atom_style full",
		     "boundary p p p"]

    amendments = ["set type 1 charge 0.71",
		  "set type 2 charge 1.2",
		  "set type 3 charge 1.48",
		  "set type 4 charge -1.69"]

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
######   Automatically generated lammps method m63279_auto   #####
@subdir_calc
def m63279_auto(atoms: ase.Atoms):

    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
	       "dielectric 1.0",
	       "pair_modify tail yes mix arithmetic",
	       "special_bonds lj/coul 0.0 0.0 1.0",
	       "kspace_style ewald 1.0e-5",
	       "pair_coeff 1 1 0.025000 2.183593",
	       "pair_coeff 2 2 0.111000 2.691405",
	       "pair_coeff 3 3 0.016000 2.800986",
	       "pair_coeff 4 4 0.060000 3.118146"]
    
    atom_types = {'Li': 1, 'Mg': 2, 'V': 3, 'O': 4}

    atom_type_masses = {'Li': 6.941, 'Mg': 24.305, 'V': 50.9415, 'O': 15.9994}

    log_file = "m63279_auto.log"

    lammps_header = ["units real",
		     "atom_style full",
		     "boundary p p p"]

    amendments = ["set type 1 charge 0.88",
		  "set type 2 charge 1.84",
		  "set type 3 charge 1.93",
		  "set type 4 charge -1.16"]

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
######   Automatically generated lammps method m63600_auto   #####
@subdir_calc
def m63600_auto(atoms: ase.Atoms):

    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
	       "dielectric 1.0",
	       "pair_modify tail yes mix arithmetic",
	       "special_bonds lj/coul 0.0 0.0 1.0",
	       "kspace_style ewald 1.0e-5",
	       "pair_coeff 1 1 0.025000 2.183593",
	       "pair_coeff 2 2 0.379000 3.813047",
	       "pair_coeff 3 3 0.398000 3.982317"]
    
    atom_types = {'Li': 1, 'Ge': 2, 'Te': 3}

    atom_type_masses = {'Li': 6.941, 'Ge': 72.64, 'Te': 127.6}

    log_file = "m63600_auto.log"

    lammps_header = ["units real",
		     "atom_style full",
		     "boundary p p p"]

    amendments = ["set type 1 charge 0.56",
		  "set type 2 charge -0.01",
		  "set type 3 charge -0.28"]

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
######   Automatically generated lammps method m63663_auto   #####
@subdir_calc
def m63663_auto(atoms: ase.Atoms):

    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
	       "dielectric 1.0",
	       "special_bonds lj/coul 0.0 0.0 1.0",
	       "pair_modify tail yes mix arithmetic",
	       "kspace_style ewald 1.0e-5",
	       "pair_coeff 1 1 0.025000 2.183593",
	       "pair_coeff 2 2 0.228000 2.537280",
	       "pair_coeff 3 3 0.014000 2.558661",
	       "pair_coeff 4 4 0.050000 2.996983"]
    
    atom_types = {'Li': 1, 'Cd': 2, 'Co': 3, 'F': 4}

    atom_type_masses = {'Li': 6.941, 'Cd': 112.411, 'Co': 58.933195, 'F': 18.9984032}

    log_file = "m63663_auto.log"

    lammps_header = ["units real",
		     "atom_style full",
		     "boundary p p p"]

    amendments = ["set type 1 charge 0.76",
		  "set type 2 charge 1.49",
		  "set type 3 charge 0.86",
		  "set type 4 charge -0.52"]

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
######   Automatically generated lammps method m63922_auto   #####
@subdir_calc
def m63922_auto(atoms: ase.Atoms):

    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
	       "dielectric 1.0",
	       "special_bonds lj/coul 0.0 0.0 1.0",
	       "pair_modify tail yes mix arithmetic",
	       "kspace_style ewald 1.0e-5",
	       "pair_coeff 1 1 0.025000 2.183593",
	       "pair_coeff 2 2 0.567000 3.912827",
	       "pair_coeff 3 3 0.039000 2.933729"]
    
    atom_types = {'Li': 1, 'Sn': 2, 'Au': 3}

    atom_type_masses = {'Li': 6.941, 'Sn': 118.71, 'Au': 196.966569}

    log_file = "m63922_auto.log"

    lammps_header = ["units real",
		     "atom_style full",
		     "boundary p p p"]

    amendments = ["set type 1 charge 0.81",
		  "set type 2 charge -1.11",
		  "set type 3 charge 0.3"]

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
######   Automatically generated lammps method m64383_auto   #####
@subdir_calc
def m64383_auto(atoms: ase.Atoms):

    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
	       "dielectric 1.0",
	       "pair_modify tail yes mix arithmetic",
	       "special_bonds lj/coul 0.0 0.0 1.0",
	       "kspace_style ewald 1.0e-5",
	       "pair_coeff 1 1 0.235000 3.243762",
	       "pair_coeff 2 2 0.025000 2.183593",
	       "pair_coeff 3 3 0.305000 3.694557"]
    
    atom_types = {'Sr': 1, 'Li': 2, 'P': 3}

    atom_type_masses = {'Sr': 87.62, 'Li': 6.941, 'P': 30.973762}

    log_file = "m64383_auto.log"

    lammps_header = ["units real",
		     "atom_style full",
		     "boundary p p p"]

    amendments = ["set type 1 charge 0.98",
		  "set type 1 charge 0.98",
		  "set type 2 charge 0.48",
		  "set type 2 charge 0.48",
		  "set type 3 charge -1.44",
		  "set type 3 charge -1.44"]

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
######   Automatically generated lammps method m64396_auto   #####
@subdir_calc
def m64396_auto(atoms: ase.Atoms):

    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
	       "dielectric 1.0",
	       "special_bonds lj/coul 0.0 0.0 1.0",
	       "pair_modify tail yes mix arithmetic",
	       "kspace_style ewald 1.0e-5",
	       "pair_coeff 1 1 0.364000 3.298998",
	       "pair_coeff 2 2 0.518000 3.893227",
	       "pair_coeff 3 3 0.025000 2.183593"]
    
    atom_types = {'Ba': 1, 'Bi': 2, 'Li': 3}

    atom_type_masses = {'Ba': 137.327, 'Bi': 208.9804, 'Li': 6.941}

    log_file = "m64396_auto.log"

    lammps_header = ["units real",
		     "atom_style full",
		     "boundary p p p"]

    amendments = ["set type 1 charge 1.29",
		  "set type 2 charge -1.66",
		  "set type 3 charge 0.69"]

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
######   Automatically generated lammps method m64754_auto   #####
@subdir_calc
def m64754_auto(atoms: ase.Atoms):

    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
	       "dielectric 1.0",
	       "pair_modify tail yes mix arithmetic",
	       "special_bonds lj/coul 0.0 0.0 1.0",
	       "kspace_style ewald 1.0e-5",
	       "pair_coeff 1 1 0.025000 2.183593",
	       "pair_coeff 2 2 0.599000 3.976081",
	       "pair_coeff 3 3 0.073000 2.530152"]
    
    atom_types = {'Li': 1, 'In': 2, 'Ir': 3}

    atom_type_masses = {'Li': 6.941, 'In': 114.818, 'Ir': 192.217}

    log_file = "m64754_auto.log"

    lammps_header = ["units real",
		     "atom_style full",
		     "boundary p p p"]

    amendments = ["set type 1 charge 0.74",
		  "set type 2 charge -2.91",
		  "set type 3 charge 1.43"]

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
######   Automatically generated lammps method m64867_auto   #####
@subdir_calc
def m64867_auto(atoms: ase.Atoms):

    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
	       "pair_modify tail yes mix arithmetic",
	       "special_bonds lj/coul 0.0 0.0 1.0",
	       "dielectric 1.0",
	       "kspace_style ewald 1.0e-5",
	       "pair_coeff 1 1 0.025000 2.183593",
	       "pair_coeff 2 2 0.415000 3.904809",
	       "pair_coeff 3 3 0.060000 3.118146"]
    
    atom_types = {'Li': 1, 'Ga': 2, 'O': 3}

    atom_type_masses = {'Li': 6.941, 'Ga': 69.723, 'O': 15.9994}

    log_file = "m64867_auto.log"

    lammps_header = ["units real",
		     "atom_style full",
		     "boundary p p p"]

    amendments = ["set type 1 charge 0.84",
		  "set type 2 charge 1.32",
		  "set type 3 charge -1.08"]

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
######   Automatically generated lammps method m65350_auto   #####
@subdir_calc
def m65350_auto(atoms: ase.Atoms):

    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
	       "dielectric 1.0",
	       "special_bonds lj/coul 0.0 0.0 1.0",
	       "pair_modify tail yes mix arithmetic",
	       "kspace_style ewald 1.0e-5",
	       "pair_coeff 1 1 0.040000 3.665157",
	       "pair_coeff 2 2 0.025000 2.183593",
	       "pair_coeff 3 3 0.072000 2.980056",
	       "pair_coeff 4 4 0.227000 3.516377"]
    
    atom_types = {'Rb': 1, 'Li': 2, 'Y': 3, 'Cl': 4}

    atom_type_masses = {'Rb': 85.4678, 'Li': 6.941, 'Y': 88.90585, 'Cl': 35.453}

    log_file = "m65350_auto.log"

    lammps_header = ["units real",
		     "atom_style full",
		     "boundary p p p"]

    amendments = ["set type 1 charge 0.78",
		  "set type 2 charge 1.04",
		  "set type 3 charge 1.57",
		  "set type 4 charge -0.7"]

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
######   Automatically generated lammps method m65898_auto   #####
@subdir_calc
def m65898_auto(atoms: ase.Atoms):

    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
	       "pair_modify tail yes mix arithmetic",
	       "special_bonds lj/coul 0.0 0.0 1.0",
	       "dielectric 1.0",
	       "kspace_style ewald 1.0e-5",
	       "pair_coeff 1 1 0.045000 4.024190",
	       "pair_coeff 2 2 0.025000 2.183593",
	       "pair_coeff 3 3 0.017000 2.828603",
	       "pair_coeff 4 4 0.060000 3.118146"]
    
    atom_types = {'Cs': 1, 'Li': 2, 'Ti': 3, 'O': 4}

    atom_type_masses = {'Cs': 132.9054519, 'Li': 6.941, 'Ti': 47.867, 'O': 15.9994}

    log_file = "m65898_auto.log"

    lammps_header = ["units real",
		     "atom_style full",
		     "boundary p p p"]

    amendments = ["set type 1 charge 0.82",
		  "set type 2 charge 0.84",
		  "set type 3 charge 1.66",
		  "set type 4 charge -1.25"]

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
######   Automatically generated lammps method m66685_auto   #####
@subdir_calc
def m66685_auto(atoms: ase.Atoms):

    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
	       "pair_modify tail yes mix arithmetic",
	       "dielectric 1.0",
	       "special_bonds lj/coul 0.0 0.0 1.0",
	       "kspace_style ewald 1.0e-5",
	       "pair_coeff 1 1 0.025000 2.183593",
	       "pair_coeff 2 2 0.041000 3.242871",
	       "pair_coeff 3 3 0.050000 2.996983"]
    
    atom_types = {'Li': 1, 'Lu': 2, 'F': 3}

    atom_type_masses = {'Li': 6.941, 'Lu': 174.9668, 'F': 18.9984032}

    log_file = "m66685_auto.log"

    lammps_header = ["units real",
		     "atom_style full",
		     "boundary p p p"]

    amendments = ["set type 1 charge 1.03",
		  "set type 2 charge 2.31",
		  "set type 3 charge -0.83"]

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
######   Automatically generated lammps method m66876_auto   #####
@subdir_calc
def m66876_auto(atoms: ase.Atoms):

    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
	       "pair_modify tail yes mix arithmetic",
	       "special_bonds lj/coul 0.0 0.0 1.0",
	       "dielectric 1.0",
	       "kspace_style ewald 1.0e-5",
	       "pair_coeff 1 1 0.035000 3.396106",
	       "pair_coeff 2 2 0.025000 2.183593",
	       "pair_coeff 3 3 0.056000 2.719023",
	       "pair_coeff 4 4 0.060000 3.118146"]
    
    atom_types = {'K': 1, 'Li': 2, 'Mo': 3, 'O': 4}

    atom_type_masses = {'K': 39.0983, 'Li': 6.941, 'Mo': 95.96, 'O': 15.9994}

    log_file = "m66876_auto.log"

    lammps_header = ["units real",
		     "atom_style full",
		     "boundary p p p"]

    amendments = ["set type 1 charge 0.89",
		  "set type 2 charge 1.28",
		  "set type 3 charge 2.34",
		  "set type 4 charge -1.13"]

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
######   Automatically generated lammps method m66985_auto   #####
@subdir_calc
def m66985_auto(atoms: ase.Atoms):

    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
	       "pair_modify tail yes mix arithmetic",
	       "special_bonds lj/coul 0.0 0.0 1.0",
	       "dielectric 1.0",
	       "kspace_style ewald 1.0e-5",
	       "pair_coeff 1 1 0.025000 2.183593",
	       "pair_coeff 2 2 0.072000 2.980056",
	       "pair_coeff 3 3 0.050000 2.996983"]
    
    atom_types = {'Li': 1, 'Y': 2, 'F': 3}

    atom_type_masses = {'Li': 6.941, 'Y': 88.90585, 'F': 18.9984032}

    log_file = "m66985_auto.log"

    lammps_header = ["units real",
		     "atom_style full",
		     "boundary p p p"]

    amendments = ["set type 1 charge 1.06",
		  "set type 2 charge 2.44",
		  "set type 3 charge -0.88"]

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
######   Automatically generated lammps method m67415_auto   #####
@subdir_calc
def m67415_auto(atoms: ase.Atoms):

    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
	       "dielectric 1.0",
	       "special_bonds lj/coul 0.0 0.0 1.0",
	       "pair_modify tail yes mix arithmetic",
	       "kspace_style ewald 1.0e-5",
	       "pair_coeff 1 1 0.045000 4.024190",
	       "pair_coeff 2 2 0.025000 2.183593",
	       "pair_coeff 3 3 0.044000 2.571134",
	       "pair_coeff 4 4 0.069000 3.260689"]
    
    atom_types = {'Cs': 1, 'Li': 2, 'H': 3, 'N': 4}

    atom_type_masses = {'Cs': 132.9054519, 'Li': 6.941, 'H': 1.00794, 'N': 14.0067}

    log_file = "m67415_auto.log"

    lammps_header = ["units real",
		     "atom_style full",
		     "boundary p p p"]

    amendments = ["set type 1 charge 0.59",
		  "set type 2 charge 0.83",
		  "set type 3 charge 0.32",
		  "set type 4 charge -1.34"]

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
######   Automatically generated lammps method m67418_auto   #####
@subdir_calc
def m67418_auto(atoms: ase.Atoms):

    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
	       "pair_modify tail yes mix arithmetic",
	       "dielectric 1.0",
	       "special_bonds lj/coul 0.0 0.0 1.0",
	       "kspace_style ewald 1.0e-5",
	       "pair_coeff 1 1 0.025000 2.183593",
	       "pair_coeff 2 2 0.238000 3.028165",
	       "pair_coeff 3 3 0.402000 3.826410",
	       "pair_coeff 4 4 0.069000 3.260689"]
    
    atom_types = {'Li': 1, 'Ca': 2, 'Si': 3, 'N': 4}

    atom_type_masses = {'Li': 6.941, 'Ca': 40.078, 'Si': 28.0855, 'N': 14.0067}

    log_file = "m67418_auto.log"

    lammps_header = ["units real",
		     "atom_style full",
		     "boundary p p p"]

    amendments = ["set type 1 charge 0.83",
		  "set type 2 charge 1.33",
		  "set type 3 charge 1.86",
		  "set type 4 charge -1.84"]

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
######   Automatically generated lammps method m67779_auto   #####
@subdir_calc
def m67779_auto(atoms: ase.Atoms):

    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
	       "pair_modify tail yes mix arithmetic",
	       "special_bonds lj/coul 0.0 0.0 1.0",
	       "dielectric 1.0",
	       "kspace_style ewald 1.0e-5",
	       "pair_coeff 1 1 0.025000 2.183593",
	       "pair_coeff 2 2 0.238000 3.028165",
	       "pair_coeff 3 3 0.505000 4.008153",
	       "pair_coeff 4 4 0.050000 2.996983"]
    
    atom_types = {'Li': 1, 'Ca': 2, 'Al': 3, 'F': 4}

    atom_type_masses = {'Li': 6.941, 'Ca': 40.078, 'Al': 26.9815386, 'F': 18.9984032}

    log_file = "m67779_auto.log"

    lammps_header = ["units real",
		     "atom_style full",
		     "boundary p p p"]

    amendments = ["set type 1 charge 0.77",
		  "set type 2 charge 1.63",
		  "set type 3 charge 0.89",
		  "set type 4 charge -0.55"]

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
######   Automatically generated lammps method m67785_auto   #####
@subdir_calc
def m67785_auto(atoms: ase.Atoms):

    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
	       "dielectric 1.0",
	       "special_bonds lj/coul 0.0 0.0 1.0",
	       "pair_modify tail yes mix arithmetic",
	       "kspace_style ewald 1.0e-5",
	       "pair_coeff 1 1 0.025000 2.183593",
	       "pair_coeff 2 2 0.305000 3.694557",
	       "pair_coeff 3 3 0.274000 3.594776"]
    
    atom_types = {'Li': 1, 'P': 2, 'S': 3}

    atom_type_masses = {'Li': 6.941, 'P': 30.973762, 'S': 32.065}

    log_file = "m67785_auto.log"

    lammps_header = ["units real",
		     "atom_style full",
		     "boundary p p p"]

    amendments = ["set type 1 charge 0.79",
		  "set type 2 charge 0.95",
		  "set type 3 charge -0.83"]

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
######   Automatically generated lammps method m67806_auto   #####
@subdir_calc
def m67806_auto(atoms: ase.Atoms):

    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
	       "special_bonds lj/coul 0.0 0.0 1.0",
	       "pair_modify tail yes mix arithmetic",
	       "dielectric 1.0",
	       "kspace_style ewald 1.0e-5",
	       "pair_coeff 1 1 0.025000 2.183593",
	       "pair_coeff 2 2 0.059000 2.819694",
	       "pair_coeff 3 3 0.067000 2.734168",
	       "pair_coeff 4 4 0.060000 3.118146"]
    
    atom_types = {'Li': 1, 'Nb': 2, 'W': 3, 'O': 4}

    atom_type_masses = {'Li': 6.941, 'Nb': 92.90638, 'W': 183.84, 'O': 15.9994}

    log_file = "m67806_auto.log"

    lammps_header = ["units real",
		     "atom_style full",
		     "boundary p p p"]

    amendments = ["set type 1 charge 0.82",
		  "set type 2 charge 2.58",
		  "set type 3 charge 1.86",
		  "set type 4 charge -0.88"]

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
######   Automatically generated lammps method m67999_auto   #####
@subdir_calc
def m67999_auto(atoms: ase.Atoms):

    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
	       "special_bonds lj/coul 0.0 0.0 1.0",
	       "pair_modify tail yes mix arithmetic",
	       "dielectric 1.0",
	       "kspace_style ewald 1.0e-5",
	       "pair_coeff 1 1 0.035000 3.396106",
	       "pair_coeff 2 2 0.025000 2.183593",
	       "pair_coeff 3 3 0.067000 2.734168",
	       "pair_coeff 4 4 0.060000 3.118146"]
    
    atom_types = {'K': 1, 'Li': 2, 'W': 3, 'O': 4}

    atom_type_masses = {'K': 39.0983, 'Li': 6.941, 'W': 183.84, 'O': 15.9994}

    log_file = "m67999_auto.log"

    lammps_header = ["units real",
		     "atom_style full",
		     "boundary p p p"]

    amendments = ["set type 1 charge 0.88",
		  "set type 2 charge 1.29",
		  "set type 3 charge 2.41",
		  "set type 4 charge -1.15"]

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
######   Automatically generated lammps method m68042_auto   #####
@subdir_calc
def m68042_auto(atoms: ase.Atoms):

    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
	       "special_bonds lj/coul 0.0 0.0 1.0",
	       "dielectric 1.0",
	       "pair_modify tail yes mix arithmetic",
	       "kspace_style ewald 1.0e-5",
	       "pair_coeff 1 1 0.025000 2.183593",
	       "pair_coeff 2 2 0.111000 2.691405",
	       "pair_coeff 3 3 0.069000 3.260689"]
    
    atom_types = {'Li': 1, 'Mg': 2, 'N': 3}

    atom_type_masses = {'Li': 6.941, 'Mg': 24.305, 'N': 14.0067}

    log_file = "m68042_auto.log"

    lammps_header = ["units real",
		     "atom_style full",
		     "boundary p p p"]

    amendments = ["set type 1 charge 1.27",
		  "set type 2 charge 2.07",
		  "set type 3 charge -3.34"]

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
######   Automatically generated lammps method m68103_auto   #####
@subdir_calc
def m68103_auto(atoms: ase.Atoms):

    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
	       "special_bonds lj/coul 0.0 0.0 1.0",
	       "dielectric 1.0",
	       "pair_modify tail yes mix arithmetic",
	       "kspace_style ewald 1.0e-5",
	       "pair_coeff 1 1 0.235000 3.243762",
	       "pair_coeff 2 2 0.379000 3.813047",
	       "pair_coeff 3 3 0.069000 3.260689",
	       "pair_coeff 4 4 0.025000 2.183593"]
    
    atom_types = {'Sr': 1, 'Ge': 2, 'N': 3, 'Li': 4}

    atom_type_masses = {'Sr': 87.62, 'Ge': 72.64, 'N': 14.0067, 'Li': 6.941}

    log_file = "m68103_auto.log"

    lammps_header = ["units real",
		     "atom_style full",
		     "boundary p p p"]

    amendments = ["set type 1 charge 1.47",
		  "set type 2 charge 1.66",
		  "set type 3 charge -1.89",
		  "set type 4 charge 0.9"]

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
######   Automatically generated lammps method m68120_auto   #####
@subdir_calc
def m68120_auto(atoms: ase.Atoms):

    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
	       "special_bonds lj/coul 0.0 0.0 1.0",
	       "pair_modify tail yes mix arithmetic",
	       "dielectric 1.0",
	       "kspace_style ewald 1.0e-5",
	       "pair_coeff 1 1 0.025000 2.183593",
	       "pair_coeff 2 2 0.056000 2.719023",
	       "pair_coeff 3 3 0.274000 3.594776"]
    
    atom_types = {'Li': 1, 'Mo': 2, 'S': 3}

    atom_type_masses = {'Li': 6.941, 'Mo': 95.96, 'S': 32.065}

    log_file = "m68120_auto.log"

    lammps_header = ["units real",
		     "atom_style full",
		     "boundary p p p"]

    amendments = ["set type 1 charge 0.62",
		  "set type 2 charge 0.7",
		  "set type 3 charge -0.66"]

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
######   Automatically generated lammps method m68242_auto   #####
@subdir_calc
def m68242_auto(atoms: ase.Atoms):

    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
	       "pair_modify tail yes mix arithmetic",
	       "dielectric 1.0",
	       "special_bonds lj/coul 0.0 0.0 1.0",
	       "kspace_style ewald 1.0e-5",
	       "pair_coeff 1 1 0.025000 2.183593",
	       "pair_coeff 2 2 0.019000 2.935511",
	       "pair_coeff 3 3 0.402000 3.826410",
	       "pair_coeff 4 4 0.060000 3.118146"]
    
    atom_types = {'Li': 1, 'Sc': 2, 'Si': 3, 'O': 4}

    atom_type_masses = {'Li': 6.941, 'Sc': 44.955912, 'Si': 28.0855, 'O': 15.9994}

    log_file = "m68242_auto.log"

    lammps_header = ["units real",
		     "atom_style full",
		     "boundary p p p"]

    amendments = ["set type 1 charge 0.84",
		  "set type 2 charge 1.62",
		  "set type 3 charge 1.37",
		  "set type 4 charge -0.87"]

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
######   Automatically generated lammps method m68674_auto   #####
@subdir_calc
def m68674_auto(atoms: ase.Atoms):

    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
	       "pair_modify tail yes mix arithmetic",
	       "dielectric 1.0",
	       "special_bonds lj/coul 0.0 0.0 1.0",
	       "kspace_style ewald 1.0e-5",
	       "pair_coeff 1 1 0.025000 2.183593",
	       "pair_coeff 2 2 0.017000 3.137745",
	       "pair_coeff 3 3 0.449000 3.937772",
	       "pair_coeff 4 4 0.060000 3.118146"]
    
    atom_types = {'Li': 1, 'La': 2, 'Sb': 3, 'O': 4}

    atom_type_masses = {'Li': 6.941, 'La': 138.90547, 'Sb': 121.76, 'O': 15.9994}

    log_file = "m68674_auto.log"

    lammps_header = ["units real",
		     "atom_style full",
		     "boundary p p p"]

    amendments = ["set type 1 charge 1.21",
		  "set type 2 charge 2.38",
		  "set type 3 charge 2.45",
		  "set type 4 charge -1.4"]

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
######   Automatically generated lammps method m68894_auto   #####
@subdir_calc
def m68894_auto(atoms: ase.Atoms):

    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
	       "special_bonds lj/coul 0.0 0.0 1.0",
	       "dielectric 1.0",
	       "pair_modify tail yes mix arithmetic",
	       "kspace_style ewald 1.0e-5",
	       "pair_coeff 1 1 0.025000 2.183593",
	       "pair_coeff 2 2 0.449000 3.937772",
	       "pair_coeff 3 3 0.060000 3.118146"]
    
    atom_types = {'Li': 1, 'Sb': 2, 'O': 3}

    atom_type_masses = {'Li': 6.941, 'Sb': 121.76, 'O': 15.9994}

    log_file = "m68894_auto.log"

    lammps_header = ["units real",
		     "atom_style full",
		     "boundary p p p"]

    amendments = ["set type 1 charge 0.95",
		  "set type 2 charge 2.41",
		  "set type 3 charge -1.12"]

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
######   Automatically generated lammps method m68960_auto   #####
@subdir_calc
def m68960_auto(atoms: ase.Atoms):

    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
	       "dielectric 1.0",
	       "special_bonds lj/coul 0.0 0.0 1.0",
	       "pair_modify tail yes mix arithmetic",
	       "kspace_style ewald 1.0e-5",
	       "pair_coeff 1 1 0.030000 2.657551",
	       "pair_coeff 2 2 0.025000 2.183593",
	       "pair_coeff 3 3 0.505000 4.008153",
	       "pair_coeff 4 4 0.044000 2.571134"]
    
    atom_types = {'Na': 1, 'Li': 2, 'Al': 3, 'H': 4}

    atom_type_masses = {'Na': 22.98976928, 'Li': 6.941, 'Al': 26.9815386, 'H': 1.00794}

    log_file = "m68960_auto.log"

    lammps_header = ["units real",
		     "atom_style full",
		     "boundary p p p"]

    amendments = ["set type 1 charge 0.89",
		  "set type 1 charge 0.89",
		  "set type 2 charge 0.85",
		  "set type 2 charge 0.85",
		  "set type 3 charge 0.39",
		  "set type 3 charge 0.39",
		  "set type 4 charge -0.5",
		  "set type 4 charge -0.5"]

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
######   Automatically generated lammps method m69116_auto   #####
@subdir_calc
def m69116_auto(atoms: ase.Atoms):

    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
	       "pair_modify tail yes mix arithmetic",
	       "dielectric 1.0",
	       "special_bonds lj/coul 0.0 0.0 1.0",
	       "kspace_style ewald 1.0e-5",
	       "pair_coeff 1 1 0.025000 2.183593",
	       "pair_coeff 2 2 0.072000 2.980056",
	       "pair_coeff 3 3 0.291000 3.746229"]
    
    atom_types = {'Li': 1, 'Y': 2, 'Se': 3}

    atom_type_masses = {'Li': 6.941, 'Y': 88.90585, 'Se': 78.96}

    log_file = "m69116_auto.log"

    lammps_header = ["units real",
		     "atom_style full",
		     "boundary p p p"]

    amendments = ["set type 1 charge 0.79",
		  "set type 1 charge 0.79",
		  "set type 2 charge 0.92",
		  "set type 2 charge 0.92",
		  "set type 3 charge -0.85",
		  "set type 3 charge -0.85"]

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
######   Automatically generated lammps method m69642_auto   #####
@subdir_calc
def m69642_auto(atoms: ase.Atoms):

    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
	       "pair_modify tail yes mix arithmetic",
	       "dielectric 1.0",
	       "special_bonds lj/coul 0.0 0.0 1.0",
	       "kspace_style ewald 1.0e-5",
	       "pair_coeff 1 1 0.025000 2.183593",
	       "pair_coeff 2 2 0.599000 3.976081",
	       "pair_coeff 3 3 0.402000 3.826410",
	       "pair_coeff 4 4 0.060000 3.118146"]
    
    atom_types = {'Li': 1, 'In': 2, 'Si': 3, 'O': 4}

    atom_type_masses = {'Li': 6.941, 'In': 114.818, 'Si': 28.0855, 'O': 15.9994}

    log_file = "m69642_auto.log"

    lammps_header = ["units real",
		     "atom_style full",
		     "boundary p p p"]

    amendments = ["set type 1 charge 0.86",
		  "set type 2 charge 1.37",
		  "set type 3 charge 1.53",
		  "set type 4 charge -0.88"]

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
######   Automatically generated lammps method m70003_auto   #####
@subdir_calc
def m70003_auto(atoms: ase.Atoms):

    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
	       "special_bonds lj/coul 0.0 0.0 1.0",
	       "dielectric 1.0",
	       "pair_modify tail yes mix arithmetic",
	       "kspace_style ewald 1.0e-5",
	       "pair_coeff 1 1 0.025000 2.183593",
	       "pair_coeff 2 2 0.105000 3.430851",
	       "pair_coeff 3 3 0.274000 3.594776",
	       "pair_coeff 4 4 0.069000 3.260689"]
    
    atom_types = {'Li': 1, 'C': 2, 'S': 3, 'N': 4}

    atom_type_masses = {'Li': 6.941, 'C': 12.0107, 'S': 32.065, 'N': 14.0067}

    log_file = "m70003_auto.log"

    lammps_header = ["units real",
		     "atom_style full",
		     "boundary p p p"]

    amendments = ["set type 1 charge 0.84",
		  "set type 2 charge 0.54",
		  "set type 3 charge -0.4",
		  "set type 4 charge -0.98"]

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
######   Automatically generated lammps method m70442_auto   #####
@subdir_calc
def m70442_auto(atoms: ase.Atoms):

    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
	       "dielectric 1.0",
	       "pair_modify tail yes mix arithmetic",
	       "special_bonds lj/coul 0.0 0.0 1.0",
	       "kspace_style ewald 1.0e-5",
	       "pair_coeff 1 1 0.235000 3.243762",
	       "pair_coeff 2 2 0.066000 2.631715",
	       "pair_coeff 3 3 0.060000 3.118146",
	       "pair_coeff 4 4 0.025000 2.183593"]
    
    atom_types = {'Sr': 1, 'Re': 2, 'O': 3, 'Li': 4}

    atom_type_masses = {'Sr': 87.62, 'Re': 186.207, 'O': 15.9994, 'Li': 6.941}

    log_file = "m70442_auto.log"

    lammps_header = ["units real",
		     "atom_style full",
		     "boundary p p p"]

    amendments = ["set type 1 charge 1.69",
		  "set type 2 charge 3.26",
		  "set type 3 charge -1.45",
		  "set type 4 charge 2.07"]

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
######   Automatically generated lammps method m70681_auto   #####
@subdir_calc
def m70681_auto(atoms: ase.Atoms):

    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
	       "dielectric 1.0",
	       "special_bonds lj/coul 0.0 0.0 1.0",
	       "pair_modify tail yes mix arithmetic",
	       "kspace_style ewald 1.0e-5",
	       "pair_coeff 1 1 0.025000 2.183593",
	       "pair_coeff 2 2 0.044000 2.571134",
	       "pair_coeff 3 3 0.291000 3.746229",
	       "pair_coeff 4 4 0.060000 3.118146",
	       "pair_coeff 5 5 0.060000 3.118146"]
    
    atom_types = {'Li': 1, 'H': 2, 'Se': 3, 'O': 5}

    atom_type_masses = {'Li': 6.941, 'H': 1.00794, 'Se': 78.96, 'O': 15.9994}

    log_file = "m70681_auto.log"

    lammps_header = ["units real",
		     "atom_style full",
		     "boundary p p p"]

    amendments = ["set type 1 charge 0.87",
		  "set type 2 charge 0.41",
		  "set type 3 charge 0.53",
		  "set type 4 charge -0.6",
		  "set type 5 charge -0.6"]

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
######   Automatically generated lammps method m71368_auto   #####
@subdir_calc
def m71368_auto(atoms: ase.Atoms):

    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
	       "special_bonds lj/coul 0.0 0.0 1.0",
	       "pair_modify tail yes mix arithmetic",
	       "dielectric 1.0",
	       "kspace_style ewald 1.0e-5",
	       "pair_coeff 1 1 0.025000 2.183593",
	       "pair_coeff 2 2 0.072000 2.980056",
	       "pair_coeff 3 3 0.050000 2.996983"]
    
    atom_types = {'Li': 1, 'Y': 2, 'F': 3}

    atom_type_masses = {'Li': 6.941, 'Y': 88.90585, 'F': 18.9984032}

    log_file = "m71368_auto.log"

    lammps_header = ["units real",
		     "atom_style full",
		     "boundary p p p"]

    amendments = ["set type 1 charge 0.59",
		  "set type 2 charge 0.41",
		  "set type 3 charge -0.5"]

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
######   Automatically generated lammps method m71681_auto   #####
@subdir_calc
def m71681_auto(atoms: ase.Atoms):

    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
	       "special_bonds lj/coul 0.0 0.0 1.0",
	       "dielectric 1.0",
	       "pair_modify tail yes mix arithmetic",
	       "kspace_style ewald 1.0e-5",
	       "pair_coeff 1 1 0.025000 2.183593",
	       "pair_coeff 2 2 0.019000 2.935511",
	       "pair_coeff 3 3 0.067000 2.734168",
	       "pair_coeff 4 4 0.060000 3.118146"]
    
    atom_types = {'Li': 1, 'Sc': 2, 'W': 3, 'O': 4}

    atom_type_masses = {'Li': 6.941, 'Sc': 44.955912, 'W': 183.84, 'O': 15.9994}

    log_file = "m71681_auto.log"

    lammps_header = ["units real",
		     "atom_style full",
		     "boundary p p p"]

    amendments = ["set type 1 charge 1.02",
		  "set type 2 charge 2.16",
		  "set type 3 charge 2.38",
		  "set type 4 charge -0.99"]

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
######   Automatically generated lammps method m71925_auto   #####
@subdir_calc
def m71925_auto(atoms: ase.Atoms):

    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
	       "pair_modify tail yes mix arithmetic",
	       "special_bonds lj/coul 0.0 0.0 1.0",
	       "dielectric 1.0",
	       "kspace_style ewald 1.0e-5",
	       "pair_coeff 1 1 0.025000 2.183593",
	       "pair_coeff 2 2 0.041000 3.242871",
	       "pair_coeff 3 3 0.067000 2.734168",
	       "pair_coeff 4 4 0.060000 3.118146",
	       "pair_coeff 5 5 0.060000 3.118146"]
    
    atom_types = {'Li': 1, 'Lu': 2, 'W': 3, 'O': 5}

    atom_type_masses = {'Li': 6.941, 'Lu': 174.9668, 'W': 183.84, 'O': 15.9994}

    log_file = "m71925_auto.log"

    lammps_header = ["units real",
		     "atom_style full",
		     "boundary p p p"]

    amendments = ["set type 1 charge 1.05",
		  "set type 2 charge 2.02",
		  "set type 3 charge 2.28",
		  "set type 4 charge -0.95",
		  "set type 5 charge -0.95"]

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
######   Automatically generated lammps method m71990_auto   #####
@subdir_calc
def m71990_auto(atoms: ase.Atoms):

    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
	       "dielectric 1.0",
	       "pair_modify tail yes mix arithmetic",
	       "special_bonds lj/coul 0.0 0.0 1.0",
	       "kspace_style ewald 1.0e-5",
	       "pair_coeff 1 1 0.235000 3.243762",
	       "pair_coeff 2 2 0.069000 3.260689",
	       "pair_coeff 3 3 0.025000 2.183593",
	       "pair_coeff 4 4 0.044000 2.571134"]
    
    atom_types = {'Sr': 1, 'N': 2, 'Li': 3, 'H': 4}

    atom_type_masses = {'Sr': 87.62, 'N': 14.0067, 'Li': 6.941, 'H': 1.00794}

    log_file = "m71990_auto.log"

    lammps_header = ["units real",
		     "atom_style full",
		     "boundary p p p"]

    amendments = ["set type 1 charge 0.49",
		  "set type 2 charge -0.55",
		  "set type 3 charge 0.18",
		  "set type 4 charge -0.31"]

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
######   Automatically generated lammps method m72018_auto   #####
@subdir_calc
def m72018_auto(atoms: ase.Atoms):

    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
	       "special_bonds lj/coul 0.0 0.0 1.0",
	       "dielectric 1.0",
	       "pair_modify tail yes mix arithmetic",
	       "kspace_style ewald 1.0e-5",
	       "pair_coeff 1 1 0.045000 4.024190",
	       "pair_coeff 2 2 0.025000 2.183593",
	       "pair_coeff 3 3 0.072000 2.980056",
	       "pair_coeff 4 4 0.227000 3.516377"]
    
    atom_types = {'Cs': 1, 'Li': 2, 'Y': 3, 'Cl': 4}

    atom_type_masses = {'Cs': 132.9054519, 'Li': 6.941, 'Y': 88.90585, 'Cl': 35.453}

    log_file = "m72018_auto.log"

    lammps_header = ["units real",
		     "atom_style full",
		     "boundary p p p"]

    amendments = ["set type 1 charge 0.68",
		  "set type 2 charge 0.91",
		  "set type 3 charge 1.48",
		  "set type 4 charge -0.62"]

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
######   Automatically generated lammps method m72372_auto   #####
@subdir_calc
def m72372_auto(atoms: ase.Atoms):

    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
	       "dielectric 1.0",
	       "pair_modify tail yes mix arithmetic",
	       "special_bonds lj/coul 0.0 0.0 1.0",
	       "kspace_style ewald 1.0e-5",
	       "pair_coeff 1 1 0.025000 2.183593",
	       "pair_coeff 2 2 0.041000 3.242871",
	       "pair_coeff 3 3 0.379000 3.813047"]
    
    atom_types = {'Li': 1, 'Lu': 2, 'Ge': 3}

    atom_type_masses = {'Li': 6.941, 'Lu': 174.9668, 'Ge': 72.64}

    log_file = "m72372_auto.log"

    lammps_header = ["units real",
		     "atom_style full",
		     "boundary p p p"]

    amendments = ["set type 1 charge 0.77",
		  "set type 2 charge 1.35",
		  "set type 3 charge -2.12"]

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
######   Automatically generated lammps method m72631_auto   #####
@subdir_calc
def m72631_auto(atoms: ase.Atoms):

    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
	       "dielectric 1.0",
	       "pair_modify tail yes mix arithmetic",
	       "special_bonds lj/coul 0.0 0.0 1.0",
	       "kspace_style ewald 1.0e-5",
	       "pair_coeff 1 1 0.025000 2.183593",
	       "pair_coeff 2 2 0.402000 3.826410",
	       "pair_coeff 3 3 0.274000 3.594776",
	       "pair_coeff 4 4 0.274000 3.594776"]
    
    atom_types = {'Li': 1, 'Si': 2, 'S': 4}

    atom_type_masses = {'Li': 6.941, 'Si': 28.0855, 'S': 32.065}

    log_file = "m72631_auto.log"

    lammps_header = ["units real",
		     "atom_style full",
		     "boundary p p p"]

    amendments = ["set type 1 charge 0.67",
		  "set type 2 charge -0.58",
		  "set type 2 charge 0.38",
		  "set type 3 charge -0.58",
		  "set type 4 charge -0.58"]

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
######   Automatically generated lammps method m72652_auto   #####
@subdir_calc
def m72652_auto(atoms: ase.Atoms):

    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
	       "pair_modify tail yes mix arithmetic",
	       "special_bonds lj/coul 0.0 0.0 1.0",
	       "dielectric 1.0",
	       "kspace_style ewald 1.0e-5",
	       "pair_coeff 1 1 0.025000 2.183593",
	       "pair_coeff 2 2 0.379000 3.813047",
	       "pair_coeff 3 3 0.060000 3.118146"]
    
    atom_types = {'Li': 1, 'Ge': 2, 'O': 3}

    atom_type_masses = {'Li': 6.941, 'Ge': 72.64, 'O': 15.9994}

    log_file = "m72652_auto.log"

    lammps_header = ["units real",
		     "atom_style full",
		     "boundary p p p"]

    amendments = ["set type 1 charge 0.98",
		  "set type 2 charge 1.47",
		  "set type 3 charge -0.98"]

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
######   Automatically generated lammps method m73298_auto   #####
@subdir_calc
def m73298_auto(atoms: ase.Atoms):

    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
	       "dielectric 1.0",
	       "special_bonds lj/coul 0.0 0.0 1.0",
	       "pair_modify tail yes mix arithmetic",
	       "kspace_style ewald 1.0e-5",
	       "pair_coeff 1 1 0.045000 4.024190",
	       "pair_coeff 2 2 0.025000 2.183593",
	       "pair_coeff 3 3 0.251000 3.731975"]
    
    atom_types = {'Cs': 1, 'Li': 2, 'Br': 3}

    atom_type_masses = {'Cs': 132.9054519, 'Li': 6.941, 'Br': 79.904}

    log_file = "m73298_auto.log"

    lammps_header = ["units real",
		     "atom_style full",
		     "boundary p p p"]

    amendments = ["set type 1 charge 1.19",
		  "set type 1 charge 1.19",
		  "set type 2 charge 0.96",
		  "set type 2 charge 0.96",
		  "set type 3 charge -1.07",
		  "set type 3 charge -1.07"]

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
######   Automatically generated lammps method m73679_auto   #####
@subdir_calc
def m73679_auto(atoms: ase.Atoms):

    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
	       "special_bonds lj/coul 0.0 0.0 1.0",
	       "dielectric 1.0",
	       "pair_modify tail yes mix arithmetic",
	       "kspace_style ewald 1.0e-5",
	       "pair_coeff 1 1 0.025000 2.183593",
	       "pair_coeff 2 2 0.505000 4.008153",
	       "pair_coeff 3 3 0.053000 2.609442"]
    
    atom_types = {'Li': 1, 'Al': 2, 'Rh': 3}

    atom_type_masses = {'Li': 6.941, 'Al': 26.9815386, 'Rh': 102.9055}

    log_file = "m73679_auto.log"

    lammps_header = ["units real",
		     "atom_style full",
		     "boundary p p p"]

    amendments = ["set type 1 charge 1.23",
		  "set type 2 charge -4.11",
		  "set type 3 charge 0.82"]

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
######   Automatically generated lammps method m73766_auto   #####
@subdir_calc
def m73766_auto(atoms: ase.Atoms):

    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
	       "dielectric 1.0",
	       "pair_modify tail yes mix arithmetic",
	       "special_bonds lj/coul 0.0 0.0 1.0",
	       "kspace_style ewald 1.0e-5",
	       "pair_coeff 1 1 0.025000 2.183593",
	       "pair_coeff 2 2 0.238000 3.028165",
	       "pair_coeff 3 3 0.309000 3.768502"]
    
    atom_types = {'Li': 1, 'Ca': 2, 'As': 3}

    atom_type_masses = {'Li': 6.941, 'Ca': 40.078, 'As': 74.9216}

    log_file = "m73766_auto.log"

    lammps_header = ["units real",
		     "atom_style full",
		     "boundary p p p"]

    amendments = ["set type 1 charge 0.84",
		  "set type 2 charge 1.32",
		  "set type 3 charge -2.16"]

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
######   Automatically generated lammps method m74395_auto   #####
@subdir_calc
def m74395_auto(atoms: ase.Atoms):

    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
	       "special_bonds lj/coul 0.0 0.0 1.0",
	       "dielectric 1.0",
	       "pair_modify tail yes mix arithmetic",
	       "kspace_style ewald 1.0e-5",
	       "pair_coeff 1 1 0.025000 2.183593",
	       "pair_coeff 2 2 0.379000 3.813047",
	       "pair_coeff 3 3 0.060000 3.118146"]
    
    atom_types = {'Li': 1, 'Ge': 2, 'O': 3}

    atom_type_masses = {'Li': 6.941, 'Ge': 72.64, 'O': 15.9994}

    log_file = "m74395_auto.log"

    lammps_header = ["units real",
		     "atom_style full",
		     "boundary p p p"]

    amendments = ["set type 1 charge 0.85",
		  "set type 2 charge 1.27",
		  "set type 3 charge -1.17"]

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
######   Automatically generated lammps method m74865_auto   #####
@subdir_calc
def m74865_auto(atoms: ase.Atoms):

    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
	       "special_bonds lj/coul 0.0 0.0 1.0",
	       "dielectric 1.0",
	       "pair_modify tail yes mix arithmetic",
	       "kspace_style ewald 1.0e-5",
	       "pair_coeff 1 1 0.025000 2.183593",
	       "pair_coeff 2 2 0.505000 4.008153",
	       "pair_coeff 3 3 0.402000 3.826410",
	       "pair_coeff 4 4 0.060000 3.118146"]
    
    atom_types = {'Li': 1, 'Al': 2, 'Si': 3, 'O': 4}

    atom_type_masses = {'Li': 6.941, 'Al': 26.9815386, 'Si': 28.0855, 'O': 15.9994}

    log_file = "m74865_auto.log"

    lammps_header = ["units real",
		     "atom_style full",
		     "boundary p p p"]

    amendments = ["set type 1 charge 0.91",
		  "set type 2 charge 1.67",
		  "set type 3 charge 1.8",
		  "set type 4 charge -1.09"]

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
######   Automatically generated lammps method m75073_auto   #####
@subdir_calc
def m75073_auto(atoms: ase.Atoms):

    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
	       "pair_modify tail yes mix arithmetic",
	       "special_bonds lj/coul 0.0 0.0 1.0",
	       "dielectric 1.0",
	       "kspace_style ewald 1.0e-5",
	       "pair_coeff 1 1 0.025000 2.183593",
	       "pair_coeff 2 2 0.228000 2.537280",
	       "pair_coeff 3 3 0.402000 3.826410",
	       "pair_coeff 4 4 0.060000 3.118146"]
    
    atom_types = {'Li': 1, 'Cd': 2, 'Si': 3, 'O': 4}

    atom_type_masses = {'Li': 6.941, 'Cd': 112.411, 'Si': 28.0855, 'O': 15.9994}

    log_file = "m75073_auto.log"

    lammps_header = ["units real",
		     "atom_style full",
		     "boundary p p p"]

    amendments = ["set type 1 charge 0.86",
		  "set type 2 charge 1.34",
		  "set type 3 charge 1.31",
		  "set type 4 charge -1.1"]

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
######   Automatically generated lammps method m75138_auto   #####
@subdir_calc
def m75138_auto(atoms: ase.Atoms):

    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
	       "dielectric 1.0",
	       "special_bonds lj/coul 0.0 0.0 1.0",
	       "pair_modify tail yes mix arithmetic",
	       "kspace_style ewald 1.0e-5",
	       "pair_coeff 1 1 0.025000 2.183593",
	       "pair_coeff 2 2 0.111000 2.691405",
	       "pair_coeff 3 3 0.402000 3.826410"]
    
    atom_types = {'Li': 1, 'Mg': 2, 'Si': 3}

    atom_type_masses = {'Li': 6.941, 'Mg': 24.305, 'Si': 28.0855}

    log_file = "m75138_auto.log"

    lammps_header = ["units real",
		     "atom_style full",
		     "boundary p p p"]

    amendments = ["set type 1 charge 0.62",
		  "set type 2 charge 0.96",
		  "set type 3 charge -2.2"]

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
######   Automatically generated lammps method m75630_auto   #####
@subdir_calc
def m75630_auto(atoms: ase.Atoms):

    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
	       "dielectric 1.0",
	       "pair_modify tail yes mix arithmetic",
	       "special_bonds lj/coul 0.0 0.0 1.0",
	       "kspace_style ewald 1.0e-5",
	       "pair_coeff 1 1 0.235000 3.243762",
	       "pair_coeff 2 2 0.449000 3.937772",
	       "pair_coeff 3 3 0.025000 2.183593"]
    
    atom_types = {'Sr': 1, 'Sb': 2, 'Li': 3}

    atom_type_masses = {'Sr': 87.62, 'Sb': 121.76, 'Li': 6.941}

    log_file = "m75630_auto.log"

    lammps_header = ["units real",
		     "atom_style full",
		     "boundary p p p"]

    amendments = ["set type 1 charge 1.36",
		  "set type 2 charge -2.22",
		  "set type 3 charge 0.86"]

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
######   Automatically generated lammps method m75961_auto   #####
@subdir_calc
def m75961_auto(atoms: ase.Atoms):

    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
	       "pair_modify tail yes mix arithmetic",
	       "dielectric 1.0",
	       "special_bonds lj/coul 0.0 0.0 1.0",
	       "kspace_style ewald 1.0e-5",
	       "pair_coeff 1 1 0.025000 2.183593",
	       "pair_coeff 2 2 0.398000 3.982317",
	       "pair_coeff 3 3 0.069000 3.260689"]
    
    atom_types = {'Li': 1, 'Te': 2, 'N': 3}

    atom_type_masses = {'Li': 6.941, 'Te': 127.6, 'N': 14.0067}

    log_file = "m75961_auto.log"

    lammps_header = ["units real",
		     "atom_style full",
		     "boundary p p p"]

    amendments = ["set type 1 charge 0.64",
		  "set type 1 charge 0.64",
		  "set type 2 charge -1.4",
		  "set type 2 charge -1.4",
		  "set type 3 charge -1.86",
		  "set type 3 charge -1.86"]

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
######   Automatically generated lammps method m76276_auto   #####
@subdir_calc
def m76276_auto(atoms: ase.Atoms):

    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
	       "pair_modify tail yes mix arithmetic",
	       "dielectric 1.0",
	       "special_bonds lj/coul 0.0 0.0 1.0",
	       "kspace_style ewald 1.0e-5",
	       "pair_coeff 1 1 0.025000 2.183593",
	       "pair_coeff 2 2 0.238000 3.028165",
	       "pair_coeff 3 3 0.013000 2.594297",
	       "pair_coeff 4 4 0.050000 2.996983"]
    
    atom_types = {'Li': 1, 'Ca': 2, 'Fe': 3, 'F': 4}

    atom_type_masses = {'Li': 6.941, 'Ca': 40.078, 'Fe': 55.845, 'F': 18.9984032}

    log_file = "m76276_auto.log"

    lammps_header = ["units real",
		     "atom_style full",
		     "boundary p p p"]

    amendments = ["set type 1 charge 0.82",
		  "set type 2 charge 1.62",
		  "set type 3 charge -0.59",
		  "set type 3 charge 1.07",
		  "set type 4 charge -0.59"]

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
######   Automatically generated lammps method m76595_auto   #####
@subdir_calc
def m76595_auto(atoms: ase.Atoms):

    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
	       "special_bonds lj/coul 0.0 0.0 1.0",
	       "dielectric 1.0",
	       "pair_modify tail yes mix arithmetic",
	       "kspace_style ewald 1.0e-5",
	       "pair_coeff 1 1 0.025000 2.183593",
	       "pair_coeff 2 2 0.069000 2.783168",
	       "pair_coeff 3 3 0.050000 2.996983"]
    
    atom_types = {'Li': 1, 'Zr': 2, 'F': 3}

    atom_type_masses = {'Li': 6.941, 'Zr': 91.224, 'F': 18.9984032}

    log_file = "m76595_auto.log"

    lammps_header = ["units real",
		     "atom_style full",
		     "boundary p p p"]

    amendments = ["set type 1 charge 0.92",
		  "set type 2 charge 2.31",
		  "set type 3 charge -0.69"]

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
######   Automatically generated lammps method m77893_auto   #####
@subdir_calc
def m77893_auto(atoms: ase.Atoms):

    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
	       "dielectric 1.0",
	       "pair_modify tail yes mix arithmetic",
	       "special_bonds lj/coul 0.0 0.0 1.0",
	       "kspace_style ewald 1.0e-5",
	       "pair_coeff 1 1 0.025000 2.183593",
	       "pair_coeff 2 2 0.019000 2.935511",
	       "pair_coeff 3 3 0.305000 3.694557",
	       "pair_coeff 4 4 0.060000 3.118146"]
    
    atom_types = {'Li': 1, 'Sc': 2, 'P': 3, 'O': 4}

    atom_type_masses = {'Li': 6.941, 'Sc': 44.955912, 'P': 30.973762, 'O': 15.9994}

    log_file = "m77893_auto.log"

    lammps_header = ["units real",
		     "atom_style full",
		     "boundary p p p"]

    amendments = ["set type 1 charge 0.91",
		  "set type 2 charge 1.94",
		  "set type 3 charge 1.64",
		  "set type 4 charge -0.88"]

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
######   Automatically generated lammps method m77899_auto   #####
@subdir_calc
def m77899_auto(atoms: ase.Atoms):

    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
	       "special_bonds lj/coul 0.0 0.0 1.0",
	       "pair_modify tail yes mix arithmetic",
	       "dielectric 1.0",
	       "kspace_style ewald 1.0e-5",
	       "pair_coeff 1 1 0.025000 2.183593",
	       "pair_coeff 2 2 0.072000 2.980056",
	       "pair_coeff 3 3 0.056000 2.719023",
	       "pair_coeff 4 4 0.060000 3.118146"]
    
    atom_types = {'Li': 1, 'Y': 2, 'Mo': 3, 'O': 4}

    atom_type_masses = {'Li': 6.941, 'Y': 88.90585, 'Mo': 95.96, 'O': 15.9994}

    log_file = "m77899_auto.log"

    lammps_header = ["units real",
		     "atom_style full",
		     "boundary p p p"]

    amendments = ["set type 1 charge 0.78",
		  "set type 2 charge 2.14",
		  "set type 3 charge 1.43",
		  "set type 4 charge -0.9"]

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
######   Automatically generated lammps method m78086_auto   #####
@subdir_calc
def m78086_auto(atoms: ase.Atoms):

    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
	       "pair_modify tail yes mix arithmetic",
	       "dielectric 1.0",
	       "special_bonds lj/coul 0.0 0.0 1.0",
	       "kspace_style ewald 1.0e-5",
	       "pair_coeff 1 1 0.025000 2.183593",
	       "pair_coeff 2 2 0.059000 2.819694",
	       "pair_coeff 3 3 0.291000 3.746229"]
    
    atom_types = {'Li': 1, 'Nb': 2, 'Se': 3}

    atom_type_masses = {'Li': 6.941, 'Nb': 92.90638, 'Se': 78.96}

    log_file = "m78086_auto.log"

    lammps_header = ["units real",
		     "atom_style full",
		     "boundary p p p"]

    amendments = ["set type 1 charge 0.53",
		  "set type 2 charge 0.57",
		  "set type 3 charge -0.55"]

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
######   Automatically generated lammps method m78099_auto   #####
@subdir_calc
def m78099_auto(atoms: ase.Atoms):

    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
	       "dielectric 1.0",
	       "pair_modify tail yes mix arithmetic",
	       "special_bonds lj/coul 0.0 0.0 1.0",
	       "kspace_style ewald 1.0e-5",
	       "pair_coeff 1 1 0.040000 3.665157",
	       "pair_coeff 2 2 0.025000 2.183593",
	       "pair_coeff 3 3 0.016000 2.800986",
	       "pair_coeff 4 4 0.274000 3.594776"]
    
    atom_types = {'Rb': 1, 'Li': 2, 'V': 3, 'S': 4}

    atom_type_masses = {'Rb': 85.4678, 'Li': 6.941, 'V': 50.9415, 'S': 32.065}

    log_file = "m78099_auto.log"

    lammps_header = ["units real",
		     "atom_style full",
		     "boundary p p p"]

    amendments = ["set type 1 charge 0.58",
		  "set type 2 charge 0.29",
		  "set type 3 charge 0.16",
		  "set type 4 charge -0.8"]

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
######   Automatically generated lammps method m78355_auto   #####
@subdir_calc
def m78355_auto(atoms: ase.Atoms):

    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
	       "dielectric 1.0",
	       "pair_modify tail yes mix arithmetic",
	       "special_bonds lj/coul 0.0 0.0 1.0",
	       "kspace_style ewald 1.0e-5",
	       "pair_coeff 1 1 0.025000 2.183593",
	       "pair_coeff 2 2 0.238000 3.028165",
	       "pair_coeff 3 3 0.015000 2.693187",
	       "pair_coeff 4 4 0.050000 2.996983"]
    
    atom_types = {'Li': 1, 'Ca': 2, 'Cr': 3, 'F': 4}

    atom_type_masses = {'Li': 6.941, 'Ca': 40.078, 'Cr': 51.9961, 'F': 18.9984032}

    log_file = "m78355_auto.log"

    lammps_header = ["units real",
		     "atom_style full",
		     "boundary p p p"]

    amendments = ["set type 1 charge 0.8",
		  "set type 2 charge 1.62",
		  "set type 3 charge 1.2",
		  "set type 4 charge -0.6"]

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
