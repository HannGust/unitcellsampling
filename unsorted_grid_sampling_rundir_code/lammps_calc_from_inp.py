import ase
from ase.calculators.lammpslib import LAMMPSlib
from decorators import subdir_calc
from read_atom_types import read_lammpsdata_atom_info

# unsure if needed
import os
import sys
import numpy as np
import math as ma

##############################################################################
# Creates a lammps energy calculator from given lammps input to use with UCS #
# i.e. a function with signature ase.Atoms -> float64.                       #
#                                                                            #
# Should also make it produce a text version of the method and print it.     #
##############################################################################


##############################################################################
# The structure of the methods:     #
"""@subdir_calc
def struct_55319_ff(atoms: ase.Atoms):
    
    lmpcmds = ["pair_style lj/cut/coul/long 12.500",      #
               "pair_modify     tail yes mix arithmetic", #
               "special_bonds   lj/coul 0.0 0.0 1.0",     #   THESE FIVE ARE FROM lammps "in."-file
               "dielectric      1.0",                     #
               "kspace_style ewald 1.0e-5",               #
               "pair_coeff 1 1 0.025000 2.183593",        ##
               "pair_coeff 2 2 0.015000 2.524807",        ##  PAIR POTENTIALS FROM lammps "data."-file
               "pair_coeff 3 3 0.069000 3.260689"]        ##
    
    atom_types = {"Li":1, "Ni":2, "N":3}                                  # FROM DATA FILE - Have read_atom_types fcn for these 2
    atom_type_masses = {"Li":6.9410, "Ni":58.693400000, "N":14.006700000} # FROM DATA FILE
    log_file = "55319_job.log" # -> should be constructred from input!
    
    lammps_header = ["units real",      #
                     "atom_style full", #   THESE ARE ALL FROM lammps "in."-file
                     "boundary p p p",  #
                     "box tilt large"]  # 

    amendments = ["set type 1 charge 0.59",   #
                  "set type 2 charge 0.03",   # THESE ARE FROM lammps "set_charges"-file
                  "set type 3 charge -0.63"]  #

    #post_changebox_cmds = []

    calc = LAMMPSlib(lmpcmds=lmpcmds,
                     atom_types=atom_types,
                     atom_type_masses=atom_type_masses,
                     log_file=log_file,
                     lammps_header=lammps_header,
                     amendments=amendments)
    
    atoms.set_calculator(calc)
    return atoms.get_potential_energy()"""
##############################################################################


def reduce_whitespace(string: str) -> str:
    """Removes superfluous (>1) whitespaces from string."""
    reduced_str = " ".join(string.split())
    return reduced_str

# Has not been fully tested
def lammps_method_from_data(lmpcmds=None,
                            atom_types=None,
                            atom_type_masses=None,
                            log_file=None,
                            lammps_header=None,
                            amendments=None,
                            post_changebox_cmds=None):
    """
    Constructs a LAMMPSlib-based energy calculator method to be used with UCS
    """
    # Section to define/read the lammps input
    # lmpcmds =
    # atom_types =
    # atom_type_masses =
    # log_file =
    # lammps_header =
    # amendments =
    # post_changebox_cmds =

    init_input_list = [atom_types,
                       atom_type_masses,
                       log_file,
                       lammps_header,
                       amendments,
                       post_changebox_cmds]

    input_names = ["atom_types",
                   "atom_type_masses",
                   "log_file",
                   "lammps_header",
                   "amendments",
                   "post_changebox_cmds"]

    indices = [i for i in range(len(init_input_list))
               if init_input_list[i] is not None]

    input_list = [i for i in init_input_list if (i is not None)]

    assert len(indices) == len(input_list), "Mismatch in length. "
    assert [init_input_list.copy()[i]
            for i in indices] == input_list, "List mismatch."

    input_dict = dict()
    for i in indices:
        input_dict[input_names[i]] = init_input_list[i]

    # Debug:
    print(input_dict)

    @subdir_calc
    def lammps_method(atoms: ase.Atoms):
        calc = LAMMPSlib(lmpcmds, **input_dict)

        atoms.set_calculator(calc)
        return atoms.get_potential_energy()

    return lammps_method


def lammps_method_from_data_to_text(method_name,
                                    lmpcmds=None,
                                    atom_types=None,
                                    atom_type_masses=None,
                                    log_file=None,
                                    lammps_header=None,
                                    amendments=None,
                                    post_changebox_cmds=None):
    """
    Constructs a LAMMPSlib-based energy calculator method in plain text.
    """

    lammps_method_text = """
@subdir_calc
def {method_name}(atoms: ase.Atoms):

    lmpcmds = {lmpcmds}
    
    atom_types = {atom_types}

    atom_type_masses = {atom_type_masses}

    log_file = {log_file}

    lammps_header = {lammps_header}

    amendments = {amendments}

    post_changebox_cmds = {post_changebox_cmds}

    calc = LAMMPSlib(lmpcmds=lmpcmds,
                     atom_types=atom_types,
                     atom_type_masses=atom_type_masses,
                     log_file=log_file,
                     lammps_header=lammps_header,
                     amendments=amendments,
                     post_changebox_cmds=post_changebox_cmds)
    
    atoms.set_calculator(calc)

    return atoms.get_potential_energy()
""".format(method_name=method_name,
               lmpcmds="[\""+"\",\n\t       \"".join(lmpcmds)+"\"]",
               atom_types=atom_types,
               atom_type_masses=atom_type_masses,
               log_file=log_file,
               lammps_header="[\""+"\",\n\t\t     \"".join(lammps_header)+"\"]",
               amendments="[\""+"\",\n\t\t  \"".join(amendments)+"\"]",
               post_changebox_cmds=post_changebox_cmds)

    return lammps_method_text


# A function to parse the data specifically for data from the dataset from Melania
# Lets try it as this, simple, first

def mel_read_pp_from_lmp_data_file(datafile):
    """Reads pair potentials from a lammps data file."""
    with open(datafile) as f_lmp:
        lmp_content = f_lmp.readlines()

    for i, line in enumerate(lmp_content):
        lmp_content[i] = line.strip("\n")

    line_no_1 = lmp_content.index("Pair Coeffs")
    line_no_2 = lmp_content.index("Atoms")

    pp_content = lmp_content[line_no_1+1:line_no_2]

    for i, line in enumerate(pp_content):
        if line.strip() == "":
            del pp_content[i]

    print(pp_content)

    pp_list = []
    for line in pp_content:
        pp_list.append(line.split()[0:3])

    lmpcmds_pp = []
    for i, line in enumerate(pp_list):
        lmpcmds_pp.append("pair_coeff {0} {0} {1} {2}".format(*line))

    return lmpcmds_pp


def mel_lmp_read_in_file(infile):
    """Reads relevant lines from an lammps in file."""
    with open(infile) as f_lmp:
        lmp_content = f_lmp.readlines()

    for i, line in enumerate(lmp_content):
        lmp_content[i] = line.strip("\n")

    lammps_header = lmp_content[1:4]

    lmpscmds_in = lmp_content[4:11]

    # delete empty lines and strip superlfuous whitespaces
    for i, line in enumerate(lammps_header):
        if line.strip() == "":
            del lammps_header[i]

    for i, line in enumerate(lammps_header):
        lammps_header[i] = reduce_whitespace(lammps_header[i])

    for i, line in enumerate(lmpscmds_in):
        if line.strip() == "":
            del lmpscmds_in[i]

    for i, line in enumerate(lmpscmds_in):
        lmpscmds_in[i] = reduce_whitespace(lmpscmds_in[i])

    return lammps_header, lmpscmds_in


def mel_lmp_read_charges_file(charge_file):
    """Read an parse the charge file."""
    with open(charge_file) as cf:
        charges = cf.readlines()
        charges.sort()

    for i, line in enumerate(charges):
        charges[i] = line.strip("\n")

    for i, line in enumerate(charges):
        if line.strip() == "":
            del charges[i]

    for i, line in enumerate(charges):
        charges[i] = reduce_whitespace(line)

    return charges


def parser_lammps_mel_inp(datafile, infile, charge_file):
    """Read lammps input from data, input and charges file."""
    lammps_header, lmpcmds_in = mel_lmp_read_in_file(infile)
    lmpcmds_pp = mel_read_pp_from_lmp_data_file(datafile)
    charges = mel_lmp_read_charges_file(charge_file)
    atom_labels, atom_masses, chem_symbs = read_lammpsdata_atom_info(datafile)

    lmpcmds = lmpcmds_in
    lmpcmds.extend(lmpcmds_pp)
    # Debug
    print("lmpcmds")
    print(lmpcmds)
    atom_types = dict(zip(chem_symbs, atom_labels))
    atom_type_masses = dict(zip(chem_symbs, atom_masses))
    amendments = charges

    # Remove log_file and post_changebox_cmds from this dict - redundant
    lmps_dict = {"lmpcmds": lmpcmds,
                 "atom_types": atom_types,
                 "atom_type_masses": atom_type_masses,
                 "lammps_header": lammps_header,
                 "amendments": amendments,
                 }

    return lmps_dict
