
import os, sys, ase, ase.io, argparse, numpy as np
from unitcellsampling.lammps_calc_from_inp import parser_lammps_mel_inp,\
                                 mel_read_pp_from_lmp_data_file,\
                                 mel_lmp_read_in_file,\
                                 mel_lmp_read_charges_file,\
                                 lammps_method_from_data_to_text
                                

from pathlib import Path

parser = argparse.ArgumentParser()

parser.add_argument("txtfile", type=str, action="store",
    help="Textfile to write (append) method to.")
parser.add_argument("method_name", type=str, action="store",
    help="Naming of the method, i.e. of the python function.") # Also affects logfile if logfile is none. NOPE, not currently! Changed this behaviour, because we don't want a logfile now.")
parser.add_argument("infile", type=str, action="store",
    help="Lammps input file to parse for info to construct method.")
parser.add_argument("datafile", type=str, action="store",
    help="Lammps data file to parse for info to construct method.")
parser.add_argument("chargefile", type=str, action="store",
    help="Lammps set_charges file to parse for info to construct method.")
parser.add_argument("--logfile", type=str, action="store", default=None,
    help="Explicit naming of the logfile. Optional. Default is None, in which case logilfe is set to None, so that no logfile is produced.")


args = parser.parse_args()

#
method_name = args.method_name
txtfile = args.txtfile
infile = args.infile
datafile = args.datafile
chargefile = args.chargefile

if args.logfile:
    if args.logfile[-4:] == ".log":
        logfile = args.logfile
    else:
        logfile = args.logfile + ".log"
else:
    #logfile=".".join((method_name,"log"))
    logfile = None
#

lmps_dict = parser_lammps_mel_inp(datafile, infile, chargefile, atomic_mass_diff_threshold=0.1)
lmps_dict["log_file"] = logfile


method_txt = lammps_method_from_data_to_text(method_name, **lmps_dict)

# Method comment:
top_comment = """
####################################################################################
#                                                                                  #
# This file contains automatically generated LAMMPSlib-energy calculators for UCS. #
#                                                                                  #
####################################################################################
"""
import_statements = """
import ase
from ase.calculators.lammpsrun import LAMMPS
from ase.calculators.lammpslib import LAMMPSlib
from decorators import subdir_calc

"""
header_comment = "#####   Automatically generated lammps method {}   #####".format(method_name)

if Path(txtfile).exists() and Path(txtfile).is_file():
    with open(txtfile, 'a') as file:
        file.write("#")
        file.write("\n")
        file.write(header_comment)
        file.write(method_txt)
        file.write("\n")
else:
    with open(txtfile, 'w') as file:
        file.write(top_comment)
        file.write("\n")
        file.write(header_comment)
        file.write(method_txt)
        file.write("\n")
