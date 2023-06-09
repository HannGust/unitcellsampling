"""Program that reads a certain lammpsinterface output and can create a dictionary of atom types: atomic charges. 
Works for a certain specially formatted version (the format the lammps-parsers in ucs can read). Also, currently assumes only
one charge per atom type."""


import os, sys, re, ase, ase.io, argparse, numpy as np

from unitcellsampling.lammps_calc_from_inp import parser_lammps_mel_inp,\
                                 mel_lmp_read_charges_file
                                
from unitcellsampling.read_atom_types import read_lammpsdata_atom_info 
from pathlib import Path


def label_charge_dict_from_parsed_charge_input(charges_input):
    """Extracts the numeric charges from the output of mel_lmp_read_charges_file
    and returns a label:charge dictionary, that can be used as a mapping from
    the numneric atom labels to charges."""
    
    label_charge_dict = {}
    for i, line in enumerate(charges_input):
        label_charge_match = re.findall("[-0-9.]*[0-9]+", line)
        print(label_charge_match)   
        if (not label_charge_match) or len(label_charge_match) != 2:
            raise Exception("ERROR: Processed a line and found unexpected number of matches (!= 2)!")

        label_charge_dict[int(label_charge_match[0])] = float(label_charge_match[1])
    # Debug:
    print(label_charge_dict)
    #
    return label_charge_dict
    

# Define the argument parser
parser = argparse.ArgumentParser()

parser.add_argument("txtfile", type=str, action="store",
    help="Textfile to write charge dictionary to.")
parser.add_argument("datafile", type=str, action="store",
    help="Lammps data file to parse for atomic labels.")
parser.add_argument("chargefile", type=str, action="store",
    help="Lammps set_charges file to read charges from.")


args = parser.parse_args()

#
txtfile = args.txtfile
datafile = args.datafile
chargefile = args.chargefile
#

# Create label-to-charge mapping from charge file:
labels_charges_dict = label_charge_dict_from_parsed_charge_input(mel_lmp_read_charges_file(chargefile))

# Create a label-to-atom symbol mapping from datafile:
atom_labels, atom_masses, chem_symbs = read_lammpsdata_atom_info(datafile) 
atom_label_to_symbol_dict = dict(zip(atom_labels, chem_symbs))
#atom_symbol_to_label_dict = dict(zip(chem_symbs, atom_labels))

# LETS SAY WE HAVE MULTIPLE LABELS PER ATOM. THEN HOW DO WE DO IT?
# WHAT HAPPENS RATHER? WELL IF WE MAKE A LABEL TO CHARGE MAPPING 
# AND A LABEL TO SYMBOL MAPPINGS, we get the full duplicate thing.
# This is good because we have te full correspondence. When we would construct a 
# new dictionary mapping symbols to charges, we get say "Na" : +1
# FROM ALL LABELS CORREPSONDING TO THE ATOM Na AND TS CHARGE,
# so multiplicity is taken care of by re-overwriting the same 
# dictionary entry, but since the pair are consistent they do always
# define he same mapping, i.e. all pairs dictate that Na has charge +1.
# unless we overwrite it into a dictionary, that works. Nut we have to 
# make a dictionary first because otherwise we get multiple of the 
# same pair if we loop over labels.

# Sanity check: Check that labels are the same:
print(labels_charges_dict.keys())
print(atom_label_to_symbol_dict.keys())
for key in labels_charges_dict.keys():
    assert key in atom_label_to_symbol_dict.keys(), "ERROR: Dictionaries with charges and atom symbols contain different labels!" 
#assert labels_charges_dict.keys()) == atom_label_to_symbol_dict.keys(), "ERROR: Dictionaries with charges and atom symbols contain different number of labels!"
assert len(set(labels_charges_dict.values())) == len(set(atom_label_to_symbol_dict.values())), "ERROR: Different number of unique atomic symbols as unique charges!"
#

# Now covert into symbol-charge dict
charge_dict = {}
for label,symbol in atom_label_to_symbol_dict.items():
    charge_dict[symbol] = labels_charges_dict[label]

##

# Then write to txt-file:



if Path(txtfile).exists() and Path(txtfile).is_file():
    with open(txtfile, 'a') as file:
        for symbol,charge in charge_dict.items():
            file.write(" ".join([symbol, str(charge), "\n"]))
else:
    with open(txtfile, 'w') as file:
        for symbol,charge in charge_dict.items():
            file.write(" ".join([symbol, str(charge), "\n"]))


