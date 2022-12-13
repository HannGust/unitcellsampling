
import sys
import numpy as np
from pathlib import Path
import re

import ase # only needed for sanity check

#############################################################################
# Function and script to read atomic symbols and masses from a lammps       #
# data file, in order to pass this information to ase and ucs. Turns out    #
# that ase read_lammps_input has problems identifying atomic species,       #
# hence the need for this function. It is tailored to input in the          #
# 'Melania_data' folder, and might not work as well for other input         #
# formats.                                                                  #
#############################################################################


def read_lammpsdata_atom_info(data_file, atomic_mass_diff_threshold=0.05):
    with open(data_file) as fp:
        data_content = fp.readlines()
    
    print("Atomic mass threshold, for issuing a warning, is "+str(atomic_mass_diff_threshold))

    for i in range(len(data_content)):
        data_content[i] = data_content[i].strip('\n ')
    
    indx1 = data_content.index("Masses") 
    indx2 = data_content.index("Pair Coeffs")
    
    atoms_masses_content = data_content[indx1:indx2]
    
    atom_type_lines = []
    for line in atoms_masses_content:
        if len(line) > 0:
            if (line[0].isdigit()):
                atom_type_lines.append(line)
    
    atom_labels = []
    atom_masses = []
    chem_symbs = []
    
    for line in atom_type_lines:
        atom_label_match = re.search('^[0-9]+ ', line)
        atom_mass_match = re.search(' [0-9]+.[0-9]+ ', line)
        chem_symb_match = re.search('# ([A-Z][a-z]*)', line)
        
        atom_labels.append(atom_label_match.group().strip())
        atom_masses.append(atom_mass_match.group().strip())
        chem_symbs.append(chem_symb_match.group(1).strip())
   
    for a,m in zip(chem_symbs,atom_masses):
        mass_diff = np.abs(ase.Atom(a).mass - float(m))
        if mass_diff >= atomic_mass_diff_threshold:
            raise Warning("ASE mass and read mass differs with "+str(mass_diff)+", i.e. more than "+str(atomic_mass_diff_threshold)+" units.")

    return (atom_labels, atom_masses, chem_symbs)
    


if __name__ == "__main__":
    
        
    data_file = sys.argv[1]
    
    with open(data_file) as fp:
        data_content = fp.readlines()
    
    
    for i in range(len(data_content)):
        data_content[i] = data_content[i].strip('\n ')
    
    indx1 = data_content.index("Masses")
    
    indx2 = data_content.index("Pair Coeffs")
    
    atoms_masses_content = data_content[indx1:indx2]
    
    atom_type_lines = []
    for line in atoms_masses_content:
        if len(line) > 0:
            if (line[0].isdigit()):
                atom_type_lines.append(line)
    
    atom_labels = []
    atom_masses = []
    chem_symbs = []
    
    for line in atom_type_lines:
        # Line structure, atom label no, atom mass, # chemical_symbol
        # (+some extra characters)
        atom_label_match = re.search('^[0-9]+ ', line)
        atom_mass_match = re.search(' [0-9]+.[0-9]+ ', line)
        chem_symb_match = re.search('# ([A-Z][a-z]*)', line)
        
        #print(atom_label_match, atom_mass_match, chem_symb_match, sep="\n")
        #print("atom_label: ", atom_label_match.group().strip())
        #print("atom_mass: ", atom_mass_match.group().strip())
        #print("chem_symb: ", chem_symb_match.group(1).strip())
    
        atom_labels.append(atom_label_match.group().strip())
        atom_masses.append(atom_mass_match.group().strip())
        chem_symbs.append(chem_symb_match.group(1).strip())
    
    
    
    # Sanity check like so. Check that the ase mass of the atom is close
    # to mass read from file.
    atomic_mass_diff_threshold = 0.05
    print("Atomic mass threshold, for issuing a warning, is "+str(atomic_mass_diff_threshold))

    for i,a in enumerate(chem_symbs):
        tmp_mass = ase.Atom(a).mass
        mass = float(atom_masses[i])
        print(mass, "  ", tmp_mass)
        print("Diff: ", abs(mass-tmp_mass))
        print("Less than ", atomic_mass_diff_threshold, ": ", np.abs(mass-tmp_mass) < atomic_mass_diff_threshold)
    
    print("\n"*3)
    print("atom_labels: ", atom_labels)
    print("atom_masses: ", atom_masses)
    print("chem_symbs: ", chem_symbs)

    # Test 
    (a_l_2, a_m_2, c_s_2) = read_lammpsdata_atom_info(data_file)
    print("atom_labels: ", a_l_2)
    print("atom_masses: ", a_m_2)
    print("chem_symbs: ", c_s_2)

