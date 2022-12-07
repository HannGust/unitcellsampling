# Short python script taking a cif file, reading it,
# calculating the spacegroup and updates that information
# finally printing an updated cif file containing the 
# symmetry information.

import sys
import os
from pathlib import Path
import ase
from ase.io import read, write
from ase.spacegroup import get_spacegroup


if (len(sys.argv) > 1):
    file = sys.argv[1]
else:
    print("Need at least one arg (file).")
    sys.exit()

if (len(sys.argv) > 2):
    tol = float(sys.argv[2])
else:
    tol = False

atoms = read(file)

if tol:
    sg = get_spacegroup(atoms,tol)
else:
    sg = get_spacegroup(atoms)

print("Found spacegroup: ",sg,sg.symbol)
atoms.info["spacegroup"] = sg


new_file = str(Path(file).stem)+"_symmetrized"+str(Path(file).suffix)
print("New filename: ", new_file)


print("New spacegroup: ",atoms.info["spacegroup"].symbol)
write(filename=new_file, images=atoms, format="cif")
