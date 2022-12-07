# This is just a script to cut out one eigth of the unit cell of the zeolite LTA4A 
# to see whether DFT calculations works on this.

import numpy as np
import ase
from ase.io import read, write
import sys




file = sys.argv[1]

atoms = read(file)

cell = atoms.get_cell()
one_eight_cell = cell * 0.5

a, b, c, alfa, beta, gamma = atoms.get_cell_lengths_and_angles()

print(a, b, c, alfa, beta, gamma) 
print("original cell: ", cell)
print("reduced cell: ", one_eight_cell)

atom_pos = atoms.get_positions()

reduced_atoms_obj = atoms.copy()
indx_atoms_to_remove = []
for i,coords in enumerate(atom_pos):
    assert (coords == atom_pos[i]).all(), "coords not equal to atom_pos[i], what!??"
    print(atom_pos[i])
    print(coords)
    x,y,z = coords
    if not ((x < a/2) and (y < b/2) and (z < c/2)):
        indx_atoms_to_remove.append(i)
        print("Atom ",i," should be removed.")

indx_atoms_to_remove.reverse()
print(indx_atoms_to_remove)
for j in indx_atoms_to_remove:
    reduced_atoms_obj.pop(j)


reduced_atoms_obj.set_cell(one_eight_cell)


write("LTA4A_reduced.cif",reduced_atoms_obj, "cif")
