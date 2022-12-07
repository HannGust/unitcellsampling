### TEST OF ATOM REMOVING FUNCTION 
from preparatory_fcns import remove_nonframework_caions_fancy, remove_nonframework_cations
from ase import Atom
from ase.io import read, write
from ase.io.cube import write_cube
from pathlib import Path

parent_path = Path("Melania_data")
data_path = Path(parent_path, "data_lmp")
cif_path =  Path(parent_path, 'init_cif', '54189.cif')

#with open(cif_path) as fp:
#test_atoms = read(fp)
test_atoms = read(cif_path)


cation = Atom('Li')

catremoved_atoms_1 = remove_nonframework_cations(test_atoms, cation)
catremoved_atoms_2 = remove_nonframework_caions_fancy(test_atoms, cation)

print(test_atoms)
print(catremoved_atoms_1)
print(catremoved_atoms_2)
