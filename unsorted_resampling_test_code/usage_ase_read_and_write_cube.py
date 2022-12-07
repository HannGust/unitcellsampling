import ase
import ase.io
from ase.io import read
from ase.io.cube import write_cube
from ase.io.cube import read_cube

import numpy as np
import sys
from pathlib import Path

#
use_read_cube_init = False
use_read_cube_final= True
#

inputfile = sys.argv[1]
# First read a cube or xsf file
if use_read_cube_init:
    with open(inputfile) as fp:
       cube_content = read_cube(fp)
       atoms = cube_content['atoms']
       data = cube_content['data']
else:
    #with open(inputfile) as fp:
    atoms = read(inputfile)
    data = read(inputfile, read_data=True)

print("data type and shape: ", type(data), data.shape)

# Then write a cube file
outputfile = Path(inputfile).stem + "_NEW.cube"

with open(outputfile, 'w') as fp:
    write_cube(fp, atoms, data=data)

# Then try to read the cubefile again with ase
if use_read_cube_final:
    with open(outputfile, 'r') as fp:
        new_content = read_cube(fp)
        new_atoms = new_content['atoms']
        new_data = new_content['data']
else:
    new_atoms = read(outputfile)
    new_data = read(outputfile, read_data=True)

print("data type and shape: ", type(new_data), new_data.shape)

print("Are atoms objects equal: ", (atoms == new_atoms))
print("Are grid objects equal: ", np.isclose(data, new_data))


