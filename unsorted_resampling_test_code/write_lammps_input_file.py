import ase
import ase.io
from ase.io.cube import read_cube
from ase.io import read
import ase.calculators
from ase.calculators.lammpsrun import LAMMPS

import sys
import os
from pathlib import Path


os.environ['ASE_LAMMPSRUN_COMMAND'] = "/home/hannes/lammps-23Jun2022/build_benbolenv/lmp"

sample_atom = ase.Atom("Na", [0.0,0.0,0.0])

if len(sys.argv) == 1:
    raise Exception("Need atom structure input.")

inp = sys.argv[1]

with open(inp, 'r') as fp:
    atoms = read(inp)
    assert isinstance(atoms,ase.Atoms), "atoms is not ase.Atoms object!!!"

atoms = atoms + sample_atom
out = Path(inp).stem + "_with_Na.lmpdata"
#else:
#    out = Path(inp).stem + ".lmpdata"

atoms.calc = LAMMPS(keep_tmp_files=True)

energy = atoms.get_potential_energy()
print("Enrgy of system with LAMMPS() genereic calculator: ", energy)

#with open(out, 'w') as fp:
#    write_lammps_data(out, atoms)
