# sampling a set of coordinates with unitcellsampler

import numpy as np
import ase
import ase.io
import ase.build
from unitcellsampling import sample
from pathlib import Path
from special_methods import structure_54189_ff_manually_coded
import os
from decorators import subdir_calc
from preparatory_fcns import remove_nonframework_cations_fancy


os.environ["UCS_CALCULATION_DIR"] = "./ucs_dir"

atom_file ="../Melania_data/init_cif/54189.cif"
coord_file ="54189_wierd_coords.txt"

atoms_init = ase.io.read(atom_file)
supcell_size = 8
P = supcell_size * np.diag(np.ones(3))

atoms = remove_nonframework_cations_fancy(ase.build.make_supercell(atoms_init, P), ase.Atom("Li"))

coords = np.loadtxt(coord_file)


sampler = sample.UnitCellSampler(atoms)
#sampler.grid_vectors = coords
grid = sampler.calculate_energies(coords, atom='Li', method=structure_54189_ff_manually_coded, exploit_symmetry=False)

np.savetxt("ase_lammps_calc_54189_wierd_coord.txt", grid)

