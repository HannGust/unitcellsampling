"""Test keep_alive option for LAMMPSlib methods"""

import numpy as np
import ase, ase.io, ase.build
from ase.io.cube import write_cube
from unitcellsampling import sample
import os
from pathlib import Path
from preparatory_fcns import remove_nonframework_cations_fancy,\
                             compute_req_supercells
from decorators import subdir_calc
from ase.calculators.lammpsrun import LAMMPS
from ase.calculators.lammpslib import LAMMPSlib

# Idea of structure: 
# 1. Define two methods for a small structure that are identical except for
#    the keep_alive parameter.
# 2. Sequentially use both methods to sample the structure, in designated 
#    directory.
# 3. Record the time taken, and compare the resulting grids.

# Define the methods

@subdir_calc
def struct_55319_ff(atoms: ase.Atoms):
    """Universal forcefield method for structure with keep_alive unset (default setting, should be False)."""
    
    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
               "pair_modify     tail yes mix arithmetic",
               "special_bonds   lj/coul 0.0 0.0 1.0",
               "dielectric      1.0",
               "kspace_style ewald 1.0e-5",
               "pair_coeff 1 1 0.025000 2.183593",
               "pair_coeff 2 2 0.015000 2.524807",
               "pair_coeff 3 3 0.069000 3.260689"] 
    
    atom_types = {"Li":1, "Ni":2, "N":3}
    atom_type_masses = {"Li":6.9410, "Ni":58.693400000, "N":14.006700000}
    log_file = "55319_job.log"
    
    lammps_header = ["units real",
                     "atom_style full",
                     "boundary p p p",
                     "box tilt large"]

    amendments = ["set type 1 charge 0.59",
                  "set type 2 charge 0.03",
                  "set type 3 charge -0.63"]

    #post_changebox_cmds = []

    calc = LAMMPSlib(lmpcmds=lmpcmds,
                     atom_types=atom_types,
                     atom_type_masses=atom_type_masses,
                     log_file=log_file,
                     lammps_header=lammps_header,
                     amendments=amendments)
    
    atoms.set_calculator(calc)
    return atoms.get_potential_energy()



@subdir_calc
def struct_55319_ff_keep_alive(atoms: ase.Atoms):
    """Universal forcefield method for structure with keep_alive set to True."""
    
    lmpcmds = ["pair_style lj/cut/coul/long 12.500",
               "pair_modify     tail yes mix arithmetic",
               "special_bonds   lj/coul 0.0 0.0 1.0",
               "dielectric      1.0",
               "kspace_style ewald 1.0e-5",
               "pair_coeff 1 1 0.025000 2.183593",
               "pair_coeff 2 2 0.015000 2.524807",
               "pair_coeff 3 3 0.069000 3.260689"] 
    
    atom_types = {"Li":1, "Ni":2, "N":3}
    atom_type_masses = {"Li":6.9410, "Ni":58.693400000, "N":14.006700000}
    log_file = "55319_job.log"
    
    lammps_header = ["units real",
                     "atom_style full",
                     "boundary p p p",
                     "box tilt large"]

    amendments = ["set type 1 charge 0.59",
                  "set type 2 charge 0.03",
                  "set type 3 charge -0.63"]

    #post_changebox_cmds = []

    calc = LAMMPSlib(lmpcmds=lmpcmds,
                     atom_types=atom_types,
                     atom_type_masses=atom_type_masses,
                     log_file=log_file,
                     lammps_header=lammps_header,
                     amendments=amendments,
                     keep_alive=True)
    
    atoms.set_calculator(calc)
    return atoms.get_potential_energy()

#


# Set the indata
structure_no = "55319"
#atoms_file = "../Melania_data/init_cif/" + structure_no + ".cif"
atoms_file = "/home/hannes/Desktop/Python_Coding/BenjaminBolbrinker_Code/grid_sampling_rundir/Melania_data/init_cif/" + structure_no + ".cif"
#/home/hannes/Desktop/Python_Coding/BenjaminBolbrinker_Code/grid_sampling_rundir/Melania_data/init_cif/


# Set the outdata names
cube_out_1 = "_".join((structure_no, "default", ".cube"))
cube_out_2 = "_".join((structure_no, "keep_alive_true", ".cube"))
#

# Read indata
atoms_init = ase.io.read(atoms_file)
#

# Prepare atoms object
atoms = remove_nonframework_cations_fancy(atoms=atoms_init, cation=ase.Atom('Li'))

cutoff = 12.5 # Force cutoff in Å, used to determine supercell
num_cells = compute_req_supercells(atoms, cutoff)
print("num cells: ", num_cells)
supercell = ase.build.make_supercell(atoms, np.diag((num_cells, num_cells, num_cells)))
#


# Define samplers for supercells
ucs_sampler = sample.UnitCellSampler(supercell)
#ucs_2 = sample.UnitCellSampler(supercell)

n_frac = (40,40,40) # Corresponds to an approximate spacing of X.Y Å CHOOSE THIS
vdw_scale = 0.0

grid_coords, included_points = ucs_sampler.generate_grid_vectors(n_frac, vdw_scale=vdw_scale)

# Despite vdv = 0.0, it still exlcudes points...
print(len(grid_coords.shape))
#print(np.product(included_points.shape), np.count_nonzero(included_points))
#assert included_points.all(), "Not all points inlcuded, something went wrong!"



## Preform calculations, method 1, with keep_alive left as default
os.environ["UCS_CALCULATION_DIR"] = "ucs_calcs/" + structure_no + "_default"
#

energies_1 = ucs_sampler.calculate_energies(grid_coords.reshape(-1,3), method=struct_55319_ff, atom='Li', exploit_symmetry=False)

##

## Preform calculations, method 2, with keep_alive=True
# Set envronment variables (working dir etc.)
os.environ["UCS_CALCULATION_DIR"] = "ucs_calcs/" + structure_no + "_keep_alive_true"
#

energies_2 = ucs_sampler.calculate_energies(grid_points=grid_coords.reshape(-1,3), method=struct_55319_ff_keep_alive, atom='Li', exploit_symmetry=False)

##


# Write resulting grids to cubefiles
write_cube(cube_out_1, supercell, data=energies_1.reshape(n_frac))

write_cube(cube_out_2, supercell, data=energies_2.reshape(n_frac))

diff_cube_name = "_".join(("diff", structure_no)) + ".cube"
write_cube(diff_cube_name, supercell, data=np.abs(energies_1 - energies_2).reshape(n_frac))


max_diff = np.max(np.abs(energies_1 - energies_2))
rmsd = np.sqrt(np.mean(np.square(np.subtract(energies_1, energies_2))))

with open("grid_comparison.txt",'w') as f:
    f.write("Maximum absolute difference: "+str(max_diff))
    f.write("RMSD: "+str(rmsd))

#
