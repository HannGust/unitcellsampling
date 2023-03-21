""" Test of new midvox option in unitcellsampling.sample.UnitCelSampler.generate_grid_vectors """

import sample
import ase
import numpy as np

atoms = ase.Atoms(symbols=["H","H","O"], positions=[[0,0,0], [1,2,3], [3,1,2]], cell=[[6,0,0],[0,8,0],[0,0,10]], pbc=True)

sampler = sample.UnitCellSampler(atoms)

grid_points_midvox, included_midvox= sampler.generate_grid_vectors((10,10,10), midvox=True)
print(grid_points_midvox)

sampler = sample.UnitCellSampler(atoms)

grid_points, included= sampler.generate_grid_vectors((10,10,10), midvox=False)

print(grid_points_midvox)
