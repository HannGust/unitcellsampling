import numpy as np, matplotlib, matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
fig = plt.figure()
ax = plt.axes(projection="3d")
energies = np.loadtxt("ase_lammps_calc_54189_wierd_coord.txt")
coords = np.loadtxt("54189_wierd_coords.txt")
threshold = 1
energy_mask = energies < threshold
energies_plot = energies[energy_mask]
coords_plot = coords[energy_mask]
ax.scatter3D(coords_plot[:,1], coords_plot[:,2], energies_plot[:], c=energies_plot[:], cmap="brg")
plt.show()

