import numpy as np


def create_supercell_grid(unit_cell: np.array, grid: np.array, dim=(1, 1, 1)):
    supercell_grid = np.empty([g*d for (g, d) in zip(grid.shape, dim)])

    return supercell_grid
