import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import ase
import ase.io
from ase.io import read
import glob
import sys
import os

from pathlib import Path
from plotting import plot_grid

assert len(sys.argv) > 2, "Need at least 2 inputs."

atoms = read(sys.argv[1])
file = sys.argv[2]



with open(file, 'r') as f:
    coords = f.readlines()


for i,line in enumerate(coords):
    new_line = line.split('_')
    new_line = list(map(float, new_line))
    coords[i] = new_line

coords = np.array(coords)
#print(coords)

data = np.ones(coords.shape[0])


plot_grid(atoms, data, coords)