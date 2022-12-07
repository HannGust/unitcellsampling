
# Check if gridpoints are inside cell

import numpy as np
import ase, ase.io
import argparse
from pathlib import Path

# Have not been tested:
def coords_within_cell_boundaries(coords, atoms, cart_coord=False):
    """Checks whether the coordinates are inside a unitcell."""

    cell = atoms.get_cell()

    if cart_coord:
        frac_coords = np.linalg.solve(cell.T, coords.T).T
    else:
        frac_coords = coords
    
    is_in_cell =  ((0.0 <= frac_coords[:,0]).all() and (frac_coords[:,0] <= 1.0).all() and\
                   (0.0 <= frac_coords[:,1]).all() and (frac_coords[:,1] <= 1.0).all() and\
                   (0.0 <= frac_coords[:,2]).all() and (frac_coords[:,2] <= 1.0).all())

    return is_in_cell


def find_min_shift(shift_incr, coords, atoms, cart_coord=False):
    """Finds the minimum shift that makes sure 
    coordinates are within the unitcell."""
    
    a, b, c = atoms.get_cell()[0,:], atoms.get_cell()[1,:], atoms.get_cell()[2,:]

    shift = 0.0
    tmp_coords = coords
    while not coords_within_cell_boundaries(tmp_coords, atoms, cart_coord=cart_coord):
        shift += shift_incr
        tmp_coords = coords + shift * (a + b + c)
        # Debug
        print("inside func find_min_shift: current shift = ", shift)

    return shift

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("coords", action="store", type=str, help="Text file containing coordinates of the grid.")
    parser.add_argument("atoms", action="store", type=str, help="File containing the unitcell.")
    parser.add_argument("-s", "--shift", action="store", type=float, default=0.0,
                        help="Specifies that a certain shift should applied\
                             to the coordinates before the test.\
                            \"shift\" here means a displacement of the form\
                            s * (a + b + c), where a, b, c are the cell\
                            vectors, and the shift value is the value of\
                            the scaling factor s.")
    parser.add_argument("--ms", "--minshift", action="store", type=float,
                        default=0.0,
                        help="Indicates that it should calculate the minmum\
                        shift making sure the gridpoints are within the unit\
                        cell. The shift increment should be given.")

    args = parser.parse_args()

    coords_file = args.coords
    atoms_file = args.atoms
    
    # Read coordinates and atoms object
    coords = np.loadtxt(coords_file)
    atoms = ase.io.read(atoms_file)

    if args.ms != 0.0:
        min_shift = find_min_shift(args.ms, coords, atoms, cart_coord=True)
        print("Minimum shift for ", atoms_file, " and ", coords_file, ": ", min_shift) 

    elif args.shift != 0.0:
        a, b, c = atoms.get_cell()[0,:], atoms.get_cell()[1,:], atoms.get_cell()[2,:]
        coords = coords + args.shift * (a + b + c)
        print(coords_within_cell_boundaries(coords, atoms, True))

    else:
        print(coords_within_cell_boundaries(coords, atoms, True))
        


