#!/usr/bin/env python

"""Test of grid point coordinate (vector) generation
method, generate_grid_vectors, in the UnitCellSampler class.
Made after update of the radial/van der Waals cutoff routines."""

### Example, calculation of distances and scaling factor
### for vdw radii for a 3,3,3 grid in te test structure generated 
### with the function below assuming a vdw radius of O of 1.52 Å.
#
#In [1]: import numpy as np
#
#In [2]: rO = np.array([3.0, 2.4, 5.0])
#
#In [3]: p1_incl = np.array([4.0,4.0,4.0])
#
#In [4]: p2_excl = np.array([2.0,2.0,2.0])
#
#In [5]: np.linalg.norm(np.subtract(rO,p1_incl))
#Out[5]: 2.1354156504062622
#
#In [6]: np.linalg.norm(np.subtract(rO,p2_excl))
#Out[6]: 3.1874754901018454
#
#In [7]: min_k = np.linalg.norm(np.subtract(rO,p1_incl))/1.52
#
#In [8]: max_k = np.linalg.norm(np.subtract(rO,p2_excl))/1.52
#
#In [9]: print(min_k)
#1.404878717372541
#In [10]: print(max_k)
#2.097023348751214
#
#In [12]: np.mean([min_k, max_k])
#Out[12]: 1.7509510330618776
#
######

import ase
import numpy as np
from mendeleev import element
import unitcellsampling.sample
from unitcellsampling.sample import UnitCellSampler

def generate_test_atoms():
    """Generates a test atoms object for the grid vector generation.
    This is a single O atom at fractional coordinates [1/2, 2/5, 5/6] = [0.5, 0.4, 5/6]
    within an orthogonal cell with dimensions 6.0 Å x 6.0 Å x 6.0 Å. That is, the 
    cartesian coordinates of the oxygen atom is [3.0 Å, 2.4 Å, 5.0 Å]."""
    symbols = "O"
    cell = np.eye(3) * 6.0
    frac_coord = [[1.0/2.0, 2.0/5.0, 5.0/6.0]]
    cart_coord = [[3.0, 2.4, 5.0]] # in Å
    atoms = ase.Atoms(symbols=symbols, scaled_positions=frac_coord, cell=cell, pbc=True)

    # make sanity checks tthat this is correct:
    assert np.allclose(np.array(cart_coord), atoms.get_positions())
    assert [8] == atoms.get_atomic_numbers()
    
    return atoms



def test_1_generate_grid_vectors():
    atoms = generate_test_atoms()
    
    # Assert that we assume the correct radius of 1.52 Å of Oxygen here
    assert (element(atoms[0].symbol).vdw_radius)/100.0 == 1.52
    
    sampler = UnitCellSampler(atoms=atoms)
    
    n_frac = (3,3,3)
    cutoff_radii = np.sqrt(3) * 2.0 / 2.0
    vdw_scale = 1.75 # See calculation above
    grid_vectors, included = sampler.generate_grid_vectors(n_frac=n_frac, cutoff_radii=cutoff_radii, vdw_scale=vdw_scale)
    
    
    #############################################
    ####        Define correct values        ####
    #############################################
    
    # Explicit definition of the 3x3x3 grid fractional coordinates
    correct_scaled_coords = np.array([[0.0, 0.0, 0.0],
                               [0.0, 0.0, 1.0/3.0],
                               [0.0, 0.0, 2.0/3.0],
                               [0.0, 1.0/3.0, 0.0],
                               [0.0, 1.0/3.0, 1.0/3.0],
                               [0.0, 1.0/3.0, 2.0/3.0],
                               [0.0, 2.0/3.0, 0.0],
                               [0.0, 2.0/3.0, 1.0/3.0],
                               [0.0, 2.0/3.0, 2.0/3.0],
                               [1.0/3.0, 0.0, 0.0],
                               [1.0/3.0, 0.0, 1.0/3.0],
                               [1.0/3.0, 0.0, 2.0/3.0],
                               [1.0/3.0, 1.0/3.0, 0.0],
                               [1.0/3.0, 1.0/3.0, 1.0/3.0],
                               [1.0/3.0, 1.0/3.0, 2.0/3.0],
                               [1.0/3.0, 2.0/3.0, 0.0],
                               [1.0/3.0, 2.0/3.0, 1.0/3.0],
                               [1.0/3.0, 2.0/3.0, 2.0/3.0],
                               [2.0/3.0, 0.0, 0.0],
                               [2.0/3.0, 0.0, 1.0/3.0],
                               [2.0/3.0, 0.0, 2.0/3.0],
                               [2.0/3.0, 1.0/3.0, 0.0],
                               [2.0/3.0, 1.0/3.0, 1.0/3.0],
                               [2.0/3.0, 1.0/3.0, 2.0/3.0],
                               [2.0/3.0, 2.0/3.0, 0.0],
                               [2.0/3.0, 2.0/3.0, 1.0/3.0],
                               [2.0/3.0, 2.0/3.0, 2.0/3.0]]
                               )
    
    # Explicit definition of the 3x3x3 grid cartesian coordinates (in Å)
    correct_cartesian_coords_1 = np.array([[0.0, 0.0, 0.0],
                                           [0.0, 0.0, 2.0],
                                           [0.0, 0.0, 4.0],
                                           [0.0, 2.0, 0.0],
                                           [0.0, 2.0, 2.0],
                                           [0.0, 2.0, 4.0],
                                           [0.0, 4.0, 0.0],
                                           [0.0, 4.0, 2.0],
                                           [0.0, 4.0, 4.0],
                                           [2.0, 0.0, 0.0],
                                           [2.0, 0.0, 2.0],
                                           [2.0, 0.0, 4.0],
                                           [2.0, 2.0, 0.0],
                                           [2.0, 2.0, 2.0],
                                           [2.0, 2.0, 4.0],
                                           [2.0, 4.0, 0.0],
                                           [2.0, 4.0, 2.0],
                                           [2.0, 4.0, 4.0],
                                           [4.0, 0.0, 0.0],
                                           [4.0, 0.0, 2.0],
                                           [4.0, 0.0, 4.0],
                                           [4.0, 2.0, 0.0],
                                           [4.0, 2.0, 2.0],
                                           [4.0, 2.0, 4.0],
                                           [4.0, 4.0, 0.0],
                                           [4.0, 4.0, 2.0],
                                           [4.0, 4.0, 4.0]]
                                           )
    
    # infer correct cartesian coordinates from correct scaled
    correct_cartesian_coords_2 = correct_scaled_coords @ np.array(atoms.cell)
    
    assert np.all(correct_cartesian_coords_1 == correct_cartesian_coords_2), "Correct Cartesian coords not equal."
    
    correct_vdw_included = np.full(27, True, dtype=bool)
    # Loop over indices corresponding to the points within the cutoff:
    # Here we have four more points, since vdw cutoff is larger
    # It corresponds to the closest encapslating cube formed by the grid points
    for idx in [12,14,21,23,15,17,24,26]:
        correct_vdw_included[idx] = False
    
    correct_cutoff_included = np.full(27, True, dtype=bool)
    # Loop over indices corresponding to the points within the cutoff:
    # This corresponds to the closest square of grid points
    for idx in [12,14,21,23]:
        correct_cutoff_included[idx] = False
    
    
    correct_included = np.logical_and(correct_cutoff_included, correct_vdw_included)
    
    correct_cartesian_coords_1 = correct_cartesian_coords_1.reshape(n_frac + (3,))
    correct_included = correct_included.reshape(n_frac)
    correct_cutoff_included = correct_cutoff_included.reshape(n_frac)
    correct_vdw_included = correct_vdw_included.reshape(n_frac)

    
    correct_values = {"grid_vectors":correct_cartesian_coords_1,
                      "included": correct_included,
                      "cutoff_included":correct_cutoff_included,
                      "vdw_included":correct_vdw_included
                     }
    
    #############################################
    ####          End of definition          ####
    #############################################
    
    ##### Now do the comparison to test: ####
    
    assert np.all(correct_values["grid_vectors"] == grid_vectors)
    assert np.all(correct_values["included"] == included)
    assert np.all(correct_values["cutoff_included"] == sampler.cutoff_included)
    assert np.all(correct_values["vdw_included"] == sampler.vdw_included)
    
    # Here we test the printing of the settings, and access of the
    # excluders that are set in generate_grid_vectors method.
    print("Radial excluder settings:")
    sampler.radial_cutoff_excluder.print_settings()

    print()
    print("vdW excluder settings:")
    sampler.scaled_vdw_excluder.print_settings()

    return True



def main():
    if test_1_generate_grid_vectors():
        print("Successful test!")


if __name__ == "__main__":
    main()