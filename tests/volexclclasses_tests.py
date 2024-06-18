#!/usr/bin/env python

"""Test of volume exlcuder classes.
Originally taken and adapted from the volume_exclusion development repo.
Can be generalised and formatted better, but tests the basics of the
functions and methods."""

import ase
import ase.visualize
import numpy as np
from unitcellsampling.volume_exclusion import RadialExcluder, ScaledVdWExcluder

def define_test_atoms_obj():
    cell = np.eye(3)*3.0 # 3.0 Å cell in ech direction, orthogonal
    atom_pos = [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]
    symbols= "OH"
    atoms = ase.Atoms(symbols=symbols,scaled_positions=atom_pos, cell=cell)
    return atoms

def generate_test_grid_coord(grid_coord_base=None):
    if not grid_coord_base:
        grid_coord_base = np.array([0.25,0.75]) 

    X,Y,Z = np.meshgrid(grid_coord_base,grid_coord_base,grid_coord_base, indexing="ij")
    grid_coord = np.stack((X,Y,Z), axis=-1).reshape(-1,3)

    return grid_coord
    

def create_dummy_visualization_atom_obj(atoms, grid_coord):
    dummy_atom = ["X"]*grid_coord.shape[0]
    vis_atoms_object = atoms.copy() +  ase.Atoms(symbols=dummy_atom, scaled_positions=grid_coord, cell=atoms.get_cell())
    ase.visualize.view(vis_atoms_object, data=grid_coord)



def test_radial_excluder(visualize=False):
    test_successful = True

    atoms = define_test_atoms_obj()
    grid_coord = generate_test_grid_coord()


    print(grid_coord @ atoms.get_cell())
    
    if visualize:
        create_dummy_visualization_atom_obj(atoms, grid_coord)

    radii={"O":1.00, "default":1.38}
    vol_excl = RadialExcluder(radii)

    filters = vol_excl.construct_cutoff_filters(atoms, grid_coord, frac_input=True, periodic=True)
    vol_excl.print_settings()
    print(vol_excl.get_full_radii_mapping(atoms))
    print("filters:", filters)
    print("Filters are complete: ", np.logical_or(filters[0], filters[1]).all())
    print("Filters are disjiont: ", not np.logical_and(filters[0], filters[1]).any())

    if not np.logical_or(filters[0], filters[1]).all():
        test_successful = False
    
    if np.logical_and(filters[0], filters[1]).any():
        test_successful = False

    return test_successful

def test_vdw_excluder(visualize=False):
    test_successful = True

    atoms = define_test_atoms_obj()
    grid_coord = generate_test_grid_coord()

    vdw_scale_mapping={"O":0.7, "H":1.3}
    vdw_excl = ScaledVdWExcluder(vdw_scaling_map=vdw_scale_mapping)
    try:
        vdw_excl.print_settings()
    except Exception as exc:
        print("An exception occurred 1:", exc)
        test_successful = False
    
    excl_obj = vdw_excl.generate_excluder(atoms)
    try:
        vdw_excl.print_settings()
    except:
        print("An exception occurred 2:", exc)
        test_successful = False
    
    vdw_filters = vdw_excl.construct_cutoff_filters(atoms, grid_coord, frac_input=True)
    print("filters:", vdw_filters)

    if visualize:
        create_dummy_visualization_atom_obj(atoms, grid_coord)
    
    return test_successful

# Change the main function to what you want to test
def main():
    print("Running excluder tests")
    print("======================")
    print("1. Radial excluder test:")
    radial_excluder_test_status = test_radial_excluder()
    print("Passed: ", radial_excluder_test_status)
    print()
    print("2. van der Waals excluder test:")
    vdw_excluder_test_status = test_vdw_excluder()
    print("Passed: ", vdw_excluder_test_status)
    print()
    print("======================")
    print("Tests complete.")


if __name__ == "__main__":
    main()
