#!/usr/bin/env python

"""Tests of spacegroup matching and grid shape compatibility functions.
OBS: NOT DONE.
TODO: Implement tests.
"""
import ase.lattice
import ase.spacegroup
import ase.spectrum
import ase.structure
import gemmi
import ase
import numpy as np
import unitcellsampling.symmetry
from unitcellsampling.symmetry import is_shape_compatible, find_spacegroup_compatible_gridshape, find_grid_shape_closest_to_spacing, prepare_matching_structure_and_spacegroup


# Construct the 530 first entries of the gemmi spacegroup table 
gemmi_spgrp_table = list(gemmi.spacegroup_table_itb())[:530]

def test_find_spacegroup_compatible_gridshape():
    grid_shape = (10, 10, 10)
    a = 5
    b = 5
    c = 5
    tmp_grid = gemmi.FloatGrid()
    
    for i, spgrp in enumerate(gemmi_spgrp_table):
        #try:
        tmp_shape = find_spacegroup_compatible_gridshape(grid_shape, spgrp, a, b, c, search_extent=3)
        
        #    if "Unexpected RuntimeError occured" in str(exc):
        #        print("RuntimeError handled: Need another grid shape for spgrp entry {} : {}".format(i, spgrp) )
                #print(exc.__context__)
                #continue
        #    else:
        #        raise RuntimeError from exc
        assert isinstance(tmp_shape, list)
        assert len(tmp_shape) == 3
        assert all([isinstance(x, int) and x >= 10 for x in tmp_shape])
        print("spgrp no: {}  -   tmp_shape: {}".format(i, tmp_shape))
        tmp_grid.spacegroup = spgrp
        tmp_grid.set_size(*tmp_shape)



def test_prepare_matching_structure_and_spacegroup():
    
    # spacegroup 190 - hexagonal system
    cellpar = [5, 5, 7.5, 90, 90, 120]
    nacl_190_conv = ase.spacegroup.crystal("NaCl", basis=[[0.0,0.0,0.0], [0.5, 0.0, 0.0]], cellpar=cellpar, ab_normal=(0.3, 0.75, 0.9), spacegroup=190, pbc=True)
    nacl_190_prim = ase.spacegroup.crystal("NaCl", basis=[[0.0,0.0,0.0], [0.5, 0.0, 0.0]], cellpar=cellpar, ab_normal=(0.3, 0.75, 0.9), spacegroup=190, pbc=True, primitive_cell=True)

    # spacegroup 74 - Orthorhombic system
    cellpar = [4, 6, 9, 90, 90, 90]
    nacl_74_conv = ase.spacegroup.crystal("NaCl", basis=[[0.0,0.0,0.0], [0.5, 0.0, 0.0]], cellpar=cellpar, ab_normal=(0.3, 0.75, 0.9), spacegroup=74, pbc=True)
    nacl_74_prim = ase.spacegroup.crystal("NaCl", basis=[[0.0,0.0,0.0], [0.5, 0.0, 0.0]], cellpar=cellpar, ab_normal=(0.3, 0.75, 0.9), spacegroup=74, pbc=True, primitive_cell=True)

    
    #atoms = ase.Atoms()
    for atoms in [nacl_190_conv, nacl_190_prim, nacl_74_conv, nacl_74_prim]:
        matching_atoms, matching_spacegroup, basis_change_tuple = prepare_matching_structure_and_spacegroup(atoms, spglib_standardize="on", change_basis="on")
        print(matching_atoms, matching_spacegroup, basis_change_tuple, "\n")


def test_prepare_matching_structure_and_spacegroup_on_NaCl():
    nacl_atoms = ase.spacegroup.crystal("NaCl", basis=[[0.0, 0.0, 0.0],[0.5, 0.5, 0.5]], cellpar=[5.64], spacegroup=225, primitive_cell=True)
    matching_atoms, matching_spacegroup, basis_change_tuple = prepare_matching_structure_and_spacegroup(nacl_atoms, spglib_standardize="allow", change_basis="allow")
    print(matching_atoms, matching_spacegroup, basis_change_tuple, "\n")


def main():
    test_find_spacegroup_compatible_gridshape()
    test_prepare_matching_structure_and_spacegroup()
    test_prepare_matching_structure_and_spacegroup_on_NaCl()

if __name__ == "__main__":
    main()     
