#!/usr/bin/env python

"""Tests of spacegroup matching and grid shape compatibility functions.
OBS: NOT DONE.
TODO: Implement tests.
"""
import gemmi
import numpy as np
import unitcellsampling.symmetry
from unitcellsampling.symmetry import is_shape_compatible, find_spacegroup_compatible_gridshape, find_grid_shape_closest_to_spacing

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


def main():
    test_find_spacegroup_compatible_gridshape()

if __name__ == "__main__":
    main()     
