import ase
import numpy as np
import gemmi
from ase.spacegroup import get_spacegroup
import itertools as it

def find_grid_assymmetric_domain(atoms=None, grid=None, n_frac=None, spacegroup=None):
    
    assert atoms or spacegroup

    if atoms and not spacegroup:
        spacegroup = get_spacegroup(atoms)
    
    assert grid or n_frac

    if grid:
        n_frac = grid.shape
        assert len(n_frac) == 3
        
    bool_grid_counted = gemmi.Int8Grid(*n_frac)
    bool_grid_counted.spacegroup = gemmi.find_spacegroup_by_number(spacegroup)

    frac_grid_coords = np.array(list(it.product(np.linspace(0, 1, n_frac[0], endpoint=False),
                                                np.linspace(0, 1, n_frac[1], endpoint=False), 
                                                np.linspace(0, 1, n_frac[2], endpoint=False)))).reshape(*n_frac, 3)
    
    frac_grid_coords_list = frac_grid_coords.reshape(-1,3)
    
    cart_grid_coords_list = frac_grid_coords @ atoms.get_cell()
    
    # fast writing pseudocode here
    
    indx_list = []
    f_coord_list = []
    f_coord_list_2 = []
    num_points = 0
    print("Before loop")
    for i in range(n_frac[0]):
        print("i=",i, "out of",n_frac[0])
        for j in range(n_frac[1]):
            for k in range(n_frac[2]):
                if not bool_grid_counted.get_value(i,j,k):
                    bool_grid_counted.set_value(i,j,k, 1)
                    bool_grid_counted.symmetrize_max()
                    indx_list.append((i,j,k))
                    f_coord_list.append(frac_grid_coords[i,j,k,:])
                    f_coord_list_2.append(frac_grid_coords_list[i*n_frac[1]*n_frac[2] + j*n_frac[2] + k])
                    num_points += 1
                else:
                    continue

    c_coord_list = np.array(f_coord_list) @ np.array(atoms.get_cell())

    print("Num points:" , num_points)
    assert (np.array(f_coord_list) == np.array(f_coord_list_2)).all(), "f_coord_lists are not the same!!!"
    return f_coord_list, indx_list, c_coord_list


        
    




