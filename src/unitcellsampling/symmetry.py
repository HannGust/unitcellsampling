# Contains functions relating to spacegroup symmetry and unitcell type

import gemmi
import numpy as np
import ase
import pymatgen as pmg
import argparse
import itertools as it

import ase.io
import ase.spacegroup
import pymatgen.symmetry as pmgsym

#from ase.spacegroup import get_spacegroup
from pymatgen.core import Structure, Lattice
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.analysis.structure_matcher import StructureMatcher


import unitcellsampling.preparatory_fcns
from unitcellsampling.preparatory_fcns import remove_nonframework_cations_fancy


# Here I define the functions that can be used to find the grid shape 
# compatible with the spacegroup 


def is_shape_and_spacegroup_compatible(shape, spacegroup):
    """ Tests whether the grid shape and spacegroup are compatible according
        to gemmi's grid.set_size() function.
    """
    test_grid = gemmi.FloatGrid()
    test_grid.spacegroup = spacegroup

    try:
       test_grid.set_size(*shape)
    except:
       return False
    else:
       return True 
        

def find_spacegroup_compatible_gridshape(grid_shape, spacegroup,
                                    search_denser=True):
    """ Starting from a given grid shape, searches for the closest grid
        shape compatible with the given spacegroup. 
    """
    # Make sure spacegroup is of the right type
    assert isinstance(spacegroup, (gemmi.SpaceGroup, 
                                   ase.spacegroup.Spacegroup, int, str))

    if isinstance(spacegroup, (int, ase.spacegroup.Spacegroup)):
        sg = gemmi.find_spacegroup_by_number(spacegroup)
    elif isinstance(spacegroup, str):
        sg = gemmi.find_spacegroup_by_name(spacegroup)
    elif isinstance(spacegroup, gemmi.SpaceGroup):
        sg = spacegroup
    else:
        raise Error("Something went wrong: spacegroup argument is"
                    +" not correct type.")


    if search_denser:
        shape_variations = tuple(it.product((0,1), repeat=3))[0:-1]
        increment_factor = 1
    else:
        shape_variations = tuple(it.product((0,-1), repeat=3))[0:-1]
        increment_factor = -1

    sg_compatible = False
    current_grid_shape = grid_shape

    while not sg_compatible:

        for shape in shape_variations:
            compatible_grid_shape = tuple(np.array(current_grid_shape)
                                          + np.array(shape))
            
            sg_compatible = is_shape_and_spacegroup_compatible(
                                                       compatible_grid_shape,
                                                       sg)

            if sg_compatible:
                break
        
        current_grid_shape = tuple(np.array(current_grid_shape) 
                                   + np.array([1,1,1]) * increment_factor)
        
    return compatible_grid_shape


def is_conventional_cell(unitcell, ltol=1.0E-05, stol=1.0E-05, angle_tol=0.1):
    """Uses pymatgens structure matcher to check whether the given 
       structure (unitcell) is equivalent (up to some tolerance) to the
       conventional standard cell that pymatgen can construct from it.
    """

    assert isinstance(unitcell, (ase.Atoms, pmg.core.Structure)), \
                      "Error: unitcell must be ase.Atoms object or \
                      pymatgen.core.Structure object"

    if isinstance(unitcell, ase.Atoms):
        f_atomic_coords = unitcell.get_scaled_positions()
        atomic_symb = unitcell.get_chemical_symbols()
        cell_params = unitcell.cell.cellpar()

        lattice = Lattice.from_parameters(*cell_params)
        struct = Structure(lattice, atomic_symb, f_atomic_coords)
    else:
        struct = unitcell
    
    # Initialize the spacegroup analyzer to find conventional
    # cell of structure
    spgrp_analyzer = SpacegroupAnalyzer(struct)
    conv_cell = spgrp_analyzer.get_conventional_standard_structure()
    
    # Define the structure matcher
    structure_matcher = StructureMatcher(ltol=ltol,
                                         stol=stol,
                                         angle_tol=angle_tol,
                                         primitive_cell=False) 

    return structure_matcher.fit(struct, conv_cell, symmetric=True)


def is_primitive_cell(unitcell, ltol=1.0E-05, stol=1.0E-05, angle_tol=0.1):
    """Uses pymatgens structure matcher to check whether the given structure
       (unitcell) is equivalent (up to some tolerance) to the primitive
       cell that pymatgen can construct from it. """

    assert isinstance(unitcell, (ase.Atoms, pmg.core.Structure)), \
                      "Error: unitcell must be ase.Atoms object or \
                      pymatgen.core.Structure object" 

    if isinstance(unitcell, ase.Atoms):
        f_atomic_coords = unitcell.get_scaled_positions()
        atomic_symb = unitcell.get_chemical_symbols()
        cell_params = unitcell.cell.cellpar()

        lattice = Lattice.from_parameters(*cell_params)
        struct = Structure(lattice, atomic_symb, f_atomic_coords)
    else:
        struct = unitcell
    
    # Initialize the spacegroup analyzer to find primitive 
    # cell of structure
    spgrp_analyzer = SpacegroupAnalyzer(struct, symprec=0.01, angle_tolerance=0.1)
    #prim_cell = spgrp_analyzer.find_primitive()
    prim_cell = spgrp_analyzer.get_primitive_standard_structure()
   
    # Define the structure matcher
    structure_matcher = StructureMatcher(ltol=ltol,
                                         stol=stol,
                                         angle_tol=angle_tol,
                                         primitive_cell=False) 

    return structure_matcher.fit(struct, prim_cell, symmetric=True)


