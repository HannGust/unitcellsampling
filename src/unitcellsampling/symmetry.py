# Contains functions relating to spacegroup symmetry and unitcell type

"""Contains functions and routines relating to symmetry,
both for structure manipulation and standardization, and
for preparation for application to grids.
All custom symmetry functions are defined here."""

import spglib
import gemmi
import numpy as np
import ase
import pymatgen as pmg
import itertools as it
import copy

import ase.io
import ase.spacegroup

from pymatgen.core import Structure, Lattice
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.io.ase import AseAtomsAdaptor


# Here I define the functions that can be used to find the grid shape 
# compatible with the spacegroup 

# TODO Use the new grid shape compatibility functions to write a new one that returns the closest grid shape??

def is_shape_compatible(shape, gemmigrid):
    """Tests whether the grid shape and gemmi grid (with associated spacegroup)
    are compatible according to gemmi's grid.set_size() function.
    """
    try:
       gemmigrid.set_size(*shape)
    except RuntimeError as rt_exc:
       if "Grid not compatible with the space group" in str(rt_exc):
           return False
       else:
           print("Unexpected RuntimeError occurred:")
           raise RuntimeError("Unexpected RuntimeError occured.") from rt_exc
    else:
       return True


# NOTE: Has been updated. Now takes spacegroup an makes dummy grid.
# Also changed "sizes" in function and variable names to "shapes".
def search_compatible_grid_shapes(gemmi_spgrp, start_shape, search_extent=10, mode="up"):
    """Searches for compatible grid sizes in a range."""
    assert isinstance(gemmi_spgrp, gemmi.SpaceGroup), "Error: gemmi_spgrp must be a gemmi.SpaceGroup!"
    assert isinstance(start_shape, (tuple, list)) and len(start_shape) == 3, "Error: Start shape must be tuple/list of lenght 3!"
    assert all([isinstance(i, int) for i in start_shape]), "Error: Start shape must be tuple/list of integers!"

    gen_of_shapes = it.product(range(start_shape[0], start_shape[0] + search_extent),
                               range(start_shape[1], start_shape[1] + search_extent),
                               range(start_shape[2], start_shape[2] + search_extent)
                               )

    list_of_shapes = list(map(list, gen_of_shapes))

    # Make a dummy grid instance, and set its spacegroup
    gemmigrid = gemmi.Int8Grid()
    #gemmigrid = gemmi.FloatGrid() # It should not matter at all which grid-type
    gemmigrid.spacegroup = gemmi_spgrp

    compatible_shapes = []
    incompat_shapes = []

    for shp in list_of_shapes:
        if is_shape_compatible(shp, gemmigrid):
            compatible_shapes.append(shp)
        else:
            incompat_shapes.append(shp)

    return compatible_shapes, incompat_shapes


# TODO: DONE, but test it!
def find_grid_shape_closest_to_spacing(grid_shapes, ref_spacing, a, b, c):
    """Finds the grid shape(s) in a list that is closest 
    to a certain reference spacing.
    grid_shapes : list of lists of grid shapes
    a,b,c : cell parameters of the unitcell.
    
    OBS: DONE but should be tested."""

    assert isinstance(ref_spacing, (tuple, list)) and len(ref_spacing) == 3, "Error: ref_spacing must be tuple or list of length 3."
    assert len(np.array(grid_shapes).shape) == 2 and np.array(grid_shapes).shape[1] == 3, "Error: grid_shapes needs to be of shape (n, 3)."

    grid_spacings = [(a/s[0], b/s[1], c/s[2]) for s in grid_shapes]

    grid_spacing_norms = np.linalg.norm(np.subtract(grid_spacings, np.array(ref_spacing)).reshape(1,3), axis=1)
    best_spacing = np.min(grid_spacing_norms)
    best_indx = np.nonzero(grid_spacing_norms == best_spacing)

    closest_shape = grid_shapes[best_indx[0][0]]

    return closest_shape


# TODO: 1. Finish implementation. 2. Test it! 3. Check usage in cli scripts.
def find_spacegroup_compatible_gridshape(grid_shape:tuple, gemmi_spacegroup:gemmi.SpaceGroup,
                                         a:float, b:float, c:float, ref_spacing=None, search_extent=10, mode="up"):
    """Finds a grid shape compatible with a given spacegroup setting. 
    Based on a given, initial grid shape, and a spacegroup entry, find a grid shape
    that is compatible with the spacegroup symmetry, and has a spacing closest to
    a given reference spacing (by default hat of the initial grid shape)."""

    # Start of strong with NotImplemented check
    if not mode == "up":
        raise NotImplementedError("ERROR: mode = {} is not implemented yet. Only available option is \"up\".")
    
    # Type checks
    assert isinstance(grid_shape, tuple) and len(grid_shape) == 3 and all([isinstance(x, int) for x in grid_shape]), "ERROR: grid_shape must be len=3 tuple of integers!"
    assert isinstance(gemmi_spacegroup, gemmi.SpaceGroup), "ERROR: gemmi_spacegroup must be gemmi.SpaceGroup object!"
    assert ref_spacing is None or (isinstance(ref_spacing, (tuple, list)) and len(ref_spacing) == 3 and all([isinstance(s, (float, int)) for s in ref_spacing])), "ERROR: ref_spacing must be len=3 tuple or list of floats!"

    # Now search for list of compatible shapes
    compatible_shapes, incompat_shapes = search_compatible_grid_shapes(gemmi_spgrp=gemmi_spacegroup,
                                                                       start_shape=grid_shape,
                                                                       search_extent=search_extent,
                                                                       mode=mode)
    
    ### TODO: Could implement this. Not neccessary
    # Extend and iterate search if nothing is found:
    #max_while_iter = 5
    #while_iter = 0
    #bool_stop = False
    #start_shape = grid_shape
    #while compatible_shapes == [] and not bool_stop:
    #    new_start_shape = (i + search_extent for i in start_shape,)
    #    start_shape = new_start_shape

    #    compatible_shapes, incompat_shapes = search_compatible_grid_shapes(gemmi_spgrp=gemmi_spacegroup,
    #                                                                   start_shape=start_shape,
    #                                                                   search_extent=search_extent,
    #                                                                   mode=mode)
    #    while_iter += 1
    #    if while_iter == max_while_iter:
    #        bool_stop = True
    ###

    if compatible_shapes == []:
        raise RuntimeError("ERROR: Search for compatible grid shape failed; no shape found.") # max iterations exceeded: while_iter = {}, max_while_iter = {}".format([while_iter, max_while_iter]))


    # Now, make selection of best shape based on reference spacing
    if ref_spacing is None:
        # Default is spacing corresponding to input grid shape
        ref_spacing = [a/grid_shape[0], b/grid_shape[1], c/grid_shape[2]]
    
    closest_compatile_shape = find_grid_shape_closest_to_spacing(compatible_shapes,
                                                                 ref_spacing=ref_spacing,
                                                                 a=a, b=b, c=c)
    
    return closest_compatile_shape



#TODO: Deprecate this?
def is_shape_and_spacegroup_compatible(shape, spacegroup):
    """ Tests whether the grid shape and spacegroup are compatible according
        to gemmi's grid.set_size() function.
    """
    raise DeprecationWarning("This function has been deprecated. Please use is_shape_compatible instead!")
    test_grid = gemmi.FloatGrid()
    test_grid.spacegroup = spacegroup

    try:
       test_grid.set_size(*shape)
    except:
       return False
    else:
       return True 


# TODO: deprecate this _old_
def _old_find_spacegroup_compatible_gridshape(grid_shape, spacegroup,
                                    search_denser=True):
    """ Starting from a given grid shape, searches for the closest grid
        shape compatible with the given spacegroup. THIS IS THE OLD VERISON.
        IT IS DEPRECATED IN FAVOR OF A NEW, IMPROVED, IN DEVELOPMENT.
    """
    raise DeprecationWarning("WARNING: This function is deprecated.")
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
    # Conversion into int-type
    final_compatible_grid_shape = (int(compatible_grid_shape[0]),
                                   int(compatible_grid_shape[1]),
                                   int(compatible_grid_shape[2]))

    # And check that it is indeed the same
    assert (np.array(final_compatible_grid_shape) 
            == np.array(compatible_grid_shape)).all(), "ERROR: Compatible grid shapes not equal after type conversion!"

    return final_compatible_grid_shape


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


def get_conv_std_structure_pymatgen(atoms):
        sg_analyzer = SpacegroupAnalyzer(AseAtomsAdaptor.get_structure(atoms))
        conv_cell = sg_analyzer.get_conventional_standard_structure()
        conv_std_atoms = AseAtomsAdaptor.get_atoms(conv_cell)
        return conv_std_atoms


def is_close_atoms(atoms1:ase.Atoms, atoms2:ase.Atoms, detailed=False):
    """Makes the equality check for atoms, which consists of 
    positions, numbers, cell, and periodic boundary conditions,
    but applies np.isclose instead of strict equality for positions and
    cell, to make small differences negligible.
    Also makes the np.isclose-check symmetric by construction, i.e.
    defines it as np.isclose(a,b) AND np.isclose(b,a)."""
    eq_positions = (np.isclose(atoms1.positions, atoms2.positions).all()
                    and np.isclose(atoms2.positions, atoms1.positions).all())
    eq_numbers = (atoms1.numbers == atoms2.numbers).all()
    eq_cells = (np.isclose(atoms1.cell, atoms2.cell).all() 
                and np.isclose(atoms2.cell, atoms1.cell).all())
    eq_pbc = (atoms1.pbc == atoms2.pbc).all()

    return (eq_positions and eq_numbers and eq_cells and eq_pbc)


def atoms_to_spglib_cell_tuple(atoms:ase.Atoms):
    """Converts an atoms object into a cell tuple as used in the 
    input to spglib, consisting of a tuple of the cartesian cell matrix
    in row vector form, the list of scaled coordinates (i.e. fractional) of the
    atoms, and the list of atom number in the same order as the coordinates.
    Only the cell matrix and the position and type of the atoms are tranferred.
    Any other informaiton in the atoms object are not carried over."""
    cell = np.array(atoms.get_cell())
    positions = atoms.get_scaled_positions()
    numbers = atoms.get_atomic_numbers()
    return (cell, positions, numbers)


# Test this!
def spglib_cell_tuple_to_atoms(cell_tuple, pbc=True, **kwargs):
    """Converts a spglib cell tuple, i.e. tuple of a cartesian cell
    matrix in row vector form, a list of scaled coordinates, and list
    of atom numbers, into an atoms object.
    Can pass other information that the ase.Atoms object constructor
    can take as keyword arguments as well, although these will be external
    additions to the infomation in the cell tuple. Only the pbc is required
    to be specified (it is True by default)."""
    atoms = ase.Atoms(cell=cell_tuple[0],
                      scaled_positions=cell_tuple[1],
                      numbers=cell_tuple[2],
                      pbc=pbc,
                      **kwargs,
                      )
    return atoms


# NOTE/TODO: Test This. The old basic seems to work. Now inlcuded clean/wrap options
def spglib_dataset_to_symops_list(spglib_dataset, clean=False, wrap=False):
    """Utility function that extracts the symops in a spglib symmetry dataset
    and returns it as a list of tuples given by list(zip(rotations, translations)).
    If clean is set to True, then the translation parts of the symops are wrapped and 
    cleaned (e.g. values close to 0.0 and 1.0 will be set to 0.0.), before being returned.
    If wrap is set to True, then the translations are wrapped before being returned."""
    assert spglib_dataset is not None, "ERROR: spglib dataset is None."
    rots = spglib_dataset["rotations"]
    trans = spglib_dataset["translations"]
    symops = list(zip(rots,trans))

    if clean:
        # if clean, then return symops with cleaned translations
        return clean_symops_translations(symops)
    if wrap:
        # if wrap, then return symops with wrapped translations
        return wrap_symops_translations(symops)
    
    return symops


# NOTE: This is in row form, both rot and P!!!!!!
## NOTE THEre was an error in the order. Fixed it now. Check,.
def rot_change_of_basis_row_form(rot, P):
    """Changes the basis of a rotation matrix e.g. of a symmetry operation.
    The P matrix encodes the change of basis through x_2 = x_1 * P + p,
    where x_1 is a point represented in the current basis and x_2 is its
    representation in the new basis.
    This function assumes both the basis change matrix P AND the rotation
    matrix to be in row form, meaning it acts on coordinate vectors x in row
    format as: x * P + p, and x * rot."""
    rot = np.array(rot)
    P = np.array(P)

    new_rot = np.linalg.inv(P) @ rot @ P

    return new_rot


# NOTE: YOU NEED TO CHECK WHETHER THIS FORM ASSUMES P in row or column form. 
#       AND IF IT ASSUMES THAT THE SYMMETRY OPERATION IS IN ROW OR COLUMN FORM.
def symop_change_of_basis_row_form(rot, trans, P, p):
    """Changes the basis forward of a symmetry operation represented by
    a tuple of a rotation matrix and translation vector.
    The basis change is of the form x_2 = x_1*P + p, where this
    maps a point x_1 represented in the current basis to its
    representation in the new basis, x_2.
    The symop is transformed to be applicable in the new basis,
    that is, if x_1 and x_1~ are related through this symop in the
    current basis, x_2 = x_1 * P + p and x_2~ = x_1~ * P  + p
    will be related in the new basis. 
    Note that this function assumes both P as shown above and the
    rotation matrix of the symmetry operation rot to be in row form, i.e.
    act on row coordinate vectors as P does before, e.g. for symmetry
    operation (rot, trans) with rotation matrix rot and translation vector
    trans:
    y = x * rot + trans."""
    rot = np.array(rot)
    trans = np.array(trans)
    P = np.array(P)
    p = np.array(p)

    new_rot = rot_change_of_basis_row_form(rot, P)
    
    #P_inv = np.linalg.inv(P)
    new_trans = trans @ P + p @ (np.eye(3) - new_rot)
    
    return (new_rot, new_trans)


### NOTE: NEW! This is in column form!!!
### This is fixed and should be correct.
def rot_change_of_basis(rot, P):
    """Changes the basis of a rotation matrix e.g. of a symmetry operation.
    The P matrix encodes the change of basis through x_2.T = P*x_1.T + p.T,
    where x_1 is a point represented in the current basis and x_2 is its
    representation in the new basis. Here, note that x_2.T, x_1.T and p.T are
    column vectors (they are transposed row vectors as it is assumed coordinates
    and translation vectors are in row form).
    This function assumes both the basis change matrix P AND the rotation
    matrix to be in column form, meaning it acts on coordinate row vectors x
    in column format as:
    P * x.T + p.T, and rot * x.T."""
    rot = np.array(rot)
    P = np.array(P)

    new_rot = P @ rot @ np.linalg.inv(P)

    return new_rot


# NOTE: New! Test.
def symop_change_of_basis(rot, trans, P, p):
    """Changes the basis forward of a symmetry operation represented by
    a tuple of a rotation matrix and translation vector.
    The basis change is of the form x_2.T = P*x_1.T + p.T, where this
    maps a point x_1 represented in the current basis to its
    representation in the new basis, x_2. Here, note that x_2.T, x_1.T,
    and p.T are column vectors (they are transposed row vectors as it 
    is assumed coordinates and translation vectors are in row form).
    The symop is transformed to be applicable in the new basis,
    that is, if x_1 and x_1~ are related through this symop in the
    current basis, x_2.T = P*x_1.T + p.T and x_2~.T = P * x_1~.T + p.T
    will be related in the new basis. 
    Note that this function assumes both P as shown above and the
    rotation matrix of the symmetry operation rot to be in column form, i.e.
    act on row coordinate vectors as P does before, e.g. for symmetry
    operation (rot, trans) with rotation matrix rot and translation vector
    trans:
    y.T = rot * x.T + trans.T."""
    rot = np.array(rot)
    trans = np.array(trans)
    P = np.array(P)
    p = np.array(p)

    new_rot = rot_change_of_basis(rot, P)
    
    #P_inv = np.linalg.inv(P)
    new_trans = trans @ P.T + p @ (np.eye(3) - new_rot.T) # New, column form of matrices, but translation in row form!
    
    return (new_rot, new_trans)


# NOTE: NEW! TODO: TEST!
def change_basis_cell_tuple(cell_tuple, P, p, clean=False, wrap=False):
    """Changes the basis of a spglib cell tuple (lattice, positions, numbers),
    according to a change of basis and origin shift given by P and p, e.g.
    coordinates are transformed according to x_2 = x_1 * P.T + p,
    and cell matrix is tranformed as C_2 = (P.T)^-1 * C_1.
    This is equivalent to how the cell matrix (lattice) and positions are
    changed in standardization in spglib, from the transformation matrix P and
    origin shift p."""
    cell1 = np.array(cell_tuple[0])
    pos1 = np.array(cell_tuple[1])
    num = cell_tuple[2]

    cell2 = np.linalg.solve(P.T, cell1)
    tmp_pos2 = pos1 @ P.T + p

    # Assert that the neccessaty equations are fullfilled by the found solution
    # Do symmetric allclose comparison
    assert np.allclose(cell1, P.T @ cell2) and np.allclose(P.T @ cell2, cell1)
    assert np.allclose(pos1 @ cell1, (tmp_pos2 - p) @ cell2) and np.allclose((tmp_pos2 - p) @ cell2, pos1 @ cell1)
    
    # Test basis change matrix, assert that determinant and metric tensor are close to unity
    #assert np.linalg.det(P) == 1.0
    #assert P@P.T == np.eye(3)
    # Same but symmetric allclose
    #assert np.allclose(np.linalg.det(P), 1.0, atol=10*-7, rtol=10**-7) and np.allclose(1.0, np.linalg.det(P), atol=10*-7, rtol=10**-7)
    #assert np.allclose(P @ P.T, np.eye(3), atol=10*-7, rtol=10**-7) and np.allclose(np.eye(3), P @ P.T, atol=10*-7, rtol=10**-7)


    if clean:
        pos2 = clean_translation(tmp_pos2)
    elif wrap:
        pos2 = wrap_translation(tmp_pos2)
    else:
        pos2 = tmp_pos2

    new_cell_tuple = (cell2, pos2, num)
    return new_cell_tuple



def symop_to_gemmi_op(rot, trans):
    """Transforms a symmetry operation representated by a
    tuple of a rotation matrix and a translation vector, into
    a gemmi Op symmetry operation object.
    Asserts that the inverse operation is successful.
    Assumes rot = 3x3 rotation matrix, trans = len 3 translation vector."""
    
    # Convert rot and trans into integer form understood by gemmi, with DEN(ominator) = 24
    rot_int_form = np.rint(np.array(rot)*24).astype(int)
    trans_int_form = np.rint(np.array(trans)*24).astype(int)

    # Convert into lists, from arrays (shape of trans?):
    gemmi_rot = list(map(list, rot_int_form))
    gemmi_trans = list(trans_int_form.flatten()) 

    # Instantiate gemmi symop and 
    # then set its rotation and translation parts:
    gemmi_sym_op = gemmi.Op()

    gemmi_sym_op.rot = gemmi_rot
    gemmi_sym_op.tran = gemmi_trans

    assert np.isclose(np.array(gemmi_sym_op.rot)/24.0, np.array(rot)).all(), "ERROR: Found gemmi rotation: "+str(np.array(gemmi_sym_op.rot)/24) + " is not similar to input: " + str(rot)
    assert np.isclose(np.array(gemmi_sym_op.tran)/24.0, np.array(trans)).all(), "ERROR: Found gemmi translation: " + str(np.array(gemmi_sym_op.tran)/24) + " is not similar to input: " + str(trans)

    return gemmi_sym_op


### Test gemmi spacegroup things:
def test_gemmi_spacegroup_to_ops_functions(spgrp:gemmi.SpaceGroup):
    trans = gemmi_spacegroup_to_translation_list(spgrp)
    rot = gemmi_spacegroup_to_rotations_list(spgrp)
    full_ops1 = list(zip(rot, trans))
    full_ops2 = gemmi_spacegroup_to_symops_list(spgrp)
    #print(full_ops1)
    #print(full_ops2)
    not_equal = False
    for i,op in enumerate(full_ops1):
        r1, t1 = op
        r2, t2 = full_ops2[i][0], full_ops2[i][1]
        if not ((r1 == r2).all() and (t1 == t2).all()):
            print("ERROR: symops not equivalent???")
            no_equal = True
    print("")   
    if not_equal:
        print("List of spacegroups are NOT the same!")
    else:
        print("List of spacegroups are THE SAME! YAY")


# NOTE: Test changes to these three below from removing .sym_ops. Symops does not include pure translations! This is the  root cause of the different number of symops error!
# NOTE: Done, has updated to get all ops, not just sym_ops.
def gemmi_spacegroup_to_rotations_list(spgrp:gemmi.SpaceGroup):
    """Returns a list of rotation-parts of of the symmetry operations
    of a gemmi.SpaceGroup object."""
    #rot_list = [np.array(op.rot)/24 for op in spgrp.operations().sym_ops] # Old!
    rot_list = [np.array(op.rot)/24 for op in spgrp.operations()] # New!
    return rot_list


# NOTE: Done, has updated to get all ops, not just sym_ops.
def gemmi_spacegroup_to_translation_list(spgrp:gemmi.SpaceGroup):
    """Returns a list of translation-parts of of the symmetry operations
    of a gemmi.SpaceGroup object."""
    #trans_list = [np.array(op.tran)/24 for op in spgrp.operations().sym_ops] # Old!
    trans_list = [np.array(op.tran)/24 for op in spgrp.operations()] # New!
    return trans_list


# NOTE: Done, has updated to get all ops, not just sym_ops.
def gemmi_spacegroup_to_symops_list(spgrp:gemmi.SpaceGroup):
    """Returns list of symmetry operations in a gemmi.SpaceGroup object,
    in the form of a list of tuples (rotation, translation).
    Should be equivalent to list(zip(rotation_list, translation_list)).
    """
    #symoplist = [(np.array(op.rot)/24, np.array(op.tran)/24) for op in spgrp.operations().sym_ops] # Old!
    symoplist = [(np.array(op.rot)/24, np.array(op.tran)/24) for op in spgrp.operations()] # New!

    return symoplist


# NOTE: This needs to be changed, either the name or the functionality.
# as of now it is not symmetric, and it is easy to see that if symops1 is
# a subset of symops2, then all symops in symops1 will be found in symops2.
# Hence the eq_symops will be equal to symops1, since all elements will be found,
# and no elements will be added to diff_symops. Yet there might be elements in 
# symops2 that are not inlcuded in symops1.
# This function really checks whether symops1 is subset of symops2!!!
# THIS IS IMPORTANT BECAUSE THIS COMPLETELY AFFECTS THE ACTUAL FUNCITONALITY OF are_symops_equal, AND SHOULD BE ADDRESSED IMMEDIATELY!

# NOTE: Update to above note: Functions below (compare_symops, are_symops_equal have been updated (names changed)
#       and two new functions have been introduces. New compare symops has been written, and are_symops_equal not relies on this new version.
#       Compare symops now outputs three lists and hsould be symmetric rather than assymetric. 
#       symops_inclusion_partition is the old compare_symops, which was assymetric in the comparison.
#       The occurence of the old function names (compare_symops mainly) need to be assessed and possibly
#       altered in the other scripts. The functions also need to be tested.

# NOTE: Name has been altered (this was compare_symops before) and has added additional assertion
def symops_inclusion_partition(symops1, symops2):
    """Compares two symops-lists and partitions the first
    into two disjoint lists: one list of operations that
    are also found in the second symops-list,
    and one list of the rest of the operations, e.g. those operations not found
    in the second symops-list.
    np.isclose is used to compare symops."""
    incl_symops = []
    diff_symops = []

    for rot1, trans1 in symops1:
        was_found = False
        for rot2, trans2 in symops2:
            # Test with all close and with lower absolute tolerances
            if np.allclose(rot1, rot2, atol=10**-6) and np.allclose(trans1, trans2, atol=10**-6):
                incl_symops.append((rot1, trans1))
                was_found = True
                break

        if not was_found:
            diff_symops.append((rot1, trans1))
    
    # Sanity check: Make sure lengths add up correctly
    assert len(incl_symops) + len(diff_symops) == len(symops1)

    return incl_symops, diff_symops

# NOTE: New, need testing.
def are_symops_contained_in(symops1, symops2):
    """Checks if symops1 are subset of symops2, i.e. if all symmetry
    operations in the symops1 list are also in the symops2 list. 
    Returns True if symops1 is subset of symops2, and False otherwise."""
    incl_list, excl_list = symops_inclusion_partition(symops1, symops2)
    if incl_list == symops1 and excl_list == []:
        return True
    else:
        return False

# NOTE: Has been changed. Need to be tested. Need to find it and alter it if necessary in in other scripts
def compare_symops(symops1, symops2):
    """Compares a pair of symops-lists and returns three lists.
    The first containing the common operations as found with symops_inclusion_partition,
    the second with the operations exclusively found in 1, and the third with the
    operations exclusively found in 2."""
    symops_in_both_1, symops_only_in_1 = symops_inclusion_partition(symops1, symops2)
    symops_in_both_2, symops_only_in_2 = symops_inclusion_partition(symops2, symops1)
    
    # Sanity check 1: Make sure that the operations found in both
    # are the same in both comparisons
    # This can actually be done with two calls of symops_inclusion_partition
    #assert len(symops_in_both_1) == len(symops_in_both_2)
    if not len(symops_in_both_1) == len(symops_in_both_2):
        # Debug
        print("WARNING: len of symops list in both are not of the same length - this should not have happened.")
        print("Symops in both 1, len = ",len(symops_in_both_1),":")
        print(*symops_in_both_1, sep="\n")
        print("Symops in both 2, len = ",len(symops_in_both_2),":")
        print(*symops_in_both_2, sep="\n")
    #assert are_symops_contained_in(symops_in_both_1, symops_in_both_2) and are_symops_contained_in(symops_in_both_2, symops_in_both_1)

    # Sanity check 2: Make sure the lengths add up - this is done within the symops_inclusion_partition function
    # This is probably redundant
    #assert len(symops1) - len(symops2) == len(symops_only_in_1) - len(symops_only_in_2)

    return symops_in_both_1, symops_only_in_1, symops_only_in_2


# NOTE: This has been changed. Seems to work fine.
def are_symops_equal(symops1, symops2):
    eqlist, difflist1, difflist2 = compare_symops(symops1, symops2)
    if eqlist == symops1 and difflist1 == [] and difflist2 == []:
        return True
    else:
        return False
    

# NOTE: Seems to work fine
def get_all_gemmi_spgrp_entries_by_spgrp_number(number:int):
    """Fetches all gemmi.SpaceGroup entries from the gemmi spacegroup table
    which matches the given spacegroup number, and returns them as a list."""
    #gemmi_spgrp_table = list(list(gemmi.spacegroup_table_itb())[0:530]) # First 530 entries - although we might just grab the whole thing ...
    gemmi_spgrp_table = list(gemmi.spacegroup_table_itb()) # The whole list
    gemmi_matching_entries = [SpGrpObj for SpGrpObj in gemmi_spgrp_table if SpGrpObj.number == number]
    return gemmi_matching_entries


# NOTE: Seems to work fine
def get_all_spglib_spgrp_entries_by_spgrp_number(number:int):
    """Fetches all spacegroup symmetry dataset entries from the spglib spacegroup list
    (the listing according to hall numbers) which matches the given spacegroup number,
    and returns them as a list."""
    spglib_table_symdicts = [spglib.get_spacegroup_type(i) for i in range(1,531)]
    spglib_matching_entries = [spgrp_dict for spgrp_dict in spglib_table_symdicts if spgrp_dict["number"] == number]
    return spglib_matching_entries


# NOTE/TODO: NEW! Test it maybe, or maybe it is not really needed. Whaetver is more convenient.
#            I think a better approach is to inlcude this in the original find and match function.
def find_matching_gemmi_spacegroup_from_symops(symops, spacegroup_number=None, match_type="full", clean=True, wrap=False):
    """Searches the gemmi spacegroup table and find a match with the same symmetry
    operations as the given list of symops. This function employs a full search, unless spacegroup_number
    is given, in which case it only looks up entries matching this number.
    The wrap option turns on wrapping of the translational parts in the symops list,
    before matching. Similarily, the clean option turns on wrapping and cleaning up of
    translation vector in the symops list, e.g. values numerically close to 0.0 and 1.0 are
    set to 0.0. """
    # Clean or wrap according to settings:
    if clean:
        cmp_symops = clean_symops_translations(symops)
    elif wrap:
        cmp_symops = wrap_symops_translations(symops)
    else:
        cmp_symops = symops

    # Set method for spacegorup matching (criteria)
    if match_type == "full":
        # Considered a match if gemmi spacegroup has ALL symops of spglib dataset
        spgrp_match_function = are_symops_equal
    elif match_type == "inclusion":
        # considered a match if gemmi spacegroup is a SUBSET of symops of the spglib dataset
        spgrp_match_function = are_symops_contained_in
    else:
        raise ValueError("Unsupported match_type option. match_type must be either \"full\"or \"inclusion\".")
    
    # Get all candidates with this number
    if spacegroup_number:
        gemmi_sg_candidates =  get_all_gemmi_spgrp_entries_by_spgrp_number(spacegroup_number)
    else:
        gemmi_sg_candidates =  list(gemmi.spacegroup_table_itb())
    
    # Loop through the candidates, and try symops matching
    # To find all spacegroups that have matching symops
    matching_sgs = []
    for sg_candidate in gemmi_sg_candidates:
        tmp_gemmi_symops = gemmi_spacegroup_to_symops_list(sg_candidate)
        if spgrp_match_function(tmp_gemmi_symops, cmp_symops):
            matching_sgs.append(sg_candidate)
    

    return matching_sgs



# NOTE: Seems to work well. TODO: Update by exctracting match function from symops.
def find_matching_gemmi_spacegroup_from_spglib_dataset(spglib_symmetry_dataset, match_type="full", clean=True, wrap=False):
    """Searches the gemmi spacegroup table and find a match with the same symmetry
    operations as the given spglib dataset."""
    # Set method for spacegorup matching (criteria)
    if match_type == "full":
        # Considered a match if gemmi spacegroup has ALL symops of spglib dataset
        spgrp_match_function = are_symops_equal
    elif match_type == "inclusion":
        # considered a match if gemmi spacegroup is a SUBSET of symops of the spglib dataset
        spgrp_match_function = are_symops_contained_in
    else:
        raise ValueError("Unsupported match_type option. match_type must be either \"full\"or \"inclusion\".")
    
    # Get spacegroup number
    sg_number = spglib_symmetry_dataset["number"]
    spglib_symops = spglib_dataset_to_symops_list(spglib_symmetry_dataset, clean=clean, wrap=wrap)

    # Get all candidates with this number
    gemmi_sg_candidates =  get_all_gemmi_spgrp_entries_by_spgrp_number(sg_number)

    # Loop through the candidates, and try symops matching
    # To find all spacegroups that have matching symops
    matching_sgs = []
    for sg_candidate in gemmi_sg_candidates:
        tmp_gemmi_symops = gemmi_spacegroup_to_symops_list(sg_candidate)
        if spgrp_match_function(tmp_gemmi_symops, spglib_symops):
            matching_sgs.append(sg_candidate)
    

    return matching_sgs


def brute_force_matching_gemmi_spacegroup_from_spglib_dataset(spglib_symmetry_dataset):
    # Get spacegroup symops
    spglib_symops = spglib_dataset_to_symops_list(spglib_symmetry_dataset)

    # Get all candidates as list: This is all spacegroup entries in this case
    gemmi_sg_candidates = list(gemmi.spacegroup_table_itb())

    # Loop through the candidates, and try symops matching
    # To find all spacegroups that have mahcing symops
    matching_sgs = []
    for sg_candidate in gemmi_sg_candidates:
        tmp_gemmi_symops = gemmi_spacegroup_to_symops_list(sg_candidate)
        if are_symops_equal(spglib_symops, tmp_gemmi_symops):
            matching_sgs.append(sg_candidate)
    

    return matching_sgs


# TODO: Test this
def spglib_dataset_to_gemmi_spacegroup_by_ops(spglib_dataset:dict):
    """Constructs a gemmi.SpaceGroup object from the symmetry operations
    from a spglib symmetry dataset."""
    symop_list =  spglib_dataset_to_symops_list(spglib_dataset)
    gemmi_ops_list = [symop_to_gemmi_op(op[0], op[1]) for op in symop_list]
    gemmi_grpops = gemmi.GroupOps(gemmi_ops_list)
    gemmi_spgrp_from_ops = gemmi.find_spacegroup_by_ops(gemmi_grpops)
    return gemmi_spgrp_from_ops


# TODO: Test this
def spglib_dataset_to_gemmi_spacegroup_by_hall(spglib_dataset:dict):
    """Constructs a gemmi.SpaceGroup object from the Hall symbol
    from a spglib symmetry dataset."""
    gemmi_spgrp_from_hall = gemmi.find_spacegroup_by_ops(gemmi.symops_from_hall(spglib_dataset["hall"]))
    return gemmi_spgrp_from_hall


# TODO:Test this - This seems to work now. Corrected the spglib_dataset key to the correct one: 'hall'
def check_spacegroup_match_between_spglib_and_gemmi(spglib_dataset:dict, gemmi_spgrp:gemmi.SpaceGroup, return_list=False, print_list=False):
    """Compares an spglib symmetry dataset, which contains spacegroup information,
    with a gemmi.SpaceGroup object. This is done by passing the Hall symbol from the
    spglib dataset to gemmi, which generates a group of symmetry operations from this,
    and then constructs a corresponding spacegroup object.
    The properties of this spacegroup object derived from the spglib dataset, are compared
    with the given input gemmi SpaceGroup. The comparison consists of:
        Direct "==" comparison of the Spacegroup Objects
        list of operations,
        Hall symbol,
        Spacegroup number,
        Hermann-Mauguin symbol,
        Extended Herman-Mauguin symbol.
    """
    # First check that hall symbol - derived spacegroups fully match:
    # Perhaps a bit redundant...
    gemmi_spgrp_from_hall = gemmi.find_spacegroup_by_ops(gemmi.symops_from_hall(spglib_dataset["hall"]))
    hall_comparison_list = (
               gemmi_spgrp_from_hall == gemmi_spgrp,
               gemmi_spgrp_from_hall.operations() == gemmi_spgrp.operations(),
               gemmi_spgrp_from_hall.hall == gemmi_spgrp.hall,
               gemmi_spgrp_from_hall.number == gemmi_spgrp.number,
               gemmi_spgrp_from_hall.hm == gemmi_spgrp.hm,
               gemmi_spgrp_from_hall.xhm() == gemmi_spgrp.xhm()
                           )
    hall_comparison = (
               gemmi_spgrp_from_hall == gemmi_spgrp
               and gemmi_spgrp_from_hall.operations() == gemmi_spgrp.operations()
               and gemmi_spgrp_from_hall.hall == gemmi_spgrp.hall
               and gemmi_spgrp_from_hall.number == gemmi_spgrp.number
               and gemmi_spgrp_from_hall.hm == gemmi_spgrp.hm
               and gemmi_spgrp_from_hall.xhm() == gemmi_spgrp.xhm()
                           )
    
    assert hall_comparison == all(hall_comparison_list)

    if print_list:
        test_labels = ["Spacegroup equality: ",
                       "Operations: ",
                       "Hall symbol: ",
                       "Spacegroup number: ",
                       "Hermann-Mauguin: ",
                       "Extended HM symbol:"
                       ]
        print(*list(zip(test_labels,hall_comparison_list)), sep="\n")

    if return_list:
        return hall_comparison_list
    else:
        return hall_comparison


# How to test this?
def change_basis_symops(symops:list, P:np.ndarray, p:np.ndarray):
    """Changes the basis of symmetry operations contained in list symops,
    of the form (rot, tran), by using the given basis change matrix P and
    translation vector p."""
    new_symops = []
    for rot, trans in symops:
        new_symops.append(symop_change_of_basis(rot, trans, P, p))
    return new_symops


def wrap_translation(tran:np.array):
    """Wraps a translation vector, e.g. reduces each component mod 1 to
    be 0.0 <= x < 1.0. Note, can result in values numerically close to 1.0."""
    wrapped_tran = np.copy(tran) % 1.0
    #wrapped_tran[wrapped_tran == 1.0] = 0.0 
    return wrapped_tran


def clean_translation(tran:np.array, **kwargs):
    """Wraps a translation vector, and then cleans up values numerically
    very close to 0.0 and 1.0 and sets these to identically 0.0.
    Identification of value \"close\" to 0.0 and 1.0 are done with
    np.isclose, with default atol and rtol values. These can be changed by their
    np.isclose keyword arguments, as **kwargs are passed."""
    wrapped_tran = wrap_translation(tran)
    where_to_clean = np.logical_or(np.isclose(wrapped_tran, 0.0, **kwargs), np.isclose(wrapped_tran, 1.0, **kwargs))
    clean_tran = np.copy(wrapped_tran)
    clean_tran[where_to_clean] = 0.0
    return clean_tran


def wrap_symops_translations(symops:list):
    wrapped_trans_symops = [(op[0], wrap_translation(op[1])) for op in symops]
    return wrapped_trans_symops


def clean_symops_translations(symops:list):
    cleaned_trans_symops = [(op[0], clean_translation(op[1])) for op in symops]
    return cleaned_trans_symops



# Suggested names:
# conform_structure_and_spacegroup
# prepare_structure_and_spacegroup
# match_structure_and_spacegroup
# standardize_structure_and_spacegroup
# Want function that takes atoms and returns atoms, and gemmi.SpaceGroup
# What the spacegroup finding functions should do is to standardize
# This should return a structure (either cell tuple or atoms) and a spacegroup that matches

# This is the simple/direct function
# TODO: Test this in sampling.
def match_cell_tuple_with_gemmi_spacegroup(cell_tuple, clean=True, wrap=False):
    """From an atoms object, attempt to find the matching spacegroup
    in gemmi, and return a gemmi.SpaceGroup object that are in accord with the atoms object
    as well as the symmetry dataset of this structure from spglib.
    The clean and wrap arguments are directly passed to the spacegroup matching function
    find_matching_gemmi_spacegroup_from_spglib_dataset that is used to find the corresponding 
    spacegroup in gemmi. Also returns the symmetry dataset from spglib.
    
    Default settings involves comparisons of the spacegroup operations after \"cleaning\" (clean=True)
    of the spglib operations, which means wrapping and reducing numbers very close to
    0.0 or 1.0 to 0.0 in the translational parts of the symmetry operations.
    """
    
    sym_dataset = spglib.get_symmetry_dataset(cell_tuple)

    # Here, try matching.
    try:
        matching_gemmi_spgrps = find_matching_gemmi_spacegroup_from_spglib_dataset(sym_dataset, match_type="full", clean=clean, wrap=wrap)
    except Exception as exc:
        raise Exception("Spacegroup matching failed!") from exc

    # Now we perform checks: Check that there is one and only one unique match - or if several, it is number 68, since it has multiple entries.
    if not len(matching_gemmi_spgrps) <= 1:
        # We only allow this ambiguity for spacegroup number 68, since it has duplicate entries 
        # in the gemmi table. Hence it is not ambiguous. Moreover, we sanity check for operation equality.
        for grp in matching_gemmi_spgrps:
            assert grp.number == 68, "ERROR: Ambiguous spacegroup matching for spacegroup other than 68. Found {} matching spacegroups. This should not have happened!".format(len(matching_gemmi_spgrps))
            assert matching_gemmi_spgrps[0].operations() == grp.operations(), "ERROR: Spacegroup matches (>1) have differing operations! This should not have happened."

    else:
        assert len(matching_gemmi_spgrps) == 1, "ERROR: Spacegroup matching failed! Found {} matching spacegroups.".format(len(matching_gemmi_spgrps))
    
    gemmi_spgrp = matching_gemmi_spgrps[0]
    
    return gemmi_spgrp, sym_dataset


# Not sure if needed
def match_atoms_with_gemmi_spacegroup(atoms:ase.Atoms, clean=True, wrap=False):
    cell_tuple = atoms_to_spglib_cell_tuple(atoms)

    gemmi_spgrp, sym_dataset = match_cell_tuple_with_gemmi_spacegroup(cell_tuple, clean=clean, wrap=wrap)

    # With these checks above, we are done, and can convert back to atoms object.
    matching_atoms = spglib_cell_tuple_to_atoms(cell_tuple)

    # This is a paranoia sanity check
    assert is_close_atoms(matching_atoms, atoms), "ERROR: Atoms from cell tuple somehow not close to input atoms. This should not have happened!"

    return cell_tuple, gemmi_spgrp, sym_dataset
    

# Standardization function - unclear if useful or not.
# TODO: Figure out if remove or use
def spglib_standardize_structure(atoms:ase.Atoms, to_primitive=False, no_idealize=True, change_basis=False):
    """Converts an atoms object into a spglib tuple and standardizes it using
    the spglib.standardize_cell function, according to the input settings.
    Arguments to_primitive and no_idealize are directly passed to this function.
    The change_basis flag controls whether an additional basis change should be 
    performed after the standardization with the spglib.standardize_cell function
    has been applied. The change of basis that can be applied is based on the 
    transformation matrix and origin shift of the spglib symmetry dataset of 
    the obtained cell tuple from the standardization. If the change of basis is
    preformad, the symmetry dataset is again determined for the resulting cell 
    tuple.

    Returns the final, standardized cell tuple, its symmetry dataset, a tuple of
    the transformation matrix P_mat and origin shift p_shift, and lastly a tuple of
    the symmetry datasets from spglib of the input stucture and from the
    structure obtained in the standardization before a potential change of basis.
    """
    cell_tuple = atoms_to_spglib_cell_tuple(atoms)
    orig_symmetry_dataset = spglib.get_symmetry_dataset(cell_tuple)

    tmp_std_cell_tuple = spglib.standardize_cell(cell_tuple, to_primitive=to_primitive, no_idealize=no_idealize)
    tmp_std_symmetry_dataset = spglib.get_symmetry_dataset(tmp_std_cell_tuple)

    P_mat = tmp_std_symmetry_dataset["transformation_matrix"]
    p_shift = tmp_std_symmetry_dataset["origin_shift"]


    if change_basis:
        std_cell_tuple = change_basis_cell_tuple(tmp_std_cell_tuple, P_mat, p_shift)
    else:
        std_cell_tuple = tmp_std_cell_tuple
    
    std_symmetry_dataset = spglib.get_symmetry_dataset(std_cell_tuple)

    
    return std_cell_tuple, std_symmetry_dataset, (P_mat, p_shift), (orig_symmetry_dataset, tmp_std_symmetry_dataset, std_symmetry_dataset)


def prepare_matching_structure_and_spacegroup(atoms:ase.Atoms, clean=True, wrap=False, to_primitive=False, no_idealize=True, spglib_standardize="allow", change_basis="allow", debug_output=False):
    """From an atoms object, attempt to find a matching spacegroup in the table in gemmi,
    and return a pair of a matching atoms object and a gemmi.SpaceGroup object that are
    in accord.
    Settings allow some amount of structural processing to be performed, either prior to
    attempting the spacegorup matching, or as needed, through spglib.

    The clean and wrap arguments are directly passed to the spacegroup matching function
    find_matching_gemmi_spacegroup_from_spglib_dataset that is used to find the corresponding 
    spacegroup in gemmi.

    The spglib_standardize argument take on one of the values "off", "on" or "allow", which
    specifies whether standardization through spglib should be done prior to matching ("on"),
    as needed if no spacegroup matches are found at first ("allow"), or be completely turned off.
    The change_basis argument similarily controls whether, after said standardization in spglib,
    a basis change should be applied to the structure. "on" specifies this as part of the preprocessing
    step, while "allow" enables it only if needed (e.g. no spacegroup mathces at first attempt), and 
    "off" disables it completely.

    The symmetry search and standardization are done with spglib. to_primitive and no_idealize
    are passed to the standardize_cell function, and controls how standardization is done if enabled.
    
    Default settings are spglib_standardize = "allow" and change_basis = "allow", i.e. an attempt is
    to match the structure as is, but if neccessary, a standardization is done with spglib, and if no
    match is found after this, a final matching is attempted after a change of basis.
    
    The default spglib standardizaton settings are standardization to the conventional unitcell 
    without idealization through spglibs standardize_cell-function with to_primitive=False 
    and no_idealize=True), as well as comparisons of the spacegroup operations after \"cleaning\"
    (clean=True) of the spglib operations, which means wrapping and reducing numbers very close
    to 0.0 or 1.0 to 0.0 in the translational parts of the symmetry operations.
    
    The change of basis that is applied is contained in the spglib dataset for the standardized structure, 
    e.g. the transformation matrix and the origin shift. This does not change the structure, 
    but merely alters the coordinate system and basis, and thus may alter the representation of 
    the symmetry operations."""
    
    # Input logic
    assert spglib_standardize in ["off", "allow", "on"], "ERROR: spglib_standardize argument was {}. It must be off, allow, or on.".format(spglib_standardize)
    assert change_basis in ["off", "allow", "on"], "ERROR: change_basis argument was {}. It must be off, allow, or on.".format(change_basis)

    if spglib_standardize == "allow":
        assert change_basis in ["off", "allow"], "ERROR: change_basis must be off or allow when spglib_standardize = {}".format(spglib_standardize)
    
    if spglib_standardize == "off":
        assert change_basis == "off", "ERROR: change_basis must be off when spglib_standardize = {}".format(spglib_standardize)
    # End input logic

    # Print settings, in case:
    print("spglib_standardize = ", spglib_standardize, "  and  ", "change_basis = ", change_basis)
    print("spglib standardization settings: \n", "to_primitive: ", to_primitive, "\nno_idealize: ", no_idealize)
    #

    # Here, set the initial settings for the structure:
    # These should always be found. Mostly in case of debugging,
    init_cell_tuple = atoms_to_spglib_cell_tuple(atoms)
    init_sym_dataset = spglib.get_symmetry_dataset(init_cell_tuple)

    std_cell_tuple = spglib.standardize_cell(init_cell_tuple, to_primitive=to_primitive, no_idealize=no_idealize)
    std_sym_dataset = spglib.get_symmetry_dataset(std_cell_tuple)

    P_mat, p_shift = std_sym_dataset["transformation_matrix"], std_sym_dataset["origin_shift"]
    cb_cell_tuple = change_basis_cell_tuple(std_cell_tuple,P=P_mat, p=p_shift, clean=clean, wrap=wrap)
    cb_sym_dataset = spglib.get_symmetry_dataset(cb_cell_tuple)
    #

    ### Initialization - setting of initial attempt structure ###
    if spglib_standardize == "on":
        print("Preprocessing structure by standardizing with spglib according to settings.")
        if change_basis == "on":
            print("Applying change of basis in structure preprocessing.")
            print("Change of basis info:")
            print("P_mat = ", P_mat)
            print("p_shift = ", p_shift)

            # Set the tuple and dataset to the change-of-basis ones.
            match_cell_tuple = cb_cell_tuple
            match_sym_dataset = cb_sym_dataset

        elif change_basis in ["off", "allow"]:
            print("No change of basis applied in structure preprocessing.")
            match_cell_tuple = std_cell_tuple
            match_cell_tuple = std_sym_dataset
        else:
            raise ValueError("Unsupported setting change_basis = {}".format(change_basis))
    
    elif spglib_standardize in ["off", "allow"]:
        print("No preprocessing of the structure is done. Matching as is.")
        match_cell_tuple = init_cell_tuple
        match_sym_dataset = init_sym_dataset
    else:
        raise ValueError("Unsupported setting spglib_standardize = {}".format(spglib_standardize))
    ### END OF INITIALIZATION ###

    # Attempt to match:
    try:
        matching_gemmi_spgrps = find_matching_gemmi_spacegroup_from_spglib_dataset(match_sym_dataset, match_type="full", clean=clean, wrap=wrap)
    except Exception as exc:
        raise Exception("Spacegroup matching failed!") from exc
    

    ### Now, if we do not find a match, try again if settings allow it!
    # First: Standardize if permitted and not already done.
    if matching_gemmi_spgrps == [] and spglib_standardize == "allow":
        print("No matching spacegroups found, and spglib_standardize = allow. Standardizing, and attempting to match again...")
        
        match_cell_tuple = std_cell_tuple
        match_sym_dataset = std_sym_dataset

        # Attempt to match:
        try:
            matching_gemmi_spgrps = find_matching_gemmi_spacegroup_from_spglib_dataset(match_sym_dataset, match_type="full", clean=clean, wrap=wrap)
        except Exception as exc:
            raise Exception("Spacegroup matching failed!") from exc

    # Second: Change basis if permitted and not already done.
    if matching_gemmi_spgrps == [] and change_basis == "allow":
        print("No matching spacegroups found, and change_basis = allow. Changing basis, and attempting to match again...")
        print("Change of basis info:")
        print("P_mat = ", P_mat)
        print("p_shift = ", p_shift)

        match_cell_tuple = cb_cell_tuple
        match_sym_dataset = cb_sym_dataset

        # Attempt to match:
        try:
            matching_gemmi_spgrps = find_matching_gemmi_spacegroup_from_spglib_dataset(match_sym_dataset, match_type="full", clean=clean, wrap=wrap)
        except Exception as exc:
            raise Exception("Spacegroup matching failed!") from exc
    
    
    ### Finally: Do final checks and controls of the spacegroup match results
    #            and finalize the structure-spacegroup match results.

    # Check that there is one and only one unique match - or if several, it is spacegroup number 68, since it has multiple entries.
    if not len(matching_gemmi_spgrps) <= 1:
        # We only allow this ambiguity for spacegroup number 68, since it has duplicate entries 
        # in the gemmi table. Hence it is not ambiguous. Moreover, we sanity check for operation equality.
        for grp in matching_gemmi_spgrps:
            assert grp.number == 68, "ERROR: Ambiguous spacegroup matching for spacegroup other than 68. Found {} matching spacegroups. This should not have happened!".format(len(matching_gemmi_spgrps))
            assert matching_gemmi_spgrps[0].operations() == grp.operations(), "ERROR: Spacegroup matches (>1) have differing operations! This should not have happened."

        #raise Warning("ERROR: Ambiguous spacegroup matching. Found {} matching spacegroups. This should not have happened!".format(len(matching_gemmi_spgrps)))
    else:
        assert len(matching_gemmi_spgrps) == 1, "ERROR: Spacegroup matching failed! Found {} matching spacegroups.".format(len(matching_gemmi_spgrps))
    
    ### End of spacegroup match results checks

    # With these checks above, we are done, and can convert back to atoms object.
    match_gemmi_spgrp = matching_gemmi_spgrps[0]
    match_atoms = spglib_cell_tuple_to_atoms(match_cell_tuple)

    # What do I want to check here in the end?
    # Currently, I do the conversion atoms -> spglib tuple -> atoms
    # If I do not make any processing to the atoms object, I would like to return a copy of of the original.
    # Options: If I check to see that no operations has been done, I could check
    # that tuples are the same, and atoms are the same (or close)

    if match_cell_tuple == init_cell_tuple:
        # Sanity check
        assert is_close_atoms(atoms, match_atoms), "ERROR: No processing of structure, but input atoms and match atoms differ!"
        # This need not be true, numerical differences can cause this to be false
        print("CHECK: Are atoms == match_atoms?: ", atoms == match_atoms)
    
    if is_close_atoms(atoms, match_atoms):
        # Retain the original atoms in this case
        print("Input atoms and match atoms are close. Retaining original atoms.")
        match_atoms = atoms.copy()
    
    if not debug_output:
        return match_atoms, match_gemmi_spgrp, (P_mat, p_shift)
    else:
        return match_atoms, match_gemmi_spgrp, (P_mat, p_shift), (init_sym_dataset, std_sym_dataset, cb_sym_dataset)
    



# This is the more intricate/advanced workflow function - This is an attempt at a first version.
# Shuold probably scrap this version, since I have another more easily readable version.
# TODO: Scrap this?
def _dev_new_prepare_matching_structure_and_spacegroup(atoms:ase.Atoms, clean=True, wrap=False, to_primitive=False, no_idealize=True, spglib_standardize="allow", change_basis="allow"):
    """From an atoms object, standardize the structure as needed to find the matching spacegroup
    in gemmi, and return a pair of a standardized atoms object and a gemmi.SpaceGroup object
    that are in accord.
    The clean and wrap arguments are directly passed to the spacegroup matching function
    find_matching_gemmi_spacegroup_from_spglib_dataset that is used to find the corresponding 
    spacegroup in gemmi.
    The symmetry search and standardization are done with spglib. to_primitive and no_idealize
    are passed to the standardize_cell function, and controls how standardization is done.
    If no matching spacegroup is found in gemmi, an attempt to change the basis of the structure
    is done based on the transformation matrix and origin shift information in the spglib symmetry dataset.
    This is allowed by default but can be disabled by setting change_basis="off".
    
    Default settings involves a standardization to the conventional unitcell without idealization through spglibs
    standardize_cell-function (with to_primitive=False and no_idealize=True), as well as comparisons of
    the spacegroup operations after \"cleaning\" (clean=True) of the spglib operations, which means wrapping
    and reducing numbers very close to 0.0 or 1.0 to 0.0 in the translational parts of the symmetry operations.
    Moreover, if a match is not found for this standardized structure, a change of basis is performed and the
    matching is attempted again (allow_change_basis=True).
    The change of basis that is applied is contained in the spglib dataset for this structure, 
    e.g. the transformation matrix and the origin shift. This does not change the structure, 
    but merely alters the coordinate system and basis, and thus may alter the representation of 
    the symmetry operations."""
    raise NotImplementedError("ERROR: This function is not completely implemented and not for use. It is awaiting probable deprecation.")
    # Input logic
    assert spglib_standardize in ["off", "allow", "on"], "ERROR: spglib_standardize argument was {}. It must be off, allow, or on.".format(spglib_standardize)
    assert change_basis in ["off", "allow", "on"], "ERROR: change_basis argument was {}. It must be off, allow, or on.".format(change_basis)

    if spglib_standardize == "allow":
        assert change_basis in ["off", "allow"], "ERROR: change_basis must be off or allow when spglib_standardize = {}".format(spglib_standardize)
    
    if spglib_standardize == "off":
        assert change_basis == "off", "ERROR: change_basis must be off when spglib_standardize = {}".format(spglib_standardize)
    # End input logic

    # Print settings, in case:
    print("spglib_standardize = ", spglib_standardize, "  and  ", "change_basis = ", change_basis)
    print("spglib standardization settings: \n", "to_primitive: ", to_primitive, "\nno_idealize: ", no_idealize)
    #

    # Here, set the initial settings for the structure:
    # These should always be found. Mostly in case of debugging,
    init_cell_tuple = atoms_to_spglib_cell_tuple(atoms)
    init_sym_dataset = spglib.get_symmetry_dataset(init_cell_tuple)

    std_cell_tuple = spglib.standardize_cell(init_cell_tuple, to_primitive=to_primitive, no_idealize=no_idealize)
    std_sym_dataset = spglib.get_symmetry_dataset(std_cell_tuple)

    P_mat, p_shift = std_sym_dataset["transformation_matrix"], std_sym_dataset["origin_shift"]
    cb_cell_tuple = change_basis_cell_tuple(std_cell_tuple,P=P_mat, p=p_shift, clean=clean, wrap=wrap)
    cb_sym_dataset = spglib.get_symmetry_dataset(cb_cell_tuple)
    #

    if spglib_standardize == "off":
        # With standardization completely off, this should be practically identical to match_atoms_with_gemmi_spacegroup
        print("No standardization or change of basis performed. Matching structure as is.")
        match_atoms = atoms.copy()
        match_cell_tuple, match_gemmi_spgrp, match_sym_dataset = match_atoms_with_gemmi_spacegroup(match_atoms, clean=clean, wrap=wrap)

        return match_atoms, match_gemmi_spgrp, (None, None), (init_sym_dataset, None, match_sym_dataset)


    # If spglib_standardize is on, we standardize first. Otherwise, we prepare initial structure as is.
    # In both these case, we prepare match_sym_dataset, and either match_cell_tuple or match_atoms.
    # We should also keep track and pass forward the basis change, if applicable.
    if spglib_standardize == "on":
        bool_change_basis = (change_basis == "on")
        print("Standardizing structure through spglib before matching spacegroup.")
        if bool_change_basis:
            print("Applying change of basis during initial standardization.")
        else:
            print("No change of basis applied in initial standardization.")

        match_cell_tuple, match_sym_dataset, basis_change_transform, sym_datasets_tuple = spglib_standardize_structure(atoms,
                                                                                                                       to_primitive=to_primitive,
                                                                                                                       no_idealize=no_idealize,
                                                                                                                       change_basis=bool_change_basis
                                                                                                                       )
        
    elif spglib_standardize == "allow":
        ####################################
        ### Set init tuples and datasets ###
        ####################################
        # Init tuple and dataset, from input
        init_cell_tuple = atoms_to_spglib_cell_tuple(atoms)
        init_sym_dataset = spglib.get_symmetry_dataset(init_cell_tuple)

        match_cell_tuple = init_cell_tuple
        match_sym_dataset = init_sym_dataset
    else:
        raise RuntimeError("ERROR: Logic Error. This should not have happened.")


    # TODO: CONTINUE WORKFLOW HERE: NOW WE ATTEMPT MATCHING HERE...

    #####################################
    # Standardized tuple and dataset, according to settings
    #std_cell_tuple = spglib.standardize_cell(cell_tuple, to_primitive=to_primitive, no_idealize=no_idealize)
    #std_sym_dataset =  spglib.get_symmetry_dataset(std_cell_tuple)
    #orig_std_sym_dataset = copy.deepcopy(std_sym_dataset)

    #P_mat = std_sym_dataset["transformation_matrix"]
    #p_shift = std_sym_dataset["origin_shift"]
    ####################################
    

    # Here, try matching with standardized structure if allowed and needed
    # We have the name match_X for object type X of the structure we want to match
    #match_atoms = ??
    #match_cell_tuple = should be set
    #match_sym_dataset = should be set

    try:
        matching_gemmi_spgrps = find_matching_gemmi_spacegroup_from_spglib_dataset(match_sym_dataset, match_type="full", clean=clean, wrap=wrap)
    except Exception as exc:
        raise Exception("Spacegroup matching failed!") from exc


    # Now we check the results - and depending on settings, try again with appropriate action.
    if spglib_standardize == "allow" and matching_gemmi_spgrps == []:
        print("No matching spacegroups found, and spglib_standardize = allow. Standardizing structure according to spglib, and attempting to match again...")
        match_cell_tuple, match_sym_dataset, basis_change_transform, sym_datasets_tuple = spglib_standardize_structure(spglib_cell_tuple_to_atoms(match_cell_tuple),
                                                                                                                       to_primitive=to_primitive,
                                                                                                                       no_idealize=no_idealize,
                                                                                                                       change_basis=False
                                                                                                                       )
        try:
            matching_gemmi_spgrps = find_matching_gemmi_spacegroup_from_spglib_dataset(match_sym_dataset, match_type="full", clean=clean, wrap=wrap)
        except Exception as exc:
            raise Exception("Spacegroup matching failed!") from exc
        

    # Here: If no matches are found, then we attempt change of basis if allowed
    if change_basis == "allow" and matching_gemmi_spgrps == []:
        print("No matching spacegroups found, and change_basis = allow. Performing change of basis, and attempting to match again...")
        # Attempt basis change
        cb_cell_tuple = change_basis_cell_tuple(std_cell_tuple, P_mat, p_shift)
        cb_sym_dataset = spglib.get_symmetry_dataset(cb_cell_tuple)
        # Here, try matching again.
        try:
            matching_gemmi_spgrps = find_matching_gemmi_spacegroup_from_spglib_dataset(cb_sym_dataset, match_type="full", clean=clean, wrap=wrap)
        except Exception as exc:
            raise Exception("Spacegroup matching failed (after attempted change of basis)!") from exc
        
        # Then overwrite the std_cell and dataset
        std_cell_tuple = cb_cell_tuple
        std_sym_dataset = cb_sym_dataset



    # Now we perform checks: Check that there is one and only one unique match - or if several, it is number 68, since it has multiple entries.
    if not len(matching_gemmi_spgrps) <= 1:
        # We only allow this ambiguity for spacegroup number 68, since it has duplicate entries 
        # in the gemmi table. Hence it is not ambiguous. Moreover, we sanity check for operation equality.
        for grp in matching_gemmi_spgrps:
            assert grp.number == 68, "ERROR: Ambiguous spacegroup matching for spacegroup other than 68. Found {} matching spacegroups. This should not have happened!".format(len(matching_gemmi_spgrps))
            assert matching_gemmi_spgrps[0].operations() == grp.operations(), "ERROR: Spacegroup matches (>1) have differing operations! This should not have happened."

        #raise Warning("ERROR: Ambiguous spacegroup matching. Found {} matching spacegroups. This should not have happened!".format(len(matching_gemmi_spgrps)))
    else:
        assert len(matching_gemmi_spgrps) == 1, "ERROR: Spacegroup matching failed! Found {} matching spacegroups.".format(len(matching_gemmi_spgrps))
    
    gemmi_spgrp = matching_gemmi_spgrps[0]
    
    
    # With these checks above, we are done, and can convert back to atoms object.
    std_atoms = spglib_cell_tuple_to_atoms(std_cell_tuple)

    return std_atoms, gemmi_spgrp, (P_mat, p_shift), (init_sym_dataset, orig_std_sym_dataset, std_sym_dataset)

