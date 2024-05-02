import numpy as np
import ase
import ase.io
from ase import Atoms
from ase import Atom
import itertools as it

##############################################################
## Functions to prepare structure for gridsampling with UCS ##
##############################################################


# To circumvent the annoying np.cross NoReturn bug
def cross2(a:np.ndarray,b:np.ndarray)->np.ndarray:
    return np.cross(a, b)


def remove_nonframework_cations(atoms:Atoms, cation:Atom):
    new_atoms = atoms.copy()
    while cation.symbol in new_atoms.get_chemical_symbols():
        ind = new_atoms.get_chemical_symbols().index(cation.symbol)
        new_atoms.pop(i=ind)

    return new_atoms

def remove_nonframework_cations_fancy(atoms:Atoms, cation:Atom):
    new_atoms = atoms.copy()
    del new_atoms[[atom.index for atom in new_atoms if atom.symbol==cation.symbol]]

    return new_atoms


# Presently, I do not need this. The supercells are already included in the input.
def compute_req_supercells(atoms, cutoff, isotropic_sc=False):
    assert isinstance(atoms, (ase.Atoms, tuple, list, np.ndarray)), "atoms needs to be instance of ase.Atoms, tuple, list or ndarray."
    
    if isinstance(atoms, ase.Atoms):
        cell = atoms.get_cell()
    else:
        cell = np.array(atoms)
    
    
    a, b, c = cell
    
    ab = cross2(a, b)
    ac = cross2(a, c)
    bc = cross2(b, c)

    AB = np.vstack((a, b, ab)).T
    AC = np.vstack((a, c, ac)).T
    BC = np.vstack((b, c, bc)).T

    xi_a = np.linalg.solve(BC, a)
    xi_b = np.linalg.solve(AC, b)
    xi_c = np.linalg.solve(AB, c)
    
    
    orth_dist = (abs(xi_a[-1]) * np.linalg.norm(bc), 
                abs(xi_b[-1]) * np.linalg.norm(ac), 
                abs(xi_c[-1]) * np.linalg.norm(ab)
                )


    if isotropic_sc:                
        min_dist = np.min(orth_dist)
        print(orth_dist)
        print(min_dist)
        
        # we need that cutoff < (min_dist_in_supercell/2)
        # If we hve a supercell of n unitcell along every edge,
        # the min distance would increase n times.
        # We need to determine smallest n such that n*(min_dist_in_unitcell/2) > cutoff
        n = int((cutoff // (min_dist/2)) + 1)
       
        nx, ny, nz = n, n, n  
    else:
        nx = int((cutoff // (orth_dist[0]/2))+ 1)
        ny = int((cutoff // (orth_dist[1]/2))+ 1)
        nz = int((cutoff // (orth_dist[2]/2))+ 1)

    # Extra check/control
    #combs = np.array(list(it.product((0,1,-1), repeat=3)))
    #assert (combs[0] == np.array((0, 0, 0))).all(), "Something went wrong. First row of combs is not [0, 0, 0]."
    #combs = combs[1:]
    #vector_combs = combs @ cell
    #distances = np.linalg.norm(vector_combs, 2, axis=1)
    #min_distance = np.min(distances)
    #print("Min distance check: ", min_distance)
    return (nx, ny, nz)


def make_supercell(atoms:Atoms, size:tuple):
    supercell_atoms = atoms.repeat(size)
    return supercell_atoms

# Takes fractional coordinates in unitcell and expresses them as corrdinates w r t supercell
def unitcell_to_supercell_frac_coords(coords, size:tuple, unitcell_ind=None):
    # Given supercell of size (nx, ny, nz), we assume the unitcell of interest is always in the middle
    # So nx, ny, nz must ALWAYS be odd numbers.
    # Hence (ni + 1) / 2 is an integer (i=x,y,z). For the unitcell (j, k, l) we need to add j*a + k*b + l*c to get the coordinates (if the first unitell is (0,0,0))
    # If the first unitcell is (1,1,1), then we clearly need to add (j-1)*a + (k-1)*b + (l-1)*c instead. a, b, and c are here the cell vectors.
    # It is easier to consider nx,ny,nz to be the total number of unitcells in each direction, in the supercell
    # Then we are interested in the unitcell ((nx+1)/2, (ny+1)/2, (nz+1)/2) since it has (ni+1)/2 - 1 = (ni-1)/2 = ni - (ni+1)/2, i.e. the same number of unitcells on each side.
    assert isinstance(size, tuple), "Error: size argument must be tuple."
    
    assert isinstance(unitcell_ind, (type(None), tuple, int, str)), "Unitcell index must be None, tuple, int or str type."

    if len(size) == 1:
        size = size*3 
    
    if unitcell_ind in ["center", "central", "middle"]:
        nx,ny,nz = size
        if nx % 2 == 0:
            nx -= 1
        if ny % 2 == 0:
            ny -= 1
        if nz % 2 == 0:
            nz -= 1
        unitcell_ind = np.array(((nx-1)/2, (ny-1)/2, (nz-1)/2))

    elif unitcell_ind is None:
        unitcell_ind = np.array((0,0,0))

    elif isinstance(unitcell_ind, int):
        unitcell_ind = (unitcell_ind,) * 3
   
    # transform coords: add unitcell index to unitcell scaled coordinates
    supercell_coords = np.copy(np.array(coords)).reshape(-1,3)
    unitcell_shift = np.array(unitcell_ind).reshape(-1,3)
    supercell_coords += unitcell_shift

    # Now it is shifted to correct unitcell. Then rescale to supercell
    # scaled coordinates, by dividing with supercell size:
    supercell_coords = np.divide(supercell_coords, np.array(size).reshape(-1,3))

    return supercell_coords

# unitcell_cart_coords -> supercell_cart_coords
# via fractional coordinates
def unitcell_to_supercell_cart_coords(coords, unitcell:Atoms, supercell:Atoms, size:tuple, frac_input=False):
    coords = np.array(coords)
    assert len(coords.shape) == 2 and coords.shape[1] == 3, "coords must be array-like of shape (n, 3)"
    if frac_input:
        unitcell_frac_coords = coords
    else:
        #unitcell_frac_coords = coords @ np.linalg.inv(np.array(unitcell.get_cell()))
        unitcell_frac_coords = np.linalg.solve(np.array(unitcell.get_cell()).T, coords.T).T
    
    supercell_frac_coords = unitcell_to_supercell_frac_coords(unitcell_frac_coords, size)

    supercell_cart_coords = supercell_frac_coords @ np.array(supercell.get_cell())
    

    return supercell_cart_coords
