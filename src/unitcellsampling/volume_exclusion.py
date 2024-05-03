""" Module containing the new volume exclusion classes and routines that are
    in use and/or tested for intended use in the unitcellsampling
    package.

    (NOTE: Originally, based on the VolExclusionClass.py module in the volume_exclusion
    package/repository.)"""

### OBSERVE: ONLY INITIALLY, MINIMALLY TESTED! ###

import ase
import ase.build
import ase.cell
import numpy as np

from mendeleev import element
from sklearn.neighbors import KDTree

import unitcellsampling.sample
from unitcellsampling.sample import UnitCellSampler as UCS
from unitcellsampling.preparatory_fcns import compute_req_supercells
from unitcellsampling.preparatory_fcns import unitcell_to_supercell_frac_coords

#### IDEAS 
# I could either create these classes and apply them as functions 
# That is, the contain the strategy to exclude, but not the 
# information about anything else, like the grid, or so.
#
# I could also make the m contain the informaiton about the grid, atoms and so on.
# I think maybe the first approach is better.
####

# Not finished - might become something at a later point
# Not planned currently.
#class GridPointExcluder:
#    def __init__(self, grid=None, shape=None):
#        self.grid = None
#        self.mask = None
#        self.shape = None

### Idea for the radial excluder:
# Take a radius - either constant float/int
# Take a map that takes atom types and returns floats

# This should is meant as an overarching class for radial exclusion
## 
class RadialExcluder:
    def __init__(self, radii):
        assert isinstance(radii, (float, int, dict)), "radii must be float, int or dict."
        
        # Type check and conversions
        if isinstance(radii, int):
            radii = float(radii)
        #self.radii = radii
        
        # If dict, cast int values to floats, and check that all are non-negative floats

        if isinstance(radii, dict):
            self.atomtype2radius_map = radii
            for key, val in radii.items():
                new_val = float(val)
                assert new_val >= 0.0, "Radii mapping must return non-negative radii."
                self.atomtype2radius_map[key] = new_val
            radii = list(self.atomtype2radius_map.values())
            self.radii = radii
        elif isinstance(radii, float):
            assert radii >= 0.0, "Radii must be non-negative."
            self.atomtype2radius_map = {"default":radii}
            self.radii = [radii]
        else:
            raise TypeError("Radii must be convertable to float or be dictionary mapping returning floats.")
        
        # Sanity check
        for key, val in self.atomtype2radius_map.items():
            assert isinstance(val, float) and val >= 0.0, "Radii mapping must return non-negative radii. This should not have happened!"
    
    
    def print_settings(self):
        """Prints attributes."""
        print("Radii: ", self.radii)
        print("Atomtype2radius_map: ", self.atomtype2radius_map)


    def get_radii_list(self, atoms:ase.Atoms):
        """Constructs and returns a list of radii."""

        cutoff_radii = [self.atomtype2radius_map[atom] if atom in self.atomtype2radius_map.keys() else self.atomtype2radius_map["default"] for atom in atoms.get_chemical_symbols()]
        return cutoff_radii
    

    def get_full_radii_mapping(self, atoms:ase.Atoms):
        """Simple short way of constructing full radii mapping for an atoms object."""

        full_radii_mapping = {atom:self.atomtype2radius_map[atom]
                            if (atom in self.atomtype2radius_map.keys()) 
                            else self.atomtype2radius_map["default"] 
                            for atom in 
                            list(dict.fromkeys(atoms.get_chemical_symbols()))}
        
        return full_radii_mapping


    def _get_full_radii_mapping1(self, atoms:ase.Atoms):
        """Constructs the full atom-type to radii mapping for an atoms-object 
        by unfolding the default value, i.e. adding missing atom-types by 
        filling in the default value."""

        # Check that neccesary radii and atom types are present
        required_atom_types = set(atoms.get_chemical_symbols())
        existing_atom_types = set(self.atomtype2radius_map.keys())
        existing_atom_types.discard('default')

        assert 'default' in self.atomtype2radius_map or existing_atom_types.issuperset(required_atom_types), "Radii mapping missing atom-types! Add the missing atoms or a default radius!"

        # Init byt copying existing mapping:
        full_radii_mapping = self.atomtype2radius_map.copy()

        # Remove default if present:
        if "default" in full_radii_mapping.keys():
            _ = full_radii_mapping.pop("default")
        

        types_to_be_added = list(required_atom_types.difference(existing_atom_types))
        additional_types_dict = dict.fromkeys(types_to_be_added, self.atomtype2radius_map["default"])
        full_radii_mapping.update(additional_types_dict)

        # For extra info:
        atom_types_mapped_to_default = list(additional_types_dict.keys())
        atom_types_with_existing_mapping = list(existing_atom_types)

        return full_radii_mapping, atom_types_with_existing_mapping, atom_types_mapped_to_default


    def _get_full_radii_mapping2(self, atoms:ase.Atoms):
        """Alternative method to get full radii mapping for an atoms object."""

        # Check that neccesary radii and atom types are present
        required_atom_types = set(atoms.get_chemical_symbols())
        existing_atom_types = set(self.atomtype2radius_map.keys())
        existing_atom_types.discard('default')

        assert 'default' in self.atomtype2radius_map or existing_atom_types.issuperset(required_atom_types), "Radii mapping missing atom-types! Add the missing atoms or a default radius!"

        # Initialize new dictionary mapping
        full_radii_mapping = dict()

        # For extra info:
        atom_types_with_existing_mapping = []
        atom_types_mapped_to_default = []
        # Alternative method:
        for atomtype in list(required_atom_types):
            if atomtype in existing_atom_types:
                full_radii_mapping[atomtype] = self.atomtype2radius_map[atomtype]
                atom_types_with_existing_mapping.append(atomtype)

            else:
                full_radii_mapping[atomtype] = self.atomtype2radius_map["default"]
                atom_types_mapped_to_default.append(atomtype)

        return full_radii_mapping, atom_types_with_existing_mapping, atom_types_mapped_to_default


    def construct_cutoff_filters(self, atoms, grid_coords, frac_input=False, periodic=True):
        """Constructs masks that yields included and excluded points according to
        the raidal cutoff."""
        assert isinstance(grid_coords, np.ndarray) and len(grid_coords.shape) == 2 and grid_coords.shape[1] == 3, "grid_coords must be shape (n, 3) numpy-array."
        # What do we need to do?
        # We need to first, if perioic is True, we need to make sure we construct
        # and appropriately sized supercell. The supercell size in terms of 
        # unitcells should be the smallest odd integer in each direction
        # such that the minimum image convention is fulfilled with 
        # respect to the largest cutoff radius. It must also be at least 3,
        # in each direction, since this is the smallest case where we have a central
        # unitcell surrounded by at least one unitcell on each side.
        #
        # Then we need to map the grid coordinates to this supercell, more 
        # specifically, make the unitcell grid lie within the central
        # unitcell in the constructed supercell, in terms of supercell coordinates.
        #
        # If we do not want a periodic application, we can simply skip the
        # above steps.
        #
        # With these coordinates, we need to apply the radial cutoff taking 
        # all atoms in the supercell into account.

        if periodic:
            max_radius = max(self.radii)
            sc_nx, sc_ny, sc_nz = compute_req_supercells(atoms, max_radius)

            # x
            if sc_nx % 2 == 0:
                sc_nx += 1
            sc_nx = max(3,sc_nx)
            # y
            if sc_ny % 2 == 0:
                sc_ny += 1
            sc_ny = max(3,sc_ny)
            # z
            if sc_nz % 2 == 0:
                sc_nz += 1
            sc_nz = max(3,sc_nz)

            sc_size = (sc_nx, sc_ny, sc_nz)
            P = np.diag(np.array(sc_size))
            supercell = ase.build.make_supercell(atoms, P, wrap=True)

            # Should this be fractional coordinates? Yes? Maybe? 
            # the uc to sc cart coord function does not include the unitcell index...
            # Thus, I want to have control
            
            # Let nx be number of uc in x direction. Then the indices are 0, 1, ..., nx-1.
            # nx is odd, hence nx-1 is even. (nx-1)/2 is an integer. This has (nx-1)/2 
            # unitcells before it, since 0, 1, ... , (nx-1)/2 - 1 are the indices of unitcells
            # before it. Likewise (nx-1)/2 + 1, (nx-1)/2 + 2, ... , nx-1 = (nx-1)/2 + ((nx-1)/2)
            # i.e. it has (nx-1)/2 unitcells after it. Thus unitcell (nx-1)/2 it in the middle,
            # if nx is the total number of unitcells in this direction and it is odd.
            # 
            # alt argument: k+1 In the middle if 2 k + 1 = total, since k, 1, k. k+1 has index k if we start at 0.
            # k = (total-1)/2.

            # We round to make sure floating point errors do not make it truncate errouneously
            center_unitcell_idx = (round((sc_nx-1)/2), 
                                   round((sc_ny-1)/2), 
                                   round((sc_nz-1)/2)
                                   )
            # Sanity check
            assert center_unitcell_idx == ((sc_nx-1)//2, (sc_ny-1)//2, (sc_nz-1)//2), "Error: Center unitcell indices from // and round did not match!!!"

            if not frac_input:
                # Need to convert to fractional coord in this case
                grid_coords = np.linalg.solve(np.array(atoms.get_cell()).T, np.array(grid_coords).T).T
            
            # Compute fractional coordinates of the grid within the supercell:
            coords = unitcell_to_supercell_frac_coords(grid_coords,
                                                size=(sc_nx, sc_ny, sc_nz),
                                                unitcell_ind=center_unitcell_idx)
            # Convert to 
            coords = coords @ np.array(supercell.get_cell())

        else:
            # We do not need to do anything here
            supercell = atoms
            if frac_input:
                # Convert to cartesian coordinates
                coords = grid_coords @ np.array(supercell.get_cell())

        # At this point it should be the case that:
        # supercell is the correct supercell no matter what setting
        # coords are the correct cartesian coordinates no matter what setting
        #print("DEBUG: supercell: ", supercell)
        #print("DEBUG: supercell.get_positions(): ", supercell.get_positions())
        #print("DEBUG: coords: ", coords)

        # Construct the tree:
        kdtree = KDTree(coords, metric="euclidean")

        # This line is supposed to create a list by mapping each atom to a radius using the atom type (symbol) to radius mapping
        # We use the "default" entry in the map when the atom type does not exist
        # If the default entry does not exist, and the atom symbol does not exist either,
        # then this would result in some kind of error hopefully, since the attempted key is not in the dictionary
        
        cutoff_radii = self.get_radii_list(atoms=supercell)
        #print("DEBUG: cutoff_radii: ", cutoff_radii)

        indices = kdtree.query_radius(X=supercell.get_positions(), r=cutoff_radii)
        #print("DEBUG: indices: ", indices)
        uniq_indices = np.unique(np.hstack(indices))
        # indices is an array of index arrays, one for each point in X
        # each index array lists the points in the initial tree X data, i.e.
        # the array of point coordinates given in the initialization of the tree.
        #print("DEBUG: uniq_indices: ", uniq_indices)
        # Now, use the uniq_indices to obtain a filer for grid coords
        
        within_cutoff_mask = np.full(grid_coords.shape[0], fill_value=False, dtype=np.bool8)
        
        outside_cutoff_mask_test = np.full(grid_coords.shape[0], fill_value=True, dtype=np.bool8)

        np.put(within_cutoff_mask, uniq_indices, True, mode="raise")
        
        np.put(outside_cutoff_mask_test, uniq_indices, False, mode="raise")
        outside_cutoff_mask = np.logical_not(within_cutoff_mask)

        assert (outside_cutoff_mask == outside_cutoff_mask_test).all(), "outside_cutoff_mask differs with different methods!"

        return within_cutoff_mask, outside_cutoff_mask


# This is utility wrapper class around the Radial excluder
# Purpose to apply scaled vdw-radial cutoff exlcusion
# Constructs a radial exlcuder object with appropriately scaled
# van der Waals-radii cutoffs.
class ScaledVdWExcluder:
    def __init__(self, vdw_scaling_map):
        if not isinstance(vdw_scaling_map, (int, float, dict)):
            raise TypeError("vdw_scaling_map must be int, float or dict type.")
        
        if isinstance(vdw_scaling_map, dict):
            self.vdw_scaling_map = {key:float(val) for key,val in vdw_scaling_map.items()}
        else:
            vdw_scaling_map = float(vdw_scaling_map)
            assert vdw_scaling_map >= 0.0, "vdw_scaling_must be >= 0.0 (non-negative float or int)"
            self.vdw_scaling_map = dict(default=vdw_scaling_map)
        
        # Sanity check:
        for key,val in self.vdw_scaling_map.items():
                assert isinstance(val, float) and val >= 0.0, "vdw_scaling_map must rern non-negative float."
    

    def generate_excluder(self, atoms:ase.Atoms):
        """Generates a vdW-excluder for an atoms object from
        the settings (vdW scaling factors). Returns a RadialExlcuderObject
        with the desired settings."""
        # Note, dividing with 100 to get to Ångström from picometer
        #radii = must be dict of atom-types to radii in Ångström
        # This is obtained from scaling each atom type with its
        # desirec radii
        #required_atom_types = list(dict.fromkeys(atoms.get_chemical_symbols()))
        #existing_atom_types = [atomtype if atomtype != "default" 
        #                       for atomtype in self.vdw_scaling_map.keys()]
        scaled_vdw_radii_map = {
            atom:(element(atom).vdw_radius * self.vdw_scaling_map[atom]/100.0)
            if atom in self.vdw_scaling_map.keys()
            else
            (element(atom).vdw_radius * self.vdw_scaling_map["default"]/100.0)
            for atom in list(dict.fromkeys(atoms.get_chemical_symbols()))
            }
        
        self.radial_excluder = RadialExcluder(radii=scaled_vdw_radii_map)
        return self.radial_excluder
    

    def construct_cutoff_filters(self, atoms, grid_coords, frac_input=False, periodic=True):
        """Constructs cutoff filters for a grid from an atoms object and the grid coordinates."""

        radial_excluder = self.generate_excluder(atoms)

        assert radial_excluder is self.radial_excluder, "Radial excluder not same as self.radial_excluder!"

        filters = radial_excluder.construct_cutoff_filters(atoms,
                                                        grid_coords,
                                                        frac_input=frac_input,
                                                        periodic=periodic
                                                        )
        
        return filters
    

    def print_settings(self):
        print("vdw_scaling_map: ", self.vdw_scaling_map)
        print()
        if self.radial_excluder is not None:
            print("Radial excluder info\n---------------------")
            self.radial_excluder.print_settings()
        else:
            print("radial_excluder: None - Not generated yet.")
   
