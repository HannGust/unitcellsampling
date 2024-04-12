#!/usr/bin/env python

# Name suggestions for this script:
# ucs_batch_collect.py
# ucs_grid_from_batch_run.py
# ucs_assemble_batch_grid.py

"""Contains functions to collect and construct a grid from the individual
batch results, output and informaton from a batch run grid sampling.  """


import numpy as np
import os
import re
import ase
import ase.io
import gemmi
from ase.spacegroup import Spacegroup, get_spacegroup
from ase.io.cube import write_cube


def read_array(file, dtype=np.float64):
    """Read array from either a txt or npy file.
    """
    if not os.path.exists(file):
        raise FileNotFoundError("ERROR: File not found.")
    
    if not os.path.isfile(file):
        raise FileExistsError("ERROR: File is not an ordinary file.")

    if file[-4:] == ".txt":
        array = np.loadtxt(file)

    elif file[-4:] == ".npy":
        array = np.load(file)

    else:
        raise KeyError("ERROR: Not a supported file type. Array file need to be txt or npy.")

    return array

def read_indices(file):
    """Read indices from a batch indices file generated in a ucs batch run.
    """
    indices = np.atleast_2d(read_array(file))
    try:
        indices.astype(int, casting='safe')
    except:
        indices = np.rint(indices)
        indices = indices.astype(int)
        print("WARNING: Could not cast indices to int-type with safe casting. Array values were rounded before casting.")

    assert indices.dtype == np.int, "ERROR: Read and processed index array has dtype != int!"
    return indices


def read_energies(file, prec=np.float64):
    """Read energies from a batch energies file produced in a ucs batch run.
    """
    
    energies = read_array(file).astype(prec, casting='safe')

    return energies


def fill_grid_indices(grid, energies, indices, copy=True):
    """Fill the grid at the given indices (-1, 3) with the corresponding values
    in energies.
    """
    assert isinstance(grid, np.ndarray), "ERROR: Grid must be np.ndarray!"
    assert isinstance(energies, np.ndarray), "ERROR: Energies must be np.ndarray!"
    assert isinstance(indices, np.ndarray), "ERROR: Indices must be np.ndarray!"

    assert indices.dtype in [np.int, int], "ERROR: Index array must have dtype int!"
    
    assert indices.shape[0] == energies.shape[0], "ERROR: Indices and energies must have the same leading dimension (lenght)!"

    if copy:
        filled_grid = np.copy(grid)
    else:
        filled_grid = grid

    #indices_list_tmp = np.hsplit(indices, indices_or_sections=3) # could set indices_or_sections=indices.shape[1]
    #indices_list = [indxs.flatten() for indxs in indices_list_tmp]
    indices_list = [idxs.flatten() for idxs in np.hsplit(indices, indices_or_sections=3)]
    

    filled_grid[indices_list[0], indices_list[1], indices_list[2]] = energies
    
    #filled_grid[indices_list] = energies

    return filled_grid

def partial_grid_from_indices(energies, indices, grid_shape:tuple, fill_value=0.0):
    """Constructs a partial 3D grid from (-1) and (-1, 3)-shaped arrays of energies and 
    indices, together with the total grid shape.
    """

    grid = np.full(grid_shape, fill_value=fill_value, dtype=energies.dtype)
    indices_list = np.hsplit(indices, indices_or_sections=3) # could set indices_or_sections=indices.shape[1]
    grid[indices_list] = energies

    return grid


# TODO: Test this function, will it return list of full paths or only 
def construct_batch_dir_list(main_batch_run_wd, batch_dir_label="BATCH"):
    """Constructs and returns a sorted list of individual batch subdirectories
    given the path of the main batch working directory.
    Returns only the list of batch subdirectory names, not the full paths.
    """
    assert os.path.exists(main_batch_run_wd), "Main batch directory could not be found!"
    assert os.path.isdir(main_batch_run_wd), "Main batch directory must be directory!"
    batch_dir_list = [os.path.realpath(os.path.join(main_batch_run_wd,dir)) for dir in os.listdir(main_batch_run_wd) if (os.path.isdir(os.path.join(main_batch_run_wd,dir)) 
                             and re.fullmatch(batch_dir_label+"[0-9]+" , str(dir)))] 
    
    return sorted(batch_dir_list)



def read_batch_log(logfile):
    """Reads the contents of a main batch log to a single string.
    """
    if not os.path.exists(logfile):
        raise FileNotFoundError("Batch log \""+str(logfile)+"\" does not exists.")
    
    with open(logfile) as lg:
        log_txt = lg.read()

    return log_txt


def get_grid_shape_from_batch_log(log_txt):
    """Reads the grid size from a batch log.
    """
    grid_shape_pattern="\(nx, ny, nz\) =\(([0-9]+), *([0-9]+), *([0-9]+)\)"

    grid_shape_match = re.search(grid_shape_pattern, log_txt)

    if grid_shape_match is None:
        raise ValueError("Grid shape could not be read from batch log")
    
    grid_shape = tuple(int(grid_shape_match.group(i+1)) for i in range(3))

    return grid_shape


def get_symmetry_from_batch_log(log_txt):
    """Reads the batch log and extracts the symmetry option,
    and if needed, the spacegroup.
    """
    symmetry_setting_match_1 = re.search("nosym *: *(\w+)", log_txt)

    if symmetry_setting_match_1 is None:
        raise ValueError("symmetry_setting could not be read from batch log")
    
    symmetry_setting_txt_1 = symmetry_setting_match_1.group(1)

    if symmetry_setting_txt_1 not in ["True", "False"]:
        raise ValueError("symmetry_setting="+str(symmetry_setting_txt_1)+" read from batch log neither True nor False!")

    if symmetry_setting_txt_1 == "True":
        symmetry_setting = True
    elif symmetry_setting_txt_1 == "False":
        symmetry_setting = False
    else:
        raise ValueError("ERROR: symmetry setting neither True nor False.")

    # Read spacegroup if needed
    if symmetry_setting:
        spgrp_sc_match = re.search("Spacegroup, supercell: +([0-9]+) +([^ ].+)", log_txt)
        spgrp_uc_match = re.search("Spacegroup, unitcell: +([0-9]+) +([^ ].+)", log_txt)


        # DEBUG:
        #if spgrp_sc_match is None:
        #    print("sc spacegroup match is None.")
        #    print(spgrp_sc_match)

        #if spgrp_uc_match is None:
        #    print("uc spacegroup match is None.")
        #    print(spgrp_uc_match)

        # DEBUG:
        #print("spgrp mathed groups:")
        #print(spgrp_sc_match.groups())
        #print(spgrp_uc_match.groups())

        spgrp_sc_num = int(spgrp_sc_match.group(1))
        spgrp_uc_num = int(spgrp_uc_match.group(1))

        spgrp_sc_name = spgrp_sc_match.group(2)
        spgrp_uc_name = spgrp_uc_match.group(2)

        assert spgrp_sc_num == spgrp_uc_num, "Spacegroup numbers not consistent!"
        assert spgrp_sc_name == spgrp_uc_name, "Spacegroup names not consistent!"

        spgrp_from_num = Spacegroup(spgrp_uc_num)
        spgrp_from_name = Spacegroup(spgrp_uc_name)

        assert spgrp_from_num == spgrp_from_name, "Spacegroups are not the same!"

        spgrp = spgrp_from_num

    else:
        spgrp = None
        
    return symmetry_setting, spgrp

## Pseudocode to construct the grid:
# read batch_run log file, get relevant options (symmetry, grid size)
# Initialize an empty grid
# intialize full indices list
# initialize full energy list
# for every batch dir in batch dir list:
#    read indices
#    read energies
#    fill grid with energies on indices positions
#    store indices in full indices list

# Note that we need to think about normalization...
# the grid should be normalized in the way that the
# most negative energy is subtracted from all grid points
# However, we need to keep track so that when symmetry
# is applied, the sampled points are used and not any other
# non-sampled point, the latter which would result in
# sampled points being overwritten.

# we can find the shift by finding the minimum over all 
# sampled points
# If we initialize the grid to -np.inf, then 
# iteratively fill the grid, we can normalize
# the entire grid at once in the end. The grid points will
# be mapped to >=0 while -np.inf will remain negative. 
# Finally a symmetrize max operation can yield the correct application
# of symmetry.

# run unique on full indices list to check that it indeed is unique

# if symmetry:
#    read structure -> done
#    determine symmetry from structure -> done
#    read symmetry from batch run log -> done
#    assert that symmetry is the same from structrue and from logfile -> done
#    create gemmi grid of same size as grid 

# Perhaps a class could be useful
class UCSBatchAnalyzer():
    def __init__(self, batch_wd="."):
        self.work_dir = None
        self._set_batch_workdir(batch_wd) # This sets self.work_dir

        self.batch_dir_label = "BATCH"
        self.batch_dirs = self.get_list_of_batch_dirs()
        self.batch_size = None

        self.structure_file = os.path.join(self.work_dir, "structure.cif")
        assert os.path.exists(self.structure_file) and os.path.isfile(self.structure_file), "ERROR: structure.cif not found or not a file!"
        self.atoms = ase.io.read(self.structure_file)

        self.batch_logfile = os.path.join(self.work_dir, "ucs_batch.log")
        assert os.path.exists(self.batch_logfile) and os.path.isfile(self.batch_logfile), "ERROR; ucs_batch.log not found or not a file!"
        
        self.batch_logtxt = read_batch_log(self.batch_logfile)

        # Array lists:
        self.energy_arrays = None
        self.index_arrays = None

        # Results:        
        self.grid = None
        self.cube_file = os.path.join(self.work_dir,"energy_grid.cube")

        
    def _set_batch_workdir(self, batch_wd):
        """Checks logics and asserts existence of necessary files
        and folders in the batch run."""
        
        self._check_batch_workdir(batch_wd)

        self.work_dir = os.path.realpath(batch_wd)


    def _check_batch_workdir(self, batch_wd=None):
        """Checks existence of batch working directory."""
        if batch_wd is None:
            if self.work_dir is not None:
                batch_wd = self.work_dir
            else:
                raise ValueError("ERROR: Batch work dir and batch_wd are None.")
        
        if not (os.path.exists(batch_wd) and os.path.isdir(batch_wd)):
            raise FileNotFoundError("ERROR: Batch work dir:"+str(batch_wd))


    def read_main_log(self):
        """Reads in the main batch run log.
        """
        self.batch_logtxt = read_batch_log(self.batch_logfile)


    def get_list_of_batch_dirs(self):
        """Compiles a list of existing batch subdirectories
        belonging to the batch run.
        """
        self.batch_dirs = construct_batch_dir_list(self.work_dir, self.batch_dir_label)
        
        return self.batch_dirs


    def grid_shape_from_log(self):
        """Reads the grid shape from the log file.
        """
        grid_shape = get_grid_shape_from_batch_log(self.batch_logtxt)
        return grid_shape
    

    def symmetry_from_log(self):
        """Reads symmetry information, such as if symmetry is/was used,
        and in that case which spacegroup was employed.
        """
        use_symmetry, sp_group = get_symmetry_from_batch_log(self.batch_logtxt)
        return (use_symmetry, sp_group)


    def _collect_energy_arrays(self, fmt=".npy"):
        """Collects the energies into one array.
        """
        assert fmt in [".npy", ".txt"], "ERROR: Format not supported."
        
        energy_file = "energies" + fmt

        self.energy_arrays = tuple(
            read_energies(os.path.join(subdir, energy_file)) for subdir in self.batch_dirs
        )

        return self.energy_arrays


    def collect_energy_arrays(self):
        """Collects the energies into one array.
        fmt controls the prioritized format of the array-files.
        """
        fmts = [".npy", ".txt"]
        for fmt in fmts:
            try:
                energy_arrays = self._collect_energy_arrays(fmt=fmt)
            except FileNotFoundError:
                continue
            except Exception as exc:
                raise exc from Exception
            else:
                return self.energy_arrays
        raise Exception()
            

    def _collect_index_arrays(self, fmt=".txt"):
        """Collects the indices into one tuple of arrays."""
        assert fmt in [".npy", ".txt"], "ERROR: Format not supported."
        
        index_file = "indices" + fmt

        self.index_arrays = tuple(
                    read_indices(os.path.join(subdir, index_file)) for subdir in self.batch_dirs
                )
        
        return self.index_arrays


    def collect_index_arrays(self, fmt=".txt"):
        """Collects the indices into one tuple of arrays.
        Wrapper that first checks prioritized format, then the other
        if exception occurs.
        """
        allowed_fmts = [".npy", ".txt"]
        assert fmt in allowed_fmts, "ERROR: Format not supported."

        fmts = allowed_fmts.copy()
        fmts.remove(fmt)
        fmts.insert(0,fmt)

        for ext in fmts:
            try:
                index_arrays = self._collect_index_arrays(fmt=ext)
            except Exception as ex:
                exc = ex
                continue
            else:
                return self.index_arrays
        
        raise exc
    

    def sanity_check_indices(self):
        """Sanity check of indices: They should be unique.
        They should not overlap.
        THEY SHOULD WORK.
        """
        if not self.index_arrays:
            raise ValueError("ERROR: Cannot sanity check no-existent index arrays. Set index arrays first.")
 
        full_stacked_indices = np.vstack(self.index_arrays)
        stacked_unique_indices, inv_indx = np.unique(full_stacked_indices, axis=0, return_inverse=True)
        
        # Check uniqueness: Both check shape, and full arrays
        assert stacked_unique_indices.shape == full_stacked_indices.shape, "ERROR: Full stacked index array is not unique, shape mismatch!!!"
        assert (full_stacked_indices == stacked_unique_indices[inv_indx]).all(), "ERROR: Full index array is not unique!!!"
        

    def compile_grid(self):
        """Collects the results of the energy calculations in each individual
        and combines/compiles them into a full grid in the form of a
        3D numpy array.
        """
        # Get necessary info from logfile
        grid_shape = self.grid_shape_from_log()
        use_symmetry, spgrp = self.symmetry_from_log()

        # Determine spacegroup from structure, and
        # assert spacegroups detemrined this way and read from log
        # are identical, if using symmetry.
        # Added possibility to relax the precision when determining
        # the spacegroup a little if necessary
        symprec = 10**-5
        max_symprec = 10**-3

        if use_symmetry:
            while True:
                spgrp_from_atoms = get_spacegroup(self.atoms, symprec=symprec)
                try:
                    assert spgrp == spgrp_from_atoms, "ERROR: Spacegroups based on log vs. structure are not the same!!!"
                except AssertionError as ass_err:
                    if symprec < max_symprec:
                        print("WARNING: Spacegroup mismatch occurred at symmetry precision: ", symprec)
                        symprec = 10.0*symprec
                        continue
                    else:
                        raise AssertionError("Spacegroup mismatch at maximum symprec tolerance ", symprec) from ass_err
                else:
                    print("Spacegroup read from log vs. determined from structure equal at symprec: ", symprec)
                    break

        
        # Check indices
        idx_arrs = self.collect_index_arrays()
        self.sanity_check_indices()

        # Get list of batch directories
        batch_dirs = self.get_list_of_batch_dirs()
        
        # Initialize grid
        grid = np.full(grid_shape, -np.inf)
        
        # Loop over batch_dirs, collect things from each in turn
        for bch in batch_dirs:
            # Read batch energies:
            energies_file = os.path.join(bch, "energies.npy")
            if os.path.exists(energies_file) and os.path.isfile(energies_file):
                batch_energies = read_energies(energies_file)
            else:
                energies_file = os.path.join(bch, "energies.txt")
                batch_energies = read_energies(energies_file)

            # Read corresponding batch indices
            indices_file = os.path.join(bch, "indices.npy")
            if os.path.exists(indices_file) and os.path.isfile(indices_file):
                batch_indices = read_indices(indices_file)
            else:
                indices_file = os.path.join(bch, "indices.txt")
                batch_indices = read_indices(indices_file)
            
            # Fill grid with the energies at the indices:
            grid = fill_grid_indices(grid, batch_energies, batch_indices)

        # Then, we have filled the grid with all batches.
        # Now it remains to normalize and symmetrize.
        # We should first normalize - then it happens in np.float64 precision
            
        energy_arrs = self.collect_energy_arrays()
        full_energies_array = np.concatenate(energy_arrs, axis=0)
        min_energy = np.min(full_energies_array)
        
        # Normalize
        grid = grid - min_energy

        # symmetrize if necessary:
        if use_symmetry:
            #bool_grid = gemmi.Int8Grid(*grid_shape)
            float32_grid = grid.astype(np.float32)
            np.nan_to_num(float32_grid, copy=False)

            gemmi_float_grid = gemmi.FloatGrid(float32_grid)
            gemmi_float_grid.spacegroup = gemmi.find_spacegroup_by_number(spgrp_from_atoms.no)

            gemmi_float_grid.symmetrize_max()
            grid = np.array(gemmi_float_grid.array, dtype=np.float32)

            fill_value = np.nan_to_num(np.array(np.inf, dtype=np.float32))
            if not (grid >= 0.0).all():
                assert np.max(grid[grid < 0.0]) == -fill_value, "ERROR: Grid has negative values that are not -np.inf!!?"
                grid[np.isclose(grid, -fill_value)] = fill_value
        else:
            fill_value = np.nan_to_num(np.inf)

        # Fill the remaining unsampled points:
        np.nan_to_num(grid, copy=False, nan=fill_value, neginf=fill_value)
    
        self.grid = grid
        return self.grid
    

    def write_grid2cube(self, path=None):
        """Writes the compiled energy grid to a cube file."""
        if self.grid is None:
            raise ValueError("ERROR: There is no assembled grid to write. Compile a grid first.")
        with open(self.cube_file, 'w') as cube_file:
            write_cube(cube_file, self.atoms, data=self.grid)




