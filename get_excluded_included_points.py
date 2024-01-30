#!/usr/bin/env python

import numpy as np
import math as ma

import ase
import argparse
import unitcellsampling.sample as sample
import unitcellsampling.symmetry as symmetry
from unitcellsampling.preparatory_fcns import remove_nonframework_cations_fancy
from pathlib import Path

import gemmi
import itertools as it
from ase.spacegroup import get_spacegroup
import pymatgen
from pymatgen.core import Structure, Lattice
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer


def find_symmetrically_equivalent_indices(grid_shape, atoms, spacegroup_no=None):
    """ Finds the indices of symmetrically equivalent points, i.e. the
        equivalence classes of points due to symmetry, in erms of their
        indices. Returns a list of indices for each equivalence class."""

    tot_bol_grid = gemmi.Int8Grid()

    if not spacegroup_no:
        spgroup = gemmi.find_spacegroup_by_number(get_spacegroup(atoms))
    else:
        spgroup = gemmi.find_spacegroup_by_number(spacegroup_no)

    tot_bol_grid.spacegroup = spgroup

    tot_bol_grid.set_size(*grid_shape)


    eq_class_list = []

    for indx in it.product(range(grid_shape[0]), range(grid_shape[1]), range(grid_shape[2])):
        if tot_bol_grid.get_value(*indx) != 0:
            continue
        else:
            # Inititalize a new temporary integer grid with all 0 entries            
            tmp_bol_grid = gemmi.Int8Grid()
            tmp_bol_grid.spacegroup = spgroup
            tmp_bol_grid.set_size(*grid_shape)

            # Set the value of the current point to 1, then symmetrize
            tmp_bol_grid.set_value(*indx, 1)
            tmp_bol_grid.symmetrize_max()

            # Do the same for the total grid (to keep track of which points have been considered)
            tot_bol_grid.set_value(*indx, 1)
            tot_bol_grid.symmetrize_max()

            # Find the list (array) of indices for nonzero entries in the temporary grid, append it to the list of equivalence classes
            eq_class_list.append(np.argwhere(tmp_bol_grid.array))

    return eq_class_list




parser = argparse.ArgumentParser(description="Computes number of grid points in a structure, both total, to be calculated, and excluded, given specified grid settings.")

### Definition and parsing of arguments
parser = argparse.ArgumentParser(description='Energy sampling of a (periodic) solid system with an added atom/ion on the grid.')
### Structure parser: (identical to the arguments in grid_gen_script.py))
parser.add_argument('file', metavar='filename', type=str, action='store', help="Name of the structure file containing the unit cell structure to be sampled. Formats: cif, cube, xsf, ...")
parser.add_argument('-a', '--atom', type=str, action='store', default='Li', help="Specify the atom/ion used for sampling.")
parser.add_argument('-g', '--grid', type=int, action='store', default=[10], nargs='+', help="Specify the number of grid points in each dimension (or cubic grid) (mutually exclusive with \"--space\").")
parser.add_argument('-s', '--space', type=float, action='store', default=None, nargs='+', help="Specify the spacing between the grid points in each dimension (mutually exclusive with \"--grid\").")
parser.add_argument('--vdw', type=float, action='store', default=0.0, help="Specify the fraction of the van der Waals radius that should be excluded from the sampled volume around each atom in the host framework structure.")


parser.add_argument('--mic-cutoff', type=float, action='store', default=0.0, help="Specify the cut-off used to construct a supercell obeying the minimum image convention.")

parser.add_argument('--nosym', action='store_false', help="Turns off usage of spacegroup symmetry. Default is to apply symmetry to save the number of required calculations.")
parser.add_argument('--ra', action='store_true', help="Specify whether to remove all atoms of the type that is used for sampling from the structure, before doing the sampling.")
parser.add_argument('--sg', type=int, default=None, action='store', help="Manually specify the spacegroup to use for symmetry. Default is None, in which case spacegroup will be automatically determined from the structure.")
parser.add_argument('--conv', '--conventional-cell', action='store_true', help="If given, will check whether unit cell is conventional cell, and if not, determine the conventional cell based on the input structure and then sample the found conventional cell instead.")
parser.add_argument('--midvox', action='store_true', help="Specifies that sampling should be done in the center of voxels rather than in the corners, as is default. Corresponds to shifting the coordinates of the grid points with 0.5 * (1/na, 1/nb, 1/nc), where ni is the number of grid points in direction of unit cell vector i.")

parser.add_argument("-p", "--printlevel", type=int, default=0, action='store', help="Sets the print level. Current options 0 or 1 (default 0). 1 prints exact input settings before results. 0 or anything other than 1, prints only results.")
###

# Parse args 
args = parser.parse_args()

### Preprocessing the structure input

input_file = Path(args.file)

use_sym = args.nosym

atom = args.atom

## Check existence of inputfile, extract the name of the inputfile and directories
if not input_file.exists():
    raise FileNotFoundError('The given inputfile \'' + str(input_file) + '\' was not found.')

# Set the input filenames and directory, and read input
indir = input_file.parent
infile =  input_file.name

# TODO: Change this reading to the full given path for infile??
input_unitcell = ase.io.read(Path(indir, infile))
init_unitcell = input_unitcell.copy()

# Change the unitcell to the conventional cell if symmetry 
# is used or if conv argument is specified
if use_sym or args.conv:
    if not symmetry.is_conventional_cell(init_unitcell):
        sg_analyzer = SpacegroupAnalyzer(AseAtomsAdaptor.get_structure(init_unitcell))
        conv_cell = sg_analyzer.get_conventional_standard_structure()
        init_unitcell = AseAtomsAdaptor.get_atoms(conv_cell)

## FOR DEBUG/TEST: Is initial spacegroup the same as
# the one obtained after removing ions???
init_spgrp = get_spacegroup(init_unitcell)


if args.ra:
    print("Removing sampling atoms", str(atom),"from structure.")
    atoms_temp = remove_nonframework_cations_fancy(init_unitcell, ase.Atom(atom))
    unitcell = atoms_temp
else:
    unitcell = init_unitcell

if use_sym:
    if args.sg:
        spacegroup = args.sg
    else:
        spacegroup = get_spacegroup(unitcell)
else:
    spacegroup = None

## Set the number of points in the grid in each dimension (or equivalently, the mesh size)
if not args.space:
    if len(args.grid) == 1:
        nx,ny,nz = args.grid * 3
    elif len(args.grid) == 2:
        raise Exception("--grid must be either 1 or 3 int inputs.")
    else:
        nx,ny,nz = args.grid[0:3]

    a,b,c = np.linalg.norm(unitcell.get_cell()[0,:]),\
            np.linalg.norm(unitcell.get_cell()[1,:]),\
            np.linalg.norm(unitcell.get_cell()[2,:])
    true_spacing = (a/nx, b/ny, c/nz)

    # There should be a symmetry adaption here probably ... 

    print('True spacing: ', true_spacing)
    print('(nx, ny, nz) =', (nx, ny, nz))

else:
    # Small section to compute nx, ny, nz i.e. the # of points in each direction in case the --space option is used
    if args.space:
        if len(args.space) == 1:
            spacing_x, spacing_y, spacing_z = args.space * 3
        elif len(args.space) == 2:
            raise Exception("--space must be either 1 or 3 float inputs.")
        else:
            spacing_x, spacing_y, spacing_z =  args.space[0:3]

        a,b,c = np.linalg.norm(unitcell.get_cell()[0,:]),\
                np.linalg.norm(unitcell.get_cell()[1,:]),\
                np.linalg.norm(unitcell.get_cell()[2,:])
        nx,ny,nz = ma.ceil(a/spacing_x), ma.ceil(b/spacing_y), ma.ceil(c/spacing_z)

        if use_sym:
            nx,ny,nz = symmetry.find_spacegroup_compatible_gridshape((nx,ny,nz), spacegroup, search_denser=True)
            print("Determined grid spacing compatible with spacegroup.")

        true_spacing = (a/nx, b/ny, c/nz)
        print('Desired spacing: ', (spacing_x, spacing_y, spacing_z),' True spacing: ', true_spacing)
        print('(nx, ny, nz) =', (nx, ny, nz))

print("vdW cutoff factor: ", args.vdw)


# Now generate grid for unitcell:
unitcell_ucs = sample.UnitCellSampler(unitcell)

if (args.midvox and use_sym):
    print("WARNING: Both symmetry and midvoxel enabled. Their compatibility has NOT been properly tested!!!")

unitcell_grid, unitcell_included = unitcell_ucs.generate_grid_vectors((nx, ny, nz), vdw_scale=args.vdw, midvox=args.midvox) # DONE

unitcell_grid = unitcell_grid.reshape(-1,3)
unitcell_included = unitcell_included.reshape(-1)

if args.printlevel == 1:
    print("Computed number of grid points, total, inlcuded and excluded, \
       for the following settings:")
    print("-"*79)
    print("Structure: ", args.file)
    print("Sampling atom: ", args.atom)
    print("Grid mesh size: ", (nx, ny, nz))
    print("Grid spacing: ", true_spacing)
    print("Removed sampling atoms: ", args.ra)
    print("Excluded vdw cutoff: ", args.vdw)
    print("Using symmetry: ", use_sym, 
      int(use_sym)*("; Spacegroup number: "+str(spacegroup.no)
         + "; Spacegroup name: " + str(spacegroup.symbol)))

    ## for spacegroup debug
    print("Init spacegroup number and name: ", str(init_spgrp.no), str(init_spgrp.symbol))
    print("Inital spacegroup same as after ion removal?: ", init_spgrp == spacegroup)
    ###
    print("Using conventional cell: ", args.conv)
    print("Midvox-grid: ", args.midvox)

    print("-"*79)
    print()

print("-"*79)
print("RESULTS FOR STRUCTURE: ", infile)
print("-"*79)

#print("unitcell_grid shape: ", unitcell_grid.shape)
#print("nx * ny * nz: ", np.product((nx,ny,nz)))

assert np.product(unitcell_grid.shape[0]) == np.product((nx,ny,nz)) == len(unitcell_grid), "Total number of grid points differnt with different methods?!?!"

print("Total number of points: ", unitcell_grid.shape[0])
print("Number of INCLUDED points (i.e. calculations): ", np.sum(unitcell_included))
print("Number of EXCLUDED points: ", np.sum(np.logical_not(unitcell_included)))
assert unitcell_grid.shape[0] == (np.sum(unitcell_included) + np.sum(np.ones_like(unitcell_included) - 1*unitcell_included)), "Number of calculations and points does not add up correctly!!!" # Sanity check
if use_sym:
    sym_ind = find_symmetrically_equivalent_indices((nx, ny, nz), unitcell)
    print(type(sym_ind), type(sym_ind[0]))
    vdw_excluded_ind = np.argwhere(np.logical_not(unitcell_included))
    vdw_included_ind = np.argwhere(unitcell_included)

    num_excluded_sym_classes = 0
    excl_sym_classes = []

    num_included_sym_classes = 0
    incl_sym_classes = []
    #incl_sym_classes = []
    for indx_list in sym_ind:
        #excl_switch = False
        for excl_ind in vdw_excluded_ind:
            if list(excl_ind) in indx_list:
                num_excluded_sym_classes += 1
                excl_sym_classes.append(indx_list)
                #excl_switch = True
                break
        #if not excl_switch: 
        #    incl_sym_classes.append(indx_list)

        for incl_ind in vdw_included_ind:
            if list(incl_ind) in indx_list:
                num_included_sym_classes += 1
                incl_sym_classes.append(indx_list)
                break
        
    num_sym_classes_in_both = 0
    for i in incl_sym_classes:
        for j in excl_sym_classes:
            if i.shape == j.shape:
                if (i == j).all():
                    num_sym_classes_in_both += 1
    print("Number of equivalence classes in both excluded and inlcuded", num_sym_classes_in_both )

    print("before scaled coord")
    scaled_grid_coord = np.linalg.solve(np.array(unitcell.cell[:]).T, unitcell_grid.T).T
    #print("ase_spgrp_1")
    #ase_spgrp_unique_sites = spacegroup.unique_sites(scaled_grid_coord, map_to_unitcell=True)
    #print("ase_spgrp_2")
    #ase_spgrp_unique_sites_incl = spacegroup.unique_sites(scaled_grid_coord[unitcell_included], map_to_unitcell=True)

    assert scaled_grid_coord.shape == unitcell_grid.shape, "Scaled and cartesian corodinate arrays do not have the same shape!!!"
    print("Order of spacegroup: ", spacegroup.nsymop)
    #print("Number of symmetry-unique points in full grid (with ase.spacegroup): ", len(ase_spgrp_unique_sites), ase_spgrp_unique_sites.shape)
    print("Number of symmetry-unique points in full grid (with gemmi): ", len(sym_ind))
    #print("Number of INCLUDED symmetry-unique points (ase.spacegroup): ", len(ase_spgrp_unique_sites_incl), type(ase_spgrp_unique_sites_incl), ase_spgrp_unique_sites_incl.shape)
    print("Number of INCLUDED symmetry-unique points (gemmi): ", len(sym_ind) - num_excluded_sym_classes)
    print("Number of excluded symmetry-classes (gemmi): ", num_excluded_sym_classes)
    print("Number of excluded points (gemmi): ", "Sym Classes: ", len(excl_sym_classes), "Total excluded points (unfolding symmetry): ", np.sum([len(indxs) for indxs in excl_sym_classes]))

    print("Number of inlcuded symmetry-classes (gemmi): ", num_included_sym_classes)
    print("Number of included points (gemmi): ", "Sym Classes: ", len(incl_sym_classes), "Total included points (unfolding symmetry): ", np.sum([len(indxs) for indxs in incl_sym_classes]))
    print(sym_ind[0],sym_ind[0][0], sym_ind[1][0], sym_ind[2][0], sym_ind[3][0])
    #print(np.divide([sym_ind[0][0], sym_ind[1][0], sym_ind[2][0], sym_ind[3][0]] ,[nx, ny, nz]))
    #print(ase_spgrp_unique_sites * np.array((nx,ny,nz)))
print("-"*79)
