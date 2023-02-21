import numpy as np
import pymatgen as pmg
from pymatgen.core import Structure, Lattice
import pymatgen.symmetry as pmgsym
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.analysis.structure_matcher import StructureMatcher

import ase
import ase.io

import sys

# A routine for testing symmetry operations on a primitive cell


def pymatgen_structure_to_ase_atoms(structure: pmg.core.structure):
    pass


def test_sym_op(atoms: ase.Atoms):
    pass


# Testing the pymatgen symmetry functions
# and the structure Li2ZnGe (structure no. 54189)

def main():
    if len(sys.argv) <= 1:
        raise Exception("Need an input.")

    atoms = ase.io.read(sys.argv[1])

    assert isinstance(
        atoms, ase.Atoms), "Input need to contain ase.Atoms obj readable with ase."

    f_atomic_coords = atoms.get_scaled_positions()
    atomic_symb = atoms.get_chemical_symbols()
    cell_params = atoms.get_cell_lengths_and_angles()

    lattice = Lattice.from_parameters(*cell_params)
    struct = Structure(lattice, atomic_symb, f_atomic_coords)

    struct_from_file = Structure.from_file(sys.argv[1])

    # print(struct)
    # print()
    # print(struct_from_file)

    # -> they yield the same file so just read from cif is best
    print("Is structure contructed Structure same as read from file: ", struct == struct_from_file)

    sg_analyzer = SpacegroupAnalyzer(struct_from_file)

    symops = sg_analyzer.get_symmetry_operations()
    conventional_struct = sg_analyzer.get_conventional_standard_structure()
    found_primitive = sg_analyzer.find_primitive()
    spgrp_num = sg_analyzer.get_space_group_number()
    spgrp_symb = sg_analyzer.get_space_group_symbol()
    lattice_type = sg_analyzer.get_lattice_type()
    crys_sys = sg_analyzer.get_crystal_system()

    print()
    print()
    print(symops)
    print(len(symops))
    print(conventional_struct)
    print(found_primitive)
    print(spgrp_num)
    print(spgrp_symb)
    print(lattice_type)
    print(crys_sys)

    # Naive comparisons:
    print("Naive comparisons of the structures:")
    print("------------------------------------")
    print("Is found primitive the same as original?: ",
          found_primitive == struct_from_file)
    print("Is found conventional the same as original?: ",
          conventional_struct == struct_from_file)
    print("Is primitive same as conventional?: ",
           found_primitive == conventional_struct)
    print("------------------------------------")

    print()
    print()
    print("Pymatgen StructureMatcher comparisons: ")
    ltol=0.01
    stol=0.01
    angletol=1
    matcher=StructureMatcher(ltol=ltol, stol=stol, angle_tol=angletol, primitive_cell=False)
    print("Parameters: frac length tol, itol="+str(ltol)+" ; site tolerance, stol="+str(stol)+" ; angle_tol="+str(angletol))
    print("--------------------------------------")
    print("Is found primitive the same as original?: ",
          matcher.fit(found_primitive, 
                               struct_from_file, 
                               symmetric=True))
    print("Is found conventional the same as original?: ",
          matcher.fit(conventional_struct, 
                      struct_from_file,
                      symmetric=True))
    print("Is primitive same as conventional?: ",
           matcher.fit(found_primitive, 
                       conventional_struct,
                       symmetric=True))
    print("--------------------------------------")

    #primitive_filename = "Li2ZnGe_pmg_found_primitive.cif"
    #conventional_filename = "Li2ZnGe_conventional.cif"
    primitive_filename = "pmg_found_primitive.cif"
    conventional_filename = "pmg_conventional.cif"
    found_primitive.to(filename=primitive_filename)
    conventional_struct.to(filename=conventional_filename)

    # Atoms in the primitive unit cell seen inside the conventional unit cell

    """primitive_cell_within_conventional_cell = ase.Atoms("Li9ZnGe",
                                                        scaled_positions=[[0.0, 0.0, 0.0],
                                                                          [0.0, 0.5,
                                                                              0.5],
                                                                          [0.5, 0.0,
                                                                              0.5],
                                                                          [0.5, 0.5,
                                                                              0.0],
                                                                          [1.0, 0.5,
                                                                              0.5],
                                                                          [0.5, 1.0,
                                                                              0.5],
                                                                          [0.5, 0.5,
                                                                              1.0],
                                                                          [1.0, 1.0,
                                                                              1.0],
                                                                          [0.25, 0.25,
                                                                              0.25],
                                                                          [0.75, 0.75,
                                                                              0.75],
                                                                          [0.5, 0.5, 0.5]],
                                                        cell=np.diag(
                                                            [6.164957, 6.164957, 6.164957]),
                                                        pbc=False
                                                        )

    ase.io.write("Li2ZnGe_primitive_within_convetional.cube",
                 primitive_cell_within_conventional_cell)"""

    return


if __name__ == "__main__":
    main()
