import numpy as np
import pymatgen as pmg
from pymatgen.core import Structure, Lattice
import pymatgen.symmetry as pmgsym
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

import ase
import ase.io

import sys

# A routine for testing symmetry operations on a primitive cell


def pymatgen_structure_to_ase_atoms(structure: pymatgen.core.structure):
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
    print(struct == struct_from_file)

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

    print("Is found primitive the same as original?: ",
          found_primitive == struct_from_file)

    found_primitive.to(filename="Li2ZnGe_pmg_found_primitive.cif")
    conventional_struct.to(filename="Li2ZnGe_conventional.cif")

    # Atoms in the primitive unit cell seen inside the conventional unit cell

    primitive_cell_within_conventional_cell = ase.Atoms("Li9ZnGe",
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
                 primitive_cell_within_conventional_cell)

    return


if __name__ == "__main__":
    main()
