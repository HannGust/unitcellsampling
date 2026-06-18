#!/usr/bin/env python

"""Simple test that the imports can be done."""

import os

def test_main_package_import():
    import unitcellsampling
    return

def test_import_all():
    from unitcellsampling import batch_analyzer
    from unitcellsampling import cp2k_calculators
    from unitcellsampling import decorators
    from unitcellsampling import energy_calculator
    from unitcellsampling import lammps_calc_from_inp
    from unitcellsampling import preparatory_fcns
    from unitcellsampling import read_atom_types
    from unitcellsampling import sample
    from unitcellsampling import special_methods
    from unitcellsampling import symmetry
    from unitcellsampling import volume_exclusion
    
    return


def main():
    print("Testing main import:")
    test_main_package_import()
    print()
    print("Testing import of modules:")
    test_import_all()
    print("DONE.")


if __name__ == "__main__":
    main()
