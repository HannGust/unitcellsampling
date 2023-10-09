import ase

import os
from pathlib import Path

## Original subdir_calc
# TODO: Add argument for getting information from 'master' calculation
#def subdir_calc(calculate_energy):
#    """
#    Decorator for performing calculation in folder specified
#    by environment variable UCS_CALCULATION_DIR.
#    """
#    def calculate_in_subdir(atoms: ase.Atoms):
#        cwd = os.getcwd()
#        calc_dir = os.getenv("UCS_CALCULATION_DIR")
#        sub_dir_name = '_'.join(
#            ['{:.3f}'.format(c)
#             for c in atoms.get_scaled_positions()[-1]])
#        sub_dir = Path(calc_dir, sub_dir_name)
#
#        Path(calc_dir).mkdir(exist_ok=True, parents=True)
#        sub_dir.mkdir(exist_ok=False)
#        os.chdir(sub_dir)
#
#        ase.io.write('atoms.cif', atoms)
#        potential_energy = calculate_energy(atoms)
#
#        os.chdir(str(cwd))
#        return potential_energy
#
#    return calculate_in_subdir
##

## Modified original subdir_calc
# DONE: Added general **kwargs argument passing from decoartor to calculator
def subdir_calc(calculate_energy):
    """
    Decorator for performing calculation in folder specified
    by environment variable UCS_CALCULATION_DIR.
    """
    def calculate_in_subdir(atoms: ase.Atoms, **kwargs):
        cwd = os.getcwd()
        calc_dir = os.getenv("UCS_CALCULATION_DIR")
        sub_dir_name = '_'.join(
            ['{:.3f}'.format(c)
             for c in atoms.get_scaled_positions()[-1]])
        sub_dir = Path(calc_dir, sub_dir_name)

        Path(calc_dir).mkdir(exist_ok=True, parents=True)
        sub_dir.mkdir(exist_ok=False)
        os.chdir(sub_dir)

        ase.io.write('atoms.cif', atoms)
        potential_energy = calculate_energy(atoms, **kwargs)

        os.chdir(str(cwd))
        return potential_energy

    return calculate_in_subdir
#


# New subdir_calc : in progess
# TODO: Finish code, decide labels, add **kwargs passing if this works in prev subdir calc, and test
# TODO: Think about passing restart-files if theres any. Also think about this and parallelization
def subdir_calc_new(calculate_energy):
    """
    Decorator for controlling output and
    performing calculations in folder specified
    by environment variable UCS_CALCULATION_DIR.

    TODO: (Implement) Reads a second environment variable, UCS_OUTPUT_LEVEL,
    to determine how much output to produce. LOW is default, produces only
    output for the first and last calculations.
    FULL produces output for all calculations.
    NONE disables output altogether.
    Finally, setting it to an integer yields printing every nth calculation.

    """

    calc_dir = os.getenv("UCS_CALCULATION_DIR")
    ucs_out_level = os.getenv("UCS_OUTPUT_LEVEL")

    if ucs_out_level in {"FULL", "DEBUG"}: 

        def calculate_in_subdir(atoms: ase.Atoms):
            cwd = os.getcwd()
            #calc_dir = calc_dir
            #ucs_out_level = ucs_out_level
    
            sub_dir_name = '_'.join(
                ['{:.3f}'.format(c)
                 for c in atoms.get_scaled_positions()[-1]])
            sub_dir = Path(calc_dir, sub_dir_name)
    
            Path(calc_dir).mkdir(exist_ok=True, parents=True)
            sub_dir.mkdir(exist_ok=False)
            os.chdir(sub_dir)
    
            ase.io.write('atoms.cif', atoms)
            potential_energy = calculate_energy(atoms)
    
            os.chdir(str(cwd))
            return potential_energy

    elif ucs_out_level in {"LOW", "FIRST", "", None}:

        def calculate_in_subdir(atoms: ase.Atoms):
            cwd = os.getcwd()
            # New/Different: Access data in function attributes in __dict__
            fcn_dict = calculate_in_subdir.__dict__

            # Check the "iter" key in function attributes
            if "iter" in fcn_dict.keys():
                iter_status = fcn_dict["iter"]
            else:
                iter_status = None

            # Check the status
            if iter_status == "FIRST":

                sub_dir_name = '_'.join(
                    ['{:.3f}'.format(c)
                     for c in atoms.get_scaled_positions()[-1]] +
                     ["FIRST", "CALC"])

                sub_dir = Path(calc_dir, sub_dir_name)
                Path(calc_dir).mkdir(exist_ok=True, parents=True)
                sub_dir.mkdir(exist_ok=False)
                os.chdir(sub_dir)
            #elif iter_status = "LAST":

            #    sub_dir_name = '_'.join(
            #        ['{:.3f}'.format(c)
            #         for c in atoms.get_scaled_positions()[-1]],
            #         "LAST")

            #    sub_dir = Path(calc_dir, sub_dir_name)
            #    Path(calc_dir).mkdir(exist_ok=True, parents=True)
            #    sub_dir.mkdir(exist_ok=False)

            else:
                sub_dir_name = 'LATEST_CALC'

                sub_dir = Path(calc_dir, sub_dir_name) 
                Path(calc_dir).mkdir(exist_ok=True, parents=True)
                sub_dir.mkdir(exist_ok=True)
                os.chdir(sub_dir)
                with open("grid_point.txt",'w') as f:
                    f.write(" ".join(["Grid point sampled"] +  ['{:.3f}'.format(c)
                     for c in atoms.get_scaled_positions()[-1]]))

            ase.io.write('atoms.cif', atoms)
            potential_energy = calculate_energy(atoms)
    
            os.chdir(str(cwd))
            return potential_energy

    return calculate_in_subdir
