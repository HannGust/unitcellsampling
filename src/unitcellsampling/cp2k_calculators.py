"""This module contains cp2k-based energy calculators and associated
   wrapper functions for unitcellsampling."""

import ase
from ase.calculators.cp2k import CP2K
from unitcellsampling.decorators import subdir_calc

# Decorator to set need to pass charge attribute
# TODO: Test this, in conjuction with the dft methods
def requires_charge_attr(method):
    method.requires_charge = True
    return method

# general template for attribute setting decorator
def set_requires_attr(method, attr):
    req_attr =  "requires_" + attr
    method.__dict__[req_attr] = True
    return method


# We need to define some wrapper functions to pass info to cp2k.
# These are actually not decorators, but still functions that
# return functions with desired parameters.
# Wrapper function for passing the charge
# TODO: Test this, in conjuction with the dft calculators.
def pass_cp2k_charge(cp2k_calc, cp2k_total_charge:int):
    """Wrapper function passing the total charge to e.g. a dft calculator 
    based on CP2K."""
    # Only for debug:
    print("Passing CP2K total charge to method:", cp2k_total_charge)
    def cp2k_calculator(atoms:ase.Atoms, charge=cp2k_total_charge):
        return cp2k_calc(atoms, charge=charge)

    return cp2k_calculator


# Wrapper function for passing parsed CP2K input
# TODO: Test this.
def pass_cp2k_input(cp2k_calc, parsed_cp2k_input):
    """Wrapper function to pass parsed text input from an input file to
    an ucs calculator based on CP2K."""
    def cp2k_calculator(atoms:ase.Atoms, cp2k_input=parsed_cp2k_input):
        return cp2k_calc(atoms, cp2k_input=cp2k_input)
    return cp2k_calculator


##########################################################################
# The following template is from the original energy_calculators-module
# where it was called cp2k_dft. This is presently only meant as a template
##########################################################################
#@subdir_calc
#def cp2k_dft_simple(atoms: ase.Atoms):
#    calc = CP2K(inp='''&FORCE_EVAL
#    &DFT
#        LSD
#    &END DFT
#&END FORCE_EVAL''')
#    atoms.set_calculator(calc)
#    return atoms.get_potential_energy()
##########################################################################

### CP2K through ase info ###
# The kwargs list for CP2K()-calculator object constructor is as follows, 
# with types and defaults
# See https://wiki.fysik.dtu.dk/ase/ase/calculators/cp2k.html for reference. 
# (Accessed 2023-09-19 @ 18:12)
###
"""
complete_cp2k_kwargs = {
        "auto_write": None, False
            "basis_set": str, DZVP-MOLOPT-SR-GTH
            "basis_set_file": str, BASIS_MOLOPT
            "charge": float, 0
            "command": str, CP2K.command -> $ASE_CP2K_COMMAND -> cp2K_shell
            "cutoff": float, 400 * Rydberg
            "debug": bool, False
            "force_eval_method": str, Quickstep
            "inp": str, -
            "max_scf": int, 50
            "poisson_solver": str, auto
            "potential_file": str, POTENTIAL
            "pseudo_potential": str, auto -> GTH-PBE
            "stress_tensor": bool, True
            "uks": bool, False
            "xc": str, LDA
            "print_level": str, LOW
            }
"""
### End CP2K through ase info ###

# CP2K purely from parsed input
# TODO: Test this
# TODO: See if some keywords should be set to "None" to prevent 
# interference of the default template used by the ase CP2K-calculator
# object.
@subdir_calc
def cp2k_from_input(atoms: ase.Atoms, cp2k_input=""):
    """General energy calculator from parsed input from a cp2k input file.
    The input must be provided as a single string."""
 
    calc = CP2K(inp=cp2k_input)
    atoms.set_calculator(calc)
    return atoms.get_potential_energy()



# Ordinary cp2k_dft_pbe, restricted spin
# TODO: Test this
@requires_charge_attr
@subdir_calc
def cp2k_dft_pbe(atoms: ase.Atoms, charge=None):
    """CP2K spin-restricted DFT calculator using the PBE-functional,
    with the default pseudo-potential (auto, based on functional) 
    and basis set (DVZP-MOLOPT-SR-GTH). Most parameters are explicitly
    set to the default of the ase CP2K-calculator object, but max_scf
    is set to a higher value of 300 (instead of 50)."""

    
    cp2k_kwargs = {
            "basis_set": "DZVP-MOLOPT-SR-GTH",
            "charge": charge,
            "force_eval_method": "Quickstep",
            "max_scf": 300,
            "pseudo_potential": "auto",
            "stress_tensor": True,
            "uks": False,
            "xc": "PBE",
            "print_level": "LOW"
            }
 
    calc = CP2K(**cp2k_kwargs)

    atoms.set_calculator(calc)
    return atoms.get_potential_energy()


# cp2k_dft_uks_pbe, unrestricted spin-polarised pbe
# TODO: Test this
@requires_charge_attr
@subdir_calc
def cp2k_dft_uks_pbe(atoms: ase.Atoms, charge=None):
    """CP2K unrestricted, spin-polarised UKS-DFT calculator using the 
    PBE-functional, with the default pseudo-potential (auto, based on 
    functional) and basis set (DVZP-MOLOPT-SR-GTH).
    Equivalent to cp2k_dft_pbe, except for being spin-polarised (UKS).
    That is, again most parameters are set explicitly to the default
    of the ase CP2K-calculator object, but max_scf is set to 300."""
    
    cp2k_kwargs = {
            "basis_set": "DZVP-MOLOPT-SR-GTH",
            "charge": charge,
            "force_eval_method": "Quickstep",
            "max_scf": 300,
            "pseudo_potential": "auto",
            "stress_tensor": True,
            "uks": True,
            "xc": "PBE",
            "print_level": "LOW"
            }
 
    calc = CP2K(**cp2k_kwargs)

    atoms.set_calculator(calc)
    return atoms.get_potential_energy()



