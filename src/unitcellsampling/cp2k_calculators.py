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


# Function converting a cp2k-calculator to an ucs calculator
# (i.e. function with singature ase.Atoms -> float
def cp2k2ucs(cp2k_calc, base_calc_dir="UCS_CALCS/cp2k_calc", base_label="cp2k",
             update_mode="label", restart=True):
    """Embeds an ase CP2K calculator into a function
    with signature atoms -> float,  required for being
    a ucs calculator. Also provides an update scheme for label and 
    directories for individual calculations, and an iteration
    tracker.
    update_mode = label, dir or both. Both is default.
    Other options not implemented."""

    # TODO: Possibly implement file-existence check here
    
    # if update_mode == "label_only":
    label_switch = 1
    dir_switch = 1
    
    if update_mode == "label":
        dir_swtch = 0
    if update_mode == "dir":
        label_switch = 0
         

    def ucs_calc(atoms):
        ucs_calc.iter += 1
         
        calc_label = base_label + label_switch * ("_" + str(ucs_calc.iter))
        calc_dir = base_calc_dir + dir_switch * ("_" + str(ucs_calc.iter))

        # New simple restart functionality - to be tested:
        # Uses the wfn-file of the previous calculations as restart
        if restart and ucs_calc.iter > 1:
            prev_calc_dir = base_calc_dir + dir_switch * ("_" 
                                          + str(ucs_calc.iter - 1))
            prev_label = base_label + label_switch * ("_" 
                                    + str(ucs_calc.iter - 1))
            restart_file = prev_calc_dir + "/" + prev_label + "-RESTART.wfn"

            # The setting needs to be checked:
            # When it works, embed into function!
            # Pseudo code: get original input as str -> parse it to InputSection -> add_c2pk_restart_file to it -> transform back into input -> overwrite input.
            inp_section = ase.calculators.cp2k.parse_input(ucs_calc.init_inp)
            inp_section.add_keyword('FORCE_EVAL/DFT', 'WFN_RESTART_FILE_NAME ' + restart_file)
            inp_section.add_keyword('FORCE_EVAL/DFT/SCF', 'SCF_GUESS RESTART')
            new_inp = "\n".join(inp_section.write()) # It should be newline!
            cp2k_calc.set(inp=new_inp)

        cp2k_calc.label = calc_label
        cp2k_calc.directory = calc_dir
        
        cp2k_calc.calculate(atoms=atoms, properties=['energy'])

        energy = cp2k_calc.results['energy']

        return energy

    ucs_calc.iter = 0
    ucs_calc.init_inp = cp2k_calc.__dict__['parameters']['inp']
    return ucs_calc


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
# TODO: Switch this into the tested one
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


def cp2k_calculator_from_input(cp2k_input, **kwargs):
    """General energy calculator from parsed input from a cp2k input file.
    The input must be provided as a single string. Need to test parsing.
    Can also specify additional kwargs that are passed to the CP2K
    calculator; these are otherwise set to None, which will prohibit the ase
    caluclator from automatically generating undesired input sections."""

    accepted_kwargs = ['auto_write',
        'basis_set',
        'basis_set_file',
        'charge',
        'cutoff',
        'force_eval_method',
        'max_scf',
        'potential_file',
        'pseudo_potential',
        'stress_tensor',
        'uks',
        'poisson_solver',
        'xc',
        'print_level',
        'debug',
       ]

    cp2k_kwargs = dict(
        auto_write=False,
        basis_set=None,
        basis_set_file=None,
        charge=None,
        cutoff=None,
        force_eval_method=None,
        max_scf=None,
        potential_file=None,
        pseudo_potential=None,
        stress_tensor=None,
        uks=None,
        poisson_solver=None,
        xc=None,
        print_level=None,
        debug=False
        )

    assert accepted_kwargs == list(cp2k_kwargs.keys())
    assert set(kwargs.keys()).issubset(set(accepted_kwargs))

    cp2k_kwargs.update(kwargs)

    cp2k_calc = CP2K(inp=cp2k_input, **cp2k_kwargs)

    return cp2k_calc


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



