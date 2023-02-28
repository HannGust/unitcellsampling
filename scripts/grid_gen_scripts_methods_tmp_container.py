### THIS CONTAINS THE METHODS THAT WAS ORIGINALLY DEFINED DIRCETLY IN THE grid_gen_script.py-script
# THIS IS ATEMPORARY CONTAINER: IT contain the code as-is. It should be cleaned up and put into a 
# properly structured file if it should be used / H.

# TODO: Implement method to read lammps parameters from file - decorator maybe
##########  Define methods for computing energy here, like so /H  ########## 

## CP2K DFT (PBE functional)
@subdir_calc
def dft_pbe(atoms: ase.Atoms):
    global restart_file_name

    # Use the same parameter here as in `single_point_calculation/` ## Need to check how to structure this, the input /H
    ## Remember to alter the charge to the correct one
    inp = '''
&FORCE_EVAL
  &DFT
    CHARGE -11
    LSD
    RESTART_FILE_NAME {}
    &SCF
       EPS_SCF 5.00000000E-005
       SCF_GUESS RESTART
       &OT T
           PRECONDITIONER FULL_SINGLE_INVERSE
           MINIMIZER CG
           LINESEARCH 3PNT
       &END OT
    &END SCF
  &END DFT
&END FORCE_EVAL
    '''.format(restart_file_name)

    # Change the restart file to the previous calculation ## As a note, this makes parallellization of calculations impossible (unless using initial calculation) /H
    p = Path('..').glob('**/*')
    dirs = [x for x in p if x.is_dir()]
    restart_file_name = str(
        Path(max(dirs, key=os.path.getmtime), 'LTA4A_reduced-RESTART.wfn'))

    calc = CP2K(inp=inp, pseudo_potential='GTH-PBE',
                max_scf=600, xc='PBE', print_level='LOW')
    atoms.set_calculator(calc)
    return atoms.get_potential_energy()


## CP2K Lennard-Jones
@subdir_calc
def cp2k_lj(atoms: ase.Atoms):
    # TODO: Finish this energy evaluation method (the input is half done at this point...)
    inp = '''
&FORCE_EVAL
  &MM
    &FORCEFIELD
      &NONBONDED
