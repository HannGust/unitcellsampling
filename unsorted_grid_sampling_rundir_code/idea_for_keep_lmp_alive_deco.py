
import os, sys, numpy as np, ase, ase.io
import decorator
# Decorator to keep lamps alive or not

# Is this a good idea?
# This could potentially also be managed through a global variable in the sampling code
# Like so: 
# if method is lammps pass keep_lammps_alive
# sampling_loop for all points do
#     compute energy keeping lammps alive
# end loop
# kill lammps
#
# Yes simply try with keep_alive = true in LAMMPSlib-calculator
@decorator
def lmp_lifesupport(calc):
    keep_lmp_alive = os.getenv("UCS_KEEP_LAMMPS_ALIVE")
    return
