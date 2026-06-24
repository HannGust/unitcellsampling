#!/bin/bash

#SBATCH -A naiss2025-1-25
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH -t 00-04:00:00
#SBATCH -J LGPS_example_python

# It appears to work to run sbatch with the environment activated.
# However, I activate it here like so:
source /proj/teoroo/users/x_hagus/venv/ucs_env_venv_2024-01-25/bin/activate # NOTE: python environment where unitcellsampling is installed

# Make sure CP2K is available (load module in cluster environment)
module load CP2K/2023.1-psmp-PLUMED

## Environment Vars:
export ASE_CP2K_COMMAND="mpprun -q cp2k_shell.psmp" # 
##

# The sampler will fill in code to run below:
