#!/bin/bash -l
 
#SBATCH -A snic2020-1-41
#SBATCH -N1 --exclusive
#SBATCH -t 08:00:00
#SBATCH -J new_LGPS_PBE_GRID

module load CP2K/6.1-psmp-PLUMED
# module load CP2K/5.1-nsc1-intel-2018a-eb
rm -rf ucs_calculation/

export UCS_WORK_DIR="."
export UCS_CALCULATION_DIR="./ucs_calculation"
export OMP_NUM_THREADS=32
export ASE_CP2K_COMMAND=""

./grid_gen.py > grid_gen.out
