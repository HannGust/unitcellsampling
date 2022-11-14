#!/bin/bash -l
 
#SBATCH -A snic2020-1-41
#SBATCH -n 16
#SBATCH -t 00:20:00
#SBATCH -J LGPS_PBE_GRID

#module load CP2K/7.1-psmp-PLUMED
#module load CP2K/5.1-nsc1-intel-2018a-eb
module load CP2K/6.1-psmp-PLUMED

OMP_NUM_THREADS=32
OMP_PROC_BIND=spread

srun cp2k.psmp cp2k.inp
