#!/bin/bash

#SBATCH -A naiss2025-1-25
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH -t 00-04:00:00
#SBATCH -J LGPS_reftraj_example

module add CP2K/2024.1-psmp-PLUMED
mpprun cp2k.popt -i "cp2k.inp" > "cp2k.out"
