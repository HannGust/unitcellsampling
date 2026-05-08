#!/bin/bash

# Example of how the batch sampler can be used
# This is an example of a runscript

cif="./54189.cif" # Cif file of the input structure
method="A54189"    # This is a classical force field specifically for Li2ZnGe based on UFF lj-parameters and REPEAT charges

# Here we setup and automatically start a batch run:
# Name is given by --name (this will be applied to the calculation directory)
# -a "Li" gives Li as the sample atom
# -s 0.2: 0.2 Å grid spacing in all directions
# --ra: remove existing sampling atoms (Li) prior to sampling
# --vdw 0.5: use vdW cutoff exclusion with scaling factor 0.5

python ../../scripts/ucs_batch_run.py $cif $method \
       	                              --name "Li2GeZn_example" \
				      -a "Li" \
				      -s 0.2 \
				      --ra \
				      --vdw 0.5 \
				      --conv \
				      --batch-size 1000 \
				      --jobscript-cmd bash \
				      --jobscript-template jobscript.sh


