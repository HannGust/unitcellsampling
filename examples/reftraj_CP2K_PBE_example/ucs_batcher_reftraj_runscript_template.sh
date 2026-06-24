#!/bin/bash

# Runscript template for the ucs batch REFTRAJ grid sampling
## For the example sampling Li in LGPS with DFT (PBE)

### Change this accordingly: ###
cif="Li10Ge(PS6)2.cif"              # Cif of the structure to be sampled
method="cp2k_calculator_from_input" # cp2k calculation from input file template
nametag="EXAMPLE_REFTRAJ"           # Applied name as a tag at the end of calculation directory

CP2K_INP_TEMPLATE="cp2k_template_reftraj_singlet.inp"            # Cp2k input template
JOBSCRIPT_TEMPLATE="ucs_batcher_reftraj_jobscript_template.sh"   # Jobscript template, used to submit/start the individual batches

### Set the following appropriately: ###

ucs_batch_settings="$cif $method -n $nametag \
                   -a Li \
                   -s 0.4 \
                   --rc 0.5 \
		   --vdw 0.5 \
                   --ra \
                   --conv \
                   --cp2k_q auto \
		   --cp2k_template ${CP2K_INP_TEMPLATE} \
                   --batch-size 750 \
                   --jobscript-cmd sbatch \
                   --jobscript-template ${JOBSCRIPT_TEMPLATE}"

# Here the calculation is run:
python -u ../../scripts/ucs_batch_run_reftraj.py ${ucs_batch_settings}


