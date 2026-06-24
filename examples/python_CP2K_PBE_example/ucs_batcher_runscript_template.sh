#!/bin/bash

# Runscript template for the ucs batch grid sampling

## For the example sampling Li in Li10Ge(PS6)2 (LGPS) with DFT (PBE)
# This shows an example with --no-run applied, so that the sampling
# directory is just setup, but no calculations are started.

###### Note on ucs batch cli settings: by running python ucs_batch_run.py --help, one can see arguments

### Change this accordingly: ###
cif="Li10Ge(PS6)2.cif"              # Cif of the structure to be sampled
method="cp2k_calculator_from_input" # cp2k calculation from input file template
nametag="EXAMPLE"                   # Applied name as a tag at the end of calculation directory

CP2K_INP_TEMPLATE="cp2k_template_singlet_many_kinds_specified.inp"  # Cp2k input template
JOBSCRIPT_TEMPLATE="ucs_batcher_jobscript_template.sh"              # Jobscript template, used to submit/start the individual batches

### Set the following appropriately, this is just an example ###
ucs_batch_settings="$cif $method -n $nametag \
                   -a Li \
                   -s 0.4 \
                   --rc 0.5 \
		   --vdw 0.5 \
                   --ra \
                   --conv \
		   --cp2k_q auto \
		   --cp2k_wfn_mode on \
		   --cp2k_template ${CP2K_INP_TEMPLATE} \
                   --batch-size 750 \
                   --jobscript-template ${JOBSCRIPT_TEMPLATE} \
                   --no-run"


# Here the calculation is run:
python -u ../../scripts/ucs_batch_run.py ${ucs_batch_settings}


