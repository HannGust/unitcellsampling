#!/usr/bin/env python

"""CLI to the unicellsampling batch analyzer."""

import unitcellsampling.batch_analyzer as ucs_ba
from unitcellsampling.batch_analyzer import UCSBatchAnalyzer

import argparse

#
def init_batch_analysis_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--directory", action="store", type=str, default=".", help="The path to the batch run main working directory, of the batch run to be analyzed.")
    parser.add_argument("-g", "--grid", action="store_true", help="Compile an energy grid to from the individual batch results of the batch run, and write it to a cube file.")
    return parser
#

def args2batch_analysis_options(args):
    """Converts parsed arguments into dictionary of options for 
    """
    opt_dict = {"batch_wd":args.directory,
                "grid":args.grid}
    return opt_dict



def main():
    parser = init_batch_analysis_parser()
    args = parser.parse_args()
    options = args2batch_analysis_options(args)
    print("analyze_batch.py: Analyzing batch run with work dir:", options["batch_wd"])
    batch_analyzer = UCSBatchAnalyzer(options["batch_wd"])

    if options["grid"]:
        print("Compiling grid...")
        batch_analyzer.compile_grid()
        print("Grid compiled.")
        print("Writing grid to cube-file...")
        batch_analyzer.write_grid2cube()
        print("Cube file written.")
    


if __name__ == "__main__":
    main()