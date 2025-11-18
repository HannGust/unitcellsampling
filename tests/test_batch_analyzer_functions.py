#!/usr/bin/env python

"""Tests for functions in the batch analyser module,
e. g. for some batch log reading functions that extract 
information from batch logs."""


from unitcellsampling.batch_analyzer import read_batch_log, get_grid_shape_from_batch_log, get_symmetry_from_batch_log, get_cp2k_sampling_charge_from_log, get_sampling_atom_from_log,  get_sampling_atom_charge_from_log, get_sampling_supercell_size_from_log
import gemmi
from unitcellsampling.symmetry import compare_gemmi_spacegroups
import os

### Log files to test on
# For ucs_batch_run.py
testlog = os.path.realpath(os.path.join(os.path.dirname(__file__), "ucs_batch.log"))
true_info = {"grid_shape":(14,14,14),
             "symmetry":(True, list(gemmi.spacegroup_table_itb())[511]),
             "sampling_atom":"Li",
             "sampling_atom_q":1,
             "sampling_q":-47,
             "sampling_supercell_size":(1, 2, 3)
             }

# For ucs_batch_run_reftraj.py
testlog_reftraj = os.path.realpath(os.path.join(os.path.dirname(__file__), "ucs_batch_reftraj.log"))
true_info_reftraj = {"grid_shape":(4, 12, 9),
                     "symmetry":(False, None),
                     "sampling_atom":"Li",
                     "sampling_atom_q":2,
                     "sampling_q":-62,
                     "sampling_supercell_size":(2, 2, 1)
                     }

def test_read_batch_log(log):
    _ = read_batch_log(log)
    return

def test_get_grid_shape_from_batch_log(logtxt, ans):
    gs = get_grid_shape_from_batch_log(logtxt)
    assert gs == ans
    return

def test_get_symmetry_from_batch_log(logtxt, ans):
    spgrp, flag = get_symmetry_from_batch_log(logtxt)
    assert flag == ans[1]
    assert (spgrp is None and ans[0] is None) or compare_gemmi_spacegroups(spgrp,ans[0])
    return

def test_get_sampling_atom_from_log(logtxt, ans):
    sample_atom = get_sampling_atom_from_log(logtxt)
    assert sample_atom == ans
    return


def test_get_sampling_atom_charge_from_log(logtxt, ans):
    sample_atom_q = get_sampling_atom_charge_from_log(logtxt)
    assert sample_atom_q == ans
    return


def test_get_cp2k_sampling_charge_from_log(logtxt, ans):
    sample_q = get_cp2k_sampling_charge_from_log(logtxt)
    assert sample_q == ans
    return


def test_get_sampling_supercell_size_from_log(logtxt, ans):
    sample_ssc = get_sampling_supercell_size_from_log(logtxt)
    assert sample_ssc == ans
    return


test_map = {"grid_shape":test_get_grid_shape_from_batch_log,
            "symmetry":test_get_symmetry_from_batch_log,
            "sampling_atom":test_get_sampling_atom_from_log,
            "sampling_atom_q":test_get_sampling_atom_charge_from_log,
            "sampling_q":test_get_cp2k_sampling_charge_from_log,
            "sampling_supercell_size":test_get_sampling_supercell_size_from_log}

# test reading logs
test_read_batch_log(testlog)
test_read_batch_log(testlog_reftraj)

# Read it!
log_txt = read_batch_log(testlog)
log_txt_reftraj = read_batch_log(testlog_reftraj)


print("Explicitly printing some retrieved info (ucs_batch.log):")
print("-"*len("Explicitly printing some retrieved info (ucs_batch.log):"))
# Try to get the grid shape:
grid_shape = get_grid_shape_from_batch_log(log_txt)
print("Extracted grid shape: ", grid_shape)

# Not try to get symmetry settings and spacegroup
symmetry_setting, spgrp = get_symmetry_from_batch_log(log_txt)
print("Symmetry setting found:", symmetry_setting)
print("Spacegroup found:", spgrp, spgrp.number, spgrp.hall)

sample_atom = get_sampling_atom_from_log(log_txt)
atom_charge = get_sampling_atom_charge_from_log(log_txt)
sample_charge = get_cp2k_sampling_charge_from_log(log_txt)
sample_supercell_size = get_sampling_supercell_size_from_log(log_txt)
print(f"Sample atom, sample atom charge, and total sampling charge:\n{sample_atom} ,{atom_charge}, {sample_charge}")
print(f"Sampling supercell size: {sample_supercell_size}")
print("-"*len("Explicitly printing some retrieved info (ucs_batch.log):"))
print()


print("Explicitly printing some retrieved info (ucs_batch_reftraj.log):")
print("-"*len("Explicitly printing some retrieved info (ucs_batch_reftraj.log):"))
# Try to get the grid shape:
grid_shape = get_grid_shape_from_batch_log(log_txt_reftraj)
print("Extracted grid shape: ", grid_shape)

# Not try to get symmetry settings and spacegroup
symmetry_setting, spgrp = get_symmetry_from_batch_log(log_txt_reftraj)
print("Symmetry setting found:", symmetry_setting)
print("Spacegroup found:", spgrp)

sample_atom = get_sampling_atom_from_log(log_txt_reftraj)
atom_charge = get_sampling_atom_charge_from_log(log_txt_reftraj)
sample_charge = get_cp2k_sampling_charge_from_log(log_txt_reftraj)
sample_supercell_size = get_sampling_supercell_size_from_log(log_txt_reftraj)
print(f"Sample atom, sample atom charge, and total sampling charge:\n{sample_atom} ,{atom_charge}, {sample_charge}")
print(f"Sampling supercell size: {sample_supercell_size}")
print("-"*len("Explicitly printing some retrieved info (ucs_batch_reftraj.log):"))
print()



### Explicit testing here
def test_for_batch_log():
    logtxt = read_batch_log(testlog)
    ans_dict = true_info
    for key, test_func in test_map.items():
        test_func(logtxt, ans=ans_dict[key])

    # Should do this, but compactly written:
    #test_get_grid_shape_from_batch_log(logtxt, ans=true_info["grid_shape"])
    #test_get_symmetry_from_batch_log(logtxt, ans=true_info["symmetry"])
    #test_get_sampling_atom_from_log(logtxt, ans=true_info["sampling_atom"])
    #test_get_sampling_atom_charge_from_log(logtxt, ans=true_info["sampling_atom_q"])
    #test_get_cp2k_sampling_charge_from_log(logtxt, ans=true_info["sampling_q"])
    #test_get_sampling_supercell_size_from_log(logtxt, ans=true_info["sampling_supercell_size"])
    return

print("Performing test for ucs_batch.log, from ucs_batch_run.py")
test_for_batch_log()
print("Done testing for ucs_batch.log, from ucs_batch_run.py.")
print()

### Explicit testing here
def test_for_batch_log_reftraj():
    logtxt = read_batch_log(testlog_reftraj)
    ans_dict = true_info_reftraj
    for key, test_func in test_map.items():
        test_func(logtxt, ans=ans_dict[key])

    # Should do this, but compactly written:
    #test_get_grid_shape_from_batch_log(logtxt, ans=true_info_reftraj["grid_shape"])
    #test_get_symmetry_from_batch_log(logtxt, ans=true_info_reftraj["symmetry"])
    #test_get_sampling_atom_from_log(logtxt, ans=true_info_reftraj["sampling_atom"])
    #test_get_sampling_atom_charge_from_log(logtxt, ans=true_info_reftraj["sampling_atom_q"])
    #test_get_cp2k_sampling_charge_from_log(logtxt, ans=true_info_reftraj["sampling_q"])
    #test_get_sampling_supercell_size_from_log(logtxt, ans=true_info_reftraj["sampling_supercell_size"])
    return

print("Performing test for ucs_batch.log, from ucs_batch_run_reftraj.py")
test_for_batch_log_reftraj()
print("Done testing for ucs_batch.log, from ucs_batch_run_reftraj.py.")
print()




    

