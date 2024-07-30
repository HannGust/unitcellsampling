


from unitcellsampling.batch_analyzer import read_batch_log, get_grid_shape_from_batch_log, get_symmetry_from_batch_log



# Log file to test on
testlog = "./tests/ucs_batch.log"

# Read it!
log_txt = read_batch_log(testlog)

# If you want to test, print it:
print(log_txt)

# Try to get the grid shape:
grid_shape = get_grid_shape_from_batch_log(log_txt)

print("Extracted grid shape: ", grid_shape)


# Not try to get symmetry settings and spacegroup

symmetry_setting, spgrp = get_symmetry_from_batch_log(log_txt)


print("Symmetry setting found:", symmetry_setting)
print("Spacegroup found:", spgrp, spgrp.number, spgrp.hall)

