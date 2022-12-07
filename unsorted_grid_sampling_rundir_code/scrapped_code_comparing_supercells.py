# INFO: This code originally comes from custom_lammps_grid_gen.py, and
# server to investigate the differentce between a supercell read from
# a lammps datafile, and one constructed using ase from a corresponding
# unitcell. It did not resolve the differences in the end, but it was 
# (and is quite messy) and was removed from that script. 



##############################################################################
# Below here is just experimenting and comparing/checking the atomic
# coordinates. Lets ignore this for now-
# TODO: LATER, investigate more closely the difference between the 
# two sets of coordinates
##############################################################################
"""#print(atoms.get_positions()[0:20])
#print(wrapped_ref_atoms.get_positions()[0:20])
#print(atoms.get_chemical_symbols()[0:20])
#print(wrapped_ref_atoms.get_chemical_symbols()[0:20])
print(np.sum(np.isclose(atoms.get_positions(), wrapped_ref_atoms.get_positions(), rtol=1.0e-5)))
print(np.product(np.array(atoms.get_positions()).shape))

print((wrapped_ref_atoms.get_positions() - atoms.get_positions()) > 1.0e-4)

print(unitcell.get_positions())



### compare sorted coordinates
sorted_atoms_coords = vector_sort(np.copy(np.array(atoms.get_positions())))

sorted_wrapped_ref_atoms = vector_sort(np.copy(np.array(wrapped_ref_atoms.get_positions())))

print("Are sorted coords close (rtol = 1.0e-2)?: ", np.isclose(sorted_atoms_coords, sorted_wrapped_ref_atoms, 1.0e-2).all())
#print(((sorted_wrapped_ref_atoms - sorted_atoms_coords) < 1.0e-4))
#print(((sorted_wrapped_ref_atoms - sorted_atoms_coords) < 1.0e+0.all())
print("Maximum cooridnate error: ", np.max(np.abs(sorted_atoms_coords - sorted_wrapped_ref_atoms)))
print(wrapped_ref_atoms.get_cell_lengths_and_angles())
print(sorted_atoms_coords[0:20])
print(sorted_wrapped_ref_atoms[0:20])
#print(wrapped_ref_atoms.get_positions()[0:20])

for i in range(len(sorted_wrapped_ref_atoms)-1):
    #print("Comparing: \n",sorted_wrapped_ref_atoms[i],sorted_wrapped_ref_atoms[i+1])
    assert vector_xyz_leq(sorted_wrapped_ref_atoms[i],sorted_wrapped_ref_atoms[i+1]), "Sorted array is not sorted, check index "+str(i)+"!"

write_cubes = False
if write_cubes:
    with open("ref_atoms_test.cube", "w") as f:
        write_cube(f, ref_atoms, data=None)
    
    with open("atoms_test.cube", "w") as f:
        write_cube(f, atoms, data=None)

"""

"""points_in_atoms_not_in_ref_atoms = []
points_in_ref_atoms_not_in_atoms = []
in_ref_atoms = True
bol_ref_points_not_in_atoms = [1] * ref_atoms.get_positions().shape[0]"""

"""for point in atoms.get_positions():
    for indx, ref_point in enumerate(ref_atoms.get_positions()):
        in_ref_atoms = False
        
        if np.isclose(point, ref_point, rtol=1.0e-2).all():
            bol_ref_points_not_in_atoms[indx] = 0
            in_ref_atoms = True
            continue
        
    if not in_ref_atoms:
        points_in_atoms_not_in_ref_atoms.append(point)


points_in_ref_atoms_not_in_atoms = [ref_atoms.get_positions()[i,:] for i in range(len(bol_ref_points_not_in_atoms)) if bol_ref_points_not_in_atoms[i] == 1]

print(points_in_atoms_not_in_ref_atoms)
print(points_in_ref_atoms_not_in_atoms)"""
#set_of_positions_1 = set()
#set_of_positions_2 = set()
#for coord in atoms.get_positions():
#    #print(coord)
#    set_of_positions_1.add(tuple(np.round(coord, 4)))

#for coord in ref_atoms.get_positions():
#    #print(coord)
#    set_of_positions_2.add(tuple(np.round(coord, 4)))

#print(set_of_positions_1.symmetric_difference(set_of_positions_2))
#print(set_of_positions_1)
#print(set(list(atoms.get_positions())))
# Debug: THIS SEEMS TO WORK
#print(np.array(supcell_atoms.get_cell()))
#print(np.array(unitcell.get_cell()))
#print(C)
#print(C @ unitcell.get_cell())
#print(c)

#print("cell params A:", unitcell.get_cell_lengths_and_angles())
#print("cell params B:", supcell_atoms.get_cell_lengths_and_angles())
#print("Fraction between them: ", np.divide(supcell_atoms.get_cell_lengths_and_angles()[0:3], unitcell.get_cell_lengths_and_angles()[0:3]))
#
##############################################################################