import numpy as np
from numpy import random
import math as ma

"""A =np.array([[0., 0., 0.],
     [0., 0., 1.],
     [0., 2., 2.],
     [1., 2., 0.],
     [1., 0., 2.],
     [1., 1., 1.],
     [0., 1., 0.]])
"""

A =np.array([[0., 0., 0.],
             [1., 1., 2.],
             [0., 0., 1.],
             [1., 1., 1.],
             [0., 0., 2.],
             [2., 1., 1.],
             [0., 1., 0.],
             [0., 1., 2.],
             [0., 2., 1.],
             [2., 2., 1.],
             [0., 2., 2.],
             [1., 0., 0.],
             [1., 0., 1.],
             [2., 1., 0.],
             [1., 0., 2.],
             [1., 1., 0.],
             [1., 2., 0.],
             [0., 2., 0.],
             [1., 2., 1.],
             [1., 2., 2.],
             [0., 1., 1.],
             [2., 0., 0.],
             [2., 0., 1.],
             [2., 0., 2.],
             [2., 1., 2.],
             [2., 2., 0.],
             [2., 2., 2.]])


block_len = 3
A_shp = A.shape

indxs_0 = np.argsort(A[:,0], axis=0)
print(indxs_0)
A_sorted_0 = np.copy(A)[indxs_0]
print(A_sorted_0)
indxs_1_list = []
A_sorted_1 = np.copy(A_sorted_0)
for i in range(0, block_len**3, block_len**2):
    indxs_1 = np.argsort(A_sorted_1[i:i+block_len**2, 1])
    A_sorted_1[i:i+block_len**2] = (A_sorted_1[i:i+block_len**2])[indxs_1]

A_sorted_2 = np.copy(A_sorted_1)
for j in range(0, block_len**3, block_len):
    indxs_2 = np.argsort(A_sorted_2[j:j+block_len, 2])
    A_sorted_2[j:j+block_len] = (A_sorted_2[j:j+block_len**2])[indxs_2]


"""indxs_1_0 = np.argsort(A_sorted_0[0:4, 1])
indxs_1_1 = np.argsort(A_sorted_0[4:, 1])
A_sorted_1 = np.copy(A_sorted_0)
A_sorted_1[0:4] = A_sorted_1[0:4][indxs_1_0]
A_sorted_1[4:] = A_sorted_1[4:][indxs_1_1]
"""
print(A_sorted_2)
A_sorted = np.sort(A, axis=0)

### sort_coords

# Sorts (n, 3) arrays which are meant to be list of 3-dimensional coordinates.
# Sorts with respect to x-coordinate first, y-coordinate second, z-coordinate last.

# array : np.array of shape (n, 3), to be sorted.

# block_size : Tuple of three ints, giving the block size in that dimension,
#              i.e. how many different values of that coordinate for each
#              pair of values of the other coordinates. Note that this assumes
#              somewhat regular structure of the array. 
#              If None, takes the block sizes to be equal in each dimension,
#              and equal to the cuberoot of the number of points, as in e.g.
#              a cartesian grid with the same number of points k in each
#              dimension (k * k * k).


def sort_coords(array, block_size=None):
    assert len(array.shape) == 2 and array.shape[1] == 3, "Array must be 2-dimensional of shape (n, 3)"

    if not block_size:
        block_size = int(ma.pow(array.shape[0], 1.0/3.0))
        print(block_size, array.shape[0])
        assert block_size**3 == array.shape[0], "Error, truncation or rounding error in cuberoot calculation."
        block_size = (block_size,) * 3
    else:
        assert np.product(block_size) == np.product(array.shape), "Error, block_sizes does not add up to total number of points."

    array_sorted_0 = np.copy(array)[np.argsort(array[:, 0], axis=0)]

    array_sorted_1 = np.copy(array_sorted_0)
    for i in range(0, np.product(block_size), np.product(block_size[1:])):
        np.argsort(array_sorted_1[i:i + np.product(block_size[1:]), 1])
        
        array_sorted_1[i:i + np.product(block_size[1:])] = (
                array_sorted_1[i:i + np.product(block_size[1:])]
                )[np.argsort(array_sorted_1[i:i + np.product(block_size[1:]), 1])]

    array_sorted_2 = np.copy(array_sorted_1)
    for i in range(0, np.product(block_size), np.product(block_size[2:])):
        np.argsort(array_sorted_2[i:i + np.product(block_size[2:]), 2])
        
        array_sorted_2[i:i + np.product(block_size[2:])] = (
                array_sorted_2[i:i + np.product(block_size[2:])]
                )[np.argsort(array_sorted_1[i:i + np.product(block_size[2:]), 2])]


    return array_sorted_2


test_sort_coords_array = sort_coords(A)
print(test_sort_coords_array)
#print(A)

#print("\n\n")

#print(A_sorted)
# Compares vectors in order x, y, z and returns true if x1<=x2, y1<=y2 and z1<=z2 all holds.
def vector_xyz_leq(v1, v2):
    v1_test = np.copy(v1).squeeze()
    v2_test = np.copy(v2).squeeze()
    assert len(v1_test) == 3 and len(v2_test) == 3, "Vectors must be of length 3."
    if v1_test[0] < v2_test[0]:
        return True
    elif v1_test[0] > v2_test[0]:
        return False
    elif v1_test[0] == v2_test[0]:
        if v1_test[1] < v2_test[1]:
            return True
        elif v1_test[1] > v2_test[1]:
            return False
        elif v1_test[1] == v2_test[1]:
            if v1_test[2] < v2_test[2]:
                return True
            elif v1_test[2] > v2_test[2]:
                return False
            elif v1_test[2] == v2_test[2]:
                return True
            else:
                raise Exception("This should not have happened, check logic and input.")
        else:
            raise Exception("This should not have happened, check logic and input.")
    else:
        raise Exception("This should not have happened, check logic and input.")



def vector_sort(array):
    iterations = 0
    assert len(array.shape) == 2 and array.shape[1] == 3, "Array must be 2-dimensional of shape (n, 3)"
    sorted_array = np.round(np.copy(array),5)
    for i in range(1,len(array)):
        j = i
        while j > 0 and (not vector_xyz_leq(sorted_array[j-1], sorted_array[j])):
            if iterations < 20:
                #print("Comparing ", sorted_array[j-1], " with ", sorted_array[j], " and decided to swap.")
                iterations += 1
            sorted_array[j-1], sorted_array[j] = np.copy(sorted_array[j]), np.copy(sorted_array[j-1])
            j -= 1
        #print("Done swapping ", array[i])

    return sorted_array



def general_n3_sort(array):
    assert len(array.shape) == 2 and array.shape[1] == 3, "Array must be 2-dimensional of shape (n, 3)"
    
    xblock_sizes = []
    xblock_lims = []

    yblock_sizes = []
    zblock_sizes = []
    
    sort_x_indx = np.argsort(array[:, 0])

    x_sorted_array = np.copy(array)[sort_x_indx]
    
    first_iter = True
    num_x_blocks = 0

    for i,x in enumerate(x_sorted_array):
        if first_iter:
            first_iter = False
            curr_block_value = x[0]
            num_x_blocks += 1
            temp_start_lim = 0
            temp_block_size = 1
        else:
            if x[0] != curr_block_value:
                temp_stop_lim = i - 1
                xblock_lims.append([temp_start_lim, temp_stop_lim])
                xblock_sizes.append(temp_block_size)
                temp_block_size = 1
                num_x_blocks += 1
                temp_start_lim = i
                curr_block_value = x[0]
            else:
                temp_block_size += 1

            if (i == (len(x_sorted_array) - 1)):
                xblock_sizes.append(temp_block_size)
                xblock_lims.append([temp_start_lim, i])

    assert num_x_blocks == len(xblock_sizes) and num_x_blocks == len(xblock_lims) 

    # Now sort each x-block according to y
    y_sorted_array = np.copy(x_sorted_array)
    for i,lims in xblock_lims:
        tmp_y_block = np.copy(x_sorted_array[lims[0]:lims[1]+1])
        temp_sort_y_indx = np.argsort(tmp_y_block, 1)
        y_sorted_block = np.copy(tmp_y_block[temp_sort_y_indx])
        
        y_sorted_array[lims[0], lims[1]+1] = y_sorted_block
    
    #   for j,y in enumerate(y_sorted_block):
    pass


    #return num_x_blocks, xblock_sizes, xblock_lims, x_sorted_array


def sort_atom_coords(coords):
    coords = np.array(coords).reshape(-1,3)
    assert len(coords.shape) == 2 and (coords.shape[1] == 3), "Coordinate array must have shape (n, 3)."
    
    sort_array = np.vstack(( coords.T[1, :], coords.T[0, :]))

    sorting_indx = np.lexsort(sort_array)

    sorted_coords = np.copy(coords)[sorting_indx]

    return sorted_coords


A_sorted_by_lexsort = sort_atom_coords(A)

print("\n\n\n\n")
print("A sorted using lexsort: ")
print(A_sorted_by_lexsort)
