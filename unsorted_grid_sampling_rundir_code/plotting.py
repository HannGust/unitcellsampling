import matplotlib.pyplot as plt
import matplotlib as mpl
#import mpl_toolkits.mplot3d as 3dplt
import numpy as np
import ase
from ase.io.cube import read_cube
from compare_grids import custom_read_cube as r_cube
import sys
import itertools as it


#####
# Program to easily plot colormapped plots of grids from commandline
# Takes (currently) one argument which is the name of a cubefile containing
# the atoms object and grid data to plot
#####


plot_cart = True # Boolean value determineing if plotting will be done w.r.t. cartesian coordinates. If false, fractional coordinates are used.
# We want to plot a grid, as points, colored according to some gradient based ontheir value
# Need to read the 
# Plots the grid in 3D space, either in fractional or cartesian coordinates, and applies a colormap to the points
def plot_grid(atoms, grid, frac_coords, colorbar=None, plot_cart=True, cmap=None, norm=None, vmin=None, vmax=None, **kwargs): # Might be useful to define a fcn

    #grid -> assume read from the file
    energies = grid.reshape((-1,1))

    # Should be fractional coordinates -> can either plot in fractional or cartesian coords. Fractional is simple, but cartesian gives a realistic depiction
    cart_coords = np.array(frac_coords) @ np.array(atoms.get_cell())

    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')

    # Specify the colormap
    #cmap = 'viridis'
    #"norm = None # whether to normalize the colormap so that min = 0 and max = 1
    #"vmin, vmax = None, None # Lower and upper values used in the colormapping, in case norm is not given

    if plot_cart:
        ax_sca = ax.scatter(cart_coords[:,0], cart_coords[:,1], cart_coords[:,2], c=energies, cmap=cmap, norm=norm, s=100, **kwargs)

        cell_corners = np.array(list(it.product([0,1], repeat=3))) @ atoms.get_cell()
        print(cell_corners)
        xmin, xmax = np.min(cell_corners[:,0]) - 1.0, np.max(cell_corners[:,0]) + 1.0
        ymin, ymax = np.min(cell_corners[:,1]) - 1.0, np.max(cell_corners[:,1]) + 1.0
        zmin, zmax = np.min(cell_corners[:,2]) - 1.0, np.max(cell_corners[:,2]) + 1.0

        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
        ax.set_zlim(zmin, zmax)
    else:
        ax_sca = ax.scatter(frac_coords[:,0], frac_coords[:,1], frac_coords[:,2], c=energies, cmap=cmap, norm=norm, s=100, **kwargs)
        ax.set_xlim(-0.1, 1.1)
        ax.set_ylim(-0.1, 1.1)
        ax.set_zlim(-0.1, 1.1)
        print(ax.get_xlim())

  
    fig.colorbar(ax_sca, ax=ax)
    #fig.show()
    fig.savefig("grid_plot_from_plotting.png")
    input("Press enter button to continue")
    return

if __name__ == "__main__":
    if len(sys.argv) > 1:
        filename = sys.argv[1]
    else:
        filename = None
    
    if filename:
        #(atoms, grid) = r_cube(filename)
        with open(filename, "r") as cf:
             cube_cont = read_cube(cf)
             atoms = cube_cont["atoms"]
             grid = cube_cont["data"]

        g_inf = np.max(grid)
        mask = grid < g_inf
    
    
        lower_max = np.max(grid[mask])
        print(g_inf, lower_max)
        modified_grid = np.ones_like(grid)*2*lower_max
        modified_grid[mask] = grid[mask]
        nx,ny,nz = grid.shape
        frac_coords = np.array(list(it.product(np.linspace(0, 1, nx, endpoint=False), np.linspace(0, 1, ny, endpoint=False), np.linspace(0, 1, nz, endpoint=False))))
    else:
        #atoms = 
        #grid = 
        raise NotImplementedError("Need input.")
    
    
    #plot_grid(atoms, grid, frac_coords)
    plot_grid(atoms, modified_grid, frac_coords)
    
    grid_list = grid.flatten()
    g_min_list = []
    min_indices_list = []
    
    cart_coords = frac_coords @ atoms.get_cell()
    
    if plot_cart:
        plot_coords = cart_coords
    else:
        plot_coords = frac_coords
    
    for i in range(1):
       if i == 0:
           min_mod_grid = grid_list 
       g_min = np.min(min_mod_grid)
       print(g_min)
       g_min_list.append(g_min)
       
       min_indices = (plot_coords.reshape(-1,3))[(grid_list == g_min),:]
       
       print(min_indices.reshape(-1,3))
       min_indices_list.append((min_indices.reshape(-1,3)))
       min_mod_grid = grid_list[grid_list > g_min]
    
    print('g min list ', g_min_list)
    print('min_ind_list ', min_indices_list)
    
    
    # energy_array_list = [np.ones((min_indices_list[i].shape[0])).fill(g_min_list[i]) for i in range(len(g_min_list))] THIS DOES NOT WORK FOR SOME REASON -> BECAUSE FILL IS A METHOD AND DOESN RETUNR ANYTHING... IT MODIFIES THE OBJECT AS IS... 
    #print('energy array list: ', energy_array_list)
    energy_array_list = [ np.ones(min_indices_list[i].shape[0]) * g_min_list[i] for i in range(len(g_min_list))]
    print('energy array list: ', energy_array_list)
    
    fig2 = plt.figure()
    ax2 = fig2.add_subplot(projection='3d')
    ax2.set_xlim(np.max(atoms.get_cell()[:,0]))
    ax2.set_ylim(np.max(atoms.get_cell()[:,1]))
    ax2.set_zlim(np.max(atoms.get_cell()[:,2]))
    
    cmap2 = plt.get_cmap('viridis')
    cnorm = mpl.colors.Normalize(vmin=g_min_list[0], vmax=g_min_list[-1])
    
    #vmin = g_min_list[0]
    #vmax = g_min_list[-1]
    
    
    for l in range(len(g_min_list)):
       ax2_sca = ax2.scatter(min_indices_list[l][:,0], min_indices_list[l][:,1], min_indices_list[l][:,2], c=energy_array_list[l], cmap=cmap2, norm=cnorm, s=1)
    
    fig2.colorbar(ax2_sca, ax=ax2)
    #fig2.show()
    fig2.savefig("fig2_grid_plot_from_plotting.png")
    input('press enter to continue')    
