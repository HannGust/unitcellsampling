""" Utility plotting functions for grids. """

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import ase
from ase.io.cube import read_cube
import sys # sys for input handling for now...


default_kwargs = {'axs':"z", 'num':0, 'start_x':0, 'end_x':None,
                    'start_y':0, 'end_y':None, 'start_z':0, 'end_z':None, 'fignum':None,
                    'cap_grid':False, 'cap_level':100.0, 'cmap':'bwr', 'saveas':None}

# Define grid plotting function
def plot_grid_plane(grid, axs="z", num=0, start_x=0, end_x=None, start_y=0,
                    end_y=None,start_z=0,end_z=None,fignum=None, 
                    cap_grid=False, cap_level=100.0, cmap='bwr', saveas=None):
    """Plot the specified plane of the given grid. 
    Only planes perpendicular to an axis is currently supported."""

    if cap_grid:
        grid[grid > cap_level] = cap_level
    
    plt.figure(fignum)
    if axs=="x" or axs==1:
        plt.imshow(grid[num,start_y:end_y,start_z:end_z], cmap=cmap)
        plt.ylabel('y')
        plt.xlabel('z')
        plt.title('x = '+str(num))
        plt.colorbar()
    elif axs=="y" or axs==2:
        plt.imshow(grid[start_x:end_x,num,start_z:end_z], cmap=cmap)
        plt.ylabel('x')
        plt.xlabel('z')
        plt.title('y = '+str(num))
        plt.colorbar()
    else:
        plt.imshow(grid[start_x:end_x,start_y:end_y,num], cmap=cmap)
        plt.ylabel('x')
        plt.xlabel('y')
        plt.title('z = '+str(num))
        plt.colorbar()

    if saveas:
        plt.savefig(saveas, format='eps')
    plt.show()

def plot_grid_plane_contour(grid, axs="z", num=0, start_x=0, end_x=None,
                    start_y=0, end_y=None, start_z=0, end_z=None, fignum=None,
                    cap_grid=False, cap_level=100.0, cmap='bwr', saveas=None,
                    contourf=False):
    """Plots a contour plot of the specified plane of the given grid. 
    Only planes perpendicular to an axis is currently supported."""

    if cap_grid:
        grid[grid > cap_level] = cap_level
    
    if end_x:
        x_coord = np.arange(start_x, end_x)
    else:
        x_coord = np.arange(start_x, grid.shape[0])

    if end_y:
        y_coord = np.arange(start_y, end_y)
    else:
        y_coord = np.arange(start_y, grid.shape[1])

    if end_z:
        z_coord = np.arange(start_z, end_z)
    else:
        z_coord = np.arange(start_z, grid.shape[2])


    plt.figure(fignum)
    if axs=="x" or axs==1:
        Y, Z = np.meshgrid(y_coord, z_coord, indexing='ij')
        if contourf:
            contour_plot = plt.contourf(Y, Z, grid[num, start_y:end_y, start_z:end_z], cmap=cmap)
        else:
            contour_plot = plt.contour(Y, Z, grid[num, start_y:end_y, start_z:end_z], cmap=cmap)
        plt.colorbar(contour_plot)
        plt.ylabel('z')
        plt.xlabel('y')
        plt.title('x = '+str(num))

    elif axs=="y" or axs==2:
        X, Z = np.meshgrid(x_coord, z_coord, indexing='ij')
        if contourf:
            contour_plot = plt.contourf(X, Z, grid[start_x:end_x, num, start_z:end_z], cmap=cmap)
        else:
            contour_plot = plt.contour(X, Z, grid[start_x:end_x, num, start_z:end_z], cmap=cmap)
        plt.colorbar(contour_plot)
        plt.ylabel('z')
        plt.xlabel('x')
        plt.title('y = '+str(num))

    else:
        X, Y = np.meshgrid(x_coord, y_coord, indexing='ij')
        if contourf:
            contour_plot = plt.contourf(X, Y, grid[start_x:end_x, start_y:end_y, num], cmap=cmap)
        else:
            contour_plot = plt.contour(X, Y, grid[start_x:end_x, start_y:end_y, num], cmap=cmap)
        plt.colorbar(contour_plot)
        plt.ylabel('x')
        plt.xlabel('y')
        plt.title('z = '+str(num))

    if saveas:
        plt.savefig(saveas, format='eps')
    plt.show()

################
# Main program #
################
def main():
    print(len(sys.argv))
    gridfile = sys.argv[1]
    kwordargs = {}
    str_args = ["axs", "cmap", "saveas"]
    bool_args = ["cap_grid"]
    #int_args = []
    float_args = ["cap_level"]
    for i in range(2,len(sys.argv),2):
        if sys.argv[i] in str_args:
            kwordargs[sys.argv[i]] = str(sys.argv[i+1])
        elif sys.argv[i] in bool_args:
            kwordargs[sys.argv[i]] = bool(sys.argv[i+1])
        elif sys.argv[i] in float_args:
            kwordargs[sys.argv[i]] = float(sys.argv[i+1])
        else:
            kwordargs[sys.argv[i]] = int(sys.argv[i+1])

    print(kwordargs)

    # mode = sys.argv[2]
    with open(gridfile) as f:
        cnt = read_cube(f, read_data=True)

    grid = cnt["data"]
    grid = np.array(grid)
    print(grid.shape)

    #plot_grid_plane(grid, **kwordargs)
    plot_grid_plane_contour(grid, **kwordargs)



if __name__ == "__main__":
    main()
