import numpy as np
import matplotlib.pyplot as plt
import re

import ase
import ase.io

from ase.io import read,write
from ase.io.lammpsdata import write_lammps_data
from ase import Atom

import pathlib

import copy

import os

import argparse

### Structure names - should be the same as in ucs_batch_run_reftraj.py
batch_structure = "structure.cif"
batch_structure_unitcell = "structure_unitcell.cif"
###


def init_batch_analysis_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--directory", action="store", type=str, default=".", help="The path to the batch run main working directory, of the batch run to be analyzed.")
    parser.add_argument("-g", "--grid", action="store_true", help="Compile an energy grid to from the individual batch results of the batch run, and write it to a cube file.")
    parser.add_argument("--rm-files", action="store_true", help="Toggles removal of mostly unneccessary or redundant files: *.restart, *.RESTART.wfn, cp2k.out, and initial trajectory files.")
    return parser
#

def args2batch_analysis_options(args):
    """Converts parsed arguments into dictionary of options for batch analysis. 
    """
    opt_dict = {"batch_wd":args.directory,
                "grid":args.grid}
    return opt_dict

def main():
    parser = init_batch_analysis_parser()
    args = parser.parse_args()
    options = args2batch_analysis_options(args)
    print("analyze_batch_reftraj.py: Analyzing batch run with work dir:\t", options["batch_wd"])

    struc_file = os.path.join(options["batch_wd"], batch_structure_unitcell)
    if os.path.exists(struc_file) and os.path.isfile(struc_file):
        print("analyze_batch_reftraj.py: Using unit cell for grid compilation:\t", struc_file)
    else:
        struc_file = os.path.join(options["batch_wd"], batch_structure)
        print("analyze_batch_reftraj.py: Using sampling cell for grid compilation:\t", struc_file)


    items = os.listdir(options["batch_wd"])

    for names in items:
        if names.endswith("symInfo.cube"):
            symmetry_file = os.path.join(options["batch_wd"], names)

    print("analyze_batch_reftraj.py: Symmetry file:\t\t\t", symmetry_file)

    batch_count = 0

    current_batch_dir = os.path.join(options["batch_wd"], "BATCH"+str(batch_count), )

    energies = []
    time_info = []
    while os.path.isdir(current_batch_dir):
        #print(current_batch_dir)
        print("analyze_batch_reftraj.py: Analyzing:\t\t\t\t", current_batch_dir)
        
        energy_file = os.path.join(current_batch_dir, "Batch_reftraj-1.ener")
        energies_i, time_info_i = read_cp2k_energy_file(energy_file)
        energies.extend(energies_i)
        time_info.extend(time_info_i)

        if args.rm_files:
	        # delete unnecessary files
            print("analyze_batch_reftraj.py: rm_files = {} - Removing unneccesary/redundant files.".format(args.rm_files))
            if os.path.isfile(os.path.join(current_batch_dir, "Batch_reftraj-1.restart")):
                os.remove(os.path.join(current_batch_dir, "Batch_reftraj-1.restart"))
            if os.path.isfile(os.path.join(current_batch_dir, "Batch_reftraj-1.restart.bak-1")):
                os.remove(os.path.join(current_batch_dir, "Batch_reftraj-1.restart.bak-1"))
            if os.path.isfile(os.path.join(current_batch_dir, "Batch_reftraj-RESTART.wfn")):
                os.remove(os.path.join(current_batch_dir, "Batch_reftraj-RESTART.wfn"))
            if os.path.isfile(os.path.join(current_batch_dir, "Batch_reftraj-RESTART.wfn.bak-1")):
                os.remove(os.path.join(current_batch_dir, "Batch_reftraj-RESTART.wfn.bak-1"))
            if os.path.isfile(os.path.join(current_batch_dir, "cp2k.out")):
                os.remove(os.path.join(current_batch_dir, "cp2k.out"))
            if os.path.isfile(os.path.join(current_batch_dir, "trajectory.xyz")):
                os.remove(os.path.join(current_batch_dir, "trajectory.xyz"))
#	        if os.path.isfile(os.path.join(current_batch_dir, "Batch_reftraj-pos-1.xyz")):
#               os.remove(os.path.join(current_batch_dir, "Batch_reftraj-pos-1.xyz"))

        batch_count+=1
        current_batch_dir = os.path.join(options["batch_wd"], "BATCH"+str(batch_count), )

    energies = np.array(energies)
    time_info = np.array(time_info)

    energies = energies - np.min(energies)
    # au to eV
    energies = 27.211324570273*energies

    number_of_sampled_energies = np.size(energies)
    print("analyze_batch_reftraj.py: Number of grid points sampled:\t", number_of_sampled_energies)

    # print time information
    total_time = np.sum(time_info)
    print("analyze_batch_reftraj.py: Total simulation time:\t\t", total_time)
    time_per_grid_point = total_time/number_of_sampled_energies
    print("analyze_batch_reftraj.py: Time per grid point:\t\t\t", time_per_grid_point)

    if options["grid"]:
        print("Compiling grid...")

        frame = read(struc_file)

        Na = frame.get_global_number_of_atoms()

        print("analyze_batch_reftraj.py: Number of framework atoms:\t\t",Na)

        atom_list = np.zeros((Na,5))
        atom_list[:,0] = frame.get_atomic_numbers()
        atom_list[:,2:5] = frame.get_positions()

        # get symmmetry information
        symmetry_info, Nx, Ny, Nz, dx, dy, dz  = read_symmetry_info(symmetry_file)

        print("analyze_batch_reftraj.py: Grid dimensions:\t\t\t",(Nx,Ny,Nz))

        cell_info = frame.cell[:][:]
        #dx = cell_info[0][:]/Nx
        #dy = cell_info[1][:]/Ny
        #dz = cell_info[2][:]/Nz

        print("analyze_batch_reftraj.py: Grid spacings:\t\t\t",(dx,dy,dz))

        Np = Nx*Ny*Nz

        energy_grid = np.zeros((Np))

        energy_of_excluded_points = 1e8*np.max(energies)
        for i in range(Np):
            if symmetry_info[i]:
                energy_grid[i] = energies[symmetry_info[i]-1]
            else:
                energy_grid[i] = energy_of_excluded_points
        print("Grid compiled.")
        print("Writing grid to cube-file...")

        grid_file = os.path.join(options["batch_wd"], "energy_grid.cube")
        writeCube(grid_file,energy_grid,Nx,Ny,Nz,dx,dy,dz,Na,atom_list)

        print("Cube file written.")

        
def writeCube(filename,data,Nx,Ny,Nz,dx,dy,dz,Na,a_list,cell_vec_unit="Bohr",atom_coord_unit="Ang"):
    if cell_vec_unit == "Bohr":
        CellUnit2Bohr = 1.0 # Bohr to Bohr
    elif cell_vec_unit == "Ang":
        CellUnit2Bohr = 1.8897161646320724 # Ang to Bohr

    if atom_coord_unit == "Bohr":
        CoordUnit2Bohr = 1.0 # Bohr to Bohr
    elif atom_coord_unit == "Ang":
        CoordUnit2Bohr = 1.8897161646320724 # Ang to Bohr


    f = open(filename,'w')
    f.write('Energy grid, Energy units: eV, length units: Bohr\n')
    f.write('--------------------------------\n')
    f.write(''.join(["{:d}".format(Na),' ','0.000000 0.000000 0.000000\n']))
    f.write(''.join(["{:d}".format(Nx),' ',"{:1.6f}".format(CellUnit2Bohr*dx[0]),' ',"{:1.6f}".format(CellUnit2Bohr*dx[1]),' ',"{:1.6f}".format(CellUnit2Bohr*dx[2]),'\n']))
    f.write(''.join(["{:d}".format(Ny),' ',"{:1.6f}".format(CellUnit2Bohr*dy[0]),' ',"{:1.6f}".format(CellUnit2Bohr*dy[1]),' ',"{:1.6f}".format(CellUnit2Bohr*dy[2]),'\n']))
    f.write(''.join(["{:d}".format(Nz),' ',"{:1.6f}".format(CellUnit2Bohr*dz[0]),' ',"{:1.6f}".format(CellUnit2Bohr*dz[1]),' ',"{:1.6f}".format(CellUnit2Bohr*dz[2]),'\n']))
    for i in range(Na):
        f.write(''.join(["{:d}".format(int(a_list[i][0])),' ',"{:.6f}".format(a_list[i][1]),' ',"{:.6f}".format(CoordUnit2Bohr*a_list[i][2]),' ',"{:.6f}".format(CoordUnit2Bohr*a_list[i][3]),' ',"{:.6f}".format(CoordUnit2Bohr*a_list[i][4]),'\n']))
    
    for val in data:
        f.write(''.join(["{:.6f}".format(val),'\n']))
    f.close()


def read_symmetry_info(filename):
    f = open(filename,"r")

    f.readline()
    f.readline()
    f.readline()

    line = f.readline()
    sp = line.split()
    Nx = int(sp[0])
    dx = np.array([float(sp[1]),float(sp[2]),float(sp[3])])
    line = f.readline()
    sp = line.split()
    Ny = int(sp[0])
    dy = np.array([float(sp[1]),float(sp[2]),float(sp[3])])
    line = f.readline()
    sp = line.split()
    dz = np.array([float(sp[1]),float(sp[2]),float(sp[3])])
    Nz = int(sp[0])

    Np = Nx*Ny*Nz

    data = f.read().splitlines()

    data = np.asarray(data)
    data = data.astype(int)
    
    f.close()

    return data, Nx, Ny, Nz, dx, dy, dz

def read_cp2k_energy_file(filename):
    E = []
    time_info = []

    f = open(filename, "r")

    line = f.readline()

    line = f.readline()
    sp = line.split()
    while line:
        E.append(float(sp[4]))
        time_info.append(float(sp[6]))
        line = f.readline()
        sp = line.split()

    E = np.array(E)

    f.close()
    
    return E, time_info


main()
