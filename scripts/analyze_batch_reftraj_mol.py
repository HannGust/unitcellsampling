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
os.chdir(pathlib.Path(__file__).parent.resolve())

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
    print("analyze_batch_reftraj.py: Analyzing batch run with work dir:\t", options["batch_wd"])

    struc_file = options["batch_wd"]+"structure.cif"

    items = os.listdir(options["batch_wd"])

    for names in items:
        if names.endswith("symInfo.cube"):
            symmetry_file = options["batch_wd"]+names

    print("analyze_batch_reftraj.py: Symmetry file:\t\t\t", symmetry_file)

    batch_count = 0

    current_batch_dir = options["batch_wd"]+"BATCH"+str(batch_count)+"/"

    energies = []
    time_info = []
    while os.path.isdir(current_batch_dir):
        #print(current_batch_dir)
        print("analyze_batch_reftraj.py: Analyzing:\t\t\t\t", current_batch_dir)
        
        energy_file = current_batch_dir+"Batch_reftraj-1.ener"
        energies_i, time_info_i = read_cp2k_energy_file(energy_file)
        energies.extend(energies_i)
        time_info.extend(time_info_i)
	
	    # delete unnecessary files
        if os.path.isfile(current_batch_dir+"Batch_reftraj-1.restart"):
            os.remove(current_batch_dir+"Batch_reftraj-1.restart")
        if os.path.isfile(current_batch_dir+"Batch_reftraj-1.restart.bak-1"):
            os.remove(current_batch_dir+"Batch_reftraj-1.restart.bak-1")
        if os.path.isfile(current_batch_dir+"Batch_reftraj-RESTART.wfn"):
            os.remove(current_batch_dir+"Batch_reftraj-RESTART.wfn")
        if os.path.isfile(current_batch_dir+"Batch_reftraj-RESTART.wfn.bak-1"):
            os.remove(current_batch_dir+"Batch_reftraj-RESTART.wfn.bak-1")
        if os.path.isfile(current_batch_dir+"cp2k.out"):
            os.remove(current_batch_dir+"cp2k.out")
        if os.path.isfile(current_batch_dir+"trajectory.xyz"):
            os.remove(current_batch_dir+"trajectory.xyz")
        if os.path.isfile(current_batch_dir+"Batch_reftraj-pos-1.xyz"):
            os.remove(current_batch_dir+"Batch_reftraj-pos-1.xyz")

        batch_count+=1
        current_batch_dir = options["batch_wd"]+"BATCH"+str(batch_count)+"/"

    energies = np.array(energies)
    time_info = np.array(time_info)

    energies = energies - np.min(energies)
    # au to eV
    energies = 27.211324570273*energies

    number_of_sampled_energies = np.size(energies)
    print("analyze_batch_reftraj.py: Number of energies sampled:\t", number_of_sampled_energies)

    # print time information
    total_time = np.sum(time_info)
    print("analyze_batch_reftraj.py: Total simulation time:\t\t", total_time)
    time_per_grid_point = total_time/number_of_sampled_energies
    print("analyze_batch_reftraj.py: Time per sp calculation:\t\t\t", time_per_grid_point)

    if options["grid"]:
        print("Compiling grid...")

        frame = read(struc_file)

        Na = frame.get_global_number_of_atoms()

        print("analyze_batch_reftraj.py: Number of framework atoms:\t\t",Na)

        atom_list = np.zeros((Na,5))
        atom_list[:,0] = frame.get_atomic_numbers()
        atom_list[:,2:5] = frame.get_positions()

        # get symmmetry information
        symmetry_info, Nx, Ny, Nz = read_symmetry_info(symmetry_file)

        print("analyze_batch_reftraj.py: Grid dimensions:\t\t\t",(Nx,Ny,Nz))

        number_of_unique_grid_points = np.max(symmetry_info)
        rotations_per_grid_points = int(number_of_sampled_energies/number_of_unique_grid_points)
        print("analyze_batch_reftraj.py: Number of unique grid points:\t", number_of_unique_grid_points)
        print("analyze_batch_reftraj.py: Number of rotations per point:\t", rotations_per_grid_points)

        Temp = 300
        kBT = 1e-5*8.617333262*Temp

        energies_avg = np.zeros((number_of_unique_grid_points))
        for i in range(number_of_unique_grid_points):
            sum1=0
            sum2=0
            E_0 = np.min(energies[(i*rotations_per_grid_points):((i+1)*rotations_per_grid_points)])
            for j in range(rotations_per_grid_points):
                k = i*rotations_per_grid_points+j
                w = np.exp(-(energies[k]-E_0)/kBT)
                #if (E_lmp[k]-E_0)/kBT < 100:
                sum1 = sum1 + w*energies[k]
                sum2 = sum2 + w

            energies_avg[i] = sum1/sum2

        energies_avg= energies_avg - np.min(energies_avg)

        cell_info = frame.cell[:][:]
        dx = cell_info[0][:]/Nx
        dy = cell_info[1][:]/Ny
        dz = cell_info[2][:]/Nz

        print("analyze_batch_reftraj.py: Grid spacings:\t\t\t",(dx,dy,dz))

        Np = Nx*Ny*Nz

        energy_grid = np.zeros((Np))

        energy_of_excluded_points = 1e8*np.max(energies_avg)
        for i in range(Np):
            if symmetry_info[i]:
                energy_grid[i] = energies_avg[symmetry_info[i]-1]
            else:
                energy_grid[i] = energy_of_excluded_points
        print("Grid compiled.")
        print("Writing grid to cube-file...")

        grid_file = options["batch_wd"]+"energy_grid.cube"
        writeCube(grid_file,energy_grid,Nx,Ny,Nz,dx,dy,dz,Na,atom_list)

        print("Cube file written.")

        
def writeCube(filename,data,Nx,Ny,Nz,dx,dy,dz,Na,a_list):
    AngtoBohr = 1.8897161646320724
    f = open(filename,'w')
    f.write('Energy grid obtained via a classical potential in lammps, Energy units: eV, length units: Bohr\n')
    f.write('--------------------------------\n')
    f.write(''.join(["{:d}".format(Na),' ','0.000000 0.000000 0.000000\n']))
    f.write(''.join(["{:d}".format(Nx),' ',"{:1.6f}".format(AngtoBohr*dx[0]),' ',"{:1.6f}".format(AngtoBohr*dx[1]),' ',"{:1.6f}".format(AngtoBohr*dx[2]),'\n']))
    f.write(''.join(["{:d}".format(Ny),' ',"{:1.6f}".format(AngtoBohr*dy[0]),' ',"{:1.6f}".format(AngtoBohr*dy[1]),' ',"{:1.6f}".format(AngtoBohr*dy[2]),'\n']))
    f.write(''.join(["{:d}".format(Nz),' ',"{:1.6f}".format(AngtoBohr*dz[0]),' ',"{:1.6f}".format(AngtoBohr*dz[1]),' ',"{:1.6f}".format(AngtoBohr*dz[2]),'\n']))
    for i in range(Na):
        f.write(''.join(["{:d}".format(int(a_list[i][0])),' ',"{:.6f}".format(a_list[i][1]),' ',"{:.6f}".format(AngtoBohr*a_list[i][2]),' ',"{:.6f}".format(AngtoBohr*a_list[i][3]),' ',"{:.6f}".format(AngtoBohr*a_list[i][4]),'\n']))
    
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
    line = f.readline()
    sp = line.split()
    Ny = int(sp[0])
    line = f.readline()
    sp = line.split()
    Nz = int(sp[0])

    Np = Nx*Ny*Nz

    data = f.read().splitlines()

    data = np.asarray(data)
    data = data.astype(int)

    return data, Nx, Ny, Nz

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
    return E, time_info


main()
