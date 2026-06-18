import copy
import numpy as np
import re

import ase.io
from pymatgen.core.lattice import Lattice

from mendeleev import element


def parse_ddec_file(filename='DDEC6_even_tempered_net_atomic_charges.xyz'):
    ''' Parser for DDEC6 atomic charge output file
    '''
    with open(filename) as fp:
        n_atoms = int(fp.readline().strip())
        fp.readline()
        regex = r'(\w+)\s+[+-]?(?:[0-9]*[.])?[0-9]+\s+[+-]?(?:[0-9]*[.])?[0-9]+\s+[+-]?(?:[0-9]*[.])?[0-9]+\s+([+-]?(?:[0-9]*[.])?[0-9]+)'
        element_charge_pair = re.search(
            regex,
            fp.readline())
        elements = []
        charge = []
        for _ in range(n_atoms):
            elements.append(element_charge_pair.group(1))
            charge.append(element_charge_pair.group(2))
            element_charge_pair = re.search(
                regex,
                fp.readline())
        return elements, charge


def parse_repeat_file(filename='cp2k-RESP_CHARGES.resp'):
    ''' Parser for *.resp file originating from
    '''
    with open(filename) as fp:
        assert(fp.readline().strip() == 'RESP charges:')
        assert(fp.readline().strip() == 'Type |   Atom   |    Charge')
        fp.readline()
        regex = r'\sRESP\s+\d+\s+(\w+)\s+([+-]?([0-9]*[.])?[0-9]+)'
        element_charge_pair = re.search(
            regex,
            fp.readline())
        elements = []
        charge = []
        while element_charge_pair:
            elements.append(element_charge_pair.group(1))
            charge.append(element_charge_pair.group(2))
            element_charge_pair = re.search(
                regex,
                fp.readline())
        return elements, charge


def fill_with_ions(atoms, data, type='Li'):
    ''' Fill up potential local minima of PES with ions while considering vdW radius.
          Parameters
        ----------
        atoms
            Atoms object containing the atoms and unitcell.

        data
            (n, m, k) shaped list containing the PES.

        Returns
        -------
        ase:Atoms
            The updated atoms object.

    '''

    # times 2 and devide by 100 for convert units (pm to Å)
    r_Li = element(type).vdw_radius * 0.02

    min_indices = []
    cart_coords = []
    frac_coords = []
    temp_data = copy.deepcopy(data)

    min_index = np.array(np.unravel_index(
        np.argmin(temp_data, axis=None), temp_data.shape))
    frac_coord = min_index / np.array(temp_data.shape)
    cart_coord = np.matmul(frac_coord, atoms.get_cell()[:])
    lat = Lattice(atoms.get_cell()[:])

    min_indices.append(min_index)
    cart_coords.append(cart_coord)
    frac_coords.append(frac_coord)
    temp_data[min_index] = data.max()

    # TODO: adjust for universal case: fill up for all grid points
    for i in range(500):
        min_index = np.array(np.unravel_index(
            np.argmin(temp_data, axis=None), temp_data.shape))
        frac_coord = min_index / np.array(temp_data.shape)
        cart_coord = np.matmul(frac_coord, atoms.get_cell()[:])
        lat = Lattice(atoms.get_cell()[:])
        fcoord, dist, _, _ = lat.get_points_in_sphere(
            frac_coords, cart_coord, r_Li, zip_results=False)

        if len(dist) == 0:
            min_indices.append(min_index)
            cart_coords.append(cart_coord)
            frac_coords.append(frac_coord)
        temp_data[min_index[0]][min_index[1]][min_index[2]] = data.max()

    for cart_coord in cart_coords:
        atoms.append(ase.Atom('Li', np.array(cart_coord)))

    return atoms
