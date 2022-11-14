================
unitcellsampling
================


.. image:: https://img.shields.io/pypi/v/unitcellsampling.svg
        :target: https://pypi.python.org/pypi/unitcellsampling

.. image:: https://img.shields.io/travis/benjaminbolbrinker/unitcellsampling.svg
        :target: https://travis-ci.com/benjaminbolbrinker/unitcellsampling

.. image:: https://readthedocs.org/projects/unitcellsampling/badge/?version=latest
        :target: https://unitcellsampling.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status



The purpose of this package is to handle 3D energy fields (energy grids) within a crystallographic unit-cell.

* Free software: GNU General Public License v3
* Documentation: https://unitcellsampling.readthedocs.io.


Features
--------

* Calculate energy grids using first-principles or force-fields.
* Supports CP2K and LAMMPS. Easily extensible.
* Create corrected energy grid (incoorporating loading effect) using a Monte Carlo simulation on the energy grids. The algorithm is decriped here (see paper).
