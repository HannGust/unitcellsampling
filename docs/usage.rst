=====
Usage
=====

Command-line interface
----------------------

To use unitcellsampling using the CLI:



Run Metropolis Monte Carlo on energy grid
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To correct a single particle energy grid using the Metropolis Monte Carlo algorithm on the energy grid run

.. code:: bash

    grid-correct-mmc

In order to create a ASE trajectory run

.. code:: bash

    grid-correct-mmc-create-traj


Calculate partial charges
^^^^^^^^^^^^^^^^^^^^^^^^^

For getting accurate results for incoorporating loading effects it might be desired to get partial charges of Lithium atoms within a framework.
To calculate these run

.. code:: bash

    grid-get-partial-charge       

Adjust grid size
^^^^^^^^^^^^^^^^

Grids can be up- and down sampled using cubic or linear interpolation using 

.. code:: bash

    grid-downsample

Average energy grid using Spacegroup
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For increasing the accuracy of the resulting tunnelling barriers averaging over symmetrically identical gridpoints might be desired:

.. code:: bash

    grid-average-by-symmetry    

Within a project
----------------

To use unitcellsampling in a project

.. code:: python

    import unitcellsampling


Create an energy grid
^^^^^^^^^^^^^^^^^^^^^

To create an energy grid you first have to create a instance of the :code:`UnitCellSampler()`.
You have to provide an :code:`ase.Atoms` object where the unitcell is defined. 

.. code:: python

    sampler = unitcellsampling.sample.UnitCellSampler(lgps)

To create you can then create a 16x16x24 grid using 

.. code:: python

    sampler.generate_grid_vectors(n_frac=(16, 16, 24), vdw_scale=0)

Often we are not interested in grid points located directly to atoms.
To ignore grid points within a constant times the respective VdW radius you can set :code:`vdw_scale`.

.. code:: python

    sampler.generate_grid_vectors(n_frac=(16, 16, 24), vdw_scale=1.0)

Create a function which calculates the energy at each grid point and call the :code:`calculate_energies` method on the :code:`UnitCellSampler()` object.

.. code:: python

    energies = sampler.calculate_energies(
        method=dft_pbe, atom="Li", exploit_symmetry=True)

A full example is provided in the `examples/` directory:

.. literalinclude:: ../examples/generate_dft_energy_grid_optimized/grid_gen.py
   :language: python
