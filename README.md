# unitcellsampling

Grid-based potential energy sampling of an atom/ion in unit cells of periodic  
structures.

Intended to be used to sample the single particle potential energy as a function of  
the position of a single atom or ion inside a unit cell, on a grid. The single  
particle potential grids are meant to be used in the Ionic-TuTraSt method to generate  
lattice models for diffusion in solids.


### Features:
- Basic, facilitating structure preprocessing
- Utilization of crystallographic symmetry and radial cutoffs based on scaled   
van der Waals radii to reduce the number of calculations
- Energy calculators based on CP2K and LAMMPS through ase, highly flexible



