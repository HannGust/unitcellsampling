use ndarray::*;

use physical_constants::{ELECTRIC_CONSTANT, ELECTRON_VOLT};

use lumol::energy::{Ewald, SharedEwald};
use lumol::Vector3D;
use lumol::{Molecule, Particle, System, UnitCell};

use log::info;

fn get_lengths(unitcell: &Vec<Vec<f64>>) -> Vec<f64> {
    let a = unitcell[0].iter().map(|x| x * x).sum::<f64>().sqrt();
    let b = unitcell[1].iter().map(|x| x * x).sum::<f64>().sqrt();
    let c = unitcell[2].iter().map(|x| x * x).sum::<f64>().sqrt();
    vec![a, b, c]
}

fn get_angles(unitcell: &Vec<Vec<f64>>) -> Vec<f64> {
    let lengths = get_lengths(unitcell);
    let mut angles = Vec::new();
    for i in 0..3 {
        let j = if (i as i8 - 1) >= 0 { i - 1 } else { i + 2 };
        let k = if (i as i8 - 2) >= 0 { i - 2 } else { i + 1 };
        let ll = lengths[j] * lengths[k];

        if ll > 1e-16 {
            let x = unitcell[j]
                .iter()
                .zip(unitcell[k].iter())
                .map(|(a, b)| a * b)
                .sum::<f64>()
                / ll;
            angles.push(f64::acos(x).to_degrees());
        } else {
            angles.push(90.0);
        }
    }
    angles
}

/// EwaldSummation instance
pub struct EwaldSummation {
    system: System,
    atom_type: String,
    partial_charge: f64,
    real_space_cutoff: f64,
    kmax: usize,
    alpha: Option<f64>,
}

impl EwaldSummation {
    /// Create a new EwaldSummation instance according to gulp3.4_manual.
    /// Point charge defaults to 1
    /// Accuracy factor defaults to 12
    /// Weight defaults to 1/sqrt(2)
    ///
    /// # Arguments
    ///
    /// * `unitcell_vectors` - Vectors specifying the unit cell, shape (3, 3)
    /// * `scaled_position_vectors` - Scaled position vectors of atoms / ions in unitcell
    pub fn new(
        unitcell_vectors: Array2<f64>,
        scaled_position_vectors: Array2<f64>,
        partial_charge: f64,
        atom_type: &str,
        real_space_cutoff: f64,
        kmax: usize,
    ) -> EwaldSummation {
        let cell = vec![
            vec![
                unitcell_vectors[[0, 0]],
                unitcell_vectors[[0, 1]],
                unitcell_vectors[[0, 2]],
            ],
            vec![
                unitcell_vectors[[1, 0]],
                unitcell_vectors[[1, 1]],
                unitcell_vectors[[1, 2]],
            ],
            vec![
                unitcell_vectors[[2, 0]],
                unitcell_vectors[[2, 1]],
                unitcell_vectors[[2, 2]],
            ],
        ];
        let lengths = get_lengths(&cell);
        let angles = get_angles(&cell);

        let unitcell = UnitCell::triclinic(
            lengths[0], lengths[1], lengths[2], angles[0], angles[1], angles[2],
        );

        let mut system = System::with_cell(unitcell);

        for position_vector in scaled_position_vectors.outer_iter() {
            let position_vector = unitcell.matrix().transposed()
                * Vector3D::new(position_vector[0], position_vector[1], position_vector[2]);
            let mut mol = Molecule::new(Particle::with_position(atom_type, position_vector));
            mol.particles_mut().charge[0] = partial_charge;
            system.add_molecule(mol);
        }
        let ewald = SharedEwald::new(Ewald::new(
            /* cutoff */ real_space_cutoff,
            /* kmax */ kmax,
            /* alpha */ None,
        ));
        system.set_coulomb_potential(Box::new(ewald));
        return Self {
            system,
            atom_type: String::from(atom_type),
            partial_charge,
            real_space_cutoff,
            kmax,
            alpha: None,
        };
    }

    fn create_new_system(&self, scaled_position_vectors: Array2<f64>) -> System {
        let unitcell = self.system.cell;
        let mut system = System::with_cell(unitcell);

        for position_vector in scaled_position_vectors.outer_iter() {
            let position_vector = unitcell.matrix().transposed()
                * Vector3D::new(position_vector[0], position_vector[1], position_vector[2]);
            let mut mol = Molecule::new(Particle::with_position(
                self.atom_type.clone(),
                position_vector,
            ));
            mol.particles_mut().charge[0] = self.partial_charge;
            system.add_molecule(mol);
        }
        let ewald = SharedEwald::new(Ewald::new(
            /* cutoff */ self.real_space_cutoff,
            /* kmax */ self.kmax,
            /* alpha */ self.alpha,
        ));
        system.set_coulomb_potential(Box::new(ewald));
        return system;
    }

    pub fn set_scaled_position_vectors(&mut self, scaled_position_vectors: Array2<f64>) {
        let mut n_particles = 0;
        for mol in self.system.molecules() {
            for part in mol.particles() {
                n_particles += 1;
            }
        }
        if scaled_position_vectors.shape()[0] != n_particles {
            self.system = self.create_new_system(scaled_position_vectors.clone());
            return;
        }

        for (idx, position_vector) in scaled_position_vectors.outer_iter().enumerate() {
            let position_vector = self.system.cell.matrix().transposed()
                * Vector3D::new(position_vector[0], position_vector[1], position_vector[2]);
            self.system.molecule_mut(idx).particles_mut().position[0] = position_vector;
        }
    }
    pub fn ewald_energy(&self) -> f64 {
        let e = self.system.potential_energy();
        e
    }
}
