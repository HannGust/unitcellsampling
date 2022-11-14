use crate::array2d::Array2D;
use crate::mathutils::*;
use crate::pointsinsphere::*;
use float_cmp::approx_eq;

use itertools::izip;
use ndarray::*;
use statrs::function::erf::erfc;

use physical_constants::{ELECTRIC_CONSTANT, ELECTRON_VOLT};
use std::f64::consts::PI;

// Converts unit of q*q/r into eV
static CONV_FACT: f64 = 1e10 * ELECTRON_VOLT / (4f64 * PI * ELECTRIC_CONSTANT);

fn to_ndarray(array: Array2D) -> Array2<f64> {
    return Array::from_shape_vec((array.n_rows, array.n_cols), array.data).unwrap();
}

/// EwaldSummation instance
#[derive(Debug)]
pub struct EwaldSummation {
    /// cartesian unitcell vectors in Å
    unitcell_vectors: Array2D,
    /// reciprocal unndarray = "0.13"itcell vectors in Å^-1
    reciprocal_cell: Array2D,
    /// (unitless) scaled / fractional position vectors
    scaled_position_vectors: Array2D,
    /// Each atom is assigned this charge
    point_charge: f64,
    /// Eta see gulp3.4_manual
    eta: f64,
    /// accuracy see gulp3.4_manual
    accuracy_factor: f64,

    /// real space energy cutoff
    real_cutoff: f64,
    /// reciprocal space energy cutoff
    reciprocal_cutoff: f64,

    /// sqrt of `self.eta`
    sqrt_eta: f64,
    /// see gulp3.4_manual
    accf: f64,
    /// total charge = sum of all point charges
    total_charge: f64,
    /// unitcell volume Å^3
    volume: f64,
    /// number of atoms in unitcell
    n_atoms: u32,

    /// previous result for real part
    previous_real_space_energy: f64,
    /// previous result for reciprocal part
    previous_reciprocal_space_energy: f64,

    /// real-space supercell size
    super_cell_repeats_real: Vec<Array1<f64>>,
    /// reciprocal-space supercell size
    super_cell_repeats_reci: Vec<Array1<f64>>,
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
    ) -> EwaldSummation {
        // first unitcell vector
        let a: ArrayView1<f64> = unitcell_vectors.slice(s![0, ..]);
        // second unitcell vector
        let b: ArrayView1<f64> = unitcell_vectors.slice(s![1, ..]);
        // third unitcell vector
        let c: ArrayView1<f64> = unitcell_vectors.slice(s![2, ..]);
        // Volume of unitcell
        let volume = a.dot(&cross_product(&b, &c)).abs();

        let n_atoms: u32 = scaled_position_vectors.len_of(Axis(0)) as u32;
        let point_charge = 1f64;
        let total_charge = point_charge * (n_atoms as f64);

        let acc_factor = 12f64;

        let accf = (acc_factor * 10f64.ln()).sqrt();
        // weight
        let w = 1. / (2.0f64.sqrt());
        let eta = ((n_atoms as f64) * w / (volume * volume)).powf(1. / 3.) * PI;
        let sqrt_eta = eta.sqrt();

        let real_cutoff = accf / sqrt_eta;

        let reciprocal_cutoff = 2. * sqrt_eta * accf;

        let reciprocal_cell: Array2<f64> = invert3(&unitcell_vectors); // TODO implement inversion
        let reciprocal_cell = reciprocal_cell.t();
        let reciprocal_cell = reciprocal_cell.map(|&x| x * 2. * PI);

        // Supercell size in units of unitcell vectors for calculating
        // realspace and
        let super_cell_repeats_real = get_super_cell(&unitcell_vectors, real_cutoff);
        // reciprocal term
        let super_cell_repeats_reci = get_super_cell(&reciprocal_cell, reciprocal_cutoff);
        //TODO: Update when unitcell or r_cut changes

        EwaldSummation {
            unitcell_vectors: Array2D::from_ndarray(unitcell_vectors),
            reciprocal_cell: Array2D::from_ndarray(reciprocal_cell),
            scaled_position_vectors: Array2D::from_ndarray(scaled_position_vectors),
            point_charge: point_charge,
            eta: eta,
            accuracy_factor: acc_factor,

            real_cutoff,
            reciprocal_cutoff,

            sqrt_eta: sqrt_eta,
            accf: accf,
            total_charge: total_charge,
            volume: volume,
            n_atoms: n_atoms,

            previous_real_space_energy: 0f64,
            previous_reciprocal_space_energy: 0f64,

            super_cell_repeats_real,
            super_cell_repeats_reci,
        }
    }

    pub fn get_eta(&self) -> f64 {
        self.eta
    }

    pub fn get_accuracy_factor(&self) -> f64 {
        self.accuracy_factor
    }

    pub fn get_point_charge(&self) -> f64 {
        self.point_charge
    }

    pub fn get_real_cutoff(&self) -> f64 {
        self.real_cutoff
    }

    pub fn get_reciprocal_cutoff(&self) -> f64 {
        self.reciprocal_cutoff
    }

    pub fn get_unitcell_vectors(&self) -> &Array2D {
        &self.unitcell_vectors
    }

    pub fn set_scaled_position_vectors(&mut self, vectors: Array2<f64>) {
        self.n_atoms = vectors.len_of(Axis(0)) as u32;
        self.scaled_position_vectors = Array2D::from_ndarray(vectors);
        self.total_charge = self.point_charge * (self.n_atoms as f64);

        let w = 1. / (2.0f64.sqrt());
        self.eta = ((self.n_atoms as f64) * w / (self.volume * self.volume)).powf(1. / 3.) * PI;
        self.sqrt_eta = self.eta.sqrt();
        self.real_cutoff = self.accf / self.sqrt_eta;
        self.reciprocal_cutoff = 2. * self.sqrt_eta * self.accf;
    }

    pub fn guess_cutoffs_from_eta_and_accuracy_factor(&mut self, eta: f64, accuracy_factor: f64) {
        self.eta = eta;
        self.sqrt_eta = eta.sqrt();
        self.accuracy_factor = accuracy_factor;
        self.accf = (self.accuracy_factor * 10f64.ln()).sqrt();

        self.real_cutoff = self.accf / self.sqrt_eta;
        self.reciprocal_cutoff = 2. * self.sqrt_eta * self.accf;

        self.super_cell_repeats_real =
            get_super_cell(&to_ndarray(self.unitcell_vectors.clone()), self.real_cutoff);
        self.super_cell_repeats_reci = get_super_cell(
            &to_ndarray(self.reciprocal_cell.clone()),
            self.reciprocal_cutoff,
        );
    }

    pub fn guess_cutoffs_from_eta(&mut self, eta: f64) {
        self.eta = eta;
        self.sqrt_eta = eta.sqrt();

        self.real_cutoff = self.accf / self.sqrt_eta;
        self.reciprocal_cutoff = 2. * self.sqrt_eta * self.accf;

        self.super_cell_repeats_real =
            get_super_cell(&to_ndarray(self.unitcell_vectors.clone()), self.real_cutoff);
        self.super_cell_repeats_reci = get_super_cell(
            &to_ndarray(self.reciprocal_cell.clone()),
            self.reciprocal_cutoff,
        );
    }

    pub fn guess_cutoffs_from_accuracy_factor(&mut self, accuracy_factor: f64) {
        self.accuracy_factor = accuracy_factor;
        self.accf = (self.accuracy_factor * 10f64.ln()).sqrt();

        self.real_cutoff = self.accf / self.sqrt_eta;
        self.reciprocal_cutoff = 2. * self.sqrt_eta * self.accf;

        self.super_cell_repeats_real =
            get_super_cell(&to_ndarray(self.unitcell_vectors.clone()), self.real_cutoff);
        self.super_cell_repeats_reci = get_super_cell(
            &to_ndarray(self.reciprocal_cell.clone()),
            self.reciprocal_cutoff,
        );
    }

    pub fn set_point_charge(&mut self, charge: f64) {
        self.point_charge = charge;
        self.total_charge = charge * (self.n_atoms as f64);
    }

    pub fn set_real_cutoff(&mut self, cutoff: f64) {
        self.real_cutoff = cutoff;

        self.super_cell_repeats_real =
            get_super_cell(&to_ndarray(self.unitcell_vectors.clone()), self.real_cutoff);
    }

    pub fn set_reciprocal_cutoff(&mut self, cutoff: f64) {
        self.reciprocal_cutoff = cutoff;
        self.super_cell_repeats_reci = get_super_cell(
            &to_ndarray(self.reciprocal_cell.clone()),
            self.reciprocal_cutoff,
        );
    }

    /// Calculate real space term of Ewald Summation
    pub fn ewald_real_energy(&self) -> f64 {
        // Array::from(vec![1.0; 1000])
        //     .map(|x| erfc(x * self.sqrt_eta) * self.point_charge * self.point_charge / x)
        //     .sum()
        //     * CONV_FACT
        //     * 0.5
        let distances = distances_within(
            self.real_cutoff,
            &self.unitcell_vectors,
            &self.scaled_position_vectors,
            &self.super_cell_repeats_real,
        );
        distances
            .iter()
            .map(|x| erfc(x * self.sqrt_eta) * self.point_charge * self.point_charge / x)
            .sum::<f64>()
            * CONV_FACT
            * 0.5
    }

    /// Calculate reciprocal space term of Ewald Summation
    pub fn ewald_reciprocal_energy(&self) -> f64 {
        let prefactor = 2. * PI / self.volume;

        let cartesian_positions = self.scaled_position_vectors.dot(&self.unitcell_vectors);

        let recip_nn = coordinates_within(
            self.reciprocal_cutoff,
            &self.reciprocal_cell,
            // &Array::zeros((1, 3)),
            &Array2D::new(vec![0f64; 3], 1, 3), // TODO: check if shape correct
            &self.super_cell_repeats_reci,
        );
        let g2s: Vec<f64> = recip_nn
            .rows_iter()
            .map(|row| row.map(|&x| x * x).sum())
            .collect(); //.map(|x| ((&x * &x).sum())).collect();
        let exp_values: Vec<f64> = g2s.iter().map(|x| (-x / (4. * self.eta)).exp()).collect();

        let grs = recip_nn.dot(&cartesian_positions.t());
        let grs = Array::from_shape_vec((grs.n_rows, grs.n_cols), grs.data).unwrap();
        let charge_squared = self.point_charge * self.point_charge;

        let mut reciprocal_energy: f64 = 0.;
        for (gr, g2, exp_value) in izip!(grs.outer_iter(), g2s.iter(), exp_values.iter()) {
            let m: Array2<f64> = vec_to_array(
                gr.iter()
                    .map(|&x| {
                        gr.iter()
                            .map(|&y| (x - y + PI / 4.).sin() * (exp_value / g2))
                            .collect()
                    })
                    .collect(),
            );
            reciprocal_energy += m.sum();
        }
        reciprocal_energy *= prefactor * CONV_FACT * charge_squared * 2f64.sqrt();
        reciprocal_energy
    }

    /// Calculate self energy of Ewald Summation
    pub fn ewald_self_energy(&self) -> f64 {
        -((self.n_atoms as f64) * self.point_charge * self.point_charge * (self.eta / PI).sqrt())
            * CONV_FACT
    }

    /// Calculate charge correction of Ewald Summation
    pub fn ewald_charge_correction(&self) -> f64 {
        -CONV_FACT / 2. * PI / self.volume / self.eta * self.total_charge * self.total_charge
    }

    /// Calculate total Ewald Summation energy
    pub fn ewald_energy(&self) -> f64 {
        let energy_real = self.ewald_real_energy();

        let energy_reci = self.ewald_reciprocal_energy();

        let energy_self = self.ewald_self_energy();

        let energy_charge_correction: f64 = self.ewald_charge_correction();

        energy_real + energy_reci + energy_self + energy_charge_correction
    }

    pub fn ewald_energy_reduced(&self) -> f64 {
        let energy_real = self.ewald_real_energy();

        // let energy_reci = self.ewald_reciprocal_energy();

        energy_real //+ energy_reci
                    // 0.0
    }

    // TODO: replace f64 by i32 on grid to make more robust and efficient
    pub fn _compare_scaled_vectors(
        &self,
        vectors_1: &Array2<f64>,
        vectors_2: &Array2<f64>,
    ) -> Vec<bool> {
        let vectors_1 = vectors_1.map(|&x| {
            if x > 0.0 {
                return x % 1.;
            }
            if x < 0.0 {
                return (x % 1.) + 1.;
            }
            return 0.0;
        });

        let vectors_2 = vectors_2.map(|&x| {
            if x > 0.0 {
                return x % 1.;
            }
            if x < 0.0 {
                return (x % 1.) + 1.;
            }
            return 0.0;
        });
        let mut bool_vector: Vec<bool> = Vec::with_capacity(vectors_1.dim().0);
        for (vec_1, vec_2) in vectors_1.outer_iter().zip(vectors_2.outer_iter()) {
            bool_vector.push(vec_1.all_close(&vec_2, 1e-6));
        }
        return bool_vector;
    }

    // pub fn move_cost_real_space(&self, scaled_vectors_to: &Array2<f64>) -> f64 {
    //     assert_eq!(self.scaled_position_vectors.len(), scaled_vectors_to.len());
    //     /// Order important!
    //     let moved_walkers =
    //         self._compare_scaled_vectors(scaled_vectors_to, &self.scaled_position_vectors);

    //     let distances = distances_within(
    //         self.real_cutoff,
    //         &self.unitcell_vectors,
    //         &self.scaled_position_vectors,
    //         &self.super_cell_repeats_real,
    //     );

    //     distances
    //         .map(|x| erfc(x * self.sqrt_eta) * self.point_charge * self.point_charge / x)
    //         .sum()
    //         * CONV_FACT
    //         * 0.5
    // }

    // pub fn move_cost_reciprocal_space(&self, scaled_vectors_to: &Array2<f64>) -> f64 {
    //     assert_eq!(self.scaled_position_vectors.len(), scaled_vectors_to.len());
    //     0.0
    // }

    // pub fn move_cost_total(&self, scaled_vectors_to: &Array2<f64>) -> f64 {
    //     self.move_cost_real_space(scaled_vectors_to)
    //         + self.move_cost_reciprocal_space(scaled_vectors_to)
    // }
}

#[cfg(test)]
mod tests {
    extern crate ndarray;
    use super::*;

    use float_cmp::approx_eq;
    #[test]
    fn move_cost() {
        let cell = array![[10., 0., 0.], [0., 10., 0.], [0., 0., 10.]];
        let scaled_positions_old = array![[0., 0., 1.], [0., 0.4, 0.1]];
        let scaled_positions_new = array![[0., 0., 0.], [0.2, 0., 1.1]];

        println!("scaled positions {}", scaled_positions_old);
        let mut ewald = EwaldSummation::new(cell, scaled_positions_old.clone());
        println!(
            "Bool {:?}",
            ewald._compare_scaled_vectors(&scaled_positions_old, &scaled_positions_new)
        );
        let energy_old = ewald.ewald_energy();
        println!("Ewald old {}", energy_old);
        ewald.set_scaled_position_vectors(scaled_positions_new.clone());
        let energy_new = ewald.ewald_energy();
        println!("Ewald new {}", energy_new);
        ewald.set_scaled_position_vectors(scaled_positions_old);
        // let move_cost = ewald.move_cost_total(&scaled_positions_new);
        // println!("Ewald move cost {}", move_cost);
        // assert!(approx_eq!(f64, energy_new - energy_old, move_cost))
    }

    #[test]
    fn test_positions_within() {
        let cell = array![[10., 0., 0.], [0., 10., 0.0], [0., 0., 10.]];
        let scaled_positions = array![[0., 0., 0.], [0., 0., 0.1]];
        let dist = points_within(14.5, &cell, &scaled_positions);
        let mut dist = dist.0.to_vec();
        dist.sort_by(|a, b| a.partial_cmp(b).unwrap());
        println!("Distances {:?}", dist);
    }

    #[test]
    fn test_ewald_energy_cubic() {
        let cell = array![[10., 0., 0.], [0., 10., 0.], [0., 0., 10.]];
        let scaled_positions = array![[0., 0., 0.], [0., 0., 0.1]];
        let mut ewald = EwaldSummation::new(cell, scaled_positions);
        ewald.guess_cutoffs_from_eta(2.);
        assert!(approx_eq!(
            f64,
            ewald.ewald_self_energy(),
            -22.978509616541267,
            epsilon = 0.00001
        ));
        println!("REAL EN {}", ewald.ewald_real_energy());
        ewald.guess_cutoffs_from_eta(4.);
        ewald.set_real_cutoff(20.);
        ewald.set_reciprocal_cutoff(2.);
        assert!(approx_eq!(
            f64,
            ewald.ewald_real_energy(),
            0.06735772536911086,
            epsilon = 0.00001
        ));
        assert!(approx_eq!(
            f64,
            ewald.ewald_reciprocal_energy(),
            23.153001445996434,
            epsilon = 0.00001
        ));
    }

    #[test]
    fn test_ewald_energy_noncubic() {
        let cell = array![[10., 0., 9.], [0., 10., 9.], [10., 10., 0.]];
        let scaled_positions = array![
            [0., 0., 0.],
            [0.0, 0.0, 0.1],
            [0.5, 0.6, 0.15],
            [0.9, 0.3, 0.85]
        ];

        let mut ewald = EwaldSummation::new(cell, scaled_positions);
        ewald.guess_cutoffs_from_eta(2.);
        assert!(approx_eq!(
            f64,
            ewald.ewald_self_energy(),
            -45.957019233082534,
            epsilon = 0.00001
        ));

        ewald.guess_cutoffs_from_eta(4.);
        ewald.set_real_cutoff(20.);
        ewald.set_reciprocal_cutoff(0.6);
        assert!(approx_eq!(
            f64,
            ewald.ewald_real_energy(),
            0.0006449586772356167,
            epsilon = 0.00001
        ));
        assert!(approx_eq!(
            f64,
            ewald.ewald_reciprocal_energy(),
            7.716058193645823,
            epsilon = 0.00001
        ));
    }
}

// #[cfg(test)]
// mod tests {
//     extern crate ndarray;
//     use super::*;
//     use float_cmp::approx_eq;

//     #[test]
//     fn test_positions_within() {
//         let cell = array![[10., 0., 0.], [0., 10., 0.0], [0., 0., 10.]];
//         let scaled_positions = array![[0., 0., 0.], [0., 0., 0.1]];
//         let (dist, _) = distances_within(14.5, &cell, &scaled_positions);
//         let mut dist = dist.to_vec();
//         dist.sort_by(|a, b| a.partial_cmp(b).unwrap());
//         println!("Distances {:?}", dist);
//     }

//     #[test]
//     fn test_ewald_energy_cubic() {
//         let cell = array![[10., 0., 0.], [0., 10., 0.], [0., 0., 10.]];
//         let scaled_positions = array![[0., 0., 0.], [0., 0., 0.1]];

//         let a: ArrayView1<f64> = cell.slice(s![0, ..]);
//         let b: ArrayView1<f64> = cell.slice(s![1, ..]);
//         let c: ArrayView1<f64> = cell.slice(s![2, ..]);
//         let volume = a.dot(&cross_product(&b, &c));

//         println!(
//             "Vol: {} Rec_E: {}",
//             volume,
//             ewald_reciprocal_energy(&cell, &scaled_positions, 2., 1., volume, 4.)
//         );
//         assert!(approx_eq!(
//             f64,
//             ewald_self_energy(scaled_positions.nrows() as u32, 1., 2.),
//             -22.978509616541267,
//             epsilon = 0.00001
//         ));
//         assert!(approx_eq!(
//             f64,
//             ewald_real_energy(&cell, &scaled_positions, 1., 20., 2.),
//             0.06735772536911086,
//             epsilon = 0.00001
//         ));
//         assert!(approx_eq!(
//             f64,
//             ewald_reciprocal_energy(&cell, &scaled_positions, 2., 1., volume, 4.),
//             23.153001445996434,
//             epsilon = 0.00001
//         ));
//     }

//     #[test]
//     fn test_ewald_energy_noncubic() {
//         let cell = array![[10., 0., 9.], [0., 10., 9.], [10., 10., 0.]];
//         let scaled_positions = array![
//             [0., 0., 0.],
//             [0.0, 0.0, 0.1],
//             [0.5, 0.6, 0.15],
//             [0.9, 0.3, 0.85]
//         ];

//         let a: ArrayView1<f64> = cell.slice(s![0, ..]);
//         let b: ArrayView1<f64> = cell.slice(s![1, ..]);
//         let c: ArrayView1<f64> = cell.slice(s![2, ..]);
//         let volume = a.dot(&cross_product(&b, &c)).abs();

//         println!(
//             "Vol: {} Rec_E: {}",
//             volume,
//             ewald_reciprocal_energy(&cell, &scaled_positions, 0.6, 1., volume, 4.)
//         );
//         assert!(approx_eq!(
//             f64,
//             ewald_self_energy(scaled_positions.nrows() as u32, 1., 2.),
//             -45.957019233082534,
//             epsilon = 0.00001
//         ));

//         println!(
//             "REAL: {}",
//             ewald_real_energy(&cell, &scaled_positions, 1., 20., 2.)
//         );

//         assert!(approx_eq!(
//             f64,
//             ewald_real_energy(&cell, &scaled_positions, 1., 20., 2.),
//             0.0006449586772356167,
//             epsilon = 0.00001
//         ));
//         assert!(approx_eq!(
//             f64,
//             ewald_reciprocal_energy(&cell, &scaled_positions, 0.6, 1., volume, 4.),
//             7.716058193645823,
//             epsilon = 0.00001
//         ));
//     }
// }
