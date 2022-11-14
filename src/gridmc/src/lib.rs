#[macro_use]
extern crate getset;

<<<<<<< HEAD
pub mod ewald;
=======
pub mod array2d;
// mod ewald;
mod ewald_lumol;
>>>>>>> lumol_energy
pub mod io;
pub mod mathutils;
pub mod montecarlo;
pub mod neighbourdisplacer;
pub mod parseconfig;
pub mod restart;
pub mod walker;

mod pointsinsphere;

use ewald_lumol::*;
use mathutils::*;
use montecarlo::*;

use ndarray::*;

use pyo3::prelude::*;
use pyo3::wrap_pyfunction;

#[pyfunction]
fn correct_energy_grid_metropolis_monte_carlo(
    single_particle_grid: Vec<f64>,
    grid_shape: Vec<usize>,
    unitcell: Vec<Vec<f64>>,
    temperature: f64,
    point_charge: f64,
    n_walkers: u8,
    n_runs: u32,
    mmc_steps: u64,
    equilibrate: u64,
    verbosity_level: u8,
    initial_grid_positions: Vec<Vec<i32>>,
) -> PyResult<Vec<f64>> {
    let single_particle_grid = Array3::from_shape_vec(
        (grid_shape[0], grid_shape[1], grid_shape[2]),
        single_particle_grid,
    )
    .unwrap();
    Ok(correct_energy_grid_mmc(
        &single_particle_grid,
        &vec_to_array(unitcell),
        temperature,
        point_charge,
        n_walkers,
        n_runs,
        mmc_steps,
        equilibrate,
        verbosity_level,
        initial_grid_positions,
    )
    .into_raw_vec())
}

#[pyfunction]
fn ewald_summation(
    cell: Vec<Vec<f64>>,
    scaled_positions: Vec<Vec<f64>>,
    partial_charge: f64,
) -> PyResult<f64> {
    let mut ewald = EwaldSummation::new(
        vec_to_array(cell.clone()),
        vec_to_array(scaled_positions),
        partial_charge,
        "Li",
        /* realspace_cutoff */ 8.0,
        /* kmax */ 5,
    );
    // ewald.set_point_charge(partial_charge);
    Ok(ewald.ewald_energy())
}

// #[pyfunction]
// fn ewald_summation_self(n_atoms: u32, charge: f64, eta: f64) -> PyResult<f64> {
//     Ok(ewald_self_energy(n_atoms, charge, eta))
// }

// #[pyfunction]
// fn ewald_summation_real(
//     cell: Vec<Vec<f64>>,
//     scaled_positions: Vec<Vec<f64>>,
//     charge: f64,
//     cutoff: f64,
//     sqrt_eta: f64,
// ) -> PyResult<f64> {
//     Ok(ewald_real_energy(
//         &vec_to_array(cell),
//         &vec_to_array(scaled_positions),
//         charge,
//         cutoff,
//         sqrt_eta,
//     ))
// }

// #[pyfunction]
// fn ewald_summation_reciprocal(
//     cell: Vec<Vec<f64>>,
//     scaled_positions: Vec<Vec<f64>>,
//     cutoff: f64,
//     charge: f64,
//     volume: f64,
//     eta: f64,
// ) -> PyResult<f64> {
//     Ok(ewald_reciprocal_energy(
//         &vec_to_array(cell),
//         &vec_to_array(scaled_positions),
//         cutoff,
//         charge,
//         volume,
//         eta,
//     ))
// }

#[pymodule]
fn gridmc(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(ewald_summation, m)?)?;
    // m.add_function(wrap_pyfunction!(ewald_summation_self, m)?)?;
    // m.add_function(wrap_pyfunction!(ewald_summation_real, m)?)?;
    // m.add_function(wrap_pyfunction!(ewald_summation_reciprocal, m)?)?;
    m.add_function(wrap_pyfunction!(
        correct_energy_grid_metropolis_monte_carlo,
        m
    )?)?;
    Ok(())
}
