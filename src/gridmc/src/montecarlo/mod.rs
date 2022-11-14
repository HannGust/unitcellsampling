// use crate::ewald::EwaldSummation;
use crate::ewald_lumol::EwaldSummation;
use crate::io;
use crate::mathutils::*;
use crate::neighbourdisplacer::NeighbourDisplacerAuto;
use crate::restart::RestartMMC;
use crate::walker::Walker;

use float_cmp::approx_eq;
use itertools::izip;
use ndarray::*;
use rand::Rng;

use log::info;
use std::time::Instant;

// use num_cpus::get;
use std::sync::mpsc::channel;
use std::sync::{Arc, Mutex};
use threadpool::ThreadPool;

use physical_constants::BOLTZMANN_CONSTANT_IN_EV_PER_K;

extern crate lru;

use lru::LruCache;

#[derive(Debug)]
pub struct Trajectory {
    /// Walker coordinates
    pub coordinates: Vec<Vec<Vec<i32>>>,
    /// Energies corresponding to the geometry
    pub energy: Vec<f64>,
    /// Step index
    pub step_idx: Vec<u64>,
}

impl Trajectory {
    pub fn new() -> Self {
        Self {
            coordinates: vec![],
            energy: vec![],
            step_idx: vec![],
        }
    }
    pub fn append(&mut self, step_idx: u64, energy: f64, coordinates: Vec<Vec<i32>>) {
        self.step_idx.push(step_idx);
        self.energy.push(energy);
        self.coordinates.push(coordinates);
    }
}

#[derive(Getters, Setters)]
pub struct MonteCarloBuilder {
    /// 3D energy grid
    #[get = "pub"]
    data: Array3<f64>,
    /// Unit cell vectors (3 * 3)
    #[get]
    #[set]
    unitcell: Array2<f64>,

    /// number of walkers (= loading)
    #[get = "pub"]
    #[set]
    n_walkers: u8,
    /// simulation temperature
    #[get]
    #[set]
    temperature: f64,
    /// partial charge of walker atoms
    #[get]
    #[set]
    partial_charge: f64,

    /// Number of Monte Carlo runs
    #[set]
    #[get = "pub"]
    n_runs: u32,
    /// Number of steps per run
    #[set]
    n_steps: u64,

    #[set]
    n_equilibration_steps: u64,

    perturbate_probs: Vec<f64>,

    #[set]
    initial_positions_on_grid: Vec<Vec<i32>>,

    #[set]
    save_trajectory: bool,

    #[set]
    verbosity_level: u8,

    /// Instance of EwaldSummation for calculating Coulomb interaction
    #[set]
    ewald: EwaldSummation,
    /// Displacing method
    #[set]
    nn_displacer: NeighbourDisplacerAuto,

    #[set]
    logfile_name: &'static str,
}

impl MonteCarloBuilder {
    pub fn new(data: Array3<f64>, unitcell: Array2<f64>) -> MonteCarloBuilder {
        let partial_charge = 1f64;
        // let mut ewald = EwaldSummation::new(unitcell.clone(), Array2::zeros((1, 3)));
        // ewald.set_point_charge(partial_charge);
        let ewald = EwaldSummation::new(
            unitcell.clone(),
            Array2::zeros((1, 3)),
            partial_charge,
            "Li",
            /* realspace_cutoff */ 2.0,
            /* kmax */ 5,
        );

        let n_walkers = 1;
        let walker_limits: Vec<i32> = data.shape().iter().map(|&x| x as i32).collect();
        let nn_displacer = NeighbourDisplacerAuto::new(n_walkers as u32, 1, &walker_limits, 500);

        let perturbate_probs = vec![1.0, 0.0, 0.0];

        Self {
            data,
            unitcell,
            n_walkers: 1,
            temperature: 1000f64,
            partial_charge: 1.0,
            n_runs: 10,
            n_steps: 100,
            n_equilibration_steps: 0,
            perturbate_probs,
            initial_positions_on_grid: vec![], // Random positions
            save_trajectory: false,
            verbosity_level: 0, // 0 (minimal logging), 1 (medium logging), 2 (verbose logging)
            ewald,
            nn_displacer,
            logfile_name: "mmc-config.log",
        }
    }

    /// Create a new Metropolis Monte Carlo simulation.
    pub fn build(&self) -> MetropolisMonteCarlo {
        let mut rng = rand::thread_rng();
        let walker_limits: Vec<i32> = self.data.shape().iter().map(|&x| x as i32).collect();

        let walkers: Vec<Walker> = (0..self.n_walkers)
            .map(|_| Walker::new_random(walker_limits.clone(), &mut rng))
            .collect();

        let visited: Array3<u64> = Array3::zeros(self.data.raw_dim());
        let trajectory = Trajectory::new();

        assert_eq!(walkers.len(), self.n_walkers as usize);

        let ewald = EwaldSummation::new(
            self.unitcell.clone(),
            Array2::zeros((1, 3)),
            self.partial_charge,
            "Li",
            /* realspace_cutoff */ 2.0,
            /* kmax */ 5,
        );

        let n_accepted = 0;

        let nn_displacer =
            NeighbourDisplacerAuto::new(self.n_walkers as u32, 1, &walker_limits, 500);

        let logfile_name = "mmc-trajectory.csv";

        // TODO: cache size in builder
        MetropolisMonteCarlo {
            data: self.data.clone(),
            n_walkers: self.n_walkers,
            temperature: self.temperature,
            walkers,
            visited,
            walker_limits,
            ewald,
            save_trajectory: self.save_trajectory,
            verbosity_level: self.verbosity_level,
            logfile_name,
            trajectory,
            perturbate_probs: self.perturbate_probs.clone(),
            nn_displacer,

            n_accepted,
            n_steps: self.n_steps,
            n_equilibration_steps: self.n_equilibration_steps,
            lru_cache: LruCache::<Vec<Vec<i32>>, f64>::new(1000000),
        }

        // MetropolisMonteCarlo::new(
        //     self.data.clone(),
        //     self.unitcell.clone(),
        //     self.partial_charge,
        //     self.n_walkers,
        //     self.temperature,
        //     self.n_steps,
        //     self.n_equilibration_steps,
        //     self.perturbate_probs,
        //     self.save_trajectory,
        //     self.verbosity_level,
        // )
    }

    fn set_perturbation_probs(&mut self, probs: Vec<f64>) {
        assert!(approx_eq!(
            f64,
            probs.iter().fold(0.0, |sum, i| sum + i),
            1.0
        ));
        self.perturbate_probs = probs;
    }
    // // Setter insert setter functions here
    // fn set_n_runs(&mut self, n_runs: u32) {
    //     self.n_runs = n_runs
    // }
}

/// Instance of Metropolis Monte Carlo simulation.
#[derive(Getters)]
pub struct MetropolisMonteCarlo {
    /// 3D energy grid
    data: Array3<f64>,
    /// number of walkers (= loading)
    #[get = "pub"]
    n_walkers: u8,
    /// simulation temperature
    #[get = "pub"]
    temperature: f64,

    /// Save each walker independently here
    #[get = "pub"]
    walkers: Vec<Walker>,
    /// Keep track of already visited grid points
    visited: Array3<u64>,
    #[get = "pub"]
    walker_limits: Vec<i32>,
    /// Instance of EwaldSummation for calculating Coulomb interaction
    #[get = "pub"]
    ewald: EwaldSummation,

    #[get = "pub"]
    save_trajectory: bool,

    #[get = "pub"]
    verbosity_level: u8,

    /// Keep track of trajectory
    #[get = "pub"]
    trajectory: Trajectory,

    #[get = "pub"]
    perturbate_probs: Vec<f64>,

    #[get = "pub"]
    nn_displacer: NeighbourDisplacerAuto,

    #[get = "pub"]
    n_accepted: u64,
    #[get = "pub"]
    n_steps: u64,
    #[get = "pub"]
    n_equilibration_steps: u64,

    #[get = "pub"]
    logfile_name: &'static str,

    lru_cache: LruCache<Vec<Vec<i32>>, f64>,
}

mod logger;

impl MetropolisMonteCarlo {
    /// Create a new Metropolis Monte Carlo simulation.
    ///
    /// # Arguments
    ///
    /// * `data` - 3D energy grid
    /// * `unitcell` - Unitcell vectors
    /// * `point_charge` - Each walker is assigned this charge to calculate Coulomb interaction
    /// * `temperature` - Simulation temperature

    // pub fn get_trajectory(&self) -> (&Vec<u64>, &Vec<f64>, &Vec<Vec<Vec<i32>>>) {
    //     (
    //         &self.trajectory_steps,
    //         &self.trajectory_energy,
    //         &self.trajectory_coords,
    //     )
    // }

    pub fn set_walker_positions(&mut self, grid_coords: Vec<Vec<i32>>) {
        let n_walkers = grid_coords.len();
        assert_eq!(n_walkers, self.n_walkers as usize);
        for idx in 0..n_walkers {
            assert_eq!(grid_coords[idx].len(), 3);
        }
        let grid_coords = grid_coords.clone();
        for (walker, grid_coord) in izip!(&mut self.walkers, grid_coords) {
            walker.set_coords(grid_coord);
        }
    }

    fn perturbate(&mut self) -> Vec<Vec<i32>> {
        let rand: f64 = rand::thread_rng().gen();
        if rand < self.perturbate_probs[0] {
            self.displace_all_to_random_neighboring()
        } else if rand < self.perturbate_probs[0] + self.perturbate_probs[1] {
            self.displace_random_walker_to_random()
        } else {
            self.displace_all_to_random()
        }
    }

    /// Displace a random walker to a random gridpoint
    ///
    /// # Returns
    ///
    /// New list of walkers with respective positions.
    fn displace_random_walker_to_random(&self) -> Vec<Vec<i32>> {
        let mut positions: Vec<Vec<i32>> = self
            .walkers
            .iter()
            .map(|walker| walker.coords.clone())
            .collect();
        let pos: Vec<i32> = (0..3)
            .map(|idx| rand::thread_rng().gen_range(0..self.walkers[0].limits[idx]))
            .collect();
        let rand_walker_idx = rand::thread_rng().gen_range(0..self.walkers.len());
        positions[rand_walker_idx] = pos;

        return positions;
    }

    /// Displace each walker to a random gridpoint
    ///
    /// # Returns
    ///
    /// New list of walkers with respective positions.
    fn displace_all_to_random(&self) -> Vec<Vec<i32>> {
        let mut positions = Vec::with_capacity(self.n_walkers as usize);
        for _ in self.walkers.iter() {
            let pos: Vec<i32> = (0..3)
                .map(|idx| rand::thread_rng().gen_range(0..self.walkers[0].limits[idx]))
                .collect();
            positions.push(pos);
        }
        return positions;
    }

    /// Displace each walker to a random neighboring gridpoint
    ///
    /// # Returns
    ///
    /// New list of walkers with respective positions.
    fn displace_all_to_random_neighboring(&mut self) -> Vec<Vec<i32>> {
        let mut positions = Vec::with_capacity(self.n_walkers as usize);

        let steps = self.nn_displacer.step();
        for (walker, step) in self.walkers.iter().zip(steps.iter()) {
            positions.push(walker.step_by(step.clone()));
        }
        return positions;
    }

    // fn calculate_energy_cost(&mut self, coords: &Vec<Vec<i32>>) -> f64 {
    //     let scaled_positions: Vec<Vec<f64>> = coords
    //         .iter()
    //         .map(|coord| {
    //             coord
    //                 .iter()
    //                 .zip(self.walkers[0].limits.iter())
    //                 .map(|(&c, &limit)| (c as f64 / limit as f64))
    //                 .collect()
    //         })
    //         .collect();
    //     let scaled_positions = vec_to_array(scaled_positions);
    //     return self.ewald.move_cost_total(&scaled_positions);
    // }

    /// Calculate total energy of the walkers on the grid.
    /// The energy is calculated as a sum of two terms:
    ///     
    /// * Sum of energies of each walker on (single-particle) energy grid
    /// * Coulomb interaction calculated via Ewald summation
    ///
    /// # Arguments
    ///
    /// * `coords` - calculate energy depending on the Walker positions
    fn calculate_energy(&mut self, coords: &Vec<Vec<i32>>) -> f64 {
        let mut energy = 0f64;
        // First term: sum of single particle energies
        for coord in coords.iter() {
            energy += self.data[[coord[0] as usize, coord[1] as usize, coord[2] as usize]]
        }
        if coords.len() > 1 {
            let scaled_positions: Vec<Vec<f64>> = coords
                .iter()
                .map(|coord| {
                    coord
                        .iter()
                        .zip(self.walkers[0].limits.iter())
                        .map(|(&c, &limit)| (c as f64 / limit as f64))
                        .collect()
                })
                .collect();
            self.ewald
                .set_scaled_position_vectors(vec_to_array(scaled_positions));

            // Second term: ewald summation
            let ewald_energy = self.ewald.ewald_energy();
            energy += ewald_energy;
            // TODO: check if correct; Probably not necessary to do division every step
        }
        return energy;
    }

    /// Boltzmann probability
    ///
    /// # Arguments
    ///
    /// * `de` - energy difference between two energy states
    fn boltzmann(&self, de: f64) -> f64 {
        if self.temperature < NUMERICAL_TOLERANCE {
            return 0f64;
        }
        return (de / (BOLTZMANN_CONSTANT_IN_EV_PER_K * self.temperature)).exp();
    }

    /// Run Metropolis Monte Carlo simulation.
    ///
    /// # Arguments
    ///
    /// * `n_steps` - Number of Markov permutations
    fn run(&mut self) {
        // Initial Walker coordinates
        let coords: Vec<Vec<i32>> = self.walkers.iter().map(|x| x.coords.clone()).collect();

        // Initial energy
        let mut energy = self.calculate_energy(&coords);
        self.lru_cache.put(coords.clone(), energy);

        // Initialize trajectory
        if self.save_trajectory {
            self.trajectory.append(0, energy, coords);
        }

        let mut energy_new: f64;
        let mut in_equilibration = true;

        // if in_equilibration {
        //     // println!("Starting equilibration...");
        // }
        for step_idx in 0..self.n_steps {
            // Displace Walkers to adjacent grid points
            let pos_new = self.perturbate();
            // let pos_new: Vec<Vec<i32>> = self.walkers.iter().map(|x| x.coords.clone()).collect();
            // Calculate new energy
            let lru_energy = self.lru_cache.get(&pos_new);
            match lru_energy {
                Some(&e) => energy_new = e,
                None => {
                    energy_new = self.calculate_energy(&pos_new);
                    self.lru_cache.put(pos_new.clone(), energy_new);
                }
            }

            let de = energy - energy_new;

            // Accept new state according to Boltzmann
            if de >= 0f64 || self.boltzmann(de) > rand::thread_rng().gen() {
                self.n_accepted += 1;
                // TODO: update only when non-nn-displacement
                self.nn_displacer.n_accepted += 1;
                // Perform MMC step
                energy = energy_new;
                for (walker, p_new) in izip!(self.walkers.iter_mut(), pos_new.iter()) {
                    walker.set_coords(p_new.clone());
                }

                // Save in trajectory
                // TODO: only push when verbose
                if self.save_trajectory {
                    self.trajectory.append(step_idx, energy, pos_new);
                }
            }

            // Keep track of visited points after equilibration
            if in_equilibration && step_idx > self.n_equilibration_steps {
                in_equilibration = false;
                // info!("DONE!");
                // info!("Equlibration stopped after {} steps.", step_idx);
                // info!("Building up the statistics...");
            }
            if !in_equilibration {
                for walker in self.walkers.iter() {
                    self.visited[[
                        walker.coords[0] as usize,
                        walker.coords[1] as usize,
                        walker.coords[2] as usize,
                    ]] += 1;
                }
            }
        }
        // info!("MMC finished after {} steps.", self.n_steps);
    }
}

/// Wrapper function for concurency
fn run_montecarlo(
    single_particle_grid: Array3<f64>,
    unitcell: Array2<f64>,
    temperature: f64,
    point_charge: f64,
    n_walkers: u8,
    mmc_steps: u64,
    equilibrate: u64,
    thread_idx: u32,
    verbosity_level: u8,
    continue_calculation: bool,
    initial_positions_on_grid: Vec<Vec<i32>>,
    save_trajectory: bool,
) -> Array3<u64> {
    let mut mmc_builder = MonteCarloBuilder::new(single_particle_grid.clone(), unitcell.clone());

    // Set variables which are the same for each run
    mmc_builder.set_temperature(temperature);
    mmc_builder.set_partial_charge(point_charge);
    mmc_builder.set_n_walkers(n_walkers);
    mmc_builder.set_n_steps(mmc_steps);
    mmc_builder.set_n_equilibration_steps(equilibrate);
    mmc_builder.set_initial_positions_on_grid(initial_positions_on_grid);
    mmc_builder.set_save_trajectory(verbosity_level > 0);
    mmc_builder.set_verbosity_level(1);

    let mut mc = mmc_builder.build();

    // Set initial walker position on grid if specified
    if mmc_builder.initial_positions_on_grid.len() > 0 {
        mc.set_walker_positions(mmc_builder.initial_positions_on_grid);
    }

    // mc.set_walkers(walkers);
    {
        let mut mmc_logger = io::LoggerMMC::new(verbosity_level, &mc);
        mmc_logger.log_before_run();
    }
    mc.run();
    {
        let mut mmc_logger = io::LoggerMMC::new(verbosity_level, &mc);
        mmc_logger.log_after_run();
    }
    return mc.visited;
}

/// Correct energy grid using a Metropolis Monte Carlo algorithm
///
/// # Arguments
///
/// * `single_particle_grid` - 3D energy grid
/// * `unitcell` - Unitcell vectors
/// * `temperature` - Simulation temperature
/// * `point_charge` - Each walker is assigned this charge to calculate Coulomb interaction
/// * `n_walkers` - Number of interacting walkers during independent simulation run (= loading)
/// * `n_runs` - Number of independent simulation runs
/// * `mmc_steps` - Number of Markov permutations
pub fn correct_energy_grid_mmc(
    single_particle_grid: &Array3<f64>,
    unitcell: &Array2<f64>,
    temperature: f64,
    point_charge: f64,
    n_walkers: u8,
    n_runs: u32,
    mmc_steps: u64,
    equilibrate: u64,
    verbosity_level: u8,
    initial_positions_on_grid: Vec<Vec<i32>>,
) -> Array3<f64> {
    let start = Instant::now();
    // Create empty grid
    let probability_grid: Array3<u64> = Array3::zeros(single_particle_grid.raw_dim());

    let mut mmc_builder = MonteCarloBuilder::new(single_particle_grid.clone(), unitcell.clone());

    // Set variables which are the same for each run
    mmc_builder.set_n_runs(n_runs);
    mmc_builder.set_temperature(temperature);
    mmc_builder.set_partial_charge(point_charge);
    mmc_builder.set_n_walkers(n_walkers);
    mmc_builder.set_n_steps(mmc_steps);
    mmc_builder.set_n_equilibration_steps(equilibrate);
    mmc_builder.set_initial_positions_on_grid(initial_positions_on_grid);
    mmc_builder.set_save_trajectory(verbosity_level > 0);
    mmc_builder.set_verbosity_level(verbosity_level);

    mmc_builder.log();

    // TODO: include Ewald and NN Displacer here, also others that I forgot?

    // Restartability
    let restarter = RestartMMC::new(&mmc_builder);
    let n_remaining_runs = *restarter.n_remaining_runs();

    let restarter = Arc::new(Mutex::new(restarter));

    // Create new ThreadPools. Plattform dependent number of cores
    let n_cpus = num_cpus::get();
    let pool = ThreadPool::new(n_cpus);
    info!("Running on {} threads.", n_cpus);

    info!("Number of total runs {:?}", n_runs);
    info!("Number of remaining runs {:?}", n_remaining_runs);
    // Assign independent Monte Carlo runs to threads
    let (tx, rx) = channel();
    for run_idx in 0..n_remaining_runs {
        let tx = tx.clone();
        let restarter = restarter.clone();
        let mut mmc = mmc_builder.build();

        pool.execute(move || {
            // TODO: include logging and continue and verbosity here.
            info!("Running job {} out of {}", run_idx + 1, n_remaining_runs);
            mmc.run();

            // let visited = run_montecarlo(
            //     single_particle_grid,
            //     unitcell,
            //     temperature,
            //     point_charge,
            //     n_walkers,
            //     mmc_steps,
            //     equilibrate,
            //     run_idx,
            //     verbose,
            //     continue_calculation,
            //     initial_positions_on_grid,
            // );

            let mut restarter = restarter.lock().unwrap();
            restarter.submit(&mmc.visited);
            tx.send(mmc.visited).unwrap();
        });
    }

    // Accumulate all MC runs from this run
    let probability_grid = rx
        .iter()
        .take(n_remaining_runs as usize)
        .fold(probability_grid, |acc, prob_grid| acc + prob_grid);

    // and from the restart file
    let probability_grid = restarter
        .lock()
        .unwrap()
        .parse_restartfile()
        .iter()
        .fold(probability_grid, |acc, prob_grid| acc + prob_grid);

    // Get the number of all sampled points
    let n_sample_points = probability_grid.sum() as f64;

    // Normalize
    let probability_grid: Array3<f64> = probability_grid.map(|&x| (x as f64) / n_sample_points);

    // Calculate the energy according to boltzmann
    let corrected_energy_grid: Array3<f64> = probability_grid.map(|&x| {
        let mut x = x;
        if x < NUMERICAL_TOLERANCE {
            x += NUMERICAL_TOLERANCE;
        }
        return -(x.ln() * BOLTZMANN_CONSTANT_IN_EV_PER_K * temperature);
    });

    // Shift energy to zero
    let minimum_energy = min_3d(&corrected_energy_grid);
    let corrected_energy_grid = corrected_energy_grid.map(|&x| x - minimum_energy);
    let duration = start.elapsed();
    info!("Total time elapsed: {:.2} s", duration.as_secs_f32());

    return corrected_energy_grid;
}

#[cfg(test)]
mod tests {
    use super::*;
    use env_logger;
    use std::convert::TryInto;
    #[test]
    fn test_nn_displacer() {
        let mut nn_displacer = NeighbourDisplacerAuto::new(3, 1, &vec![1, 1, 1], 10);
        for idx in 0..20 {
            println!(
                "{}, {:?}, {}, {}, {}",
                idx,
                nn_displacer.step(),
                nn_displacer.get_distance(),
                nn_displacer.get_n_walkers_to_displace(),
                nn_displacer.get_current_acceptance_ratio()
            );
            // nn_displacer.step();
            nn_displacer.n_accepted += 1;
        }
        assert!(true);
    }

    #[test]
    fn test_mc() {
        env_logger::init();
        let kjmol_to_ev = 0.010364269574711572;
        let path = "../../tests/structures/tdc.cube";
        let cube_file = cubeio::read_cube(path).unwrap();
        let data_shape: Vec<usize> = cube_file.shape.iter().map(|&x| x as usize).collect();
        let data_shape: [usize; 3] = data_shape.try_into().unwrap();
        let data = Array3::from_shape_vec(data_shape, cube_file.data)
            .unwrap()
            .map(|x| x * kjmol_to_ev);
        let cell_vec = cube_file.cell.clone();
        let cell = vec_to_array(cube_file.cell);

        let corrected_grid =
            correct_energy_grid_mmc(&data, &cell, 1000f64, 1f64, 8u8, 4u32, 10000, 0, 1, vec![]);
        let cube_data = cubeio::CubeData {
            atoms: cube_file.atoms,
            data: corrected_grid.into_raw_vec(),
            cell: cell_vec,
            n_atoms: cube_file.n_atoms,
            origin: cube_file.origin,
            shape: cube_file.shape,
        };
        cubeio::write_cube("tdc_loading_1.cube", &cube_data).unwrap();
        assert_eq!(3, 3);
    }
}
