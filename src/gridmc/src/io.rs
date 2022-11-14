use crate::montecarlo;

use std::fs::File;
use std::fs::{remove_file, OpenOptions};
use std::io::prelude::*;
use std::path::Path;

pub struct LoggerMMC<'a> {
    verbosity_level: u8,
    log_filename: String,
    montecarlo: &'a montecarlo::MetropolisMonteCarlo,
    file: File,
}

impl<'a> LoggerMMC<'a> {
    pub fn new(verbosity_level: u8, montecarlo: &'a montecarlo::MetropolisMonteCarlo) -> Self {
        let mut outfile_idx = 0;
        let mut log_filename = format!("mmc_{}.log", outfile_idx);

        while Path::new(&format!("mmc_{}.log", outfile_idx)).exists() {
            log_filename = format!("mmc_{}.log", outfile_idx);
            outfile_idx += 1;
        }

        let file = OpenOptions::new()
            .write(true)
            .create(true)
            .open(log_filename.clone())
            .unwrap();

        // if verbosity_level > 0 {
        //     if !continue_calculation {
        //         let path = Path::new(&log_self.filename);
        //         if path.exists() {
        //             remove_self.file(&log_filename).unwrap();
        //         }
        //     } else {
        //         let mut continue_idx = 0;
        //         let mut new_self.file = Path::new(&log_filename);
        //         while new_self.file.exists() {
        //             log_self.filename = format!("mmc_{}_continue_{}.log", thread_idx, continue_idx);
        //             continue_idx = continue_idx + 1;
        //             new_self.file = Path::new(&log_filename);
        //         }
        //     }
        // }
        Self {
            verbosity_level,
            log_filename,
            montecarlo,
            file,
        }
    }

    pub fn log_before_run(&mut self) {
        if self.verbosity_level == 0 {
            return;
        }
        self.file.write_all(b"Runing MMC calculating:\n").unwrap();
        self.file.write_all(b"\n").unwrap();

        let perturbation = format!("Perturbation probabilities:\n");
        self.file.write_all(perturbation.as_bytes()).unwrap();
        let perturbation = format!(
            "Displace to all walkers to adjacent grid points: {}\n",
            self.montecarlo.perturbate_probs()[0]
        );
        self.file.write_all(perturbation.as_bytes()).unwrap();
        let perturbation = format!(
            "Displace to a random walker to an adjacent grid point: {}\n",
            self.montecarlo.perturbate_probs()[1]
        );
        self.file.write_all(perturbation.as_bytes()).unwrap();
        let perturbation = format!(
            "Displace to all walkers to random grid points: {}\n",
            self.montecarlo.perturbate_probs()[2]
        );
        self.file.write_all(perturbation.as_bytes()).unwrap();
        self.file.write_all(b"\n").unwrap();

        self.file
            .write_all(
                format!(
                    "Temperature: {}\n",
                    self.montecarlo.temperature().to_string()
                )
                .as_bytes(),
            )
            .unwrap();

        let loading = format!(
            "Loading / number of walkers: {}\n",
            self.montecarlo.n_walkers().to_string()
        );
        self.file.write_all(loading.as_bytes()).unwrap();

        let loading = format!("Grid size: {:?}\n", self.montecarlo.walker_limits());
        self.file.write_all(loading.as_bytes()).unwrap();

        let nn_displacement = format!(
            "The inital grid displacement is {}\n",
            self.montecarlo.nn_displacer().get_distance()
        );
        self.file.write_all(nn_displacement.as_bytes()).unwrap();
        self.file.write_all(b"\n").unwrap();

        let nn_displacement = format!(
            "The initial number of walkers to displace is {} out of {}\n",
            self.montecarlo.nn_displacer().get_n_walkers_to_displace(),
            self.montecarlo.n_walkers()
        );
        self.file.write_all(nn_displacement.as_bytes()).unwrap();
        self.file.write_all(b"\n").unwrap();

        let nn_displacement = format!(
            "The batch size for adjusting the acceptance ratio is {} steps.\n",
            self.montecarlo.nn_displacer().get_batch_size()
        );
        self.file.write_all(nn_displacement.as_bytes()).unwrap();
        self.file.write_all(b"\n").unwrap();

        self.file.write_all(b"\n").unwrap();
        self.file
            .write_all(b"Initial walker positions on grid:\n")
            .unwrap();
        let mut walker_str = String::new();
        for walker in self.montecarlo.walkers().iter() {
            let coords = walker.get_coords();
            let coord_str = format!("{} {} {}\n", coords[0], coords[1], coords[2]);
            walker_str.push_str(&coord_str);
        }
        self.file.write_all(walker_str.as_bytes()).unwrap();

        self.file.write_all(b"\n").unwrap();

        // self.file.write_all(b"Ewald summation:\n").unwrap();
        // self.file.write_all(b"\n").unwrap();
        // let ewald = self.montecarlo.ewald();

        // let ewald_str = format!("Point charge: {}\n", ewald.get_point_charge().to_string());
        // self.file.write_all(ewald_str.as_bytes()).unwrap();

        // let ewald_str = format!("accuracy factor: {}\n", ewald.get_accuracy_factor());
        // self.file.write_all(ewald_str.as_bytes()).unwrap();
        // let ewald_str = format!("eta: {}\n", ewald.get_eta().to_string());
        // self.file.write_all(ewald_str.as_bytes()).unwrap();

        // let ewald_str = format!("Realspace cutoff: {}\n", ewald.get_real_cutoff());
        // self.file.write_all(ewald_str.as_bytes()).unwrap();
        // let ewald_str = format!("Reciprocal cutoff: {}\n", ewald.get_reciprocal_cutoff());
        // self.file.write_all(ewald_str.as_bytes()).unwrap();

        self.file.write_all(b"\n").unwrap();
        self.file.write_all(b"Running MMC...").unwrap();
    }

    pub fn log_after_run(&mut self) {
        if self.verbosity_level == 0 {
            return;
        }
        let mut file = OpenOptions::new()
            .append(true)
            .open(self.log_filename.clone())
            .unwrap();

        self.file.write_all(b" DONE!\n").unwrap();
        self.file.write_all(b"\n").unwrap();

        let steps_str = format!("Performed {} MMC steps.\n", self.montecarlo.n_steps());
        self.file.write_all(steps_str.as_bytes()).unwrap();
        let steps_str = format!(
            "The first {} steps are excluded from generating a statistics.\n",
            self.montecarlo.n_equilibration_steps()
        );
        self.file.write_all(steps_str.as_bytes()).unwrap();

        let acceptance_rate = format!(
            "The acceptance rate was {:+.5e}.\n",
            (*self.montecarlo.n_accepted() as f64) / (*self.montecarlo.n_steps() as f64)
        );
        self.file.write_all(acceptance_rate.as_bytes()).unwrap();
        self.file.write_all(b"\n").unwrap();

        let nn_displacement = format!(
            "The final grid displacement is {}\n",
            self.montecarlo.nn_displacer().get_distance()
        );
        self.file.write_all(nn_displacement.as_bytes()).unwrap();
        self.file.write_all(b"\n").unwrap();

        let nn_displacement = format!(
            "The final number of walkers to displace is {} out of {}\n",
            self.montecarlo.nn_displacer().get_n_walkers_to_displace(),
            self.montecarlo.n_walkers()
        );
        self.file.write_all(nn_displacement.as_bytes()).unwrap();
        self.file.write_all(b"\n").unwrap();

        // TODO: switch delimiter to , without space: more efficient pandas parsing
        self.file.write_all(b"Printing trajectory:\n").unwrap();
        self.file.write_all(b"\n").unwrap();
        self.file.write_all(b"step, energy, ").unwrap();
        for i in 0..*self.montecarlo.n_walkers() {
            self.file
                .write_all(format!("grid coordinate walker_{}, ", i.to_string()).as_bytes())
                .unwrap();
        }
        self.file.write_all(b"\n").unwrap();
        // TODO: print trajectory in seperate self.file
        let traj = self.montecarlo.trajectory();
        for (step, (energy, coords)) in traj
            .step_idx
            .iter()
            .zip(traj.energy.iter().zip(traj.coordinates.iter()))
        {
            let step_energy_str = format!("{}, {:+.3e}, ", step, energy);
            self.file.write_all(step_energy_str.as_bytes()).unwrap();

            let mut walker_str = String::new();
            for coord in coords.iter() {
                let coord_str = format!("{} {} {}, ", coord[0], coord[1], coord[2]);
                walker_str.push_str(&coord_str);
            }
            self.file.write_all(walker_str.as_bytes()).unwrap();
            self.file.write_all(b"\n").unwrap();
        }

        self.file.write_all(b"End of trajectory:\n").unwrap();
        self.file.write_all(b"\n").unwrap();
    }
}
