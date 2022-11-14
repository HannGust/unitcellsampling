use super::MetropolisMonteCarlo;
use super::MonteCarloBuilder;

use std::fs::OpenOptions;
use std::io::prelude::*;

impl MonteCarloBuilder {
    pub fn log(&self) {
        if self.verbosity_level == 0 {
            return;
        }

        println!(
            "###################################################################################"
        );
        println!("  Metropolis Monte Carlo Sampling on a Discrete Potential Energy Surface");
        println!(
            "###################################################################################"
        );
        println!();
        println!("                       Benjamin Bolbrinker");
        println!("             Uppsala University - Ångström Laboratory");
        println!("                Department for Structural Chemistry");
        println!();
        println!("      This routine corrects the a single particle potential energy surface");
        println!(
            "of ions moving in a rigid periodic framework using metropolis monte carlo sampling."
        );
        println!("           You can find details in <include citation here>.");

        println!();
        println!("              For full docs see <insert link here>.");
        println!();
        println!(
            "###################################################################################"
        );
        println!();

        println!("The following parameters have been specified for the simulation:");
        println!();

        println!("temperature : {} K", self.temperature().to_string());

        println!("loading : {}", self.n_walkers().to_string());

        println!(
            "grid_displacement_initial : {}",
            self.nn_displacer.get_distance()
        );

        println!(
            "n_displace_initial : {} (out of {})",
            self.nn_displacer.get_n_walkers_to_displace(),
            self.n_walkers()
        );

        println!(
            "acceptance_ratio_batchsize : {}",
            self.nn_displacer.get_batch_size()
        );

        println!();
        println!("Perturbation probabilities:");
        println!();
        println!(
            "all_to_adjacent : {}", // "Displace to all walkers to adjacent grid points : {}",
            self.perturbate_probs[0]
        );
        println!(
            "random_to_adjacent : {}", // "Displace to a random walker to an adjacent grid point: {}",
            self.perturbate_probs[1]
        );
        println!(
            "all_to_random : {}", // "Displace to all walkers to random grid points: {}",
            self.perturbate_probs[2]
        );

        println!("Ewald summation parameters:");
        println!();

        // println!(
        //     "partial_charge : {} C",
        //     self.ewald.get_point_charge().to_string()
        // );

        // println!("accuracy_factor : {}", self.ewald.get_accuracy_factor());
        // println!("eta : {}", self.ewald.get_eta().to_string());

        // // TODO: insert units here!
        // println!("real_cutoff : {}", self.ewald.get_real_cutoff());
        // println!("reciprocal_cutoff : {}", self.ewald.get_reciprocal_cutoff());
    }
}

impl MetropolisMonteCarlo {
    fn log(&self) {
        if self.verbosity_level == 0 {
            return;
        }
        let mut file = OpenOptions::new()
            .create(true)
            .open(self.logfile_name)
            .unwrap();

        let loading = format!("Grid size: {:?}\n", self.walker_limits);
        file.write_all(loading.as_bytes()).unwrap();
    }
}
