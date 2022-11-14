// use log::info;
use rand::seq::SliceRandom;
use rand::{thread_rng, Rng};

#[derive(Copy, Clone, Debug)]
pub struct NeighbourDisplacerAuto {
    n_walkers: u32,

    pub n_accepted: u64,
    pub n_steps: u64,
    distance: u32,
    max_distance: u32,
    n_walkers_to_displace: u32,

    optimal_acceptance_ratio: f64,
    current_acceptance_ratio: f64,
    batch_counter: u32,
    batch_size: u32,
}

impl NeighbourDisplacerAuto {
    pub fn new(n_walkers: u32, distance: u32, limits: &Vec<i32>, batch_size: u32) -> Self {
        assert!(n_walkers > 0);
        assert!(distance > 0);
        Self {
            n_walkers,

            n_accepted: 0,
            n_steps: 0,
            distance,
            max_distance: (*limits.iter().min().unwrap() as u32) / 2,
            n_walkers_to_displace: (n_walkers + 1) / 2,

            optimal_acceptance_ratio: 0.5,
            current_acceptance_ratio: 0.0,
            batch_counter: 0,
            batch_size,
        }
    }

    pub fn get_batchsize(self) -> u32 {
        self.batch_size
    }

    pub fn get_distance(self) -> u32 {
        self.distance
    }

    pub fn get_n_walkers_to_displace(self) -> u32 {
        self.n_walkers_to_displace
    }

    pub fn get_current_acceptance_ratio(self) -> f64 {
        self.current_acceptance_ratio
    }

    pub fn get_batch_size(self) -> u32 {
        self.batch_size
    }

    pub fn get_n_accepted(&self) -> u64 {
        self.n_accepted
    }

    pub fn get_n_steps(&self) -> u64 {
        self.n_steps
    }

    fn update_distance(&mut self) {
        if self.current_acceptance_ratio < self.optimal_acceptance_ratio {
            if self.distance > 1 {
                self.distance -= 1;
            }

        // info!(
        //     "Change nearest-neighbor displacement distance to {}.",
        //     self.distance
        // );
        } else {
            if self.distance < self.max_distance {
                self.distance += 1;
            }
            // info!(
            //     "Change nearest-neighbor displacement distance to {}.",
            //     self.distance
            // );
        }
    }

    fn update_n_walkers_to_displace(&mut self) {
        if self.current_acceptance_ratio < self.optimal_acceptance_ratio {
            if self.n_walkers_to_displace > 1 {
                self.n_walkers_to_displace -= 1;
            }

        // info!(
        //     "Number of walkers to displace {}.",
        //     self.n_walkers_to_displace
        // );
        } else {
            if self.n_walkers_to_displace < self.n_walkers {
                self.n_walkers_to_displace += 1;
            }
            // info!(
            //     "Change nearest-neighbor displacement distance to {}.",
            //     self.n_walkers_to_displace
            // );
        }
    }

    pub fn step(&mut self) -> Vec<Vec<i32>> {
        if self.batch_counter >= self.batch_size {
            self.current_acceptance_ratio = (self.n_accepted as f64) / (self.n_steps as f64);

            // info!("N Accepted: {}", self.n_accepted);
            // info!("N Steps: {}", self.n_steps);
            // info!(
            //     "Acceptance ratio: {:+.5e} percent.",
            //     self.current_acceptance_ratio
            // );

            // Choose a random adjustment
            let r: f64 = thread_rng().gen();
            if r < 0.5 {
                self.update_distance();
            } else {
                self.update_n_walkers_to_displace();
            }

            self.n_accepted = 0;
            self.n_steps = 0;
            self.batch_counter = 0;
        }
        self.n_steps += 1;
        self.batch_counter += 1;

        let mut steps: Vec<Vec<i32>> = Vec::with_capacity(self.n_walkers as usize);

        // Create array containing indecies
        let index_arr: Vec<u32> = (0..self.n_walkers).collect();

        // Choose random elements from this array
        let indicies_to_displace =
            index_arr.choose_multiple(&mut thread_rng(), self.get_n_walkers_to_displace() as usize);

        // Create a boolean array which denotes which walkers to displace
        let mut to_displace_arr: Vec<bool> = vec![false; self.n_walkers as usize];
        for &to_displace_idx in indicies_to_displace {
            to_displace_arr[to_displace_idx as usize] = true;
        }

        for to_displace in to_displace_arr {
            if to_displace {
                let step: Vec<i32> = (0..3)
                    .map(|_| {
                        thread_rng().gen_range(-(self.distance as i32)..=(self.distance as i32))
                    })
                    .collect();
                steps.push(step);
            } else {
                steps.push(vec![0, 0, 0]);
            }
        }
        return steps;
    }
}
