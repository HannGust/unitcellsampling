use std::fs::File;
use std::fs::OpenOptions;
use std::io::prelude::*;
use std::io::{BufReader, Write};
use std::path::Path;

use linecount::count_lines;

use ndarray::Array3;

use crate::montecarlo::MonteCarloBuilder;

use log::info;

#[derive(Clone, Debug, Getters)]
pub struct RestartMMC {
    filename: String,
    counter: u32,
    #[get = "pub"]
    n_remaining_runs: u32,
    #[get = "pub"]
    shape: Vec<usize>,
}

impl RestartMMC {
    pub fn new(builder: &MonteCarloBuilder) -> Self {
        let filename = String::from("mmc-restart.json");

        let mut restart_mmc = Self {
            filename: filename.clone(),
            counter: 0,
            n_remaining_runs: *builder.n_runs(),
            shape: Vec::from(builder.data().shape()),
        };
        if Path::new(&filename).exists() {
            let lines: usize = count_lines(File::open(&restart_mmc.filename).unwrap()).unwrap();
            restart_mmc.counter = lines as u32;
        } else {
            // Create file if it does not exist
            File::create(&restart_mmc.filename).unwrap();
        }
        if *builder.n_runs() >= restart_mmc.counter {
            restart_mmc.n_remaining_runs = *builder.n_runs() - restart_mmc.counter;
        } else {
            restart_mmc.n_remaining_runs = 0;
        }
        info!(
            "Remaining Monte Carlo runs: {:?}",
            restart_mmc.n_remaining_runs
        );
        restart_mmc
    }

    pub fn submit(&mut self, visited: &Array3<u64>) {
        let mut file = OpenOptions::new()
            .write(true)
            .append(true)
            .open(&self.filename)
            .unwrap();
        self.counter += 1;

        for number in visited {
            file.write_all(format!(" {}", number).as_bytes()).unwrap();
        }
        file.write_all(b"\n").unwrap();
    }

    pub fn parse_restartfile(&self) -> Vec<Array3<u64>> {
        let file = File::open(&self.filename).unwrap();
        let reader = BufReader::new(file);

        let mut histrograms: Vec<Array3<u64>> =
            Vec::with_capacity((self.counter + self.n_remaining_runs) as usize);

        let mut index = 0;
        for line in reader.lines() {
            // Skip the lines which where calculated during this run
            index += 1;
            if index > self.counter - self.n_remaining_runs {
                continue;
            }
            let histrogram: Vec<u64> = line
                .unwrap()
                .split(" ")
                .filter_map(|s| s.parse::<u64>().ok())
                .collect::<Vec<_>>();
            let histrogram: Array3<u64> =
                Array3::from_shape_vec((self.shape[0], self.shape[1], self.shape[2]), histrogram)
                    .unwrap();
            histrograms.push(histrogram);
        }
        histrograms
    }
}
