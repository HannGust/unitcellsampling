use serde::Deserialize;
use serde_json::{Result, Value};

use std::error::Error;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;

struct MonteCarloConfig {
    /// Cube file containing 3D energy grid and Unit cell vectors (3 * 3)
    unitcell: String,

    /// number of walkers (= loading)
    n_walkers: u8,
    /// simulation temperature
    temperature: f64,
    /// partial charge of walker atoms
    partial_charge: f64,

    /// Instance of EwaldSummation for calculating Coulomb interaction
    ewald: String,
    /// Displacing method
    nn_displacer: String,
}

fn read_user_from_file<P: AsRef<Path>>(path: P) -> Result<Value> {
    // Open the file in read-only mode with buffer.
    let file = File::open(path).unwrap();
    let reader = BufReader::new(file);

    // Read the JSON contents of the file as an instance of `User`.
    let u: Value = serde_json::from_reader(reader)?;
    Ok(u)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_config() {
        let u = read_user_from_file("mmc-config-defaults.json").unwrap();
        println!("{:#?}", u);
    }
}
