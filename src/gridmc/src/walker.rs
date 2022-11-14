use itertools::izip;
use ndarray::*;
use rand::Rng;

/// Minimum inmage convention
/// Returns the minimum image vector between a pair of positions.
///
/// # Arguments
///
/// * `r1` - position of point 1
/// * `r2` - position of point 2
/// * `cell` - unitcell vectors with shape (3, 3)
pub fn mic(r1: &Array1<f64>, r2: &Array1<f64>, cell: &Array2<f64>) -> Array1<f64> {
    assert_eq!(r1.dim(), 3);
    assert_eq!(r1.dim(), r2.dim());
    assert_eq!(cell.shape(), &[3, 3]);
    let cell_lengths: Vec<f64> = cell.outer_iter().map(|x| (&x * &x).sum().sqrt()).collect();
    let mut vector: Vec<f64> = Vec::new();

    for (comp1, comp2, cl) in izip!(r1, r2, cell_lengths) {
        let mut dist = (comp1 - comp2).abs();
        while dist <= cl / 2. {
            dist = dist - cl;
        }
        vector.push(dist);
    }

    return Array1::from(vector);
}

/// A walker instance representing a point (ion, molecule or atom)
/// in the Metropolis Monte Carlo Simulation.
#[derive(Clone, Debug)]
pub struct Walker {
    /// A walker is specified by a coordinate vector
    /// which specifies the position on the energy grid
    pub coords: Vec<i32>,
    /// and the limits / shape of the grid (e.g. limits of (50, 50, 50)
    /// specify a grid with 125000 grid points in total)
    pub limits: Vec<i32>,
}

impl Walker {
    /// Create a new Walker instance on a random position on the grid.
    ///
    /// # Arguments
    ///
    /// * `limits` - the limits of the grid
    /// * `rng` - reference to random number generator
    pub fn new_random(limits: Vec<i32>, rng: &mut rand::prelude::ThreadRng) -> Self {
        Self {
            coords: vec![
                rng.gen_range(0..limits[0]),
                rng.gen_range(0..limits[1]),
                rng.gen_range(0..limits[2]),
            ],
            limits: limits,
        }
    }

    pub fn new(limits: Vec<i32>, coords: Vec<i32>) -> Self {
        Self { coords, limits }
    }

    /// Set Walker coordinates.
    ///
    /// # Arguments
    ///
    /// * `coords` - new Walker coordinates, e.g. vec![12, 4, 9]
    pub fn set_coords(&mut self, coords: Vec<i32>) {
        self.coords = coords;
    }

    pub fn get_coords(&self) -> &Vec<i32> {
        &self.coords
    }

    /// Displace Walker position by vector.
    ///
    /// # Arguments
    ///
    /// * `grid_displacement` - add this vector to `self.coords`
    pub fn step_by(&self, grid_displacement: Vec<i32>) -> Vec<i32> {
        assert!(grid_displacement < self.limits);
        let mut new_coords: Vec<i32> = vec![
            self.coords[0] + grid_displacement[0],
            self.coords[1] + grid_displacement[1],
            self.coords[2] + grid_displacement[2],
        ];
        for (idx, new_coord) in new_coords.iter_mut().enumerate() {
            if *new_coord < 0 {
                *new_coord += self.limits[idx];
                continue;
            }
            if *new_coord >= self.limits[idx] {
                *new_coord -= self.limits[idx];
            }
        }
        assert!(new_coords.iter().all(|&x| x >= 0));
        assert!(new_coords
            .iter()
            .zip(self.limits.iter())
            .all(|(&x, &limit)| x < limit));
        return new_coords;
    }
}
