use gridmc::mathutils::vec_to_array;
use ndarray::*;
use std::convert::TryInto;

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use gridmc::montecarlo::correct_energy_grid_mmc;

pub fn criterion_benchmark(c: &mut Criterion) {
    let mut group = c.benchmark_group("MMC");
    // Configure Criterion.rs to detect smaller differences and increase sample size to improve
    // precision and counteract the resulting noise.
    group.sample_size(10);
    let path = "/home/benjamin/Code/unitcellsampling/tests/structures/tdc.cube";
    let cube_file = cubeio::read_cube(path).unwrap();
    let data_shape: Vec<usize> = cube_file.shape.iter().map(|&x| x as usize).collect();
    let data_shape: [usize; 3] = data_shape.try_into().unwrap();
    let data = Array3::from_shape_vec(data_shape, cube_file.data).unwrap();
    let cell = vec_to_array(cube_file.cell);
    group.bench_function("Correct energy with grid Monte Carlo", |b| {
        b.iter(|| {
            // correct_energy_grid_mmc(
            //     black_box(&data),
            //     black_box(&cell),
            //     black_box(100f64),
            //     black_box(1f64),
            //     black_box(2u8),
            //     black_box(8u32),
            //     black_box(500),
            // );
            correct_energy_grid_mmc(
                black_box(&data),
                black_box(&cell),
                black_box(1000f64),
                black_box(1f64),
                black_box(8u8),
                black_box(8u32),
                black_box(500),
                black_box(0),
                black_box(0),
                vec![],
            );
        })
    });
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
