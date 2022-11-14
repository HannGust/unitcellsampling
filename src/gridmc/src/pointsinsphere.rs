use crate::array2d::Array2D;
use crate::mathutils::*;
use ndarray::*;
use std::f64::consts::PI;

static PI2: f64 = 2. * PI;

pub fn get_super_cell(cell: &Array2<f64>, r_cut: f64) -> Vec<Array1<f64>> {
    let reciprocal_cell: Array2<f64> = invert3(cell); // TODO implement inversion
    let reciprocal_cell = reciprocal_cell.t();
    let reciprocal_cell = reciprocal_cell.map(|&x| x * PI2);

    let reciprocal_lengths: Array1<f64> = reciprocal_cell
        .outer_iter()
        .map(|vec| (&vec * &vec).sum().sqrt())
        .collect();

    let max_repeats: Array1<i32> = reciprocal_lengths // TODO: check why
        .map(|len| r_cut * len / (PI2))
        .map(|x| x.ceil() as i32);

    let mut repeats: Vec<Array1<f64>> =
        Vec::with_capacity((8 * max_repeats[0] * max_repeats[1] * max_repeats[2] + 8) as usize);
    for a in -max_repeats[0]..=max_repeats[0] {
        for b in -max_repeats[1]..=max_repeats[1] {
            for c in -max_repeats[2]..=max_repeats[2] {
                //TODO: Filter out vectors which are bigger r_cut + max(posible distance in unitcell)
                repeats.push(array![a as f64, b as f64, c as f64]);
            }
        }
    }
    repeats
}

// Version for general use
pub fn points_within(
    r_cut: f64,
    cell: &Array2<f64>,
    scaled_positions: &Array2<f64>,
) -> (Array1<f64>, Array2<f64>) {
    get_points_in_spheres(&cell, &scaled_positions, &scaled_positions.dot(cell), r_cut)
}

// Version for general use
pub fn get_points_in_spheres(
    cell: &Array2<f64>,
    fractional_coordinates: &Array2<f64>,
    center_coordinates: &Array2<f64>,
    r_cut: f64,
) -> (Array1<f64>, Array2<f64>) {
    let fractional_coordinates = fractional_coordinates.map(|&x| {
        if x > 0.0 {
            return x % 1.;
        }
        if x < 0.0 {
            return (x % 1.) + 1.;
        }
        return 0.0;
    });

    let r_cut_squared = r_cut * r_cut;
    let super_cell_repeats = get_super_cell(cell, r_cut);
    let mut super_cell_fractional_vectors: Vec<Array1<f64>> =
        Vec::with_capacity(super_cell_repeats.len() * fractional_coordinates.nrows());
    for repeat in super_cell_repeats {
        for frac_coord in fractional_coordinates.outer_iter() {
            super_cell_fractional_vectors.push(&frac_coord + &repeat);
        }
    }

    let mut distances: Vec<f64> = Vec::new();
    let mut coordinates: Vec<Vec<f64>> = Vec::new();
    for center_coord in center_coordinates.outer_iter() {
        for frac_coord in super_cell_fractional_vectors.iter() {
            let valid_coord = frac_coord.dot(cell);
            let distance_squared: f64 = valid_coord
                .iter()
                .zip(center_coord.iter())
                .map(|(valid, center)| (valid - center).powi(2))
                .sum();
            if distance_squared < r_cut_squared + NUMERICAL_TOLERANCE
                && distance_squared > NUMERICAL_TOLERANCE
            {
                coordinates.push(valid_coord.to_vec());
                distances.push(distance_squared.sqrt());
            }
        }
    }
    (Array::from(distances), vec_to_array(coordinates))
}

// TODO include Verlet List and Cell list to improve scaling
// Optimized version for Metropolis Monte Carlo
pub fn get_coords_in_spheres(
    cell: &Array2D,
    fractional_coordinates: &Array2D,
    center_coordinates: &Array2D,
    super_cell_repeats: &Vec<Array1<f64>>,
    r_cut: f64,
) -> Array2D {
    let fractional_coordinates = fractional_coordinates.clone().map(|&x| {
        if x > 0.0 {
            return x % 1.;
        }
        if x < 0.0 {
            return (x % 1.) + 1.;
        }
        return 0.0;
    });

    let r_cut_squared = r_cut * r_cut;

    let super_cell_repeat_size = super_cell_repeats.len() * fractional_coordinates.n_rows;
    let mut super_cell_fractional_vectors: Array2D = Array2D::new(
        vec![0f64; super_cell_repeat_size * 3],
        super_cell_repeat_size,
        3,
    );

    let mut idx = 0;
<<<<<<< HEAD
    for frac_coord in fractional_coordinates.outer_iter() {
        for repeat in super_cell_repeats {
            super_cell_fractional_vectors[[idx, 0]] = repeat[0] + frac_coord[0];
            super_cell_fractional_vectors[[idx, 1]] = repeat[1] + frac_coord[1];
            super_cell_fractional_vectors[[idx, 2]] = repeat[2] + frac_coord[2];
=======
    // TODO: don't use Array but directly Vectors: should be much faster https://github.com/rust-ndarray/ndarray/issues/571
    for frac_coord in fractional_coordinates.rows_iter() {
        let frac: Vec<f64> = frac_coord.map(|&x| x).collect();
        for repeat in super_cell_repeats {
            super_cell_fractional_vectors
                .get_mut(idx, 0)
                .map(|x| *x = repeat[0] + frac[0]);
            super_cell_fractional_vectors
                .get_mut(idx, 1)
                .map(|x| *x = repeat[1] + frac[1]);
            super_cell_fractional_vectors
                .get_mut(idx, 2)
                .map(|x| *x = repeat[2] + frac[2]);
            // super_cell_fractional_vectors[[idx, 1]] = repeat[1] + frac[1];
            // super_cell_fractional_vectors[[idx, 2]] = repeat[2] + frac[2];
>>>>>>> lumol_energy

            // TODO: check if correct
            // let vec: Array1<f64> = repeat + &frac_coord;
            // super_cell_fractional_vectors[[idx, 0]] = vec[0];
            // super_cell_fractional_vectors[[idx, 1]] = vec[1];
            // super_cell_fractional_vectors[[idx, 2]] = vec[2];
            idx += 1;
        }
    }

    // let mut distances: Vec<f64> = Vec::with_capacity(fractional_coordinates.n_rows * 27);
    let mut coordinates: Vec<Vec<f64>> = Vec::new();
    // println!("{:?}", center_coordinates);
    // println!("{:?}", super_cell_fractional_vectors);

    // TODO: also do not do this with Array!
    let super_cell_cartesian_vectors = super_cell_fractional_vectors.dot(&cell);
    // TODO: check if good!
    for center_coord in center_coordinates.rows_iter() {
        let center_coord: Vec<f64> = center_coord.map(|&x| x).collect();
        for valid_coord in super_cell_cartesian_vectors.rows_iter() {
            let valid_coord: Vec<f64> = valid_coord.map(|&x| x).collect();

            let distance_squared = (valid_coord[0] - center_coord[0]).powi(2)
                + (valid_coord[1] - center_coord[1]).powi(2)
                + (valid_coord[2] - center_coord[2]).powi(2);
            if distance_squared < r_cut_squared + NUMERICAL_TOLERANCE
                && distance_squared > NUMERICAL_TOLERANCE
            {
                coordinates.push(valid_coord.to_vec());
            }
        }
    }
    Array2D::from_rows(&coordinates)
}

// fn to_ndarray(array: Array2D) -> Array2<f64> {
//     return Array::from_shape_vec((array.n_rows, array.n_cols), array.data).unwrap();
// }

// TODO include Verlet List and Cell list to improve scaling
// Optimized version for Metropolis Monte Carlo
pub fn get_distances_in_spheres(
    cell: &Array2D,
    fractional_coordinates: &Array2D,
    center_coordinates: &Array2D,
    super_cell_repeats: &Vec<Array1<f64>>,
    r_cut: f64,
) -> Vec<f64> {
    // let cell = to_ndarray(cell.clone());
    // let fractional_coordinates = to_ndarray(fractional_coordinates.clone());
    // let center_coordinates = to_ndarray(center_coordinates.clone());

    let fractional_coordinates = fractional_coordinates.clone().map(|&x| {
        if x > 0.0 {
            return x % 1.;
        }
        if x < 0.0 {
            return (x % 1.) + 1.;
        }
        return 0.0;
    });

    let r_cut_squared = r_cut * r_cut;

    let super_cell_repeat_size = super_cell_repeats.len() * fractional_coordinates.n_rows;
    let mut super_cell_fractional_vectors: Array2D = Array2D::new(
        vec![0f64; super_cell_repeat_size * 3],
        super_cell_repeat_size,
        3,
    );

    let mut idx = 0;
    // TODO: don't use Array but directly Vectors: should be much faster https://github.com/rust-ndarray/ndarray/issues/571
    for frac_coord in fractional_coordinates.rows_iter() {
        let frac: Vec<f64> = frac_coord.map(|&x| x).collect();
        for repeat in super_cell_repeats {
            super_cell_fractional_vectors
                .get_mut(idx, 0)
                .map(|x| *x = repeat[0] + frac[0]);
            super_cell_fractional_vectors
                .get_mut(idx, 1)
                .map(|x| *x = repeat[1] + frac[1]);
            super_cell_fractional_vectors
                .get_mut(idx, 2)
                .map(|x| *x = repeat[2] + frac[2]);
            // super_cell_fractional_vectors[[idx, 1]] = repeat[1] + frac[1];
            // super_cell_fractional_vectors[[idx, 2]] = repeat[2] + frac[2];

            // TODO: check if correct
            // let vec: Array1<f64> = repeat + &frac_coord;
            // super_cell_fractional_vectors[[idx, 0]] = vec[0];
            // super_cell_fractional_vectors[[idx, 1]] = vec[1];
            // super_cell_fractional_vectors[[idx, 2]] = vec[2];
            idx += 1;
        }
    }

    let mut distances: Vec<f64> = Vec::with_capacity(fractional_coordinates.n_rows * 27);
    // println!("{:?}", center_coordinates);
    // println!("{:?}", super_cell_fractional_vectors);

    // TODO: also do not do this with Array!
    let super_cell_cartesian_vectors = super_cell_fractional_vectors.dot(&cell);
    // TODO: check if good!
    for center_coord in center_coordinates.rows_iter() {
        let center_coord: Vec<f64> = center_coord.map(|&x| x).collect();
        for valid_coord in super_cell_cartesian_vectors.rows_iter() {
            let valid_coord: Vec<f64> = valid_coord.map(|&x| x).collect();

            let distance_squared = (valid_coord[0] - center_coord[0]).powi(2)
                + (valid_coord[1] - center_coord[1]).powi(2)
                + (valid_coord[2] - center_coord[2]).powi(2);
            if distance_squared < r_cut_squared + NUMERICAL_TOLERANCE
                && distance_squared > NUMERICAL_TOLERANCE
            {
                distances.push(distance_squared.sqrt());
            }
        }
    }
    distances
}

// TODO include Verlet List and Cell list to improve scaling
// Optimized version for Metropolis Monte Carlo
pub fn coordinates_within(
    r_cut: f64,
    cell: &Array2D,
    scaled_positions: &Array2D,
    super_cell_repeats: &Vec<Array1<f64>>,
) -> Array2D {
    get_coords_in_spheres(
        &cell,
        &scaled_positions,
        &scaled_positions.dot(cell),
        super_cell_repeats,
        r_cut,
    )
}

// TODO include Verlet List and Cell list to improve scaling
// Optimized version for Metropolis Monte Carlo
pub fn distances_within(
    r_cut: f64,
    cell: &Array2D,
    scaled_positions: &Array2D,
    super_cell_repeats: &Vec<Array1<f64>>,
) -> Vec<f64> {
    get_distances_in_spheres(
        &cell,
        &scaled_positions,
        &scaled_positions.dot(cell),
        super_cell_repeats,
        r_cut,
    )
}

// #[cfg(test)]
// mod tests {
//     use super::*;

//     #[test]
//     fn test_super_cell() {
//         let cell: Array2<f64> = array![[9.5, 0., 0.], [0., -12., 0.], [0., 0., 8.]];
//         // let cell: Array2<f64> = array![[9., 0., 2.3], [0., -12., 7.2], [-4., 1., 8.]];
//         let r_cut: f64 = 9.3;
//         get_super_cell(&cell, r_cut);
//         assert!(true);
//     }

//     #[test]
//     fn non_cubic_cell() {
//         let cell: Array2<f64> = array![[9., 0., 2.3], [0., -12., 7.2], [-4., 1., 8.]];
//         let fractional_coordinates: Array2<f64> = array![
//             [1., 2., 3.],
//             [3., 4., -5.6],
//             [1., 1., 1.2],
//             [-3.4, 2.2, -9.8]
//         ];
//         // let fractional_coordinates = cartesian_coordinates.dot(&invert3(&cell));
//         let center_coords_cartesian: Array2<f64> = array![[0.0, 0.0, 0.0]];
//         let r_cut: f64 = 9.3;

//         let get_points_in_sphere = get_points_in_spheres(
//             &cell,
//             &fractional_coordinates,
//             &center_coords_cartesian,
//             r_cut,
//         );

//         let get_distances_in_sphere = get_distances_in_spheres(
//             &cell,
//             &fractional_coordinates,
//             &center_coords_cartesian,
//             &get_super_cell(&cell, r_cut),
//             r_cut,
//         );

//         let get_coords_in_sphere = get_coords_in_spheres(
//             &cell,
//             &fractional_coordinates,
//             &center_coords_cartesian,
//             &get_super_cell(&cell, r_cut),
//             r_cut,
//         );
//         assert!(get_distances_in_sphere.all_close(&get_points_in_sphere.0, 0.0001));
//         assert!(get_coords_in_sphere.all_close(&get_points_in_sphere.1, 0.0001));

//         assert!(get_points_in_sphere.0.all_close(
//             &array![
//                 6.706295549705521,
//                 9.289241088485108,
//                 5.356715411518514,
//                 9.0,
//                 5.399999999999997,
//                 7.200000000000001,
//                 3.600000000000003,
//                 1.7999999999999996,
//                 6.748066389714909,
//                 9.0,
//                 9.289241088485108,
//                 9.228759396581971,
//                 9.082400563727631
//             ],
//             0.0001
//         ));

//         let r_cut: f64 = 9.3;
//         let center_coords_cartesian: Array2<f64> = array![[9., -3.1, 0.8]];

//         let get_points_in_sphere = get_points_in_spheres(
//             &cell,
//             &fractional_coordinates,
//             &center_coords_cartesian,
//             r_cut,
//         );

//         assert!(get_points_in_sphere.0.all_close(
//             &array![
//                 9.009439494219382,
//                 4.399363590338952,
//                 5.768396657651062,
//                 7.9158069708653205,
//                 4.785394445602157,
//                 6.288083968904996,
//                 8.848525300862288,
//                 3.443835071544513,
//                 6.074537019394979,
//                 4.597825573029059,
//                 7.550920473690606
//             ],
//             0.0001
//         ));
//     }
// }

// // pub fn points_in_spheres(
// //     cell: &Array2<f64>,
// //     fractional_coordinates: &Array2<f64>,
// //     center_coordinates: &Array2<f64>,
// //     r_cut: f64,
// // ) -> (Array1<f64>, Array2<f64>) {
// //     let center_coordinates_min: Array1<f64> = center_coordinates
// //         .axis_iter(Axis(1))
// //         .map(|x| min(&x))
// //         .collect();

// //     let center_coordinates_max: Array1<f64> = center_coordinates
// //         .axis_iter(Axis(1))
// //         .map(|x| max(&x))
// //         .collect();

// //     let global_min: Array1<f64> = center_coordinates_min
// //         .iter()
// //         .map(|x| x - r_cut - NUMERICAL_TOLERANCE)
// //         .collect();

// //     let global_max: Array1<f64> = center_coordinates_max
// //         .iter()
// //         .map(|x| x + r_cut + NUMERICAL_TOLERANCE)
// //         .collect();

// //     let reciprocal_cell: Array2<f64> = invert3(cell); // TODO implement inversion
// //     let reciprocal_cell = reciprocal_cell.t();
// //     let reciprocal_cell = reciprocal_cell.map(|&x| x * 2. * PI);

// //     let reciprocal_lengths: Array1<f64> = reciprocal_cell
// //         .outer_iter()
// //         .map(|vec| (&vec * &vec).sum().sqrt())
// //         .collect();

// //     let max_repeats: Array1<u32> = reciprocal_lengths // TODO: check why
// //         .map(|len| (r_cut + 0.15) * len / (2. * PI))
// //         .map(|x| x.ceil() as u32);

// //     let fractional_center_coordinates = center_coordinates.dot(&invert3(&cell));

// //     let nmin: Array1<i32> = fractional_center_coordinates
// //         .axis_iter(Axis(1))
// //         .zip(max_repeats.iter())
// //         .map(|(x, &max_r)| min(&x).floor() as i32 - max_r as i32)
// //         .collect();

// //     let nmax: Array1<i32> = fractional_center_coordinates
// //         .axis_iter(Axis(1))
// //         .zip(max_repeats.iter())
// //         .map(|(x, &max_r)| max(&x).ceil() as i32 + max_r as i32)
// //         .collect();

// //     let all_ranges: Vec<Vec<i32>> = nmin
// //         .iter()
// //         .zip(nmax.iter())
// //         .map(|(&min, &max)| (min..max).collect())
// //         .collect();

// //     let all_fcoords = fractional_coordinates.map(|&x| {
// //         if x > 0.0 {
// //             return x % 1.;
// //         }
// //         if x < 0.0 {
// //             return (x % 1.) + 1.;
// //         }
// //         return 0.0;
// //     }); // TODO: check if correcy

// //     let cartesian_coordinates_in_cell = all_fcoords.dot(cell);

// //     assert_eq!(all_ranges.len(), 3);

// //     let mut valid_coordinates: Vec<Vec<f64>> = Vec::new();
// //     for &range1 in all_ranges[0].iter() {
// //         for &range2 in all_ranges[1].iter() {
// //             for &range3 in all_ranges[2].iter() {
// //                 let image = array![range1 as f64, range2 as f64, range3 as f64];
// //                 let coords = &image.dot(cell);

// //                 // TODO:optimize with rust documentation https://docs.rs/ndarray/0.14.0/ndarray/struct.ArrayBase.html#conversions
// //                 let cartesian_coords: Vec<f64> = cartesian_coordinates_in_cell
// //                     .outer_iter()
// //                     .map(|cart_coord| (&cart_coord + coords).to_vec())
// //                     .flatten()
// //                     .collect();
// //                 let cartesian_coords = Array2::from_shape_vec(
// //                     cartesian_coordinates_in_cell.raw_dim(),
// //                     cartesian_coords,
// //                 )
// //                 .unwrap();

// //                 let bigger_global_min: Vec<bool> = cartesian_coords
// //                     .outer_iter()
// //                     .map(|cart_coord| {
// //                         cart_coord
// //                             .iter()
// //                             .zip(global_min.iter())
// //                             .all(|(coord, min)| coord > min)
// //                     })
// //                     .collect();

// //                 let smaller_global_max: Vec<bool> = cartesian_coords
// //                     .outer_iter()
// //                     .map(|cart_coord| {
// //                         cart_coord
// //                             .iter()
// //                             .zip(global_max.iter())
// //                             .all(|(coord, max)| coord < max)
// //                     })
// //                     .collect();

// //                 let valid_index: Vec<bool> = smaller_global_max
// //                     .iter()
// //                     .zip(bigger_global_min.iter())
// //                     .map(|(&smaller, &bigger)| smaller && bigger)
// //                     .collect();

// //                 if valid_index.iter().any(|&x| x) {
// //                     for (idx, &valid) in valid_index.iter().enumerate() {
// //                         if valid {
// //                             valid_coordinates.push(cartesian_coords.slice(s!(idx, ..)).to_vec());
// //                         }
// //                     }
// //                 }
// //             }
// //         }
// //     }

// //     if valid_coordinates.len() < 1 {
// //         return (array![], array![[]]); // TODO RETURN
// //     }

// //     let mut coordinates: Vec<Vec<f64>> = Vec::new();
// //     let mut distances: Vec<f64> = Vec::new();
// //     let r_cut_squared = r_cut * r_cut;
// //     for center_coord in center_coordinates.outer_iter() {
// //         for valid_coord in valid_coordinates.iter() {
// //             let distance_squared: f64 = valid_coord
// //                 .iter()
// //                 .zip(center_coord.iter())
// //                 .map(|(valid, center)| (valid - center).powi(2))
// //                 .sum();
// //             if distance_squared < r_cut_squared + NUMERICAL_TOLERANCE
// //                 && distance_squared > NUMERICAL_TOLERANCE
// //             {
// //                 distances.push(distance_squared.sqrt());
// //                 coordinates.push(valid_coord.to_vec());
// //             }
// //         }
// //     }
// //     (Array::from(distances), vec_to_array(coordinates))
// //     // return (Array1::zeros(3), Array2::zeros((3, 3))); // TODO RETURN
// // }

// // // Compute cell index from coordinates
// // fn _compute_cube_index(
// //     coordinates: &Array2<f64>,
// //     global_min: &Array1<f64>,
// //     radius: f64,
// // ) -> Array2<i32> {
// //     let vec = coordinates
// //         .outer_iter()
// //         .map(|coord| {
// //             (&coord - global_min)
// //                 .map(|x| (x / radius).floor() as i32)
// //                 .to_vec()
// //         })
// //         .flatten()
// //         .collect();
// //     Array2::from_shape_vec(coordinates.raw_dim(), vec).unwrap()
// // }
