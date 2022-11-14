use nalgebra::Matrix3;
use ndarray::prelude::*;

pub static NUMERICAL_TOLERANCE: f64 = 1e-10;

/// Calculate the cross product of two 3D-vectors.
///
/// # Arguments
///
/// * `arr1` - vector 1
/// * `arr2` - vector 2
///
/// # Returns
///
/// Cross product arr1 x arr2
pub fn cross_product(arr1: &ArrayView1<f64>, arr2: &ArrayView1<f64>) -> Array1<f64> {
    assert_eq!(arr1.dim(), 3);
    assert_eq!(arr1.dim(), arr2.dim());
    let mut cross = Array1::zeros(3);
    cross[0] = arr1[1] * arr2[2] - arr1[2] * arr2[1];
    cross[1] = arr1[2] * arr2[0] - arr1[0] * arr2[2];
    cross[2] = arr1[0] * arr2[1] - arr1[1] * arr2[0];
    cross
}

/// Get minimum value of 3D array
pub fn min_3d(vec: &Array3<f64>) -> f64 {
    let mut minimum = vec[[0, 0, 0]];
    for &v in vec {
        if v < minimum {
            minimum = v;
        }
    }
    minimum
}

/// Get minimum value of 1D array
pub fn min(vec: &ArrayView1<f64>) -> f64 {
    let mut minimum = vec[0];
    for &v in vec.iter() {
        if v < minimum {
            minimum = v;
        }
    }
    minimum
}

/// Get maximum value of 1D array
pub fn max(vec: &ArrayView1<f64>) -> f64 {
    let mut maximum = vec[0];
    for &v in vec.iter() {
        if v > maximum {
            maximum = v;
        }
    }
    maximum
}

/// Invert matrix of shape (3, 3).
pub fn invert3(arr: &Array2<f64>) -> Array2<f64> {
    let vec: Vec<f64> = arr.iter().map(|&x| x).collect();
    let m1 = Matrix3::from_vec(vec);
    let inv = m1.try_inverse().unwrap();
    let arr_2 = array![
        [inv[0], inv[1], inv[2]],
        [inv[3], inv[4], inv[5]],
        [inv[6], inv[7], inv[8]]
    ];
    arr_2
}

/// Converts nested `Vec`s to a 2-D array by cloning the elements.
///
/// **Panics** if the length of any axis overflows `isize`, if the
/// size in bytes of all the data overflows `isize`, or if not all the
/// rows have the same length.
pub fn vec_to_array<T: Clone>(v: Vec<Vec<T>>) -> Array2<T> {
    if v.is_empty() {
        return Array2::from_shape_vec((0, 0), Vec::new()).unwrap();
    }
    let nrows = v.len();
    let ncols = v[0].len();
    let mut data = Vec::with_capacity(nrows * ncols);
    for row in &v {
        assert_eq!(row.len(), ncols);
        data.extend_from_slice(&row);
    }
    Array2::from_shape_vec((nrows, ncols), data).unwrap()
}

#[cfg(test)]
mod tests {
    use super::*;
    use float_cmp::approx_eq;

    #[test]
    fn test_min() {
        let array = array![[1., -1., -3.], [2., 4., -2.], [3., 4., 5.]];
        let min: Vec<f64> = array.outer_iter().map(|x| min(&x)).collect();
        assert!(approx_eq!(f64, min[0], -3.));
        assert!(approx_eq!(f64, min[1], -2.));
        assert!(approx_eq!(f64, min[2], 3.));
    }
}
