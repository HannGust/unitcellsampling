extern crate openblas_src;
use cblas::*;
use ndarray::Array2;
use transpose::transpose;

#[derive(Debug, Clone)]
pub struct Array2D {
    pub data: Vec<f64>,
    pub n_rows: usize,
    pub n_cols: usize,
}

fn flatten<T: Clone>(nested: &[Vec<T>]) -> Vec<T> {
    nested.iter().flat_map(|row| row.clone()).collect()
}
impl Array2D {
    pub fn map(mut self, f: fn(&f64) -> f64) -> Array2D {
        self.data = self.data.iter().map(f).collect::<Vec<f64>>();
        self
    }
    pub fn from_ndarray(array: Array2<f64>) -> Array2D {
        Array2D::new(array.clone().into_raw_vec(), array.dim().0, array.dim().1)
    }

    pub fn new(vec: Vec<f64>, n_rows: usize, n_cols: usize) -> Array2D {
        assert!(vec.len() % n_rows == 0);
        assert!(vec.len() % n_cols == 0);
        assert!(n_cols * n_rows == vec.len());
        Self {
            data: vec,
            n_rows,
            n_cols,
        }
    }
    pub fn from_rows(elements: &[Vec<f64>]) -> Self {
        let row_len = elements.get(0).map(Vec::len).unwrap_or(0);
        if !elements.iter().all(|row| row.len() == row_len) {
            panic!("Rows were not all {} elements long", row_len);
        }
        Array2D {
            data: flatten(elements),
            n_rows: elements.len(),
            n_cols: row_len,
        }
    }

    pub fn from_column_major(elements: &[f64], n_rows: usize, n_cols: usize) -> Self {
        let total_len = n_rows * n_cols;
        if total_len != elements.len() {
            panic!(
                "The number of elements ({}) did not match the expected size ({})",
                elements.len(),
                total_len
            );
        }
        let indices_row_major =
            (0..n_rows).flat_map(move |row| (0..n_cols).map(move |column| (row, column)));
        let data = indices_row_major
            .map(|(row, column)| {
                let index = column * n_rows + row;
                elements[index].clone()
            })
            .collect();
        Array2D {
            data,
            n_rows,
            n_cols,
        }
    }
    pub fn get_index(&self, row: usize, column: usize) -> Option<usize> {
        if row < self.n_rows && column < self.n_cols {
            Some(row * 3 + column)
        } else {
            None
        }
    }
    pub fn get(&self, row: usize, column: usize) -> Option<&f64> {
        self.get_index(row, column).map(|index| &self.data[index])
    }

    pub fn get_mut(&mut self, row: usize, column: usize) -> Option<&mut f64> {
        self.get_index(row, column)
            .map(move |index| &mut self.data[index])
    }

    pub fn row_iter(&self, row_index: usize) -> impl Iterator<Item = &f64> {
        let start = self.get_index(row_index, 0).expect(&format!(
            "Row index, {}, was out of bounds (>= number of rows, {})",
            row_index, self.n_rows,
        ));
        let end = start + self.n_cols;
        self.data[start..end].iter()
    }

    pub fn rows_iter(&self) -> impl Iterator<Item = impl Iterator<Item = &f64>> {
        (0..self.n_rows).map(move |row_index| self.row_iter(row_index))
    }

    pub fn dot(&self, b: &Array2D) -> Array2D {
        assert_eq!(self.n_cols, b.n_rows);

        let m: usize = self.n_rows;
        let n: usize = self.n_cols;
        let k: usize = b.n_cols;

        let mut c = Array2D::new(vec![0f64; (k * m) as usize], m, k);

        let m: i32 = m as i32;
        let n: i32 = n as i32;
        let k: i32 = k as i32;
        unsafe {
            dgemm(
                Layout::RowMajor,
                Transpose::None,
                Transpose::None,
                m,
                n,
                k,
                1.0,
                &self.data,
                k,
                &b.data,
                n,
                1.0,
                &mut c.data,
                n,
            );
        }
        c
    }

    pub fn t(&self) -> Array2D {
        let mut transposed = Array2D::new(
            vec![0f64; self.n_cols * self.n_rows],
            self.n_cols,
            self.n_rows,
        );
        transpose(&self.data, &mut transposed.data, self.n_cols, self.n_rows);
        transposed
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use float_cmp::approx_eq;
    use ndarray::prelude::*;

    #[test]
    fn array2d() {
        let a = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let b = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0];

        assert!(
            Array2D::new(a.clone(), 2, 3)
                .dot(&Array2D::new(b.clone(), 3, 3))
                .data
                == vec![30., 36., 42., 66., 81., 96.]
        );
        assert!(
            Array::from_shape_vec((2, 3), a.clone())
                .unwrap()
                .dot(&Array::from_shape_vec((3, 3), b.clone()).unwrap())
                .into_raw_vec()
                == vec![30., 36., 42., 66., 81., 96.]
        );
        let arr = Array2D::new(a.clone(), 2, 3).dot(&Array2D::new(b.clone(), 3, 3));
        assert!(approx_eq!(f64, *arr.get(0, 0).unwrap(), 30.0, ulps = 2));
        assert!(approx_eq!(f64, *arr.get(0, 1).unwrap(), 36.0, ulps = 2));
        assert!(approx_eq!(f64, *arr.get(0, 2).unwrap(), 42.0, ulps = 2));
        assert!(approx_eq!(f64, *arr.get(1, 0).unwrap(), 66.0, ulps = 2));
        assert!(approx_eq!(f64, *arr.get(1, 1).unwrap(), 81.0, ulps = 2));
        assert!(approx_eq!(f64, *arr.get(1, 2).unwrap(), 96.0, ulps = 2));
        for row_iter in arr.rows_iter() {
            for element in row_iter {
                print!("{} ", element);
            }
            println!();
        }
        for (a, b) in arr.rows_iter().zip(arr.rows_iter()) {
            let mut vec: Vec<f64> = Vec::with_capacity(arr.n_cols);
            vec = a.zip(b).map(|(x, y)| x + y).collect();
        }

        println!(
            "{:?} {:?}",
            Array2D::from_ndarray(array![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]),
            Array2D::new(vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0], 2, 3)
        );

        println!("{:?}", arr.t())
    }
}
