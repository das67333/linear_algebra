extern crate ndarray;
extern crate num;

use ndarray::Array2;
use num::rational::Ratio;
use num::zero;

type T = Ratio<i128>;

pub fn reduced_row_echelon_form(x: &Array2<T>) -> (Array2<T>, Array2<T>) {
    let mut result = x.clone();
    let mut factor = Array2::<T>::eye(x.nrows());
    let mut row = 0;
    for column in 0..result.ncols() {
        if row == result.nrows() {
            break;
        }
        if result[[row, column]] == zero() {
            if let Some(i) = (row..result.nrows()).find(|&i| result[[i, column]] != zero()) {
                for j in 0..result.ncols() {
                    let t = result[[i, j]];
                    result[[row, j]] += t;
                }
                for j in 0..factor.ncols() {
                    let t = factor[[i, j]];
                    factor[[row, j]] += t;
                }
            } else {
                continue;
            }
        }
        let divider = result[[row, column]];
        for j in column..result.ncols() {
            result[[row, j]] /= divider;
        }
        for j in 0..factor.ncols() {
            factor[[row, j]] /= divider;
        }
        for i in 0..result.nrows() {
            let k = result[[i, column]];
            if i == row || k == zero() {
                continue;
            }
            for j in 0..result.ncols() {
                let t = result[[row, j]] * k;
                result[[i, j]] -= t;
            }
            for j in 0..factor.ncols() {
                let t = factor[[row, j]] * k;
                factor[[i, j]] -= t;
            }
        }
        row += 1;
    }
    assert_eq!(factor.dot(&x.view()), result);
    (result, factor)
}

pub fn canonical_form(x: &Array2<T>) -> (Array2<T>, Array2<T>) {
    // x must be symmetrical
    assert!(x.is_square());
    let n = x.nrows();
    let mut result = x.clone();
    let mut conversion = Array2::<T>::eye(n);
    let mut row = 0;
    for column in 0..n {
        if row == n {
            break;
        }
        if result[[row, column]] == zero() {
            if let Some(i) = (row..n).find(|&i| result[[i, column]] != zero()) {
                for j in 0..n {
                    let t = result[[i, j]];
                    result[[row, j]] += t;
                }
                for j in 0..n {
                    let t = result[[j, i]];
                    result[[j, row]] += t;
                }
                for j in 0..n {
                    let t = conversion[[j, i]];
                    conversion[[j, row]] += t;
                }
            } else {
                continue;
            }
        }
        for i in row + 1..n {
            let k = result[[i, column]] / result[[row, column]];
            if k == zero() {
                continue;
            }
            for j in 0..n {
                let t = result[[row, j]] * k;
                result[[i, j]] -= t;
            }
            for j in 0..n {
                let t = result[[j, row]] * k;
                result[[j, i]] -= t;
            }
            for j in 0..n {
                let t = conversion[[j, row]] * k;
                conversion[[j, i]] -= t;
            }
        }
        row += 1;
    }
    assert_eq!(conversion.t().dot(&x.view()).dot(&conversion.view()), result);
    (result, conversion)
}

pub fn inverse(x: &Array2<T>) -> Array2<T> {
    assert!(x.is_square());
    let (r, f) = reduced_row_echelon_form(&x);
    assert_eq!(r, Array2::<T>::eye(x.nrows()));
    f
}

/* pub fn determinant(&self) -> Rational {
    assert_eq!(self.size().0, self.size().1);
    let n = self.size().0;
    let mut result = self;
    let mut sign = 1;
    let zero = Rational::from(0);
    let one = Rational::from(1);
    let mut prev = one;
    for k in 0..n - 1 {
        if result[k][k] == zero {
            if let Some(i) = (k + 1..n).find(|&i| result[i][k] != zero) {
                result._data.swap(k, i);
                sign = -sign;
            } else {
                return zero;
            }
        }
        for i in k + 1..n {
            for j in k + 1..n {
                result[i][j] =
                    (result[i][j] * result[k][k] - result[i][k] * result[k][j]) / prev;
            }
        }
        prev = result[k][k];
    }
    result[n - 1][n - 1] * sign
} */
