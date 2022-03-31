extern crate num;
pub type Rational = num::rational::Ratio<i32>;

use std::fs::File;
use std::io::{Read, Result};
use std::ops::*;
use std::{mem, str, vec};

#[derive(Debug, Clone)]
pub struct Matrix<T> {
    _data: Vec<Vec<T>>,
}

impl<T> Matrix<T> {
    pub fn size(&self) -> (usize, usize) {
        (self._data.len(), self._data[0].len())
    }

    pub fn new(data: Vec<Vec<T>>) -> Matrix<T>
    where
        T: Clone,
    {
        assert!(data.len() != 0 && data[0].len() != 0);
        for i in 0..data.len() - 1 {
            assert_eq!(data[i].len(), data[i + 1].len());
        }
        Matrix { _data: data }
    }

    pub fn from_file(path: &str) -> Result<Matrix<T>>
    where
        T: Clone + str::FromStr,
        <T as str::FromStr>::Err: std::fmt::Debug,
    {
        let mut file = File::open(path)?;
        let mut contents = String::new();
        file.read_to_string(&mut contents)?;
        let mut data = Vec::new();
        for line in contents.trim().split('\n') {
            let mut row: Vec<T> = Vec::new();
            for item in line.trim().split(',') {
                row.push(item.parse().expect("parsing string failed"));
            }
            data.push(row);
        }
        Ok(Matrix::new(data))
    }

    pub fn zero(n: usize, m: usize) -> Matrix<T>
    where
        T: Clone + From<i32>,
    {
        Matrix::new(vec![vec![T::from(0); m]; n])
    }

    pub fn one(n: usize) -> Matrix<T>
    where
        T: Clone + From<i32>,
    {
        let mut result = Matrix::zero(n, n);
        for i in 0..result.size().0 {
            result[i][i] = T::from(1);
        }
        result
    }

    pub fn transposed(&self) -> Matrix<T>
    where
        T: Clone + From<i32>,
    {
        let mut result = Matrix::zero(self.size().1, self.size().0);
        for i in 0..self.size().0 {
            for j in 0..self.size().1 {
                result[j][i] = self[i][j].clone();
            }
        }
        result
    }

    pub fn reduced(&self, n: usize, m: usize) -> Matrix<T>
    where
        T: Clone + From<i32>,
    {
        assert!(n <= self.size().0 && m <= self.size().1);
        let mut result = Matrix::zero(n, m);
        for i in 0..n {
            result[i].clone_from_slice(&self[i][0..m]);
        }
        result
    }
}

#[macro_export]
macro_rules!matrix {
    [ $( [ $( $d:expr ),* ] ),* ] => {
        Matrix::new(vec![
            $(
                vec![$($d),*],
            )*
        ])
    }
}

impl Matrix<Rational> {
    pub fn gaussian_elimination(&self) -> (Matrix<Rational>, Matrix<Rational>) {
        let mut result = Matrix::<Rational>::from(self.clone());
        let mut factor = Matrix::<Rational>::one(self.size().0);
        let mut lead_row = 0;
        let zero = Rational::from(0);
        for lead_column in 0..result.size().1 {
            if lead_row == result.size().0 {
                break;
            }
            if result[lead_row][lead_column] == zero {
                if let Some(i) =
                    (lead_row..result.size().0).find(|&i| result[i][lead_column] != zero)
                {
                    result._data.swap(lead_row, i);
                    factor._data.swap(lead_row, i);
                } else {
                    continue;
                }
            }
            let divider = result[lead_row][lead_column];
            for j in lead_column..result.size().1 {
                result[lead_row][j] /= divider;
            }
            for j in 0..factor.size().1 {
                factor[lead_row][j] /= divider;
            }
            for i in 0..result.size().0 {
                if i == lead_row || result[i][lead_column] == zero {
                    continue;
                }
                for j in lead_column + 1..result.size().1 {
                    let mut t = result[lead_row][j];
                    t *= result[i][lead_column];
                    result[i][j] -= t;
                }
                for j in 0..factor.size().1 {
                    let mut t = factor[lead_row][j];
                    t *= result[i][lead_column];
                    factor[i][j] -= t;
                }
                result[i][lead_column] = zero;
            }
            lead_row += 1;
        }
        (result, factor)
    }

    pub fn inversed(&self) -> Option<Matrix<Rational>> {
        assert_eq!(self.size().0, self.size().1);
        let (m1, m2) = self.gaussian_elimination();
        if m1 == Matrix::<Rational>::one(self.size().0) {
            Some(m2)
        } else {
            None
        }
    }

    pub fn determinant(&self) -> Rational {
        assert_eq!(self.size().0, self.size().1);
        let n = self.size().0;
        let mut result = self.clone();
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
    }
}

impl<T> Index<usize> for Matrix<T> {
    type Output = [T];
    fn index(&self, i: usize) -> &Self::Output {
        &self._data[i]
    }
}

impl<T> IndexMut<usize> for Matrix<T> {
    fn index_mut(&mut self, i: usize) -> &mut Self::Output {
        &mut self._data[i]
    }
}

impl<T> std::fmt::Display for Matrix<T>
where
    T: ToString,
{
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let max_displayable_size = (8, 8);
        let mut max_width = 0;
        let (mut i, mut j) = (0, 0);
        while i < self.size().0 {
            if i == max_displayable_size.0 - 3 && max_displayable_size.0 < self.size().0 {
                i = self.size().0 - 1;
            }
            while j < self.size().1 {
                if j == max_displayable_size.1 - 3 && max_displayable_size.1 < self.size().1 {
                    j = self.size().1 - 1;
                    continue;
                }
                max_width = std::cmp::max(max_width, self[i][j].to_string().len());
                j += 1;
            }
            j = 0;
            i += 1;
        }
        let mut result = String::new();
        i = 0;
        while i < self.size().0 {
            if i == max_displayable_size.0 - 3 && max_displayable_size.0 < self.size().0 {
                result.push_str(&format!(
                    "|{:.^1$}|\n",
                    "",
                    (max_width + 1) * std::cmp::min(self.size().1, max_displayable_size.1) - 1
                ));
                i = self.size().0 - 1;
                continue;
            }
            result.push('|');
            while j < self.size().1 {
                if j == max_displayable_size.1 - 3 && max_displayable_size.1 < self.size().1 {
                    result.push_str(&format!("{:^1$}", "...", (max_width + 1) * 2));
                    j = self.size().1 - 1;
                    continue;
                }
                result.push_str(&format!("{:^1$}", self[i][j].to_string(), max_width + 1));
                j += 1;
            }
            result.pop();
            result.push_str("|\n");
            j = 0;
            i += 1;
        }
        write!(f, "Matrix: {:?}\n{}", self.size(), result)
    }
}

impl<T> PartialEq for Matrix<T>
where
    T: PartialEq,
{
    fn eq(&self, other: &Self) -> bool {
        if self.size() != other.size() {
            return false;
        }
        for i in 0..self.size().0 {
            for j in 0..self.size().1 {
                if self[i][j] != other[i][j] {
                    return false;
                }
            }
        }
        true
    }
}

impl<T> MulAssign<T> for Matrix<T>
where
    T: Clone + MulAssign<T>,
{
    fn mul_assign(&mut self, rhs: T) {
        for i in 0..self.size().0 {
            for j in 0..self.size().1 {
                self[i][j] *= rhs.clone();
            }
        }
    }
}

impl<T> Mul<T> for Matrix<T>
where
    T: Clone + MulAssign<T>,
{
    type Output = Matrix<T>;
    fn mul(self, rhs: T) -> Self::Output {
        let mut result = self.clone();
        result *= rhs;
        result
    }
}

impl<T> AddAssign for Matrix<T>
where
    T: Clone + AddAssign,
{
    fn add_assign(&mut self, rhs: Self) {
        assert_eq!(self.size(), rhs.size());
        for i in 0..self.size().0 {
            for j in 0..self.size().1 {
                self[i][j] += rhs[i][j].clone();
            }
        }
    }
}

impl<T> Add for Matrix<T>
where
    T: Clone + AddAssign,
{
    type Output = Matrix<T>;
    fn add(self, rhs: Self) -> Self::Output {
        let mut result = self.clone();
        result += rhs;
        result
    }
}

impl<T> MulAssign for Matrix<T>
where
    T: Clone + From<i32> + AddAssign<<T as Mul>::Output> + Mul,
{
    fn mul_assign(&mut self, rhs: Self) {
        assert_eq!(self.size().1, rhs.size().0);
        let mut result = Matrix::zero(self.size().0, rhs.size().1);
        for i in 0..result.size().0 {
            for j in 0..result.size().1 {
                for k in 0..self.size().1 {
                    result[i][j] += self[i][k].clone() * rhs[k][j].clone();
                }
            }
        }
        mem::swap(self, &mut result);
    }
}

impl<T> Mul for Matrix<T>
where
    T: Clone + From<i32> + AddAssign<<T as Mul>::Output> + Mul,
{
    type Output = Matrix<T>;
    fn mul(self, rhs: Self) -> Self::Output {
        let mut result = self.clone();
        result *= rhs;
        result
    }
}

impl From<Matrix<i32>> for Matrix<Rational> {
    fn from(source: Matrix<i32>) -> Self {
        let mut result = Matrix::<Rational>::zero(source.size().0, source.size().1);
        for i in 0..source.size().0 {
            for j in 0..source.size().1 {
                result[i][j] = Rational::from(source[i][j]);
            }
        }
        result
    }
}
