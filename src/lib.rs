// extern crate num;
extern crate rug;

use std::fs::File;
use std::io::prelude::*;
use std::{mem, vec};

// type Num = num::BigInt;
type Num = rug::Integer;

fn zero() -> Num {
    Num::from(0)
}

fn one() -> Num {
    Num::from(1)
}

pub fn gcd(a: Num, b: Num) -> Num {
    /* if b == zero() {
        Num::abs(a)
    } else {
        gcd(b.clone(), a % b)
    } */
    a.gcd(&b)
}

pub fn lcm(a: Num, b: Num) -> Num {
    // a.clone() * b.clone() / gcd(a, b)
    a.lcm(&b)
}

#[derive(Debug, Clone)]
pub struct Matrix {
    _size: (usize, usize),
    _data: Vec<Num>,
}

impl Matrix {
    pub fn size(&self) -> (usize, usize) {
        self._size
    }

    pub fn new(data: Vec<Vec<Num>>) -> Matrix {
        assert!(data.len() != 0 && data[0].len() != 0);
        for i in 0..data.len() - 1 {
            assert_eq!(data[i].len(), data[i + 1].len());
        }
        let size = (data.len(), data[0].len());
        let mut result: Vec<Num> = Vec::new();
        for row in data {
            result.extend(row);
        }
        Matrix {
            _size: size,
            _data: result,
        }
    }

    pub fn from_file(path: &str) -> std::io::Result<Matrix> {
        let mut file = File::open(path)?;
        let mut contents = String::new();
        file.read_to_string(&mut contents)?;
        let mut data = Vec::new();
        for line in contents.trim().split('\n') {
            let mut row: Vec<Num> = Vec::new();
            for item in line.trim().split(',') {
                row.push(item.parse().expect("parsing string failed"));
            }
            data.push(row);
        }
        Ok(Matrix::new(data))
    }

    pub fn zero(n: usize, m: usize) -> Matrix {
        Matrix::new(vec![vec![zero(); m]; n])
    }

    pub fn one(n: usize) -> Matrix {
        let mut result = Matrix::zero(n, n);
        for i in 0..result.size().0 {
            result[i][i] = one();
        }
        result
    }

    pub fn transposed(&self) -> Matrix {
        let mut result = Matrix::zero(self.size().1, self.size().0);
        for i in 0..self.size().0 {
            for j in 0..self.size().1 {
                result[j][i] = self[i][j].clone();
            }
        }
        result
    }

    pub fn reduced(&self, n: usize, m: usize) -> Matrix {
        assert!(n <= self.size().0 && m <= self.size().1);
        let mut result = Matrix::zero(n, m);
        for i in 0..n {
            result[i].clone_from_slice(&self[i][0..m]);
        }
        result
    }

    fn transformation1(&mut self, i1: usize, i2: usize, lambda: Num) {
        for j in 0..self.size().1 {
            let temp = self[i2][j].clone();
            self[i1][j] += temp * lambda.clone();
        }
    }

    fn transformation2(&mut self, i1: usize, i2: usize) {
        self._data.swap(i1, i2);
    }

    fn transformation3(&mut self, i: usize, lambda: Num) {
        for j in 0..self.size().1 {
            self[i][j] *= lambda.clone();
        }
    }

    fn extract_gcd(row: &mut [Num]) -> Num {
        let mut result = zero();
        let mut j1 = 0;
        while row[j1] == zero() {
            j1 += 1;
            if j1 == row.len() {
                return one();
            }
        }
        for j2 in j1..row.len() {
            result = gcd(row[j2].clone(), result);
        }
        result *= Num::signum(row[j1].clone());
        for j in 0..row.len() {
            row[j] /= result.clone();
        }
        result
    }

    pub fn gaussian_row_echelon_form(&self) -> (Matrix, Num, Num) {
        let mut result = self.clone();
        let mut numerator: Num = one();
        let mut denominator: Num = one();
        let mut lead_row = 0;
        for i in 0..self.size().0 {
            numerator *= Matrix::extract_gcd(&mut result[i]);
        }
        for j in 0..self.size().1 {
            if lead_row == self.size().0 {
                break;
            }
            let mut i = lead_row;
            let mut are_zeros = false;
            while result[i][j] == zero() {
                i += 1;
                if i == self.size().0 {
                    are_zeros = true;
                    break;
                }
            }
            if are_zeros {
                continue;
            }
            if i != lead_row {
                result.transformation2(lead_row, i);
                numerator *= Num::from(-1);
            }
            for i in lead_row + 1..self.size().0 {
                if result[i][j] == zero() {
                    continue;
                }
                let lcm = lcm(result[lead_row][j].clone(), result[i][j].clone());
                denominator *= lcm.clone() / result[i][j].clone();
                result.transformation3(i, lcm.clone() / result[i][j].clone());
                result.transformation1(i, lead_row, -lcm.clone() / result[lead_row][j].clone());
                numerator *= Matrix::extract_gcd(&mut result[i][j..]);
                let gcd = gcd(numerator.clone(), denominator.clone());
                numerator /= gcd.clone();
                denominator /= gcd.clone();
            }
            lead_row += 1;
        }
        (result, numerator, denominator)
    }

    pub fn determinant_gaussian(&self) -> Num {
        assert_eq!(self.size().0, self.size().1);
        let (triangular, mut numerator, denominator) = self.gaussian_row_echelon_form();
        for i in 0..self.size().0 {
            numerator *= triangular[i][i].clone();
        }
        numerator / denominator
    }

    pub fn determinant_bareiss(&self) -> Num {
        let mut result = self.clone();
        let n = result.size().0;
        let mut sign = 1;
        let mut prev = one();
        for k in 0..n - 1 {
            if result[k][k] == zero() {
                let mut i = k + 1;
                while result[i][k] == zero() {
                    i += 1;
                    if i == n {
                        return zero();
                    }
                }
                result.transformation2(k, i);
                sign *= -1;
            }
            for i in k + 1..n {
                for j in k + 1..n {
                    result[i][j] = (result[i][j].clone() * result[k][k].clone() - 
                    result[i][k].clone() * result[k][j].clone()) / prev.clone();
                }
            }
            prev = result[k][k].clone();
        }
        result[n - 1][n - 1].clone() * Num::from(sign)
    }
}

#[macro_export]
macro_rules! matrix {
    ( $( $x:expr ),* ) => {
        {
            let mut result = Vec::new();
            $(
                result.push($x.to_vec());
            )*
            Matrix::new(result)
        }
    };
}

impl std::ops::Index<usize> for Matrix {
    type Output = [Num];
    fn index(&self, i: usize) -> &Self::Output {
        let i1 = i * self.size().1;
        let i2 = (i + 1) * self.size().1;
        &self._data[i1..i2]
    }
}

impl std::ops::IndexMut<usize> for Matrix {
    fn index_mut(&mut self, i: usize) -> &mut Self::Output {
        let i1 = i * self.size().1;
        let i2 = (i + 1) * self.size().1;
        &mut self._data[i1..i2]
    }
}

impl std::fmt::Display for Matrix {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let max_displayable_size = (5, 5);
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

impl std::ops::MulAssign<Num> for Matrix {
    fn mul_assign(&mut self, rhs: Num) {
        for i in 0..self.size().0 {
            for j in 0..self.size().1 {
                self[i][j] *= rhs.clone();
            }
        }
    }
}

impl std::ops::Mul<Num> for Matrix {
    type Output = Matrix;
    fn mul(self, rhs: Num) -> Self::Output {
        let mut result = self.clone();
        result *= rhs;
        result
    }
}

impl std::ops::Mul<Matrix> for Num {
    type Output = Matrix;
    fn mul(self, rhs: Matrix) -> Self::Output {
        rhs * self
    }
}

impl std::ops::AddAssign for Matrix {
    fn add_assign(&mut self, rhs: Self) {
        assert_eq!(self.size(), rhs.size());
        for i in 0..self.size().0 {
            for j in 0..self.size().1 {
                self[i][j] += rhs[i][j].clone();
            }
        }
    }
}

impl std::ops::Add for Matrix {
    type Output = Matrix;
    fn add(self, rhs: Self) -> Self::Output {
        let mut result = self.clone();
        result += rhs;
        result
    }
}

impl std::ops::MulAssign for Matrix {
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

impl std::ops::Mul for Matrix {
    type Output = Matrix;
    fn mul(self, rhs: Self) -> Self::Output {
        let mut result = self.clone();
        result *= rhs;
        result
    }
}
