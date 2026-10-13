use std::ops::{Add, Index, IndexMut, Mul, Sub};

use super::Vector;

/// A 3x3 matrix, indexed `m[(row, col)]` from 0.
#[doc(alias = "gp_Mat")]
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Matrix3 {
    rows: [[f64; 3]; 3],
}

impl Matrix3 {
    pub const ZERO: Matrix3 = Matrix3::from_rows([[0.0; 3]; 3]);
    pub const IDENTITY: Matrix3 =
        Matrix3::from_rows([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]);

    pub const fn from_rows(rows: [[f64; 3]; 3]) -> Self {
        Self { rows }
    }

    pub fn determinant(&self) -> f64 {
        let m = &self.rows;
        m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
    }

    /// `None` if the matrix is singular.
    pub fn inverse(&self) -> Option<Matrix3> {
        let det = self.determinant();
        if det.abs() <= crate::tolerance::RESOLUTION {
            return None;
        }
        let m = &self.rows;
        let mut inv = [[0.0; 3]; 3];
        for r in 0..3 {
            for c in 0..3 {
                // Cofactor of (c, r), divided by the determinant.
                let (r1, r2) = ((c + 1) % 3, (c + 2) % 3);
                let (c1, c2) = ((r + 1) % 3, (r + 2) % 3);
                inv[r][c] = (m[r1][c1] * m[r2][c2] - m[r1][c2] * m[r2][c1]) / det;
            }
        }
        Some(Matrix3::from_rows(inv))
    }

    /// The matrix of cofactors: `C[i][j] = (-1)^(i+j)` times the minor of `(i, j)`, which
    /// is `determinant * inverse` transposed (and stays defined for a singular matrix).
    pub fn cofactors(&self) -> Matrix3 {
        let m = &self.rows;
        Matrix3::from_rows(std::array::from_fn(|i| {
            std::array::from_fn(|j| {
                let (r1, r2) = ((i + 1) % 3, (i + 2) % 3);
                let (c1, c2) = ((j + 1) % 3, (j + 2) % 3);
                m[r1][c1] * m[r2][c2] - m[r1][c2] * m[r2][c1]
            })
        }))
    }

    pub fn transpose(&self) -> Matrix3 {
        let m = &self.rows;
        Matrix3::from_rows(std::array::from_fn(|r| std::array::from_fn(|c| m[c][r])))
    }
}

impl Index<(usize, usize)> for Matrix3 {
    type Output = f64;
    fn index(&self, (row, col): (usize, usize)) -> &f64 {
        &self.rows[row][col]
    }
}

impl IndexMut<(usize, usize)> for Matrix3 {
    fn index_mut(&mut self, (row, col): (usize, usize)) -> &mut f64 {
        &mut self.rows[row][col]
    }
}

impl Add for Matrix3 {
    type Output = Matrix3;
    fn add(self, rhs: Matrix3) -> Matrix3 {
        Matrix3::from_rows(std::array::from_fn(|r| std::array::from_fn(|c| self.rows[r][c] + rhs.rows[r][c])))
    }
}

impl Sub for Matrix3 {
    type Output = Matrix3;
    fn sub(self, rhs: Matrix3) -> Matrix3 {
        Matrix3::from_rows(std::array::from_fn(|r| std::array::from_fn(|c| self.rows[r][c] - rhs.rows[r][c])))
    }
}

impl Mul for Matrix3 {
    type Output = Matrix3;
    fn mul(self, rhs: Matrix3) -> Matrix3 {
        Matrix3::from_rows(std::array::from_fn(|r| {
            std::array::from_fn(|c| (0..3).map(|k| self.rows[r][k] * rhs.rows[k][c]).sum())
        }))
    }
}

impl Mul<Vector> for Matrix3 {
    type Output = Vector;
    fn mul(self, v: Vector) -> Vector {
        let m = &self.rows;
        Vector::new(
            m[0][0] * v.x + m[0][1] * v.y + m[0][2] * v.z,
            m[1][0] * v.x + m[1][1] * v.y + m[1][2] * v.z,
            m[2][0] * v.x + m[2][1] * v.y + m[2][2] * v.z,
        )
    }
}
