//! OCCT's `math_FunctionSetRoot` (the root of a system of equations near a starting point,
//! between bounds: Newton's method with a line search, falling back to the steepest descent,
//! sliding along the bounds it meets) and the dense solvers it calls, `math_Gauss` (LU with
//! partial pivoting, `LU_Decompose` / `LU_Solve` of math_Recipes), `math_SVD`
//! (`SVD_Decompose` / `SVD_Solve`) and `math_GaussLeastSquare` (LU of the normal equations).
//!
//! The port is written for square systems but keeps OCCT's general shapes: `DF` has
//! `NbEquations()` rows and `NbVariables()` columns. A square system still reaches the
//! non-square solvers: with some (not all) unknowns at a bound the constrained
//! `SearchDirection` solves the sub-problem of the free unknowns, which has more equations
//! than unknowns and goes to `math_GaussLeastSquare`. (The wide branch, `math_SVD` of a
//! matrix with more columns than rows, is ported too; a square system never reaches it.)
//!
//! Left out: the `theStopOnDivergent` mode. Every return of OCCT's `Perform` is written as
//! `Done = false; if (!theStopOnDivergent || !myIsDivergent) { Done = ...; State = ...; }`,
//! and with the default `theStopOnDivergent = false` the condition always holds, so
//! `myIsDivergent` is never read; the port keeps only that path. The one-variable case of
//! the same solver is [`super::function_root`].

use super::adaptor::epsilon;
use super::brent_minimum::BrentMinimum;
use super::surface_analysis::is_infinite;
use std::ops::{Index, IndexMut};

/// `Precision::SquareConfusion()`.
const SQUARE_CONFUSION: f64 = 1e-14;

/// `Precision::Infinite()`.
const INFINITE: f64 = 2e100;

/// `math_FunctionSetWithDerivatives`.
#[doc(alias = "math_FunctionSetWithDerivatives")]
pub(crate) trait FunctionSetWithDerivatives {
    /// `NbVariables()`.
    fn nb_variables(&self) -> usize;
    /// `NbEquations()`.
    fn nb_equations(&self) -> usize;
    /// `Values(X, F, D)`: F (n) and the Jacobian D (row-major, n*n, D[i*n + j] = dF_i/dX_j); false where it cannot be evaluated.
    /// (For a system that is not square, D has `nb_equations()` rows of `nb_variables()` entries.)
    fn values(&mut self, x: &[f64], f: &mut [f64], d: &mut [f64]) -> bool;
    /// `Value(X, F)`: F alone; by default through [`values`](Self::values) with the
    /// Jacobian thrown away.
    fn value(&mut self, x: &[f64], f: &mut [f64]) -> bool {
        let mut d = vec![0.0; self.nb_equations() * self.nb_variables()];
        self.values(x, f, &mut d)
    }
    /// `GetStateNumber()` (default: nothing).
    fn get_state_number(&mut self) {}
}

/// A `math_Matrix` indexed from 1 as OCCT's (row, column), stored row-major: the recipes
/// below are translated with their own indices.
#[derive(Clone, Debug)]
struct Matrix {
    rows: usize,
    cols: usize,
    data: Vec<f64>,
}

impl Matrix {
    fn zeros(rows: usize, cols: usize) -> Self {
        Matrix { rows, cols, data: vec![0.0; rows * cols] }
    }

    /// A matrix from row-major entries.
    fn from_row_major(rows: usize, cols: usize, entries: &[f64]) -> Self {
        Matrix { rows, cols, data: entries[..rows * cols].to_vec() }
    }
}

impl Index<(usize, usize)> for Matrix {
    type Output = f64;
    fn index(&self, (i, j): (usize, usize)) -> &f64 {
        &self.data[(i - 1) * self.cols + (j - 1)]
    }
}

impl IndexMut<(usize, usize)> for Matrix {
    fn index_mut(&mut self, (i, j): (usize, usize)) -> &mut f64 {
        &mut self.data[(i - 1) * self.cols + (j - 1)]
    }
}

/// `SIGN(a, b)` of math_Recipes: `|a|` with the sign of `b` (`b >= 0` counting as positive).
fn sign(a: f64, b: f64) -> f64 {
    if b >= 0.0 { a.abs() } else { -a.abs() }
}

/// `PYTHAG(a, b)` of math_Recipes: `sqrt(a^2 + b^2)` without overflow.
fn pythag(a: f64, b: f64) -> f64 {
    let at = a.abs();
    let bt = b.abs();
    let mut ct = 0.0;
    if at > bt {
        ct = bt / at;
        ct = at * (1.0 + ct * ct).sqrt();
    } else if bt != 0.0 {
        ct = at / bt;
        ct = bt * (1.0 + ct * ct).sqrt();
    }
    ct
}

/// `LU_Decompose(a, indx, d, TINY)` of math_Recipes: Crout's LU decomposition with partial
/// pivoting in place, `indx` (from 1) recording the row permutation and `d` its parity.
/// False (`math_Status_SingularMatrix`) where a row's largest entry or a pivot is not above
/// `tiny` in size.
#[doc(alias = "LU_Decompose")]
fn lu_decompose(a: &mut Matrix, indx: &mut [usize], d: &mut f64, tiny: f64) -> bool {
    let n = a.rows;
    let mut imax = 0;
    // (`vv` from 1, the 0 entry unused.)
    let mut vv = vec![0.0; n + 1];
    *d = 1.0;
    for i in 1..=n {
        let mut big: f64 = 0.0;
        for j in 1..=n {
            let temp = a[(i, j)].abs();
            if temp > big {
                big = temp;
            }
        }
        if big <= tiny {
            return false;
        }
        vv[i] = 1.0 / big;
    }
    // (OCCT's progress scope can abort the loop; without one it never does.)
    for j in 1..=n {
        for i in 1..j {
            let mut sum = a[(i, j)];
            for k in 1..i {
                sum -= a[(i, k)] * a[(k, j)];
            }
            a[(i, j)] = sum;
        }
        let mut big = 0.0;
        for i in j..=n {
            let mut sum = a[(i, j)];
            for k in 1..j {
                sum -= a[(i, k)] * a[(k, j)];
            }
            a[(i, j)] = sum;
            // The comparison updates imax even for a NaN or infinite argument (OCCT #25559).
            let dum = vv[i] * sum.abs();
            if dum < big {
                continue;
            }
            big = dum;
            imax = i;
        }
        if j != imax {
            for k in 1..=n {
                let dum = a[(imax, k)];
                a[(imax, k)] = a[(j, k)];
                a[(j, k)] = dum;
            }
            *d = -*d;
            vv[imax] = vv[j];
        }
        indx[j - 1] = imax;
        if a[(j, j)].abs() <= tiny {
            return false;
        }
        if j != n {
            let dum = 1.0 / a[(j, j)];
            for i in j + 1..=n {
                a[(i, j)] *= dum;
            }
        }
    }
    true
}

/// `LU_Solve(a, indx, b)` of math_Recipes: forward and back substitution through the
/// decomposition, `b` (from 0 here) replaced by the solution.
#[doc(alias = "LU_Solve")]
fn lu_solve(a: &Matrix, indx: &[usize], b: &mut [f64]) {
    let n = a.rows;
    // `ii`: the first nonzero entry of b, 0 for none yet.
    let mut ii = 0;
    for i in 1..=n {
        let ip = indx[i - 1];
        let mut sum = b[ip - 1];
        b[ip - 1] = b[i - 1];
        if ii != 0 {
            for j in ii..i {
                sum -= a[(i, j)] * b[j - 1];
            }
        } else if sum != 0.0 {
            // (C's `if (sum)`: NaN counts as nonzero.)
            ii = i;
        }
        b[i - 1] = sum;
    }
    for i in (1..=n).rev() {
        let mut sum = b[i - 1];
        for j in i + 1..=n {
            sum -= a[(i, j)] * b[j - 1];
        }
        b[i - 1] = sum / a[(i, i)];
    }
}

/// `SVD_Decompose(a, w, v)` of math_Recipes (Golub and Reinsch, as in Numerical Recipes'
/// `svdcmp`): `a` (m x n, m >= n) is replaced by U, `w` (from 1) gets the singular values
/// and `v` (n x n) V, with A = U diag(w) V^T. False (`math_Status_NoConvergence`) where an
/// implicit QR sweep takes more than 30 iterations.
#[doc(alias = "SVD_Decompose")]
// (Some assignments are overwritten before being read, as in OCCT.)
#[allow(unused_assignments)]
fn svd_decompose(a: &mut Matrix, w: &mut [f64], v: &mut Matrix) -> bool {
    let m = a.rows;
    let n = a.cols;
    // (`rv1` and `w` from 1, the 0 entries unused.)
    let mut rv1 = vec![0.0; n + 1];
    let mut l = 0;
    let mut nm = 0;
    let (mut f, mut h, mut s, mut x, mut y, mut z);
    let mut anorm: f64 = 0.0;
    let mut g: f64 = 0.0;
    let mut scale: f64 = 0.0;

    // Householder reduction to bidiagonal form.
    for i in 1..=n {
        l = i + 1;
        rv1[i] = scale * g;
        g = 0.0;
        s = 0.0;
        scale = 0.0;
        if i <= m {
            for k in i..=m {
                let aki = a[(k, i)];
                if aki > 0.0 {
                    scale += aki;
                } else {
                    scale -= aki;
                }
            }
            if scale != 0.0 {
                for k in i..=m {
                    a[(k, i)] /= scale;
                    s += a[(k, i)] * a[(k, i)];
                }
                f = a[(i, i)];
                g = -sign(s.sqrt(), f);
                h = f * g - s;
                a[(i, i)] = f - g;
                if i != n {
                    for j in l..=n {
                        s = 0.0;
                        for k in i..=m {
                            s += a[(k, i)] * a[(k, j)];
                        }
                        f = s / h;
                        for k in i..=m {
                            a[(k, j)] += f * a[(k, i)];
                        }
                    }
                }
                for k in i..=m {
                    a[(k, i)] *= scale;
                }
            }
        }
        w[i] = scale * g;
        g = 0.0;
        s = 0.0;
        scale = 0.0;
        if i <= m && i != n {
            for k in l..=n {
                let aik = a[(i, k)];
                if aik > 0.0 {
                    scale += aik;
                } else {
                    scale -= aik;
                }
            }
            if scale != 0.0 {
                for k in l..=n {
                    a[(i, k)] /= scale;
                    s += a[(i, k)] * a[(i, k)];
                }
                f = a[(i, l)];
                g = -sign(s.sqrt(), f);
                h = f * g - s;
                a[(i, l)] = f - g;
                for k in l..=n {
                    rv1[k] = a[(i, k)] / h;
                }
                if i != m {
                    for j in l..=m {
                        s = 0.0;
                        for k in l..=n {
                            s += a[(j, k)] * a[(i, k)];
                        }
                        for k in l..=n {
                            a[(j, k)] += s * rv1[k];
                        }
                    }
                }
                for k in l..=n {
                    a[(i, k)] *= scale;
                }
            }
        }
        let mut aw = w[i];
        if aw < 0.0 {
            aw = -aw;
        }
        let mut ar = rv1[i];
        if ar > 0.0 {
            ar += aw;
        } else {
            ar = aw - ar;
        }
        if anorm < ar {
            anorm = ar;
        }
    }
    // Accumulation of the right-hand transformations.
    for i in (1..=n).rev() {
        if i < n {
            if g != 0.0 {
                for j in l..=n {
                    // (Double division, to avoid a possible underflow.)
                    v[(j, i)] = (a[(i, j)] / a[(i, l)]) / g;
                }
                for j in l..=n {
                    s = 0.0;
                    for k in l..=n {
                        s += a[(i, k)] * v[(k, j)];
                    }
                    for k in l..=n {
                        v[(k, j)] += s * v[(k, i)];
                    }
                }
            }
            for j in l..=n {
                v[(i, j)] = 0.0;
                v[(j, i)] = 0.0;
            }
        }
        v[(i, i)] = 1.0;
        g = rv1[i];
        l = i;
    }
    // Accumulation of the left-hand transformations.
    for i in (1..=n).rev() {
        l = i + 1;
        g = w[i];
        if i < n {
            for j in l..=n {
                a[(i, j)] = 0.0;
            }
        }
        if g != 0.0 {
            g = 1.0 / g;
            if i != n {
                for j in l..=n {
                    s = 0.0;
                    for k in l..=m {
                        s += a[(k, i)] * a[(k, j)];
                    }
                    f = (s / a[(i, i)]) * g;
                    for k in i..=m {
                        a[(k, j)] += f * a[(k, i)];
                    }
                }
            }
            for j in i..=m {
                a[(j, i)] *= g;
            }
        } else {
            for j in i..=m {
                a[(j, i)] = 0.0;
            }
        }
        a[(i, i)] += 1.0;
    }
    // Diagonalization of the bidiagonal form.
    for k in (1..=n).rev() {
        for its in 1..=30 {
            let mut flag = true;
            // Test for splitting (rv1[1] is always zero, so the loop ends by l = 1 and
            // w[nm] is never read at nm = 0).
            l = k;
            loop {
                nm = l - 1;
                if rv1[l].abs() + anorm == anorm {
                    flag = false;
                    break;
                }
                if w[nm].abs() + anorm == anorm {
                    break;
                }
                if l == 1 {
                    // (Only with a NaN anorm does the loop get past l = 1; OCCT then reads
                    // out of range. Reported as not converged.)
                    return false;
                }
                l -= 1;
            }
            if flag {
                // Cancellation of rv1[l] if l > 1.
                s = 1.0;
                for i in l..=k {
                    f = s * rv1[i];
                    if f.abs() + anorm != anorm {
                        g = w[i];
                        h = pythag(f, g);
                        w[i] = h;
                        h = 1.0 / h;
                        let c = g * h;
                        s = -f * h;
                        for j in 1..=m {
                            y = a[(j, nm)];
                            z = a[(j, i)];
                            a[(j, nm)] = y * c + z * s;
                            a[(j, i)] = z * c - y * s;
                        }
                    }
                }
            }
            z = w[k];
            if l == k {
                // Convergence: the singular value made nonnegative.
                if z < 0.0 {
                    w[k] = -z;
                    for j in 1..=n {
                        v[(j, k)] = -v[(j, k)];
                    }
                }
                break;
            }
            if its == 30 {
                return false;
            }
            // Shift from the bottom 2 x 2 minor.
            x = w[l];
            nm = k - 1;
            y = w[nm];
            g = rv1[nm];
            h = rv1[k];
            f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
            g = pythag(f, 1.0);
            f = ((x - z) * (x + z) + h * ((y / (f + sign(g, f))) - h)) / x;
            // The next QR transformation.
            let mut c = 1.0;
            s = 1.0;
            for j in l..=nm {
                let i = j + 1;
                g = rv1[i];
                y = w[i];
                h = s * g;
                g *= c;
                z = pythag(f, h);
                rv1[j] = z;
                c = f / z;
                s = h / z;
                f = x * c + g * s;
                g = g * c - x * s;
                h = y * s;
                y *= c;
                for jj in 1..=n {
                    x = v[(jj, j)];
                    z = v[(jj, i)];
                    v[(jj, j)] = x * c + z * s;
                    v[(jj, i)] = z * c - x * s;
                }
                z = pythag(f, h);
                w[j] = z;
                if z != 0.0 {
                    z = 1.0 / z;
                    c = f * z;
                    s = h * z;
                }
                f = (c * g) + (s * y);
                x = (c * y) - (s * g);
                for jj in 1..=m {
                    y = a[(jj, j)];
                    z = a[(jj, i)];
                    a[(jj, j)] = y * c + z * s;
                    a[(jj, i)] = z * c - y * s;
                }
            }
            rv1[l] = 0.0;
            rv1[k] = f;
            w[k] = x;
        }
    }
    true
}

/// `SVD_Solve(u, w, v, b, x)` of math_Recipes: x = V diag(1/w) U^T b, the zero singular
/// values left out; `w` from 1, `b` (`u.rows` entries) and `x` (`u.cols`) from 0.
#[doc(alias = "SVD_Solve")]
fn svd_solve(u: &Matrix, w: &[f64], v: &Matrix, b: &[f64], x: &mut [f64]) {
    let m = u.rows;
    let n = u.cols;
    let mut tmp = vec![0.0; n + 1];
    for j in 1..=n {
        let mut s = 0.0;
        if w[j] != 0.0 {
            for i in 1..=m {
                s += u[(i, j)] * b[i - 1];
            }
            s /= w[j];
        }
        tmp[j] = s;
    }
    for j in 1..=n {
        let mut s = 0.0;
        for jj in 1..=n {
            s += v[(j, jj)] * tmp[jj];
        }
        x[j - 1] = s;
    }
}

/// `math_Gauss(A, MinPivot)`: the LU decomposition of a square matrix.
#[doc(alias = "math_Gauss")]
struct Gauss {
    lu: Matrix,
    index: Vec<usize>,
    done: bool,
}

impl Gauss {
    /// `math_Gauss(A, MinPivot)`: `a` row-major, `n` x `n`.
    fn new(a: &[f64], n: usize, min_pivot: f64) -> Self {
        let mut lu = Matrix::from_row_major(n, n, a);
        let mut index = vec![0; n];
        let mut d = 0.0;
        let done = lu_decompose(&mut lu, &mut index, &mut d, min_pivot);
        Gauss { lu, index, done }
    }

    /// `math_Gauss::Solve(X)`: `x` holds the right-hand side and receives the solution.
    fn solve(&self, x: &mut [f64]) {
        debug_assert!(self.done);
        lu_solve(&self.lu, &self.index, x);
    }
}

/// `math_SVD(A)`: the singular value decomposition of a matrix.
#[doc(alias = "math_SVD")]
struct Svd {
    /// `U`: max(rows, cols) x cols, A padded with zero rows.
    u: Matrix,
    v: Matrix,
    /// `Diag`, from 1.
    diag: Vec<f64>,
    /// `RowA`.
    row_a: usize,
    done: bool,
}

impl Svd {
    /// `math_SVD(A)`: `a` row-major, `rows` x `cols`.
    fn new(a: &[f64], rows: usize, cols: usize) -> Self {
        let mut u = Matrix::zeros(rows.max(cols), cols);
        u.data[..rows * cols].copy_from_slice(&a[..rows * cols]);
        let mut v = Matrix::zeros(cols, cols);
        let mut diag = vec![0.0; cols + 1];
        let done = svd_decompose(&mut u, &mut diag, &mut v);
        Svd { u, v, diag, row_a: rows, done }
    }

    /// `math_SVD::Solve(B, X, Eps = 1.0e-6)`: the singular values below `eps` times the
    /// largest are zeroed (for good: `Diag` is changed), then `SVD_Solve`.
    fn solve(&mut self, b: &[f64], x: &mut [f64]) {
        let eps = 1.0e-6;
        debug_assert!(self.done && b.len() == self.row_a);
        let mut bb = vec![0.0; self.u.rows];
        bb[..b.len()].copy_from_slice(b);
        // `Diag(Diag.Max())`: the largest value (the first of equal ones).
        let mut largest = self.diag[1];
        for &value in &self.diag[2..] {
            if value > largest {
                largest = value;
            }
        }
        let wmin = eps * largest;
        for value in &mut self.diag[1..] {
            if *value < wmin {
                *value = 0.0;
            }
        }
        svd_solve(&self.u, &self.diag, &self.v, &bb, x);
    }
}

/// `math_GaussLeastSquare(A, MinPivot = 1.0e-20)`: the LU decomposition of A^T A, for the
/// least-squares solution of an overdetermined system.
#[doc(alias = "math_GaussLeastSquare")]
struct GaussLeastSquare {
    lu: Matrix,
    /// `A2`: A^T.
    a2: Matrix,
    index: Vec<usize>,
    done: bool,
}

impl GaussLeastSquare {
    /// `math_GaussLeastSquare(A)`: `a` row-major, `rows` x `cols`.
    fn new(a: &[f64], rows: usize, cols: usize) -> Self {
        let min_pivot = 1.0e-20;
        let a = Matrix::from_row_major(rows, cols, a);
        let mut a2 = Matrix::zeros(cols, rows);
        for i in 1..=rows {
            for j in 1..=cols {
                a2[(j, i)] = a[(i, j)];
            }
        }
        // LU.Multiply(A2, A).
        let mut lu = Matrix::zeros(cols, cols);
        for i in 1..=cols {
            for j in 1..=cols {
                let mut som = 0.0;
                for k in 1..=rows {
                    som += a2[(i, k)] * a[(k, j)];
                }
                lu[(i, j)] = som;
            }
        }
        let mut index = vec![0; cols];
        let mut d = 0.0;
        let done = lu_decompose(&mut lu, &mut index, &mut d, min_pivot);
        GaussLeastSquare { lu, a2, index, done }
    }

    /// `math_GaussLeastSquare::Solve(B, X)`: X = A^T B, then `LU_Solve`.
    fn solve(&self, b: &[f64], x: &mut [f64]) {
        debug_assert!(self.done);
        for i in 1..=self.a2.rows {
            let mut som = 0.0;
            for k in 1..=self.a2.cols {
                som += self.a2[(i, k)] * b[k - 1];
            }
            x[i - 1] = som;
        }
        lu_solve(&self.lu, &self.index, x);
    }
}

/// `MyDirFunction`: the system along the line `p0 + x * dir`, as half its square norm.
/// `p` and `fv` are OCCT's scratch vectors (`Temp3`, `Temp4`).
struct DirFunction {
    p0: Vec<f64>,
    dir: Vec<f64>,
    p: Vec<f64>,
    fv: Vec<f64>,
}

impl DirFunction {
    fn new(nb_variables: usize, nb_equations: usize) -> Self {
        DirFunction { p0: vec![0.0; nb_variables], dir: vec![0.0; nb_variables], p: vec![0.0; nb_variables], fv: vec![0.0; nb_equations] }
    }

    /// `MyDirFunction::Initialize(p0, dir)`.
    fn initialize(&mut self, p0: &[f64], dir: &[f64]) {
        self.p0.copy_from_slice(p0);
        self.dir.copy_from_slice(dir);
    }

    /// `MyDirFunction::Value(x, fval)`.
    fn value(&mut self, function: &mut dyn FunctionSetWithDerivatives, x: f64) -> Option<f64> {
        for i in 0..self.p.len() {
            self.p[i] = self.dir[i] * x + self.p0[i];
        }
        if !function.value(&self.p, &mut self.fv) {
            return None;
        }
        for &value in &self.fv {
            // (Precision::HalfInfinite later, says OCCT.)
            if value <= -1.0e100 || value >= 1.0e100 {
                return None;
            }
        }
        Some(0.5 * norm2(&self.fv))
    }
}

/// `math_Vector::Norm2()`: the square norm.
fn norm2(v: &[f64]) -> f64 {
    v.iter().map(|x| x * x).sum()
}

/// `math_Vector::Norm()`.
fn norm(v: &[f64]) -> f64 {
    norm2(v).sqrt()
}

/// `u * v` of math_Vector: the dot product.
fn dot(u: &[f64], v: &[f64]) -> f64 {
    u.iter().zip(v).map(|(a, b)| a * b).sum()
}

/// `MyDirFunction::Value(Sol, FF, DF, GH, F2, Gnr1)`: the system and its Jacobian at `sol`,
/// the gradient `GH = DF^T FF` of half its square norm, and returns `(F2, Gnr1)`: that half
/// square norm and the gradient's square norm.
fn evaluate(function: &mut dyn FunctionSetWithDerivatives, sol: &[f64], ff: &mut [f64], df: &mut [f64], gh: &mut [f64]) -> Option<(f64, f64)> {
    if !function.values(sol, ff, df) {
        return None;
    }
    for &value in ff.iter() {
        if value < 0.0 {
            if value <= -1.0e100 {
                return None;
            }
        } else if value >= 1.0e100 {
            return None;
        }
    }
    let f2 = 0.5 * norm2(ff);
    // GH.TMultiply(DF, FF).
    let ninc = gh.len();
    for j in 0..ninc {
        let mut som = 0.0;
        for i in 0..ff.len() {
            som += df[i * ninc + j] * ff[i];
        }
        gh[j] = som;
    }
    if gh.iter().any(|&value| is_infinite(value)) {
        return None;
    }
    Some((f2, norm2(gh)))
}

/// `MinimizeDirection(P0, P1, P2, F1, Delta, Tol, F)`: a Brent search along `p1 - p0`
/// through `p1` for a value below `f1`; on success `delta` becomes the step from `p1`.
#[allow(clippy::too_many_arguments)]
fn minimize_direction_brent(function: &mut dyn FunctionSetWithDerivatives, f_dir: &mut DirFunction, p0: &[f64], p1: &[f64], p2: &[f64], f1: f64, delta: &mut [f64], tol: &[f64]) -> bool {
    // (1) A one-dimensional parametric tolerance.
    let mut tol1d: f64 = 2.1;
    let eps = 1.0e-16;
    for ii in 0..tol.len() {
        let invnorme = delta[ii].abs();
        if invnorme > eps {
            tol1d = tol1d.min(tol[ii] / invnorme);
        }
    }
    if tol1d > 1.9 {
        // Not worth the trouble.
        return false;
    }
    tol1d /= 3.0;

    for i in 0..delta.len() {
        delta[i] = p1[i] - p0[i];
    }
    let mut invnorme = norm(delta);
    if invnorme <= eps {
        return false;
    }
    invnorme = 1.0 / invnorme;

    f_dir.initialize(p1, delta);

    // (2) Minimize.
    let ax = -1.0;
    let bx = 0.0;
    let p2_p1: Vec<f64> = p2.iter().zip(p1).map(|(a, b)| a - b).collect();
    let cx = norm(&p2_p1) * invnorme;
    if cx < 1.0e-2 {
        return false;
    }
    let mut solver = BrentMinimum::new(tol1d, 100, tol1d);
    solver.perform(&mut |x| f_dir.value(function, x), ax, bx, cx);
    if solver.is_done() {
        let tsol = solver.location();
        if solver.minimum() < f1 {
            for value in delta.iter_mut() {
                *value *= tsol;
            }
            return true;
        }
    }
    false
}

/// `MinimizeDirection(P, Dir, PValue, PDirValue, Gradient, DGradient, Tol, F)`: a quadratic
/// interpolation along `dir` from `p`, then a Brent search if it gained too little; on
/// success `dir` is scaled to the step taken.
#[allow(clippy::too_many_arguments)]
fn minimize_direction_interpolated(function: &mut dyn FunctionSetWithDerivatives, f_dir: &mut DirFunction, p: &[f64], dir: &mut [f64], p_value: f64, p_dir_value: f64, gradient: &[f64], d_gradient: &[f64], tol: &[f64]) -> bool {
    if is_infinite(p_value) || is_infinite(p_dir_value) {
        return false;
    }
    // (0) A one-dimensional parametric tolerance.
    let mut good = false;
    let eps = 1.0e-20;
    let mut tol1d: f64 = 1.1;
    let mut result = p_value;
    for ii in 0..tol.len() {
        let absdir = dir[ii].abs();
        if absdir > eps {
            tol1d = tol1d.min(tol[ii] / absdir);
        }
    }
    if tol1d > 0.9 {
        return false;
    }

    // (1) A first quadratic interpolation.
    let (mut ax, mut bx, mut cx);
    let mut tsol;
    let df1 = dot(gradient, dir);
    let df2 = dot(d_gradient, dir);
    if df1 < -eps && df2 > eps {
        // A basin.
        tsol = -df1 / (df2 - df1);
    } else {
        cx = p_value;
        bx = df1;
        ax = p_dir_value - (bx + cx);
        if ax.abs() <= eps {
            // Linear.
            tsol = if bx.abs() >= eps { -cx / bx } else { 0.0 };
        } else {
            // Quadratic.
            let mut delta = bx * bx - 4.0 * ax * cx;
            if delta > 1.0e-9 {
                // Roots: take the one nearest 0.
                delta = delta.sqrt();
                tsol = -(bx + delta);
                let tsolbis = delta - bx;
                if tsolbis.abs() < tsol.abs() {
                    tsol = tsolbis;
                }
                tsol /= 2.0 * ax;
            } else {
                // No root or hardly any: take the extremum.
                tsol = -(0.5 * bx) / ax;
            }
        }
    }
    if tsol.abs() >= 1.0 {
        // A result of no interest.
        return false;
    }

    f_dir.initialize(p, dir);
    // (OCCT ignores whether this evaluation succeeds; its value is then uninitialized. The
    // port reads it as infinite, which keeps the point from being taken.)
    let fsol = f_dir.value(function, tsol).unwrap_or(f64::INFINITY);
    if fsol < p_value {
        good = true;
        result = fsol;
    }

    // (2) Too little progress: a proper search from the same start.
    if fsol > 0.2 * p_value && tol1d < 0.5 {
        if tsol < 0.0 {
            ax = tsol;
            bx = 0.0;
            cx = 1.0;
        } else {
            ax = 0.0;
            bx = tsol;
            cx = 1.0;
        }
        let mut solver = BrentMinimum::new(tol1d, 100, tol1d);
        solver.perform(&mut |x| f_dir.value(function, x), ax, bx, cx);
        if solver.is_done() && solver.minimum() <= result {
            tsol = solver.location();
            good = true;
            result = solver.minimum();
            // The objective changes too fast: search both halves again (the solution being
            // inside (ax, cx)).
            if norm2(gradient) > 1.0 / SQUARE_CONFUSION && tsol > ax && tsol < cx {
                solver.perform(&mut |x| f_dir.value(function, x), ax, (ax + tsol) / 2.0, tsol);
                if solver.is_done() && solver.minimum() <= result {
                    tsol = solver.location();
                    good = true;
                    result = solver.minimum();
                }
                solver.perform(&mut |x| f_dir.value(function, x), tsol, (cx + tsol) / 2.0, cx);
                if solver.is_done() && solver.minimum() <= result {
                    tsol = solver.location();
                    good = true;
                }
            }
        }
    }
    if good {
        // Update Delta.
        for value in dir.iter_mut() {
            *value *= tsol;
        }
    }
    good
}

/// `SearchDirection(DF, GH, FF, ChangeDirection, InvLengthMax, Direction, Dy)` (the
/// unconstrained one): Newton's step (`math_Gauss`, or `math_SVD` where the matrix is
/// singular, wide or tall as it is), shortened to a quarter of the domain, or the steepest
/// descent where Newton's does not descend. `df` is `neq` x `direction.len()` row-major;
/// `change_direction` is taken by value, as in OCCT. Returns `Dy`.
#[allow(clippy::too_many_arguments)]
fn search_direction(df: &[f64], neq: usize, gh: &[f64], ff: &[f64], mut change_direction: bool, inv_length_max: &[f64], direction: &mut [f64]) -> f64 {
    let ninc = direction.len();
    let eps = 1.0e-32;
    let minus_ff: Vec<f64> = ff.iter().map(|value| -value).collect();
    if !change_direction {
        if ninc == neq {
            direction.copy_from_slice(&minus_ff);
            let solut = Gauss::new(df, ninc, 1.0e-9);
            if solut.done {
                solut.solve(direction);
            } else {
                // We have to "forget" singular directions.
                let mut solve_by_svd = Svd::new(df, neq, ninc);
                if solve_by_svd.done {
                    solve_by_svd.solve(&minus_ff, direction);
                } else {
                    change_direction = true;
                }
            }
        } else if ninc > neq {
            let mut solut = Svd::new(df, neq, ninc);
            if solut.done {
                solut.solve(&minus_ff, direction);
            } else {
                change_direction = true;
            }
        } else {
            // By GaussLeastSquare.
            let solut = GaussLeastSquare::new(df, neq, ninc);
            if solut.done {
                solut.solve(&minus_ff, direction);
            } else {
                change_direction = true;
            }
        }
    }
    // Directions too long are better forbidden, to guard the badly conditioned cases (PMN
    // 12/05/97: singularities in fillets on periodic surfaces).
    let mut ratio = (direction[0] * inv_length_max[0]).abs();
    for i in 1..ninc {
        ratio = ratio.max((direction[i] * inv_length_max[i]).abs());
    }
    if ratio > 1.0 {
        for value in direction.iter_mut() {
            *value /= ratio;
        }
    }

    let mut dy = dot(direction, gh);
    if dy >= -eps {
        // Newton "does not descend": take the gradient.
        change_direction = true;
    }
    if change_direction {
        for i in 0..ninc {
            direction[i] = -gh[i];
        }
        dy = -norm2(gh);
    }
    dy
}

/// The constrained `SearchDirection(DF, GH, FF, Constraints, X, ChangeDirection,
/// InvLengthMax, Direction, Dy)`: a direction (and a step, where Newton works) along a
/// boundary, solving the sub-problem of the unknowns not at a bound. Returns `Dy`.
#[allow(clippy::too_many_arguments)]
fn search_direction_constrained(df: &[f64], neq: usize, gh: &[f64], ff: &[f64], constraints: &[bool], change_direction: bool, inv_length_max: &[f64], direction: &mut [f64]) -> f64 {
    let ninc = direction.len();
    // Check the imposed bounds (otherwise the system to solve does not change).
    let cons = constraints.iter().filter(|&&c| c).count();
    if cons == 0 {
        search_direction(df, neq, gh, ff, change_direction, inv_length_max, direction)
    } else if cons == ninc {
        // Nothing left to do.
        for value in direction.iter_mut() {
            *value = 0.0;
        }
        0.0
    } else {
        // (1) The general case: define a sub-problem.
        let free = ninc - cons;
        let mut df2 = vec![0.0; neq * free];
        let mut my_gh = vec![0.0; free];
        let mut my_direction = vec![0.0; free];
        let mut my_inv_length_max = vec![0.0; free];
        let mut k = 0;
        for i in 0..ninc {
            if !constraints[i] {
                my_gh[k] = gh[i];
                my_inv_length_max[k] = inv_length_max[i];
                my_direction[k] = direction[i];
                for j in 0..neq {
                    df2[j * free + k] = df[j * ninc + i];
                }
                // On to the next unknown.
                k += 1;
            }
        }
        // (2) Solve it.
        let dy = search_direction(&df2, neq, &my_gh, ff, change_direction, &my_inv_length_max, &mut my_direction);
        // (3) Interpret it: rebuild Direction.
        let mut k = 0;
        for i in 0..ninc {
            if !constraints[i] {
                if !change_direction {
                    direction[i] = my_direction[k];
                } else {
                    direction[i] = -gh[i];
                }
                k += 1;
            } else {
                direction[i] = 0.0;
            }
        }
        dy
    }
}

/// `Bounds(InfBound, SupBound, Tol, Sol, SolSave, Constraints, Delta, isNewSol)`: trims
/// `sol` into `[inf, sup]` along the step from `sol_save`, `delta` getting the trimmed step
/// and `constraints` the unknowns at a bound. Returns `(Out, isNewSol)`: whether the step
/// left the domain, and false where `sol` stays `sol_save` (it was on a bound and the step
/// went straight beyond it).
#[allow(clippy::too_many_arguments)]
fn bounds(inf: &[f64], sup: &[f64], tol: &[f64], sol: &mut [f64], sol_save: &[f64], constraints: &mut [bool], delta: &mut [f64]) -> (bool, bool) {
    let mut out = false;
    let ninc = sol.len();
    let mut monratio: f64 = 1.0;
    let mut is_new_sol = true;

    // The ratio to bring the step back by.
    for i in 0..ninc {
        constraints[i] = false;
        delta[i] = sol[i] - sol_save[i];
        if inf[i] == sup[i] {
            constraints[i] = true;
            // ("Ok, but this should be avoided.")
            out = true;
        } else if sol[i] < inf[i] {
            constraints[i] = true;
            out = true;
            // (Delta is negative.) Avoid zero ratios for nothing.
            if -delta[i] > tol[i] {
                monratio = monratio.min((inf[i] - sol_save[i]) / delta[i]);
            }
        } else if sol[i] > sup[i] {
            constraints[i] = true;
            out = true;
            // (Delta is positive.)
            if delta[i] > tol[i] {
                monratio = monratio.min((sup[i] - sol_save[i]) / delta[i]);
            }
        }
    }

    if out {
        // Truncate, and clamp once more against numerical trouble.
        if monratio == 0.0 {
            is_new_sol = false;
            sol.copy_from_slice(sol_save);
            for value in delta.iter_mut() {
                *value = 0.0;
            }
        } else {
            for i in 0..ninc {
                delta[i] *= monratio;
                sol[i] = sol_save[i] + delta[i];
            }
            for i in 0..ninc {
                if sol[i] < inf[i] {
                    sol[i] = inf[i];
                    delta[i] = sol[i] - sol_save[i];
                } else if sol[i] > sup[i] {
                    sol[i] = sup[i];
                    delta[i] = sol[i] - sol_save[i];
                }
            }
        }
    }
    (out, is_new_sol)
}

/// `math_FunctionSetRoot`: the root of a system of equations near a starting point,
/// between bounds, to within a tolerance on each unknown.
#[doc(alias = "math_FunctionSetRoot")]
pub(crate) struct FunctionSetRoot {
    /// `Delta`: the last step (`FunctionSetErrors()`).
    delta: Vec<f64>,
    /// `Sol`.
    sol: Vec<f64>,
    /// `DF`: the Jacobian, `NbEquations()` x `NbVariables()` row-major.
    df: Vec<f64>,
    /// `Tol`.
    tol: Vec<f64>,
    done: bool,
    /// `Kount`: the iteration reached.
    kount: usize,
    /// `Itermax`.
    iteration_limit: usize,
    sol_save: Vec<f64>,
    gh: Vec<f64>,
    dh: Vec<f64>,
    dh_save: Vec<f64>,
    ff: Vec<f64>,
    previous_solution: Vec<f64>,
    /// `Save(0..Itermax)`: the half square norm of the system after each iteration.
    save: Vec<f64>,
}

/// `Sol = A + B` of math_Vector, into `sol`.
fn set_sum(sol: &mut [f64], a: &[f64], b: &[f64]) {
    for i in 0..sol.len() {
        sol[i] = a[i] + b[i];
    }
}

/// `Sol(i) = Sol(i) + Ambda * DH(i)`.
fn step(sol: &mut [f64], ambda: f64, dh: &[f64]) {
    for i in 0..sol.len() {
        sol[i] += ambda * dh[i];
    }
}

impl FunctionSetRoot {
    /// `math_FunctionSetRoot(F, NbIterations = 100)` then `SetTolerance(tolerance)`.
    /// (`n` unknowns and as many equations; [`perform`](Self::perform) resizes for a system
    /// with another number of equations.)
    pub(crate) fn new(n: usize, iteration_limit: usize) -> Self {
        FunctionSetRoot {
            delta: vec![0.0; n],
            sol: vec![0.0; n],
            df: vec![0.0; n * n],
            tol: vec![0.0; n],
            done: false,
            kount: 0,
            iteration_limit,
            sol_save: vec![0.0; n],
            gh: vec![0.0; n],
            dh: vec![0.0; n],
            dh_save: vec![0.0; n],
            ff: vec![0.0; n],
            previous_solution: vec![0.0; n],
            save: vec![0.0; iteration_limit + 1],
        }
    }

    /// `SetTolerance(Tolerance)`: the tolerance on each unknown.
    pub(crate) fn set_tolerance(&mut self, tolerance: &[f64]) {
        let n = self.tol.len();
        self.tol.copy_from_slice(&tolerance[..n]);
    }

    /// `IsDone()`.
    pub(crate) fn is_done(&self) -> bool {
        self.done
    }

    /// `Root()`: the solution (OCCT raises `StdFail_NotDone` where not [`is_done`](Self::is_done);
    /// here it is the last point reached).
    pub(crate) fn root(&self) -> &[f64] {
        &self.sol
    }

    /// `NbIterations()`.
    pub(crate) fn nb_iterations(&self) -> usize {
        self.kount
    }

    /// `IsSolutionReached`: OCCT's default test (subclasses may redefine it), the last step
    /// within the tolerance on every unknown: `|Delta(i)| <= Tol(i)`.
    #[doc(alias = "math_FunctionSetRoot::IsSolutionReached")]
    fn is_solution_reached(&self) -> bool {
        self.delta.iter().zip(&self.tol).all(|(delta, tol)| delta.abs() <= *tol)
    }

    /// `MyDirFunction::Value(Sol, FF, DF, GH, F2, Gnr1)` at the current `Sol`.
    fn evaluate_at_sol(&mut self, f: &mut dyn FunctionSetWithDerivatives) -> Option<(f64, f64)> {
        evaluate(f, &self.sol, &mut self.ff, &mut self.df, &mut self.gh)
    }

    /// The returns where OCCT sets `Done = true`: `F.Value(Sol, FF)` ("update F before
    /// GetStateNumber", jgv 31.08.2011), then `GetStateNumber()`.
    fn finish_done(&mut self, f: &mut dyn FunctionSetWithDerivatives) {
        self.done = true;
        f.value(&self.sol, &mut self.ff);
        f.get_state_number();
    }

    /// `Perform(F, StartingPoint, InfBound, SupBound)` (theStopOnDivergent = false).
    #[doc(alias = "math_FunctionSetRoot::Perform")]
    // (Some assignments are overwritten before being read, as in OCCT.)
    #[allow(unused_assignments)]
    pub(crate) fn perform(&mut self, f: &mut dyn FunctionSetWithDerivatives, start: &[f64], inf: &[f64], sup: &[f64]) {
        let ninc = f.nb_variables();
        let neq = f.nb_equations();
        // (Standard_DimensionError.)
        assert!(neq > 0 && start.len() == ninc && inf.len() == ninc && sup.len() == ninc && self.sol.len() == ninc);
        self.ff.resize(neq, 0.0);
        self.df.resize(neq * ninc, 0.0);

        let mut change_direction = false;
        let mut sort;
        let mut is_new_sol;
        let mut good;
        let mut verif;
        let mut stop;
        const EPS_SQRT: f64 = 1.0e-16;
        const EPS: f64 = 1.0e-32;
        const EPS2: f64 = 1.0e-64;
        const PROGRES: f64 = 0.005;
        let (mut f2, mut gnr1);
        let mut previous_minimum;
        let mut dy;
        let mut old_f;
        let mut ambda;
        let mut ambda2;
        let mut oldgr;
        // To keep steps to a quarter of the domain.
        let mut inv_length_max = vec![0.0; ninc];
        for i in 0..ninc {
            let a_sup_bound = sup[i].min(INFINITE);
            let an_inf_bound = inf[i].max(-INFINITE);
            inv_length_max[i] = 1.0 / ((a_sup_bound - an_inf_bound) / 4.0).max(1.0e-9);
        }
        // `aConstraints`: which bounds the unknowns are on.
        let mut constraints = vec![false; ninc];
        let mut f_dir = DirFunction::new(ninc, neq);
        let mut descente_iter;

        self.done = false;
        self.sol.copy_from_slice(start);
        self.kount = 0;

        // (`myIsDivergent`, whether the start or a step left the bounds, is read only in the
        // theStopOnDivergent mode.)

        // Check the unknowns against the bounds; bring them onto the bounds if not valid.
        for i in 0..ninc {
            if self.sol[i] <= inf[i] {
                self.sol[i] = inf[i];
            } else if self.sol[i] > sup[i] {
                self.sol[i] = sup[i];
            }
        }

        // The first value of F and its gradient.
        match self.evaluate_at_sol(f) {
            Some(e) => (f2, gnr1) = e,
            None => {
                self.done = false;
                f.get_state_number();
                return;
            }
        }
        ambda2 = gnr1;
        // Save(0) is used only by the accelerating test at the end of the loop; on a
        // solution already, it must not force a second iteration.
        self.save[0] = f2.max(EPS_SQRT);
        let a_tol_func = epsilon(f2);

        if f2 <= EPS || gnr1 <= EPS2 {
            // (No F.Value update here in OCCT.)
            self.done = true;
            f.get_state_number();
            return;
        }

        self.kount = 1;
        while self.kount <= self.iteration_limit {
            let kount = self.kount;
            previous_minimum = f2;
            oldgr = gnr1;
            self.previous_solution.copy_from_slice(&self.sol);
            self.sol_save.copy_from_slice(&self.sol);

            dy = search_direction(&self.df, neq, &self.gh, &self.ff, change_direction, &inv_length_max, &mut self.dh);
            if dy.abs() <= EPS {
                self.finish_done(f);
                return;
            }
            if change_direction {
                ambda = ambda2 / dy.abs().sqrt();
                if ambda > 1.0 {
                    ambda = 1.0;
                }
            } else {
                ambda = 1.0;
                ambda2 = 0.5 * ambda / norm(&self.dh);
            }

            step(&mut self.sol, ambda, &self.dh);
            (sort, is_new_sol) = bounds(inf, sup, &self.tol, &mut self.sol, &self.sol_save, &mut constraints, &mut self.delta);

            self.dh_save.copy_from_slice(&self.gh);
            if is_new_sol {
                match self.evaluate_at_sol(f) {
                    Some(e) => (f2, gnr1) = e,
                    None => {
                        self.done = false;
                        f.get_state_number();
                        return;
                    }
                }
            }

            if f2 <= EPS || gnr1 <= EPS2 {
                self.finish_done(f);
                return;
            }

            if sort || f2 / previous_minimum > PROGRES {
                dy = dot(&self.gh, &self.dh);
                old_f = previous_minimum;
                stop = false;
                good = false;
                descente_iter = 0;
                let mut sortbis;

                // Standard processing, without the bounds.
                if !sort {
                    // Not out: try to go further.
                    while f2 / previous_minimum > PROGRES && !stop {
                        if f2 < old_f && dy < 0.0 {
                            // Try to progress in this direction.
                            descente_iter += 1;
                            self.sol_save.copy_from_slice(&self.sol);
                            old_f = f2;
                            step(&mut self.sol, ambda, &self.dh);
                            (stop, is_new_sol) = bounds(inf, sup, &self.tol, &mut self.sol, &self.sol_save, &mut constraints, &mut self.delta);
                            ambda *= 1.7;
                        } else {
                            if f2 >= old_f || f2 >= previous_minimum {
                                good = false;
                                if descente_iter == 0 {
                                    // The first step fails: interpolate, and minimize if need be.
                                    descente_iter += 1;
                                    good = minimize_direction_interpolated(f, &mut f_dir, &self.sol_save, &mut self.delta, old_f, f2, &self.dh_save, &self.gh, &self.tol);
                                } else if change_direction || descente_iter > 1 || old_f > previous_minimum {
                                    // The progress was useful: minimize.
                                    descente_iter += 1;
                                    good = minimize_direction_brent(f, &mut f_dir, &self.previous_solution, &self.sol_save, &self.sol, old_f, &mut self.delta, &self.tol);
                                }
                                if !good {
                                    self.sol.copy_from_slice(&self.sol_save);
                                    f2 = old_f;
                                } else {
                                    set_sum(&mut self.sol, &self.sol_save, &self.delta);
                                    (sort, is_new_sol) = bounds(inf, sup, &self.tol, &mut self.sol, &self.sol_save, &mut constraints, &mut self.delta);
                                }
                                // The point on the boundary is rejected.
                                sort = false;
                            }
                            // And leave in any case.
                            stop = true;
                        }
                        self.dh_save.copy_from_slice(&self.gh);
                        if is_new_sol {
                            match self.evaluate_at_sol(f) {
                                Some(e) => (f2, gnr1) = e,
                                None => {
                                    self.done = false;
                                    f.get_state_number();
                                    return;
                                }
                            }
                        }
                        dy = dot(&self.gh, &self.dh);
                        if dy.abs() <= EPS {
                            if f2 > old_f {
                                self.sol.copy_from_slice(&self.sol_save);
                            }
                            self.finish_done(f);
                            return;
                        }
                        if descente_iter >= 100 {
                            stop = true;
                        }
                    }
                }

                // Processing on the bounds.
                if sort {
                    // Not to progress along the boundary.
                    stop = f2 > 1.001 * old_f;
                    sortbis = sort;
                    descente_iter = 0;
                    while sortbis && (f2 < old_f || descente_iter == 0) && !stop {
                        descente_iter += 1;
                        // Try to progress along the boundary.
                        self.sol_save.copy_from_slice(&self.sol);
                        old_f = f2;
                        dy = search_direction_constrained(&self.df, neq, &self.gh, &self.ff, &constraints, change_direction, &inv_length_max, &mut self.dh);
                        if dy < -EPS {
                            // (Avoids useless computations and division by zero.)
                            if change_direction {
                                ambda = ambda2 / (-dy).sqrt();
                                if ambda > 1.0 {
                                    ambda = 1.0;
                                }
                            } else {
                                ambda = 1.0;
                                ambda2 = 0.5 * ambda / norm(&self.dh);
                            }
                            step(&mut self.sol, ambda, &self.dh);
                            (sortbis, is_new_sol) = bounds(inf, sup, &self.tol, &mut self.sol, &self.sol_save, &mut constraints, &mut self.delta);

                            self.dh_save.copy_from_slice(&self.gh);
                            if is_new_sol {
                                match self.evaluate_at_sol(f) {
                                    Some(e) => (f2, gnr1) = e,
                                    None => {
                                        self.done = false;
                                        f.get_state_number();
                                        return;
                                    }
                                }
                            }
                            ambda2 = gnr1;
                        } else {
                            stop = true;
                        }

                        while f2 / previous_minimum > PROGRES && f2 < old_f && !stop {
                            descente_iter += 1;
                            if f2 < old_f && dy < 0.0 {
                                // Try to progress in this direction.
                                self.sol_save.copy_from_slice(&self.sol);
                                old_f = f2;
                                step(&mut self.sol, ambda, &self.dh);
                                (sortbis, is_new_sol) = bounds(inf, sup, &self.tol, &mut self.sol, &self.sol_save, &mut constraints, &mut self.delta);
                            }
                            self.dh_save.copy_from_slice(&self.gh);
                            if is_new_sol {
                                match self.evaluate_at_sol(f) {
                                    Some(e) => (f2, gnr1) = e,
                                    None => {
                                        self.done = false;
                                        f.get_state_number();
                                        return;
                                    }
                                }
                            }
                            ambda2 = gnr1;
                            dy = dot(&self.gh, &self.dh);
                            stop = dy >= 0.0 || descente_iter >= 10 || sortbis;
                        }
                        stop = dy >= 0.0 || descente_iter >= 10;
                    }
                    if (f2 / previous_minimum > PROGRES && f2 >= old_f) || f2 >= previous_minimum {
                        // Minimize by Brent.
                        descente_iter += 1;
                        good = minimize_direction_interpolated(f, &mut f_dir, &self.sol_save, &mut self.delta, old_f, f2, &self.dh_save, &self.gh, &self.tol);
                        if !good {
                            self.sol.copy_from_slice(&self.sol_save);
                            sort = false;
                        } else {
                            set_sum(&mut self.sol, &self.sol_save, &self.delta);
                            (sort, is_new_sol) = bounds(inf, sup, &self.tol, &mut self.sol, &self.sol_save, &mut constraints, &mut self.delta);
                            if is_new_sol {
                                match self.evaluate_at_sol(f) {
                                    Some(e) => (f2, gnr1) = e,
                                    None => {
                                        self.done = false;
                                        f.get_state_number();
                                        return;
                                    }
                                }
                            }
                        }
                        dy = dot(&self.gh, &self.dh);
                    }
                }
            }

            // The stopping tests.
            self.save[kount] = f2;
            // Is this the solution?
            if change_direction {
                // Gradient: avoid looping.
                verif = true;
            } else {
                verif = false;
                if kount > 1 {
                    // Speeds up the nearly quadratic cases.
                    if self.save[kount - 1] < 1.0e-4 * self.save[kount - 2] {
                        verif = true;
                    }
                } else {
                    // For the cases already solved.
                    verif = f2 < 1.0e-6 * self.save[0];
                }
            }
            if verif {
                for i in 0..ninc {
                    self.delta[i] = self.previous_solution[i] - self.sol[i];
                }
                if self.is_solution_reached() {
                    if previous_minimum < f2 {
                        self.sol.copy_from_slice(&self.sol_save);
                    }
                    self.finish_done(f);
                    return;
                }
            }

            // The progress: the current minimum against the previous one.
            if f2 - previous_minimum <= a_tol_func {
                if kount > 5 {
                    // Is the history good?
                    if f2 >= 0.95 * self.save[kount - 5] {
                        if !change_direction {
                            change_direction = true;
                        } else {
                            // Less than 5% gained: leave. (No F.Value update here in OCCT.)
                            self.done = true;
                            f.get_state_number();
                            return;
                        }
                    } else {
                        // It is: start again.
                        change_direction = false;
                    }
                } else {
                    // No history: go on.
                    change_direction = false;
                }
                // If Newton does not decrease the gradient enough, try the gradient method,
                // unless f decreases (with Newton the gradient of f can grow while f
                // decreases: Newton must then be kept).
                if gnr1 > 0.9 * oldgr && f2 > 0.5 * previous_minimum {
                    change_direction = true;
                }
                // Not changing strategy: check, if not done already.
                if !change_direction && !verif {
                    for i in 0..ninc {
                        self.delta[i] = self.previous_solution[i] - self.sol[i];
                    }
                    if self.is_solution_reached() {
                        self.finish_done(f);
                        return;
                    }
                }
            } else {
                // A regression.
                if !change_direction {
                    // Switch to the gradient.
                    change_direction = true;
                    self.sol.copy_from_slice(&self.previous_solution);
                    match self.evaluate_at_sol(f) {
                        Some(e) => (f2, gnr1) = e,
                        None => {
                            self.done = false;
                            f.get_state_number();
                            return;
                        }
                    }
                } else {
                    // No way out left.
                    f.get_state_number();
                    return;
                }
            }
            self.kount += 1;
        }
        f.get_state_number();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A system given by closures for F and its Jacobian.
    struct System<V: FnMut(&[f64], &mut [f64], &mut [f64])> {
        nb_variables: usize,
        nb_equations: usize,
        values: V,
        evaluations: usize,
        states: usize,
    }

    impl<V: FnMut(&[f64], &mut [f64], &mut [f64])> FunctionSetWithDerivatives for System<V> {
        fn nb_variables(&self) -> usize {
            self.nb_variables
        }
        fn nb_equations(&self) -> usize {
            self.nb_equations
        }
        fn values(&mut self, x: &[f64], f: &mut [f64], d: &mut [f64]) -> bool {
            self.evaluations += 1;
            (self.values)(x, f, d);
            true
        }
        fn get_state_number(&mut self) {
            self.states += 1;
        }
    }

    fn system<V: FnMut(&[f64], &mut [f64], &mut [f64])>(nb_variables: usize, nb_equations: usize, values: V) -> System<V> {
        System { nb_variables, nb_equations, values, evaluations: 0, states: 0 }
    }

    fn solve(f: &mut dyn FunctionSetWithDerivatives, start: &[f64], inf: &[f64], sup: &[f64], tolerance: f64) -> FunctionSetRoot {
        let n = start.len();
        let mut solver = FunctionSetRoot::new(n, 100);
        solver.set_tolerance(&vec![tolerance; n]);
        solver.perform(f, start, inf, sup);
        solver
    }

    fn assert_close(actual: &[f64], expected: &[f64], tolerance: f64) {
        for (a, e) in actual.iter().zip(expected) {
            assert!((a - e).abs() <= tolerance, "{actual:?} against {expected:?}");
        }
    }

    const A: [f64; 9] = [4.0, 1.0, -2.0, 1.0, 3.0, 0.5, -1.0, 2.0, 5.0];

    #[test]
    fn lu_solves_a_linear_system() {
        // A x = b for x = (1, -2, 3).
        let x = [1.0, -2.0, 3.0];
        let mut b: Vec<f64> = (0..3).map(|i| (0..3).map(|j| A[i * 3 + j] * x[j]).sum()).collect();
        let gauss = Gauss::new(&A, 3, 1.0e-9);
        assert!(gauss.done);
        gauss.solve(&mut b);
        assert_close(&b, &x, 1e-14);
    }

    #[test]
    fn svd_reconstructs_its_matrix() {
        let mut svd = Svd::new(&A, 3, 3);
        assert!(svd.done);
        for i in 1..=3 {
            for j in 1..=3 {
                let entry: f64 = (1..=3).map(|k| svd.u[(i, k)] * svd.diag[k] * svd.v[(j, k)]).sum();
                assert!((entry - A[(i - 1) * 3 + j - 1]).abs() < 1e-13);
            }
        }
        assert!(svd.diag[1..].iter().all(|&w| w >= 0.0));
        // And solves the system.
        let x = [0.5, 1.5, -1.0];
        let b: Vec<f64> = (0..3).map(|i| (0..3).map(|j| A[i * 3 + j] * x[j]).sum()).collect();
        let mut solution = [0.0; 3];
        svd.solve(&b, &mut solution);
        assert_close(&solution, &x, 1e-13);
    }

    #[test]
    fn singular_matrix_goes_to_the_least_norm_solution() {
        // Rank 2: LU fails at the 1e-9 pivot; the SVD gives the least-norm solution.
        let a = [1.0, 1.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 2.0];
        assert!(!Gauss::new(&a, 3, 1.0e-9).done);
        let mut svd = Svd::new(&a, 3, 3);
        assert!(svd.done);
        let mut x = [0.0; 3];
        svd.solve(&[2.0, 2.0, 4.0], &mut x);
        assert_close(&x, &[1.0, 1.0, 2.0], 1e-14);
    }

    #[test]
    fn least_squares_of_a_tall_system() {
        // Fitting y = c0 + c1 t through (0, 1), (1, 3), (2, 5), (3, 7.5): the normal equations
        // give c1 = 2.15, c0 = 4.125 - 1.5 c1 = 0.9.
        let a = [1.0, 0.0, 1.0, 1.0, 1.0, 2.0, 1.0, 3.0];
        let solver = GaussLeastSquare::new(&a, 4, 2);
        assert!(solver.done);
        let mut x = [0.0; 2];
        solver.solve(&[1.0, 3.0, 5.0, 7.5], &mut x);
        assert_close(&x, &[0.9, 2.15], 1e-13);
    }

    #[test]
    fn linear_system() {
        let x = [1.0, -2.0, 3.0];
        let b: Vec<f64> = (0..3).map(|i| (0..3).map(|j| A[i * 3 + j] * x[j]).sum()).collect();
        let mut f = system(3, 3, |p: &[f64], f: &mut [f64], d: &mut [f64]| {
            for i in 0..3 {
                f[i] = (0..3).map(|j| A[i * 3 + j] * p[j]).sum::<f64>() - b[i];
            }
            d.copy_from_slice(&A);
        });
        let solver = solve(&mut f, &[0.0; 3], &[-10.0; 3], &[10.0; 3], 1e-10);
        assert!(solver.is_done());
        assert_close(solver.root(), &x, 1e-10);
        assert!(solver.nb_iterations() <= 2);
        assert!(f.states >= 1);
    }

    /// A sphere of radius `r` about the origin: S(u, v) = r (cos u cos v, sin u cos v, sin v).
    fn sphere(r: f64, u: f64, v: f64) -> [f64; 3] {
        [r * u.cos() * v.cos(), r * u.sin() * v.cos(), r * v.sin()]
    }

    #[test]
    fn line_meets_sphere() {
        // F(u, v, w) = S(u, v) - L(w), the line through S(0.3, 0.2) at w = 1.5.
        let r = 2.0;
        let (u0, v0, w0) = (0.3, 0.2, 1.5);
        let dir = [-0.6, 0.3, 0.2];
        let hit = sphere(r, u0, v0);
        let origin = [hit[0] - w0 * dir[0], hit[1] - w0 * dir[1], hit[2] - w0 * dir[2]];
        let mut f = system(3, 3, |p: &[f64], f: &mut [f64], d: &mut [f64]| {
            let (u, v, w) = (p[0], p[1], p[2]);
            let s = sphere(r, u, v);
            for i in 0..3 {
                f[i] = s[i] - (origin[i] + w * dir[i]);
            }
            let su = [-r * u.sin() * v.cos(), r * u.cos() * v.cos(), 0.0];
            let sv = [-r * u.cos() * v.sin(), -r * u.sin() * v.sin(), r * v.cos()];
            for i in 0..3 {
                d[i * 3] = su[i];
                d[i * 3 + 1] = sv[i];
                d[i * 3 + 2] = -dir[i];
            }
        });
        let inf = [-std::f64::consts::PI, -std::f64::consts::FRAC_PI_2, -10.0];
        let sup = [std::f64::consts::PI, std::f64::consts::FRAC_PI_2, 10.0];
        let solver = solve(&mut f, &[0.5, 0.0, 1.0], &inf, &sup, 1e-10);
        assert!(solver.is_done());
        assert_close(solver.root(), &[u0, v0, w0], 1e-9);
    }

    #[test]
    fn root_on_a_bound_from_outside() {
        // F = (x0^2 - 1, x1 - x0 / 2, x2 + x0 x1 - 1): the root (1, 1/2, 1/2) has x0 on its
        // upper bound; the start, beyond it, is brought back onto it.
        let mut f = system(3, 3, |p: &[f64], f: &mut [f64], d: &mut [f64]| {
            f[0] = p[0] * p[0] - 1.0;
            f[1] = p[1] - 0.5 * p[0];
            f[2] = p[2] + p[0] * p[1] - 1.0;
            d.copy_from_slice(&[2.0 * p[0], 0.0, 0.0, -0.5, 1.0, 0.0, p[1], p[0], 1.0]);
        });
        let solver = solve(&mut f, &[3.0, -1.0, 0.0], &[0.0, -2.0, -2.0], &[1.0, 2.0, 2.0], 1e-10);
        assert!(solver.is_done());
        assert_close(solver.root(), &[1.0, 0.5, 0.5], 1e-9);
    }

    #[test]
    fn least_squares_along_a_bound() {
        // F = x - c with c outside the box [0, 1]^3: the minimum of |F| over the box is c
        // clamped, (1, 0.5, 0.25), reached by sliding along the face x0 = 1 (the
        // constrained sub-problem, solved by math_GaussLeastSquare).
        let c = [2.0, 0.5, 0.25];
        let mut f = system(3, 3, |p: &[f64], f: &mut [f64], d: &mut [f64]| {
            for i in 0..3 {
                f[i] = p[i] - c[i];
            }
            d.copy_from_slice(&[1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]);
        });
        let solver = solve(&mut f, &[0.5, 0.9, 0.9], &[0.0; 3], &[1.0; 3], 1e-10);
        assert_close(solver.root(), &[1.0, 0.5, 0.25], 1e-9);
    }

    #[test]
    fn singular_jacobian_at_the_start() {
        // F = (x0^2 + x1 - 2, x0 + x1 - 2, x2 - 0.5): roots (0, 2, 0.5) and (1, 1, 0.5). The
        // Jacobian is singular where x0 = 1/2, as at the start: the first direction is the
        // SVD's.
        let mut f = system(3, 3, |p: &[f64], f: &mut [f64], d: &mut [f64]| {
            f[0] = p[0] * p[0] + p[1] - 2.0;
            f[1] = p[0] + p[1] - 2.0;
            f[2] = p[2] - 0.5;
            d.copy_from_slice(&[2.0 * p[0], 1.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0]);
        });
        let solver = solve(&mut f, &[0.5, 0.3, 0.0], &[-5.0; 3], &[5.0; 3], 1e-10);
        assert!(solver.is_done());
        let root = solver.root();
        let near = |expected: &[f64]| root.iter().zip(expected).all(|(a, e)| (a - e).abs() <= 1e-9);
        assert!(near(&[0.0, 2.0, 0.5]) || near(&[1.0, 1.0, 0.5]), "{root:?}");
    }

    #[test]
    fn already_a_root() {
        let mut f = system(3, 3, |p: &[f64], f: &mut [f64], d: &mut [f64]| {
            for i in 0..3 {
                f[i] = p[i] - 1.0;
            }
            d.copy_from_slice(&[1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]);
        });
        let solver = solve(&mut f, &[1.0; 3], &[-5.0; 3], &[5.0; 3], 1e-10);
        assert!(solver.is_done());
        assert_eq!(solver.nb_iterations(), 0);
        assert_eq!(f.evaluations, 1);
    }

    #[test]
    fn wide_and_tall_systems() {
        // Two equations in three unknowns (the plane x0 + x1 + x2 = 3 meeting x0 = x1):
        // the wide SVD branch, any point of the line is a root.
        let mut wide = system(3, 2, |p: &[f64], f: &mut [f64], d: &mut [f64]| {
            f[0] = p[0] + p[1] + p[2] - 3.0;
            f[1] = p[0] - p[1];
            d.copy_from_slice(&[1.0, 1.0, 1.0, 1.0, -1.0, 0.0]);
        });
        let solver = solve(&mut wide, &[0.0, 2.0, 0.0], &[-5.0; 3], &[5.0; 3], 1e-10);
        assert!(solver.is_done());
        let p = solver.root();
        assert!((p[0] + p[1] + p[2] - 3.0).abs() < 1e-9 && (p[0] - p[1]).abs() < 1e-9, "{p:?}");
        // Three consistent equations in two unknowns: math_GaussLeastSquare.
        let mut tall = system(2, 3, |p: &[f64], f: &mut [f64], d: &mut [f64]| {
            f[0] = p[0] - 1.0;
            f[1] = p[1] - 2.0;
            f[2] = p[0] + p[1] - 3.0;
            d.copy_from_slice(&[1.0, 0.0, 0.0, 1.0, 1.0, 1.0]);
        });
        let mut solver = FunctionSetRoot::new(2, 100);
        solver.set_tolerance(&[1e-10; 2]);
        solver.perform(&mut tall, &[0.0, 0.0], &[-5.0; 2], &[5.0; 2]);
        assert!(solver.is_done());
        assert_close(solver.root(), &[1.0, 2.0], 1e-9);
    }
}
