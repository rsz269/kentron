//! B-spline curves through or near given points.

use crate::error::{Error, Result};

use super::bspline::BSplineCurve;
use super::control::ControlPoint;
use super::parametric::ParametricCurve;
use super::Continuity;

/// Solves the square system `rows · x = rhs` (three right-hand sides) by Gaussian
/// elimination with partial pivoting; `None` if it is singular.
fn solve_dense(mut rows: Vec<Vec<f64>>, mut rhs: Vec<[f64; 3]>) -> Option<Vec<[f64; 3]>> {
    let n = rows.len();
    for col in 0..n {
        let pivot = (col..n).max_by(|&a, &b| rows[a][col].abs().total_cmp(&rows[b][col].abs()))?;
        if rows[pivot][col].abs() < 1e-300 {
            return None;
        }
        rows.swap(col, pivot);
        rhs.swap(col, pivot);
        for r in col + 1..n {
            let f = rows[r][col] / rows[col][col];
            if f == 0.0 {
                continue;
            }
            for k in col..n {
                rows[r][k] -= f * rows[col][k];
            }
            for d in 0..3 {
                rhs[r][d] -= f * rhs[col][d];
            }
        }
    }
    let mut x = vec![[0.0; 3]; n];
    for r in (0..n).rev() {
        for d in 0..3 {
            let s: f64 = (r + 1..n).map(|k| rows[r][k] * x[k][d]).sum();
            x[r][d] = (rhs[r][d] - s) / rows[r][r];
        }
    }
    Some(x)
}

/// The least-squares B-spline of `degree` with `m` poles through the first and last of
/// `q` (at `params`, from 0 to 1) and near the others (The NURBS Book, A9.7).
fn least_squares<P: ControlPoint>(q: &[[f64; 3]], params: &[f64], degree: usize, m: usize) -> Option<BSplineCurve<P>> {
    least_squares_with_knots(q, params, degree, knots_for(params, degree, m))
}

/// The least-squares B-spline of `degree` on flat knots `flat` (clamped, spanning
/// `params`) through the first and last of `q` and near the others.
pub(crate) fn least_squares_with_knots<P: ControlPoint>(q: &[[f64; 3]], params: &[f64], degree: usize, flat: Vec<f64>) -> Option<BSplineCurve<P>> {
    faired_least_squares(q, params, degree, flat, &[])
}

/// The least-squares B-spline of [`least_squares_with_knots`], faired: plus, for each
/// `(order, weight)` of `fairing`, `weight` times the integral of the square of the curve's
/// derivative of that order (as a share of the fit's own size: 1 weighs them alike).
fn faired_least_squares<P: ControlPoint>(q: &[[f64; 3]], params: &[f64], degree: usize, flat: Vec<f64>, fairing: &[(usize, f64)]) -> Option<BSplineCurve<P>> {
    let (n, p) = (q.len(), degree);
    let m = flat.len() - p - 1;
    if flat.windows(2).any(|w| w[1] < w[0]) || m < 2 {
        return None;
    }
    let basis = |t: f64| -> Vec<f64> {
        let span = (p..m).rev().find(|&k| flat[k] <= t && flat[k] < flat[k + 1]).unwrap_or(p).min(m - 1);
        let b = crate::geometry::bspline::basis_derivatives(&flat, span, t, p, 0);
        let mut row = vec![0.0; m];
        for (j, &v) in b[0].iter().enumerate() {
            row[span - p + j] = v;
        }
        row
    };
    // Interior poles from the normal equations, the end ones fixed at the end points.
    let inner = m - 2;
    let (first, last) = (q[0], q[n - 1]);
    let mut poles = vec![first];
    if inner > 0 {
        let mut a = vec![vec![0.0; inner]; inner];
        let mut b = vec![[0.0; 3]; inner];
        for (k, &t) in params.iter().enumerate().take(n - 1).skip(1) {
            let row = basis(t);
            let r: [f64; 3] = std::array::from_fn(|c| q[k][c] - row[0] * first[c] - row[m - 1] * last[c]);
            for i in 0..inner {
                let ni = row[i + 1];
                if ni == 0.0 {
                    continue;
                }
                for j in 0..inner {
                    a[i][j] += ni * row[j + 1];
                }
                for c in 0..3 {
                    b[i][c] += ni * r[c];
                }
            }
        }
        let fit_size: f64 = (0..inner).map(|i| a[i][i]).sum();
        for &(order, weight) in fairing {
            if weight <= 0.0 || order > p {
                continue;
            }
            let gram = derivative_gram(&flat, p, order);
            let size: f64 = (1..m - 1).map(|i| gram[i][i]).sum();
            if !(size > 0.0) {
                continue;
            }
            let w = weight * fit_size / size;
            for i in 0..inner {
                for j in 0..inner {
                    a[i][j] += w * gram[i + 1][j + 1];
                }
                for c in 0..3 {
                    b[i][c] -= w * (gram[i + 1][0] * first[c] + gram[i + 1][m - 1] * last[c]);
                }
            }
        }
        poles.extend(solve_dense(a, b)?);
    }
    poles.push(last);
    let (knots, multiplicities) = super::knots::group(&flat);
    BSplineCurve::try_new(poles.into_iter().map(P::from_coords).collect(), knots, multiplicities, p).ok()
}

/// The integrals over the knots' range of the products of the B-spline basis functions'
/// derivatives of `order` (degree `p`, flat knots `flat`).
fn derivative_gram(flat: &[f64], p: usize, order: usize) -> Vec<Vec<f64>> {
    let m = flat.len() - p - 1;
    let mut gram = vec![vec![0.0; m]; m];
    for span in p..m {
        let (a, b) = (flat[span], flat[span + 1]);
        if b <= a {
            continue;
        }
        let (mid, half) = (0.5 * (a + b), 0.5 * (b - a));
        for &(x, w) in &super::measure::GAUSS {
            let d = crate::geometry::bspline::basis_derivatives(flat, span, mid + half * x, p, order);
            for i in 0..=p {
                for j in 0..=p {
                    gram[span - p + i][span - p + j] += w * half * d[order][i] * d[order][j];
                }
            }
        }
    }
    gram
}

/// The B-spline curve passing exactly through `points`, in order (periodic: closing back
/// to the first point smoothly). Points closer than `tolerance` to each other make the
/// problem fail.
#[doc(alias = "GeomAPI_Interpolate")]
#[doc(alias = "Geom2dAPI_Interpolate")]
pub fn interpolate<P: ControlPoint>(points: &[P], periodic: bool, tolerance: f64) -> Result<BSplineCurve<P>> {
    // Chord-length parameters, as OCCT's GeomAPI_Interpolate.
    let mut params = Vec::with_capacity(points.len());
    let mut t = 0.0;
    for (i, p) in points.iter().enumerate() {
        if i > 0 {
            let step = distance(points[i - 1], *p);
            if step <= tolerance {
                return Err(Error::InvalidInput(format!("interpolation points {} and {i} coincide", i - 1)));
            }
            t += step;
        }
        params.push(t);
    }
    if periodic {
        let closing = distance(points[points.len() - 1], points[0]);
        if closing <= tolerance {
            return Err(Error::InvalidInput("the last interpolation point repeats the first".into()));
        }
        return interpolate_periodic_with_parameters(points, &params, t + closing);
    }
    interpolate_with_parameters(points, &params)
}

/// The periodic cubic B-spline curve passing through `points[i]` at parameter
/// `params[i]` and closing back to the first point smoothly (C2) after `period`; the
/// parameters must increase within less than a period, and the first point is not
/// repeated at the end. The knots are the parameters.
#[doc(alias = "GeomAPI_Interpolate")]
pub fn interpolate_periodic_with_parameters<P: ControlPoint>(points: &[P], params: &[f64], period: f64) -> Result<BSplineCurve<P>> {
    let n = points.len();
    if n < 4 || params.len() != n || params.windows(2).any(|w| !(w[0] < w[1])) || !(params[n - 1] < params[0] + period) {
        return Err(Error::InvalidInput("periodic interpolation needs four or more points with increasing parameters within a period".into()));
    }
    let p = 3;
    let mut knots = params.to_vec();
    knots.push(params[0] + period);
    let multiplicities = vec![1; n + 1];
    let flat = crate::geometry::bspline::periodic_flatten(&knots, &multiplicities, p);
    // At knot k the non-zero basis functions are those of poles k, k + 1 and k + 2
    // (cyclically): with unknown x[c] for pole c + 1 the system is cyclic tridiagonal.
    let (mut below, mut diagonal, mut above) = (vec![0.0; n], vec![0.0; n], vec![0.0; n]);
    for k in 0..n {
        let b = crate::geometry::bspline::basis_derivatives(&flat, p + k, params[k], p, 0);
        (below[k], diagonal[k], above[k]) = (b[0][0], b[0][1], b[0][2]);
    }
    let rhs: Vec<[f64; 3]> = points.iter().map(|q| q.coords()).collect();
    let x = solve_cyclic_tridiagonal(&below, &diagonal, &above, rhs).ok_or_else(|| Error::InvalidInput("singular interpolation problem".into()))?;
    let poles = (0..n).map(|j| P::from_coords(x[(j + n - 1) % n])).collect();
    BSplineCurve::try_new_periodic(poles, knots, multiplicities, p)
}

/// Solves the cyclic tridiagonal system whose row `i` is `below[i] x[i-1] + diagonal[i]
/// x[i] + above[i] x[i+1] = rhs[i]` (indices modulo the size) by the Sherman–Morrison
/// correction of the plain tridiagonal one; `None` if it is singular.
fn solve_cyclic_tridiagonal(below: &[f64], diagonal: &[f64], above: &[f64], rhs: Vec<[f64; 3]>) -> Option<Vec<[f64; 3]>> {
    let n = diagonal.len();
    let (alpha, beta) = (above[n - 1], below[0]);
    let gamma = -diagonal[0];
    let mut modified = diagonal.to_vec();
    modified[0] -= gamma;
    modified[n - 1] -= alpha * beta / gamma;
    // Thomas algorithm on several right-hand sides.
    let solve = |rhs: Vec<[f64; 3]>| -> Option<Vec<[f64; 3]>> {
        let mut c = vec![0.0; n];
        let mut d = rhs;
        let mut m = modified[0];
        if m.abs() < 1e-300 {
            return None;
        }
        c[0] = above[0] / m;
        d[0] = d[0].map(|v| v / m);
        for i in 1..n {
            m = modified[i] - below[i] * c[i - 1];
            if m.abs() < 1e-300 {
                return None;
            }
            c[i] = above[i] / m;
            let prev = d[i - 1];
            d[i] = std::array::from_fn(|k| (d[i][k] - below[i] * prev[k]) / m);
        }
        for i in (0..n - 1).rev() {
            let next = d[i + 1];
            d[i] = std::array::from_fn(|k| d[i][k] - c[i] * next[k]);
        }
        Some(d)
    };
    let x = solve(rhs)?;
    let mut u = vec![[0.0; 3]; n];
    u[0] = [gamma; 3];
    u[n - 1] = [alpha; 3];
    let z = solve(u)?;
    let denominator = 1.0 + z[0][0] + beta * z[n - 1][0] / gamma;
    if denominator.abs() < 1e-300 {
        return None;
    }
    Some(
        (0..n)
            .map(|i| {
                std::array::from_fn(|k| {
                    let factor = (x[0][k] + beta * x[n - 1][k] / gamma) / denominator;
                    x[i][k] - factor * z[i][0]
                })
            })
            .collect(),
    )
}

/// The cubic (lower degree for fewer than four points) B-spline curve passing through
/// `points[i]` at parameter `params[i]`; parameters must increase. Knots are placed by
/// averaging the parameters, so the curve is C2.
#[doc(alias = "GeomAPI_Interpolate")]
pub fn interpolate_with_parameters<P: ControlPoint>(points: &[P], params: &[f64]) -> Result<BSplineCurve<P>> {
    let n = points.len();
    if n < 2 || params.len() != n || params.windows(2).any(|w| !(w[0] < w[1])) {
        return Err(Error::InvalidInput("interpolation needs two or more points with increasing parameters".into()));
    }
    let p = (n - 1).min(3);
    // Flat knots by averaging (The NURBS Book, eq. 9.8).
    let mut flat = vec![params[0]; p + 1];
    for j in 1..n - p {
        flat.push(params[j..j + p].iter().sum::<f64>() / p as f64);
    }
    flat.extend(std::iter::repeat(params[n - 1]).take(p + 1));
    // The collocation matrix: row i holds the basis functions at params[i], within p of
    // the diagonal (with knots averaged, params[i] lies in the support of the i-th one).
    let mut matrix = Band::new(n, p);
    for (i, &t) in params.iter().enumerate() {
        let mut span = p;
        while span + 1 < n && flat[span + 1] <= t {
            span += 1;
        }
        // (Degree three at most: the basis functions on the stack.)
        let (basis, _) = crate::geometry::bspline::basis_and_first::<4>(&flat, 0, span, t, p, false);
        for (j, &b) in basis[..=p].iter().enumerate() {
            if b != 0.0 {
                *matrix.at(i, span - p + j).ok_or_else(|| Error::InvalidInput("interpolation parameters out of order with the knots".into()))? = b;
            }
        }
    }
    let rhs: Vec<[f64; 3]> = points.iter().map(|p| p.coords()).collect();
    let poles = matrix.solve(rhs).ok_or_else(|| Error::InvalidInput("singular interpolation problem".into()))?;
    let mut knots: Vec<f64> = Vec::new();
    let mut multiplicities: Vec<usize> = Vec::new();
    for k in flat {
        match knots.last() {
            Some(&last) if last == k => *multiplicities.last_mut().expect("multiplicity") += 1,
            _ => {
                knots.push(k);
                multiplicities.push(1);
            }
        }
    }
    BSplineCurve::try_new(poles.into_iter().map(P::from_coords).collect(), knots, multiplicities, p)
}

fn distance<P: ControlPoint>(a: P, b: P) -> f64 {
    let (a, b) = (a.coords(), b.coords());
    (0..3).map(|i| (a[i] - b[i]).powi(2)).sum::<f64>().sqrt()
}

/// A square matrix whose non-zeros lie within `band` of the diagonal, stored by rows of
/// `2 band + 1`.
struct Band {
    n: usize,
    band: usize,
    values: Vec<f64>,
}

impl Band {
    fn new(n: usize, band: usize) -> Band {
        Band { n, band, values: vec![0.0; n * (2 * band + 1)] }
    }

    /// The entry at (`row`, `col`); `None` outside the band.
    fn at(&mut self, row: usize, col: usize) -> Option<&mut f64> {
        let offset = (col + self.band).checked_sub(row).filter(|&o| o <= 2 * self.band && col < self.n)?;
        self.values.get_mut(row * (2 * self.band + 1) + offset)
    }

    fn get(&self, row: usize, col: usize) -> f64 {
        match (col + self.band).checked_sub(row).filter(|&o| o <= 2 * self.band && col < self.n) {
            Some(offset) => self.values[row * (2 * self.band + 1) + offset],
            None => 0.0,
        }
    }

    /// Solves `self · x = rhs` (three right-hand sides) by Gaussian elimination without
    /// pivoting, which stays within the band and is stable for B-spline collocation
    /// matrices (totally positive: de Boor, A Practical Guide to Splines, ch. XIII);
    /// `None` if a pivot vanishes.
    fn solve(mut self, mut rhs: Vec<[f64; 3]>) -> Option<Vec<[f64; 3]>> {
        let (n, band) = (self.n, self.band);
        for col in 0..n {
            let pivot = self.get(col, col);
            if pivot.abs() < 1e-300 {
                return None;
            }
            for row in col + 1..(col + band + 1).min(n) {
                let f = self.get(row, col) / pivot;
                if f == 0.0 {
                    continue;
                }
                for k in col..(col + band + 1).min(n) {
                    let above = self.get(col, k);
                    *self.at(row, k).expect("within the band") -= f * above;
                }
                for d in 0..3 {
                    rhs[row][d] -= f * rhs[col][d];
                }
            }
        }
        let mut x = vec![[0.0; 3]; n];
        for row in (0..n).rev() {
            for d in 0..3 {
                let s: f64 = (row + 1..(row + band + 1).min(n)).map(|k| self.get(row, k) * x[k][d]).sum();
                x[row][d] = (rhs[row][d] - s) / self.get(row, row);
            }
        }
        Some(x)
    }
}

/// Like [`interpolate`] (non-periodic) with prescribed tangents at the first and last
/// points. If `scale`, the tangents only give the directions and are rescaled to the
/// parametrization.
#[doc(alias = "GeomAPI_Interpolate::Load")]
pub fn interpolate_with_end_tangents<P: ControlPoint>(
    points: &[P],
    first_tangent: P::Vector,
    last_tangent: P::Vector,
    scale: bool,
    tolerance: f64,
) -> Result<BSplineCurve<P>> {
    let n = points.len();
    let mut params = Vec::with_capacity(n);
    let mut t = 0.0;
    for (i, p) in points.iter().enumerate() {
        if i > 0 {
            let step = distance(points[i - 1], *p);
            if step <= tolerance {
                return Err(Error::InvalidInput(format!("interpolation points {} and {i} coincide", i - 1)));
            }
            t += step;
        }
        params.push(t);
    }
    if n < 2 {
        return Err(Error::InvalidInput("interpolation needs two or more points".into()));
    }
    // With chord-length parameters the curve runs at about unit speed: scaled tangents are
    // unit ones.
    let scaled = |v: P::Vector| {
        let c = P::vector_coords(v);
        let length = (c[0] * c[0] + c[1] * c[1] + c[2] * c[2]).sqrt();
        if scale && length > 0.0 { [c[0] / length, c[1] / length, c[2] / length] } else { c }
    };
    // A cubic through the points with the end derivatives: knots at the parameters, ends
    // of multiplicity 4; n + 2 poles from n points and 2 derivatives.
    let p = 3;
    let mut flat = vec![params[0]; p + 1];
    flat.extend_from_slice(&params[1..n - 1]);
    flat.extend(std::iter::repeat(params[n - 1]).take(p + 1));
    let m = n + 2;
    let mut rows: Vec<Vec<f64>> = Vec::with_capacity(m);
    let mut rhs: Vec<[f64; 3]> = Vec::with_capacity(m);
    let basis_row = |t: f64, order: usize| -> Vec<f64> {
        let span = (p..m).rev().find(|&k| flat[k] <= t).unwrap_or(p).min(m - 1);
        let b = crate::geometry::bspline::basis_derivatives(&flat, span, t, p, order);
        let mut row = vec![0.0; m];
        for (j, &v) in b[order].iter().enumerate() {
            row[span - p + j] = v;
        }
        row
    };
    for (i, (&t, q)) in params.iter().zip(points).enumerate() {
        rows.push(basis_row(t, 0));
        rhs.push(q.coords());
        if i == 0 {
            rows.push(basis_row(t, 1));
            rhs.push(scaled(first_tangent));
        }
    }
    rows.insert(rows.len() - 1, basis_row(params[n - 1], 1));
    rhs.insert(rhs.len() - 1, scaled(last_tangent));
    let poles = solve_dense(rows, rhs).ok_or_else(|| Error::InvalidInput("singular interpolation problem".into()))?;
    let (knots, multiplicities) = super::knots::group(&flat);
    BSplineCurve::try_new(poles.into_iter().map(P::from_coords).collect(), knots, multiplicities, p)
}

/// The B-spline curve through `points` laid out as OCCT's GeomAPI_Interpolate lays it out:
/// knots at the parameters (`params`, else chord lengths), cubic, with the tangents given
/// (`tangents`: none, or one per point, `None` where free) as conditions at their points
/// (doubling the knot there). The end tangents not given are those of the Lagrange
/// polynomial through the first (last) four points; with fewer points and no tangents the
/// curve is the line, or the parabola through three. If `scale`, the given tangents only
/// give directions: each is scaled to the derivative there of the parabola through
/// neighbouring points. Periodic curves close back to the first point (C1 there, its
/// tangent being the Lagrange one or given); their `params` have one more value, where the
/// curve returns to the first point. Points (and tangents) shorter than `tolerance` apart
/// (long) make the problem fail.
#[doc(alias = "GeomAPI_Interpolate::Load")]
pub fn interpolate_with_tangents<P: ControlPoint>(
    points: &[P],
    params: Option<&[f64]>,
    tangents: &[Option<P::Vector>],
    periodic: bool,
    scale: bool,
    tolerance: f64,
) -> Result<BSplineCurve<P>> {
    let invalid = |msg: &str| Err(Error::InvalidInput(msg.into()));
    let n = points.len();
    if n < 2 {
        return invalid("interpolation needs two or more points");
    }
    let q: Vec<[f64; 3]> = points.iter().map(|p| p.coords()).collect();
    if (1..n).any(|i| dist3(q[i - 1], q[i]) < tolerance) {
        return invalid("interpolation points closer than the tolerance");
    }
    let t: Vec<f64> = match params {
        Some(params) => params.to_vec(),
        None => {
            let mut t = vec![0.0];
            for i in 1..n {
                t.push(t[i - 1] + dist3(q[i - 1], q[i]));
            }
            if periodic {
                t.push(t[n - 1] + dist3(q[n - 1], q[0]));
            }
            t
        }
    };
    if t.len() != n + usize::from(periodic) {
        return invalid("interpolation parameters must number the points (and one more, periodic)");
    }
    if t.windows(2).any(|w| !(w[1] - w[0] >= f64::MIN_POSITIVE)) {
        return invalid("interpolation parameters must increase");
    }
    let requested = !tangents.is_empty();
    if requested && tangents.len() != n {
        return invalid("interpolation tangents must number the points");
    }
    let mut d: Vec<Option<[f64; 3]>> = tangents.iter().map(|v| v.map(P::vector_coords)).collect();
    if d.iter().flatten().any(|v| norm3(*v) < tolerance) {
        return invalid("an interpolation tangent shorter than the tolerance");
    }
    if requested && scale {
        // (OCCT's ScaleTangents: L1 norms, parabolas through a sliding window of points.)
        let degree = if n == 2 { 1 } else { 2 };
        let mut index = 0;
        for i in 0..n {
            if let Some(v) = d[i] {
                let derivative = lagrange_derivative(&q[index..=index + degree], &t[index..=index + degree], t[i]);
                let ratio = l1(derivative) / l1(v);
                d[i] = Some(v.map(|c| c * ratio));
                if i != 0 {
                    index += 1;
                }
                index = index.min(n - 1 - degree);
            }
        }
    }
    if !requested {
        d = vec![None; n];
        if periodic && n == 2 {
            let knots = t.clone();
            return BSplineCurve::try_new_periodic(points.to_vec(), knots, vec![1; n + 1], 1);
        }
        if !periodic && n == 2 {
            return BSplineCurve::try_new(points.to_vec(), vec![t[0], t[1]], vec![2, 2], 1);
        }
        if !periodic && n == 3 {
            let flat = [t[0], t[0], t[0], t[2], t[2], t[2]];
            let rows = t.iter().map(|&u| basis_row(&flat, 2, u, 0)).collect();
            let poles = solve_dense(rows, q.clone()).ok_or_else(|| Error::InvalidInput("singular interpolation problem".into()))?;
            return BSplineCurve::try_new(poles.into_iter().map(P::from_coords).collect(), vec![t[0], t[2]], vec![3, 3], 2);
        }
    }
    // The free end tangents (the first only, periodic): Lagrange cubics (parabolas through
    // three points) through the end points.
    let degree = if n == 3 { 2 } else { 3 };
    if d[0].is_none() || (!periodic && d[n - 1].is_none()) {
        if n < 3 {
            return invalid("interpolation through two points needs both end tangents");
        }
    }
    if d[0].is_none() {
        d[0] = Some(lagrange_derivative(&q[..=degree], &t[..=degree], t[0]));
    }
    if !periodic && d[n - 1].is_none() {
        d[n - 1] = Some(lagrange_derivative(&q[n - 1 - degree..], &t[n - 1 - degree..n], t[n - 1]));
    }
    // The conditions in order (parameter, derivative order, value), with the flat knots.
    let p = 3;
    let mut conditions: Vec<(f64, usize, [f64; 3])> = vec![(t[0], 0, q[0]), (t[0], 1, d[0].expect("first tangent"))];
    let mut inner = Vec::new();
    let last_interior = if periodic { n } else { n - 1 };
    for i in 1..last_interior {
        conditions.push((t[i], 0, q[i]));
        inner.push(t[i]);
        if let Some(v) = d[i] {
            conditions.push((t[i], 1, v));
            inner.push(t[i]);
        }
    }
    let (end, closing) = if periodic { (t[n], q[0]) } else { (t[n - 1], q[n - 1]) };
    let end_tangent = if periodic { d[0] } else { d[n - 1] }.expect("last tangent");
    conditions.push((end, 1, end_tangent));
    conditions.push((end, 0, closing));
    let m = conditions.len();
    let mut flat = Vec::with_capacity(m + p + 1);
    if periodic {
        let period = t[n] - t[0];
        flat.extend([t[n - 1] - period, t[n - 1] - period, t[0], t[0]]);
        flat.extend(&inner);
        flat.extend([t[n], t[n], t[1] + period, t[1] + period]);
    } else {
        flat.extend([t[0]; 4]);
        flat.extend(&inner);
        flat.extend([t[n - 1]; 4]);
    }
    let rows = conditions.iter().map(|&(u, order, _)| basis_row(&flat, p, u, order)).collect();
    let rhs = conditions.iter().map(|&(_, _, v)| v).collect();
    let poles = solve_dense(rows, rhs).ok_or_else(|| Error::InvalidInput("singular interpolation problem".into()))?;
    let poles: Vec<P> = poles.into_iter().map(P::from_coords).collect();
    if periodic {
        // Knots at the parameters, doubled at the ends and where tangents are given.
        let mut multiplicities = vec![2];
        multiplicities.extend((1..n).map(|i| if d[i].is_some() { 2 } else { 1 }));
        multiplicities.push(2);
        return BSplineCurve::try_new_periodic(poles[..m - 2].to_vec(), t, multiplicities, p);
    }
    let (knots, multiplicities) = super::knots::group(&flat);
    BSplineCurve::try_new(poles, knots, multiplicities, p)
}

/// The row of `degree` B-spline basis functions' derivatives of `order` at `t` over flat
/// knots `flat` (the last non-empty span at the end).
fn basis_row(flat: &[f64], degree: usize, t: f64, order: usize) -> Vec<f64> {
    let m = flat.len() - degree - 1;
    let span = (degree..m).rev().find(|&k| flat[k] <= t && flat[k] < flat[k + 1]).unwrap_or(degree);
    let b = crate::geometry::bspline::basis_derivatives(flat, span, t, degree, order);
    let mut row = vec![0.0; m];
    for (j, &v) in b[order].iter().enumerate() {
        row[span - degree + j] = v;
    }
    row
}

/// The derivative at `x` of the polynomial through `q[i]` at `t[i]` (Lagrange's form).
fn lagrange_derivative(q: &[[f64; 3]], t: &[f64], x: f64) -> [f64; 3] {
    let k = q.len();
    let mut out = [0.0; 3];
    for j in 0..k {
        // l_j'(x) = sum over i != j of 1 / (t_j - t_i) prod over m != i, j of (x - t_m) / (t_j - t_m).
        let mut w = 0.0;
        for i in (0..k).filter(|&i| i != j) {
            let mut term = 1.0 / (t[j] - t[i]);
            for m in (0..k).filter(|&m| m != i && m != j) {
                term *= (x - t[m]) / (t[j] - t[m]);
            }
            w += term;
        }
        for c in 0..3 {
            out[c] += w * q[j][c];
        }
    }
    out
}

fn dist3(a: [f64; 3], b: [f64; 3]) -> f64 {
    norm3([a[0] - b[0], a[1] - b[1], a[2] - b[2]])
}

fn norm3(v: [f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

fn l1(v: [f64; 3]) -> f64 {
    v[0].abs() + v[1].abs() + v[2].abs()
}

/// Options for [`approximate_with`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ApproximationOptions {
    pub min_degree: usize,
    pub max_degree: usize,
    pub continuity: Continuity,
    /// Largest allowed distance from the points.
    pub tolerance: f64,
    /// Weights of the curve's length, bending and twisting to keep small besides (its first,
    /// second and third derivatives' squares integrated), for a fair curve: then of the
    /// greatest degree, and as fair as it may be while near enough the points.
    pub smoothing: Option<[f64; 3]>,
}

impl Default for ApproximationOptions {
    fn default() -> Self {
        Self { min_degree: 3, max_degree: 8, continuity: Continuity::C2, tolerance: 1.0e-3, smoothing: None }
    }
}

/// A smooth B-spline curve passing near `points` (within 1e-3, degree 3 to 8, C2).
#[doc(alias = "GeomAPI_PointsToBSpline")]
#[doc(alias = "Geom2dAPI_PointsToBSpline")]
pub fn approximate<P: ControlPoint>(points: &[P]) -> Result<BSplineCurve<P>> {
    approximate_with(points, ApproximationOptions::default())
}

pub fn approximate_with<P: ControlPoint>(points: &[P], options: ApproximationOptions) -> Result<BSplineCurve<P>> {
    let n = points.len();
    if n < 2 {
        return Err(Error::InvalidInput("approximation needs two or more points".into()));
    }
    // Chord-length parameters over [0, 1].
    let mut params = vec![0.0];
    for i in 1..n {
        params.push(params[i - 1] + distance(points[i - 1], points[i]));
    }
    let total = params[n - 1];
    if total <= 0.0 {
        return Err(Error::InvalidInput("approximation points coincide".into()));
    }
    params.iter_mut().for_each(|t| *t /= total);
    let q: Vec<[f64; 3]> = points.iter().map(|p| p.coords()).collect();
    if let Some(weights) = options.smoothing {
        return smoothed(points, &q, &params, weights, options);
    }
    let degree = options.min_degree.max(1).min(n - 1).min(options.max_degree.max(1));
    // More poles until the curve passes within the tolerance (with n poles it interpolates).
    for m in (degree + 1)..=n {
        let Some(curve) = least_squares::<P>(&q, &params, degree, m) else { continue };
        let worst = params.iter().zip(&q).map(|(&t, target)| {
            let c = curve.point(t).coords();
            (0..3).map(|k| (c[k] - target[k]).powi(2)).sum::<f64>().sqrt()
        }).fold(0.0, f64::max);
        if worst <= options.tolerance || m == n {
            return Ok(curve);
        }
    }
    interpolate_with_parameters(points, &params)
}

/// [`approximate_with`] with smoothing: of the greatest degree, with the fewest poles (as
/// near as a search finds) that, faired as much as may be, pass near enough the points.
fn smoothed<P: ControlPoint>(points: &[P], q: &[[f64; 3]], params: &[f64], weights: [f64; 3], options: ApproximationOptions) -> Result<BSplineCurve<P>> {
    let n = q.len();
    let least = match options.continuity {
        Continuity::C0 => 1,
        Continuity::C1 | Continuity::G1 => 2,
        _ => 3,
    };
    let degree = options.max_degree.max(least).min(n - 1).max(1);
    let near = |curve: &BSplineCurve<P>| {
        params.iter().zip(q).all(|(&t, target)| {
            let c = curve.point(t).coords();
            (0..3).map(|k| (c[k] - target[k]).powi(2)).sum::<f64>().sqrt() <= options.tolerance
        })
    };
    // (Fairing ever less, down to none: the fairest near enough.)
    let with = |m: usize| -> Option<BSplineCurve<P>> {
        for share in [1e-2, 1e-4, 1e-6, 1e-8, 0.0] {
            let fairing: Vec<(usize, f64)> = weights.iter().enumerate().map(|(k, &w)| (k + 1, w * share)).collect();
            let flat = knots_for(params, degree, m);
            if let Some(curve) = faired_least_squares::<P>(q, params, degree, flat, &fairing) {
                if near(&curve) {
                    return Some(curve);
                }
            }
        }
        None
    };
    // Poles growing by half until near enough, then the fewest between.
    let (mut lo, mut hi) = (degree, None);
    let mut m = degree + 1;
    while m < n {
        if let Some(curve) = with(m) {
            hi = Some((m, curve));
            break;
        }
        lo = m;
        m = (m + m / 2).max(m + 1).min(n);
    }
    let Some((mut m_hi, mut best)) = hi else {
        return with(n).ok_or(()).or_else(|_| interpolate_with_parameters(points, params));
    };
    while m_hi - lo > 1 {
        let mid = (lo + m_hi) / 2;
        match with(mid) {
            Some(curve) => (m_hi, best) = (mid, curve),
            None => lo = mid,
        }
    }
    Ok(best)
}

/// A B-spline surface passing near a grid of points (`points[i][j]`, `i` along u and `j`
/// along v), within `options.tolerance`: of degree 3 in each direction where the options
/// allow it (else the nearest they allow, and less where there are too few points), with
/// as few poles as a search finds. Each column of points fitted along u on common knots,
/// then each row of those curves' poles along v.
#[doc(alias = "GeomAPI_PointsToBSplineSurface")]
pub fn approximate_surface(points: &[Vec<crate::math::Point>], options: ApproximationOptions) -> Result<super::bspline::BSplineSurface> {
    use crate::math::Point;
    let (nu, nv) = (points.len(), points.first().map_or(0, Vec::len));
    if nu < 2 || nv < 2 || points.iter().any(|r| r.len() != nv) {
        return Err(Error::InvalidInput("surface approximation needs a grid of two or more points each way".into()));
    }
    // Chord-length parameters, averaged across the grid, over [0, 1].
    let average = |count: usize, other: usize, at: &dyn Fn(usize, usize) -> Point| -> Result<Vec<f64>> {
        let mut params = vec![0.0; count];
        for k in 0..other {
            let mut t = vec![0.0];
            for i in 1..count {
                t.push(t[i - 1] + at(i - 1, k).distance(at(i, k)));
            }
            let total = t[count - 1];
            if total > 0.0 {
                for i in 0..count {
                    params[i] += t[i] / total / other as f64;
                }
            } else {
                for i in 0..count {
                    params[i] += i as f64 / (count - 1) as f64 / other as f64;
                }
            }
        }
        params[0] = 0.0;
        params[count - 1] = 1.0;
        if params.windows(2).any(|w| !(w[1] > w[0])) {
            return Err(Error::InvalidInput("surface approximation points coincide".into()));
        }
        Ok(params)
    };
    let us = average(nu, nv, &|i, j| points[i][j])?;
    let vs = average(nv, nu, &|j, i| points[i][j])?;
    // (A way with no more points than a span of the greatest degree has poles: that one
    // span through all of them, exactly: as OCCT's, so that faces fitted to grids sharing a
    // row of points share the curve through it.)
    let whole = |n: usize| n <= options.max_degree.max(1) + 1;
    let degree = |n: usize| if whole(n) { n - 1 } else { 3usize.clamp(options.min_degree.max(1), options.max_degree.max(1)).min(n - 1) };
    let (p, q) = (degree(nu), degree(nv));
    // The poles of every column fitted along u with `mu` poles, then every row of those along
    // v with `mv`: the surface, or `None` where a fit fails.
    let fit = |mu: usize, mv: usize| -> Option<super::bspline::BSplineSurface> {
        let (ku, kv) = (knots_for(&us, p, mu), knots_for(&vs, q, mv));
        let mut columns: Vec<Vec<Point>> = Vec::with_capacity(nv);
        for j in 0..nv {
            let q3: Vec<[f64; 3]> = (0..nu).map(|i| points[i][j].coords()).collect();
            let curve: BSplineCurve<Point> = least_squares_with_knots(&q3, &us, p, ku.clone())?;
            columns.push(curve.poles().to_vec());
        }
        let mut rows: Vec<Vec<Point>> = Vec::with_capacity(mu);
        for i in 0..mu {
            let q3: Vec<[f64; 3]> = (0..nv).map(|j| columns[j][i].coords()).collect();
            let curve: BSplineCurve<Point> = least_squares_with_knots(&q3, &vs, q, kv.clone())?;
            rows.push(curve.poles().to_vec());
        }
        let (ku_knots, ku_mults) = super::knots::group(&ku);
        let (kv_knots, kv_mults) = super::knots::group(&kv);
        super::bspline::BSplineSurface::try_new(
            rows,
            super::bspline::KnotVector::new(ku_knots, ku_mults, p),
            super::bspline::KnotVector::new(kv_knots, kv_mults, q),
        )
        .ok()
    };
    let near = |surface: &super::bspline::BSplineSurface| {
        use super::ParametricSurface;
        (0..nu).all(|i| (0..nv).all(|j| surface.point(us[i], vs[j]).distance(points[i][j]) <= options.tolerance))
    };
    // Poles each way growing by half until near enough (or as many as points), then the
    // fewest each way between.
    let grow = |m: usize, n: usize| (m + m / 2).max(m + 1).min(n);
    let (mut mu, mut mv) = (p + 1, q + 1);
    let mut found = None;
    loop {
        if let Some(surface) = fit(mu, mv).filter(|s| near(s)) {
            found = Some((mu, mv, surface));
            break;
        }
        if mu == nu && mv == nv {
            break;
        }
        (mu, mv) = (grow(mu, nu), grow(mv, nv));
    }
    let Some((mut mu, mut mv, mut best)) = found else {
        return fit(nu, nv).ok_or_else(|| Error::Failed("no surface through the points".into()));
    };
    for side in 0..2 {
        let mut lo = if side == 0 { p } else { q };
        loop {
            let hi = if side == 0 { mu } else { mv };
            if hi - lo <= 1 {
                break;
            }
            let mid = (lo + hi) / 2;
            let (a, b) = if side == 0 { (mid, mv) } else { (mu, mid) };
            match fit(a, b).filter(|s| near(s)) {
                Some(surface) => {
                    best = surface;
                    if side == 0 { mu = mid } else { mv = mid }
                }
                None => lo = mid,
            }
        }
    }
    Ok(best)
}

/// Clamped flat knots for `m` poles of `degree`, spread so that every span holds parameters
/// (The NURBS Book, eq. 9.69).
fn knots_for(params: &[f64], degree: usize, m: usize) -> Vec<f64> {
    let (n, p) = (params.len(), degree);
    let mut flat = vec![0.0; p + 1];
    let d = n as f64 / (m - p) as f64;
    for j in 1..m - p {
        let i = (j as f64 * d).floor() as usize;
        let alpha = j as f64 * d - i as f64;
        flat.push((1.0 - alpha) * params[i.saturating_sub(1).min(n - 1)] + alpha * params[i.min(n - 1)]);
    }
    flat.extend(std::iter::repeat(1.0).take(p + 1));
    flat
}

#[cfg(test)]
mod tests {
    use std::f64::consts::PI;

    use super::*;
    use crate::geometry::surface::Surface;
    use crate::geometry::ParametricCurve;
    use crate::math::Point;

    // CadQuery's spline tests (test_cadquery.py), whose expected points are OCCT's.
    #[test]
    fn interpolation_with_tangents_matches_occt() {
        use crate::math::Vector;
        let points = [Point::new(0.0, 0.0, 0.0), Point::new(1.0, 1.0, 0.0), Point::new(2.0, 0.0, 0.0), Point::new(1.0, -1.0, 0.0)];
        let all = [Vector::new(0.0, 1.0, 0.0), Vector::new(1.0, 0.0, 0.0), Vector::new(0.0, -1.0, 0.0), Vector::new(-1.0, 0.0, 0.0)].map(Some);
        let params = [0.0, 1.0, 2.0, 3.0];
        let at = |c: &BSplineCurve, t: f64, expected: (f64, f64)| {
            let p = c.point(t);
            assert!((p.x - expected.0).abs() < 1e-12 && (p.y - expected.1).abs() < 1e-12 && p.z.abs() < 1e-12, "{p:?} at {t}, not {expected:?}");
        };
        let scaled = interpolate_with_tangents(&points, Some(&params), &all, false, true, 1e-6).unwrap();
        at(&scaled, 0.5, (0.375, 0.875));
        at(&scaled, 2.5, (1.875, -0.625));
        let unscaled = interpolate_with_tangents(&points, Some(&params), &all, false, false, 1e-6).unwrap();
        at(&unscaled, 0.5, (0.375, 0.625));
        let mut some = all;
        some[1] = None;
        let partial = interpolate_with_tangents(&points, Some(&params), &some, false, true, 1e-6).unwrap();
        at(&partial, 1.5, (1.6875, 0.875));
        let tenfold = interpolate_with_tangents(&points, Some(&[0.0, 10.0, 20.0, 30.0]), &all, false, true, 1e-6).unwrap();
        at(&tenfold, 15.0, (1.625, 0.625));
        at(&scaled, 1.5, (1.625, 0.625));
        // Periodic, with parameters (one more: the return to the first point).
        let closed = interpolate_with_tangents(&points, Some(&[0.0, 1.0, 2.0, 3.0, 4.0]), &[], true, true, 1e-6).unwrap();
        at(&closed, 4.0, (0.0, 0.0));
        // A tangent shorter than the tolerance fails.
        let mut short = all;
        short[0] = Some(Vector::new(0.0, 0.5, 0.0));
        assert!(interpolate_with_tangents(&points, Some(&params), &short, false, true, 1.0).is_err());
    }

    #[test]
    fn periodic_interpolation_closes_smoothly() {
        // Uneven spacing round a whole circle.
        let angles: Vec<f64> = (0..24).map(|i| 2.0 * std::f64::consts::PI * (i as f64 + 0.3 * (i % 3) as f64) / 24.0).collect();
        let points: Vec<Point> = angles.iter().map(|a| Point::new(a.cos(), a.sin(), 0.0)).collect();
        let curve = interpolate_periodic_with_parameters(&points, &angles, 2.0 * std::f64::consts::PI).unwrap();
        assert!(curve.is_periodic());
        for (a, q) in angles.iter().zip(&points) {
            assert!(curve.point(*a).distance(*q) < 1e-12);
        }
        for i in 0..=300 {
            let p = curve.point(2.0 * std::f64::consts::PI * i as f64 / 300.0);
            assert!(((p - Point::ORIGIN).length() - 1.0).abs() < 1e-3, "{p:?}");
        }
        let (t0, t1) = (curve.first_parameter(), curve.last_parameter());
        let ((p0, v0, a0), (p1, v1, a1)) = (curve.d2(t0), curve.d2(t1));
        assert!(p0.distance(p1) < 1e-12 && (v0 - v1).length() < 1e-9 && (a0 - a1).length() < 1e-9, "a jump at the seam");
        let chord = interpolate(&points, true, 1e-9).unwrap();
        assert!(chord.point(chord.first_parameter()).distance(points[0]) < 1e-12);
    }

    #[test]
    fn periodic_curves_and_surfaces_reverse() {
        use crate::geometry::bspline::{BSplineSurface, KnotVector};
        use crate::geometry::{ParametricSurface, curve::Curve};
        let angles: Vec<f64> = (0..7).map(|i| 0.9 * i as f64).collect();
        let points: Vec<Point> = angles.iter().map(|a| Point::new(a.cos(), (2.0 * a).sin(), 0.1 * a)).collect();
        let fitted = interpolate_periodic_with_parameters(&points, &angles, 6.5).unwrap();
        // (A whole circle: rational, double knots.)
        let circle = crate::geometry::convert::curve_to_bspline(&Curve::Circle(crate::geometry::curve::Circle::new(crate::math::Frame::WORLD, 2.0))).unwrap();
        for curve in [fitted.clone(), circle] {
            let mut reversed = curve.clone();
            reversed.reverse();
            for i in 0..=20 {
                let t = curve.first_parameter() + (curve.last_parameter() - curve.first_parameter()) * i as f64 / 20.0;
                assert!(reversed.point(curve.reversed_parameter(t)).distance(curve.point(t)) < 1e-12);
            }
        }
        let rows: Vec<Vec<Point>> = fitted.poles().iter().map(|&q| vec![q, q + crate::math::Vector::new(0.0, 0.0, 1.0)]).collect();
        let u = KnotVector::periodic(fitted.knots().to_vec(), fitted.multiplicities().to_vec(), 3);
        let surface = BSplineSurface::try_new(rows, u, KnotVector::new(vec![0.0, 1.0], vec![2, 2], 1)).unwrap();
        // The same surface with u and v swapped.
        let columns: Vec<Vec<Point>> = (0..2).map(|k| (0..fitted.poles().len()).map(|i| surface.poles()[i][k]).collect()).collect();
        let swapped = BSplineSurface::try_new(columns, KnotVector::new(vec![0.0, 1.0], vec![2, 2], 1), KnotVector::periodic(fitted.knots().to_vec(), fitted.multiplicities().to_vec(), 3)).unwrap();
        let mut reversed = surface.clone();
        reversed.u_reverse();
        for i in 0..=20 {
            let u = 6.5 * i as f64 / 20.0;
            assert!(reversed.point(surface.u_reversed_parameter(u), 0.3).distance(surface.point(u, 0.3)) < 1e-12);
        }
        let mut reversed = swapped.clone();
        reversed.v_reverse();
        for i in 0..=20 {
            let v = 6.5 * i as f64 / 20.0;
            assert!(reversed.point(0.3, swapped.v_reversed_parameter(v)).distance(swapped.point(0.3, v)) < 1e-12);
        }
    }

    #[test]
    fn smoothed_approximation_passes_near_the_points() {
        // Three turns of a helix (5 across, 4 high) sampled densely: within the tolerance of
        // every point, of the greatest degree, with far fewer poles than points; its length
        // the helix's (OCCT's, from the same points: 94.332624).
        let at = |t: f64| Point::new(5.0 * (6.0 * PI * t).cos(), 5.0 * (6.0 * PI * t).sin(), 4.0 * t);
        let points: Vec<Point> = (0..=400).map(|i| at(i as f64 / 400.0)).collect();
        let options = ApproximationOptions { min_degree: 1, max_degree: 6, tolerance: 1e-6, smoothing: Some([1.0, 1.0, 1.0]), ..Default::default() };
        let curve = approximate_with(&points, options).unwrap();
        assert_eq!(curve.degree(), 6);
        assert!(curve.poles().len() < 100, "{} poles", curve.poles().len());
        for q in &points {
            let found = crate::geometry::extrema::nearest_on_curve(*q, &crate::geometry::curve::Curve::BSpline(curve.clone()), 0.0, 1.0);
            assert!(found.distance <= 1e-6, "{q:?} {}", found.distance);
        }
        let length = crate::geometry::measure::length(&crate::geometry::curve::Curve::BSpline(curve));
        assert!((length - 94.332624).abs() < 1e-5, "{length}");
    }

    /// A grid of five points one way: the surface passes through every point exactly (one
    /// span of degree four that way), as OCCT's does, so that faces fitted to grids sharing
    /// a row of points meet along it (gear teeth's faces, sewn).
    #[test]
    fn short_ways_of_grids_are_interpolated() {
        let at = |u: f64, v: f64| Point::new(u, 0.5 * (0.6 * v).sin() + 0.3 * u, v + 0.2 * (u * 1.3).cos());
        let grid: Vec<Vec<Point>> = (0..12).map(|i| (0..5).map(|j| at(-4.0 + 4.0 * i as f64 / 11.0, j as f64)).collect()).collect();
        let options = ApproximationOptions { min_degree: 3, max_degree: 8, tolerance: 1e-2, ..Default::default() };
        let surface = Surface::BSpline(approximate_surface(&grid, options).unwrap());
        for row in [&grid[0], &grid[11]] {
            for p in row {
                let (_, _, d) = crate::geometry::extrema::nearest_on_surface(*p, &surface, None);
                assert!(d < 1e-9, "{d} from {p:?}");
            }
        }
    }

    #[test]
    fn surfaces_pass_near_their_grids() {
        use crate::geometry::ParametricSurface;
        // A saddle sampled on a 30 x 20 grid: within the tolerance at every point, between
        // them too, with fewer poles than points.
        let at = |x: f64, y: f64| Point::new(x, y, 0.3 * x * y - 0.1 * x * x + (y * 0.7).sin());
        let points: Vec<Vec<Point>> = (0..30).map(|i| (0..20).map(|j| at(i as f64 * 0.2, j as f64 * 0.25)).collect()).collect();
        let options = ApproximationOptions { min_degree: 3, max_degree: 8, tolerance: 1e-3, ..Default::default() };
        let surface = approximate_surface(&points, options).unwrap();
        assert!(surface.poles().len() < 30 && surface.poles()[0].len() < 20, "{} x {} poles", surface.poles().len(), surface.poles()[0].len());
        for row in &points {
            for &p in row {
                let (u, v, _) = crate::geometry::extrema::nearest_on_surface(p, &Surface::BSpline(surface.clone()), None);
                assert!(surface.point(u, v).distance(p) <= 1e-3, "{p:?}");
            }
        }
        let mid = at(2.9, 2.1);
        let (u, v, _) = crate::geometry::extrema::nearest_on_surface(mid, &Surface::BSpline(surface.clone()), None);
        assert!(surface.point(u, v).distance(mid) <= 1e-2);
    }

    #[test]
    fn interpolated_circle_stays_close() {
        let points: Vec<Point> = (0..=32)
            .map(|i| {
                let a = std::f64::consts::PI * i as f64 / 32.0;
                Point::new(a.cos(), a.sin(), 0.0)
            })
            .collect();
        let curve = interpolate(&points, false, 1e-9).unwrap();
        assert_eq!(curve.degree(), 3);
        let (t0, t1) = (curve.first_parameter(), curve.last_parameter());
        for i in 0..=200 {
            let p = curve.point(t0 + (t1 - t0) * i as f64 / 200.0);
            assert!(((p - Point::ORIGIN).length() - 1.0).abs() < 1e-5, "{p:?}");
        }
        for (p, q) in [(curve.point(t0), points[0]), (curve.point(t1), points[32])] {
            assert!(p.distance(q) < 1e-12);
        }
        // Derivatives against finite differences.
        let t = 0.3 * (t1 - t0);
        let (_, d1, d2) = curve.d2(t);
        let h = 1e-5;
        let fd1 = (curve.point(t + h) - curve.point(t - h)) * (0.5 / h);
        assert!((d1 - fd1).length() < 1e-6);
        let (_, d1h) = curve.d1(t + h);
        let (_, d1m) = curve.d1(t - h);
        assert!((d2 - (d1h - d1m) * (0.5 / h)).length() < 1e-5);
    }

    #[test]
    fn rational_quarter_circle_is_exact() {
        // The classic quadratic rational arc.
        let w = std::f64::consts::FRAC_1_SQRT_2;
        let curve = BSplineCurve::new(
            vec![Point::new(1.0, 0.0, 0.0), Point::new(1.0, 1.0, 0.0), Point::new(0.0, 1.0, 0.0)],
            vec![0.0, 1.0],
            vec![3, 3],
            2,
        )
        .with_weights(vec![1.0, w, 1.0]);
        for i in 0..=10 {
            let p = curve.point(i as f64 / 10.0);
            assert!(((p - Point::ORIGIN).length() - 1.0).abs() < 1e-14);
        }
    }
}
