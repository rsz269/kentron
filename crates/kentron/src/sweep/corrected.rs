//! OCCT's corrected Frenet trihedron (`GeomFill_CorrectedFrenet`), the frame its pipes carry
//! a profile in by default: the spine's Frenet frame turned about the tangent by an angle
//! law, which takes the Frenet frame's turning about the tangent out.
//!
//! The law is found as OCCT finds it, and is not the rotation-minimizing frame it stands
//! for: along each interval where the spine is C2, in steps (of at most a tenth of the
//! spine, a third of the interval, and the step the spine's curvature allows), each step's
//! angle the Frenet normal's turn against the previous one carried by the least rotation
//! taking the previous tangent to this one; the angles summed, interpolated by a cubic
//! B-spline in the parameter (`Law_Interpolate`), and smoothed by removing knots within a
//! tenth of a radian (`smoothlaw`). Over a helix those steps are some way from parallel
//! transport, and OCCT's pipes turn their profiles accordingly; kentron's follow.
//!
//! Only spines OCCT corrects: curves not analytic, not planar, with no point where the
//! curvature vanishes (there OCCT's Frenet frame has its own rules, not ported); elsewhere
//! the rotation-minimizing frame, which equals OCCT's there.

use crate::geometry::curve::Curve;
use crate::geometry::ParametricCurve;
use crate::math::{Point, Vector};

/// `Precision::PConfusion()` and `Precision::Confusion()`.
const P_CONFUSION: f64 = 1e-9;
const CONFUSION: f64 = 1e-7;
/// `Precision::Angular()`.
const ANGULAR: f64 = 1e-12;

/// `gp_Dir::Angle`: the angle between two vectors (normalized first), in `[0, π]`.
fn angle(a: Vector, b: Vector) -> f64 {
    let (a, b) = (a * (1.0 / a.length()), b * (1.0 / b.length()));
    let cos = a.dot(b);
    if cos > -std::f64::consts::FRAC_1_SQRT_2 && cos < std::f64::consts::FRAC_1_SQRT_2 {
        cos.acos()
    } else {
        let sin = a.cross(b).length().min(1.0);
        if cos < 0.0 { std::f64::consts::PI - sin.asin() } else { sin.asin() }
    }
}

/// `v` turned by `angle` about the unit vector `axis` (OCCT's `SetLinearForm` of it).
fn turned(v: Vector, axis: Vector, angle: f64) -> Vector {
    let cross = axis.cross(v);
    v + cross * angle.sin() + axis.cross(cross) * (1.0 - angle.cos())
}

/// A curve as the law sees it: run over `[0, length]` from `first` toward `last`.
struct Spine<'a> {
    curve: &'a Curve,
    first: f64,
    sign: f64,
}

impl Spine<'_> {
    fn u(&self, s: f64) -> f64 {
        self.first + self.sign * s
    }

    /// Derivatives up to the third along `s`.
    fn d3(&self, s: f64) -> (Vector, Vector, Vector) {
        let (_, d1, d2, d3) = self.curve.d3(self.u(s));
        (d1 * self.sign, d2, d3 * self.sign)
    }

    /// The Frenet frame's tangent and normal (`GeomFill_Frenet::D0`, away from singular
    /// points).
    fn frenet(&self, s: f64) -> (Vector, Vector) {
        let (d1, d2, _) = self.d3(s);
        let t = d1 * (1.0 / d1.length());
        let b = t.cross(d2);
        let b = if b.length() <= 1e-290 { default_y(t) } else { b * (1.0 / b.length()) };
        (t, b.cross(t))
    }
}

/// `gp_Ax2(P, V).YDirection()`: the Y direction OCCT gives an axis system about `v`.
fn default_y(v: Vector) -> Vector {
    let (a, b, c) = (v.x.abs(), v.y.abs(), v.z.abs());
    // (gp_Ax2's X direction: the axis of the smallest component dropped.)
    let x = if b <= a && b <= c {
        if a > c { Vector::new(-v.z, 0.0, v.x) } else { Vector::new(v.z, 0.0, -v.x) }
    } else if a <= b && a <= c {
        if b > c { Vector::new(0.0, -v.z, v.y) } else { Vector::new(0.0, v.z, -v.y) }
    } else if a > b {
        Vector::new(-v.y, v.x, 0.0)
    } else {
        Vector::new(v.y, -v.x, 0.0)
    };
    let x = x * (1.0 / x.length());
    let z = v * (1.0 / v.length());
    z.cross(x)
}

/// `GeomFill_CorrectedFrenet::CalcAngleAT`: the turn about the tangent from the frame
/// `(pt, pn)` to `(t, n)`, `n` first carried back by the least rotation taking `t` to `pt`.
fn turn_between(t: Vector, n: Vector, pt: Vector, pn: Vector) -> f64 {
    let a = angle(t, pt);
    let rotated = if a > ANGULAR && a < std::f64::consts::PI - ANGULAR {
        let axis = t.cross(pt);
        let axis = axis * (1.0 / axis.length());
        turned(n, axis, a)
    } else {
        n
    };
    let mut at = angle(rotated, pn);
    if at > ANGULAR && std::f64::consts::PI - at > ANGULAR {
        let c = rotated.cross(pn);
        if c.length() > 0.0 && std::f64::consts::PI - angle(c, pt) <= ANGULAR {
            at = -at;
        }
    }
    at
}

/// A scalar B-spline (OCCT's `Law_BSpline`, not rational, not periodic).
#[derive(Clone, Debug)]
struct Law {
    degree: usize,
    knots: Vec<f64>,
    mults: Vec<usize>,
    poles: Vec<f64>,
}

impl Law {
    fn flat(&self) -> Vec<f64> {
        self.knots.iter().zip(&self.mults).flat_map(|(&k, &m)| std::iter::repeat_n(k, m)).collect()
    }

    /// The value at `t` (de Boor, in the span holding `t`, the last for the end).
    fn value(&self, t: f64) -> f64 {
        let flat = self.flat();
        let p = self.degree;
        let n = self.poles.len();
        let t = t.clamp(flat[p], flat[n]);
        let mut k = p;
        while k + 1 < n && flat[k + 1] <= t {
            k += 1;
        }
        let mut d: Vec<f64> = (0..=p).map(|j| self.poles[j + k - p]).collect();
        for r in 1..=p {
            for j in (r..=p).rev() {
                let i = j + k - p;
                let denom = flat[i + p + 1 - r] - flat[i];
                let alpha = if denom != 0.0 { (t - flat[i]) / denom } else { 0.0 };
                d[j] = (1.0 - alpha) * d[j - 1] + alpha * d[j];
            }
        }
        d[p]
    }

    /// The basis functions' values (`order` 0) or first derivatives (`order` 1) at `t`, all
    /// `n` of them, over the flat knots.
    fn basis(flat: &[f64], degree: usize, n: usize, t: f64, order: usize) -> Vec<f64> {
        // Cox-de Boor over every function (small systems: fine).
        let last = t >= flat[n];
        let mut b: Vec<f64> = (0..flat.len() - 1)
            .map(|i| {
                let inside = flat[i] <= t && t < flat[i + 1];
                // (At the end, the last span not empty.)
                let at_end = last && flat[i] < flat[i + 1] && flat[i + 1] == flat[n];
                if inside || at_end { 1.0 } else { 0.0 }
            })
            .collect();
        let mut prev = b.clone();
        for p in 1..=degree {
            prev = b.clone();
            for i in 0..flat.len() - 1 - p {
                let left = if flat[i + p] > flat[i] { (t - flat[i]) / (flat[i + p] - flat[i]) * prev[i] } else { 0.0 };
                let right = if flat[i + p + 1] > flat[i + 1] { (flat[i + p + 1] - t) / (flat[i + p + 1] - flat[i + 1]) * prev[i + 1] } else { 0.0 };
                b[i] = left + right;
            }
        }
        if order == 0 {
            return b[..n].to_vec();
        }
        // Derivatives from the degree - 1 functions.
        (0..n)
            .map(|i| {
                let p = degree as f64;
                let a = if flat[i + degree] > flat[i] { p / (flat[i + degree] - flat[i]) * prev[i] } else { 0.0 };
                let c = if flat[i + degree + 1] > flat[i + 1] { p / (flat[i + degree + 1] - flat[i + 1]) * prev[i + 1] } else { 0.0 };
                a - c
            })
            .collect()
    }

    /// `Law_Interpolate` (not periodic, no tangents given): through `points` at `params`;
    /// cubic with the end slopes of the Lagrange cubics through the first and last four.
    fn interpolate(points: &[f64], params: &[f64]) -> Option<Law> {
        let n = points.len();
        match n {
            0 | 1 => None,
            2 => Some(Law { degree: 1, knots: params.to_vec(), mults: vec![2, 2], poles: points.to_vec() }),
            3 => {
                let flat = vec![params[0], params[0], params[0], params[2], params[2], params[2]];
                let rows: Vec<Vec<f64>> = params.iter().map(|&t| Law::basis(&flat, 2, 3, t, 0)).collect();
                let poles = solve(rows, points.to_vec())?;
                Some(Law { degree: 2, knots: vec![params[0], params[2]], mults: vec![3, 3], poles })
            }
            _ => {
                let slope = |ts: &[f64], ys: &[f64], at: f64| -> f64 {
                    // (The derivative at `at` of the Lagrange polynomial through the points.)
                    let mut d = 0.0;
                    for i in 0..ts.len() {
                        let mut term = 0.0;
                        for j in 0..ts.len() {
                            if j == i {
                                continue;
                            }
                            let mut prod = 1.0 / (ts[i] - ts[j]);
                            for k in 0..ts.len() {
                                if k != i && k != j {
                                    prod *= (at - ts[k]) / (ts[i] - ts[k]);
                                }
                            }
                            term += prod;
                        }
                        d += ys[i] * term;
                    }
                    d
                };
                let start = slope(&params[..4], &points[..4], params[0]);
                let end = slope(&params[n - 4..], &points[n - 4..], params[n - 1]);
                let m = n + 2;
                let mut flat = vec![params[0]; 4];
                flat.extend_from_slice(&params[1..n - 1]);
                flat.extend(std::iter::repeat_n(params[n - 1], 4));
                let mut rows = Vec::with_capacity(m);
                let mut rhs = Vec::with_capacity(m);
                rows.push(Law::basis(&flat, 3, m, params[0], 0));
                rhs.push(points[0]);
                rows.push(Law::basis(&flat, 3, m, params[0], 1));
                rhs.push(start);
                for i in 1..n - 1 {
                    rows.push(Law::basis(&flat, 3, m, params[i], 0));
                    rhs.push(points[i]);
                }
                rows.push(Law::basis(&flat, 3, m, params[n - 1], 1));
                rhs.push(end);
                rows.push(Law::basis(&flat, 3, m, params[n - 1], 0));
                rhs.push(points[n - 1]);
                let poles = solve(rows, rhs)?;
                let mut mults = vec![1; n];
                mults[0] = 4;
                mults[n - 1] = 4;
                Some(Law { degree: 3, knots: params.to_vec(), mults, poles })
            }
        }
    }

    /// `Law_BSpline::RemoveKnot(index, 0, tolerance)` (`BSplCLib::RemoveKnot` and its
    /// anti-Boor scheme): knot `index` (from 0, interior) taken out if the curve without it
    /// stays within `tolerance`.
    fn remove_knot(&mut self, index: usize, tolerance: f64) -> bool {
        let degree = self.degree;
        if index == 0 || index + 1 >= self.knots.len() || self.mults[index] == 0 {
            return false;
        }
        let depth = self.mults[index];
        let length = degree;
        if depth != 1 {
            return false;
        }
        // The knots around: the `degree` ending at knots[index - 1], then the `degree`
        // starting at knots[index + 1] (BuildKnots of the spans either side).
        let flat = self.flat();
        let last_of = |k: usize| -> usize { self.mults[..=k].iter().sum::<usize>() - 1 };
        let a = last_of(index - 1);
        let b = last_of(index) + 1;
        let mut knots = Vec::with_capacity(2 * degree);
        for i in 0..degree {
            knots.push(flat[a + 1 + i - degree]);
        }
        for i in 0..degree {
            knots.push(flat[b + i]);
        }
        // The poles of the anti-Boor scheme, at BoorIndex places.
        let pole_index = self.mults[..index].iter().sum::<usize>() - (degree + 1);
        let boor = |i: usize| -> usize {
            if i <= depth {
                i
            } else if i <= length {
                2 * i - depth
            } else {
                length + i - depth
            }
        };
        let mut local = vec![0.0; 2 * degree + 2];
        for i in 0..=length + depth {
            local[boor(i)] = self.poles[pole_index + i];
        }
        let u = self.knots[index];
        // (Offsets into `local` kept signed, as OCCT's pointers.)
        let at = |i: isize| -> usize { i as usize };
        if length == 1 {
            let x = (knots[degree] - u) / (knots[degree] - knots[0]);
            let y = 1.0 - x;
            let first = (depth - 1) as isize;
            let z = x * local[at(first)] + y * local[at(first + 2)];
            if (z - local[at(first + 1)]).abs() > tolerance {
                return false;
            }
        } else {
            let mut firstpole = (depth - 1) as isize;
            for step in (0..depth).rev() {
                firstpole -= 1;
                let mut pole = firstpole;
                for i in step..length - 1 {
                    pole += 2;
                    let x = (knots[i + degree - step] - u) / (knots[i + degree - step] - knots[i]);
                    let y = 1.0 - x;
                    local[at(pole + 1)] = (local[at(pole)] - x * local[at(pole - 1)]) / y;
                }
                pole += 4;
                let half = (length - 1 + step) / 2;
                let mut i = length - 1;
                while i > half {
                    pole -= 2;
                    let x = (knots[i + degree - step] - u) / (knots[i + degree - step] - knots[i]);
                    let y = 1.0 - x;
                    let z = (local[at(pole)] - y * local[at(pole + 1)]) / x;
                    if (z - local[at(pole - 1)]).abs() > tolerance {
                        return false;
                    }
                    local[at(pole - 1)] = (local[at(pole - 1)] + z) / 2.0;
                    i -= 1;
                }
            }
        }
        // The new poles: those before, the scheme's (at BoorIndex(i, length, 0) = 2 i), those
        // after.
        let mut poles = self.poles[..pole_index + 1].to_vec();
        for i in 1..=length {
            poles.push(local[2 * i]);
        }
        poles.extend_from_slice(&self.poles[pole_index + 1 + length + depth..]);
        self.poles = poles;
        self.knots.remove(index);
        self.mults.remove(index);
        true
    }

    /// `smoothlaw`: knots removed in two passes (within a tenth of `tolerance`, then what is
    /// left of it), each kept only if the law stays within `tolerance` of the points.
    fn smoothed(self, points: &[f64], params: &[f64], tolerance: f64) -> Law {
        let mut law = self.clone();
        let mut bs = self;
        let mut tol = tolerance / 10.0;
        let mut ok = false;
        for ii in (1..bs.knots.len().saturating_sub(1)).rev() {
            if bs.remove_knot(ii, tol) {
                ok = true;
            }
        }
        let worst = |bs: &Law| params.iter().zip(points).map(|(&t, &y)| (bs.value(t) - y).abs()).fold(0.0, f64::max);
        if ok {
            let w = worst(&bs);
            if w <= tolerance {
                tol = (tolerance - w) / 2.0;
                law = bs.clone();
            } else {
                return law;
            }
        } else {
            tol = tolerance / 2.0;
        }
        let mut ok = false;
        for ii in (1..bs.knots.len().saturating_sub(1)).rev() {
            if bs.remove_knot(ii, tol) {
                ok = true;
            }
        }
        if ok && worst(&bs) <= tolerance {
            law = bs;
        }
        law
    }
}

/// Gaussian elimination with partial pivoting.
fn solve(mut a: Vec<Vec<f64>>, mut b: Vec<f64>) -> Option<Vec<f64>> {
    let n = b.len();
    for c in 0..n {
        let p = (c..n).max_by(|&i, &j| a[i][c].abs().total_cmp(&a[j][c].abs()))?;
        if a[p][c].abs() < 1e-300 {
            return None;
        }
        a.swap(c, p);
        b.swap(c, p);
        for r in c + 1..n {
            let f = a[r][c] / a[c][c];
            if f != 0.0 {
                let (above, below) = a.split_at_mut(r);
                for (x, y) in below[0][c..].iter_mut().zip(&above[c][c..]) {
                    *x -= f * y;
                }
                b[r] -= f * b[c];
            }
        }
    }
    let mut x = vec![0.0; n];
    for r in (0..n).rev() {
        let s: f64 = (r + 1..n).map(|k| a[r][k] * x[k]).sum();
        x[r] = (b[r] - s) / a[r][r];
    }
    Some(x)
}

/// The angle law over one interval (ending at its first field): constant, or a B-spline.
enum Piece {
    Constant(f64, f64),
    Spline(f64, Law),
}

/// OCCT's corrected Frenet frames along one spine piece.
pub(super) struct CorrectedFrenet {
    curve: Curve,
    first: f64,
    sign: f64,
    pieces: Vec<Piece>,
    /// Where each step started: its parameter, the angle there, the Frenet tangent and
    /// normal (the arrays OCCT added for its bug OCC78, for `angle_at`).
    poles: Vec<(f64, f64, Vector, Vector)>,
}

impl CorrectedFrenet {
    /// The law along `curve` from `first` to `last` (either way); `None` where OCCT keeps the
    /// Frenet frame (analytic or planar spines: equal to the rotation-minimizing one) or
    /// where the curvature vanishes (not ported).
    pub(super) fn new(curve: &Curve, first: f64, last: f64) -> Option<CorrectedFrenet> {
        if !matches!(curve, Curve::BSpline(_) | Curve::Bezier(_) | Curve::Trimmed(_) | Curve::Offset(_)) {
            return None;
        }
        let length = (last - first).abs();
        if length.is_nan() || length <= 0.0 {
            return None;
        }
        let sign = if last >= first { 1.0 } else { -1.0 };
        let spine = Spine { curve, first, sign };
        if planar(curve, first, last) || vanishing_curvature(&spine, length) {
            return None;
        }
        // The intervals where the spine is C2 (GeomFill_Frenet's for C0).
        let mut breaks = vec![0.0];
        if let Curve::BSpline(b) = curve {
            let degree = b.degree();
            for (&k, &m) in b.knots().iter().zip(b.multiplicities()) {
                if m + 2 > degree {
                    let s = (k - first) * sign;
                    if s > P_CONFUSION && s < length - P_CONFUSION {
                        breaks.push(s);
                    }
                }
            }
            breaks.sort_by(f64::total_cmp);
        }
        breaks.push(length);
        let (mut prev_t, mut prev_n) = spine.frenet(0.0);
        let mut start_angle = 0.0;
        let average_step = length / 10.0;
        let mut pieces = Vec::with_capacity(breaks.len() - 1);
        let mut poles = Vec::new();
        let mut zero = true;
        for w in breaks.windows(2) {
            let (a, b) = (w[0], w[1]);
            let steps = (((b - a) / average_step) as usize).max(3);
            let step = (b - a) / steps as f64;
            let (piece, is_zero) = interval(&spine, a, b, step, &mut start_angle, &mut prev_t, &mut prev_n, &mut poles)?;
            zero &= is_zero;
            pieces.push(piece);
        }
        if zero {
            return None;
        }
        Some(CorrectedFrenet { curve: curve.clone(), first, sign, pieces, poles })
    }

    fn spine(&self) -> Spine<'_> {
        Spine { curve: &self.curve, first: self.first, sign: self.sign }
    }

    /// The angle law's value at `s` (`Law_Composite`: the interval holding it, the last
    /// past the end).
    fn law(&self, s: f64) -> f64 {
        let end = |p: &Piece| match p {
            Piece::Constant(b, _) | Piece::Spline(b, _) => *b,
        };
        let piece = self.pieces.iter().find(|p| s <= end(p)).or(self.pieces.last());
        match piece {
            Some(Piece::Constant(_, v)) => *v,
            Some(Piece::Spline(_, law)) => law.value(s),
            None => 0.0,
        }
    }

    /// `GetAngleAT`: the law's angle, unless it strays more than a quarter turn from the
    /// turn measured from the step before (OCCT's fix for its bug OCC78).
    fn angle_at(&self, s: f64) -> f64 {
        let poles = &self.poles;
        let n = poles.len();
        if n == 0 || s == poles[0].0 {
            return self.law(s);
        }
        let (mut ib, mut ie) = (0usize, n - 1);
        let mut ic = (ib + ie) / 2;
        if s > poles[ie].0 {
            ic = ie;
        }
        if ic < ie {
            while poles[ic].0 > s || s > poles[ic + 1].0 {
                if poles[ic].0 < s {
                    ib = ic;
                } else {
                    ie = ic;
                }
                let next = (ie + ib) / 2;
                if next == ic {
                    break;
                }
                ic = next;
            }
            if poles[ic].0 == s || s == poles[ic + 1].0 {
                return self.law(s);
            }
        }
        let mut angle = self.law(s);
        let (_, angle_o, pt, pn) = poles[ic];
        let d_ang = angle - angle_o;
        let (t, nrm) = self.spine().frenet(s);
        let measured = turn_between(t, nrm, pt, pn);
        if diff_angle(measured, d_ang).abs() > std::f64::consts::FRAC_PI_2 {
            angle = angle_o + measured;
        }
        angle
    }

    /// The frame's normal at fraction `fraction` of the piece (from its `first` toward its
    /// `last`): the Frenet normal turned by the law's angle about the tangent.
    pub(super) fn normal(&self, fraction: f64, length: f64) -> Vector {
        let s = fraction * length;
        let (t, n) = self.spine().frenet(s);
        turned(n, t, self.angle_at(s))
    }
}

/// An angle in `[0, 2π)` brought into `[-π, π)`.
fn corr(a: f64) -> f64 {
    if a < std::f64::consts::PI { a } else { a - 2.0 * std::f64::consts::PI }
}

/// `diffAng`: `a - o` brought into `(-π, π]`.
fn diff_angle(a: f64, o: f64) -> f64 {
    let tau = 2.0 * std::f64::consts::PI;
    let d = (a - o) - ((a - o) / tau).floor() * tau;
    if d >= 0.0 { corr(d) } else { -corr(-d) }
}

/// `InitInterval`: the angles over one interval, from `a` to `b` in steps of at most
/// `step`, continuing from the frame `(prev_t, prev_n)` and angle `start_angle` the last one
/// ended with; the law over it, and whether every angle was nought.
#[allow(clippy::too_many_arguments)]
fn interval(spine: &Spine, a: f64, b: f64, step: f64, start_angle: &mut f64, prev_t: &mut Vector, prev_n: &mut Vector, poles: &mut Vec<(f64, f64, Vector, Vector)>) -> Option<(Piece, bool)> {
    const LENGTH_MIN: f64 = 1e-2 * 1e-4;
    let mut params = Vec::new();
    let mut angles: Vec<f64> = Vec::new();
    let (mut zero, mut constant) = (true, true);
    let mut param = a;
    let mut current = param;
    let mut current_step;
    let d_last = b - P_CONFUSION;
    let mut angle_at = 0.0;
    let mut i = 1;
    let mut guard = 0;
    while param < b {
        guard += 1;
        if guard > 1_000_000 {
            return None;
        }
        if current > d_last {
            if (d_last - param).abs() < CONFUSION * CONFUSION {
                param = current;
            }
            current = b;
        }
        let (t, n) = spine.frenet(current);
        if angle(*prev_t, t) < std::f64::consts::FRAC_PI_3 || i == 1 {
            params.push(current);
            poles.push((param, if i > 1 { angles[i - 2] } else { *start_angle }, *prev_t, *prev_n));
            angle_at = turn_between(t, n, *prev_t, *prev_n);
            if constant && i > 1 && angle_at.abs() > P_CONFUSION {
                constant = false;
            }
            angle_at += if i > 1 { angles[i - 2] } else { *start_angle };
            angles.push(angle_at);
            *prev_n = n;
            if zero && angle_at.abs() > P_CONFUSION {
                zero = false;
            }
            *prev_t = t;
            param = current;
            i += 1;
            // The next step: as the curvature (the singular function D1 × D2) changes.
            let (d1, d2, d3) = spine.d3(param);
            let l = (d1.cross(d2).length() / 2.0).max(LENGTH_MIN);
            let norm = d1.cross(d3).length().max(CONFUSION);
            current_step = (l / norm).min(step);
        } else {
            current_step = (current - param) / 2.0;
        }
        current = param + current_step;
    }
    *start_angle = angle_at;
    let piece = if constant {
        Piece::Constant(b, angle_at)
    } else {
        let law = Law::interpolate(&angles, &params)?;
        Piece::Spline(b, law.smoothed(&angles, &params, 0.1))
    };
    Some((piece, zero))
}

/// `FindPlane`: the spine lies in a plane (its poles, a B-spline's or Bezier's, else 15
/// samples, within `Precision::Confusion` of their plane of inertia).
fn planar(curve: &Curve, first: f64, last: f64) -> bool {
    let points: Vec<Point> = match curve {
        Curve::BSpline(b) => b.poles().to_vec(),
        Curve::Bezier(b) => b.poles().to_vec(),
        _ => (0..15).map(|i| curve.point(first + (last - first) * i as f64 / 14.0)).collect(),
    };
    if points.len() < 3 {
        return false;
    }
    let n = points.len() as f64;
    let c = points.iter().fold(Vector::new(0.0, 0.0, 0.0), |s, p| s + Vector::new(p.x, p.y, p.z)) * (1.0 / n);
    let mut m = [[0.0; 3]; 3];
    for p in &points {
        let d = [p.x - c.x, p.y - c.y, p.z - c.z];
        for i in 0..3 {
            for j in 0..3 {
                m[i][j] += d[i] * d[j];
            }
        }
    }
    let matrix = crate::math::Matrix3::from_rows(m);
    let (values, vectors) = crate::properties::symmetric_eigen(&matrix);
    // (All the points on one line: no plane.)
    if values[1] <= 1e-30 * values[2].max(1e-300) {
        return false;
    }
    let normal = vectors[0];
    points.iter().all(|p| ((*p - Point::new(c.x, c.y, c.z)).dot(normal)).abs() <= CONFUSION)
}

/// Points where the spine's curvature vanishes (D1 × D2 nought), found from samples.
fn vanishing_curvature(spine: &Spine, length: f64) -> bool {
    let values: Vec<f64> = (0..=256).map(|i| {
        let (d1, d2, _) = spine.d3(length * i as f64 / 256.0);
        d1.cross(d2).length()
    }).collect();
    let average = values.iter().sum::<f64>() / values.len() as f64;
    values.iter().any(|&v| v <= 1e-6 * average)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Knot removal leaves a cubic through points on a line (all its knots removable) the
    /// line, and interpolation passes through its points.
    #[test]
    fn laws_interpolate_and_smooth() {
        let params = [0.0, 0.3, 0.7, 1.2, 2.0, 2.4];
        let points: Vec<f64> = params.iter().map(|t| 2.0 * t - 1.0).collect();
        let law = Law::interpolate(&points, &params).unwrap();
        for (t, y) in params.iter().zip(&points) {
            assert!((law.value(*t) - y).abs() < 1e-12);
        }
        let smooth = law.smoothed(&points, &params, 0.1);
        assert_eq!(smooth.knots.len(), 2);
        assert!((smooth.value(1.0) - 1.0).abs() < 1e-12);
        let wavy: Vec<f64> = params.iter().map(|t: &f64| t.sin()).collect();
        let law = Law::interpolate(&wavy, &params).unwrap();
        for (t, y) in params.iter().zip(&wavy) {
            assert!((law.value(*t) - y).abs() < 1e-12);
        }
    }
}
