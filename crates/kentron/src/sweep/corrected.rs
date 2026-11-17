//! OCCT's corrected Frenet trihedron (`GeomFill_CorrectedFrenet`), ported: the frame its pipes
//! carry a profile in by default, the spine's Frenet frame (`GeomFill_Frenet`, see `frenet`)
//! turned about the tangent by an angle law that takes the Frenet frame's turning about the
//! tangent out.
//!
//! As OCCT: lines and conics keep their Frenet frame. Other curves get the law, found
//! interval by interval (the Frenet frame's intervals: where the curve is C2, split at its
//! singular points) in steps (of at most a tenth of the curve's range, a third of the
//! interval, and the step the curvature allows), each step's angle the Frenet normal's turn
//! against the previous one carried by the least rotation taking the previous tangent to this
//! one (`CalcAngleAT`); the angles summed and interpolated by a cubic B-spline in the
//! parameter (`Law_Interpolate`), smoothed by removing knots within a tenth of a radian
//! (`smoothlaw`), or held constant over an interval where they do not change or the curve is
//! planar (`FindPlane`, one step then). A curve whose angles all stay nought keeps the
//! Frenet frame. The law's value is checked against the turn measured from the step before
//! (`GetAngleAT`, OCCT's fix for its bug OCC78).

use crate::error::{Error, Result};
use crate::geometry::adaptor::{CurveAdaptor, CurveKind, GeomAdaptor};
use crate::geometry::curve::Curve;
use crate::geometry::{Continuity, ParametricCurve};
use crate::math::{Point, Vector};

use super::frenet::{angle, Frenet, SingularFunction, ANGULAR, CONFUSION, P_CONFUSION};

/// `v` turned by `angle` about the unit vector `axis` (OCCT's `SetLinearForm` of it).
fn turned(v: Vector, axis: Vector, angle: f64) -> Vector {
    let cross = axis.cross(v);
    v + cross * angle.sin() + axis.cross(cross) * (1.0 - angle.cos())
}

/// `GeomFill_CorrectedFrenet::CalcAngleAT`: the turn about the tangent from the frame
/// `(pt, pn)` to `(t, n)`, `n` first carried back by the least rotation taking `t` to `pt`.
fn turn_between(t: Vector, n: Vector, pt: Vector, pn: Vector) -> f64 {
    let a = angle(t, pt);
    let rotated = if a.abs() > ANGULAR && a.abs() < std::f64::consts::PI - ANGULAR {
        // (`Normalized()`: divided by its length.)
        let axis = t.cross(pt);
        let axis = axis / axis.length();
        turned(n, axis, a)
    } else {
        n
    };
    let mut at = angle(rotated, pn);
    // (`IsOpposite(prevTangent, Angular)`: the cross product is not nought here.)
    if at > ANGULAR && std::f64::consts::PI - at > ANGULAR && std::f64::consts::PI - angle(rotated.cross(pn), pt) <= ANGULAR {
        at = -at;
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

/// One interval's law (`Law_Constant` or `Law_BSpFunc`), with its bounds.
#[derive(Clone, Debug)]
enum IntervalLaw {
    Constant(f64),
    Spline(Law),
}

/// OCCT's corrected Frenet trihedron along one curve (`SetCurve` done).
#[derive(Clone, Debug)]
pub(super) struct CorrectedFrenet {
    frenet: Frenet,
    /// `isFrenet`: the Frenet frame as it is.
    is_frenet: bool,
    /// `EvolAroundT`, a `Law_Composite`: each interval's bounds and law.
    laws: Vec<(f64, f64, IntervalLaw)>,
    periodic: bool,
    /// `HArrPoles`, `HArrAngle`, `HArrTangent`, `HArrNormal` (OCC78): where each step
    /// started, the angle there, the Frenet tangent and normal there.
    poles: Vec<(f64, f64, Vector, Vector)>,
    /// The Frenet frame's intervals (`Intervals(T, GeomAbs_C0)`), the law's.
    intervals: Vec<f64>,
}

impl CorrectedFrenet {
    /// The trihedron law along `curve` run from `first` to `last` (`last` the smaller for an
    /// edge run backwards: the curve reversed then, as `BRepFill_Edge3DLaw` reverses a
    /// reversed edge's).
    pub(super) fn new(curve: &Curve, first: f64, last: f64) -> Result<CorrectedFrenet> {
        CorrectedFrenet::of(edge_adaptor(curve, first, last))
    }

    /// `SetCurve`.
    fn of(curve: GeomAdaptor) -> Result<CorrectedFrenet> {
        let frenet = Frenet::new(curve);
        let periodic = frenet.curve().is_periodic();
        let intervals = frenet.intervals(Continuity::C0);
        let mut law = CorrectedFrenet { frenet, is_frenet: true, laws: Vec::new(), periodic, poles: Vec::new(), intervals };
        if !matches!(law.frenet.curve().kind(), CurveKind::Circle | CurveKind::Ellipse | CurveKind::Hyperbola | CurveKind::Parabola | CurveKind::Line) {
            law.init()?;
        }
        Ok(law)
    }

    fn curve(&self) -> &GeomAdaptor {
        self.frenet.curve()
    }

    /// `Init`.
    fn init(&mut self) -> Result<()> {
        let curve = self.curve().clone();
        let (first, last) = (curve.first_parameter(), curve.last_parameter());
        let (_, (mut tangent, mut normal, _)) = self.frenet.d0(&curve, first);
        let mut start_angle = 0.0;
        let average_step = (last - first) / 10.0;
        let intervals = self.intervals.clone();
        for w in intervals.windows(2) {
            // (C++'s `int()`: toward nought.)
            let steps = (((w[1] - w[0]) / average_step) as i64).max(3);
            let step = (w[1] - w[0]) / steps as f64;
            let (law, is_zero) = self.init_interval(w[0], w[1], step, &mut start_angle, &mut tangent, &mut normal)?;
            if !is_zero {
                self.is_frenet = false;
            }
            self.laws.push((w[0], w[1], law));
        }
        Ok(())
    }

    /// `InitInterval`: the angles over `[first, last]`, in steps of at most `step`, carried on
    /// from the frame `(prev_tangent, prev_normal)` and angle `start_angle` the interval
    /// before ended with; the law over it, and whether every angle was nought.
    fn init_interval(&mut self, first: f64, last: f64, step: f64, start_angle: &mut f64, prev_tangent: &mut Vector, prev_normal: &mut Vector) -> Result<(IntervalLaw, bool)> {
        // (`frenet->SetInterval(First, Last)`.)
        let curve = self.curve().clone();
        let trimmed = curve.trim(first, last);
        let singular = SingularFunction { curve: &curve, ratio: 1.0 };
        // `aCurveBoundBox.GetGap() * 1e-4`: the box of the singular function (an "other"
        // curve to BndLib, its sampled bounds kept as values) enlarged by the tolerance
        // 1e-2 given, the gap that.
        const LENGTH_MIN: f64 = 1e-2 * 1e-4;
        let mut params: Vec<f64> = Vec::new();
        let mut angles: Vec<f64> = Vec::new();
        let mut param = first;
        let (mut is_zero, mut is_const) = (true, true);
        let mut angle_at = 0.0;
        let mut curr_step = step;
        let is_planar = find_plane(&curve);
        let mut i = 1;
        let mut curr_param = param;
        let d_last = last - P_CONFUSION;
        while param < last {
            if curr_param > d_last {
                if (d_last - param).abs() < CONFUSION * CONFUSION {
                    param = curr_param;
                }
                curr_step = d_last - param;
                curr_param = last;
            }
            if is_planar {
                curr_param = last;
            }
            let (_, (t, n, _)) = self.frenet.d0(&trimmed, curr_param);
            if angle(*prev_tangent, t) < std::f64::consts::FRAC_PI_3 || i == 1 {
                params.push(curr_param);
                let before = if i > 1 { angles[i - 2] } else { *start_angle };
                self.poles.push((param, before, *prev_tangent, *prev_normal));
                angle_at = turn_between(t, n, *prev_tangent, *prev_normal);
                if is_const && i > 1 && angle_at.abs() > P_CONFUSION {
                    is_const = false;
                }
                angle_at += before;
                angles.push(angle_at);
                *prev_normal = n;
                if is_zero && angle_at.abs() > P_CONFUSION {
                    is_zero = false;
                }
                *prev_tangent = t;
                param = curr_param;
                i += 1;
                // The next step, as the singular function changes.
                let (p, d1) = singular.d1(param);
                let l = (Vector::new(p.x, p.y, p.z).length() / 2.0).max(LENGTH_MIN);
                let norm = d1.length().max(CONFUSION);
                curr_step = (l / norm).min(step);
            } else {
                // (Step too long.)
                curr_step /= 2.0;
            }
            curr_param = param + curr_step;
        }
        *start_angle = angle_at;
        if is_const || is_planar {
            return Ok((IntervalLaw::Constant(angle_at), is_zero));
        }
        // `Law_Interpolate` (its parameters checked increasing: else it throws).
        if params.windows(2).any(|w| w[1] - w[0] < f64::MIN_POSITIVE || (w[1] - w[0]).is_nan()) {
            return Err(Error::Failed("a pipe's angle law at parameters not increasing".into()));
        }
        let law = Law::interpolate(&angles, &params).ok_or_else(|| Error::Failed("a pipe's angle law not interpolated".into()))?;
        Ok((IntervalLaw::Spline(law.smoothed(&angles, &params, 0.1)), is_zero))
    }

    /// `EvolAroundT->Value(W)` (`Law_Composite::Value`): the interval law holding `w`, the
    /// first of two at a bound between them (the one a sweep run forward holds there, as
    /// OCCT's `Prepare` keeps its current law while it holds the value).
    fn law_value(&self, w: f64) -> f64 {
        let (first, last) = (self.laws[0].0, self.laws[self.laws.len() - 1].1);
        let w = if self.periodic { crate::geometry::adaptor::in_period(w, first, last) } else { w };
        let index = if w <= first {
            0
        } else if w >= last {
            self.laws.len() - 1
        } else {
            self.laws.iter().position(|(f, l, _)| *f <= w && w <= *l).unwrap_or(0)
        };
        match &self.laws[index].2 {
            IntervalLaw::Constant(v) => *v,
            IntervalLaw::Spline(law) => law.value(w),
        }
    }

    /// `GetAngleAT`: the law's angle, unless it strays more than a quarter turn from the
    /// turn measured from the step before (OCC78).
    fn angle_at(&self, trimmed: &GeomAdaptor, param: f64) -> f64 {
        let poles = &self.poles;
        let (mut ib, mut ie) = (0usize, poles.len() - 1);
        let mut ic = (ie + ib) / 2;
        if param == poles[ib].0 {
            return self.law_value(param);
        }
        if param > poles[ie].0 {
            ic = ie;
        }
        if ic < ie {
            while poles[ic].0 > param || param > poles[ic + 1].0 {
                if poles[ic].0 < param {
                    ib = ic;
                } else {
                    ie = ic;
                }
                let next = (ie + ib) / 2;
                // (Where OCCT's search would go round for ever: not reached for the
                // parameters a pipe asks for.)
                if next == ic {
                    break;
                }
                ic = next;
            }
            if poles[ic].0 == param || param == poles[ic + 1].0 {
                return self.law_value(param);
            }
        }
        let mut ang_p = self.law_value(param);
        let (_, ang_po, pt, pn) = poles[ic];
        let d_ang_law = ang_p - ang_po;
        let (_, (t, n, _)) = self.frenet.d0(trimmed, param);
        let d_ang = turn_between(t, n, pt, pn);
        if diff_angle(d_ang, d_ang_law).abs() > std::f64::consts::FRAC_PI_2 {
            ang_p = ang_po + d_ang;
        }
        ang_p
    }

    /// `D0` at `param`, evaluated as a sweep does after `SetInterval` to the law's interval
    /// holding it (the first of two at a bound): tangent and normal.
    pub(super) fn d0(&self, param: f64) -> (Vector, Vector) {
        let k = self.intervals.windows(2).position(|w| w[0] <= param && param <= w[1]).unwrap_or(if param < self.intervals[0] { 0 } else { self.intervals.len() - 2 });
        let trimmed = self.curve().trim(self.intervals[k], self.intervals[k + 1]);
        let (_, (t, n, _)) = self.frenet.d0(&trimmed, param);
        if self.is_frenet {
            return (t, n);
        }
        // (`Normal.SetLinearForm(sin, T × N, 1 - cos, T × (T × N), Normal)`: summed in that
        // order.)
        let angle_at = self.angle_at(&trimmed, param);
        let cross = t.cross(n);
        let (sin, one_less_cos, twice) = (angle_at.sin(), 1.0 - angle_at.cos(), t.cross(cross));
        let normal = Vector::new(sin * cross.x + one_less_cos * twice.x + n.x, sin * cross.y + one_less_cos * twice.y + n.y, sin * cross.z + one_less_cos * twice.z + n.z);
        (t, normal)
    }

    /// The frame's tangent and normal at fraction `fraction` of the way along.
    pub(super) fn at_fraction(&self, fraction: f64) -> (Vector, Vector) {
        let curve = self.curve();
        let (first, last) = (curve.first_parameter(), curve.last_parameter());
        self.d0(first + (last - first) * fraction)
    }
}

/// `GeomFill_Frenet` alone (a pipe's Frenet mode), along `curve` run from `first` to `last`.
#[derive(Clone, Debug)]
pub(super) struct FrenetLaw {
    frenet: Frenet,
}

impl FrenetLaw {
    pub(super) fn new(curve: &Curve, first: f64, last: f64) -> FrenetLaw {
        FrenetLaw { frenet: Frenet::new(edge_adaptor(curve, first, last)) }
    }

    /// The frame's tangent and normal at fraction `fraction` of the way along (evaluated
    /// in the Frenet frame's interval holding it, as a sweep sets it).
    pub(super) fn at_fraction(&self, fraction: f64) -> (Vector, Vector) {
        let curve = self.frenet.curve();
        let (first, last) = (curve.first_parameter(), curve.last_parameter());
        let param = first + (last - first) * fraction;
        let intervals = self.frenet.intervals(Continuity::C0);
        let k = intervals.windows(2).position(|w| w[0] <= param && param <= w[1]).unwrap_or(if param < intervals[0] { 0 } else { intervals.len() - 2 });
        let (_, (t, n, _)) = self.frenet.d0(&curve.trim(intervals[k], intervals[k + 1]), param);
        (t, n)
    }
}

/// An edge's curve as `BRepFill_Edge3DLaw` adapts it: over `[first, last]`, or (`last` the
/// smaller, the edge reversed) a trimmed copy reversed (`Geom_TrimmedCurve::Reverse`).
pub(super) fn edge_adaptor(curve: &Curve, first: f64, last: f64) -> GeomAdaptor {
    if last >= first {
        return GeomAdaptor::new(curve, first, last);
    }
    let mut basis = curve;
    while let Curve::Trimmed(t) = basis {
        basis = t.basis();
    }
    let mut reversed = basis.clone();
    let (u1, u2) = (basis.reversed_parameter(first), basis.reversed_parameter(last));
    match &mut reversed {
        Curve::BSpline(b) => b.reverse_as_occt(),
        // (`Geom_OffsetCurve::Reverse`: its basis reversed, the offset negated.)
        Curve::Offset(o) if matches!(o.basis(), Curve::BSpline(_)) => {
            let mut inner = o.basis().clone();
            if let Curve::BSpline(b) = &mut inner {
                b.reverse_as_occt();
            }
            *o = crate::geometry::curve::OffsetCurve::new(inner, -o.offset(), o.direction());
        }
        other => other.reverse(),
    }
    GeomAdaptor::new(&reversed, u1, u2)
}

/// `corr2PI_PI`: an angle in `[0, 2π)` brought into `[-π, π)`.
fn corr(a: f64) -> f64 {
    if a < std::f64::consts::PI { a } else { a - 2.0 * std::f64::consts::PI }
}

/// `diffAng`.
fn diff_angle(a: f64, o: f64) -> f64 {
    let tau = 2.0 * std::f64::consts::PI;
    let d = (a - o) - ((a - o) / 2.0 / std::f64::consts::PI).floor() * tau;
    if d >= 0.0 { corr(d) } else { -corr(-d) }
}

/// `FindPlane`: whether the curve lies in a plane (a conic's, else the plane of inertia of
/// its poles, or of samples of it, all within `Precision::Confusion` of it).
fn find_plane(curve: &GeomAdaptor) -> bool {
    let points: Vec<Point> = match (curve.kind(), curve.curve()) {
        (CurveKind::Line, _) => return false,
        (CurveKind::Circle | CurveKind::Ellipse | CurveKind::Hyperbola | CurveKind::Parabola, _) => return true,
        (CurveKind::Bezier, Curve::Bezier(b)) => {
            if b.poles().len() <= 2 {
                return false;
            }
            b.poles().to_vec()
        }
        (CurveKind::BSpline, Curve::BSpline(b)) => {
            if b.poles().len() <= 2 {
                return false;
            }
            b.poles().to_vec()
        }
        _ => {
            let nbp = 15 + curve.intervals(Continuity::C3).len() - 1;
            let (f, l) = (curve.first_parameter(), curve.last_parameter());
            let inv = 1.0 / (nbp - 1) as f64;
            (1..=nbp).map(|ii| curve.value((f * (nbp - ii) as f64 + l * (ii - 1) as f64) * inv)).collect()
        }
    };
    let Some((origin, normal)) = axe_of_inertia(&points, 1.0e-7) else { return false };
    let d = -normal.dot(Vector::new(origin.x, origin.y, origin.z));
    points.iter().all(|p| (normal.x * p.x + normal.y * p.y + normal.z * p.z + d).abs() <= CONFUSION)
}

/// `GeomLib::AxeOfInertia` (with `GeomLib::Inertia`): the barycentre and the normal of the
/// plane of inertia of `points` (the largest axis crossed with the middle one); `None` where
/// it is singular (the middle spread, times the count, within `tolerance`). The
/// eigen-decomposition is kentron's Jacobi rotations (OCCT's `math_Jacobi`).
fn axe_of_inertia(points: &[Point], tolerance: f64) -> Option<(Point, Vector)> {
    let nb = points.len() as f64;
    let g = points.iter().fold(Vector::new(0.0, 0.0, 0.0), |s, p| s + Vector::new(p.x, p.y, p.z)) * (1.0 / nb);
    let mut m = [[0.0; 3]; 3];
    for p in points {
        let d = [p.x - g.x, p.y - g.y, p.z - g.z];
        for i in 0..3 {
            for j in 0..3 {
                m[i][j] += d[i] * d[j];
            }
        }
    }
    for row in &mut m {
        for x in row.iter_mut() {
            *x /= nb;
        }
    }
    let (values, vectors) = crate::properties::symmetric_eigen(&crate::math::Matrix3::from_rows(m));
    // (Ascending: the largest is OCCT's X direction, the middle its Y.)
    let y_gap = values[1].abs().sqrt();
    if y_gap * nb <= tolerance {
        return None;
    }
    let z = vectors[2].cross(vectors[1]);
    Some((Point::new(g.x, g.y, g.z), z * (1.0 / z.length())))
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
