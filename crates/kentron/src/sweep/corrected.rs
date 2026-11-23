//! The default pipe frame: a rotation-minimizing frame (also called a Bishop or
//! parallel-transport frame) along each spine piece.
//!
//! The normal is carried along the spine by parallel transport in the normal plane: it never
//! turns about the tangent. It starts as the piece's Frenet normal at its start (the default
//! normal of the tangent where the piece starts straight).
//!
//! The transport is integrated by the double-reflection method (Wang, Jüttler, Zheng and Liu,
//! "Computation of rotation minimizing frames", ACM TOG 2008), whose step is accurate to the
//! fourth order: each step reflects the frame in the plane bisecting the chord, then in the
//! plane taking the reflected tangent onto the new one. Steps are chosen adaptively, each
//! compared with two half steps; the stored frames take the extrapolated (Richardson) value.
//! Between stored frames the law makes one double-reflection step from the frame below, then
//! turns by its share of the small difference to the stored frame above, so that the frames
//! are continuous and smooth in the parameter.
//!
//! At a tangent break inside a piece the normal is carried across by the least rotation
//! between the two tangents. Along a closed piece (ending where it starts, with the same
//! tangent) the transported normal comes back turned by the holonomy angle; a twist growing
//! evenly with arc length takes it out, so the end frame is the start frame.

use crate::error::{Error, Result};
use crate::geometry::abscissa::signed_length;
use crate::geometry::adaptor::{CurveAdaptor, CurveKind, GeomAdaptor};
use crate::geometry::curve::Curve;
use crate::geometry::Continuity;
use crate::math::{Point, Vector};

use super::frenet::{least_rotation, rotate_about, signed_angle, travelled, unit, unit_tangent, Frenet};
pub(super) use super::frenet::FrenetLaw;

/// The largest angle (radians) a step may differ from its two half steps by.
const STEP_TOLERANCE: f64 = 1e-10;

/// A transported frame at a parameter.
#[derive(Clone, Copy, Debug)]
struct Sample {
    u: f64,
    point: Point,
    /// The derivatives of the curve there (for chords to nearby parameters).
    d1: Vector,
    d2: Vector,
    d3: Vector,
    tangent: Vector,
    normal: Vector,
    /// The turn about the tangent taking a single step from this sample to the next one
    /// onto the stored next frame.
    correction: f64,
    /// The arc length from the piece's start (closed pieces only).
    length: f64,
}

/// The transported frames along one C2 stretch of the piece.
#[derive(Clone, Debug)]
struct Stretch {
    view: GeomAdaptor,
    samples: Vec<Sample>,
}

/// The even twist taking out a closed piece's holonomy.
#[derive(Clone, Copy, Debug)]
struct Closure {
    /// The turn at the end (minus the holonomy).
    angle: f64,
    length: f64,
}

/// The rotation-minimizing frame along one spine piece.
#[derive(Clone, Debug)]
pub(super) struct CorrectedFrenet {
    first: f64,
    last: f64,
    /// A line: the start frame everywhere.
    constant: Option<(Vector, Vector)>,
    stretches: Vec<Stretch>,
    closure: Option<Closure>,
}

/// The normal `n0` (at `p0` with tangent `t0`) carried to `p1` with tangent `t1` by two
/// reflections: in the plane bisecting the chord `chord = p1 − p0`, then in the plane
/// taking the reflected tangent onto `t1`.
fn double_reflection(chord: Vector, t0: Vector, n0: Vector, t1: Vector) -> Vector {
    let c1 = chord.dot(chord);
    let (n_l, t_l) = if c1 > 0.0 {
        (n0 - chord * (2.0 * chord.dot(n0) / c1), t0 - chord * (2.0 * chord.dot(t0) / c1))
    } else {
        (n0, -t0)
    };
    let v2 = t1 - t_l;
    let c2 = v2.dot(v2);
    let n = if c2 > 0.0 { n_l - v2 * (2.0 * v2.dot(n_l) / c2) } else { n_l };
    // (Square to the tangent and unit, to rounding.)
    unit(n - t1 * t1.dot(n))
}

impl Stretch {
    /// The chord from sample `s` to the curve point at `u` (`p` there): by Taylor's
    /// expansion close to the sample, where the difference of points would lose digits.
    fn chord(s: &Sample, u: f64, p: Point, step: f64) -> Vector {
        let h = u - s.u;
        if h.abs() < 1e-3 * step.abs() {
            s.d1 * h + s.d2 * (0.5 * h * h) + s.d3 * (h * h * h / 6.0)
        } else {
            p - s.point
        }
    }

    fn sample(view: &GeomAdaptor, u: f64, normal: Vector) -> Sample {
        let (point, d1, d2, d3) = view.d3(u);
        let tangent = unit_tangent(view, u);
        Sample { u, point, d1, d2, d3, tangent, normal, correction: 0.0, length: 0.0 }
    }

    /// One step from sample `s` to `u`.
    fn step(view: &GeomAdaptor, s: &Sample, u: f64, step: f64) -> (Point, Vector, Vector) {
        let p = view.value(u);
        let t = unit_tangent(view, u);
        let n = double_reflection(Self::chord(s, u, p, step), s.tangent, s.normal, t);
        (p, t, n)
    }

    /// The frames along `view` from the normal `normal` at its start.
    fn transport(view: GeomAdaptor, normal: Vector) -> Result<Stretch> {
        let (a, b) = (view.first_parameter(), view.last_parameter());
        let width = b - a;
        let start = Self::sample(&view, a, normal);
        if !start.tangent.length().is_finite() {
            return Err(Error::Failed("a pipe's spine with no tangent".into()));
        }
        let mut samples = vec![start];
        let mut h = width / 8.0;
        let smallest = 1e-12 * width.abs().max(1e-300);
        while samples.last().expect("a sample").u < b {
            let s = *samples.last().expect("a sample");
            let u1 = if s.u + h >= b - 1e-3 * h { b } else { s.u + h };
            let step = u1 - s.u;
            let (_, t1, full) = Self::step(&view, &s, u1, step);
            let um = s.u + 0.5 * step;
            let (_, _, mid_normal) = Self::step(&view, &s, um, step);
            let mid = Sample { normal: mid_normal, ..Self::sample(&view, um, mid_normal) };
            let (_, _, two) = Self::step(&view, &mid, u1, 0.5 * step);
            let difference = signed_angle(full, two, t1);
            if !difference.is_finite() || !t1.length().is_finite() {
                return Err(Error::Failed("a pipe's spine with no tangent".into()));
            }
            let error = difference.abs();
            if error > STEP_TOLERANCE && step > smallest {
                h = step * (0.9 * (STEP_TOLERANCE / error).powf(0.2)).max(0.2);
                continue;
            }
            // (Extrapolated: the error of a step is 16/15 of its difference from two.)
            let correction = difference * 16.0 / 15.0;
            let normal = rotate_about(full, t1, correction);
            samples.last_mut().expect("a sample").correction = correction;
            samples.push(Self::sample(&view, u1, unit(normal - t1 * t1.dot(normal))));
            let grow = if error > 0.0 { 0.9 * (STEP_TOLERANCE / error).powf(0.2) } else { 2.0 };
            h = step * grow.clamp(0.2, 2.0);
        }
        Ok(Stretch { view, samples })
    }

    /// The tangent and normal at `u` inside the stretch, and the arc length from the
    /// piece's start (when `lengths`).
    fn at(&self, u: f64, lengths: bool) -> (Vector, Vector, f64) {
        let samples = &self.samples;
        let k = samples.partition_point(|s| s.u <= u).clamp(1, samples.len() - 1) - 1;
        let (s, next) = (&samples[k], &samples[k + 1]);
        let length = |s: &Sample| if lengths { s.length + signed_length(&self.view, s.u, u, 1e-13 * (1.0 + s.length)) } else { 0.0 };
        if u == s.u {
            return (s.tangent, s.normal, if lengths { s.length } else { 0.0 });
        }
        if u == next.u {
            return (next.tangent, next.normal, if lengths { next.length } else { 0.0 });
        }
        let step = next.u - s.u;
        let (_, t, n) = Self::step(&self.view, s, u, step);
        let n = rotate_about(n, t, s.correction * (u - s.u) / step);
        (t, n, length(s))
    }
}

impl CorrectedFrenet {
    /// The law along `curve` travelled from `first` to `last` (`last < first`: backwards).
    pub(super) fn new(curve: &Curve, first: f64, last: f64) -> Result<CorrectedFrenet> {
        if !(first != last) || !first.is_finite() || !last.is_finite() {
            return Err(Error::InvalidInput(format!("a pipe's spine piece over the empty range [{first}, {last}]")));
        }
        let (view, a, b) = travelled(curve, first, last);
        let frenet = Frenet::without_singular_points(view.clone());
        let start = frenet.frame_in(a, 0);
        if !(start.t.length().is_finite() && start.n.length().is_finite()) {
            return Err(Error::Failed("a pipe's spine with no tangent at its start".into()));
        }
        if view.kind() == CurveKind::Line {
            return Ok(CorrectedFrenet { first: a, last: b, constant: Some((start.t, start.n)), stretches: Vec::new(), closure: None });
        }
        let breaks = view.intervals(Continuity::C2);
        let mut stretches: Vec<Stretch> = Vec::with_capacity(breaks.len() - 1);
        let mut normal = start.n;
        for w in breaks.windows(2) {
            let piece = if breaks.len() == 2 { view.clone() } else { view.trim(w[0], w[1]) };
            if let Some(before) = stretches.last() {
                // (Across a break: the least rotation between the tangents either side.)
                let end = before.samples.last().expect("a sample");
                let t = unit_tangent(&piece, w[0]);
                let n = least_rotation(end.normal, end.tangent, t);
                normal = unit(n - t * t.dot(n));
            }
            let stretch = Stretch::transport(piece, normal)?;
            normal = stretch.samples.last().expect("a sample").normal;
            stretches.push(stretch);
        }
        let mut law = CorrectedFrenet { first: a, last: b, constant: None, stretches, closure: None };
        law.close(&view, start.n);
        Ok(law)
    }

    /// The twist taking out the holonomy, where the piece is closed.
    fn close(&mut self, view: &GeomAdaptor, start_normal: Vector) {
        let first = &self.stretches[0].samples[0];
        let end = *self.stretches.last().expect("a stretch").samples.last().expect("a sample");
        let size = crate::geometry::abscissa::length(view, None);
        let tangent_gap = first.tangent.cross(end.tangent).length().atan2(first.tangent.dot(end.tangent));
        if !(first.point.distance(end.point) <= 1e-7 * size && tangent_gap <= 1e-9) {
            return;
        }
        let holonomy = signed_angle(start_normal, end.normal, end.tangent);
        // (Arc lengths at the samples.)
        let mut total = 0.0;
        for stretch in &mut self.stretches {
            let mut previous: Option<f64> = None;
            for i in 0..stretch.samples.len() {
                if let Some(u) = previous {
                    let s = stretch.samples[i].u;
                    total += signed_length(&stretch.view, u, s, 1e-13 * (1.0 + total));
                }
                stretch.samples[i].length = total;
                previous = Some(stretch.samples[i].u);
            }
        }
        if total > 0.0 {
            self.closure = Some(Closure { angle: -holonomy, length: total });
        }
    }

    /// (tangent, normal) at fraction `fraction` ∈ [0, 1] of the way.
    pub(super) fn at_fraction(&self, fraction: f64) -> (Vector, Vector) {
        if let Some(frame) = self.constant {
            return frame;
        }
        let u = if fraction == 1.0 { self.last } else { self.first + (self.last - self.first) * fraction };
        // (The stretch holding u: the one before where it is a bound.)
        let i = self.stretches.partition_point(|s| s.view.last_parameter() < u).min(self.stretches.len() - 1);
        let (t, n, length) = self.stretches[i].at(u, self.closure.is_some());
        match self.closure {
            Some(c) => (t, rotate_about(n, t, c.angle * (length / c.length))),
            None => (t, n),
        }
    }
}

#[cfg(test)]
#[path = "corrected_tests.rs"]
pub(super) mod tests;
