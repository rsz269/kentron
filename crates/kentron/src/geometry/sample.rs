//! Sampling points along curves.

use crate::math::{Point, Vector};

use super::curve::Curve;
use super::ParametricCurve;

/// Points sampled along a curve.
#[derive(Debug, Clone, PartialEq)]
pub struct CurveSampling {
    /// Increasing parameters of the samples, including both ends.
    pub parameters: Vec<f64>,
    pub points: Vec<Point>,
    /// The largest distance between the curve and the polyline through the points.
    pub deflection: f64,
}

/// Samples a bounded curve so that the polyline through the points stays within
/// `deflection` of it, refining each span by halving until it is flat enough. `None` if the
/// curve is unbounded or `deflection` is not positive.
#[doc(alias = "GCPnts_QuasiUniformDeflection")]
pub fn quasi_uniform_deflection(curve: &Curve, deflection: f64) -> Option<CurveSampling> {
    let (u1, u2) = (curve.first_parameter(), curve.last_parameter());
    let finite = |t: f64| t.is_finite() && t.abs() < 1e50;
    if !(deflection > 0.0) || !finite(u1) || !finite(u2) {
        return None;
    }
    // As OCCT's GCPnts_QuasiUniformDeflection (for C1 curves).
    let basis = match curve {
        Curve::Trimmed(t) => t.basis(),
        c => c,
    };
    let mut sampler = Sampler::new(curve, deflection);
    match basis {
        Curve::Line(_) => sampler.straight(u1, u2),
        Curve::BSpline(b) if b.poles().len() == 2 => sampler.straight(u1, u2),
        Curve::Bezier(b) if b.poles().len() == 2 => sampler.straight(u1, u2),
        Curve::Circle(c) => {
            let step = 2.0 * (1.0 - deflection / c.radius()).max(0.0).acos();
            let n = ((u2 - u1) / step) as usize + 2;
            let step = (u2 - u1) / (n - 1) as f64;
            for i in 0..n {
                let u = u1 + step * i as f64;
                sampler.parameters.push(u);
                sampler.points.push(curve.point(u));
            }
        }
        Curve::BSpline(b) => {
            // A curve with kinks is sampled piece by piece between its knots of continuity
            // below C2.
            let p = b.degree();
            let inner = |m: usize| b.knots().iter().zip(b.multiplicities()).filter(move |&(&k, &mult)| mult >= m && u1 < k && k < u2);
            if inner(p).next().is_some() {
                let mut cuts: Vec<f64> = std::iter::once(u1).chain(inner(p.saturating_sub(1).max(1)).map(|(&k, _)| k)).collect();
                cuts.push(u2);
                for w in cuts.windows(2) {
                    if sampler.points.len() > 0 {
                        sampler.parameters.pop();
                        sampler.points.pop();
                    }
                    sampler.curved(w[0], w[1]);
                }
            } else {
                sampler.curved(u1, u2);
            }
        }
        _ => sampler.curved(u1, u2),
    }
    // The deflection actually reached.
    let mut worst = 0.0f64;
    for i in 0..sampler.parameters.len() - 1 {
        let (a, b) = (sampler.points[i], sampler.points[i + 1]);
        let (s, t) = (sampler.parameters[i], sampler.parameters[i + 1]);
        for k in 1..16 {
            let q = curve.point(s + (t - s) * k as f64 / 16.0);
            worst = worst.max(distance_to_segment(q, a, b));
        }
    }
    Some(CurveSampling { parameters: sampler.parameters, points: sampler.points, deflection: worst })
}

/// The gap between `x` and the next larger float.
fn ulp(x: f64) -> f64 {
    let a = x.abs();
    if a == 0.0 { f64::MIN_POSITIVE } else { f64::from_bits(a.to_bits() + 1) - a }
}

fn distance_to_segment(q: Point, a: Point, b: Point) -> f64 {
    let ab = b - a;
    let len2 = ab.dot(ab);
    let t = if len2 > 0.0 { ((q - a).dot(ab) / len2).clamp(0.0, 1.0) } else { 0.0 };
    q.distance(a + t * ab)
}

/// The most refinement steps one sampling may take.
const MAX_CALLS: usize = 2000;

struct Sampler<'a> {
    curve: &'a Curve,
    deflection2: f64,
    /// Squared lengths below this count as zero.
    eps2: f64,
    calls: usize,
    parameters: Vec<f64>,
    points: Vec<Point>,
}

impl<'a> Sampler<'a> {
    fn new(curve: &'a Curve, deflection: f64) -> Self {
        // The parameter step that moves the curve by the confusion tolerance, estimated
        // from its fastest speed.
        let (u1, u2) = (curve.first_parameter(), curve.last_parameter());
        let speed = (0..=32).map(|i| curve.d1(u1 + (u2 - u1) * i as f64 / 32.0).1.length()).fold(0.0, f64::max);
        let eps = if speed > 0.0 { crate::tolerance::LINEAR / speed } else { crate::tolerance::LINEAR };
        Self { curve, deflection2: deflection * deflection, eps2: eps * eps, calls: 0, parameters: Vec::new(), points: Vec::new() }
    }

    fn straight(&mut self, u1: f64, u2: f64) {
        for u in [u1, u2] {
            self.parameters.push(u);
            self.points.push(self.curve.point(u));
        }
    }

    fn curved(&mut self, u1: f64, u2: f64) {
        let (p1, d1) = self.curve.d1(u1);
        self.parameters.push(u1);
        self.points.push(p1);
        let (p2, d2) = self.curve.d1(u2 - 10.0 * ulp(u2));
        self.calls = 0;
        self.refine((u1, p1, d1), (u2, p2, d2), 2);
    }

    /// Appends samples after `start` up to `end`, the span being cut into at least
    /// `pieces - 1` parts.
    fn refine(&mut self, start: (f64, Point, Vector), end: (f64, Point, Vector), pieces: usize) {
        self.calls += 1;
        let before = self.points.len();
        if self.calls >= MAX_CALLS || (self.calls > 100 && before < 2) {
            return;
        }
        let (u_start, p_start, v_start) = start;
        let mut step = end.0 - u_start;
        let (p_step, v_step) = if pieces > 2 {
            step /= (pieces - 1) as f64;
            self.curve.d1(u_start + step)
        } else {
            (end.1, end.2)
        };
        // The sag from how much the tangent turns over the chord, unless the chord is short.
        let chord2 = (p_step - p_start).length_squared();
        let mut turn_sag = None;
        if chord2 > self.eps2 && chord2 > 16.0 * self.deflection2 {
            let (n1, n2) = (v_start.length_squared(), v_step.length_squared());
            if n1 > self.eps2 && n2 > self.eps2 {
                let turn = (v_start.normalized() - v_step.normalized()).length_squared();
                if turn > self.eps2 {
                    turn_sag = Some(turn * chord2 / 64.0);
                }
            }
        }
        // The midpoint's distance from the chord's, which non-uniform speed can inflate.
        let mid = p_start + 0.5 * (p_step - p_start);
        let mid_sag = mid.distance_squared(self.curve.point(u_start + 0.5 * step));
        let sag = match turn_sag {
            Some(s) if mid_sag <= chord2 / 4.0 + self.deflection2 => s,
            _ => mid_sag,
        };
        let u_step = u_start + step;
        if sag < self.deflection2 {
            self.parameters.push(u_step);
            self.points.push(p_step);
        } else {
            self.refine(start, (u_step, p_step, v_step), 3);
        }
        if pieces > 2 {
            let added = self.points.len() - before;
            self.refine((u_step, p_step, v_step), end, (pieces as isize - added as isize).max(0) as usize);
        }
        self.calls -= 1;
    }
}
