//! What a B-spline strip that a straight segment sweeps along a screw is made from: the
//! segment and the screw. Kept with the strip (with its clones, until it changes), so that a
//! point's foot on it, and its sections with other strips of the same screw, are found from
//! the exact surface (the helicoid of the segment) rather than from the B-spline's samples.

use crate::math::{Axis, Point, Transform, Transformable, Vector};

/// A curve that is an arc of a screw's helix (one point's path as the screw turns it): the
/// screw's axis and how far it rises along it for each radian turned (right-handed about
/// the axis's direction). The B-spline it is kept with lies within about 1e-10 of its size of
/// the helix.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct Helix {
    pub(crate) axis: Axis,
    pub(crate) rise_per_radian: f64,
}

impl Helix {
    /// The helix moved rigidly by `t`; `None` where `t` mirrors.
    pub(crate) fn transformed(&self, t: &Transform) -> Option<Helix> {
        if t.is_negative() || (t.scale_factor() - 1.0).abs() > 1e-12 {
            return None;
        }
        Some(Helix { axis: self.axis.transformed(t), rise_per_radian: self.rise_per_radian })
    }

    /// Whether `other` is a helix of the same screw: about the same line, rising as much a
    /// radian (about either direction along it: a screw turned about the line the other way
    /// rises the other way along it too), within `tolerance` of a size of `scale`.
    pub(crate) fn same_screw(&self, axis: Axis, rise_per_radian: f64, scale: f64) -> bool {
        let z = Vector::from(self.axis.direction);
        let d = axis.origin - self.axis.origin;
        Vector::from(axis.direction).cross(z).length() <= 1e-12
            && (d - z * d.dot(z)).length() <= 1e-9 * (1.0 + scale)
            && (rise_per_radian - self.rise_per_radian).abs() <= 1e-12 * self.rise_per_radian.abs()
    }
}

/// The surface `(u, v) -> M(v)·S(u)`: `S` the segment from `start` (at `u0`) to `end` (at
/// `u1`), `M(v)` the screw about `axis`, right-handed, turning `turn·(v - v0)` and rising
/// `rise·(v - v0)` along it. The B-spline it is kept with lies within about 1e-10 of its
/// size of it, at the same parameters.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct ScrewStrip {
    pub(crate) axis: Axis,
    pub(crate) turn: f64,
    pub(crate) rise: f64,
    pub(crate) v0: f64,
    pub(crate) start: Point,
    pub(crate) end: Point,
    pub(crate) u0: f64,
    pub(crate) u1: f64,
}

impl ScrewStrip {
    /// The screw at `v`.
    pub(crate) fn place(&self, v: f64) -> Transform {
        let t = v - self.v0;
        Transform::translation(Vector::from(self.axis.direction) * (self.rise * t)) * Transform::rotation(self.axis, self.turn * t)
    }

    /// The segment's point at `u`.
    pub(crate) fn segment_point(&self, u: f64) -> Point {
        self.start + (self.end - self.start) * ((u - self.u0) / (self.u1 - self.u0))
    }

    pub(crate) fn point(&self, u: f64, v: f64) -> Point {
        self.segment_point(u).transformed(&self.place(v))
    }

    /// How far the screw rises along its axis for a turn of one radian (the same for every
    /// strip of one screw, whatever its parameters).
    pub(crate) fn rise_per_radian(&self) -> f64 {
        self.rise / self.turn
    }

    /// The strip moved rigidly by `t`; `None` where `t` mirrors (the screw would turn the
    /// other way).
    pub(crate) fn transformed(&self, t: &Transform) -> Option<ScrewStrip> {
        if t.is_negative() || (t.scale_factor() - 1.0).abs() > 1e-12 {
            return None;
        }
        Some(ScrewStrip { axis: self.axis.transformed(t), start: self.start.transformed(t), end: self.end.transformed(t), ..self.clone() })
    }

    /// The parameters within `[u_lo, u_hi] × [v_lo, v_hi]` of the point of the exact surface
    /// nearest `p`, and its distance. (Screwed back by `v`, `p` runs along its helix; the
    /// segment's point nearest it is had in closed form, which leaves a search along `v`:
    /// samples close enough to tell the helix's passes by the segment apart, then, about the
    /// nearest, the root of the distance's derivative along `v` (regula falsi, kept in its
    /// bracket by halving the far end's value), else that sample.)
    pub(crate) fn nearest(&self, p: Point, (u_lo, u_hi): (f64, f64), (v_lo, v_hi): (f64, f64)) -> (f64, f64, f64) {
        let o = self.axis.origin;
        let z = Vector::from(self.axis.direction);
        let d = self.end - self.start;
        let length2 = d.dot(d);
        let (s_lo, s_hi) = ((u_lo - self.u0) / (self.u1 - self.u0), (u_hi - self.u0) / (self.u1 - self.u0));
        let (s_lo, s_hi) = (s_lo.min(s_hi), s_lo.max(s_hi));
        // p screwed back by v: rotated by -turn·(v - v0) about the axis, lowered by rise·(v - v0).
        let rel = p - o;
        let (h, radial) = (rel.dot(z), rel - z * rel.dot(z));
        let across = z.cross(radial);
        // (The point screwed back from its angle's cosine and sine; the squared distance to
        // the segment, the segment's parameter, and the squared distance's derivative in v.)
        let measure = |t: f64, cos: f64, sin: f64| -> (f64, f64, f64) {
            let q = o + radial * cos + across * sin + z * (h - self.rise * t);
            let s = if length2 > 0.0 { ((q - self.start).dot(d) / length2).clamp(s_lo, s_hi) } else { s_lo };
            let r = q - (self.start + d * s);
            let dq = z.cross(q - o) * (-self.turn) - z * self.rise;
            (r.dot(r), s, 2.0 * r.dot(dq))
        };
        let at = |v: f64| {
            let t = v - self.v0;
            let (sin, cos) = (-self.turn * t).sin_cos();
            measure(t, cos, sin)
        };
        // (Each sample a 64th of a turn apart at most: the helix passes the segment's line at
        // most twice a turn. Turned from one to the next by the step's rotation.)
        let span = (self.turn * (v_hi - v_lo)).abs();
        let n = 16usize.max((span / (std::f64::consts::TAU / 64.0)).ceil() as usize).min(4096);
        let step = (v_hi - v_lo) / n as f64;
        let (sin0, cos0) = (-self.turn * (v_lo - self.v0)).sin_cos();
        let (sin_step, cos_step) = (-self.turn * step).sin_cos();
        let (mut cos, mut sin) = (cos0, sin0);
        let mut samples = Vec::with_capacity(n + 1);
        for i in 0..=n {
            let v = if i == n { v_hi } else { v_lo + step * i as f64 };
            // (Exactly at the ends and every 16 steps, against the rounding gathered.)
            if i % 16 == 0 || i == n {
                let (s, c) = (-self.turn * (v - self.v0)).sin_cos();
                (cos, sin) = (c, s);
            }
            let (f, s, df) = measure(v - self.v0, cos, sin);
            samples.push((v, f, s, df));
            (cos, sin) = (cos * cos_step - sin * sin_step, sin * cos_step + cos * sin_step);
        }
        let best = (0..=n).min_by(|&i, &j| samples[i].1.total_cmp(&samples[j].1)).expect("samples");
        let (mut v, mut f, mut s) = (samples[best].0, samples[best].1, samples[best].2);
        // The bracket about the nearest sample where the derivative goes from - to +.
        let bracket = if best > 0 && samples[best - 1].3 < 0.0 && samples[best].3 > 0.0 {
            Some((best - 1, best))
        } else if best < n && samples[best].3 < 0.0 && samples[best + 1].3 > 0.0 {
            Some((best, best + 1))
        } else {
            None
        };
        if let Some((i, j)) = bracket {
            let (mut a, mut fa) = (samples[i].0, samples[i].3);
            let (mut b, mut fb) = (samples[j].0, samples[j].3);
            let mut side = 0;
            for _ in 0..60 {
                let m = (a * fb - b * fa) / (fb - fa);
                let m = if m > a && m < b { m } else { 0.5 * (a + b) };
                let (fm, sm, dm) = at(m);
                if fm < f {
                    (v, f, s) = (m, fm, sm);
                }
                if dm == 0.0 || (b - a).abs() <= 1e-15 * (1.0 + a.abs() + b.abs()) {
                    break;
                }
                if dm < 0.0 {
                    (a, fa) = (m, dm);
                    if side == -1 {
                        fb *= 0.5;
                    }
                    side = -1;
                } else {
                    (b, fb) = (m, dm);
                    if side == 1 {
                        fa *= 0.5;
                    }
                    side = 1;
                }
            }
        }
        (self.u0 + s * (self.u1 - self.u0), v, f.max(0.0).sqrt())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::Direction;

    /// The nearest point of a strip (a segment slanting across a screw's axis, swept half a
    /// turn) against the nearest of a fine grid of its points, for points around it.
    #[test]
    fn nearest_point_of_a_strip() {
        let axis = Axis::new(Point::new(0.3, -0.2, 0.1), Direction::new(0.1, 0.2, 1.0));
        let strip = ScrewStrip { axis, turn: std::f64::consts::PI, rise: 0.7, v0: 0.2, start: Point::new(2.0, 0.5, -0.3), end: Point::new(4.5, -1.0, 0.6), u0: 1.0, u1: 3.0 };
        let (us, vs) = ((1.0, 3.0), (0.2, 1.2));
        for i in 0..40 {
            let f = i as f64;
            let p = Point::new(3.0 * (0.7 * f).cos(), 3.0 * (1.3 * f).sin(), 0.8 * (0.4 * f).sin());
            let (u, v, d) = strip.nearest(p, us, vs);
            assert!((strip.point(u, v).distance(p) - d).abs() <= 1e-12, "distance at the parameters found");
            let n = 400;
            let grid = (0..=n).flat_map(|a| (0..=n).map(move |b| (a, b))).map(|(a, b)| strip.point(1.0 + 2.0 * a as f64 / n as f64, 0.2 + b as f64 / n as f64).distance(p)).fold(f64::INFINITY, f64::min);
            assert!(d <= grid + 1e-12, "{d} farther than the grid's {grid}");
            assert!(d >= grid - 0.02, "{d} far nearer than the grid's {grid}");
        }
    }
}
