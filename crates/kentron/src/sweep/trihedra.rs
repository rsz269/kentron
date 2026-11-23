//! Two more pipe frame laws, per spine piece: a frame keeping a fixed binormal, and a
//! frame carried by least rotations between tangents at fixed samples (a discrete
//! rotation-method frame, after Klok 1986 and Bloomenthal 1990).

use crate::error::{Error, Result};
use crate::geometry::adaptor::{CurveAdaptor, CurveKind, GeomAdaptor};
use crate::geometry::curve::Curve;
use crate::geometry::Continuity;
use crate::math::Vector;

use super::frenet::{default_normal, least_rotation, travelled, unit, Frenet};

/// Parameters closer than this to a sample are the sample.
const AT_SAMPLE: f64 = 1e-9;

/// The frame whose binormal is a fixed direction: the spine's tangent with its component
/// along the direction taken out.
#[derive(Clone, Debug)]
pub(super) struct ConstantBinormalLaw {
    frenet: Frenet,
    binormal: Vector,
}

impl ConstantBinormalLaw {
    /// The law along `curve` travelled from `first` to `last`, keeping `binormal` (made unit).
    pub(super) fn new(curve: &Curve, first: f64, last: f64, binormal: Vector) -> ConstantBinormalLaw {
        let (view, _, _) = travelled(curve, first, last);
        ConstantBinormalLaw { frenet: Frenet::new(view), binormal: unit(binormal) }
    }

    /// (tangent, normal) at fraction `fraction` of the way: `N = D × T_F` and `T = N × D`
    /// normalised, or, where the Frenet tangent `T_F` runs along `D`, `T` along `N_F × D`.
    pub(super) fn at_fraction(&self, fraction: f64) -> (Vector, Vector) {
        let d = self.binormal;
        let frame = self.frenet.frame(self.frenet.at_fraction(fraction));
        let across = d.cross(frame.t);
        if across.length() > 1e-7 {
            let n = unit(across);
            (n.cross(d), n)
        } else {
            let t = unit(frame.n.cross(d));
            (t, d.cross(t))
        }
    }
}

/// A frame turned from sample to sample by the least rotation between their tangents, and
/// from the sample below to any parameter between them.
#[derive(Clone, Debug)]
pub(super) struct DiscreteLaw {
    frenet: Frenet,
    /// Lines and conics: the Frenet frame.
    analytic: bool,
    /// The curve's smooth spans (`CN` intervals).
    spans: Vec<f64>,
    /// Sample parameters, tangents and normals.
    samples: Vec<(f64, Vector, Vector)>,
}

impl DiscreteLaw {
    /// The law along `curve` travelled from `first` to `last`. Fails where a sample has no
    /// tangent and no chord to give one.
    pub(super) fn new(curve: &Curve, first: f64, last: f64) -> Result<DiscreteLaw> {
        let (view, _, b) = travelled(curve, first, last);
        let analytic = matches!(view.kind(), CurveKind::Line | CurveKind::Circle | CurveKind::Ellipse | CurveKind::Hyperbola | CurveKind::Parabola);
        let frenet = Frenet::without_singular_points(view.clone());
        let spans = view.intervals(Continuity::CN);
        if analytic {
            return Ok(DiscreteLaw { frenet, analytic, spans, samples: Vec::new() });
        }
        let mut us: Vec<f64> = Vec::with_capacity(10 * spans.len());
        for w in spans.windows(2) {
            us.extend((0..10).map(|k| w[0] + (w[1] - w[0]) * k as f64 / 10.0));
        }
        us.push(b);
        let tangent = |k: usize, us: &[f64]| -> Result<Vector> {
            let u = us[k];
            let (_, d1) = view.d1(u);
            if d1.length() >= 1e-7 {
                return Ok(unit(d1));
            }
            let Some(&next) = us.get(k + 1) else { return Err(Error::Failed("a discrete frame's last sample has no tangent".into())) };
            chord_tangent(&view, u, next)
        };
        let mut samples: Vec<(f64, Vector, Vector)> = Vec::with_capacity(us.len());
        let t0 = tangent(0, &us)?;
        samples.push((us[0], t0, default_normal(t0).0));
        let mut k = 1;
        while k < us.len() {
            let t = tangent(k, &us)?;
            let &(u_prev, t_prev, x_prev) = samples.last().expect("a sample");
            let axis = t_prev.cross(t);
            if axis == Vector::new(0.0, 0.0, 0.0) && t_prev.dot(t) <= 0.0 {
                // (Opposite tangents: a sample midway first.)
                let mid = 0.5 * (u_prev + us[k]);
                if !(us[k] - u_prev >= f64::MIN_POSITIVE) || mid <= u_prev || mid >= us[k] {
                    return Err(Error::Failed("a discrete frame between opposite tangents".into()));
                }
                us.insert(k, mid);
                continue;
            }
            samples.push((us[k], t, carry(x_prev, t_prev, t)));
            k += 1;
        }
        Ok(DiscreteLaw { frenet, analytic, spans, samples })
    }

    /// (tangent, normal) at fraction `fraction` of the way.
    pub(super) fn at_fraction(&self, fraction: f64) -> (Vector, Vector) {
        let u = self.frenet.at_fraction(fraction);
        if self.analytic {
            let frame = self.frenet.frame(u);
            return (frame.t, frame.n);
        }
        let samples = &self.samples;
        let last = samples.len() - 1;
        let mut i = if u <= samples[0].0 { 0 } else { samples.partition_point(|s| s.0 < u).clamp(1, last) - 1 };
        if i < last && (u - samples[i + 1].0).abs() < AT_SAMPLE {
            i += 1;
        }
        let (ui, ti, xi) = samples[i];
        if (u - ui).abs() < AT_SAMPLE {
            return (ti, xi);
        }
        // (The tangent in the smooth span holding u: the one before where it is a bound.)
        let j = self.spans[1..].partition_point(|&b| b < u).min(self.spans.len() - 2);
        let view = self.frenet.curve().trim(self.spans[j], self.spans[j + 1]);
        let (_, d1) = view.d1(u);
        let t = if d1.length() >= 1e-7 {
            unit(d1)
        } else {
            let next = samples.get(i + 1).map_or(self.frenet.last(), |s| s.0);
            match chord_tangent(&view, u, next) {
                Ok(t) => t,
                Err(_) => ti,
            }
        };
        (t, carry(xi, ti, t))
    }
}

/// `x` (square to the unit `from`) turned by the least rotation onto the unit `to`, made
/// square to it again.
fn carry(x: Vector, from: Vector, to: Vector) -> Vector {
    let turned = least_rotation(x, from, to);
    unit(to.cross(turned.cross(to)))
}

/// The direction of the chord from `u` a tenth of the way towards `next` (all the way where
/// a tenth is under `1e-9`).
fn chord_tangent(view: &GeomAdaptor, u: f64, next: f64) -> Result<Vector> {
    let mut e = (next - u) / 10.0;
    if e.abs() < 1e-9 {
        e = next - u;
    }
    let chord = view.value(u + e) - view.value(u);
    if !(chord.length() > 0.0) {
        return Err(Error::Failed("a discrete frame's sample has no tangent".into()));
    }
    Ok(unit(chord))
}

#[cfg(test)]
#[path = "trihedra_tests.rs"]
mod tests;
