//! Arc length along a curve view, and its inverse: the parameter reached by travelling a given
//! arc length from a parameter.
//!
//! Lengths integrate the speed `|C′(t)|` by adaptive Gauss–Kronrod quadrature (the 7-point
//! Gauss rule embedded in the 15-point Kronrod rule; the difference of the two estimates the
//! error), starting from the curve's smooth spans so that no rule straddles a knot, and
//! splitting the interval of largest estimated error until the total is below the requested
//! tolerance. Curves of constant speed (lines, circles, straight two-pole polynomial
//! curves) are measured by the formula.
//!
//! The inverse walks the spans from the start parameter in the direction of travel,
//! accumulating their lengths, and solves for the parameter inside the span that holds the
//! target by Newton's method on the arc length (whose derivative is the speed), safeguarded
//! by bisection of the span, integrating only the increments between iterates.

use crate::error::{Error, Result};
use crate::geometry::adaptor::{CurveAdaptor, CurveKind};
use crate::geometry::curve::Curve;
use crate::geometry::Continuity;

/// Nodes of the 15-point Kronrod rule on `[-1, 1]` (non-negative half, largest first); the
/// odd-indexed ones (and the centre) are the 7-point Gauss rule's.
const KRONROD_NODES: [f64; 8] = [
    0.991_455_371_120_812_6,
    0.949_107_912_342_758_5,
    0.864_864_423_359_769_1,
    0.741_531_185_599_394_4,
    0.586_087_235_467_691_1,
    0.405_845_151_377_397_2,
    0.207_784_955_007_898_5,
    0.0,
];

const KRONROD_WEIGHTS: [f64; 8] = [
    0.022_935_322_010_529_225,
    0.063_092_092_629_978_55,
    0.104_790_010_322_250_18,
    0.140_653_259_715_525_92,
    0.169_004_726_639_267_9,
    0.190_350_578_064_785_4,
    0.204_432_940_075_298_9,
    0.209_482_141_084_727_83,
];

/// Weights of the 7-point Gauss rule at `KRONROD_NODES[1]`, `[3]`, `[5]` and the centre.
const GAUSS_WEIGHTS: [f64; 4] = [0.129_484_966_168_869_7, 0.279_705_391_489_276_7, 0.381_830_050_505_118_94, 0.417_959_183_673_469_4];

/// The most intervals one adaptive integration splits into.
const MAX_PIECES: usize = 4000;

/// A whole curve shorter than this cannot be travelled by [`AbscissaPoint::new`].
const SHORTEST: f64 = 1e-7;

/// The 15-point Kronrod estimate of `∫_a^b f` and its difference from the 7-point Gauss one.
fn kronrod(f: &dyn Fn(f64) -> f64, a: f64, b: f64) -> (f64, f64) {
    let (c, h) = (0.5 * (a + b), 0.5 * (b - a));
    let fc = f(c);
    let mut k = KRONROD_WEIGHTS[7] * fc;
    let mut g = GAUSS_WEIGHTS[3] * fc;
    for i in 0..7 {
        let x = h * KRONROD_NODES[i];
        let pair = f(c - x) + f(c + x);
        k += KRONROD_WEIGHTS[i] * pair;
        if i % 2 == 1 {
            g += GAUSS_WEIGHTS[i / 2] * pair;
        }
    }
    (k * h, ((k - g) * h).abs())
}

/// `∫ f` over consecutive intervals `breaks[i]..breaks[i + 1]` (a signed integral when they
/// decrease) to an absolute accuracy `tolerance` (or to rounding, or after `MAX_PIECES`
/// intervals, whichever comes first).
fn integrate(f: &dyn Fn(f64) -> f64, breaks: &[f64], tolerance: f64) -> f64 {
    let mut pieces: Vec<(f64, f64, f64, f64)> = Vec::new();
    for w in breaks.windows(2) {
        if w[0] != w[1] {
            let (v, e) = kronrod(f, w[0], w[1]);
            pieces.push((w[0], w[1], v, e));
        }
    }
    loop {
        let total: f64 = pieces.iter().map(|p| p.2).sum();
        let error: f64 = pieces.iter().map(|p| p.3).sum();
        let floor = 64.0 * f64::EPSILON * pieces.iter().map(|p| p.2.abs()).sum::<f64>();
        if error <= tolerance.max(floor) || pieces.len() >= MAX_PIECES || !total.is_finite() {
            return total;
        }
        let worst = (0..pieces.len()).max_by(|&i, &j| pieces[i].3.total_cmp(&pieces[j].3)).unwrap();
        let (a, b, _, _) = pieces[worst];
        let m = 0.5 * (a + b);
        if m == a || m == b {
            // (Cannot split further: settle for what this piece gives.)
            pieces[worst].3 = 0.0;
            continue;
        }
        let (v1, e1) = kronrod(f, a, m);
        let (v2, e2) = kronrod(f, m, b);
        pieces[worst] = (a, m, v1, e1);
        pieces.push((m, b, v2, e2));
    }
}

/// The speed `|C′(t)|`.
fn speed(curve: &dyn CurveAdaptor, t: f64) -> f64 {
    curve.d1(t).1.length()
}

/// The constant speed of a line, a circle or a straight non-rational two-pole curve.
fn constant_speed(curve: &dyn CurveAdaptor) -> Option<f64> {
    match (curve.curve_type()?, curve.basis_curve()?) {
        (CurveKind::Line, Curve::Line(_)) => Some(speed(curve, curve.first_parameter())),
        (CurveKind::Circle, Curve::Circle(c)) => Some(c.radius()),
        (CurveKind::BSpline, Curve::BSpline(b)) if b.poles().len() == 2 && !b.is_rational() => Some(speed(curve, curve.first_parameter())),
        (CurveKind::Bezier, Curve::Bezier(b)) if b.poles().len() == 2 && !b.is_rational() => Some(speed(curve, curve.first_parameter())),
        _ => None,
    }
}

/// `u1`, the curve's smooth-span breakpoints strictly between, and `u2` (in the order from
/// `u1` to `u2`).
fn span_breaks(curve: &dyn CurveAdaptor, u1: f64, u2: f64) -> Vec<f64> {
    let (lo, hi) = (u1.min(u2), u1.max(u2));
    let mut out = vec![lo];
    out.extend(curve.intervals(Continuity::CN).into_iter().filter(|&k| k > lo && k < hi));
    out.push(hi);
    if u1 > u2 {
        out.reverse();
    }
    out
}

/// The signed arc length `∫_{u1}^{u2} |C′|` to the absolute accuracy `tolerance`.
pub(crate) fn signed_length(curve: &dyn CurveAdaptor, u1: f64, u2: f64, tolerance: f64) -> f64 {
    if u1 == u2 {
        return 0.0;
    }
    if let Some(s) = constant_speed(curve) {
        return (u2 - u1) * s;
    }
    let f = |t: f64| speed(curve, t);
    integrate(&f, &span_breaks(curve, u1, u2), tolerance)
}

/// The default absolute accuracy for a length near `estimate`: `1e-10` relative (and
/// `1e-13` absolute for short lengths).
fn default_tolerance(estimate: f64) -> f64 {
    (1e-10 * estimate.abs()).max(1e-13)
}

/// A first, rough length (one rule per span) to scale the default tolerance by.
fn rough_length(curve: &dyn CurveAdaptor, u1: f64, u2: f64) -> f64 {
    let f = |t: f64| speed(curve, t);
    span_breaks(curve, u1, u2).windows(2).map(|w| kronrod(&f, w[0], w[1]).0).sum::<f64>().abs()
}

/// The arc length of the whole curve view, to `tolerance` (a length) or, without one, to a
/// relative accuracy of about `1e-10`.
pub(crate) fn length(curve: &dyn CurveAdaptor, tolerance: Option<f64>) -> f64 {
    length_between(curve, curve.first_parameter(), curve.last_parameter(), tolerance)
}

/// The arc length between the parameters `u1` and `u2` (either order), to `tolerance` (a
/// length) or, without one, to a relative accuracy of about `1e-10`.
pub(crate) fn length_between(curve: &dyn CurveAdaptor, u1: f64, u2: f64, tolerance: Option<f64>) -> f64 {
    if u1 == u2 {
        return 0.0;
    }
    if let Some(s) = constant_speed(curve) {
        return (u2 - u1).abs() * s;
    }
    let tolerance = match tolerance {
        Some(t) => t,
        None => default_tolerance(rough_length(curve, u1, u2)),
    };
    signed_length(curve, u1, u2, tolerance).abs()
}

/// The parameter reached by travelling the arc length `distance ≥ 0` from `u0` in the
/// direction `forward` (increasing parameters) or back, the lengths of the spans integrated
/// to `tolerance` each and the final one solved to it. Beyond the curve's range the walk goes
/// on for another range's width (the view extrapolates the curve); `None` past that.
pub(crate) fn travel(curve: &dyn CurveAdaptor, u0: f64, distance: f64, forward: bool, tolerance: f64) -> Option<f64> {
    if distance == 0.0 {
        return Some(u0);
    }
    let (first, last) = (curve.first_parameter(), curve.last_parameter());
    if let Some(s) = constant_speed(curve) {
        if s > 0.0 {
            return Some(if forward { u0 + distance / s } else { u0 - distance / s });
        }
        return None;
    }
    let width = (last - first).abs().max(f64::MIN_POSITIVE);
    let mut stops: Vec<f64> = curve.intervals(Continuity::CN);
    if forward {
        stops.retain(|&k| k > u0);
        if last <= u0 {
            stops.clear();
        }
        stops.push(last.max(u0) + width);
    } else {
        stops.retain(|&k| k < u0);
        stops.reverse();
        if first >= u0 {
            stops.clear();
        }
        stops.push(first.min(u0) - width);
    }
    let piece_tolerance = 0.25 * tolerance;
    let mut current = u0;
    let mut travelled = 0.0;
    for next in stops {
        let piece = signed_length(curve, current, next, piece_tolerance).abs();
        if travelled + piece >= distance {
            return solve_within(curve, current, next, distance - travelled, piece, tolerance);
        }
        travelled += piece;
        current = next;
    }
    None
}

/// The parameter between `start` and `end` (a smooth piece of length `piece`) at the arc
/// length `remaining ≤ piece` from `start`: Newton's method on the arc length, kept inside
/// a shrinking bracket by bisection.
fn solve_within(curve: &dyn CurveAdaptor, start: f64, end: f64, remaining: f64, piece: f64, tolerance: f64) -> Option<f64> {
    let dir = if end > start { 1.0 } else { -1.0 };
    let quad = 0.25 * tolerance;
    let target_error = 0.5 * tolerance;
    // Bracket: `behind` short of the target, `ahead` beyond it (as positions along `dir`).
    let (mut behind, mut ahead) = (start, end);
    let mut u = if piece > 0.0 { start + (end - start) * (remaining / piece) } else { end };
    let mut along = signed_length(curve, start, u, quad) * dir;
    for _ in 0..100 {
        let h = along - remaining;
        if h.abs() <= target_error {
            return Some(u);
        }
        if h < 0.0 {
            behind = u;
        } else {
            ahead = u;
        }
        let s = speed(curve, u);
        let mut next = if s > 0.0 && s.is_finite() { u - dir * h / s } else { f64::NAN };
        let inside = |x: f64| (x - behind) * dir > 0.0 && (ahead - x) * dir > 0.0;
        if !next.is_finite() || !inside(next) {
            next = 0.5 * (behind + ahead);
        }
        if next == u || next == behind || next == ahead {
            return Some(u);
        }
        along += signed_length(curve, u, next, quad) * dir;
        u = next;
    }
    Some(u)
}

/// The point at a signed arc length from a parameter: forwards for a positive abscissa,
/// backwards for a negative one. Computed on construction; [`parameter`](Self::parameter) is
/// `None` when no such point was found.
#[derive(Clone, Copy)]
pub(crate) struct AbscissaPoint<'a> {
    curve: &'a dyn CurveAdaptor,
    parameter: Option<f64>,
}

impl<'a> AbscissaPoint<'a> {
    /// To the default accuracy. Fails for a curve shorter than `1e-7` altogether; an abscissa
    /// of at most `1e-7` gives `u0`.
    pub(crate) fn new(curve: &'a dyn CurveAdaptor, abscissa: f64, u0: f64) -> Result<Self> {
        let total = length(curve, None);
        if total < SHORTEST {
            return Err(Error::InvalidInput(format!("a curve of length {total} has no point at an abscissa")));
        }
        if abscissa.abs() <= SHORTEST {
            return Ok(AbscissaPoint { curve, parameter: Some(u0) });
        }
        Ok(Self::solved(curve, abscissa, u0, default_tolerance(abscissa)))
    }

    /// Lengths and the result to the absolute length `tolerance`.
    pub(crate) fn with_tolerance(tolerance: f64, curve: &'a dyn CurveAdaptor, abscissa: f64, u0: f64) -> Self {
        Self::solved(curve, abscissa, u0, tolerance)
    }

    /// To the default accuracy; `ui` (a guess) does not change the result.
    pub(crate) fn with_guess(curve: &'a dyn CurveAdaptor, abscissa: f64, u0: f64, ui: f64) -> Self {
        let _ = ui;
        Self::solved(curve, abscissa, u0, default_tolerance(abscissa))
    }

    /// Lengths and the result to the absolute length `tolerance`; `ui` (a guess) does not
    /// change the result.
    pub(crate) fn with_guess_and_tolerance(curve: &'a dyn CurveAdaptor, abscissa: f64, u0: f64, ui: f64, tolerance: f64) -> Self {
        let _ = ui;
        Self::solved(curve, abscissa, u0, tolerance)
    }

    fn solved(curve: &'a dyn CurveAdaptor, abscissa: f64, u0: f64, tolerance: f64) -> Self {
        let parameter = if abscissa == 0.0 || !abscissa.is_finite() {
            abscissa.is_finite().then_some(u0)
        } else {
            let tolerance = if tolerance > 0.0 { tolerance } else { default_tolerance(abscissa) };
            travel(curve, u0, abscissa.abs(), abscissa > 0.0, tolerance)
        };
        AbscissaPoint { curve, parameter }
    }

    pub(crate) fn is_done(&self) -> bool {
        self.parameter.is_some()
    }

    /// The parameter found, `None` when not done.
    pub(crate) fn parameter(&self) -> Option<f64> {
        let _ = self.curve;
        self.parameter
    }
}

#[cfg(test)]
mod tests;
