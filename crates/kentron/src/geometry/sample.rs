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

/// The most samples [`quasi_uniform_deflection`] takes before it stops refining (the
/// deflection it reports may then exceed the one asked for).
const MOST_SAMPLES: usize = 8192;

/// How many equal steps of a span the chord's distance from the curve is measured at.
const STEPS: usize = 16;

/// Samples of a bounded `curve` whose polyline stays within `deflection` of it: both ends,
/// plus for a straight curve (a line, a two-pole B-spline or Bézier curve) nothing else,
/// for a circle equal steps each within the deflection by its sagitta (and at most a half
/// turn), and for other curves the halving of spans until every chord is within the
/// deflection of the curve, measured at [`STEPS`] equal steps of its span. A B-spline's
/// inner knots where it may have a corner (multiplicity at least the degree) are samples
/// too. The `deflection` reported is the largest distance so measured.
///
/// `None` for a deflection not above zero, or an unbounded or empty parameter range.
#[doc(alias = "GCPnts_QuasiUniformDeflection")]
pub fn quasi_uniform_deflection(curve: &Curve, deflection: f64) -> Option<CurveSampling> {
    let (a, b) = (curve.first_parameter(), curve.last_parameter());
    let bounded = |t: f64| t.is_finite() && t.abs() < 1e50;
    if !(deflection > 0.0) || !bounded(a) || !bounded(b) || !(b > a) {
        return None;
    }
    let basis = match curve {
        Curve::Trimmed(t) => t.basis(),
        c => c,
    };
    let mut parameters = match basis {
        Curve::Line(_) => vec![a, b],
        Curve::BSpline(s) if s.poles().len() == 2 => vec![a, b],
        Curve::Bezier(s) if s.poles().len() == 2 => vec![a, b],
        Curve::Circle(c) => {
            let r = c.radius();
            let step = if deflection >= r { std::f64::consts::PI } else { (2.0 * (1.0 - deflection / r).acos()).min(std::f64::consts::PI) };
            let mut n = ((b - a) / step).ceil().max(1.0) as usize;
            // (Rounding at the boundary: one more step rather than a sagitta just over.)
            if r * (1.0 - (0.5 * (b - a) / n as f64).cos()) > deflection {
                n += 1;
            }
            (0..=n).map(|i| if i == n { b } else { a + (b - a) * i as f64 / n as f64 }).collect()
        }
        _ => {
            let mut breaks = vec![a];
            if let Curve::BSpline(s) = basis {
                let last = s.knots().len() - 1;
                for i in 1..last {
                    let k = s.knots()[i];
                    if s.multiplicities()[i] >= s.degree() && k > a && k < b {
                        breaks.push(k);
                    }
                }
            }
            breaks.push(b);
            refine(curve, &breaks, deflection)
        }
    };
    parameters.dedup();
    let points: Vec<Point> = parameters.iter().map(|&t| curve.point(t)).collect();
    let measured = (1..parameters.len()).map(|i| chord_deviation(curve, parameters[i - 1], parameters[i], points[i - 1], points[i])).fold(0.0, f64::max);
    Some(CurveSampling { parameters, points, deflection: measured })
}

/// Halves the spans between consecutive `breaks` until each chord is within `deflection`
/// of the curve (or the sample budget or rounding stops it): the sample parameters.
fn refine(curve: &Curve, breaks: &[f64], deflection: f64) -> Vec<f64> {
    let mut out = vec![breaks[0]];
    for w in breaks.windows(2) {
        // Spans still to check, the next one on top, so the output stays in order.
        let mut stack = vec![(w[0], w[1], curve.point(w[0]), curve.point(w[1]))];
        while let Some((s, e, ps, pe)) = stack.pop() {
            let mid = 0.5 * (s + e);
            let unsplittable = !(mid > s && mid < e);
            if unsplittable || out.len() + stack.len() >= MOST_SAMPLES || chord_deviation(curve, s, e, ps, pe) <= deflection {
                out.push(e);
                continue;
            }
            let pm = curve.point(mid);
            stack.push((mid, e, pm, pe));
            stack.push((s, mid, ps, pm));
        }
    }
    out
}

/// The largest distance from the curve between `s` and `e` (measured at [`STEPS`] equal
/// steps) to the chord from `ps` to `pe`.
fn chord_deviation(curve: &Curve, s: f64, e: f64, ps: Point, pe: Point) -> f64 {
    let chord = pe - ps;
    let length2 = chord.dot(chord);
    (1..STEPS)
        .map(|i| {
            let q = curve.point(s + (e - s) * i as f64 / STEPS as f64);
            let w = q - ps;
            let along = if length2 > 0.0 { (w.dot(chord) / length2).clamp(0.0, 1.0) } else { 0.0 };
            let off: Vector = w - chord * along;
            off.length()
        })
        .fold(0.0, f64::max)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::curve::{BSplineCurve, Circle, Ellipse, Line, TrimmedCurve};
    use crate::math::{Direction, Frame};

    #[test]
    fn straight_curves_take_their_ends() {
        let line = Curve::Trimmed(TrimmedCurve::new(Curve::Line(Line::new(Point::ORIGIN, Direction::X)), -1.0, 3.0));
        assert_eq!(quasi_uniform_deflection(&line, 1e-6).unwrap().parameters, vec![-1.0, 3.0]);
        let spline = Curve::BSpline(BSplineCurve::new(vec![Point::ORIGIN, Point::new(1.0, 2.0, 3.0)], vec![0.0, 2.0], vec![2, 2], 1));
        assert_eq!(quasi_uniform_deflection(&spline, 1e-6).unwrap().parameters, vec![0.0, 2.0]);
        assert!(quasi_uniform_deflection(&spline, 0.0).is_none());
        assert!(quasi_uniform_deflection(&Curve::Line(Line::new(Point::ORIGIN, Direction::X)), 1.0).is_none());
    }

    #[test]
    fn circles_take_equal_steps_by_their_sagitta() {
        let r = 3.0;
        let circle = Curve::Circle(Circle::new(Frame::WORLD, r));
        for d in [1e-1, 1e-3, 1e-6, 10.0] {
            let s = quasi_uniform_deflection(&circle, d).unwrap();
            let fewest = (std::f64::consts::TAU / (2.0 * (1.0 - (d / r).min(1.0)).acos()).min(std::f64::consts::PI)).ceil() as usize + 1;
            assert!(s.parameters.len() == fewest || s.parameters.len() == fewest + 1, "{d}: {} vs {fewest}", s.parameters.len());
            assert!(s.deflection <= d);
            let step = s.parameters[1] - s.parameters[0];
            assert!(s.parameters.windows(2).all(|w| ((w[1] - w[0]) - step).abs() < 1e-12));
        }
    }

    #[test]
    fn ellipses_keep_within_the_deflection_and_need_more_points_for_less() {
        let ellipse = Curve::Ellipse(Ellipse::new(Frame::WORLD, 10.0, 1.0));
        let mut last = 0;
        for k in 1..=6 {
            let d = 10f64.powi(-k);
            let s = quasi_uniform_deflection(&ellipse, d).unwrap();
            assert!(s.deflection <= d);
            assert!(s.parameters.len() >= last);
            assert!(s.parameters.windows(2).all(|w| w[0] < w[1]));
            for (t, p) in s.parameters.iter().zip(&s.points) {
                assert_eq!(*p, ellipse.point(*t));
            }
            last = s.parameters.len();
        }
    }

    #[test]
    fn corners_of_b_splines_are_samples() {
        let poles = vec![Point::ORIGIN, Point::new(1.0, 1.0, 0.0), Point::new(2.0, 0.0, 0.0), Point::new(3.0, 1.0, 0.0), Point::new(4.0, 0.0, 0.0)];
        let corner = Curve::BSpline(BSplineCurve::new(poles, vec![0.0, 0.37, 1.0], vec![3, 2, 3], 2));
        let s = quasi_uniform_deflection(&corner, 0.5).unwrap();
        assert!(s.parameters.contains(&0.37));
        // All poles in one place: a point.
        let dot = Curve::BSpline(BSplineCurve::new(vec![Point::ORIGIN; 4], vec![0.0, 1.0], vec![4, 4], 3));
        assert_eq!(quasi_uniform_deflection(&dot, 1e-9).unwrap().parameters, vec![0.0, 1.0]);
    }
}
