//! Converting curves and surfaces to B-splines, and joining curves into one B-spline.

use crate::math::{Point, Point2};

use super::bspline::BSplineCurve;
use super::control::ControlPoint;
use super::parametric::{Continuity, ParametricCurve, ParametricSurface};
use super::curve::Curve;
use super::curve2::{BSplineCurve2, Curve2};
use super::surface::{BSplineSurface, Surface};

/// The B-spline (possibly rational) representation of a curve. Bounded curves keep their
/// parametrization where the conversion allows it; a full circle or ellipse becomes a
/// periodic curve. Returns `None` for an unbounded curve (an untrimmed line, parabola or
/// hyperbola).
#[doc(alias = "GeomConvert::CurveToBSplineCurve")]
pub fn curve_to_bspline(curve: &Curve) -> Option<BSplineCurve> {
    let form = match curve {
        Curve::Line(_) => Form::Line,
        Curve::Circle(c) => Form::Ellipse(c.center()),
        Curve::Ellipse(c) => Form::Ellipse(c.center()),
        Curve::Hyperbola(c) => Form::Hyperbola(c.center()),
        Curve::Parabola(_) => Form::Parabola,
        Curve::Bezier(b) => return Some(b.to_bspline()),
        Curve::BSpline(b) => return Some(b.clone()),
        Curve::Trimmed(t) => {
            let (a, b) = (t.first_parameter(), t.last_parameter());
            return match t.basis() {
                Curve::Bezier(_) | Curve::BSpline(_) | Curve::Trimmed(_) => {
                    let mut whole = curve_to_bspline(t.basis())?;
                    whole.segment(a, b);
                    Some(whole)
                }
                basis => conic(basis, form_of(basis)?, a, b),
            };
        }
        Curve::Offset(_) => return None,
    };
    conic(curve, form, curve.first_parameter(), curve.last_parameter())
}

/// The kind of a curve [`conic`] converts, with its center where it has one.
#[derive(Clone, Copy)]
enum Form<P> {
    Line,
    Ellipse(P),
    Hyperbola(P),
    Parabola,
}

fn form_of(curve: &Curve) -> Option<Form<Point>> {
    Some(match curve {
        Curve::Line(_) => Form::Line,
        Curve::Circle(c) => Form::Ellipse(c.center()),
        Curve::Ellipse(c) => Form::Ellipse(c.center()),
        Curve::Hyperbola(c) => Form::Hyperbola(c.center()),
        Curve::Parabola(_) => Form::Parabola,
        _ => return None,
    })
}

fn form_of2(curve: &Curve2) -> Option<Form<Point2>> {
    Some(match curve {
        Curve2::Line(_) => Form::Line,
        Curve2::Circle(c) => Form::Ellipse(c.center()),
        Curve2::Ellipse(c) => Form::Ellipse(c.center()),
        Curve2::Hyperbola(c) => Form::Hyperbola(c.center()),
        Curve2::Parabola(_) => Form::Parabola,
        _ => return None,
    })
}

/// A line, conic arc or parabola over `[a, b]` as an exact B-spline: a segment, rational
/// quadratic arcs of at most 120° (a whole ellipse a periodic curve), one rational
/// quadratic piece of hyperbola, a polynomial quadratic piece of parabola. `None` if the
/// range is unbounded.
fn conic<P: ControlPoint, C: ParametricCurve<Point = P, Vector = P::Vector>>(curve: &C, form: Form<P>, a: f64, b: f64) -> Option<BSplineCurve<P>> {
    if !(a.is_finite() && b.is_finite() && a < b) {
        return None;
    }
    let point = |t: f64| curve.point(t).coords();
    let combine = |x: [f64; 3], s: f64, y: [f64; 3]| P::from_coords([x[0] + s * y[0], x[1] + s * y[1], x[2] + s * y[2]]);
    let diff = |x: [f64; 3], y: [f64; 3]| [x[0] - y[0], x[1] - y[1], x[2] - y[2]];
    match form {
        Form::Line => BSplineCurve::try_new(vec![curve.point(a), curve.point(b)], vec![a, b], vec![2, 2], 1).ok(),
        Form::Parabola => {
            let tangent = P::vector_coords(curve.d1(a).1);
            let middle = combine(point(a), 0.5 * (b - a), tangent);
            BSplineCurve::try_new(vec![curve.point(a), middle, curve.point(b)], vec![a, b], vec![3, 3], 2).ok()
        }
        Form::Hyperbola(center) => {
            let half = 0.5 * (b - a);
            let w = half.cosh();
            let middle = combine(center.coords(), 1.0 / w, diff(point(0.5 * (a + b)), center.coords()));
            let curve = BSplineCurve::try_new(vec![curve.point(a), middle, curve.point(b)], vec![a, b], vec![3, 3], 2).ok()?;
            Some(curve.with_weights(vec![1.0, w, 1.0]))
        }
        Form::Ellipse(center) => {
            let full = (b - a - std::f64::consts::TAU).abs() <= crate::tolerance::PARAMETRIC;
            let spans = (((b - a) / (std::f64::consts::TAU / 3.0)) - 1e-9).ceil().max(1.0) as usize;
            let step = (b - a) / spans as f64;
            let w = (0.5 * step).cos();
            let knots: Vec<f64> = (0..=spans).map(|i| a + i as f64 * step).collect();
            let mut poles = Vec::new();
            let mut weights = Vec::new();
            for i in 0..spans {
                let middle = combine(center.coords(), 1.0 / w, diff(point(knots[i] + 0.5 * step), center.coords()));
                poles.extend([curve.point(knots[i]), middle]);
                weights.extend([1.0, w]);
            }
            let curve = if full {
                BSplineCurve::try_new_periodic(poles, knots, vec![2; spans + 1], 2).ok()?
            } else {
                poles.push(curve.point(b));
                weights.push(1.0);
                let mut multiplicities = vec![2; spans + 1];
                multiplicities[0] = 3;
                multiplicities[spans] = 3;
                BSplineCurve::try_new(poles, knots, multiplicities, 2).ok()?
            };
            Some(curve.with_weights(weights))
        }
    }
}

/// [`curve_to_bspline`] for 2D curves.
#[doc(alias = "Geom2dConvert::CurveToBSplineCurve")]
pub fn curve2_to_bspline(curve: &Curve2) -> Option<BSplineCurve2> {
    match curve {
        Curve2::Bezier(b) => Some(b.to_bspline()),
        Curve2::BSpline(b) => Some(b.clone()),
        Curve2::Offset(_) => None,
        Curve2::Trimmed(t) => {
            let (a, b) = (t.first_parameter(), t.last_parameter());
            match t.basis() {
                Curve2::Bezier(_) | Curve2::BSpline(_) | Curve2::Trimmed(_) => {
                    let mut whole = curve2_to_bspline(t.basis())?;
                    whole.segment(a, b);
                    Some(whole)
                }
                basis => conic(basis, form_of2(basis)?, a, b),
            }
        }
        _ => conic(curve, form_of2(curve)?, curve.first_parameter(), curve.last_parameter()),
    }
}

/// The B-spline representation of a bounded surface (e.g. a trimmed plane); `None` if the
/// surface is unbounded.
#[doc(alias = "GeomConvert::SurfaceToBSplineSurface")]
pub fn surface_to_bspline(surface: &Surface) -> Option<BSplineSurface> {
    match surface {
        Surface::BSpline(b) => Some(b.clone()),
        Surface::Bezier(b) => {
            // A Bezier surface is a B-spline with one span.
            let (p, q) = (b.u_degree(), b.v_degree());
            let u = super::surface::KnotVector::new(vec![0.0, 1.0], vec![p + 1, p + 1], p);
            let v = super::surface::KnotVector::new(vec![0.0, 1.0], vec![q + 1, q + 1], q);
            let surface = BSplineSurface::try_new(b.poles().to_vec(), u, v).ok()?;
            Some(match b.weights() {
                Some(w) => surface.with_weights(w.to_vec()),
                None => surface,
            })
        }
        Surface::Trimmed(t) => {
            let (u1, u2, v1, v2) = t.bounds();
            surface_patch_to_bspline(t.basis(), (u1, u2), (v1, v2))
        }
        other => {
            let (u1, u2, v1, v2) = other.bounds();
            surface_patch_to_bspline(other, (u1, u2), (v1, v2))
        }
    }
}

/// The exact B-spline form of `surface` over a bounded parameter rectangle (`None` if it
/// is unbounded or of a kind without one), keeping its parameters at the rectangle's
/// corners: a tensor product of a circle arc (rational quadratic) in the angular direction
/// with the profile in the other (a line, a circle arc, the basis curve).
fn surface_patch_to_bspline(surface: &Surface, (u1, u2): (f64, f64), (v1, v2): (f64, f64)) -> Option<BSplineSurface> {
    use crate::geometry::curve::{Circle, Line, TrimmedCurve};
    use crate::math::{Frame, Vector};
    use super::surface::KnotVector;
    if ![u1, u2, v1, v2].iter().all(|x| x.is_finite() && x.abs() < 1e50) {
        return None;
    }
    // A unit circle arc over [a, b] (poles (cos-like, sin-like) with weights), and a
    // segment over [a, b] (poles a and b).
    let arc = |a: f64, b: f64| curve_to_bspline(&Curve::Trimmed(TrimmedCurve::new(Curve::Circle(Circle::new(Frame::WORLD, 1.0)), a, b))).map(|c| c.clamped());
    let segment = |a: f64, b: f64| {
        let line = Curve::Line(Line::new(Point::new(0.0, 0.0, 0.0), crate::math::Direction::X));
        curve_to_bspline(&Curve::Trimmed(TrimmedCurve::new(line, a, b)))
    };
    let knots = |c: &BSplineCurve| KnotVector::new(c.knots().to_vec(), c.multiplicities().to_vec(), c.degree());
    // The tensor product of `u` and `v` with poles `f(u pole, v pole)`.
    let product = |u: &BSplineCurve, v: &BSplineCurve, f: &dyn Fn(Point, Point) -> Point| -> Option<BSplineSurface> {
        let poles: Vec<Vec<Point>> = u.poles().iter().map(|&a| v.poles().iter().map(|&b| f(a, b)).collect()).collect();
        let weights: Vec<Vec<f64>> = (0..u.poles().len()).map(|i| (0..v.poles().len()).map(|j| u.weight(i) * v.weight(j)).collect()).collect();
        let surface = BSplineSurface::try_new(poles, knots(u), knots(v)).ok()?;
        Some(if u.is_rational() || v.is_rational() { surface.with_weights(weights) } else { surface })
    };
    let frame_vectors = |f: &Frame| (f.origin(), Vector::from(f.x()), Vector::from(f.y()), Vector::from(f.z()));
    match surface {
        Surface::Plane(_) => {
            let (lu, lv) = (segment(u1, u2)?, segment(v1, v2)?);
            product(&lu, &lv, &|a, b| surface.point(a.x, b.x))
        }
        Surface::Cylinder(c) => {
            let (o, x, y, z) = frame_vectors(c.frame());
            let r = c.radius();
            product(&arc(u1, u2)?, &segment(v1, v2)?, &|a, b| o + (x * a.x + y * a.y) * r + z * b.x)
        }
        Surface::Cone(c) => {
            let (o, x, y, z) = frame_vectors(c.frame());
            let (sin, cos) = c.semi_angle_sin_cos();
            let r = c.radius();
            product(&arc(u1, u2)?, &segment(v1, v2)?, &|a, b| o + (x * a.x + y * a.y) * (r + b.x * sin) + z * (b.x * cos))
        }
        Surface::Sphere(sp) => {
            let (o, x, y, z) = frame_vectors(sp.frame());
            let r = sp.radius();
            product(&arc(u1, u2)?, &arc(v1, v2)?, &|a, b| o + (x * a.x + y * a.y) * (r * b.x) + z * (r * b.y))
        }
        Surface::Torus(t) => {
            let (o, x, y, z) = frame_vectors(t.frame());
            let (big, small) = (t.major_radius(), t.minor_radius());
            product(&arc(u1, u2)?, &arc(v1, v2)?, &|a, b| o + (x * a.x + y * a.y) * (big + small * b.x) + z * (small * b.y))
        }
        Surface::Extrusion(e) => {
            let basis = curve_to_bspline(&Curve::Trimmed(TrimmedCurve::new(e.basis().clone(), u1, u2)))?.clamped();
            let d = Vector::from(e.direction());
            product(&basis, &segment(v1, v2)?, &|a, b| a + d * b.x)
        }
        Surface::Revolution(r) => {
            let basis = curve_to_bspline(&Curve::Trimmed(TrimmedCurve::new(r.basis().clone(), v1, v2)))?.clamped();
            let axis = r.axis();
            let (o, z) = (axis.origin, Vector::from(axis.direction));
            // A point turned by the angle whose (cos, sin) is (a.x, a.y) about the axis.
            product(&arc(u1, u2)?, &basis, &|a, b| {
                let d = b - o;
                let along = z * d.dot(z);
                let radial = d - along;
                o + along + radial * a.x + z.cross(radial) * a.y
            })
        }
        Surface::BSpline(b) => {
            let mut b = b.clone();
            b.segment(u1, u2, v1, v2);
            Some(b)
        }
        Surface::Bezier(_) => surface_to_bspline(surface),
        Surface::Trimmed(t) => surface_patch_to_bspline(t.basis(), (u1, u2), (v1, v2)),
        Surface::Offset(_) => None,
    }
}

/// Settings for [`approximate_curve`] and [`approximate_curve2`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CurveApproximationOptions {
    /// The largest allowed distance between the curve and its approximation.
    pub tolerance: f64,
    /// The continuity required at the approximation's inner knots.
    pub continuity: Continuity,
    /// The largest number of polynomial spans.
    pub max_segments: usize,
    pub max_degree: usize,
}

impl Default for CurveApproximationOptions {
    fn default() -> Self {
        Self { tolerance: 1.0e-4, continuity: Continuity::C2, max_segments: 100, max_degree: 14 }
    }
}

/// A B-spline approximating another curve, and how far apart they are.
#[derive(Debug, Clone, PartialEq)]
pub struct ApproximatedCurve<P: ControlPoint = Point> {
    /// Parametrized over the same range as the approximated curve.
    pub curve: BSplineCurve<P>,
    pub max_error: f64,
}

/// A polynomial B-spline approximating a bounded curve (e.g. an offset curve, which has
/// no exact B-spline form). `None` if the curve is unbounded or the tolerance cannot be
/// met within the options' limits.
#[doc(alias = "GeomConvert_ApproxCurve")]
pub fn approximate_curve(curve: &Curve, options: &CurveApproximationOptions) -> Option<ApproximatedCurve> {
    approximate_parametric(curve, options)
}

/// [`approximate_curve`] for 2D curves.
#[doc(alias = "Geom2dConvert_ApproxCurve")]
pub fn approximate_curve2(
    curve: &Curve2,
    options: &CurveApproximationOptions,
) -> Option<ApproximatedCurve<Point2>> {
    approximate_parametric(curve, options)
}

/// Least-squares B-splines on more and more uniform spans (kept apart at the curve's
/// breaks) until one is within the tolerance.
fn approximate_parametric<P, C>(curve: &C, options: &CurveApproximationOptions) -> Option<ApproximatedCurve<P>>
where
    P: ControlPoint,
    C: ParametricCurve<Point = P, Vector = P::Vector>,
{
    let (a, b) = (curve.first_parameter(), curve.last_parameter());
    if !(a.is_finite() && b.is_finite() && a.abs() < 1e50 && b.abs() < 1e50 && a < b) {
        return None;
    }
    let order = match options.continuity {
        Continuity::C0 | Continuity::G1 => 0,
        Continuity::C1 | Continuity::G2 => 1,
        Continuity::C2 => 2,
        _ => 3,
    };
    let degree = options.max_degree.min(5).max(order + 1).max(1);
    if degree > options.max_degree.max(1) {
        return None;
    }
    // The curve's own breaks stay knots of full multiplicity (no continuity claimed there).
    let mut breaks: Vec<f64> = curve.breaks().into_iter().filter(|&k| a < k && k < b && !curve.is_cn(order as u32)).collect();
    breaks.sort_by(f64::total_cmp);
    breaks.dedup();
    let pieces = breaks.len() + 1;
    let mut best: Option<ApproximatedCurve<P>> = None;
    let mut per_piece = 1;
    while per_piece * pieces <= options.max_segments.max(1) {
        let cuts: Vec<f64> = std::iter::once(a).chain(breaks.iter().copied()).chain(std::iter::once(b)).collect();
        let mut flat = vec![a; degree + 1];
        for (w, piece) in cuts.windows(2).zip(0..) {
            if piece > 0 {
                flat.extend(std::iter::repeat(w[0]).take(degree));
            }
            for i in 1..per_piece {
                flat.push(w[0] + (w[1] - w[0]) * i as f64 / per_piece as f64);
            }
        }
        flat.extend(std::iter::repeat(b).take(degree + 1));
        // Samples: several per span, including the span ends.
        let spans: Vec<f64> = { let mut k = flat.clone(); k.dedup(); k };
        let mut params = Vec::new();
        for w in spans.windows(2) {
            for i in 0..4 * (degree + 1) {
                params.push(w[0] + (w[1] - w[0]) * i as f64 / (4 * (degree + 1)) as f64);
            }
        }
        params.push(b);
        let q: Vec<[f64; 3]> = params.iter().map(|&t| curve.point(t).coords()).collect();
        if let Some(fitted) = super::fit::least_squares_with_knots::<P>(&q, &params, degree, flat) {
            // The error, checked between the samples.
            let mut max_error = 0.0f64;
            for w in spans.windows(2) {
                for i in 0..=8 * (degree + 1) {
                    let t = w[0] + (w[1] - w[0]) * i as f64 / (8 * (degree + 1)) as f64;
                    max_error = max_error.max(distance(fitted.point(t).coords(), curve.point(t).coords()));
                }
            }
            let better = best.as_ref().map_or(true, |x| max_error < x.max_error);
            if better {
                best = Some(ApproximatedCurve { curve: fitted, max_error });
            }
            if max_error <= options.tolerance {
                break;
            }
        }
        per_piece *= 2;
    }
    best.filter(|x| x.max_error <= options.tolerance)
}

fn distance(a: [f64; 3], b: [f64; 3]) -> f64 {
    ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)).sqrt()
}

/// Joins curves end to end into a single B-spline.
///
/// Each added curve must touch the current result at one of its ends (within the given
/// tolerance); it is reversed and/or prepended as needed.
#[doc(alias = "GeomConvert_CompCurveToBSplineCurve")]
#[doc(alias = "Geom2dConvert_CompCurveToBSplineCurve")]
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Concatenation<P: ControlPoint = Point> {
    curve: Option<BSplineCurve<P>>,
}

impl<P: ControlPoint> Concatenation<P> {
    /// An empty concatenation.
    pub fn new() -> Self {
        Self { curve: None }
    }

    /// The joined curve so far; `None` if nothing has been added.
    pub fn result(&self) -> Option<&BSplineCurve<P>> {
        self.curve.as_ref()
    }

    pub fn into_result(self) -> Option<BSplineCurve<P>> {
        self.curve
    }

    pub fn clear(&mut self) {
        self.curve = None;
    }

    /// Starts with `curve` (converted to a B-spline).
    pub fn from_curve<C: ToBSpline<Point = P>>(curve: &C) -> Self {
        let mut chain = Self::new();
        chain.add(curve, 0.0);
        chain
    }

    /// Converts `curve` to a B-spline and joins it. False (leaving the result unchanged)
    /// if it cannot be converted or does not touch the result.
    pub fn add<C: ToBSpline<Point = P>>(&mut self, curve: &C, tolerance: f64) -> bool {
        match curve.to_bspline() {
            Some(bspline) => self.add_bspline(&bspline, tolerance),
            None => false,
        }
    }

    /// Joins `curve`; false (leaving the result unchanged) if it does not touch it.
    pub fn add_bspline(&mut self, curve: &BSplineCurve<P>, tolerance: f64) -> bool {
        let curve = curve.clamped();
        let Some(current) = &self.curve else {
            self.curve = Some(curve);
            return true;
        };
        let near = |a: P, b: P| {
            let (a, b) = (a.coords(), b.coords());
            (0..3).map(|k| (a[k] - b[k]).powi(2)).sum::<f64>().sqrt() <= tolerance.max(crate::tolerance::LINEAR)
        };
        let reversed = || {
            let mut r = curve.clone();
            r.reverse();
            r
        };
        // Which ends meet: after the result (as it is or reversed), else before it.
        let (append, next) = if near(current.end_point(), curve.start_point()) {
            (true, curve.clone())
        } else if near(current.end_point(), curve.end_point()) {
            (true, reversed())
        } else if near(current.start_point(), curve.end_point()) {
            (false, curve.clone())
        } else if near(current.start_point(), curve.start_point()) {
            (false, reversed())
        } else {
            return false;
        };
        let (mut a, mut b) = if append { (current.clone(), next) } else { (next, current.clone()) };
        let degree = a.degree().max(b.degree());
        a.increase_degree(degree);
        b.increase_degree(degree);
        self.curve = Some(a.joined(&b));
        true
    }
}

/// Curves with a B-spline representation (see [`curve_to_bspline`]).
pub trait ToBSpline {
    type Point: ControlPoint;
    fn to_bspline(&self) -> Option<BSplineCurve<Self::Point>>;
}

impl ToBSpline for Curve {
    type Point = Point;
    fn to_bspline(&self) -> Option<BSplineCurve> {
        curve_to_bspline(self)
    }
}

impl ToBSpline for Curve2 {
    type Point = Point2;
    fn to_bspline(&self) -> Option<BSplineCurve2> {
        curve2_to_bspline(self)
    }
}
