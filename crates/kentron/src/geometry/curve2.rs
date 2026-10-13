//! 2D curves, used for planar sketches and as curves in a surface's (u, v) parameter space.

use crate::math::{Axis2, Direction2, Frame, Frame2, Point2, Transform2, Transformable2, Vector2};

use super::bezier::BezierCurve;
use super::bspline::BSplineCurve;
use super::curve::Curve;
use super::curve::TrimmedCurve;

/// A 2D Bezier curve.
#[doc(alias = "Geom2d_BezierCurve")]
pub type BezierCurve2 = BezierCurve<Point2>;

/// A 2D B-spline curve.
#[doc(alias = "Geom2d_BSplineCurve")]
pub type BSplineCurve2 = BSplineCurve<Point2>;

/// A trimmed 2D curve.
#[doc(alias = "Geom2d_TrimmedCurve")]
pub type TrimmedCurve2 = TrimmedCurve<Curve2>;

/// Any 2D curve.
#[doc(alias = "Geom2d_Curve")]
#[derive(Debug, Clone, PartialEq)]
pub enum Curve2 {
    Line(Line2),
    Circle(Circle2),
    Ellipse(Ellipse2),
    Hyperbola(Hyperbola2),
    Parabola(Parabola2),
    Bezier(BezierCurve2),
    BSpline(BSplineCurve2),
    Trimmed(TrimmedCurve2),
    Offset(OffsetCurve2),
}

// The `ParametricCurve` impl of `Curve2` is in eval.rs.

impl Curve2 {
    /// The curve placed in the XY plane of `frame` (lines, conics, Bezier and B-spline
    /// curves keep their kind; offsets and trimmed curves are mapped recursively).
    #[doc(alias = "GeomLib::To3d")]
    pub fn to_3d(&self, frame: &Frame) -> Curve {
        use super::curve::{Circle, Ellipse, Line};
        use crate::math::{Direction, Vector};
        let (x, y) = (Vector::from(frame.x()), Vector::from(frame.y()));
        let point = |p: Point2| frame.origin() + x * p.x + y * p.y;
        let direction = |d: Direction2| Direction::from_vector(x * d.x() + y * d.y());
        // A conic's frame: its X axis mapped, its normal the plane's (or the opposite one
        // for a clockwise conic).
        let conic_frame = |f: &Frame2| {
            let ccw = f.x().x() * f.y().y() - f.x().y() * f.y().x() > 0.0;
            let normal = if ccw { frame.z() } else { -frame.z() };
            Frame::new(point(f.origin()), normal, direction(f.x()))
        };
        match self {
            Curve2::Line(l) => Curve::Line(Line { origin: point(l.origin), direction: direction(l.direction) }),
            Curve2::Circle(c) => Curve::Circle(Circle::new(conic_frame(c.frame()), c.radius())),
            Curve2::Ellipse(e) => Curve::Ellipse(Ellipse::new(conic_frame(e.frame()), e.major_radius(), e.minor_radius())),
            Curve2::BSpline(b) => {
                let poles = b.poles().iter().map(|&p| point(p)).collect();
                let (knots, mults) = (b.knots().to_vec(), b.multiplicities().to_vec());
                let curve = if super::ParametricCurve::is_periodic(b) {
                    BSplineCurve::new_periodic(poles, knots, mults, b.degree())
                } else {
                    BSplineCurve::new(poles, knots, mults, b.degree())
                };
                Curve::BSpline(match b.weights() {
                    Some(w) => curve.with_weights(w.to_vec()),
                    None => curve,
                })
            }
            Curve2::Bezier(b) => {
                let poles = b.poles().iter().map(|&p| point(p)).collect();
                Curve::Bezier(match b.weights() {
                    Some(w) => BezierCurve::rational(poles, w.to_vec()),
                    None => BezierCurve::new(poles),
                })
            }
            other => todo!("placing a {other:?} in 3D"),
        }
    }
}

impl Transformable2 for Curve2 {
    fn transformed(&self, t: &Transform2) -> Self {
        match self {
            Curve2::Line(c) => Curve2::Line(c.transformed(t)),
            Curve2::Circle(c) => Curve2::Circle(c.transformed(t)),
            Curve2::Ellipse(c) => Curve2::Ellipse(c.transformed(t)),
            Curve2::Hyperbola(c) => Curve2::Hyperbola(c.transformed(t)),
            Curve2::Parabola(c) => Curve2::Parabola(c.transformed(t)),
            Curve2::Bezier(c) => Curve2::Bezier(c.transformed(t)),
            Curve2::BSpline(c) => Curve2::BSpline(c.transformed(t)),
            Curve2::Trimmed(c) => Curve2::Trimmed(c.transformed(t)),
            Curve2::Offset(c) => Curve2::Offset(c.transformed(t)),
        }
    }
}

impl Curve2 {
    /// The factor the curve's parameters are multiplied by when it is transformed by `t`
    /// (see [`Curve::parameter_scale`]).
    pub fn parameter_scale(&self, t: &Transform2) -> f64 {
        match self {
            Curve2::Line(_) | Curve2::Parabola(_) => t.scale_factor().abs(),
            Curve2::Trimmed(c) => c.basis().parameter_scale(t),
            Curve2::Offset(c) => c.basis().parameter_scale(t),
            _ => 1.0,
        }
    }
}

impl Transformable2 for TrimmedCurve2 {
    fn transformed(&self, t: &Transform2) -> Self {
        let k = self.basis().parameter_scale(t);
        let (first, last) = self.range();
        let mut c = self.clone();
        *c.basis_mut() = self.basis().transformed(t);
        c.set_range_unchecked(first * k, last * k);
        c
    }
}

/// An infinite straight line in 2D, parametrized by signed distance from `origin`.
#[doc(alias = "gp_Lin2d")]
#[doc(alias = "Geom2d_Line")]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Line2 {
    pub origin: Point2,
    pub direction: Direction2,
}

impl Line2 {
    pub const fn new(origin: Point2, direction: Direction2) -> Self {
        Self { origin, direction }
    }

    /// The line `a*x + b*y + c = 0`, with its origin at the point closest to (0, 0); `None`
    /// if `(a, b)` is (near) zero.
    pub fn from_equation(a: f64, b: f64, c: f64) -> Option<Self> {
        let n2 = a * a + b * b;
        if n2.sqrt() <= crate::tolerance::ANGULAR {
            return None;
        }
        let origin = Point2::new(-a * c / n2, -b * c / n2);
        Some(Line2::new(origin, Direction2::new(-b, a)))
    }

    pub fn axis(&self) -> Axis2 {
        Axis2::new(self.origin, self.direction)
    }

    pub fn distance_to_point(&self, p: Point2) -> f64 {
        (p - self.origin).cross(Vector2::from(self.direction)).abs()
    }
}

impl From<Axis2> for Line2 {
    fn from(axis: Axis2) -> Self {
        Line2::new(axis.origin, axis.direction)
    }
}

impl Transformable2 for Line2 {
    fn transformed(&self, t: &Transform2) -> Self {
        Line2::new(self.origin.transformed(t), self.direction.transformed(t))
    }
}

// The `ParametricCurve` impl of `Line2` is in eval.rs.

/// A circle centered at the origin of `frame`, parametrized by angle from the frame's X
/// axis, counter-clockwise for a counter-clockwise frame (period 2π).
#[doc(alias = "gp_Circ2d")]
#[doc(alias = "Geom2d_Circle")]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Circle2 {
    frame: Frame2,
    radius: f64,
}

impl Circle2 {
    /// Panics if `radius` is negative.
    pub fn new(frame: Frame2, radius: f64) -> Self {
        assert!(radius >= 0.0, "negative circle radius {radius}");
        Self { frame, radius }
    }

    pub fn frame(&self) -> &Frame2 {
        &self.frame
    }

    pub fn center(&self) -> Point2 {
        self.frame.origin()
    }

    pub fn radius(&self) -> f64 {
        self.radius
    }

    /// Panics if `radius` is negative.
    pub fn set_radius(&mut self, radius: f64) {
        assert!(radius >= 0.0, "negative circle radius {radius}");
        self.radius = radius;
    }
}

impl Transformable2 for Circle2 {
    fn transformed(&self, t: &Transform2) -> Self {
        Circle2::new(self.frame.transformed(t), self.radius * t.scale_factor().abs())
    }
}

// The `ParametricCurve` impl of `Circle2` is in eval.rs.

/// An ellipse centered at the origin of `frame`, major axis along the frame's X axis.
#[doc(alias = "gp_Elips2d")]
#[doc(alias = "Geom2d_Ellipse")]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Ellipse2 {
    frame: Frame2,
    major_radius: f64,
    minor_radius: f64,
}

impl Ellipse2 {
    /// Panics unless `0 <= minor_radius <= major_radius`.
    pub fn new(frame: Frame2, major_radius: f64, minor_radius: f64) -> Self {
        assert!(
            0.0 <= minor_radius && minor_radius <= major_radius,
            "invalid ellipse radii {major_radius}, {minor_radius}"
        );
        Self { frame, major_radius, minor_radius }
    }

    pub fn frame(&self) -> &Frame2 {
        &self.frame
    }

    pub fn center(&self) -> Point2 {
        self.frame.origin()
    }

    pub fn major_radius(&self) -> f64 {
        self.major_radius
    }

    pub fn minor_radius(&self) -> f64 {
        self.minor_radius
    }

    /// Panics if the radius would be smaller than the minor radius.
    pub fn set_major_radius(&mut self, radius: f64) {
        assert!(radius >= self.minor_radius, "major radius {radius} < minor radius");
        self.major_radius = radius;
    }

    /// Panics if the radius is negative or larger than the major radius.
    pub fn set_minor_radius(&mut self, radius: f64) {
        assert!((0.0..=self.major_radius).contains(&radius), "invalid minor radius {radius}");
        self.minor_radius = radius;
    }

    pub fn eccentricity(&self) -> f64 {
        if self.major_radius == 0.0 { 0.0 } else { self.focal_distance() / (2.0 * self.major_radius) }
    }

    pub fn focal_distance(&self) -> f64 {
        2.0 * (self.major_radius.powi(2) - self.minor_radius.powi(2)).sqrt()
    }

    /// The foci on the positive and negative side of the major axis.
    pub fn foci(&self) -> (Point2, Point2) {
        let offset = Vector2::from(self.frame.x()) * (self.focal_distance() / 2.0);
        (self.center() + offset, self.center() + offset * -1.0)
    }

    #[doc(alias = "Parameter")]
    pub fn semi_latus_rectum(&self) -> f64 {
        if self.major_radius == 0.0 { 0.0 } else { self.minor_radius.powi(2) / self.major_radius }
    }
}

impl Transformable2 for Ellipse2 {
    fn transformed(&self, t: &Transform2) -> Self {
        let k = t.scale_factor().abs();
        Ellipse2::new(self.frame.transformed(t), self.major_radius * k, self.minor_radius * k)
    }
}

// The `ParametricCurve` impl of `Ellipse2` is in eval.rs.

/// The branch of a hyperbola centered at the origin of `frame` crossing its positive X
/// axis: `t -> major·cosh t·x + minor·sinh t·y`.
#[doc(alias = "gp_Hypr2d")]
#[doc(alias = "Geom2d_Hyperbola")]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Hyperbola2 {
    frame: Frame2,
    major_radius: f64,
    minor_radius: f64,
}

impl Hyperbola2 {
    /// Panics if a radius is negative.
    pub fn new(frame: Frame2, major_radius: f64, minor_radius: f64) -> Self {
        assert!(
            major_radius >= 0.0 && minor_radius >= 0.0,
            "invalid hyperbola radii {major_radius}, {minor_radius}"
        );
        Self { frame, major_radius, minor_radius }
    }

    pub fn frame(&self) -> &Frame2 {
        &self.frame
    }

    pub fn center(&self) -> Point2 {
        self.frame.origin()
    }

    pub fn major_radius(&self) -> f64 {
        self.major_radius
    }

    pub fn minor_radius(&self) -> f64 {
        self.minor_radius
    }

    /// Panics if negative.
    pub fn set_major_radius(&mut self, radius: f64) {
        assert!(radius >= 0.0, "negative radius {radius}");
        self.major_radius = radius;
    }

    /// Panics if negative.
    pub fn set_minor_radius(&mut self, radius: f64) {
        assert!(radius >= 0.0, "negative radius {radius}");
        self.minor_radius = radius;
    }

    pub fn eccentricity(&self) -> f64 {
        self.major_radius.hypot(self.minor_radius) / self.major_radius
    }

    /// The distance between the foci.
    pub fn focal_distance(&self) -> f64 {
        2.0 * self.major_radius.hypot(self.minor_radius)
    }

    pub fn foci(&self) -> (Point2, Point2) {
        let c = Vector2::from(self.frame.x()) * self.major_radius.hypot(self.minor_radius);
        (self.center() + c, self.center() + c * -1.0)
    }

    #[doc(alias = "Parameter")]
    pub fn semi_latus_rectum(&self) -> f64 {
        self.minor_radius * self.minor_radius / self.major_radius
    }

    /// The two asymptotes, through the center.
    pub fn asymptotes(&self) -> (Axis2, Axis2) {
        let (x, y) = (Vector2::from(self.frame.x()) * self.major_radius, Vector2::from(self.frame.y()) * self.minor_radius);
        let dir = |v: Vector2| Direction2::new(v.x, v.y);
        (Axis2::new(self.center(), dir(x + y)), Axis2::new(self.center(), dir(x - y)))
    }
}

impl Transformable2 for Hyperbola2 {
    fn transformed(&self, t: &Transform2) -> Self {
        let s = t.scale_factor().abs();
        Hyperbola2::new(self.frame.transformed(t), self.major_radius * s, self.minor_radius * s)
    }
}


/// A parabola with apex at the origin of `frame`, opening along its X axis, focus at
/// distance `focal`: `t -> t²/(4·focal)·x + t·y`.
///
/// A zero focal length gives a degenerate parabola: the line along the frame's X axis,
/// `t -> t·x`.
#[doc(alias = "gp_Parab2d")]
#[doc(alias = "Geom2d_Parabola")]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Parabola2 {
    frame: Frame2,
    focal: f64,
}

impl Parabola2 {
    /// Panics if `focal` is negative.
    pub fn new(frame: Frame2, focal: f64) -> Self {
        assert!(focal >= 0.0, "negative parabola focal length {focal}");
        Self { frame, focal }
    }

    pub fn frame(&self) -> &Frame2 {
        &self.frame
    }

    pub fn apex(&self) -> Point2 {
        self.frame.origin()
    }

    pub fn focal(&self) -> f64 {
        self.focal
    }

    /// Panics if negative.
    pub fn set_focal(&mut self, focal: f64) {
        assert!(focal >= 0.0, "negative parabola focal length {focal}");
        self.focal = focal;
    }

    pub fn focus(&self) -> Point2 {
        self.apex() + Vector2::from(self.frame.x()) * self.focal
    }

    pub fn directrix(&self) -> Axis2 {
        Axis2::new(self.apex() + Vector2::from(self.frame.x()) * -self.focal, self.frame.y())
    }

    #[doc(alias = "Parameter")]
    pub fn semi_latus_rectum(&self) -> f64 {
        2.0 * self.focal
    }

    /// Coefficients `[a, b, c, d, e, f]` of the implicit equation
    /// `a·x² + b·y² + 2c·xy + 2d·x + 2e·y + f = 0`, scaled so that `a + b = 1`.
    #[doc(alias = "Coefficients")]
    pub fn implicit_coefficients(&self) -> [f64; 6] {
        // In local coordinates y'² - 4f x' = 0, with x' = X·(p - o), y' = Y·(p - o).
        let (o, x, y) = (self.frame.origin(), self.frame.x(), self.frame.y());
        let (xx, xy, yx, yy) = (x.x(), x.y(), y.x(), y.y());
        let k = yx * o.x + yy * o.y;
        let m = xx * o.x + xy * o.y;
        let f = self.focal;
        [yx * yx, yy * yy, yx * yy, -k * yx - 2.0 * f * xx, -k * yy - 2.0 * f * xy, k * k + 4.0 * f * m]
    }
}

impl Transformable2 for Parabola2 {
    fn transformed(&self, t: &Transform2) -> Self {
        Parabola2::new(self.frame.transformed(t), self.focal * t.scale_factor().abs())
    }
}


/// The curve at signed distance `offset` from a basis curve, to the right of the direction
/// of travel (along the tangent turned clockwise). Undefined where the basis tangent
/// vanishes.
#[doc(alias = "Geom2d_OffsetCurve")]
#[doc(alias = "Adaptor2d_OffsetCurve")]
#[derive(Debug, Clone, PartialEq)]
pub struct OffsetCurve2 {
    basis: Box<Curve2>,
    offset: f64,
}

impl OffsetCurve2 {
    pub fn new(basis: Curve2, offset: f64) -> Self {
        Self { basis: Box::new(basis), offset }
    }

    pub fn basis(&self) -> &Curve2 {
        &self.basis
    }

    pub fn offset(&self) -> f64 {
        self.offset
    }

    pub fn set_offset(&mut self, offset: f64) {
        self.offset = offset;
    }
}

impl Transformable2 for OffsetCurve2 {
    fn transformed(&self, t: &Transform2) -> Self {
        // A reflection turns the offset side over.
        let sign = if t.is_negative() { -1.0 } else { 1.0 };
        OffsetCurve2::new(self.basis.transformed(t), sign * self.offset * t.scale_factor().abs())
    }
}

