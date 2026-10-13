//! 3D curves.

use crate::math::{Axis, Direction, Frame, Point, Transform, Transformable, Vector};

use super::ParametricCurve;

pub use super::bezier::BezierCurve;
pub use super::bspline::BSplineCurve;

/// Any 3D curve.
#[doc(alias = "Geom_Curve")]
#[derive(Debug, Clone, PartialEq)]
pub enum Curve {
    Line(Line),
    Circle(Circle),
    Ellipse(Ellipse),
    Hyperbola(Hyperbola),
    Parabola(Parabola),
    Bezier(BezierCurve),
    BSpline(BSplineCurve),
    Trimmed(TrimmedCurve),
    Offset(OffsetCurve),
}

// The `ParametricCurve` impl of `Curve` is in eval.rs.

impl Transformable for Curve {
    fn transformed(&self, t: &Transform) -> Self {
        match self {
            Curve::Line(c) => Curve::Line(c.transformed(t)),
            Curve::Circle(c) => Curve::Circle(c.transformed(t)),
            Curve::Ellipse(c) => Curve::Ellipse(c.transformed(t)),
            Curve::Hyperbola(c) => Curve::Hyperbola(c.transformed(t)),
            Curve::Parabola(c) => Curve::Parabola(c.transformed(t)),
            Curve::Bezier(c) => Curve::Bezier(c.transformed(t)),
            Curve::BSpline(c) => Curve::BSpline(c.transformed(t)),
            Curve::Trimmed(c) => Curve::Trimmed(c.transformed(t)),
            Curve::Offset(c) => Curve::Offset(c.transformed(t)),
        }
    }
}

impl Curve {
    /// The factor a curve's parameters are multiplied by when it is transformed by `t`
    /// (OCCT's `TransformedParameter`): a distance along a line or a parabola's axis
    /// scales, an angle or a B-spline's parameter does not.
    #[doc(alias = "Geom_Curve::TransformedParameter")]
    pub fn parameter_scale(&self, t: &Transform) -> f64 {
        match self {
            Curve::Line(_) | Curve::Parabola(_) => t.scale_factor().abs(),
            Curve::Trimmed(c) => c.basis().parameter_scale(t),
            Curve::Offset(c) => c.basis().parameter_scale(t),
            _ => 1.0,
        }
    }
}

/// The part of a curve between two parameters, keeping the basis' parametrization.
/// Generic over the basis so that it serves 3D ([`Curve`]) and 2D
/// ([`Curve2`](super::curve2::Curve2)) curves.
#[doc(alias = "Geom_TrimmedCurve")]
#[doc(alias = "Geom2d_TrimmedCurve")]
#[derive(Debug, Clone, PartialEq)]
pub struct TrimmedCurve<C = Curve> {
    basis: Box<C>,
    first: f64,
    last: f64,
}

impl<C: ParametricCurve> TrimmedCurve<C> {
    /// `basis` restricted to `[first, last]`. For a periodic basis, `last` is shifted by
    /// periods to lie in `(first, first + period]`. For a non-periodic basis the range may
    /// be a single parameter (a range reversed by at most `tolerance::LINEAR` counts as
    /// one); panics if it is reversed by more, or outside the basis' domain. To run a curve
    /// backwards, reverse it first (unlike OCCT's `Geom_TrimmedCurve`, a reversed range is
    /// not swapped).
    pub fn new(basis: C, first: f64, last: f64) -> Self {
        let (first, last) = Self::checked_range(&basis, first, last);
        Self { basis: Box::new(basis), first, last }
    }

    pub fn basis(&self) -> &C {
        &self.basis
    }

    /// Changes the trimming range, with the same rules as [`TrimmedCurve::new`].
    pub fn set_trim(&mut self, first: f64, last: f64) {
        (self.first, self.last) = Self::checked_range(&self.basis, first, last);
    }

    /// The trimming range `(first, last)`.
    pub fn range(&self) -> (f64, f64) {
        (self.first, self.last)
    }

    pub(super) fn basis_mut(&mut self) -> &mut C {
        &mut self.basis
    }

    pub(super) fn set_range_unchecked(&mut self, first: f64, last: f64) {
        (self.first, self.last) = (first, last);
    }
}

// The `ParametricCurve` impl of `TrimmedCurve` is in eval.rs.

impl Transformable for TrimmedCurve {
    fn transformed(&self, t: &Transform) -> Self {
        // The basis' parameters change with a scaling for lines and parabolas.
        let k = self.basis.parameter_scale(t);
        Self { basis: Box::new(self.basis.transformed(t)), first: self.first * k, last: self.last * k }
    }
}

/// An infinite straight line, parametrized by signed distance from `origin`.
#[doc(alias = "gp_Lin")]
#[doc(alias = "Geom_Line")]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Line {
    pub origin: Point,
    pub direction: Direction,
}

impl Line {
    pub const fn new(origin: Point, direction: Direction) -> Self {
        Self { origin, direction }
    }

    pub fn axis(&self) -> Axis {
        Axis::new(self.origin, self.direction)
    }

    /// Angle in [0, PI/2] between the lines.
    pub fn angle(&self, other: &Line) -> f64 {
        let (a, b) = (Vector::from(self.direction), Vector::from(other.direction));
        a.cross(b).length().atan2(a.dot(b).abs())
    }

    pub fn distance_to_point(&self, p: Point) -> f64 {
        (p - self.origin).cross(Vector::from(self.direction)).length()
    }

    pub fn distance_to_line(&self, other: &Line) -> f64 {
        let normal = Vector::from(self.direction).cross(Vector::from(other.direction));
        let len = normal.length();
        if len <= crate::tolerance::ANGULAR {
            return self.distance_to_point(other.origin);
        }
        (other.origin - self.origin).dot(normal).abs() / len
    }

    pub fn contains(&self, p: Point, tolerance: f64) -> bool {
        self.distance_to_point(p) <= tolerance
    }
}

impl From<Axis> for Line {
    fn from(axis: Axis) -> Self {
        Line::new(axis.origin, axis.direction)
    }
}

/// The image line, with the parameter scaled like distances (OCCT's convention).
impl Transformable for Line {
    fn transformed(&self, t: &Transform) -> Self {
        Line::new(self.origin.transformed(t), self.direction.transformed(t))
    }
}

// The `ParametricCurve` impl of `Line` is in eval.rs.

/// A circle in the XY plane of `frame`, centered at its origin, parametrized by angle from
/// the frame's X axis (period 2π).
#[doc(alias = "gp_Circ")]
#[doc(alias = "Geom_Circle")]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Circle {
    frame: Frame,
    radius: f64,
}

impl Circle {
    /// Panics if `radius` is negative.
    pub fn new(frame: Frame, radius: f64) -> Self {
        assert!(radius >= 0.0, "negative circle radius {radius}");
        Self { frame, radius }
    }

    pub fn frame(&self) -> &Frame {
        &self.frame
    }

    pub fn center(&self) -> Point {
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

    pub fn area(&self) -> f64 {
        std::f64::consts::PI * self.radius * self.radius
    }

    pub fn circumference(&self) -> f64 {
        std::f64::consts::TAU * self.radius
    }

    /// The arc from angle `start` to angle `end`, counter-clockwise about the frame's Z
    /// axis (`ccw`), or clockwise.
    #[doc(alias = "GC_MakeArcOfCircle")]
    pub fn arc(&self, start: f64, end: f64, ccw: bool) -> TrimmedCurve {
        if ccw {
            return TrimmedCurve::new(Curve::Circle(*self), start, end);
        }
        // Clockwise: on the circle with its Y axis flipped, angle a is at 2π - a.
        let mut frame = self.frame;
        frame.flip_y();
        let tau = std::f64::consts::TAU;
        TrimmedCurve::new(Curve::Circle(Circle::new(frame, self.radius)), tau - start, tau - end)
    }
}

/// `frame` with its Z axis `X × Y`: a conic's frame moved by a reflection keeps its X and Y
/// axes (and so its parametrization) and turns its axis to stay the one it runs
/// counter-clockwise about (as OCCT's gp_Ax2 does).
fn right_handed(frame: Frame) -> Frame {
    let z = Vector::from(frame.x()).cross(Vector::from(frame.y()));
    if z.dot(Vector::from(frame.z())) >= 0.0 {
        return frame;
    }
    Frame::from_axes_unchecked(frame.origin(), frame.x(), frame.y(), Direction::from_vector(z))
}

impl Transformable for Circle {
    fn transformed(&self, t: &Transform) -> Self {
        Circle::new(right_handed(self.frame.transformed(t)), self.radius * t.scale_factor().abs())
    }
}

// The `ParametricCurve` impl of `Circle` is in eval.rs.

/// An ellipse in the XY plane of `frame`, centered at its origin, with its major axis
/// along the frame's X axis: `t -> major·cos t·x + minor·sin t·y` (period 2π).
#[doc(alias = "gp_Elips")]
#[doc(alias = "Geom_Ellipse")]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Ellipse {
    frame: Frame,
    major_radius: f64,
    minor_radius: f64,
}

impl Ellipse {
    /// Panics unless `0 <= minor_radius <= major_radius`.
    pub fn new(frame: Frame, major_radius: f64, minor_radius: f64) -> Self {
        assert!(
            0.0 <= minor_radius && minor_radius <= major_radius,
            "invalid ellipse radii {major_radius}, {minor_radius}"
        );
        Self { frame, major_radius, minor_radius }
    }

    pub fn frame(&self) -> &Frame {
        &self.frame
    }

    pub fn center(&self) -> Point {
        self.frame.origin()
    }

    pub fn major_radius(&self) -> f64 {
        self.major_radius
    }

    pub fn minor_radius(&self) -> f64 {
        self.minor_radius
    }

    pub fn eccentricity(&self) -> f64 {
        if self.major_radius == 0.0 { 0.0 } else { self.focal_distance() / (2.0 * self.major_radius) }
    }

    /// Distance between the foci.
    pub fn focal_distance(&self) -> f64 {
        2.0 * (self.major_radius.powi(2) - self.minor_radius.powi(2)).sqrt()
    }

    /// The foci on the positive and negative side of the major axis.
    pub fn foci(&self) -> (Point, Point) {
        let offset = Vector::from(self.frame.x()) * (self.focal_distance() / 2.0);
        (self.center() + offset, self.center() - offset)
    }

    /// The semi-latus rectum `minor² / major`.
    #[doc(alias = "Parameter")]
    pub fn semi_latus_rectum(&self) -> f64 {
        if self.major_radius == 0.0 { 0.0 } else { self.minor_radius.powi(2) / self.major_radius }
    }
}

impl Transformable for Ellipse {
    fn transformed(&self, t: &Transform) -> Self {
        let k = t.scale_factor().abs();
        Ellipse::new(right_handed(self.frame.transformed(t)), self.major_radius * k, self.minor_radius * k)
    }
}

// The `ParametricCurve` impl of `Ellipse` is in eval.rs.

/// The branch of a hyperbola in the XY plane of `frame` that crosses the positive X axis:
/// `t -> major·cosh t·x + minor·sinh t·y` (unbounded).
#[doc(alias = "gp_Hypr")]
#[doc(alias = "Geom_Hyperbola")]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Hyperbola {
    frame: Frame,
    major_radius: f64,
    minor_radius: f64,
}

impl Hyperbola {
    /// Panics if a radius is negative.
    pub fn new(frame: Frame, major_radius: f64, minor_radius: f64) -> Self {
        assert!(
            major_radius >= 0.0 && minor_radius >= 0.0,
            "invalid hyperbola radii {major_radius}, {minor_radius}"
        );
        Self { frame, major_radius, minor_radius }
    }

    pub fn frame(&self) -> &Frame {
        &self.frame
    }

    pub fn center(&self) -> Point {
        self.frame.origin()
    }

    pub fn major_radius(&self) -> f64 {
        self.major_radius
    }

    pub fn minor_radius(&self) -> f64 {
        self.minor_radius
    }

    pub fn eccentricity(&self) -> f64 {
        self.major_radius.hypot(self.minor_radius) / self.major_radius
    }

    /// The distance between the foci.
    pub fn focal_distance(&self) -> f64 {
        2.0 * self.major_radius.hypot(self.minor_radius)
    }

    pub fn foci(&self) -> (Point, Point) {
        let c = Vector::from(self.frame.x()) * self.major_radius.hypot(self.minor_radius);
        (self.center() + c, self.center() - c)
    }

    #[doc(alias = "Parameter")]
    pub fn semi_latus_rectum(&self) -> f64 {
        self.minor_radius * self.minor_radius / self.major_radius
    }

    /// The two asymptotes, through the center.
    pub fn asymptotes(&self) -> (Axis, Axis) {
        let (x, y) = (Vector::from(self.frame.x()) * self.major_radius, Vector::from(self.frame.y()) * self.minor_radius);
        (Axis::new(self.center(), Direction::from_vector(x + y)), Axis::new(self.center(), Direction::from_vector(x - y)))
    }
}

impl Transformable for Hyperbola {
    fn transformed(&self, t: &Transform) -> Self {
        let s = t.scale_factor().abs();
        Hyperbola::new(right_handed(self.frame.transformed(t)), self.major_radius * s, self.minor_radius * s)
    }
}


/// A parabola in the XY plane of `frame` with its apex at the origin, opening along the
/// frame's X axis, focus at distance `focal` from the apex:
/// `t -> t²/(4·focal)·x + t·y` (unbounded).
///
/// A zero focal length gives a degenerate parabola: the line along the frame's X axis,
/// `t -> t·x`.
#[doc(alias = "gp_Parab")]
#[doc(alias = "Geom_Parabola")]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Parabola {
    frame: Frame,
    focal: f64,
}

impl Parabola {
    /// Panics if `focal` is negative.
    pub fn new(frame: Frame, focal: f64) -> Self {
        assert!(focal >= 0.0, "negative parabola focal length {focal}");
        Self { frame, focal }
    }

    pub fn frame(&self) -> &Frame {
        &self.frame
    }

    pub fn apex(&self) -> Point {
        self.frame.origin()
    }

    /// Distance from the apex to the focus.
    pub fn focal(&self) -> f64 {
        self.focal
    }

    pub fn focus(&self) -> Point {
        self.apex() + Vector::from(self.frame.x()) * self.focal
    }

    /// The directrix (perpendicular to the X axis, at `-focal`).
    pub fn directrix(&self) -> Axis {
        Axis::new(self.apex() - Vector::from(self.frame.x()) * self.focal, self.frame.y())
    }

    /// `2·focal`.
    #[doc(alias = "Parameter")]
    pub fn semi_latus_rectum(&self) -> f64 {
        2.0 * self.focal
    }
}

impl Transformable for Parabola {
    fn transformed(&self, t: &Transform) -> Self {
        Parabola::new(right_handed(self.frame.transformed(t)), self.focal * t.scale_factor().abs())
    }
}


/// The curve at constant distance `offset` from a basis curve, measured along
/// `basis' × direction` (to the right of the basis seen from `direction`). Where the basis
/// tangent is parallel to `direction` the offset is undefined.
#[doc(alias = "Geom_OffsetCurve")]
#[derive(Debug, Clone, PartialEq)]
pub struct OffsetCurve {
    basis: Box<Curve>,
    offset: f64,
    direction: Direction,
}

impl OffsetCurve {
    pub fn new(basis: Curve, offset: f64, direction: Direction) -> Self {
        Self { basis: Box::new(basis), offset, direction }
    }

    pub fn basis(&self) -> &Curve {
        &self.basis
    }

    pub fn offset(&self) -> f64 {
        self.offset
    }

    pub fn set_offset(&mut self, offset: f64) {
        self.offset = offset;
    }

    pub fn direction(&self) -> Direction {
        self.direction
    }
}

impl Transformable for OffsetCurve {
    fn transformed(&self, t: &Transform) -> Self {
        // A reflection turns the offset side over.
        let sign = if t.is_negative() { -1.0 } else { 1.0 };
        OffsetCurve::new(self.basis.transformed(t), sign * self.offset * t.scale_factor().abs(), self.direction.transformed(t))
    }
}

