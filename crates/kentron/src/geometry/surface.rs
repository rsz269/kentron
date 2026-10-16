//! Surfaces.

use crate::math::{Axis, Direction, Frame, Point, Transform, Transformable, Vector};

use super::curve::Curve;

pub use super::bezier::BezierSurface;
pub use super::bspline::{BSplineSurface, KnotVector};

/// Any surface.
#[doc(alias = "Geom_Surface")]
#[derive(Debug, Clone, PartialEq)]
pub enum Surface {
    Plane(Plane),
    Cylinder(Cylinder),
    Cone(Cone),
    Sphere(Sphere),
    Torus(Torus),
    Bezier(BezierSurface),
    BSpline(BSplineSurface),
    Extrusion(ExtrusionSurface),
    Revolution(RevolutionSurface),
    Offset(OffsetSurface),
    Trimmed(TrimmedSurface),
}

// The `ParametricSurface` impl of `Surface` is in eval.rs.

impl Surface {
    /// What identifies a B-spline or Bézier surface and its unchanged clones (for caches
    /// keyed by surface); `None` for others.
    pub(crate) fn shared_id(&self) -> Option<super::SurfaceId> {
        match self {
            Surface::BSpline(b) => Some(b.id()),
            Surface::Bezier(b) => Some(b.id()),
            _ => None,
        }
    }

    /// The point and first derivatives at `(u, v)` for each `u` of `us` (best ascending), as
    /// `d1` gives them one by one; at once (for less work) on a B-spline surface.
    pub fn d1_along_u(&self, us: &[f64], v: f64) -> Vec<(Point, Vector, Vector)> {
        use super::ParametricSurface;
        match self {
            Surface::BSpline(b) => b.d1_along_u(us, v),
            _ => us.iter().map(|&u| self.d1(u, v)).collect(),
        }
    }

    /// The points at `(u, v)` for each `u` of `us` (best ascending), as `point` gives them
    /// one by one.
    pub fn points_along_u(&self, us: &[f64], v: f64) -> Vec<Point> {
        use super::{ParametricCurve, ParametricSurface};
        match self {
            Surface::BSpline(b) => b.points_along_u(us, v),
            // (The basis point at v once, turned for each u.)
            Surface::Revolution(r) => {
                let at = r.basis().point(v);
                us.iter().map(|&u| at.transformed(&Transform::rotation(r.axis(), u))).collect()
            }
            _ => us.iter().map(|&u| self.point(u, v)).collect(),
        }
    }

    /// The points at `(u, v)` for each `v` of `vs` (a row each) and each `u` of `us`, as
    /// [`Self::points_along_u`] gives them row by row; an extrusion's basis points found once
    /// for all rows.
    pub fn points_on_grid(&self, us: &[f64], vs: &[f64]) -> Vec<Vec<Point>> {
        use super::ParametricCurve;
        match self {
            Surface::Extrusion(e) => {
                let at: Vec<Point> = us.iter().map(|&u| e.basis().point(u)).collect();
                let direction = Vector::from(e.direction());
                vs.iter().map(|&v| at.iter().map(|&p| p + direction * v).collect()).collect()
            }
            _ => vs.iter().map(|&v| self.points_along_u(us, v)).collect(),
        }
    }

    /// The parameters `(u, v)` of a point on (or near) an elementary surface; `None` for
    /// other surfaces. Angular parameters are in `[0, 2π)`.
    pub fn parameters(&self, p: Point) -> Option<(f64, f64)> {
        Some(match self {
            Surface::Plane(s) => s.parameters(p),
            Surface::Cylinder(s) => s.parameters(p),
            Surface::Cone(s) => s.parameters(p),
            Surface::Sphere(s) => s.parameters(p),
            Surface::Torus(s) => s.parameters(p),
            Surface::Revolution(s) => s.parameters(p),
            Surface::Extrusion(s) => s.parameters(p),
            _ => return None,
        })
    }
}

impl Transformable for Surface {
    fn transformed(&self, t: &Transform) -> Self {
        match self {
            Surface::Plane(s) => Surface::Plane(s.transformed(t)),
            Surface::Cylinder(s) => Surface::Cylinder(s.transformed(t)),
            Surface::Cone(s) => Surface::Cone(s.transformed(t)),
            Surface::Sphere(s) => Surface::Sphere(s.transformed(t)),
            Surface::Torus(s) => Surface::Torus(s.transformed(t)),
            Surface::Revolution(s) => Surface::Revolution(s.transformed(t)),
            Surface::BSpline(s) => Surface::BSpline(s.transformed(t)),
            Surface::Extrusion(s) => Surface::Extrusion(s.transformed(t)),
            Surface::Bezier(s) => Surface::Bezier(s.transformed(t)),
            Surface::Offset(s) => Surface::Offset(s.transformed(t)),
            Surface::Trimmed(s) => Surface::Trimmed(s.transformed(t)),
        }
    }
}

/// The coordinates of `p` in `frame`.
fn local_coordinates(frame: &Frame, p: Point) -> (f64, f64, f64) {
    let d = p - frame.origin();
    (d.dot(Vector::from(frame.x())), d.dot(Vector::from(frame.y())), d.dot(Vector::from(frame.z())))
}

/// An angle from `atan2` (in `[-π, π]`) moved into `[0, 2π)`.
fn angle_in_period(angle: f64) -> f64 {
    // As OCCT's ElSLib: an angle a rounding error below 0 is 0, not 2π.
    if angle < -f64::EPSILON {
        angle + std::f64::consts::TAU
    } else {
        angle.max(0.0)
    }
}

/// An infinite plane: the XY plane of `frame`, parametrized by the frame's X and Y axes.
#[doc(alias = "gp_Pln")]
#[doc(alias = "Geom_Plane")]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Plane {
    frame: Frame,
}

impl Plane {
    /// The XY plane of the world frame.
    #[doc(alias = "gp::XOY")]
    pub const XY: Plane = Plane::from_frame(Frame::WORLD);

    /// The plane through `origin` with normal `normal`.
    pub fn new(origin: Point, normal: Direction) -> Self {
        Self::from_frame(Frame::from_z(origin, normal))
    }

    pub const fn from_frame(frame: Frame) -> Self {
        Self { frame }
    }

    /// The plane `a*x + b*y + c*z + d = 0`; `None` if `(a, b, c)` is (near) zero.
    pub fn from_equation(a: f64, b: f64, c: f64, d: f64) -> Option<Self> {
        // As OCCT's gp_Pln: the origin on the axis of the largest coefficient but one.
        let normal = Direction::try_from_vector(Vector::new(a, b, c))?;
        let (aa, ab, ac) = (a.abs(), b.abs(), c.abs());
        let (origin, x) = if ab <= aa && ab <= ac {
            if aa > ac { (Point::new(-d / a, 0.0, 0.0), Vector::new(-c, 0.0, a)) } else { (Point::new(0.0, 0.0, -d / c), Vector::new(c, 0.0, -a)) }
        } else if aa <= ab && aa <= ac {
            if ab > ac { (Point::new(0.0, -d / b, 0.0), Vector::new(0.0, -c, b)) } else { (Point::new(0.0, 0.0, -d / c), Vector::new(0.0, c, -b)) }
        } else if aa > ab {
            (Point::new(-d / a, 0.0, 0.0), Vector::new(-b, a, 0.0))
        } else {
            (Point::new(0.0, -d / b, 0.0), Vector::new(b, -a, 0.0))
        };
        Some(Plane::from_frame(Frame::new(origin, normal, Direction::try_from_vector(x)?)))
    }

    /// Coefficients `(a, b, c, d)` of the plane's equation, with `(a, b, c)` the unit normal.
    pub fn equation(&self) -> (f64, f64, f64, f64) {
        let n = self.frame().z();
        let o = self.frame().origin();
        (n.x(), n.y(), n.z(), -(n.x() * o.x + n.y() * o.y + n.z() * o.z))
    }

    pub fn frame(&self) -> &Frame {
        &self.frame
    }

    pub fn origin(&self) -> Point {
        self.frame.origin()
    }

    pub fn normal(&self) -> Direction {
        self.frame.z()
    }

    /// Unsigned distance from `p` to the plane.
    pub fn distance(&self, p: Point) -> f64 {
        (p - self.origin()).dot(Vector::from(self.normal())).abs()
    }

    pub fn contains(&self, p: Point, tolerance: f64) -> bool {
        self.distance(p) <= tolerance
    }

    /// The parameters `(u, v)` of `p`, which should lie on the surface (for other points,
    /// those of a nearby surface point; `u` in `[0, 2π)` for an angular parameter).
    #[doc(alias = "ElSLib::Parameters")]
    pub fn parameters(&self, p: Point) -> (f64, f64) {
        let d = p - self.origin();
        (d.dot(Vector::from(self.frame.x())), d.dot(Vector::from(self.frame.y())))
    }
}

impl Transformable for Plane {
    fn transformed(&self, t: &Transform) -> Self {
        Plane::from_frame(self.frame.transformed(t))
    }
}

// The `ParametricSurface` impl of `Plane` is in eval.rs.

/// An infinite cylinder of the given radius around the Z axis of `frame`, parametrized by
/// `u` (angle from the frame's X axis) and `v` (height along Z).
#[doc(alias = "gp_Cylinder")]
#[doc(alias = "Geom_CylindricalSurface")]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Cylinder {
    frame: Frame,
    radius: f64,
}

impl Cylinder {
    /// Panics if `radius` is negative.
    pub fn new(frame: Frame, radius: f64) -> Self {
        assert!(radius >= 0.0, "negative cylinder radius {radius}");
        Self { frame, radius }
    }

    pub fn frame(&self) -> &Frame {
        &self.frame
    }

    pub fn radius(&self) -> f64 {
        self.radius
    }

    /// The parameters `(u, v)` of `p`, which should lie on the surface (for other points,
    /// those of a nearby surface point; `u` in `[0, 2π)` for an angular parameter).
    #[doc(alias = "ElSLib::Parameters")]
    pub fn parameters(&self, p: Point) -> (f64, f64) {
        let (x, y, z) = local_coordinates(&self.frame, p);
        (angle_in_period(y.atan2(x)), z)
    }
}

impl Transformable for Cylinder {
    fn transformed(&self, t: &Transform) -> Self {
        Cylinder::new(self.frame.transformed(t), self.radius * t.scale_factor().abs())
    }
}

// The `ParametricSurface` impl of `Cylinder` is in eval.rs.

/// An infinite cone around the Z axis of `frame`, with radius `radius` in the frame's XY
/// plane and half-angle `semi_angle` (radians, in (-PI/2, PI/2), non-zero). Parametrized
/// by `u` (angle) and `v` (distance along a generating line).
#[doc(alias = "gp_Cone")]
#[doc(alias = "Geom_ConicalSurface")]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Cone {
    frame: Frame,
    semi_angle: f64,
    radius: f64,
    /// The semi-angle's sine and cosine (every evaluation needs them).
    semi_angle_sin_cos: (f64, f64),
}

impl Cone {
    /// Panics on an invalid angle or a negative radius.
    pub fn new(frame: Frame, semi_angle: f64, radius: f64) -> Self {
        assert!(
            semi_angle != 0.0 && semi_angle.abs() < std::f64::consts::FRAC_PI_2,
            "invalid cone semi-angle {semi_angle}"
        );
        assert!(radius >= 0.0, "negative cone radius {radius}");
        Self { frame, semi_angle, radius, semi_angle_sin_cos: semi_angle.sin_cos() }
    }

    pub fn frame(&self) -> &Frame {
        &self.frame
    }

    /// The semi-angle's sine and cosine.
    pub fn semi_angle_sin_cos(&self) -> (f64, f64) {
        self.semi_angle_sin_cos
    }

    pub fn semi_angle(&self) -> f64 {
        self.semi_angle
    }

    pub fn radius(&self) -> f64 {
        self.radius
    }

    /// The parameters `(u, v)` of `p`, which should lie on the surface (for other points,
    /// those of a nearby surface point; `u` in `[0, 2π)` for an angular parameter).
    #[doc(alias = "ElSLib::Parameters")]
    pub fn parameters(&self, p: Point) -> (f64, f64) {
        let (x, y, z) = local_coordinates(&self.frame, p);
        let radial = x.hypot(y);
        let u = if radial > 0.0 { angle_in_period(y.atan2(x)) } else { 0.0 };
        let (sin, cos) = self.semi_angle.sin_cos();
        // On the generator at `u`, or past the apex on the one opposite (`v` beyond the
        // apex's, the radius negative): whichever the point is nearer.
        let off = |signed_radial: f64| ((signed_radial - self.radius) * cos - z * sin).abs();
        if radial > 0.0 && off(-radial) < off(radial) {
            (angle_in_period(u + std::f64::consts::PI), (-radial - self.radius) * sin + z * cos)
        } else {
            (u, (radial - self.radius) * sin + z * cos)
        }
    }

    /// The apex: where the generating lines meet, at `v = -radius / sin(semi_angle)`.
    pub fn apex(&self) -> Point {
        let v = -self.radius / self.semi_angle.sin();
        self.frame.origin() + Vector::from(self.frame.z()) * (v * self.semi_angle.cos())
    }
}

impl Transformable for Cone {
    fn transformed(&self, t: &Transform) -> Self {
        Cone::new(self.frame.transformed(t), self.semi_angle, self.radius * t.scale_factor().abs())
    }
}

// The `ParametricSurface` impl of `Cone` is in eval.rs.

/// A sphere centered at the origin of `frame`, parametrized by `u` (longitude from the X
/// axis) and `v` (latitude from the XY plane, in [-PI/2, PI/2]).
#[doc(alias = "gp_Sphere")]
#[doc(alias = "Geom_SphericalSurface")]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Sphere {
    frame: Frame,
    radius: f64,
}

impl Sphere {
    /// Panics if `radius` is negative.
    pub fn new(frame: Frame, radius: f64) -> Self {
        assert!(radius >= 0.0, "negative sphere radius {radius}");
        Self { frame, radius }
    }

    pub fn frame(&self) -> &Frame {
        &self.frame
    }

    pub fn radius(&self) -> f64 {
        self.radius
    }

    /// The parameters `(u, v)` of `p`, which should lie on the surface (for other points,
    /// those of a nearby surface point; `u` in `[0, 2π)` for an angular parameter).
    #[doc(alias = "ElSLib::Parameters")]
    pub fn parameters(&self, p: Point) -> (f64, f64) {
        let (x, y, z) = local_coordinates(&self.frame, p);
        let radial = x.hypot(y);
        let u = if radial > 0.0 { angle_in_period(y.atan2(x)) } else { 0.0 };
        (u, z.atan2(radial))
    }
}

impl Transformable for Sphere {
    fn transformed(&self, t: &Transform) -> Self {
        Sphere::new(self.frame.transformed(t), self.radius * t.scale_factor().abs())
    }
}

// The `ParametricSurface` impl of `Sphere` is in eval.rs.

/// The part of a surface over the parameter rectangle `[u1, u2] x [v1, v2]`.
#[doc(alias = "Geom_RectangularTrimmedSurface")]
#[derive(Debug, Clone, PartialEq)]
pub struct TrimmedSurface {
    basis: Box<Surface>,
    u1: f64,
    u2: f64,
    v1: f64,
    v2: f64,
}

impl TrimmedSurface {
    /// Panics if a range is empty or outside the basis' domain in a non-periodic direction.
    /// As OCCT's `Geom_RectangularTrimmedSurface`: in a periodic direction the range's end is
    /// shifted by periods to follow its start; in another, a reversed range trims the basis
    /// reversed in that direction.
    pub fn new(basis: Surface, u1: f64, u2: f64, v1: f64, v2: f64) -> Self {
        use super::ParametricSurface;
        let mut basis = basis;
        let (bu1, bu2, bv1, bv2) = basis.bounds();
        let (u1, u2) = Self::range(&mut basis, true, u1, u2, bu1, bu2);
        let (v1, v2) = Self::range(&mut basis, false, v1, v2, bv1, bv2);
        Self { basis: Box::new(basis), u1, u2, v1, v2 }
    }

    /// One direction's range, normalized (reversing the basis if needed).
    fn range(basis: &mut Surface, u: bool, t1: f64, t2: f64, first: f64, last: f64) -> (f64, f64) {
        use super::ParametricSurface;
        assert!((t1 - t2).abs() > crate::tolerance::PARAMETRIC, "an empty trimming range [{t1}, {t2}]");
        let period = if u { basis.u_period() } else { basis.v_period() };
        if let Some(period) = period {
            let mut t2 = t2;
            while t2 <= t1 {
                t2 += period;
            }
            while t2 > t1 + period + crate::tolerance::PARAMETRIC {
                t2 -= period;
            }
            return (t1, t2);
        }
        let outside = |t: f64| t < first - crate::tolerance::PARAMETRIC || t > last + crate::tolerance::PARAMETRIC;
        assert!(!outside(t1) && !outside(t2), "trimming range [{t1}, {t2}] outside [{first}, {last}]");
        if t1 < t2 {
            return (t1, t2);
        }
        if u {
            let r = (basis.u_reversed_parameter(t1), basis.u_reversed_parameter(t2));
            basis.u_reverse();
            r
        } else {
            let r = (basis.v_reversed_parameter(t1), basis.v_reversed_parameter(t2));
            basis.v_reverse();
            r
        }
    }

    pub fn basis(&self) -> &Surface {
        &self.basis
    }

    /// The trimming rectangle `(u1, u2, v1, v2)`.
    pub fn trim(&self) -> (f64, f64, f64, f64) {
        (self.u1, self.u2, self.v1, self.v2)
    }
}

impl super::ParametricSurface for TrimmedSurface {
    fn bounds(&self) -> (f64, f64, f64, f64) {
        (self.u1, self.u2, self.v1, self.v2)
    }
    // As OCCT: a trimmed direction is not closed, but periodicity is the basis'.
    fn is_u_closed(&self) -> bool {
        false
    }
    fn is_v_closed(&self) -> bool {
        false
    }
    fn is_u_periodic(&self) -> bool {
        self.basis.is_u_periodic()
    }
    fn is_v_periodic(&self) -> bool {
        self.basis.is_v_periodic()
    }
    fn u_period(&self) -> Option<f64> {
        self.basis.u_period()
    }
    fn v_period(&self) -> Option<f64> {
        self.basis.v_period()
    }
    fn continuity(&self) -> super::Continuity {
        self.basis.continuity()
    }
    fn is_cn_u(&self, n: u32) -> bool {
        self.basis.is_cn_u(n)
    }
    fn is_cn_v(&self, n: u32) -> bool {
        self.basis.is_cn_v(n)
    }
    fn point(&self, u: f64, v: f64) -> Point {
        self.basis.point(u, v)
    }
    fn d1(&self, u: f64, v: f64) -> (Point, Vector, Vector) {
        self.basis.d1(u, v)
    }
    fn d2(&self, u: f64, v: f64) -> super::SurfaceD2 {
        self.basis.d2(u, v)
    }
    fn d3(&self, u: f64, v: f64) -> super::SurfaceD3 {
        self.basis.d3(u, v)
    }
    fn dn(&self, u: f64, v: f64, nu: u32, nv: u32) -> Vector {
        self.basis.dn(u, v, nu, nv)
    }
    fn u_iso(&self, u: f64) -> Curve {
        Curve::Trimmed(super::curve::TrimmedCurve::new(self.basis.u_iso(u), self.v1, self.v2))
    }
    fn v_iso(&self, v: f64) -> Curve {
        Curve::Trimmed(super::curve::TrimmedCurve::new(self.basis.v_iso(v), self.u1, self.u2))
    }
    fn u_reversed_parameter(&self, u: f64) -> f64 {
        self.basis.u_reversed_parameter(u)
    }
    fn v_reversed_parameter(&self, v: f64) -> f64 {
        self.basis.v_reversed_parameter(v)
    }
    fn u_reverse(&mut self) {
        (self.u1, self.u2) = (self.basis.u_reversed_parameter(self.u2), self.basis.u_reversed_parameter(self.u1));
        self.basis.u_reverse();
    }
    fn v_reverse(&mut self) {
        (self.v1, self.v2) = (self.basis.v_reversed_parameter(self.v2), self.basis.v_reversed_parameter(self.v1));
        self.basis.v_reverse();
    }
    fn resolution(&self, tolerance: f64) -> (f64, f64) {
        self.basis.resolution(tolerance)
    }
}

impl Transformable for TrimmedSurface {
    fn transformed(&self, t: &Transform) -> Self {
        TrimmedSurface { basis: Box::new(self.basis.transformed(t)), ..*self }
    }
}

/// A torus around the Z axis of `frame`: the circle of radius `minor_radius` whose center
/// runs along the circle of radius `major_radius` in the frame's XY plane. Parametrized by
/// `u` (angle around the main axis) and `v` (angle around the tube).
#[doc(alias = "gp_Torus")]
#[doc(alias = "Geom_ToroidalSurface")]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Torus {
    frame: Frame,
    major_radius: f64,
    minor_radius: f64,
}

impl Torus {
    /// Panics if a radius is negative. A minor radius larger than the major one makes a
    /// self-intersecting (spindle) torus.
    pub fn new(frame: Frame, major_radius: f64, minor_radius: f64) -> Self {
        assert!(
            0.0 <= minor_radius && 0.0 <= major_radius,
            "invalid torus radii {major_radius}, {minor_radius}"
        );
        Self { frame, major_radius, minor_radius }
    }

    pub fn frame(&self) -> &Frame {
        &self.frame
    }

    pub fn major_radius(&self) -> f64 {
        self.major_radius
    }

    pub fn minor_radius(&self) -> f64 {
        self.minor_radius
    }

    /// The parameters `(u, v)` of `p`, which should lie on the surface (for other points,
    /// those of a nearby surface point; `u` in `[0, 2π)` for an angular parameter).
    #[doc(alias = "ElSLib::Parameters")]
    pub fn parameters(&self, p: Point) -> (f64, f64) {
        let (x, y, z) = local_coordinates(&self.frame, p);
        let radial = x.hypot(y);
        let u = if radial > 0.0 { angle_in_period(y.atan2(x)) } else { 0.0 };
        (u, angle_in_period(z.atan2(radial - self.major_radius)))
    }


    /// The 35 coefficients of the torus' implicit quartic equation `F(x, y, z) = 0`, in the
    /// order x^4, y^4, z^4, x^3y, x^3z, y^3x, y^3z, z^3x, z^3y, x^2y^2, x^2z^2, y^2z^2,
    /// x^2yz, xy^2z, xyz^2, x^3, y^3, z^3, x^2y, x^2z, y^2x, y^2z, z^2x, z^2y, xyz, x^2, y^2,
    /// z^2, xy, xz, yz, x, y, z, 1.
    pub fn implicit_coefficients(&self) -> [f64; 35] {
        // As OCCT's gp_Torus::Coefficients: the quartic of the torus in its frame, in
        // global coordinates through T, the map into the frame.
        let f = self.frame();
        let (x, y, z, o) = (Vector::from(f.x()), Vector::from(f.y()), Vector::from(f.z()), f.origin() - Point::ORIGIN);
        let t = [[x.x, x.y, x.z, -x.dot(o)], [y.x, y.y, y.z, -y.dot(o)], [z.x, z.y, z.z, -z.dot(o)]];
        let col = |i: usize, j: usize| t[0][i] * t[0][j] + t[1][i] * t[1][j] + t[2][i] * t[2][j];
        let (c11, c22, c33, c44) = (col(0, 0), col(1, 1), col(2, 2), col(3, 3));
        let (c12, c13, c23) = (col(0, 1), col(0, 2), col(1, 2));
        let (c14, c24, c34) = (col(0, 3), col(1, 3), col(2, 3));
        let (r, q) = (self.major_radius(), self.minor_radius());
        let (sum, sub) = (r * r + q * q, r * r - q * q);
        let [t11, t12, t13, t14] = t[0];
        let [t21, t22, t23, t24] = t[1];
        let [t31, t32, t33, t34] = t[2];
        [
            c11 * c11,
            c22 * c22,
            c33 * c33,
            4.0 * c11 * c12,
            4.0 * c11 * c13,
            4.0 * c22 * c12,
            4.0 * c22 * c23,
            4.0 * c33 * c13,
            4.0 * c33 * c23,
            2.0 * (c11 * c22 + 2.0 * c12 * c12),
            2.0 * (c11 * c33 + 2.0 * c13 * c13),
            2.0 * (c22 * c33 + 2.0 * c23 * c23),
            4.0 * (c11 * c23 + 2.0 * c12 * c13),
            4.0 * (c22 * c13 + 2.0 * c12 * c23),
            4.0 * (c33 * c12 + 2.0 * c13 * c23),
            4.0 * c11 * c14,
            4.0 * c22 * c24,
            4.0 * c33 * c34,
            4.0 * (c11 * c24 + 2.0 * c14 * c12),
            4.0 * (c11 * c34 + 2.0 * c14 * c13),
            4.0 * (c22 * c14 + 2.0 * c24 * c12),
            4.0 * (c22 * c34 + 2.0 * c24 * c23),
            4.0 * (c33 * c14 + 2.0 * c34 * c13),
            4.0 * (c33 * c24 + 2.0 * c34 * c23),
            8.0 * (c12 * c34 + c23 * c14 + c24 * c13),
            2.0 * (sub * t31 * t31 - sum * (t11 * t11 + t21 * t21) + c44 * c11 + 2.0 * c14 * c14),
            2.0 * (sub * t32 * t32 - sum * (t12 * t12 + t22 * t22) + c44 * c22 + 2.0 * c24 * c24),
            2.0 * (sub * t33 * t33 - sum * (t13 * t13 + t23 * t23) + c44 * c33 + 2.0 * c34 * c34),
            4.0 * (sub * t31 * t32 - sum * (t11 * t12 + t21 * t22) + c44 * c12 + 2.0 * c14 * c24),
            4.0 * (sub * t31 * t33 - sum * (t11 * t13 + t21 * t23) + c44 * c13 + 2.0 * c14 * c34),
            4.0 * (sub * t32 * t33 - sum * (t12 * t13 + t22 * t23) + c44 * c23 + 2.0 * c24 * c34),
            4.0 * (c44 * c14 + sub * t31 * t34 - sum * (t11 * t14 + t21 * t24)),
            4.0 * (c44 * c24 + sub * t32 * t34 - sum * (t12 * t14 + t22 * t24)),
            4.0 * (c44 * c34 + sub * t33 * t34 - sum * (t13 * t14 + t23 * t24)),
            2.0 * sub * t34 * t34 - 2.0 * sum * (t14 * t14 + t24 * t24) + c44 * c44 + sub * sub,
        ]
    }
}

impl Transformable for Torus {
    fn transformed(&self, t: &Transform) -> Self {
        let scale = t.scale_factor().abs();
        Torus::new(self.frame.transformed(t), self.major_radius * scale, self.minor_radius * scale)
    }
}

// The `ParametricSurface` impl of `Torus` is in eval.rs.

/// The surface swept by a curve translated along a direction:
/// `(u, v) -> basis(u) + v·direction`.
#[doc(alias = "Geom_SurfaceOfLinearExtrusion")]
#[derive(Debug, Clone, PartialEq)]
pub struct ExtrusionSurface {
    basis: Box<Curve>,
    direction: Direction,
    /// The points of a bounded basis that [`Self::parameters`] scans (made on first use).
    scan: super::patches::Shared<Vec<Point>>,
}

impl ExtrusionSurface {
    pub fn new(basis: Curve, direction: Direction) -> Self {
        Self { basis: Box::new(basis), direction, scan: Default::default() }
    }

    pub fn basis(&self) -> &Curve {
        &self.basis
    }

    pub fn direction(&self) -> Direction {
        self.direction
    }
}

// The `ParametricSurface` impl of `ExtrusionSurface` is in eval.rs.

impl ExtrusionSurface {
    /// The parameters `(u, v)` of `p` (on or near the surface): `u` that of the basis
    /// curve's point closest to `p` seen along the direction, `v` the distance from that
    /// point to `p` along the direction.
    #[doc(alias = "ElSLib::Parameters")]
    pub fn parameters(&self, p: Point) -> (f64, f64) {
        use super::ParametricCurve;
        let d = Vector::from(self.direction);
        let flat = |v: Vector| v - d * v.dot(d);
        let gap = |u: f64| flat(self.basis.point(u) - p).length_squared();
        let (mut a, mut b) = (self.basis.first_parameter(), self.basis.last_parameter());
        let bounded = !(a <= -crate::tolerance::INFINITE || b >= crate::tolerance::INFINITE);
        if !bounded {
            // An unbounded basis: around the parameter of its point nearest to p's shadow.
            let t0 = match &*self.basis {
                super::curve::Curve::Line(l) => (p - l.origin).dot(Vector::from(l.direction)),
                _ => 0.0,
            };
            let reach = (p - self.basis.point(t0)).length() * 4.0 + 1.0;
            (a, b) = (t0 - reach, t0 + reach);
        }
        // A coarse scan, then Newton steps on the derivative of the squared gap.
        let n = 64;
        let step = (b - a) / n as f64;
        // (A bounded basis's points scanned are the same each time: kept.)
        let scan = bounded.then(|| self.scan.get(|| Some((0..=n).map(|i| self.basis.point(a + step * i as f64)).collect()))).flatten();
        let (i, scanned) = match scan {
            Some(points) => points.iter().enumerate().map(|(i, &q)| (i, flat(q - p).length_squared())).min_by(|x, y| x.1.total_cmp(&y.1)),
            None => (0..=n).map(|i| (i, gap(a + step * i as f64))).min_by(|x, y| x.1.total_cmp(&y.1)),
        }
        .expect("samples");
        let start = a + step * i as f64;
        let mut u = start;
        let (lo, hi) = (u - step, u + step);
        for _ in 0..30 {
            let (c, d1, d2) = self.basis.d2(u);
            let (w, w1, w2) = (flat(c - p), flat(d1), flat(d2));
            let (g1, g2) = (w.dot(w1), w1.dot(w1) + w.dot(w2));
            if g2 <= 0.0 {
                break;
            }
            let next = (u - g1 / g2).clamp(lo, hi);
            if (next - u).abs() <= 1e-15 * (1.0 + u.abs()) {
                u = next;
                break;
            }
            u = next;
        }
        if gap(u) > scanned {
            u = start;
        }
        if let Some(period) = self.basis.period() {
            u = a + (u - a).rem_euclid(period);
        }
        (u, (p - self.basis.point(u)).dot(d))
    }
}

impl Transformable for ExtrusionSurface {
    fn transformed(&self, t: &Transform) -> Self {
        ExtrusionSurface::new(self.basis.transformed(t), self.direction.transformed(t))
    }
}

/// The surface swept by a curve rotated about an axis: `(u, v) -> basis(v)` rotated by
/// angle `u` (u in [0, 2π), v along the basis).
#[doc(alias = "Geom_SurfaceOfRevolution")]
#[derive(Debug, Clone, PartialEq)]
pub struct RevolutionSurface {
    basis: Box<Curve>,
    axis: Axis,
    /// The meridian points (radius, height) that [`Self::parameters`] scans over the
    /// basis's whole range (made on first use).
    scan: super::patches::Shared<Vec<(f64, f64)>>,
}

impl RevolutionSurface {
    pub fn new(basis: Curve, axis: Axis) -> Self {
        Self { basis: Box::new(basis), axis, scan: Default::default() }
    }

    pub fn basis(&self) -> &Curve {
        &self.basis
    }

    pub fn axis(&self) -> Axis {
        self.axis
    }

    /// The parameters `(u, v)` of `p` (on or near the surface): `v` that of the basis
    /// curve's point closest to `p` in the meridian half-plane (radius from the axis against
    /// height along it), `u` in `[0, 2π)` the angle from that point to `p` about the axis.
    #[doc(alias = "ElSLib::Parameters")]
    pub fn parameters(&self, p: Point) -> (f64, f64) {
        use super::ParametricCurve;
        let o = self.axis.origin;
        let z = Vector::from(self.axis.direction);
        let meridian = |q: Point| {
            let d = q - o;
            let h = d.dot(z);
            ((d - z * h).length(), h)
        };
        let (rp, hp) = meridian(p);
        // A circle in a plane through the axis (a torus's meridian): in closed form, the
        // point's angle on the circle seen in its meridian half-plane.
        if let super::curve::Curve::Circle(circle) = &*self.basis {
            let f = circle.frame();
            let (x, y) = (Vector::from(f.x()), Vector::from(f.y()));
            let normal = x.cross(y);
            let center = circle.center() - o;
            let in_meridian = normal.dot(z).abs() <= 1e-12 && normal.dot(center).abs() <= 1e-12 * (1.0 + center.length());
            let radial_c = center - z * center.dot(z);
            if in_meridian && radial_c.length() > 1e-12 {
                // p turned back into the circle's half-plane.
                let e = radial_c.normalized();
                let q = o + e * rp + z * hp;
                let d = q - circle.center();
                let v = d.dot(y).atan2(d.dot(x));
                let (first, last) = (self.basis.first_parameter(), self.basis.last_parameter());
                let v = if v < first - 1e-12 && last - first >= std::f64::consts::TAU - 1e-12 { v + std::f64::consts::TAU } else { v };
                let v = if v < first - 1e-12 { v + std::f64::consts::TAU } else { v };
                let c = self.basis.point(v) - o;
                let (rc, rq) = (c - z * c.dot(z), (p - o) - z * (p - o).dot(z));
                let u = z.dot(rc.cross(rq)).atan2(rc.dot(rq));
                return (angle_in_period(u), v);
            }
        }
        self.parameters_between(p, self.basis.first_parameter(), self.basis.last_parameter())
    }

    /// [`parameters`](Self::parameters) with `v` sought within `[first, last]` of the basis
    /// only (a face's range: where the basis crosses the axis, the point's nearest point may
    /// lie on its half across the axis, with the normal the other way).
    pub fn parameters_between(&self, p: Point, first: f64, last: f64) -> (f64, f64) {
        use super::ParametricCurve;
        let o = self.axis.origin;
        let z = Vector::from(self.axis.direction);
        let meridian = |q: Point| {
            let d = q - o;
            let h = d.dot(z);
            ((d - z * h).length(), h)
        };
        let (rp, hp) = meridian(p);
        let gap = |v: f64| {
            let (r, h) = meridian(self.basis.point(v));
            (r - rp).powi(2) + (h - hp).powi(2)
        };
        let (mut a, mut b) = (first, last);
        if a <= -crate::tolerance::INFINITE || b >= crate::tolerance::INFINITE {
            // An unbounded basis (a line): around the parameter closest to p.
            let t0 = match &*self.basis {
                super::curve::Curve::Line(l) => (p - l.origin).dot(Vector::from(l.direction)),
                _ => 0.0,
            };
            let reach = (p - o).length() * 4.0 + 1.0;
            (a, b) = (t0 - reach, t0 + reach);
        }
        // A coarse scan, then Newton steps on the derivative of the squared gap.
        let n = 48;
        let step = (b - a) / n as f64;
        // (Over the basis's whole range, the points scanned are the same each time: kept.)
        let whole = a == self.basis.first_parameter() && b == self.basis.last_parameter();
        let scan = whole.then(|| self.scan.get(|| Some((0..=n).map(|i| meridian(self.basis.point(a + step * i as f64))).collect()))).flatten();
        let (i, scanned) = match scan {
            Some(points) => points.iter().enumerate().map(|(i, &(r, h))| (i, (r - rp).powi(2) + (h - hp).powi(2))).min_by(|x, y| x.1.total_cmp(&y.1)),
            None => (0..=n).map(|i| (i, gap(a + step * i as f64))).min_by(|x, y| x.1.total_cmp(&y.1)),
        }
        .expect("samples");
        let mut v = a + step * i as f64;
        let (lo, hi) = (v - step, v + step);
        let mut converged = false;
        for _ in 0..30 {
            let (c, d1, d2) = self.basis.d2(v);
            let d = c - o;
            let h = d.dot(z);
            let radial = d - z * h;
            let r = radial.length().max(1e-300);
            let (h1, h2) = (d1.dot(z), d2.dot(z));
            let (rad1, rad2) = (d1 - z * h1, d2 - z * h2);
            let r1 = radial.dot(rad1) / r;
            let r2 = (rad1.dot(rad1) + radial.dot(rad2) - r1 * r1) / r;
            // With g = (r - rp)² + (h - hp)²: g' / 2 and g'' / 2.
            let g1 = (r - rp) * r1 + (h - hp) * h1;
            let g2 = r1 * r1 + (r - rp) * r2 + h1 * h1 + (h - hp) * h2;
            if g2 <= 0.0 {
                break;
            }
            let next = (v - g1 / g2).clamp(lo, hi);
            if (next - v).abs() <= 1e-15 * (1.0 + v.abs()) {
                v = next;
                converged = true;
                break;
            }
            v = next;
        }
        // Newton may stall (at the axis, where the meridian's radius has a kink): keep the
        // best of it, the scan and a golden-section search around the scan's best.
        let scan_best = a + step * i as f64;
        if gap(v) > scanned {
            v = scan_best;
        }
        let (mut x, mut y) = ((scan_best - step).max(a.min(b)), (scan_best + step).min(a.max(b)));
        let ratio = 0.618_033_988_749_895;
        for _ in 0..if converged { 0 } else { 80 } {
            let (p1, p2) = (y - ratio * (y - x), x + ratio * (y - x));
            if gap(p1) < gap(p2) {
                y = p2;
            } else {
                x = p1;
            }
        }
        let golden = 0.5 * (x + y);
        if gap(golden) < gap(v) {
            v = golden;
        }
        let c = self.basis.point(v) - o;
        let (rc, rq) = (c - z * c.dot(z), (p - o) - z * (p - o).dot(z));
        let u = z.dot(rc.cross(rq)).atan2(rc.dot(rq));
        (angle_in_period(u), v)
    }
}

// The `ParametricSurface` impl of `RevolutionSurface` is in eval.rs.

impl Transformable for RevolutionSurface {
    fn transformed(&self, t: &Transform) -> Self {
        RevolutionSurface::new(self.basis.transformed(t), self.axis.transformed(t))
    }
}

/// The surface at signed distance `offset` from a basis surface, along its normal.
#[doc(alias = "Geom_OffsetSurface")]
#[derive(Debug, Clone, PartialEq)]
pub struct OffsetSurface {
    basis: Box<Surface>,
    offset: f64,
    /// The equivalent non-offset surface, if there is one (made on first use), which
    /// evaluation then goes through.
    equivalent: super::patches::Shared<Surface>,
}

impl OffsetSurface {
    /// An offset of an offset surface is made on its basis, the offsets added (as OCCT's
    /// Geom_OffsetSurface does).
    pub fn new(basis: Surface, offset: f64) -> Self {
        match basis {
            Surface::Offset(inner) => Self { basis: inner.basis, offset: inner.offset + offset, equivalent: Default::default() },
            basis => Self { basis: Box::new(basis), offset, equivalent: Default::default() },
        }
    }

    /// The equivalent surface, kept.
    fn equivalent(&self) -> Option<&Surface> {
        self.equivalent.get(|| self.equivalent_surface())
    }

    pub fn basis(&self) -> &Surface {
        &self.basis
    }

    pub fn offset(&self) -> f64 {
        self.offset
    }

    pub fn set_offset(&mut self, offset: f64) {
        self.offset = offset;
        self.equivalent = Default::default();
    }

    /// The same surface as a non-offset one where that exists: a parallel plane, a
    /// cylinder, cone, sphere or torus with adjusted radius (flipped if the offset goes
    /// through the axis), the basis for a zero offset, and the same for a trimmed basis.
    /// `None` otherwise, or when the result would degenerate (e.g. zero radius).
    #[doc(alias = "Geom_OffsetSurface::Surface")]
    pub fn equivalent_surface(&self) -> Option<Surface> {
        use super::ParametricSurface;
        let d = self.offset;
        let basis = &*self.basis;
        if d == 0.0 {
            return Some(basis.clone());
        }
        // Whether the normal ∂u × ∂v points away from the axis at a sample point (frames
        // may be left-handed).
        let normal_at = |u: f64, v: f64| {
            let (p, du, dv) = basis.d1(u, v);
            (p, du.cross(dv))
        };
        let away = |frame: &Frame, p: Point, n: Vector| {
            let z = Vector::from(frame.z());
            let axis_point = frame.origin() + z * (p - frame.origin()).dot(z);
            n.dot(p - axis_point) > 0.0
        };
        // A radius that would vanish gives no surface.
        let usable = |r: f64| (r.abs() > crate::tolerance::LINEAR).then_some(r);
        Some(match basis {
            Surface::Plane(plane) => {
                let shift = Vector::from(plane.normal()) * d;
                Surface::Plane(Plane::from_frame(plane.frame().transformed(&Transform::translation(shift))))
            }
            Surface::Cylinder(c) => {
                let (p, n) = normal_at(0.0, 0.0);
                let r = usable(c.radius() + if away(c.frame(), p, n) { d } else { -d })?;
                // Through the axis: the cylinder on the other side, its frame turned round so
                // that the parameters still match.
                let mut frame = *c.frame();
                if r < 0.0 {
                    frame.flip_x();
                    frame.flip_y();
                }
                Surface::Cylinder(Cylinder::new(frame, r.abs()))
            }
            Surface::Sphere(s) => {
                let (p, n) = normal_at(0.0, 0.0);
                let r = usable(s.radius() + if n.dot(p - s.frame().origin()) > 0.0 { d } else { -d })?;
                let mut frame = *s.frame();
                if r < 0.0 {
                    frame.flip_x();
                    frame.flip_y();
                    frame.flip_z();
                }
                Surface::Sphere(Sphere::new(frame, r.abs()))
            }
            Surface::Torus(t) => {
                // At the outer equator (elsewhere the tube may reach the axis).
                let (p, n) = normal_at(0.0, 0.0);
                let f = t.frame();
                let centre = f.origin() + Vector::from(f.x()) * t.major_radius();
                let minor = usable(t.minor_radius() + if n.dot(p - centre) > 0.0 { d } else { -d })?;
                Surface::Torus(Torus::new(*f, t.major_radius(), minor.abs()))
            }
            Surface::Cone(c) => {
                // Points move cos α radially and -sin α along the axis (for an outward
                // normal): the cone with its reference radius grown by d·cos α, moved back by
                // d·sin α.
                let (p, n) = normal_at(0.0, 0.0);
                let d = if away(c.frame(), p, n) { d } else { -d };
                let (sin, cos) = c.semi_angle_sin_cos();
                let frame = c.frame().transformed(&Transform::translation(Vector::from(c.frame().z()) * (-d * sin)));
                Surface::Cone(Cone::new(frame, c.semi_angle(), c.radius() + d * cos))
            }
            Surface::Trimmed(t) => {
                let (u1, u2, v1, v2) = t.bounds();
                let offset = OffsetSurface::new(t.basis().clone(), d).equivalent_surface()?;
                Surface::Trimmed(TrimmedSurface::new(offset, u1, u2, v1, v2))
            }
            // The normal of an extrusion is the basis's tangent across the direction: the
            // extrusion of the basis offset across it.
            Surface::Extrusion(e) => {
                let direction = e.direction();
                let candidate = Surface::Extrusion(ExtrusionSurface::new(offset_across(e.basis(), d, direction)?, direction));
                return self.reproduced_by(candidate);
            }
            // The normal of a surface of revolution lies in the meridian plane: the revolution
            // of the meridian offset within its plane.
            Surface::Revolution(r) => {
                use super::ParametricCurve;
                let axis = r.axis();
                let (a, o) = (Vector::from(axis.direction), axis.origin);
                let meridian = r.basis();
                // (A window of unbounded parameters where points stay finite: hyperbolas grow
                // exponentially.)
                let (lo, hi) = (meridian.first_parameter().max(-10.0), meridian.last_parameter().min(10.0));
                let samples: Vec<Point> = (0..=16).map(|i| meridian.point(lo + (hi - lo) * i as f64 / 16.0)).collect();
                let plane_normal = samples.iter().filter(|p| (**p - Point::ORIGIN).length().is_finite()).find_map(|&p| {
                    let radial = (p - o) - a * (p - o).dot(a);
                    (radial.length() > crate::tolerance::LINEAR).then(|| Direction::try_from_vector(a.cross(radial))).flatten()
                })?;
                // ∂u × ∂v = (a × r) × C' = -(C' × m) |a × r|.
                let offset_meridian = offset_across(meridian, -d, plane_normal)?;
                return self.reproduced_by(Surface::Revolution(RevolutionSurface::new(offset_meridian, axis)));
            }
            _ => return None,
        })
    }

    /// `candidate` if it matches this offset surface at sample points (a check for
    /// equivalents derived under assumptions, such as a planar meridian).
    fn reproduced_by(&self, candidate: Surface) -> Option<Surface> {
        use super::ParametricSurface;
        let (u1, u2, v1, v2) = self.basis.bounds();
        let clip = |x: f64| x.clamp(-10.0, 10.0);
        let (u1, u2, v1, v2) = (clip(u1), clip(u2), clip(v1), clip(v2));
        let scale = self.offset.abs().max(1.0);
        for i in 0..=4 {
            for j in 0..=4 {
                let (u, v) = (u1 + (u2 - u1) * (0.1 + 0.2 * i as f64), v1 + (v2 - v1) * (0.1 + 0.2 * j as f64));
                let b = self.basis.d1(u, v);
                let n = b.1.cross(b.2);
                if n.length() <= 1e-300 {
                    continue;
                }
                let expected = b.0 + n * (self.offset / n.length());
                if !(candidate.point(u, v).distance(expected) <= 1e-9 * scale * (expected - Point::ORIGIN).length().max(1.0)) {
                    return None;
                }
            }
        }
        Some(candidate)
    }

    /// An iso-curve without closed form: a cubic B-spline through points of it, at their
    /// own parameters (finite parts of infinite ranges), with the points doubled until it
    /// stays within 1e-10 (relative) of the iso-curve between them.
    fn iso(&self, at: &dyn Fn(f64) -> Point, lo: f64, hi: f64) -> Curve {
        let (lo, hi) = (lo.max(-1e3), hi.min(1e3));
        let mut n = 64;
        loop {
            let params: Vec<f64> = (0..=n).map(|i| lo + (hi - lo) * i as f64 / n as f64).collect();
            let points: Vec<Point> = params.iter().map(|&t| at(t)).collect();
            let Ok(curve) = super::fit::interpolate_with_parameters(&points, &params) else {
                return Curve::Line(super::curve::Line::new(points[0], Direction::X));
            };
            let size = points.iter().fold(1.0f64, |m, p| m.max((*p - Point::ORIGIN).length()));
            let error = params.windows(2).map(|w| {
                let t = 0.5 * (w[0] + w[1]);
                use super::ParametricCurve;
                curve.point(t).distance(at(t))
            }).fold(0.0, f64::max);
            if error <= 1e-10 * size || n >= 4096 {
                return Curve::BSpline(curve);
            }
            n *= 2;
        }
    }

    /// The unit normal `∂u × ∂v` of the basis and its derivatives at `(u, v)` (zero where the
    /// basis is singular).
    fn normal_d1(&self, u: f64, v: f64) -> (super::SurfaceD2, Vector, Vector, Vector) {
        use super::ParametricSurface;
        let b = self.basis.d2(u, v);
        let n = b.du.cross(b.dv);
        let length = n.length();
        if length <= 1e-300 {
            return (b, Vector::ZERO, Vector::ZERO, Vector::ZERO);
        }
        let unit = n * (1.0 / length);
        // d(n/|n|) = (dn - unit (unit · dn)) / |n|
        let dn_u = b.duu.cross(b.dv) + b.du.cross(b.duv);
        let dn_v = b.duv.cross(b.dv) + b.du.cross(b.dvv);
        let unit_u = (dn_u - unit * unit.dot(dn_u)) * (1.0 / length);
        let unit_v = (dn_v - unit * unit.dot(dn_v)) * (1.0 / length);
        (b, unit, unit_u, unit_v)
    }
}

/// `curve` moved by `distance` along `curve' × direction` (normalized): a line or a circle
/// about `direction` stays one, others become offset curves.
fn offset_across(curve: &Curve, distance: f64, direction: Direction) -> Option<Curve> {
    use super::curve::{Circle, Line, OffsetCurve};
    let z = Vector::from(direction);
    Some(match curve {
        Curve::Line(l) => {
            let n = Direction::try_from_vector(Vector::from(l.direction).cross(z))?;
            Curve::Line(Line::new(l.origin + Vector::from(n) * distance, l.direction))
        }
        Curve::Circle(c) if Vector::from(c.frame().z()).cross(z).length() <= crate::tolerance::ANGULAR => {
            // At t = 0 the tangent is the frame's y: outwards if y × direction is x.
            let f = c.frame();
            let outwards = Vector::from(f.y()).cross(z).dot(Vector::from(f.x())) > 0.0;
            let r = c.radius() + if outwards { distance } else { -distance };
            if r.abs() <= crate::tolerance::LINEAR {
                return None;
            }
            let mut frame = *f;
            if r < 0.0 {
                // Through the centre: the opposite point at each angle.
                frame.flip_x();
                frame.flip_y();
            }
            Curve::Circle(Circle::new(frame, r.abs()))
        }
        Curve::Trimmed(t) => {
            let (a, b) = t.range();
            Curve::Trimmed(super::curve::TrimmedCurve::new(offset_across(t.basis(), distance, direction)?, a, b))
        }
        other => Curve::Offset(OffsetCurve::new(other.clone(), distance, direction)),
    })
}

impl super::ParametricSurface for OffsetSurface {
    fn bounds(&self) -> (f64, f64, f64, f64) {
        self.basis.bounds()
    }
    fn is_u_closed(&self) -> bool {
        self.basis.is_u_closed()
    }
    fn is_v_closed(&self) -> bool {
        self.basis.is_v_closed()
    }
    fn is_u_periodic(&self) -> bool {
        self.basis.is_u_periodic()
    }
    fn is_v_periodic(&self) -> bool {
        self.basis.is_v_periodic()
    }
    fn u_period(&self) -> Option<f64> {
        self.basis.u_period()
    }
    fn v_period(&self) -> Option<f64> {
        self.basis.v_period()
    }
    fn continuity(&self) -> super::Continuity {
        // One order less than the basis's (the normal takes a derivative).
        match self.basis.continuity() {
            super::Continuity::C0 | super::Continuity::C1 => super::Continuity::C0,
            super::Continuity::G1 | super::Continuity::G2 => super::Continuity::C0,
            super::Continuity::C2 => super::Continuity::C1,
            super::Continuity::C3 => super::Continuity::C2,
            super::Continuity::CN => super::Continuity::CN,
        }
    }
    fn is_cn_u(&self, n: u32) -> bool {
        self.basis.is_cn_u(n + 1)
    }
    fn is_cn_v(&self, n: u32) -> bool {
        self.basis.is_cn_v(n + 1)
    }
    fn point(&self, u: f64, v: f64) -> Point {
        if let Some(e) = self.equivalent() {
            return e.point(u, v);
        }
        let (b, unit, _, _) = self.normal_d1(u, v);
        b.point + unit * self.offset
    }
    fn d1(&self, u: f64, v: f64) -> (Point, Vector, Vector) {
        if let Some(e) = self.equivalent() {
            return e.d1(u, v);
        }
        let (b, unit, unit_u, unit_v) = self.normal_d1(u, v);
        (b.point + unit * self.offset, b.du + unit_u * self.offset, b.dv + unit_v * self.offset)
    }
    fn d2(&self, u: f64, v: f64) -> super::SurfaceD2 {
        if let Some(e) = self.equivalent() {
            return e.d2(u, v);
        }
        // N = n / |n| with n = Bu × Bv: N_x = f_x / L with f_x = n_x - N (N·n_x), L = |n|,
        // and N_xy = (∂y f_x L - f_x L_y) / L².
        let b = self.basis.d3(u, v);
        let n = b.du.cross(b.dv);
        let l = n.length();
        if l <= 1e-300 {
            let (p, du, dv) = self.d1(u, v);
            return super::SurfaceD2 { point: p, du, dv, duu: Vector::ZERO, dvv: Vector::ZERO, duv: Vector::ZERO };
        }
        let unit = n * (1.0 / l);
        let n_u = b.duu.cross(b.dv) + b.du.cross(b.duv);
        let n_v = b.duv.cross(b.dv) + b.du.cross(b.dvv);
        let n_uu = b.duuu.cross(b.dv) + b.duu.cross(b.duv) * 2.0 + b.du.cross(b.duuv);
        let n_uv = b.duuv.cross(b.dv) + b.duu.cross(b.dvv) + b.du.cross(b.duvv);
        let n_vv = b.duvv.cross(b.dv) + b.duv.cross(b.dvv) * 2.0 + b.du.cross(b.dvvv);
        let (l_u, l_v) = (unit.dot(n_u), unit.dot(n_v));
        let (f_u, f_v) = (n_u - unit * l_u, n_v - unit * l_v);
        let (unit_u, unit_v) = (f_u * (1.0 / l), f_v * (1.0 / l));
        let second = |f: Vector, n_x: Vector, l_x: f64, unit_y: Vector, n_xy: Vector, l_y: f64| {
            // ∂y f_x = n_xy - N_y (N·n_x) - N (N_y·n_x + N·n_xy)
            let df = n_xy - unit_y * l_x - unit * (unit_y.dot(n_x) + unit.dot(n_xy));
            (df * l - f * l_y) * (1.0 / (l * l))
        };
        let d = self.offset;
        super::SurfaceD2 {
            point: b.point + unit * d,
            du: b.du + unit_u * d,
            dv: b.dv + unit_v * d,
            duu: b.duu + second(f_u, n_u, l_u, unit_u, n_uu, l_u) * d,
            dvv: b.dvv + second(f_v, n_v, l_v, unit_v, n_vv, l_v) * d,
            duv: b.duv + second(f_u, n_u, l_u, unit_v, n_uv, l_v) * d,
        }
    }
    fn d3(&self, u: f64, v: f64) -> super::SurfaceD3 {
        if let Some(e) = self.equivalent() {
            return e.d3(u, v);
        }
        // Third derivatives from the second by central differences, Richardson-extrapolated.
        let h = 1e-3 * (1.0 + u.abs().max(v.abs()));
        let diff = |f: &dyn Fn(f64) -> super::SurfaceD2, h: f64| {
            let (a, b) = (f(h), f(-h));
            [(a.duu - b.duu) * (0.5 / h), (a.dvv - b.dvv) * (0.5 / h), (a.duv - b.duv) * (0.5 / h)]
        };
        let along_u = |t: f64| self.d2(u + t, v);
        let along_v = |t: f64| self.d2(u, v + t);
        let rich = |f: &dyn Fn(f64) -> super::SurfaceD2| {
            let (coarse, fine) = (diff(f, h), diff(f, 0.5 * h));
            [0, 1, 2].map(|k| (fine[k] * 4.0 - coarse[k]) * (1.0 / 3.0))
        };
        let (du, dv) = (rich(&along_u), rich(&along_v));
        let s = self.d2(u, v);
        super::SurfaceD3 {
            point: s.point,
            du: s.du,
            dv: s.dv,
            duu: s.duu,
            dvv: s.dvv,
            duv: s.duv,
            duuu: du[0],
            dvvv: dv[1],
            duuv: du[2],
            duvv: dv[2],
        }
    }
    fn dn(&self, u: f64, v: f64, nu: u32, nv: u32) -> Vector {
        if let Some(e) = self.equivalent() {
            return e.dn(u, v, nu, nv);
        }
        match (nu, nv) {
            (0, 0) => self.point(u, v) - Point::ORIGIN,
            (1, 0) => self.d1(u, v).1,
            (0, 1) => self.d1(u, v).2,
            (2, 0) => self.d2(u, v).duu,
            (0, 2) => self.d2(u, v).dvv,
            (1, 1) => self.d2(u, v).duv,
            (3, 0) => self.d3(u, v).duuu,
            (0, 3) => self.d3(u, v).dvvv,
            (2, 1) => self.d3(u, v).duuv,
            (1, 2) => self.d3(u, v).duvv,
            // Higher orders by central differences of the order below.
            (nu, nv) => {
                let h = 1e-3 * (1.0 + u.abs().max(v.abs()));
                if nu > 0 {
                    (self.dn(u + h, v, nu - 1, nv) - self.dn(u - h, v, nu - 1, nv)) * (0.5 / h)
                } else {
                    (self.dn(u, v + h, nu, nv - 1) - self.dn(u, v - h, nu, nv - 1)) * (0.5 / h)
                }
            }
        }
    }
    fn u_iso(&self, u: f64) -> Curve {
        if let Some(e) = self.equivalent() {
            return e.u_iso(u);
        }
        let (_, _, v0, v1) = self.bounds();
        self.iso(&|t| self.point(u, t), v0, v1)
    }
    fn v_iso(&self, v: f64) -> Curve {
        if let Some(e) = self.equivalent() {
            return e.v_iso(v);
        }
        let (u0, u1, _, _) = self.bounds();
        self.iso(&|t| self.point(t, v), u0, u1)
    }
    fn u_reversed_parameter(&self, u: f64) -> f64 {
        self.basis.u_reversed_parameter(u)
    }
    fn v_reversed_parameter(&self, v: f64) -> f64 {
        self.basis.v_reversed_parameter(v)
    }
    fn u_reverse(&mut self) {
        // The normal turns round with the parameter: the offset changes side to stay put.
        self.basis.u_reverse();
        self.offset = -self.offset;
        self.equivalent = Default::default();
    }
    fn v_reverse(&mut self) {
        self.basis.v_reverse();
        self.offset = -self.offset;
        self.equivalent = Default::default();
    }
    fn resolution(&self, tolerance: f64) -> (f64, f64) {
        self.basis.resolution(tolerance)
    }
}

impl Transformable for OffsetSurface {
    fn transformed(&self, t: &Transform) -> Self {
        // As OCCT: a negative scale (a reflection) moves the offset to the other side.
        OffsetSurface::new(self.basis.transformed(t), self.offset * t.scale_factor())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::ParametricSurface;

    #[test]
    fn cone_parameters_round_trip_on_both_sides_of_the_apex() {
        let cone = Cone::new(Frame::from_z(Point::new(1.0, 2.0, 3.0), Direction::new(0.2, -0.1, 1.0)), 0.4, 1.5);
        let surface = Surface::Cone(cone);
        for (u, v) in [(0.3, 2.0), (2.5, -1.0), (0.7, -9.0), (4.0, -12.5), (5.9, 0.0)] {
            let p = surface.point(u, v);
            let (pu, pv) = cone.parameters(p);
            assert!(surface.point(pu, pv).distance(p) < 1e-9, "({u}, {v}) came back as ({pu}, {pv})");
            assert!((pv - v).abs() < 1e-9);
        }
    }
}
