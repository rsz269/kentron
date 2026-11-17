//! OCCT's intersection of a curve with an infinite plane, `IntCurveSurface_HInter::Perform`
//! for a plane (`GeomAdaptor_Surface` of a `Geom_Plane`), ported line by line.
//!
//! The dispatch is OCCT's `PerformBounds`: lines, circles, ellipses, parabolas and hyperbolas
//! are solved in closed form (`IntAna_IntConicQuad` against the plane, with
//! `IntAna_QuadQuadGeo` for the plane of a circle and `IntAna2d_AnaIntersection` for the line
//! met there, `math_TrigonometricFunctionRoots` for ellipses and `math_DirectPolynomialRoots`
//! for parabolas and hyperbolas); any other curve goes to
//! `IntCurveSurface_TheQuadCurvExactHInter`, which samples the plane's signed distance along
//! each C1 span of the curve and solves it with `math_FunctionAllRoots`. Each solution is then
//! filtered against the curve's bounds (`ComputeAppendPoint`), given the plane's parameters
//! (`ElSLib::Parameters`) and a transition, and appended unless a point with the same
//! parameters and transition is there already (`IntCurveSurface_Intersection::Append`).
//!
//! The plane is OCCT's whole infinite plane: its parameter bounds are
//! `±Precision::Infinite()`, which `PerformBounds` clamps to `±1.0e50` (`ClampUVParameters`)
//! for the polyhedral paths a plane never takes; the parameter filter tests against the
//! unclamped bounds, as OCCT's does.

use std::f64::consts::PI;

use crate::math::{Direction, Frame, Point, Vector};

use super::adaptor::{CurveAdaptor, CurveKind, epsilon, in_period};
use super::curve::Curve;
use super::function_roots::{FunctionSample, FunctionWithDerivative, NumericError, function_all_roots};
use super::polynomial_roots::{PolynomialRoots, trigonometric_roots};
use super::surface::Plane;
use super::Continuity;

/// `THE_TOLTANGENCY`: the tolerance on the curve's and the plane's parameter bounds.
const TOLTANGENCY: f64 = 0.00000001;
/// `THE_TOLERANCE_ANGULAIRE`.
const TOLERANCE_ANGULAIRE: f64 = 1.0e-12;
/// `THE_TOLERANCE`: the distance tolerance handed to the circle and ellipse cases.
const TOLERANCE: f64 = 0.00000001;
/// `Precision::Infinite()`: the plane's parameter bounds.
const INFINITE: f64 = 2.0e100;
/// `THE_PARAM_LIMIT` of `ClampUVParameters`.
const PARAM_LIMIT: f64 = 1.0e50;
/// `Precision::PConfusion()`.
const P_CONFUSION: f64 = 1.0e-9;
/// `EPSX`, `EPSDIST`, `EPSNUL` of `IntCurveSurface_QuadricCurveExactInterUtils`.
const EPSX: f64 = 0.00000000000001;
const EPSDIST: f64 = 0.00000001;
const EPSNUL: f64 = 0.00000001;
/// `PARAMEQUAL`'s tolerance in `IntCurveSurface_Intersection::Append`.
const PARAM_EQUAL: f64 = 1e-8;

/// How the curve crosses the plane at a point, against the plane's normal.
#[doc(alias = "IntCurveSurface_TransitionOnCurve")]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Transition {
    /// The curve's tangent is (nearly) in the plane.
    Tangent,
    /// The curve goes against the normal.
    In,
    /// The curve goes along the normal.
    Out,
}

/// One intersection point.
#[doc(alias = "IntCurveSurface_IntersectionPoint")]
#[derive(Clone, Copy, Debug)]
pub(crate) struct CurvePlanePoint {
    /// `Pnt()`: the curve's point at `w`.
    pub point: Point,
    /// `U()`, `V()`: the plane's parameters of the point the solver found.
    pub u: f64,
    pub v: f64,
    /// `W()`: the curve's parameter.
    pub w: f64,
    /// `Transition()`.
    pub transition: Transition,
}

/// The result of [`intersect_curve_plane`].
#[derive(Clone, Debug, Default)]
pub(crate) struct CurvePlaneIntersection {
    /// `Point(1..=NbPoints())`, in OCCT's order.
    pub points: Vec<CurvePlanePoint>,
    /// `IsParallel()`: a conic parallel to the plane or in it (no points then).
    pub parallel: bool,
}

/// The plane OCCT's `Geom_Plane(P, V)` makes: `gp_Ax3(P, V)`, its X direction chosen by
/// zeroing the smallest coordinate of `V` and swapping the other two, its Y direction
/// `V × X`.
#[doc(alias = "gp_Ax3")]
pub(crate) fn plane_through(p: Point, v: Direction) -> Plane {
    let (a, b, c) = (v.x(), v.y(), v.z());
    let (aabs, babs, cabs) = (a.abs(), b.abs(), c.abs());
    let d = if babs <= aabs && babs <= cabs {
        if aabs > cabs { Direction::new(-c, 0.0, a) } else { Direction::new(c, 0.0, -a) }
    } else if aabs <= babs && aabs <= cabs {
        if babs > cabs { Direction::new(0.0, -c, b) } else { Direction::new(0.0, c, -b) }
    } else if aabs > babs {
        Direction::new(-b, a, 0.0)
    } else {
        Direction::new(b, -a, 0.0)
    };
    let y = Direction::from_vector(v.cross(d));
    Plane::from_frame(Frame::from_axes_unchecked(p, d, y, v))
}

/// A `gp_Ax3`: location and three directions (left-handed allowed).
#[derive(Clone, Copy, Debug)]
struct Ax3 {
    location: [f64; 3],
    x: [f64; 3],
    y: [f64; 3],
    z: [f64; 3],
}

fn coords(d: Direction) -> [f64; 3] {
    [d.x(), d.y(), d.z()]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}

fn magnitude(a: [f64; 3]) -> f64 {
    (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt()
}

/// `gp_Dir(x, y, z)`: normalized (left as is when its squared length is exactly 1).
fn gp_dir(a: [f64; 3]) -> [f64; 3] {
    let sq = a[0] * a[0] + a[1] * a[1] + a[2] * a[2];
    if (1.0 - f64::MIN_POSITIVE..=1.0 + f64::MIN_POSITIVE).contains(&sq) {
        return a;
    }
    let d = sq.sqrt();
    [a[0] / d, a[1] / d, a[2] / d]
}

impl Ax3 {
    fn from_frame(frame: &Frame) -> Ax3 {
        let o = frame.origin();
        Ax3 { location: [o.x, o.y, o.z], x: coords(frame.x()), y: coords(frame.y()), z: coords(frame.z()) }
    }

    /// `Direct()`.
    fn direct(&self) -> bool {
        dot(cross(self.x, self.y), self.z) > 0.0
    }

    /// `gp_Pln::Coefficients(A, B, C, D)`: the plane `A x + B y + C z + D = 0`, its normal
    /// reversed for a left-handed position.
    fn coefficients(&self) -> (f64, f64, f64, f64) {
        let (a, b, c) = if self.direct() { (self.z[0], self.z[1], self.z[2]) } else { (-self.z[0], -self.z[1], -self.z[2]) };
        let p = self.location;
        let d = -(a * p[0] + b * p[1] + c * p[2]);
        (a, b, c, d)
    }

    /// `gp_Pln::Distance(P)`: the unsigned distance along the main direction.
    fn distance(&self, p: [f64; 3]) -> f64 {
        let loc = self.location;
        let dir = self.z;
        (dir[0] * (p[0] - loc[0]) + dir[1] * (p[1] - loc[1]) + dir[2] * (p[2] - loc[2])).abs()
    }

    /// `gp_Trsf::SetTransformation(A3)`: the rows of the matrix (the axes) and the
    /// translation `-(M · Location)`.
    fn transformation(&self) -> ([[f64; 3]; 3], [f64; 3]) {
        let m = [self.x, self.y, self.z];
        let o = self.location;
        let mut loc = [0.0; 3];
        for i in 0..3 {
            loc[i] = m[i][0] * o[0] + m[i][1] * o[1] + m[i][2] * o[2];
        }
        (m, [-loc[0], -loc[1], -loc[2]])
    }

    /// `ElSLib::PlaneParameters(Pos, P, U, V)`: `P` transformed into the plane's frame.
    fn plane_parameters(&self, p: [f64; 3]) -> (f64, f64) {
        let (m, loc) = self.transformation();
        let u = m[0][0] * p[0] + m[0][1] * p[1] + m[0][2] * p[2];
        let v = m[1][0] * p[0] + m[1][1] * p[1] + m[1][2] * p[2];
        (u + loc[0], v + loc[1])
    }
}

/// `IntAna_Quadric`: `Cxx x² + Cyy y² + Czz z² + 2(Cxy xy + Cxz xz + Cyz yz) + 2(Cx x + Cy y
/// + Cz z) + Ccte`.
#[doc(alias = "IntAna_Quadric")]
#[derive(Clone, Copy, Debug)]
struct Quadric {
    cxx: f64,
    cyy: f64,
    czz: f64,
    cxy: f64,
    cxz: f64,
    cyz: f64,
    cx: f64,
    cy: f64,
    cz: f64,
    ccte: f64,
}

impl Quadric {
    /// `IntAna_Quadric(gp_Pln)`.
    fn plane(p: &Ax3) -> Quadric {
        let (a, b, c, d) = p.coefficients();
        Quadric { cxx: 0.0, cyy: 0.0, czz: 0.0, cxy: 0.0, cxz: 0.0, cyz: 0.0, cx: a * 0.5, cy: b * 0.5, cz: c * 0.5, ccte: d }
    }

    /// `NewCoefficients(..., Axis)`: the coefficients in the frame `axis`.
    fn new_coefficients(&self, axis: &Ax3) -> Quadric {
        // The inverse of `SetTransformation(Axis)`: the transposed matrix, and the
        // translation `-(Mᵀ · loc)`.
        let (m, loc) = axis.transformation();
        let mut inv_loc = [0.0; 3];
        for i in 0..3 {
            inv_loc[i] = m[0][i] * loc[0] + m[1][i] * loc[1] + m[2][i] * loc[2];
        }
        let inv_loc = [-inv_loc[0], -inv_loc[1], -inv_loc[2]];
        // x = t11 X + t12 Y + t13 Z + t14, and so on.
        let (t11, t12, t13, t14) = (m[0][0], m[1][0], m[2][0], inv_loc[0]);
        let (t21, t22, t23, t24) = (m[0][1], m[1][1], m[2][1], inv_loc[1]);
        let (t31, t32, t33, t34) = (m[0][2], m[1][2], m[2][2], inv_loc[2]);
        let Quadric { cxx, cyy, czz, cxy, cxz, cyz, cx, cy, cz, ccte } = *self;

        let t11_p2 = t11 * t11;
        let t21_p2 = t21 * t21;
        let t31_p2 = t31 * t31;
        let t12_p2 = t12 * t12;
        let t22_p2 = t22 * t22;
        let t32_p2 = t32 * t32;
        let t13_p2 = t13 * t13;
        let t23_p2 = t23 * t23;
        let t33_p2 = t33 * t33;
        let t14_p2 = t14 * t14;
        let t24_p2 = t24 * t24;
        let t34_p2 = t34 * t34;

        let new_ccte = ccte + t14_p2 * cxx + t24_p2 * cyy + t34_p2 * czz + 2.0 * (t14 * (cx + t24 * cxy + t34 * cxz) + t24 * (cy + t34 * cyz) + t34 * cz);
        let new_cxx = t11_p2 * cxx + t21_p2 * cyy + t31_p2 * czz + 2.0 * (t11 * (t21 * cxy + t31 * cxz) + t21 * t31 * cyz);
        let new_cyy = t12_p2 * cxx + t22_p2 * cyy + t32_p2 * czz + 2.0 * (t12 * (t22 * cxy + t32 * cxz) + t22 * t32 * cyz);
        let new_czz = t13_p2 * cxx + t33_p2 * czz + t23_p2 * cyy + 2.0 * (t13 * (t23 * cxy + t33 * cxz) + t23 * t33 * cyz);
        let new_cz = t13 * cx + t13 * (t14 * cxx + t24 * cxy + t34 * cxz) + t14 * (t23 * cxy + t33 * cxz) + t23 * (cy + t24 * cyy + t34 * cyz) + t33 * (t24 * cyz + cz + t34 * czz);
        let new_cx = t11 * (cx + t14 * cxx + t24 * cxy + t34 * cxz) + t14 * (t21 * cxy + t31 * cxz) + t21 * (cy + t24 * cyy + t34 * cyz) + t31 * (t24 * cyz + cz + t34 * czz);
        let new_cxy = t11 * (t12 * cxx + t22 * cxy + t32 * cxz) + t12 * (t21 * cxy + t31 * cxz) + t21 * (t22 * cyy + t32 * cyz) + t31 * (t22 * cyz + t32 * czz);
        let new_cxz = t11 * (t13 * cxx + t23 * cxy + t33 * cxz) + t13 * (t21 * cxy + t31 * cxz) + t21 * (t23 * cyy + t33 * cyz) + t31 * (t23 * cyz + t33 * czz);
        let new_cy = t12 * (cx + t14 * cxx + t24 * cxy + t34 * cxz) + t14 * (t22 * cxy + t32 * cxz) + t22 * (cy + t24 * cyy + t34 * cyz) + t32 * (cz + t24 * cyz + t34 * czz);
        let new_cyz = t12 * (t13 * cxx + t23 * cxy + t33 * cxz) + t13 * (t22 * cxy + t32 * cxz) + t22 * (t23 * cyy + t33 * cyz) + t32 * (t23 * cyz + t33 * czz);
        Quadric { cxx: new_cxx, cyy: new_cyy, czz: new_czz, cxy: new_cxy, cxz: new_cxz, cyz: new_cyz, cx: new_cx, cy: new_cy, cz: new_cz, ccte: new_ccte }
    }
}

/// The result of `IntAna_IntConicQuad`.
#[doc(alias = "IntAna_IntConicQuad")]
#[derive(Clone, Debug, Default)]
struct ConicQuad {
    done: bool,
    parallel: bool,
    inquadric: bool,
    /// `Point(i)`, `ParamOnConic(i)`.
    points: Vec<([f64; 3], f64)>,
}

/// `IntAna_IntConicQuad::Perform(gp_Lin, gp_Pln, Tolang, Tol, Len)`.
fn line_plane(origin: [f64; 3], direction: [f64; 3], plane: &Ax3, tolang: f64, tol: f64, len: f64) -> ConicQuad {
    // Tolang is the angular tolerance under which the angle between two vectors counts as
    // zero (reasoning on its cosine, cos(t) close to t near Pi/2).
    let mut out = ConicQuad::default();
    let (a, b, c, d) = plane.coefficients();
    let (al, bl, cl) = (direction[0], direction[1], direction[2]);
    let direc = a * al + b * bl + c * cl;
    let dis = a * origin[0] + b * origin[1] + c * origin[2] + d;
    out.parallel = false;
    if direc.abs() < tolang {
        out.parallel = true;
        if len != 0.0 && direc != 0.0 {
            // The distance from the bounding point of the line to the plane.
            let p1 = [origin[0] - dis * a, origin[1] - dis * b, origin[2] - dis * c];
            let p2 = [p1[0] + len * al, p1[1] + len * bl, p1[2] + len * cl];
            if plane.distance(p2) > tol {
                out.parallel = false;
            }
        }
    }
    if out.parallel {
        out.inquadric = dis.abs() < tolang;
    } else {
        out.inquadric = false;
        let t = -dis / direc;
        out.points.push(([origin[0] + t * al, origin[1] + t * bl, origin[2] + t * cl], t));
    }
    out.done = true;
    out
}

/// The result of `IntAna_QuadQuadGeo` for two planes.
enum PlanePlane {
    NotDone,
    Empty,
    Same,
    /// `Line(1)`: location and direction.
    Line([f64; 3], [f64; 3]),
}

/// `IntAna_QuadQuadGeo(gp_Pln P1, gp_Pln P2, TolAng, Tol)`.
fn plane_plane(p1: &Ax3, p2: &Ax3, tol_ang: f64, tol: f64) -> PlanePlane {
    let (a1, b1, c1, d1) = p1.coefficients();
    let (a2, b2, c2, d2) = p2.coefficients();
    let vn1 = [a1, b1, c1];
    let vn2 = [a2, b2, c2];
    let vd = cross(vn1, vn2);
    let loc1 = p1.location;
    let loc2 = p2.location;
    let dist1 = a2 * loc1[0] + b2 * loc1[1] + c2 * loc1[2] + d2;
    let dist2 = a1 * loc2[0] + b1 * loc2[1] + c1 * loc2[2] + d1;
    let mvd = magnitude(vd);
    if mvd <= tol_ang {
        // The normals are collinear: the planes are the same or parallel.
        return if dist1.abs() <= tol && dist2.abs() <= tol { PlanePlane::Same } else { PlanePlane::Empty };
    }
    let eps = 1.0e-16;
    let mut denom = a1 * a2 + b1 * b2 + c1 * c2;
    let denom2 = denom * denom;
    let ddenom = 1.0 - denom2;
    denom = if ddenom.abs() <= eps { eps } else { ddenom };
    let par1 = dist1 / denom;
    let par2 = -dist2 / denom;
    let inter1 = cross(vn1, vd);
    let inter2 = cross(vn2, vd);
    let x1 = loc1[0] + par1 * inter1[0];
    let y1 = loc1[1] + par1 * inter1[1];
    let z1 = loc1[2] + par1 * inter1[2];
    let x2 = loc2[0] + par2 * inter2[0];
    let y2 = loc2[1] + par2 * inter2[1];
    let z2 = loc2[2] + par2 * inter2[2];
    let mut pt1 = [(x1 + x2) * 0.5, (y1 + y2) * 0.5, (z1 + z2) * 0.5];
    let dir1 = gp_dir(vd);
    // When the angle between the planes is small, the origin of the line is computed with
    // an error that cannot be considered small: it is refined where possible.
    let tresh_ang = 2.0e-6;
    let tresh_dist = 1.0e-12;
    if mvd < tresh_ang {
        let dist_1 = a1 * pt1[0] + b1 * pt1[1] + c1 * pt1[2] + d1;
        let dist_2 = a2 * pt1[0] + b2 * pt1[1] + c2 * pt1[2] + d2;
        if dist_1.abs() > tresh_dist || dist_2.abs() > tresh_dist {
            // 1.
            let dn1 = gp_dir(vn1);
            let icq = line_plane(pt1, dn1, p1, tol_ang, tol, 0.0);
            if !icq.done {
                return PlanePlane::NotDone;
            }
            // (`Point(1)` throws on a parallel line: a line along the plane's normal is not.)
            let pnt1 = icq.points[0].0;
            // 2.
            let dl2 = gp_dir(cross(dir1, dn1));
            let icq = line_plane(pnt1, dl2, p2, tol_ang, tol, 0.0);
            if !icq.done || icq.parallel {
                return PlanePlane::NotDone;
            }
            pt1 = icq.points[0].0;
        }
    }
    PlanePlane::Line(pt1, dir1)
}

/// `gp_Dir2d::Angle(Other)` (with OCCT's literal for cos 45°).
#[allow(clippy::approx_constant)]
fn dir2d_angle(a: [f64; 2], b: [f64; 2]) -> f64 {
    let cosinus = a[0] * b[0] + a[1] * b[1];
    let sinus = a[0] * b[1] - a[1] * b[0];
    if cosinus > -0.70710678118655 && cosinus < 0.70710678118655 {
        if sinus > 0.0 { cosinus.acos() } else { -cosinus.acos() }
    } else if cosinus > 0.0 {
        sinus.asin()
    } else if sinus > 0.0 {
        PI - sinus.asin()
    } else {
        -PI - sinus.asin()
    }
}

/// `gp_Dir2d(x, y)`.
fn gp_dir2d(x: f64, y: f64) -> [f64; 2] {
    let sq = x * x + y * y;
    if (1.0 - f64::MIN_POSITIVE..=1.0 + f64::MIN_POSITIVE).contains(&sq) {
        return [x, y];
    }
    let d = sq.sqrt();
    [x / d, y / d]
}

/// `ElCLib::CircleParameter` for the circle at the origin with the standard axes:
/// the angle of `p`, normalized into `[0, 2π]` (`normalizeAngle`).
fn circle_parameter(p: [f64; 2]) -> f64 {
    let mut teta = dir2d_angle([1.0, 0.0], gp_dir2d(p[0] - 0.0, p[1] - 0.0));
    // (X ^ Y = 1: a direct circle keeps the angle.)
    let pipi = PI + PI;
    while teta < -f64::EPSILON {
        teta += pipi;
    }
    while teta > pipi * (1.0 + f64::MIN_POSITIVE) {
        teta -= pipi;
    }
    if teta < 0.0 {
        teta = 0.0;
    }
    teta
}

/// `IntAna2d_AnaIntersection(gp_Lin2d L, gp_Circ2d C)` for a circle of radius `radius` at
/// the origin with the standard axes: `(x, y, parameter on the line, angle on the circle)`.
fn line2d_circle(origin: [f64; 2], direction: [f64; 2], radius: f64) -> Vec<(f64, f64, f64, f64)> {
    // gp_Lin2d::Coefficients.
    let a = direction[1];
    let b = -direction[0];
    let c0 = -(a * origin[0] + b * origin[1]);
    let d = a * 0.0 + b * 0.0 + c0;
    let line_parameter = |x: f64, y: f64| (x - origin[0]) * direction[0] + (y - origin[1]) * direction[1];
    if d.abs() - radius > epsilon(radius) {
        return Vec::new();
    }
    if (d.abs() - radius).abs() <= epsilon(radius) {
        // The tangent case.
        let xs = 0.0 - d * a;
        let ys = 0.0 - d * b;
        return vec![(xs, ys, line_parameter(xs, ys), circle_parameter([xs, ys]))];
    }
    // Two intersection points.
    let h = (radius * radius - d * d).sqrt();
    let xs1 = 0.0 - d * a - h * b;
    let ys1 = 0.0 - d * b + h * a;
    let xs2 = 0.0 - d * a + h * b;
    let ys2 = 0.0 - d * b - h * a;
    vec![
        (xs1, ys1, line_parameter(xs1, ys1), circle_parameter([xs1, ys1])),
        (xs2, ys2, line_parameter(xs2, ys2), circle_parameter([xs2, ys2])),
    ]
}

/// `IntAna_IntConicQuad::Perform(gp_Circ, gp_Pln, Tolang, Tol)`.
fn circle_plane(position: &Ax3, radius: f64, plane: &Ax3, tolang: f64, tol: f64) -> ConicQuad {
    let mut out = ConicQuad::default();
    let plconic = *position;
    match plane_plane(&plconic, plane, tolang, tol) {
        PlanePlane::NotDone => out,
        PlanePlane::Empty => {
            out.parallel = true;
            let distmax = plane.distance(position.location) + radius * tolang;
            out.inquadric = distmax < tol;
            out.done = true;
            out
        }
        PlanePlane::Same => {
            out.inquadric = true;
            out.done = true;
            out
        }
        PlanePlane::Line(location, direction) => {
            let o = plconic.location;
            let v0 = [location[0] - o[0], location[1] - o[1], location[2] - o[2]];
            let axex = plconic.x;
            let axey = plconic.y;
            let orig = [dot(axex, v0), dot(axey, v0)];
            let dire = gp_dir2d(dot(axex, direction), dot(axey, direction));
            for (x, y, _, ang) in line2d_circle(orig, dire, radius) {
                let p = [o[0] + x * axex[0] + y * axey[0], o[1] + x * axex[1] + y * axey[1], o[2] + x * axex[2] + y * axey[2]];
                out.points.push((p, ang));
            }
            out.done = true;
            out
        }
    }
}

/// `ElCLib::EllipseValue` (and `CircleValue` with two equal radii).
fn ellipse_value(u: f64, pos: &Ax3, major: f64, minor: f64) -> [f64; 3] {
    let a1 = major * u.cos();
    let a2 = minor * u.sin();
    let (x, y, l) = (pos.x, pos.y, pos.location);
    [a1 * x[0] + a2 * y[0] + l[0], a1 * x[1] + a2 * y[1] + l[1], a1 * x[2] + a2 * y[2] + l[2]]
}

/// `ElCLib::HyperbolaValue`.
fn hyperbola_value(u: f64, pos: &Ax3, major: f64, minor: f64) -> [f64; 3] {
    let a1 = major * u.cosh();
    let a2 = minor * u.sinh();
    let (x, y, l) = (pos.x, pos.y, pos.location);
    [a1 * x[0] + a2 * y[0] + l[0], a1 * x[1] + a2 * y[1] + l[1], a1 * x[2] + a2 * y[2] + l[2]]
}

/// `ElCLib::ParabolaValue`.
fn parabola_value(u: f64, pos: &Ax3, focal: f64) -> [f64; 3] {
    let (x, y, l) = (pos.x, pos.y, pos.location);
    if focal.abs() <= f64::MIN_POSITIVE {
        return [u * x[0] + l[0], u * x[1] + l[1], u * x[2] + l[2]];
    }
    let a1 = u * u / (4.0 * focal);
    [a1 * x[0] + u * y[0] + l[0], a1 * x[1] + u * y[1] + l[1], a1 * x[2] + u * y[2] + l[2]]
}

/// `IntAna_IntConicQuad::Perform(gp_Elips, IntAna_Quadric)`: in the ellipse's frame,
/// `x = R cos t`, `y = r sin t`, `z = 0` turn the quadric into a trigonometric equation.
fn ellipse_quadric(pos: &Ax3, major: f64, minor: f64, quad: &Quadric) -> ConicQuad {
    let mut out = ConicQuad::default();
    let q = quad.new_coefficients(pos);
    let (r_big, r) = (major, minor);
    let p_cos_cos = r_big * r_big * q.cxx;
    let p_sin_sin = r * r * q.cyy;
    let p_sin = r * q.cy;
    let p_cos = r_big * q.cx;
    let p_cos_sin = r_big * r * q.cxy;
    let p_cte = q.ccte;
    let pol = trigonometric_roots(p_cos_cos - p_sin_sin, p_cos_sin, p_cos + p_cos, p_sin + p_sin, p_cte + p_sin_sin, 0.0, PI + PI);
    if pol.done {
        out.done = true;
        if pol.infinite {
            out.inquadric = true;
        } else {
            for &t in &pol.roots {
                out.points.push((ellipse_value(t, pos, r_big, r), t));
            }
        }
    }
    out
}

/// `IntAna_IntConicQuad::Perform(gp_Parab, IntAna_Quadric)`: in the parabola's frame,
/// `x = y² / 2p`, which makes a quartic in `y`.
fn parabola_quadric(pos: &Ax3, focal: f64, quad: &Quadric) -> ConicQuad {
    let mut out = ConicQuad::default();
    let q = quad.new_coefficients(pos);
    let f = focal;
    let un_sur_2p = 0.25 / f;
    let a4 = q.cxx * un_sur_2p * un_sur_2p;
    let a3 = (q.cxy + q.cxy) * un_sur_2p;
    let a2 = q.cyy + (q.cx + q.cx) * un_sur_2p;
    let a1 = q.cy + q.cy;
    let a0 = q.ccte;
    let pol = PolynomialRoots::quartic(a4, a3, a2, a1, a0);
    if pol.done {
        out.done = true;
        if pol.infinite {
            out.inquadric = true;
        } else {
            for &t in &pol.roots {
                out.points.push((parabola_value(t, pos, f), t));
            }
        }
    }
    out
}

/// `IntAna_IntConicQuad::Perform(gp_Hypr, IntAna_Quadric)`: in the hyperbola's frame,
/// `x = R ch t`, `y = r sh t`, which makes a quartic in `exp t`.
fn hyperbola_quadric(pos: &Ax3, major: f64, minor: f64, quad: &Quadric) -> ConicQuad {
    let mut out = ConicQuad::default();
    let q = quad.new_coefficients(pos);
    let (r_big, r) = (major, minor);
    let rr_big = r_big * r_big;
    let rr = r * r;
    let r_r = r_big * r;
    let a4 = rr_big * q.cxx + r_r * (q.cxy + q.cxy) + rr * q.cyy;
    let a3 = 4.0 * (r_big * q.cx + r * q.cy);
    let a2 = 2.0 * ((q.ccte + q.ccte) + q.cxx * rr_big - q.cyy * rr);
    let a1 = 4.0 * (r_big * q.cx - r * q.cy);
    let a0 = q.cxx * rr_big - r_r * (q.cxy + q.cxy) + q.cyy * rr;
    let pol = PolynomialRoots::quartic(a4, a3, a2, a1, a0);
    if pol.done {
        out.done = true;
        if pol.infinite {
            out.inquadric = true;
        } else {
            for &t in &pol.roots {
                if t >= f64::EPSILON {
                    let lnt = t.ln();
                    out.points.push((hyperbola_value(lnt, pos, r_big, r), lnt));
                }
            }
        }
    }
    out
}

/// `IntSurf_Quadric` of a plane: its signed distance `A x + B y + C z + D` and gradient.
#[doc(alias = "IntSurf_Quadric")]
#[derive(Clone, Copy, Debug)]
struct PlaneQuadric {
    prm1: f64,
    prm2: f64,
    prm3: f64,
    prm4: f64,
}

impl PlaneQuadric {
    fn distance(&self, p: Point) -> f64 {
        self.prm1 * p.x + self.prm2 * p.y + self.prm3 * p.z + self.prm4
    }

    fn gradient(&self) -> [f64; 3] {
        [self.prm1, self.prm2, self.prm3]
    }
}

/// `IntCurveSurface_TheQuadCurvFuncOfTheQuadCurvExactHInter`: the quadric's distance along
/// the curve, with its derivative.
#[doc(alias = "IntCurveSurface_TheQuadCurvFuncOfTheQuadCurvExactHInter")]
struct QuadCurvFunction<'a> {
    quadric: PlaneQuadric,
    curve: &'a dyn CurveAdaptor,
}

impl FunctionWithDerivative for QuadCurvFunction<'_> {
    fn value(&mut self, param: f64, f: &mut f64) -> bool {
        *f = self.quadric.distance(self.curve.value(param));
        true
    }

    fn derivative(&mut self, param: f64, d: &mut f64) -> bool {
        let (_, t) = self.curve.d1(param);
        let g = self.quadric.gradient();
        *d = t.x * g[0] + t.y * g[1] + t.z * g[2];
        true
    }

    fn values(&mut self, param: f64, f: &mut f64, d: &mut f64) -> bool {
        let (p, t) = self.curve.d1(param);
        // (ValAndGrad.)
        *f = self.quadric.distance(p);
        let g = self.quadric.gradient();
        *d = t.x * g[0] + t.y * g[1] + t.z * g[2];
        true
    }

    fn get_state_number(&mut self) -> i32 {
        0
    }
}

/// `IntCurveSurface_TheHCurveTool::NbSamples(C, U0, U1)`.
#[doc(alias = "IntCurveSurface_TheHCurveTool::NbSamples")]
fn nb_samples(curve: &dyn CurveAdaptor, u0: f64, u1: f64) -> i32 {
    let nbs_other = 10.0;
    let mut nbs: f64 = nbs_other;
    match (curve.curve_type(), curve.basis_curve()) {
        (Some(CurveKind::Line), _) => nbs = 2.0,
        (Some(CurveKind::Bezier), Some(Curve::Bezier(b))) => nbs = 3.0 + b.poles().len() as f64,
        (Some(CurveKind::BSpline), Some(Curve::BSpline(b))) => {
            nbs = b.knots().len() as f64;
            nbs *= b.degree() as f64;
            nbs *= curve.last_parameter() - curve.first_parameter();
            let range = u1 - u0;
            if range.abs() > P_CONFUSION {
                nbs /= range;
            }
            if nbs < 2.0 {
                nbs = 2.0;
            }
        }
        _ => {}
    }
    if nbs > 50.0 {
        nbs = 50.0;
    }
    nbs as i32
}

/// `IntCurveSurface_TheQuadCurvExactHInter(S, C)` for a plane: the roots of the plane's
/// distance along each C1 span of the curve (`PerformIntersection`). The null intervals
/// it also finds are not reported (OCCT's `IntCurveSurface_HInter` reads only the roots).
#[doc(alias = "IntCurveSurface_TheQuadCurvExactHInter")]
fn quad_curv_exact(plane: &Ax3, curve: &dyn CurveAdaptor) -> Result<Vec<f64>, NumericError> {
    let (prm1, prm2, prm3, prm4) = plane.coefficients();
    let quadric = PlaneQuadric { prm1, prm2, prm3, prm4 };
    let intervals = curve.intervals(Continuity::C1);
    let mut pnts = Vec::new();
    for span in intervals.windows(2) {
        let (u1, u2) = (span[0], span[1]);
        let sample = FunctionSample { a: u1, b: u2, n: nb_samples(curve, u1, u2) };
        let mut function = QuadCurvFunction { quadric, curve };
        // (`IsDone()` is always true when no `Standard_NumericError` is raised.)
        let roots = function_all_roots(&mut function, &sample, EPSX, EPSDIST, EPSNUL)?;
        pnts.extend_from_slice(&roots.points);
    }
    Ok(pnts)
}

/// `IntCurveSurface_HInter::Perform(curve, surface)` for an infinite plane: the points where
/// `curve` crosses or touches `plane`, in OCCT's order. Errs where OCCT raises
/// `Standard_NumericError` from `math_FunctionAllRoots` (in debug builds; release builds read
/// past an empty result there).
#[doc(alias = "IntCurveSurface_HInter")]
pub(crate) fn intersect_curve_plane(curve: &dyn CurveAdaptor, plane: &Plane) -> Result<CurvePlaneIntersection, NumericError> {
    let mut out = CurvePlaneIntersection::default();
    let pln = Ax3::from_frame(plane.frame());
    // `DecomposeSurfaceIntervals`: a plane is one C2 interval over its whole (infinite)
    // parameter range, which `PerformBounds` clamps against overflow; the plane's paths do
    // not use the clamped bounds.
    let (mut uu1, mut uu2, mut vv1, mut vv2) = (-INFINITE, INFINITE, -INFINITE, INFINITE);
    clamp_uv_parameters(&mut uu1, &mut uu2, &mut vv1, &mut vv2);

    let conic = |kind: CurveKind| curve.basis_curve().unwrap_or_else(|| panic!("a {kind:?} adaptor without its curve"));
    match curve.curve_type() {
        Some(CurveKind::Line) => {
            let Curve::Line(line) = conic(CurveKind::Line) else { unreachable!("a line adaptor's curve is a line") };
            let o = line.origin;
            let ana = line_plane([o.x, o.y, o.z], coords(line.direction), &pln, TOLERANCE_ANGULAIRE, 0.0, 0.0);
            append_int_ana(curve, &pln, &ana, &mut out);
        }
        Some(CurveKind::Circle) => {
            let Curve::Circle(circle) = conic(CurveKind::Circle) else { unreachable!("a circle adaptor's curve is a circle") };
            let ana = circle_plane(&Ax3::from_frame(circle.frame()), circle.radius(), &pln, TOLERANCE_ANGULAIRE, TOLERANCE);
            append_int_ana(curve, &pln, &ana, &mut out);
        }
        Some(CurveKind::Ellipse) => {
            let Curve::Ellipse(ellipse) = conic(CurveKind::Ellipse) else { unreachable!("an ellipse adaptor's curve is an ellipse") };
            let ana = ellipse_quadric(&Ax3::from_frame(ellipse.frame()), ellipse.major_radius(), ellipse.minor_radius(), &Quadric::plane(&pln));
            append_int_ana(curve, &pln, &ana, &mut out);
        }
        Some(CurveKind::Parabola) => {
            let Curve::Parabola(parabola) = conic(CurveKind::Parabola) else { unreachable!("a parabola adaptor's curve is a parabola") };
            let ana = parabola_quadric(&Ax3::from_frame(parabola.frame()), parabola.focal(), &Quadric::plane(&pln));
            append_int_ana(curve, &pln, &ana, &mut out);
        }
        Some(CurveKind::Hyperbola) => {
            let Curve::Hyperbola(hyperbola) = conic(CurveKind::Hyperbola) else { unreachable!("a hyperbola adaptor's curve is a hyperbola") };
            let ana = hyperbola_quadric(&Ax3::from_frame(hyperbola.frame()), hyperbola.major_radius(), hyperbola.minor_radius(), &Quadric::plane(&pln));
            append_int_ana(curve, &pln, &ana, &mut out);
        }
        _ => {
            // `InternalPerformCurveQuadric` (a plane is a quadric).
            for w in quad_curv_exact(&pln, curve)? {
                let p = curve.value(w);
                let (u, v) = pln.plane_parameters([p.x, p.y, p.z]);
                if let Some(point) = compute_append_point(curve, w, &pln, u, v) {
                    append(&mut out.points, point);
                }
            }
        }
    }
    Ok(out)
}

/// `ClampUVParameters`: protection from double overflow (bug26525).
fn clamp_uv_parameters(u1: &mut f64, u2: &mut f64, v1: &mut f64, v2: &mut f64) {
    if *u1 < -PARAM_LIMIT {
        *u1 = -PARAM_LIMIT;
    }
    if *u2 > PARAM_LIMIT {
        *u2 = PARAM_LIMIT;
    }
    if *v1 < -PARAM_LIMIT {
        *v1 = -PARAM_LIMIT;
    }
    if *v2 > PARAM_LIMIT {
        *v2 = PARAM_LIMIT;
    }
}

/// `AppendIntAna` (`ProcessIntAna`): a conic parallel to the plane or in it only sets
/// `IsParallel`; otherwise each solution is filtered and appended.
fn append_int_ana(curve: &dyn CurveAdaptor, plane: &Ax3, ana: &ConicQuad, out: &mut CurvePlaneIntersection) {
    if !ana.done {
        return;
    }
    if ana.inquadric || ana.parallel {
        out.parallel = true;
        return;
    }
    let mut points = Vec::new();
    for &(p, w) in &ana.points {
        let (u, v) = plane.plane_parameters(p);
        if let Some(point) = compute_append_point(curve, w, plane, u, v) {
            points.push(point);
        }
    }
    for point in points {
        append(&mut out.points, point);
    }
}

/// `ComputeAppendPoint`: the solution brought into the curve's period, rejected outside the
/// curve's or the plane's bounds, given its transition and the curve's point.
fn compute_append_point(curve: &dyn CurveAdaptor, lw: f64, plane: &Ax3, su: f64, sv: f64) -> Option<CurvePlanePoint> {
    let w0 = curve.first_parameter();
    let w1 = curve.last_parameter();
    // (The plane's bounds: `GeomAdaptor_Surface` of an infinite `Geom_Plane`.)
    let (u0, u1, v0, v1) = (-INFINITE, INFINITE, -INFINITE, INFINITE);
    let mut w = lw;
    let (u, v) = (su, sv);
    let kind = curve.curve_type();
    if curve.is_periodic() || kind == Some(CurveKind::Circle) || kind == Some(CurveKind::Ellipse) {
        w = in_period(w, w0, w0 + curve.period());
    }
    if w0 - w >= TOLTANGENCY || w - w1 >= TOLTANGENCY {
        return None;
    }
    // (A plane is neither periodic nor one of the surfaces whose U is an angle.)
    if u0 - u >= TOLTANGENCY || u - u1 >= TOLTANGENCY {
        return None;
    }
    if v0 - v >= TOLTANGENCY || v - v1 >= TOLTANGENCY {
        return None;
    }
    let transition = compute_transition(curve, w, plane);
    let p = curve.value(w);
    Some(CurvePlanePoint { point: p, u, v, w, transition })
}

/// `ComputeTransitions`: the curve's unit tangent against the plane's normal `D1U ^ D1V`.
fn compute_transition(curve: &dyn CurveAdaptor, w: f64, plane: &Ax3) -> Transition {
    // (ElSLib::PlaneD1: the derivatives are the plane's X and Y directions.)
    let nsurf = cross(plane.x, plane.y);
    let (_, d1u) = curve.d1(w);
    let d1u: Vector = d1u;
    let norm = magnitude(nsurf);
    if norm > TOLERANCE_ANGULAIRE && d1u.length_squared() > TOLERANCE_ANGULAIRE {
        let d = magnitude([d1u.x, d1u.y, d1u.z]);
        let t = [d1u.x / d, d1u.y / d, d1u.z / d];
        let mut cos_dir = dot(nsurf, t);
        cos_dir /= norm;
        if -cos_dir > TOLERANCE_ANGULAIRE {
            // --Curve--->    <----Surface----
            Transition::In
        } else if cos_dir > TOLERANCE_ANGULAIRE {
            // --Curve--->  ----Surface-->
            Transition::Out
        } else {
            Transition::Tangent
        }
    } else {
        Transition::Tangent
    }
}

/// `IntCurveSurface_Intersection::Append(Point)`: skipped when a point with the same
/// parameters (within 1e-8) and transition is there already.
fn append(points: &mut Vec<CurvePlanePoint>, point: CurvePlanePoint) {
    let equal = |a: f64, b: f64| (a - b).abs() < PARAM_EQUAL;
    for other in points.iter() {
        if equal(point.u, other.u) && equal(point.v, other.v) && equal(point.w, other.w) && other.transition == point.transition {
            return;
        }
    }
    points.push(point);
}

#[cfg(test)]
pub(crate) mod tests;
