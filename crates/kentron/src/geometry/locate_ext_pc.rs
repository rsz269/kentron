//! OCCT's local point-to-curve extremum (`Extrema_LocateExtPC`, that is
//! `Extrema_GLocateExtPC` over `Adaptor3d_Curve`), ported: the extremum of the distance
//! from a point to a curve nearest a starting parameter. B-spline, offset and other curves
//! go by Newton (`Extrema_GenLocateExtPC` on `Extrema_GFuncExtPC`) on the `C2` interval of
//! the starting parameter, then on the intervals either side; Bézier curves by Newton over
//! the whole range; lines and conics take the global extrema (`Extrema_GGExtPC` with
//! `Extrema_ExtPElC`) and keep the one whose parameter is nearest the start.

use std::f64::consts::PI;

use crate::math::{Direction, Point, Vector};

use super::ParametricCurve;
use super::Continuity;
use super::adaptor::{CurveAdaptor, CurveKind, adjust_periodic, in_period};
use super::curve::{Circle, Curve, Ellipse, Hyperbola, Line, Parabola};
use super::ext_pc::{DistanceFunction, ExtPCPoint};
use super::function_root::function_root;
use super::function_roots::FunctionWithDerivative;
use super::polynomial_roots::{PolynomialRoots, trigonometric_roots};

/// `Precision::Confusion()`.
const CONFUSION: f64 = 1.0e-7;
/// `Precision::SquareConfusion()`.
const SQUARE_CONFUSION: f64 = CONFUSION * CONFUSION;
/// `Precision::Angular()`.
const ANGULAR: f64 = 1.0e-12;
/// `Precision::Infinite()`.
const INFINITE: f64 = 2.0e100;
/// `RealLast()`.
const REAL_LAST: f64 = f64::MAX;
/// `RealEpsilon()`.
const REAL_EPSILON: f64 = f64::EPSILON;
/// `gp::Resolution()`.
const RESOLUTION: f64 = f64::MIN_POSITIVE;

/// `Precision::IsInfinite`.
fn is_infinite(value: f64) -> bool {
    value.abs() >= 0.5 * INFINITE
}

/// `std::max(a, b)`: `b` only when `a < b`.
fn std_max(a: f64, b: f64) -> f64 {
    if a < b { b } else { a }
}

/// `std::min(a, b)`: `b` only when `b < a`.
fn std_min(a: f64, b: f64) -> f64 {
    if b < a { b } else { a }
}

/// The extremum found (`Point()`, as `Extrema_POnCurv`, with `SquareDistance()` and
/// `IsMin()`).
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct LocatedExtremum {
    pub parameter: f64,
    pub point: Point,
    pub square_distance: f64,
    pub is_min: bool,
}

impl From<ExtPCPoint> for LocatedExtremum {
    fn from(extremum: ExtPCPoint) -> Self {
        LocatedExtremum { parameter: extremum.parameter, point: extremum.point, square_distance: extremum.square_distance, is_min: extremum.is_min }
    }
}

/// `Extrema_LocateExtPC(P, C, U0, Umin, Usup, TolF)`: the extremum of the distance from `p`
/// to `curve` nearest the parameter `u0`, searched within `[umin, usup]`; `None` where
/// `IsDone()` is false. `tolf` stops the Newton iterations of B-spline, offset and other
/// curves (`|Un - Un-1| < TolF`); Bézier curves, lines and conics use the curve's
/// `Resolution(Precision::Confusion())` instead, as OCCT.
///
/// The kind of curve is `curve.curve_type()` (`GetType()`; `None` is `GeomAbs_OtherCurve`),
/// a line or conic read from `curve.basis_curve()`.
///
/// For OCCT's `Extrema_LocateExtPC(P, C, U0, TolF)` (the whole curve), pass the curve's
/// first and last parameters.
#[doc(alias = "Extrema_LocateExtPC")]
#[doc(alias = "Extrema_GLocateExtPC")]
pub(crate) fn locate_ext_pc(p: Point, curve: &dyn CurveAdaptor, u0: f64, umin: f64, usup: f64, tolf: f64) -> Option<LocatedExtremum> {
    // `Initialize`.
    let tolu = curve.resolution(CONFUSION);
    match curve.curve_type() {
        None | Some(CurveKind::BSpline) | Some(CurveKind::Offset) => locate_by_intervals(p, curve, u0, umin, usup, tolf),
        Some(CurveKind::Bezier) => {
            // (`myLocExtPC` as `Initialize` left it: over `[Umin, Usup]` to `tolu`, from `U0`
            // unclamped.)
            let mut function = DistanceFunction::new(curve);
            gen_locate_ext_pc(&mut function, p, u0, umin, usup, tolu).map(LocatedExtremum::from)
        }
        Some(kind) => {
            let basis = curve.basis_curve().expect("GeomAdaptor_Curve of a line or conic");
            // `myExtremPC.Initialize(C, Umin, Usup, tolu)` (`TolF` is `tolu`, unused by
            // lines and conics), `Perform(P)`.
            let extrema = ext_pc_elementary(p, curve, kind, basis, umin, usup)?;
            // The extremum whose parameter is nearest `U0` (the last of equals).
            let mut nearest = None;
            let mut nearest_gap = REAL_LAST;
            for extremum in extrema {
                let gap = (extremum.parameter - u0).abs();
                if gap <= nearest_gap {
                    nearest_gap = gap;
                    nearest = Some(extremum);
                }
            }
            nearest.map(LocatedExtremum::from)
        }
    }
}

/// `Perform` for a B-spline, offset or other curve: Newton on the `C2` interval holding
/// `u0` (clamped into `[umin, usup]`), else interval by interval outwards from it, taking
/// an interval bound where `(C(u) - P).C'(u)` changes sign across it.
fn locate_by_intervals(p: Point, curve: &dyn CurveAdaptor, u0: f64, umin: f64, usup: f64, tol: f64) -> Option<LocatedExtremum> {
    let mut function = DistanceFunction::new(curve);
    // Search for extremum is done interval by continuous C2 interval.
    let inter_bounds = curve.intervals(Continuity::C2);
    let n = inter_bounds.len() as i32 - 1;
    // `theInter(i)`, 1-based.
    let the_inter = |i: i32| inter_bounds[i as usize - 1];
    // Be gentle with the caller.
    let mut local_u0 = u0;
    if local_u0 < umin {
        local_u0 = umin;
    } else if local_u0 > usup {
        local_u0 = usup;
    }
    // Search for interval containing U0.
    let mut found = false;
    let mut inter = 1;
    let mut int_uinf = 0.0;
    let mut int_usup = 0.0;
    while !found && inter <= n {
        int_uinf = std_max(the_inter(inter), umin);
        int_usup = std_min(the_inter(inter + 1), usup);
        if local_u0 >= int_uinf && local_u0 < int_usup {
            found = true;
        }
        inter += 1;
    }
    if found {
        // IFV 16.06.00 - inter is increased after found!
        inter -= 1;
    }

    // Try on found interval.
    if let Some(extremum) = gen_locate_ext_pc(&mut function, p, local_u0, int_uinf, int_usup, tol) {
        return Some(extremum.into());
    }

    // Try on neighboring intervals.
    let mut k = 1;
    let mut i1 = inter;
    let mut i2 = inter;
    let (p1, v1) = curve.d1(int_uinf);
    let mut s2inf = (p1 - p).dot(v1);
    let (p1, v1) = curve.d1(int_usup);
    let mut s1sup = (p1 - p).dot(v1);
    while i2 > 0 && i1 <= n {
        i1 = inter + k;
        i2 = inter - k;
        if i1 <= n {
            int_uinf = std_max(the_inter(i1), umin);
            int_usup = std_min(the_inter(i1 + 1), usup);
            if int_uinf < int_usup {
                let (p1, v1) = curve.d1(int_uinf);
                let s2sup = (p1 - p).dot(v1);
                if is_infinite(s2sup) || is_infinite(s1sup) {
                    break;
                }
                if s1sup * s2sup <= REAL_EPSILON {
                    // Extremum.
                    return Some(LocatedExtremum { parameter: int_uinf, point: p1, square_distance: p.distance_squared(p1), is_min: s1sup <= 0.0 });
                }
                let (p1, v1) = curve.d1(int_usup);
                s1sup = (p1 - p).dot(v1);
                if let Some(extremum) = gen_locate_ext_pc(&mut function, p, (int_uinf + int_usup) * 0.5, int_uinf, int_usup, tol) {
                    return Some(extremum.into());
                }
            }
        }

        if i2 > 0 {
            int_uinf = std_max(the_inter(i2), umin);
            int_usup = std_min(the_inter(i2 + 1), usup);
            if int_uinf < int_usup {
                let (p1, v1) = curve.d1(int_usup);
                let s1inf = (p1 - p).dot(v1);
                if is_infinite(s2inf) || is_infinite(s1inf) {
                    break;
                }
                if s1inf * s2inf <= REAL_EPSILON {
                    // Extremum.
                    return Some(LocatedExtremum { parameter: int_usup, point: p1, square_distance: p.distance_squared(p1), is_min: s1inf <= 0.0 });
                }
                let (p1, v1) = curve.d1(int_uinf);
                s2inf = (p1 - p).dot(v1);
                if let Some(extremum) = gen_locate_ext_pc(&mut function, p, (int_uinf + int_usup) * 0.5, int_uinf, int_usup, tol) {
                    return Some(extremum.into());
                }
            }
        }

        k += 1;
    }
    None
}

/// `Extrema_GenLocateExtPC` (`Extrema_LocEPCOfLocateExtPC`): `Initialize(C, Umin, Usup,
/// TolU)` then `Perform(P, U0)`. The root of the distance function by `math_FunctionRoot`
/// from `u0` within `[umin, usup]`, rejected where the function is not below `1e-7` there.
/// (`Initialize(C)` of the function only depends on the curve, so one function serves every
/// call.)
#[doc(alias = "Extrema_GenLocateExtPC")]
#[doc(alias = "Extrema_LocEPCOfLocateExtPC")]
fn gen_locate_ext_pc(function: &mut DistanceFunction<'_>, p: Point, u0: f64, umin: f64, usup: f64, tolu: f64) -> Option<ExtPCPoint> {
    function.set_point(p);
    let root = function_root(function, u0, tolu, umin, usup, 100);
    if !root.done {
        return None;
    }
    // `Point()`: `myF.Point(1)`, saved by `math_FunctionRoot` (`GetStateNumber`).
    let extremum = function.extrema[0];
    let mut value = 0.0;
    if !function.value(extremum.parameter, &mut value) || value.abs() >= 1.0e-7 {
        return None;
    }
    Some(extremum)
}

/// `Extrema_CurveTool::IsPeriodic`: circles and ellipses always, else the adaptor's.
fn curve_tool_is_periodic(curve: &dyn CurveAdaptor, kind: CurveKind) -> bool {
    matches!(kind, CurveKind::Circle | CurveKind::Ellipse) || curve.is_periodic()
}

/// `Extrema_ELPCOfLocateExtPC` (`Extrema_GGExtPC`) `Initialize(C, Uinf, Usup)`, `Perform(P)`
/// for a line or conic: the extrema of `Extrema_ExtPElC` brought into the period of a
/// periodic curve and kept within `mytolu` of `[uinf, usup]`; `None` where not done. (The
/// end points' distances `Perform` computes are only read for other kinds of curve.)
#[doc(alias = "Extrema_ELPCOfLocateExtPC")]
#[doc(alias = "Extrema_GGExtPC")]
fn ext_pc_elementary(p: Point, curve: &dyn CurveAdaptor, kind: CurveKind, basis: &Curve, uinf: f64, usup: f64) -> Option<Vec<ExtPCPoint>> {
    let tolu = curve.resolution(CONFUSION);
    let t3d = CONFUSION;
    let elementary = match (kind, basis) {
        (CurveKind::Circle, Curve::Circle(c)) => ext_p_elc_circle(p, c, t3d, uinf, usup),
        (CurveKind::Ellipse, Curve::Ellipse(c)) => ext_p_elc_ellipse(p, c, t3d, uinf, usup),
        (CurveKind::Parabola, Curve::Parabola(c)) => ext_p_elc_parabola(p, c, uinf, usup),
        (CurveKind::Hyperbola, Curve::Hyperbola(c)) => ext_p_elc_hyperbola(p, c, t3d, uinf, usup),
        (CurveKind::Line, Curve::Line(l)) => ext_p_elc_line(p, l, t3d, uinf, usup),
        _ => unreachable!("GetType() names the curve"),
    }?;
    let is_periodic = curve_tool_is_periodic(curve, kind);
    let mut extrema = Vec::new();
    for extremum in elementary {
        let mut u = extremum.parameter;
        if is_periodic {
            u = in_period(u, uinf, uinf + curve.period());
        }
        if u >= uinf - tolu && u <= usup + tolu {
            extrema.push(ExtPCPoint { parameter: u, ..extremum });
        }
    }
    Some(extrema)
}

/// `P` projected into the plane of a conic (`Pp`).
fn projected(p: Point, origin: Point, axis: Vector) -> Point {
    let translation = axis * -((p - origin).dot(axis));
    p + translation
}

/// `Extrema_ExtPElC::Perform(P, gp_Lin, Tol, Uinf, Usup)`.
#[doc(alias = "Extrema_ExtPElC")]
fn ext_p_elc_line(p: Point, line: &Line, tol: f64, uinf: f64, usup: f64) -> Option<Vec<ExtPCPoint>> {
    let v1 = Vector::from(line.direction);
    let origin = line.origin;
    let distance = v1.dot(p - origin);
    if distance >= uinf - tol && distance <= usup + tol {
        let foot = origin + v1 * distance;
        Some(vec![ExtPCPoint { parameter: distance, point: foot, square_distance: p.distance_squared(foot), is_min: true }])
    } else {
        // (OCCT leaves the line's extrema not done.)
        None
    }
}

/// `gp_Dir::AngleWithRef`: the angle from `a` to `b` (unit vectors), in `[-π, π]`,
/// positive where `a × b` points the way `reference` does.
#[allow(clippy::approx_constant)]
fn angle_with_ref(a: Vector, b: Vector, reference: Vector) -> f64 {
    let cross = a.cross(b);
    let cos = a.dot(b);
    let sin = cross.length();
    // (OCCT's constant, not quite 1/√2.)
    let angle = if cos > -0.70710678118655 && cos < 0.70710678118655 {
        cos.acos()
    } else if cos < 0.0 {
        PI - sin.asin()
    } else {
        sin.asin()
    };
    if cross.dot(reference) >= 0.0 { angle } else { -angle }
}

/// `Extrema_ExtPElC::Perform(P, gp_Circ, Tol, Uinf, Usup)`: the minimum at the angle of `P`
/// projected into the circle's plane, the maximum opposite; not done where that projection
/// is within `tol` of the centre.
#[doc(alias = "Extrema_ExtPElC")]
fn ext_p_elc_circle(p: Point, circle: &Circle, tol: f64, uinf: f64, usup: f64) -> Option<Vec<ExtPCPoint>> {
    let frame = circle.frame();
    let o = circle.center();
    let axe = Vector::from(frame.z());
    let pp = projected(p, o, axe);
    let opp = pp - o;
    if opp.length() < tol {
        return None;
    }
    let mut usol = [0.0; 2];
    // -PI < U1 < PI.
    usol[0] = angle_with_ref(Vector::from(frame.x()), Vector::from(Direction::from_vector(opp)), Vector::from(Direction::from_vector(axe)));
    if usol[0] + PI < ANGULAR {
        usol[0] = -PI;
    } else if usol[0] - PI > -ANGULAR {
        usol[0] = PI;
    }
    usol[1] = usol[0] + PI;

    let radius = circle.radius();
    let mut tol_u = INFINITE;
    if radius > RESOLUTION {
        tol_u = tol / radius;
    }
    let mut my_uinf = uinf;
    (my_uinf, usol[0]) = adjust_periodic(uinf, uinf + 2.0 * PI, tol_u, my_uinf, usol[0]);
    (_, usol[1]) = adjust_periodic(uinf, uinf + 2.0 * PI, tol_u, my_uinf, usol[1]);
    for u in &mut usol {
        if (*u - 2.0 * PI - uinf) < tol_u && (*u - 2.0 * PI - uinf) > -tol_u {
            *u = uinf;
        }
    }

    let mut extrema = Vec::new();
    for (index, &us) in usol.iter().enumerate() {
        if (uinf - us) < tol_u && (us - usup) < tol_u {
            let cu = circle.point(us);
            extrema.push(ExtPCPoint { parameter: us, point: cu, square_distance: cu.distance_squared(p), is_min: index == 0 });
        }
    }
    Some(extrema)
}

/// `Extrema_ExtPElC::Perform(P, gp_Elips, Tol, Uinf, Usup)`: the roots in `[uinf, usup]` of
/// `(B² - A²) cos u sin u - B Y cos u + A X sin u` (`math_TrigonometricFunctionRoots`).
/// Not done for a circular ellipse with `P` on its axis; nor where the equation vanishes
/// identically (where OCCT's `NbSolutions` throws `StdFail_InfiniteSolutions`).
#[doc(alias = "Extrema_ExtPElC")]
fn ext_p_elc_ellipse(p: Point, ellipse: &Ellipse, tol: f64, uinf: f64, usup: f64) -> Option<Vec<ExtPCPoint>> {
    let frame = ellipse.frame();
    let o = ellipse.center();
    let axe = Vector::from(frame.z());
    let pp = projected(p, o, axe);

    let a = ellipse.major_radius();
    let b = ellipse.minor_radius();
    let opp = pp - o;
    let opp_magnitude = opp.length();
    if opp_magnitude < tol && (a - b).abs() < tol {
        return None;
    }
    let x = opp.dot(Vector::from(frame.x()));
    let y = opp.dot(Vector::from(frame.y()));

    let ko2 = (b * b - a * a) / 2.0;
    let mut ko3 = -b * y;
    let ko4 = a * x;
    if ko3.abs() < 1.0e-16 * std_max(ko2.abs(), ko3.abs()) {
        ko3 = 0.0;
    }

    let solution = trigonometric_roots(0.0, ko2, ko3, ko4, 0.0, uinf, usup);
    if !solution.done || solution.infinite {
        return None;
    }
    let mut extrema = Vec::new();
    for &us in &solution.roots {
        let cu = ellipse.point(us);
        let square_distance = cu.distance_squared(p);
        let next = ellipse.point(us + 0.1);
        extrema.push(ExtPCPoint { parameter: us, point: cu, square_distance, is_min: square_distance < next.distance_squared(p) });
    }
    Some(extrema)
}

/// `Extrema_ExtPElC::Perform(P, gp_Hypr, Tol, Uinf, Usup)`: `u = ln v` for the positive
/// roots `v` of `(R² + r²)/4 v⁴ - (XR + Yr)/2 v³ + (XR - Yr)/2 v - (R² + r²)/4`
/// (`math_DirectPolynomialRoots`) in `[uinf, usup]`, points within `tol` of one kept once.
#[doc(alias = "Extrema_ExtPElC")]
fn ext_p_elc_hyperbola(p: Point, hyperbola: &Hyperbola, tol: f64, uinf: f64, usup: f64) -> Option<Vec<ExtPCPoint>> {
    let frame = hyperbola.frame();
    let o = hyperbola.center();
    let axe = Vector::from(frame.z());
    let pp = projected(p, o, axe);

    let tol2 = tol * tol;
    let big_r = hyperbola.major_radius();
    let r = hyperbola.minor_radius();
    let opp = pp - o;
    let x = opp.dot(Vector::from(frame.x()));
    let y = opp.dot(Vector::from(frame.y()));

    let c1 = (big_r * big_r + r * r) / 4.0;
    let solution = PolynomialRoots::quartic(c1, -(x * big_r + y * r) / 2.0, 0.0, (x * big_r - y * r) / 2.0, -c1);
    if !solution.done || solution.infinite {
        return None;
    }
    let mut extrema: Vec<ExtPCPoint> = Vec::new();
    for &vs in &solution.roots {
        if vs > 0.0 {
            let us = vs.ln();
            if us >= uinf && us <= usup {
                let cu = hyperbola.point(us);
                // (`TbExt` holds the points already recorded.)
                if !extrema.iter().any(|extremum| extremum.point.distance_squared(cu) < tol2) {
                    let square_distance = cu.distance_squared(p);
                    let is_min = square_distance < p.distance_squared(hyperbola.point(us + 1.0));
                    extrema.push(ExtPCPoint { parameter: us, point: cu, square_distance, is_min });
                }
            }
        }
    }
    Some(extrema)
}

/// `Extrema_ExtPElC::Perform(P, gp_Parab, Tol, Uinf, Usup)`: the roots in `[uinf, usup]` of
/// `u³/(4F) + (2F - X) u - 2FY` (`math_DirectPolynomialRoots`), points within
/// `Precision::Confusion()` of one kept once (the tolerance given is unused, as OCCT's).
#[doc(alias = "Extrema_ExtPElC")]
fn ext_p_elc_parabola(p: Point, parabola: &Parabola, uinf: f64, usup: f64) -> Option<Vec<ExtPCPoint>> {
    let frame = parabola.frame();
    let o = frame.origin();
    let axe = Vector::from(frame.z());
    let pp = projected(p, o, axe);

    let f = parabola.focal();
    let opp = pp - o;
    let x = opp.dot(Vector::from(frame.x()));
    let y = opp.dot(Vector::from(frame.y()));
    let solution = PolynomialRoots::cubic(1.0 / (4.0 * f), 0.0, 2.0 * f - x, -2.0 * f * y);
    if !solution.done || solution.infinite {
        return None;
    }
    let mut extrema: Vec<ExtPCPoint> = Vec::new();
    for &us in &solution.roots {
        if us >= uinf && us <= usup {
            let cu = parabola.point(us);
            // (`TbExt` holds the points already recorded.)
            if !extrema.iter().any(|extremum| extremum.point.distance_squared(cu) < SQUARE_CONFUSION) {
                let square_distance = cu.distance_squared(p);
                let is_min = square_distance < p.distance_squared(parabola.point(us + 1.0));
                extrema.push(ExtPCPoint { parameter: us, point: cu, square_distance, is_min });
            }
        }
    }
    Some(extrema)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::adaptor::GeomAdaptor;
    use crate::geometry::curve::{BSplineCurve, BezierCurve, OffsetCurve};
    use crate::math::Frame;

    /// The B-spline of the OCCT references: degree 3, knots 0 to 4 of multiplicities
    /// 4 2 3 1 4 (C1 at 1, C0 at 2, C2 at 3), so the `C2` intervals [0, 1], [1, 2], [2, 4].
    fn bspline_poles() -> Vec<Point> {
        [(0.0, 0.0, 0.0), (1.0, 2.0, 0.0), (2.0, -1.0, 0.5), (3.0, 1.0, 1.0), (4.0, 3.0, 0.0), (5.0, 0.0, -1.0), (6.0, 1.0, 0.0), (7.0, 2.0, 1.0), (8.0, 0.0, 0.5), (9.0, 1.0, 0.0)]
            .iter()
            .map(|&(x, y, z)| Point::new(x, y, z))
            .collect()
    }

    /// The curve of a reference case by name, as `GeomAdaptor_Curve(C, first, last)`.
    fn adaptor(name: &str) -> GeomAdaptor {
        let frame = Frame::new(Point::new(1.0, 2.0, 3.0), Direction::new(0.0, 0.6, 0.8), Direction::X);
        let knots = vec![0.0, 1.0, 2.0, 3.0, 4.0];
        let (curve, first, last) = match name {
            "bspline" | "bspline_sub" => (Curve::BSpline(BSplineCurve::new(bspline_poles(), knots, vec![4, 2, 3, 1, 4], 3)), 0.0, 4.0),
            "bezier" => {
                let poles = vec![Point::new(0.0, 0.0, 0.0), Point::new(1.0, 3.0, 0.0), Point::new(2.0, -2.0, 1.0), Point::new(3.0, 3.0, 0.0), Point::new(4.0, 0.0, 0.0)];
                (Curve::Bezier(BezierCurve::new(poles)), 0.0, 1.0)
            }
            "circle" => (Curve::Circle(Circle::new(frame, 2.0)), 0.0, 2.0 * PI),
            "circle_arc" => (Curve::Circle(Circle::new(frame, 2.0)), 0.5, 4.0),
            "line" => (Curve::Line(Line::new(Point::new(1.0, 0.0, 2.0), Direction::new(0.6, 0.0, 0.8))), -5.0, 5.0),
            "ellipse" => (Curve::Ellipse(Ellipse::new(frame, 3.0, 1.0)), 0.0, 2.0 * PI),
            "hyperbola" => (Curve::Hyperbola(Hyperbola::new(frame, 2.0, 1.0)), -2.0, 2.0),
            "parabola" => (Curve::Parabola(Parabola::new(frame, 0.5)), -3.0, 3.0),
            "offset" => {
                let poles = bspline_poles().iter().map(|p| Point::new(p.x, p.y, 0.0)).collect();
                let basis = Curve::BSpline(BSplineCurve::new(poles, knots, vec![4, 2, 2, 2, 4], 3));
                (Curve::Offset(OffsetCurve::new(basis, 0.3, Direction::Z)), 0.0, 4.0)
            }
            _ => unreachable!(),
        };
        GeomAdaptor::new(&curve, first, last)
    }

    type Expected = Option<(f64, [f64; 3], f64, bool)>;
    /// The curve's name, `P`, `U0`, `Umin`, `Usup` and OCCT's result.
    type Case = (&'static str, [f64; 3], f64, f64, f64, Expected);

    /// `Extrema_LocateExtPC(P, GeomAdaptor_Curve(C, first, last), U0, Umin, Usup, 1e-9)`
    /// from OCCT 8.1: `IsDone()`, then `Point().Parameter()`, `Point().Value()`,
    /// `SquareDistance()` and `IsMin()`.
    // (Values as OCCT printed them, π and 2π among them.)
    #[rustfmt::skip]
    #[allow(clippy::approx_constant)]
    const OCCT_CASES: &[Case] = &[
        ("bspline", [3.0, 3.0, 0.0], 0.5, 0.0, 4.0, Some((0.42510371023737975, [1.2369002123623105, 0.5313213248574051, 0.21345336916163593], 9.248457603078316, true))),
        ("bspline", [3.0, 3.0, 0.0], 1.0, 0.0, 4.0, Some((2.2343024951660184, [5.626992325905995, 0.6173453320017858, -0.37783117104610964], 12.720888341116744, true))),
        ("bspline", [3.0, 3.0, 0.0], 1.5, 0.0, 4.0, Some((1.5465239226725231, [3.6861983036083803, 1.5561948324122827, 0.24386189300540395], 2.6149100966881944, true))),
        ("bspline", [3.0, 3.0, 0.0], 1.9999999, 0.0, 4.0, Some((1.5465239226725231, [3.6861983036083803, 1.5561948324122827, 0.24386189300540395], 2.6149100966881944, true))),
        ("bspline", [3.0, 3.0, 0.0], 2.0, 0.0, 4.0, Some((2.2343024951660184, [5.626992325905995, 0.6173453320017858, -0.37783117104610964], 12.720888341116744, true))),
        ("bspline", [3.0, 3.0, 0.0], 2.0000001, 0.0, 4.0, Some((2.2343024951660184, [5.626992325905995, 0.6173453320017858, -0.37783117104610964], 12.720888341116744, true))),
        ("bspline", [3.0, 3.0, 0.0], 2.7, 0.0, 4.0, Some((2.2343024951660184, [5.626992325905995, 0.6173453320017858, -0.37783117104610964], 12.720888341116744, true))),
        ("bspline", [3.0, 3.0, 0.0], 3.0, 0.0, 4.0, Some((2.2343024951660184, [5.626992325905995, 0.6173453320017858, -0.37783117104610964], 12.720888341116744, true))),
        ("bspline", [3.0, 3.0, 0.0], 3.6, 0.0, 4.0, Some((1.5465239226725231, [3.6861983036083803, 1.5561948324122827, 0.24386189300540395], 2.6149100966881944, true))),
        ("bspline", [3.0, 3.0, 0.0], -1.0, 0.0, 4.0, Some((0.42510371023737975, [1.2369002123623105, 0.5313213248574051, 0.21345336916163593], 9.248457603078316, true))),
        ("bspline", [3.0, 3.0, 0.0], 5.0, 0.0, 4.0, None),
        ("bspline", [1.0, 0.0, 0.0], 0.5, 0.0, 4.0, Some((0.4444727247247247, [1.2895140466992234, 0.49377059753414704, 0.2304778133239053], 0.38074780865995766, true))),
        ("bspline", [1.0, 0.0, 0.0], 1.0, 0.0, 4.0, Some((0.4444727247247247, [1.2895140466992234, 0.49377059753414704, 0.2304778133239053], 0.38074780865995766, true))),
        ("bspline", [1.0, 0.0, 0.0], 1.5, 0.0, 4.0, Some((0.4444727247247247, [1.2895140466992234, 0.49377059753414704, 0.2304778133239053], 0.38074780865995766, true))),
        ("bspline", [1.0, 0.0, 0.0], 1.9999999, 0.0, 4.0, Some((0.4444727247247247, [1.2895140466992234, 0.49377059753414704, 0.2304778133239053], 0.38074780865995766, true))),
        ("bspline", [1.0, 0.0, 0.0], 2.0, 0.0, 4.0, None),
        ("bspline", [1.0, 0.0, 0.0], 2.0000001, 0.0, 4.0, None),
        ("bspline", [1.0, 0.0, 0.0], 2.7, 0.0, 4.0, None),
        ("bspline", [1.0, 0.0, 0.0], 3.0, 0.0, 4.0, None),
        ("bspline", [1.0, 0.0, 0.0], 3.6, 0.0, 4.0, None),
        ("bspline", [1.0, 0.0, 0.0], -1.0, 0.0, 4.0, None),
        ("bspline", [1.0, 0.0, 0.0], 5.0, 0.0, 4.0, None),
        ("bspline", [5.0, 2.0, 0.0], 0.5, 0.0, 4.0, Some((1.7457952228683342, [4.2455990336881655, 1.4171023444345503, -0.2579190813129102], 0.9754127473412229, true))),
        ("bspline", [5.0, 2.0, 0.0], 1.0, 0.0, 4.0, Some((2.385199738270018, [5.961608702116293, 0.9187420848395315, -0.059824606522087764], 2.0973889586284677, true))),
        ("bspline", [5.0, 2.0, 0.0], 1.5, 0.0, 4.0, Some((1.7457952228683342, [4.2455990336881655, 1.4171023444345503, -0.2579190813129102], 0.9754127473412229, true))),
        ("bspline", [5.0, 2.0, 0.0], 1.9999999, 0.0, 4.0, Some((1.7457952228683342, [4.2455990336881655, 1.4171023444345503, -0.2579190813129102], 0.9754127473412229, true))),
        ("bspline", [5.0, 2.0, 0.0], 2.0, 0.0, 4.0, Some((2.385199738270018, [5.961608702116293, 0.9187420848395315, -0.059824606522087764], 2.0973889586284677, true))),
        ("bspline", [5.0, 2.0, 0.0], 2.0000001, 0.0, 4.0, Some((2.385199738270018, [5.961608702116293, 0.9187420848395315, -0.059824606522087764], 2.0973889586284677, true))),
        ("bspline", [5.0, 2.0, 0.0], 2.7, 0.0, 4.0, Some((2.385199738270018, [5.961608702116293, 0.9187420848395315, -0.059824606522087764], 2.0973889586284677, true))),
        ("bspline", [5.0, 2.0, 0.0], 3.0, 0.0, 4.0, Some((2.385199738270018, [5.961608702116293, 0.9187420848395315, -0.059824606522087764], 2.0973889586284677, true))),
        ("bspline", [5.0, 2.0, 0.0], 3.6, 0.0, 4.0, Some((1.7457952228683342, [4.2455990336881655, 1.4171023444345503, -0.2579190813129102], 0.9754127473412229, true))),
        ("bspline", [5.0, 2.0, 0.0], -1.0, 0.0, 4.0, Some((1.7457952228683342, [4.2455990336881655, 1.4171023444345503, -0.2579190813129102], 0.9754127473412229, true))),
        ("bspline", [5.0, 2.0, 0.0], 5.0, 0.0, 4.0, None),
        ("bspline", [8.0, -1.0, 0.0], 0.5, 0.0, 4.0, None),
        ("bspline", [8.0, -1.0, 0.0], 1.0, 0.0, 4.0, None),
        ("bspline", [8.0, -1.0, 0.0], 1.5, 0.0, 4.0, None),
        ("bspline", [8.0, -1.0, 0.0], 1.9999999, 0.0, 4.0, None),
        ("bspline", [8.0, -1.0, 0.0], 2.0, 0.0, 4.0, None),
        ("bspline", [8.0, -1.0, 0.0], 2.0000001, 0.0, 4.0, None),
        ("bspline", [8.0, -1.0, 0.0], 2.7, 0.0, 4.0, None),
        ("bspline", [8.0, -1.0, 0.0], 3.0, 0.0, 4.0, None),
        ("bspline", [8.0, -1.0, 0.0], 3.6, 0.0, 4.0, Some((3.679718518509398, [8.17659862300753, 0.564286915046461, 0.39938023329669553], 2.6376851969818524, true))),
        ("bspline", [8.0, -1.0, 0.0], -1.0, 0.0, 4.0, None),
        ("bspline", [8.0, -1.0, 0.0], 5.0, 0.0, 4.0, Some((3.679718518509398, [8.17659862300753, 0.564286915046461, 0.39938023329669553], 2.6376851969818524, true))),
        ("bspline", [2.5, 0.5, 2.0], 0.5, 0.0, 4.0, Some((1.1151582181354716, [2.691865869536981, 0.37609593434964694, 0.7885473077671693], 1.5197823548960423, true))),
        ("bspline", [2.5, 0.5, 2.0], 1.0, 0.0, 4.0, Some((1.1151582181354716, [2.691865869536981, 0.37609593434964694, 0.7885473077671693], 1.5197823548960423, true))),
        ("bspline", [2.5, 0.5, 2.0], 1.5, 0.0, 4.0, Some((1.1151582181354716, [2.691865869536981, 0.37609593434964694, 0.7885473077671693], 1.5197823548960423, true))),
        ("bspline", [2.5, 0.5, 2.0], 1.9999999, 0.0, 4.0, Some((2.115863381188173, [5.32823135264627, 0.327064812047843, -0.6723519176529434], 15.170263935106762, true))),
        ("bspline", [2.5, 0.5, 2.0], 2.0, 0.0, 4.0, Some((2.115863381188173, [5.32823135264627, 0.327064812047843, -0.6723519176529434], 15.170263935106762, true))),
        ("bspline", [2.5, 0.5, 2.0], 2.0000001, 0.0, 4.0, Some((2.115863381188173, [5.32823135264627, 0.327064812047843, -0.6723519176529434], 15.170263935106762, true))),
        ("bspline", [2.5, 0.5, 2.0], 2.7, 0.0, 4.0, Some((2.115863381188173, [5.32823135264627, 0.327064812047843, -0.6723519176529434], 15.170263935106762, true))),
        ("bspline", [2.5, 0.5, 2.0], 3.0, 0.0, 4.0, Some((2.115863381188173, [5.32823135264627, 0.327064812047843, -0.6723519176529434], 15.170263935106762, true))),
        ("bspline", [2.5, 0.5, 2.0], 3.6, 0.0, 4.0, Some((2.115863381188173, [5.32823135264627, 0.327064812047843, -0.6723519176529434], 15.170263935106762, true))),
        ("bspline", [2.5, 0.5, 2.0], -1.0, 0.0, 4.0, Some((1.1151582181354716, [2.691865869536981, 0.37609593434964694, 0.7885473077671693], 1.5197823548960423, true))),
        ("bspline", [2.5, 0.5, 2.0], 5.0, 0.0, 4.0, None),
        ("bspline", [6.0, 5.0, 0.0], 0.5, 0.0, 4.0, Some((0.5460825525852128, [1.5568250688041898, 0.2690096030676219, 0.32517534792979996], 42.22981261197454, true))),
        ("bspline", [6.0, 5.0, 0.0], 1.0, 0.0, 4.0, Some((2.7378049407954657, [6.597694947770533, 1.2964734657043777, 0.4470842067374552], 14.27323232873622, true))),
        ("bspline", [6.0, 5.0, 0.0], 1.5, 0.0, 4.0, Some((1.692713457420534, [4.092648141110836, 1.523296820649371, -0.11440979438468701], 15.738545711965866, true))),
        ("bspline", [6.0, 5.0, 0.0], 1.9999999, 0.0, 4.0, Some((1.692713457420534, [4.092648141110836, 1.523296820649371, -0.11440979438468701], 15.738545711965866, true))),
        ("bspline", [6.0, 5.0, 0.0], 2.0, 0.0, 4.0, Some((2.7378049407954657, [6.597694947770533, 1.2964734657043777, 0.4470842067374552], 14.27323232873622, true))),
        ("bspline", [6.0, 5.0, 0.0], 2.0000001, 0.0, 4.0, Some((2.7378049407954657, [6.597694947770533, 1.2964734657043777, 0.4470842067374552], 14.27323232873622, true))),
        ("bspline", [6.0, 5.0, 0.0], 2.7, 0.0, 4.0, Some((2.737804940795466, [6.597694947770534, 1.2964734657043777, 0.44708420673745586], 14.273232328736222, true))),
        ("bspline", [6.0, 5.0, 0.0], 3.0, 0.0, 4.0, Some((2.7378049407954657, [6.597694947770533, 1.2964734657043777, 0.4470842067374552], 14.27323232873622, true))),
        ("bspline", [6.0, 5.0, 0.0], 3.6, 0.0, 4.0, Some((3.866860948896683, [8.625991845211917, 0.7004488229707726, 0.18611906849066473], 25.416613802668742, false))),
        ("bspline", [6.0, 5.0, 0.0], -1.0, 0.0, 4.0, Some((0.5460825525852128, [1.5568250688041898, 0.2690096030676219, 0.32517534792979996], 42.22981261197454, true))),
        ("bspline", [6.0, 5.0, 0.0], 5.0, 0.0, 4.0, Some((3.866860948896683, [8.625991845211917, 0.7004488229707726, 0.18611906849066473], 25.416613802668742, false))),
        ("bspline", [0.0, -2.0, 0.0], 0.5, 0.0, 4.0, Some((0.6229623282028294, [1.748006732057752, 0.09238757973188344, 0.4008027147137302], 7.594256135257366, true))),
        ("bspline", [0.0, -2.0, 0.0], 1.0, 0.0, 4.0, Some((0.6229623282028294, [1.748006732057752, 0.09238757973188344, 0.4008027147137302], 7.594256135257366, true))),
        ("bspline", [0.0, -2.0, 0.0], 1.5, 0.0, 4.0, Some((1.8581848522849924, [4.575980617087895, 0.9917760104086443, -0.5781197074372724], 30.224544700548137, false))),
        ("bspline", [0.0, -2.0, 0.0], 1.9999999, 0.0, 4.0, Some((1.9999999988713764, [4.999999996614129, 1.0157612961496683e-08, -0.9999999966141291], 30.0, true))),
        ("bspline", [0.0, -2.0, 0.0], 2.0, 0.0, 4.0, Some((1.8581848522849924, [4.575980617087895, 0.9917760104086443, -0.5781197074372724], 30.224544700548137, false))),
        ("bspline", [0.0, -2.0, 0.0], 2.0000001, 0.0, 4.0, Some((1.8581848522849924, [4.575980617087895, 0.9917760104086443, -0.5781197074372724], 30.224544700548137, false))),
        ("bspline", [0.0, -2.0, 0.0], 2.7, 0.0, 4.0, Some((1.8581848522849924, [4.575980617087895, 0.9917760104086443, -0.5781197074372724], 30.224544700548137, false))),
        ("bspline", [0.0, -2.0, 0.0], 3.0, 0.0, 4.0, Some((1.8581848522849924, [4.575980617087895, 0.9917760104086443, -0.5781197074372724], 30.224544700548137, false))),
        ("bspline", [0.0, -2.0, 0.0], 3.6, 0.0, 4.0, Some((1.8581848522849924, [4.575980617087895, 0.9917760104086443, -0.5781197074372724], 30.224544700548137, false))),
        ("bspline", [0.0, -2.0, 0.0], -1.0, 0.0, 4.0, Some((1.8581848522849924, [4.575980617087895, 0.9917760104086443, -0.5781197074372724], 30.224544700548137, false))),
        ("bspline", [0.0, -2.0, 0.0], 5.0, 0.0, 4.0, None),
        ("bspline", [9.5, 2.0, 0.0], 0.5, 0.0, 4.0, Some((1.9450445086499437, [4.835216511655222, 0.45029366017230593, -0.8353409902133087], 24.859589302766775, true))),
        ("bspline", [9.5, 2.0, 0.0], 1.0, 0.0, 4.0, None),
        ("bspline", [9.5, 2.0, 0.0], 1.5, 0.0, 4.0, Some((1.9450445086499437, [4.835216511655222, 0.45029366017230593, -0.8353409902133087], 24.859589302766775, true))),
        ("bspline", [9.5, 2.0, 0.0], 1.9999999, 0.0, 4.0, Some((1.9450445086499435, [4.835216511655221, 0.4502936601723084, -0.8353409902133078], 24.85958930276677, true))),
        ("bspline", [9.5, 2.0, 0.0], 2.0, 0.0, 4.0, Some((1.9450445086499437, [4.835216511655222, 0.45029366017230593, -0.8353409902133087], 24.859589302766775, true))),
        ("bspline", [9.5, 2.0, 0.0], 2.0000001, 0.0, 4.0, Some((1.9450445086499437, [4.835216511655222, 0.45029366017230593, -0.8353409902133087], 24.859589302766775, true))),
        ("bspline", [9.5, 2.0, 0.0], 2.7, 0.0, 4.0, Some((1.9450445086499437, [4.835216511655222, 0.45029366017230593, -0.8353409902133087], 24.859589302766775, true))),
        ("bspline", [9.5, 2.0, 0.0], 3.0, 0.0, 4.0, Some((1.9450445086499437, [4.835216511655222, 0.45029366017230593, -0.8353409902133087], 24.859589302766775, true))),
        ("bspline", [9.5, 2.0, 0.0], 3.6, 0.0, 4.0, Some((1.9450445086499437, [4.835216511655222, 0.45029366017230593, -0.8353409902133087], 24.859589302766775, true))),
        ("bspline", [9.5, 2.0, 0.0], -1.0, 0.0, 4.0, Some((1.9450445086499437, [4.835216511655222, 0.45029366017230593, -0.8353409902133087], 24.859589302766775, true))),
        ("bspline", [9.5, 2.0, 0.0], 5.0, 0.0, 4.0, None),
        ("bspline_sub", [3.0, 3.0, 0.0], 0.7, 0.6, 3.1, Some((0.8291292929829928, [2.202393180502324, -0.20715043729943325, 0.6036910290537352], 11.286433424539101, false))),
        ("bspline_sub", [3.0, 3.0, 0.0], 1.5, 0.6, 3.1, Some((1.5465239226725231, [3.6861983036083803, 1.5561948324122827, 0.24386189300540395], 2.6149100966881944, true))),
        ("bspline_sub", [3.0, 3.0, 0.0], 2.2, 0.6, 3.1, Some((2.2343024951660184, [5.626992325905995, 0.6173453320017858, -0.37783117104610964], 12.720888341116744, true))),
        ("bspline_sub", [3.0, 3.0, 0.0], 2.9, 0.6, 3.1, Some((2.2343024951660184, [5.626992325905995, 0.6173453320017858, -0.37783117104610964], 12.720888341116744, true))),
        ("bspline_sub", [1.0, 0.0, 0.0], 0.7, 0.6, 3.1, None),
        ("bspline_sub", [1.0, 0.0, 0.0], 1.5, 0.6, 3.1, None),
        ("bspline_sub", [1.0, 0.0, 0.0], 2.2, 0.6, 3.1, None),
        ("bspline_sub", [1.0, 0.0, 0.0], 2.9, 0.6, 3.1, None),
        ("bspline_sub", [5.0, 2.0, 0.0], 0.7, 0.6, 3.1, Some((1.7457952228683342, [4.2455990336881655, 1.4171023444345503, -0.2579190813129102], 0.9754127473412229, true))),
        ("bspline_sub", [5.0, 2.0, 0.0], 1.5, 0.6, 3.1, Some((1.7457952228683342, [4.2455990336881655, 1.4171023444345503, -0.2579190813129102], 0.9754127473412229, true))),
        ("bspline_sub", [5.0, 2.0, 0.0], 2.2, 0.6, 3.1, Some((2.385199738270018, [5.961608702116293, 0.9187420848395315, -0.059824606522087764], 2.0973889586284677, true))),
        ("bspline_sub", [5.0, 2.0, 0.0], 2.9, 0.6, 3.1, Some((2.385199738270018, [5.961608702116293, 0.9187420848395315, -0.059824606522087764], 2.0973889586284677, true))),
        ("bspline_sub", [8.0, -1.0, 0.0], 0.7, 0.6, 3.1, None),
        ("bspline_sub", [8.0, -1.0, 0.0], 1.5, 0.6, 3.1, None),
        ("bspline_sub", [8.0, -1.0, 0.0], 2.2, 0.6, 3.1, None),
        ("bspline_sub", [8.0, -1.0, 0.0], 2.9, 0.6, 3.1, None),
        ("bspline_sub", [2.5, 0.5, 2.0], 0.7, 0.6, 3.1, Some((1.1151582181354716, [2.691865869536981, 0.37609593434964694, 0.7885473077671693], 1.5197823548960423, true))),
        ("bspline_sub", [2.5, 0.5, 2.0], 1.5, 0.6, 3.1, Some((1.1151582181354716, [2.691865869536981, 0.37609593434964694, 0.7885473077671693], 1.5197823548960423, true))),
        ("bspline_sub", [2.5, 0.5, 2.0], 2.2, 0.6, 3.1, Some((2.115863381188173, [5.32823135264627, 0.327064812047843, -0.6723519176529434], 15.170263935106762, true))),
        ("bspline_sub", [2.5, 0.5, 2.0], 2.9, 0.6, 3.1, Some((2.115863381188173, [5.32823135264627, 0.327064812047843, -0.6723519176529434], 15.170263935106762, true))),
        ("bspline_sub", [6.0, 5.0, 0.0], 0.7, 0.6, 3.1, Some((0.7092981126148645, [1.9494690454958699, -0.0791132504189194, 0.48701778031542786], 42.44137874232056, false))),
        ("bspline_sub", [6.0, 5.0, 0.0], 1.5, 0.6, 3.1, Some((1.692713457420534, [4.092648141110836, 1.523296820649371, -0.11440979438468701], 15.738545711965866, true))),
        ("bspline_sub", [6.0, 5.0, 0.0], 2.2, 0.6, 3.1, Some((2.737804940795466, [6.597694947770534, 1.2964734657043777, 0.44708420673745586], 14.273232328736222, true))),
        ("bspline_sub", [6.0, 5.0, 0.0], 2.9, 0.6, 3.1, Some((2.737804940795466, [6.597694947770534, 1.2964734657043777, 0.44708420673745586], 14.273232328736222, true))),
        ("bspline_sub", [0.0, -2.0, 0.0], 0.7, 0.6, 3.1, Some((0.6229623282028293, [1.7480067320577517, 0.09238757973188397, 0.4008027147137301], 7.594256135257367, true))),
        ("bspline_sub", [0.0, -2.0, 0.0], 1.5, 0.6, 3.1, Some((1.8581848522849924, [4.575980617087895, 0.9917760104086443, -0.5781197074372724], 30.224544700548137, false))),
        ("bspline_sub", [0.0, -2.0, 0.0], 2.2, 0.6, 3.1, Some((1.8581848522849924, [4.575980617087895, 0.9917760104086443, -0.5781197074372724], 30.224544700548137, false))),
        ("bspline_sub", [0.0, -2.0, 0.0], 2.9, 0.6, 3.1, Some((1.8581848522849924, [4.575980617087895, 0.9917760104086443, -0.5781197074372724], 30.224544700548137, false))),
        ("bspline_sub", [9.5, 2.0, 0.0], 0.7, 0.6, 3.1, Some((1.9450445086499437, [4.835216511655222, 0.45029366017230593, -0.8353409902133087], 24.859589302766775, true))),
        ("bspline_sub", [9.5, 2.0, 0.0], 1.5, 0.6, 3.1, Some((1.9450445086499437, [4.835216511655222, 0.45029366017230593, -0.8353409902133087], 24.859589302766775, true))),
        ("bspline_sub", [9.5, 2.0, 0.0], 2.2, 0.6, 3.1, Some((1.9450445086499437, [4.835216511655222, 0.45029366017230593, -0.8353409902133087], 24.859589302766775, true))),
        ("bspline_sub", [9.5, 2.0, 0.0], 2.9, 0.6, 3.1, Some((1.9450445086499437, [4.835216511655222, 0.45029366017230593, -0.8353409902133087], 24.859589302766775, true))),
        ("bezier", [1.5, 3.0, 0.0], 0.1, 0.0, 1.0, Some((0.2553812390675034, [1.0215249562700135, 0.9801274673329007, 0.21696874571088334], 4.3558988523107685, true))),
        ("bezier", [1.5, 3.0, 0.0], 0.35, 0.0, 1.0, Some((0.2553812390675033, [1.0215249562700133, 0.9801274673329009, 0.21696874571088326], 4.3558988523107685, true))),
        ("bezier", [1.5, 3.0, 0.0], 0.5, 0.0, 1.0, None),
        ("bezier", [1.5, 3.0, 0.0], 0.8, 0.0, 1.0, None),
        ("bezier", [1.5, 3.0, 0.0], 1.2, 0.0, 1.0, None),
        ("bezier", [0.5, 0.0, 0.0], 0.1, 0.0, 1.0, Some((0.014669071506337755, [0.05867628602535102, 0.1659257406277257, 0.001253489590338461], 0.22229954315539016, true))),
        ("bezier", [0.5, 0.0, 0.0], 0.35, 0.0, 1.0, Some((0.014669071506337755, [0.05867628602535102, 0.1659257406277257, 0.001253489590338461], 0.22229954315539016, true))),
        ("bezier", [0.5, 0.0, 0.0], 0.5, 0.0, 1.0, Some((0.014669071506337755, [0.05867628602535102, 0.1659257406277257, 0.001253489590338461], 0.22229954315539016, true))),
        ("bezier", [0.5, 0.0, 0.0], 0.8, 0.0, 1.0, None),
        ("bezier", [0.5, 0.0, 0.0], 1.2, 0.0, 1.0, None),
        ("bezier", [2.5, 2.0, 0.0], 0.1, 0.0, 1.0, Some((0.7160601390653454, [2.8642405562613815, 0.9516404285033431, 0.2480292945302362], 1.2932475049194214, true))),
        ("bezier", [2.5, 2.0, 0.0], 0.35, 0.0, 1.0, Some((0.7160601390653453, [2.864240556261381, 0.9516404285033442, 0.24802929453023656], 1.2932475049194192, true))),
        ("bezier", [2.5, 2.0, 0.0], 0.5, 0.0, 1.0, Some((0.7160601390653455, [2.864240556261382, 0.9516404285033457, 0.24802929453023628], 1.293247504919416, true))),
        ("bezier", [2.5, 2.0, 0.0], 0.8, 0.0, 1.0, Some((0.7160601390653455, [2.864240556261382, 0.9516404285033457, 0.24802929453023628], 1.293247504919416, true))),
        ("bezier", [2.5, 2.0, 0.0], 1.2, 0.0, 1.0, Some((0.7160601390653453, [2.864240556261381, 0.9516404285033442, 0.24802929453023656], 1.2932475049194192, true))),
        ("bezier", [4.0, -1.0, 0.0], 0.1, 0.0, 1.0, None),
        ("bezier", [4.0, -1.0, 0.0], 0.35, 0.0, 1.0, None),
        ("bezier", [4.0, -1.0, 0.0], 0.5, 0.0, 1.0, None),
        ("bezier", [4.0, -1.0, 0.0], 0.8, 0.0, 1.0, None),
        ("bezier", [4.0, -1.0, 0.0], 1.2, 0.0, 1.0, None),
        ("bezier", [1.25, 0.5, 2.0], 0.1, 0.0, 1.0, Some((0.395665919930308, [1.582663679721232, 0.8110477350320685, 0.34305417696667734], 2.9528852777418253, true))),
        ("bezier", [1.25, 0.5, 2.0], 0.35, 0.0, 1.0, Some((0.39566591993030803, [1.5826636797212321, 0.8110477350320686, 0.3430541769666775], 2.9528852777418253, true))),
        ("bezier", [1.25, 0.5, 2.0], 0.5, 0.0, 1.0, Some((0.395665919930308, [1.582663679721232, 0.8110477350320685, 0.34305417696667734], 2.9528852777418253, true))),
        ("bezier", [1.25, 0.5, 2.0], 0.8, 0.0, 1.0, None),
        ("bezier", [1.25, 0.5, 2.0], 1.2, 0.0, 1.0, None),
        ("bezier", [3.0, 5.0, 0.0], 0.1, 0.0, 1.0, None),
        ("bezier", [3.0, 5.0, 0.0], 0.35, 0.0, 1.0, None),
        ("bezier", [3.0, 5.0, 0.0], 0.5, 0.0, 1.0, None),
        ("bezier", [3.0, 5.0, 0.0], 0.8, 0.0, 1.0, Some((0.7848487222815443, [3.1393948891261774, 0.9998266239886231, 0.17108464016437336], 16.050087927364928, true))),
        ("bezier", [3.0, 5.0, 0.0], 1.2, 0.0, 1.0, Some((0.7848487222815442, [3.139394889126177, 0.999826623988623, 0.17108464016437333], 16.050087927364928, true))),
        ("bezier", [0.0, -2.0, 0.0], 0.1, 0.0, 1.0, None),
        ("bezier", [0.0, -2.0, 0.0], 0.35, 0.0, 1.0, None),
        ("bezier", [0.0, -2.0, 0.0], 0.5, 0.0, 1.0, None),
        ("bezier", [0.0, -2.0, 0.0], 0.8, 0.0, 1.0, Some((0.9164993135773571, [3.6659972543094286, 0.7075028361815638, 0.03513950421558657], 20.771342261291995, false))),
        ("bezier", [0.0, -2.0, 0.0], 1.2, 0.0, 1.0, None),
        ("bezier", [4.75, 2.0, 0.0], 0.1, 0.0, 1.0, Some((0.8930412872638701, [3.5721651490554804, 0.8177660199135166, 0.05474275635966745], 2.78796888914448, true))),
        ("bezier", [4.75, 2.0, 0.0], 0.35, 0.0, 1.0, Some((0.8930412872638701, [3.5721651490554804, 0.8177660199135166, 0.05474275635966745], 2.78796888914448, true))),
        ("bezier", [4.75, 2.0, 0.0], 0.5, 0.0, 1.0, Some((0.8930412872638701, [3.5721651490554804, 0.8177660199135166, 0.05474275635966745], 2.78796888914448, true))),
        ("bezier", [4.75, 2.0, 0.0], 0.8, 0.0, 1.0, Some((0.8930412872638699, [3.5721651490554795, 0.8177660199135132, 0.054742756359667424], 2.78796888914449, true))),
        ("bezier", [4.75, 2.0, 0.0], 1.2, 0.0, 1.0, Some((0.8930412872638701, [3.5721651490554804, 0.8177660199135166, 0.05474275635966745], 2.78796888914448, true))),
        ("circle", [4.0, 2.0, 3.0], 0.0, 0.0, 6.283185307179586, Some((0.0, [3.0, 2.0, 3.0], 1.0, true))),
        ("circle_arc", [4.0, 2.0, 3.0], 0.0, 0.5, 4.0, Some((3.141592653589793, [-1.0, 2.0, 3.0], 25.0, false))),
        ("circle", [4.0, 2.0, 3.0], 1.0, 0.0, 6.283185307179586, Some((0.0, [3.0, 2.0, 3.0], 1.0, true))),
        ("circle_arc", [4.0, 2.0, 3.0], 1.0, 0.5, 4.0, Some((3.141592653589793, [-1.0, 2.0, 3.0], 25.0, false))),
        ("circle", [4.0, 2.0, 3.0], 3.0, 0.0, 6.283185307179586, Some((3.141592653589793, [-1.0, 2.0, 3.0], 25.0, false))),
        ("circle_arc", [4.0, 2.0, 3.0], 3.0, 0.5, 4.0, Some((3.141592653589793, [-1.0, 2.0, 3.0], 25.0, false))),
        ("circle", [4.0, 2.0, 3.0], 5.5, 0.0, 6.283185307179586, Some((3.141592653589793, [-1.0, 2.0, 3.0], 25.0, false))),
        ("circle_arc", [4.0, 2.0, 3.0], 5.5, 0.5, 4.0, Some((3.141592653589793, [-1.0, 2.0, 3.0], 25.0, false))),
        ("circle", [1.0, 5.0, 1.0], 0.0, 0.0, 6.283185307179586, Some((1.5707963267948966, [1.0000000000000002, 3.6, 1.8], 2.5999999999999996, true))),
        ("circle_arc", [1.0, 5.0, 1.0], 0.0, 0.5, 4.0, Some((1.5707963267948966, [1.0000000000000002, 3.6, 1.8], 2.5999999999999996, true))),
        ("circle", [1.0, 5.0, 1.0], 1.0, 0.0, 6.283185307179586, Some((1.5707963267948966, [1.0000000000000002, 3.6, 1.8], 2.5999999999999996, true))),
        ("circle_arc", [1.0, 5.0, 1.0], 1.0, 0.5, 4.0, Some((1.5707963267948966, [1.0000000000000002, 3.6, 1.8], 2.5999999999999996, true))),
        ("circle", [1.0, 5.0, 1.0], 3.0, 0.0, 6.283185307179586, Some((1.5707963267948966, [1.0000000000000002, 3.6, 1.8], 2.5999999999999996, true))),
        ("circle_arc", [1.0, 5.0, 1.0], 3.0, 0.5, 4.0, Some((1.5707963267948966, [1.0000000000000002, 3.6, 1.8], 2.5999999999999996, true))),
        ("circle", [1.0, 5.0, 1.0], 5.5, 0.0, 6.283185307179586, Some((4.71238898038469, [0.9999999999999997, 0.3999999999999999, 4.2], 31.4, false))),
        ("circle_arc", [1.0, 5.0, 1.0], 5.5, 0.5, 4.0, Some((1.5707963267948966, [1.0000000000000002, 3.6, 1.8], 2.5999999999999996, true))),
        ("circle", [-3.0, 1.0, 4.0], 0.0, 0.0, 6.283185307179586, Some((0.3366748193867273, [2.8877167127320345, 2.5285606795649698, 2.603579490326273], 38.95169608033367, false))),
        ("circle_arc", [-3.0, 1.0, 4.0], 0.0, 0.5, 4.0, Some((3.4782674729765204, [-0.8877167127320349, 1.4714393204350302, 3.396420509673727], 5.048303919666328, true))),
        ("circle", [-3.0, 1.0, 4.0], 1.0, 0.0, 6.283185307179586, Some((0.3366748193867273, [2.8877167127320345, 2.5285606795649698, 2.603579490326273], 38.95169608033367, false))),
        ("circle_arc", [-3.0, 1.0, 4.0], 1.0, 0.5, 4.0, Some((3.4782674729765204, [-0.8877167127320349, 1.4714393204350302, 3.396420509673727], 5.048303919666328, true))),
        ("circle", [-3.0, 1.0, 4.0], 3.0, 0.0, 6.283185307179586, Some((3.4782674729765204, [-0.8877167127320349, 1.4714393204350302, 3.396420509673727], 5.048303919666328, true))),
        ("circle_arc", [-3.0, 1.0, 4.0], 3.0, 0.5, 4.0, Some((3.4782674729765204, [-0.8877167127320349, 1.4714393204350302, 3.396420509673727], 5.048303919666328, true))),
        ("circle", [-3.0, 1.0, 4.0], 5.5, 0.0, 6.283185307179586, Some((3.4782674729765204, [-0.8877167127320349, 1.4714393204350302, 3.396420509673727], 5.048303919666328, true))),
        ("circle_arc", [-3.0, 1.0, 4.0], 5.5, 0.5, 4.0, Some((3.4782674729765204, [-0.8877167127320349, 1.4714393204350302, 3.396420509673727], 5.048303919666328, true))),
        ("circle", [1.0, 2.0, 3.0000001], 0.0, 0.0, 6.283185307179586, None),
        ("circle_arc", [1.0, 2.0, 3.0000001], 0.0, 0.5, 4.0, None),
        ("circle", [1.0, 2.0, 3.0000001], 1.0, 0.0, 6.283185307179586, None),
        ("circle_arc", [1.0, 2.0, 3.0000001], 1.0, 0.5, 4.0, None),
        ("circle", [1.0, 2.0, 3.0000001], 3.0, 0.0, 6.283185307179586, None),
        ("circle_arc", [1.0, 2.0, 3.0000001], 3.0, 0.5, 4.0, None),
        ("circle", [1.0, 2.0, 3.0000001], 5.5, 0.0, 6.283185307179586, None),
        ("circle_arc", [1.0, 2.0, 3.0000001], 5.5, 0.5, 4.0, None),
        ("circle", [2.0, 2.0, 10.0], 0.0, 0.0, 6.283185307179586, Some((1.804539507663798, [0.5367589453879207, 3.556489943496587, 1.8326325423775598], 71.26962651593831, false))),
        ("circle_arc", [2.0, 2.0, 10.0], 0.0, 0.5, 4.0, Some((1.804539507663798, [0.5367589453879207, 3.556489943496587, 1.8326325423775598], 71.26962651593831, false))),
        ("circle", [2.0, 2.0, 10.0], 1.0, 0.0, 6.283185307179586, Some((1.804539507663798, [0.5367589453879207, 3.556489943496587, 1.8326325423775598], 71.26962651593831, false))),
        ("circle_arc", [2.0, 2.0, 10.0], 1.0, 0.5, 4.0, Some((1.804539507663798, [0.5367589453879207, 3.556489943496587, 1.8326325423775598], 71.26962651593831, false))),
        ("circle", [2.0, 2.0, 10.0], 3.0, 0.0, 6.283185307179586, Some((1.804539507663798, [0.5367589453879207, 3.556489943496587, 1.8326325423775598], 71.26962651593831, false))),
        ("circle_arc", [2.0, 2.0, 10.0], 3.0, 0.5, 4.0, Some((1.804539507663798, [0.5367589453879207, 3.556489943496587, 1.8326325423775598], 71.26962651593831, false))),
        ("circle", [2.0, 2.0, 10.0], 5.5, 0.0, 6.283185307179586, Some((4.946132161253591, [1.4632410546120795, 0.4435100565034129, 4.16736745762244], 36.73037348406168, true))),
        ("circle_arc", [2.0, 2.0, 10.0], 5.5, 0.5, 4.0, Some((1.804539507663798, [0.5367589453879207, 3.556489943496587, 1.8326325423775598], 71.26962651593831, false))),
        ("line", [4.0, 2.0, 3.0], -1.0, -5.0, 5.0, Some((2.5999999999999996, [2.5599999999999996, 0.0, 4.08], 7.240000000000001, true))),
        ("line", [4.0, 2.0, 3.0], 2.0, -5.0, 5.0, Some((2.5999999999999996, [2.5599999999999996, 0.0, 4.08], 7.240000000000001, true))),
        ("line", [1.0, 5.0, 1.0], -1.0, -5.0, 5.0, Some((-0.8, [0.52, 0.0, 1.3599999999999999], 25.36, true))),
        ("line", [1.0, 5.0, 1.0], 2.0, -5.0, 5.0, Some((-0.8, [0.52, 0.0, 1.3599999999999999], 25.36, true))),
        ("line", [-3.0, 1.0, 4.0], -1.0, -5.0, 5.0, Some((-0.7999999999999998, [0.5200000000000001, 0.0, 1.36], 20.36, true))),
        ("line", [-3.0, 1.0, 4.0], 2.0, -5.0, 5.0, Some((-0.7999999999999998, [0.5200000000000001, 0.0, 1.36], 20.36, true))),
        ("line", [1.0, 2.0, 3.0000001], -1.0, -5.0, 5.0, Some((0.8000000799999999, [1.480000048, 0.0, 2.640000064], 4.360000072000004, true))),
        ("line", [1.0, 2.0, 3.0000001], 2.0, -5.0, 5.0, Some((0.8000000799999999, [1.480000048, 0.0, 2.640000064], 4.360000072000004, true))),
        ("line", [2.0, 2.0, 10.0], -1.0, -5.0, 5.0, None),
        ("line", [2.0, 2.0, 10.0], 2.0, -5.0, 5.0, None),
        ("line", [30.0, 0.0, 0.0], 0.0, -5.0, 5.0, None),
        ("ellipse", [4.0, 2.0, 3.0], 0.0, 0.0, 6.283185307179586, Some((0.0, [4.0, 2.0, 3.0], 0.0, true))),
        ("ellipse", [4.0, 2.0, 3.0], 1.0, 0.0, 6.283185307179586, Some((0.0, [4.0, 2.0, 3.0], 0.0, true))),
        ("ellipse", [4.0, 2.0, 3.0], 3.0, 0.0, 6.283185307179586, Some((3.141592653589793, [-2.0, 2.0, 3.0], 36.0, false))),
        ("ellipse", [4.0, 2.0, 3.0], 5.5, 0.0, 6.283185307179586, Some((3.141592653589793, [-2.0, 2.0, 3.0], 36.0, false))),
        ("ellipse", [1.0, 5.0, 1.0], 0.0, 0.0, 6.283185307179586, Some((1.5707963267948966, [1.0000000000000002, 2.8, 2.4], 6.800000000000001, true))),
        ("ellipse", [1.0, 5.0, 1.0], 1.0, 0.0, 6.283185307179586, Some((1.5707963267948966, [1.0000000000000002, 2.8, 2.4], 6.800000000000001, true))),
        ("ellipse", [1.0, 5.0, 1.0], 3.0, 0.0, 6.283185307179586, Some((3.6083579926370897, [-1.6790856649237629, 1.64, 3.27], 23.62, false))),
        ("ellipse", [1.0, 5.0, 1.0], 5.5, 0.0, 6.283185307179586, Some((5.81641996813229, [3.6790856649237624, 1.6399999999999997, 3.2700000000000005], 23.620000000000005, false))),
        ("ellipse", [-3.0, 1.0, 4.0], 0.0, 0.0, 6.283185307179586, Some((0.06995421647745646, [3.992662604313123, 2.0559177406924283, 2.9580616944806786], 51.097928005376666, false))),
        ("ellipse", [-3.0, 1.0, 4.0], 1.0, 0.0, 6.283185307179586, Some((0.06995421647745646, [3.992662604313123, 2.0559177406924283, 2.9580616944806786], 51.097928005376666, false))),
        ("ellipse", [-3.0, 1.0, 4.0], 3.0, 0.0, 6.283185307179586, Some((3.4509756461041174, [-1.8575648361805293, 1.7564232148227155, 3.1826825888829635], 2.5453419339690075, true))),
        ("ellipse", [-3.0, 1.0, 4.0], 5.5, 0.0, 6.283185307179586, Some((3.4509756461041174, [-1.8575648361805293, 1.7564232148227155, 3.1826825888829635], 2.5453419339690075, true))),
        ("ellipse", [1.0, 2.0, 3.0000001], 0.0, 0.0, 6.283185307179586, Some((7.49999998772566e-09, [4.0, 2.000000006, 2.9999999955], 9.00000000000001, false))),
        ("ellipse", [1.0, 2.0, 3.0000001], 1.0, 0.0, 6.283185307179586, Some((1.5707963267948966, [1.0000000000000002, 2.8, 2.4], 1.0000001200000095, true))),
        ("ellipse", [1.0, 2.0, 3.0000001], 3.0, 0.0, 6.283185307179586, Some((3.141592646089793, [-2.0, 2.000000006, 2.9999999954999996], 9.00000000000001, false))),
        ("ellipse", [1.0, 2.0, 3.0000001], 5.5, 0.0, 6.283185307179586, Some((4.71238898038469, [0.9999999999999994, 1.2, 3.6], 0.9999998800000105, true))),
        ("ellipse", [2.0, 2.0, 10.0], 0.0, 0.0, 6.283185307179586, Some((2.758439069884682, [-1.7824708728325103, 2.299077793443029, 2.775691654917728], 66.58716449518244, false))),
        ("ellipse", [2.0, 2.0, 10.0], 1.0, 0.0, 6.283185307179586, Some((2.758439069884682, [-1.7824708728325103, 2.299077793443029, 2.775691654917728], 66.58716449518244, false))),
        ("ellipse", [2.0, 2.0, 10.0], 3.0, 0.0, 6.283185307179586, Some((2.758439069884682, [-1.7824708728325103, 2.299077793443029, 2.775691654917728], 66.58716449518244, false))),
        ("ellipse", [2.0, 2.0, 10.0], 5.5, 0.0, 6.283185307179586, Some((4.9581624073934005, [1.7299197477272727, 1.224040451827419, 3.5819696611294356], 41.86616999373091, true))),
        ("hyperbola", [4.0, 2.0, 3.0], -1.0, -2.0, 2.0, Some((-0.6223625037147787, [3.4000000000000004, 1.469340033543136, 3.397994974842648], 0.7999999999999998, true))),
        ("parabola", [4.0, 2.0, 3.0], -1.0, -3.0, 3.0, Some((0.0, [1.0, 2.0, 3.0], 9.0, false))),
        ("hyperbola", [4.0, 2.0, 3.0], 0.0, -2.0, 2.0, Some((0.0, [3.0, 2.0, 3.0], 1.0, true))),
        ("parabola", [4.0, 2.0, 3.0], 0.0, -3.0, 3.0, Some((0.0, [1.0, 2.0, 3.0], 9.0, false))),
        ("hyperbola", [4.0, 2.0, 3.0], 1.5, -2.0, 2.0, Some((0.6223625037147789, [3.4000000000000004, 2.5306599664568643, 2.602005025157352], 0.8, true))),
        ("parabola", [4.0, 2.0, 3.0], 1.5, -3.0, 3.0, Some((2.0, [3.0, 3.6, 1.8], 5.0, true))),
        ("hyperbola", [1.0, 5.0, 1.0], -1.0, -2.0, 2.0, Some((0.6689742267542799, [3.4644674881198982, 2.576, 2.568], 14.408000000000001, true))),
        ("parabola", [1.0, 5.0, 1.0], -1.0, -3.0, 3.0, Some((1.5900334917155636, [2.264103252388594, 3.2720267933724507, 2.045979904970662], 5.677922397124552, true))),
        ("hyperbola", [1.0, 5.0, 1.0], 0.0, -2.0, 2.0, Some((0.6689742267542799, [3.4644674881198982, 2.576, 2.568], 14.408000000000001, true))),
        ("parabola", [1.0, 5.0, 1.0], 0.0, -3.0, 3.0, Some((1.5900334917155636, [2.264103252388594, 3.2720267933724507, 2.045979904970662], 5.677922397124552, true))),
        ("hyperbola", [1.0, 5.0, 1.0], 1.5, -2.0, 2.0, Some((0.6689742267542799, [3.4644674881198982, 2.576, 2.568], 14.408000000000001, true))),
        ("parabola", [1.0, 5.0, 1.0], 1.5, -3.0, 3.0, Some((1.5900334917155636, [2.264103252388594, 3.2720267933724507, 2.045979904970662], 5.677922397124552, true))),
        ("hyperbola", [-3.0, 1.0, 4.0], -1.0, -2.0, 2.0, Some((-0.10786806967403892, [3.0116468069432423, 1.9135381000579395, 3.0648464249565452], 37.8489614006049, true))),
        ("parabola", [-3.0, 1.0, 4.0], -1.0, -3.0, 3.0, Some((-0.27785486795902214, [1.0386016638242628, 1.7777161056327824, 3.166712920775413], 17.60951309640737, true))),
        ("hyperbola", [-3.0, 1.0, 4.0], 0.0, -2.0, 2.0, Some((-0.10786806967403892, [3.0116468069432423, 1.9135381000579395, 3.0648464249565452], 37.8489614006049, true))),
        ("parabola", [-3.0, 1.0, 4.0], 0.0, -3.0, 3.0, Some((-0.27785486795902214, [1.0386016638242628, 1.7777161056327824, 3.166712920775413], 17.60951309640737, true))),
        ("hyperbola", [-3.0, 1.0, 4.0], 1.5, -2.0, 2.0, Some((-0.10786806967403892, [3.0116468069432423, 1.9135381000579395, 3.0648464249565452], 37.8489614006049, true))),
        ("parabola", [-3.0, 1.0, 4.0], 1.5, -3.0, 3.0, Some((-0.27785486795902214, [1.0386016638242628, 1.7777161056327824, 3.166712920775413], 17.60951309640737, true))),
        ("hyperbola", [1.0, 2.0, 3.0000001], -1.0, -2.0, 2.0, Some((-1.1999999843639122e-08, [3.0, 1.9999999904, 3.0000000071999997], 4.000000000000009, true))),
        ("parabola", [1.0, 2.0, 3.0000001], -1.0, -3.0, 3.0, Some((-5.999999990180517e-08, [1.0000000000000018, 1.999999952, 3.000000036], 6.399999993262648e-15, true))),
        ("hyperbola", [1.0, 2.0, 3.0000001], 0.0, -2.0, 2.0, Some((-1.1999999843639122e-08, [3.0, 1.9999999904, 3.0000000071999997], 4.000000000000009, true))),
        ("parabola", [1.0, 2.0, 3.0000001], 0.0, -3.0, 3.0, Some((-5.999999990180517e-08, [1.0000000000000018, 1.999999952, 3.000000036], 6.399999993262648e-15, true))),
        ("hyperbola", [1.0, 2.0, 3.0000001], 1.5, -2.0, 2.0, Some((-1.1999999843639122e-08, [3.0, 1.9999999904, 3.0000000071999997], 4.000000000000009, true))),
        ("parabola", [1.0, 2.0, 3.0000001], 1.5, -3.0, 3.0, Some((-5.999999990180517e-08, [1.0000000000000018, 1.999999952, 3.000000036], 6.399999993262648e-15, true))),
        ("hyperbola", [2.0, 2.0, 10.0], -1.0, -2.0, 2.0, Some((-0.9775120749499221, [4.034081521097536, 1.087364197356642, 3.6844768519825184], 44.85622437588159, true))),
        ("parabola", [2.0, 2.0, 10.0], -1.0, -3.0, 3.0, Some((-2.0327927136297066, [3.0661231082930134, 0.3737658290962347, 4.219675628177824], 37.193405904132845, true))),
        ("hyperbola", [2.0, 2.0, 10.0], 0.0, -2.0, 2.0, Some((-0.9775120749499221, [4.034081521097536, 1.087364197356642, 3.6844768519825184], 44.85622437588159, true))),
        ("parabola", [2.0, 2.0, 10.0], 0.0, -3.0, 3.0, Some((-2.0327927136297066, [3.0661231082930134, 0.3737658290962347, 4.219675628177824], 37.193405904132845, true))),
        ("hyperbola", [2.0, 2.0, 10.0], 1.5, -2.0, 2.0, Some((-0.9775120749499221, [4.034081521097536, 1.087364197356642, 3.6844768519825184], 44.85622437588159, true))),
        ("parabola", [2.0, 2.0, 10.0], 1.5, -3.0, 3.0, Some((-2.0327927136297066, [3.0661231082930134, 0.3737658290962347, 4.219675628177824], 37.193405904132845, true))),
        ("offset", [3.0, 3.0, 0.0], 0.5, 0.0, 4.0, Some((0.431224895553831, [1.0808638128527923, 0.2744276008708022, 0.0], 11.111828607712813, true))),
        ("offset", [3.0, 3.0, 0.0], 1.5, 0.0, 4.0, Some((1.6159059397943798, [3.9333088956620528, 1.69003828240809, 0.0], 2.5870651962782674, true))),
        ("offset", [3.0, 3.0, 0.0], 2.5, 0.0, 4.0, Some((1.6159059397943798, [3.9333088956620528, 1.69003828240809, 0.0], 2.5870651962782674, true))),
        ("offset", [3.0, 3.0, 0.0], 3.5, 0.0, 4.0, None),
        ("offset", [1.0, 0.0, 0.0], 0.5, 0.0, 4.0, Some((0.46846665742363863, [1.1670202415620776, 0.2095692182692181, 0.0], 0.07181501833742594, true))),
        ("offset", [1.0, 0.0, 0.0], 1.5, 0.0, 4.0, Some((2.0287190968310034, [4.264651924004473, 1.2682448796582388, 0.0], 12.266397259685444, true))),
        ("offset", [1.0, 0.0, 0.0], 2.5, 0.0, 4.0, Some((2.0287190968310034, [4.264651924004473, 1.2682448796582388, 0.0], 12.266397259685444, true))),
        ("offset", [1.0, 0.0, 0.0], 3.5, 0.0, 4.0, Some((2.0287190968310034, [4.264651924004473, 1.2682448796582388, 0.0], 12.266397259685444, true))),
        ("offset", [5.0, 2.0, 0.0], 0.5, 0.0, 4.0, Some((1.9455563054237872, [4.144341358388834, 1.5832320651967868, 0.0], 0.9058472224440011, true))),
        ("offset", [5.0, 2.0, 0.0], 1.5, 0.0, 4.0, Some((1.9455563054237872, [4.144341358388834, 1.5832320651967868, 0.0], 0.9058472224440011, true))),
        ("offset", [5.0, 2.0, 0.0], 2.5, 0.0, 4.0, Some((1.9455563054237872, [4.144341358388834, 1.5832320651967868, 0.0], 0.9058472224440011, true))),
        ("offset", [5.0, 2.0, 0.0], 3.5, 0.0, 4.0, None),
        ("offset", [8.0, -1.0, 0.0], 0.5, 0.0, 4.0, None),
        ("offset", [8.0, -1.0, 0.0], 1.5, 0.0, 4.0, Some((2.7916088371554575, [6.328153205738166, 0.9230205079007626, 0.0], 6.493079577290477, true))),
        ("offset", [8.0, -1.0, 0.0], 2.5, 0.0, 4.0, Some((2.7916088371554575, [6.328153205738166, 0.9230205079007626, 0.0], 6.493079577290477, true))),
        ("offset", [8.0, -1.0, 0.0], 3.5, 0.0, 4.0, Some((3.752332308135213, [8.219121390115964, 0.4289619243103159, 0.0], 2.0899463647349936, true))),
        ("offset", [2.5, 0.5, 0.0], 0.5, 0.0, 4.0, Some((1.1203499281404088, [2.96717308041857, 0.2592181588974158, 0.0], 0.27622658207252565, true))),
        ("offset", [2.5, 0.5, 0.0], 1.5, 0.0, 4.0, Some((1.1203499281404088, [2.96717308041857, 0.2592181588974158, 0.0], 0.27622658207252565, true))),
        ("offset", [2.5, 0.5, 0.0], 2.5, 0.0, 4.0, Some((2.0378972008629734, [4.281176502059924, 1.2267381475856918, 0.0], 3.700738066646708, true))),
        ("offset", [2.5, 0.5, 0.0], 3.5, 0.0, 4.0, Some((2.0378972008629734, [4.281176502059924, 1.2267381475856918, 0.0], 3.700738066646708, true))),
        ("offset", [6.0, 5.0, 0.0], 0.5, 0.0, 4.0, Some((0.5621202721085621, [1.3940464914592143, 0.011176934057826104, 0.0], 46.103163306115846, true))),
        ("offset", [6.0, 5.0, 0.0], 1.5, 0.0, 4.0, Some((1.8207888158183472, [4.03389293443508, 1.738192196587849, 0.0], 14.504967139664702, true))),
        ("offset", [6.0, 5.0, 0.0], 2.5, 0.0, 4.0, Some((2.340745425621043, [5.087011148706932, 0.42456205469932173, 0.0], 21.768181033882733, false))),
        ("offset", [6.0, 5.0, 0.0], 3.5, 0.0, 4.0, Some((3.0864065579959945, [6.695665142804487, 1.2822804396101226, 0.0], 14.30538872061869, false))),
        ("offset", [0.0, -2.0, 0.0], 0.5, 0.0, 4.0, Some((0.6480389893789016, [1.6089635069453194, -0.18618811815801822, 0.0], 5.878677109392932, true))),
        ("offset", [0.0, -2.0, 0.0], 1.5, 0.0, 4.0, Some((2.122089799871008, [4.452129461210507, 0.8934513963179848, 0.0], 28.193517722233057, true))),
        ("offset", [0.0, -2.0, 0.0], 2.5, 0.0, 4.0, Some((2.122089799871008, [4.452129461210507, 0.8934513963179848, 0.0], 28.193517722233057, true))),
        ("offset", [0.0, -2.0, 0.0], 3.5, 0.0, 4.0, Some((2.122089799871008, [4.452129461210507, 0.8934513963179848, 0.0], 28.193517722233057, true))),
        ("offset", [9.5, 2.0, 0.0], 0.5, 0.0, 4.0, None),
        ("offset", [9.5, 2.0, 0.0], 1.5, 0.0, 4.0, None),
        ("offset", [9.5, 2.0, 0.0], 2.5, 0.0, 4.0, None),
        ("offset", [9.5, 2.0, 0.0], 3.5, 0.0, 4.0, None),
        // Points off the B-spline on normals at its knots (extrema at interval bounds, found
        // by the sign of (C(u) - P).C'(u) either side). Left out: P = (2.20722997811544,
        // -0.5855400437691198, 3.6777002188455996) from U0 = 3 and 3.5, where Newton on
        // [2, 4] crawls along a minimum of |F| without a root (F ~ 2.669, F' ~ 1e-8) and
        // succeeds or fails on rounding of the curve's evaluation (OCCT 2.15344..., kentron
        // the knot 1).
        ("bspline", [2.20722997811544, -0.5855400437691198, 3.6777002188455996], 0.0, 0.0, 4.0, Some((1.0, [2.5, 0.0, 0.75], 9.0, true))),
        ("bspline", [2.20722997811544, -0.5855400437691198, 3.6777002188455996], 1.0, 0.0, 4.0, Some((1.0, [2.5, 0.0, 0.75], 9.0, true))),
        ("bspline", [2.20722997811544, -0.5855400437691198, 3.6777002188455996], 1.5, 0.0, 4.0, Some((1.0, [2.5, 0.0, 0.75], 9.0, true))),
        ("bspline", [2.4024099927051465, -0.19518001458970663, 1.7259000729485332], 0.0, 0.0, 4.0, Some((1.0, [2.5, 0.0, 0.75], 1.0, false))),
        ("bspline", [2.4024099927051465, -0.19518001458970663, 1.7259000729485332], 1.0, 0.0, 4.0, Some((1.0, [2.5, 0.0, 0.75], 1.0, true))),
        ("bspline", [2.4024099927051465, -0.19518001458970663, 1.7259000729485332], 1.5, 0.0, 4.0, Some((1.0, [2.5, 0.0, 0.75], 1.0, true))),
        ("bspline", [2.4024099927051465, -0.19518001458970663, 1.7259000729485332], 3.0, 0.0, 4.0, Some((1.0, [2.5, 0.0, 0.75], 1.0, true))),
        ("bspline", [2.4024099927051465, -0.19518001458970663, 1.7259000729485332], 3.5, 0.0, 4.0, Some((1.0, [2.5, 0.0, 0.75], 1.0, true))),
        ("bspline", [2.2316718427000253, 0.13416407864998736, 0.75], 0.0, 0.0, 4.0, Some((1.0, [2.5, 0.0, 0.75], 0.08999999999999997, true))),
        ("bspline", [2.2316718427000253, 0.13416407864998736, 0.75], 1.0, 0.0, 4.0, Some((1.0, [2.5, 0.0, 0.75], 0.08999999999999997, true))),
        ("bspline", [2.2316718427000253, 0.13416407864998736, 0.75], 1.5, 0.0, 4.0, Some((1.0, [2.5, 0.0, 0.75], 0.08999999999999997, true))),
        ("bspline", [2.2316718427000253, 0.13416407864998736, 0.75], 3.0, 0.0, 4.0, Some((1.0, [2.5, 0.0, 0.75], 0.08999999999999997, true))),
        ("bspline", [2.2316718427000253, 0.13416407864998736, 0.75], 3.5, 0.0, 4.0, Some((1.0, [2.5, 0.0, 0.75], 0.08999999999999997, true))),
        ("bspline", [2.470722997811544, -0.05855400437691199, 1.04277002188456], 0.0, 0.0, 4.0, Some((1.0, [2.5, 0.0, 0.75], 0.08999999999999998, false))),
        ("bspline", [2.470722997811544, -0.05855400437691199, 1.04277002188456], 1.0, 0.0, 4.0, Some((1.0, [2.5, 0.0, 0.75], 0.08999999999999998, true))),
        ("bspline", [2.470722997811544, -0.05855400437691199, 1.04277002188456], 1.5, 0.0, 4.0, Some((1.0, [2.5, 0.0, 0.75], 0.08999999999999998, true))),
        ("bspline", [2.470722997811544, -0.05855400437691199, 1.04277002188456], 3.0, 0.0, 4.0, Some((1.0, [2.5, 0.0, 0.75], 0.08999999999999998, true))),
        ("bspline", [2.470722997811544, -0.05855400437691199, 1.04277002188456], 3.5, 0.0, 4.0, Some((1.0, [2.5, 0.0, 0.75], 0.08999999999999998, true))),
        ("bspline", [2.7683281572999747, -0.13416407864998736, 0.75], 0.0, 0.0, 4.0, Some((1.0, [2.5, 0.0, 0.75], 0.08999999999999997, true))),
        ("bspline", [2.7683281572999747, -0.13416407864998736, 0.75], 1.0, 0.0, 4.0, Some((1.0, [2.5, 0.0, 0.75], 0.08999999999999997, true))),
        ("bspline", [2.7683281572999747, -0.13416407864998736, 0.75], 1.5, 0.0, 4.0, Some((1.0, [2.5, 0.0, 0.75], 0.08999999999999997, true))),
        ("bspline", [2.7683281572999747, -0.13416407864998736, 0.75], 3.0, 0.0, 4.0, Some((1.0, [2.5, 0.0, 0.75], 0.08999999999999997, true))),
        ("bspline", [2.7683281572999747, -0.13416407864998736, 0.75], 3.5, 0.0, 4.0, Some((1.0, [2.5, 0.0, 0.75], 0.08999999999999997, true))),
        ("bspline", [2.529277002188456, 0.05855400437691199, 0.45722997811544003], 0.0, 0.0, 4.0, Some((1.0, [2.5, 0.0, 0.75], 0.09000000000000002, true))),
        ("bspline", [2.529277002188456, 0.05855400437691199, 0.45722997811544003], 1.0, 0.0, 4.0, Some((1.0, [2.5, 0.0, 0.75], 0.09000000000000002, true))),
        ("bspline", [2.529277002188456, 0.05855400437691199, 0.45722997811544003], 1.5, 0.0, 4.0, Some((1.0, [2.5, 0.0, 0.75], 0.09000000000000002, true))),
        ("bspline", [2.529277002188456, 0.05855400437691199, 0.45722997811544003], 3.0, 0.0, 4.0, Some((1.0, [2.5, 0.0, 0.75], 0.09000000000000002, true))),
        ("bspline", [2.529277002188456, 0.05855400437691199, 0.45722997811544003], 3.5, 0.0, 4.0, Some((1.0, [2.5, 0.0, 0.75], 0.09000000000000002, true))),
        ("bspline", [3.3944271909999157, -0.4472135954999579, 0.75], 0.0, 0.0, 4.0, Some((1.0, [2.5, 0.0, 0.75], 0.9999999999999997, false))),
        ("bspline", [3.3944271909999157, -0.4472135954999579, 0.75], 1.0, 0.0, 4.0, Some((1.0, [2.5, 0.0, 0.75], 0.9999999999999997, true))),
        ("bspline", [3.3944271909999157, -0.4472135954999579, 0.75], 1.5, 0.0, 4.0, Some((1.794836954136244, [4.388828711089101, 1.2669034661168335, -0.3953054841096555], 5.238756336013116, false))),
        ("bspline", [3.3944271909999157, -0.4472135954999579, 0.75], 3.0, 0.0, 4.0, Some((1.794836954136244, [4.388828711089101, 1.2669034661168335, -0.3953054841096555], 5.238756336013116, false))),
        ("bspline", [3.3944271909999157, -0.4472135954999579, 0.75], 3.5, 0.0, 4.0, Some((1.794836954136244, [4.388828711089101, 1.2669034661168335, -0.3953054841096555], 5.238756336013116, false))),
        ("bspline", [2.5975900072948535, 0.19518001458970663, -0.2259000729485332], 0.0, 0.0, 4.0, Some((1.0, [2.5, 0.0, 0.75], 1.0, true))),
        ("bspline", [2.5975900072948535, 0.19518001458970663, -0.2259000729485332], 1.0, 0.0, 4.0, Some((1.0, [2.5, 0.0, 0.75], 1.0, true))),
        ("bspline", [2.5975900072948535, 0.19518001458970663, -0.2259000729485332], 1.5, 0.0, 4.0, Some((1.0, [2.5, 0.0, 0.75], 1.0, true))),
        ("bspline", [2.5975900072948535, 0.19518001458970663, -0.2259000729485332], 3.0, 0.0, 4.0, None),
        ("bspline", [2.5975900072948535, 0.19518001458970663, -0.2259000729485332], 3.5, 0.0, 4.0, None),
        ("bspline", [5.183281572999748, -1.3416407864998738, 0.75], 0.0, 0.0, 4.0, Some((1.0, [2.5, 0.0, 0.75], 9.000000000000004, true))),
        ("bspline", [5.183281572999748, -1.3416407864998738, 0.75], 1.0, 0.0, 4.0, Some((1.0, [2.5, 0.0, 0.75], 9.000000000000004, true))),
        ("bspline", [5.183281572999748, -1.3416407864998738, 0.75], 1.5, 0.0, 4.0, Some((1.5356870766498376, [3.6571110268362608, 1.5456165109428284, 0.2678142778336176], 10.89795430827421, false))),
        ("bspline", [5.183281572999748, -1.3416407864998738, 0.75], 3.0, 0.0, 4.0, Some((1.5356870766498376, [3.6571110268362608, 1.5456165109428284, 0.2678142778336176], 10.89795430827421, false))),
        ("bspline", [5.183281572999748, -1.3416407864998738, 0.75], 3.5, 0.0, 4.0, Some((1.5356870766498376, [3.6571110268362608, 1.5456165109428284, 0.2678142778336176], 10.89795430827421, false))),
        ("bspline", [2.79277002188456, 0.5855400437691198, -2.1777002188455996], 0.0, 0.0, 4.0, Some((0.5458507193634121, [1.556233226260109, 0.2695499716794474, 0.32495111399914073], 7.892136666381977, true))),
        ("bspline", [2.79277002188456, 0.5855400437691198, -2.1777002188455996], 1.0, 0.0, 4.0, Some((1.0, [2.5, 0.0, 0.75], 9.0, false))),
        ("bspline", [2.79277002188456, 0.5855400437691198, -2.1777002188455996], 1.5, 0.0, 4.0, Some((1.0, [2.5, 0.0, 0.75], 9.0, false))),
        ("bspline", [2.79277002188456, 0.5855400437691198, -2.1777002188455996], 3.0, 0.0, 4.0, None),
        ("bspline", [2.79277002188456, 0.5855400437691198, -2.1777002188455996], 3.5, 0.0, 4.0, None),
        ("bspline", [2.878679656440357, 2.121320343559643, -1.0], 0.0, 0.0, 4.0, Some((0.4423918649368796, [1.2838852145192197, 0.4979196024318685, 0.2286302728063667], 6.6883316255330065, true))),
        ("bspline", [2.878679656440357, 2.121320343559643, -1.0], 1.0, 0.0, 4.0, Some((2.0, [5.0, 0.0, -1.0], 9.000000000000002, true))),
        ("bspline", [2.878679656440357, 2.121320343559643, -1.0], 1.5, 0.0, 4.0, Some((1.606152294580062, [3.8490029270437143, 1.5844416424535341, 0.10517800800099208], 2.451184418544819, true))),
        ("bspline", [2.878679656440357, 2.121320343559643, -1.0], 3.0, 0.0, 4.0, Some((2.0, [5.0, 0.0, -1.0], 9.000000000000002, true))),
        ("bspline", [2.878679656440357, 2.121320343559643, -1.0], 3.5, 0.0, 4.0, Some((2.0, [5.0, 0.0, -1.0], 9.000000000000002, true))),
        ("bspline", [3.775255128608411, -1.224744871391589, 1.4494897427831779], 0.0, 0.0, 4.0, Some((1.7245927248642694, [4.184222880822423, 1.4662644430116336, -0.19988994016684525], 10.129239091084203, false))),
        ("bspline", [3.775255128608411, -1.224744871391589, 1.4494897427831779], 1.0, 0.0, 4.0, Some((2.0, [5.0, 0.0, -1.0], 8.999999999999998, true))),
        ("bspline", [3.775255128608411, -1.224744871391589, 1.4494897427831779], 1.5, 0.0, 4.0, Some((1.7245927248642694, [4.184222880822423, 1.4662644430116336, -0.19988994016684525], 10.129239091084203, false))),
        ("bspline", [3.775255128608411, -1.224744871391589, 1.4494897427831779], 3.0, 0.0, 4.0, Some((2.0, [5.0, 0.0, -1.0], 8.999999999999998, true))),
        ("bspline", [3.775255128608411, -1.224744871391589, 1.4494897427831779], 3.5, 0.0, 4.0, Some((2.0, [5.0, 0.0, -1.0], 8.999999999999998, true))),
        ("bspline", [4.292893218813452, 0.7071067811865476, -1.0], 0.0, 0.0, 4.0, None),
        ("bspline", [4.292893218813452, 0.7071067811865476, -1.0], 1.0, 0.0, 4.0, Some((2.0, [5.0, 0.0, -1.0], 1.0000000000000007, false))),
        ("bspline", [4.292893218813452, 0.7071067811865476, -1.0], 1.5, 0.0, 4.0, Some((2.0, [5.0, 0.0, -1.0], 1.0000000000000007, false))),
        ("bspline", [4.292893218813452, 0.7071067811865476, -1.0], 3.0, 0.0, 4.0, Some((2.0, [5.0, 0.0, -1.0], 1.0000000000000007, false))),
        ("bspline", [4.292893218813452, 0.7071067811865476, -1.0], 3.5, 0.0, 4.0, Some((2.0, [5.0, 0.0, -1.0], 1.0000000000000007, false))),
        ("bspline", [4.591751709536137, -0.40824829046386296, -0.18350341907227408], 0.0, 0.0, 4.0, Some((1.3906452534183982, [3.2850664928874456, 1.2720634030985052, 0.5452374081641781], 5.06193683621194, false))),
        ("bspline", [4.591751709536137, -0.40824829046386296, -0.18350341907227408], 1.0, 0.0, 4.0, Some((1.2755850382521468, [3.016833283277753, 0.9290171259242771, 0.6980494639402778], 5.045782328783673, true))),
        ("bspline", [4.591751709536137, -0.40824829046386296, -0.18350341907227408], 1.5, 0.0, 4.0, Some((1.3906452534183982, [3.2850664928874456, 1.2720634030985052, 0.5452374081641781], 5.06193683621194, false))),
        ("bspline", [4.591751709536137, -0.40824829046386296, -0.18350341907227408], 3.0, 0.0, 4.0, Some((2.0, [5.0, 0.0, -1.0], 0.9999999999999999, false))),
        ("bspline", [4.591751709536137, -0.40824829046386296, -0.18350341907227408], 3.5, 0.0, 4.0, Some((2.0, [5.0, 0.0, -1.0], 0.9999999999999999, false))),
        ("bspline", [4.787867965644036, 0.21213203435596426, -1.0], 0.0, 0.0, 4.0, None),
        ("bspline", [4.787867965644036, 0.21213203435596426, -1.0], 1.0, 0.0, 4.0, Some((2.0, [5.0, 0.0, -1.0], 0.08999999999999983, true))),
        ("bspline", [4.787867965644036, 0.21213203435596426, -1.0], 1.5, 0.0, 4.0, Some((2.0, [5.0, 0.0, -1.0], 0.08999999999999983, true))),
        ("bspline", [4.787867965644036, 0.21213203435596426, -1.0], 3.0, 0.0, 4.0, Some((2.0, [5.0, 0.0, -1.0], 0.08999999999999983, true))),
        ("bspline", [4.787867965644036, 0.21213203435596426, -1.0], 3.5, 0.0, 4.0, Some((2.0, [5.0, 0.0, -1.0], 0.08999999999999983, true))),
        ("bspline", [4.877525512860841, -0.12247448713915889, -0.7550510257216823], 0.0, 0.0, 4.0, None),
        ("bspline", [4.877525512860841, -0.12247448713915889, -0.7550510257216823], 1.0, 0.0, 4.0, Some((2.0, [5.0, 0.0, -1.0], 0.08999999999999993, true))),
        ("bspline", [4.877525512860841, -0.12247448713915889, -0.7550510257216823], 1.5, 0.0, 4.0, Some((2.0, [5.0, 0.0, -1.0], 0.08999999999999993, true))),
        ("bspline", [4.877525512860841, -0.12247448713915889, -0.7550510257216823], 3.0, 0.0, 4.0, Some((2.0, [5.0, 0.0, -1.0], 0.08999999999999993, true))),
        ("bspline", [4.877525512860841, -0.12247448713915889, -0.7550510257216823], 3.5, 0.0, 4.0, Some((2.0, [5.0, 0.0, -1.0], 0.08999999999999993, true))),
        ("bspline", [5.212132034355964, -0.21213203435596426, -1.0], 0.0, 0.0, 4.0, None),
        ("bspline", [5.212132034355964, -0.21213203435596426, -1.0], 1.0, 0.0, 4.0, Some((2.0, [5.0, 0.0, -1.0], 0.08999999999999983, false))),
        ("bspline", [5.212132034355964, -0.21213203435596426, -1.0], 1.5, 0.0, 4.0, Some((2.0, [5.0, 0.0, -1.0], 0.08999999999999983, false))),
        ("bspline", [5.212132034355964, -0.21213203435596426, -1.0], 3.0, 0.0, 4.0, Some((2.0, [5.0, 0.0, -1.0], 0.08999999999999983, false))),
        ("bspline", [5.212132034355964, -0.21213203435596426, -1.0], 3.5, 0.0, 4.0, Some((2.0, [5.0, 0.0, -1.0], 0.08999999999999983, false))),
        ("bspline", [5.122474487139159, 0.12247448713915889, -1.2449489742783177], 0.0, 0.0, 4.0, None),
        ("bspline", [5.122474487139159, 0.12247448713915889, -1.2449489742783177], 1.0, 0.0, 4.0, Some((2.0, [5.0, 0.0, -1.0], 0.08999999999999993, false))),
        ("bspline", [5.122474487139159, 0.12247448713915889, -1.2449489742783177], 1.5, 0.0, 4.0, Some((2.0, [5.0, 0.0, -1.0], 0.08999999999999993, false))),
        ("bspline", [5.122474487139159, 0.12247448713915889, -1.2449489742783177], 3.0, 0.0, 4.0, Some((2.0, [5.0, 0.0, -1.0], 0.08999999999999993, false))),
        ("bspline", [5.122474487139159, 0.12247448713915889, -1.2449489742783177], 3.5, 0.0, 4.0, Some((2.0, [5.0, 0.0, -1.0], 0.08999999999999993, false))),
        ("bspline", [5.707106781186548, -0.7071067811865476, -1.0], 0.0, 0.0, 4.0, None),
        ("bspline", [5.707106781186548, -0.7071067811865476, -1.0], 1.0, 0.0, 4.0, Some((2.0, [5.0, 0.0, -1.0], 1.0000000000000007, true))),
        ("bspline", [5.707106781186548, -0.7071067811865476, -1.0], 1.5, 0.0, 4.0, Some((2.0, [5.0, 0.0, -1.0], 1.0000000000000007, true))),
        ("bspline", [5.707106781186548, -0.7071067811865476, -1.0], 3.0, 0.0, 4.0, Some((2.0, [5.0, 0.0, -1.0], 1.0000000000000007, true))),
        ("bspline", [5.707106781186548, -0.7071067811865476, -1.0], 3.5, 0.0, 4.0, Some((2.0, [5.0, 0.0, -1.0], 1.0000000000000007, true))),
        ("bspline", [5.408248290463863, 0.40824829046386296, -1.8164965809277258], 0.0, 0.0, 4.0, None),
        ("bspline", [5.408248290463863, 0.40824829046386296, -1.8164965809277258], 1.0, 0.0, 4.0, Some((2.0, [5.0, 0.0, -1.0], 0.9999999999999997, true))),
        ("bspline", [5.408248290463863, 0.40824829046386296, -1.8164965809277258], 1.5, 0.0, 4.0, Some((2.0, [5.0, 0.0, -1.0], 0.9999999999999997, true))),
        ("bspline", [5.408248290463863, 0.40824829046386296, -1.8164965809277258], 3.0, 0.0, 4.0, Some((2.0, [5.0, 0.0, -1.0], 0.9999999999999997, true))),
        ("bspline", [5.408248290463863, 0.40824829046386296, -1.8164965809277258], 3.5, 0.0, 4.0, Some((2.0, [5.0, 0.0, -1.0], 0.9999999999999997, true))),
        ("bspline", [6.224744871391589, 1.224744871391589, -3.449489742783178], 0.0, 0.0, 4.0, None),
        ("bspline", [6.224744871391589, 1.224744871391589, -3.449489742783178], 1.0, 0.0, 4.0, Some((2.0, [5.0, 0.0, -1.0], 9.0, true))),
        ("bspline", [6.224744871391589, 1.224744871391589, -3.449489742783178], 1.5, 0.0, 4.0, Some((2.0, [5.0, 0.0, -1.0], 9.0, true))),
        ("bspline", [6.224744871391589, 1.224744871391589, -3.449489742783178], 3.0, 0.0, 4.0, Some((2.0, [5.0, 0.0, -1.0], 9.0, true))),
        ("bspline", [6.224744871391589, 1.224744871391589, -3.449489742783178], 3.5, 0.0, 4.0, Some((2.0, [5.0, 0.0, -1.0], 9.0, true))),
    ];

    #[test]
    fn locate_ext_pc_matches_occt() {
        let mut failures = Vec::new();
        for &(name, p, u0, umin, usup, expected) in OCCT_CASES {
            let curve = adaptor(name);
            let found = locate_ext_pc(Point::new(p[0], p[1], p[2]), &curve, u0, umin, usup, 1.0e-9);
            let matches = match (found, expected) {
                (None, None) => true,
                (Some(found), Some((parameter, point, square_distance, is_min))) => {
                    (found.parameter - parameter).abs() <= 1.0e-12 * parameter.abs().max(1.0)
                        && found.point.distance(Point::new(point[0], point[1], point[2])) <= 1.0e-12
                        && (found.square_distance - square_distance).abs() <= 1.0e-12 * square_distance.max(1.0)
                        && found.is_min == is_min
                }
                _ => false,
            };
            if !matches {
                failures.push(format!("{name} {p:?} u0={u0} [{umin}, {usup}]: {found:?}, OCCT {expected:?}"));
            }
        }
        assert!(failures.is_empty(), "{} of {} differ:\n{}", failures.len(), OCCT_CASES.len(), failures.join("\n"));
    }
}
