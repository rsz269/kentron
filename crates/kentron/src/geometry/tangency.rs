//! Circles satisfying tangency constraints in the plane.

use std::f64::consts::TAU;

use crate::error::{Error, Result};
use crate::math::{Direction2, Frame2, Point2, Vector2};

use super::curve2::{Circle2, Curve2, Line2, OffsetCurve2};
use super::{extrema, intersect, ParametricCurve};

/// Which side of a curve a tangent circle must lie on.
#[doc(alias = "GccEnt_Position")]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub enum Qualifier {
    /// Any side.
    #[default]
    Unqualified,
    /// The solution encloses the curve.
    Enclosing,
    /// The solution is enclosed by the curve.
    Enclosed,
    /// The solution and the curve are outside each other.
    Outside,
}

/// Something a circle must touch.
#[doc(alias = "Geom2dGcc_QualifiedCurve")]
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Tangent<'a> {
    Curve(&'a Curve2, Qualifier),
    /// The circle passes through the point.
    Point(Point2),
}

/// Where a solution circle touches one of its arguments.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Tangency {
    pub point: Point2,
    pub parameter_on_circle: f64,
    /// Parameter on the argument curve (0 for a point).
    pub parameter_on_argument: f64,
}

/// A solution: the circle and its tangency with each argument, in argument order.
#[derive(Debug, Clone, PartialEq)]
pub struct TangentCircle {
    pub circle: Circle2,
    pub tangencies: Vec<Tangency>,
}

/// The circles of radius `radius` tangent to `a` and `b`. An empty result means there is
/// no solution; an error means the computation failed.
#[doc(alias = "Geom2dGcc_Circ2d2TanRad")]
#[doc(alias = "GccAna_Circ2d2TanRad")]
pub fn with_radius(a: Tangent, b: Tangent, radius: f64, tolerance: f64) -> Result<Vec<TangentCircle>> {
    if !(radius > 0.0) || !radius.is_finite() {
        return Err(Error::InvalidInput(format!("tangent circle radius {radius}")));
    }
    // The centers lie on the curves at distance `radius` from each argument.
    let mut solutions: Vec<TangentCircle> = Vec::new();
    for ca in center_loci(a, radius) {
        for cb in center_loci(b, radius) {
            for x in intersect::curves2(&ca, &cb, tolerance) {
                let center = x.point;
                if solutions.iter().any(|s| s.circle.center().distance(center) <= tolerance) {
                    continue;
                }
                let circle = Circle2::new(Frame2::new(center, Direction2::X, true), radius);
                let tangencies = vec![touch(&circle, a, x.t1), touch(&circle, b, x.t2)];
                if qualified(&circle, a, &tangencies[0]) && qualified(&circle, b, &tangencies[1]) {
                    solutions.push(TangentCircle { circle, tangencies });
                }
            }
        }
    }
    Ok(solutions)
}

/// The circles centered at `center` tangent to `a`.
#[doc(alias = "Geom2dGcc_Circ2dTanCen")]
pub fn with_center(a: Tangent, center: Point2, tolerance: f64) -> Result<Vec<TangentCircle>> {
    let feet: Vec<(f64, Point2)> = match a {
        Tangent::Point(p) => vec![(0.0, p)],
        Tangent::Curve(curve, _) => extrema::point_curve2(center, curve)
            .into_iter()
            .filter(|x| {
                // Only true feet of perpendiculars, not domain ends.
                let (p, d) = curve.d1(x.parameter);
                let r = p - center;
                r.dot(d).abs() <= 1e-7 * r.length().max(tolerance) * d.length()
            })
            .map(|x| (x.parameter, x.point))
            .collect(),
    };
    let mut solutions: Vec<TangentCircle> = Vec::new();
    for (t, foot) in feet {
        let radius = foot.distance(center);
        if radius <= tolerance || solutions.iter().any(|s| (s.circle.radius() - radius).abs() <= tolerance) {
            continue;
        }
        let circle = Circle2::new(Frame2::new(center, Direction2::X, true), radius);
        let tangency = touch(&circle, a, t);
        if qualified(&circle, a, &tangency) {
            solutions.push(TangentCircle { circle, tangencies: vec![tangency] });
        }
    }
    Ok(solutions)
}

/// The circles tangent to `a` and `b` whose center lies on `on`, found iteratively from
/// initial parameters on `a`, `b` and `on`.
#[doc(alias = "Geom2dGcc_Circ2d2TanOn")]
pub fn centered_on(
    a: Tangent,
    b: Tangent,
    on: &Curve2,
    tolerance: f64,
    initial: (f64, f64, f64),
) -> Result<Vec<TangentCircle>> {
    // Unknowns: the parameters on a, b and `on`. Equations: the center is on the normal of
    // each curve argument, and equally far from both.
    let foot = |arg: Tangent, t: f64| match arg {
        Tangent::Point(p) => (p, Vector2::new(0.0, 0.0)),
        Tangent::Curve(c, _) => c.d1(t),
    };
    let residual = |x: [f64; 3]| -> [f64; 3] {
        let center = on.point(x[2]);
        let (pa, da) = foot(a, x[0]);
        let (pb, db) = foot(b, x[1]);
        let (ra, rb) = (center - pa, center - pb);
        let unit = |d: Vector2| if d.length() > 0.0 { d * (1.0 / d.length()) } else { d };
        // A point argument fixes its parameter.
        let e0 = if matches!(a, Tangent::Point(_)) { x[0] - initial.0 } else { ra.dot(unit(da)) };
        let e1 = if matches!(b, Tangent::Point(_)) { x[1] - initial.1 } else { rb.dot(unit(db)) };
        [e0, e1, ra.length() - rb.length()]
    };
    let size = |r: [f64; 3]| r.iter().map(|v| v.abs()).fold(0.0, f64::max);
    let mut x = [initial.0, initial.1, initial.2];
    let mut r = residual(x);
    for _ in 0..100 {
        if size(r) <= tolerance {
            break;
        }
        // A finite-difference Jacobian.
        let mut jac = [[0.0; 3]; 3];
        for k in 0..3 {
            let h = 1e-7 * (1.0 + x[k].abs());
            let mut y = x;
            y[k] += h;
            let ry = residual(y);
            for i in 0..3 {
                jac[i][k] = (ry[i] - r[i]) / h;
            }
        }
        let Some(step) = solve3(jac, [-r[0], -r[1], -r[2]]) else { break };
        // Halve the step until the residual drops.
        let mut lambda = 1.0;
        let mut improved = false;
        for _ in 0..30 {
            let y = [x[0] + lambda * step[0], x[1] + lambda * step[1], x[2] + lambda * step[2]];
            let ry = residual(y);
            if ry.iter().all(|v| v.is_finite()) && size(ry) < size(r) {
                (x, r, improved) = (y, ry, true);
                break;
            }
            lambda *= 0.5;
        }
        if !improved {
            break;
        }
    }
    if !(size(r) <= tolerance.max(1e-9)) {
        return Ok(Vec::new());
    }
    let center = on.point(x[2]);
    let radius = (center - foot(a, x[0]).0).length();
    if !(radius > tolerance) {
        return Ok(Vec::new());
    }
    let circle = Circle2::new(Frame2::new(center, Direction2::X, true), radius);
    let tangencies = vec![touch(&circle, a, x[0]), touch(&circle, b, x[1])];
    if !(qualified(&circle, a, &tangencies[0]) && qualified(&circle, b, &tangencies[1])) {
        return Ok(Vec::new());
    }
    Ok(vec![TangentCircle { circle, tangencies }])
}

/// The curves the center of a circle of radius `radius` tangent to `arg` lies on: offsets
/// to the sides the qualifier allows, or the circle round a point.
fn center_loci(arg: Tangent, radius: f64) -> Vec<Curve2> {
    match arg {
        Tangent::Point(p) => vec![Curve2::Circle(Circle2::new(Frame2::new(p, Direction2::X, true), radius))],
        Tangent::Curve(curve, qualifier) => {
            // Offsets go to the right; the inside of a curve is on its left.
            let sides: &[f64] = match qualifier {
                Qualifier::Unqualified => &[1.0, -1.0],
                Qualifier::Outside => &[1.0],
                Qualifier::Enclosed | Qualifier::Enclosing => &[-1.0],
            };
            sides.iter().map(|&s| offset_of(curve, s * radius)).collect()
        }
    }
}

/// `curve` offset by `distance` to the right, exactly where the offset is a conic.
fn offset_of(curve: &Curve2, distance: f64) -> Curve2 {
    match curve {
        Curve2::Line(line) => {
            let axis = line.axis();
            let d = axis.direction;
            Curve2::Line(Line2::new(axis.origin + distance * Vector2::new(d.y(), -d.x()), d))
        }
        _ => Curve2::Offset(OffsetCurve2::new(curve.clone(), distance)),
    }
}

/// Where `circle` touches `arg`, the touching point being the foot on the argument.
fn touch(circle: &Circle2, arg: Tangent, t: f64) -> Tangency {
    let (point, parameter_on_argument) = match arg {
        Tangent::Point(p) => (p, 0.0),
        Tangent::Curve(c, _) => (c.point(t), t),
    };
    let v = point - circle.center();
    let parameter_on_circle = v.y.atan2(v.x).rem_euclid(TAU);
    Tangency { point, parameter_on_circle, parameter_on_argument }
}

/// Whether the solution lies where the argument's qualifier asks.
fn qualified(circle: &Circle2, arg: Tangent, tangency: &Tangency) -> bool {
    let Tangent::Curve(curve, qualifier) = arg else { return true };
    if qualifier == Qualifier::Unqualified {
        return true;
    }
    let (_, d, dd) = curve.d2(tangency.parameter_on_argument);
    let left = Vector2::new(-d.y, d.x);
    // The center's side: left is inside.
    let inside = (circle.center() - tangency.point).dot(left) > 0.0;
    // Curvature towards the left, and the radius it bends with.
    let speed = d.length();
    let bend = if speed > 0.0 { d.cross(dd) / (speed * speed * speed) } else { 0.0 };
    match qualifier {
        Qualifier::Outside => !inside,
        // Enclosed: the circle bends at least as much as the curve; enclosing: at most.
        Qualifier::Enclosed => inside && (bend <= 0.0 || circle.radius() <= 1.0 / bend * (1.0 + 1e-9)),
        Qualifier::Enclosing => inside && bend > 0.0 && circle.radius() >= 1.0 / bend * (1.0 - 1e-9),
        Qualifier::Unqualified => true,
    }
}

/// The solution of a 3×3 linear system, by Cramer's rule.
fn solve3(m: [[f64; 3]; 3], b: [f64; 3]) -> Option<[f64; 3]> {
    let det = |m: [[f64; 3]; 3]| {
        m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
    };
    let d = det(m);
    if !(d.abs() > 1e-300) || !d.is_finite() {
        return None;
    }
    let mut x = [0.0; 3];
    for (k, xk) in x.iter_mut().enumerate() {
        let mut mk = m;
        for i in 0..3 {
            mk[i][k] = b[i];
        }
        *xk = det(mk) / d;
    }
    Some(x)
}
