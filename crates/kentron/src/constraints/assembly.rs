//! Assembly constraints: the placements of parts that best satisfy constraints between
//! points, axes, planes and lines on them (CadQuery's `Assembly.solve`).
//!
//! The problem is CadQuery's, solved as IPOPT (through casadi) solves it, by Newton's method
//! with IPOPT's regularization, line search and stopping rules ([`newton`]), so that the
//! parts end where CadQuery's do also where constraints leave them free to move (as a part
//! constrained only by an axis and a point can still turn about the axis). Each part not
//! locked has six variables, a move of its translation
//! (in units of `scale`, the assembly's size) and of its rotation, from where it starts.
//! Rotations are in CadQuery's three parameters R (modified Rodrigues parameters, the axis
//! times tan(angle / 4)), so the rotation's quaternion is ((1 − |R|²), 2R) / (1 + |R|²).
//! Each constraint's cost is CadQuery's, written as a sum of squares of residuals (with its
//! exact Hessian, from second-order duals); their sum is made least from CadQuery's start (translations unmoved, rotations moved by 0.01 about
//! each axis), with its penalty of 1e-16 times the squared moves, so that underdetermined
//! problems settle where CadQuery's do.

use super::dual::Jet;
use super::newton::{self, Evaluation, Settings, Stop};
use crate::math::{Point, Quaternion, Vector};

/// A part's placement: its translation and rotation, the rotation as CadQuery's three
/// parameters (the rotation axis times tan(angle / 4)).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Placement {
    pub translation: Vector,
    pub rotation: Vector,
}

impl Placement {
    /// The placement rotating by `rotation` (a unit quaternion) and then translating.
    pub fn new(translation: Vector, rotation: Quaternion) -> Self {
        // (The parameters of q and of −q are the same rotation; CadQuery takes q as given.)
        let scale = 1.0 / (1.0 + rotation.w);
        Placement { translation, rotation: Vector::new(rotation.x * scale, rotation.y * scale, rotation.z * scale) }
    }

    /// The rotation as a unit quaternion.
    pub fn quaternion(&self) -> Quaternion {
        let r = self.rotation;
        let m = r.x * r.x + r.y * r.y + r.z * r.z;
        let k = 2.0 / (1.0 + m);
        Quaternion::new(r.x * k, r.y * k, r.z * k, (1.0 - m) / (1.0 + m))
    }
}

/// How two axes are to lie (CadQuery's `Axis` parameter).
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum AxisAngle {
    /// Pointing the same way.
    Same,
    /// Pointing opposite ways (CadQuery's default).
    Opposite,
    /// At this angle (radians) to one another.
    Angle(f64),
}

/// A constraint on the parts numbered as in [`solve`]'s `starts`, its points and directions
/// in each part's own coordinates (where the part's placement is the identity).
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Constraint {
    /// Points on two parts at `distance` apart (0: together).
    Point { first: usize, on_first: Point, second: usize, on_second: Point, distance: f64 },
    /// Directions on two parts at an angle.
    Axis { first: usize, on_first: Vector, second: usize, on_second: Vector, angle: AxisAngle },
    /// A point on the first part in the second's plane (through `origin`, of unit normal
    /// `normal`), moved `offset` along the normal.
    PointInPlane { first: usize, point: Point, second: usize, origin: Point, normal: Vector, offset: f64 },
    /// A point on the first part on the second's line (through `origin` along the unit
    /// `direction`); `distance` nonzero: CadQuery's cost |e|² − distance, e the point's
    /// offset from the line in units of `scale` (distance not squared, as CadQuery has it).
    PointOnLine { first: usize, point: Point, second: usize, origin: Point, direction: Vector, distance: f64 },
    /// A point on a part at a place.
    FixedPoint { part: usize, point: Point, target: Point },
    /// A direction on a part along a direction.
    FixedAxis { part: usize, along: Vector, target: Vector },
    /// A part's rotation (a unit quaternion).
    FixedRotation { part: usize, rotation: Quaternion },
}

/// What [`solve`] found.
#[derive(Debug, Clone)]
pub struct Solution {
    pub placements: Vec<Placement>,
    /// The sum of the constraints' costs there (0 where all hold).
    pub cost: f64,
    /// Newton iterations taken and objective evaluations made.
    pub steps: usize,
    pub evaluations: usize,
    /// Whether the minimization converged: to the tolerance (IPOPT's "Solve_Succeeded"), or to
    /// the acceptable level where the objective stopped changing ("Solved_To_Acceptable_Level").
    pub converged: bool,
    /// Whether to the tolerance (rather than to the acceptable level).
    pub optimal: bool,
    /// Each iteration's Hessian regularization (IPOPT's "rg").
    pub regularizations: Vec<f64>,
}

/// The variables of one constraint: the two parts' six each.
const N: usize = 12;
type D = Jet<N>;

/// The placements of the parts starting at `starts` (those `locked` kept there) that best
/// satisfy `constraints`; `scale` the assembly's size (CadQuery: its bounding box diagonal),
/// in which distances are measured.
pub fn solve(starts: &[Placement], locked: &[bool], constraints: &[Constraint], scale: f64) -> Solution {
    // Each free part's variables: its translation's move / scale, then its rotation's move.
    let mut first_variable = vec![None; starts.len()];
    let mut count = 0;
    for (i, &l) in locked.iter().enumerate() {
        if !l {
            first_variable[i] = Some(count);
            count += 6;
        }
    }
    let mut start = vec![0.0; count];
    for k in first_variable.iter().flatten() {
        start[k + 3..k + 6].fill(1e-2);
    }
    let evaluate = |x: &[f64]| {
        // CadQuery's objective, the sum of the residuals' squares, with its gradient and
        // Hessian: 2 Σ r ∇r and 2 Σ (∇r ∇rᵀ + r ∇²r).
        let mut objective = 0.0;
        let mut gradient = vec![0.0; count];
        let mut hessian = vec![0.0; count * count];
        // A part's translation and rotation as duals over the constraint's variables
        // (`slot` 0 or 1: which of its parts).
        let part = |i: usize, slot: usize| -> ([D; 3], [D; 3]) {
            let s = starts[i];
            let (t0, r0) = ([s.translation.x, s.translation.y, s.translation.z], [s.rotation.x, s.rotation.y, s.rotation.z]);
            match first_variable[i] {
                Some(k) => (
                    std::array::from_fn(|j| D::variable(x[k + j], 6 * slot + j) * scale + t0[j]),
                    std::array::from_fn(|j| D::variable(x[k + 3 + j], 6 * slot + 3 + j) + r0[j]),
                ),
                None => (std::array::from_fn(|j| D::constant(t0[j])), std::array::from_fn(|j| D::constant(r0[j]))),
            }
        };
        let mut push = |values: &[D], parts: &[usize]| {
            // The constraint's variables (local index, global index).
            let mut map = Vec::with_capacity(N);
            for (slot, &p) in parts.iter().enumerate() {
                if let Some(k) = first_variable[p] {
                    for j in 0..6 {
                        map.push((6 * slot + j, k + j));
                    }
                }
            }
            for v in values {
                objective += v.value * v.value;
                for &(a, ga) in &map {
                    gradient[ga] += 2.0 * v.value * v.gradient[a];
                    for &(b, gb) in &map {
                        hessian[ga * count + gb] += 2.0 * (v.gradient[a] * v.gradient[b] + v.value * v.hessian[a][b]);
                    }
                }
            }
        };
        for c in constraints {
            match *c {
                Constraint::Point { first, on_first, second, on_second, distance } => {
                    let ((t1, r1), (t2, r2)) = (part(first, 0), part(second, 1));
                    let d = sub(transform(on_first, t1, r1), transform(on_second, t2, r2)).map(|v| v / scale);
                    if distance == 0.0 {
                        push(&d, &[first, second]);
                    } else {
                        push(&[dot(d, d) - (distance / scale).powi(2)], &[first, second]);
                    }
                }
                Constraint::Axis { first, on_first, second, on_second, angle } => {
                    let ((_, r1), (_, r2)) = (part(first, 0), part(second, 1));
                    let (d1, d2) = (rotate(array(on_first), r1), rotate(array(on_second), r2));
                    match angle {
                        AxisAngle::Same => push(&sub(d1, d2), &[first, second]),
                        AxisAngle::Opposite => push(&add(d1, d2), &[first, second]),
                        AxisAngle::Angle(a) => push(&[dot(d1, d2) - a.cos()], &[first, second]),
                    }
                }
                Constraint::PointInPlane { first, point, second, origin, normal, offset } => {
                    let ((t1, r1), (t2, r2)) = (part(first, 0), part(second, 1));
                    let through = origin + normal * offset;
                    let gap = sub(transform(through, t2, r2), transform(point, t1, r1));
                    push(&[dot(rotate(array(normal), r2), gap) / scale], &[first, second]);
                }
                Constraint::PointOnLine { first, point, second, origin, direction, distance } => {
                    let ((t1, r1), (t2, r2)) = (part(first, 0), part(second, 1));
                    let d = sub(transform(point, t1, r1), transform(origin, t2, r2));
                    let n = rotate(array(direction), r2);
                    let along = dot(d, n);
                    let e: [D; 3] = std::array::from_fn(|j| (d[j] - n[j] * along) / scale);
                    if distance == 0.0 {
                        push(&e, &[first, second]);
                    } else {
                        push(&[dot(e, e) - distance], &[first, second]);
                    }
                }
                Constraint::FixedPoint { part: p, point, target } => {
                    let (t, r) = part(p, 0);
                    let at = transform(point, t, r);
                    push(&[(at[0] - target.x) / scale, (at[1] - target.y) / scale, (at[2] - target.z) / scale], &[p]);
                }
                Constraint::FixedAxis { part: p, along, target } => {
                    let (_, r) = part(p, 0);
                    let target = target * (1.0 / target.length());
                    let d = rotate(array(along), r);
                    push(&[d[0] - target.x, d[1] - target.y, d[2] - target.z], &[p]);
                }
                Constraint::FixedRotation { part: p, rotation } => {
                    // CadQuery's cost 1 − (q·q₀)², for unit quaternions the sum of the
                    // squares of their 2×2 minors (Lagrange's identity).
                    let (_, r) = part(p, 0);
                    let (s, u) = quaternion(r);
                    let q = [s, u[0], u[1], u[2]];
                    let q0 = [rotation.w, rotation.x, rotation.y, rotation.z];
                    let mut minors = Vec::with_capacity(6);
                    for i in 0..4 {
                        for j in i + 1..4 {
                            minors.push(q[i] * q0[j] - q[j] * q0[i]);
                        }
                    }
                    push(&minors, &[p]);
                }
            }
        }
        // CadQuery's penalty, 1e-16 times the squared moves.
        for (i, &v) in x.iter().enumerate() {
            objective += 1e-16 * v * v;
            gradient[i] += 2e-16 * v;
            hessian[i * count + i] += 2e-16;
        }
        Evaluation { objective, gradient, hessian }
    };
    // (IPOPT's limit, 3000 iterations, and CadQuery's tolerances.)
    let settings = Settings { tolerance: 1e-14, acceptable: 1e-6, acceptable_change: 1e-12, max_iterations: 3000 };
    let out = newton::minimize(&start, evaluate, settings);
    let placements = starts
        .iter()
        .enumerate()
        .map(|(i, s)| match first_variable[i] {
            Some(k) => {
                let x = &out.variables[k..k + 6];
                Placement {
                    translation: s.translation + Vector::new(x[0], x[1], x[2]) * scale,
                    rotation: s.rotation + Vector::new(x[3], x[4], x[5]),
                }
            }
            None => *s,
        })
        .collect();
    // (The cost without the penalty.)
    let penalty: f64 = out.variables.iter().map(|v| 1e-16 * v * v).sum();
    Solution {
        placements,
        cost: (out.objective - penalty).max(0.0),
        steps: out.iterations,
        evaluations: out.evaluations,
        converged: matches!(out.stop, Stop::Optimal | Stop::Acceptable),
        optimal: out.stop == Stop::Optimal,
        regularizations: out.regularizations,
    }
}

fn add(a: [D; 3], b: [D; 3]) -> [D; 3] {
    std::array::from_fn(|i| a[i] + b[i])
}

fn sub(a: [D; 3], b: [D; 3]) -> [D; 3] {
    std::array::from_fn(|i| a[i] - b[i])
}

fn dot(a: [D; 3], b: [D; 3]) -> D {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// The quaternion (s, u) of rotation parameters `r`.
fn quaternion(r: [D; 3]) -> (D, [D; 3]) {
    let m = dot(r, r);
    let inverse = D::constant(1.0) / (m + 1.0);
    (-(m - 1.0) * inverse, r.map(|v| v * inverse * 2.0))
}

/// `v` rotated by the rotation of parameters `r`: 2(u·v)u + (s² − u·u)v + 2s(u × v).
fn rotate(v: [f64; 3], r: [D; 3]) -> [D; 3] {
    let (s, u) = quaternion(r);
    let uv = u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
    let k = s * s - dot(u, u);
    let cross = [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]];
    std::array::from_fn(|i| uv * u[i] * 2.0 + k * v[i] + s * cross[i] * 2.0)
}

/// `p` rotated by the rotation of parameters `r` and moved by `t`.
fn transform(p: Point, t: [D; 3], r: [D; 3]) -> [D; 3] {
    let rotated = rotate([p.x, p.y, p.z], r);
    std::array::from_fn(|i| rotated[i] + t[i])
}

fn array(v: Vector) -> [f64; 3] {
    [v.x, v.y, v.z]
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    fn unmoved() -> Placement {
        Placement::new(Vector::new(0.0, 0.0, 0.0), Quaternion::IDENTITY)
    }

    /// A box's top face centre (0, 0, 5) mated to another's bottom face centre (0, 0, −5),
    /// their normals opposed (CadQuery's "Plane"): the second, starting at the origin, ends
    /// 10 up, unturned.
    #[test]
    fn a_box_stacked_on_another() {
        let constraints = [
            Constraint::Axis { first: 0, on_first: Vector::new(0.0, 0.0, 1.0), second: 1, on_second: Vector::new(0.0, 0.0, -1.0), angle: AxisAngle::Opposite },
            Constraint::Point { first: 0, on_first: Point::new(0.0, 0.0, 5.0), second: 1, on_second: Point::new(0.0, 0.0, -5.0), distance: 0.0 },
        ];
        let s = solve(&[unmoved(), unmoved()], &[true, false], &constraints, 17.32);
        let p = s.placements[1];
        // (To IPOPT's acceptable level, as CadQuery's solves end: within about 1e-8.)
        assert!(s.converged && s.cost < 1e-15, "{s:?}");
        assert!((p.translation - Vector::new(0.0, 0.0, 10.0)).length() < 1e-7, "{p:?}");
        let q = p.quaternion();
        // (Free to turn about Z: it keeps the start's 0.01 there, as CadQuery's does.)
        assert!(q.x.abs() < 1e-9 && q.y.abs() < 1e-9, "{q:?}");
    }

    /// A part turned to a given rotation (a quarter turn about Z) and a point of it put at
    /// a place: it ends there, turned so.
    #[test]
    fn a_part_turned_and_placed() {
        let turn = Quaternion::from_axis_angle(Vector::new(0.0, 0.0, 1.0), PI / 2.0);
        let constraints = [
            Constraint::FixedRotation { part: 0, rotation: turn },
            Constraint::FixedPoint { part: 0, point: Point::new(1.0, 0.0, 0.0), target: Point::new(3.0, 4.0, 5.0) },
        ];
        let s = solve(&[unmoved()], &[false], &constraints, 10.0);
        let p = s.placements[0];
        // (1, 0, 0) turned a quarter about Z is (0, 1, 0): the translation is (3, 3, 5).
        assert!((p.translation - Vector::new(3.0, 3.0, 5.0)).length() < 1e-7, "{s:?}");
        let q = p.quaternion();
        assert!((q.w * turn.w + q.x * turn.x + q.y * turn.y + q.z * turn.z).abs() > 1.0 - 1e-10, "{q:?}");
    }

    /// The rotation parameters of a quaternion give the quaternion back.
    #[test]
    fn rotation_parameters_round_trip() {
        let q = Quaternion::from_axis_angle(Vector::new(1.0, 2.0, 3.0).normalized(), 2.0);
        let back = Placement::new(Vector::new(0.0, 0.0, 0.0), q).quaternion();
        assert!((back.x - q.x).abs() < 1e-15 && (back.w - q.w).abs() < 1e-15, "{back:?} {q:?}");
    }
}
