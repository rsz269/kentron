//! Sketch constraints: the lines and arcs of a planar sketch moved to best satisfy
//! constraints on them (CadQuery's `Sketch.solve`).
//!
//! The problem is CadQuery's, solved by least squares instead of its SLSQP solver (nlopt):
//! a line's variables are its start and end, (x₁, y₁, x₂, y₂); an arc's its centre, radius,
//! start angle and sweep, (x, y, r, a, da), its point at angle t being
//! (x + r sin t, y + r cos t) (angles from +y towards +x, as CadQuery has them), the radius
//! kept at or above 0. Each constraint's residual is CadQuery's cost, their squares summed
//! and made least from where the sketch's edges are.

use super::dual::Dual;
use super::least_squares::{self, Evaluation, Settings, Stop};
use crate::error::{Error, Result};

/// A sketch edge and its variables.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Entity {
    /// From `start` to `end`.
    Line { start: [f64; 2], end: [f64; 2] },
    /// About `centre`, of `radius`, from angle `start` through `sweep` (radians, from +y
    /// towards +x).
    Arc { centre: [f64; 2], radius: f64, start: f64, sweep: f64 },
}

impl Entity {
    fn variables(&self) -> Vec<f64> {
        match *self {
            Entity::Line { start, end } => vec![start[0], start[1], end[0], end[1]],
            Entity::Arc { centre, radius, start, sweep } => vec![centre[0], centre[1], radius, start, sweep],
        }
    }

    fn is_arc(&self) -> bool {
        matches!(self, Entity::Arc { .. })
    }

    fn from_variables(&self, x: &[f64]) -> Entity {
        match self {
            Entity::Line { .. } => Entity::Line { start: [x[0], x[1]], end: [x[2], x[3]] },
            Entity::Arc { .. } => Entity::Arc { centre: [x[0], x[1]], radius: x[2], start: x[3], sweep: x[4] },
        }
    }
}

/// A constraint on the entities numbered as in [`solve`]'s `entities`. A point "at" a
/// fraction is the point that far along the edge (by parameter); `None` is an arc's centre.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Constraint {
    /// The entity stays where it starts.
    Fixed(usize),
    /// Its point at the fraction stays where it starts.
    FixedPoint(usize, Option<f64>),
    /// The first's end at the second's start.
    Coincident(usize, usize),
    /// The second's direction at its start turned to the first's at its end by this angle
    /// (radians, anticlockwise positive).
    Angle(usize, usize, f64),
    /// A line's length, or an arc's (|r · da|).
    Length(usize, f64),
    /// The first's point at a fraction this far from the second's at a fraction.
    Distance(usize, Option<f64>, usize, Option<f64>, f64),
    /// An arc's radius.
    Radius(usize, f64),
    /// A line's direction along this vector (the angle from it to the vector made 0).
    Orientation(usize, [f64; 2]),
    /// An arc's sweep (radians).
    ArcAngle(usize, f64),
}

/// What [`solve`] found.
#[derive(Debug, Clone)]
pub struct Solution {
    pub entities: Vec<Entity>,
    /// The sum of the squares of the constraints' residuals there.
    pub cost: f64,
    pub steps: usize,
    pub evaluations: usize,
    /// Whether the minimization converged (not stopped for taking too many steps).
    pub converged: bool,
}

/// Two entities' variables (five each, at most).
const N: usize = 10;
type D = Dual<N>;

/// `entities` moved to best satisfy `constraints`.
pub fn solve(entities: &[Entity], constraints: &[Constraint]) -> Result<Solution> {
    let mut first = Vec::with_capacity(entities.len());
    let mut start = Vec::new();
    for e in entities {
        first.push(start.len());
        start.extend(e.variables());
    }
    let mut lower = vec![f64::NEG_INFINITY; start.len()];
    let upper = vec![f64::INFINITY; start.len()];
    for (e, &k) in entities.iter().zip(&first) {
        if e.is_arc() {
            lower[k + 2] = 0.0;
        }
    }
    // (Points at no fraction are arcs' centres: a line has none.)
    for c in constraints {
        let lines_without_fraction = match *c {
            Constraint::FixedPoint(e, None) => !entities[e].is_arc(),
            Constraint::Distance(a, fa, b, fb, _) => (fa.is_none() && !entities[a].is_arc()) || (fb.is_none() && !entities[b].is_arc()),
            _ => false,
        };
        if lines_without_fraction {
            return Err(Error::InvalidInput("a point on a line needs a fraction along it".into()));
        }
        let wrong_kind = match *c {
            Constraint::Radius(e, _) | Constraint::ArcAngle(e, _) => !entities[e].is_arc(),
            Constraint::Orientation(e, _) => entities[e].is_arc(),
            _ => false,
        };
        if wrong_kind {
            return Err(Error::InvalidInput("a constraint on the wrong kind of edge".into()));
        }
    }
    let evaluate = |x: &[f64]| {
        let mut residuals = Vec::new();
        let mut rows = Vec::new();
        // An entity's variables as duals (`slot` 0 or 1: which of the constraint's).
        let vars = |e: usize, slot: usize| -> Vec<D> {
            let k = first[e];
            let n = if entities[e].is_arc() { 5 } else { 4 };
            (0..n).map(|j| D::variable(x[k + j], 5 * slot + j)).collect()
        };
        let mut push = |values: &[D], parts: &[usize]| {
            for v in values {
                let mut row = Vec::new();
                for (slot, &e) in parts.iter().enumerate() {
                    let n = if entities[e].is_arc() { 5 } else { 4 };
                    for j in 0..n {
                        let d = v.derivatives[5 * slot + j];
                        if d != 0.0 {
                            row.push((first[e] + j, d));
                        }
                    }
                }
                residuals.push(v.value);
                rows.push(row);
            }
        };
        let starting = |e: usize| -> Vec<D> { entities[e].variables().into_iter().map(D::constant).collect() };
        for c in constraints {
            match *c {
                Constraint::Fixed(e) => {
                    let (v, s) = (vars(e, 0), starting(e));
                    let r: Vec<D> = v.iter().zip(&s).map(|(a, b)| *a - *b).collect();
                    push(&r, &[e]);
                }
                Constraint::FixedPoint(e, at) => {
                    let arc = entities[e].is_arc();
                    let (p, p0) = (point(&vars(e, 0), arc, at), point(&starting(e), arc, at));
                    push(&[p[0] - p0[0], p[1] - p0[1]], &[e]);
                }
                Constraint::Coincident(a, b) => {
                    let (p, q) = (end_point(&vars(a, 0), entities[a].is_arc()), start_point(&vars(b, 1), entities[b].is_arc()));
                    push(&[p[0] - q[0], p[1] - q[1]], &[a, b]);
                }
                Constraint::Angle(a, b, angle) => {
                    let v1 = end_direction(&vars(a, 0), entities[a].is_arc());
                    let v2 = start_direction(&vars(b, 1), entities[b].is_arc());
                    // gp_Vec2d::Angle: from v2 to v1, anticlockwise positive.
                    push(&[angle_between(v2, v1) - angle], &[a, b]);
                }
                Constraint::Length(e, length) => {
                    let v = vars(e, 0);
                    let r = if entities[e].is_arc() { (v[2] * v[4]).abs() } else { distance([v[0], v[1]], [v[2], v[3]]) };
                    push(&[r - length], &[e]);
                }
                Constraint::Distance(a, fa, b, fb, d) => {
                    let p = point(&vars(a, 0), entities[a].is_arc(), fa);
                    let q = point(&vars(b, 1), entities[b].is_arc(), fb);
                    push(&[distance(p, q) - d], &[a, b]);
                }
                Constraint::Radius(e, radius) => {
                    let v = vars(e, 0);
                    push(&[v[2] - radius], &[e]);
                }
                Constraint::Orientation(e, along) => {
                    let v = vars(e, 0);
                    let direction = [v[2] - v[0], v[3] - v[1]];
                    push(&[angle_between(direction, [D::constant(along[0]), D::constant(along[1])])], &[e]);
                }
                Constraint::ArcAngle(e, sweep) => {
                    let v = vars(e, 0);
                    push(&[v[4] - sweep], &[e]);
                }
            }
        }
        Evaluation { residuals, rows }
    };
    let settings = Settings { max_steps: 10_000, step_relative: 1e-13, step_absolute: 1e-15, cost: 0.0, gradient: 1e-24 };
    let out = least_squares::minimize(&start, &lower, &upper, evaluate, settings);
    let solved = entities.iter().zip(&first).map(|(e, &k)| e.from_variables(&out.variables[k..])).collect();
    Ok(Solution { entities: solved, cost: out.cost, steps: out.steps, evaluations: out.evaluations, converged: out.stop != Stop::Steps })
}

/// An entity's point at fraction `at` of its parameter range (`None`: an arc's centre).
fn point(v: &[D], arc: bool, at: Option<f64>) -> [D; 2] {
    match (arc, at) {
        (true, None) => [v[0], v[1]],
        (true, Some(t)) => {
            let a = v[3] + v[4] * t;
            [v[0] + v[2] * a.sin(), v[1] + v[2] * a.cos()]
        }
        (false, Some(t)) => [v[0] + (v[2] - v[0]) * t, v[1] + (v[3] - v[1]) * t],
        (false, None) => unreachable!("checked before solving"),
    }
}

fn start_point(v: &[D], arc: bool) -> [D; 2] {
    if arc { point(v, true, Some(0.0)) } else { [v[0], v[1]] }
}

fn end_point(v: &[D], arc: bool) -> [D; 2] {
    if arc { point(v, true, Some(1.0)) } else { [v[2], v[3]] }
}

/// An arc's direction at angle `a`, the way it runs (CadQuery: sign(da)(cos a, −sin a)).
fn arc_direction(sweep: D, a: D) -> [D; 2] {
    let s = sweep.value.signum() * if sweep.value == 0.0 { 0.0 } else { 1.0 };
    [a.cos() * s, -(a.sin() * s)]
}

fn start_direction(v: &[D], arc: bool) -> [D; 2] {
    if arc { arc_direction(v[4], v[3]) } else { [v[2] - v[0], v[3] - v[1]] }
}

fn end_direction(v: &[D], arc: bool) -> [D; 2] {
    if arc { arc_direction(v[4], v[3] + v[4]) } else { [v[2] - v[0], v[3] - v[1]] }
}

/// The angle from `a` to `b`, anticlockwise positive, in (−π, π].
fn angle_between(a: [D; 2], b: [D; 2]) -> D {
    let cross = a[0] * b[1] - a[1] * b[0];
    let dot = a[0] * b[0] + a[1] * b[1];
    cross.atan2(dot)
}

fn distance(p: [D; 2], q: [D; 2]) -> D {
    let (dx, dy) = (p[0] - q[0], p[1] - q[1]);
    (dx * dx + dy * dy).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    /// A triangle's three lines, drawn roughly: the first fixed, each joined to the next, the
    /// second at 120° to the first and 2 long, the third 2 long: an equilateral triangle on
    /// the first line, (0, 0)–(2, 0).
    #[test]
    fn a_rough_triangle_made_equilateral() {
        let entities = [
            Entity::Line { start: [0.0, 0.0], end: [2.0, 0.0] },
            Entity::Line { start: [2.1, 0.1], end: [1.1, 1.6] },
            Entity::Line { start: [0.9, 1.8], end: [0.1, -0.1] },
        ];
        let constraints = [
            Constraint::Fixed(0),
            Constraint::Coincident(0, 1),
            Constraint::Coincident(1, 2),
            Constraint::Coincident(2, 0),
            Constraint::Angle(0, 1, -2.0 * PI / 3.0),
            Constraint::Length(1, 2.0),
            Constraint::Length(2, 2.0),
        ];
        let s = solve(&entities, &constraints).unwrap();
        assert!(s.converged && s.cost < 1e-20, "{s:?}");
        let Entity::Line { end, .. } = s.entities[1] else { panic!() };
        assert!((end[0] - 1.0).abs() < 1e-9 && (end[1] - 3f64.sqrt()).abs() < 1e-9, "{:?}", s.entities);
    }

    /// An arc's radius set and its centre held: the radius changes, the centre stays.
    #[test]
    fn an_arc_given_its_radius() {
        let entities = [Entity::Arc { centre: [1.0, 2.0], radius: 3.0, start: 0.0, sweep: PI / 2.0 }];
        let s = solve(&entities, &[Constraint::FixedPoint(0, None), Constraint::Radius(0, 5.0)]).unwrap();
        let Entity::Arc { centre, radius, .. } = s.entities[0] else { panic!() };
        assert!((radius - 5.0).abs() < 1e-12 && (centre[0] - 1.0).abs() < 1e-12 && (centre[1] - 2.0).abs() < 1e-12, "{s:?}");
    }
}
