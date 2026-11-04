//! CadQuery's sketch constraint costs (`cadquery.occ_impl.sketch_solver`), computed as it
//! computes them, and its objective and gradient: the sum of the costs' squares, and that
//! sum's derivatives by forward differences of each cost (each entity's variables moved by
//! 1e-10 in turn, in a copy of that entity's variables). SLSQP's steps depend on these
//! numbers, finite-difference error and all, so they are CadQuery's to the last operation
//! where that is cheap to keep (gp_Vec2d::Angle as OCCT computes it, CadQuery's order of
//! sums).

use super::sketch::{Constraint, Entity};

/// CadQuery's DIFF_EPS.
const STEP: f64 = 1e-10;

/// A constraint's entities' variables, as CadQuery passes them to a cost: the first's, then
/// the second's (if any), with each one's starting variables.
struct Arguments<'a> {
    first: &'a [f64],
    first_start: &'a [f64],
    first_arc: bool,
    second: &'a [f64],
    second_start: &'a [f64],
    second_arc: bool,
}

fn norm2(x: f64, y: f64) -> f64 {
    (x * x + y * y).sqrt()
}

/// OCCT's gp_Vec2d::Angle: the angle from `a` to `b`, anticlockwise positive.
fn angle(a: [f64; 2], b: [f64; 2]) -> f64 {
    let (na, nb) = (norm2(a[0], a[1]), norm2(b[0], b[1]));
    let d = na * nb;
    let cosine = (a[0] * b[0] + a[1] * b[1]) / d;
    let sine = (a[0] * b[1] - a[1] * b[0]) / d;
    if cosine > -0.70710678118655 && cosine < 0.70710678118655 {
        if sine > 0.0 { cosine.acos() } else { -cosine.acos() }
    } else if cosine > 0.0 {
        sine.asin()
    } else if sine > 0.0 {
        std::f64::consts::PI - sine.asin()
    } else {
        -std::f64::consts::PI - sine.asin()
    }
}

/// numpy's sign.
fn sign(v: f64) -> f64 {
    if v > 0.0 { 1.0 } else if v < 0.0 { -1.0 } else { 0.0 }
}

fn arc_point(x: &[f64], at: Option<f64>) -> [f64; 2] {
    match at {
        None => [x[0], x[1]],
        Some(t) => {
            let a = x[3] + t * x[4];
            [x[0] + x[2] * a.sin(), x[1] + x[2] * a.cos()]
        }
    }
}

fn arc_first(x: &[f64]) -> [f64; 2] {
    [x[0] + x[2] * x[3].sin(), x[1] + x[2] * x[3].cos()]
}

fn arc_last(x: &[f64]) -> [f64; 2] {
    [x[0] + x[2] * (x[3] + x[4]).sin(), x[1] + x[2] * (x[3] + x[4]).cos()]
}

fn line_point(x: &[f64], t: f64) -> [f64; 2] {
    [x[0] + t * (x[2] - x[0]), x[1] + t * (x[3] - x[1])]
}

fn arc_first_tangent(x: &[f64]) -> [f64; 2] {
    [sign(x[4]) * x[3].cos(), -sign(x[4]) * x[3].sin()]
}

fn arc_last_tangent(x: &[f64]) -> [f64; 2] {
    [sign(x[4]) * (x[3] + x[4]).cos(), -sign(x[4]) * (x[3] + x[4]).sin()]
}

/// A point at a fraction (`None`: an arc's centre; checked before for lines).
fn point(x: &[f64], arc: bool, at: Option<f64>) -> [f64; 2] {
    if arc { arc_point(x, at) } else { line_point(x, at.expect("a fraction along a line")) }
}

/// The cost of constraint `c` (CadQuery's cost functions; their sign and form as there).
fn cost(c: &Constraint, a: &Arguments) -> f64 {
    let x1 = a.first;
    match *c {
        Constraint::Fixed(_) => {
            let s: f64 = x1.iter().zip(a.first_start).map(|(v, w)| (v - w) * (v - w)).sum();
            s.sqrt()
        }
        Constraint::FixedPoint(_, at) => {
            let (p, q) = (point(x1, a.first_arc, at), point(a.first_start, a.first_arc, at));
            norm2(p[0] - q[0], p[1] - q[1])
        }
        Constraint::Coincident(..) => {
            let x2 = a.second;
            let v1 = if a.first_arc { arc_last(x1) } else { [x1[2], x1[3]] };
            let v2 = if a.second_arc { arc_first(x2) } else { [x2[0], x2[1]] };
            norm2(v1[0] - v2[0], v1[1] - v2[1])
        }
        Constraint::Angle(_, _, value) => {
            let x2 = a.second;
            let v1 = if a.first_arc { arc_last_tangent(x1) } else { [x1[2] - x1[0], x1[3] - x1[1]] };
            let v2 = if a.second_arc { arc_first_tangent(x2) } else { [x2[2] - x2[0], x2[3] - x2[1]] };
            angle(v2, v1) - value
        }
        Constraint::Length(_, value) => {
            if a.first_arc { (x1[2] * x1[4]).abs() - value } else { norm2(x1[2] - x1[0], x1[3] - x1[1]) - value }
        }
        Constraint::Distance(_, at1, _, at2, d) => {
            let (v1, v2) = (point(x1, a.first_arc, at1), point(a.second, a.second_arc, at2));
            norm2(v1[0] - v2[0], v1[1] - v2[1]) - d
        }
        Constraint::Radius(_, value) => x1[2] - value,
        Constraint::Orientation(_, along) => angle([x1[2] - x1[0], x1[3] - x1[1]], along),
        Constraint::ArcAngle(_, value) => x1[4] - value,
    }
}

/// The entities a constraint is on (the second `None` for one on one entity).
fn entities_of(c: &Constraint) -> (usize, Option<usize>) {
    match *c {
        Constraint::Fixed(e) | Constraint::FixedPoint(e, _) | Constraint::Length(e, _) | Constraint::Radius(e, _) | Constraint::Orientation(e, _) | Constraint::ArcAngle(e, _) => (e, None),
        Constraint::Coincident(a, b) | Constraint::Angle(a, b, _) | Constraint::Distance(a, _, b, _, _) => (a, Some(b)),
    }
}

/// CadQuery's sketch objective: the variables of `entities` (each's run `first[e]..`), their
/// starting values `start`, and the constraints.
pub(super) struct Objective<'a> {
    pub entities: &'a [Entity],
    pub constraints: &'a [Constraint],
    pub first: &'a [usize],
    pub start: &'a [f64],
}

impl Objective<'_> {
    fn size(&self, e: usize) -> usize {
        if matches!(self.entities[e], Entity::Arc { .. }) { 5 } else { 4 }
    }

    fn arc(&self, e: usize) -> bool {
        matches!(self.entities[e], Entity::Arc { .. })
    }

    fn slice<'b>(&self, x: &'b [f64], e: usize) -> &'b [f64] {
        &x[self.first[e]..self.first[e] + self.size(e)]
    }

    /// The sum of the costs' squares at `x`, in CadQuery's order; with its gradient (by
    /// CadQuery's forward differences) into `gradient` when given.
    pub fn evaluate(&self, x: &[f64], gradient: Option<&mut [f64]>) -> f64 {
        let mut value = 0.0;
        for c in self.constraints {
            let (e1, e2) = entities_of(c);
            let e2v = e2.unwrap_or(e1);
            let a = Arguments {
                first: self.slice(x, e1),
                first_start: self.slice(self.start, e1),
                first_arc: self.arc(e1),
                second: self.slice(x, e2v),
                second_start: self.slice(self.start, e2v),
                second_arc: self.arc(e2v),
            };
            let r = cost(c, &a);
            value += r * r;
        }
        if let Some(g) = gradient {
            g.fill(0.0);
            for c in self.constraints {
                let (e1, e2) = entities_of(c);
                let e2v = e2.unwrap_or(e1);
                let mut x1 = self.slice(x, e1).to_vec();
                let mut x2 = self.slice(x, e2v).to_vec();
                let base = |x1: &[f64], x2: &[f64]| {
                    cost(c, &Arguments {
                        first: x1,
                        first_start: self.slice(self.start, e1),
                        first_arc: self.arc(e1),
                        second: x2,
                        second_start: self.slice(self.start, e2v),
                        second_arc: self.arc(e2v),
                    })
                };
                let r = base(&x1, &x2);
                for j in 0..x1.len() {
                    let saved = x1[j];
                    x1[j] += STEP;
                    let moved = base(&x1, &x2);
                    g[self.first[e1] + j] += 2.0 * r * (moved - r) / STEP;
                    x1[j] = saved;
                }
                if e2.is_some() {
                    for j in 0..x2.len() {
                        let saved = x2[j];
                        x2[j] += STEP;
                        let moved = base(&x1, &x2);
                        g[self.first[e2v] + j] += 2.0 * r * (moved - r) / STEP;
                        x2[j] = saved;
                    }
                }
            }
        }
        value
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// OCCT's angle from (1, 0) to (−1, 1): 135°, through its arcsine branch.
    #[test]
    fn angle_as_occt_computes_it() {
        assert!((angle([1.0, 0.0], [-1.0, 1.0]) - 0.75 * std::f64::consts::PI).abs() < 1e-15);
        assert!((angle([1.0, 0.0], [0.0, -1.0]) + std::f64::consts::FRAC_PI_2).abs() < 1e-15);
    }
}
