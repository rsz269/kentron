//! Constraint solving for CadQuery's assemblies and sketches (`pykentron.constraint_solvers`):
//! the problems in plain numbers, solved by `kentron::constraints`.

use kentron::constraints::{assembly, sketch};
use kentron::math::{Point, Quaternion, Vector};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

type Triple = [f64; 3];

/// What a solve reports: the cost left, steps tried, evaluations made, whether it converged.
type Report = (f64, usize, usize, bool);

/// The placements of an assembly's parts. `starts`: each part's translation and rotation
/// parameters (CadQuery's: the axis times tan(angle / 4)); `locked`: parts kept in place;
/// `constraints`: (kind, parts, numbers), the numbers by kind:
/// Point (p₁, p₂, distance), Axis (d₁, d₂, mode 0 same / 1 opposite / 2 at an angle, angle),
/// PointInPlane (p, origin, normal, offset), PointOnLine (p, origin, direction, distance),
/// FixedPoint (p, target), FixedAxis (d, target), FixedRotation (w, x, y, z). With the
/// report, whether it converged to the tolerance (not only to the acceptable level), and each
/// iteration's Hessian regularization (IPOPT's "rg").
#[pyfunction]
#[allow(clippy::type_complexity)]
pub fn solve_assembly(starts: Vec<(Triple, Triple)>, locked: Vec<bool>, constraints: Vec<(String, Vec<usize>, Vec<f64>)>, scale: f64) -> PyResult<(Vec<(Triple, Triple)>, Report, bool, Vec<f64>)> {
    if locked.len() != starts.len() {
        return Err(PyValueError::new_err("one locked flag per part"));
    }
    let v = |d: &[f64], i: usize| Vector::new(d[i], d[i + 1], d[i + 2]);
    let p = |d: &[f64], i: usize| Point::new(d[i], d[i + 1], d[i + 2]);
    let mut parsed = Vec::with_capacity(constraints.len());
    for (kind, parts, d) in &constraints {
        let need = |numbers: usize, count: usize| -> PyResult<()> {
            if d.len() != numbers || parts.len() != count || parts.iter().any(|&i| i >= starts.len()) {
                return Err(PyValueError::new_err(format!("{kind}: {count} parts and {numbers} numbers expected")));
            }
            Ok(())
        };
        parsed.push(match kind.as_str() {
            "Point" => {
                need(7, 2)?;
                assembly::Constraint::Point { first: parts[0], on_first: p(d, 0), second: parts[1], on_second: p(d, 3), distance: d[6] }
            }
            "Axis" => {
                need(8, 2)?;
                let angle = match d[6] as i64 {
                    0 => assembly::AxisAngle::Same,
                    1 => assembly::AxisAngle::Opposite,
                    _ => assembly::AxisAngle::Angle(d[7]),
                };
                assembly::Constraint::Axis { first: parts[0], on_first: v(d, 0), second: parts[1], on_second: v(d, 3), angle }
            }
            "PointInPlane" => {
                need(10, 2)?;
                assembly::Constraint::PointInPlane { first: parts[0], point: p(d, 0), second: parts[1], origin: p(d, 3), normal: v(d, 6), offset: d[9] }
            }
            "PointOnLine" => {
                need(10, 2)?;
                assembly::Constraint::PointOnLine { first: parts[0], point: p(d, 0), second: parts[1], origin: p(d, 3), direction: v(d, 6), distance: d[9] }
            }
            "FixedPoint" => {
                need(6, 1)?;
                assembly::Constraint::FixedPoint { part: parts[0], point: p(d, 0), target: p(d, 3) }
            }
            "FixedAxis" => {
                need(6, 1)?;
                assembly::Constraint::FixedAxis { part: parts[0], along: v(d, 0), target: v(d, 3) }
            }
            "FixedRotation" => {
                need(4, 1)?;
                assembly::Constraint::FixedRotation { part: parts[0], rotation: Quaternion::new(d[1], d[2], d[3], d[0]) }
            }
            other => return Err(PyValueError::new_err(format!("unknown assembly constraint {other}"))),
        });
    }
    let starts: Vec<assembly::Placement> = starts.iter().map(|&(t, r)| assembly::Placement { translation: Vector::new(t[0], t[1], t[2]), rotation: Vector::new(r[0], r[1], r[2]) }).collect();
    let s = assembly::solve(&starts, &locked, &parsed, scale);
    let placements = s.placements.iter().map(|q| ([q.translation.x, q.translation.y, q.translation.z], [q.rotation.x, q.rotation.y, q.rotation.z])).collect();
    Ok((placements, (s.cost, s.steps, s.evaluations, s.converged), s.optimal, s.regularizations))
}

/// A sketch's edges moved to satisfy constraints. `entities`: ("LINE", [x₁, y₁, x₂, y₂]) or
/// ("CIRCLE", [x, y, r, a, da]); `constraints`: (kind, entity, second entity or None,
/// parameters) as CadQuery's sketch takes them (angles in radians; a missing fraction None).
#[pyfunction]
#[allow(clippy::type_complexity)]
pub fn solve_sketch(entities: Vec<(String, Vec<f64>)>, constraints: Vec<(String, usize, Option<usize>, Vec<Option<f64>>)>) -> PyResult<(Vec<Vec<f64>>, Report)> {
    let mut parsed_entities = Vec::with_capacity(entities.len());
    for (geom, x) in &entities {
        parsed_entities.push(match (geom.as_str(), x.len()) {
            ("LINE", 4) => sketch::Entity::Line { start: [x[0], x[1]], end: [x[2], x[3]] },
            ("CIRCLE", 5) => sketch::Entity::Arc { centre: [x[0], x[1]], radius: x[2], start: x[3], sweep: x[4] },
            _ => return Err(PyValueError::new_err(format!("a sketch entity {geom} of {} numbers", x.len()))),
        });
    }
    let mut parsed = Vec::with_capacity(constraints.len());
    for (kind, a, b, params) in &constraints {
        let n = parsed_entities.len();
        let second = || -> PyResult<usize> { b.filter(|&b| b < n).ok_or_else(|| PyValueError::new_err(format!("{kind}: a second entity expected"))) };
        let number = |i: usize| -> PyResult<f64> { params.get(i).copied().flatten().ok_or_else(|| PyValueError::new_err(format!("{kind}: a parameter expected"))) };
        if *a >= n {
            return Err(PyValueError::new_err(format!("{kind}: no entity {a}")));
        }
        parsed.push(match kind.as_str() {
            "Fixed" => sketch::Constraint::Fixed(*a),
            "FixedPoint" => sketch::Constraint::FixedPoint(*a, params.first().copied().flatten()),
            "Coincident" => sketch::Constraint::Coincident(*a, second()?),
            "Angle" => sketch::Constraint::Angle(*a, second()?, number(0)?),
            "Length" => sketch::Constraint::Length(*a, number(0)?),
            "Distance" => sketch::Constraint::Distance(*a, params.first().copied().flatten(), second()?, params.get(1).copied().flatten(), number(2)?),
            "Radius" => sketch::Constraint::Radius(*a, number(0)?),
            "Orientation" => sketch::Constraint::Orientation(*a, [number(0)?, number(1)?]),
            "ArcAngle" => sketch::Constraint::ArcAngle(*a, number(0)?),
            other => return Err(PyValueError::new_err(format!("unknown sketch constraint {other}"))),
        });
    }
    let s = sketch::solve(&parsed_entities, &parsed).map_err(|e| PyValueError::new_err(e.to_string()))?;
    let out = s
        .entities
        .iter()
        .map(|e| match *e {
            sketch::Entity::Line { start, end } => vec![start[0], start[1], end[0], end[1]],
            sketch::Entity::Arc { centre, radius, start, sweep } => vec![centre[0], centre[1], radius, start, sweep],
        })
        .collect();
    Ok((out, (s.cost, s.steps, s.evaluations, s.converged)))
}
