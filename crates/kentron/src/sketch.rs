//! Planar profiles drawn by a moving pen: a point and a heading in a plane, moved by
//! straight lines and arcs, then closed into a face (or left as a wire).
//!
//! The pen starts at the plane's origin heading along its X axis. Lengths and coordinates
//! are in the plane's coordinates; angles are in degrees, positive turning left (toward Y).

use crate::error::{Error, Result};
use crate::geometry::curve::{Circle, Curve, Line};
use crate::geometry::surface::{Plane, Surface};
use crate::make;
use crate::math::{Direction, Frame, Point, Vector};
use crate::tolerance;
use crate::topology::{Model, Shape};

/// One step of a profile.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Step {
    /// Puts the pen at `(x, y)` without drawing; only before anything is drawn. The profile
    /// closes back to this point.
    Start(f64, f64),
    /// Moves the plane's origin to a point, keeping its axes.
    Origin(Point),
    /// Turns the plane: a new normal and X axis, keeping its origin.
    Plane { normal: Vector, x: Vector },
    /// Replaces the plane.
    Frame(Frame),
    /// A line along X by a length (negative: toward -X); the pen then heads along X.
    AlongX(f64),
    /// A line along X to an X coordinate.
    ToX(f64),
    /// A line along Y by a length; the pen then heads along Y.
    AlongY(f64),
    /// A line along Y to a Y coordinate.
    ToY(f64),
    /// A line ahead by a length (negative: backward, the pen then heading back).
    Forward(f64),
    /// A line ahead until the pen reaches an X coordinate.
    ForwardToX(f64),
    /// A line ahead until the pen reaches a Y coordinate.
    ForwardToY(f64),
    /// A line by a vector, which becomes the heading.
    By(f64, f64),
    /// A line to a point.
    To(f64, f64),
    /// Turns the heading by an angle.
    Turn(f64),
    /// Sets the heading to an angle from X.
    Heading(f64),
    /// Sets the heading to a direction.
    Direction(f64, f64),
    /// An arc tangent to the heading, of a radius (the centre on the left; negative: on the
    /// right), through an angle (negative: the other way round).
    Arc { radius: f64, angle: f64 },
}

/// What a profile makes of its edges.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Outcome {
    /// Closed back to the start, a face on the plane.
    #[default]
    Face,
    /// Closed back to the start, a wire.
    ClosedWire,
    /// The edges as drawn, a wire.
    OpenWire,
}

/// The profile drawn by `steps` on `plane`.
#[doc(alias = "DRAW profile")]
pub fn profile(model: &mut Model, plane: Frame, steps: &[Step], outcome: Outcome) -> Result<Shape> {
    let mut plane = plane;
    let (mut x0, mut y0, mut x, mut y, mut dx, mut dy) = (0.0f64, 0.0f64, 0.0f64, 0.0f64, 1.0f64, 0.0f64);
    let mut drawing = false;
    let mut edges: Vec<Shape> = Vec::new();
    let at = |plane: &Frame, u: f64, v: f64| plane.origin() + Vector::from(plane.x()) * u + Vector::from(plane.y()) * v;
    let direction = |v: Vector| Direction::try_from_vector(v).ok_or_else(|| Error::InvalidInput("profile: a null direction".into()));
    enum Move {
        Line(f64),
        Arc(f64, f64),
        None,
    }
    let closing = matches!(outcome, Outcome::Face | Outcome::ClosedWire);
    for (k, step) in steps.iter().copied().chain(closing.then_some(Step::To(f64::NAN, f64::NAN))).enumerate() {
        let setup = matches!(step, Step::Start(..) | Step::Origin(_) | Step::Plane { .. } | Step::Frame(_));
        let mv = match step {
            _ if k == steps.len() => {
                // Back to the start.
                let (vx, vy) = (x0 - x, y0 - y);
                let length = vx.hypot(vy);
                if length <= tolerance::LINEAR {
                    break;
                }
                (dx, dy) = (vx / length, vy / length);
                Move::Line(length)
            }
            Step::Start(u, v) => {
                if drawing {
                    return Err(Error::InvalidInput("profile: the start must precede all moves".into()));
                }
                (x0, y0, x, y) = (u, v, u, v);
                Move::None
            }
            Step::Origin(p) => {
                plane = Frame::new(p, plane.z(), plane.x());
                Move::None
            }
            Step::Plane { normal, x: x_axis } => {
                plane = Frame::new(plane.origin(), direction(normal)?, direction(x_axis)?);
                Move::None
            }
            Step::Frame(f) => {
                plane = f;
                Move::None
            }
            Step::AlongX(l) | Step::ToX(l) => {
                (dx, dy) = (1.0, 0.0);
                Move::Line(if matches!(step, Step::ToX(_)) { l - x } else { l })
            }
            Step::AlongY(l) | Step::ToY(l) => {
                (dx, dy) = (0.0, 1.0);
                Move::Line(if matches!(step, Step::ToY(_)) { l - y } else { l })
            }
            Step::Forward(l) => Move::Line(l),
            Step::ForwardToX(target) | Step::ForwardToY(target) => {
                let (d, from) = if matches!(step, Step::ForwardToX(_)) { (dx, x) } else { (dy, y) };
                if d.abs() < tolerance::LINEAR {
                    return Err(Error::InvalidInput(format!("profile: the heading never reaches {step:?}")));
                }
                Move::Line((target - from) / d)
            }
            Step::By(u, v) | Step::To(u, v) => {
                let (vx, vy) = if matches!(step, Step::To(..)) { (u - x, v - y) } else { (u, v) };
                let length = vx.hypot(vy);
                if length > tolerance::LINEAR {
                    (dx, dy) = (vx / length, vy / length);
                    Move::Line(length)
                } else {
                    Move::None
                }
            }
            Step::Turn(a) => {
                let (c, s) = (a.to_radians().cos(), a.to_radians().sin());
                (dx, dy) = (c * dx - s * dy, s * dx + c * dy);
                Move::None
            }
            Step::Heading(a) => {
                (dx, dy) = (a.to_radians().cos(), a.to_radians().sin());
                Move::None
            }
            Step::Direction(u, v) => {
                let length = u.hypot(v);
                if length > tolerance::LINEAR {
                    (dx, dy) = (u / length, v / length);
                }
                Move::None
            }
            Step::Arc { radius, angle } => {
                if radius.abs() > tolerance::LINEAR { Move::Arc(radius, angle.to_radians()) } else { Move::None }
            }
        };
        drawing |= !setup;
        match mv {
            Move::Line(mut length) => {
                if length < 0.0 {
                    length = -length;
                    (dx, dy) = (-dx, -dy);
                }
                let heading = Vector::from(plane.x()) * dx + Vector::from(plane.y()) * dy;
                let line = Line::new(at(&plane, x, y), direction(heading)?);
                edges.push(make::edge_from_curve_range(model, Curve::Line(line), 0.0, length)?);
                (x, y) = (x + length * dx, y + length * dy);
            }
            Move::Arc(mut radius, mut angle) => {
                let mut left = true;
                if radius < 0.0 {
                    radius = -radius;
                    left = !left;
                    (dx, dy) = (-dx, -dy);
                }
                let (cx, cy) = (x - radius * dy, y + radius * dx);
                // The circle's X axis points from its centre to the pen.
                let (ax, ay) = (dy, -dx);
                if angle < 0.0 {
                    angle = -angle;
                    left = !left;
                }
                let normal = if left { plane.z() } else { -plane.z() };
                let x_axis = direction(Vector::from(plane.x()) * ax + Vector::from(plane.y()) * ay)?;
                let circle = Circle::new(Frame::new(at(&plane, cx, cy), normal, x_axis), radius);
                edges.push(make::edge_from_curve_range(model, Curve::Circle(circle), 0.0, angle)?);
                let turn = if left { angle } else { -angle };
                let (c, s) = (turn.cos(), turn.sin());
                let (ux, uy) = (c * ax - s * ay, s * ax + c * ay);
                (x, y) = (cx + radius * ux, cy + radius * uy);
                (dx, dy) = if left { (-uy, ux) } else { (uy, -ux) };
            }
            Move::None => {}
        }
    }
    let wire = make::wire(model, &edges)?;
    match outcome {
        Outcome::Face => make::face_on(model, Surface::Plane(Plane::from_frame(plane)), wire),
        _ => Ok(wire),
    }
}
