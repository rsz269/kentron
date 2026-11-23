//! Everything a solid passes through as it moves rigidly: the union of `M(t)·body` for `t` in
//! `[0, 1]`, with an exact boundary (not a union of sampled places).
//!
//! # The envelope
//!
//! Seen from the moving body, each point `p` of it moves at `V(t, p)`. Where `V·n > 0` (`n` the
//! outward normal) the body's surface moves out into space it did not hold; call that region
//! of the surface `W(t)`. A point the moving body covers at some `t` lies in the body at the
//! start or was reached through `W`. So the swept set is the body at `t = 0` together with what
//! `W` sweeps, and its boundary lies among (Wang and Wang 1986; Martin and Stephenson 1990;
//! Weld and Leu 1990):
//!
//! * the body's surface outside `W(0)` at the start, and `W(1)` at the end;
//! * what the boundary of `W` sweeps: its **characteristic curves**, where `V·n = 0` on a
//!   smooth face, sweep the envelope, tangent to the moving body all along them; its **sharp
//!   edges**, between a face in `W` and one outside it, sweep surfaces of their own.
//!
//! For a slide, a turn or a screw (each a one-parameter group of motions) `V` is the same field
//! at every `t`, so `W` and its curves are fixed on the body, and the curves' sweeps are the
//! products of fixed curves and the motion: a slide sweeps a curve into its surface of linear
//! extrusion, a turn into its surface of revolution (both exact), a screw into the B-spline
//! surface whose columns follow its poles' helices (as [`super::screw_extrude`] sweeps an edge,
//! within about 1e-10 of the body's size).
//!
//! # Winding
//!
//! Those faces, oriented as the boundary of the 3-chain `body(0) + {M(t)·p : p ∈ W, t ∈ [0,
//! 1]}`, make a closed shell: the "product". Every preimage of a point under `(t, p) ↦ M(t)·p`
//! on `W` counts positively (the map's Jacobian there has the sign of `V·n`), so the product's
//! winding number about a point is the body at the start (1 inside) plus how many times `W`
//! passes over it: at least 1 exactly where the point is swept. Where the product does not run
//! into itself it is the swept solid. Where it does (a fold, where the motion turns past the
//! line it grazes along; a sheet a concave edge sweeps into the body; a turn coming back over
//! where it began) its faces are split where they cross ([`crate::boolean::enclosed_solids`])
//! and the pieces of space it winds round at least once are the swept solid, joined.
//!
//! Nothing is tuned: the characteristic curves are followed within a tolerance taken from the
//! body's size, and refusals are stated geometric conditions.

mod carried;
mod contour;
mod product;

use std::f64::consts::TAU;

use crate::error::{Error, Result};
use crate::math::{Axis, Point, Transform, Vector};
use crate::topology::{Model, NodeKind, Shape};

/// A rigid motion given by where it puts the body at each `t` in `[0, 1]` (a frame carried
/// along a curve, a screw whose rise is any function of its turn): its places, smooth in `t`.
#[derive(Clone)]
pub struct Carriage(pub std::sync::Arc<dyn Fn(f64) -> Transform + Send + Sync>);

impl std::fmt::Debug for Carriage {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str("Carriage")
    }
}

/// How a solid is carried, `t` running from 0 to 1, from where it is drawn.
#[derive(Debug, Clone)]
pub enum RigidMotion {
    /// Carried straight by `t·by`.
    Slide { by: Vector },
    /// Turned right-handed about `axis` from angle `from` to `to` (radians; the body is drawn
    /// at angle 0).
    Turn { axis: Axis, from: f64, to: f64 },
    /// Turned right-handed about `axis` from `from` to `to` (radians) while rising `lead` along
    /// it for every whole turn (the body is drawn at angle 0, unrisen).
    Screw { axis: Axis, lead: f64, from: f64, to: f64 },
    /// Placed by the given transforms, at each `t` (the body is drawn where the identity
    /// leaves it).
    Carried(Carriage),
}

impl RigidMotion {
    /// Where the motion puts the body at `t`.
    pub fn at(&self, t: f64) -> Transform {
        match *self {
            RigidMotion::Carried(ref c) => (c.0)(t),
            RigidMotion::Slide { by } => Transform::translation(by * t),
            RigidMotion::Turn { axis, from, to } => Transform::rotation(axis, from + (to - from) * t),
            RigidMotion::Screw { axis, lead, from, to } => {
                let angle = from + (to - from) * t;
                Transform::translation(Vector::from(axis.direction) * (lead * angle / TAU)) * Transform::rotation(axis, angle)
            }
        }
    }

    /// The place at the start and the steady motion from there (turning positively).
    fn steady(&self) -> Result<(Transform, Steady)> {
        let finite = |x: f64| x.is_finite();
        match *self {
            RigidMotion::Carried(_) => Err(Error::InvalidInput("a carried motion is not steady".into())),
            RigidMotion::Slide { by } => {
                if !(by.length() > 0.0 && finite(by.length())) {
                    return Err(Error::InvalidInput("a slide by no distance".into()));
                }
                Ok((Transform::IDENTITY, Steady::Slide { by }))
            }
            RigidMotion::Turn { axis, from, to } => {
                if !(finite(from) && finite(to) && from != to) {
                    return Err(Error::InvalidInput(format!("a turn from {from} to {to}")));
                }
                let turning = if to > from { axis } else { Axis::new(axis.origin, -axis.direction) };
                Ok((Transform::rotation(axis, from), Steady::Turn { axis: turning, angle: (to - from).abs() }))
            }
            RigidMotion::Screw { axis, lead, from, to } => {
                if !(finite(from) && finite(to) && from != to && finite(lead)) {
                    return Err(Error::InvalidInput(format!("a screw from {from} to {to} of lead {lead}")));
                }
                let start = self.at(0.0);
                let rise = lead * (to - from) / TAU;
                let axis_turning = if to > from { axis } else { Axis::new(axis.origin, -axis.direction) };
                // (Along the turning axis: the rise as signed along it.)
                let rise = if to > from { rise } else { -rise };
                Ok((start, Steady::Screw { axis: axis_turning, angle: (to - from).abs(), rise }))
            }
        }
    }
}

/// A motion whose velocity field is the same at every instant (seen from the body, which here
/// stands at its start): `t` from 0 to 1.
#[derive(Debug, Clone)]
pub(super) enum Steady {
    /// A stretch of a general motion over which `W` stays put (its field the one midway).
    Carried(std::sync::Arc<carried::Stage>),
    Slide { by: Vector },
    /// Turning by `angle` (positive) right-handed about `axis`.
    Turn { axis: Axis, angle: f64 },
    /// Turning by `angle` (positive) right-handed about `axis`, rising `rise` along it.
    Screw { axis: Axis, angle: f64, rise: f64 },
}

impl Steady {
    /// The velocity (in `t`) of the point at `p` at the start.
    pub(super) fn velocity(&self, p: Point) -> Vector {
        match *self {
            Steady::Carried(ref stage) => stage.velocity(p),
            Steady::Slide { by } => by,
            Steady::Turn { axis, angle } => Vector::from(axis.direction).cross(p - axis.origin) * angle,
            Steady::Screw { axis, angle, rise } => Vector::from(axis.direction).cross(p - axis.origin) * angle + Vector::from(axis.direction) * rise,
        }
    }

    /// Where the motion puts the start at `t`.
    pub(super) fn at(&self, t: f64) -> Transform {
        match *self {
            Steady::Carried(ref stage) => stage.at(t),
            Steady::Slide { by } => Transform::translation(by * t),
            Steady::Turn { axis, angle } => Transform::rotation(axis, angle * t),
            Steady::Screw { axis, angle, rise } => Transform::translation(Vector::from(axis.direction) * (rise * t)) * Transform::rotation(axis, angle * t),
        }
    }

    /// A scale of the speeds, for telling a vanishing `V·n`: the slide's, or the turn's at a
    /// unit's distance plus the rise.
    pub(super) fn speed_scale(&self) -> f64 {
        match *self {
            Steady::Carried(ref stage) => stage.velocity(Point::ORIGIN).length() + stage.velocity(Point::new(1.0, 0.0, 0.0)).length(),
            Steady::Slide { by } => by.length(),
            Steady::Turn { angle, .. } => angle,
            Steady::Screw { angle, rise, .. } => angle + rise.abs(),
        }
    }
}

/// The solid `tool` passes through as `motion` carries it, `t` from 0 to 1: the union of
/// `motion.at(t)·tool`.
///
/// Taken: one solid of one shell, under a slide, a turn or a screw. A turn about an axis
/// meeting the tool is refused (a point on the axis stands still there, and the curves where
/// the motion grazes the tool meet at it), as is a turn of a whole turn or more (its ends would
/// be one).
pub fn solid(model: &mut Model, tool: Shape, motion: &RigidMotion) -> Result<Shape> {
    let solids = model.unique_subshapes(tool, NodeKind::Solid);
    if solids.len() != 1 || model.unique_subshapes(solids[0], NodeKind::Shell).len() != 1 {
        return Err(Error::InvalidInput("a sweep carries one solid of one shell".into()));
    }
    if let RigidMotion::Carried(carriage) = motion {
        return carried::sweep(model, solids[0], carriage);
    }
    let (start, steady) = motion.steady()?;
    if let Steady::Turn { angle, .. } = steady {
        if angle >= TAU {
            return Err(Error::InvalidInput("a turn of a whole turn or more ends where it starts".into()));
        }
    }
    let placed = crate::modify::transform(model, solids[0], &start, Default::default())?;
    if let Steady::Turn { axis, .. } = steady {
        if meets_axis(model, placed, axis)? {
            return Err(Error::InvalidInput("a turn about an axis meeting the tool".into()));
        }
    }
    product::sweep(model, placed, &steady)
}

/// Whether the line `axis` meets the solid (within the tolerance).
fn meets_axis(model: &mut Model, solid: Shape, axis: Axis) -> Result<bool> {
    let bounds = crate::properties::bounding_box(model, solid);
    let reach = (bounds.max - bounds.min).length() + (axis.origin - bounds.min).length();
    let d = Vector::from(axis.direction);
    let line = crate::make::edge_between(model, axis.origin - d * reach, axis.origin + d * reach)?;
    let distance = crate::distance::distance(model, solid, line)?.value;
    Ok(distance <= crate::tolerance::LINEAR)
}

#[cfg(test)]
mod tests;
