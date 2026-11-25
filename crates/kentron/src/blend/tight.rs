//! Fillets along chains of edges with rounds tighter than the ball.
//!
//! A convex edge between a plane and a cylinder square to it (the rim of a prism with rounded
//! corners), rounded by a ball wider than the cylinder: no ball of that radius rolls round the
//! round touching both faces, its centre would pass the cylinder's axis and its fillet fold
//! over itself. The largest ball touching both that rolls round it is as wide as the round,
//! centred on its axis: it stands still there, and its fillet is a piece of its sphere, from
//! the plane (where it touches it at one point) down to the cylinder (round the circle it
//! touches it in).
//!
//! The straight edges carrying the chain on from the round's ends take the ball from that
//! width back up to the fillet's radius, a radius changing smoothly along each (flat at both
//! ends of the change): the fillet then meets its faces tangentially everywhere and bends
//! smoothly from the round's sphere into each straight edge's fillet. The change takes twice
//! the length the ball, wedged between the two sides meeting past the round, would need to
//! grow from the round's width to the radius (where it touches the round's plane, the side
//! and the other side): at any point of it the ball is narrower than that wedged ball, so it
//! never reaches the other side; and no less than twice the growth itself, so that the radius
//! changes no faster than three quarters of the length along the edge.

use std::f64::consts::PI;

use super::{Corner, Path};
use crate::geometry::surface::Surface;
use crate::math::{Point, Vector};
use crate::topology::{Model, Orientation, Shape};
use crate::tolerance;

/// A round of the chain tighter than the ball: its radius and the angle it turns through.
#[derive(Clone, Copy)]
struct Round {
    radius: f64,
    angle: f64,
}

/// The rounds' change on one side of a straight edge: the radius at its end there, and the
/// length along it over which the ball grows back to the fillet's radius.
#[derive(Clone, Copy)]
struct Growth {
    from: f64,
    length: f64,
}

/// `smoothstep`: 0 at 0, 1 at 1 and beyond, flat at both.
fn smooth(t: f64) -> f64 {
    let t = t.clamp(0.0, 1.0);
    t * t * (3.0 - 2.0 * t)
}

/// How much longer than the wedged ball's growth a change of the ball must be at the least
/// for the ball to stay narrower than the wedged one all along it: past the round's end by
/// `x`, the changing ball's growth is `smooth(x / length) (r − ρ)`, the wedged ball's
/// `x tan(angle / 2)`, and `smooth(t) / t` peaks at 9/8 (at t = 3/4).
const WEDGED: f64 = 9.0 / 8.0;

/// The ball's radius along each of `edges` (one tangent chain of `shape`, all filleted by
/// `radius`) as the fraction of the edge's length from its start, where some of its convex
/// edges are arcs round cylinders narrower than the ball (see the module's notes). `None`
/// where none is, or where the chain is not such rounds between straight edges carrying it
/// on from their ends, or a straight edge is too short for the ball to change along it.
#[allow(clippy::type_complexity)]
pub(super) fn tight_laws(model: &Model, shape: Shape, edges: &[Shape], radius: f64) -> Option<Vec<Box<dyn Fn(f64) -> f64>>> {
    let corners: Vec<Corner> = edges.iter().map(|&e| Corner::of(model, shape, e, false).ok().flatten()).collect::<Option<_>>()?;
    let rounds: Vec<Option<Round>> = corners.iter().map(|c| tight_round(c, radius)).collect();
    if rounds.iter().all(Option::is_none) {
        return None;
    }
    // The edges' ends (their forward start and end) and which edge meets each there.
    let ends: Vec<[Point; 2]> = edges
        .iter()
        .map(|&e| {
            let (a, b) = model.edge_vertices(Shape { orientation: Orientation::Forward, ..e });
            [model.point(a), model.point(b)]
        })
        .collect();
    let near = |p: Point, q: Point| p.distance(q) <= 10.0 * tolerance::LINEAR;
    let neighbour = |i: usize, p: Point| -> Option<usize> {
        let at: Vec<usize> = (0..edges.len()).filter(|&j| j != i && ends[j].iter().any(|&q| near(p, q))).collect();
        match at[..] {
            [j] => Some(j),
            _ => None,
        }
    };
    // Each round between straight convex edges carrying it on tangentially.
    for (i, round) in rounds.iter().enumerate() {
        if round.is_none() {
            continue;
        }
        for &p in &ends[i] {
            let j = neighbour(i, p)?;
            let line = &corners[j];
            if !matches!(line.path, Path::Line { .. }) || !line.convex {
                return None;
            }
            let (a, b) = (corners[i].direction_at(p), line.direction_at(p));
            if a.cross(b).length() > 1e-6 * a.length() * b.length() {
                return None;
            }
        }
    }
    let mut laws: Vec<Box<dyn Fn(f64) -> f64>> = Vec::with_capacity(edges.len());
    for (i, corner) in corners.iter().enumerate() {
        if let Some(round) = rounds[i] {
            let r = round.radius;
            laws.push(Box::new(move |_| r));
            continue;
        }
        let Path::Line { length } = corner.path else {
            laws.push(Box::new(move |_| radius));
            continue;
        };
        // The growth at each end next to a round.
        let mut growth: [Option<Growth>; 2] = [None, None];
        for (k, &p) in ends[i].iter().enumerate() {
            let Some(round) = neighbour(i, p).and_then(|j| rounds[j]) else { continue };
            let wedged = (radius - round.radius) * (0.5 * round.angle).tan();
            growth[k] = Some(Growth { from: round.radius, length: 2.0 * wedged.max(radius - round.radius) });
        }
        // (Both growths not fitting along the edge: each shortened alike, while each is still
        // a twentieth longer than the shortest keeping the ball narrower than the wedged one.)
        let total: f64 = growth.iter().flatten().map(|g| g.length).sum();
        if total > length {
            let scale = length / total;
            for (k, g) in growth.iter_mut().enumerate() {
                let Some(g) = g else { continue };
                let Some(round) = neighbour(i, ends[i][k]).and_then(|j| rounds[j]) else { continue };
                g.length *= scale;
                if g.length < 1.05 * WEDGED * (radius - round.radius) * (0.5 * round.angle).tan() {
                    return None;
                }
            }
        }
        laws.push(Box::new(move |u: f64| {
            let x = u.clamp(0.0, 1.0) * length;
            let mut r = radius;
            if let Some(g) = growth[0] {
                r -= (radius - g.from) * (1.0 - smooth(x / g.length));
            }
            if let Some(g) = growth[1] {
                r -= (radius - g.from) * (1.0 - smooth((length - x) / g.length));
            }
            r
        }));
    }
    Some(laws)
}

/// The round `corner` is, where it is one tighter than a ball of `radius`: a convex arc where
/// a plane meets a cylinder square to it, about the arc's axis, the cylinder narrower than the
/// ball (and turning through less than a half turn).
fn tight_round(corner: &Corner, radius: f64) -> Option<Round> {
    let Path::Arc { axis, radius: arc, angle } = corner.path else { return None };
    if !corner.convex || !(arc > tolerance::LINEAR) || radius <= arc + tolerance::LINEAR || angle >= PI - 1e-6 {
        return None;
    }
    let z = Vector::from(axis.direction);
    let (mut plane, mut cylinder) = (false, false);
    for (_, surface, _) in &corner.faces {
        match surface {
            Surface::Plane(p) => plane |= Vector::from(p.normal()).cross(z).length() <= 1e-9,
            Surface::Cylinder(c) => {
                let w = Vector::from(c.frame().z());
                let off = c.frame().origin() - axis.origin;
                cylinder |= w.cross(z).length() <= 1e-9 && (off - z * off.dot(z)).length() <= tolerance::LINEAR && (c.radius() - arc).abs() <= tolerance::LINEAR;
            }
            _ => {}
        }
    }
    (plane && cylinder).then_some(Round { radius: arc, angle })
}
