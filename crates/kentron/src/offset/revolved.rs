//! Offsets of partial solids of revolution (a face of a plane through an axis revolved less
//! than a whole turn: a cone, a sphere or a revolved profile cut to 90° or 270°) whose side
//! planes move where the profile reaches the axis, at a pole or an apex. Moved, a side plane
//! leaves the axis, and the pole is no longer on the offset: the curved faces' offsets run
//! over it, cut short by the moved plane, which moving the faces one by one (keeping the
//! solid's topology, its degenerate edge at the pole) cannot give.
//!
//! The offset is the whole solid of revolution's (its poles staying on the axis), cut by the
//! side planes' half-spaces, moved: joined by intersection, the curved faces' offsets meet
//! the moved planes where they are extended round the axis, as planes through their axis
//! are square to them.

use std::collections::HashMap;

use crate::error::Result;
use crate::geometry::extrema::nearest_on_surface;
use crate::geometry::surface::Surface;
use crate::math::{Axis, Direction, Frame, Point, Transform, Transformable, Vector};
use crate::primitives::Cuboid;
use crate::topology::{Model, NodeId, NodeKind, Orientation, Shape};

use super::Join;

/// `solid` with each face moved `distance_of` it along its outward normal, joined by
/// intersection, where it is a partial solid of revolution reaching its axis at a pole or
/// apex and a side plane moves; `None` otherwise (or where the result is not sound).
pub(super) fn offset(model: &mut Model, solid: Shape, distance_of: &dyn Fn(NodeId) -> f64) -> Result<Option<Shape>> {
    let Some(rev) = partial_revolution(model, solid) else { return Ok(None) };
    let moves = [distance_of(rev.sides[0].face.node), distance_of(rev.sides[1].face.node)];
    if moves == [0.0, 0.0] || !rev.pole {
        return Ok(None);
    }
    // The whole solid of revolution: a side face revolved a turn (its seam on that side's
    // plane), the faces its edges sweep moved as the faces of the solid beside those edges;
    // from the other side where the boolean cutting it there misses (the seam lying along
    // the cut).
    for k in 0..2 {
        if let Some(found) = cut_revolution(model, &rev, k, moves, distance_of)? {
            return Ok(Some(found));
        }
    }
    Ok(None)
}

/// [`offset`] with the whole solid of revolution swept from side `k`: `None` where the cut
/// leaves nothing, all of it, or not one sound solid.
fn cut_revolution(model: &mut Model, rev: &Revolution, k: usize, moves: [f64; 2], distance_of: &dyn Fn(NodeId) -> f64) -> Result<Option<Shape>> {
    let profile = Shape { orientation: Orientation::Forward, ..rev.sides[k].face };
    let beside = &rev.beside[k];
    let sweep = crate::sweep::revolve(model, profile, rev.axis, std::f64::consts::TAU)?;
    let mut own: Vec<(Shape, f64)> = Vec::new();
    for edge in model.unique_subshapes(profile, NodeKind::Edge) {
        let Some(&beside) = beside.get(&edge.node) else { continue };
        for swept in sweep.generated(edge) {
            if swept.kind() == NodeKind::Face {
                own.push((swept, distance_of(beside)));
            }
        }
    }
    let Some(&(_, first)) = own.first() else { return Ok(None) };
    let moved = if own.iter().all(|&(_, d)| d == first) {
        super::shape::shape(model, sweep.shape, first, Join::Intersection)
    } else {
        super::shape::shape_with_faces(model, sweep.shape, first, &own, Join::Intersection)
    };
    let Ok(moved) = moved else { return Ok(None) };
    // The side planes' half-spaces, moved: met by the solid where it turns less than half a
    // turn, joined where more.
    let bounds = crate::properties::bounding_box(model, moved);
    let size = 4.0 * bounds.size().length().max(1.0);
    let middle = bounds.center();
    // (Each half-space a box with a face on the moved plane, inside it or outside.)
    let half = |model: &mut Model, side: &Side, d: f64, outside: bool| -> Result<Shape> {
        let on = side.origin + side.normal * d;
        let normal = Direction::from_vector(if outside { -side.normal } else { side.normal });
        let n = Vector::from(normal);
        // (The middle of the box's face on the plane: the bounds' middle, moved onto it.)
        let top = middle - n * n.dot(middle - on);
        let frame = Frame::new(Point::ORIGIN, normal, if normal.dot(Direction::X).abs() < 0.9 { Direction::X } else { Direction::Y });
        let corner = top - Vector::from(frame.x()) * (0.5 * size) - Vector::from(frame.y()) * (0.5 * size) - n * size;
        Cuboid::new(size, size, size).with_frame(Frame::new(corner, normal, frame.x())).build(model)
    };
    let result = if rev.convex {
        // Within both half-spaces.
        let a = half(model, &rev.sides[0], moves[0], false)?;
        let b = half(model, &rev.sides[1], moves[1], false)?;
        crate::boolean::intersection(model, a, b).and_then(|w| crate::boolean::intersection(model, moved, w))
    } else {
        // Within either: the wedge outside both cut away.
        let a = half(model, &rev.sides[0], moves[0], true)?;
        let b = half(model, &rev.sides[1], moves[1], true)?;
        crate::boolean::intersection(model, a, b).and_then(|w| crate::boolean::difference(model, moved, w))
    };
    let Ok(result) = result else { return Ok(None) };
    let (whole, left) = (crate::properties::volume(model, moved).mass, crate::properties::volume(model, result).mass);
    let cut = left > 1e-9 * whole && left < whole * (1.0 - 1e-9);
    let sound = cut && model.unique_subshapes(result, NodeKind::Solid).len() == 1 && crate::validity::check(model, result).is_valid();
    Ok(sound.then_some(result))
}

/// A side face of a partial solid of revolution: a face of a plane through the axis.
struct Side {
    face: Shape,
    origin: Point,
    /// Its plane's normal, outwards.
    normal: Vector,
}

/// A solid of a face revolved about an axis less than a whole turn.
struct Revolution {
    axis: Axis,
    sides: [Side; 2],
    /// Whether it turns less than half a turn (it lies within both sides' half-spaces, not
    /// either's).
    convex: bool,
    /// Whether a curved face reaches the axis (a vertex on it: a pole, an apex).
    pole: bool,
    /// Each edge of either side face with the other face of the solid beside it.
    beside: [HashMap<NodeId, NodeId>; 2],
}

/// `solid` as a partial solid of revolution: two planar faces whose planes meet in a line
/// about which every other face's surface turns into itself.
fn partial_revolution(model: &Model, solid: Shape) -> Option<Revolution> {
    let solid = match solid.kind() {
        NodeKind::Solid => solid,
        NodeKind::Compound => match model.unique_subshapes(solid, NodeKind::Solid)[..] {
            [one] => one,
            _ => return None,
        },
        _ => return None,
    };
    if model.unique_subshapes(solid, NodeKind::Shell).len() != 1 {
        return None;
    }
    let faces = model.unique_subshapes(solid, NodeKind::Face);
    let surface_of = |f: Shape| model.surface(f).map(|s| s.transformed(&f.placement));
    let size = crate::properties::bounding_box(model, solid).size().length().max(1e-9);
    let near = 1e-7 * size;
    let planes: Vec<(usize, Point, Vector)> = faces
        .iter()
        .enumerate()
        .filter_map(|(i, &f)| match surface_of(f)? {
            Surface::Plane(p) => {
                let n = Vector::from(p.normal());
                Some((i, p.origin(), if f.orientation == Orientation::Reversed { -n } else { n }))
            }
            _ => None,
        })
        .collect();
    // Points along a face's edges.
    let samples = |f: Shape| -> Vec<Point> {
        let mut out = Vec::new();
        for e in model.unique_subshapes(f, NodeKind::Edge) {
            if model.is_degenerate(e) {
                continue;
            }
            let Some((curve, a, b)) = model.edge_curve(e) else { continue };
            for t in [0.13, 0.5, 0.87] {
                use crate::geometry::ParametricCurve;
                out.push(curve.point(a + t * (b - a)));
            }
        }
        out
    };
    for (x, &(i, oi, ni)) in planes.iter().enumerate() {
        for &(j, oj, nj) in &planes[x + 1..] {
            let along = ni.cross(nj);
            if along.length() < 1e-6 {
                continue;
            }
            // The line where the planes meet: n·p = h for both.
            let (hi, hj) = (ni.dot(oi - Point::ORIGIN), nj.dot(oj - Point::ORIGIN));
            let on = Point::ORIGIN + (nj.cross(along) * hi + along.cross(ni) * hj) * (1.0 / along.length_squared());
            let axis = Axis::new(on, Direction::from_vector(along));
            // Every other face turning into itself about the line.
            let turns = |f: Shape| {
                let Some(surface) = surface_of(f) else { return false };
                let points = samples(f);
                !points.is_empty()
                    && points.iter().all(|&p| {
                        [0.3, -0.3].iter().any(|&angle| {
                            let q = p.transformed(&Transform::rotation(axis, angle));
                            nearest_on_surface(q, &surface, None).2 <= 10.0 * near
                        })
                    })
            };
            if !faces.iter().enumerate().all(|(k, &f)| k == i || k == j || turns(f)) {
                continue;
            }
            let off_axis = |p: Point| {
                let v = p - on;
                v - Vector::from(axis.direction) * v.dot(Vector::from(axis.direction))
            };
            // The sides' directions away from the axis, and the angle the solid turns from the
            // first to the second, turning into it (against the first's outward normal).
            let away = |f: Shape| samples(f).into_iter().map(off_axis).max_by(|a, b| a.length().total_cmp(&b.length()));
            let (Some(ui), Some(uj)) = (away(faces[i]), away(faces[j])) else { return None };
            if ui.length() <= near || uj.length() <= near {
                return None;
            }
            let (ui, uj) = (ui.normalized(), uj.normalized());
            let mut turn = Vector::from(axis.direction);
            if turn.cross(ui).dot(-ni) < 0.0 {
                turn = -turn;
            }
            let mut angle = ui.cross(uj).dot(turn).atan2(ui.dot(uj));
            if angle <= 0.0 {
                angle += std::f64::consts::TAU;
            }
            let curved_on_axis = faces.iter().enumerate().any(|(k, &f)| {
                k != i && k != j && !matches!(surface_of(f), Some(Surface::Plane(_))) && model.unique_subshapes(f, NodeKind::Vertex).into_iter().any(|v| off_axis(model.point(v)).length() <= 10.0 * near)
            });
            let beside_of = |side: usize| {
                let edges: Vec<NodeId> = model.unique_subshapes(faces[side], NodeKind::Edge).into_iter().map(|e| e.node).collect();
                let mut beside = HashMap::new();
                for (k, &f) in faces.iter().enumerate() {
                    if k == side {
                        continue;
                    }
                    for e in model.unique_subshapes(f, NodeKind::Edge) {
                        if edges.contains(&e.node) {
                            beside.insert(e.node, f.node);
                        }
                    }
                }
                beside
            };
            let beside = [beside_of(i), beside_of(j)];
            return Some(Revolution {
                axis,
                sides: [Side { face: faces[i], origin: oi, normal: ni }, Side { face: faces[j], origin: oj, normal: nj }],
                convex: angle < std::f64::consts::PI,
                pole: curved_on_axis,
                beside,
            });
        }
    }
    None
}
