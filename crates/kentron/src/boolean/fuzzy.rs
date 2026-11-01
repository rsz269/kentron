//! Fuzzy booleans: operands meant to meet that miss by a little.
//!
//! OCCT's fuzzy value (`BRepAlgoAPI_BooleanOperation::SetFuzzyValue`) widens the tolerance of
//! every intersection, so faces within it of each other count as coincident. Kentron's
//! operations build exact geometry instead, so the tool is first made to meet: each of its
//! plane faces lying within the fuzzy value of a plane face of the other operand (and
//! overlapping it) is taken onto that plane, and its vertices moved to where their faces'
//! planes now meet. The operation then runs on the snapped tool as on any other.

use std::collections::HashMap;

use crate::error::{Error, Result};
use crate::geometry::curve::Curve;
use crate::geometry::surface::{Plane, Surface};
use crate::math::{BoundingBox, Point, Transformable, Vector};
use crate::topology::{Model, NodeId, NodeKind, Orientation, Shape};

/// `tool` with its plane faces lying within `tolerance` of plane faces of `reference` (and
/// overlapping them) taken onto those planes, its vertices moved with them: a new shape,
/// or `tool` itself where nothing is that close, or where it is not bounded by planes and
/// straight edges only (shapes other than polyhedra are left as they are).
#[doc(alias = "SetFuzzyValue")]
pub fn snapped(model: &mut Model, tool: Shape, reference: Shape, tolerance: f64) -> Result<Shape> {
    if tolerance <= crate::tolerance::LINEAR {
        return Ok(tool);
    }
    Ok(snap(model, tool, reference, tolerance)?.unwrap_or(tool))
}

/// A plane as `n · x = d`.
#[derive(Clone, Copy)]
struct Equation {
    normal: Vector,
    d: f64,
}

impl Equation {
    fn of(plane: &Plane) -> Self {
        let normal = Vector::from(plane.normal());
        Equation { normal, d: normal.dot(plane.origin() - Point::ORIGIN) }
    }

    fn distance(&self, p: Point) -> f64 {
        self.normal.dot(p - Point::ORIGIN) - self.d
    }
}

fn snap(model: &mut Model, tool: Shape, reference: Shape, tolerance: f64) -> Result<Option<Shape>> {
    // The reference's planes, where they are, with their faces' extents.
    let targets: Vec<(Equation, BoundingBox)> = model
        .explore(reference, NodeKind::Face)
        .collect::<Vec<_>>()
        .into_iter()
        .filter_map(|f| match model.surface(f)? {
            Surface::Plane(p) => Some((Equation::of(&p.transformed(&f.placement)), crate::properties::bounding_box(model, f).enlarged(tolerance))),
            _ => None,
        })
        .collect();
    if targets.is_empty() {
        return Ok(None);
    }
    let (copy, _) = crate::modify::copy_unplaced(model, tool, &[])?;
    let solids: Vec<Shape> = model.explore(copy, NodeKind::Solid).collect();
    if solids.is_empty() {
        return Ok(None);
    }
    // Only polyhedra: plane faces and straight edges.
    for f in model.explore(copy, NodeKind::Face).collect::<Vec<_>>() {
        if !matches!(model.surface(f), Some(Surface::Plane(_))) {
            return Ok(None);
        }
    }
    for e in model.explore(copy, NodeKind::Edge).collect::<Vec<_>>() {
        if !matches!(model.edge_curve(e), Some((Curve::Line(_), _, _))) {
            return Ok(None);
        }
    }
    // Each face's plane: the reference's it lies within `tolerance` of, if any, else its own.
    let mut planes: HashMap<NodeId, Equation> = HashMap::new();
    let mut moved = false;
    for f in model.explore(copy, NodeKind::Face).collect::<Vec<_>>() {
        if planes.contains_key(&f.node) {
            continue;
        }
        let Some(Surface::Plane(plane)) = model.surface(f) else { unreachable!() };
        let own = Equation::of(plane);
        let points: Vec<Point> = model.explore(f, NodeKind::Vertex).map(|v| model.point(v)).collect();
        let extent = BoundingBox::from_points(points.iter().copied());
        let best = targets
            .iter()
            .filter(|(_, b)| b.overlaps(&extent))
            .map(|(t, _)| (t, points.iter().map(|&p| t.distance(p).abs()).fold(0.0, f64::max)))
            .filter(|&(_, gap)| gap <= tolerance * (1.0 + 1e-6) + crate::tolerance::LINEAR)
            .min_by(|a, b| a.1.total_cmp(&b.1));
        let equation = match best {
            Some((t, gap)) if gap > crate::tolerance::LINEAR => {
                moved = true;
                *t
            }
            _ => own,
        };
        planes.insert(f.node, equation);
    }
    if !moved {
        return Ok(None);
    }
    // Each vertex where its faces' planes meet, as near where it was as they allow.
    let mut points: HashMap<NodeId, Point> = HashMap::new();
    for f in model.explore(copy, NodeKind::Face).collect::<Vec<_>>() {
        for v in model.explore(f, NodeKind::Vertex).collect::<Vec<_>>() {
            points.entry(v.node).or_insert_with(|| model.point(v));
        }
    }
    let mut around: HashMap<NodeId, Vec<Equation>> = HashMap::new();
    for f in model.explore(copy, NodeKind::Face).collect::<Vec<_>>() {
        let equation = planes[&f.node];
        let mut seen = Vec::new();
        for v in model.explore(f, NodeKind::Vertex).collect::<Vec<_>>() {
            if !seen.contains(&v.node) {
                seen.push(v.node);
                around.entry(v.node).or_default().push(equation);
            }
        }
    }
    for (v, equations) in &around {
        let p = points[v];
        let Some(q) = meeting_point(p, equations) else { return Ok(None) };
        if equations.iter().any(|e| e.distance(q).abs() > 1e-9 * (1.0 + q.distance(Point::ORIGIN))) || q.distance(p) > 2.0 * tolerance {
            return Ok(None);
        }
        points.insert(*v, q);
    }
    // The polyhedron again from the moved vertices: shared vertices and edges, plane faces.
    let attempt = (|| -> Result<Shape> {
        let mut vertices: HashMap<NodeId, Shape> = HashMap::new();
        for (&v, &p) in &points {
            vertices.insert(v, crate::make::vertex(model, p));
        }
        let mut edges: HashMap<NodeId, Shape> = HashMap::new();
        let mut rebuilt = Vec::new();
        for solid in &solids {
            let new_solid = model.new_empty(NodeKind::Solid);
            for shell in model.children(*solid).collect::<Vec<_>>() {
                let new_shell = model.new_empty(NodeKind::Shell);
                for f in model.explore(shell, NodeKind::Face).collect::<Vec<_>>() {
                    let face = rebuilt_face(model, f, &vertices, &mut edges)?;
                    model.add(new_shell, face)?;
                }
                model.add(new_solid, if shell.orientation == Orientation::Reversed { new_shell.reversed() } else { new_shell })?;
            }
            rebuilt.push(new_solid);
        }
        if let [one] = rebuilt[..] {
            return Ok(one);
        }
        let compound = model.new_empty(NodeKind::Compound);
        for s in rebuilt {
            model.add(compound, s)?;
        }
        Ok(compound)
    })();
    Ok(attempt.ok())
}

/// The face `f` again (on the plane its moved vertices span), through the new `vertices`, its edges shared
/// through `edges` (new ones, forward as the old, by the old edge's node).
fn rebuilt_face(model: &mut Model, f: Shape, vertices: &HashMap<NodeId, Shape>, edges: &mut HashMap<NodeId, Shape>) -> Result<Shape> {
    let id = crate::topology::FaceId::try_from(f.node).map_err(|n| Error::InvalidInput(format!("{n:?} is not a face")))?;
    let outer = model.outer_wire(id);
    let mut outer_wire = None;
    let mut holes = Vec::new();
    for w in model.children(Shape { orientation: Orientation::Forward, ..f }).collect::<Vec<_>>() {
        let mut new_edges = Vec::new();
        for e in model.ordered_edges(w) {
            let new = match edges.get(&e.node) {
                Some(&new) => new,
                None => {
                    let (start, end) = model.edge_vertices(Shape { orientation: Orientation::Forward, ..e });
                    let new = crate::make::edge(model, vertices[&start.node], vertices[&end.node])?;
                    edges.insert(e.node, new);
                    new
                }
            };
            new_edges.push(if e.orientation == Orientation::Reversed { new.reversed() } else { new });
        }
        let wire = crate::make::wire(model, &new_edges)?;
        if outer.is_some_and(|o| crate::topology::NodeId::from(o) == w.node) {
            outer_wire = Some(wire);
        } else {
            holes.push(wire);
        }
    }
    let outer_wire = outer_wire.ok_or_else(|| Error::InvalidInput("a face without its outer wire".into()))?;
    let mut face = crate::make::planar_face(model, outer_wire)?;
    if !holes.is_empty() {
        face = crate::make::face_with_holes(model, face, &holes)?;
    }
    // Facing as the old face did (outward, in its solid).
    let Some(Surface::Plane(new_plane)) = model.surface(face) else { return Err(Error::InvalidInput("a face off its plane".into())) };
    let Some(Surface::Plane(old_plane)) = model.surface(f) else { unreachable!() };
    let old_out = Vector::from(old_plane.normal()) * if f.orientation == Orientation::Reversed { -1.0 } else { 1.0 };
    Ok(if Vector::from(new_plane.normal()).dot(old_out) < 0.0 { face.reversed() } else { face })
}

/// The point nearest `p` on all the planes `equations` (which meet, if the snapping is
/// sound): through the most independent three of them, two or one.
fn meeting_point(p: Point, equations: &[Equation]) -> Option<Point> {
    let mut chosen: Vec<Equation> = Vec::new();
    for e in equations {
        let independent = match chosen.len() {
            0 => true,
            1 => chosen[0].normal.cross(e.normal).length() > 1e-6,
            2 => chosen[0].normal.cross(chosen[1].normal).dot(e.normal).abs() > 1e-6,
            _ => false,
        };
        if independent {
            chosen.push(*e);
        }
    }
    let x = Vector::from(p - Point::ORIGIN);
    let q = match chosen[..] {
        [a] => x - a.normal * a.distance(p),
        [a, b] => {
            // Nearest point on the line where two planes meet: x + α a + β b.
            let (aa, ab, bb) = (a.normal.dot(a.normal), a.normal.dot(b.normal), b.normal.dot(b.normal));
            let (ra, rb) = (-a.distance(p), -b.distance(p));
            let det = aa * bb - ab * ab;
            let (alpha, beta) = ((ra * bb - rb * ab) / det, (rb * aa - ra * ab) / det);
            x + a.normal * alpha + b.normal * beta
        }
        [a, b, c] => {
            let det = a.normal.cross(b.normal).dot(c.normal);
            (b.normal.cross(c.normal) * a.d + c.normal.cross(a.normal) * b.d + a.normal.cross(b.normal) * c.d) * (1.0 / det)
        }
        _ => return None,
    };
    Some(Point::ORIGIN + q)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::primitives::Cuboid;

    fn cuboid(model: &mut Model, low: [f64; 3], high: [f64; 3]) -> Shape {
        Cuboid::from_corners(Point::new(low[0], low[1], low[2]), Point::new(high[0], high[1], high[2])).build(model).unwrap()
    }

    #[test]
    fn boxes_apart_by_less_than_the_fuzzy_value_fuse_into_one() {
        let mut model = Model::new();
        let a = cuboid(&mut model, [0.0; 3], [1.0; 3]);
        let b = cuboid(&mut model, [1.001, 0.0, 0.0], [2.001, 1.0, 1.0]);
        let tool = snapped(&mut model, b, a, 1e-3).unwrap();
        let fused = crate::boolean::union(&mut model, a, tool).unwrap();
        assert_eq!(model.explore(fused, NodeKind::Solid).count(), 1);
        assert!((crate::properties::volume(&model, fused).mass - 2.001).abs() < 1e-9);
        // Farther apart than the value: left as they are.
        let c = cuboid(&mut model, [1.01, 0.0, 0.0], [2.01, 1.0, 1.0]);
        assert_eq!(snapped(&mut model, c, a, 1e-3).unwrap(), c);
    }

    #[test]
    fn a_box_nearly_on_another_cut_by_it_leaves_nothing() {
        let mut model = Model::new();
        let a = cuboid(&mut model, [0.0; 3], [1.0; 3]);
        let b = cuboid(&mut model, [0.001, 0.0, 0.0], [1.001, 1.0, 1.0]);
        let tool = snapped(&mut model, b, a, 1e-3).unwrap();
        let cut = crate::boolean::difference(&mut model, a, tool).unwrap();
        assert!(crate::properties::volume(&model, cut).mass.abs() < 1e-9);
    }
}
