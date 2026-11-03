//! Offsets outwards joined by intersection of solids with planar faces and straight edges,
//! found as the space they cover: the solid, with a prism over each face (the face swept
//! the distance along its outward normal), a wedge along each convex edge (between the two
//! faces' prisms, out to where their moved planes meet) and a corner at each vertex all of
//! whose edges are convex (out to where its faces' moved planes meet), all united. Moving the
//! solid's faces one by one fails where the offset closes a notch or a face vanishes (its
//! neighbours' offsets meeting past it); the pieces' union closes it as it covers it.

use crate::error::Result;
use crate::geometry::curve::Curve;
use crate::geometry::surface::Surface;
use crate::math::{Point, Transformable, Vector};
use crate::topology::{Model, NodeId, NodeKind, Orientation, Shape};

/// `solid` offset outwards by `distance` (positive), joined by intersection, as the union of
/// its pieces; `None` where it is not one solid of planar faces and straight edges.
pub(super) fn offset(model: &mut Model, solid: Shape, distance: f64) -> Result<Option<Shape>> {
    let solid = match (solid.kind(), &model.unique_subshapes(solid, NodeKind::Solid)[..]) {
        (NodeKind::Solid, _) => solid,
        (NodeKind::Compound, [one]) => *one,
        _ => return Ok(None),
    };
    if !(distance > 0.0) {
        return Ok(None);
    }
    let d = distance;
    // The faces' planes, by outward normal.
    let faces = model.unique_subshapes(solid, NodeKind::Face);
    let mut normals: Vec<Vector> = Vec::with_capacity(faces.len());
    let mut origins: Vec<Point> = Vec::with_capacity(faces.len());
    for &f in &faces {
        let Some(Surface::Plane(p)) = model.surface(f).map(|s| s.transformed(&f.placement)) else { return Ok(None) };
        let n = Vector::from(p.normal());
        normals.push(if f.orientation == Orientation::Reversed { -n } else { n });
        origins.push(p.origin());
    }
    let size = crate::properties::bounding_box(model, solid).size().length().max(1e-9);
    let edges_of: Vec<Vec<NodeId>> = faces.iter().map(|&f| model.unique_subshapes(f, NodeKind::Edge).into_iter().map(|e| e.node).collect()).collect();
    let inside = crate::boolean::PointClassifier::new(model, solid)?;
    let mut pieces: Vec<Shape> = Vec::new();
    // Over each face, its prism.
    for (i, &f) in faces.iter().enumerate() {
        let face = Shape { orientation: Orientation::Forward, ..f };
        pieces.push(crate::sweep::extrude(model, face, normals[i] * d)?.shape);
    }
    // Along each convex edge, its wedge; and which vertices have a concave edge.
    let mut at_vertex: crate::hashing::FastMap<NodeId, (Point, Vec<usize>, bool)> = Default::default();
    for edge in model.unique_subshapes(solid, NodeKind::Edge) {
        let using: Vec<usize> = (0..faces.len()).filter(|&i| edges_of[i].contains(&edge.node)).collect();
        let [i, j] = using[..] else { return Ok(None) };
        let Some((Curve::Line(_), _, _)) = model.edge_curve(edge) else { return Ok(None) };
        let (a, b) = model.edge_vertices(edge);
        let (pa, pb) = (model.point(a), model.point(b));
        let (ni, nj) = (normals[i], normals[j]);
        let flat = ni.cross(nj).length() <= 1e-9;
        // (Convex where the point just under one face's plane and out past the other's, by
        // the edge, lies outside the solid: locally its material is the half-spaces' common
        // part, at a concave edge their union.)
        let convex = !flat && {
            let middle = pa + (pb - pa) * 0.5;
            let step = 1e-4 * size;
            inside.locate(middle - ni * step + nj * step) == crate::boolean::PointLocation::Outside
        };
        for (v, p) in [(a.node, pa), (b.node, pb)] {
            let entry = at_vertex.entry(v).or_insert((p, Vec::new(), false));
            for k in [i, j] {
                if !entry.1.contains(&k) {
                    entry.1.push(k);
                }
            }
            entry.2 |= !convex && !flat;
        }
        if !convex {
            continue;
        }
        let miter = (ni + nj) * (d / (1.0 + ni.dot(nj)));
        let points: Vec<Point> = [pa, pb].iter().flat_map(|&e| [e, e + ni * d, e + miter, e + nj * d]).collect();
        if let Some(wedge) = crate::blend::convex_solid(model, &points)? {
            pieces.push(wedge);
        }
    }
    // (A convex solid's faces all stay, moved: nothing to close.)
    if !at_vertex.values().any(|(_, _, concave)| *concave) {
        return Ok(None);
    }
    // At each vertex all of whose edges are convex, its corner: out to where its faces' moved
    // planes meet (where three or more meet in one point).
    for (_, (v, around, concave)) in at_vertex {
        if concave || around.len() < 3 {
            continue;
        }
        let mut points = vec![v];
        points.extend(around.iter().map(|&k| v + normals[k] * d));
        for (x, &k) in around.iter().enumerate() {
            for &l in &around[x + 1..] {
                let (nk, nl) = (normals[k], normals[l]);
                if nk.cross(nl).length() > 1e-9 {
                    points.push(v + (nk + nl) * (d / (1.0 + nk.dot(nl))));
                }
            }
        }
        // (The point on all the moved planes: by the first three independent normals,
        // where the others pass through it too.)
        let n: Vec<Vector> = around.iter().map(|&k| normals[k]).collect();
        let found = (0..n.len()).find_map(|x| {
            (x + 1..n.len()).find_map(|y| {
                (y + 1..n.len()).find_map(|z| {
                    let det = n[x].dot(n[y].cross(n[z]));
                    (det.abs() > 1e-9).then(|| (n[y].cross(n[z]) * d + n[z].cross(n[x]) * d + n[x].cross(n[y]) * d) * (1.0 / det))
                })
            })
        });
        if let Some(m) = found {
            if n.iter().all(|k| (k.dot(m) - d).abs() <= 1e-9 * size.max(d)) {
                points.push(v + m);
            }
        }
        if let Some(corner) = crate::blend::convex_solid(model, &points)? {
            pieces.push(corner);
        }
    }
    // All united with the solid.
    let mut result = solid;
    for piece in pieces {
        result = crate::boolean::union(model, result, piece)?;
    }
    crate::unify::unify_same_domain(model, result).map(Some)
}
