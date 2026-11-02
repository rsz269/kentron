//! Removing features of a solid: faces (a hole's, a boss's, a step's) taken away, the faces
//! around them carried on over where they were to close the solid again.
//!
//! The faces next to the removed ones are extended on their surfaces (a plane far across
//! the features, a cylinder or cone all round and along, a sphere or torus whole), and
//! with the removed faces they split space into cells. The gaps are the cells against the
//! removed faces (reaching as far as the faces carried on close them: a sloping face
//! carried on up to the wall across, where a top is removed): one outside the solid is a
//! hole or pocket, filled; one inside is a boss, cut away. The faces lying on one surface are merged after, so that each face
//! around a feature is carried on over it as one face.

use crate::error::{Error, Result};
use crate::geometry::surface::Surface;
use crate::geometry::ParametricSurface;
use crate::math::{BoundingBox, Point, Transformable};
use crate::topology::{History, Model, NodeId, NodeKind, Shape};

/// A solid with features removed: the result, and what became of the solid's faces, edges
/// and vertices (those of the features deleted, faces around them modified).
pub struct Removed {
    pub shape: Shape,
    pub history: History,
}

/// `solid` (a solid, or a compound of one) without the faces `faces`, the faces around them
/// carried on over where they were.
#[doc(alias = "BRepAlgoAPI_Defeaturing")]
#[doc(alias = "removefeatures")]
pub fn remove_faces(model: &mut Model, solid: Shape, faces: &[Shape]) -> Result<Removed> {
    let solids = model.unique_subshapes(solid, NodeKind::Solid);
    let [solid] = solids[..] else {
        return Err(Error::InvalidInput("features are removed from one solid".into()));
    };
    let all: Vec<Shape> = model.unique_subshapes(solid, NodeKind::Face);
    let removed: Vec<Shape> = all.iter().copied().filter(|f| faces.iter().any(|g| g.node == f.node)).collect();
    if removed.is_empty() {
        return Err(Error::InvalidInput("no faces of the solid to remove".into()));
    }
    // The faces around the features: sharing an edge with one of them.
    let removed_edges: std::collections::HashSet<NodeId> = removed.iter().flat_map(|&f| model.unique_subshapes(f, NodeKind::Edge)).map(|e| e.node).collect();
    let around: Vec<Shape> = all
        .iter()
        .copied()
        .filter(|f| !removed.iter().any(|g| g.node == f.node))
        .filter(|&f| model.unique_subshapes(f, NodeKind::Edge).iter().any(|e| removed_edges.contains(&e.node)))
        .collect();
    // The neighbourhood: the box of the features and the faces around them.
    let mut near = BoundingBox::EMPTY;
    for &f in removed.iter().chain(&around) {
        near.add_box(&crate::properties::bounding_box(model, f));
    }
    // The faces around carried on, and the features: the cells they bound.
    let mut sheets: Vec<Shape> = removed.clone();
    for &f in &around {
        sheets.push(extended(model, f, &near)?);
    }
    let cells = crate::boolean::enclosed_solids(model, &sheets)?;
    let cells = model.unique_subshapes(cells, NodeKind::Solid);
    let classifier = crate::boolean::PointClassifier::new(model, solid)?;
    let (mut fill, mut cut) = (Vec::new(), Vec::new());
    for &cell in &cells {
        if !touches(model, cell, &removed) {
            continue;
        }
        let Some(p) = crate::boolean::point_inside(model, cell)? else { continue };
        match classifier.locate(p) {
            crate::boolean::PointLocation::Outside => fill.push(cell),
            crate::boolean::PointLocation::Inside => cut.push(cell),
            crate::boolean::PointLocation::Boundary => {}
        }
    }
    if fill.is_empty() && cut.is_empty() {
        return Err(Error::Failed("removing features: the faces around them do not close over them".into()));
    }
    let mut result = solid;
    for cell in fill {
        result = crate::boolean::union(model, result, cell)?;
    }
    for cell in cut {
        result = crate::boolean::difference(model, result, cell)?;
    }
    let options = crate::unify::UnifyOptions { same_surface_only: true, ..Default::default() };
    let result = crate::unify::unify_same_domain_with(model, result, options)?;
    // (Then, as OCCT's removal, the faces it made or changed merged where their surfaces
    // are alike (a face carried on and the piece of another on the same sphere); the
    // solid's other faces kept as they are.)
    let untouched: crate::hashing::FastSet<NodeId> = model
        .unique_subshapes(solid, NodeKind::Face)
        .into_iter()
        .map(|f| f.node)
        .collect::<crate::hashing::FastSet<NodeId>>()
        .intersection(&model.unique_subshapes(result, NodeKind::Face).into_iter().map(|f| f.node).collect())
        .copied()
        .filter(|n| !around.iter().any(|a| a.node == *n))
        .collect();
    let options = crate::unify::UnifyOptions { faces: true, edges: false, same_surface_only: false };
    let result = crate::unify::unify_same_domain_keeping(model, result, options, &untouched)?;
    // (Then edges running on one curve merged, wherever they came from: an edge of a face
    // carried on meets the one it continues end to end.)
    let options = crate::unify::UnifyOptions { faces: false, edges: true, same_surface_only: false };
    let result = crate::unify::unify_same_domain_with(model, result, options)?;
    let history = history(model, solid, &removed, result);
    Ok(Removed { shape: result, history })
}

/// `face` carried on over the neighbourhood `near` on its surface (a plane across it, a
/// cylinder or cone round its axis and along it, a sphere or torus whole); the face itself
/// on other surfaces.
fn extended(model: &mut Model, face: Shape, near: &BoundingBox) -> Result<Shape> {
    let Some(surface) = model.surface(face).cloned() else { return Ok(face) };
    let corners: Vec<Point> = [near.min.x, near.max.x]
        .into_iter()
        .flat_map(|x| [near.min.y, near.max.y].into_iter().flat_map(move |y| [near.min.z, near.max.z].into_iter().map(move |z| Point::new(x, y, z))))
        .collect();
    let size = (near.max - near.min).length().max(1e-9);
    let range = |values: &mut dyn Iterator<Item = f64>| {
        let (lo, hi) = values.fold((f64::INFINITY, f64::NEG_INFINITY), |(a, b), v| (a.min(v), b.max(v)));
        (lo - size, hi + size)
    };
    let tau = std::f64::consts::TAU;
    let (u0, u1, v0, v1) = match &surface {
        Surface::Plane(_) => {
            let uvs: Vec<(f64, f64)> = corners.iter().filter_map(|&p| surface.parameters(p)).collect();
            let (u0, u1) = range(&mut uvs.iter().map(|uv| uv.0));
            let (v0, v1) = range(&mut uvs.iter().map(|uv| uv.1));
            (u0, u1, v0, v1)
        }
        Surface::Cylinder(_) => {
            let (v0, v1) = range(&mut corners.iter().filter_map(|&p| surface.parameters(p)).map(|uv| uv.1));
            (0.0, tau, v0, v1)
        }
        Surface::Cone(c) => {
            // (Along the generators on the face's side of the apex.)
            let (v0, v1) = range(&mut corners.iter().filter_map(|&p| surface.parameters(p)).map(|uv| uv.1));
            let (_, _, fv0, fv1) = model.face_uv_bounds(face);
            let apex = -c.radius() / c.semi_angle_sin_cos().0;
            let (v0, v1) = if 0.5 * (fv0 + fv1) > apex { (v0.max(apex + 1e-6 * size), v1) } else { (v0, v1.min(apex - 1e-6 * size)) };
            (0.0, tau, v0, v1)
        }
        Surface::Sphere(_) => (0.0, tau, -std::f64::consts::FRAC_PI_2, std::f64::consts::FRAC_PI_2),
        Surface::Torus(_) => (0.0, tau, 0.0, tau),
        _ => return Ok(face),
    };
    let big = crate::make::face_on_domain(model, surface, u0, u1, v0, v1)?;
    Ok(Shape { orientation: face.orientation, ..big })
}

/// Whether `cell` lies against one of `faces`: a point inside one of its faces is on one
/// of them.
fn touches(model: &Model, cell: Shape, faces: &[Shape]) -> bool {
    model.unique_subshapes(cell, NodeKind::Face).iter().any(|&f| {
        let Some(p) = face_point(model, f) else { return false };
        faces.iter().any(|&g| on_face(model, g, p))
    })
}

/// Whether `p` lies on `face`.
fn on_face(model: &Model, face: Shape, p: Point) -> bool {
    crate::boolean::PointClassifier::new(model, face).is_ok_and(|c| c.locate(p) == crate::boolean::PointLocation::Boundary)
}

/// A point inside `face`.
fn face_point(model: &Model, face: Shape) -> Option<Point> {
    let surface = model.surface(face)?;
    let loops = crate::heal::face_loops(model, Shape { orientation: crate::topology::Orientation::Forward, ..face })?;
    let u = loops.first()?.first()?;
    let ahead = u.at(0.5 + 1e-3) - u.at(0.5 - 1e-3);
    let uv = u.at(0.5);
    Some(surface.point(uv.x - 0.5 * ahead.y, uv.y + 0.5 * ahead.x))
}

/// What became of `solid`'s faces, edges and vertices in `result`: kept where they are in
/// it, modified into its faces on their surfaces (faces around the features), else deleted.
fn history(model: &Model, solid: Shape, removed: &[Shape], result: Shape) -> History {
    let mut history = History::new();
    let in_result: std::collections::HashSet<NodeId> = [NodeKind::Face, NodeKind::Edge, NodeKind::Vertex]
        .iter()
        .flat_map(|&k| model.unique_subshapes(result, k))
        .map(|s| s.node)
        .collect();
    let result_faces = model.unique_subshapes(result, NodeKind::Face);
    // (Edges below.)
    for kind in [NodeKind::Face, NodeKind::Vertex] {
        for s in model.unique_subshapes(solid, kind) {
            if in_result.contains(&s.node) {
                continue;
            }
            let images: Vec<NodeId> = if kind == NodeKind::Face && !removed.iter().any(|r| r.node == s.node) {
                let surface = model.surface(s);
                result_faces
                    .iter()
                    .filter(|f| surface.is_some() && model.surface(**f) == surface && face_point(model, s).is_some_and(|p| on_face(model, **f, p)))
                    .map(|f| f.node)
                    .collect()
            } else {
                Vec::new()
            };
            if images.is_empty() {
                history.record_deleted("removefeatures", &[s.node]);
            } else {
                history.record("removefeatures", s.node, &images);
            }
        }
    }
    // Edges: an edge of the solid's that is gone, modified into the result's edges along it
    // (the parts of it left, carried on); and the edges of a modified face reaching off the
    // solid's edges (made where features were, or carried on past them) generated from that
    // face (as OCCT's).
    let originals: Vec<Shape> = model.unique_subshapes(solid, NodeKind::Edge);
    let boxes: Vec<BoundingBox> = originals.iter().map(|&o| crate::properties::bounding_box(model, o).enlarged(1e-6)).collect();
    let result_edges: Vec<Shape> = model.unique_subshapes(result, NodeKind::Edge).into_iter().filter(|e| !model.is_degenerate(*e)).collect();
    let samples = |model: &Model, e: Shape| -> Vec<Point> {
        let Some((curve, first, last)) = model.edge_curve(e) else { return Vec::new() };
        (1..8).map(|k| crate::geometry::ParametricCurve::point(&curve, first + (last - first) * k as f64 / 8.0).transformed(&e.placement)).collect()
    };
    for &o in &originals {
        if in_result.contains(&o.node) {
            continue;
        }
        let along = samples(model, o);
        let near = BoundingBox::from_points(along.iter().copied()).enlarged(1e-6);
        let images: Vec<NodeId> = result_edges
            .iter()
            .filter(|e| !in_result_before(model, solid, e.node))
            .filter(|&&e| crate::properties::bounding_box(model, e).enlarged(1e-6).overlaps(&near))
            .filter(|&&e| along.iter().any(|&p| on_edge(model, p, e)) || samples(model, e).iter().any(|&p| on_edge(model, p, o)))
            .filter(|&&e| samples(model, e).iter().all(|&p| on_carrier(model, p, o)))
            .map(|e| e.node)
            .collect();
        if images.is_empty() {
            history.record_deleted("removefeatures", &[o.node]);
        } else {
            history.record("removefeatures", o.node, &images);
        }
    }
    // (Each result face's original, where it is a modified face's image.)
    let modified_from: Vec<(NodeId, NodeId)> = model
        .unique_subshapes(solid, NodeKind::Face)
        .into_iter()
        .filter_map(|f| history.modified(f.node).map(|images| images.iter().map(|&i| (i, f.node)).collect::<Vec<_>>()))
        .flatten()
        .collect();
    for f in model.unique_subshapes(solid, NodeKind::Face) {
        let Some(images) = history.modified(f.node).map(|i| i.to_vec()) else { continue };
        let mut made: Vec<NodeId> = Vec::new();
        for image in result_faces.iter().filter(|g| images.contains(&g.node)) {
            for e in model.unique_subshapes(*image, NodeKind::Edge) {
                if in_result_before(model, solid, e.node) || model.is_degenerate(e) || made.contains(&e.node) {
                    continue;
                }
                // (Between this face and another modified one: not a seam of its own.)
                let others = result_faces
                    .iter()
                    .filter(|g| !images.contains(&g.node) && model.unique_subshapes(**g, NodeKind::Edge).iter().any(|x| x.node == e.node))
                    .filter(|g| modified_from.iter().any(|(i, o)| *i == g.node && *o != f.node))
                    .count();
                if others == 0 {
                    continue;
                }
                if samples(model, e).iter().any(|&p| !originals.iter().zip(&boxes).any(|(&o, b)| b.contains(p) && on_edge(model, p, o))) {
                    made.push(e.node);
                }
            }
        }
        if !made.is_empty() {
            history.record_as("removefeatures", f.node, &made, crate::topology::ChangeKind::Generated);
        }
    }
    history
}

/// Whether `p` lies on `edge`.
fn on_edge(model: &Model, p: Point, edge: Shape) -> bool {
    let Some((curve, first, last)) = model.edge_curve(edge) else { return false };
    let local = p.transformed(&edge.placement.inverse());
    crate::geometry::extrema::nearest_on_curve(local, &curve, first.min(last), first.max(last)).distance <= 1e-6
}

/// Whether `p` lies on the curve `edge` runs along (anywhere along it).
fn on_carrier(model: &Model, p: Point, edge: Shape) -> bool {
    use crate::geometry::ParametricCurve;
    let Some((curve, first, last)) = model.edge_curve(edge) else { return false };
    // (A trimmed curve as the curve it trims.)
    let curve = match curve {
        crate::geometry::curve::Curve::Trimmed(t) => t.basis().clone(),
        c => c,
    };
    let local = p.transformed(&edge.placement.inverse());
    let span = (last - first).abs().max(1.0);
    let (lo, hi) = match (curve.is_periodic(), &curve) {
        (true, _) => (curve.first_parameter(), curve.last_parameter()),
        (false, crate::geometry::curve::Curve::Line(_)) => (first.min(last) - 1e6 * span, first.max(last) + 1e6 * span),
        _ => (curve.first_parameter().max(first.min(last) - span), curve.last_parameter().min(first.max(last) + span)),
    };
    crate::geometry::extrema::nearest_on_curve(local, &curve, lo, hi).distance <= 1e-6
}

/// Whether `node` is one of `solid`'s own (kept as it was).
fn in_result_before(model: &Model, solid: Shape, node: NodeId) -> bool {
    model.unique_subshapes(solid, node.kind).iter().any(|s| s.node == node)
}
