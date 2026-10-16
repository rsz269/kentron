//! Removing features of a solid: faces (a hole's, a boss's, a step's) taken away, the faces
//! around them carried on over where they were to close the solid again.
//!
//! The faces next to the removed ones are extended on their surfaces (a plane far across
//! the features, a cylinder or cone all round and along, a sphere or torus whole), and
//! with the removed faces they split space into cells. The gaps are the cells against the
//! removed faces lying within the features' neighbourhood (no farther than the faces
//! around them reach): one outside the solid is a hole or pocket, filled; one inside is a
//! boss, cut away. The faces lying on one surface are merged after, so that each face
//! around a feature is carried on over it as one face.

use crate::error::{Error, Result};
use crate::geometry::surface::Surface;
use crate::geometry::ParametricSurface;
use crate::math::{BoundingBox, Point};
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
    let size = (near.max - near.min).length().max(1e-9);
    // The faces around carried on, and the features: the cells they bound.
    let mut sheets: Vec<Shape> = removed.clone();
    for &f in &around {
        sheets.push(extended(model, f, &near)?);
    }
    let cells = crate::boolean::enclosed_solids(model, &sheets)?;
    let cells = model.unique_subshapes(cells, NodeKind::Solid);
    let within = near.enlarged(1e-3 * size);
    let classifier = crate::boolean::PointClassifier::new(model, solid)?;
    let (mut fill, mut cut) = (Vec::new(), Vec::new());
    for &cell in &cells {
        let b = crate::properties::bounding_box(model, cell);
        if !within.contains_box(&b) || !touches(model, cell, &removed) {
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
    for kind in [NodeKind::Face, NodeKind::Edge, NodeKind::Vertex] {
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
    history
}
