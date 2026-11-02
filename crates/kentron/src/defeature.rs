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

/// `shape` (a solid, or a compound or compsolid of solids) without the faces `faces`, the
/// faces around them carried on over where they were; several solids sharing faces where
/// they meet still.
#[doc(alias = "BRepAlgoAPI_Defeaturing")]
#[doc(alias = "removefeatures")]
pub fn remove_faces(model: &mut Model, shape: Shape, faces: &[Shape]) -> Result<Removed> {
    let solids = model.unique_subshapes(shape, NodeKind::Solid);
    if solids.is_empty() {
        return Err(Error::InvalidInput("features are removed from solids".into()));
    }
    let faces_of: Vec<Vec<Shape>> = solids.iter().map(|&s| model.unique_subshapes(s, NodeKind::Face)).collect();
    let all: Vec<Shape> = model.unique_subshapes(shape, NodeKind::Face);
    let removed: Vec<Shape> = all.iter().copied().filter(|f| faces.iter().any(|g| g.node == f.node)).collect();
    if removed.is_empty() {
        return Err(Error::InvalidInput("no faces of the solid to remove".into()));
    }
    // Each feature (faces joined by edges) on its own: the faces around it, sharing an edge
    // with it, carried on over it, and the cells they and it bound; the gaps of all filled
    // or cut together after (as OCCT's, which removes features one by one; carried on over
    // all the features at once, the faces around one cut through others). A gap inside a
    // solid is cut from it; one outside them all filled into the solid whose removed faces
    // it is against.
    let classifiers: Vec<crate::boolean::PointClassifier> = solids.iter().map(|&s| crate::boolean::PointClassifier::new(model, s)).collect::<Result<_>>()?;
    let mut fill: Vec<Vec<Shape>> = vec![Vec::new(); solids.len()];
    let mut cut: Vec<Vec<Shape>> = vec![Vec::new(); solids.len()];
    let mut around: Vec<Shape> = Vec::new();
    for feature in features(model, &removed) {
        let edges: std::collections::HashSet<NodeId> = feature.iter().flat_map(|&f| model.unique_subshapes(f, NodeKind::Edge)).map(|e| e.node).collect();
        let next: Vec<Shape> = all
            .iter()
            .copied()
            .filter(|f| !removed.iter().any(|g| g.node == f.node))
            .filter(|&f| model.unique_subshapes(f, NodeKind::Edge).iter().any(|e| edges.contains(&e.node)))
            .collect();
        // The neighbourhood: the box of the feature and the faces around it.
        let mut near = BoundingBox::EMPTY;
        for &f in feature.iter().chain(&next) {
            near.add_box(&crate::properties::bounding_box(model, f));
        }
        let mut sheets: Vec<Shape> = feature.clone();
        for &f in &next {
            sheets.push(extended(model, f, &near)?);
        }
        let cells = crate::boolean::enclosed_solids(model, &sheets)?;
        let cells = model.unique_subshapes(cells, NodeKind::Solid);
        let mut found = false;
        for &cell in &cells {
            if !touches(model, cell, &feature) {
                continue;
            }
            let Some(p) = crate::boolean::point_inside(model, cell)? else { continue };
            let locations: Vec<crate::boolean::PointLocation> = classifiers.iter().map(|c| c.locate(p)).collect();
            if let Some(k) = locations.iter().position(|&l| l == crate::boolean::PointLocation::Inside) {
                cut[k].push(cell);
                found = true;
            } else if locations.iter().all(|&l| l == crate::boolean::PointLocation::Outside) {
                let mine = |k: usize| feature.iter().copied().filter(|f| faces_of[k].iter().any(|g| g.node == f.node)).collect::<Vec<_>>();
                let k = (0..solids.len()).find(|&k| touches(model, cell, &mine(k))).unwrap_or(0);
                fill[k].push(cell);
                found = true;
            }
        }
        if !found {
            return Err(Error::Failed("removing features: the faces around them do not close over them".into()));
        }
        for f in next {
            if !around.iter().any(|a| a.node == f.node) {
                around.push(f);
            }
        }
    }
    let mut results = Vec::with_capacity(solids.len());
    for (k, &solid) in solids.iter().enumerate() {
        if fill[k].is_empty() && cut[k].is_empty() {
            results.push(solid);
            continue;
        }
        results.push(closed_over(model, solid, &fill[k], &cut[k], &around)?);
    }
    // Several solids: made to share their faces again where they meet (split where one
    // meets part of another's), in a shape of the input's kind.
    let (result, solid_images) = if let [one] = results[..] {
        (one, vec![one])
    } else {
        let cells = crate::boolean::cells(model, &results)?;
        let mut images = Vec::with_capacity(results.len());
        for k in 0..results.len() {
            let mine: Vec<Shape> = cells.iter().filter(|c| c.inside[k] && c.inside.iter().filter(|&&i| i).count() == 1).map(|c| c.shape).collect();
            let [solid] = mine[..] else {
                return Err(Error::Failed("removing features: the solids overlap after".into()));
            };
            images.push(model.unique_subshapes(solid, NodeKind::Solid).first().copied().unwrap_or(solid));
        }
        let kind = if shape.kind() == NodeKind::CompSolid { NodeKind::CompSolid } else { NodeKind::Compound };
        let out = model.new_empty(kind);
        for &solid in &images {
            model.add(out, solid)?;
        }
        (out, images)
    };
    let mut history = history(model, shape, &removed, result);
    if solids.len() > 1 {
        for (&solid, &image) in solids.iter().zip(&solid_images) {
            if image.node != solid.node {
                history.record("removefeatures", solid.node, &[image.node]);
            }
        }
    }
    Ok(Removed { shape: result, history })
}

/// `solid` with the gaps `fill` filled and `cut` cut away, the faces on one surface merged as
/// OCCT's removal merges them.
fn closed_over(model: &mut Model, solid: Shape, fill: &[Shape], cut: &[Shape], around: &[Shape]) -> Result<Shape> {
    let mut result = solid;
    for &cell in fill {
        result = crate::boolean::union(model, result, cell)?;
    }
    for &cell in cut {
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
    crate::unify::unify_same_domain_with(model, result, options)
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

/// `faces` in features: those joined by edges together.
fn features(model: &Model, faces: &[Shape]) -> Vec<Vec<Shape>> {
    let edges: Vec<std::collections::HashSet<NodeId>> = faces.iter().map(|&f| model.unique_subshapes(f, NodeKind::Edge).into_iter().map(|e| e.node).collect()).collect();
    let mut feature_of: Vec<Option<usize>> = vec![None; faces.len()];
    let mut out: Vec<Vec<Shape>> = Vec::new();
    for start in 0..faces.len() {
        if feature_of[start].is_some() {
            continue;
        }
        let id = out.len();
        feature_of[start] = Some(id);
        let mut stack = vec![start];
        let mut members = Vec::new();
        while let Some(i) = stack.pop() {
            members.push(faces[i]);
            for j in 0..faces.len() {
                if feature_of[j].is_none() && !edges[i].is_disjoint(&edges[j]) {
                    feature_of[j] = Some(id);
                    stack.push(j);
                }
            }
        }
        out.push(members);
    }
    out
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
    let original_faces = model.unique_subshapes(solid, NodeKind::Face);
    // (Edges below.)
    for kind in [NodeKind::Face, NodeKind::Vertex] {
        for s in model.unique_subshapes(solid, kind) {
            if in_result.contains(&s.node) {
                continue;
            }
            let images: Vec<NodeId> = if kind == NodeKind::Face && !removed.iter().any(|r| r.node == s.node) {
                let surface = model.surface(s);
                let alike: Vec<Shape> = result_faces.iter().copied().filter(|f| surface.is_some() && model.surface(*f) == surface).collect();
                let mut images: Vec<Shape> = alike.iter().copied().filter(|&f| face_point(model, s).is_some_and(|p| on_face(model, f, p))).collect();
                // (And the pieces of it carried on that another solid split off, where it
                // meets that one in part: new faces on its surface holding none of the
                // originals, next to its images.)
                let fresh: Vec<Shape> = alike
                    .iter()
                    .copied()
                    .filter(|&f| !images.iter().any(|i| i.node == f.node) && !in_result_before(model, solid, f.node))
                    .filter(|&f| !original_faces.iter().any(|&o| model.surface(o) == surface && face_point(model, o).is_some_and(|p| on_face(model, f, p))))
                    .collect();
                if !images.is_empty() {
                    let mut grown = true;
                    while grown {
                        grown = false;
                        for &f in &fresh {
                            if images.iter().any(|i| i.node == f.node) {
                                continue;
                            }
                            let edges: Vec<NodeId> = model.unique_subshapes(f, NodeKind::Edge).into_iter().map(|e| e.node).collect();
                            if images.iter().any(|&i| model.unique_subshapes(i, NodeKind::Edge).iter().any(|e| edges.contains(&e.node))) {
                                images.push(f);
                                grown = true;
                            }
                        }
                    }
                }
                images.into_iter().map(|f| f.node).collect()
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
