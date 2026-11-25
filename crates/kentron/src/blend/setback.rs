//! Fillets set back from the corners where three of them meet.
//!
//! OCCT's fillet has no set-back corner (`ChFi3d` closes three fillets meeting at a vertex
//! where their contact lines cross, with no stop of its own), so this is built from the
//! definition with kentron's own rolling-ball fillet.
//!
//! The rolling-ball fillet of radius r of a convex corner is the union of the balls of r
//! inside the solid: the core the ball's centre may reach (the solid shrunk by r) swelled
//! again by r. Round a set-back corner the core is cut by one plane more, through the three
//! points on its edges whose feet on the solid's edges lie `setback` from the corner. The
//! core stays convex, so the core swelled by r has a tangent plane everywhere: the boundary
//! leaves no crease. Past the cut, the core is the ordinary one, so each fillet is the
//! ordinary fillet farther than `setback` from the corner, and it stops on its section there
//! (where the ball's centre reaches the cut). On each face, the flat stops along the cut's
//! trace on the core's face, short of the point where its two fillets' contact lines would
//! cross. The corner closes with a flat (the cut swelled by r), a round along each edge of the
//! cut and a sphere at each of its corners, all meeting with their tangent planes.
//!
//! The core cut so is the core of the solid cut by the plane swelled out from the core's cut
//! by r (the solid with its corner chamfered, the chamfer's sides on its edges), so the result
//! is that chamfered solid with all its edges there rounded by the ordinary fillet of r.

use crate::error::{Error, Result};
use crate::math::{Point, Vector};
use crate::topology::{Model, NodeId, NodeKind, Orientation, Shape};
use crate::{boolean, tolerance};

use super::{edge_turn, fillet, line_of, tangent_chains, volume_of};

/// `shape` with `edges` rounded by a rolling ball of `radius`, as [`fillet`] rounds them,
/// except at every corner where three of the edges (with the edges continuing them
/// tangentially) meet: there each fillet stops `setback` short of the corner, and the corner
/// is closed by faces meeting the three fillets and the three faces with their tangent planes
/// (see the module's notes).
///
/// `radius` must be positive and `setback` greater than `radius`. Supported at corners of
/// three planes meeting square to one another along straight convex edges, each edge long
/// enough to hold its set-back stops (and the fillets otherwise as [`fillet`] supports them).
pub fn fillet_setback(model: &mut Model, shape: Shape, edges: &[Shape], radius: f64, setback: f64) -> Result<Shape> {
    if !(radius > 0.0) || !radius.is_finite() {
        return Err(Error::InvalidInput(format!("a fillet radius of {radius}")));
    }
    if !(setback > radius) || !setback.is_finite() {
        return Err(Error::InvalidInput(format!("a setback of {setback} not greater than the radius {radius}")));
    }
    if edges.is_empty() {
        return Err(Error::InvalidInput("no edges to fillet".into()));
    }
    // The stop on each edge: how far along it from the corner the chamfer meets it.
    let reach = setback + (2.0 - 3f64.sqrt()) * radius;
    let corners = set_back_corners(model, shape, edges, radius, reach)?;
    if corners.is_empty() {
        return fillet(model, shape, edges, radius);
    }
    // Each corner chamfered: cut by the plane through the points `reach` along its edges,
    // with a tetrahedron whose other faces stand off the solid's faces outside it.
    let mut cut = shape;
    let before = volume_of(model, shape);
    for &(vertex, along) in &corners {
        let out = reach;
        let apex = vertex - (along[0] + along[1] + along[2]) * out;
        let base = along.map(|e| apex + e * (reach + 3.0 * out));
        let tool = super::tetrahedron(model, [apex, base[0], base[1], base[2]])?;
        cut = boolean::difference(model, cut, tool)?;
    }
    // The corners' tetrahedra, and nothing else, taken away.
    let removed = before - volume_of(model, cut);
    let expected = corners.len() as f64 * reach.powi(3) / 6.0;
    if (removed - expected).abs() > 1e-6 * expected.max(before) {
        return Err(Error::Failed("a set-back corner's chamfer cut more than the corner".into()));
    }
    // The named edges, as the chamfers left them, and the chamfers' edges.
    let mut rounded: Vec<Shape> = Vec::new();
    for &edge in edges {
        let (curve, first, last) = model.edge_curve(Shape { orientation: Orientation::Forward, ..edge }).ok_or_else(|| Error::InvalidInput("an edge without a curve".into()))?;
        let found = super::find_edge(model, cut, &curve, first, last).ok_or_else(|| Error::Failed("an edge lost to a set-back corner's chamfer".into()))?;
        rounded.push(found);
    }
    for &(vertex, along) in &corners {
        let stops = along.map(|e| vertex + e * reach);
        for i in 0..3 {
            let (a, b) = (stops[i], stops[(i + 1) % 3]);
            let tol = 1e-6 * reach;
            let side = model
                .explore(cut, NodeKind::Edge)
                .find(|&e| {
                    let (p, q) = model.edge_vertices(e);
                    let (p, q) = (model.point(p), model.point(q));
                    (p.distance(a) <= tol && q.distance(b) <= tol) || (p.distance(b) <= tol && q.distance(a) <= tol)
                })
                .ok_or_else(|| Error::Failed("a set-back corner's chamfer edge not found".into()))?;
            rounded.push(Shape { orientation: Orientation::Forward, ..side });
        }
    }
    fillet(model, cut, &rounded, radius)
}

/// The corners where three of `edges` (with the edges continuing them tangentially) meet,
/// each its vertex and the unit directions of its three edges away from it: refused unless
/// every such corner is of three planes square to one another along straight convex edges,
/// each edge long enough to hold the stops `reach` from the corners at its ends.
pub(super) fn set_back_corners(model: &Model, shape: Shape, edges: &[Shape], radius: f64, reach: f64) -> Result<Vec<(Point, [Vector; 3])>> {
    let pairs: Vec<(Shape, f64)> = edges.iter().map(|&e| (e, radius)).collect();
    let chains: Vec<Shape> = tangent_chains(model, shape, &pairs, 0.0).into_iter().map(|e| e.0).collect();
    // The vertices where three edges of the chains meet, with those edges.
    let mut at: Vec<(NodeId, Point, Vec<Shape>)> = Vec::new();
    for &edge in &chains {
        let (a, b) = model.edge_vertices(edge);
        for v in [a, b] {
            let p = model.point(v);
            match at.iter_mut().find(|(n, _, _)| *n == v.node) {
                Some((_, _, list)) => {
                    if !list.iter().any(|e| e.node == edge.node) {
                        list.push(edge);
                    }
                }
                None => at.push((v.node, p, vec![edge])),
            }
        }
    }
    let mut corners: Vec<(Point, [Vector; 3])> = Vec::new();
    // How much of each edge the corners at its ends take.
    let mut taken: Vec<(NodeId, f64, f64)> = Vec::new();
    for (_, vertex, list) in at.iter().filter(|c| c.2.len() >= 3) {
        if list.len() > 3 {
            return Err(Error::Failed("a set-back corner of more than three edges".into()));
        }
        let mut along = [Vector::new(0.0, 0.0, 0.0); 3];
        for (k, &edge) in list.iter().enumerate() {
            let (curve, first, last) = model.edge_curve(Shape { orientation: Orientation::Forward, ..edge }).ok_or_else(|| Error::InvalidInput("an edge without a curve".into()))?;
            if line_of(&curve).is_none() {
                return Err(Error::Failed("a set-back corner on an edge that is not straight".into()));
            }
            use crate::geometry::ParametricCurve;
            let (p, q) = (curve.point(first), curve.point(last));
            let (near, far) = if p.distance(*vertex) <= q.distance(*vertex) { (p, q) } else { (q, p) };
            let length = near.distance(far);
            along[k] = (far - near).normalized();
            match taken.iter_mut().find(|t| t.0 == edge.node) {
                Some(t) => t.2 += reach,
                None => taken.push((edge.node, length, reach)),
            }
            if edge_turn(model, shape, edge) != Some(true) {
                return Err(Error::Failed("a set-back corner on an edge that is not convex".into()));
            }
        }
        for i in 0..3 {
            if along[i].dot(along[(i + 1) % 3]).abs() > 1e-9 {
                return Err(Error::Failed("a set-back corner whose edges are not square to one another".into()));
            }
        }
        // Its three faces, planes.
        let mut faces: Vec<NodeId> = Vec::new();
        for f in model.explore(shape, NodeKind::Face) {
            if model.explore(f, NodeKind::Edge).any(|e| list.iter().any(|l| l.node == e.node)) && !faces.contains(&f.node) {
                if !matches!(model.surface(f), Some(crate::geometry::surface::Surface::Plane(_))) {
                    return Err(Error::Failed("a set-back corner on a face that is not a plane".into()));
                }
                faces.push(f.node);
            }
        }
        if faces.len() != 3 {
            return Err(Error::Failed("a set-back corner not of three faces".into()));
        }
        corners.push((*vertex, along));
    }
    for &(_, length, stops) in &taken {
        if !(length > stops + 10.0 * tolerance::LINEAR) {
            return Err(Error::InvalidInput("an edge too short for its set-back stops".into()));
        }
    }
    Ok(corners)
}
