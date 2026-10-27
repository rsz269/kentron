//! Unifying faces and edges that lie on the same geometry: what is left of a model built by
//! booleans, whose faces get split along the other operand's edges.
//!
//! In each shell, neighbouring faces on the same plane or quadric (and facing the same way)
//! become one face, the edges between them dropped; then edges meeting at a vertex of no
//! other edge, on one line or circle, become one edge. Faces that would close round a
//! periodic surface without a seam between them are left apart.

use crate::hashing::{FastMap as HashMap, FastSet as HashSet};

use crate::error::Result;
use crate::geometry::curve::{Circle, Curve, Line};
use crate::geometry::curve2::Curve2;
use crate::geometry::project::analytic_pcurve;
use crate::geometry::surface::Surface;
use crate::geometry::{ParametricCurve, ParametricSurface};
use crate::math::{Direction, Point, Point2, Transform, Transform2, Transformable2, Vector, Vector2};
use crate::tolerance;
use crate::topology::{EdgeId, Model, NodeKind, Orientation, Shape, VertexId, WireId};

/// What [`unify_same_domain_with`] merges.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct UnifyOptions {
    #[doc(alias = "UnifyFaces")]
    pub faces: bool,
    #[doc(alias = "UnifyEdges")]
    pub edges: bool,
    /// Faces merged only where they lie on the very same surface (pieces of one face split
    /// apart), not on others with the same geometry.
    pub same_surface_only: bool,
}

impl Default for UnifyOptions {
    fn default() -> Self {
        Self { faces: true, edges: true, same_surface_only: false }
    }
}

/// `shape` with faces on the same surfaces, and edges on the same curves, merged. Should
/// the merged shape come out invalid where `shape` is valid, `shape` is returned.
#[doc(alias = "ShapeUpgrade_UnifySameDomain")]
#[doc(alias = "unifysamedom")]
#[doc(alias = "clean")]
pub fn unify_same_domain(model: &mut Model, shape: Shape) -> Result<Shape> {
    unify_same_domain_with(model, shape, UnifyOptions::default())
}

/// [`unify_same_domain`], merging only faces or only edges if `options` say so.
pub fn unify_same_domain_with(model: &mut Model, shape: Shape, options: UnifyOptions) -> Result<Shape> {
    let local = Shape { orientation: Orientation::Forward, placement: Transform::IDENTITY, ..shape };
    let result = rebuild(model, local, options)?;
    // A valid shape is not traded for an invalid one (OCCT's safe mode).
    if result.node != local.node && !crate::validity::check(model, result).is_valid() && crate::validity::check(model, local).is_valid() {
        return Ok(shape);
    }
    Ok(Shape { orientation: shape.orientation, placement: shape.placement, ..result })
}

fn rebuild(model: &mut Model, shape: Shape, options: UnifyOptions) -> Result<Shape> {
    match shape.kind() {
        NodeKind::Shell => unify_shell(model, shape, options),
        NodeKind::Compound | NodeKind::CompSolid | NodeKind::Solid => {
            let mut children: Vec<Shape> = model.children_local(shape).collect();
            // A compound's own faces are unified together, as if in one shell.
            let faces: Vec<Shape> = children.iter().copied().filter(|c| c.kind() == NodeKind::Face && c.placement == Transform::IDENTITY).collect();
            if shape.kind() == NodeKind::Compound && !faces.is_empty() {
                let together = model.new_empty(NodeKind::Shell);
                for &f in &faces {
                    model.add(together, f)?;
                }
                let unified = unify_shell(model, together, options)?;
                if unified.node != together.node {
                    let out = model.new_empty(NodeKind::Compound);
                    for c in children.iter().filter(|c| !faces.contains(c)) {
                        model.add(out, *c)?;
                    }
                    for f in model.children_local(unified).collect::<Vec<_>>() {
                        model.add(out, f)?;
                    }
                    return rebuild(model, out, options);
                }
                children = model.children_local(shape).collect();
            }
            let mut rebuilt = Vec::with_capacity(children.len());
            for child in &children {
                let local = Shape { orientation: Orientation::Forward, placement: Transform::IDENTITY, ..*child };
                rebuilt.push(rebuild(model, local, options)?);
            }
            if rebuilt.iter().zip(&children).all(|(r, c)| r.node == c.node) {
                return Ok(shape);
            }
            let out = model.new_empty(shape.kind());
            for (r, c) in rebuilt.into_iter().zip(&children) {
                model.add(out, Shape { orientation: c.orientation, placement: c.placement, ..r })?;
            }
            Ok(out)
        }
        // A face alone: its edges merged, as in a shell of it alone.
        NodeKind::Face => {
            let alone = model.new_empty(NodeKind::Shell);
            model.add(alone, shape)?;
            let unified = unify_shell(model, alone, options)?;
            match model.children_local(unified).collect::<Vec<_>>()[..] {
                [face] if unified.node != alone.node => Ok(face),
                _ => Ok(shape),
            }
        }
        NodeKind::Wire if options.edges => unify_wire(model, shape),
        _ => Ok(shape),
    }
}

/// A wire alone with its edges merged where consecutive ones lie on one line or one circle
/// and meet at a vertex of no other edge of it: lines into the segment between the group's
/// ends, arcs into the arc of the circle they run round together.
fn unify_wire(model: &mut Model, wire: Shape) -> Result<Shape> {
    let edges = model.ordered_edges(wire);
    if edges.len() < 2 || edges.iter().any(|e| e.placement != Transform::IDENTITY) {
        return Ok(wire);
    }
    let untrimmed = |(c, a, b): (Curve, f64, f64)| match c {
        Curve::Trimmed(t) if matches!(t.basis(), Curve::Line(_) | Curve::Circle(_)) => (t.basis().clone(), a, b),
        c => (c, a, b),
    };
    let Some(curves) = edges.iter().map(|&e| model.edge_curve(Shape { orientation: Orientation::Forward, ..e }).map(untrimmed)).collect::<Option<Vec<_>>>() else { return Ok(wire) };
    let ends: Vec<(Shape, Shape)> = edges.iter().map(|&e| model.edge_vertices(e)).collect();
    let mut uses: HashMap<crate::topology::NodeId, usize> = HashMap::default();
    for (a, b) in &ends {
        *uses.entry(a.node).or_default() += 1;
        *uses.entry(b.node).or_default() += 1;
    }
    // (Each edge joins the one before it where they share a carrier and the vertex between
    // them is theirs alone.)
    let joins = |i: usize, j: usize| ends[i].1.node == ends[j].0.node && uses[&ends[i].1.node] == 2 && same_carrier(&curves[i], &curves[j]);
    let n = edges.len();
    let mut groups: Vec<Vec<usize>> = Vec::new();
    for i in 0..n {
        match groups.last_mut() {
            Some(g) if joins(*g.last().expect("a group"), i) => g.push(i),
            _ => groups.push(vec![i]),
        }
    }
    // (Round a closed wire: the last group goes on into the first.)
    if groups.len() > 1 && joins(n - 1, 0) {
        let last = groups.pop().expect("groups");
        groups[0].splice(0..0, last);
    }
    if groups.iter().all(|g| g.len() == 1) {
        return Ok(wire);
    }
    let vertex = |s: Shape| VertexId::try_from(s.node).map_err(|n| crate::error::Error::InvalidInput(format!("{n:?} is not a vertex")));
    let mut out = Vec::with_capacity(groups.len());
    for g in &groups {
        let (first, last) = (g[0], g[g.len() - 1]);
        if g.len() == 1 {
            out.push(edges[first]);
            continue;
        }
        let (start, end) = (ends[first].0, ends[last].1);
        let (p, q) = (model.point(start), model.point(end));
        let tolerance = g.iter().map(|&i| model.tolerance(edges[i])).fold(tolerance::LINEAR, f64::max);
        let (curve, from, to, forward) = match &curves[first].0 {
            Curve::Line(_) => {
                let Some(direction) = Direction::try_from_vector(q - p) else { return Ok(wire) };
                (Curve::Line(Line::new(p, direction)), 0.0, p.distance(q), true)
            }
            Curve::Circle(c) => {
                // The arcs' angles together, from the group's start the way the first runs.
                let angle = |x: Point| {
                    let d = x - c.center();
                    d.dot(Vector::from(c.frame().y())).atan2(d.dot(Vector::from(c.frame().x())))
                };
                let sweep: f64 = g.iter().map(|&i| (curves[i].2 - curves[i].1).abs()).sum();
                let along = edges[first].orientation == Orientation::Forward;
                let a = angle(p);
                if along { (curves[first].0.clone(), a, a + sweep, true) } else { (curves[first].0.clone(), a - sweep, a, false) }
            }
            _ => return Ok(wire),
        };
        let (v0, v1) = if forward { (vertex(start)?, vertex(end)?) } else { (vertex(end)?, vertex(start)?) };
        let id = model.add_edge(v0, v1, tolerance)?;
        model.set_edge_curve(id.into(), curve, tolerance)?;
        model.set_edge_range(id.into(), from, to)?;
        let merged = Shape::from(id);
        out.push(if forward { merged } else { merged.reversed() });
    }
    crate::make::wire(model, &out)
}

/// An edge's use in a face's loop, relative to the face.
#[derive(Clone)]
struct Use {
    edge: EdgeId,
    orientation: Orientation,
    pcurve: Option<Curve2>,
    first: f64,
    last: f64,
}

/// A face to be: the original one while unchanged.
struct Spec {
    original: Option<Shape>,
    /// In the shell.
    orientation: Orientation,
    surface: Surface,
    tolerance: f64,
    loops: Vec<Vec<Use>>,
}

fn unify_shell(model: &mut Model, shell: Shape, options: UnifyOptions) -> Result<Shape> {
    let mut specs: Vec<Spec> = Vec::new();
    for face in model.children_local(shell).collect::<Vec<_>>() {
        let node = Shape { orientation: Orientation::Forward, ..face };
        let Some(surface) = model.surface(node).cloned() else { return Ok(shell) };
        let mut loops = Vec::new();
        for wire in model.children_local(node).collect::<Vec<_>>() {
            let id = WireId::try_from(wire.node).expect("a face holds wires");
            let mut uses: Vec<Use> = model.get(id).expect("a wire").coedges.iter().map(|&c| {
                let record = model.get(c).expect("a coedge");
                let (pcurve, first, last) = match model.effective_pcurve(c) {
                    Some((p, a, b)) => (Some(p), a, b),
                    None => (None, 0.0, 0.0),
                };
                Use { edge: record.edge, orientation: record.orientation, pcurve, first, last }
            }).collect();
            if wire.orientation == Orientation::Reversed {
                uses.reverse();
                for u in &mut uses {
                    u.orientation = u.orientation.reversed();
                }
            }
            loops.push(uses);
        }
        specs.push(Spec { original: Some(node), orientation: face.orientation, tolerance: model.tolerance(node), surface, loops });
    }
    let mut changed = options.faces && merge_faces(model, &mut specs, options.same_surface_only);
    if options.edges {
        changed |= merge_edges(model, &mut specs)?;
    }
    if !changed {
        return Ok(shell);
    }
    let out = model.new_empty(NodeKind::Shell);
    for spec in specs {
        let face = match spec.original {
            Some(face) => face,
            None => build_face(model, &spec)?,
        };
        model.add(out, Shape { orientation: spec.orientation, ..face })?;
    }
    Ok(out)
}

fn build_face(model: &mut Model, spec: &Spec) -> Result<Shape> {
    let face = model.new_empty(NodeKind::Face);
    model.set_face_surface(face, spec.surface.clone(), spec.tolerance)?;
    // The outer loop (the largest) first.
    let mut order: Vec<usize> = (0..spec.loops.len()).collect();
    order.sort_by(|&a, &b| loop_area(&spec.loops[b]).abs().total_cmp(&loop_area(&spec.loops[a]).abs()));
    for i in order {
        let mut coedges = Vec::new();
        for u in &spec.loops[i] {
            let coedge = model.add_coedge(u.edge, u.orientation)?;
            if let Some(p) = &u.pcurve {
                model.set_pcurve_with_range(coedge, p.clone(), u.first, u.last)?;
            }
            coedges.push(coedge);
        }
        let wire = model.add_wire(&coedges)?;
        model.add(face, wire.into())?;
    }
    Ok(face)
}

/// The ends of a use, in the direction it runs.
fn ends(model: &Model, edge: EdgeId, orientation: Orientation) -> (VertexId, VertexId) {
    let (a, b) = model.edge_vertices(Shape::from(edge));
    let (a, b) = (VertexId::try_from(a.node).expect("a vertex"), VertexId::try_from(b.node).expect("a vertex"));
    if orientation == Orientation::Reversed { (b, a) } else { (a, b) }
}

/// A use's (u, v) at its start and end.
fn uv_ends(u: &Use) -> Option<(Point2, Point2)> {
    let p = u.pcurve.as_ref()?;
    let (s, e) = if u.orientation == Orientation::Reversed { (u.last, u.first) } else { (u.first, u.last) };
    Some((p.point(s), p.point(e)))
}

fn loop_area(uses: &[Use]) -> f64 {
    let mut points = Vec::new();
    for u in uses {
        let Some(p) = &u.pcurve else { continue };
        let (a, b) = if u.orientation == Orientation::Reversed { (u.last, u.first) } else { (u.first, u.last) };
        for j in 0..8 {
            points.push(p.point(a + (b - a) * j as f64 / 8.0));
        }
    }
    let n = points.len();
    (0..n).map(|i| points[i].x * points[(i + 1) % n].y - points[(i + 1) % n].x * points[i].y).sum::<f64>() / 2.0
}

/// True for planes and quadrics: surfaces whose edges have closed-form pcurves.
fn elementary(s: &Surface) -> bool {
    matches!(s, Surface::Plane(_) | Surface::Cylinder(_) | Surface::Cone(_) | Surface::Sphere(_))
}

/// True if two elementary surfaces are the same point set.
fn same_domain(a: &Surface, b: &Surface) -> bool {
    let tol = tolerance::LINEAR;
    let parallel = |x: Direction, y: Direction| Vector::from(x).cross(Vector::from(y)).length() <= 1e-9;
    let on_axis = |origin: Point, z: Direction, p: Point| {
        let d = p - origin;
        let z = Vector::from(z);
        (d - z * d.dot(z)).length() <= tol
    };
    match (a, b) {
        (Surface::Plane(p), Surface::Plane(q)) => {
            parallel(p.normal(), q.normal()) && (q.origin() - p.origin()).dot(Vector::from(p.normal())).abs() <= tol
        }
        (Surface::Cylinder(p), Surface::Cylinder(q)) => {
            parallel(p.frame().z(), q.frame().z())
                && on_axis(p.frame().origin(), p.frame().z(), q.frame().origin())
                && (p.radius() - q.radius()).abs() <= tol
        }
        (Surface::Cone(p), Surface::Cone(q)) => {
            let (zp, zq) = (Vector::from(p.frame().z()), Vector::from(q.frame().z()));
            let same_way = zp.dot(zq) > 0.0;
            let angle = if same_way { q.semi_angle() } else { -q.semi_angle() };
            parallel(p.frame().z(), q.frame().z())
                && p.apex().distance(q.apex()) <= tol
                && (p.semi_angle() - angle).abs() <= 1e-9
        }
        (Surface::Sphere(p), Surface::Sphere(q)) => {
            p.frame().origin().distance(q.frame().origin()) <= tol && (p.radius() - q.radius()).abs() <= tol
        }
        _ => false,
    }
}

/// The unit normal `∂u × ∂v` of an elementary surface at a point on it.
fn normal(surface: &Surface, p: Point) -> Option<Vector> {
    let (u, v) = surface.parameters(p)?;
    let (_, du, dv) = surface.d1(u, v);
    let n = du.cross(dv);
    (n.length() > 1e-12).then(|| n.normalized())
}

fn sign(o: Orientation) -> f64 {
    if o == Orientation::Reversed { -1.0 } else { 1.0 }
}

/// Merges groups of neighbouring faces on the same surface. True if any were.
fn merge_faces(model: &Model, specs: &mut Vec<Spec>, same_surface_only: bool) -> bool {
    // Where each edge is used.
    let mut uses_of: HashMap<EdgeId, Vec<usize>> = HashMap::default();
    for (fi, spec) in specs.iter().enumerate() {
        for l in &spec.loops {
            for u in l {
                uses_of.entry(u.edge).or_default().push(fi);
            }
        }
    }
    let mut parent: Vec<usize> = (0..specs.len()).collect();
    fn find(parent: &mut [usize], i: usize) -> usize {
        let mut r = i;
        while parent[r] != r {
            r = parent[r];
        }
        parent[i] = r;
        r
    }
    let mut internal: HashSet<EdgeId> = HashSet::default();
    for (&edge, faces) in &uses_of {
        let [f, g] = faces[..] else { continue };
        // (Pieces of one face, on the very same surface, merge whatever the surface: their
        // pcurves are on one parametrization.)
        let pieces = same_surface_only && specs[f].surface == specs[g].surface;
        if f == g || (!pieces && (!elementary(&specs[f].surface) || !same_domain(&specs[f].surface, &specs[g].surface))) {
            continue;
        }
        if same_surface_only && !pieces {
            continue;
        }
        let Some((curve, first, last)) = model.edge_curve(Shape::from(edge)) else { continue };
        let p = curve.point(0.5 * (first + last));
        let (Some(nf), Some(ng)) = (normal(&specs[f].surface, p), normal(&specs[g].surface, p)) else { continue };
        if (nf * sign(specs[f].orientation)).dot(ng * sign(specs[g].orientation)) <= 0.0 {
            continue;
        }
        internal.insert(edge);
        let (a, b) = (find(&mut parent, f), find(&mut parent, g));
        parent[a] = b;
    }
    let mut groups: HashMap<usize, Vec<usize>> = HashMap::default();
    for i in 0..specs.len() {
        let r = find(&mut parent, i);
        groups.entry(r).or_default().push(i);
    }
    let mut merged: Vec<(Vec<usize>, Spec)> = Vec::new();
    for group in groups.into_values().filter(|g| g.len() > 1) {
        if let Some(spec) = merge_group(model, specs, &group, &internal) {
            merged.push((group, spec));
        }
    }
    if merged.is_empty() {
        return false;
    }
    let gone: HashSet<usize> = merged.iter().flat_map(|(g, _)| g.iter().copied()).collect();
    let mut kept: Vec<Spec> = std::mem::take(specs).into_iter().enumerate().filter(|(i, _)| !gone.contains(i)).map(|(_, s)| s).collect();
    kept.extend(merged.into_iter().map(|(_, s)| s));
    *specs = kept;
    true
}

/// One face for a group, on its first face's surface; `None` where that does not work out
/// (degenerate edges, no closed-form pcurves, a loop closing round a period).
fn merge_group(model: &Model, specs: &[Spec], group: &[usize], internal: &HashSet<EdgeId>) -> Option<Spec> {
    let base = &specs[group[0]];
    let surface = base.surface.clone();
    let orientation = base.orientation;
    // The new face's loops: the group's loops spliced where they meet, each followed in its
    // own order (backwards in a face turned the other way in the shell) until an edge between
    // the group's faces, then on from that edge's other use (seams, used twice by one face,
    // are kept).
    let after = |s: usize, l: usize, i: usize| {
        let n = specs[s].loops[l].len();
        (s, l, if specs[s].orientation == orientation { (i + 1) % n } else { (i + n - 1) % n })
    };
    let mut twin: HashMap<(usize, usize, usize), (usize, usize, usize)> = HashMap::default();
    let mut by_edge: HashMap<EdgeId, Vec<(usize, usize, usize)>> = HashMap::default();
    for &s in group {
        for (l, uses) in specs[s].loops.iter().enumerate() {
            for (i, u) in uses.iter().enumerate() {
                if internal.contains(&u.edge) {
                    by_edge.entry(u.edge).or_default().push((s, l, i));
                }
            }
        }
    }
    for places in by_edge.values() {
        let [a, b] = places[..] else { return None };
        twin.insert(a, b);
        twin.insert(b, a);
    }
    let mut done: HashSet<(usize, usize, usize)> = HashSet::default();
    let mut loops = Vec::new();
    for &s in group {
        for (l, uses) in specs[s].loops.iter().enumerate() {
            for i in 0..uses.len() {
                let start = (s, l, i);
                if done.contains(&start) || twin.contains_key(&start) {
                    continue;
                }
                let mut chain: Vec<Use> = Vec::new();
                let mut at = start;
                loop {
                    if !done.insert(at) {
                        // Back where it started, or (inconsistent loops) elsewhere.
                        if at != start {
                            return None;
                        }
                        break;
                    }
                    let (fs, fl, fi) = at;
                    let spec = &specs[fs];
                    let u = &spec.loops[fl][fi];
                    let (curve, first, last) = model.edge_curve(Shape::from(u.edge))?;
                    let in_shell = u.orientation.compose(spec.orientation);
                    // A closed-form pcurve, else (a section traced between faces of boolean
                    // operands) the face's own, on the very surface.
                    let (pcurve, first, last) = match analytic_pcurve(&curve, &surface) {
                        Some(p) => (p, first, last),
                        None if spec.surface == surface => (u.pcurve.clone()?, u.first, u.last),
                        None => (u.pcurve.as_ref()?.transformed(&reparametrization(&spec.surface, &surface)?), u.first, u.last),
                    };
                    chain.push(Use { edge: u.edge, orientation: in_shell.compose(orientation), pcurve: Some(pcurve), first, last });
                    // The next use, across edges between the group's faces.
                    let mut next = after(fs, fl, fi);
                    let mut steps = 0;
                    while let Some(&other) = twin.get(&next) {
                        done.insert(next);
                        done.insert(other);
                        next = after(other.0, other.1, other.2);
                        steps += 1;
                        if steps > 4 * group.len() + 64 {
                            return None;
                        }
                    }
                    at = next;
                }
                if !align(&mut chain, &surface) {
                    return None;
                }
                loops.push(chain);
            }
        }
    }
    gather(&mut loops, &surface);
    let tolerance = group.iter().map(|&i| specs[i].tolerance).fold(0.0, f64::max);
    Some(Spec { original: None, orientation, surface, tolerance, loops })
}

/// The map from `from`'s parameters to `to`'s, for a plane or cylinder on the very same
/// point set, framed otherwise.
fn reparametrization(from: &Surface, to: &Surface) -> Option<Transform2> {
    match (from, to) {
        (Surface::Plane(a), Surface::Plane(b)) => {
            let (a, b) = (a.frame(), b.frame());
            let (x, y) = (Vector::from(b.x()), Vector::from(b.y()));
            let (ax, ay) = (Vector::from(a.x()), Vector::from(a.y()));
            let o = a.origin() - b.origin();
            Some(Transform2::orthogonal([[ax.dot(x), ay.dot(x)], [ax.dot(y), ay.dot(y)]], Vector2::new(o.dot(x), o.dot(y))))
        }
        (Surface::Cylinder(a), Surface::Cylinder(b)) => {
            let (a, b) = (a.frame(), b.frame());
            let (x, y, z) = (Vector::from(b.x()), Vector::from(b.y()), Vector::from(b.z()));
            // (Round the same way about the axis, or the other way where it points the other
            // way, from the angle of the one's start round the other.)
            let s = Vector::from(a.z()).dot(z).signum();
            let turn = Vector::from(a.x()).dot(y).atan2(Vector::from(a.x()).dot(x));
            Some(Transform2::orthogonal([[s, 0.0], [0.0, s]], Vector2::new(turn, (a.origin() - b.origin()).dot(z))))
        }
        _ => None,
    }
}

/// Shifts loops by whole periods into the span of the largest one (each aligned from its own
/// first pcurve, a hole may come out a period away from the loop around it).
fn gather(loops: &mut [Vec<Use>], surface: &Surface) {
    let periods = [surface.u_period(), surface.v_period()];
    if loops.len() < 2 || periods.iter().all(Option::is_none) {
        return;
    }
    let span = |uses: &[Use]| -> ([f64; 2], [f64; 2]) {
        let mut lo = [f64::INFINITY; 2];
        let mut hi = [f64::NEG_INFINITY; 2];
        for (a, b) in uses.iter().filter_map(uv_ends) {
            for p in [a, b] {
                lo = [lo[0].min(p.x), lo[1].min(p.y)];
                hi = [hi[0].max(p.x), hi[1].max(p.y)];
            }
        }
        (lo, hi)
    };
    let largest = (0..loops.len()).max_by(|&a, &b| loop_area(&loops[a]).abs().total_cmp(&loop_area(&loops[b]).abs())).expect("loops");
    let (lo, hi) = span(&loops[largest]);
    for (k, uses) in loops.iter_mut().enumerate() {
        if k == largest {
            continue;
        }
        let (a, b) = span(uses);
        let mut shift = [0.0; 2];
        for d in 0..2 {
            let Some(p) = periods[d] else { continue };
            // (The middle of the loop brought nearest the middle of the largest.)
            shift[d] = ((0.5 * (lo[d] + hi[d]) - 0.5 * (a[d] + b[d])) / p).round() * p;
        }
        if shift != [0.0; 2] {
            let t = Transform2::translation(Vector2::new(shift[0], shift[1]));
            for u in uses.iter_mut() {
                u.pcurve = u.pcurve.as_ref().map(|p| p.transformed(&t));
            }
        }
    }
}

/// Shifts the pcurves of a loop by whole periods so that each starts where the one before
/// ends; false if the loop then does not close (it runs round a period).
fn align(uses: &mut [Use], surface: &Surface) -> bool {
    let periods = [surface.u_period(), surface.v_period()];
    let shift = |period: Option<f64>, target: f64, at: f64| period.map_or(0.0, |p| ((target - at) / p).round() * p);
    let mut previous: Option<Point2> = None;
    for u in uses.iter_mut() {
        let Some((start, end)) = uv_ends(u) else { return false };
        if let Some(prev) = previous {
            let (du, dv) = (shift(periods[0], prev.x, start.x), shift(periods[1], prev.y, start.y));
            if du != 0.0 || dv != 0.0 {
                let t = Transform2::translation(Vector2::new(du, dv));
                u.pcurve = u.pcurve.as_ref().map(|p| p.transformed(&t));
            }
            previous = Some(Point2::new(end.x + du, end.y + dv));
        } else {
            previous = Some(end);
        }
    }
    let (Some(first), Some(last)) = (uses.first().and_then(uv_ends), previous) else { return false };
    (first.0.x - last.x).abs() <= 1e-6 && (first.0.y - last.y).abs() <= 1e-6
}

/// Merges chains of edges joined at vertices of no other edge, along one line or circle,
/// into single edges. True if any were.
fn merge_edges(model: &mut Model, specs: &mut [Spec]) -> Result<bool> {
    // The edges at each vertex, and each edge's ends.
    let mut at_vertex: HashMap<VertexId, Vec<EdgeId>> = HashMap::default();
    let mut edges: Vec<EdgeId> = Vec::new();
    let mut seen: HashSet<EdgeId> = HashSet::default();
    for spec in specs.iter() {
        for l in &spec.loops {
            for u in l {
                if seen.insert(u.edge) {
                    edges.push(u.edge);
                    let (a, b) = ends(model, u.edge, Orientation::Forward);
                    at_vertex.entry(a).or_default().push(u.edge);
                    if b != a {
                        at_vertex.entry(b).or_default().push(u.edge);
                    }
                }
            }
        }
    }
    // (A trimmed line or circle as that line or circle, over the same parameters.)
    let untrimmed = |(c, a, b): (Curve, f64, f64)| match c {
        Curve::Trimmed(t) if matches!(t.basis(), Curve::Line(_) | Curve::Circle(_)) => (t.basis().clone(), a, b),
        c => (c, a, b),
    };
    let curves: HashMap<EdgeId, (Curve, f64, f64)> = edges.iter().filter_map(|&e| Some((e, untrimmed(model.edge_curve(Shape::from(e))?)))).collect();
    // The faces using each edge.
    let mut faces_of: HashMap<EdgeId, Vec<usize>> = HashMap::default();
    for (si, spec) in specs.iter().enumerate() {
        for u in spec.loops.iter().flatten() {
            faces_of.entry(u.edge).or_default().push(si);
        }
    }
    for faces in faces_of.values_mut() {
        faces.sort_unstable();
    }
    // The other edge at a vertex of two, if the two lie on one line or circle, or are other
    // curves (pieces of one section) between the same two faces.
    let general = |c: &(Curve, f64, f64)| !matches!(c.0, Curve::Line(_) | Curve::Circle(_));
    let continues = |e: EdgeId, v: VertexId| -> Option<EdgeId> {
        let list = at_vertex.get(&v)?;
        let [x, y] = list[..] else { return None };
        let f = if x == e { y } else { x };
        let (ce, cf) = (curves.get(&e)?, curves.get(&f)?);
        let one_section = general(ce) && general(cf) && faces_of.get(&e)?.len() == 2 && faces_of.get(&e) == faces_of.get(&f);
        (f != e && (same_carrier(ce, cf) || one_section)).then_some(f)
    };
    let mut in_chain: HashSet<EdgeId> = HashSet::default();
    let mut changed = false;
    for &e in &edges {
        if in_chain.contains(&e) || !curves.contains_key(&e) {
            continue;
        }
        // Grown both ways from e: edges in order, with the vertices between and at the ends.
        let (a, b) = ends(model, e, Orientation::Forward);
        if a == b {
            continue;
        }
        let mut chain: Vec<EdgeId> = vec![e];
        let mut members: HashSet<EdgeId> = [e].into_iter().collect();
        let (mut head, mut tail) = (a, b);
        let mut cyclic = false;
        for forwards in [true, false] {
            loop {
                let (last, at) = if forwards { (*chain.last().expect("edge"), tail) } else { (chain[0], head) };
                let Some(next) = continues(last, at) else { break };
                if members.contains(&next) {
                    cyclic = true;
                    break;
                }
                let (p, q) = ends(model, next, Orientation::Forward);
                if p == q {
                    break;
                }
                let far = if p == at { q } else { p };
                members.insert(next);
                if forwards {
                    chain.push(next);
                    tail = far;
                } else {
                    chain.insert(0, next);
                    head = far;
                }
            }
            if cyclic {
                break;
            }
        }
        in_chain.extend(chain.iter().copied());
        if cyclic || chain.len() < 2 || head == tail {
            continue;
        }
        changed |= merge_chain(model, specs, &chain, head, tail, &curves)?;
    }
    Ok(changed)
}

/// True if two edges' curves are one line or one circle.
fn same_carrier(a: &(Curve, f64, f64), b: &(Curve, f64, f64)) -> bool {
    let tol = tolerance::LINEAR;
    match (&a.0, &b.0) {
        (Curve::Line(x), Curve::Line(y)) => {
            let dx = Vector::from(x.direction);
            let off = |p: Point| {
                let d = p - x.origin;
                (d - dx * d.dot(dx)).length()
            };
            dx.cross(Vector::from(y.direction)).length() <= 1e-9 && off(y.origin) <= tol
        }
        (Curve::Circle(x), Curve::Circle(y)) => {
            let axis = |c: &Circle| Vector::from(c.frame().x()).cross(Vector::from(c.frame().y()));
            x.center().distance(y.center()) <= tol && (x.radius() - y.radius()).abs() <= tol && axis(x).cross(axis(y)).length() <= 1e-9
        }
        _ => false,
    }
}

/// [`merge_general_chain`] for pieces of one curve over abutting ranges, with one pcurve in
/// each face (as a face split and merged back leaves its edges): the curve itself over their
/// union, and the pcurve likewise. `None` where the chain is not so.
fn merge_pieces_of_one_curve(
    model: &mut Model,
    specs: &mut [Spec],
    chain: &[EdgeId],
    head: VertexId,
    tail: VertexId,
    curves: &HashMap<EdgeId, (Curve, f64, f64)>,
) -> Result<Option<bool>> {
    let curve = &curves[&chain[0]].0;
    if chain.iter().any(|e| &curves[e].0 != curve) {
        return Ok(None);
    }
    let mut ranges: Vec<(f64, f64)> = chain.iter().map(|e| { let (_, a, b) = curves[e]; (a.min(b), a.max(b)) }).collect();
    ranges.sort_by(|x, y| x.0.total_cmp(&y.0));
    if ranges.windows(2).any(|w| (w[1].0 - w[0].1).abs() > 1e-9 * (1.0 + w[0].1.abs())) {
        return Ok(None);
    }
    let (first, last) = (ranges[0].0, ranges[ranges.len() - 1].1);
    // The edge runs along the curve from the vertex at its start.
    let (pa, pb) = (model.point(Shape::from(head)), model.point(Shape::from(tail)));
    let (start, end) = if curve.point(first).distance(pa) <= curve.point(first).distance(pb) { (head, tail) } else { (tail, head) };
    let members: HashSet<EdgeId> = chain.iter().copied().collect();
    let mut replacements: Vec<(usize, usize, usize, usize, Use)> = Vec::new();
    for (si, spec) in specs.iter().enumerate() {
        for (li, l) in spec.loops.iter().enumerate() {
            let n = l.len();
            if !l.iter().any(|u| members.contains(&u.edge)) {
                continue;
            }
            if n <= chain.len() {
                return Ok(None);
            }
            let Some(begin) = (0..n).find(|&i| members.contains(&l[i].edge) && !members.contains(&l[(i + n - 1) % n].edge)) else {
                return Ok(None);
            };
            let count = (0..n).take_while(|&k| members.contains(&l[(begin + k) % n].edge)).count();
            if count != chain.len() || l.iter().filter(|u| members.contains(&u.edge)).count() != count {
                return Ok(None);
            }
            // One pcurve for all the pieces, used all the same way.
            let uses: Vec<&Use> = (0..count).map(|k| &l[(begin + k) % n]).collect();
            let (Some(pcurve), orientation) = (uses[0].pcurve.clone(), uses[0].orientation) else { return Ok(None) };
            if uses.iter().any(|u| u.pcurve.as_ref() != Some(&pcurve) || u.orientation != orientation) {
                return Ok(None);
            }
            let u = Use { edge: chain[0], orientation, pcurve: Some(pcurve), first, last };
            replacements.push((si, li, begin, count, u));
        }
    }
    if replacements.is_empty() {
        return Ok(None);
    }
    let tolerance = chain.iter().map(|&e| model.tolerance(Shape::from(e))).fold(tolerance::LINEAR, f64::max);
    let edge = model.add_edge(start, end, tolerance)?;
    model.set_edge_curve(edge.into(), curve.clone(), tolerance)?;
    model.set_edge_range(edge.into(), first, last)?;
    for (si, li, begin, count, mut u) in replacements {
        u.edge = edge;
        let spec = &mut specs[si];
        spec.original = None;
        let l = &mut spec.loops[li];
        l.rotate_left(begin);
        l.splice(0..count, std::iter::once(u));
    }
    Ok(Some(true))
}

/// Replaces a chain of edges (in order, from vertex `head` to `tail`) by one edge in every
/// loop using them, if each such face takes a closed-form pcurve for it. True if done.
fn merge_chain(
    model: &mut Model,
    specs: &mut [Spec],
    chain: &[EdgeId],
    head: VertexId,
    tail: VertexId,
    curves: &HashMap<EdgeId, (Curve, f64, f64)>,
) -> Result<bool> {
    let (pa, pb) = (model.point(Shape::from(head)), model.point(Shape::from(tail)));
    let (first_curve, f0, l0) = &curves[&chain[0]];
    let (curve, first, last) = match first_curve {
        Curve::Line(_) => {
            let Some(direction) = Direction::try_from_vector(pb - pa) else { return Ok(false) };
            (Curve::Line(Line::new(pa, direction)), 0.0, pa.distance(pb))
        }
        Curve::Circle(_) => {
            let span: f64 = chain.iter().map(|e| (curves[e].2 - curves[e].1).abs()).sum();
            if span >= std::f64::consts::TAU - 1e-9 {
                return Ok(false);
            }
            // Along the chain from head: the first edge's circle, turned to run that way.
            let (s, _) = ends(model, chain[0], Orientation::Forward);
            let mut circle = first_curve.clone();
            if s != head {
                circle.reverse();
            }
            let _ = (f0, l0);
            let tau = std::f64::consts::TAU;
            let start = crate::geometry::extrema::nearest_on_curve(pa, &circle, 0.0, tau).parameter;
            // The end at the tail vertex itself, near where the spans add up to.
            let guess = start + span;
            let at_tail = crate::geometry::extrema::nearest_on_curve(pb, &circle, 0.0, tau).parameter;
            let end = at_tail + ((guess - at_tail) / tau).round() * tau;
            if !(end > start) {
                return Ok(false);
            }
            (circle, start, end)
        }
        _ => return merge_general_chain(model, specs, chain, head, tail, curves),
    };
    let members: HashSet<EdgeId> = chain.iter().copied().collect();
    // Each loop's run of the chain's edges: where it starts, how long, and the new use.
    let mut replacements: Vec<(usize, usize, usize, usize, Use)> = Vec::new();
    for (si, spec) in specs.iter().enumerate() {
        for (li, l) in spec.loops.iter().enumerate() {
            let n = l.len();
            if !l.iter().any(|u| members.contains(&u.edge)) {
                continue;
            }
            if n <= chain.len() {
                return Ok(false);
            }
            // The run starts after a use not in the chain.
            let Some(begin) = (0..n).find(|&i| members.contains(&l[i].edge) && !members.contains(&l[(i + n - 1) % n].edge)) else {
                return Ok(false);
            };
            let count = (0..n).take_while(|&k| members.contains(&l[(begin + k) % n].edge)).count();
            if count != chain.len() || l.iter().filter(|u| members.contains(&u.edge)).count() != count {
                return Ok(false);
            }
            let first_use = &l[begin];
            let orientation = if ends(model, first_use.edge, first_use.orientation).0 == head { Orientation::Forward } else { Orientation::Reversed };
            let Some(pcurve) = analytic_pcurve(&curve, &spec.surface) else { return Ok(false) };
            let mut u = Use { edge: chain[0], orientation, pcurve: Some(pcurve), first, last };
            let (Some((start, _)), Some((new_start, _))) = (uv_ends(first_use), uv_ends(&u)) else { return Ok(false) };
            let shift = |period: Option<f64>, target: f64, at: f64| period.map_or(0.0, |p| ((target - at) / p).round() * p);
            let (du, dv) = (shift(spec.surface.u_period(), start.x, new_start.x), shift(spec.surface.v_period(), start.y, new_start.y));
            if du != 0.0 || dv != 0.0 {
                u.pcurve = u.pcurve.map(|p| p.transformed(&Transform2::translation(Vector2::new(du, dv))));
            }
            replacements.push((si, li, begin, count, u));
        }
    }
    if replacements.is_empty() {
        return Ok(false);
    }
    let tolerance = chain.iter().map(|&e| model.tolerance(Shape::from(e))).fold(0.0, f64::max);
    let edge = model.add_edge(head, tail, tolerance)?;
    model.set_edge_curve(edge.into(), curve, tolerance)?;
    model.set_edge_range(edge.into(), first, last)?;
    for (si, li, begin, count, mut u) in replacements {
        u.edge = edge;
        let spec = &mut specs[si];
        spec.original = None;
        let l = &mut spec.loops[li];
        let n = l.len();
        // The run's uses replaced by the one (rotating so the run does not wrap round).
        l.rotate_left(begin);
        l.splice(0..count, std::iter::once(u));
        let _ = n;
    }
    Ok(true)
}

/// [`merge_chain`] for a chain of other curves (pieces of one section between two faces):
/// one curve fitted through points along them, and in each face a pcurve fitted through its
/// own pcurves' points there, both by the chord length along the chain.
fn merge_general_chain(
    model: &mut Model,
    specs: &mut [Spec],
    chain: &[EdgeId],
    head: VertexId,
    tail: VertexId,
    curves: &HashMap<EdgeId, (Curve, f64, f64)>,
) -> Result<bool> {
    if let Some(done) = merge_pieces_of_one_curve(model, specs, chain, head, tail, curves)? {
        return Ok(done);
    }
    const SAMPLES: usize = 32;
    // Each edge's parameters along the chain, from head to tail.
    let mut at = head;
    let mut runs: Vec<(EdgeId, Vec<f64>)> = Vec::new();
    for (k, &e) in chain.iter().enumerate() {
        let (s, f) = ends(model, e, Orientation::Forward);
        let (_, a, b) = curves[&e];
        let (from, to, next) = if s == at { (a, b, f) } else if f == at { (b, a, s) } else { return Ok(false) };
        let skip = usize::from(k > 0);
        runs.push((e, (skip..=SAMPLES).map(|i| from + (to - from) * i as f64 / SAMPLES as f64).collect()));
        at = next;
    }
    if at != tail {
        return Ok(false);
    }
    let points: Vec<Point> = runs.iter().flat_map(|(e, ts)| ts.iter().map(|&t| curves[e].0.point(t))).collect();
    let mut params = vec![0.0];
    for w in points.windows(2) {
        params.push(params[params.len() - 1] + w[0].distance(w[1]));
    }
    let length = params[params.len() - 1];
    if !(length > 0.0) {
        return Ok(false);
    }
    let Ok(fitted) = crate::geometry::fit::interpolate_with_parameters(&points, &params) else { return Ok(false) };
    // Within the chain's tolerance of the curves, halfway between the samples.
    let tolerance = chain.iter().map(|&e| model.tolerance(Shape::from(e))).fold(tolerance::LINEAR, f64::max);
    let mut k = 0;
    for (e, ts) in &runs {
        for w in ts.windows(2) {
            let s = 0.5 * (params[k] + params[k + 1]);
            if fitted.point(s).distance(curves[e].0.point(0.5 * (w[0] + w[1]))) > tolerance {
                return Ok(false);
            }
            k += 1;
        }
        k += 1;
    }
    let members: HashSet<EdgeId> = chain.iter().copied().collect();
    let mut replacements: Vec<(usize, usize, usize, usize, Use)> = Vec::new();
    for (si, spec) in specs.iter().enumerate() {
        for (li, l) in spec.loops.iter().enumerate() {
            let n = l.len();
            if !l.iter().any(|u| members.contains(&u.edge)) {
                continue;
            }
            if n <= chain.len() {
                return Ok(false);
            }
            let Some(begin) = (0..n).find(|&i| members.contains(&l[i].edge) && !members.contains(&l[(i + n - 1) % n].edge)) else {
                return Ok(false);
            };
            let count = (0..n).take_while(|&k| members.contains(&l[(begin + k) % n].edge)).count();
            if count != chain.len() || l.iter().filter(|u| members.contains(&u.edge)).count() != count {
                return Ok(false);
            }
            // The face's own pcurves at the chain's points, unwrapped across its periods.
            let mut uv: Vec<Point2> = Vec::with_capacity(points.len());
            for (e, ts) in &runs {
                let Some(pcurve) = l.iter().find(|u| u.edge == *e).and_then(|u| u.pcurve.as_ref()) else { return Ok(false) };
                for &t in ts {
                    let mut q = pcurve.point(t);
                    if let Some(&last) = uv.last() {
                        let unwrap = |x: f64, near: f64, period: Option<f64>| period.map_or(x, |p| x + ((near - x) / p).round() * p);
                        q = Point2::new(unwrap(q.x, last.x, spec.surface.u_period()), unwrap(q.y, last.y, spec.surface.v_period()));
                    }
                    uv.push(q);
                }
            }
            let Ok(pcurve) = crate::geometry::fit::interpolate_with_parameters(&uv, &params) else { return Ok(false) };
            let first_use = &l[begin];
            let orientation = if ends(model, first_use.edge, first_use.orientation).0 == head { Orientation::Forward } else { Orientation::Reversed };
            let u = Use { edge: chain[0], orientation, pcurve: Some(Curve2::BSpline(pcurve)), first: 0.0, last: length };
            replacements.push((si, li, begin, count, u));
        }
    }
    if replacements.is_empty() {
        return Ok(false);
    }
    let edge = model.add_edge(head, tail, tolerance)?;
    model.set_edge_curve(edge.into(), Curve::BSpline(fitted), tolerance)?;
    model.set_edge_range(edge.into(), 0.0, length)?;
    for (si, li, begin, count, mut u) in replacements {
        u.edge = edge;
        let spec = &mut specs[si];
        spec.original = None;
        let l = &mut spec.loops[li];
        l.rotate_left(begin);
        l.splice(0..count, std::iter::once(u));
    }
    Ok(true)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::boolean;
    use crate::math::Frame;
    use crate::primitives::{Cuboid, Cylinder};
    use crate::{properties, validity};

    /// The offset all round a segment (each side cut in two pieces, an arc round each end)
    /// unified: one line each side, the arcs kept (as OCCT's clean leaves a stadium).
    #[test]
    fn a_wire_alone_has_its_pieces_of_one_line_merged() {
        let mut model = Model::new();
        let segment = crate::make::edge_between(&mut model, Point::new(0.0, 0.0, 0.0), Point::new(10.0, 0.0, 0.0)).unwrap();
        let segment = crate::make::wire(&mut model, &[segment]).unwrap();
        let around = crate::offset::planar_wire(&mut model, segment, 1.0, crate::offset::Join::Arc).unwrap();
        assert_eq!(model.ordered_edges(around).len(), 6);
        let merged = unify_same_domain(&mut model, around).unwrap();
        let edges = model.ordered_edges(merged);
        assert_eq!(edges.len(), 4);
        assert_eq!(edges.iter().filter(|&&e| matches!(model.edge_curve(e), Some((Curve::Line(_), _, _)))).count(), 2);
        assert!(model.is_closed_wire(WireId::try_from(merged.node).unwrap()));
        let length = properties::length(&model, merged).mass;
        assert!((length - (20.0 + 2.0 * std::f64::consts::PI)).abs() < 1e-9, "{length}");
    }

    fn cuboid(model: &mut Model, corner: (f64, f64, f64), size: (f64, f64, f64)) -> Shape {
        let frame = Frame::new(Point::new(corner.0, corner.1, corner.2), Direction::Z, Direction::X);
        Cuboid::new(size.0, size.1, size.2).with_frame(frame).build(model).unwrap()
    }

    fn counts(model: &Model, shape: Shape) -> (usize, usize) {
        (model.unique_subshapes(shape, NodeKind::Face).len(), model.unique_subshapes(shape, NodeKind::Edge).len())
    }

    #[test]
    fn boxes_side_by_side_become_one() {
        let mut model = Model::new();
        let a = cuboid(&mut model, (0.0, 0.0, 0.0), (1.0, 1.0, 1.0));
        let b = cuboid(&mut model, (1.0, 0.0, 0.0), (1.0, 1.0, 1.0));
        let fused = boolean::union(&mut model, a, b).unwrap();
        let clean = unify_same_domain(&mut model, fused).unwrap();
        assert!(validity::check(&model, clean).is_valid(), "{:?}", validity::check(&model, clean));
        assert_eq!(counts(&model, clean), (6, 12));
        assert!((properties::volume(&model, clean).mass - 2.0).abs() < 1e-12);
    }

    #[test]
    fn a_face_alone_has_its_edges_in_line_merged() {
        // A square with a vertex partway along a side, alone and in a compound of it alone.
        let mut model = Model::new();
        let corners = [(0.0, 5.0), (0.0, 10.0), (10.0, 10.0), (10.0, 0.0), (0.0, 0.0)];
        let vertices: Vec<Shape> = corners.iter().map(|&(x, y)| crate::make::vertex(&mut model, Point::new(x, y, 0.0))).collect();
        let edges: Vec<Shape> = (0..5).map(|i| crate::make::edge(&mut model, vertices[i], vertices[(i + 1) % 5]).unwrap()).collect();
        let wire = crate::make::wire(&mut model, &edges).unwrap();
        let face = crate::make::planar_face(&mut model, wire).unwrap();
        let clean = unify_same_domain(&mut model, face).unwrap();
        assert_eq!(clean.kind(), NodeKind::Face);
        assert_eq!(counts(&model, clean), (1, 4));
        assert!((properties::area(&model, clean).mass - 100.0).abs() < 1e-12);
        let compound = model.new_empty(NodeKind::Compound);
        model.add(compound, face).unwrap();
        let clean = unify_same_domain(&mut model, compound).unwrap();
        assert_eq!(counts(&model, clean), (1, 4));
    }

    #[test]
    fn holes_kept_within_the_period_of_their_outer_loop() {
        // A cylinder with its bottom edge rounded and a hole in its side, and a cylinder on
        // top of it: their sides join, the hole's loop in the period of the side's outer one
        // (or a cut through the top then goes wrong).
        let mut model = Model::new();
        // (Each an extruded disc, as CadQuery's are.)
        let disc = |model: &mut Model, centre: (f64, f64, f64), normal: Direction, x: Direction, radius: f64, length: f64| {
            let frame = Frame::new(Point::new(centre.0, centre.1, centre.2), normal, x);
            let circle = crate::make::edge_from_curve(model, crate::geometry::curve::Curve::Circle(crate::geometry::curve::Circle::new(frame, radius))).unwrap();
            let wire = crate::make::wire(model, &[circle]).unwrap();
            let face = crate::make::planar_face(model, wire).unwrap();
            crate::sweep::extrude(model, face, Vector::from(normal) * length).unwrap().shape
        };
        let mut base = disc(&mut model, (0.0, 0.0, 0.0), Direction::Z, Direction::X, 12.5, 16.0);
        let bottom: Vec<Shape> = model.explore(base, NodeKind::Edge).filter(|&e| crate::properties::bounding_box(&model, e).max.z < 1.0).collect();
        base = crate::blend::fillet(&mut model, base, &bottom, 1.25).unwrap();
        let hole = disc(&mut model, (-20.0, -5.0, 8.0), Direction::X, Direction::Y, 2.0, 10.0);
        base = boolean::difference(&mut model, base, hole).unwrap();
        let top = disc(&mut model, (0.0, 0.0, 16.0), Direction::Z, Direction::X, 12.5, 4.75);
        let fused = boolean::union(&mut model, base, top).unwrap();
        let clean = unify_same_domain(&mut model, fused).unwrap();
        assert!(validity::check(&model, clean).is_valid(), "{:?}", validity::check(&model, clean));
        let volume = properties::volume(&model, clean).mass;
        let pocket = cuboid(&mut model, (-3.0, -3.0, 16.0), (6.0, 6.0, 4.75));
        let cut = boolean::difference(&mut model, clean, pocket).unwrap();
        assert!(validity::check(&model, cut).is_valid(), "{:?}", validity::check(&model, cut));
        let removed = volume - properties::volume(&model, cut).mass;
        assert!((removed - 36.0 * 4.75).abs() < 1e-6, "{removed}");
    }

    #[test]
    fn a_slot_of_a_box_and_two_cylinders() {
        // A box with half cylinders at its ends: their halves inside the box go, the planar
        // sides join across the cylinders' tangent lines, as in OCCT (6 faces).
        let mut model = Model::new();
        let mut slot = cuboid(&mut model, (0.0, -1.0, 0.0), (4.0, 2.0, 1.0));
        for x in [0.0, 4.0] {
            let frame = Frame::new(Point::new(x, 0.0, 0.0), Direction::Z, Direction::X);
            let c = Cylinder::new(1.0, 1.0).with_frame(frame).build(&mut model).unwrap();
            slot = boolean::union(&mut model, slot, c).unwrap();
        }
        let volume = properties::volume(&model, slot).mass;
        let clean = unify_same_domain(&mut model, slot).unwrap();
        assert!(validity::check(&model, clean).is_valid(), "{:?}", validity::check(&model, clean));
        assert!((properties::volume(&model, clean).mass - volume).abs() < 1e-9);
        assert_eq!(counts(&model, clean).0, 6);
    }
}
