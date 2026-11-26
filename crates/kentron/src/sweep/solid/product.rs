//! The swept solid under a steady motion, built as the product of the boundary of `W` and the
//! motion, capped by the rest of the body at the start and `W` at the end (see the module
//! above).

use std::collections::HashMap;

use crate::error::{Error, Result};
use crate::geometry::curve::BSplineCurve;
use crate::geometry::surface::{BSplineSurface, KnotVector};
use crate::geometry::curve::{Circle, Curve, Line, TrimmedCurve};
use crate::geometry::curve2::{Curve2, Line2};
use crate::geometry::surface::{Cylinder, Plane, RevolutionSurface, Surface};
use crate::geometry::{ParametricCurve, ParametricSurface};
use crate::math::{Direction, Direction2, Frame, Point, Point2, Transform, Transformable, Vector};
use crate::sweep::pipe::{path, PieceMotion};
use crate::tolerance;
use crate::topology::{EdgeId, FaceId, Model, NodeKind, Orientation, Shape, VertexId};

use super::{contour, Steady};

/// Which way the motion carries a face (or a piece of one) relative to its outward normal.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Side {
    /// Outwards somewhere (`V·n > 0`): in `W`.
    Front,
    /// Inwards, or along the face everywhere.
    Back,
}

/// The swept solid of `solid` (unplaced; changed: its faces are divided where the motion
/// grazes them) under `steady`, a screw built in stages of at most a quarter turn (so that no
/// face winds round its axis: each is a strip short enough that a point's place on it is found
/// from any start on it).
pub(super) fn sweep(model: &mut Model, solid: Shape, steady: &Steady) -> Result<Shape> {
    let (stage, n) = match steady.clone() {
        Steady::Screw { axis, angle, rise } => {
            let n = (angle / std::f64::consts::FRAC_PI_2).ceil().max(1.0);
            (Steady::Screw { axis, angle: angle / n, rise: rise / n }, n as usize)
        }
        other => (other, 1),
    };
    let junctions = (0..=n).map(|k| steady.at(k as f64 / n as f64)).collect();
    sweep_in_stages(model, solid, &vec![stage; n], junctions)
}

/// Whether `W` passes back into the body: a point inside one of a stage's faces in `W`,
/// carried along by the motion over the stage, found inside the body where it stood at the
/// start of the motion or of a stage before. Then the product runs into itself: its winding is
/// 2 where the body and `W` passing by, or `W` twice, cover a point (as a torus slid square to
/// its axis passes its near side and then its far side over its hole, or a flat face turned
/// about an axis passes its chords' points), and the product is resolved. (Under stages of
/// one steady motion the body's place at the start alone: one stage after another is the
/// same motion.) Found at points on a grid of each face's parameters (inside it, 16 a side),
/// each at places along the whole motion as close as the grid's points: at least 32, and
/// more where the fastest point of the body moves farther than a sixteenth of its size
/// between them (an overlap may be as thin as the window over which a point is back inside
/// the body: a pin's end turned about a far axis covers the bulge of its rim twice, each
/// point of it for a few degrees of a long turn). An overlap thinner than those steps may go
/// unnoticed. (Where a point is found outside the body, the ball round it clear of the
/// body's surfaces is kept, and points carried into one are outside without asking again:
/// see [`Clearance`].)
fn passes_back_into(model: &Model, solid: Shape, faces: &[Shape], sides: &[Vec<Side>], stages: &[Steady], junctions: &[Transform]) -> Result<bool> {
    const N: usize = 16;
    const STEPS: usize = 32;
    let n = stages.len();
    let one_motion = !matches!(stages[0], Steady::Carried(_)) && stages.iter().all(|s| same(s, &stages[0]));
    // (As far apart along the motion as the grid's points are on the body.)
    let bounds = crate::properties::bounding_box(model, solid);
    let size = (bounds.max - bounds.min).length();
    let travel: f64 = stages.iter().map(|s| speed_bound(model, solid, s)).sum();
    let steps = if size > 0.0 && travel.is_finite() { STEPS.max((travel * N as f64 / size).ceil() as usize) } else { STEPS };
    let per_stage = steps.div_ceil(n).max(1);
    let classifier = crate::boolean::PointClassifier::new(model, solid)?;
    let mut clear = Clearance::new(model, solid);
    let starts: Vec<Transform> = junctions.iter().take(n).map(|j| j.inverse()).collect();
    for (k, steady) in stages.iter().enumerate() {
        if one_motion && k > 0 {
            break;
        }
        // The places over the stage (over the whole motion, for one steady motion), from the
        // body's place at the start of the motion or of each stage up to this one.
        let places: Vec<Transform> = if one_motion {
            (1..=per_stage * n).map(|j| {
                let whole = (j - 1) / per_stage;
                steady.at((j - whole * per_stage) as f64 / per_stage as f64) * junctions[whole]
            }).collect()
        } else {
            (1..=per_stage).map(|j| steady.at(j as f64 / per_stage as f64) * junctions[k]).collect()
        };
        let references: &[Transform] = if one_motion { &starts[..1] } else { &starts[..=k] };
        for (i, &face) in faces.iter().enumerate() {
            if sides[k][i] != Side::Front {
                continue;
            }
            let Some(surface) = model.surface(face) else { continue };
            let Some(polygons) = contour::loop_polygons(model, face) else { continue };
            let (u0, u1, v0, v1) = model.face_uv_bounds(face);
            for a in 0..N {
                for b in 0..N {
                    let q = Point2::new(u0 + (u1 - u0) * (a as f64 + 0.5) / N as f64, v0 + (v1 - v0) * (b as f64 + 0.5) / N as f64);
                    if !contour::within(&polygons, q) {
                        continue;
                    }
                    let p = surface.point(q.x, q.y).transformed(&face.placement);
                    for place in &places {
                        let carried = p.transformed(place);
                        if references.iter().any(|r| clear.inside(&classifier, carried.transformed(r))) {
                            return Ok(true);
                        }
                    }
                }
            }
        }
    }
    Ok(false)
}

/// Balls of space known to be outside a body, round points found outside it: as far from
/// each as the body's faces' surfaces are (a bound below that distance), less twice the
/// distance within which a point counts as on a face. No point of the body's boundary lies
/// within one, so the whole ball is on the side its centre is on, and a point in it is told
/// outside, as [`crate::boolean::PointClassifier::locate`] tells it, without casting rays
/// through the faces. Kept in the cells of a grid an eighth of the body's size a side, each
/// ball in those its box meets, its radius at most a cell's side.
struct Clearance {
    /// The faces' surfaces, placed (or the boxes holding them); `None` where a face's
    /// surface is one whose distance from a point is not bounded here: then no ball is kept.
    surfaces: Option<Vec<Bounded>>,
    /// Points within this of a face's surface may be on it.
    margin: f64,
    cell: f64,
    balls: Vec<(Point, f64)>,
    grid: HashMap<(i64, i64, i64), Vec<u32>>,
}

/// A surface whose distance from a point is bounded below: exactly, or by a box holding it.
enum Bounded {
    Exact(Surface),
    Within(crate::math::BoundingBox),
}

impl Clearance {
    fn new(model: &Model, solid: Shape) -> Clearance {
        let surfaces: Option<Vec<Bounded>> = model
            .unique_subshapes(solid, NodeKind::Face)
            .into_iter()
            .map(|face| bounded(&model.surface(face)?.transformed(&face.placement)))
            .collect();
        let b = crate::properties::bounding_box(model, solid);
        let size = (b.max - b.min).length();
        Clearance { surfaces, margin: 2.0 * crate::boolean::on_tol(), cell: size / 8.0, balls: Vec::new(), grid: HashMap::new() }
    }

    fn key(&self, p: Point) -> (i64, i64, i64) {
        ((p.x / self.cell).floor() as i64, (p.y / self.cell).floor() as i64, (p.z / self.cell).floor() as i64)
    }

    /// Whether `p` is inside the body: `classifier` asked where it is in no ball known
    /// outside, and a ball kept round it where it is found outside.
    fn inside(&mut self, classifier: &crate::boolean::PointClassifier, p: Point) -> bool {
        let located = |p: Point| classifier.locate(p);
        let Some(surfaces) = &self.surfaces else {
            return located(p) == crate::boolean::PointLocation::Inside;
        };
        // (Cells keyed by integers: none for a point too far out, or a body of no size.)
        if !(self.cell > 0.0) || ![p.x, p.y, p.z].iter().all(|c| (c / self.cell).abs() < 1e15) {
            return located(p) == crate::boolean::PointLocation::Inside;
        }
        if self.grid.get(&self.key(p)).is_some_and(|near| near.iter().any(|&k| self.balls[k as usize].0.distance(p) < self.balls[k as usize].1)) {
            return false;
        }
        let at = located(p);
        if at == crate::boolean::PointLocation::Outside {
            let clear = surfaces.iter().map(|s| s.distance_bound(p)).fold(f64::INFINITY, f64::min);
            let r = (clear - self.margin).min(self.cell);
            if r > 0.0 {
                let k = self.balls.len() as u32;
                self.balls.push((p, r));
                let (lo, hi) = (self.key(p - Vector::new(r, r, r)), self.key(p + Vector::new(r, r, r)));
                for i in lo.0..=hi.0 {
                    for j in lo.1..=hi.1 {
                        for l in lo.2..=hi.2 {
                            self.grid.entry((i, j, l)).or_default().push(k);
                        }
                    }
                }
            }
        }
        at == crate::boolean::PointLocation::Inside
    }
}

/// `surface` as [`Bounded`]: a plane, cylinder, cone, sphere or torus exactly, a free-form
/// surface by the box of its poles; `None` for others.
fn bounded(surface: &Surface) -> Option<Bounded> {
    match surface {
        Surface::Plane(_) | Surface::Cylinder(_) | Surface::Cone(_) | Surface::Sphere(_) | Surface::Torus(_) => Some(Bounded::Exact(surface.clone())),
        Surface::BSpline(_) | Surface::Bezier(_) => crate::properties::control_box(surface).map(Bounded::Within),
        Surface::Trimmed(t) => bounded(t.basis()),
        _ => None,
    }
}

impl Bounded {
    /// A bound below the distance from `p` to the surface (the whole of it, untrimmed):
    /// exact for a plane, cylinder, sphere and torus; for a cone, the distance to the two
    /// lines it cuts in the plane through `p` and its axis; for a box, to the box.
    fn distance_bound(&self, p: Point) -> f64 {
        // `p` in a frame's coordinates: how far from its z axis, and how far along it.
        let local = |f: &Frame| {
            let d = p - f.origin();
            let z = d.dot(Vector::from(f.z()));
            ((d - Vector::from(f.z()) * z).length(), z)
        };
        match self {
            Bounded::Within(b) => {
                let gap = |x: f64, lo: f64, hi: f64| (lo - x).max(x - hi).max(0.0);
                Vector::new(gap(p.x, b.min.x, b.max.x), gap(p.y, b.min.y, b.max.y), gap(p.z, b.min.z, b.max.z)).length()
            }
            Bounded::Exact(Surface::Plane(s)) => (p - s.origin()).dot(Vector::from(s.normal())).abs(),
            Bounded::Exact(Surface::Cylinder(s)) => (local(s.frame()).0 - s.radius()).abs(),
            Bounded::Exact(Surface::Sphere(s)) => (p.distance(s.frame().origin()) - s.radius()).abs(),
            Bounded::Exact(Surface::Torus(s)) => {
                // (The tube's circles are centred on its centre circle, as near `p` as its
                // nearest point and as far as its farthest: the tube passes through `p` where
                // its radius lies between those distances.)
                let (rho, z) = local(s.frame());
                let (near, far) = ((rho - s.major_radius()).hypot(z), (rho + s.major_radius()).hypot(z));
                let r = s.minor_radius();
                if r < near { near - r } else if r > far { r - far } else { 0.0 }
            }
            Bounded::Exact(Surface::Cone(s)) => {
                let (rho, z) = local(s.frame());
                let (sin, cos) = s.semi_angle_sin_cos();
                let off_line = |x: f64| ((x - s.radius()) * cos - z * sin).abs();
                off_line(rho).min(off_line(-rho))
            }
            Bounded::Exact(_) => 0.0,
        }
    }
}

/// A point inside `cell`: just inside a face of it, off a point well inside the face (in its
/// parameters), as near it as needed for a thin cell; `None` if none is found.
fn point_within(model: &Model, cell: Shape) -> Result<Option<Point>> {
    let classifier = crate::boolean::PointClassifier::new(model, cell)?;
    let b = crate::properties::bounding_box(model, cell);
    let size = (b.max - b.min).length();
    for face in model.unique_subshapes(cell, NodeKind::Face) {
        let Some(surface) = model.surface(face) else { continue };
        let Some(polygons) = contour::loop_polygons(model, Shape { orientation: Orientation::Forward, ..face }) else { continue };
        let closed: Vec<Vec<[f64; 2]>> = polygons.iter().filter(|p| !p.is_empty()).map(|p| p.iter().chain(&p[..1]).map(|q| [q.x, q.y]).collect()).collect();
        if closed.is_empty() {
            continue;
        }
        let [u, v] = crate::boolean::interior_uv(&closed);
        let (p, du, dv) = surface.d1(u, v);
        let n = du.cross(dv);
        if !(n.length() > 0.0) {
            continue;
        }
        let n = n.normalized() * if face.orientation == Orientation::Reversed { -1.0 } else { 1.0 };
        let (p, n) = (p.transformed(&face.placement), n.transformed(&face.placement));
        for step in [1e-3, 1e-4, 1e-5] {
            let q = p - n * (step * size);
            if classifier.locate(q) == crate::boolean::PointLocation::Inside {
                return Ok(Some(q));
            }
        }
    }
    Ok(None)
}

/// The volume a solid bounds, from its faces' volume terms, which the model keeps (the
/// pieces' faces' found as the pieces were enclosed).
fn volume_of(model: &Model, solid: Shape) -> f64 {
    crate::properties::solid_volume(model, solid).map_or_else(|| crate::properties::volume(model, solid).mass, |(v, _)| v)
}

/// The product where it runs into itself (see the module's notes): its faces split where they
/// cross, the pieces of space it winds round at least once, joined. An error where the
/// pieces' volumes, each counted as often as the product winds round it, do not add up to the
/// product's (pieces missed or misjudged).
fn resolved(model: &mut Model, product: Shape) -> Result<Shape> {
    let faces = model.unique_subshapes(product, NodeKind::Face);
    // (Split all at once: the volumes adding up below check the pieces.)
    let cells = crate::boolean::enclosed_solids_at_once(model, &faces)?;
    let cells = model.unique_subshapes(cells, NodeKind::Solid);
    if cells.is_empty() {
        return Err(Error::Failed("a sweep running into itself, not split into pieces".into()));
    }
    let classifier = crate::boolean::PointClassifier::new(model, product)?;
    let mut kept = Vec::new();
    let mut counted = 0.0;
    for &cell in &cells {
        // (A sliver with no point found inside it left out: the volumes adding up below tell
        // whether it held anything.)
        let found = match crate::boolean::point_inside(model, cell)? {
            Some(p) => Some(p),
            None => point_within(model, cell)?,
        };
        let Some(p) = found else { continue };
        let winding = classifier.winding(p).ok_or_else(|| Error::Failed("how often a sweep winds round a point of it".into()))?;
        counted += winding as f64 * volume_of(model, cell);
        if winding >= 1 {
            kept.push(cell);
        }
    }
    let whole = volume_of(model, product);
    if !((counted - whole).abs() <= 1e-6 * whole.abs()) {
        return Err(Error::Failed(format!("a sweep running into itself: its pieces, as often as it winds round them, hold {counted}, not its {whole}")));
    }
    super::super::joined(model, &kept)
}

/// The swept solid of `solid` (unplaced; changed) under motions one after another, each steady
/// as `W` goes (its field fixed on the body, or its curves where it grazes the body fixed),
/// stage `k` carrying the body from where `junctions[k]` puts it to where `junctions[k + 1]`
/// does. Where a face of the body passes from `W` to outside it, or back, at a junction, it
/// closes the stages' sweeps there, facing the way the stage before carried it.
pub(super) fn sweep_in_stages(model: &mut Model, solid: Shape, stages: &[Steady], junctions: Vec<Transform>) -> Result<Shape> {
    let bounds = crate::properties::bounding_box(model, solid);
    let size = (bounds.max - bounds.min).length();
    // The characteristic curves within a billionth of the body's size.
    let tolerance = 1e-9 * size;
    // Faces closed round a seam are split, so that each curve is followed within one period.
    for _ in 0..2 {
        crate::heal::split_closed_faces(model, solid)?;
    }
    // (Each distinct stage's curves once.)
    let mut distinct: Vec<&Steady> = Vec::new();
    for s in stages {
        if !distinct.iter().any(|d| same(d, s)) {
            distinct.push(s);
        }
    }
    for steady in &distinct {
        for face in model.unique_subshapes(solid, NodeKind::Face) {
            let cuts = contour::cuts(model, face, steady, tolerance)?;
            if !cuts.is_empty() {
                crate::divide::divide_face(model, face, &cuts)?;
            }
        }
    }
    // (Under a turn, the points on its axis stand still: edges are split where they cross
    // it, so that each piece's path ends there, at a point.)
    if let Some(axis) = still_axis(stages) {
        split_at_axis(model, solid, axis, tolerance::LINEAR.max(tolerance))?;
    }
    let faces = model.unique_subshapes(solid, NodeKind::Face);
    let mut sides: Vec<Vec<Side>> = Vec::with_capacity(stages.len());
    let mut folds = 0;
    for steady in stages {
        let speed = speed_bound(model, solid, steady);
        let s: Vec<Side> = faces.iter().map(|&f| side(model, f, steady, speed)).collect::<Result<_>>()?;
        folds += split_folds(model, &faces, &s, steady)?;
        sides.push(s);
    }
    // (Where the surface an edge of `W`'s boundary sweeps folds, the motion turning past the
    // direction along it, the product runs into itself: its sheets on either side of the fold
    // meet along the fold's path, tangent there, and cross beyond it. Splitting the product
    // where it crosses itself does not part sheets meeting so, and finds no piece of space
    // after all the work of trying; so a fold is refused here, before the product is made.)
    if folds > 0 {
        return Err(Error::Failed(format!("a sweep folding over itself ({folds} fold{} where the motion turns past a sharp edge or characteristic curve of the tool): its sheets meet tangent along the fold, and are not parted", if folds == 1 { "" } else { "s" })));
    }
    // (How far the body reaches from a turn's axis: the farthest corner of its box.)
    let reach = (0..8)
        .map(|i| {
            let p = Point::new(if i & 1 == 0 { bounds.min.x } else { bounds.max.x }, if i & 2 == 0 { bounds.min.y } else { bounds.max.y }, if i & 4 == 0 { bounds.min.z } else { bounds.max.z });
            match stages[0] {
                Steady::Turn { axis, .. } | Steady::Screw { axis, .. } => {
                    let r = p - axis.origin;
                    (r - Vector::from(axis.direction) * r.dot(Vector::from(axis.direction))).length()
                }
                Steady::Slide { .. } | Steady::Carried(_) => 0.0,
            }
        })
        .fold(0.0, f64::max);
    // (A whole turn ends where it starts: its product holds the whole body, which `W` passes
    // over, and always runs into itself.)
    let whole = is_whole_turn(stages);
    let runs_into_itself = whole || passes_back_into(model, solid, &faces, &sides, stages, &junctions)?;
    let mut builder = Builder::new(model, stages, junctions, size, reach)?;
    let product = builder.product(&faces, &sides)?;
    // (Where the split fails, the product as it is: a solid running into itself, its volume
    // counting what is covered twice twice, which a boolean operation with it may yet take
    // right. Not for a turn: its product, as it is, is a wrong solid (a ring turned nearly
    // whole holds 1208 where it sweeps 761), and a whole turn's holds the body's faces within
    // it.)
    if runs_into_itself {
        match resolved(model, product) {
            Ok(resolved) => return Ok(resolved),
            Err(e) if still_axis(stages).is_some() => return Err(e),
            Err(_) => {}
        }
    }
    Ok(product)
}

/// The axis of a turn (one stage), whose points stand still.
fn still_axis(stages: &[Steady]) -> Option<crate::math::Axis> {
    match stages {
        [Steady::Turn { axis, .. }] => Some(*axis),
        _ => None,
    }
}

/// Whether the stages are one whole turn.
fn is_whole_turn(stages: &[Steady]) -> bool {
    matches!(stages, [Steady::Turn { angle, .. }] if *angle >= std::f64::consts::TAU)
}

/// How far `p` lies from the line `axis`.
fn off_axis(axis: crate::math::Axis, p: Point) -> f64 {
    let r = p - axis.origin;
    let d = Vector::from(axis.direction);
    (r - d * r.dot(d)).length()
}

/// Splits the solid's edges where they cross `axis` (pass within `tolerance` of it), but those
/// lying along it. Found as the nearest points to the axis between samples along each edge,
/// each narrowed by golden sections.
fn split_at_axis(model: &mut Model, solid: Shape, axis: crate::math::Axis, tolerance: f64) -> Result<()> {
    for edge in model.unique_subshapes(solid, NodeKind::Edge) {
        let e = EdgeId::from_index(edge.node.index);
        if model.is_degenerate(edge) {
            continue;
        }
        let Some((curve, a, b)) = model.edge_curve(edge) else { continue };
        let samples = match &curve {
            Curve::BSpline(bs) => 64.max(4 * bs.poles().len()),
            _ => 64,
        };
        let at = |i: usize| a + (b - a) * i as f64 / samples as f64;
        let d: Vec<f64> = (0..=samples).map(|i| off_axis(axis, curve.point(at(i)))).collect();
        if d.iter().all(|&x| x <= tolerance) {
            continue;
        }
        let mut roots = Vec::new();
        for i in 1..samples {
            if !(d[i] <= d[i - 1] && d[i] <= d[i + 1]) {
                continue;
            }
            // (Golden sections on the distance between the neighbours.)
            let (mut lo, mut hi) = (at(i - 1), at(i + 1));
            let g = 0.5 * (5f64.sqrt() - 1.0);
            let f = |s: f64| off_axis(axis, curve.point(s));
            let (mut x1, mut x2) = (hi - g * (hi - lo), lo + g * (hi - lo));
            let (mut f1, mut f2) = (f(x1), f(x2));
            for _ in 0..100 {
                if f1 <= f2 {
                    hi = x2;
                    (x2, f2) = (x1, f1);
                    x1 = hi - g * (hi - lo);
                    f1 = f(x1);
                } else {
                    lo = x1;
                    (x1, f1) = (x2, f2);
                    x2 = lo + g * (hi - lo);
                    f2 = f(x2);
                }
                if hi - lo <= 1e-15 * (b - a) {
                    break;
                }
            }
            let s = 0.5 * (lo + hi);
            let p = curve.point(s);
            let (start, end) = model.edge_end_vertices(e);
            let near_end = [start, end].iter().any(|&v| model[v].point.distance(p) <= tolerance.max(model[v].tolerance));
            if f(s) <= tolerance && !near_end && s > a + 1e-9 * (b - a) && s < b - 1e-9 * (b - a) && roots.last().is_none_or(|&r: &f64| s - r > 1e-9 * (b - a)) {
                roots.push(s);
            }
        }
        let mut piece = e;
        for s in roots {
            let tol = model[piece].tolerance;
            // (The vertex on the axis, where the point stands still.)
            let p = curve.point(s);
            let d = Vector::from(axis.direction);
            let foot = axis.origin + d * (p - axis.origin).dot(d);
            let v = model.add_vertex(foot, tol.max(p.distance(foot)));
            piece = model.split_edge(piece, v, s)?.1;
        }
    }
    Ok(())
}

/// Whether two stages run along the body alike (the same steady motion).
fn same(a: &Steady, b: &Steady) -> bool {
    match (a, b) {
        (Steady::Carried(x), Steady::Carried(y)) => std::sync::Arc::ptr_eq(x, y),
        (Steady::Slide { by: x }, Steady::Slide { by: y }) => x == y,
        (Steady::Turn { axis: x, angle: p }, Steady::Turn { axis: y, angle: q }) => x == y && p == q,
        (Steady::Screw { axis: x, angle: p, rise: r }, Steady::Screw { axis: y, angle: q, rise: w }) => x == y && p == q && r == w,
        _ => false,
    }
}

/// The edges between a face in `W` and one outside it (each with whether the face in `W` runs
/// along the edge's curve), but degenerate ones: the boundary of `W`.
fn boundary_edges(model: &Model, faces: &[Shape], sides: &[Side]) -> Result<Vec<(EdgeId, bool)>> {
    // Each edge's uses: the face (by index) and whether it runs along the edge's curve.
    let mut uses: HashMap<EdgeId, Vec<(usize, bool)>> = HashMap::new();
    for (k, &face) in faces.iter().enumerate() {
        let loops = crate::heal::face_loops(model, face).ok_or_else(|| Error::Failed("a face without pcurves".into()))?;
        let reversed_face = face.orientation == Orientation::Reversed;
        for u in loops.iter().flatten() {
            uses.entry(u.edge).or_default().push((k, u.reversed == reversed_face));
        }
    }
    let mut edges: Vec<EdgeId> = uses.keys().copied().collect();
    edges.sort_by_key(|e| e.index());
    let mut out = Vec::new();
    for e in edges {
        let list = &uses[&e];
        let fronts: Vec<&(usize, bool)> = list.iter().filter(|(k, _)| sides[*k] == Side::Front).collect();
        if list.len() != 2 || fronts.len() != 1 || list[0].0 == list[1].0 || model.is_degenerate(Shape::from(e)) {
            continue;
        }
        out.push((e, fronts[0].1));
    }
    Ok(out)
}

/// Splits the edges bounding `W` where the surface they sweep folds: where the motion runs
/// along the edge (its tangent crossed with the velocity, the swept surface's normal,
/// vanishes and turns over), as a characteristic curve does where it turns past the direction
/// of a slide. Each piece then sweeps a surface without a fold, meeting the next along the
/// fold's path. Returns how many folds were found.
fn split_folds(model: &mut Model, faces: &[Shape], sides: &[Side], steady: &Steady) -> Result<usize> {
    let mut count = 0;
    for (e, _) in boundary_edges(model, faces, sides)? {
        let Some((curve, a, b)) = model.edge_curve(Shape::from(e)) else { continue };
        let normal = |s: f64| {
            let (p, d) = curve.d1(s);
            d.cross(steady.velocity(p))
        };
        let samples = match &curve {
            Curve::BSpline(bs) => 64.max(4 * bs.poles().len()),
            _ => 64,
        };
        let inside = |at: f64| at - a > 1e-9 * (b - a) && b - at > 1e-9 * (b - a);
        let mut roots = Vec::new();
        let mut previous = (a, normal(a));
        for i in 1..=samples {
            let s = a + (b - a) * i as f64 / samples as f64;
            let n = normal(s);
            if n.dot(previous.1) < 0.0 {
                let reference = previous.1;
                let f = |x: f64| normal(previous.0 + (s - previous.0) * x).dot(reference);
                let x = super::contour::root(&f, f(0.0), f(1.0));
                let at = previous.0 + (s - previous.0) * x;
                if inside(at) {
                    roots.push(at);
                    count += 1;
                }
            }
            previous = (s, n);
        }
        roots.sort_by(f64::total_cmp);
        let mut piece = e;
        for at in roots {
            let tolerance = model[piece].tolerance;
            let v = model.add_vertex(curve.point(at), tolerance);
            piece = model.split_edge(piece, v, at)?.1;
        }
    }
    Ok(count)
}

/// The largest speed over the body's box.
fn speed_bound(model: &Model, solid: Shape, steady: &Steady) -> f64 {
    let b = crate::properties::bounding_box(model, solid);
    let mut speed = 0.0f64;
    for i in 0..8 {
        let p = Point::new(if i & 1 == 0 { b.min.x } else { b.max.x }, if i & 2 == 0 { b.min.y } else { b.max.y }, if i & 4 == 0 { b.min.z } else { b.max.z });
        speed = speed.max(steady.velocity(p).length());
    }
    speed
}

/// Which way the motion carries `face` (as its shell uses it): from `V·n` at points inside
/// it, the largest either way.
fn side(model: &Model, face: Shape, steady: &Steady, speed: f64) -> Result<Side> {
    let surface = model.surface(face).ok_or_else(|| Error::InvalidInput("a face without a surface".into()))?;
    let polygons = contour::loop_polygons(model, face).ok_or_else(|| Error::Failed("a face without pcurves".into()))?;
    let gauge = contour::Gauge { surface, steady };
    let sign = if face.orientation == Orientation::Reversed { -1.0 } else { 1.0 };
    let (u0, u1, v0, v1) = model.face_uv_bounds(face);
    const N: usize = 24;
    let (mut most, mut least) = (0.0f64, 0.0f64);
    for i in 0..N {
        for j in 0..N {
            let q = Point2::new(u0 + (u1 - u0) * (i as f64 + 0.5) / N as f64, v0 + (v1 - v0) * (j as f64 + 0.5) / N as f64);
            if !contour::within(&polygons, q) {
                continue;
            }
            let g = sign * gauge.g(q.x, q.y);
            most = most.max(g);
            least = least.min(g);
        }
    }
    // (V·n of a face the motion runs along is zero but for rounding. Points near a piece's
    // boundary may fall just across the curve bounding it, as the loops are sampled: the
    // piece goes by the larger.)
    let still = 1e-9 * speed;
    Ok(if most > still && most >= -least { Side::Front } else { Side::Back })
}

/// Builds the product's faces: the paths of vertices, the copies at the end, the faces edges
/// sweep.
///
/// A screw is built in stages of at most a quarter turn (each the same screw, carrying the
/// copies the stage before ended at), so that no face winds round its axis: each is a strip
/// short enough that a point's place on it is found from any start on it.
struct Builder<'m> {
    model: &'m mut Model,
    /// The motion of each stage.
    stages: Vec<Steady>,
    /// Where the stages start and the last ends: the motion's places.
    junctions: Vec<Transform>,
    /// Each stage, where a screw or a carried motion, as the pipes' motions sample it (its
    /// points' paths quintics through them).
    motions: Vec<Option<PieceMotion>>,
    size: f64,
    /// A turn's axis, whose points stand still.
    still: Option<crate::math::Axis>,
    /// A whole turn: its end is its start, the copies there the body's own.
    whole: bool,
    vertex_copies: HashMap<(VertexId, usize), VertexId>,
    edge_copies: HashMap<(EdgeId, usize), EdgeId>,
    paths: HashMap<(VertexId, usize), EdgeId>,
}

impl<'m> Builder<'m> {
    fn new(model: &'m mut Model, stages: &[Steady], junctions: Vec<Transform>, size: f64, reach: f64) -> Result<Builder<'m>> {
        let mut motions: Vec<Option<PieceMotion>> = Vec::with_capacity(stages.len());
        for stage in stages {
            motions.push(match stage {
                Steady::Screw { axis, angle, rise } => Some(screw_motion(*axis, *angle, *rise, size, reach)?),
                Steady::Carried(stage) => Some(carried_motion(stage, size)?),
                _ => None,
            });
        }
        Ok(Builder { model, stages: stages.to_vec(), junctions, motions, size, still: still_axis(stages), whole: is_whole_turn(stages), vertex_copies: HashMap::new(), edge_copies: HashMap::new(), paths: HashMap::new() })
    }

    /// Whether the vertex lies on the turn's axis, standing still.
    fn stands_still(&self, v: VertexId) -> bool {
        self.still.is_some_and(|axis| off_axis(axis, self.model[v].point) <= tolerance::LINEAR.max(self.model[v].tolerance))
    }

    /// Whether the edge lies along the turn's axis (it sweeps nothing).
    fn along_axis(&self, e: EdgeId) -> bool {
        let Some(axis) = self.still else { return false };
        let (a, b) = self.model.edge_end_vertices(e);
        if !(self.stands_still(a) && self.stands_still(b)) || self.model.is_degenerate(Shape::from(e)) {
            return false;
        }
        let Some((curve, first, last)) = self.model.edge_curve(Shape::from(e)) else { return false };
        (0..=8).all(|i| off_axis(axis, curve.point(first + (last - first) * i as f64 / 8.0)) <= tolerance::LINEAR.max(self.model[e].tolerance))
    }

    fn product(&mut self, faces: &[Shape], sides: &[Vec<Side>]) -> Result<Shape> {
        let shell = self.model.new_empty(NodeKind::Shell);
        let n = self.stages.len();
        // The faces outside W where the motion starts.
        for (i, &face) in faces.iter().enumerate() {
            if sides[0][i] == Side::Back {
                self.model.add(shell, face)?;
            }
        }
        for k in 0..n {
            // What the boundary of W sweeps over the stage.
            for (e, along) in boundary_edges(self.model, faces, &sides[k])? {
                if self.along_axis(e) {
                    continue;
                }
                let face = self.side_face(e, along, k)?;
                self.model.add(shell, face)?;
            }
            // At the junction after it: the faces passing out of W (facing out) or into it
            // (facing in); at the end, W.
            for (i, &face) in faces.iter().enumerate() {
                let next = if k + 1 < n { Some(sides[k + 1][i]) } else { None };
                match (sides[k][i], next) {
                    (Side::Front, None) | (Side::Front, Some(Side::Back)) => {
                        let copy = self.face_copy(face, k + 1)?;
                        self.model.add(shell, copy)?;
                    }
                    (Side::Back, Some(Side::Front)) => {
                        let copy = self.face_copy(face, k + 1)?;
                        self.model.add(shell, copy.reversed())?;
                    }
                    _ => {}
                }
            }
        }
        let volume = crate::properties::closed_shell_volume(self.model, shell);
        let solid = self.model.new_empty(NodeKind::Solid);
        self.model.add(solid, if volume < 0.0 { shell.reversed() } else { shell })?;
        Ok(solid)
    }

    /// A vertex's copy at junction `k` (itself at the start).
    fn vertex_copy(&mut self, v: VertexId, k: usize) -> VertexId {
        if k == 0 || self.whole || self.stands_still(v) {
            return v;
        }
        if let Some(&c) = self.vertex_copies.get(&(v, k)) {
            return c;
        }
        let (p, tol) = (self.model[v].point, self.model[v].tolerance);
        let c = self.model.add_vertex(p.transformed(&self.junctions[k]), tol);
        self.vertex_copies.insert((v, k), c);
        c
    }

    /// An edge's copy at junction `k` (itself at the start).
    fn edge_copy(&mut self, e: EdgeId, k: usize) -> Result<EdgeId> {
        if k == 0 || self.whole || self.along_axis(e) {
            return Ok(e);
        }
        if let Some(&c) = self.edge_copies.get(&(e, k)) {
            return Ok(c);
        }
        let (a, b) = self.model.edge_end_vertices(e);
        let (a, b) = (self.vertex_copy(a, k), self.vertex_copy(b, k));
        let tol = self.model[e].tolerance;
        let c = self.model.add_edge(a, b, tol)?;
        if let Some((curve, first, last)) = self.model.edge_curve(Shape::from(e)) {
            self.model.set_edge_curve(c.into(), curve.transformed(&self.junctions[k]), tol)?;
            self.model.set_edge_range(c.into(), first, last)?;
        }
        self.edge_copies.insert((e, k), c);
        Ok(c)
    }

    /// The copy at junction `j` of a face, used as the body uses it (its edges the copies of
    /// its edges, its pcurves the same on the moved surface).
    fn face_copy(&mut self, face: Shape, j: usize) -> Result<Shape> {
        if self.whole {
            return Ok(face);
        }
        let id = FaceId::from_index(face.node.index);
        let tol = self.model[id].tolerance;
        let end = self.junctions[j];
        let stages = j;
        let surface = self.model.face_surface(id).ok_or_else(|| Error::InvalidInput("a face without a surface".into()))?.transformed(&end);
        let copy = self.model.new_empty(NodeKind::Face);
        self.model.set_face_surface(copy, surface, tol)?;
        for wire_use in self.model[id].wires.clone() {
            let (wire, wire_orientation) = (self.model[wire_use].wire, self.model[wire_use].orientation);
            let mut coedges = Vec::new();
            for coedge in self.model[wire].coedges.clone() {
                let (edge, orientation) = (self.model[coedge].edge, self.model[coedge].orientation);
                let moved = self.edge_copy(edge, stages)?;
                let c = self.model.add_coedge(moved, orientation)?;
                if let Some((pcurve, first, last)) = self.model.effective_pcurve(coedge) {
                    self.model.set_pcurve_with_range(c, pcurve, first, last)?;
                }
                coedges.push(c);
            }
            let w = self.model.add_wire(&coedges)?;
            self.model.add(copy, Shape { orientation: wire_orientation, ..Shape::from(w) })?;
        }
        Ok(Shape { orientation: face.orientation, ..copy })
    }

    /// The edge a vertex's path makes over stage `k`, from its copy at junction `k` to that at
    /// the next.
    fn vertex_path(&mut self, v: VertexId, k: usize) -> Result<EdgeId> {
        let tol = self.model[v].tolerance.max(tolerance::LINEAR);
        // (A vertex standing still on a turn's axis: a degenerate edge, each face its own.)
        if self.stands_still(v) {
            return self.model.add_edge(v, v, tol);
        }
        if let Some(&e) = self.paths.get(&(v, k)) {
            return Ok(e);
        }
        let (start, end) = (self.vertex_copy(v, k), self.vertex_copy(v, k + 1));
        let p = self.model[start].point;
        let edge = self.model.add_edge(start, end, tol)?;
        let (curve, first, last) = match self.stages[k] {
            Steady::Slide { by } => (Curve::Line(Line::new(p, Direction::from_vector(by))), 0.0, by.length()),
            Steady::Turn { axis, angle } => (Curve::Circle(turn_circle(axis, p)?), 0.0, angle),
            Steady::Screw { .. } | Steady::Carried(_) => (Curve::BSpline(path(self.motions[k].as_ref().expect("a screw's motion"), p)?), 0.0, 1.0),
        };
        self.model.set_edge_curve(edge.into(), curve, tol)?;
        self.model.set_edge_range(edge.into(), first, last)?;
        self.paths.insert((v, k), edge);
        Ok(edge)
    }

    /// The face edge `e` sweeps over stage `k`, used so that it bounds the sweep of the face
    /// in `W` beside it: that face runs along the edge's curve if `along`.
    fn side_face(&mut self, e: EdgeId, along: bool, k: usize) -> Result<Shape> {
        let (start, finish) = self.model.edge_end_vertices(e);
        let (start_path, end_path) = (self.vertex_path(start, k)?, self.vertex_path(finish, k)?);
        let top = self.edge_copy(e, k + 1)?;
        let e = self.edge_copy(e, k)?;
        let (curve, a, b) = self.model.edge_curve(Shape::from(e)).ok_or_else(|| Error::Failed("an edge without a curve".into()))?;
        let line = |x: f64, y: f64, dx: f64, dy: f64| Curve2::Line(Line2::new(Point2::new(x, y), Direction2::new(dx, dy)));
        if let Some(face) = self.planar_side_face(e, along, k, top, start_path, end_path)? {
            return Ok(face);
        }
        // The surface; the (u, v) of the curve's point at s and moved to t; and the uses
        // round the face, counter-clockwise in (u, v).
        type Uses = Vec<(EdgeId, Orientation, Curve2, f64, f64)>;
        let (surface, uv, uses): (Surface, Box<dyn Fn(f64, f64) -> Point2>, Uses) = match self.stages[k].clone() {
            Steady::Slide { by } => {
                let length = by.length();
                match &curve {
                    Curve::Line(l) => {
                        let x = Vector::from(l.direction);
                        let across = by - x * by.dot(x);
                        if across.length() <= tolerance::LINEAR {
                            return Err(Error::Failed("a straight edge slid along itself bounding W".into()));
                        }
                        let plane = Plane::from_frame(Frame::new(l.origin, Direction::from_vector(x.cross(across)), l.direction));
                        let (du, dv) = (by.dot(x), across.length());
                        let path = |s: f64| (line(s, 0.0, du / length, dv / length), 0.0, length);
                        let (p0, p1) = (path(a), path(b));
                        (
                            Surface::Plane(plane),
                            Box::new(move |s, t| Point2::new(s + du * t, dv * t)),
                            vec![
                                (e, Orientation::Forward, line(0.0, 0.0, 1.0, 0.0), a, b),
                                (end_path, Orientation::Forward, p1.0, p1.1, p1.2),
                                (top, Orientation::Reversed, line(du, dv, 1.0, 0.0), a, b),
                                (start_path, Orientation::Reversed, p0.0, p0.1, p0.2),
                            ],
                        )
                    }
                    Curve::Circle(c) if Vector::from(c.frame().z()).cross(by).length() <= 1e-12 * length => {
                        let height = by.dot(Vector::from(c.frame().z()));
                        let uses = rectangle_uses(e, top, start_path, end_path, a, b, height, true);
                        (Surface::Cylinder(Cylinder::new(*c.frame(), c.radius())), Box::new(move |s, t| Point2::new(s, height * t)), uses)
                    }
                    other => {
                        // The extrusion as the B-spline surface it is: the curve's poles and
                        // their copies at the end, joined straight.
                        let profile = Profile::of(other, a, b, self.size)?;
                        let surface = profiled_surface(&profile.spline, &|q| Ok(vec![q, q + by]), KnotVector::new(vec![0.0, length], vec![2, 2], 1))?;
                        let uses = profile.uses(e, top, start_path, end_path, a, b, length);
                        (surface, Box::new(move |s, t| Point2::new(profile.at(s), length * t)), uses)
                    }
                }
            }
            Steady::Turn { axis, angle } if turned_on_sphere(&curve, a, b, axis, self.size).is_some() => {
                // A circle whose axis meets the turn's lies on the sphere about where they meet
                // through it, and turns on it (a pin's end turned about a far axis square to it:
                // as a surface of revolution of the circle, its two halves would sweep the same
                // points, and its parameters not be had from a point). u the angle turned and
                // v the latitude: the curve's pcurve, moved along u as it turns; its ends'
                // paths are parallels.
                let (sphere, pcurve) = turned_on_sphere(&curve, a, b, axis, self.size).expect("checked");
                let (p0, p1) = (pcurve.point(a), pcurve.point(b));
                let shifted = crate::math::Transformable2::transformed(&pcurve, &crate::math::Transform2::translation(crate::math::Vector2::new(angle, 0.0)));
                let (f, r) = (Orientation::Forward, Orientation::Reversed);
                let mut uses = vec![
                    (e, f, Curve2::BSpline(pcurve.clone()), a, b),
                    (end_path, f, line(p1.x, p1.y, 1.0, 0.0), 0.0, angle),
                    (top, r, Curve2::BSpline(shifted), a, b),
                    (start_path, r, line(p0.x, p0.y, 1.0, 0.0), 0.0, angle),
                ];
                // (Counter-clockwise in (s, t); the map to (u, v) turns it where the latitude
                // grows along the curve.)
                let (_, d) = pcurve.d1(0.5 * (a + b));
                if d.y > 0.0 {
                    uses.reverse();
                    for u in &mut uses {
                        u.1 = u.1.reversed();
                    }
                }
                let at = pcurve.clone();
                (Surface::Sphere(sphere), Box::new(move |s, t| {
                    let q = at.point(s);
                    Point2::new(q.x + angle * t, q.y)
                }), uses)
            }
            Steady::Turn { axis, angle } => match crate::sweep::turned_surface(axis, &curve, a, b) {
                // The elementary surface the curve turns into where it is one (a cylinder,
                // cone, sphere or torus about the axis), else its surface of revolution (an
                // unbounded curve, as a line, trimmed to the edge): u the angle turned from
                // the curve and v = v0 + sigma·s.
                crate::sweep::SweptSurface::Ruled { surface, v0, sigma } => {
                    let (f, r) = (Orientation::Forward, Orientation::Reversed);
                    let v = move |s: f64| v0 + sigma * s;
                    let mut uses = vec![
                        (e, f, line(0.0, v0, 0.0, sigma), a, b),
                        (end_path, f, line(0.0, v(b), 1.0, 0.0), 0.0, angle),
                        (top, r, line(angle, v0, 0.0, sigma), a, b),
                        (start_path, r, line(0.0, v(a), 1.0, 0.0), 0.0, angle),
                    ];
                    // (Clockwise in (u, v) where v grows along the curve.)
                    if sigma > 0.0 {
                        uses.reverse();
                        for u in &mut uses {
                            u.1 = u.1.reversed();
                        }
                    }
                    (surface, Box::new(move |s, t| Point2::new(angle * t, v(s))), uses)
                }
                _ => {
                    let revolution = RevolutionSurface::new(curve.clone(), axis);
                    let uses = rectangle_uses(e, top, start_path, end_path, a, b, angle, false);
                    (Surface::Revolution(revolution), Box::new(move |s, t| Point2::new(angle * t, s)), uses)
                }
            },
            Steady::Screw { axis, angle, rise } if screwed_along_itself(&curve, axis) => {
                // A line along the axis screwed slides along itself as it turns: the cylinder
                // about the axis through it, u the angle turned from the line and v the height
                // above the line's origin. The paths of its ends are helices, straight there:
                // (angle·t, v + rise·t).
                let Curve::Line(l) = &curve else { unreachable!() };
                let z = Vector::from(axis.direction);
                let foot = axis.origin + z * (l.origin - axis.origin).dot(z);
                let radial = l.origin - foot;
                let cylinder = Cylinder::new(Frame::new(foot, axis.direction, Direction::from_vector(radial)), radial.length());
                let sigma = Vector::from(l.direction).dot(z).signum();
                let path = |s: f64| {
                    let from = Point2::new(0.0, sigma * s);
                    Curve2::BSpline(crate::geometry::curve2::BSplineCurve2::try_new(vec![from, Point2::new(angle, sigma * s + rise)], vec![0.0, 1.0], vec![2, 2], 1).expect("a segment"))
                };
                let (f, r) = (Orientation::Forward, Orientation::Reversed);
                let mut uses = vec![
                    (e, f, line(0.0, 0.0, 0.0, sigma), a, b),
                    (end_path, f, path(b), 0.0, 1.0),
                    (top, r, line(angle, rise, 0.0, sigma), a, b),
                    (start_path, r, path(a), 0.0, 1.0),
                ];
                // (Counter-clockwise in (s, t); the map to (u, v) turns it where -sigma·angle
                // is negative.)
                if sigma > 0.0 {
                    uses.reverse();
                    for u in &mut uses {
                        u.1 = u.1.reversed();
                    }
                }
                (Surface::Cylinder(cylinder), Box::new(move |s, t| Point2::new(angle * t, sigma * s + rise * t)), uses)
            }
            Steady::Screw { .. } | Steady::Carried(_) => {
                let profile = Profile::of(&curve, a, b, self.size)?;
                let motion = self.motions[k].as_ref().expect("a screw's motion");
                let knots = path(motion, Point::ORIGIN)?;
                let knots = KnotVector::new(knots.knots().to_vec(), knots.multiplicities().to_vec(), knots.degree());
                let swept = profiled_surface(&profile.spline, &|q| Ok(path(motion, q)?.poles().to_vec()), knots)?;
                let uses = profile.uses(e, top, start_path, end_path, a, b, 1.0);
                (swept, Box::new(move |s, t| Point2::new(profile.at(s), t)), uses)
            }
        };
        // Facing as the boundary of the sweep of the face beside it: along (the way that face
        // runs the edge) crossed with the motion.
        let (sm, tm) = (0.5 * (a + b), 0.5);
        let ds = 1e-6 * (b - a);
        let moved = |s: f64, t: f64| curve.point(s).transformed(&self.stages[k].at(t));
        let along_s = (moved(sm + ds, tm) - moved(sm - ds, tm)) * if along { 1.0 } else { -1.0 };
        let along_t = moved(sm, tm + 1e-6) - moved(sm, tm - 1e-6);
        let q = uv(sm, tm);
        let (_, su, sv) = surface.d1(q.x, q.y);
        let orientation = if su.cross(sv).dot(along_s.cross(along_t)) >= 0.0 { Orientation::Forward } else { Orientation::Reversed };
        let face = self.model.new_empty(NodeKind::Face);
        self.model.set_face_surface(face, surface, tolerance::LINEAR)?;
        let mut coedges = Vec::with_capacity(uses.len());
        for (edge, o, pcurve, first, last) in uses {
            let coedge = self.model.add_coedge(edge, o)?;
            self.model.set_pcurve_with_range(coedge, pcurve, first, last)?;
            coedges.push(coedge);
        }
        let wire = self.model.add_wire(&coedges)?;
        self.model.add(face, wire.into())?;
        Ok(Shape { orientation, ..face })
    }

    /// A side face of a turn on a plane (square to the axis), bounded by `loops`, each put
    /// counter-clockwise (the outer) or clockwise (a hole) in the plane's parameters; facing
    /// as [`Builder::side_face`] faces it.
    #[allow(clippy::too_many_arguments)]
    fn planar_loops_face(&mut self, surface: Surface, mut loops: Vec<Vec<(EdgeId, Orientation, Curve2, f64, f64)>>, uv: Box<dyn Fn(f64, f64) -> Point2 + '_>, curve: &Curve, along: bool, k: usize, a: f64, b: f64) -> Result<Shape> {
        // Each loop's signed area in the parameters, from points along its pcurves.
        let area = |uses: &[(EdgeId, Orientation, Curve2, f64, f64)]| {
            let mut points = Vec::new();
            for (_, o, pc, first, last) in uses {
                for i in 0..32 {
                    let f = i as f64 / 32.0;
                    let t = if *o == Orientation::Forward { first + (last - first) * f } else { last - (last - first) * f };
                    points.push(pc.point(t));
                }
            }
            let n = points.len();
            (0..n).map(|i| points[i].x * points[(i + 1) % n].y - points[(i + 1) % n].x * points[i].y).sum::<f64>() * 0.5
        };
        let areas: Vec<f64> = loops.iter().map(|l| area(l)).collect();
        let outer = (0..loops.len()).max_by(|&i, &j| areas[i].abs().total_cmp(&areas[j].abs())).ok_or_else(|| Error::Failed("a planar side face without a boundary".into()))?;
        for (i, uses) in loops.iter_mut().enumerate() {
            if (areas[i] > 0.0) != (i == outer) {
                uses.reverse();
                for u in uses.iter_mut() {
                    u.1 = u.1.reversed();
                }
            }
        }
        let (sm, tm) = (0.5 * (a + b), 0.5);
        let ds = 1e-6 * (b - a);
        let moved = |s: f64, t: f64| curve.point(s).transformed(&self.stages[k].at(t));
        let along_s = (moved(sm + ds, tm) - moved(sm - ds, tm)) * if along { 1.0 } else { -1.0 };
        let along_t = moved(sm, tm + 1e-6) - moved(sm, tm - 1e-6);
        let q = uv(sm, tm);
        let (_, su, sv) = surface.d1(q.x, q.y);
        let orientation = if su.cross(sv).dot(along_s.cross(along_t)) >= 0.0 { Orientation::Forward } else { Orientation::Reversed };
        let face = self.model.new_empty(NodeKind::Face);
        self.model.set_face_surface(face, surface, tolerance::LINEAR)?;
        for uses in loops {
            let mut coedges = Vec::with_capacity(uses.len());
            for (edge, o, pcurve, first, last) in uses {
                let coedge = self.model.add_coedge(edge, o)?;
                self.model.set_pcurve_with_range(coedge, pcurve, first, last)?;
                coedges.push(coedge);
            }
            let wire = self.model.add_wire(&coedges)?;
            self.model.add(face, wire.into())?;
        }
        Ok(Shape { orientation, ..face })
    }
}

impl Builder<'_> {
    /// The face edge `e` (its copy at junction `k`) sweeps over stage `k` where that is a
    /// plane: the edge's curve lies in a plane the motion keeps to (square to a turn's axis, or
    /// along a slide). Such a sweep as a surface of revolution or an extrusion is singular
    /// where the curve runs along the motion (as it does at the ends of an edge bounding `W`),
    /// and a plane is what it is. `None` where the sweep is not planar.
    fn planar_side_face(&mut self, e: EdgeId, along: bool, k: usize, top: EdgeId, start_path: EdgeId, end_path: EdgeId) -> Result<Option<Shape>> {
        let (curve, a, b) = self.model.edge_curve(Shape::from(e)).ok_or_else(|| Error::Failed("an edge without a curve".into()))?;
        let stage = self.stages[k].clone();
        let tol = 1e-9 * self.size.max(1.0);
        let samples: Vec<Point> = (0..=16).map(|i| curve.point(a + (b - a) * i as f64 / 16.0)).collect();
        let normal = match stage {
            Steady::Turn { axis, .. } => Vector::from(axis.direction),
            Steady::Slide { by } => {
                if matches!(curve, Curve::Line(_)) {
                    return Ok(None);
                }
                // Square to the slide and to the chord reaching farthest across it.
                let d = by / by.length();
                let across = samples.iter().map(|&p| (p - samples[0]) - d * (p - samples[0]).dot(d)).fold(Vector::ZERO, |m, v| if v.length() > m.length() { v } else { m });
                if across.length() <= tol {
                    return Ok(None);
                }
                d.cross(across)
            }
            Steady::Screw { .. } | Steady::Carried(_) => return Ok(None),
        };
        let Some(n) = Direction::try_from_vector(normal) else { return Ok(None) };
        let nv = Vector::from(n);
        if samples.iter().any(|&p| (p - samples[0]).dot(nv).abs() > tol) {
            return Ok(None);
        }
        let plane = Plane::from_frame(Frame::from_z(samples[0], n));
        let surface = Surface::Plane(plane.clone());
        let (f, r) = (Orientation::Forward, Orientation::Reversed);
        // (A whole turn: the annulus, or the disc, between the circles its ends' paths are.)
        if self.whole {
            let mut loops = Vec::new();
            for edge in [start_path, end_path] {
                let Some((c, first, last)) = self.model.edge_curve(Shape::from(edge)) else { continue };
                let Some(pcurve) = crate::geometry::project::analytic_pcurve(&c, &surface) else { return Ok(None) };
                loops.push(vec![(edge, f, pcurve, first, last)]);
            }
            let moving = curve.clone();
            let uv = move |s: f64, t: f64| {
                let (u, v) = plane.parameters(moving.point(s).transformed(&stage.at(t)));
                Point2::new(u, v)
            };
            return self.planar_loops_face(surface, loops, Box::new(uv), &curve, along, k, a, b).map(Some);
        }
        // Counter-clockwise in (s, t): the edge, the end's path, the copy back, the start's
        // path back (but the path of an end standing still on a turn's axis: a point).
        let mut uses = Vec::with_capacity(4);
        for (edge, o) in [(e, f), (end_path, f), (top, r), (start_path, r)] {
            if self.model.is_degenerate(Shape::from(edge)) {
                continue;
            }
            let Some((c, first, last)) = self.model.edge_curve(Shape::from(edge)) else { return Ok(None) };
            let Some(pcurve) = crate::geometry::project::analytic_pcurve(&c, &surface) else { return Ok(None) };
            uses.push((edge, o, pcurve, first, last));
        }
        // Counter-clockwise about the plane's normal: (s, t) maps to the plane keeping its
        // turn where ∂s × ∂t runs along the normal (it does not change sign within an edge
        // bounding W: the curve runs along the motion only at its ends).
        let (sm, tm) = (0.5 * (a + b), 0.5);
        let ds = 1e-6 * (b - a);
        let moved = |s: f64, t: f64| curve.point(s).transformed(&stage.at(t));
        let along_s = moved(sm + ds, tm) - moved(sm - ds, tm);
        let along_t = moved(sm, tm + 1e-6) - moved(sm, tm - 1e-6);
        let turning = along_s.cross(along_t).dot(nv);
        if turning == 0.0 {
            return Ok(None);
        }
        if turning < 0.0 {
            uses.reverse();
            for u in &mut uses {
                u.1 = u.1.reversed();
            }
        }
        // Facing as the boundary of the sweep of the face beside it (see `side_face`).
        let facing = turning * if along { 1.0 } else { -1.0 };
        let orientation = if facing >= 0.0 { Orientation::Forward } else { Orientation::Reversed };
        let face = self.model.new_empty(NodeKind::Face);
        self.model.set_face_surface(face, surface, tolerance::LINEAR)?;
        let mut coedges = Vec::with_capacity(uses.len());
        for (edge, o, pcurve, first, last) in uses {
            let coedge = self.model.add_coedge(edge, o)?;
            self.model.set_pcurve_with_range(coedge, pcurve, first, last)?;
            coedges.push(coedge);
        }
        let wire = self.model.add_wire(&coedges)?;
        self.model.add(face, wire.into())?;
        Ok(Some(Shape { orientation, ..face }))
    }
}

/// The sphere a circle turned about `axis` stays on where the circle's own axis meets it (not
/// along it: a circle square to the axis turns within its plane), with the circle's pcurve on
/// it over `[a, b]`; `None` otherwise.
fn turned_on_sphere(curve: &Curve, a: f64, b: f64, axis: crate::math::Axis, size: f64) -> Option<(crate::geometry::surface::Sphere, crate::geometry::curve2::BSplineCurve2)> {
    let Curve::Circle(circle) = curve else { return None };
    let tol = 1e-9 * size.max(1.0);
    let z = Vector::from(axis.direction);
    let (c, n) = (circle.center(), Vector::from(circle.frame().z()));
    // The circle's axis, c + λn, against the turn's, o + μz: where they meet.
    let across = n.cross(z);
    if across.length() <= 1e-9 {
        return None;
    }
    let w = c - axis.origin;
    if w.dot(across).abs() > tol * across.length() {
        return None;
    }
    // (c + λn - o - μz ⊥ n and ⊥ z.)
    let nz = n.dot(z);
    let det = nz * nz - 1.0;
    let lambda = (w.dot(n) - nz * w.dot(z)) / det;
    let o = c + n * lambda;
    let radius = ((o - c).length().powi(2) + circle.radius().powi(2)).sqrt();
    let samples = (0..=8).map(|i| circle.point(a + (b - a) * i as f64 / 8.0));
    let first = circle.point(a) - o;
    let radial = first - z * first.dot(z);
    if radial.length() <= tol {
        return None;
    }
    let mut ok = true;
    for p in samples {
        let d = p - o;
        ok &= ((d.length() - radius).abs() <= tol) && (d - z * d.dot(z)).length() > tol;
    }
    if !ok {
        return None;
    }
    let sphere = crate::geometry::surface::Sphere::new(Frame::new(o, axis.direction, Direction::from_vector(radial)), radius);
    let pcurve = crate::geometry::project::pcurve_on_range(curve, &Surface::Sphere(sphere), a, b, tol)?;
    Some((sphere, pcurve))
}

/// Whether `curve` is a line along a screw's `axis`, off it (screwed, it slides along itself).
fn screwed_along_itself(curve: &Curve, axis: crate::math::Axis) -> bool {
    let Curve::Line(l) = curve else { return false };
    let z = Vector::from(axis.direction);
    let d = l.origin - axis.origin;
    Vector::from(l.direction).cross(z).length() <= 1e-12 && (d - z * d.dot(z)).length() > tolerance::LINEAR
}

/// The surface through the paths of `profile`'s poles (each pole's path's poles over `knots`,
/// the same for all): u the profile's parameter, v the motion's. The motion is rigid, so the
/// surface's iso-curves are the profile moved, exactly where the paths are.
fn profiled_surface(profile: &BSplineCurve, path_poles: &dyn Fn(Point) -> Result<Vec<Point>>, knots: KnotVector) -> Result<Surface> {
    let poles: Vec<Vec<Point>> = profile.poles().iter().map(|&q| path_poles(q)).collect::<Result<_>>()?;
    let u = if profile.is_periodic() {
        KnotVector::periodic(profile.knots().to_vec(), profile.multiplicities().to_vec(), profile.degree())
    } else {
        KnotVector::new(profile.knots().to_vec(), profile.multiplicities().to_vec(), profile.degree())
    };
    let rows = poles[0].len();
    let mut swept = BSplineSurface::try_new(poles, u, knots)?;
    if let Some(w) = profile.weights() {
        swept = swept.with_weights(w.iter().map(|&wi| vec![wi; rows]).collect());
    }
    Ok(Surface::BSpline(swept))
}

/// The uses round a side face whose surface has the curve's parameter `s` along one direction
/// and the motion (`0..extent` over the motion) along the other: `s_first` says whether `s` is
/// u. Counter-clockwise in (u, v).
#[allow(clippy::too_many_arguments)]
fn rectangle_uses(e: EdgeId, top: EdgeId, start_path: EdgeId, end_path: EdgeId, a: f64, b: f64, extent: f64, s_first: bool) -> Vec<(EdgeId, Orientation, Curve2, f64, f64)> {
    let line = |x: f64, y: f64, dx: f64, dy: f64| Curve2::Line(Line2::new(Point2::new(x, y), Direction2::new(dx, dy)));
    let (f, r) = (Orientation::Forward, Orientation::Reversed);
    // In (s, m): the edge at m = 0, the end path at s = b, the copy at m = extent, the start
    // path at s = a; counter-clockwise when s is u and the extent positive.
    let mut uses = if s_first {
        vec![
            (e, f, line(0.0, 0.0, 1.0, 0.0), a, b),
            (end_path, f, line(b, 0.0, 0.0, extent.signum()), 0.0, path_range(extent)),
            (top, r, line(0.0, extent, 1.0, 0.0), a, b),
            (start_path, r, line(a, 0.0, 0.0, extent.signum()), 0.0, path_range(extent)),
        ]
    } else {
        // (u the motion, v the curve's parameter.)
        vec![
            (e, f, line(0.0, 0.0, 0.0, 1.0), a, b),
            (end_path, f, line(0.0, b, extent.signum(), 0.0), 0.0, path_range(extent)),
            (top, r, line(extent, 0.0, 0.0, 1.0), a, b),
            (start_path, r, line(0.0, a, extent.signum(), 0.0), 0.0, path_range(extent)),
        ]
    };
    // Counter-clockwise: reversed where the rectangle is mirrored.
    let mirrored = (extent < 0.0) != !s_first;
    if mirrored {
        uses.reverse();
        for u in &mut uses {
            u.1 = u.1.reversed();
        }
    }
    uses
}

/// The parameter range of a vertex's path: its length along a slide (signed extents taken as
/// their size), the angle of a turn, the screw's fraction.
fn path_range(extent: f64) -> f64 {
    extent.abs()
}

/// The circle a turn about `axis` carries `p` round, its parameter the angle turned.
fn turn_circle(axis: crate::math::Axis, p: Point) -> Result<Circle> {
    let d = Vector::from(axis.direction);
    let foot = axis.origin + d * (p - axis.origin).dot(d);
    let radial = p - foot;
    let x = Direction::try_from_vector(radial).filter(|_| radial.length() > tolerance::LINEAR).ok_or_else(|| Error::InvalidInput("a point on the turn's axis".into()))?;
    Ok(Circle::new(Frame::new(foot, axis.direction, x), radial.length()))
}

/// The screw sampled as the pipes' motions are, finely enough that the path of each point
/// within `reach` of the axis (the quintics through its places, velocities and accelerations)
/// lies within 1e-10 of `size` of its helix.
fn screw_motion(axis: crate::math::Axis, angle: f64, rise: f64, size: f64, reach: f64) -> Result<PieceMotion> {
    let along = Vector::from(axis.direction);
    let at = move |f: f64| Transform::translation(along * (rise * f)) * Transform::rotation(axis, angle * f);
    let rates = move |f: f64, p: Point| -> (Vector, Vector) {
        let q = p.transformed(&at(f));
        let r = (q - axis.origin) - along * (q - axis.origin).dot(along);
        let turning = along.cross(r);
        (turning * angle + along * rise, along.cross(turning) * (angle * angle))
    };
    // (Probed at the farthest the body reaches from the axis.)
    let far = axis.origin + Vector::from(perpendicular(axis.direction)) * reach;
    let mut samples = 4usize.max((angle / (std::f64::consts::PI / 8.0)).ceil() as usize);
    loop {
        let params: Vec<f64> = (0..=samples).map(|i| i as f64 / samples as f64).collect();
        let transforms: Vec<Transform> = params.iter().map(|&f| at(f)).collect();
        let motion = PieceMotion { params, transforms, rates: Some(Box::new(rates)) };
        let curve = path(&motion, far)?;
        let error = (0..4 * samples).map(|i| (i as f64 + 0.5) / (4 * samples) as f64).map(|f| curve.point(f).distance(far.transformed(&at(f)))).fold(0.0, f64::max);
        if error <= 1e-10 * size || samples >= 4096 {
            return Ok(motion);
        }
        samples *= 2;
    }
}

/// A stretch of a carried motion sampled as the pipes' motions are (its points' velocities and
/// accelerations by central differences), finely enough that the paths of its box's corners lie
/// within 1e-9 of `size` of where the motion puts them.
fn carried_motion(stage: &std::sync::Arc<super::carried::Stage>, size: f64) -> Result<PieceMotion> {
    // (The probes: points as far from the body as its size, round it.)
    let probes: Vec<Point> = (0..8).map(|i| Point::new(if i & 1 == 0 { -size } else { size }, if i & 2 == 0 { -size } else { size }, if i & 4 == 0 { -size } else { size })).collect();
    let mut samples = 8usize;
    loop {
        let params: Vec<f64> = (0..=samples).map(|i| i as f64 / samples as f64).collect();
        let transforms: Vec<Transform> = params.iter().map(|&f| stage.at(f)).collect();
        let places: Vec<Transform> = transforms.clone();
        let n = samples;
        let at = {
            let stage = std::sync::Arc::clone(stage);
            move |f: f64| stage.at(f)
        };
        // (Velocity and acceleration from the places a little either way.)
        let h = 1e-4 / n as f64;
        let rates = move |f: f64, p: Point| -> (Vector, Vector) {
            let (a, b) = (f - h, f + h);
            let (pa, p0, pb) = (p.transformed(&at(a)), p.transformed(&at(f)), p.transformed(&at(b)));
            ((pb - pa) * (0.5 / h), ((pb - p0) - (p0 - pa)) * (1.0 / (h * h)))
        };
        let _ = places;
        let motion = PieceMotion { params, transforms, rates: Some(Box::new(rates)) };
        let mut error = 0.0f64;
        for &q in &probes {
            let curve = path(&motion, q)?;
            for i in 0..4 * n {
                let f = (i as f64 + 0.5) / (4 * n) as f64;
                error = error.max(curve.point(f).distance(q.transformed(&stage.at(f))));
            }
        }
        if error <= 1e-9 * size || samples >= 2048 {
            return Ok(motion);
        }
        samples *= 2;
    }
}

fn perpendicular(d: Direction) -> Direction {
    let v = Vector::from(d);
    let other = if v.x.abs() < 0.9 { Vector::new(1.0, 0.0, 0.0) } else { Vector::new(0.0, 1.0, 0.0) };
    Direction::from_vector(v.cross(other))
}

/// An edge's curve over `[a, b]` as a B-spline: exactly where it converts so (a circle's
/// arcs as rational quadratics), else cubics through its points within 1e-10 of the size. Where
/// the B-spline's parameter is not the edge's, `sigma` maps the edge's to it.
struct Profile {
    spline: BSplineCurve,
    /// The B-spline's parameter (as the first coordinate) against the edge's, within 1e-12 of
    /// the size in space; `None` where they are one.
    sigma: Option<crate::geometry::curve2::BSplineCurve2>,
}

impl Profile {
    fn of(curve: &Curve, a: f64, b: f64, size: f64) -> Result<Profile> {
        let exact = crate::geometry::convert::curve_to_bspline(&Curve::Trimmed(TrimmedCurve::new(curve.clone(), a, b)));
        if let Some(bs) = exact {
            let same = (bs.first_parameter() - a).abs() <= 1e-12 * (1.0 + a.abs()) && (bs.last_parameter() - b).abs() <= 1e-12 * (1.0 + b.abs());
            let keeps = (0..=16).all(|i| {
                let t = a + (b - a) * i as f64 / 16.0;
                curve.point(t).distance(bs.point(t)) <= 1e-12 * (1.0 + size)
            });
            if same && keeps {
                return Ok(Profile { spline: bs, sigma: None });
            }
            if let Some(sigma) = reparametrization(curve, a, b, &bs, size) {
                return Ok(Profile { spline: bs, sigma: Some(sigma) });
            }
        }
        let mut n = 32;
        loop {
            let params: Vec<f64> = (0..=n).map(|i| a + (b - a) * i as f64 / n as f64).collect();
            let points: Vec<Point> = params.iter().map(|&t| curve.point(t)).collect();
            let fit = crate::geometry::fit::interpolate_with_parameters(&points, &params)?;
            let error = (0..n).map(|i| a + (b - a) * (i as f64 + 0.5) / n as f64).map(|t| curve.point(t).distance(fit.point(t))).fold(0.0, f64::max);
            if error <= 1e-10 * size || n >= 8192 {
                return Ok(Profile { spline: fit, sigma: None });
            }
            n *= 2;
        }
    }

    /// The B-spline's parameter at the edge's `s`.
    fn at(&self, s: f64) -> f64 {
        match &self.sigma {
            Some(sigma) => sigma.point(s).x,
            None => s,
        }
    }

    /// The uses round the side face of a surface with the B-spline's parameter along u and the
    /// motion (`0..extent`) along v ([`rectangle_uses`], the edges' pcurves mapped).
    #[allow(clippy::too_many_arguments)]
    fn uses(&self, e: EdgeId, top: EdgeId, start_path: EdgeId, end_path: EdgeId, a: f64, b: f64, extent: f64) -> Vec<(EdgeId, Orientation, Curve2, f64, f64)> {
        let mut uses = rectangle_uses(e, top, start_path, end_path, a, b, extent, true);
        if let Some(sigma) = &self.sigma {
            let (ua, ub) = (self.at(a), self.at(b));
            let raised = {
                let poles = sigma.poles().iter().map(|q| Point2::new(q.x, q.y + extent)).collect();
                crate::geometry::curve2::BSplineCurve2::new(poles, sigma.knots().to_vec(), sigma.multiplicities().to_vec(), sigma.degree())
            };
            let line = |x: f64, dy: f64| Curve2::Line(Line2::new(Point2::new(x, 0.0), Direction2::new(0.0, dy)));
            for u in &mut uses {
                if u.0 == e {
                    u.2 = Curve2::BSpline(sigma.clone());
                } else if u.0 == top {
                    u.2 = Curve2::BSpline(raised.clone());
                } else if u.0 == start_path {
                    u.2 = line(ua, extent.signum());
                } else if u.0 == end_path {
                    u.2 = line(ub, extent.signum());
                }
            }
        }
        uses
    }
}

/// The map from `curve`'s parameter over `[a, b]` to that of `spline` (the same points), as
/// cubics through samples, within 1e-12 of `size` in space, in pieces between the points of
/// the spline's knots (the map is smooth only between them); `None` if it does not settle.
fn reparametrization(curve: &Curve, a: f64, b: f64, spline: &BSplineCurve, size: f64) -> Option<crate::geometry::curve2::BSplineCurve2> {
    let (u0, u1) = (spline.first_parameter(), spline.last_parameter());
    let hair = 1e-12 * (1.0 + size);
    // Newton on the distance: the parameter on one curve of the other's point.
    let newton = |target: Point, guess: f64, lo: f64, hi: f64, on: &dyn Fn(f64) -> (Point, Vector)| -> Option<f64> {
        let mut x = guess.clamp(lo, hi);
        for _ in 0..50 {
            let (p, d) = on(x);
            let dd = d.dot(d);
            if !(dd > 0.0) {
                return None;
            }
            let step = (target - p).dot(d) / dd;
            x = (x + step).clamp(lo, hi);
            if step.abs() <= 1e-15 * (hi - lo) {
                break;
            }
        }
        (on(x).0.distance(target) <= hair).then_some(x)
    };
    let on_spline = |u: f64| spline.d1(u);
    let on_curve = |s: f64| curve.d1(s);
    // The pieces: between the edge's parameters of the spline's knots.
    let mut breaks = vec![(a, u0)];
    for &k in spline.knots() {
        if k > u0 + 1e-12 * (u1 - u0) && k < u1 - 1e-12 * (u1 - u0) {
            let guess = a + (b - a) * (k - u0) / (u1 - u0);
            breaks.push((newton(spline.point(k), guess, a, b, &on_curve)?, k));
        }
    }
    breaks.push((b, u1));
    let mut pieces: Vec<crate::geometry::curve2::BSplineCurve2> = Vec::new();
    for w in breaks.windows(2) {
        let ((sa, ua), (sb, ub)) = (w[0], w[1]);
        let mut n = 8;
        let piece = loop {
            if n > 4096 {
                return None;
            }
            let params: Vec<f64> = (0..=n).map(|i| sa + (sb - sa) * i as f64 / n as f64).collect();
            let mut us: Vec<f64> = Vec::with_capacity(params.len());
            for (i, &x) in params.iter().enumerate() {
                let u = if i == 0 {
                    ua
                } else if i == n {
                    ub
                } else {
                    newton(curve.point(x), us[i - 1] + (ub - ua) / n as f64, ua, ub, &on_spline)?
                };
                us.push(u);
            }
            let points: Vec<Point2> = us.iter().map(|&u| Point2::new(u, 0.0)).collect();
            let sigma = crate::geometry::fit::interpolate_with_parameters(&points, &params).ok()?;
            let error = (0..n).map(|i| sa + (sb - sa) * (i as f64 + 0.5) / n as f64).map(|x| curve.point(x).distance(spline.point(sigma.point(x).x.clamp(u0, u1)))).fold(0.0, f64::max);
            if error <= hair {
                break sigma;
            }
            n *= 2;
        };
        pieces.push(piece);
    }
    joined_pieces(&pieces)
}

/// Cubic B-splines end to end (each clamped, each starting where the one before ends) as one,
/// continuous across the joins.
fn joined_pieces(pieces: &[crate::geometry::curve2::BSplineCurve2]) -> Option<crate::geometry::curve2::BSplineCurve2> {
    let degree = pieces.first()?.degree();
    if pieces.iter().any(|p| p.degree() != degree) {
        return None;
    }
    let (mut poles, mut knots, mut multiplicities): (Vec<Point2>, Vec<f64>, Vec<usize>) = (Vec::new(), Vec::new(), Vec::new());
    for (k, piece) in pieces.iter().enumerate() {
        let (pk, pm) = (piece.knots(), piece.multiplicities());
        let skip = if k == 0 { 0 } else { 1 };
        poles.extend_from_slice(&piece.poles()[skip..]);
        for (i, (&x, &m)) in pk.iter().zip(pm).enumerate() {
            if i == 0 && k > 0 {
                // (The join: as often as the degree.)
                *multiplicities.last_mut()? = degree;
                continue;
            }
            knots.push(x);
            multiplicities.push(m);
        }
    }
    crate::geometry::curve2::BSplineCurve2::try_new(poles, knots, multiplicities, degree).ok()
}
