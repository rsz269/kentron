//! Offsets whose faces pass through one another: an offset closing a notch narrower than
//! twice it moves the notch's faces through each other, and the faces around it (moved
//! and trimmed by their neighbours as if nothing closed) bound the space twice there, or
//! not at all. The offset is the space they wind round at least once.
//!
//! A face whose loops cross themselves is cut into simple pieces (a face of its plane
//! divided along its own edges), each covering the face as many times as its loops wind
//! round it; a face whose loops run the other way round (moved past its neighbours) covers
//! its plane the other way. All the pieces are split by one another into the cells they bound, and the
//! cells they wind round (along a ray from inside, the pieces left less those entered,
//! each counted as often as it covers) are joined.

use crate::error::Result;
use crate::geometry::surface::Surface;
use crate::heal::face_loops;
use crate::math::Point2;
use crate::topology::{Model, NodeKind, Orientation, Shape};

/// `solid`, an offset, where its faces pass through one another: the space they wind round
/// at least once; `solid` itself where no face's loops cross themselves.
pub(super) fn regularized(model: &mut Model, solid: Shape) -> Result<Shape> {
    if solid.kind() != NodeKind::Solid {
        // (A compound of one solid: that solid, regularized.)
        return match model.unique_subshapes(solid, NodeKind::Solid)[..] {
            [one] if solid.kind() == NodeKind::Compound => {
                let done = regularized(model, one)?;
                if done == one {
                    return Ok(solid);
                }
                let compound = model.new_empty(NodeKind::Compound);
                model.add(compound, done)?;
                Ok(compound)
            }
            _ => Ok(solid),
        };
    }
    let faces: Vec<Shape> = model.explore(solid, NodeKind::Face).collect();
    let mut found: Vec<Vec<(Shape, i32)>> = Vec::with_capacity(faces.len());
    for &face in &faces {
        found.push(simple_pieces(model, face)?);
    }
    let mut pieces: Vec<(Shape, i32)> = Vec::new();
    let mut split = false;
    for f in found {
        split |= f.len() != 1 || f[0].1 != 1;
        pieces.extend(f);
    }
    if !split {
        return Ok(solid);
    }
    // The cells the pieces bound.
    let cell_faces: Vec<Shape> = pieces.iter().map(|&(f, _)| f).collect();
    let cells = crate::boolean::enclosed_solids(model, &cell_faces)?;
    let cells = model.unique_subshapes(cells, NodeKind::Solid);
    // The pieces as often as they cover (a copy for each time more), turned where they
    // cover the other way.
    let counted = model.new_empty(NodeKind::Compound);
    for &(face, times) in &pieces {
        let face = if times < 0 { face.reversed() } else { face };
        model.add(counted, face)?;
        for _ in 1..times.unsigned_abs() {
            let copy = crate::modify::copy(model, face, false)?;
            model.add(counted, copy)?;
        }
    }
    let classifier = crate::boolean::PointClassifier::new(model, counted)?;
    let mut kept = Vec::new();
    for &cell in &cells {
        let Some(p) = crate::boolean::point_inside(model, cell)? else { return Ok(solid) };
        match classifier.winding(p) {
            Some(w) if w >= 1 => kept.push(cell),
            Some(_) => {}
            None => return Ok(solid),
        }
    }
    if kept.is_empty() {
        return Err(crate::error::Error::InvalidInput("an offset leaving nothing (the offset is too large)".into()));
    }
    joined(model, &kept, true)
}

/// The cells joined: the faces two of them share are inside. Faces left in pieces are
/// merged where they lie on one surface (`same_surface_only`), or on the same geometry.
pub(super) fn joined(model: &mut Model, cells: &[Shape], same_surface_only: bool) -> Result<Shape> {
    let mut uses: std::collections::HashMap<crate::topology::NodeId, (Shape, usize)> = std::collections::HashMap::new();
    let mut order = Vec::new();
    for &cell in cells {
        for face in model.explore(cell, NodeKind::Face).collect::<Vec<_>>() {
            let entry = uses.entry(face.node).or_insert_with(|| {
                order.push(face.node);
                (face, 0)
            });
            entry.1 += 1;
        }
    }
    let shell = model.new_empty(NodeKind::Shell);
    for node in order {
        let (face, count) = uses[&node];
        if count == 1 {
            model.add(shell, face)?;
        }
    }
    let joined = model.new_empty(NodeKind::Solid);
    model.add(joined, shell)?;
    let options = crate::unify::UnifyOptions { same_surface_only, ..Default::default() };
    crate::unify::unify_same_domain_with(model, joined, options)
}

/// `face` in simple pieces where its loops cross themselves (on a plane), each with how many
/// times the loops wind round it (pieces they do not cover left out); the face once
/// otherwise.
fn simple_pieces(model: &mut Model, face: Shape) -> Result<Vec<(Shape, i32)>> {
    let once = Ok(vec![(face, 1)]);
    let Some(surface @ Surface::Plane(_)) = model.surface(face).cloned() else { return once };
    let Some(loops) = face_loops(model, face) else { return once };
    let polygons: Vec<Vec<Point2>> = loops.iter().map(|l| l.iter().flat_map(|u| (0..32).map(move |i| u.at(i as f64 / 32.0))).collect()).collect();
    if !crosses_itself(&polygons) {
        // (Its loops turned round, as where an offset moved the face past its neighbours:
        // it covers its surface the other way.)
        let area: f64 = polygons.iter().map(|p| signed_area(p)).sum();
        return Ok(vec![(face, if area < 0.0 { -1 } else { 1 })]);
    }
    // A face of the plane covering the loops, divided along their edges.
    let (mut lo, mut hi) = (Point2::new(f64::INFINITY, f64::INFINITY), Point2::new(f64::NEG_INFINITY, f64::NEG_INFINITY));
    for p in polygons.iter().flatten() {
        lo = Point2::new(lo.x.min(p.x), lo.y.min(p.y));
        hi = Point2::new(hi.x.max(p.x), hi.y.max(p.y));
    }
    let margin = 0.1 * (hi - lo).length().max(1e-9);
    let cover = crate::make::face_on_domain(model, surface, lo.x - margin, hi.x + margin, lo.y - margin, hi.y + margin)?;
    let mut cuts = Vec::new();
    for u in loops.iter().flatten() {
        let Some((curve, a, b)) = model.edge_curve(Shape::from(u.edge)) else { continue };
        let (first, last) = (u.first.min(u.last), u.first.max(u.last));
        if (first - a.min(b)).abs() > 1e-9 * (1.0 + first.abs()) || (last - a.max(b)).abs() > 1e-9 * (1.0 + last.abs()) {
            // (A pcurve over another range than its edge's curve: not cut along.)
            return once;
        }
        cuts.push(crate::divide::Cut { pcurve: u.pcurve.clone(), curve, first, last });
    }
    let divided = crate::divide::divide_face(model, cover, &cuts)?;
    // (The pieces lie forward on the face's plane, and are used as the face is.)
    let mut out = Vec::new();
    for piece in divided {
        let Some(q) = interior_point(model, piece) else { return once };
        let times = winding(&polygons, q);
        if times != 0 {
            out.push((Shape { orientation: face.orientation, ..piece }, times));
        }
    }
    Ok(out)
}

/// Whether the closed polygons' sides cross (other than where neighbours meet).
fn crosses_itself(polygons: &[Vec<Point2>]) -> bool {
    let sides: Vec<(usize, usize, Point2, Point2)> = polygons
        .iter()
        .enumerate()
        .flat_map(|(k, p)| (0..p.len()).map(move |i| (k, i, p[i], p[(i + 1) % p.len()])))
        .collect();
    let cross = |a: Point2, b: Point2, c: Point2| (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
    for (i, &(k1, j1, a, b)) in sides.iter().enumerate() {
        for &(k2, j2, c, d) in &sides[i + 1..] {
            let n = polygons[k1].len();
            if k1 == k2 && (j1 + 1) % n == j2 || k1 == k2 && (j2 + 1) % n == j1 {
                continue;
            }
            let (d1, d2, d3, d4) = (cross(c, d, a), cross(c, d, b), cross(a, b, c), cross(a, b, d));
            if d1 * d2 < 0.0 && d3 * d4 < 0.0 {
                return true;
            }
        }
    }
    false
}

/// The polygon's area, positive where it runs anticlockwise.
fn signed_area(p: &[Point2]) -> f64 {
    0.5 * (0..p.len()).map(|i| {
        let (a, b) = (p[i], p[(i + 1) % p.len()]);
        a.x * b.y - b.x * a.y
    }).sum::<f64>()
}

/// How many times the closed polygons wind round `q`.
pub(super) fn winding(polygons: &[Vec<Point2>], q: Point2) -> i32 {
    let mut w = 0;
    for p in polygons {
        for i in 0..p.len() {
            let (a, b) = (p[i], p[(i + 1) % p.len()]);
            let side = (b.x - a.x) * (q.y - a.y) - (b.y - a.y) * (q.x - a.x);
            if a.y <= q.y && b.y > q.y && side > 0.0 {
                w += 1;
            } else if a.y > q.y && b.y <= q.y && side < 0.0 {
                w -= 1;
            }
        }
    }
    w
}

/// A point inside `face` in its parameters: just inside the middle of one of its outer
/// loop's sides (to the left, where a face's loops have it).
fn interior_point(model: &Model, face: Shape) -> Option<Point2> {
    let loops = face_loops(model, Shape { orientation: Orientation::Forward, ..face })?;
    let polygons: Vec<Vec<Point2>> = loops.iter().map(|l| l.iter().flat_map(|u| (0..32).map(move |i| u.at(i as f64 / 32.0))).collect()).collect();
    let outer = polygons.first()?;
    let size = outer.iter().fold(0.0f64, |m, p| m.max((*p - outer[0]).length()));
    for i in 0..outer.len() {
        let (a, b) = (outer[i], outer[(i + 1) % outer.len()]);
        let along = b - a;
        if along.length() <= 1e-12 * size {
            continue;
        }
        let k = 1e-4 * size / along.length();
        let q = Point2::new(0.5 * (a.x + b.x) - along.y * k, 0.5 * (a.y + b.y) + along.x * k);
        if winding(&polygons, q) != 0 {
            return Some(q);
        }
    }
    None
}
