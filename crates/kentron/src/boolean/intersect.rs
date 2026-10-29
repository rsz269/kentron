//! Intersecting the operands: edges with faces (split points, "paves"), then faces with
//! faces (section edges), as OCCT's pave filler does.

use crate::hashing::FastMap as HashMap;
use std::f64::consts::{PI, TAU};

use crate::geometry::curve::{Circle, Curve, Line};
use crate::geometry::fit;
use crate::geometry::intersect::{self, QuadricCurve, QuadricIntersection};
use crate::geometry::surface::{RevolutionSurface, Surface};
use crate::geometry::ParametricCurve;
use crate::math::{Direction, Frame, Point, Vector};

use super::geom::{self, on_tol, TOL};
use super::graph::{Chords, Edge, Face, Graph, Location};

/// The sine of the angle between an edge and a surface below which the edge lies nearly
/// along it (where they meet is found only roughly along the surface).
const TANGENT_SINE: f64 = 0.05;

/// How a curve meets a surface.
enum Contact {
    /// The curve lies on the surface.
    On,
    /// Parameters where it crosses or touches it.
    At(Vec<f64>),
}

/// Whether `curve` over `[first, last]` keeps clear of `face`'s surface (over the face's
/// range, and as far past it as its edges stray): no stretch of the curve between samples
/// (grown by its length and `tolerance`) meets a box of the face's grid (each grown by the
/// grid's spacing).
fn clear_of_face(face: &super::graph::Face, curve: &Curve, first: f64, last: f64, tolerance: f64) -> bool {
    use crate::math::BoundingBox;
    let cells = face_cells(face);
    const M: usize = 32;
    let points: Vec<Point> = (0..=M).map(|k| curve.point(first + (last - first) * k as f64 / M as f64)).collect();
    points.windows(2).all(|w| {
        let mut b = BoundingBox::EMPTY;
        b.add_point(w[0]);
        b.add_point(w[1]);
        !cells.meets(&b.enlarged(w[0].distance(w[1]) + tolerance + on_tol()))
    })
}

/// Whether two faces keep clear of each other: no box of one's grid (see [`face_cells`])
/// meets a box of the other's.
fn faces_clear(a: &Face, b: &Face) -> bool {
    !face_cells(a).overlaps(face_cells(b))
}

/// The boxes of the cells of a grid over `face`'s range, each grown by the largest cell's
/// size (the surface between a cell's corners stays within that of them) and the face's
/// tolerance, indexed (made on first use).
fn face_cells(face: &Face) -> &CellPyramid {
    use crate::math::BoundingBox;
    face.cells.get_or_init(|| {
        const N: usize = 32;
        let ((u0, u1), (v0, v1)) = (face.u_range, face.v_range);
        let us: Vec<f64> = (0..=N).map(|i| u0 + (u1 - u0) * i as f64 / N as f64).collect();
        let vs: Vec<f64> = (0..=N).map(|j| v0 + (v1 - v0) * j as f64 / N as f64).collect();
        let rows = face.surface.points_on_grid(&us, &vs);
        let mut boxes = Vec::with_capacity(N * N);
        for j in 0..N {
            for i in 0..N {
                let mut b = BoundingBox::EMPTY;
                for p in [rows[j][i], rows[j][i + 1], rows[j + 1][i], rows[j + 1][i + 1]] {
                    b.add_point(p);
                }
                boxes.push(b);
            }
        }
        // (Each grown by the largest cell's size: the surface between the corners stays
        // within that of them.)
        let spacing = boxes.iter().map(|b| b.min.distance(b.max)).fold(0.0, f64::max);
        let margin = spacing + face.drift + on_tol();
        let boxes: Vec<BoundingBox> = boxes.into_iter().map(|b| b.enlarged(margin)).collect();
        CellPyramid::new(N, boxes)
    })
}

/// The boxes of a square grid's cells, row by row, and above them level by level those of
/// blocks of two by two of the level below (fewer at an odd width's end): whether one of
/// the cells meets a box, or one of another grid's, looked for only within the blocks that
/// do. (The same answers as any tree over the cells' boxes.)
pub(super) struct CellPyramid {
    /// Each level's width and boxes, the cells' first, the whole grid's last.
    levels: Vec<(usize, Vec<crate::math::BoundingBox>)>,
}

impl CellPyramid {
    fn new(width: usize, cells: Vec<crate::math::BoundingBox>) -> CellPyramid {
        use crate::math::BoundingBox;
        let mut levels = vec![(width, cells)];
        while levels.last().unwrap().0 > 1 {
            let (w, below) = levels.last().unwrap();
            let (w, up) = (*w, (*w + 1) / 2);
            let mut boxes = vec![BoundingBox::EMPTY; up * up];
            for j in 0..w {
                for i in 0..w {
                    boxes[(j / 2) * up + i / 2].add_box(&below[j * w + i]);
                }
            }
            levels.push((up, boxes));
        }
        CellPyramid { levels }
    }

    /// The blocks of `level` - 1 within block `(i, j)` of `level`.
    fn below(&self, level: usize, (i, j): (usize, usize)) -> impl Iterator<Item = (usize, usize)> {
        let w = self.levels[level - 1].0;
        (2 * j..(2 * j + 2).min(w)).flat_map(move |y| (2 * i..(2 * i + 2).min(w)).map(move |x| (x, y)))
    }

    fn at(&self, level: usize, (i, j): (usize, usize)) -> &crate::math::BoundingBox {
        let (w, boxes) = &self.levels[level];
        &boxes[j * w + i]
    }

    /// Whether a cell's box meets `query`.
    fn meets(&self, query: &crate::math::BoundingBox) -> bool {
        let mut stack = vec![(self.levels.len() - 1, (0, 0))];
        while let Some((level, block)) = stack.pop() {
            if !boxes_meet(self.at(level, block), query) {
                continue;
            }
            if level == 0 {
                return true;
            }
            stack.extend(self.below(level, block).map(|b| (level - 1, b)));
        }
        false
    }

    /// Whether a cell's box meets one of `other`'s.
    fn overlaps(&self, other: &CellPyramid) -> bool {
        let mut stack = vec![((self.levels.len() - 1, (0, 0)), (other.levels.len() - 1, (0, 0)))];
        while let Some(((la, a), (lb, b))) = stack.pop() {
            if !boxes_meet(self.at(la, a), other.at(lb, b)) {
                continue;
            }
            match (la, lb) {
                (0, 0) => return true,
                // (Down the higher of the two.)
                _ if la >= lb => stack.extend(self.below(la, a).map(|a| ((la - 1, a), (lb, b)))),
                _ => stack.extend(other.below(lb, b).map(|b| ((la, a), (lb - 1, b)))),
            }
        }
        false
    }
}

/// Whether two boxes meet (touching counts).
fn boxes_meet(a: &crate::math::BoundingBox, b: &crate::math::BoundingBox) -> bool {
    a.min.x <= b.max.x && b.min.x <= a.max.x && a.min.y <= b.max.y && b.min.y <= a.max.y && a.min.z <= b.max.z && b.min.z <= a.max.z
}

/// The signed distances found lately from points of curves to free-form surfaces (by the
/// curve's parameter's bits), with the surface, the curve and its range: a boolean operation
/// run again at the looser tolerances asks for the same ones (shared by the threads of
/// operations).
static CONTACT_DISTANCES: std::sync::Mutex<Vec<(ContactKey, std::sync::Arc<HashMap<u64, f64>>)>> = std::sync::Mutex::new(Vec::new());

/// What [`CONTACT_DISTANCES`] keeps distances for: a surface by its identity, the curve whole.
struct ContactKey {
    surface: crate::geometry::extrema::SurfaceKey,
    curve: Curve,
    range: (u64, u64),
}

fn surface_contact(curve: &Curve, first: f64, last: f64, surface: &Surface, tolerance: f64) -> Contact {
    if let Some(contact) = plane_contact(curve, first, last, surface).or_else(|| coaxial_contact(curve, first, last, surface)) {
        return contact;
    }
    // (As many as the surface's sampling has, at least: one winding round many times may be
    // crossed as often.)
    let n = if matches!(curve, Curve::Line(_)) { 32 } else { 128 }.max(geom::winding_samples(surface));
    let ts: Vec<f64> = (0..=n).map(|i| first + (last - first) * i as f64 / n as f64).collect();
    // (On a free-form surface, known by its identity, the distances found before for the same
    // curve kept: they depend on nothing else, the tolerances only on what is made of them.)
    let lock = || CONTACT_DISTANCES.lock().unwrap_or_else(|poisoned| poisoned.into_inner());
    let key = surface.shared_id().map(|_| ContactKey { surface: crate::geometry::extrema::SurfaceKey::of(surface), curve: curve.clone(), range: (first.to_bits(), last.to_bits()) });
    let same = |k: &ContactKey, key: &ContactKey| k.surface == key.surface && k.range == key.range && k.curve == key.curve;
    let known = key.as_ref().and_then(|key| lock().iter().find(|(k, _)| same(k, key)).map(|(_, d)| d.clone()));
    let found: std::cell::RefCell<HashMap<u64, f64>> = Default::default();
    let g = |t: f64| {
        if let Some(&d) = known.as_ref().and_then(|k| k.get(&t.to_bits())) {
            return d;
        }
        let d = geom::signed_distance(surface, curve.point(t));
        if key.is_some() {
            found.borrow_mut().insert(t.to_bits(), d);
        }
        d
    };
    let contact = contact_from(curve, first, last, surface, tolerance, n, &ts, &g);
    if let Some(key) = key {
        let found = found.into_inner();
        if !found.is_empty() {
            let mut all = known.map(|k| (*k).clone()).unwrap_or_default();
            all.extend(found);
            let mut cache = lock();
            cache.retain(|(k, _)| !same(k, &key));
            if cache.len() >= 1024 {
                cache.remove(0);
            }
            cache.push((key, std::sync::Arc::new(all)));
        }
    }
    contact
}

/// [`surface_contact`] from the distances `g` gives.
#[allow(clippy::too_many_arguments)]
fn contact_from(curve: &Curve, first: f64, last: f64, surface: &Surface, tolerance: f64, n: usize, ts: &[f64], g: &dyn Fn(f64) -> f64) -> Contact {
    let values: Vec<f64> = ts.iter().map(|&t| g(t)).collect();
    // Lying on the surface within the tolerance of meeting points (shapes made by earlier
    // operations put edges on surfaces only that closely).
    let within = |bound: f64| values.iter().all(|v| v.abs() <= bound) && ts.windows(2).all(|w| g(0.5 * (w[0] + w[1])).abs() <= bound);
    // (Or within its own tolerance all along, crossing the surface to and fro: an edge made
    // only that close to a surface it lies on, as a fillet's along the face it rolls on.
    // Not one merely near it, as an edge of a neighbouring face may be.)
    let weaving = values.windows(2).filter(|w| w[0] * w[1] < 0.0).count() >= 4;
    let on = within(on_tol()) || (tolerance > on_tol() && weaving && within(tolerance));
    if on {
        return Contact::On;
    }
    let mut roots = Vec::new();
    for i in 0..=n {
        if values[i].abs() <= TOL {
            roots.push(ts[i]);
        }
        // (Where a sample is itself a root, the interval beside it may hold another: the
        // sign just inside from it tells, as a line through a cylinder whose one crossing
        // falls on a sample.)
        let bracket = if i == n {
            None
        } else if values[i].abs() > TOL && values[i + 1].abs() > TOL {
            (values[i] * values[i + 1] < 0.0).then_some((ts[i], ts[i + 1], values[i]))
        } else if (values[i].abs() > TOL) != (values[i + 1].abs() > TOL) {
            let nudge = 1e-3 * (ts[i + 1] - ts[i]);
            let (outer, root) = if values[i].abs() > TOL { (i, i + 1) } else { (i + 1, i) };
            let inside = ts[root] + if root > outer { -nudge } else { nudge };
            let g_inside = g(inside);
            (g_inside.abs() > TOL && g_inside * values[outer] < 0.0).then(|| if outer < root { (ts[outer], inside, values[outer]) } else { (inside, ts[outer], g_inside) })
        } else {
            None
        };
        if let Some((lo, hi, glo)) = bracket {
            let (mut lo, mut hi, mut glo) = (lo, hi, glo);
            for _ in 0..100 {
                let mid = 0.5 * (lo + hi);
                if mid == lo || mid == hi {
                    break;
                }
                let gm = g(mid);
                if gm == 0.0 {
                    lo = mid;
                    hi = mid;
                    break;
                }
                if (gm < 0.0) == (glo < 0.0) {
                    lo = mid;
                    glo = gm;
                } else {
                    hi = mid;
                }
            }
            let t = 0.5 * (lo + hi);
            // A jump of the signed distance (on surfaces without a closed form) is no root.
            if g(t).abs() <= on_tol() {
                roots.push(t);
            }
        }
        // A touching point: |g| has a small local minimum without a sign change. It is
        // where g's derivative along the curve vanishes, found precisely by bisection
        // (minimizing |g| itself is imprecise there, g being quadratic).
        if i > 0 && i < n {
            let (a, b, c) = (values[i - 1], values[i], values[i + 1]);
            // (A strict extremum: where g is flat, as along a circle about a cylinder's axis,
            // there is none to find.)
            let flat = 1e-12 * (1.0 + a.abs().max(c.abs()));
            if b.abs() > TOL && b.abs() <= a.abs() && b.abs() <= c.abs() && (a.abs() - b.abs() > flat || c.abs() - b.abs() > flat) && a * b > 0.0 && b * c > 0.0 {
                let slope = |t: f64| {
                    let (p, dp) = curve.d1(t);
                    // (On a free-form surface the gradient comes with the distance, from one
                    // foot: the unit normal there, the way the distance grows.)
                    if let (_, Some(n)) = geom::signed_distance_and_normal(surface, p) {
                        return n.dot(dp);
                    }
                    let h = 1e-7 * (1.0 + (p - Point::ORIGIN).length());
                    let n = geom::normal_at_point(surface, p);
                    let n = if geom::signed_distance(surface, p + n * h) < geom::signed_distance(surface, p - n * h) { -n } else { n };
                    n.dot(dp)
                };
                let (mut lo, mut hi) = (ts[i - 1], ts[i + 1]);
                let (slo, shi) = (slope(lo), slope(hi));
                let t = if slo * shi < 0.0 {
                    for _ in 0..100 {
                        let mid = 0.5 * (lo + hi);
                        if mid == lo || mid == hi {
                            break;
                        }
                        if (slope(mid) < 0.0) == (slo < 0.0) {
                            lo = mid;
                        } else {
                            hi = mid;
                        }
                    }
                    0.5 * (lo + hi)
                } else {
                    ts[i]
                };
                let gt = g(t);
                if gt.abs() <= TOL {
                    roots.push(t);
                } else if gt * b < 0.0 {
                    // The curve dips through the surface between samples (a short stretch
                    // across a cylinder's edge): a root on either side of the extremum.
                    for (lo, hi) in [(ts[i - 1], t), (t, ts[i + 1])] {
                        let (mut lo, mut hi, glo) = (lo, hi, g(lo));
                        for _ in 0..100 {
                            let mid = 0.5 * (lo + hi);
                            if mid == lo || mid == hi {
                                break;
                            }
                            if (g(mid) < 0.0) == (glo < 0.0) {
                                lo = mid;
                            } else {
                                hi = mid;
                            }
                        }
                        roots.push(0.5 * (lo + hi));
                    }
                }
            }
        }
    }
    roots.sort_by(f64::total_cmp);
    let resolution = 1e-9 * (last - first).abs().max(1e-12);
    roots.dedup_by(|b, a| (*b - *a).abs() <= resolution);
    Contact::At(roots)
}

/// [`surface_contact`] for a circle about a cylinder's or cone's axis (or a sphere's
/// centre): at a constant distance from it, on it or nowhere.
fn coaxial_contact(curve: &Curve, first: f64, last: f64, surface: &Surface) -> Option<Contact> {
    let Curve::Circle(circle) = curve else { return None };
    let f = circle.frame();
    let normal = Vector::from(f.x()).cross(Vector::from(f.y()));
    let around = |origin: Point, axis: Vector| {
        let d = circle.center() - origin;
        normal.cross(axis).length() <= 1e-12 && (d - axis * d.dot(axis)).length() <= 1e-12 * (1.0 + d.length())
    };
    let coaxial = match surface {
        Surface::Cylinder(c) => around(c.frame().origin(), Vector::from(c.frame().z())),
        Surface::Cone(c) => around(c.frame().origin(), Vector::from(c.frame().z())),
        Surface::Sphere(s) => circle.center().distance(s.frame().origin()) <= 1e-12 * (1.0 + s.radius()) || {
            let d = circle.center() - s.frame().origin();
            normal.cross(d).length() <= 1e-12 * (1.0 + d.length())
        },
        _ => false,
    };
    if !coaxial {
        return None;
    }
    let g = geom::signed_distance(surface, curve.point(first));
    Some(if g.abs() <= on_tol() && geom::signed_distance(surface, curve.point(0.5 * (first + last))).abs() <= on_tol() { Contact::On } else { Contact::At(Vec::new()) })
}

/// [`surface_contact`] in closed form for a line or a circle against a plane.
fn plane_contact(curve: &Curve, first: f64, last: f64, surface: &Surface) -> Option<Contact> {
    let Surface::Plane(plane) = surface else { return None };
    let g = |t: f64| geom::signed_distance(surface, curve.point(t));
    let normal = geom::plane_normal(plane);
    let (lo, hi) = (first.min(last), first.max(last));
    let mut roots: Vec<f64> = [first, last].into_iter().filter(|&t| g(t).abs() <= TOL).collect();
    let within = |t: f64| t >= lo - 1e-12 * (1.0 + t.abs()) && t <= hi + 1e-12 * (1.0 + t.abs());
    match curve {
        Curve::Line(line) => {
            // Within the tolerance at both ends: all along.
            if g(lo).abs() <= on_tol() && g(hi).abs() <= on_tol() {
                return Some(Contact::On);
            }
            let slope = normal.dot(Vector::from(line.direction));
            if slope.abs() <= 1e-12 {
                return Some(Contact::At(roots));
            }
            let t = first - g(first) / slope;
            if within(t) {
                roots.push(t.clamp(lo, hi));
            }
        }
        Curve::Circle(circle) => {
            let f = circle.frame();
            let (a, b) = (normal.dot(Vector::from(f.x())) * circle.radius(), normal.dot(Vector::from(f.y())) * circle.radius());
            let c = (circle.center() - plane.origin()).dot(normal);
            let amplitude = a.hypot(b);
            if amplitude <= 1e-12 * circle.radius().max(1.0) {
                // In a plane parallel to the surface.
                return Some(if c.abs() <= on_tol() { Contact::On } else { Contact::At(roots) });
            }
            // g(t) = c + amplitude cos(t - phase).
            let phase = b.atan2(a);
            // Within the tolerance at its ends and extremes in the range: all along.
            let extremes = [phase, phase + std::f64::consts::PI].map(|t| lo + (t - lo).rem_euclid(std::f64::consts::TAU));
            if [lo, hi].into_iter().chain(extremes.into_iter().filter(|&t| t <= hi)).all(|t| g(t).abs() <= on_tol()) {
                return Some(Contact::On);
            }
            let candidates: Vec<f64> = if c.abs() >= amplitude - TOL {
                // Touching (or missing) at the nearest point.
                let t = if c > 0.0 { phase + std::f64::consts::PI } else { phase };
                if (c.abs() - amplitude).abs() <= TOL { vec![t] } else { Vec::new() }
            } else {
                let spread = (-c / amplitude).clamp(-1.0, 1.0).acos();
                vec![phase + spread, phase - spread]
            };
            let period = std::f64::consts::TAU;
            for t in candidates {
                // Into the edge's range.
                let t = lo + (t - lo).rem_euclid(period);
                if within(t) {
                    roots.push(t.clamp(lo, hi));
                } else if within(t - period) {
                    roots.push((t - period).clamp(lo, hi));
                }
            }
        }
        _ => return None,
    }
    roots.sort_by(f64::total_cmp);
    let resolution = 1e-9 * (last - first).abs().max(1e-12);
    roots.dedup_by(|b, a| (*b - *a).abs() <= resolution);
    Some(Contact::At(roots))
}

/// Of the parameters where a curve meets a surface, those not inside a stretch lying on it:
/// a stretch's ends, found precisely.
fn stretch_ends(curve: &Curve, first: f64, last: f64, surface: &Surface, ts: Vec<f64>) -> Vec<f64> {
    let on = |t: f64| geom::signed_distance(surface, curve.point(t)).abs() <= TOL;
    let lying = |a: f64, b: f64| on(0.5 * (a + b)) && on(a + 0.25 * (b - a)) && on(a + 0.75 * (b - a));
    // The point between `inside` (on the surface) and `outside` where the curve leaves it.
    let leave = |mut inside: f64, mut outside: f64| {
        for _ in 0..60 {
            let mid = 0.5 * (inside + outside);
            if on(mid) {
                inside = mid;
            } else {
                outside = mid;
            }
        }
        inside
    };
    let step = (last - first) / 128.0;
    // A stretch where the curve only grazes the surface within the tolerance (touching it,
    // nearly tangent) is the touching point: on a stretch it lies on, the curve stays far
    // closer than the tolerance.
    let grazes = |a: f64, b: f64| {
        [0.25, 0.75].iter().any(|&s| geom::signed_distance(surface, curve.point(a + s * (b - a))).abs() > 0.01 * TOL)
    };
    let mut out = Vec::new();
    let mut run_start: Option<f64> = None;
    for i in 0..ts.len() {
        let before = i > 0 && lying(ts[i - 1], ts[i]);
        let after = i + 1 < ts.len() && lying(ts[i], ts[i + 1]);
        match (before, after) {
            (true, true) => {}
            (false, true) => run_start = Some(if ts[i] > first { leave(ts[i], (ts[i] - step).max(first)) } else { ts[i] }),
            (true, false) => {
                let end = if ts[i] < last { leave(ts[i], (ts[i] + step).min(last)) } else { ts[i] };
                let start = run_start.take().unwrap_or(end);
                if grazes(start, end) {
                    // The point of the stretch closest to the surface.
                    let closest = (0..=32)
                        .map(|k| start + (end - start) * k as f64 / 32.0)
                        .min_by(|&x, &y| {
                            let d = |t: f64| geom::signed_distance(surface, curve.point(t)).abs();
                            d(x).total_cmp(&d(y))
                        })
                        .unwrap_or(start);
                    out.push(closest);
                } else {
                    out.push(start);
                    out.push(end);
                }
            }
            _ => out.push(ts[i]),
        }
    }
    out
}

/// Where a plane across the direction of a surface of extrusion meets it: the extruded
/// curve projected along the direction onto the plane (exactly, for curves defined by poles
/// or planar across the direction).
fn plane_across_extrusion(plane: &Surface, extrusion: &Surface) -> Option<Section> {
    use crate::math::{AffineTransform, AffineTransformable, Transformable};
    let (Surface::Plane(plane), Surface::Extrusion(extrusion)) = (plane, extrusion) else { return None };
    let (n, d) = (Vector::from(plane.normal()), Vector::from(extrusion.direction()));
    if n.cross(d).length() > crate::tolerance::ANGULAR.sqrt() {
        return None;
    }
    let basis = extrusion.basis();
    let (first, last) = (basis.first_parameter(), basis.last_parameter());
    let projection = AffineTransform::plane_affinity(plane.origin(), plane.normal(), 0.0);
    let curve = match basis {
        Curve::BSpline(c) => Curve::BSpline(c.affine_transformed(&projection)),
        Curve::Bezier(c) => Curve::Bezier(c.affine_transformed(&projection)),
        c => {
            // Planar across the direction: moved onto the plane.
            let height = |t: f64| n.dot(c.point(t) - plane.origin());
            let h = height(first);
            if (0..=16).any(|i| (height(first + (last - first) * i as f64 / 16.0) - h).abs() > TOL) {
                return None;
            }
            c.transformed(&crate::math::Transform::translation(-(n * h)))
        }
    };
    Some(Section::Open(curve, first, last))
}

/// Where a plane along the direction of a surface of extrusion meets it: lines along the
/// direction through the points of the extruded curve on the plane, the roots of its
/// height above the plane, found between samples and either side of the turns of the
/// height (two roots close together, the curve dipping just through the plane, or one
/// where it touches the plane). `None` unless the plane is along the direction.
fn plane_along_extrusion(fa: &Face, fb: &Face) -> Option<Vec<Section>> {
    let (plane, extrusion, extruded) = match (&fa.surface, &fb.surface) {
        (Surface::Plane(p), Surface::Extrusion(e)) => (p, e, fb),
        (Surface::Extrusion(e), Surface::Plane(p)) => (p, e, fa),
        _ => return None,
    };
    let (n, d) = (Vector::from(plane.normal()), Vector::from(extrusion.direction()));
    if n.dot(d).abs() > crate::tolerance::ANGULAR {
        return None;
    }
    let curve = extrusion.basis();
    let (b0, b1) = (curve.first_parameter(), curve.last_parameter());
    // (Over the extruded face's range, a little beyond it, within the curve's.)
    let (u0, u1) = extruded.u_range;
    let span = u1 - u0;
    let (t0, t1) = ((u0 - 0.05 * span).max(b0), (u1 + 0.05 * span).min(b1));
    if !(t1 > t0) {
        return Some(Vec::new());
    }
    let origin = plane.origin();
    let height = |t: f64| n.dot(curve.point(t) - origin);
    let slope = |t: f64| n.dot(curve.d1(t).1);
    const N: usize = 512;
    let ts: Vec<f64> = (0..=N).map(|i| t0 + (t1 - t0) * i as f64 / N as f64).collect();
    let hs: Vec<f64> = ts.iter().map(|&t| height(t)).collect();
    let ds: Vec<f64> = ts.iter().map(|&t| slope(t)).collect();
    // (The extruded curve in the plane: the surfaces are one there, met along no line.)
    if hs.iter().all(|h| h.abs() <= TOL) {
        return Some(Vec::new());
    }
    let root = |f: &dyn Fn(f64) -> f64, mut lo: f64, mut hi: f64| {
        let flo = f(lo);
        for _ in 0..100 {
            let mid = 0.5 * (lo + hi);
            if (f(mid) < 0.0) == (flo < 0.0) {
                lo = mid;
            } else {
                hi = mid;
            }
            if hi - lo <= 1e-15 * (t1 - t0) {
                break;
            }
        }
        0.5 * (lo + hi)
    };
    let mut roots = Vec::new();
    for i in 0..N {
        let (a, b) = (ts[i], ts[i + 1]);
        if hs[i] == 0.0 {
            roots.push(a);
            continue;
        }
        if hs[i] * hs[i + 1] < 0.0 {
            roots.push(root(&height, a, b));
        } else if ds[i] * ds[i + 1] < 0.0 && hs[i + 1] != 0.0 {
            // A turn of the height: through the plane and back, or touching it.
            let turn = root(&slope, a, b);
            let h = height(turn);
            if h * hs[i] < 0.0 {
                roots.push(root(&height, a, turn));
                roots.push(root(&height, turn, b));
            } else if h.abs() <= TOL {
                roots.push(turn);
            }
        }
    }
    if hs[N] == 0.0 {
        roots.push(ts[N]);
    }
    let boxes = [&fa.bbox, &fb.bbox];
    Some(
        roots
            .into_iter()
            .map(|t| {
                let line = Line::new(curve.point(t) - n * height(t), extrusion.direction());
                let along = |bx: &(Point, Point)| {
                    let (lo, hi) = *bx;
                    let ts = (0..8).map(|i| {
                        let corner = Point::new(if i & 1 == 0 { lo.x } else { hi.x }, if i & 2 == 0 { lo.y } else { hi.y }, if i & 4 == 0 { lo.z } else { hi.z });
                        (corner - line.origin).dot(d)
                    });
                    ts.fold((f64::INFINITY, f64::NEG_INFINITY), |(a, b), t| (a.min(t), b.max(t)))
                };
                let ((lo_a, hi_a), (lo_b, hi_b)) = (along(boxes[0]), along(boxes[1]));
                Section::Open(Curve::Line(line), lo_a.max(lo_b) - 1.0, hi_a.min(hi_b) + 1.0)
            })
            .collect(),
    )
}

/// A curve where two faces' surfaces meet.
enum Section {
    /// A line, restricted to a range covering both faces.
    Open(Curve, f64, f64),
    /// A circle or ellipse (period 2π), with its samples for projections onto it (found
    /// once).
    Closed(Curve, std::sync::OnceLock<Vec<Point>>),
    /// A chain of branches of a cylinder's or cone's intersection with another quadric,
    /// each run forwards or backwards; closed, or open where it leaves the faces' range.
    Chain(Vec<(QuadricCurve, bool)>, bool, ChainSamples),
    /// A curve traced on a parametrized surface (`(u, v)` samples) where another,
    /// elementary surface meets it, closed or not; points in between are interpolated in
    /// `(u, v)` and pulled onto the other surface.
    Sampled(Vec<Point>, bool, Box<Traced>),
}

/// A chain's points at even steps of its parameter, found once (slow to find, and looked for
/// near each point the section may pass): 512 a branch, as [`Section::passes`] steps, and
/// 256, as [`Section::project`] does.
#[derive(Default)]
struct ChainSamples {
    passes: std::sync::OnceLock<Vec<Point>>,
    project: std::sync::OnceLock<Vec<Point>>,
}

/// A traced curve: its samples in the base surface's parameters, and the two surfaces.
struct Traced {
    uv: Vec<(f64, f64)>,
    /// The periods of the base's `u` and `v`, if the curve may cross their seams.
    u_period: Option<f64>,
    v_period: Option<f64>,
    /// The other surface's parameters at the samples, where known (marched along both).
    st: Option<Vec<(f64, f64)>>,
    base: Surface,
    other: Surface,
    /// The point last found by [`Section::point`] and its parameter's bits (asked for again
    /// as often as not: the point a projection onto the curve ends at).
    last: std::sync::Mutex<Option<(u64, Point)>>,
}

impl Traced {
    /// The point of the curve near `(u, v)` on the base surface.
    fn point(&self, uv: (f64, f64)) -> Point {
        use crate::geometry::ParametricSurface;
        let (u, v) = self.uv_at(uv);
        self.base.point(u, v)
    }

    /// The base's parameters of the curve near `(u, v)`: Newton steps in `(u, v)` on the
    /// other surface's signed distance, along its gradient.
    fn uv_at(&self, (mut u, mut v): (f64, f64)) -> (f64, f64) {
        use crate::geometry::ParametricSurface;
        if geom::freeform(&self.other) {
            if let Some(found) = self.on_both((u, v)) {
                return found;
            }
            // The gradient in closed form: the other surface's normal at the foot point
            // along the base's derivatives (one inversion per step).
            let mut last = None;
            for _ in 0..12 {
                let (p, du, dv) = self.base.d1(u, v);
                let (value, foot) = geom::signed_distance_from(&self.other, p, last);
                last = foot;
                if value.abs() <= 1e-12 {
                    break;
                }
                // (The normal at the foot just found: the gradient of that distance.)
                let n = match &last {
                    Some(foot) => foot.normal(),
                    None => geom::normal_at_point(&self.other, p),
                };
                let (fu, fv) = (n.dot(du), n.dot(dv));
                let norm = fu * fu + fv * fv;
                if norm <= 1e-300 {
                    break;
                }
                u -= value * fu / norm;
                v -= value * fv / norm;
            }
            return (u, v);
        }
        if geom::closed_form_gradient(&self.other) {
            // Analytic: the gradient in closed form too.
            for _ in 0..12 {
                let (p, du, dv) = self.base.d1(u, v);
                let Some((value, gradient)) = geom::signed_distance_with_gradient(&self.other, p) else { break };
                if value.abs() <= 1e-12 {
                    break;
                }
                let (fu, fv) = (gradient.dot(du), gradient.dot(dv));
                let norm = fu * fu + fv * fv;
                if norm <= 1e-300 {
                    break;
                }
                u -= value * fu / norm;
                v -= value * fv / norm;
            }
            return (u, v);
        }
        let last_foot = std::cell::Cell::new(None);
        let f = |u: f64, v: f64| {
            let (d, foot) = geom::signed_distance_from(&self.other, self.base.point(u, v), last_foot.get());
            last_foot.set(foot);
            d
        };
        for _ in 0..12 {
            let value = f(u, v);
            if value.abs() <= 1e-12 {
                break;
            }
            let h = 1e-7;
            let (fu, fv) = ((f(u + h, v) - f(u - h, v)) / (2.0 * h), (f(u, v + h) - f(u, v - h)) / (2.0 * h));
            let norm = fu * fu + fv * fv;
            if norm <= 1e-300 {
                break;
            }
            u -= value * fu / norm;
            v -= value * fv / norm;
        }
        (u, v)
    }

    /// [`Self::uv_at`] on a free-form other surface by Newton steps on the point both
    /// surfaces share, `base(u, v) = other(s, t)` (each step moving `(u, v)` the least,
    /// `(s, t)` taking up the rest), from the other's foot of the start; `None` where they
    /// do not converge.
    fn on_both(&self, (u, v): (f64, f64)) -> Option<(f64, f64)> {
        use crate::geometry::ParametricSurface;
        let foot = geom::foot_on(&self.other, self.base.point(u, v));
        self.on_both_from((u, v), (foot.u, foot.v))
    }

    /// Whether the curve at parameter `t` (by its samples, as [`Section::point`]) was
    /// followed past either surface's bounds (where it was marched).
    fn beyond_bounds(&self, t: f64) -> bool {
        use crate::geometry::ParametricSurface;
        let Some(st) = &self.st else { return false };
        let n = st.len();
        let (i, f) = ((t.floor() as usize) % n, t - t.floor());
        let at = |values: &[(f64, f64)]| {
            let (c, d) = (values[i], values[(i + 1) % n]);
            (c.0 + (d.0 - c.0) * f, c.1 + (d.1 - c.1) * f)
        };
        let past = |x: f64, lo: f64, hi: f64, periodic: bool| {
            let hair = 1e-7 * (hi - lo).abs().min(1e100).max(1.0);
            !periodic && (x < lo - hair || x > hi + hair)
        };
        let off = |surface: &Surface, (x, y): (f64, f64)| {
            let (x0, x1, y0, y1) = surface.bounds();
            past(x, x0, x1, surface.is_u_periodic()) || past(y, y0, y1, surface.is_v_periodic())
        };
        let (uv, st) = (at(&self.uv), at(st));
        if !(off(&self.base, uv) || off(&self.other, st)) {
            return false;
        }
        // (Between samples past the same bound, the curve may come back inside: a section
        // leaving a face's edge and meeting it again, its samples, only at its ends, each a
        // hair past it. Where the point found on both surfaces from there is not.)
        match self.on_both_with(uv, st) {
            Some((uv, st)) => off(&self.base, uv) || off(&self.other, st),
            None => true,
        }
    }

    /// The other surface's point and derivatives at `(s, t)`; beyond its bounds carried on
    /// across its edge in its tangent plane (square to the edge): where its distance, from
    /// the nearest point of the edge along the normal there, is 0.
    fn other_extended(&self, s: f64, t: f64) -> (Point, Vector, Vector) {
        use crate::geometry::ParametricSurface;
        let (s0, s1, t0, t1) = self.other.bounds();
        let (sc, tc) = (if self.other.is_u_periodic() { s } else { s.clamp(s0, s1) }, if self.other.is_v_periodic() { t } else { t.clamp(t0, t1) });
        let (q, qs, qt) = self.other.d1(sc, tc);
        let (ds, dt) = (s - sc, t - tc);
        if ds == 0.0 && dt == 0.0 {
            return (q, qs, qt);
        }
        let across = |a: Vector, b: Vector| {
            let bb = b.dot(b);
            if bb > 0.0 { a - b * (a.dot(b) / bb) } else { a }
        };
        let (es, et) = (if ds != 0.0 { across(qs, qt) } else { qs }, if dt != 0.0 { across(qt, qs) } else { qt });
        (q + es * ds + et * dt, es, et)
    }

    /// The base's parameters of the curve's point in the plane through `origin` square to
    /// `dir`, by Newton steps on both surfaces from `(u, v)` and `(s, t)` (as the curve was
    /// marched); `None` where they do not converge.
    fn across_chord(&self, (u, v): (f64, f64), (s, t): (f64, f64), origin: Point, dir: Vector) -> Option<(f64, f64)> {
        use crate::geometry::ParametricSurface;
        let residual = |x: [f64; 4]| {
            let (p, pu, pv) = self.base.d1(x[0], x[1]);
            let (q, qs, qt) = self.other_extended(x[2], x[3]);
            (p - q, (p - origin).dot(dir), p, [pu, pv, qs, qt])
        };
        let mut x = [u, v, s, t];
        let (mut f, mut g, mut p, mut d) = residual(x);
        for _ in 0..30 {
            let scale = 1e-12 * (1.0 + (p - Point::ORIGIN).length());
            if f.length() <= scale && g.abs() <= scale {
                return Some((x[0], x[1]));
            }
            let [pu, pv, qs, qt] = d;
            let m = [
                [pu.x, pv.x, -qs.x, -qt.x],
                [pu.y, pv.y, -qs.y, -qt.y],
                [pu.z, pv.z, -qs.z, -qt.z],
                [pu.dot(dir), pv.dot(dir), 0.0, 0.0],
            ];
            let step = super::march::solve4(m, [-f.x, -f.y, -f.z, -g])?;
            // (The step shortened until it comes nearer: from a guess off the curve, where the
            // surfaces' parameters bend, a whole step may lead far away.)
            let size = (f.length_squared() + g * g).sqrt();
            let mut k = 1.0;
            loop {
                let y = [x[0] + k * step[0], x[1] + k * step[1], x[2] + k * step[2], x[3] + k * step[3]];
                if !y.iter().all(|v| v.is_finite()) {
                    return None;
                }
                let (f2, g2, p2, d2) = residual(y);
                if (f2.length_squared() + g2 * g2).sqrt() < size || k < 1e-3 {
                    (x, f, g, p, d) = (y, f2, g2, p2, d2);
                    break;
                }
                k *= 0.5;
            }
        }
        None
    }

    /// [`Self::on_both`] from `(s, t)` on the other surface.
    fn on_both_from(&self, uv: (f64, f64), st: (f64, f64)) -> Option<(f64, f64)> {
        self.on_both_with(uv, st).map(|(uv, _)| uv)
    }

    /// [`Self::on_both_from`], with the other surface's parameters there too.
    fn on_both_with(&self, (mut u, mut v): (f64, f64), (mut s, mut t): (f64, f64)) -> Option<((f64, f64), (f64, f64))> {
        use crate::geometry::ParametricSurface;
        let (bu0, bu1, bv0, bv1) = self.base.bounds();
        // (u, v) weighed against (s, t): steps of the other's parameters a million times
        // cheaper.
        const FREE: f64 = 1e6;
        for _ in 0..20 {
            let (p, pu, pv) = self.base.d1(u, v);
            let (q, qs, qt) = self.other_extended(s, t);
            let f = p - q;
            if f.length() <= 1e-12 * (1.0 + (p - Point::ORIGIN).length()) {
                return Some(((u, v), (s, t)));
            }
            // J = [pu pv -qs -qt], the step -W J^T (J W J^T)^-1 f with W = diag(1, 1, FREE, FREE).
            let m = |a: Vector, b: Vector| a.x * b.x + a.y * b.y + a.z * b.z;
            let col = [pu, pv, qs * -1.0, qt * -1.0];
            let w = [1.0, 1.0, FREE, FREE];
            let mut mm = [[0.0; 3]; 3];
            for k in 0..4 {
                let c = [col[k].x, col[k].y, col[k].z];
                for r in 0..3 {
                    for cc in 0..3 {
                        mm[r][cc] += w[k] * c[r] * c[cc];
                    }
                }
            }
            let det = mm[0][0] * (mm[1][1] * mm[2][2] - mm[1][2] * mm[2][1]) - mm[0][1] * (mm[1][0] * mm[2][2] - mm[1][2] * mm[2][0]) + mm[0][2] * (mm[1][0] * mm[2][1] - mm[1][1] * mm[2][0]);
            let scale = mm[0][0].abs().max(mm[1][1].abs()).max(mm[2][2].abs());
            if !(det.abs() > 1e-24 * scale * scale * scale) {
                return None;
            }
            // y = M^-1 f by Cramer's rule.
            let fv = [f.x, f.y, f.z];
            let solve = |k: usize| {
                let mut a = mm;
                for r in 0..3 {
                    a[r][k] = fv[r];
                }
                (a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1]) - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0]) + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0])) / det
            };
            let y = Vector::new(solve(0), solve(1), solve(2));
            let step: Vec<f64> = (0..4).map(|k| -w[k] * m(col[k], y)).collect();
            (u, v, s, t) = (u + step[0], v + step[1], s + step[2], t + step[3]);
            if !(u.is_finite() && v.is_finite() && s.is_finite() && t.is_finite()) {
                return None;
            }
            // (Not off the base, where it is extrapolated; but a little, where the section was
            // marched, which follows them a little past its bounds.)
            let past = if self.st.is_some() { 0.05 } else { 0.0 };
            let off = |x: f64, lo: f64, hi: f64, periodic: bool| !periodic && (x < lo - past * (hi - lo) || x > hi + past * (hi - lo));
            if off(u, bu0, bu1, self.base.is_u_periodic()) || off(v, bv0, bv1, self.base.is_v_periodic()) {
                return None;
            }
        }
        None
    }

    /// The parameters a fraction `s` of the way from `a` to `b` (the short way round a
    /// periodic base's seam).
    fn between(&self, a: (f64, f64), b: (f64, f64), s: f64) -> (f64, f64) {
        let across = |d: f64, period: Option<f64>| match period {
            Some(period) => (d + 0.5 * period).rem_euclid(period) - 0.5 * period,
            None => d,
        };
        (a.0 + across(b.0 - a.0, self.u_period) * s, a.1 + across(b.1 - a.1, self.v_period) * s)
    }

    /// The samples in `uv` with more between them where the curve strays from the chords
    /// between them by more than a hundredth of their length and than `floor` (the
    /// tracing's grid may be coarse for the curve).
    fn densified(&self, closed: bool, floor: f64) -> Vec<(f64, f64)> {
        use crate::geometry::ParametricSurface;
        let at = |(u, v): (f64, f64)| self.base.point(u, v);
        fn split(traced: &Traced, at: &dyn Fn((f64, f64)) -> Point, a: (f64, f64), b: (f64, f64), depth: usize, floor: f64, out: &mut Vec<(f64, f64)>) {
            let (pa, pb) = (at(a), at(b));
            let m = traced.uv_at(traced.between(a, b, 0.5));
            let pm = at(m);
            let chord = pa.distance(pb);
            if depth < 6 && pm.distance(pa + (pb - pa) * 0.5) > (1e-2 * chord).max(floor) + on_tol() && pm.distance(pa) < chord && pm.distance(pb) < chord {
                split(traced, at, a, m, depth + 1, floor, out);
                out.push(m);
                split(traced, at, m, b, depth + 1, floor, out);
            }
        }
        let n = self.uv.len();
        let segments = if closed { n } else { n - 1 };
        let mut out = Vec::with_capacity(n);
        for k in 0..segments {
            out.push(self.uv[k]);
            split(self, &at, self.uv[k], self.uv[(k + 1) % n], 0, floor, &mut out);
        }
        if !closed {
            out.push(self.uv[n - 1]);
        }
        out
    }
}


impl Section {
    fn range(&self) -> (f64, f64) {
        match self {
            Section::Open(_, a, b) => (*a, *b),
            Section::Closed(..) => (0.0, TAU),
            Section::Chain(branches, ..) => (0.0, branches.len() as f64),
            Section::Sampled(points, closed, _) => (0.0, if *closed { points.len() } else { points.len() - 1 } as f64),
        }
    }

    fn closed(&self) -> bool {
        match self {
            Section::Open(..) => false,
            Section::Closed(..) => true,
            Section::Chain(_, closed, _) => *closed,
            Section::Sampled(_, closed, _) => *closed,
        }
    }

    fn point(&self, t: f64) -> Point {
        match self {
            Section::Open(c, ..) | Section::Closed(c, _) => c.point(t),
            Section::Chain(branches, closed, _) => {
                let n = branches.len();
                let t = if *closed { t.rem_euclid(n as f64) } else { t.clamp(0.0, n as f64) };
                let k = (t.floor() as usize).min(n - 1);
                let (branch, reversed) = &branches[k];
                let s = t - k as f64;
                let s = if *reversed { 1.0 - s } else { s };
                let (a, b) = branch.domain();
                let full = b - a >= TAU - 1e-12;
                // Away from full circles the branch ends are square-root singular: spacing
                // the angle like a cosine keeps points evenly spread along the curve.
                let angle = if full { a + s * (b - a) } else { a + (b - a) * 0.5 * (1.0 - (PI * s).cos()) };
                branch.point_clamped(angle)
            }
            Section::Sampled(points, closed, traced) => {
                let last = || traced.last.lock().unwrap_or_else(|poisoned| poisoned.into_inner());
                if let Some((bits, p)) = *last() {
                    if bits == t.to_bits() {
                        return p;
                    }
                }
                let n = points.len();
                let (_, len) = self.range();
                let t = if *closed { t.rem_euclid(len) } else { t.clamp(0.0, len) };
                let i = (t.floor() as usize).min(if *closed { n - 1 } else { n - 2 });
                let s = t - i as f64;
                let (a, b) = (traced.uv[i], traced.uv[(i + 1) % n]);
                // Across the seam of a periodic base, interpolate the short way round.
                let across = |d: f64, period: Option<f64>| match period {
                    Some(period) => (d + 0.5 * period).rem_euclid(period) - 0.5 * period,
                    None => d,
                };
                let (du, dv) = (across(b.0 - a.0, traced.u_period), across(b.1 - a.1, traced.v_period));
                let uv = (a.0 + du * s, a.1 + dv * s);
                // (Where the curve was marched along both surfaces, as it was marched: on both,
                // in the plane square to the chord between the samples through the point
                // interpolated, from the parameters on both interpolated likewise. Elsewhere,
                // or where that does not converge, from the other's foot point.)
                let across_chord = || {
                    use crate::geometry::ParametricSurface;
                    let st = traced.st.as_ref()?;
                    let (c, d) = (st[i], st[(i + 1) % n]);
                    let (ds, dt) = (across(d.0 - c.0, traced.other.u_period()), across(d.1 - c.1, traced.other.v_period()));
                    let (pa, pb) = (points[i], points[(i + 1) % n]);
                    let chord = pb - pa;
                    if !(chord.length() > 0.0) {
                        return None;
                    }
                    let dir = chord * (1.0 / chord.length());
                    // (Not on another branch crossing the plane.)
                    let (u, v) = traced.across_chord(uv, (c.0 + ds * s, c.1 + dt * s), pa + chord * s, dir)?;
                    // (The curve strays from its chords by far less: the samples are spaced
                    // for that.)
                    (traced.base.point(u, v).distance(pa + chord * s) <= 0.25 * chord.length()).then_some((u, v))
                };
                let found = if traced.st.is_some() { across_chord().or_else(|| traced.on_both(uv)) } else { None };
                let p = match found {
                    Some((u, v)) => {
                        use crate::geometry::ParametricSurface;
                        traced.base.point(u, v)
                    }
                    None => traced.point(uv),
                };
                *last() = Some((t.to_bits(), p));
                p
            }
        }
    }

    /// The parameter of the point closest to `p`, and the distance.
    fn project(&self, p: Point) -> (f64, f64) {
        match self {
            Section::Open(c, a, b) => geom::project(c, *a, *b, p),
            // (A circle's projection looks at few samples; an ellipse's at all.)
            Section::Closed(c @ Curve::Circle(_), _) => geom::project(c, 0.0, TAU, p),
            Section::Closed(c, samples) => geom::project_among(c, 0.0, TAU, Some(samples.get_or_init(|| geom::projection_samples(c, 0.0, TAU))), p),
            Section::Sampled(points, closed, _) => {
                // The nearest stretch of the polyline, then the curve near it.
                let n = points.len();
                let segments = if *closed { n } else { n - 1 };
                let mut best = (0.0, f64::INFINITY);
                for i in 0..segments {
                    let (d, s) = segment_gap(points[i], points[(i + 1) % n], p);
                    if d < best.1 {
                        best = (i as f64 + s, d);
                    }
                }
                self.project_near(p, best.0)
            }
            Section::Chain(_, _, samples) => {
                let len = self.range().1;
                let n = (len as usize) * 256;
                let points = samples.project.get_or_init(|| (0..n).map(|i| self.point(len * i as f64 / n as f64)).collect());
                let mut best = (0.0, f64::INFINITY);
                for (i, q) in points.iter().enumerate() {
                    let t = len * i as f64 / n as f64;
                    let d = q.distance(p);
                    if d < best.1 {
                        best = (t, d);
                    }
                }
                let step = len / n as f64;
                let (mut lo, mut hi) = (best.0 - step, best.0 + step);
                let ratio = 0.618_033_988_749_895;
                for _ in 0..80 {
                    let (x1, x2) = (hi - ratio * (hi - lo), lo + ratio * (hi - lo));
                    if self.point(x1).distance(p) < self.point(x2).distance(p) {
                        hi = x2;
                    } else {
                        lo = x1;
                    }
                }
                let t = 0.5 * (lo + hi);
                let t = if self.closed() { t.rem_euclid(len) } else { t.clamp(0.0, len) };
                let d = self.point(t).distance(p);
                if d < best.1 { (t, d) } else { best }
            }
        }
    }

    /// Every parameter where the curve passes within [`on_tol`] of `p` (a loop through a
    /// singular point passes it more than once), each with whether it was found only near a
    /// point where the surfaces touch crossing (see [`touching_across`]).
    fn passes(&self, p: Point) -> Vec<(f64, bool)> {
        self.passes_within(p, on_tol())
    }

    /// [`Self::passes`], `p` within `tolerance` of the curve (on both surfaces as nearly).
    fn passes_within(&self, p: Point, tolerance: f64) -> Vec<(f64, bool)> {
        if let Section::Sampled(points, closed, traced) = self {
            // On the curve only if on both surfaces (one inversion each, before refining
            // along the curve).
            if geom::signed_distance(&traced.base, p).abs() > tolerance || geom::signed_distance(&traced.other, p).abs() > tolerance {
                return Vec::new();
            }
            // Each stretch of the polyline passing near p (a local minimum of the gap), refined
            // once.
            let n = points.len();
            let segments = if *closed { n } else { n - 1 };
            let gaps: Vec<(f64, f64, f64)> = (0..segments)
                .map(|i| {
                    let (a, b) = (points[i], points[(i + 1) % n]);
                    let (d, s) = segment_gap(a, b, p);
                    (d, s, a.distance(b))
                })
                .collect();
            let gap_at = |i: isize| -> f64 {
                if *closed {
                    gaps[i.rem_euclid(segments as isize) as usize].0
                } else if i < 0 || i >= segments as isize {
                    f64::INFINITY
                } else {
                    gaps[i as usize].0
                }
            };
            let mut out: Vec<(f64, bool)> = Vec::new();
            for i in 0..segments {
                let (d, s, length) = gaps[i];
                if d > 2.0 * length + tolerance || d > gap_at(i as isize - 1) || d >= gap_at(i as isize + 1) {
                    continue;
                }
                let (t, dist) = self.project_near(p, i as f64 + s);
                // Where the surfaces touch at a point and cross there (as two tubes turning a
                // corner do, above and below it) the section's branches cross in it and
                // refining onto both is singular: a point on both surfaces within the samples'
                // spacing (two samples', where the tracing's cells are coarse) of the polyline
                // may be on it there.
                let found = if dist <= tolerance {
                    Some((t, false))
                } else {
                    (d <= 2.0 * length + tolerance && touching_across(&traced.base, &traced.other, p)).then_some((i as f64 + s, true))
                };
                if let Some(found) = found.filter(|&(t, _)| !out.iter().any(|&(u, _)| (u - t).abs() <= 1.5)) {
                    out.push(found);
                }
            }
            return out;
        }
        let closed = &self.closed();
        let (n, cached) = match self {
            Section::Chain(branches, _, samples) => (branches.len() * 512, samples),
            _ => {
                let (t, d) = self.project(p);
                return if d <= tolerance { vec![(t, false)] } else { Vec::new() };
            }
        };
        let len = self.range().1;
        let step = len / n as f64;
        let samples = cached.passes.get_or_init(|| (0..n).map(|i| self.point(step * i as f64)).collect());
        let d: Vec<f64> = samples.iter().map(|q| q.distance(p)).collect();
        let mut out: Vec<f64> = Vec::new();
        for i in 0..n {
            if !closed && (i == 0 || i == n - 1) {
                // The ends of an open chain: check them directly below.
                continue;
            }
            let (a, b, c) = (d[(i + n - 1) % n], d[i], d[(i + 1) % n]);
            // A local minimum no farther than the neighbouring samples are apart (the
            // parameter is not proportional to arc length, so compare in space).
            let spacing = samples[(i + n - 1) % n].distance(samples[(i + 1) % n]);
            if b > a || b > c || b > spacing + tolerance {
                continue;
            }
            let (mut lo, mut hi) = (step * (i as f64 - 1.0), step * (i as f64 + 1.0));
            let ratio = 0.618_033_988_749_895;
            for _ in 0..80 {
                let (x1, x2) = (hi - ratio * (hi - lo), lo + ratio * (hi - lo));
                if self.point(x1).distance(p) < self.point(x2).distance(p) {
                    hi = x2;
                } else {
                    lo = x1;
                }
            }
            let t = 0.5 * (lo + hi);
            let t = if *closed { t.rem_euclid(len) } else { t.clamp(0.0, len) };
            if self.point(t).distance(p) <= tolerance && !out.iter().any(|&u| (u - t).abs() <= 2.0 * step) {
                out.push(t);
            }
        }
        if !closed {
            for t in [0.0, len] {
                if self.point(t).distance(p) <= tolerance && !out.iter().any(|&u| (u - t).abs() <= 2.0 * step) {
                    out.push(t);
                }
            }
        }
        // The nearest point overall, where the sampled minima miss it (at a junction of
        // branches, where the parametrization kinks).
        let (t, dist) = self.project(p);
        if dist <= on_tol() && !out.iter().any(|&u| (u - t).abs() <= 2.0 * step) {
            out.push(t);
        }
        out.into_iter().map(|t| (t, false)).collect()
    }

    /// The parameter of the point of a line or circle section closest to `p` by Newton steps
    /// from `seed` (others: the closest overall), and the distance.
    fn project_seeded(&self, p: Point, seed: f64) -> (f64, f64) {
        let (curve, lo, hi, periodic) = match self {
            Section::Open(c, a, b) => (c, a.min(*b), a.max(*b), false),
            Section::Closed(c, _) => (c, 0.0, TAU, true),
            _ => return self.project(p),
        };
        let mut t = seed;
        for _ in 0..12 {
            let (c, d1, d2) = curve.d2(t);
            let f = (c - p).dot(d1);
            let df = d1.dot(d1) + (c - p).dot(d2);
            if df.abs() <= 1e-300 {
                break;
            }
            let next = t - f / df;
            let next = if periodic { next } else { next.clamp(lo, hi) };
            let done = (next - t).abs() <= 1e-15 * (1.0 + t.abs());
            t = next;
            if done {
                break;
            }
        }
        if periodic {
            t = lo + (t - lo).rem_euclid(hi - lo);
        }
        (t, curve.point(t).distance(p))
    }

    /// The parameter of the point closest to `p` within one unit of `t`, and the distance.
    fn project_near(&self, p: Point, t: f64) -> (f64, f64) {
        let len = self.range().1;
        let (lo, hi) = if self.closed() { (t - 1.0, t + 1.0) } else { ((t - 1.0).max(0.0), (t + 1.0).min(len)) };
        // Newton along the curve: move by the gap's component along the tangent (no further
        // where the curve's points cannot be found).
        if t.is_nan() {
            return (lo, f64::INFINITY);
        }
        let mut t = t.clamp(lo, hi);
        let h = 1e-4;
        for _ in 0..12 {
            let q = self.point(t);
            let tangent = (self.point((t + h).min(hi)) - self.point((t - h).max(lo))) * (1.0 / (((t + h).min(hi) - (t - h).max(lo)).max(1e-300)));
            let step = (p - q).dot(tangent) / tangent.length_squared().max(1e-300);
            if !step.is_finite() {
                break;
            }
            let next = (t + step).clamp(lo, hi);
            if (next - t).abs() <= 1e-12 {
                t = next;
                break;
            }
            t = next;
        }
        let t = if self.closed() { t.rem_euclid(len) } else { t };
        (t, self.point(t).distance(p))
    }

    /// The edge geometry from `t0` to `t1` (wrapped past `t0` for closed curves), between
    /// the given end points.
    /// For a marched section, the parameter a fraction of the way along it from `t0` to
    /// `t1` (by the lengths between its samples).
    fn even_parameters(&self, t0: f64, t1: f64) -> Option<impl Fn(f64) -> f64> {
        let Section::Sampled(points, _, traced) = self else { return None };
        traced.st.as_ref()?;
        let n = points.len();
        let at = |t: f64| self.point(t);
        // The breaks: the ends and the samples between them, with the lengths to each.
        let mut ts = vec![t0];
        let mut k = t0.floor() + 1.0;
        while k < t1 {
            ts.push(k);
            k += 1.0;
        }
        ts.push(t1);
        let ps: Vec<Point> = ts.iter().map(|&t| if t.fract() == 0.0 { points[(t as usize) % n] } else { at(t) }).collect();
        let mut lengths = vec![0.0];
        for w in ps.windows(2) {
            lengths.push(lengths.last().expect("a length") + w[0].distance(w[1]));
        }
        let total = *lengths.last().expect("a length");
        if !(total > 0.0) {
            return None;
        }
        Some(move |f: f64| {
            let target = f * total;
            let i = lengths.partition_point(|&l| l < target).clamp(1, lengths.len() - 1);
            let (l0, l1) = (lengths[i - 1], lengths[i]);
            let w = if l1 > l0 { (target - l0) / (l1 - l0) } else { 0.0 };
            ts[i - 1] + (ts[i] - ts[i - 1]) * w
        })
    }

    fn edge_curve(&self, t0: f64, t1: f64, a: Point, b: Point) -> (Curve, f64, f64) {
        match self {
            Section::Open(c, ..) | Section::Closed(c, _) => (c.clone(), t0, t1),
            Section::Chain(..) | Section::Sampled(..) => {
                // Interpolated at the section's own parameters, the samples doubled until the
                // curve follows the section between them (to 1e-9 of its size).
                let size = 1.0 + a.distance(b) + a.distance(self.point(0.5 * (t0 + t1)));
                let mut n = 48;
                // (A marched section's samples are spaced as its steps were, unevenly: its
                // points taken evenly along it instead, by the lengths between the samples.)
                let along = self.even_parameters(t0, t1);
                let curve = loop {
                    let mut params: Vec<f64> = (0..=n).map(|i| along.as_ref().map_or(t0 + (t1 - t0) * i as f64 / n as f64, |f| f(i as f64 / n as f64))).collect();
                    let mut points: Vec<Point> = params.iter().map(|&t| self.point(t)).collect();
                    points[0] = a;
                    points[n] = b;
                    // Points that coincide (a stretch where the samples bunch up, or a sample
                    // next to the end vertex it was merged into) are one.
                    let spacing = on_tol().max(1e-9 * (1.0 + a.distance(b)));
                    let mut kept: Vec<Point> = Vec::with_capacity(points.len());
                    let mut kept_params: Vec<f64> = Vec::with_capacity(points.len());
                    for (p, t) in points.into_iter().zip(params.drain(..)) {
                        if kept.last().is_none_or(|q: &Point| q.distance(p) > spacing) {
                            kept.push(p);
                            kept_params.push(t);
                        }
                    }
                    if kept.last().is_some_and(|q| q.distance(b) > spacing) {
                        kept.push(b);
                        kept_params.push(t1);
                    } else if let Some(last) = kept.last_mut() {
                        *last = b;
                        *kept_params.last_mut().expect("a point") = t1;
                    }
                    // Each point advancing along the section's way there (from the points
                    // some way before it to some way after it), else dropped: where the
                    // surfaces nearly touch, points pulled onto the section may slide along
                    // it, back over its stretch before.
                    if matches!(self, Section::Sampled(..)) && kept.len() > 2 {
                        let last = kept.len() - 1;
                        let (mut points, mut ts) = (vec![kept[0]], vec![kept_params[0]]);
                        for i in 1..=last {
                            let way = kept[(i + 4).min(last)] - kept[i.saturating_sub(4)];
                            if i == last || (kept[i] - *points.last().expect("a point")).dot(way) > 0.0 {
                                points.push(kept[i]);
                                ts.push(kept_params[i]);
                            }
                        }
                        (kept, kept_params) = (points, ts);
                    }
                    // (Straight from end to end where nothing else can be: a section all but a
                    // point, its ends maybe one.)
                    let straight = || crate::geometry::curve::BSplineCurve::new(vec![a, b], vec![0.0, 1.0], vec![2, 2], 1);
                    if kept.len() < 2 {
                        break fit::interpolate(&[a, b], false, 0.0).unwrap_or_else(|_| straight());
                    }
                    // Interpolated by chord length (a chain's own parameter slows to a stop at
                    // the ends of its branches); its accuracy measured from the section's points
                    // between the samples.
                    // (Where chord lengths do not increase, points all but coinciding: by the
                    // section's own parameter, else straight.)
                    let curve = match fit::interpolate(&kept, false, 0.0).or_else(|_| fit::interpolate_with_parameters(&kept, &kept_params)) {
                        Ok(curve) => curve,
                        Err(_) => break fit::interpolate(&[a, b], false, 0.0).unwrap_or_else(|_| straight()),
                    };
                    if matches!(self, Section::Sampled(..)) {
                        break curve;
                    }
                    let as_curve = Curve::BSpline(curve.clone());
                    let (c0, c1) = (curve.first_parameter(), curve.last_parameter());
                    let samples = geom::projection_samples(&as_curve, c0, c1);
                    let error = kept_params
                        .windows(2)
                        .map(|w| geom::project_among(&as_curve, c0, c1, Some(&samples), self.point(0.5 * (w[0] + w[1]))).1)
                        .fold(0.0, f64::max);
                    // (A traced section is a polyline in the base's parameters between its
                    // samples, not an exact curve to converge to: its first fit stands.)
                    if error <= 1e-9 * size || n >= 3072 || matches!(self, Section::Sampled(..)) {
                        break curve;
                    }
                    n *= 2;
                };
                let (first, last) = (curve.first_parameter(), curve.last_parameter());
                (Curve::BSpline(curve), first, last)
            }
        }
    }
}

/// The curves where two faces' surfaces meet; lines are restricted to the faces' boxes,
/// branches on a cylinder or cone to near the face's range of `v`.
/// The size of a face: its box's diagonal.
fn size(face: &Face) -> f64 {
    face.bbox.0.distance(face.bbox.1)
}

fn sections(fa: &Face, fb: &Face, crossings: &[Point]) -> Vec<Section> {
    let (a, b) = (&fa.surface, &fb.surface);
    let boxes = [&fa.bbox, &fb.bbox];
    if let Some(section) = plane_across_extrusion(a, b).or_else(|| plane_across_extrusion(b, a)) {
        return vec![section];
    }
    if let Some(sections) = plane_along_extrusion(fa, fb) {
        return sections;
    }
    if let Some(sections) = meridians(a, b).or_else(|| meridians(b, a)) {
        return sections;
    }
    match (a, b) {
        (Surface::Revolution(r), _) if coaxial(r, b) => return parallels(fa, fb),
        (_, Surface::Revolution(r)) if coaxial(r, a) => return parallels(fb, fa),
        // (A free-form face smaller than the swept one it meets: the free-form one traced,
        // against the swept one's distance, a projection onto a curve rather than onto a
        // patch. The smaller face, for the grid it is traced on to be fine enough.)
        (Surface::Revolution(_) | Surface::Extrusion(_), s) if geom::freeform(s) && size(fb) < size(fa) => return sliced(fb, fa, crossings),
        (s, Surface::Revolution(_) | Surface::Extrusion(_)) if geom::freeform(s) && size(fa) < size(fb) => return sliced(fa, fb, crossings),
        (Surface::Revolution(_) | Surface::Extrusion(_), _) => return sliced(fa, fb, crossings),
        (_, Surface::Revolution(_) | Surface::Extrusion(_)) => return sliced(fb, fa, crossings),
        // The same free-form surface twice (pieces of one loft made twice) meets itself
        // everywhere: no section curves, like coincident quadrics.
        (s, t) if geom::freeform(s) && geom::freeform(t) && same_freeform(s, t) => return Vec::new(),
        // A free-form surface is traced in its parameters against the other's signed
        // distance (the cheaper one when both are free-form).
        (s, _) if geom::freeform(s) => return sliced(fa, fb, crossings),
        (_, s) if geom::freeform(s) => return sliced(fb, fa, crossings),
        _ => {}
    }
    // The stretch of a line within both faces' boxes (where the section can be).
    let line_range = |line: &Line| {
        let along = |bx: &(Point, Point)| {
            let (lo, hi) = *bx;
            let ts: Vec<f64> = (0..8)
                .map(|i| {
                    let corner = Point::new(
                        if i & 1 == 0 { lo.x } else { hi.x },
                        if i & 2 == 0 { lo.y } else { hi.y },
                        if i & 4 == 0 { lo.z } else { hi.z },
                    );
                    (corner - line.origin).dot(Vector::from(line.direction))
                })
                .collect();
            (ts.iter().copied().fold(f64::INFINITY, f64::min), ts.iter().copied().fold(f64::NEG_INFINITY, f64::max))
        };
        let ((lo_a, hi_a), (lo_b, hi_b)) = (along(boxes[0]), along(boxes[1]));
        (lo_a.max(lo_b) - 1.0, hi_a.min(hi_b) + 1.0)
    };
    // Where surfaces are nearly but not exactly parallel to a line (a plane tilted within
    // the angular tolerance from a cylinder's axis), the intersection is a long, nearly
    // straight curve: the line through two of its points, found near the faces, follows it
    // closely there.
    let through_curve = |l: Line| {
        let (lo, hi) = line_range(&l);
        // Newton steps on both surfaces' signed distances, along their normals.
        let on_both = |mut p: Point| {
            for _ in 0..20 {
                let (da, db) = (geom::signed_distance(a, p), geom::signed_distance(b, p));
                let (na, nb) = (geom::normal_at_point(a, p), geom::normal_at_point(b, p));
                let (sa, sb) = (
                    if geom::signed_distance(a, p + na * 1e-6) > da { 1.0 } else { -1.0 },
                    if geom::signed_distance(b, p + nb * 1e-6) > db { 1.0 } else { -1.0 },
                );
                let (na, nb) = (na * sa, nb * sb);
                let (m11, m12, m22) = (na.dot(na), na.dot(nb), nb.dot(nb));
                let det = m11 * m22 - m12 * m12;
                if det.abs() <= 1e-18 {
                    break;
                }
                let (alpha, beta) = ((-da * m22 + db * m12) / det, (-db * m11 + da * m12) / det);
                p = p + na * alpha + nb * beta;
                if da.abs().max(db.abs()) <= 1e-15 {
                    break;
                }
            }
            p
        };
        // Lines of a plane and a cylinder only (others, through a cone's apex or of planes,
        // are exact), when they stray.
        if !matches!((a, b), (Surface::Plane(_), Surface::Cylinder(_)) | (Surface::Cylinder(_), Surface::Plane(_))) {
            return l;
        }
        let off = |p: Point| geom::signed_distance(a, p).abs().max(geom::signed_distance(b, p).abs());
        if off(l.point(lo)).max(off(l.point(hi))) <= 0.1 * on_tol() {
            return l;
        }
        let (p, q) = (on_both(l.point(lo)), on_both(l.point(hi)));
        match crate::math::Direction::try_from_vector(q - p) {
            Some(d) if Vector::from(d).dot(Vector::from(l.direction)) > 0.0 => Line::new(p, d),
            _ => l,
        }
    };
    // Parallel planes within the meeting tolerance are one plane (shapes made earlier put
    // faces that should meet only that close).
    if let (Surface::Plane(p), Surface::Plane(q)) = (a, b) {
        let (np, nq) = (geom::plane_normal(p), geom::plane_normal(q));
        if np.cross(nq).length() <= 1e-12 && (q.origin() - p.origin()).dot(np).abs() <= on_tol() {
            return Vec::new();
        }
    }
    // (A cone's conics other than circles are traced as before, by angle about its axis.)
    let found = match intersect::quadrics(a, b, TOL) {
        QuadricIntersection::Ellipses(_) | QuadricIntersection::Parabolas(_) | QuadricIntersection::Hyperbolas(_) if matches!(a, Surface::Cone(_)) || matches!(b, Surface::Cone(_)) => QuadricIntersection::NotAnalytic,
        found => found,
    };
    match found {
        QuadricIntersection::Lines(lines) => lines
            .into_iter()
            .map(|l| {
                let l = through_curve(l);
                let (lo, hi) = line_range(&l);
                Section::Open(Curve::Line(l), lo, hi)
            })
            .collect(),
        QuadricIntersection::Circles(cs) => cs.into_iter().map(|c| Section::Closed(Curve::Circle(c), Default::default())).collect(),
        QuadricIntersection::Ellipses(es) => es.into_iter().map(|e| Section::Closed(Curve::Ellipse(e), Default::default())).collect(),
        QuadricIntersection::NotAnalytic => {
            // Parametrize on a cylinder, else a cone, by angle around its axis.
            let rank = |f: &Face| match f.surface {
                Surface::Cylinder(_) => 0,
                Surface::Cone(_) => 1,
                _ => 2,
            };
            let (base, other) = if rank(fa) <= rank(fb) { (fa, fb) } else { (fb, fa) };
            let (lo, hi) = base.v_range;
            let margin = 0.5 * (hi - lo) + 1.0;
            let Some(curves) = intersect::ruled_quadric(&base.surface, &other.surface, TOL, Some((lo - margin, hi + margin))) else {
                todo!("intersecting these surfaces")
            };
            chains(&curves)
        }
        QuadricIntersection::Empty | QuadricIntersection::Same | QuadricIntersection::Points(_) => Vec::new(),
        other => todo!("section curves {other:?}"),
    }
}

/// True if `other` is a surface of revolution about the same axis as `r` (a plane square to
/// it counts as one).
fn coaxial(r: &RevolutionSurface, other: &Surface) -> bool {
    let axis = r.axis();
    let parallel = |d: Direction| Vector::from(d).cross(Vector::from(axis.direction)).length() <= 1e-9;
    let on_axis = |p: Point| {
        let d = p - axis.origin;
        (d - Vector::from(axis.direction) * d.dot(Vector::from(axis.direction))).length() <= on_tol()
    };
    match other {
        Surface::Plane(p) => parallel(p.normal()),
        Surface::Cylinder(c) => parallel(c.frame().z()) && on_axis(c.frame().origin()),
        Surface::Cone(c) => parallel(c.frame().z()) && on_axis(c.frame().origin()),
        Surface::Sphere(s) => on_axis(s.frame().origin()),
        Surface::Revolution(o) => parallel(o.axis().direction) && on_axis(o.axis().origin),
        _ => false,
    }
}

/// Where a surface of revolution whose meridian lies in a plane through its axis (as a
/// torus's) meets a plane through that axis: its meridian turned into the plane, on either
/// side of the axis. `None` for other pairs.
fn meridians(revolution: &Surface, plane: &Surface) -> Option<Vec<Section>> {
    use crate::math::{Axis, Transform, Transformable};
    let (Surface::Revolution(r), Surface::Plane(plane)) = (revolution, plane) else { return None };
    let axis = r.axis();
    let z = Vector::from(axis.direction);
    let n = Vector::from(plane.normal());
    if n.dot(z).abs() > 1e-9 || n.dot(axis.origin - plane.origin()).abs() > on_tol() {
        return None;
    }
    let basis = r.basis();
    let (first, last) = (basis.first_parameter(), basis.last_parameter());
    if !(first.is_finite() && last.is_finite()) {
        return None;
    }
    // The meridian's side of the axis: all its points off the axis there.
    let radial = |p: Point| {
        let d = p - axis.origin;
        d - z * d.dot(z)
    };
    let samples: Vec<Vector> = (0..=32).map(|i| radial(basis.point(first + (last - first) * i as f64 / 32.0))).collect();
    let side = samples.iter().copied().max_by(|a, b| a.length().total_cmp(&b.length()))?;
    if !(side.length() > on_tol()) {
        return None;
    }
    let side = side.normalized();
    if samples.iter().any(|v| v.cross(side).length() > on_tol() || v.dot(side) < -on_tol()) {
        return None;
    }
    let across = n.cross(z).normalized();
    let closed = matches!(basis, Curve::Circle(_) | Curve::Ellipse(_)) && last - first >= TAU - 1e-12;
    Some(
        [across, -across]
            .iter()
            .map(|&to| {
                let angle = z.dot(side.cross(to)).atan2(side.dot(to));
                let curve = basis.transformed(&Transform::rotation(Axis::new(axis.origin, axis.direction), angle));
                if closed { Section::Closed(curve, Default::default()) } else { Section::Open(curve, first, last) }
            })
            .collect(),
    )
}

/// The circles where a surface of revolution (`base`) meets a coaxial surface: the
/// parallels through the points where its meridian crosses or touches the other surface.
fn parallels(base: &Face, other: &Face) -> Vec<Section> {
    use crate::geometry::ParametricSurface;
    let Surface::Revolution(r) = &base.surface else { unreachable!() };
    let (s, q) = (&base.surface, &other.surface);
    let (_, _, bv0, bv1) = s.bounds();
    let (v0, v1) = match base.v_period {
        Some(period) => (base.v_range.0, base.v_range.0 + period),
        None => {
            let m = 0.1 * (base.v_range.1 - base.v_range.0);
            let (lo, hi) = (base.v_range.0 - m, base.v_range.1 + m);
            if base.v_period.is_some() { (lo, hi) } else { (lo.max(bv0), hi.min(bv1)) }
        }
    };
    let u = 0.5 * (base.u_range.0 + base.u_range.1);
    let f = |v: f64| geom::signed_distance(q, s.point(u, v));
    let n = 2000;
    let step = (v1 - v0) / n as f64;
    let values: Vec<f64> = (0..=n).map(|i| f(v0 + step * i as f64)).collect();
    // The same surface (coaxial surfaces of revolution meeting along a stretch of meridian
    // coincide): no section curves, like coincident quadrics.
    if values.iter().filter(|d| d.abs() <= on_tol()).count() > n / 2 {
        return Vec::new();
    }
    // Bisection on [a, b] for a function changing sign there.
    let bisect = |g: &dyn Fn(f64) -> f64, mut a: f64, mut b: f64| {
        let ga = g(a) >= 0.0;
        for _ in 0..80 {
            let m = 0.5 * (a + b);
            if (g(m) >= 0.0) == ga {
                a = m;
            } else {
                b = m;
            }
        }
        0.5 * (a + b)
    };
    let mut roots: Vec<f64> = Vec::new();
    for i in 0..n {
        let (a, b) = (values[i], values[i + 1]);
        let (va, vb) = (v0 + step * i as f64, v0 + step * (i + 1) as f64);
        if a == 0.0 {
            roots.push(va);
        } else if (a > 0.0) != (b > 0.0) && b != 0.0 {
            roots.push(bisect(&f, va, vb));
        }
    }
    // Touching: |f| has a local minimum near zero without a sign change; it is a root of
    // the derivative.
    let h = 1e-7 * (v1 - v0);
    let df = |v: f64| (f(v + h) - f(v - h)) / (2.0 * h);
    for i in 1..n {
        let (a, b, c) = (values[i - 1].abs(), values[i].abs(), values[i + 1].abs());
        if b > a || b > c || (values[i - 1] > 0.0) != (values[i + 1] > 0.0) {
            continue;
        }
        let (va, vc) = (v0 + step * (i - 1) as f64, v0 + step * (i + 1) as f64);
        let v = if (df(va) > 0.0) != (df(vc) > 0.0) { bisect(&df, va, vc) } else { v0 + step * i as f64 };
        if f(v).abs() <= on_tol() {
            roots.push(v);
        }
    }
    roots.sort_by(f64::total_cmp);
    roots.dedup_by(|x, y| (*x - *y).abs() <= 2.0 * step);
    let axis = r.axis();
    let z = Vector::from(axis.direction);
    roots
        .into_iter()
        .filter_map(|v| {
            let p = s.point(u, v);
            let center = axis.origin + z * (p - axis.origin).dot(z);
            let radial = p - center;
            (radial.length() > on_tol()).then(|| {
                Section::Closed(Curve::Circle(Circle::new(Frame::new(center, axis.direction, Direction::from_vector(radial)), radial.length())), Default::default())
            })
        })
        .collect()
}

/// Whether `a`'s surface lies on `b`'s over `a`'s range (the same surface, or one made from
/// it): a grid of its points all on `b`'s surface, to far within the tolerance (not merely
/// near it, as planes a hair apart).
fn coincide(a: &Face, b: &Face) -> bool {
    use crate::geometry::ParametricSurface;
    if a.operand == b.operand {
        return false;
    }
    (0..=4).all(|i| {
        (0..=4).all(|j| {
            let u = a.u_range.0 + (a.u_range.1 - a.u_range.0) * (0.1 + 0.2 * i as f64);
            let v = a.v_range.0 + (a.v_range.1 - a.v_range.0) * (0.1 + 0.2 * j as f64);
            geom::signed_distance(&b.surface, a.surface.point(u, v)).abs() <= 1e-3 * on_tol()
        })
    })
}

/// How far the curve over `[first, last]` strays from `surface` (sampled).
fn farthest(curve: &Curve, first: f64, last: f64, surface: &Surface) -> f64 {
    (0..=32).map(|i| geom::signed_distance(surface, curve.point(first + (last - first) * i as f64 / 32.0)).abs()).fold(0.0, f64::max)
}

/// The curves where a face's surface without a closed-form intersection (`base`, swept by
/// its parameters) meets an elementary one: the zero contour of the other surface's signed
/// distance over a grid in the base's `(u, v)` (marching squares, wrapping round in a
/// periodic `u`), with crossing points found exactly on the grid lines.
/// True if two free-form surfaces are the same: equal degrees, knots and weights, poles
/// within the meeting tolerance.
fn same_freeform(a: &Surface, b: &Surface) -> bool {
    let close = |p: &[Vec<Point>], q: &[Vec<Point>]| {
        p.len() == q.len() && p.iter().zip(q).all(|(r, s)| r.len() == s.len() && r.iter().zip(s).all(|(x, y)| x.distance(*y) <= on_tol()))
    };
    let near = |x: &[f64], y: &[f64]| x.len() == y.len() && x.iter().zip(y).all(|(a, b)| (a - b).abs() <= 1e-12 * (1.0 + a.abs()));
    match (a, b) {
        (Surface::BSpline(p), Surface::BSpline(q)) => {
            let (pu, pv, qu, qv) = (p.u_knot_vector(), p.v_knot_vector(), q.u_knot_vector(), q.v_knot_vector());
            pu.degree == qu.degree
                && pv.degree == qv.degree
                && pu.multiplicities == qu.multiplicities
                && pv.multiplicities == qv.multiplicities
                && near(&pu.knots, &qu.knots)
                && near(&pv.knots, &qv.knots)
                && close(p.poles(), q.poles())
                && match (p.weights(), q.weights()) {
                    (None, None) => true,
                    (Some(x), Some(y)) => x.len() == y.len() && x.iter().zip(y).all(|(r, s)| near(r, s)),
                    _ => false,
                }
        }
        (Surface::Bezier(p), Surface::Bezier(q)) => {
            close(p.poles(), q.poles())
                && match (p.weights(), q.weights()) {
                    (None, None) => true,
                    (Some(x), Some(y)) => x.len() == y.len() && x.iter().zip(y).all(|(r, s)| near(r, s)),
                    _ => false,
                }
        }
        _ => false,
    }
}

/// True if a free-form surface surely misses a plane: its poles (whose hull holds it) all
/// lie on one side, beyond the tolerance.
fn misses(s: &Surface, q: &Surface) -> bool {
    let Surface::Plane(_) = q else { return false };
    let poles: Vec<Point> = match s {
        Surface::BSpline(b) => b.poles().iter().flatten().copied().collect(),
        Surface::Bezier(b) => b.poles().iter().flatten().copied().collect(),
        _ => return false,
    };
    let d: Vec<f64> = poles.iter().map(|&p| geom::signed_distance(q, p)).collect();
    d.iter().all(|&x| x > TOL) || d.iter().all(|&x| x < -TOL)
}

/// Grid cells across `u` and `v` for tracing sections on a surface: fine for surfaces of
/// revolution and extrusion (whose sections may wind), by knot spans and degrees for
/// B-spline and Bezier patches.
fn grid_size(s: &Surface) -> (usize, usize) {
    let by = |spans: usize, degree: usize| (4 * spans * degree.max(1)).clamp(16, 240);
    match s {
        Surface::BSpline(b) => (by(b.u_knots().len() - 1, b.u_degree()), by(b.v_knots().len() - 1, b.v_degree())),
        Surface::Bezier(b) => (by(1, b.u_degree()), by(1, b.v_degree())),
        _ => (720, 120),
    }
}

fn sliced(base: &Face, other: &Face, crossings: &[Point]) -> Vec<Section> {
    use crate::geometry::ParametricSurface;
    let (s, q) = (&base.surface, &other.surface);
    if misses(s, q) {
        return Vec::new();
    }
    let (bu0, bu1, bv0, bv1) = s.bounds();
    let (nu, nv) = grid_size(s);
    // The grid is offset by an irrational fraction of a step, so that a section along an
    // iso-line of the base (a meridian of a torus cut by a plane through its axis) does
    // not run along a grid line.
    let shift = 0.381_966_011_250_105;
    let periodic = base.u_period.is_some() && base.u_range.1 - base.u_range.0 >= TAU - 1e-9;
    let (u0, u1) = if periodic {
        let u0 = base.u_range.0 + shift * TAU / nu as f64;
        (u0, u0 + TAU)
    } else {
        // Within the surface's bounds, unless u is periodic (the face's range may lie
        // across them).
        let span = base.u_range.1 - base.u_range.0;
        let (lo, hi) = (base.u_range.0 - 0.05 * span - shift * 1.1 * span / nu as f64, base.u_range.1 + 0.05 * span);
        if base.u_period.is_some() { (lo, hi) } else { (lo.max(bu0), hi.min(bu1)) }
    };
    let v_periodic = base.v_period.filter(|&period| base.v_range.1 - base.v_range.0 >= period - 1e-9);
    let (v0, v1) = match v_periodic {
        Some(period) => {
            let v0 = base.v_range.0 + shift * period / nv as f64;
            (v0, v0 + period)
        }
        None => {
            let span = base.v_range.1 - base.v_range.0;
            let (lo, hi) = (base.v_range.0 - 0.1 * span - shift * 1.2 * span / nv as f64, base.v_range.1 + 0.1 * span);
            if base.v_period.is_some() { (lo, hi) } else { (lo.max(bv0), hi.min(bv1)) }
        }
    };
    // Where the base is much larger than the other face, only the part of it near the
    // other's box (found from samples, widened by their spacing) is traced: the grid is
    // fine enough there for sections small on the base.
    let (u0, u1, v0, v1) = {
        const N: usize = 64;
        let (lo, hi) = other.extent;
        let size = lo.distance(hi);
        // (The same for every other face: kept with the base. Evaluated row by row along u,
        // as the B-spline batches them, the same points.)
        let (samples, spacing) = base.samples.get_or_init(|| {
            let us: Vec<f64> = (0..=N).map(|i| u0 + (u1 - u0) * i as f64 / N as f64).collect();
            let vs: Vec<f64> = (0..=N).map(|j| v0 + (v1 - v0) * j as f64 / N as f64).collect();
            let rows = s.points_on_grid(&us, &vs);
            let samples: Vec<Vec<Point>> = (0..=N).map(|i| (0..=N).map(|j| rows[j][i]).collect()).collect();
            let spacing = (0..N).flat_map(|i| (0..N).map(move |j| (i, j))).map(|(i, j)| samples[i][j].distance(samples[i + 1][j]).max(samples[i][j].distance(samples[i][j + 1]))).fold(0.0, f64::max);
            (samples, spacing)
        });
        let spacing = *spacing;
        let margin = spacing + 0.05 * size;
        let near = |p: Point| p.x >= lo.x - margin && p.x <= hi.x + margin && p.y >= lo.y - margin && p.y <= hi.y + margin && p.z >= lo.z - margin && p.z <= hi.z + margin;
        let (mut imin, mut imax, mut jmin, mut jmax) = (usize::MAX, 0, usize::MAX, 0);
        for i in 0..=N {
            for j in 0..=N {
                if near(samples[i][j]) {
                    (imin, imax, jmin, jmax) = (imin.min(i), imax.max(i), jmin.min(j), jmax.max(j));
                }
            }
        }
        if imin == usize::MAX {
            (u0, u1, v0, v1)
        } else {
            let at = |a: f64, b: f64, k: usize| a + (b - a) * k as f64 / N as f64;
            let (iu, ju) = if periodic { (0, N) } else { (imin.saturating_sub(1), (imax + 1).min(N)) };
            let (iv, jv) = if v_periodic.is_some() { (0, N) } else { (jmin.saturating_sub(1), (jmax + 1).min(N)) };
            (at(u0, u1, iu), at(u0, u1, ju), at(v0, v1, iv), at(v0, v1, jv))
        }
    };
    // Traced by marching along the curves, where that is done: seeds sought on the grid's
    // lattice, steps no longer than its cells.
    let marched = || {
        let spacing = base.samples.get().map_or(0.0, |s| s.1);
        let range = super::march::Range { u: (u0, u1), v: (v0, v1), u_period: if periodic { base.u_period } else { None }, v_period: v_periodic, u_wrap: base.u_period, v_wrap: base.v_period };
        super::march::trace(s, q, range, (nu, nv), spacing * 64.0 / nu.max(nv) as f64, crossings, 10.0 * (base.drift.max(other.drift) + on_tol()), (other.u_range, other.v_range))
    };
    // (Where either surface is free-form, or a surface of extrusion: on its grid of signed
    // distances, as fine as a quadric's, sections are slow to find.)
    let extrusion = |x: &Surface| matches!(x, Surface::Extrusion(_));
    if geom::freeform(s) || geom::freeform(q) || extrusion(s) || extrusion(q) {
        if let Some(chains) = marched() {
            let mut out = Vec::new();
            for super::march::Chain { uv, st, closed } in chains {
                let traced = Traced { uv, u_period: base.u_period, v_period: v_periodic, st: Some(st), base: s.clone(), other: q.clone(), last: Default::default() };
                let points: Vec<Point> = traced.uv.iter().map(|&(u, v)| s.point(u, v)).collect();
                out.push(Section::Sampled(points, closed, Box::new(traced)));
            }
            return out;
        }
    }
    // Along the grid, each foot point on a free-form other surface starts the next search.
    let last_foot = std::cell::Cell::new(None);
    let f = |u: f64, v: f64| {
        let (d, foot) = geom::signed_distance_from(q, s.point(u, v), last_foot.get());
        last_foot.set(foot);
        d
    };
    let columns = if periodic { nu } else { nu + 1 };
    let rows = if v_periodic.is_some() { nv } else { nv + 1 };
    let u_at = |i: usize| u0 + (u1 - u0) * i as f64 / nu as f64;
    let v_at = |j: usize| v0 + (v1 - v0) * j as f64 / nv as f64;
    // The grid's points evaluated row by row along u, then the distances in the order the
    // foot points chain.
    let us: Vec<f64> = (0..columns).map(u_at).collect();
    let grid = s.points_on_grid(&us, &(0..rows).map(v_at).collect::<Vec<f64>>());
    // (Only the signs matter, and a distance moves by no more than the point does: a point
    // farther from the other surface than from a neighbour has the neighbour's sign, a
    // distance at least as far as it less the gap kept for it. Not on a free-form other,
    // whose feet chain.)
    let lipschitz = !geom::freeform(q);
    let margin = 1e3 * on_tol();
    // (Blocks of nodes all far from the other surface, by the distance at the block's middle
    // less how far its points (and their neighbours beyond) reach from it, and farther than
    // the block across: one sign, none of its nodes next to a node of the other.)
    const BLOCK: usize = 8;
    let filled: Vec<Vec<Option<f64>>> = if lipschitz && columns >= 2 * BLOCK && rows >= 2 * BLOCK {
        let mut filled = vec![vec![None; rows]; columns];
        for bi in (0..columns).step_by(BLOCK) {
            for bj in (0..rows).step_by(BLOCK) {
                let (ie, je) = ((bi + BLOCK).min(columns), (bj + BLOCK).min(rows));
                let mut b = crate::math::BoundingBox::EMPTY;
                for i in bi.saturating_sub(1)..(ie + 1).min(columns) {
                    for j in bj.saturating_sub(1)..(je + 1).min(rows) {
                        b.add_point(grid[j][i]);
                    }
                }
                let c = grid[(bj + je) / 2][(bi + ie) / 2];
                let (d, _) = geom::signed_distance_from(q, c, None);
                let corner = |k: usize| Point::new(if k & 1 == 0 { b.min.x } else { b.max.x }, if k & 2 == 0 { b.min.y } else { b.max.y }, if k & 4 == 0 { b.min.z } else { b.max.z });
                let reach = d.abs() - (0..8).map(|k| c.distance(corner(k))).fold(0.0, f64::max);
                if reach > margin + (b.max - b.min).length() {
                    for column in &mut filled[bi..ie] {
                        for value in &mut column[bj..je] {
                            *value = Some(reach.copysign(d));
                        }
                    }
                }
            }
        }
        filled
    } else {
        Vec::new()
    };
    let mut values: Vec<Vec<f64>> = Vec::with_capacity(columns);
    for i in 0..columns {
        let mut column = Vec::with_capacity(rows);
        for j in 0..rows {
            if let Some(value) = filled.get(i).and_then(|c| c[j]) {
                column.push(value);
                continue;
            }
            let p = grid[j][i];
            if lipschitz {
                let below = (j > 0).then(|| (column[j - 1], grid[j - 1][i]));
                let left = (i > 0).then(|| (values[i - 1][j], grid[j][i - 1]));
                let (reach, d) = [below, left].into_iter().flatten().map(|(d, q): (f64, Point)| (d.abs() - p.distance(q), d)).fold((0.0, 0.0), |a, b| if b.0 > a.0 { b } else { a });
                if reach > margin {
                    column.push(reach.copysign(d));
                    continue;
                }
            }
            let (d, foot) = geom::signed_distance_from(q, p, last_foot.get());
            last_foot.set(foot);
            column.push(d);
        }
        values.push(column);
    }
    let sign = |i: usize, j: usize| values[i % columns][j % rows] >= 0.0;
    // Whether the distance crosses zero between two neighbouring nodes: of different signs,
    // and near enough the other surface for that. (A distance moves by no more than the
    // point does, here along the base between them, at most a few times their chord: where
    // it is farther at both than that, the signs differ only as feet on the other's bounds
    // give them, the point not square to it there.)
    let near = 10.0 * (base.drift.max(other.drift) + on_tol());
    let crosses = |a: (usize, usize), b: (usize, usize)| {
        let (da, db) = (values[a.0 % columns][a.1 % rows], values[b.0 % columns][b.1 % rows]);
        (da >= 0.0) != (db >= 0.0) && da.abs() + db.abs() <= 2.0 * grid[a.1 % rows][a.0 % columns].distance(grid[b.1 % rows][b.0 % columns]) + near
    };
    // A root of g on [a, b] (g(a) and g(b) of different signs, zero counting as positive).
    // (Done once a and b are neighbouring numbers, where the steps would not move them,
    // unless the other surface is free-form: each of its distances starts the next.)
    let stateless = !geom::freeform(q);
    let bisect = |g: &dyn Fn(f64) -> f64, mut a: f64, mut b: f64| {
        let ga = g(a) >= 0.0;
        for _ in 0..60 {
            let m = 0.5 * (a + b);
            if stateless && m == a {
                break;
            }
            if (g(m) >= 0.0) == ga {
                a = m;
            } else {
                if stateless && m == b {
                    break;
                }
                b = m;
            }
        }
        0.5 * (a + b)
    };
    // Crossing points, keyed by grid edge: (i, j, along u?).
    let mut crossing: HashMap<(usize, usize, bool), (f64, f64)> = HashMap::default();
    let mut point_on = |i: usize, j: usize, along_u: bool| -> (f64, f64) {
        let key = (i % columns, j % rows, along_u);
        *crossing.entry(key).or_insert_with(|| {
            if along_u {
                let v = v_at(j);
                let a = u_at(i);
                (bisect(&|u| f(u, v), a, a + (u1 - u0) / nu as f64), v)
            } else {
                let u = u_at(i);
                (u, bisect(&|v| f(u, v), v_at(j), v_at(j + 1)))
            }
        })
    };
    // Segments between crossing points, cell by cell.
    let cells_u = if periodic { nu } else { nu };
    let mut links: HashMap<(usize, usize, bool), Vec<(usize, usize, bool)>> = HashMap::default();
    let mut link = |a: (usize, usize, bool), b: (usize, usize, bool)| {
        let norm = |k: (usize, usize, bool)| (k.0 % columns, k.1 % rows, k.2);
        let (a, b) = (norm(a), norm(b));
        links.entry(a).or_default().push(b);
        links.entry(b).or_default().push(a);
    };
    for i in 0..cells_u {
        let (here, next) = (&values[i], &values[(i + 1) % columns]);
        for j in 0..nv {
            let s00 = sign(i, j);
            // (No edge crossed where the corners' signs are all the same.)
            let up = if j + 1 == rows { 0 } else { j + 1 };
            if (next[j] >= 0.0) == s00 && (here[up] >= 0.0) == s00 && (next[up] >= 0.0) == s00 {
                continue;
            }
            // Cell edges: bottom, right, top, left.
            let edges = [((i, j, true), crosses((i, j), (i + 1, j))), ((i + 1, j, false), crosses((i + 1, j), (i + 1, j + 1))), ((i, j + 1, true), crosses((i, j + 1), (i + 1, j + 1))), ((i, j, false), crosses((i, j), (i, j + 1)))];
            let crossed: Vec<(usize, usize, bool)> = edges.iter().filter(|e| e.1).map(|e| e.0).collect();
            match crossed.len() {
                2 => link(crossed[0], crossed[1]),
                4 => {
                    // A saddle: the centre's sign tells which corners connect.
                    let centre = f(0.5 * (u_at(i) + u_at(i + 1)), 0.5 * (v_at(j) + v_at(j + 1))) >= 0.0;
                    if centre == s00 {
                        link(crossed[0], crossed[1]);
                        link(crossed[2], crossed[3]);
                    } else {
                        link(crossed[0], crossed[3]);
                        link(crossed[1], crossed[2]);
                    }
                }
                _ => {}
            }
        }
    }
    // Walk the links into polylines, from free ends first, then loops.
    let mut keys: Vec<(usize, usize, bool)> = links.keys().copied().collect();
    keys.sort_unstable();
    let mut used: HashMap<(usize, usize, bool), bool> = HashMap::default();
    let mut out = Vec::new();
    let starts: Vec<(usize, usize, bool)> =
        keys.iter().copied().filter(|k| links[k].len() == 1).chain(keys.iter().copied()).collect();
    for start in starts {
        if used.get(&start).copied().unwrap_or(false) {
            continue;
        }
        let mut chain = vec![start];
        used.insert(start, true);
        let mut current = start;
        let closed = loop {
            let next = links[&current].iter().copied().find(|k| !used.get(k).copied().unwrap_or(false));
            match next {
                Some(k) => {
                    used.insert(k, true);
                    chain.push(k);
                    current = k;
                }
                None => break chain.len() > 2 && links[&current].contains(&start),
            }
        };
        let mut uv: Vec<(f64, f64)> = chain.iter().map(|&(i, j, along)| point_on(i, j, along)).collect();
        // Unwrap u (and a periodic v) along the curve.
        let u_period = base.u_period.unwrap_or(TAU);
        for k in 1..uv.len() {
            let du = uv[k].0 - uv[k - 1].0;
            uv[k].0 -= (du / u_period).round() * u_period * if periodic { 1.0 } else { 0.0 };
            if let Some(period) = v_periodic {
                let dv = uv[k].1 - uv[k - 1].1;
                uv[k].1 -= (dv / period).round() * period;
            }
        }
        if uv.len() >= 2 {
            let mut traced = Traced { uv, u_period: base.u_period, v_period: v_periodic, st: None, base: s.clone(), other: q.clone(), last: Default::default() };
            // (Where it strays from the samples' chords by more than a few ten-thousandths of
            // the other face. The factor is tuned on the tests: finer, a section running along
            // an edge two faces share moves across it, where telling which face it lies in
            // is ill-conditioned.)
            // (On a free-form base, whose grid follows its poles and may be coarse.)
            if geom::freeform(s) {
                traced.uv = traced.densified(closed, 3e-4 * other.extent.0.distance(other.extent.1));
            }
            let points: Vec<Point> = traced.uv.iter().map(|&(u, v)| s.point(u, v)).collect();
            out.push(Section::Sampled(points, closed, Box::new(traced)));
        }
    }
    out
}

/// Whether surfaces `a` and `b` touch at `p` (on both, their normals parallel) crossing each
/// other there: the difference of their curvatures is indefinite, so that their section is
/// two branches crossing in `p` (not a curve along which they touch, as a fillet and its
/// face, nor `p` alone).
fn touching_across(a: &Surface, b: &Surface, p: Point) -> bool {
    use crate::geometry::ParametricSurface;
    let n = geom::normal_at_point(a, p);
    if n.cross(geom::normal_at_point(b, p)).length() > 1e-3 {
        return false;
    }
    // The normal curvature of a surface along a tangent direction (as seen along `n`).
    let curvature = |s: &Surface, w: Vector| -> Option<f64> {
        let [u, v] = geom::uv_of(s, p);
        let d = s.d2(u, v);
        let (e, f, g) = (d.du.dot(d.du), d.du.dot(d.dv), d.dv.dot(d.dv));
        let det = e * g - f * f;
        if !(det.abs() > 1e-300) {
            return None;
        }
        // w = α du + β dv (w is tangent).
        let (wu, wv) = (w.dot(d.du), w.dot(d.dv));
        let (alpha, beta) = ((g * wu - f * wv) / det, (e * wv - f * wu) / det);
        let second = (d.duu * (alpha * alpha) + d.duv * (2.0 * alpha * beta) + d.dvv * (beta * beta)).dot(n);
        let first = (d.du * alpha + d.dv * beta).length_squared();
        (first > 1e-300).then(|| second / first)
    };
    let e1 = {
        let x = if n.x.abs() < 0.9 { Vector::new(1.0, 0.0, 0.0) } else { Vector::new(0.0, 1.0, 0.0) };
        let x = x - n * x.dot(n);
        x * (1.0 / x.length())
    };
    let e2 = n.cross(e1);
    let diff = |w: Vector| -> Option<f64> { Some(curvature(a, w)? - curvature(b, w)?) };
    let (Some(h11), Some(h22), Some(hd)) = (diff(e1), diff(e2), diff((e1 + e2) * std::f64::consts::FRAC_1_SQRT_2)) else {
        return false;
    };
    // κ(d) for d = (e1 + e2)/√2 is (h11 + 2 h12 + h22) / 2.
    let h12 = hd - 0.5 * (h11 + h22);
    let size = h11 * h11 + h22 * h22 + 2.0 * h12 * h12;
    // (Not the same surface twice, where the difference is rounding.)
    size.sqrt() > 1e-9 && h11 * h22 - h12 * h12 < -1e-3 * size
}

/// The distance from `p` to the segment `ab`, and the parameter of the closest point.
fn segment_gap(a: Point, b: Point, p: Point) -> (f64, f64) {
    let ab = b - a;
    let s = ((p - a).dot(ab) / ab.length_squared().max(1e-300)).clamp(0.0, 1.0);
    ((a + ab * s).distance(p), s)
}

/// Branches joined end to end into chains: open ones from a free end, then closed loops.
fn chains(curves: &intersect::QuadricCurves) -> Vec<Section> {
    let n = curves.curves.len();
    let mut used = vec![false; n];
    let mut out = Vec::new();
    // Walk from branch `start`, entered at its `start_end` (0 = start, 1 = end).
    let walk = |start: usize, from_end: usize, used: &mut Vec<bool>| {
        let mut chain = Vec::new();
        let (mut current, mut entered) = (start, from_end);
        loop {
            used[current] = true;
            chain.push((curves.curves[current].clone(), entered == 1));
            let exit = 1 - entered;
            match curves.joins[current][exit] {
                Some(next) if !used[next] => {
                    // Joined branches meet at the same angle: entered at the matching end.
                    let at = |i: usize, e: usize| if e == 0 { curves.curves[i].domain().0 } else { curves.curves[i].domain().1 };
                    let angle = at(current, exit);
                    let near = |x: f64| ((x - angle).rem_euclid(TAU)).min((angle - x).rem_euclid(TAU));
                    entered = if near(at(next, 0)) <= near(at(next, 1)) { 0 } else { 1 };
                    current = next;
                }
                Some(next) => return (chain, next == start),
                None => return (chain, false),
            }
        }
    };
    for i in 0..n {
        if used[i] {
            continue;
        }
        let free = curves.joins[i].iter().position(|j| j.is_none());
        let full = {
            let (a, b) = curves.curves[i].domain();
            b - a >= TAU - 1e-9
        };
        if full {
            used[i] = true;
            out.push(Section::Chain(vec![(curves.curves[i].clone(), false)], true, ChainSamples::default()));
            continue;
        }
        if let Some(e) = free {
            let (chain, closed) = walk(i, e, &mut used);
            out.push(Section::Chain(chain, closed, ChainSamples::default()));
        }
    }
    for i in 0..n {
        if !used[i] {
            let (chain, closed) = walk(i, 0, &mut used);
            out.push(Section::Chain(chain, closed, ChainSamples::default()));
        }
    }
    out
}

impl Graph {
    /// Splits edges where they meet faces of the other operand, recording the meeting
    /// points on those faces, and the edges lying on their surfaces.
    pub fn intersect_edges_with_faces(&mut self) {
        /// Where an edge meets a face of the other operand.
        enum Meeting {
            /// A degenerate edge's vertex lies on the face.
            Vertex(usize),
            /// The edge lies on the face's surface.
            On,
            /// The edge crosses the face at these parameters and points (and, where it meets
            /// the face at its boundary, within the tolerances of a boundary edge, that edge
            /// and its parameter: both split at one vertex).
            At(Vec<(f64, Point, Option<(usize, f64)>)>),
        }
        // Each edge's meetings, found for all edges at once (locating points reads the
        // operands' own edges and vertices, which this does not change); then recorded in
        // order.
        let edges: Vec<usize> = (0..self.edges.len()).filter(|&ei| self.edges[ei].operand.is_some()).collect();
        // Each edge's own faces (of its operand).
        let mut own_faces: Vec<Vec<usize>> = vec![Vec::new(); self.edges.len()];
        for (fi, face) in self.faces.iter().enumerate() {
            for u in &face.boundary {
                if !own_faces[u.edge].contains(&fi) {
                    own_faces[u.edge].push(fi);
                }
            }
        }
        let found = crate::parallel::map(&edges, |&ei| {
            let edge = &self.edges[ei];
            let operand = edge.operand.expect("an operand's edge");
            let mut meetings = Vec::new();
            for fi in self.index().faces_overlapping_others(operand, &edge.bbox) {
                let face = &self.faces[fi];
                let Some(curve) = &edge.curve else {
                    let p = self.vertices.points[edge.start];
                    if geom::signed_distance(&face.surface, p).abs() <= on_tol() && face.locate(p, &self.edges, &self.vertices) != Location::Outside {
                        meetings.push((fi, Meeting::Vertex(edge.start)));
                    }
                    continue;
                };
                // (Clear of a free-form face's grid of boxes: no contact, without the search.)
                if geom::freeform(&face.surface) && clear_of_face(face, curve, edge.first, edge.last, edge.tolerance) {
                    continue;
                }
                match surface_contact(curve, edge.first, edge.last, &face.surface, edge.tolerance) {
                    Contact::On => meetings.push((fi, Meeting::On)),
                    // (An edge of a face on the very surface of this one lies on it as far as
                    // it strays from its own face's: a fitted section bounding a face cut to
                    // a surface strays from it between its fit points by more than its
                    // tolerance says.)
                    Contact::At(_)
                        if own_faces[ei].iter().any(|&own| {
                            let own = &self.faces[own];
                            let strays = farthest(curve, edge.first, edge.last, &own.surface);
                            strays <= 1e3 * on_tol() && coincide(own, face) && farthest(curve, edge.first, edge.last, &face.surface) <= strays + on_tol()
                        }) =>
                    {
                        meetings.push((fi, Meeting::On))
                    }
                    Contact::At(ts) => {
                        let ts = if edge.free { stretch_ends(curve, edge.first, edge.last, &face.surface, ts) } else { ts };
                        // (An edge may meet a face just past its boundary, by as far as the
                        // face's edges stray from it: tolerant edges read from a file.)
                        let at = ts
                            .into_iter()
                            .map(|t| (t, curve.point(t)))
                            .filter_map(|(t, p)| {
                                let location = face.locate(p, &self.edges, &self.vertices);
                                // (Inside, but met at a shallow angle (the edge lying nearly
                                // along the surface, as where a face tangent to the edge's own
                                // meets it), so found only roughly along the edge: where a
                                // boundary edge passes within the tolerances that near, the
                                // edge meets the face there, on its boundary (one vertex with
                                // its meeting the face across that edge).)
                                if location == Location::Inside {
                                    let snapped = self.meeting_at_boundary(curve, edge.first, edge.last, t, face, 0.0).filter(|m| m.2.is_none() && m.1.distance(p) <= 1e3 * on_tol());
                                    return Some(snapped.unwrap_or((t, p, None)));
                                }
                                // (On its boundary, or near it: a boundary edge passing within
                                // the tolerances is met there too; on a free-form face at any
                                // angle, within the face's own tolerance where that is wider
                                // (at most a thousand times the classifying distance), as
                                // where the rims of two filled faces' walls cross.)
                                if location == Location::Boundary || (face.drift > on_tol() && face.near_boundary(p, face.drift, &self.edges, &self.vertices)) {
                                    let within = if geom::freeform(&face.surface) && face.drift > on_tol() { face.drift.min(1e3 * on_tol()) } else { 0.0 };
                                    return match self.meeting_at_boundary(curve, edge.first, edge.last, t, face, within) {
                                        Some(crossing @ (_, _, Some(_))) => Some(crossing),
                                        _ => Some((t, p, None)),
                                    };
                                }
                                self.meeting_at_boundary(curve, edge.first, edge.last, t, face, 0.0)
                            })
                            .collect();
                        meetings.push((fi, Meeting::At(at)));
                    }
                }
            }
            meetings
        });
        for (&ei, meetings) in edges.iter().zip(found) {
            for (fi, meeting) in meetings {
                match meeting {
                    Meeting::Vertex(v) => self.faces[fi].vertices_on.push(v),
                    Meeting::On => self.faces[fi].on_edges.push(ei),
                    Meeting::At(at) => {
                        for (t, p, crossing) in at {
                            match crossing {
                                Some((other, s)) => {
                                    // (As far off the surfaces as the edges are apart there.)
                                    let gap = self.edges[other].point(s, &self.vertices).distance(p);
                                    let within = 2.0 * gap + on_tol();
                                    // (The vertex either edge has there already, from the
                                    // other's meeting this one's face: one for both.)
                                    let known = [ei, other].iter().flat_map(|&e| self.edges[e].paves.iter().map(|&(_, v)| v)).find(|&v| self.vertices.points[v].distance(p) <= within);
                                    let v = match known {
                                        Some(v) => {
                                            let q = self.vertices.points[v];
                                            let (a, b) = (&self.edges[ei], &self.edges[other]);
                                            let (t, s) = (geom::project(a.curve.as_ref().unwrap(), a.first, a.last, q).0, geom::project(b.curve.as_ref().unwrap(), b.first, b.last, q).0);
                                            self.add_pave(ei, t, v);
                                            self.add_pave(other, s, v);
                                            v
                                        }
                                        None => self.add_crossing(ei, t, other, s),
                                    };
                                    self.faces[fi].vertices_on.push(v);
                                    let allowance = self.loose_vertices.entry(v).or_insert(0.0);
                                    *allowance = allowance.max(within);
                                }
                                None => {
                                    let v = self.vertices.add(p);
                                    self.add_pave(ei, t, v);
                                    self.faces[fi].vertices_on.push(v);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    /// Splits two edges where they cross (at `t` on `e1`, `s` on `e2`). Where they cross at
    /// a small angle the point is found only roughly along them (by the tolerance over the
    /// angle's sine): an end vertex of either within that is where they meet.
    fn add_crossing(&mut self, e1: usize, t: f64, e2: usize, s: f64) -> usize {
        let p = self.edges[e1].point(t, &self.vertices);
        let (Some(c1), Some(c2)) = (self.edges[e1].curve.clone(), self.edges[e2].curve.clone()) else {
            let v = self.vertices.add(p);
            self.add_pave(e1, t, v);
            self.add_pave(e2, s, v);
            return v;
        };
        let (d1, d2) = (c1.d1(t).1, c2.d1(s).1);
        let sine = d1.cross(d2).length() / (d1.length() * d2.length()).max(1e-300);
        let reach = (on_tol() / sine.max(1e-12)).min(1e2 * on_tol());
        let ends = [self.edges[e1].start, self.edges[e1].end, self.edges[e2].start, self.edges[e2].end];
        let snapped = ends
            .into_iter()
            .map(|v| (v, self.vertices.points[v].distance(p)))
            .filter(|&(_, d)| d <= reach)
            .min_by(|a, b| a.1.total_cmp(&b.1));
        match snapped {
            Some((v, _)) => {
                let q = self.vertices.points[v];
                let (t, s) = {
                    let (a, b) = (&self.edges[e1], &self.edges[e2]);
                    (geom::project(&c1, a.first, a.last, q).0, geom::project(&c2, b.first, b.last, q).0)
                };
                self.add_pave(e1, t, v);
                self.add_pave(e2, s, v);
                v
            }
            None => {
                let v = self.vertices.add(p);
                self.add_pave(e1, t, v);
                self.add_pave(e2, s, v);
                v
            }
        }
    }

    /// Records a split point on an edge (ignoring its own end vertices: an edge does not
    /// pass through them again).
    /// Vertex `v`, put on edge `ei` by a section ending on it, put on the edges of the other
    /// operand lying along `ei` (on the same line or circle) through it too: where the two
    /// operands' faces meet along one edge each has (two blocks side by side, their faces
    /// touching there), a section met only one of them, and the other, left whole, did not
    /// match its pieces. (Not on an edge a section ended on near there itself: each has its
    /// own vertex there, and both on each would leave a sliver between them.)
    fn pave_coincident(&mut self, ei: usize, v: usize) {
        let (Some(curve), operand) = (self.edges[ei].curve.clone(), self.edges[ei].operand) else { return };
        let p = self.vertices.points[v];
        for j in 0..self.edges.len() {
            let e = &self.edges[j];
            if j == ei || e.operand.is_none() || e.operand == operand || v == e.start || v == e.end {
                continue;
            }
            let Some(c) = &e.curve else { continue };
            if !same_carrier(&curve, c) {
                continue;
            }
            let (s, d) = geom::project(c, e.first, e.last, p);
            let inside = s > e.first.min(e.last) + 1e-9 && s < e.first.max(e.last) - 1e-9;
            let reach = 1e-3 * c.point(e.first).distance(c.point(e.last));
            if d <= on_tol() && inside && !e.paves.iter().any(|&(_, w)| self.vertices.points[w].distance(p) <= reach) {
                self.add_pave(j, s, v);
            }
        }
    }

    pub fn add_pave(&mut self, ei: usize, t: f64, v: usize) {
        let e = &mut self.edges[ei];
        if v == e.start || v == e.end {
            return;
        }
        if !e.paves.iter().any(|&(_, w)| w == v) {
            e.paves.push((t, v));
        }
    }

    /// Splits edges of the two operands where they cross, when one of them is not on a face
    /// (crossings of edges on faces are found from the faces).
    pub fn intersect_free_edges(&mut self) {
        let mut crossings: Vec<(usize, usize, f64, f64)> = Vec::new();
        for (i, j) in self.index().edge_pairs() {
            let (a, b) = (&self.edges[i], &self.edges[j]);
            {
                if !(a.free || b.free) {
                    continue;
                }
                // Within the edges' default tolerances (as OCCT, whose edges carry 1e-7 each):
                // curves passing closer than on_tol but further than that only come near.
                for (t, s) in edge_crossings(a, b, 2.0 * TOL) {
                    crossings.push((i, j, t, s));
                }
            }
        }
        for (e1, e2, t, s) in crossings {
            self.add_crossing(e1, t, e2, s);
        }
    }

    /// Where an edge meeting a face's surface at a shallow angle (lying nearly along it)
    /// meets the face itself, found just outside the face: the meeting point is found only
    /// within the tolerance over the sine of that angle along the surface, so where the face's
    /// boundary is that near, the edge meets it there: at the edge's point nearest the
    /// boundary edge, where the two pass within the tolerance (as where two faces lying
    /// nearly on one another have their edges cross), the boundary edge given to be split
    /// there too. So too at any angle where the boundary edge passes `within` (the face's
    /// own tolerance, met on or near its boundary): the edges cross there, one vertex
    /// splitting both. `None` otherwise.
    fn meeting_at_boundary(&self, curve: &Curve, first: f64, last: f64, t: f64, face: &Face, within: f64) -> Option<(f64, Point, Option<(usize, f64)>)> {
        let p = curve.point(t);
        let tangent = curve.d1(t).1;
        let n = face.outward_normal(p);
        let (a, b) = (tangent.length(), n.length());
        if !(a > 0.0 && b > 0.0) {
            return None;
        }
        let sine = (tangent.dot(n) / (a * b)).abs();
        let tangent = sine < TANGENT_SINE;
        if !tangent && within <= 0.0 {
            return None;
        }
        let reach = if tangent { on_tol() / sine.max(1e-6) } else { 0.0 }.max(within);
        // The nearest boundary edge of the face, and the closest points of the two curves
        // (alternate projections from `p`).
        let nearest = face
            .boundary
            .iter()
            .filter(|b| b.bounds)
            .filter_map(|b| {
                let e = &self.edges[b.edge];
                let c = e.curve.as_ref()?;
                let (s, d) = geom::project(c, e.first, e.last, p);
                Some((d, b.edge, s))
            })
            .min_by(|x, y| x.0.total_cmp(&y.0))?;
        if nearest.0 > reach {
            return None;
        }
        let other = &self.edges[nearest.1];
        let other_curve = other.curve.as_ref()?;
        let (mut s, mut u) = (nearest.2, t);
        for _ in 0..16 {
            u = geom::project(curve, first, last, other_curve.point(s)).0;
            s = geom::project(other_curve, other.first, other.last, curve.point(u)).0;
        }
        let (x, y) = (curve.point(u), other_curve.point(s));
        let gap = x.distance(y);
        // (Within the boundary edge's tolerance, or the tolerance over the sine of the angle
        // the edge meets the face's surface at (as the rims of two faces lying nearly on one
        // another, passing a hair apart where they would cross): one vertex splits both. At
        // most a thousandth, and only on a free-form face (an edge lying along an exact one
        // is found on it as such).)
        // (Within the boundary edge's own tolerance, or where the edge comes nearest it at an
        // end of it (passing a corner of the face, not crossing an edge): the edge's point
        // only.)
        let at_end = [other.start, other.end].iter().any(|&v| self.vertices.points[v].distance(y) <= gap.max(on_tol()));
        if gap <= on_tol().max(other.tolerance) {
            return Some((u, x, None));
        }
        if at_end {
            return None;
        }
        let loose = if tangent && geom::freeform(&face.surface) { (on_tol() / sine.max(1e-6)).min(1e-3) } else { 0.0 };
        (gap <= loose.max(other.tolerance).max(within)).then_some((u, x, Some((nearest.1, s))))
    }

    /// Edges lying on a face of the other operand cross its boundary edges: splits both
    /// there.
    pub fn intersect_edges_on_faces(&mut self) {
        let mut crossings: Vec<(usize, usize, f64, f64)> = Vec::new();
        for face in &self.faces {
            // (With many edges both ways, the boundary's edges in a tree over their boxes;
            // those met taken in the boundary's order all the same.)
            let tree = (face.on_edges.len() * face.boundary.len() > 1024)
                .then(|| crate::math::Bvh::new(&face.boundary.iter().map(|b| super::graph::to_box(&self.edges[b.edge].bbox)).collect::<Vec<_>>()));
            let mut near: Vec<usize> = Vec::new();
            for &ei in &face.on_edges {
                near.clear();
                match &tree {
                    Some(tree) => {
                        near.extend(tree.overlapping(&super::graph::to_box(&self.edges[ei].bbox)));
                        near.sort_unstable();
                    }
                    None => near.extend(0..face.boundary.len()),
                }
                for &k in &near {
                    let b = &face.boundary[k];
                    let (e1, e2) = (&self.edges[ei], &self.edges[b.edge]);
                    if !geom::boxes_overlap(&e1.bbox, &e2.bbox) {
                        continue;
                    }
                    for (t, s) in edge_crossings(e1, e2, on_tol()) {
                        crossings.push((ei, b.edge, t, s));
                    }
                }
            }
        }
        for (e1, e2, t, s) in crossings {
            self.add_crossing(e1, t, e2, s);
        }
    }

    /// Every vertex lying on an edge (away from its ends) splits it, so that both
    /// operands agree on their meeting points.
    pub fn share_vertices(&mut self) {
        // The vertices as points in a tree, to find those in each edge's box.
        let points: Vec<crate::math::BoundingBox> = self.vertices.points.iter().map(|&p| crate::math::BoundingBox::new(p, p)).collect();
        let tree = crate::math::Bvh::new(&points);
        for ei in 0..self.edges.len() {
            let e = &self.edges[ei];
            let Some(curve) = &e.curve else { continue };
            let mut found = Vec::new();
            let mut near = tree.overlapping(&super::graph::to_box(&e.bbox));
            near.sort_unstable();
            for v in near {
                let p = self.vertices.points[v];
                if v == e.start || v == e.end {
                    continue;
                }
                let (t, d) = geom::project(curve, e.first, e.last, p);
                if d <= on_tol() {
                    found.push((t, v));
                }
            }
            for (t, v) in found {
                self.add_pave(ei, t, v);
            }
        }
    }

    /// Section edges on a face, where it meets faces of two other operands, cross where all
    /// three meet: splits them there (with more than two operands; with two, such points
    /// lie on an operand's edge, found from it).
    pub fn intersect_sections(&mut self) {
        let mut owners: Vec<Vec<usize>> = vec![Vec::new(); self.edges.len()];
        for (fi, face) in self.faces.iter().enumerate() {
            for &e in &face.sections {
                owners[e].push(fi);
            }
        }
        let partner = |e: usize, fi: usize| owners[e].iter().copied().find(|&g| g != fi);
        let mut crossings: Vec<(usize, usize, f64, f64)> = Vec::new();
        for (fi, face) in self.faces.iter().enumerate() {
            for (a, &e1) in face.sections.iter().enumerate() {
                for &e2 in &face.sections[a + 1..] {
                    let (Some(g1), Some(g2)) = (partner(e1, fi), partner(e2, fi)) else { continue };
                    if self.faces[g1].operand == self.faces[g2].operand {
                        continue;
                    }
                    let (x, y) = (&self.edges[e1], &self.edges[e2]);
                    if !geom::boxes_overlap(&x.bbox, &y.bbox) {
                        continue;
                    }
                    for (t, u) in edge_crossings(x, y, on_tol()) {
                        crossings.push((e1, e2, t, u));
                    }
                }
            }
        }
        for (e1, e2, t, u) in crossings {
            self.add_crossing(e1, t, e2, u);
        }
    }

    /// Section edges: where faces of the two operands meet, inside both.
    pub fn intersect_faces(&mut self) {
        let pairs = self.index().face_pairs();
        // Each pair's section curves, and where they cross each other (two ellipses of equal
        // cylinders), found for all pairs at once; then added in order.
        let found = crate::parallel::map(&pairs, |&(fa, fb)| {
            // (Where edges of either meet the other, in both boxes: points the sections pass.)
            let (x, y) = (&self.faces[fa], &self.faces[fb]);
            let crossings: Vec<Point> = x.vertices_on.iter().chain(&y.vertices_on).map(|&v| self.vertices.points[v]).filter(|&p| geom::box_contains(&x.bbox, p) && geom::box_contains(&y.bbox, p)).collect();
            // (Free-form faces near each other by their boxes, clear of each other by their
            // grids': no section, without the search.)
            if (geom::freeform(&x.surface) || geom::freeform(&y.surface)) && faces_clear(x, y) {
                return (Vec::new(), Vec::new());
            }
            let curves = sections(x, y, &crossings);
            // Each curve's box (from samples, grown by their spacing): curves whose boxes are
            // apart do not cross (nearly tangent surfaces can give many short curves).
            let boxes: Vec<crate::math::BoundingBox> = curves
                .iter()
                .map(|c| {
                    let (lo, hi) = c.range();
                    let points: Vec<Point> = (0..=64).map(|i| c.point(lo + (hi - lo) * i as f64 / 64.0)).collect();
                    let spacing = points.windows(2).map(|w| w[0].distance(w[1])).fold(0.0, f64::max);
                    let mut b = crate::math::BoundingBox::EMPTY;
                    for p in points {
                        b.add_point(p);
                    }
                    b.enlarged(spacing + on_tol())
                })
                .collect();
            // (A plane meets a cone clear of its apex in one conic, whose pieces and branches
            // do not cross: through it, in lines crossing there, which are found below.)
            let one_conic = |p: &Surface, c: &Surface| match (p, c) {
                (Surface::Plane(p), Surface::Cone(c)) => {
                    let n = Vector::from(geom::plane_normal(p));
                    let apex = c.apex();
                    n.dot(apex - p.origin()).abs() > 1e3 * on_tol() * (1.0 + (apex - Point::ORIGIN).length())
                }
                _ => false,
            };
            let conic = one_conic(&x.surface, &y.surface) || one_conic(&y.surface, &x.surface);
            let mut crossings = Vec::new();
            for i in 0..curves.len() {
                for j in i + 1..curves.len() {
                    if boxes[i].distance_to_box(&boxes[j]) > 0.0 {
                        continue;
                    }
                    if conic && matches!((&curves[i], &curves[j]), (Section::Chain(..), Section::Chain(..))) {
                        continue;
                    }
                    crossings.extend(section_crossings(&curves[i], &curves[j]));
                }
            }
            (curves, crossings)
        });
        // (The edges sections ended on, with the vertices put there.)
        let mut ended: Vec<(usize, usize)> = Vec::new();
        for (&(fa, fb), (curves, crossings)) in pairs.iter().zip(found) {
            // The crossings are vertices of both sections.
            for p in crossings {
                let v = self.vertices.add(p);
                self.faces[fa].vertices_on.push(v);
                self.faces[fb].vertices_on.push(v);
            }
            // A vertex the section passes only near (where the surfaces touch crossing, its
            // branches crossing in it) is on it if passed so twice: by both arcs the tracing
            // makes of the crossing (one alone may be another branch's, passing by).
            // (A traced section ending on an edge of either face that lies on the other face
            // too, as the edge two faces lying nearly on one another share, meets that edge
            // where the edge passes it: no edge crosses the other face there to tell.)
            let mut ends: Vec<Vec<(f64, usize, bool)>> = Vec::with_capacity(curves.len());
            for section in &curves {
                let mut met: Vec<(f64, usize, bool)> = Vec::new();
                for (edge, s, t, q) in self.section_end_meetings(fa, fb, section) {
                    // (One vertex where it meets coincident edges, or meets an edge where it
                    // has one already, or the sections cross: a vertex there lying on the
                    // edge. Found along the edge only roughly, where the surfaces touch.)
                    let e = &self.edges[edge];
                    let on_edge = |v: usize| {
                        let x = self.vertices.points[v];
                        x.distance(q) <= 100.0 * on_tol() && e.curve.as_ref().is_some_and(|c| geom::project(c, e.first, e.last, x).1 <= 10.0 * on_tol())
                    };
                    let known = met.iter().map(|m| m.1).chain(e.paves.iter().map(|p| p.1)).chain(self.faces[fa].vertices_on.iter().chain(&self.faces[fb].vertices_on).copied());
                    if let Some(v) = known.filter(|&v| on_edge(v)).min_by(|&v, &w| self.vertices.points[v].distance(q).total_cmp(&self.vertices.points[w].distance(q))) {
                        let s = geom::project(e.curve.as_ref().unwrap(), e.first, e.last, self.vertices.points[v]).0;
                        self.add_pave(edge, s, v);
                        ended.push((edge, v));
                        if !met.iter().any(|m| m.1 == v) {
                            met.push((t, v, false));
                        }
                        continue;
                    }
                    let v = self.vertices.add(q);
                    self.add_pave(edge, s, v);
                    ended.push((edge, v));
                    self.faces[fa].vertices_on.push(v);
                    self.faces[fb].vertices_on.push(v);
                    met.push((t, v, false));
                }
                ends.push(met);
            }
            let points: Vec<Vec<(f64, usize, bool)>> = curves
                .iter()
                .zip(ends)
                .map(|(section, met)| {
                    let mut found = self.section_points(fa, fb, section);
                    // (Met exactly: on it, not only passed near.)
                    for m in met {
                        match found.iter_mut().find(|f| f.1 == m.1) {
                            Some(f) => f.2 = false,
                            None => found.push(m),
                        }
                    }
                    found
                })
                .collect();
            let mut near: HashMap<usize, usize> = HashMap::default();
            for &(_, v, _) in points.iter().flatten().filter(|p| p.2) {
                *near.entry(v).or_default() += 1;
            }
            for (section, found) in curves.iter().zip(points) {
                let on = found.into_iter().filter(|&(_, v, only_near)| !only_near || near[&v] >= 2).map(|(t, v, _)| (t, v)).collect();
                self.add_section(fa, fb, section, on);
            }
        }
        for (edge, v) in ended {
            self.pave_coincident(edge, v);
        }
    }

    /// Where an open traced section, near either end, passes a boundary edge of either face
    /// lying on the other face's surface too, where the surfaces touch nearly tangent
    /// (within ten times the classifying distance): the edge, its parameter there, the
    /// section's, and the point on the edge.
    fn section_end_meetings(&self, fa: usize, fb: usize, section: &Section) -> Vec<(usize, f64, f64, Point)> {
        let Section::Sampled(points, false, traced) = section else { return Vec::new() };
        // (Between free-form faces: two lying nearly on one another as offsets and ruled
        // walls of a filling can.)
        if traced.st.is_none() || points.len() < 2 || !geom::freeform(&self.faces[fa].surface) || !geom::freeform(&self.faces[fb].surface) {
            return Vec::new();
        }
        let (lo, hi) = section.range();
        let span = hi - lo;
        let reach = 10.0 * on_tol();
        // (Where they touch tangent, a trace may stop anywhere in the band in which the
        // surfaces are within the tolerances of each other (as wide as the root of the
        // tolerance times their radius): some way short of the edge, met all the same.)
        let tangent_reach = 1e3 * on_tol();
        let mut out: Vec<(usize, f64, f64, Point)> = Vec::new();
        for (end, a, b) in [(lo, lo, lo + 0.1 * span), (hi, hi - 0.1 * span, hi)] {
            let at_end = section.point(end);
            let far = section.point(if end == lo { b } else { a }).distance(at_end);
            for (f, g) in [(fa, fb), (fb, fa)] {
                for bound in self.faces[f].boundary.iter().filter(|b| b.bounds) {
                    let e = &self.edges[bound.edge];
                    let Some(c) = &e.curve else { continue };
                    if geom::project(c, e.first, e.last, at_end).1 > far + tangent_reach {
                        continue;
                    }
                    // (Only an edge lying on the other face's surface too.)
                    let lying = (0..=4).all(|k| {
                        let x = c.point(e.first + (e.last - e.first) * k as f64 / 4.0);
                        geom::signed_distance(&self.faces[g].surface, x).abs() <= reach.max(self.faces[g].drift)
                    });
                    if !lying {
                        continue;
                    }
                    let gap = |t: f64| geom::project(c, e.first, e.last, section.point(t)).1;
                    // (Its nearest approach over the stretch: sampled, then narrowed.)
                    let n = 32;
                    let best = (0..=n).map(|i| a + (b - a) * i as f64 / n as f64).min_by(|x, y| gap(*x).total_cmp(&gap(*y))).unwrap_or(end);
                    let step = (b - a) / n as f64;
                    let (mut x, mut y) = ((best - step).max(a), (best + step).min(b));
                    for _ in 0..60 {
                        let (p, q) = (x + (y - x) / 3.0, y - (y - x) / 3.0);
                        if gap(p) < gap(q) { y = q } else { x = p }
                    }
                    let t = 0.5 * (x + y);
                    if gap(t) > tangent_reach {
                        continue;
                    }
                    let (s, _) = geom::project(c, e.first, e.last, section.point(t));
                    // (Not at the edge's own ends: those are vertices already.)
                    let q = c.point(s);
                    // (Only where the surfaces touch there, nearly tangent: elsewhere a section
                    // crossing the edge is found where the edge crosses the face.)
                    let (n1, n2) = (geom::normal_at_point(&traced.base, q), geom::normal_at_point(&traced.other, q));
                    let sine = n1.cross(n2).length() / (n1.length() * n2.length()).max(1e-300);
                    if sine >= TANGENT_SINE || (gap(t) > reach && sine > 1e-3) {
                        continue;
                    }
                    if q.distance(self.vertices.points[e.start]) <= reach || q.distance(self.vertices.points[e.end]) <= reach {
                        continue;
                    }
                    // (Once on each edge: on edges of both faces where they coincide, as the
                    // edge two faces lying nearly on one another each have.)
                    if !out.iter().any(|&(other, _, _, p)| other == bound.edge && p.distance(q) <= reach) {
                        out.push((bound.edge, s, t, q));
                    }
                }
            }
        }
        out
    }

    /// The known vertices `section` passes, with its parameter there and whether it passes
    /// only near a point where the surfaces touch crossing.
    fn section_points(&self, fa: usize, fb: usize, section: &Section) -> Vec<(f64, usize, bool)> {
        // Points of the curve known to be on either face's boundary or met by an edge.
        let mut candidates: Vec<usize> = Vec::new();
        for &f in &[fa, fb] {
            candidates.extend(&self.faces[f].vertices_on);
            for b in &self.faces[f].boundary {
                let e = &self.edges[b.edge];
                candidates.push(e.start);
                candidates.push(e.end);
                candidates.extend(e.paves.iter().map(|&(_, v)| v));
            }
        }
        // Only points in both faces' boxes can be on the section. (Kept before sorting: a
        // face with many edges has many.)
        let (ba, bb) = (&self.faces[fa].bbox, &self.faces[fb].bbox);
        candidates.retain(|&v| {
            let p = self.vertices.points[v];
            geom::box_contains(ba, p) && geom::box_contains(bb, p)
        });
        candidates.sort_unstable();
        candidates.dedup();
        let (lo, hi) = section.range();
        let period = hi - lo;
        // Only points near the curve can be on it: within its samples' box, widened by the
        // samples' spacing (the curve strays from them by less).
        let samples: Vec<Point> = match section {
            Section::Sampled(points, ..) => points.clone(),
            // A line's ends bound it (exactly: it strays nowhere from them).
            Section::Open(Curve::Line(_), ..) => vec![section.point(lo), section.point(hi)],
            _ => (0..=64).map(|i| section.point(lo + period * i as f64 / 64.0)).collect(),
        };
        let straight = matches!(section, Section::Open(Curve::Line(_), ..));
        let reach = if straight { 0.0 } else { samples.windows(2).map(|w| w[0].distance(w[1])).fold(0.0, f64::max) } + on_tol();
        let (near_lo, near_hi) = samples.iter().fold(
            ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]),
            |(lo, hi), p| ([lo[0].min(p.x), lo[1].min(p.y), lo[2].min(p.z)], [hi[0].max(p.x), hi[1].max(p.y), hi[2].max(p.z)]),
        );
        candidates.retain(|&v| {
            let p = self.vertices.points[v];
            [p.x, p.y, p.z].iter().enumerate().all(|(k, &x)| x >= near_lo[k] - reach && x <= near_hi[k] + reach)
        });
        // On a chain of quadric branches (slow to find points on), only points near a
        // polyline of it, sampled once, can be on it.
        if let Section::Chain(branches, ..) = section {
            let n = branches.len() * 256;
            let polyline: Vec<Point> = (0..=n).map(|i| section.point(hi * i as f64 / n as f64)).collect();
            let spacing = polyline.windows(2).map(|w| w[0].distance(w[1])).fold(0.0, f64::max) + on_tol();
            candidates.retain(|&v| {
                let p = self.vertices.points[v];
                polyline.windows(2).any(|w| segment_gap(w[0], w[1], p).0 <= spacing)
            });
        }
        // (A vertex where an edge meets a face may lie off the surfaces by as far as the edges
        // stray from their faces' (tolerant edges read from a file): passed by the section
        // where its point on both surfaces is.)
        let drift = self.faces[fa].drift.max(self.faces[fb].drift);
        let (a, b) = (&self.faces[fa].surface, &self.faces[fb].surface);
        // (Only of the points where an edge of one meets the other: a corner of a face near
        // the section is not on it.)
        let crossings: Vec<usize> = self.faces[fa].vertices_on.iter().chain(&self.faces[fb].vertices_on).copied().collect();
        // (A vertex where edges were taken to meet a hair apart is on the section within
        // that.)
        let exact: Vec<(f64, usize, bool)> = candidates
            .iter()
            .flat_map(|&v| {
                let p = self.vertices.points[v];
                let found = match self.loose_vertices.get(&v) {
                    Some(&allowance) => section.passes_within(p, allowance),
                    None => section.passes(p),
                };
                found.into_iter().map(move |(t, near)| (t, v, near))
            })
            .collect();
        // (Only where the section would otherwise have too few to make a piece of: elsewhere
        // the points found exactly bound it. Points within a few tolerances of each other
        // count once: the same point met by the edges of both faces, each off its face.)
        let mut distinct: Vec<Point> = Vec::new();
        for &(_, v, _) in &exact {
            let p = self.vertices.points[v];
            if !distinct.iter().any(|q| q.distance(p) <= 10.0 * drift) {
                distinct.push(p);
            }
        }
        let enough = if section.closed() { 1 } else { 2 };
        if drift <= on_tol() || distinct.len() >= enough {
            return exact;
        }
        // Points where an edge meets a face not found exactly on the section: their point on
        // both surfaces may be (an edge strays from its face by its tolerance, which where the
        // surfaces meet at a shallow angle moves the point along them by more).
        let mut found = exact;
        for &v in candidates.iter().filter(|&&v| crossings.contains(&v)) {
            if found.iter().any(|&(_, w, _)| w == v) {
                continue;
            }
            let p = self.vertices.points[v];
            let q = onto_both(a, b, p);
            let sine = geom::normal_at_point(a, q).cross(geom::normal_at_point(b, q)).length();
            // (Not where the surfaces meet nearly tangent: there edges of one lie as near the
            // other all along.)
            let crossing = sine > 0.05;
            let reach = if crossing { (drift / sine).min(10.0 * drift) } else { drift };
            // (A point found already far closer than the tolerance is the same point, met by
            // another edge.)
            if found.iter().any(|&(_, w, _)| self.vertices.points[w].distance(p) <= 0.01 * drift) {
                continue;
            }
            let mut passes = if crossing && q.distance(p) <= reach { section.passes(q) } else { Vec::new() };
            // (Or the point itself, where the section was carried on to it: a traced branch
            // ending just short of it ends on it, off the surfaces as far as it is.)
            if passes.is_empty() && crossing && q.distance(p) <= reach {
                let (t, d) = section.project(p);
                if d <= on_tol() {
                    passes.push((t, false));
                }
            }
            // (Or at an end of an open section, where the curve stops at a surface's bounds,
            // which such a vertex may lie past.)
            if passes.is_empty() && !section.closed() {
                let (lo, hi) = section.range();
                passes = [lo, hi].into_iter().filter(|&t| section.point(t).distance(p) <= reach).map(|t| (t, false)).collect();
            }
            found.extend(passes.into_iter().map(|(t, near)| (t, v, near)));
        }
        found
    }

    /// The pieces of `section` between the vertices `on` it (at their parameters) inside
    /// both faces, as section edges of both.
    /// The edges of either face lying on the other, and for each of their vertices (ends and
    /// split points) those it is on.
    fn edges_lying_on_both(&self, fa: usize, fb: usize) -> (Vec<usize>, HashMap<usize, Vec<usize>>) {
        let (a, b) = (&self.faces[fa], &self.faces[fb]);
        let lying = |f: &Face, other: &Face| {
            let on: crate::hashing::FastSet<usize> = other.on_edges.iter().copied().collect();
            f.boundary.iter().map(|u| u.edge).filter(move |e| on.contains(e)).collect::<Vec<usize>>()
        };
        let mut edges = lying(a, b);
        edges.extend(lying(b, a));
        let mut by_vertex: HashMap<usize, Vec<usize>> = HashMap::default();
        for &e in &edges {
            let edge = &self.edges[e];
            for v in [edge.start, edge.end].into_iter().chain(edge.paves.iter().map(|&(_, v)| v)) {
                by_vertex.entry(v).or_default().push(e);
            }
        }
        (edges, by_vertex)
    }

    /// Whether the piece of `section` from `t0` to `t1` (from vertex `v0` to `v1`) runs along
    /// `lying`, edges of one face lying on the other (see [`Self::edges_lying_on_both`]),
    /// within the faces' tolerance over the sine of the angle between their surfaces: points
    /// that far from where they meet are on both (where they meet nearly tangent, as a fillet
    /// the face it rolls on, far), and the section traced among them is that edge (made of
    /// it, the faces would be split twice, by two curves a little apart).
    fn along_edge_on_both(&self, fa: usize, fb: usize, section: &Section, (t0, v0): (f64, usize), (t1, v1): (f64, usize), lying: &(Vec<usize>, HashMap<usize, Vec<usize>>)) -> bool {
        let (edges, by_vertex) = lying;
        if edges.is_empty() {
            return false;
        }
        let (a, b) = (&self.faces[fa], &self.faces[fb]);
        let tolerance = a.drift.max(b.drift) + 10.0 * on_tol();
        let scale = a.bbox.0.distance(a.bbox.1).min(b.bbox.0.distance(b.bbox.1));
        // (Those through the piece's ends first: usually the one it runs along.)
        let mut first: Vec<usize> = [v0, v1].iter().filter_map(|v| by_vertex.get(v)).flatten().copied().collect();
        first.sort_unstable();
        first.dedup();
        (0..=8).all(|i| {
            let p = section.point(t0 + (t1 - t0) * i as f64 / 8.0);
            let sine = geom::normal_at_point(&a.surface, p).cross(geom::normal_at_point(&b.surface, p)).length();
            let reach = (tolerance / sine.max(1e-12)).min(1e-3 * scale).max(tolerance);
            let near = |&e: &usize| {
                let edge = &self.edges[e];
                let Some(curve) = &edge.curve else { return false };
                geom::box_contains(&(edge.bbox.0 - Vector::new(reach, reach, reach), edge.bbox.1 + Vector::new(reach, reach, reach)), p) && geom::project(curve, edge.first, edge.last, p).1 <= reach
            };
            first.iter().any(near) || edges.iter().any(near)
        })
    }

    /// Whether the piece of `section` from `t0` to `t1` runs along the boundary of either
    /// face past its surface's bound: in that face's parameters on its boundary all the way
    /// (a point past a patch's bound found on it), near the face in space, and somewhere off
    /// the surface by more than the faces' tolerance (followed on past its bound). Where the
    /// other surface grazes a face along an edge on its patch's bound, the section is that
    /// edge: made of it, it would lie on the face's own boundary arcs without being one with
    /// them (as a section along a boundary on the surface is), cutting a sliver of no area
    /// off the face.
    fn along_boundary(&self, fa: usize, fb: usize, section: &Section, (t0, v0): (f64, usize), (t1, v1): (f64, usize)) -> bool {
        let tolerance = self.faces[fa].drift.max(self.faces[fb].drift) + 10.0 * on_tol();
        [fa, fb].iter().any(|&f| {
            let face = &self.faces[f];
            // (Its ends on the face's boundary edges first.)
            let on_boundary = |v: usize| {
                face.boundary.iter().any(|u| {
                    let edge = &self.edges[u.edge];
                    edge.start == v || edge.end == v || edge.paves.iter().any(|&(_, w)| w == v)
                })
            };
            if !(on_boundary(v0) && on_boundary(v1)) {
                return false;
            }
            let (du, dv) = (face.u_range.1 - face.u_range.0, face.v_range.1 - face.v_range.0);
            let near_uv = 1e-7 * (du.abs() + dv.abs()).max(1e-300);
            let near = 1e-3 * face.bbox.0.distance(face.bbox.1);
            let gap = |q: geom::Uv| {
                face.polylines
                    .iter()
                    .flat_map(|uv| uv.windows(2))
                    .map(|w| {
                        let (a, b) = (w[0], w[1]);
                        let d = [b[0] - a[0], b[1] - a[1]];
                        let len2 = d[0] * d[0] + d[1] * d[1];
                        let s = if len2 > 0.0 { (((q[0] - a[0]) * d[0] + (q[1] - a[1]) * d[1]) / len2).clamp(0.0, 1.0) } else { 0.0 };
                        (q[0] - a[0] - d[0] * s).hypot(q[1] - a[1] - d[1] * s)
                    })
                    .fold(f64::INFINITY, f64::min)
            };
            let mut past = false;
            // (The middle first.)
            [4, 0, 1, 2, 3, 5, 6, 7, 8].iter().all(|&i| {
                let p = section.point(t0 + (t1 - t0) * i as f64 / 8.0);
                let uv = geom::uv_of(&face.surface, p);
                let off = geom::surface_point(&face.surface, uv).distance(p);
                past |= off > tolerance;
                gap(uv) <= near_uv && off <= near
            }) && past
        })
    }

    fn add_section(&mut self, fa: usize, fb: usize, section: &Section, mut on: Vec<(f64, usize)>) {
        let (lo, hi) = section.range();
        let period = hi - lo;
        on.sort_by(|a, b| a.0.total_cmp(&b.0));
        let inside_both = |graph: &Graph, p: Point| {
            graph.faces[fa].locate(p, &graph.edges, &graph.vertices) == Location::Inside
                && graph.faces[fb].locate(p, &graph.edges, &graph.vertices) == Location::Inside
        };
        let mut pieces: Vec<(f64, usize, f64, usize)> = Vec::new();
        if on.is_empty() {
            if section.closed() && inside_both(self, section.point(lo + 0.37 * period)) {
                let v = self.vertices.add(section.point(lo));
                pieces.push((lo, v, hi, v));
            }
        } else {
            for w in on.windows(2) {
                pieces.push((w[0].0, w[0].1, w[1].0, w[1].1));
            }
            if section.closed() {
                let (first, last) = (on[0], on[on.len() - 1]);
                pieces.push((last.0, last.1, first.0 + period, first.1));
            }
        }
        // (The edges either face has lying on the other, found once for all the pieces.)
        let lying = self.edges_lying_on_both(fa, fb);
        for (t0, v0, t1, v1) in pieces {
            if t1 - t0 <= 1e-12 * period.max(1.0) {
                continue;
            }
            if self.along_edge_on_both(fa, fb, section, (t0, v0), (t1, v1), &lying) {
                continue;
            }
            // (A marched section, followed a little past the surfaces' bounds: there off them,
            // not inside the faces even where locating it on them puts it so.)
            let marched = matches!(section, Section::Sampled(_, _, traced) if traced.st.is_some());
            let near = 10.0 * (self.faces[fa].drift.max(self.faces[fb].drift) + on_tol());
            let off_surfaces = |graph: &Graph, p: Point| marched && [fa, fb].iter().any(|&f| geom::signed_distance(&graph.faces[f].surface, p).abs() > near);
            // (Nor where it was followed past either surface's bounds: there the surfaces do
            // not meet. Unless only as far as the faces' tolerance: across a corner of a face
            // whose edges stray by that much.)
            if let Section::Sampled(_, _, traced) = section {
                let (a, b) = (self.vertices.points[v0], self.vertices.points[v1]);
                if a.distance(b) > near && [0.25, 0.5, 0.75].iter().all(|&f| traced.beyond_bounds(t0 + (t1 - t0) * f)) {
                    continue;
                }
            }
            let middle = section.point(0.5 * (t0 + t1));
            if off_surfaces(self, middle) || !inside_both(self, middle) {
                continue;
            }
            if marched && self.along_boundary(fa, fb, section, (t0, v0), (t1, v1)) {
                continue;
            }
            let (a, b) = (self.vertices.points[v0], self.vertices.points[v1]);
            let (curve, first, last) = section.edge_curve(t0, t1, a, b);
            let edge = Edge::new(Some(curve), first, last, v0, v1, &self.vertices, None);
            self.edges.push(edge);
            let i = self.edges.len() - 1;
            self.faces[fa].sections.push(i);
            self.faces[fb].sections.push(i);
        }
    }
}

/// The point on both surfaces nearest `p` (Newton steps along their normals).
fn onto_both(a: &Surface, b: &Surface, mut p: Point) -> Point {
    for _ in 0..20 {
        // (On a free-form surface the distance and its normal from one foot: that normal is
        // the distance's gradient, the way it grows.)
        let ((da, foot_na), (db, foot_nb)) = (geom::signed_distance_and_normal(a, p), geom::signed_distance_and_normal(b, p));
        if da.abs().max(db.abs()) <= 1e-12 {
            break;
        }
        // (The normals as the distances grow along them.)
        let toward = |s: &Surface, d: f64| {
            let n = geom::normal_at_point(s, p);
            if geom::signed_distance(s, p + n * 1e-6) > d { n } else { -n }
        };
        let na = foot_na.unwrap_or_else(|| toward(a, da));
        let nb = foot_nb.unwrap_or_else(|| toward(b, db));
        let (m11, m12, m22) = (na.dot(na), na.dot(nb), nb.dot(nb));
        let det = m11 * m22 - m12 * m12;
        if det.abs() <= 1e-18 {
            break;
        }
        let (alpha, beta) = ((-da * m22 + db * m12) / det, (-db * m11 + da * m12) / det);
        p = p + na * alpha + nb * beta;
    }
    p
}

/// True if two curves are the same circle, or lie on the same line.
fn same_carrier(a: &Curve, b: &Curve) -> bool {
    let tol = on_tol();
    match (a, b) {
        (Curve::Line(x), Curve::Line(y)) => {
            let dx = Vector::from(x.direction);
            let d = y.origin - x.origin;
            dx.cross(Vector::from(y.direction)).length() <= 1e-9 && (d - dx * d.dot(dx)).length() <= tol
        }
        (Curve::Circle(x), Curve::Circle(y)) => {
            let axis = |c: &Circle| Vector::from(c.frame().x()).cross(Vector::from(c.frame().y()));
            x.center().distance(y.center()) <= tol && (x.radius() - y.radius()).abs() <= tol && axis(x).cross(axis(y)).length() <= 1e-9
        }
        _ => false,
    }
}

/// The smallest distance between segments `x0 x1` and `y0 y1`.
fn segments_gap(x0: Point, x1: Point, y0: Point, y1: Point) -> f64 {
    let (d1, d2, r) = (x1 - x0, y1 - y0, x0 - y0);
    let (a, e, f) = (d1.dot(d1), d2.dot(d2), d2.dot(r));
    let (s, t) = if a <= 1e-300 && e <= 1e-300 {
        (0.0, 0.0)
    } else if a <= 1e-300 {
        (0.0, (f / e).clamp(0.0, 1.0))
    } else {
        let c = d1.dot(r);
        if e <= 1e-300 {
            ((-c / a).clamp(0.0, 1.0), 0.0)
        } else {
            let b = d1.dot(d2);
            let denom = a * e - b * b;
            let mut s = if denom > 0.0 { ((b * f - c * e) / denom).clamp(0.0, 1.0) } else { 0.0 };
            let mut t = (b * s + f) / e;
            if t < 0.0 {
                t = 0.0;
                s = (-c / a).clamp(0.0, 1.0);
            } else if t > 1.0 {
                t = 1.0;
                s = ((b - c) / a).clamp(0.0, 1.0);
            }
            (s, t)
        }
    };
    (x0 + d1 * s).distance(y0 + d2 * t)
}

/// Whether two sections on circles (whole or arcs) cannot meet within the tolerance: a
/// point on one within it of the other lies (nearly) on both circles' spheres and planes, so
/// the spheres meet and each circle reaches the other's plane.
fn circles_apart(a: &Section, b: &Section) -> bool {
    let circle = |s: &Section| match s {
        Section::Open(Curve::Circle(c), ..) | Section::Closed(Curve::Circle(c), _) => Some(c.clone()),
        _ => None,
    };
    let (Some(p), Some(q)) = (circle(a), circle(b)) else { return false };
    let normal = |c: &Circle| {
        let f = c.frame();
        Vector::from(f.x()).cross(Vector::from(f.y()))
    };
    let (np, nq) = (normal(&p), normal(&q));
    let between = q.center() - p.center();
    let d = between.length();
    let (rp, rq) = (p.radius(), q.radius());
    // (Twice the tolerance, and a little for rounding.)
    let slack = 2.0 * on_tol() + 1e-9 * (1.0 + d + rp + rq);
    if d > rp + rq + slack || d < (rp - rq).abs() - slack {
        return true;
    }
    // How far a circle reaches from its centre along the other's normal.
    let reach = |r: f64, n: Vector, other: Vector| r * (1.0 - n.dot(other).powi(2)).max(0.0).sqrt();
    between.dot(nq).abs() > reach(rp, np, nq) + slack || between.dot(np).abs() > reach(rq, nq, np) + slack
}

/// Points where two section curves meet.
fn section_crossings(a: &Section, b: &Section) -> Vec<Point> {
    if circles_apart(a, b) {
        return Vec::new();
    }
    if let (Section::Sampled(pa, ca, _), Section::Sampled(pb, cb, _)) = (a, b) {
        // Two traced curves: where their polylines come close, refined near there.
        let segments = |p: &Vec<Point>, closed: bool| if closed { p.len() } else { p.len() - 1 };
        let (na, nb) = (segments(pa, *ca), segments(pb, *cb));
        let length = |p: &Vec<Point>, i: usize| p[i].distance(p[(i + 1) % p.len()]);
        let mut out: Vec<Point> = Vec::new();
        // Refined from `(ta, tb)` (stretches `i` of `a` and `j` of `b`) by alternate projections.
        let refine = |i: usize, j: usize, ta: f64, tb: f64, out: &mut Vec<Point>| {
            let (mut q, mut sa, mut sb) = (a.point(ta), ta, tb);
            // (Done where a round leaves everything as it was: so would the next.)
            let state = |q: Point, sa: f64, sb: f64| [q.x, q.y, q.z, sa, sb].map(f64::to_bits);
            // (Also done where the gap between the curves, shrinking as it does, would not
            // come within the tolerance in the rounds left: curves running side by side
            // without crossing.)
            const ROUNDS: usize = 8;
            let mut gap = f64::INFINITY;
            for round in 0..ROUNDS {
                let before = state(q, sa, sb);
                (sb, _) = b.project_near(q, sb);
                let qb = b.point(sb);
                (sa, _) = a.project_near(qb, sa);
                q = a.point(sa);
                if state(q, sa, sb) == before {
                    break;
                }
                let now = q.distance(qb);
                let ratio = now / gap;
                let left = (ROUNDS - 1 - round) as i32;
                if now > on_tol() && (ratio >= 1.0 || now * ratio.powi(left) > on_tol()) {
                    break;
                }
                gap = now;
            }
            // Only crossings within this pair of stretches.
            let within = |s: f64, k: usize| s >= k as f64 - 1e-9 && s <= k as f64 + 1.0 + 1e-9;
            if within(sa, i) && within(sb, j) && b.point(sb).distance(q) <= on_tol() && !out.iter().any(|x| x.distance(q) <= on_tol()) {
                out.push(q);
            }
        };
        if let (Section::Sampled(_, _, ta), Section::Sampled(_, _, tb)) = (a, b) {
            if ta.base == tb.base && ta.uv.len() == pa.len() && tb.uv.len() == pb.len() {
                // Traced on one surface: where their polylines cross in its parameters
                // (curves running side by side there do not cross).
                let period = |x: Option<f64>| x.filter(|p| *p > 0.0);
                let (pu, pv) = (period(ta.u_period.or(tb.u_period)), period(ta.v_period.or(tb.v_period)));
                // `b`'s point moved by whole periods to be nearest `near`.
                let shifted = |p: (f64, f64), near: (f64, f64)| {
                    let by = |x: f64, to: f64, period: Option<f64>| period.map_or(x, |p| x + ((to - x) / p).round() * p);
                    (by(p.0, near.0, pu), by(p.1, near.1, pv))
                };
                for i in 0..na {
                    let (a0, a1) = (ta.uv[i], ta.uv[(i + 1) % ta.uv.len()]);
                    let a1 = shifted(a1, a0);
                    for j in 0..nb {
                        let b0 = shifted(tb.uv[j], a0);
                        let b1 = shifted(tb.uv[(j + 1) % tb.uv.len()], b0);
                        // Where a0 + s (a1 - a0) = b0 + r (b1 - b0), both within [0, 1].
                        let (dx, dy, ex, ey) = (a1.0 - a0.0, a1.1 - a0.1, b1.0 - b0.0, b1.1 - b0.1);
                        let den = dx * ey - dy * ex;
                        if den == 0.0 {
                            continue;
                        }
                        let (fx, fy) = (b0.0 - a0.0, b0.1 - a0.1);
                        let (s, r) = ((fx * ey - fy * ex) / den, (fx * dy - fy * dx) / den);
                        let slack = 1e-9;
                        if s < -slack || s > 1.0 + slack || r < -slack || r > 1.0 + slack {
                            continue;
                        }
                        refine(i, j, i as f64 + s.clamp(0.0, 1.0), j as f64 + r.clamp(0.0, 1.0), &mut out);
                    }
                }
                return out;
            }
        }
        for i in 0..na {
            for j in 0..nb {
                let (x0, x1, y0, y1) = (pa[i], pa[(i + 1) % pa.len()], pb[j], pb[(j + 1) % pb.len()]);
                let reach = length(pa, i).max(length(pb, j));
                // Quick rejection by the segments' middles, then by how close they come (the
                // polylines of crossing curves cross within their sag).
                if (x0 + (x1 - x0) * 0.5).distance(y0 + (y1 - y0) * 0.5) > 2.0 * reach + on_tol() {
                    continue;
                }
                if segments_gap(x0, x1, y0, y1) > 0.05 * reach + on_tol() {
                    continue;
                }
                refine(i, j, i as f64 + 0.5, j as f64 + 0.5, &mut out);
            }
        }
        return out;
    }
    if let Section::Sampled(points, closed, _) = b {
        // Against a traced curve: compare with its polyline first.
        let n = points.len();
        let segments = if *closed { n } else { n - 1 };
        let near_b = |q: Point| (0..segments).map(|i| segment_gap(points[i], points[(i + 1) % n], q).0).fold(f64::INFINITY, f64::min);
        let (lo, hi) = a.range();
        let m = 720;
        let step = (hi - lo) / m as f64;
        let d: Vec<f64> = (0..=m).map(|i| near_b(a.point(lo + step * i as f64))).collect();
        let reach = (0..segments).map(|i| points[i].distance(points[(i + 1) % n])).fold(0.0, f64::max);
        let mut out: Vec<Point> = Vec::new();
        for i in 1..m {
            if d[i] > d[i - 1] || d[i] > d[i + 1] || d[i] > reach {
                continue;
            }
            // Alternate projections between the curves.
            let mut q = a.point(lo + step * i as f64);
            for _ in 0..30 {
                q = b.point(b.project(q).0);
                q = a.point(a.project(q).0);
            }
            if b.project(q).1 <= on_tol() && !out.iter().any(|x| x.distance(q) <= on_tol()) {
                out.push(q);
            }
        }
        return out;
    }
    if let Section::Sampled(..) = a {
        return section_crossings(b, a);
    }
    // A line or circle projects onto in closed form: the distances from samples of the
    // other to it, refined at their minima.
    let exact = |s: &Section| matches!(s, Section::Open(..) | Section::Closed(..));
    if exact(a) && !exact(b) {
        return section_crossings(b, a);
    }
    if exact(b) {
        let (lo, hi) = a.range();
        let n = 720;
        let step = (hi - lo) / n as f64;
        // Distances from samples of `a`, with the feet on `b`.
        let samples: Vec<(f64, f64, f64)> = (0..=n + 1)
            .map(|i| {
                let t = lo + (hi - lo) * i as f64 / n as f64;
                let (s, d) = b.project(a.point(t));
                (t, s, d)
            })
            .collect();
        // Local minima of the distance; where the curves (nearly) coincide there are many,
        // and the deepest are enough.
        let mut minima: Vec<usize> = (1..samples.len() - 1)
            .filter(|&i| samples[i].2 < samples[i - 1].2 && samples[i].2 <= samples[i + 1].2)
            .collect();
        minima.sort_by(|&x, &y| samples[x].2.total_cmp(&samples[y].2));
        minima.truncate(64);
        let mut out: Vec<Point> = Vec::new();
        for i in minima {
            let (t, seed, _) = samples[i];
            // The distance to `b` near the sample's foot (Newton from it).
            let gap = |x: f64| b.project_seeded(a.point(x), seed).1;
            let (mut x, mut y) = (t - step, t + step);
            let ratio = 0.618_033_988_749_895;
            for _ in 0..80 {
                let (x1, x2) = (y - ratio * (y - x), x + ratio * (y - x));
                if gap(x1) < gap(x2) {
                    y = x2;
                } else {
                    x = x1;
                }
            }
            let p = a.point(0.5 * (x + y));
            if b.project(p).1 <= on_tol() && !out.iter().any(|q| q.distance(p) <= on_tol()) {
                out.push(p);
            }
        }
        return out;
    }
    // Two chains (slow to project onto): distances from samples of `a` to a polyline of
    // `b`, then, where they are least, the nearest points of the curves from there (see
    // `crossing_near`).
    let (blo, bhi) = b.range();
    let m = 1024;
    let bt: Vec<f64> = (0..=m).map(|j| blo + (bhi - blo) * j as f64 / m as f64).collect();
    let bp: Vec<Point> = bt.iter().map(|&t| b.point(t)).collect();
    let reach_b = bp.windows(2).map(|w| w[0].distance(w[1])).fold(0.0, f64::max);
    let (lo, hi) = a.range();
    let n = 720;
    let ta: Vec<f64> = (0..=n).map(|i| lo + (hi - lo) * i as f64 / n as f64).collect();
    let pa: Vec<Point> = ta.iter().map(|&t| a.point(t)).collect();
    let reach = reach_b + pa.windows(2).map(|w| w[0].distance(w[1])).fold(0.0, f64::max);
    // Polylines whose boxes (grown by the spacing of their points) are apart do not meet.
    let bbox = |ps: &[Point]| {
        let mut b = crate::math::BoundingBox::EMPTY;
        for &p in ps {
            b.add_point(p);
        }
        b
    };
    if bbox(&pa).enlarged(reach + on_tol()).distance_to_box(&bbox(&bp)) > 0.0 {
        return Vec::new();
    }
    // The nearest point of `b`'s polyline: blocks of segments with their boxes, skipped when
    // farther than the best so far.
    const BLOCK: usize = 32;
    let blocks: Vec<(usize, crate::math::BoundingBox)> = (0..m).step_by(BLOCK).map(|j| (j, bbox(&bp[j..(j + BLOCK + 1).min(m + 1)]))).collect();
    // (From the segment nearest the sample before, whose foot is near: a small best so far
    // skips most blocks.)
    let hint = std::cell::Cell::new(0usize);
    let near_b = |q: Point| {
        let mut best = (f64::INFINITY, blo);
        let mut at = hint.get();
        for j in hint.get().saturating_sub(2)..(hint.get() + 3).min(m) {
            let (d, s) = segment_gap(bp[j], bp[j + 1], q);
            if d < best.0 {
                (best, at) = ((d, bt[j] + s * (bt[j + 1] - bt[j])), j);
            }
        }
        for &(start, ref bx) in &blocks {
            if bx.distance_to_point(q) >= best.0 {
                continue;
            }
            for j in start..(start + BLOCK).min(m) {
                let (d, s) = segment_gap(bp[j], bp[j + 1], q);
                if d < best.0 {
                    (best, at) = ((d, bt[j] + s * (bt[j + 1] - bt[j])), j);
                }
            }
        }
        hint.set(at);
        best
    };
    let gaps: Vec<(f64, f64)> = pa.iter().map(|&q| near_b(q)).collect();
    // Local minima of the gap (where the curves nearly coincide there are many: the deepest
    // are enough).
    let mut minima: Vec<usize> = (0..=n)
        .filter(|&i| {
            let d = gaps[i].0;
            let (before, after) = (if i > 0 { gaps[i - 1].0 } else { f64::INFINITY }, if i < n { gaps[i + 1].0 } else { f64::INFINITY });
            d <= before && d <= after && d <= reach
        })
        .collect();
    minima.sort_by(|&x, &y| gaps[x].0.total_cmp(&gaps[y].0));
    minima.truncate(64);
    minima.sort_unstable();
    let mut out: Vec<Point> = Vec::new();
    for i in minima {
        let (sa, sb) = crossing_near(a, b, ta[i], gaps[i].1, 2.0 * (hi - lo).abs() / n as f64, 2.0 * (bhi - blo).abs() / m as f64);
        let q = a.point(sa);
        // (A point found already first: the projection onto all of `b` costs most.)
        if b.point(sb).distance(q) <= on_tol() && !out.iter().any(|x| x.distance(q) <= on_tol()) && b.project(q).1 <= on_tol() {
            out.push(q);
        }
    }
    out
}

/// Where curves `a` and `b` come nearest from `(sa, sb)` on them: Gauss–Newton on both
/// parameters at once (damped where a step would not bring them nearer), which converges
/// fast where they cross at however small an angle; alternate projections between them,
/// converging slowly there, where it stalls short or strays farther than `within_a` and
/// `within_b` (in their parameters) from where it started, nearer another crossing.
fn crossing_near(a: &Section, b: &Section, start_a: f64, start_b: f64, within_a: f64, within_b: f64) -> (f64, f64) {
    let (mut sa, mut sb) = (start_a, start_b);
    let ((alo, ahi), (blo, bhi)) = (a.range(), b.range());
    let clamp = |s: f64, lo: f64, hi: f64| s.clamp(lo.min(hi), lo.max(hi));
    let tangent = |c: &Section, s: f64, lo: f64, hi: f64| {
        let h = 1e-7 * (hi - lo).abs().max(1e-12);
        let (x, y) = (clamp(s - h, lo, hi), clamp(s + h, lo, hi));
        (c.point(y) - c.point(x)) * (1.0 / (y - x).max(1e-300))
    };
    let mut gap = a.point(sa).distance(b.point(sb));
    let mut damping = 0.0;
    for _ in 0..40 {
        if gap <= 1e-3 * on_tol() {
            break;
        }
        let (da, db) = (tangent(a, sa, alo, ahi), tangent(b, sb, blo, bhi));
        let r = a.point(sa) - b.point(sb);
        // (J = [da, -db]: Jᵀ J + λ diag, Jᵀ r.)
        let (m11, m12, m22) = (da.dot(da), -da.dot(db), db.dot(db));
        let (g1, g2) = (da.dot(r), -db.dot(r));
        let mut stepped = false;
        for _ in 0..8 {
            let (k11, k22) = (m11 * (1.0 + damping), m22 * (1.0 + damping));
            let det = k11 * k22 - m12 * m12;
            if !(det.abs() > 1e-300) {
                damping = damping.max(1e-6) * 10.0;
                continue;
            }
            let (x1, x2) = ((k22 * g1 - m12 * g2) / det, (k11 * g2 - m12 * g1) / det);
            let (na, nb) = (clamp(sa - x1, alo, ahi), clamp(sb - x2, blo, bhi));
            let now = a.point(na).distance(b.point(nb));
            if now < gap {
                (sa, sb, gap) = (na, nb, now);
                damping *= 0.1;
                stepped = true;
                break;
            }
            damping = damping.max(1e-6) * 10.0;
        }
        if !stepped {
            break;
        }
    }
    // (Where it stalls short of meeting, or strays, alternate projections from the start.)
    if gap > on_tol() || (sa - start_a).abs() > within_a || (sb - start_b).abs() > within_b {
        (sa, sb) = (start_a, start_b);
        let mut q = a.point(sa);
        for _ in 0..60 {
            let previous = q;
            (sb, _) = b.project_near(q, sb);
            (sa, _) = a.project_near(b.point(sb), sa);
            q = a.point(sa);
            if q.distance(previous) <= 1e-12 {
                break;
            }
        }
    }
    (sa, sb)
}

/// Points where two edges cross (pass within `tolerance`) away from their ends: `(t on e1,
/// s on e2)`.
fn edge_crossings(e1: &Edge, e2: &Edge, tolerance: f64) -> Vec<(f64, f64)> {
    let (Some(c1), Some(c2)) = (&e1.curve, &e2.curve) else { return Vec::new() };
    // On one circle or line (pieces of coincident edges of the two operands) the edges
    // cross nowhere: their ends split each other as vertices lying on them.
    if same_carrier(c1, c2) {
        return Vec::new();
    }
    // The edges' samples (kept with them: an edge meets many others).
    let (k1, k2) = (e1.chords.get_or_init(|| Chords::of(c1, e1.first, e1.last)), e2.chords.get_or_init(|| Chords::of(c2, e2.first, e2.last)));
    let (t1, p1, t2, p2) = (&k1.params, &k1.points, &k2.params, &k2.points);
    let mut out: Vec<(f64, f64)> = Vec::new();
    // Chords' boxes widened by half their length (the reach below): pairs whose boxes miss
    // each other are passed over, a run of the second edge's at once where its box misses.
    let misses = |(lo1, hi1): &([f64; 3], [f64; 3]), (lo2, hi2): &([f64; 3], [f64; 3])| (0..3).any(|k| lo1[k] > hi2[k] || lo2[k] > hi1[k]);
    for i in 0..p1.len() - 1 {
        let box1 = k1.boxes[i];
        for (run, run_box) in k2.runs.iter().enumerate() {
            if misses(&box1, run_box) {
                continue;
            }
            for j in run * Chords::RUN..((run + 1) * Chords::RUN).min(p2.len() - 1) {
                if misses(&box1, &k2.boxes[j]) {
                    continue;
                }
                // Closest points of the two chords.
                let (a, b, c, d) = (p1[i], p1[i + 1], p2[j], p2[j + 1]);
                let (u, v, w) = (b - a, d - c, a - c);
                let (aa, bb, cc, dd, ee) = (u.dot(u), u.dot(v), v.dot(v), u.dot(w), v.dot(w));
                let den = aa * cc - bb * bb;
                if den <= 1e-24 * aa * cc {
                    continue;
                }
                let s = ((bb * ee - cc * dd) / den).clamp(0.0, 1.0);
                let r = ((aa * ee - bb * dd) / den).clamp(0.0, 1.0);
                let gap = (a + u * s).distance(c + v * r);
                // (Curves stray from their chords by less than half a chord; a line's chord is
                // the line.)
                let straight = |c: &Curve| matches!(c, Curve::Line(_));
                let half = |c: &Curve, chord: Vector| if straight(c) { 0.0 } else { 0.5 * chord.length() };
                let reach = if straight(c1) || straight(c2) { half(c1, u) + half(c2, v) + tolerance } else { 0.5 * (u.length() + v.length()) };
                if gap > reach {
                    continue;
                }
                // Newton on |c1(t) - c2(s)|².
                let (mut t, mut q) = (t1[i] + (t1[i + 1] - t1[i]) * s, t2[j] + (t2[j + 1] - t2[j]) * r);
                for _ in 0..40 {
                    let (x, dx) = c1.d1(t);
                    let (y, dy) = c2.d1(q);
                    let f = x - y;
                    // Gauss-Newton: [dx -dy] [δt δq]ᵀ = -f.
                    let (m11, m12, m22) = (dx.dot(dx), -dx.dot(dy), dy.dot(dy));
                    let (r1, r2) = (-f.dot(dx), f.dot(dy));
                    let det = m11 * m22 - m12 * m12;
                    if det.abs() <= 1e-300 {
                        break;
                    }
                    let (dt, dq) = ((r1 * m22 - m12 * r2) / det, (m11 * r2 - m12 * r1) / det);
                    t += dt;
                    q += dq;
                    // Converged: the steps are down to rounding.
                    if dt.abs() <= 1e-15 * (1.0 + t.abs()) && dq.abs() <= 1e-15 * (1.0 + q.abs()) {
                        break;
                    }
                }
                let inside = |x: f64, a: f64, b: f64| x > a.min(b) + 1e-9 && x < a.max(b) - 1e-9;
                if !inside(t, e1.first, e1.last) || !inside(q, e2.first, e2.last) {
                    continue;
                }
                // Where a free edge, or a free-form one (a traced section's curve, shared by faces
                // of both operands), coincides with the other (a common stretch), they do not
                // cross: the stretch's ends split the edges as vertices lying on them.
                let step = 1e-3 * (e1.last - e1.first).abs();
                let along = [t - step, t + step].into_iter().filter(|&x| inside(x, e1.first, e1.last));
                let free_form = |c: &Curve| matches!(c, Curve::BSpline(_) | Curve::Bezier(_));
                if (e1.free || e2.free || (free_form(c1) && free_form(c2)))
                    && along.clone().count() > 0
                    && along.into_iter().all(|x| geom::project(c2, e2.first, e2.last, c1.point(x)).1 <= on_tol())
                {
                    continue;
                }
                if c1.point(t).distance(c2.point(q)) <= tolerance && !out.iter().any(|&(x, _)| (x - t).abs() <= 1e-9) {
                    out.push((t, q));
                }
            }
        }
    }
    out
}
