//! Where a line meets a face's surface, over the face's parameter box: the ray-casting core.
//!
//! Planes, cylinders, cones, spheres and tori are solved in closed form: the line, written
//! in the surface's own frame (so that a surface far from the world origin is solved as
//! accurately as one at it), substituted into the surface's implicit equation gives a
//! linear, quadratic or (for a torus, with the line re-based at its point nearest the
//! centre) quartic equation in the line parameter.
//!
//! Any other surface is cut, once per face, into small patches (along a B-spline's knot
//! spans first, then halved until each patch is nearly flat and its normals nearly
//! parallel), each bounded by a box that holds the whole patch; the boxes go into a bounding
//! volume tree. A line is tried only on the patches whose boxes it crosses: from where it
//! meets the patch's two corner triangles, and from the patch's centre, Newton steps solve
//! `S(u, v) = L(w)` (three equations in `u`, `v`, `w`). A line nearly tangent to a patch
//! has the patch cut further on the fly, so that two crossings close together each get a
//! start of their own.

use crate::math::{BoundingBox, Bvh, Point, Vector};

use super::adaptor::in_period;
use super::int_curve_plane::Transition;
use super::parametric::ParametricSurface;
use super::polynomial_roots::PolynomialRoots;
use super::surface::Surface;

/// The bounds of the infinite line.
pub(crate) const INFINITE: f64 = 2.0e100;

/// A point where the line meets the surface.
#[derive(Clone, Copy, Debug)]
pub(crate) struct LinePoint {
    /// The line's point, `origin + w·direction`.
    pub point: Point,
    /// The surface parameters there.
    pub u: f64,
    pub v: f64,
    /// The line parameter.
    pub w: f64,
    /// How the line crosses the surface, against its normal `S_u × S_v` (`In`) or along it.
    pub transition: Transition,
}

/// Solutions closer than this in `u`, `v` and `w` (and crossing alike) are one.
const SAME: f64 = 1e-8;

/// How far outside the face's box a solution may lie and be kept.
const BOX_SLACK: f64 = 1e-8;

/// The most patches a face is cut into ahead of the lines.
const MAX_PATCHES: usize = 4096;

/// Halvings of a starting patch at most (ahead of the lines).
const MAX_DEPTH: usize = 6;

/// Further halvings of a patch a line nearly touches.
const MAX_TANGENT_DEPTH: usize = 5;

/// A face's surface made ready for lines.
pub(crate) struct LineFace {
    surface: Surface,
    bounds: (f64, f64, f64, f64),
    method: Method,
}

enum Method {
    /// A plane, cylinder, cone, sphere or torus (the basis of the surface).
    ClosedForm(Surface),
    Sampled(Sampled),
    /// A box without extent in one direction.
    NotReady,
}

/// The patches of a free-form face.
struct Sampled {
    patches: Vec<Patch>,
    tree: Bvh,
    /// The size of the face (for tolerances on residuals).
    size: f64,
}

/// A parameter rectangle of the surface, its box, corners and the range of its normals.
#[derive(Clone, Debug)]
struct Patch {
    u0: f64,
    u1: f64,
    v0: f64,
    v1: f64,
    bbox: BoundingBox,
    /// The corners at `(u0, v0)`, `(u1, v0)`, `(u1, v1)`, `(u0, v1)`.
    corners: [Point; 4],
    /// The unit normal at the centre (zero where degenerate) and the largest angle (radians)
    /// between it and the normals at the corners.
    normal: Vector,
    spread: f64,
}

fn untrimmed(surface: &Surface) -> &Surface {
    match surface {
        Surface::Trimmed(t) => untrimmed(t.basis()),
        s => s,
    }
}

fn infinite(x: f64) -> bool {
    !(x.abs() < 1e100)
}

impl LineFace {
    /// The face's `surface` (placed) over its parameter box `bounds = (u0, u1, v0, v1)`.
    pub(crate) fn new(surface: Surface, bounds: (f64, f64, f64, f64)) -> LineFace {
        let basis = untrimmed(&surface);
        if matches!(basis, Surface::Plane(_) | Surface::Cylinder(_) | Surface::Cone(_) | Surface::Sphere(_) | Surface::Torus(_)) {
            let basis = basis.clone();
            return LineFace { surface, bounds, method: Method::ClosedForm(basis) };
        }
        // Infinite directions (an extrusion's) are cut to a large window.
        let finite = |a: f64, b: f64| if infinite(a) || infinite(b) { (a.max(-1e3), b.min(1e3)) } else { (a, b) };
        let ((u0, u1), (v0, v1)) = (finite(bounds.0, bounds.1), finite(bounds.2, bounds.3));
        let box_ = (u0, u1, v0, v1);
        let extent = |width: f64, res: f64| if res > 0.0 && res.is_finite() { width / res } else { f64::INFINITY };
        let ready = extent(u1 - u0, super::resolution::u_resolution(&surface, box_, 1.0)) > 1e-7
            && extent(v1 - v0, super::resolution::v_resolution(&surface, box_, 1.0)) > 1e-7
            && u1 > u0
            && v1 > v0;
        if !ready {
            return LineFace { surface, bounds, method: Method::NotReady };
        }
        let sampled = sample(&surface, box_);
        LineFace { surface, bounds, method: Method::Sampled(sampled) }
    }

    /// Every point where the line `origin + w·direction` meets the surface within the face's
    /// box (and `w` in `[par_min, par_max]`; points beyond may be listed too), with the range
    /// of `w` to keep. `None` for a face without extent along one of its directions.
    pub(crate) fn perform(&self, origin: Point, direction: Vector, par_min: f64, par_max: f64) -> Option<(Vec<LinePoint>, f64, f64)> {
        let mut found: Vec<LinePoint> = Vec::new();
        if !(par_min <= par_max) || !(direction.length_squared() > 0.0) {
            return match self.method {
                Method::NotReady => None,
                _ => Some((found, par_min, par_max)),
            };
        }
        match &self.method {
            Method::NotReady => return None,
            Method::ClosedForm(basis) => {
                for w in closed_form_roots(basis, origin, direction) {
                    let point = origin + direction * w;
                    let Some((u, v)) = basis.parameters(point) else { continue };
                    self.keep(&mut found, point, u, v, w, direction, false);
                }
            }
            Method::Sampled(sampled) => self.sampled(sampled, origin, direction, par_min, par_max, &mut found),
        }
        Some((found, par_min, par_max))
    }

    /// Adds the solution if it lies in the box (periodic parameters normalised) and is not
    /// one found already (`by_point`: also one at the same line parameter and point, as on a
    /// seam or pole of a sampled surface).
    fn keep(&self, found: &mut Vec<LinePoint>, point: Point, u: f64, v: f64, w: f64, d: Vector, by_point: bool) {
        let s = &self.surface;
        let (u0, u1, v0, v1) = self.bounds;
        let angular_u = matches!(untrimmed(s), Surface::Cylinder(_) | Surface::Cone(_) | Surface::Sphere(_));
        let u = if s.is_u_periodic() || angular_u {
            let p = s.u_period().unwrap_or(std::f64::consts::TAU);
            if infinite(u0) { u } else { in_period(u, u0, u0 + p) }
        } else {
            u
        };
        let v = if s.is_v_periodic() {
            let p = s.v_period().unwrap_or(std::f64::consts::TAU);
            if infinite(v0) { v } else { in_period(v, v0, v0 + p) }
        } else {
            v
        };
        if !(u.is_finite() && v.is_finite() && w.is_finite()) {
            return;
        }
        if u0 - u >= BOX_SLACK || u - u1 >= BOX_SLACK || v0 - v >= BOX_SLACK || v - v1 >= BOX_SLACK {
            return;
        }
        let transition = crossing(s, u, v, d);
        let same = |p: &LinePoint| {
            p.transition == transition
                && (p.w - w).abs() < SAME
                && (((p.u - u).abs() < SAME && (p.v - v).abs() < SAME) || (by_point && p.point.distance(point) < SAME))
        };
        if !found.iter().any(same) {
            found.push(LinePoint { point, u, v, w, transition });
        }
    }

    fn sampled(&self, sampled: &Sampled, origin: Point, d: Vector, par_min: f64, par_max: f64, found: &mut Vec<LinePoint>) {
        // The stretch of the line within the face's box.
        let Some((ta, tb)) = line_in_box(&sampled.tree.bbox(), origin, d) else { return };
        let (ta, tb) = (ta.max(par_min), tb.min(par_max));
        if !(ta <= tb) {
            return;
        }
        let start = origin + d * ta;
        let mut solver = Solver { surface: &self.surface, origin, d, bounds: self.patch_bounds(sampled), size: sampled.size };
        for i in sampled.tree.along_ray(start, d) {
            let patch = &sampled.patches[i];
            let Some((pa, _)) = line_in_box(&patch.bbox, origin, d) else { continue };
            if pa > tb {
                continue;
            }
            solver.patch(patch, 0, &mut |u, v, w| {
                let point = origin + d * w;
                self.keep(found, point, u, v, w, d, true);
            });
        }
    }

    fn patch_bounds(&self, sampled: &Sampled) -> (f64, f64, f64, f64) {
        let (mut u0, mut u1, mut v0, mut v1) = (f64::INFINITY, f64::NEG_INFINITY, f64::INFINITY, f64::NEG_INFINITY);
        for p in &sampled.patches {
            u0 = u0.min(p.u0);
            u1 = u1.max(p.u1);
            v0 = v0.min(p.v0);
            v1 = v1.max(p.v1);
        }
        (u0, u1, v0, v1)
    }
}

/// How the line along `d` crosses the surface at `(u, v)`.
fn crossing(s: &Surface, u: f64, v: f64, d: Vector) -> Transition {
    let (_, du, dv) = s.d1(u, v);
    let n = du.cross(dv);
    let (nl, dl2) = (n.length(), d.length_squared());
    if !(nl > 1e-12) || !(dl2 > 1e-12) {
        return Transition::Tangent;
    }
    let c = n.dot(d) / (nl * dl2.sqrt());
    if c < -1e-12 {
        Transition::In
    } else if c > 1e-12 {
        Transition::Out
    } else {
        Transition::Tangent
    }
}

/// The line parameters where the line meets the elementary surface `basis` (all its
/// sheets, the whole surface).
fn closed_form_roots(basis: &Surface, origin: Point, d: Vector) -> Vec<f64> {
    let local = |frame: &crate::math::Frame, v: Vector| (v.dot(Vector::from(frame.x())), v.dot(Vector::from(frame.y())), v.dot(Vector::from(frame.z())));
    let quadratic = |a: f64, b: f64, c: f64| -> Vec<f64> {
        let scale = a.abs().max(b.abs()).max(c.abs());
        // (Identically zero: a line lying in the surface.)
        if !(scale > 0.0) || (a.abs() <= 1e-24 * scale.max(1.0) && b.abs() <= 1e-12 * scale.max(1.0)) {
            return Vec::new();
        }
        let r = PolynomialRoots::quadratic(a, b, c);
        if r.done && !r.infinite { r.roots } else { Vec::new() }
    };
    match basis {
        Surface::Plane(plane) => {
            let n = Vector::from(plane.normal());
            let denom = d.dot(n);
            if denom.abs() < 1e-12 * d.length() {
                return Vec::new();
            }
            vec![-(origin - plane.origin()).dot(n) / denom]
        }
        Surface::Cylinder(c) => {
            let f = c.frame();
            let (qx, qy, _) = local(f, origin - f.origin());
            let (dx, dy, _) = local(f, d);
            let a = dx * dx + dy * dy;
            if a <= 1e-24 * d.length_squared() {
                return Vec::new();
            }
            quadratic(a, 2.0 * (qx * dx + qy * dy), qx * qx + qy * qy - c.radius() * c.radius())
        }
        Surface::Sphere(s) => {
            let q = origin - s.frame().origin();
            quadratic(d.dot(d), 2.0 * q.dot(d), q.dot(q) - s.radius() * s.radius())
        }
        Surface::Cone(c) => {
            let f = c.frame();
            let (qx, qy, qz) = local(f, origin - f.origin());
            let (dx, dy, dz) = local(f, d);
            let k = c.semi_angle().tan();
            let r0 = c.radius() + k * qz;
            quadratic(dx * dx + dy * dy - k * k * dz * dz, 2.0 * (qx * dx + qy * dy - k * dz * r0), qx * qx + qy * qy - r0 * r0)
        }
        Surface::Torus(t) => torus_roots(t, origin, d),
        _ => Vec::new(),
    }
}

/// The real roots of the torus's quartic along the line, re-based at the line's point
/// nearest the centre, polished and checked on the torus.
fn torus_roots(t: &super::surface::Torus, origin: Point, d: Vector) -> Vec<f64> {
    let f = t.frame();
    let local = |v: Vector| (v.dot(Vector::from(f.x())), v.dot(Vector::from(f.y())), v.dot(Vector::from(f.z())));
    let (big, small) = (t.major_radius(), t.minor_radius());
    let dd = d.dot(d);
    let q = origin - f.origin();
    let w0 = -q.dot(d) / dd;
    let (qx, qy, qz) = local(q + d * w0);
    let (dx, dy, dz) = local(d);
    let _ = (dx, dy);
    let a = qx * qx + qy * qy + qz * qz;
    let k = a + big * big - small * small;
    let r2 = 4.0 * big * big;
    let c4 = dd * dd;
    let c2 = 2.0 * dd * k - r2 * (dd - dz * dz);
    let c1 = 2.0 * r2 * qz * dz;
    let c0 = k * k - r2 * (a - qz * qz);
    let quartic = |s: f64| (((c4 * s) * s + c2) * s + c1) * s + c0;
    let slope = |s: f64| ((4.0 * c4 * s) * s + 2.0 * c2) * s + c1;
    let roots = PolynomialRoots::quartic(c4, 0.0, c2, c1, c0);
    if !roots.done || roots.infinite {
        return Vec::new();
    }
    let at = |s: f64| (qx + dx * s, qy + dy * s, qz + dz * s);
    let off_torus = |s: f64| {
        let (x, y, z) = at(s);
        ((x.hypot(y) - big).hypot(z) - small).abs()
    };
    let guard = 1e-6 * (big + small);
    let mut out: Vec<f64> = Vec::new();
    for mut s in roots.roots {
        // (Newton on the quartic while it gains.)
        for _ in 0..4 {
            let (fv, df) = (quartic(s), slope(s));
            if !(df.abs() > 0.0) {
                break;
            }
            let next = s - fv / df;
            if !next.is_finite() || quartic(next).abs() >= fv.abs() {
                break;
            }
            s = next;
        }
        // (And on the distance to the tube, where the quartic is flat at a double root.)
        if off_torus(s) <= guard && !out.iter().any(|&o| (o - s).abs() <= 1e-12 * (1.0 + s.abs())) {
            out.push(s);
        }
    }
    out.into_iter().map(|s| w0 + s).collect()
}

/// Where the line `origin + t·d` lies within the box: `(t_enter, t_leave)`, both ways along
/// the line; `None` if it misses the box.
fn line_in_box(b: &BoundingBox, origin: Point, d: Vector) -> Option<(f64, f64)> {
    if b.is_empty() {
        return None;
    }
    let (mut enter, mut leave) = (f64::NEG_INFINITY, f64::INFINITY);
    for (o, d, lo, hi) in [(origin.x, d.x, b.min.x, b.max.x), (origin.y, d.y, b.min.y, b.max.y), (origin.z, d.z, b.min.z, b.max.z)] {
        if d == 0.0 {
            if o < lo || o > hi {
                return None;
            }
            continue;
        }
        let (a, c) = ((lo - o) / d, (hi - o) / d);
        enter = enter.max(a.min(c));
        leave = leave.min(a.max(c));
    }
    (enter <= leave).then_some((enter, leave))
}

/// The surface cut into patches over the box, each nearly flat with nearly parallel
/// normals (within the limits on their number).
fn sample(surface: &Surface, (u0, u1, v0, v1): (f64, f64, f64, f64)) -> Sampled {
    // Starting cuts: a B-spline's knots within the box.
    let cuts = |lo: f64, hi: f64, knots: &[f64]| -> Vec<f64> {
        let mut out = vec![lo];
        out.extend(knots.iter().copied().filter(|&k| k > lo + 1e-12 * (hi - lo) && k < hi - 1e-12 * (hi - lo)));
        out.push(hi);
        out
    };
    let (us, vs) = match untrimmed(surface) {
        Surface::BSpline(b) => {
            // (Periodic knots shifted into the box.)
            let shifted = |knots: &[f64], period: Option<f64>, lo: f64, hi: f64| -> Vec<f64> {
                match period {
                    Some(p) if p > 0.0 => {
                        let mut out = Vec::new();
                        for &k in knots {
                            let mut x = in_period(k, lo, lo + p);
                            while x < hi {
                                out.push(x);
                                x += p;
                            }
                        }
                        out.sort_by(f64::total_cmp);
                        out
                    }
                    _ => knots.to_vec(),
                }
            };
            let up = b.u_period().filter(|_| b.is_u_periodic());
            let vp = b.v_period().filter(|_| b.is_v_periodic());
            (cuts(u0, u1, &shifted(b.u_knots(), up, u0, u1)), cuts(v0, v1, &shifted(b.v_knots(), vp, v0, v1)))
        }
        _ => {
            let even = |lo: f64, hi: f64| (0..=4).map(|i| lo + (hi - lo) * i as f64 / 4.0).collect::<Vec<f64>>();
            (even(u0, u1), even(v0, v1))
        }
    };
    // (Thinned where there are too many knots for the budget.)
    let thin = |cuts: Vec<f64>, most: usize| -> Vec<f64> {
        if cuts.len() <= most + 1 {
            return cuts;
        }
        let step = (cuts.len() - 1).div_ceil(most);
        let mut out: Vec<f64> = cuts.iter().copied().step_by(step).collect();
        if *out.last().expect("cuts") != *cuts.last().expect("cuts") {
            out.push(*cuts.last().expect("cuts"));
        }
        out
    };
    let (us, vs) = (thin(us, 48), thin(vs, 48));
    let mut queue: Vec<(f64, f64, f64, f64, usize)> = Vec::new();
    for wu in us.windows(2) {
        for wv in vs.windows(2) {
            queue.push((wu[0], wu[1], wv[0], wv[1], 0));
        }
    }
    let mut patches: Vec<Patch> = Vec::new();
    let mut size_box = BoundingBox::from_points(std::iter::empty());
    // Breadth first, so that the budget is spread evenly.
    let mut k = 0;
    while k < queue.len() {
        let (a, b, c, d, depth) = queue[k];
        k += 1;
        let (patch, flat) = make_patch(surface, a, b, c, d);
        let room = patches.len() + (queue.len() - k) + 4 <= MAX_PATCHES;
        if !flat && depth < MAX_DEPTH && room {
            let (um, vm) = (0.5 * (a + b), 0.5 * (c + d));
            queue.push((a, um, c, vm, depth + 1));
            queue.push((um, b, c, vm, depth + 1));
            queue.push((a, um, vm, d, depth + 1));
            queue.push((um, b, vm, d, depth + 1));
            continue;
        }
        size_box.add_box(&patch.bbox);
        patches.push(patch);
    }
    let boxes: Vec<BoundingBox> = patches.iter().map(|p| p.bbox).collect();
    let size = if size_box.is_empty() { 1.0 } else { size_box.size().length() };
    Sampled { tree: Bvh::new(&boxes), patches, size }
}

/// The patch over `[a, b] × [c, d]`, and whether it is flat enough to stop cutting.
fn make_patch(surface: &Surface, a: f64, b: f64, c: f64, d: f64) -> (Patch, bool) {
    const N: usize = 4;
    let mut grid = [[Point::ORIGIN; N + 1]; N + 1];
    for (i, row) in grid.iter_mut().enumerate() {
        for (j, q) in row.iter_mut().enumerate() {
            *q = surface.point(a + (b - a) * i as f64 / N as f64, c + (d - c) * j as f64 / N as f64);
        }
    }
    let corners = [grid[0][0], grid[N][0], grid[N][N], grid[0][N]];
    // How far the surface strays from the bilinear patch of its corners, and between
    // neighbouring samples (from second differences along both directions).
    let mut deviation: f64 = 0.0;
    let mut between: f64 = 0.0;
    for i in 0..=N {
        for j in 0..=N {
            let (s, t) = (i as f64 / N as f64, j as f64 / N as f64);
            let bilinear = Point::ORIGIN
                + ((corners[0] - Point::ORIGIN) * ((1.0 - s) * (1.0 - t))
                    + (corners[1] - Point::ORIGIN) * (s * (1.0 - t))
                    + (corners[2] - Point::ORIGIN) * (s * t)
                    + (corners[3] - Point::ORIGIN) * ((1.0 - s) * t));
            deviation = deviation.max(grid[i][j].distance(bilinear));
            if i > 0 && i < N {
                let second = (grid[i - 1][j] - grid[i][j]) + (grid[i + 1][j] - grid[i][j]);
                between = between.max(second.length());
            }
            if j > 0 && j < N {
                let second = (grid[i][j - 1] - grid[i][j]) + (grid[i][j + 1] - grid[i][j]);
                between = between.max(second.length());
            }
        }
    }
    let mut bbox = BoundingBox::from_points(grid.iter().flat_map(|r| r.iter().copied()));
    let diagonal = bbox.size().length();
    // (The surface between samples bulges by about a quarter of the second difference;
    // twice that, and a little more, is room enough.)
    let margin = 0.5 * between + 0.05 * deviation / (N * N) as f64 + 1e-7 + 1e-12 * (Point::ORIGIN - bbox.center()).length();
    bbox = bbox.enlarged(margin);
    let normal_at = |u: f64, v: f64| {
        let (_, du, dv) = surface.d1(u, v);
        let n = du.cross(dv);
        let l = n.length();
        if l > 1e-300 && l.is_finite() { n * (1.0 / l) } else { Vector::new(0.0, 0.0, 0.0) }
    };
    let (um, vm) = (0.5 * (a + b), 0.5 * (c + d));
    let normal = normal_at(um, vm);
    let mut spread: f64 = 0.0;
    for (u, v) in [(a, c), (b, c), (b, d), (a, d), (um, c), (um, d), (a, vm), (b, vm)] {
        let n = normal_at(u, v);
        if n.length_squared() > 0.5 && normal.length_squared() > 0.5 {
            spread = spread.max(n.dot(normal).clamp(-1.0, 1.0).acos());
        }
    }
    if normal.length_squared() < 0.5 {
        spread = std::f64::consts::PI;
    }
    let flat = spread <= 0.35 && deviation <= 0.05 * diagonal.max(1e-300);
    (Patch { u0: a, u1: b, v0: c, v1: d, bbox, corners, normal, spread }, flat)
}

/// Newton's method on `S(u, v) = origin + w·d` from starts in the patches.
struct Solver<'a> {
    surface: &'a Surface,
    origin: Point,
    d: Vector,
    /// The parameter box the steps keep to (the patches' union).
    bounds: (f64, f64, f64, f64),
    size: f64,
}

impl Solver<'_> {
    /// The solutions started from within `patch` (cut further where the line nearly touches
    /// it), each passed to `found`.
    fn patch(&mut self, patch: &Patch, depth: usize, found: &mut dyn FnMut(f64, f64, f64)) {
        let dl = self.d.length();
        let grazing = patch.normal.length_squared() > 0.5 && (patch.normal.dot(self.d) / dl).abs() < (patch.spread + 0.2).sin().min(1.0);
        let wide = patch.spread > 0.35;
        if (grazing || wide) && depth < MAX_TANGENT_DEPTH {
            let (um, vm) = (0.5 * (patch.u0 + patch.u1), 0.5 * (patch.v0 + patch.v1));
            for (a, b, c, d) in [(patch.u0, um, patch.v0, vm), (um, patch.u1, patch.v0, vm), (patch.u0, um, vm, patch.v1), (um, patch.u1, vm, patch.v1)] {
                let (sub, _) = make_patch(self.surface, a, b, c, d);
                if line_in_box(&sub.bbox, self.origin, self.d).is_some() {
                    self.patch(&sub, depth + 1, found);
                }
            }
            return;
        }
        for (u, v, w) in self.starts(patch) {
            if let Some((u, v, w)) = self.newton(u, v, w) {
                found(u, v, w);
            }
        }
    }

    /// Starts in the patch: where the line meets its two corner triangles, and its centre
    /// with the line's point nearest it.
    fn starts(&self, patch: &Patch) -> Vec<(f64, f64, f64)> {
        let mut out = Vec::with_capacity(3);
        let c = patch.corners;
        let (a, b, cc, d) = (patch.u0, patch.u1, patch.v0, patch.v1);
        // Triangle (00, 10, 11) with barycentrics mapped to (u, v); then (00, 11, 01).
        for (p0, p1, p2, uv0, uv1, uv2) in [(c[0], c[1], c[2], (a, cc), (b, cc), (b, d)), (c[0], c[2], c[3], (a, cc), (b, d), (a, d))] {
            if let Some((w, s, t)) = line_triangle(self.origin, self.d, p0, p1, p2) {
                let u = uv0.0 + s * (uv1.0 - uv0.0) + t * (uv2.0 - uv0.0);
                let v = uv0.1 + s * (uv1.1 - uv0.1) + t * (uv2.1 - uv0.1);
                out.push((u, v, w));
            }
        }
        let (um, vm) = (0.5 * (a + b), 0.5 * (cc + d));
        let q = self.surface.point(um, vm);
        let w = (q - self.origin).dot(self.d) / self.d.length_squared();
        out.push((um, vm, w));
        out
    }

    /// Newton steps on `S(u, v) − (origin + w·d) = 0` from `(u, v, w)`; the solution where
    /// they converge within the box (a little beyond).
    fn newton(&self, mut u: f64, mut v: f64, mut w: f64) -> Option<(f64, f64, f64)> {
        let s = self.surface;
        let (u0, u1, v0, v1) = self.bounds;
        let (mu, mv) = (0.05 * (u1 - u0), 0.05 * (v1 - v0));
        let (up, vp) = (s.is_u_periodic(), s.is_v_periodic());
        let tolerance = 1e-12 * (1.0 + self.size + (self.origin - Point::ORIGIN).length());
        let accept = 1e-9 * (1.0 + self.size);
        let mut residual = f64::INFINITY;
        for _ in 0..40 {
            let (p, su, sv) = s.d1(u, v);
            let f = p - (self.origin + self.d * w);
            residual = f.length();
            if !residual.is_finite() {
                return None;
            }
            let c = -self.d;
            let bc = sv.cross(c);
            let det = su.dot(bc);
            let scale = su.length() * sv.length() * c.length();
            if !(det.abs() > 1e-14 * scale) {
                break;
            }
            let r = -f;
            let du = r.dot(bc) / det;
            let dv = su.dot(r.cross(c)) / det;
            let dw = su.dot(sv.cross(r)) / det;
            u += du;
            v += dv;
            w += dw;
            if !up {
                u = u.clamp(u0 - mu, u1 + mu);
            }
            if !vp {
                v = v.clamp(v0 - mv, v1 + mv);
            }
            let small = du.abs() <= 1e-15 * (1.0 + u.abs()) && dv.abs() <= 1e-15 * (1.0 + v.abs()) && dw.abs() <= 1e-15 * (1.0 + w.abs());
            if small || residual <= tolerance * 1e-3 {
                let (p, _, _) = s.d1(u, v);
                residual = p.distance(self.origin + self.d * w);
                break;
            }
        }
        (residual <= accept).then_some((u, v, w))
    }
}

/// Where the line meets the triangle `(p0, p1, p2)`: the line parameter and the barycentric
/// coordinates of `p1` and `p2` (Möller–Trumbore), a little slack round the edges.
fn line_triangle(o: Point, d: Vector, p0: Point, p1: Point, p2: Point) -> Option<(f64, f64, f64)> {
    let (e1, e2) = (p1 - p0, p2 - p0);
    let h = d.cross(e2);
    let det = e1.dot(h);
    if !(det.abs() > 1e-300) {
        return None;
    }
    let inv = 1.0 / det;
    let s = o - p0;
    let a = s.dot(h) * inv;
    let q = s.cross(e1);
    let b = d.dot(q) * inv;
    let slack = 0.05;
    if a < -slack || b < -slack || a + b > 1.0 + slack {
        return None;
    }
    let w = e2.dot(q) * inv;
    w.is_finite().then_some((w, a, b))
}

#[cfg(test)]
mod tests;
