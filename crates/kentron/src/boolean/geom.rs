//! Geometric helpers of the boolean operations: distances to elementary surfaces, rays,
//! projections onto curves, polylines in parameter space.

use crate::geometry::curve::Curve;
use crate::geometry::curve2::Curve2;
use crate::geometry::surface::Surface;
use crate::geometry::{ParametricCurve, ParametricSurface};
use crate::math::{BoundingBox, Bvh, Point, Vector};
use crate::tolerance;

/// Distance below which points coincide (vertices merge).
pub(super) const TOL: f64 = tolerance::LINEAR;

/// Distance below which a point counts as lying on a curve or surface when classifying.
pub(super) const ON_TOL: f64 = 1.0e-6;

/// [`ON_TOL`] with the looser tolerances on (see [`crate::tolerance::loosened`]): for shapes
/// whose points meant to coincide (a cylinder through a box's corners) are only that close.
pub(super) const LOOSE_ON_TOL: f64 = 1.0e-5;

/// The distance below which a point counts as lying on a curve or surface when
/// classifying: [`ON_TOL`], or [`LOOSE_ON_TOL`] with the looser tolerances on.
pub(super) fn on_tol() -> f64 {
    if crate::tolerance::loosened() { LOOSE_ON_TOL } else { ON_TOL }
}

/// A point in a surface's (u, v) parameter space.
pub(super) type Uv = [f64; 2];

/// Signed distance from `p` to an elementary surface: along the normal for a plane, from
/// the axis (centre) outwards for a cylinder (sphere).
/// Whether [`signed_distance_with_gradient`] handles the surface.
pub(super) fn closed_form_gradient(surface: &Surface) -> bool {
    matches!(surface, Surface::Plane(_) | Surface::Cylinder(_) | Surface::Sphere(_) | Surface::Cone(_) | Surface::Revolution(_) | Surface::Extrusion(_))
}

/// [`signed_distance`] and its gradient in `p`, for the analytic surfaces it handles in
/// closed form (`None` for others).
pub(super) fn signed_distance_with_gradient(surface: &Surface, p: Point) -> Option<(f64, Vector)> {
    let unit = |v: Vector| {
        let len = v.length();
        if len > 0.0 { v * (1.0 / len) } else { v }
    };
    Some(match surface {
        Surface::Plane(plane) => {
            let n = plane_normal(plane);
            ((p - plane.origin()).dot(n), n)
        }
        Surface::Cylinder(c) => {
            let f = c.frame();
            let d = p - f.origin();
            let z = Vector::from(f.z());
            let radial = d - z * d.dot(z);
            (radial.length() - c.radius(), unit(radial))
        }
        Surface::Sphere(s) => (p.distance(s.frame().origin()) - s.radius(), unit(p - s.frame().origin())),
        Surface::Cone(c) => {
            let f = c.frame();
            let d = p - f.origin();
            let z = Vector::from(f.z());
            let h = d.dot(z);
            let radial = d - z * h;
            let (sin, cos) = c.semi_angle_sin_cos();
            // To the nearer nappe: past the apex the surface faces the other way.
            let (near, beyond) = ((radial.length() - c.radius()) * cos - h * sin, (radial.length() + c.radius()) * cos + h * sin);
            if beyond.abs() < near.abs() {
                (beyond, unit(radial) * cos + z * sin)
            } else {
                (near, unit(radial) * cos - z * sin)
            }
        }
        Surface::Revolution(_) | Surface::Extrusion(_) => {
            // Along the normal at the foot point, which is the gradient there.
            let uv = uv_of(surface, p);
            let n = normal_at(surface, uv);
            ((p - surface_point(surface, uv)).dot(n), n)
        }
        _ => return None,
    })
}

pub(super) fn signed_distance(surface: &Surface, p: Point) -> f64 {
    match surface {
        Surface::Plane(plane) => (p - plane.origin()).dot(plane_normal(plane)),
        Surface::Cylinder(c) => {
            let f = c.frame();
            let d = p - f.origin();
            let z = Vector::from(f.z());
            (d - z * d.dot(z)).length() - c.radius()
        }
        Surface::Sphere(s) => p.distance(s.frame().origin()) - s.radius(),
        Surface::Cone(c) => {
            // In the half-plane through the axis and p: the distance to the generator, or to
            // the one opposite past the apex if nearer (the surface facing the other way
            // there).
            let f = c.frame();
            let d = p - f.origin();
            let z = Vector::from(f.z());
            let h = d.dot(z);
            let radial = (d - z * h).length();
            let (sin, cos) = c.semi_angle_sin_cos();
            let (near, beyond) = ((radial - c.radius()) * cos - h * sin, (radial + c.radius()) * cos + h * sin);
            if beyond.abs() < near.abs() { beyond } else { near }
        }
        Surface::Revolution(_) | Surface::Extrusion(_) => {
            // Along the normal at the foot point; or, on an extrusion where there is no normal
            // (the foot at an end of its curve where that has no tangent, as font outlines'
            // curves start), the whole way to it: along nothing, the point would seem on the
            // surface. (Not at a surface of revolution's pole, where the distance changes sign
            // through the axis.)
            let uv = uv_of(surface, p);
            let (r, n) = (p - surface_point(surface, uv), normal_at(surface, uv));
            if n.length() > 0.5 || matches!(surface, Surface::Revolution(_)) { r.dot(n) } else { r.length() }
        }
        other if freeform(other) => {
            let foot = foot_on(surface, p);
            (p - foot.point).dot(unit_normal(foot.du, foot.dv))
        }
        other => unsupported(other),
    }
}

pub(super) fn unsupported(surface: &Surface) -> ! {
    let kind = format!("{surface:?}");
    todo!("boolean operations on {} faces", kind.split('(').next().unwrap_or_default().to_lowercase())
}

/// The (u, v) parameters of a point on the surface (`u` in `[0, 2π)` when angular).
pub(super) fn uv_of(surface: &Surface, p: Point) -> Uv {
    if let Some((u, v)) = surface.parameters(p) {
        return [u, v];
    }
    if !freeform(surface) {
        unsupported(surface);
    }
    let foot = foot_on(surface, p);
    [foot.u, foot.v]
}

/// The parameters of `p` on the surface; `None` on a free-form surface where its foot lies
/// farther than `reach` (beyond the surface's edge, the foot found on a bound: a point in
/// the surface's tangent extension has no place on it, though its signed distance is 0).
pub(super) fn uv_within(surface: &Surface, p: Point, reach: f64) -> Option<Uv> {
    if !freeform(surface) {
        return Some(uv_of(surface, p));
    }
    let foot = foot_on(surface, p);
    (foot.point.distance(p) <= reach).then_some([foot.u, foot.v])
}

/// The foot of `p` on a free-form surface: Newton steps from the nearest of a grid of the
/// surface's points. Where that leads astray (to a foot on the surface's bounds, or farther
/// than the sample it started from) and another part of the surface passes nearly as near
/// (the next turn of one winding round), from the nearest sample there too: the nearer foot.
pub(super) fn foot_on(surface: &Surface, p: Point) -> crate::geometry::extrema::SurfaceFoot {
    let (best, start) = with_samples(surface, |grid| {
        let best = grid.nearest(p, 0, grid.samples.len());
        (best, (grid.samples[best.1].0, grid.samples[best.1].1))
    });
    let foot = crate::geometry::extrema::newton_foot_from_sample(p, surface, start);
    let (u0, u1, v0, v1) = surface.bounds();
    let on_bound = |x: f64, a: f64, b: f64, periodic: bool| !periodic && ((x - a).abs() <= 1e-12 * (b - a).abs().max(1.0) || (x - b).abs() <= 1e-12 * (b - a).abs().max(1.0));
    let distance = foot.point.distance(p);
    let astray = distance > best.0.sqrt() * (1.0 + 1e-9) + 1e-12 || on_bound(foot.u, u0, u1, surface.u_period().is_some()) || on_bound(foot.v, v0, v1, surface.v_period().is_some());
    if !astray || winding_samples(surface) == 0 {
        return foot;
    }
    // (Away from the nearest sample by more than a cell either way.)
    let (_, nv) = grid_steps(surface);
    let other = with_samples(surface, |grid| {
        // (The stretches of the grid before, between and after the rows of the 3 × 3 cells
        // around the nearest sample, in order.)
        let row = nv + 1;
        let (bi, bj) = (best.1 / row, best.1 % row);
        let (first, last) = (bj.saturating_sub(1), (bj + 1).min(nv));
        let mut other = (f64::INFINITY, 0);
        let mut from = 0;
        for i in bi.saturating_sub(1)..=(bi + 1).min(grid.samples.len() / row - 1) {
            let found = grid.nearest(p, from, i * row + first);
            if found.0 < other.0 {
                other = found;
            }
            from = i * row + last + 1;
        }
        let found = grid.nearest(p, from, grid.samples.len());
        if found.0 < other.0 {
            other = found;
        }
        (other.0 < 4.0 * best.0).then(|| (grid.samples[other.1].0, grid.samples[other.1].1))
    });
    let Some(other) = other else { return foot };
    let second = crate::geometry::extrema::newton_foot_from_sample(p, surface, other);
    if second.point.distance(p) < distance { second } else { foot }
}

/// Where the last point's foot was found on a free-form surface: its parameters, the point
/// and its distance, and the unit normal there.
#[derive(Clone, Copy)]
pub(super) struct Foot {
    uv: Uv,
    point: Point,
    distance: f64,
    normal: Vector,
}

impl Foot {
    /// The unit normal at the foot: the gradient of the signed distance there.
    pub(super) fn normal(&self) -> Vector {
        self.normal
    }
}

/// The signed distance from `p` to a surface; on a free-form surface numerically, from the
/// foot of a nearby point (as along a grid or a curve) when given. The distance moves by no
/// more than the point does, so a start that leads farther than that went astray, and the
/// search starts afresh.
pub(super) fn signed_distance_from(surface: &Surface, p: Point, last: Option<Foot>) -> (f64, Option<Foot>) {
    if !freeform(surface) {
        return (signed_distance(surface, p), None);
    }
    let fresh = || foot_on(surface, p);
    let foot = match last {
        Some(f) => {
            let foot = crate::geometry::extrema::newton_foot_to_the_last(p, surface, (f.uv[0], f.uv[1]));
            if foot.point.distance(p) <= f.distance + p.distance(f.point) + 1e-9 { foot } else { fresh() }
        }
        None => fresh(),
    };
    let normal = unit_normal(foot.du, foot.dv);
    let signed = (p - foot.point).dot(normal);
    (signed, Some(Foot { uv: [foot.u, foot.v], point: p, distance: foot.point.distance(p), normal }))
}

/// The signed distance from `p` to a surface and, on a free-form surface, the unit normal
/// at the same foot (the distance's gradient there): both from one search, as
/// [`signed_distance`] and [`normal_at_point`] would give them. (`None` for the normal where
/// the point's parameters have a closed form, which `normal_at_point` takes.)
pub(super) fn signed_distance_and_normal(surface: &Surface, p: Point) -> (f64, Option<Vector>) {
    if !freeform(surface) || surface.parameters(p).is_some() {
        return (signed_distance(surface, p), None);
    }
    let foot = foot_on(surface, p);
    let normal = unit_normal(foot.du, foot.dv);
    ((p - foot.point).dot(normal), Some(normal))
}

/// Surfaces without closed-form parameters, handled numerically (B-spline and Bezier
/// patches, offsets and trimmed surfaces of them).
pub(super) fn freeform(surface: &Surface) -> bool {
    match surface {
        Surface::BSpline(_) | Surface::Bezier(_) => true,
        Surface::Offset(o) => freeform(o.basis()) || matches!(o.basis(), Surface::Plane(_) | Surface::Cylinder(_) | Surface::Cone(_) | Surface::Sphere(_)),
        Surface::Trimmed(t) => freeform(t.basis()),
        _ => false,
    }
}

/// Grids of points of the free-form surfaces met lately, to start point inversions from
/// (shared by the threads of an operation).
static SAMPLES: std::sync::Mutex<Vec<(crate::geometry::extrema::SurfaceKey, std::sync::Arc<SampleGrid>)>> = std::sync::Mutex::new(Vec::new());

/// A surface's sample grid: its points `(u, v, p)`, row by row along u, and their
/// coordinates again in arrays of their own, for scans of all of them.
struct SampleGrid {
    samples: Vec<(f64, f64, Point)>,
    x: Vec<f64>,
    y: Vec<f64>,
    z: Vec<f64>,
}

impl SampleGrid {
    fn new(samples: Vec<(f64, f64, Point)>) -> SampleGrid {
        let x = samples.iter().map(|s| s.2.x).collect();
        let y = samples.iter().map(|s| s.2.y).collect();
        let z = samples.iter().map(|s| s.2.z).collect();
        SampleGrid { samples, x, y, z }
    }

    /// The sample among `from..to` nearest to `p` (the first of equally near ones) and its
    /// squared distance, `(INFINITY, 0)` if none is nearer than that. (Four at a time, each
    /// lane keeping the first nearest of its own: four chains of compares the processor
    /// runs side by side, where one running minimum makes every sample wait for the one
    /// before. The squared distances are the same sums as `(q - p).dot(q - p)`.)
    fn nearest(&self, p: Point, from: usize, to: usize) -> (f64, usize) {
        let (x, y, z) = (&self.x[from..to], &self.y[from..to], &self.z[from..to]);
        let mut best = [f64::INFINITY; 4];
        let mut at = [0usize; 4];
        let mut offer = |k: usize, lane: usize| {
            let (dx, dy, dz) = (x[k] - p.x, y[k] - p.y, z[k] - p.z);
            let d = dx * dx + dy * dy + dz * dz;
            if d < best[lane] {
                best[lane] = d;
                at[lane] = k;
            }
        };
        let whole = (to - from) / 4 * 4;
        for k in (0..whole).step_by(4) {
            for lane in 0..4 {
                offer(k + lane, lane);
            }
        }
        for k in whole..to - from {
            offer(k, k % 4);
        }
        let mut out = (f64::INFINITY, 0);
        for lane in 0..4 {
            if best[lane] < out.0 || (best[lane] == out.0 && best[lane] < f64::INFINITY && from + at[lane] < out.1) {
                out = (best[lane], from + at[lane]);
            }
        }
        out
    }
}

/// The steps of a surface's sample grid along u and v.
fn grid_steps(surface: &Surface) -> (usize, usize) {
    // A few samples per knot span and degree (a bilinear patch needs few); more on one
    // winding round twice or more, its turns close together
    // ([`crate::geometry::extrema::windings_from`]).
    let by = |spans: usize, degree: usize| (4 * spans * degree.max(1)).clamp(4, 32);
    match surface {
        Surface::BSpline(b) => {
            let (wu, wv) = crate::geometry::extrema::windings_from(surface, 2.0);
            (by(b.u_knots().len() - 1, b.u_degree()).max(wu), by(b.v_knots().len() - 1, b.v_degree()).max(wv))
        }
        Surface::Bezier(b) => (by(1, b.u_degree()), by(1, b.v_degree())),
        _ => (32, 32),
    }
}


/// The samples of a surface winding round twice or more, its turns close together
/// ([`crate::geometry::extrema::windings_from`]), along its more winding direction; 0 for
/// others.
pub(super) fn winding_samples(surface: &Surface) -> usize {
    let (u, v) = crate::geometry::extrema::windings_from(surface, 2.0);
    u.max(v)
}

/// Runs `f` on the surface's sample grid (points `(u, v, p)`, row by row along u), cached.
fn with_samples<R>(surface: &Surface, f: impl FnOnce(&SampleGrid) -> R) -> R {
    let lock = || SAMPLES.lock().unwrap_or_else(|poisoned| poisoned.into_inner());
    let key = crate::geometry::extrema::SurfaceKey::of(surface);
    let found = lock().iter().find(|(k, _)| *k == key).map(|(_, grid)| grid.clone());
    let grid = match found {
        Some(grid) => grid,
        None => {
            // (Made outside the lock; another thread may make it too.)
            let grid = std::sync::Arc::new(SampleGrid::new(sample_grid(surface)));
            let mut cache = lock();
            if !cache.iter().any(|(k, _)| *k == key) {
                if cache.len() >= 64 {
                    cache.remove(0);
                }
                cache.push((key, grid.clone()));
            }
            grid
        }
    };
    f(&grid)
}

fn sample_grid(surface: &Surface) -> Vec<(f64, f64, Point)> {
    let (u0, u1, v0, v1) = surface.bounds();
    let (nu, nv) = grid_steps(surface);
    let mut out = Vec::with_capacity((nu + 1) * (nv + 1));
    for i in 0..=nu {
        for j in 0..=nv {
            let (u, v) = (u0 + (u1 - u0) * i as f64 / nu as f64, v0 + (v1 - v0) * j as f64 / nv as f64);
            out.push((u, v, surface.point(u, v)));
        }
    }
    out
}

/// The point of the surface at (u, v).
pub(super) fn surface_point(surface: &Surface, uv: Uv) -> Point {
    surface.point(uv[0], uv[1])
}

/// True at a pole of a sphere or the apex of a cone (where `u` is undefined).
pub(super) fn is_pole(surface: &Surface, p: Point) -> bool {
    match surface {
        Surface::Cone(c) => p.distance(c.apex()) <= on_tol(),
        Surface::Sphere(s) => {
            let z = Vector::from(s.frame().z()) * s.radius();
            let c = s.frame().origin();
            p.distance(c + z) <= on_tol() || p.distance(c - z) <= on_tol()
        }
        _ => false,
    }
}

/// The unit normal `∂u × ∂v` at (u, v).
pub(super) fn normal_at(surface: &Surface, uv: Uv) -> Vector {
    let (_, su, sv) = surface.d1(uv[0], uv[1]);
    unit_normal(su, sv)
}

/// `su × sv` made unit (zero where they are parallel).
fn unit_normal(su: Vector, sv: Vector) -> Vector {
    let n = su.cross(sv);
    let len = n.length();
    if len > 0.0 { n * (1.0 / len) } else { n }
}

/// A plane's `∂u × ∂v`: its frame's Z, or the opposite for a left-handed frame (a
/// reflected plane).
pub(super) fn plane_normal(plane: &crate::geometry::surface::Plane) -> Vector {
    let f = plane.frame();
    Vector::from(f.x()).cross(Vector::from(f.y()))
}

/// The unit normal `∂u × ∂v` at a point of the surface.
pub(super) fn normal_at_point(surface: &Surface, p: Point) -> Vector {
    match surface {
        Surface::Plane(plane) => plane_normal(plane),
        _ => normal_at(surface, uv_of(surface, p)),
    }
}

/// The parameters `t > 0` where the ray `origin + t·direction` meets the surface, or
/// `None` if it grazes it (tangent, ambiguous for counting crossings).
pub(super) fn ray_hits(surface: &Surface, origin: Point, direction: Vector, bbox: &(Point, Point)) -> Option<Vec<f64>> {
    let quadratic = |a: f64, b: f64, c: f64| -> Option<Vec<f64>> {
        let disc = b * b - 4.0 * a * c;
        if disc < 0.0 {
            return Some(Vec::new());
        }
        if disc <= 1e-12 * (b * b).max(1e-300) {
            return None;
        }
        let root = disc.sqrt();
        Some(vec![(-b - root) / (2.0 * a), (-b + root) / (2.0 * a)])
    };
    let hits = match surface {
        Surface::Plane(plane) => {
            let n = Vector::from(plane.normal());
            let denom = n.dot(direction);
            if denom.abs() <= 1e-12 {
                return Some(Vec::new());
            }
            vec![n.dot(plane.origin() - origin) / denom]
        }
        Surface::Cylinder(c) => {
            let f = c.frame();
            let z = Vector::from(f.z());
            let o = origin - f.origin();
            let (o, d) = (o - z * o.dot(z), direction - z * direction.dot(z));
            if d.length_squared() <= 1e-24 {
                return Some(Vec::new());
            }
            quadratic(d.dot(d), 2.0 * o.dot(d), o.dot(o) - c.radius() * c.radius())?
        }
        Surface::Sphere(s) => {
            let o = origin - s.frame().origin();
            quadratic(direction.dot(direction), 2.0 * o.dot(direction), o.dot(o) - s.radius() * s.radius())?
        }
        Surface::Cone(c) => {
            // radial² = (R + h tan α)², both nappes (hits off the face are located away).
            let f = c.frame();
            let z = Vector::from(f.z());
            let o = origin - f.origin();
            let (oh, dh) = (o.dot(z), direction.dot(z));
            let (op, dp) = (o - z * oh, direction - z * dh);
            let k = c.semi_angle().tan();
            let r0 = c.radius() + k * oh;
            let a = dp.dot(dp) - k * k * dh * dh;
            let b = 2.0 * (op.dot(dp) - k * dh * r0);
            let cc = op.dot(op) - r0 * r0;
            if a.abs() <= 1e-12 {
                if b.abs() <= 1e-12 {
                    return Some(Vec::new());
                }
                vec![-cc / b]
            } else {
                quadratic(a, b, cc)?
            }
        }
        Surface::BSpline(_) | Surface::Bezier(_) => patch_ray_hits(surface, origin, direction)?,
        Surface::Revolution(_) | Surface::Extrusion(_) => numeric_ray_hits(surface, origin, direction, bbox)?,
        other if freeform(other) => numeric_ray_hits(surface, origin, direction, bbox)?,
        other => unsupported(other),
    };
    Some(hits.into_iter().filter(|&t| t > on_tol()).collect())
}

/// Ray hits on a B-spline or Bezier patch: Newton steps on S(u, v) = origin + t·direction
/// from the cells of the sample grid the ray passes by; `None` if the ray grazes it.
fn patch_ray_hits(surface: &Surface, origin: Point, direction: Vector) -> Option<Vec<f64>> {
    let d = direction.normalized();
    // Two directions across the ray: grid points are seen on the plane across it.
    let helper = if d.x.abs() < 0.6 { Vector::new(1.0, 0.0, 0.0) } else { Vector::new(0.0, 1.0, 0.0) };
    let e1 = d.cross(helper).normalized();
    let e2 = d.cross(e1);
    let (nu, nv) = grid_steps(surface);
    let (u0, u1, v0, v1) = surface.bounds();
    let starts: Vec<(f64, f64, f64)> = with_samples(surface, |grid| {
        let at = |i: usize, j: usize| grid.samples[i * (nv + 1) + j];
        let mut starts = Vec::new();
        for i in 0..nu {
            for j in 0..nv {
                let corners = [at(i, j), at(i + 1, j), at(i, j + 1), at(i + 1, j + 1)];
                let (mut lo, mut hi) = ([f64::INFINITY; 2], [f64::NEG_INFINITY; 2]);
                let mut size = 0.0f64;
                for (k, c) in corners.iter().enumerate() {
                    let q = c.2 - origin;
                    let (a, b) = (q.dot(e1), q.dot(e2));
                    lo = [lo[0].min(a), lo[1].min(b)];
                    hi = [hi[0].max(a), hi[1].max(b)];
                    size = size.max(c.2.distance(corners[(k + 1) % 4].2)).max(c.2.distance(corners[3 - k].2));
                }
                // The patch strays from its corners' quad by a fraction of the cell.
                let margin = 0.5 * size + on_tol();
                if lo[0] - margin <= 0.0 && hi[0] + margin >= 0.0 && lo[1] - margin <= 0.0 && hi[1] + margin >= 0.0 {
                    let t = corners.iter().map(|c| (c.2 - origin).dot(d)).sum::<f64>() / 4.0;
                    starts.push((0.5 * (corners[0].0 + corners[3].0), 0.5 * (corners[0].1 + corners[3].1), t));
                }
            }
        }
        starts
    });
    let mut hits: Vec<f64> = Vec::new();
    for (mut u, mut v, mut t) in starts {
        let mut converged = false;
        let mut residual = f64::INFINITY;
        for _ in 0..30 {
            let (p, su, sv) = surface.d1(u, v);
            let f = p - (origin + d * t);
            residual = f.length();
            if residual <= 1e-11 * (1.0 + p.distance(origin)) {
                converged = true;
                break;
            }
            // [su sv -d] · (du, dv, dt) = -f, by Cramer's rule.
            let det = su.dot(sv.cross(-d));
            if det.abs() <= 1e-300 {
                break;
            }
            let du = -f.dot(sv.cross(-d)) / det;
            let dv = -su.dot(f.cross(-d)) / det;
            let dt = -su.dot(sv.cross(f)) / det;
            u = (u + du).clamp(u0, u1);
            v = (v + dv).clamp(v0, v1);
            t += dt;
        }
        if !converged {
            if residual <= 100.0 * on_tol() {
                return None; // nearly touching: a graze, or too close to tell
            }
            continue;
        }
        let normal = normal_at(surface, [u, v]);
        if normal.dot(d).abs() <= 1e-3 {
            return None;
        }
        if !hits.iter().any(|&h| (h - t).abs() <= 1e-7 * (1.0 + t.abs())) {
            hits.push(t);
        }
    }
    let scale = direction.length();
    Some(hits.into_iter().map(|t| t / scale).collect())
}

/// Ray hits on a surface without a closed form, by root finding on the signed distance
/// where the ray crosses `bbox`; `None` if the ray grazes the surface.
fn numeric_ray_hits(surface: &Surface, origin: Point, direction: Vector, bbox: &(Point, Point)) -> Option<Vec<f64>> {
    // The stretch of the ray inside the box (slabs).
    let (mut t0, mut t1) = (0.0f64, f64::INFINITY);
    for (o, d, lo, hi) in [
        (origin.x, direction.x, bbox.0.x, bbox.1.x),
        (origin.y, direction.y, bbox.0.y, bbox.1.y),
        (origin.z, direction.z, bbox.0.z, bbox.1.z),
    ] {
        if d.abs() <= 1e-300 {
            if o < lo || o > hi {
                return Some(Vec::new());
            }
            continue;
        }
        let (a, b) = ((lo - o) / d, (hi - o) / d);
        t0 = t0.max(a.min(b));
        t1 = t1.min(a.max(b));
    }
    if t0 >= t1 {
        return Some(Vec::new());
    }
    let f = |t: f64| signed_distance(surface, origin + direction * t);
    let n = 400;
    let ts: Vec<f64> = (0..=n).map(|i| t0 + (t1 - t0) * i as f64 / n as f64).collect();
    let fs: Vec<f64> = ts.iter().map(|&t| f(t)).collect();
    let mut hits = Vec::new();
    for i in 0..n {
        if fs[i] * fs[i + 1] < 0.0 {
            let (mut lo, mut hi, mut flo) = (ts[i], ts[i + 1], fs[i]);
            for _ in 0..80 {
                let mid = 0.5 * (lo + hi);
                let fm = f(mid);
                if (fm < 0.0) == (flo < 0.0) {
                    lo = mid;
                    flo = fm;
                } else {
                    hi = mid;
                }
            }
            let t = 0.5 * (lo + hi);
            // A jump of the signed distance (another sheet's foot point) is not a hit, nor is
            // passing beside a bounded patch (where the foot point stays on its border): the
            // point must lie on the surface.
            let q = origin + direction * t;
            if f(t).abs() <= on_tol() && surface_point(surface, uv_of(surface, q)).distance(q) <= 10.0 * on_tol() {
                hits.push(t);
            }
        } else if i > 0 && fs[i].abs() <= on_tol() && fs[i].abs() <= fs[i - 1].abs() && fs[i].abs() <= fs[i + 1].abs() {
            return None;
        }
    }
    Some(hits.into_iter().filter(|&t| t > on_tol()).collect())
}

/// The parameter in `[first, last]` of the point of `curve` closest to `p`, and the
/// distance to it.
pub(super) fn project(curve: &Curve, first: f64, last: f64, p: Point) -> (f64, f64) {
    project_among(curve, first, last, None, p)
}

/// The points [`project`] samples `curve` at over `first..last` (to project many points
/// onto it with [`project_among`]).
pub(super) fn projection_samples(curve: &Curve, first: f64, last: f64) -> Vec<Point> {
    (0..=PROJECTION_SAMPLES).map(|i| curve.point(first + (last - first) * i as f64 / PROJECTION_SAMPLES as f64)).collect()
}

const PROJECTION_SAMPLES: usize = 64;

/// [`project`], with the curve's samples found before ([`projection_samples`]) when given.
pub(super) fn project_among(curve: &Curve, first: f64, last: f64, samples: Option<&[Point]>, p: Point) -> (f64, f64) {
    let clamp = |t: f64| t.clamp(first.min(last), first.max(last));
    if let Curve::Line(line) = curve {
        let t = clamp((p - line.origin).dot(Vector::from(line.direction)));
        return (t, curve.point(t).distance(p));
    }
    // Coarse sampling, then Newton steps on (C(t) - p)·C'(t) = 0.
    let n = PROJECTION_SAMPLES;
    let mut best = (first, f64::INFINITY);
    let mut sample = |i: usize| {
        let t = first + (last - first) * i as f64 / n as f64;
        let d = samples.map_or_else(|| curve.point(t), |s| s[i]).distance(p);
        if d < best.1 {
            best = (t, d);
        }
    };
    let near = match curve {
        Curve::Circle(circle) => circle_samples_near(circle, first, last, n, p),
        _ => None,
    };
    match near {
        Some(near) => near.into_iter().for_each(&mut sample),
        None => (0..=n).for_each(&mut sample),
    }
    // Over a whole period of a periodic curve the parameter wraps round (a point just
    // before the curve's origin is not stuck at the range's start).
    let (lo, hi) = (first.min(last), first.max(last));
    let period = curve.period().filter(|&q| hi - lo >= q - 1e-12);
    let mut t = best.0;
    for _ in 0..30 {
        let (c, d1, d2) = curve.d2(t);
        let f = (c - p).dot(d1);
        let df = d1.dot(d1) + (c - p).dot(d2);
        if df.abs() <= 1e-300 {
            break;
        }
        let step = t - f / df;
        let next = if period.is_some() { step } else { clamp(step) };
        if (next - t).abs() <= 1e-15 * (1.0 + t.abs()) {
            t = next;
            break;
        }
        t = next;
    }
    if let Some(q) = period {
        t = lo + (t - lo).rem_euclid(q);
    }
    let d = curve.point(t).distance(p);
    if d < best.1 { (t, d) } else { best }
}

/// The samples `i` of `0..=n` (at `first + (last - first) i / n`) of a circle that can be
/// nearest to `p`, in order: those next to the angle of `p` (and the ends), when all others
/// are farther by much more than rounding (so that the samples' nearest is the same as
/// among all); `None` where that is not clear (`p` near the axis, short steps).
fn circle_samples_near(circle: &crate::geometry::curve::Circle, first: f64, last: f64, n: usize, p: Point) -> Option<Vec<usize>> {
    let f = circle.frame();
    let q = p - f.origin();
    let (x, y) = (q.dot(Vector::from(f.x())), q.dot(Vector::from(f.y())));
    let (r, rho) = (circle.radius(), x.hypot(y));
    let step = (last - first) / n as f64;
    // Samples 2.5 steps from the angle of p are farther than one half a step from it by
    // about 2.5 r ρ step² / distance; required to dwarf the rounding of the distances.
    let distance = (q.dot(q) + r * r).sqrt();
    let scale = (f.origin() - Point::ORIGIN).length() + r + (p - Point::ORIGIN).length();
    let gap = 2.5 * r * rho * step * step / distance.max(1e-300);
    if !(gap > 1e-9 * scale) || (last - first).abs() > 4.0 * std::f64::consts::PI || n < 8 {
        return None;
    }
    let angle = y.atan2(x);
    let tau = std::f64::consts::TAU;
    let mut near = vec![0, 1, 2, n - 2, n - 1, n];
    let lowest = ((first.min(last) - angle) / tau).floor() as i64;
    for k in lowest - 1..=lowest + 3 {
        let at = (angle + tau * k as f64 - first) / step;
        if at > -4.0 && at < n as f64 + 4.0 {
            let centre = at.round() as i64;
            near.extend((centre - 3..=centre + 3).filter(|&i| i >= 0 && i <= n as i64).map(|i| i as usize));
        }
    }
    near.sort_unstable();
    near.dedup();
    Some(near)
}

/// Sample parameters along `[first, last]` of a curve, fine enough to approximate it in
/// parameter space (two for a line).
pub(super) fn sample_parameters(curve: Option<&Curve>, first: f64, last: f64) -> Vec<f64> {
    let n = match curve {
        Some(Curve::Line(_)) => 1,
        Some(Curve::Circle(_)) | Some(Curve::Ellipse(_)) => {
            // About 2° steps.
            (((last - first).abs() / (std::f64::consts::PI / 90.0)).ceil() as usize).max(4)
        }
        _ => 64,
    };
    (0..=n).map(|i| first + (last - first) * i as f64 / n as f64).collect()
}

/// A pcurve sampled over `[first, last]` (its two ends for a line).
pub(super) fn sample_pcurve(pcurve: &Curve2, first: f64, last: f64) -> Vec<Uv> {
    let n = match pcurve {
        Curve2::Line(_) => 1,
        Curve2::Circle(_) | Curve2::Ellipse(_) => (((last - first).abs() / (std::f64::consts::PI / 90.0)).ceil() as usize).max(4),
        _ => 64,
    };
    let params: Vec<f64> = (0..=n).map(|i| first + (last - first) * i as f64 / n as f64).collect();
    pcurve.points_along(&params).into_iter().map(|p| [p.x, p.y]).collect()
}

/// [`sample_pcurve`], with free-form spans halved while their midpoints stray from their
/// chords: for telling points close to the boundary on which side of it they lie.
pub(super) fn sample_pcurve_closely(pcurve: &Curve2, first: f64, last: f64) -> Vec<Uv> {
    if matches!(pcurve, Curve2::Line(_) | Curve2::Circle(_) | Curve2::Ellipse(_)) {
        return sample_pcurve(pcurve, first, last);
    }
    let n = 64;
    let mut params: Vec<f64> = (0..=n).map(|i| first + (last - first) * i as f64 / n as f64).collect();
    let mut points: Vec<Uv> = pcurve.points_along(&params).into_iter().map(|p| [p.x, p.y]).collect();
    let (lo, hi) = polyline_box(&points);
    let tol = 1e-7 * (hi[0] - lo[0]).hypot(hi[1] - lo[1]);
    // (Each round halves the stretches whose middles stray; a stretch whose middle did not
    // stays as it is, and is not looked at again.)
    let mut settled = vec![false; params.len() - 1];
    for _ in 0..6 {
        let open: Vec<usize> = (0..settled.len()).filter(|&i| !settled[i]).collect();
        let mids: Vec<f64> = open.iter().map(|&i| 0.5 * (params[i] + params[i + 1])).collect();
        let at: Vec<Uv> = pcurve.points_along(&mids).into_iter().map(|p| [p.x, p.y]).collect();
        let stray = |k: usize| {
            let i = open[k];
            let (a, b, m) = (points[i], points[i + 1], at[k]);
            let (dx, dy) = (b[0] - a[0], b[1] - a[1]);
            let len = dx.hypot(dy);
            let off = if len > 0.0 { ((m[0] - a[0]) * dy - (m[1] - a[1]) * dx).abs() / len } else { (m[0] - a[0]).hypot(m[1] - a[1]) };
            off > tol
        };
        let strays: Vec<bool> = (0..open.len()).map(stray).collect();
        if !strays.iter().any(|&s| s) {
            break;
        }
        let (mut p2, mut q2, mut s2) = (vec![params[0]], vec![points[0]], Vec::with_capacity(settled.len() * 2));
        let mut k = 0;
        for i in 0..settled.len() {
            if !settled[i] && strays[k] {
                p2.push(mids[k]);
                q2.push(at[k]);
                s2.extend([false, false]);
            } else {
                s2.push(true);
            }
            if !settled[i] {
                k += 1;
            }
            p2.push(params[i + 1]);
            q2.push(points[i + 1]);
        }
        (params, points, settled) = (p2, q2, s2);
    }
    points
}

/// Twice the signed area enclosed by a closed polyline (positive counter-clockwise).
pub(super) fn signed_area(points: &[Uv]) -> f64 {
    let n = points.len();
    (0..n).map(|i| points[i][0] * points[(i + 1) % n][1] - points[(i + 1) % n][0] * points[i][1]).sum::<f64>() / 2.0
}

/// Even-odd point location against closed polylines: true if inside. The ray is cast at
/// an odd angle, so that it does not pass exactly through polyline vertices at "round"
/// coordinates (where nearly coincident ends of adjacent polylines would count twice).
pub(super) fn inside_polylines(polylines: &[Vec<Uv>], p: Uv) -> bool {
    let p = ray_turn(p);
    let mut inside = false;
    for poly in polylines {
        for w in poly.windows(2) {
            if crosses_ray(ray_turn(w[0]), ray_turn(w[1]), p) {
                inside = !inside;
            }
        }
    }
    inside
}

/// A point turned so that [`inside_polylines`]' ray, along the turned u, is not along u or v
/// (where boundaries run).
fn ray_turn(q: Uv) -> Uv {
    let (sin, cos) = 0.618_033_988_7_f64.sin_cos();
    [q[0] * cos + q[1] * sin, q[1] * cos - q[0] * sin]
}

/// Whether the segment `a`-`b` crosses the ray from `p` along u (all turned).
fn crosses_ray(a: Uv, b: Uv, p: Uv) -> bool {
    (a[1] > p[1]) != (b[1] > p[1]) && p[0] < a[0] + (p[1] - a[1]) * (b[0] - a[0]) / (b[1] - a[1])
}

/// [`inside_polylines`] for many points of the same polylines: their segments (turned) in
/// strips across the ray, each point testing only those of its strip (the same segments
/// that can cross its ray, so the same answers).
pub(super) struct InsidePolylines {
    low: f64,
    width: f64,
    strips: Vec<Vec<(Uv, Uv)>>,
}

impl InsidePolylines {
    pub fn new(polylines: &[Vec<Uv>]) -> Self {
        let segments: Vec<(Uv, Uv)> = polylines.iter().flat_map(|poly| poly.windows(2).map(|w| (ray_turn(w[0]), ray_turn(w[1])))).collect();
        let (low, high) = segments.iter().fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), (a, b)| (lo.min(a[1]).min(b[1]), hi.max(a[1]).max(b[1])));
        let count = (segments.len() / 2).clamp(1, 4096);
        let width = if high > low { (high - low) / count as f64 } else { 1.0 };
        let mut strips = vec![Vec::new(); count];
        let strip = |y: f64| (((y - low) / width).floor().max(0.0) as usize).min(count - 1);
        for &(a, b) in &segments {
            for k in strip(a[1].min(b[1]))..=strip(a[1].max(b[1])) {
                strips[k].push((a, b));
            }
        }
        InsidePolylines { low, width, strips }
    }

    pub fn contains(&self, p: Uv) -> bool {
        let p = ray_turn(p);
        if self.strips.is_empty() || !(p[1] >= self.low) {
            return false;
        }
        let k = (((p[1] - self.low) / self.width).floor() as usize).min(self.strips.len() - 1);
        self.strips[k].iter().filter(|(a, b)| crosses_ray(*a, *b, p)).count() % 2 == 1
    }
}

/// A point well inside the region bounded by closed polylines: among the middles of the
/// spans inside the region along horizontal lines halfway between consecutive vertex
/// heights, the one farthest from the boundary.
pub(super) fn interior_point(polylines: &[Vec<Uv>]) -> Uv {
    interior_point_avoiding(polylines, &[])
}

/// Like [`interior_point`], also keeping clear of `obstacles` (polylines inside the region,
/// such as arcs that bound nothing there).
pub(super) fn interior_point_avoiding(polylines: &[Vec<Uv>], obstacles: &[&[Uv]]) -> Uv {
    interior_point_among(polylines, &Obstacles::new(obstacles))
}

/// The box of a polyline in parameter space.
pub(super) fn polyline_box(points: &[Uv]) -> (Uv, Uv) {
    points.iter().fold(([f64::INFINITY; 2], [f64::NEG_INFINITY; 2]), |(lo, hi), p| {
        ([lo[0].min(p[0]), lo[1].min(p[1])], [hi[0].max(p[0]), hi[1].max(p[1])])
    })
}

/// Polylines in parameter space, with a tree over their boxes to find those meeting a
/// region.
pub(super) struct Obstacles<'a> {
    polylines: Vec<&'a [Uv]>,
    tree: Bvh,
}

impl<'a> Obstacles<'a> {
    pub fn new(polylines: &[&'a [Uv]]) -> Self {
        let boxes: Vec<BoundingBox> = polylines
            .iter()
            .map(|o| {
                let (lo, hi) = polyline_box(o);
                BoundingBox::new(Point::new(lo[0], lo[1], 0.0), Point::new(hi[0], hi[1], 0.0))
            })
            .collect();
        Obstacles { polylines: polylines.to_vec(), tree: Bvh::new(&boxes) }
    }

    /// Those meeting the box `(lo, hi)`, in their order.
    fn meeting(&self, (lo, hi): (Uv, Uv)) -> Vec<&'a [Uv]> {
        let query = BoundingBox::new(Point::new(lo[0], lo[1], 0.0), Point::new(hi[0], hi[1], 0.0));
        let mut found = self.tree.overlapping(&query);
        found.sort_unstable();
        found.into_iter().map(|i| self.polylines[i]).collect()
    }
}

/// [`interior_point_avoiding`], with the obstacles indexed.
pub(super) fn interior_point_among(polylines: &[Vec<Uv>], obstacles: &Obstacles) -> Uv {
    // Only obstacles meeting the region's box can cut it.
    let region_box = polylines.iter().fold(([f64::INFINITY; 2], [f64::NEG_INFINITY; 2]), |(lo, hi), p| {
        let (a, b) = polyline_box(p);
        ([lo[0].min(a[0]), lo[1].min(a[1])], [hi[0].max(b[0]), hi[1].max(b[1])])
    });
    let obstacles = obstacles.meeting(region_box);
    // The boundary's segments, then the obstacles': one list, seen in its two parts.
    let mut segments: Vec<(Uv, Uv)> = polylines.iter().flat_map(|p| p.windows(2).map(|w| (w[0], w[1]))).collect();
    let boundary_count = segments.len();
    segments.extend(obstacles.iter().flat_map(|o| o.windows(2).map(|w| (w[0], w[1]))));
    let (boundary, blocking) = segments.split_at(boundary_count);
    // Distances to the segments: by going through them for the first few points, through a
    // grid (worth building for many) after.
    let grid: std::cell::OnceCell<SegmentGrid> = std::cell::OnceCell::new();
    let queries = std::cell::Cell::new(0usize);
    let clearance = |p: Uv| {
        queries.set(queries.get() + 1);
        if queries.get() <= 16 {
            return segments.iter().map(|&s| segment_distance(s, p)).fold(f64::INFINITY, f64::min);
        }
        grid.get_or_init(|| SegmentGrid::new(&segments)).distance(&segments, p)
    };
    // Where a horizontal line crosses segments, sorted (into a buffer used line after line).
    let crossings = |list: &[(Uv, Uv)], y: f64, xs: &mut Vec<f64>| {
        xs.clear();
        xs.extend(list.iter().filter(|(a, b)| (a[1] > y) != (b[1] > y)).map(|(a, b)| a[0] + (y - a[1]) * (b[0] - a[0]) / (b[1] - a[1])));
        xs.sort_by(f64::total_cmp);
    };
    let (lo, hi) = (region_box.0[1].min(region_box.1[1]), region_box.1[1]);
    let size = (region_box.1[0] - region_box.0[0]).hypot(region_box.1[1] - region_box.0[1]);
    let mut best = ([0.0, lo], -1.0);
    let (mut xs, mut cuts_all, mut cuts) = (Vec::new(), Vec::new(), Vec::new());
    let mut try_line = |y: f64, best: &mut (Uv, f64)| {
        crossings(boundary, y, &mut xs);
        crossings(blocking, y, &mut cuts_all);
        for c in xs.chunks_exact(2) {
            if 0.5 * (c[1] - c[0]) <= best.1 {
                continue;
            }
            // Obstacles may cut the span: try the middles of its pieces.
            cuts.clear();
            cuts.push(c[0]);
            cuts.extend(cuts_all.iter().copied().filter(|&x| x > c[0] && x < c[1]));
            cuts.push(c[1]);
            for w in cuts.windows(2) {
                let p = [0.5 * (w[0] + w[1]), y];
                if 0.5 * (w[1] - w[0]) <= best.1 {
                    continue;
                }
                let d = clearance(p);
                if d > best.1 {
                    *best = (p, d);
                }
            }
        }
    };
    // A few evenly spaced lines first, middle ones first (the order of the bits reversed),
    // stopping at a point well clear of the boundary (any such point will do): on most
    // regions they find one; the lines through every slab (between consecutive vertex
    // heights) are needed only for narrow ones.
    for k in 1..32u32 {
        let i = k.reverse_bits() >> 27;
        try_line(lo + (hi - lo) * (i as f64 + 0.123_456_7) / 32.0, &mut best);
        if best.1 >= 0.05 * size {
            return best.0;
        }
    }
    if best.1 >= 1e-3 * size {
        return best.0;
    }
    let mut ys: Vec<f64> = polylines.iter().flatten().map(|p| p[1]).collect();
    ys.extend(obstacles.iter().flat_map(|o| o.iter().map(|p| p[1])));
    ys.sort_by(f64::total_cmp);
    ys.dedup();
    for w in ys.windows(2) {
        try_line(0.5 * (w[0] + w[1]), &mut best);
    }
    best.0
}

/// The distance from `p` to the segment `ab`.
fn segment_distance((a, b): (Uv, Uv), p: Uv) -> f64 {
    let (ab, ap) = ([b[0] - a[0], b[1] - a[1]], [p[0] - a[0], p[1] - a[1]]);
    let len2 = ab[0] * ab[0] + ab[1] * ab[1];
    let t = if len2 > 0.0 { ((ap[0] * ab[0] + ap[1] * ab[1]) / len2).clamp(0.0, 1.0) } else { 0.0 };
    let (dx, dy) = (ap[0] - ab[0] * t, ap[1] - ab[1] * t);
    (dx * dx + dy * dy).sqrt()
}

/// Segments bucketed in a uniform grid of square cells, for nearest-segment distances (or
/// tested one by one when they are few).
struct SegmentGrid {
    origin: Uv,
    cell: f64,
    /// Cells per side; none for few segments.
    n: [i64; 2],
    /// Cell `c`'s segments are `segments[starts[c]..starts[c + 1]]` (indices into the
    /// segment list).
    starts: Vec<usize>,
    segments: Vec<usize>,
}

impl SegmentGrid {
    fn new(segments: &[(Uv, Uv)]) -> Self {
        if segments.len() < 64 {
            return SegmentGrid { origin: [0.0, 0.0], cell: 0.0, n: [0, 0], starts: Vec::new(), segments: Vec::new() };
        }
        let (lo, hi) = segments.iter().fold(([f64::INFINITY; 2], [f64::NEG_INFINITY; 2]), |(lo, hi), (a, b)| {
            ([lo[0].min(a[0]).min(b[0]), lo[1].min(a[1]).min(b[1])], [hi[0].max(a[0]).max(b[0]), hi[1].max(a[1]).max(b[1])])
        });
        let extent = (hi[0] - lo[0]).max(hi[1] - lo[1]).max(1e-12);
        // About two segments per cell along a side.
        let per_side = ((segments.len() as f64).sqrt() * 0.5).ceil().clamp(1.0, 512.0);
        let cell = extent / per_side;
        let n = [((hi[0] - lo[0]) / cell).floor() as i64 + 1, ((hi[1] - lo[1]) / cell).floor() as i64 + 1];
        let at = |x: f64, o: f64, m: i64| (((x - o) / cell).floor() as i64).clamp(0, m - 1);
        let cells_of = |&(a, b): &(Uv, Uv)| {
            let (i0, i1) = (at(a[0].min(b[0]), lo[0], n[0]), at(a[0].max(b[0]), lo[0], n[0]));
            let (j0, j1) = (at(a[1].min(b[1]), lo[1], n[1]), at(a[1].max(b[1]), lo[1], n[1]));
            (i0..=i1).flat_map(move |i| (j0..=j1).map(move |j| (i * n[1] + j) as usize))
        };
        // Counted, then filled: two arrays instead of a list per cell.
        let count = (n[0] * n[1]) as usize;
        let mut starts = vec![0usize; count + 1];
        for s in segments {
            for c in cells_of(s) {
                starts[c + 1] += 1;
            }
        }
        for c in 0..count {
            starts[c + 1] += starts[c];
        }
        let mut next = starts.clone();
        let mut members = vec![0usize; starts[count]];
        for (k, s) in segments.iter().enumerate() {
            for c in cells_of(s) {
                members[next[c]] = k;
                next[c] += 1;
            }
        }
        SegmentGrid { origin: lo, cell, n, starts, segments: members }
    }

    /// The distance from `p` to the nearest segment, searching rings of cells outwards.
    fn distance(&self, segments: &[(Uv, Uv)], p: Uv) -> f64 {
        if self.starts.is_empty() {
            return segments.iter().map(|&s| segment_distance(s, p)).fold(f64::INFINITY, f64::min);
        }
        let (ci, cj) = (((p[0] - self.origin[0]) / self.cell).floor() as i64, ((p[1] - self.origin[1]) / self.cell).floor() as i64);
        let mut best = f64::INFINITY;
        let reach = self.n[0].max(self.n[1]) + ci.abs().max(cj.abs()) + 1;
        for r in 0..=reach {
            // Cells beyond ring r are at least r - 1 cells away.
            if best <= (r as f64 - 1.0).max(0.0) * self.cell {
                break;
            }
            for i in ci - r..=ci + r {
                if i < 0 || i >= self.n[0] {
                    continue;
                }
                for j in cj - r..=cj + r {
                    if j < 0 || j >= self.n[1] || ((i - ci).abs() != r && (j - cj).abs() != r) {
                        continue;
                    }
                    let c = (i * self.n[1] + j) as usize;
                    for &k in &self.segments[self.starts[c]..self.starts[c + 1]] {
                        best = best.min(segment_distance(segments[k], p));
                    }
                }
            }
        }
        best
    }
}

/// The bounding box of points, enlarged by `margin`.
pub(super) fn bounding_box(points: impl IntoIterator<Item = Point>, margin: f64) -> (Point, Point) {
    let mut min = Point::new(f64::INFINITY, f64::INFINITY, f64::INFINITY);
    let mut max = Point::new(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY);
    for p in points {
        min = Point::new(min.x.min(p.x), min.y.min(p.y), min.z.min(p.z));
        max = Point::new(max.x.max(p.x), max.y.max(p.y), max.z.max(p.z));
    }
    (Point::new(min.x - margin, min.y - margin, min.z - margin), Point::new(max.x + margin, max.y + margin, max.z + margin))
}

pub(super) fn boxes_overlap(a: &(Point, Point), b: &(Point, Point)) -> bool {
    !(a.0.x > b.1.x || b.0.x > a.1.x || a.0.y > b.1.y || b.0.y > a.1.y || a.0.z > b.1.z || b.0.z > a.1.z)
}

pub(super) fn box_contains(b: &(Point, Point), p: Point) -> bool {
    p.x >= b.0.x && p.x <= b.1.x && p.y >= b.0.y && p.y <= b.1.y && p.z >= b.0.z && p.z <= b.1.z
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::curve2::BezierCurve2;
    use crate::math::Point2;

    /// A point between an arch and the chord of its 64-part polyline is inside the arch's
    /// closely sampled outline.
    #[test]
    fn points_just_inside_curved_boundaries_are_inside() {
        let arch = Curve2::Bezier(BezierCurve2::new(vec![Point2::new(0.0, 0.0), Point2::new(5.0, 10.0), Point2::new(10.0, 0.0)]));
        let at = |t: f64| {
            let p = arch.point(t);
            [p.x, p.y]
        };
        let (a, b, c) = (at(0.5), at(33.0 / 64.0), at(65.0 / 128.0));
        let q = [c[0] + 0.5 * (0.5 * (a[0] + b[0]) - c[0]), c[1] + 0.5 * (0.5 * (a[1] + b[1]) - c[1])];
        let outline = |arc: Vec<Uv>| vec![arc, vec![[10.0, 0.0], [0.0, 0.0]]];
        assert!(!inside_polylines(&outline(sample_pcurve(&arch, 0.0, 1.0)), q));
        assert!(inside_polylines(&outline(sample_pcurve_closely(&arch, 0.0, 1.0)), q));
    }
}
