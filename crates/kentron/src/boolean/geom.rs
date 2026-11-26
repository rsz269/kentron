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
pub(crate) fn on_tol() -> f64 {
    if crate::tolerance::loosened() { LOOSE_ON_TOL } else { ON_TOL }
}

/// A point in a surface's (u, v) parameter space.
pub(super) type Uv = [f64; 2];

/// A polyline in parameter space, shared (a face's boundary samples, kept from one
/// operation to the next).
pub(super) type Polyline = std::sync::Arc<[Uv]>;

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
            let (foot, su, sv) = surface.d1(uv[0], uv[1]);
            let n = unit_normal(su, sv);
            ((p - foot).dot(n), n)
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
        Surface::Torus(t) => {
            // From the tube's centre circle, less the tube's radius (outward from the tube,
            // as its normal points).
            let f = t.frame();
            let d = p - f.origin();
            let z = Vector::from(f.z());
            let h = d.dot(z);
            let radial = (d - z * h).length();
            (radial - t.major_radius()).hypot(h) - t.minor_radius()
        }
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
            let (foot, su, sv) = surface.d1(uv[0], uv[1]);
            let (r, n) = (p - foot, unit_normal(su, sv));
            if n.length() > 0.5 || matches!(surface, Surface::Revolution(_)) { r.dot(n) } else { r.length() }
        }
        other if freeform(other) => {
            let foot = foot_on(surface, p);
            (p - foot.point).dot(unit_normal(foot.du, foot.dv))
        }
        other => unsupported(other),
    }
}

/// How far the magnitude of [`signed_distance`] to `surface` may stray from a function
/// changing no faster than its point moves (a distance to a set: the surface itself, or
/// curves in the half-plane through the axis): nothing for the elementary surfaces, the
/// tie between a torus's two candidate feet for a surface of revolution of a circle in a
/// plane through its axis clear of it (its closed form; the foot's normal never
/// degenerate). `None` for others, found numerically (a foot may be a local one).
pub(super) fn distance_slack(surface: &Surface) -> Option<f64> {
    match surface {
        Surface::Plane(_) | Surface::Cylinder(_) | Surface::Sphere(_) | Surface::Cone(_) | Surface::Torus(_) => Some(0.0),
        Surface::Revolution(r) => {
            let crate::geometry::curve::Curve::Circle(circle) = r.basis() else { return None };
            let axis = r.axis();
            let (o, z) = (axis.origin, Vector::from(axis.direction));
            let f = circle.frame();
            let normal = Vector::from(f.x()).cross(Vector::from(f.y()));
            let center = circle.center() - o;
            let radial = center - z * center.dot(z);
            // (Where RevolutionSurface::parameters takes it in closed form; within the tie it
            // allows between the two feet.)
            let in_meridian = normal.dot(z).abs() <= 1e-12 && normal.dot(center).abs() <= 1e-12 * (1.0 + center.length());
            (in_meridian && radial.length() > circle.radius()).then(|| 1e-9 * (1.0 + circle.radius()))
        }
        _ => None,
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
    let mut foot = foot_on(surface, p);
    // (Farther than the reach, on a surface coming back near itself, from samples nearer the
    // point too: the foot found may be on the next turn of one winding round, nearer than its
    // samples on the point's own turn, with no sample of that turn near enough to lead there.)
    // (Or on any surface whose sheets pass close to each other, as a ring squashed across
    // its hole, where the nearest sample may lie on the other sheet: the point, taken to
    // lie on the surface, is looked for from the nearest samples before it is put off it.)
    if foot.point.distance(p) > reach {
        foot = foot_from_nearest_samples(surface, p, foot, 32);
    }
    (foot.point.distance(p) <= reach).then_some([foot.u, foot.v])
}

/// The nearest of `foot` and the feet found from the `n` samples nearest `p`.
fn foot_from_nearest_samples(surface: &Surface, p: Point, foot: crate::geometry::extrema::SurfaceFoot, n: usize) -> crate::geometry::extrema::SurfaceFoot {
    let starts: Vec<(f64, f64)> = with_samples(surface, |grid| {
        let mut ranked: Vec<(f64, usize)> = grid.samples.iter().enumerate().map(|(i, s)| (s.2.distance(p), i)).collect();
        ranked.sort_by(|a, b| a.0.total_cmp(&b.0));
        ranked.iter().take(n).map(|&(_, i)| (grid.samples[i].0, grid.samples[i].1)).collect()
    });
    let mut nearest = foot;
    for start in starts {
        let other = crate::geometry::extrema::newton_foot_from_sample(p, surface, start);
        if other.point.distance(p) < nearest.point.distance(p) {
            nearest = other;
        }
    }
    nearest
}

/// The foot of `p` on a B-spline strip kept with the segment and screw it is swept from
/// ([`crate::geometry::screw_strip`]): Newton steps on the strip from the point of their
/// helicoid nearest `p` (within the strip's parameters), which lies within about 1e-10 of the
/// size of the strip's foot. `None` for other surfaces, or where the steps come out farther
/// than that point (then the general search is taken).
fn screw_strip_foot(surface: &Surface, p: Point) -> Option<crate::geometry::extrema::SurfaceFoot> {
    let Surface::BSpline(b) = surface else { return None };
    let strip = b.screw_strip()?;
    let (u0, u1, v0, v1) = surface.bounds();
    let (u, v, d) = strip.nearest(p, (u0, u1), (v0, v1));
    // (One Newton step on the strip from there, about 1e-10 of the size off its foot: then
    // as near as rounding allows.)
    let at = surface.d2(u, v);
    let r = at.point - p;
    let (gu, gv) = (r.dot(at.du), r.dot(at.dv));
    let (huu, huv, hvv) = (at.du.dot(at.du) + r.dot(at.duu), at.du.dot(at.dv) + r.dot(at.duv), at.dv.dot(at.dv) + r.dot(at.dvv));
    let det = huu * hvv - huv * huv;
    let (mut nu, mut nv) = (u, v);
    if det > 0.0 && huu > 0.0 {
        nu = (u + (-gu * hvv + gv * huv) / det).clamp(u0, u1);
        nv = (v + (-gv * huu + gu * huv) / det).clamp(v0, v1);
    }
    let (point, du, dv) = surface.d1(nu, nv);
    let foot = if point.distance(p) <= at.point.distance(p) {
        crate::geometry::extrema::SurfaceFoot { u: nu, v: nv, point, du, dv }
    } else {
        crate::geometry::extrema::SurfaceFoot { u, v, point: at.point, du: at.du, dv: at.dv }
    };
    let o = strip.axis.origin;
    let slack = 1e-8 * (1.0 + strip.start.distance(o) + strip.end.distance(o));
    (foot.point.distance(p) <= d + slack).then_some(foot)
}

/// The foot of `p` on a free-form surface: Newton steps from the nearest of a grid of the
/// surface's points. Where that leads astray (to a foot on the surface's bounds, or farther
/// than the sample it started from) and another part of the surface passes nearly as near
/// (the next turn of one winding round), from the nearest sample there too: the nearer foot.
pub(super) fn foot_on(surface: &Surface, p: Point) -> crate::geometry::extrema::SurfaceFoot {
    // (A strip a segment sweeps along a screw: from the exact surface's nearest point.)
    if let Some(foot) = screw_strip_foot(surface, p) {
        return foot;
    }
    let (best, start) = with_samples(surface, |grid| {
        let best = grid.nearest_of_all(p);
        (best, (grid.samples[best.1].0, grid.samples[best.1].1))
    });
    let mut foot = crate::geometry::extrema::newton_foot_from_sample(p, surface, start);
    // (A strip ruled one way (a B-spline of one span of degree 1 that way, of several spans
    // the other), whose rulings turn along it, as a straight edge carried along a screw
    // sweeps: its few samples across the rulings may lie far from the foot, and the steps
    // from the nearest settle on another, farther one, off a point lying on the strip. Where
    // the foot found is nearer than a ruling is long, from the point of the rulings, taken
    // closely along the strip, nearest it (each ruling's exactly, a segment's): the nearer
    // foot.)
    if foot.point.distance(p) > on_tol() {
        if let Some(rulings) = rulings(surface) {
            let d = foot.point.distance(p);
            if d < rulings.length {
                let nearest = rulings.nearest(p);
                if nearest.0 < d {
                    let second = crate::geometry::extrema::newton_foot_from_sample(p, surface, nearest.1);
                    if second.point.distance(p) < d {
                        foot = second;
                    }
                }
            }
        }
    }
    // (Off the surface at a collapsed row, a pole or an apex where a free-form surface's
    // every u (or v) is one point, or from a sample there: the steps may stay at it, the
    // distance not changing along the row, or run into it from the side away from the foot.
    // From the nearest of a finer row of samples just beside it too.)
    if freeform(surface) && foot.point.distance(p) > on_tol() {
        let (u0, u1, v0, v1) = surface.bounds();
        let collapsed = |u: f64, v: f64| {
            let (_, su, sv) = surface.d1(u, v);
            let (lu, lv) = (su.length(), sv.length());
            (lu <= 1e-9 * lv, lv <= 1e-9 * lu)
        };
        let (at_start, at_foot) = (collapsed(start.0, start.1), collapsed(foot.u, foot.v));
        let (along_u, along_v) = (at_start.0 || at_foot.0, at_start.1 || at_foot.1);
        if (along_u || along_v) && [u0, u1, v0, v1].iter().all(|x| x.is_finite()) {
            let (nu, nv) = grid_steps(surface);
            // The collapsed row's parameter, and rows a fraction of a grid step in from it.
            let (row, (lo, hi), steps, from) = if along_u {
                (if at_foot.0 { foot.v } else { start.1 }, (v0, v1), nv, (u0, u1))
            } else {
                (if at_foot.1 { foot.u } else { start.0 }, (u0, u1), nu, (v0, v1))
            };
            let step = (hi - lo) / steps.max(1) as f64;
            let mut nearest = (f64::INFINITY, start);
            for k in [1.0, 0.5, 0.25, 0.125] {
                for side in [-1.0, 1.0] {
                    let w = row + side * k * step;
                    if w < lo || w > hi {
                        continue;
                    }
                    for i in 0..=64 {
                        let t = from.0 + (from.1 - from.0) * i as f64 / 64.0;
                        let (u, v) = if along_u { (t, w) } else { (w, t) };
                        let d = surface.point(u, v).distance(p);
                        if d < nearest.0 {
                            nearest = (d, (u, v));
                        }
                    }
                }
            }
            if nearest.0.is_finite() {
                let steps = if along_u { ((from.1 - from.0) / 64.0, step) } else { (step, (from.1 - from.0) / 64.0) };
                let second = narrowed(surface, p, nearest.1, steps);
                if second.point.distance(p) < foot.point.distance(p) {
                    foot = second;
                }
            }
        }
    }
    // (No nearer than the nearest sample itself, the steps having run off or stalled, as
    // from beside a pole, or along a bound to a farther point: from that sample by steps
    // each bringing the point nearer.)
    if freeform(surface) && foot.point.distance(p) > on_tol() && foot.point.distance(p) > best.0.sqrt() * (1.0 + 1e-9) + 1e-12 {
        let (u0, u1, v0, v1) = surface.bounds();
        if [u0, u1, v0, v1].iter().all(|x| x.is_finite()) {
            let second = descended(surface, p, start);
            if second.point.distance(p) < foot.point.distance(p) {
                foot = second;
            }
        }
    }
    // (From a sample far off compared with the foot found, as on a long narrow strip
    // sampled only by its spans, the steps may settle off the true foot: the general search
    // too, kept where it finds the point on the surface.)
    if foot.point.distance(p) > on_tol() && best.0.sqrt() > 2.0 * foot.point.distance(p) {
        let (u, v, d) = crate::geometry::extrema::nearest_on_surface(p, surface, Some(start));
        if d <= on_tol() {
            let (point, du, dv) = surface.d1(u, v);
            foot = crate::geometry::extrema::SurfaceFoot { u, v, point, du, dv };
        }
    }
    let (u0, u1, v0, v1) = surface.bounds();
    let on_bound = |x: f64, a: f64, b: f64, periodic: bool| !periodic && ((x - a).abs() <= 1e-12 * (b - a).abs().max(1.0) || (x - b).abs() <= 1e-12 * (b - a).abs().max(1.0));
    let distance = foot.point.distance(p);
    // (Or a foot just off the point, far nearer than the samples are: where a surface winds
    // back onto itself, as a thread whose flanks touch from turn to turn, the point may lie
    // on the next turn, its nearest sample on this one.)
    let grazing = distance > on_tol() && distance < 1e-3 * best.0.sqrt();
    let (u_bound, v_bound) = (on_bound(foot.u, u0, u1, surface.u_period().is_some()), on_bound(foot.v, v0, v1, surface.v_period().is_some()));
    // On a bound of a surface closed there without a period (a B-spline round a whole
    // turn), the steps may have run into the seam from the side away from the foot: from
    // the seam's other copy too.
    if (u_bound && surface.is_u_closed()) || (v_bound && surface.is_v_closed()) {
        let across = |x: f64, a: f64, b: f64| if (x - a).abs() < (x - b).abs() { b } else { a };
        let start = (if u_bound && surface.is_u_closed() { across(foot.u, u0, u1) } else { foot.u }, if v_bound && surface.is_v_closed() { across(foot.v, v0, v1) } else { foot.v });
        let (nu, nv) = grid_steps(surface);
        let second = narrowed(surface, p, start, ((u1 - u0) / nu.max(1) as f64, (v1 - v0) / nv.max(1) as f64));
        if second.point.distance(p) < distance {
            return second;
        }
    }
    let astray = distance > best.0.sqrt() * (1.0 + 1e-9) + 1e-12 || grazing || u_bound || v_bound;
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
    let nearer = match other {
        Some(other) => {
            let second = crate::geometry::extrema::newton_foot_from_sample(p, surface, other);
            if second.point.distance(p) < distance { second } else { foot }
        }
        None => foot,
    };
    // (Grazing still, from more of the samples nearest the point: the next turn's may lie
    // farther than those of the turn the point is not on, as [`uv_within`] finds them.)
    if grazing && nearer.point.distance(p) > on_tol() {
        return foot_from_nearest_samples(surface, p, nearer, 32);
    }
    nearer
}

/// The point of `surface` nearest `p` about `start`, narrowed down by ever finer grids of
/// 5 × 5 points round the nearest found (from steps of `steps`, halved each time), then
/// Newton's steps from there where they come nearer still: slow, but sure where the steps
/// alone run off (near a pole, where u barely moves the point).
fn narrowed(surface: &Surface, p: Point, start: (f64, f64), steps: (f64, f64)) -> crate::geometry::extrema::SurfaceFoot {
    let (u0, u1, v0, v1) = surface.bounds();
    let (mut du, mut dv) = steps;
    let mut nearest = (surface.point(start.0, start.1).distance(p), start);
    for _ in 0..40 {
        let (cu, cv) = nearest.1;
        for i in -2..=2 {
            for j in -2..=2 {
                let (u, v) = ((cu + i as f64 * du).clamp(u0, u1), (cv + j as f64 * dv).clamp(v0, v1));
                let d = surface.point(u, v).distance(p);
                if d < nearest.0 {
                    nearest = (d, (u, v));
                }
            }
        }
        du *= 0.5;
        dv *= 0.5;
    }
    let (u, v) = nearest.1;
    let (point, su, sv) = surface.d1(u, v);
    let found = crate::geometry::extrema::SurfaceFoot { u, v, point, du: su, dv: sv };
    let polished = crate::geometry::extrema::newton_foot_from_sample(p, surface, nearest.1);
    if polished.point.distance(p) < found.point.distance(p) { polished } else { found }
}

/// The point of `surface` nearest `p` about `start` (on a surface with finite bounds), by
/// steps that each bring the point nearer: Newton's where the distance is convex, else down
/// the gradient, either halved until it comes nearer. A parameter on its bound that the
/// distance would take past it is held there, the other alone moving (the foot sought
/// along that edge); where both are held, or no step comes nearer, the foot is there. (As
/// sure as [`narrowed`] from a sample where Newton's own steps run off, a few evaluations
/// a step instead of a thousand.)
fn descended(surface: &Surface, p: Point, start: (f64, f64)) -> crate::geometry::extrema::SurfaceFoot {
    let (u0, u1, v0, v1) = surface.bounds();
    let (u_periodic, v_periodic) = (surface.is_u_periodic(), surface.is_v_periodic());
    let clamp = |u: f64, v: f64| (if u_periodic { u } else { u.clamp(u0, u1) }, if v_periodic { v } else { v.clamp(v0, v1) });
    let (mut u, mut v) = clamp(start.0, start.1);
    let (u_small, v_small) = (1e-12 * (u1 - u0).abs().max(1e-300), 1e-12 * (v1 - v0).abs().max(1e-300));
    let mut f = surface.point(u, v).distance_squared(p);
    for _ in 0..100 {
        let d = surface.d2(u, v);
        let r = d.point - p;
        let (gu, gv) = (r.dot(d.du), r.dot(d.dv));
        let (huu, huv, hvv) = (d.du.dot(d.du) + r.dot(d.duu), d.du.dot(d.dv) + r.dot(d.duv), d.dv.dot(d.dv) + r.dot(d.dvv));
        let held = |x: f64, g: f64, lo: f64, hi: f64, periodic: bool| !periodic && ((x <= lo && g > 0.0) || (x >= hi && g < 0.0));
        let (hu, hv) = (held(u, gu, u0, u1, u_periodic), held(v, gv, v0, v1, v_periodic));
        let (lu, lv) = (d.du.dot(d.du).max(1e-300), d.dv.dot(d.dv).max(1e-300));
        let det = huu * hvv - huv * huv;
        let newton = match (hu, hv) {
            (true, true) => break,
            (true, false) => (0.0, if hvv > 0.0 { -gv / hvv } else { -gv / lv }),
            (false, true) => (if huu > 0.0 { -gu / huu } else { -gu / lu }, 0.0),
            (false, false) if huu > 0.0 && det > 0.0 => ((-gu * hvv + gv * huv) / det, (-gv * huu + gu * huv) / det),
            (false, false) => (-gu / lu.max(lv), -gv / lu.max(lv)),
        };
        let gradient = (if hu { 0.0 } else { -gu / lu.max(lv) }, if hv { 0.0 } else { -gv / lu.max(lv) });
        let mut moved = None;
        'directions: for (su, sv) in [newton, gradient] {
            let mut t = 1.0;
            for _ in 0..40 {
                let (nu, nv) = clamp(u + t * su, v + t * sv);
                if (nu - u).abs() <= u_small && (nv - v).abs() <= v_small {
                    break;
                }
                let fc = surface.point(nu, nv).distance_squared(p);
                if fc < f {
                    moved = Some((nu, nv, fc));
                    break 'directions;
                }
                t *= 0.5;
            }
        }
        // (Steps of a trillionth of the range are not taken: the foot is there.)
        let Some((nu, nv, fc)) = moved else { break };
        (u, v, f) = (nu, nv, fc);
    }
    let (point, du, dv) = surface.d1(u, v);
    crate::geometry::extrema::SurfaceFoot { u, v, point, du, dv }
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
        // (Any offset surface (one not the elementary surface it equals): numerically, as
        // free-form; a surface of revolution's or extrusion's too.)
        Surface::Offset(_) => true,
        Surface::Trimmed(t) => freeform(t.basis()),
        _ => false,
    }
}

/// Grids of points of the free-form surfaces met lately, to start point inversions from
/// (shared by the threads of an operation).
static SAMPLES: std::sync::RwLock<Option<crate::hashing::FastMap<crate::geometry::extrema::SurfaceKey, std::sync::Arc<SampleGrid>>>> = std::sync::RwLock::new(None);

/// A surface's sample grid: its points `(u, v, p)`, row by row along u, and their
/// coordinates again in arrays of their own, for scans of all of them.
struct SampleGrid {
    samples: Vec<(f64, f64, Point)>,
    x: Vec<f64>,
    y: Vec<f64>,
    z: Vec<f64>,
    /// The boxes of the samples in blocks of [`BLOCK`], in order.
    boxes: Vec<(Point, Point)>,
    /// The block the last scan found the nearest sample in, scanned first.
    last: std::sync::atomic::AtomicUsize,
}

/// The samples in a block of [`SampleGrid::boxes`].
const BLOCK: usize = 32;

impl SampleGrid {
    fn new(samples: Vec<(f64, f64, Point)>) -> SampleGrid {
        let x = samples.iter().map(|s| s.2.x).collect();
        let y = samples.iter().map(|s| s.2.y).collect();
        let z = samples.iter().map(|s| s.2.z).collect();
        let boxes = samples
            .chunks(BLOCK)
            .map(|block| {
                block.iter().fold((Point::new(f64::INFINITY, f64::INFINITY, f64::INFINITY), Point::new(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY)), |(lo, hi), s| {
                    let q = s.2;
                    (Point::new(lo.x.min(q.x), lo.y.min(q.y), lo.z.min(q.z)), Point::new(hi.x.max(q.x), hi.y.max(q.y), hi.z.max(q.z)))
                })
            })
            .collect();
        SampleGrid { samples, x, y, z, boxes, last: std::sync::atomic::AtomicUsize::new(0) }
    }

    /// [`Self::nearest`] among all the samples, the same one found: the block of the last
    /// one found scanned first, then the others in order, those whose box lies farther than
    /// the nearest sample so far passed over (none of theirs could be as near). (A box
    /// with a NaN corner is never passed over.)
    fn nearest_of_all(&self, p: Point) -> (f64, usize) {
        use std::sync::atomic::Ordering;
        let blocks = self.boxes.len();
        let first = self.last.load(Ordering::Relaxed).min(blocks.saturating_sub(1));
        let scan = |b: usize| self.nearest(p, b * BLOCK, ((b + 1) * BLOCK).min(self.samples.len()));
        let mut best = if blocks > 0 { scan(first) } else { (f64::INFINITY, 0) };
        let mut found = first;
        for b in (0..blocks).filter(|&b| b != first) {
            let (lo, hi) = self.boxes[b];
            let gap = |a: f64, l: f64, h: f64| if a < l { l - a } else if a > h { a - h } else { 0.0 };
            let (gx, gy, gz) = (gap(p.x, lo.x, hi.x), gap(p.y, lo.y, hi.y), gap(p.z, lo.z, hi.z));
            if gx * gx + gy * gy + gz * gz > best.0 {
                continue;
            }
            let here = scan(b);
            if here.0 < best.0 || (here.0 == best.0 && here.0 < f64::INFINITY && here.1 < best.1) {
                (best, found) = (here, b);
            }
        }
        if best.0 < f64::INFINITY {
            self.last.store(found, Ordering::Relaxed);
        }
        best
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
    // ([`crate::geometry::extrema::windings_from`]), or coming back near itself, two in each
    // gap between its turns ([`crate::geometry::extrema::close_samples`]).
    let by = |spans: usize, degree: usize| (4 * spans * degree.max(1)).clamp(4, 32);
    match surface {
        Surface::BSpline(b) => {
            let (wu, wv) = crate::geometry::extrema::windings_from(surface, 2.0);
            let (cu, cv) = crate::geometry::extrema::close_samples(surface);
            (by(b.u_knots().len() - 1, b.u_degree()).max(wu).max(cu), by(b.v_knots().len() - 1, b.v_degree()).max(wv).max(cv))
        }
        Surface::Bezier(b) => (by(1, b.u_degree()), by(1, b.v_degree())),
        _ => (32, 32),
    }
}


/// The rulings of a B-spline strip ruled one way (one span of degree 1 that way, several spans
/// the other way), taken closely along it (a single patch, as a bilinear one, may fold over
/// itself, and a foot found on the far fold is not the one its sections follow); `None` for
/// other surfaces. Kept with the surface.
fn rulings(surface: &Surface) -> Option<&crate::geometry::Rulings> {
    let Surface::BSpline(b) = surface else { return None };
    b.shared()
        .rulings
        .get_or_init(|| {
            let spans = |knots: &[f64]| knots.len().saturating_sub(1);
            let (u0, u1, v0, v1) = surface.bounds();
            if ![u0, u1, v0, v1].iter().all(|x| x.is_finite()) {
                return None;
            }
            let along_u = b.u_degree() == 1 && spans(b.u_knots()) == 1 && spans(b.v_knots()) >= 2;
            let along_v = b.v_degree() == 1 && spans(b.v_knots()) == 1 && spans(b.u_knots()) >= 2;
            if along_u == along_v {
                return None;
            }
            // (Eight to each span and degree of the other way, at most 1024.)
            let n = if along_u { 8 * spans(b.v_knots()) * b.v_degree().max(1) } else { 8 * spans(b.u_knots()) * b.u_degree().max(1) }.clamp(16, 1024);
            let lines: Vec<(f64, Point, Point)> = (0..=n)
                .map(|k| {
                    if along_u {
                        let v = v0 + (v1 - v0) * k as f64 / n as f64;
                        (v, surface.point(u0, v), surface.point(u1, v))
                    } else {
                        let u = u0 + (u1 - u0) * k as f64 / n as f64;
                        (u, surface.point(u, v0), surface.point(u, v1))
                    }
                })
                .collect();
            Some(crate::geometry::Rulings::new(along_u, if along_u { (u0, u1) } else { (v0, v1) }, lines))
        })
        .as_ref()
}

/// The samples of a surface winding round twice or more, its turns close together
/// ([`crate::geometry::extrema::windings_from`]), or coming back near itself
/// ([`crate::geometry::extrema::close_samples`]), along its more winding direction; 0 for
/// others.
pub(super) fn winding_samples(surface: &Surface) -> usize {
    let (u, v) = crate::geometry::extrema::windings_from(surface, 2.0);
    let (cu, cv) = crate::geometry::extrema::close_samples(surface);
    u.max(v).max(cu).max(cv)
}

/// Runs `f` on the surface's sample grid (points `(u, v, p)`, row by row along u), cached.
fn with_samples<R>(surface: &Surface, f: impl FnOnce(&SampleGrid) -> R) -> R {
    let key = crate::geometry::extrema::SurfaceKey::of(surface);
    let found = SAMPLES.read().unwrap_or_else(|poisoned| poisoned.into_inner()).as_ref().and_then(|cache| cache.get(&key).cloned());
    let grid = match found {
        Some(grid) => grid,
        None => {
            // (Made outside the lock; another thread may make it too.)
            let grid = std::sync::Arc::new(SampleGrid::new(sample_grid(surface)));
            let mut cache = SAMPLES.write().unwrap_or_else(|poisoned| poisoned.into_inner());
            let cache = cache.get_or_insert_with(Default::default);
            if cache.len() >= 512 {
                cache.clear();
            }
            cache.entry(key).or_insert_with(|| grid.clone());
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

/// True at a pole of a sphere, the apex of a cone or where a surface of revolution's
/// meridian meets its axis (where `u` is undefined).
pub(super) fn is_pole(surface: &Surface, p: Point) -> bool {
    match surface {
        Surface::Revolution(r) => {
            let axis = r.axis();
            let d = p - axis.origin;
            let z = Vector::from(axis.direction);
            (d - z * d.dot(z)).length() <= on_tol()
        }
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

/// [`normal_at_point`] at a point of the surface found by other means (a ray's hit on it):
/// on a free-form surface whose foot of the point is found off it (on another turn of a
/// surface winding round, the steps from the nearest of its samples led astray), the foot
/// from the samples nearest the point too, the nearer taken.
pub(super) fn normal_at_point_on(surface: &Surface, p: Point) -> Vector {
    if !freeform(surface) || surface.parameters(p).is_some() {
        return normal_at_point(surface, p);
    }
    let mut foot = foot_on(surface, p);
    if foot.point.distance(p) > on_tol() {
        foot = foot_from_nearest_samples(surface, p, foot, 32);
    }
    unit_normal(foot.du, foot.dv)
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
                    // (From the middle of the corners' reach along the ray, and from either end
                    // of it: a ray in and out again across one cell, as by a rim, meets it
                    // twice.)
                    let along: Vec<f64> = corners.iter().map(|c| (c.2 - origin).dot(d)).collect();
                    let (first, last) = along.iter().fold((f64::INFINITY, f64::NEG_INFINITY), |(a, b), &t| (a.min(t), b.max(t)));
                    let (u, v) = (0.5 * (corners[0].0 + corners[3].0), 0.5 * (corners[0].1 + corners[3].1));
                    starts.push((u, v, along.iter().sum::<f64>() / 4.0));
                    starts.push((u, v, first - margin));
                    starts.push((u, v, last + margin));
                    // (And from its corners: two hits in one cell, Newton's steps from its
                    // middle all may run to the one.)
                    for (c, &t) in corners.iter().zip(&along) {
                        starts.push((c.0, c.1, t));
                    }
                }
            }
        }
        starts
    });
    let mut hits: Vec<f64> = Vec::new();
    let cell = ((u1 - u0) / nu.max(1) as f64, (v1 - v0) / nv.max(1) as f64);
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
            let (du, dv, dt) = (-f.dot(sv.cross(-d)) / det, -su.dot(f.cross(-d)) / det, -su.dot(sv.cross(f)) / det);
            // (Each step at most a cell's size: from a cell's middle, the full step may leap
            // to another cell's hit, the one in this cell never found.)
            let k = (cell.0 / du.abs().max(1e-300)).min(cell.1 / dv.abs().max(1e-300)).min(1.0);
            let (du, dv, dt) = (k * du, k * dv, k * dt);
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
    // (Sampled along the ray, past the samples that cannot be on the surface: the distance
    // shrinks by no more than the way along the ray, half of it kept for a foot found short
    // of the nearest.)
    let reach = (t1 - t0) / n as f64 * direction.length();
    let (mut ts, mut fs) = (Vec::new(), Vec::new());
    let mut k = 0;
    while k <= n {
        let t = t0 + (t1 - t0) * k as f64 / n as f64;
        let d = f(t);
        ts.push(t);
        fs.push(d);
        let clear = if reach > 0.0 { (0.5 * d.abs() / reach).floor() } else { 0.0 };
        k = if k == n { n + 1 } else { (k + (clear.min(n as f64) as usize).max(1)).min(n) };
    }
    let n = ts.len() - 1;
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
            // (A sample on the surface, or next to it, between samples on either side of
            // it: a crossing, not a graze. One on it exactly is a hit there; one next to it
            // is found by the sign change on its other side.)
            if fs[i - 1] * fs[i + 1] < 0.0 {
                let q = origin + direction * ts[i];
                if fs[i] == 0.0 && surface_point(surface, uv_of(surface, q)).distance(q) <= 10.0 * on_tol() {
                    hits.push(ts[i]);
                }
                continue;
            }
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
    // Conics in closed form: the nearest point within the range is one of the curve's
    // critical points of the distance there, or an end.
    if let Some(found) = conic_foot(curve, first, last, p) {
        return found;
    }
    if let Some(found) = helix_foot(curve, first, last, p) {
        return found;
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
    (0..=n).for_each(&mut sample);
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

/// [`project_among`] on a B-spline kept with the screw it is a helix of
/// ([`crate::geometry::screw_strip::Helix`]): from the point of the exact helix nearest `p`
/// (its parameter taken as turning evenly, as the screw's paths do, within about 1e-10 of
/// the size), Newton steps on the curve. `None` for other curves, or where the steps come out
/// farther than that point (then the general search is taken).
fn helix_foot(curve: &Curve, first: f64, last: f64, p: Point) -> Option<(f64, f64)> {
    let Curve::BSpline(b) = curve else { return None };
    let helix = b.helix()?;
    let k = helix.rise_per_radian;
    let (f0, f1) = (b.first_parameter(), b.last_parameter());
    let (c0, c1) = (b.point(f0), b.point(f1));
    let z = Vector::from(helix.axis.direction);
    if !(k.is_finite() && k != 0.0 && f1 > f0) {
        return None;
    }
    // The turn over the curve's whole range, from the rise between its ends.
    let angle = (c1 - c0).dot(z) / k;
    let rate = angle / (f1 - f0);
    let exact = crate::geometry::screw_strip::ScrewStrip { axis: helix.axis, turn: rate, rise: k * rate, v0: f0, start: c0, end: c0, u0: 0.0, u1: 1.0 };
    let (lo, hi) = (first.min(last), first.max(last));
    let (_, mut t, d) = exact.nearest(p, (0.0, 1.0), (lo, hi));
    for _ in 0..3 {
        let (c, d1, d2) = curve.d2(t);
        let f = (c - p).dot(d1);
        let df = d1.dot(d1) + (c - p).dot(d2);
        if !(df > 1e-300) {
            break;
        }
        let next = (t - f / df).clamp(lo, hi);
        let done = (next - t).abs() <= 1e-15 * (1.0 + t.abs());
        t = next;
        if done {
            break;
        }
    }
    let found = curve.point(t).distance(p);
    let o = helix.axis.origin;
    let slack = 1e-8 * (1.0 + c0.distance(o) + c1.distance(o));
    (found <= d + slack).then_some((t, found))
}

/// [`project_among`] on a circle or an ellipse in closed form: the parameter in
/// `[first, last]` of the point of `curve` nearest to `p`, and the distance; `None` for other
/// curves and degenerate conics.
///
/// Measured in the conic's plane (the height above it adds the same to every point). The
/// nearest point of the whole conic is taken where it lies in the range; elsewhere the
/// nearest is at an end or at another minimum of the distance within the range (a circle has
/// none; an ellipse one where the point is inside its evolute, see [`ellipse_feet`]).
pub(super) fn conic_foot(curve: &Curve, first: f64, last: f64, p: Point) -> Option<(f64, f64)> {
    let (lo, hi) = (first.min(last), first.max(last));
    let tau = std::f64::consts::TAU;
    let (frame, a, b) = match curve {
        Curve::Circle(c) => (c.frame(), c.radius(), c.radius()),
        Curve::Ellipse(e) => (e.frame(), e.major_radius(), e.minor_radius()),
        _ => return None,
    };
    if !(b > 0.0) || !(a.is_finite() && lo.is_finite() && hi.is_finite()) {
        return None;
    }
    let q = p - frame.origin();
    let (x, y) = (q.dot(Vector::from(frame.x())), q.dot(Vector::from(frame.y())));
    // The first turn of an angle at or past the range's start, if within the range.
    let wrap = |angle: f64| {
        let t = lo + (angle - lo).rem_euclid(tau);
        if t <= hi {
            Some(t)
        } else {
            // (Rounding in the wrap.)
            (t - tau >= lo).then_some(t - tau)
        }
    };
    let at = |t: f64| (t, curve.point(t).distance(p));
    let mut feet = [0.0f64; 6];
    let (nearest, count) = if a == b {
        // (On the axis every point is as near: the range's start is taken.)
        if x == 0.0 && y == 0.0 {
            return Some(at(lo));
        }
        feet[0] = y.atan2(x);
        (1, 1)
    } else {
        ellipse_feet(a, b, x, y, false, &mut feet)
    };
    // The nearest of the whole conic, where within the range.
    if let Some(t) = feet[..nearest].iter().find_map(|&angle| wrap(angle)) {
        return Some(at(t));
    }
    let (nearest, count) = if a == b { (nearest, count) } else { ellipse_feet(a, b, x, y, true, &mut feet) };
    let mut best = at(lo);
    for t in feet[nearest..count].iter().filter_map(|&angle| wrap(angle)).chain([hi]) {
        let found = at(t);
        if found.1 < best.1 {
            best = found;
        }
    }
    Some(best)
}

/// The angles `t` of the ellipse `(a cos t, b sin t)` (`a > b > 0`) where the distance to
/// the point `(x, y)` of its plane can be least, written to `out`: first those of the nearest
/// point (two where `y = 0` puts it on either side of the major axis), then the others (a
/// local minimum and the maximum beside it); returns how many are nearest and how many there
/// are in all (only the nearest unless `all`).
///
/// At a foot point `(X, Y)` the gap to the point is along the normal: with a multiplier `s`,
/// `X = a² x / (s + d)` and `Y = b² y / s` where `d = a² - b²`, the foot lying on the ellipse
/// where `F(s) = (a x / (s + d))² + (b y / s)² = 1`. `F` falls from +∞ to 0 over `s > 0` (one
/// root: the nearest point), is convex over `-d < s < 0` (two roots where the point is inside
/// the evolute, `(a x)^⅔ + (b y)^⅔ < d^⅔`, else none) and rises over `s < -d` (the farthest
/// point, not wanted). Each root is reached by Newton steps from a side where `F > 1`, which
/// approach it without overshooting (`F` is convex). Where `y = 0` the roots near `s = 0` run
/// into the pole there: the feet are then `X = a² x / d` with either sign of `Y`.
pub(super) fn ellipse_feet(a: f64, b: f64, x: f64, y: f64, all: bool, out: &mut [f64; 6]) -> (usize, usize) {
    let d = (a - b) * (a + b);
    let (ax, by) = (a * x, b * y);
    let (aa, bb) = (ax * ax, by * by);
    // Newton steps on F(s) = 1 from `s` (where F(s) > 1), towards the root.
    let root = |mut s: f64| {
        for _ in 0..200 {
            let (u, v) = (ax / (s + d), by / s);
            let f = u * u + v * v - 1.0;
            let df = -2.0 * (u * u / (s + d) + v * v / s);
            if !(f > 0.0) || !(df != 0.0) || !df.is_finite() {
                break;
            }
            let next = s - f / df;
            if !next.is_finite() || next == s {
                break;
            }
            s = next;
        }
        s
    };
    let mut n = 0;
    let push = |xx: f64, yy: f64, n: &mut usize, out: &mut [f64; 6]| {
        if xx.is_finite() && yy.is_finite() && *n < out.len() {
            out[*n] = (yy / b).atan2(xx / a);
            *n += 1;
        }
    };
    let foot = |s: f64| (a * a * x / (s + d), b * b * y / s);
    // Where F <= 1 at the pole (y = 0 inside the evolute's span, or the centre): the feet at
    // the pole.
    let start = by.abs().max(ax.abs() - d);
    if start > 0.0 {
        let (xx, yy) = foot(root(start));
        push(xx, yy, &mut n, out);
    } else {
        let xx = (a * a * x / d).clamp(-a, a);
        let yy = b * (1.0 - (xx / a).powi(2)).max(0.0).sqrt();
        push(xx, yy, &mut n, out);
        push(xx, -yy, &mut n, out);
    }
    let nearest = n;
    // Inside the evolute: the other minimum and the maximum beside it.
    if all && aa.cbrt() + bb.cbrt() < d.cbrt().powi(2) {
        if aa > 0.0 {
            // From the left of F's least value (at `s + d = a|x| / 2`, F >= 4).
            let (xx, yy) = foot(root(-d + 0.5 * ax.abs()));
            push(xx, yy, &mut n, out);
        }
        if bb > 0.0 {
            // From its right (at `s = -b|y| / 2`).
            let (xx, yy) = foot(root(-0.5 * by.abs()));
            push(xx, yy, &mut n, out);
        }
    }
    // Polish the others (found less precisely near the poles): Newton steps on
    // g(t) = (E(t) - p)·E'(t).
    let k = b * b - a * a;
    for t in &mut out[nearest..n] {
        for _ in 0..4 {
            let (sn, cs) = t.sin_cos();
            let g = k * sn * cs + ax * sn - by * cs;
            let dg = k * (cs * cs - sn * sn) + ax * cs + by * sn;
            let step = g / dg;
            if !(step.abs() < 0.1) {
                break;
            }
            *t -= step;
            if step.abs() <= 1e-15 {
                break;
            }
        }
    }
    (nearest, n)
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
/// The angle between the samples [`sample_pcurve`] takes along a circle or an ellipse.
const ARC_STEP: f64 = std::f64::consts::PI / 90.0;

/// How far the chords between [`sample_pcurve`]'s samples of an arc of `radius` stray from it.
pub(super) fn arc_sag(radius: f64) -> f64 {
    radius.abs() * (1.0 - (0.5 * ARC_STEP).cos())
}

pub(super) fn sample_pcurve(pcurve: &Curve2, first: f64, last: f64) -> Vec<Uv> {
    let n = match pcurve {
        Curve2::Line(_) => 1,
        Curve2::Circle(_) | Curve2::Ellipse(_) => (((last - first).abs() / ARC_STEP).ceil() as usize).max(4),
        _ => 64,
    };
    let params: Vec<f64> = (0..=n).map(|i| first + (last - first) * i as f64 / n as f64).collect();
    pcurve.points_along(&params).into_iter().map(|p| [p.x, p.y]).collect()
}

/// [`sample_pcurve`], with free-form spans halved while their midpoints stray from their
/// chords: for telling points close to the boundary on which side of it they lie.
/// (A free-form pcurve's samples are kept for its curve and range: a face's boundary is
/// sampled again by each operation on it, unchanged where the operation does not reach.)
pub(super) fn sample_pcurve_closely(pcurve: &Curve2, first: f64, last: f64) -> Polyline {
    sample_pcurve_closely_boxed(pcurve, first, last).0
}

/// [`sample_pcurve_closely`], with the samples' box ([`polyline_box`]).
pub(super) fn sample_pcurve_closely_boxed(pcurve: &Curve2, first: f64, last: f64) -> (Polyline, (Uv, Uv)) {
    let boxed = |points: Vec<Uv>| {
        let b = polyline_box(&points);
        (Polyline::from(points), b)
    };
    if matches!(pcurve, Curve2::Line(_) | Curve2::Circle(_) | Curve2::Ellipse(_)) {
        return boxed(sample_pcurve(pcurve, first, last));
    }
    let Some(key) = closely_key(pcurve, first, last) else { return boxed(sample_closely(pcurve, first, last)) };
    let lock = || CLOSE_SAMPLES.lock().unwrap_or_else(|poisoned| poisoned.into_inner());
    let same = |e: &CloseSamples| e.key == key && e.range == [first.to_bits(), last.to_bits()] && e.pcurve == *pcurve;
    if let Some(found) = lock().entries.iter().find(|e| same(e)) {
        return (found.points.clone(), found.bounds);
    }
    // (Sampled outside the lock; another thread may sample it too.)
    let (points, bounds) = boxed(sample_closely(pcurve, first, last));
    let mut cache = lock();
    if !cache.entries.iter().any(|e| same(e)) {
        cache.points += points.len();
        cache.entries.push_back(CloseSamples { key, range: [first.to_bits(), last.to_bits()], pcurve: pcurve.clone(), points: points.clone(), bounds });
        while cache.points > CLOSE_SAMPLES_POINTS || cache.entries.len() > CLOSE_SAMPLES_ENTRIES {
            let Some(old) = cache.entries.pop_front() else { break };
            cache.points -= old.points.len();
        }
    }
    (points, bounds)
}

/// Free-form pcurves' samples kept lately ([`sample_pcurve_closely`]), shared by the threads
/// of operations: oldest first.
struct CloseSampleCache {
    entries: std::collections::VecDeque<CloseSamples>,
    /// The points they hold, in all.
    points: usize,
}

struct CloseSamples {
    /// A hash of the pcurve (its definition, by bits) and range, to tell entries apart
    /// quickly; then compared whole.
    key: u64,
    range: [u64; 2],
    pcurve: Curve2,
    points: Polyline,
    bounds: (Uv, Uv),
}

static CLOSE_SAMPLES: std::sync::Mutex<CloseSampleCache> = std::sync::Mutex::new(CloseSampleCache { entries: std::collections::VecDeque::new(), points: 0 });

/// Points kept at most in all (16 bytes each), and pcurves (each looked through on a lookup).
const CLOSE_SAMPLES_POINTS: usize = 1 << 21;
const CLOSE_SAMPLES_ENTRIES: usize = 4096;

/// [`CloseSamples::key`] for a B-spline or Bézier pcurve (others are not kept).
fn closely_key(pcurve: &Curve2, first: f64, last: f64) -> Option<u64> {
    use std::hash::Hasher;
    let mut h = crate::hashing::FastHasher::default();
    let mut put = |x: f64| h.write_u64(x.to_bits());
    put(first);
    put(last);
    let (poles, weights) = match pcurve {
        Curve2::BSpline(b) => {
            for &k in b.knots() {
                put(k);
            }
            for &m in b.multiplicities() {
                put(m as f64);
            }
            put(b.degree() as f64);
            (b.poles(), b.weights())
        }
        Curve2::Bezier(b) => {
            put(b.degree() as f64);
            (b.poles(), b.weights())
        }
        _ => return None,
    };
    for p in poles {
        put(p.x);
        put(p.y);
    }
    for &w in weights.unwrap_or(&[]) {
        put(w);
    }
    Some(h.finish())
}

fn sample_closely(pcurve: &Curve2, first: f64, last: f64) -> Vec<Uv> {
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
pub(super) fn inside_polylines<P: AsRef<[Uv]>>(polylines: &[P], p: Uv) -> bool {
    let p = ray_turn(p);
    let mut inside = false;
    for poly in polylines {
        for w in poly.as_ref().windows(2) {
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

/// [`inside_polylines`] for many points of the same polylines (given again with each
/// point): the polylines in runs of consecutive segments, each with its box (turned), and,
/// once a few points' rays may have crossed it, its segments (turned) in strips across the
/// ray, each point testing only those of its strip (the same segments that can cross its
/// ray, so the same answers). A face's boundary may be long (many holes, finely sampled)
/// and the points asked about all near one part of it.
pub(super) struct InsidePolylines {
    runs: Vec<Run>,
}

/// Segments of [`InsidePolylines`] at most in a run.
const RUN: usize = 256;

/// A run of consecutive segments of one polyline: those between its points `from` and
/// `to` (included).
struct Run {
    polyline: usize,
    from: usize,
    to: usize,
    /// The range of its points (turned) across the ray, and along it.
    across: (f64, f64),
    along: (f64, f64),
    /// Points whose rays may cross it, so far: its strips are made once a few have been
    /// tested segment by segment.
    asked: std::sync::atomic::AtomicUsize,
    strips: std::sync::OnceLock<Strips>,
}

/// A polyline's segments (by the index of their first point) in strips across the ray
/// (turned): strip `k`'s are `segments[starts[k]..starts[k + 1]]`.
struct Strips {
    low: f64,
    width: f64,
    starts: Vec<u32>,
    segments: Vec<u32>,
}

impl Strips {
    fn new(poly: &[Uv]) -> Self {
        let across = |q: Uv| ray_turn(q)[1];
        let (low, high) = poly.iter().fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), &q| (lo.min(across(q)), hi.max(across(q))));
        let count = (poly.len() / 4).max(1);
        let width = if high > low { (high - low) / count as f64 } else { 1.0 };
        // (A cast rounds down past zero, and up to it from below, as `floor` and `max(0.0)`
        // would, with no call.)
        let strip = |y: f64| (((y - low) / width) as usize).min(count - 1);
        let spans: Vec<(u32, u32)> = poly
            .windows(2)
            .map(|w| {
                let (a, b) = (across(w[0]), across(w[1]));
                (strip(a.min(b)) as u32, strip(a.max(b)) as u32)
            })
            .collect();
        // (Counted first, then placed.)
        let mut starts = vec![0u32; count + 1];
        for &(lo, hi) in &spans {
            for k in lo..=hi {
                starts[k as usize + 1] += 1;
            }
        }
        for k in 0..count {
            starts[k + 1] += starts[k];
        }
        let mut next = starts.clone();
        let mut segments = vec![0u32; starts[count] as usize];
        for (i, &(lo, hi)) in spans.iter().enumerate() {
            for k in lo..=hi {
                segments[next[k as usize] as usize] = i as u32;
                next[k as usize] += 1;
            }
        }
        Strips { low, width, starts, segments }
    }

    /// How many of the segments of `poly` (this was made for) cross the ray from `p`
    /// (turned).
    fn crossings(&self, poly: &[Uv], p: Uv) -> usize {
        if !(p[1] >= self.low) {
            return 0;
        }
        let k = (((p[1] - self.low) / self.width) as usize).min(self.starts.len() - 2);
        let (from, to) = (self.starts[k] as usize, self.starts[k + 1] as usize);
        self.segments[from..to].iter().filter(|&&i| crosses_ray(ray_turn(poly[i as usize]), ray_turn(poly[i as usize + 1]), p)).count()
    }
}

/// Points tested against a polyline segment by segment before its strips are made.
const ASKED_BEFORE_STRIPS: usize = 8;

impl InsidePolylines {
    pub fn new<P: AsRef<[Uv]>>(polylines: &[P]) -> Self {
        let mut runs = Vec::new();
        for (polyline, poly) in polylines.iter().map(AsRef::as_ref).enumerate() {
            let mut from = 0;
            while from + 1 < poly.len() {
                let to = (from + RUN).min(poly.len() - 1);
                let (mut across, mut along) = ((f64::INFINITY, f64::NEG_INFINITY), (f64::INFINITY, f64::NEG_INFINITY));
                for &q in &poly[from..=to] {
                    let t = ray_turn(q);
                    along = (along.0.min(t[0]), along.1.max(t[0]));
                    across = (across.0.min(t[1]), across.1.max(t[1]));
                }
                runs.push(Run { polyline, from, to, across, along, asked: Default::default(), strips: std::sync::OnceLock::new() });
                from = to;
            }
        }
        InsidePolylines { runs }
    }

    /// Whether `p` lies inside `polylines` (those this was made for).
    pub fn contains<P: AsRef<[Uv]>>(&self, polylines: &[P], p: Uv) -> bool {
        let p = ray_turn(p);
        let mut crossings = 0;
        for run in &self.runs {
            // (A segment crosses the ray only between its ends across it, and ahead of `p`
            // only if some of it is: none of a run below, above or behind it, a margin left
            // for the crossing's rounding.)
            let (lo, hi) = run.across;
            if !(p[1] >= lo && p[1] < hi) {
                continue;
            }
            let margin = 1e-9 * (run.along.1.abs() + run.along.0.abs() + 1.0);
            if p[0] > run.along.1 + margin {
                continue;
            }
            let points = &polylines[run.polyline].as_ref()[run.from..=run.to];
            let strips = match run.strips.get() {
                Some(strips) => Some(strips),
                None if run.asked.fetch_add(1, std::sync::atomic::Ordering::Relaxed) >= ASKED_BEFORE_STRIPS => Some(run.strips.get_or_init(|| Strips::new(points))),
                None => None,
            };
            crossings += match strips {
                Some(strips) => strips.crossings(points, p),
                None => points.windows(2).filter(|w| crosses_ray(ray_turn(w[0]), ray_turn(w[1]), p)).count(),
            };
        }
        crossings % 2 == 1
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
pub(super) fn interior_point_among<P: AsRef<[Uv]>>(polylines: &[P], obstacles: &Obstacles) -> Uv {
    // Only obstacles meeting the region's box can cut it.
    let region_box = polylines.iter().fold(([f64::INFINITY; 2], [f64::NEG_INFINITY; 2]), |(lo, hi), p| {
        let (a, b) = polyline_box(p.as_ref());
        ([lo[0].min(a[0]), lo[1].min(a[1])], [hi[0].max(b[0]), hi[1].max(b[1])])
    });
    let obstacles = obstacles.meeting(region_box);
    // The boundary's segments, then the obstacles': one list, seen in its two parts.
    let mut segments: Vec<(Uv, Uv)> = polylines.iter().flat_map(|p| p.as_ref().windows(2).map(|w| (w[0], w[1]))).collect();
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
    // (On a strip far longer than it is wide, as on a cylinder of a small solid, where an
    // angle runs along a tiny length, no point is clear by a part of the diagonal: there,
    // nearly half the width, which only a straight strip's middle is. A curved sliver's box
    // is wider than it is, and it is searched as any other.)
    let (width, height) = (region_box.1[0] - region_box.0[0], region_box.1[1] - region_box.0[1]);
    let narrow = if width.min(height) < 0.02 * width.max(height) { 0.45 * width.min(height) } else { f64::INFINITY };
    let mut best = ([0.0, lo], -1.0);
    let (mut xs, mut cuts_all, mut cuts) = (Vec::new(), Vec::new(), Vec::new());
    let mut try_line = |boundary: &[(Uv, Uv)], blocking: &[(Uv, Uv)], y: f64, best: &mut (Uv, f64)| {
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
        try_line(boundary, blocking, lo + (hi - lo) * (i as f64 + 0.123_456_7) / 32.0, &mut best);
        if best.1 >= (0.05 * size).min(narrow) {
            return best.0;
        }
    }
    if best.1 >= (1e-3 * size).min(narrow) {
        return best.0;
    }
    let mut ys: Vec<f64> = polylines.iter().flat_map(|p| p.as_ref()).map(|p| p[1]).collect();
    ys.extend(obstacles.iter().flat_map(|o| o.iter().map(|p| p[1])));
    ys.sort_by(f64::total_cmp);
    ys.dedup();
    // (Line after line upwards, each against only the segments across its height: those
    // starting below it, kept while they end above. A narrow region of many points has
    // few at any height.)
    let low = |s: &(Uv, Uv)| s.0[1].min(s.1[1]);
    let high = |s: &(Uv, Uv)| s.0[1].max(s.1[1]);
    let by_low = |list: &[(Uv, Uv)]| {
        let mut list = list.to_vec();
        list.sort_by(|p, q| low(p).total_cmp(&low(q)));
        list
    };
    let (boundary, blocking) = (by_low(boundary), by_low(blocking));
    let (mut next, mut across) = ([0usize; 2], [Vec::new(), Vec::new()]);
    for w in ys.windows(2) {
        let y = 0.5 * (w[0] + w[1]);
        for (k, list) in [&boundary, &blocking].into_iter().enumerate() {
            while next[k] < list.len() && low(&list[next[k]]) <= y {
                across[k].push(list[next[k]]);
                next[k] += 1;
            }
            across[k].retain(|s| high(s) > y);
        }
        try_line(&across[0], &across[1], y, &mut best);
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
    /// A tree over the segments' boxes instead, where a cell would hold many (segments
    /// crowded along a thin strip across the cells).
    tree: Option<Bvh>,
}

impl SegmentGrid {
    fn new(segments: &[(Uv, Uv)]) -> Self {
        if segments.len() < 64 {
            return SegmentGrid { origin: [0.0, 0.0], cell: 0.0, n: [0, 0], starts: Vec::new(), segments: Vec::new(), tree: None };
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
        if starts.iter().max().is_some_and(|&most| most > 64) {
            let boxes: Vec<BoundingBox> = segments.iter().map(|(a, b)| BoundingBox::new(Point::new(a[0].min(b[0]), a[1].min(b[1]), 0.0), Point::new(a[0].max(b[0]), a[1].max(b[1]), 0.0))).collect();
            return SegmentGrid { origin: lo, cell, n, starts: Vec::new(), segments: Vec::new(), tree: Some(Bvh::new(&boxes)) };
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
        SegmentGrid { origin: lo, cell, n, starts, segments: members, tree: None }
    }

    /// The distance from `p` to the nearest segment, searching rings of cells outwards.
    fn distance(&self, segments: &[(Uv, Uv)], p: Uv) -> f64 {
        if let Some(tree) = &self.tree {
            return tree.nearest(Point::new(p[0], p[1], 0.0), |k| segment_distance(segments[k], p)).map_or(f64::INFINITY, |(_, d)| d);
        }
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

    /// The nearest point of a curve over a range by brute force: dense samples, then golden
    /// section steps about the best.
    fn nearest_by_samples(curve: &Curve, lo: f64, hi: f64, p: Point) -> (f64, f64) {
        let n = 20000;
        let step = (hi - lo) / n as f64;
        let mut best = (lo, f64::INFINITY);
        for i in 0..=n {
            let t = lo + step * i as f64;
            let d = curve.point(t).distance(p);
            if d < best.1 {
                best = (t, d);
            }
        }
        let (mut x, mut y) = ((best.0 - step).max(lo), (best.0 + step).min(hi));
        let ratio = 0.618_033_988_749_895;
        for _ in 0..100 {
            let (x1, x2) = (y - ratio * (y - x), x + ratio * (y - x));
            if curve.point(x1).distance(p) < curve.point(x2).distance(p) {
                y = x2;
            } else {
                x = x1;
            }
        }
        let t = 0.5 * (x + y);
        let d = curve.point(t).distance(p);
        if d < best.1 { (t, d) } else { best }
    }

    /// Circles and ellipses (round, slender, flat to a sliver) projected onto in closed form
    /// find the nearest point found by dense sampling, over whole turns and arcs (across the
    /// angle's seam too), for points about them, on their axes (and the centre), on them,
    /// inside the ellipse's evolute and on its major and minor axes there, and off the plane.
    #[test]
    fn conics_project_in_closed_form_as_dense_sampling_does() {
        use crate::geometry::curve::{Circle, Ellipse};
        use crate::math::{Direction, Frame};
        let mut seed = 977u64;
        let mut next = || {
            seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            (seed >> 11) as f64 / (1u64 << 53) as f64
        };
        let frame = Frame::new(Point::new(1.0, -2.0, 0.5), Direction::new(0.3, -0.4, 0.85), Direction::new(1.0, 0.2, 0.0));
        let (x, y, z) = (Vector::from(frame.x()), Vector::from(frame.y()), Vector::from(frame.z()));
        let curves = [
            Curve::Circle(Circle::new(frame, 2.5)),
            Curve::Ellipse(Ellipse::new(frame, 3.0, 2.0)),
            Curve::Ellipse(Ellipse::new(frame, 5.0, 0.7)),
            Curve::Ellipse(Ellipse::new(frame, 4.0, 0.01)),
            Curve::Ellipse(Ellipse::new(frame, 1.0, 1.0 - 1e-9)),
        ];
        let tau = std::f64::consts::TAU;
        let ranges = [(0.0, tau), (-1.0, tau - 1.0), (0.3, 2.0), (5.5, 7.0), (-3.0, -2.5), (2.0, 4.5)];
        for curve in &curves {
            let (a, b) = match curve {
                Curve::Circle(c) => (c.radius(), c.radius()),
                Curve::Ellipse(e) => (e.major_radius(), e.minor_radius()),
                _ => unreachable!(),
            };
            let at = |u: f64, v: f64, h: f64| frame.origin() + x * u + y * v + z * h;
            let mut points = vec![at(0.0, 0.0, 0.0), at(0.0, 0.0, 1.5), at(a, 0.0, 0.0), at(-a, 0.0, 0.3), at(0.0, b, 0.0), at(0.0, -b, 0.0)];
            // Inside the evolute, on and beside the axes.
            let d = (a * a - b * b) / a;
            for f in [0.1, 0.5, 0.9, 0.999] {
                points.extend([at(f * d, 0.0, 0.0), at(-f * d, 0.0, 0.2), at(f * d, 1e-12, 0.0), at(f * d, -1e-7, 0.0)]);
                let e = ((a * a - b * b) / b).min(2.0 * a);
                points.extend([at(0.0, f * e, 0.0), at(1e-12, f * e, 0.0), at(0.3 * f * d, 0.2 * f * e, 0.0)]);
            }
            for _ in 0..300 {
                let (u, v, h) = ((2.0 * next() - 1.0) * 1.6 * a, (2.0 * next() - 1.0) * 1.6 * a, (2.0 * next() - 1.0) * a);
                points.push(at(u, v, if next() < 0.5 { 0.0 } else { h }));
            }
            // On the curve.
            for k in 0..40 {
                points.push(curve.point(k as f64 * 0.1618 * tau));
            }
            for &(lo, hi) in &ranges {
                for &p in &points {
                    let (t, dist) = conic_foot(curve, lo, hi, p).unwrap();
                    assert!(t >= lo && t <= hi, "{curve:?} {p:?}: {t} outside {lo}..{hi}");
                    // (A distance to a point of the curve within the range: none nearer.)
                    assert!((curve.point(t).distance(p) - dist).abs() <= 1e-12 * (1.0 + dist));
                    let (_, brute) = nearest_by_samples(curve, lo, hi, p);
                    assert!(dist <= brute + 1e-12 * (1.0 + a), "{curve:?} {p:?} over {lo}..{hi}: {dist} farther than sampled {brute}");
                }
            }
        }
    }

    /// The feet of an ellipse's critical points: the distance's derivative vanishes at each,
    /// and the nearest is listed first.
    #[test]
    fn ellipse_feet_are_critical_points() {
        let (a, b) = (3.0, 1.0);
        let mut out = [0.0; 6];
        for &(x, y) in &[(0.5, 0.2), (2.0, 0.1), (0.0, 0.3), (0.7, 0.0), (5.0, 4.0), (-0.3, -0.1), (2.5, 0.0), (0.0, 0.0), (1e-9, 1e-9)] {
            let (nearest, all) = ellipse_feet(a, b, x, y, true, &mut out);
            assert!(nearest >= 1 && all >= nearest);
            let dist = |t: f64| (a * t.cos() - x).hypot(b * t.sin() - y);
            for &t in &out[..all] {
                let g = (b * b - a * a) * t.sin() * t.cos() + a * x * t.sin() - b * y * t.cos();
                assert!(g.abs() <= 1e-9, "({x}, {y}): g({t}) = {g}");
            }
            let least = (0..100000).map(|i| dist(std::f64::consts::TAU * i as f64 / 100000.0)).fold(f64::INFINITY, f64::min);
            assert!(dist(out[0]) <= least + 1e-12, "({x}, {y}): {} > {least}", dist(out[0]));
        }
    }
    use crate::geometry::curve2::BezierCurve2;
    use crate::math::Point2;

    /// The index answers as testing every segment does, for points asked many times over
    /// (its strips made) and few, near and on long polylines and short ones.
    #[test]
    fn indexed_polylines_answer_as_every_segment_does() {
        let circle = |c: Uv, r: f64, n: usize, wobble: f64| -> Vec<Uv> {
            (0..=n).map(|i| {
                let t = std::f64::consts::TAU * (i % n) as f64 / n as f64;
                let r = r * (1.0 + wobble * (7.0 * t).sin());
                [c[0] + r * t.cos(), c[1] + r * t.sin()]
            }).collect()
        };
        let mut polylines = vec![circle([0.0, 0.0], 10.0, 3000, 0.1)];
        for k in 0..12 {
            let a = k as f64 * 0.5;
            polylines.push(circle([5.0 * a.cos(), 5.0 * a.sin()], 0.8, 700 + 37 * k, 0.05));
        }
        // An open polyline and a single segment, as stretches of boundary are.
        polylines.push(vec![[-3.0, -9.0], [3.0, -8.5], [4.0, -7.0]]);
        polylines.push(vec![[-3.0, -9.0], [4.0, -7.0]]);
        let index = InsidePolylines::new(&polylines);
        let mut seed = 12345u64;
        let mut next = || {
            seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            (seed >> 11) as f64 / (1u64 << 53) as f64
        };
        for round in 0..4 {
            for _ in 0..2000 {
                let p = [-12.0 + 24.0 * next(), -12.0 + 24.0 * next()];
                assert_eq!(index.contains(&polylines, p), inside_polylines(&polylines, p), "{p:?} in round {round}");
            }
            // The polylines' own points, on the boundary.
            for poly in &polylines {
                for &q in poly.iter().step_by(17) {
                    assert_eq!(index.contains(&polylines, q), inside_polylines(&polylines, q), "{q:?}");
                }
            }
        }
    }

    /// A band winding a turn and a half round an axis, as a thread's flank (its turns 8 apart,
    /// 88 round): every point on it found on it, on its own turn, not held on the next turn
    /// nearer than the band's samples on its own.
    #[test]
    fn points_on_a_band_winding_round_are_found_on_their_own_turn() {
        use crate::geometry::fit;
        use crate::geometry::surface::{BSplineSurface, KnotVector};
        let (turns, pitch) = (1.5, 8.0);
        let n = 96;
        let params: Vec<f64> = (0..=n).map(|i| i as f64 / n as f64).collect();
        let helix = |r: f64, dz: f64| {
            let points: Vec<Point> = params
                .iter()
                .map(|&t| {
                    let a = std::f64::consts::TAU * turns * t;
                    Point::new(r * a.cos(), r * a.sin(), pitch * turns * t + dz)
                })
                .collect();
            fit::interpolate_with_parameters(&points, &params).unwrap()
        };
        let (inner, outer) = (helix(84.0, 0.0), helix(88.0, 2.5));
        let u = KnotVector::new(inner.knots().to_vec(), inner.multiplicities().to_vec(), inner.degree());
        let v = KnotVector::new(vec![0.0, 1.0], vec![2, 2], 1);
        let poles: Vec<Vec<Point>> = inner.poles().iter().zip(outer.poles()).map(|(&a, &b)| vec![a, b]).collect();
        let band = Surface::BSpline(BSplineSurface::new(poles, u, v));
        for i in 0..=40 {
            for j in [0.2, 0.5, 0.8] {
                let p = surface_point(&band, [i as f64 / 40.0, j]);
                let uv = uv_within(&band, p, 1e-6).unwrap_or_else(|| panic!("{p:?} not found on the band"));
                assert!(surface_point(&band, uv).distance(p) <= 1e-6);
            }
        }
    }

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
        assert!(inside_polylines(&outline(sample_pcurve_closely(&arch, 0.0, 1.0).to_vec()), q));
    }
}

