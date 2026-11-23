//! Where a steady motion runs along a face: its characteristic curves, `n · V = 0`, found in
//! the face's parameters.
//!
//! On each face `g(u, v) = n(u, v) · V(S(u, v))` (the unit normal, the velocity of the point
//! the surface has there) is sampled on a grid over the face's parameter box (a little wider,
//! so that each curve crosses the face's boundary rather than ending on it). The grid is as
//! fine as the normal turns: no cell spans more than a sixty-fourth of a half turn of it, so
//! a curve's branches, which part where the normal turns, fall in different cells. The zero
//! set is followed through the cells (marching squares, a saddle cell decided by its centre),
//! each place it crosses a grid line found exactly (by regula falsi on `g` along the line).
//! Between those points the curve is fitted as cubics in the parameters and in space through
//! the same parameters, and points are added (each projected onto `g = 0` by Newton steps
//! along the gradient) until the fits stay within the given tolerance of the exact curve
//! everywhere between them.

use crate::divide::Cut;
use crate::error::{Error, Result};
use crate::geometry::curve::Curve;
use crate::geometry::curve2::Curve2;
use crate::geometry::surface::Surface;
use crate::geometry::{ParametricCurve, ParametricSurface};
use crate::math::{Point, Point2, Vector};
use crate::topology::{Model, Shape};

use super::Steady;

/// `g` on one surface.
pub(super) struct Gauge<'a> {
    pub surface: &'a Surface,
    pub steady: &'a Steady,
}

impl Gauge<'_> {
    /// The unit normal (as the surface's parameters turn it) and the point at `(u, v)`.
    pub fn normal(&self, u: f64, v: f64) -> (Point, Vector) {
        let (p, su, sv) = self.surface.d1(u, v);
        let n = su.cross(sv);
        let length = n.length();
        if length > 1e-12 * su.length() * sv.length() && length > 0.0 {
            return (p, n * (1.0 / length));
        }
        let n = crate::geometry::props::surface_normal(self.surface, u, v, 0.0).map(|n| Vector::from(n.direction())).unwrap_or(Vector::ZERO);
        (p, n)
    }

    pub fn g(&self, u: f64, v: f64) -> f64 {
        let (p, n) = self.normal(u, v);
        n.dot(self.steady.velocity(p))
    }
}

/// The parameter box `(u0, u1, v0, v1)` a face's curves are followed over: its own, widened
/// by a fiftieth each way, kept within the surface's parameters where they end.
fn search_box(surface: &Surface, face: (f64, f64, f64, f64)) -> ((f64, f64, f64, f64), [bool; 4]) {
    let (u0, u1, v0, v1) = face;
    let (su0, su1, sv0, sv1) = surface.bounds();
    let (mu, mv) = (0.02 * (u1 - u0).max(1e-9), 0.02 * (v1 - v0).max(1e-9));
    let clamp = |x: f64, lo: f64, hi: f64, periodic: bool| if periodic { x } else { x.clamp(lo, hi) };
    let (pu, pv) = (surface.u_period().is_some(), surface.v_period().is_some());
    let b = (clamp(u0 - mu, su0, su1, pu), clamp(u1 + mu, su0, su1, pu), clamp(v0 - mv, sv0, sv1, pv), clamp(v1 + mv, sv0, sv1, pv));
    // (Which sides could not be widened: the surface ends there.)
    let held = [b.0 > u0 - 0.5 * mu, b.1 < u1 + 0.5 * mu, b.2 > v0 - 0.5 * mv, b.3 < v1 + 0.5 * mv];
    (b, held)
}

/// How many cells along u and along v: enough that the normal turns by at most a
/// sixty-fourth of a half turn across a cell, at least four (a plane's `g` is linear).
fn grid_size(gauge: &Gauge, b: (f64, f64, f64, f64)) -> (usize, usize) {
    const PROBES: usize = 16;
    let (u0, u1, v0, v1) = b;
    let normals: Vec<Vec<Vector>> = (0..=PROBES)
        .map(|i| (0..=PROBES).map(|j| gauge.normal(u0 + (u1 - u0) * i as f64 / PROBES as f64, v0 + (v1 - v0) * j as f64 / PROBES as f64).1).collect())
        .collect();
    let angle = |a: Vector, b: Vector| a.dot(b).clamp(-1.0, 1.0).acos();
    let (mut turn_u, mut turn_v) = (0.0f64, 0.0f64);
    for k in 0..=PROBES {
        turn_u = turn_u.max((0..PROBES).map(|i| angle(normals[i][k], normals[i + 1][k])).sum());
        turn_v = turn_v.max((0..PROBES).map(|j| angle(normals[k][j], normals[k][j + 1])).sum());
    }
    let cells = |turn: f64| ((turn / (std::f64::consts::PI / 64.0)).ceil() as usize).clamp(4, 512);
    (cells(turn_u), cells(turn_v))
}

/// The root of `f` on `[0, 1]` where it changes sign between `fa = f(0)` and `fb = f(1)`
/// (regula falsi, Illinois).
fn root(f: &dyn Fn(f64) -> f64, mut fa: f64, mut fb: f64) -> f64 {
    let (mut a, mut b) = (0.0f64, 1.0f64);
    let mut side = 0;
    for _ in 0..100 {
        let x = if fb != fa { (a * fb - b * fa) / (fb - fa) } else { 0.5 * (a + b) };
        let x = if x > a && x < b { x } else { 0.5 * (a + b) };
        let fx = f(x);
        if fx == 0.0 || (b - a) <= 1e-15 {
            return x;
        }
        if (fx > 0.0) == (fa > 0.0) {
            a = x;
            fa = fx;
            if side == -1 {
                fb *= 0.5;
            }
            side = -1;
        } else {
            b = x;
            fb = fx;
            if side == 1 {
                fa *= 0.5;
            }
            side = 1;
        }
    }
    0.5 * (a + b)
}

/// One curve of the zero set as the grid follows it: its points in the parameters, in order,
/// and whether it closes.
struct Chain {
    points: Vec<Point2>,
    closed: bool,
}

/// The zero set of `g` over the box, as chains of the points where it crosses the grid lines.
fn chains(gauge: &Gauge, b: (f64, f64, f64, f64), nu: usize, nv: usize) -> Vec<Chain> {
    let (u0, u1, v0, v1) = b;
    let us: Vec<f64> = (0..=nu).map(|i| u0 + (u1 - u0) * i as f64 / nu as f64).collect();
    let vs: Vec<f64> = (0..=nv).map(|j| v0 + (v1 - v0) * j as f64 / nv as f64).collect();
    let values: Vec<Vec<f64>> = us.iter().map(|&u| vs.iter().map(|&v| gauge.g(u, v)).collect()).collect();
    let positive = |i: usize, j: usize| values[i][j] >= 0.0;
    // Crossing points, by grid line: along u (i, j)–(i+1, j), and along v (i, j)–(i, j+1).
    let along_u = |i: usize, j: usize| j * nu + i;
    let along_v = |i: usize, j: usize| nu * (nv + 1) + i * nv + j;
    let count = nu * (nv + 1) + (nu + 1) * nv;
    let mut points: Vec<Option<Point2>> = vec![None; count];
    for j in 0..=nv {
        for i in 0..nu {
            if positive(i, j) != positive(i + 1, j) {
                let (ua, ub, v) = (us[i], us[i + 1], vs[j]);
                let s = root(&|s| gauge.g(ua + (ub - ua) * s, v), values[i][j], values[i + 1][j]);
                points[along_u(i, j)] = Some(Point2::new(ua + (ub - ua) * s, v));
            }
        }
    }
    for i in 0..=nu {
        for j in 0..nv {
            if positive(i, j) != positive(i, j + 1) {
                let (va, vb, u) = (vs[j], vs[j + 1], us[i]);
                let s = root(&|s| gauge.g(u, va + (vb - va) * s), values[i][j], values[i][j + 1]);
                points[along_v(i, j)] = Some(Point2::new(u, va + (vb - va) * s));
            }
        }
    }
    // The segments through each cell.
    let mut links: Vec<Vec<usize>> = vec![Vec::new(); count];
    let mut link = |a: usize, b: usize| {
        links[a].push(b);
        links[b].push(a);
    };
    for i in 0..nu {
        for j in 0..nv {
            let (bottom, top, left, right) = (along_u(i, j), along_u(i, j + 1), along_v(i, j), along_v(i + 1, j));
            let crossed: Vec<usize> = [bottom, right, top, left].into_iter().filter(|&k| points[k].is_some()).collect();
            match crossed.len() {
                2 => link(crossed[0], crossed[1]),
                4 => {
                    let centre = gauge.g(0.5 * (us[i] + us[i + 1]), 0.5 * (vs[j] + vs[j + 1])) >= 0.0;
                    if centre == positive(i, j) {
                        link(bottom, right);
                        link(top, left);
                    } else {
                        link(bottom, left);
                        link(top, right);
                    }
                }
                _ => {}
            }
        }
    }
    // Walked: from the ends first, then round what is left (closed).
    let mut used = vec![false; count];
    let mut out = Vec::new();
    let mut walk = |start: usize, used: &mut Vec<bool>| -> Chain {
        let mut chain = vec![start];
        used[start] = true;
        let mut at = start;
        loop {
            let Some(&next) = links[at].iter().find(|&&k| !used[k]) else { break };
            used[next] = true;
            chain.push(next);
            at = next;
        }
        let closed = chain.len() > 2 && links[at].contains(&start);
        Chain { points: chain.iter().map(|&k| points[k].expect("a crossing")).collect(), closed }
    };
    for k in 0..count {
        if points[k].is_some() && !used[k] && links[k].len() == 1 {
            out.push(walk(k, &mut used));
        }
    }
    for k in 0..count {
        if points[k].is_some() && !used[k] && !links[k].is_empty() {
            out.push(walk(k, &mut used));
        }
    }
    out.retain(|c| c.points.len() >= 2);
    out
}

/// `uv` moved onto `g = 0` by Newton steps along the gradient; `None` if they do not settle.
fn project(gauge: &Gauge, uv: Point2, h: (f64, f64)) -> Option<Point2> {
    let mut q = uv;
    for _ in 0..30 {
        let g = gauge.g(q.x, q.y);
        let gu = (gauge.g(q.x + h.0, q.y) - gauge.g(q.x - h.0, q.y)) / (2.0 * h.0);
        let gv = (gauge.g(q.x, q.y + h.1) - gauge.g(q.x, q.y - h.1)) / (2.0 * h.1);
        let norm = gu * gu + gv * gv;
        if !(norm > 0.0) {
            return None;
        }
        let step = Point2::new(-g * gu / norm, -g * gv / norm);
        q = Point2::new(q.x + step.x, q.y + step.y);
        if step.x.abs() <= 1e-14 * h.0 * 1e6 && step.y.abs() <= 1e-14 * h.1 * 1e6 {
            return Some(q);
        }
    }
    let g = gauge.g(q.x, q.y);
    (g.abs() <= 1e-12 * gauge.steady.speed_scale()).then_some(q)
}

/// A curve of the zero set fitted: in the parameters and in space, over one parameter.
pub(super) struct Fitted {
    pub pcurve: crate::geometry::curve2::BSplineCurve2,
    pub curve: crate::geometry::curve::BSplineCurve,
    pub first: f64,
    pub last: f64,
}

/// Fits the points (on the zero set, in order; the first not repeated at the end if
/// `closed`) as cubics within `tolerance` in space of the zero set, adding points where
/// they stray.
fn fit(gauge: &Gauge, mut points: Vec<Point2>, closed: bool, extend: (bool, bool), h: (f64, f64), tolerance: f64) -> Result<Fitted> {
    let failed = |what: &str| Error::Failed(format!("following where a motion grazes a face: {what}"));
    // (Points closer than a hair in space are one.)
    let hair = 1e-3 * tolerance;
    for _round in 0..16 {
        let mut kept: Vec<Point2> = Vec::with_capacity(points.len());
        for &p in &points {
            if kept.last().is_none_or(|&q: &Point2| gauge.surface.point(q.x, q.y).distance(gauge.surface.point(p.x, p.y)) > hair) {
                kept.push(p);
            }
        }
        if closed && kept.len() > 2 && gauge.surface.point(kept[0].x, kept[0].y).distance(gauge.surface.point(kept[kept.len() - 1].x, kept[kept.len() - 1].y)) <= hair {
            kept.pop();
        }
        points = kept;
        if points.len() < 2 || (closed && points.len() < 3) {
            return Err(failed("a curve of fewer than two points"));
        }
        let space: Vec<Point> = points.iter().map(|q| gauge.surface.point(q.x, q.y)).collect();
        let mut params = vec![0.0];
        for w in space.windows(2) {
            params.push(params.last().unwrap() + w[0].distance(w[1]));
        }
        let (pcurve, curve, first, last) = if closed {
            let period = params.last().unwrap() + space.last().unwrap().distance(space[0]);
            let pc = crate::geometry::fit::interpolate_periodic_with_parameters(&points, &params, period)?;
            let c = crate::geometry::fit::interpolate_periodic_with_parameters(&space, &params, period)?;
            (pc, c, 0.0, period)
        } else {
            // (Carried on straight past an end the box holds, a fiftieth of its length.)
            let (mut uv, mut xyz, mut at) = (points.clone(), space.clone(), params.clone());
            let reach = 0.02 * params.last().unwrap();
            let n = points.len();
            if extend.0 {
                let (du, dx) = (points[0] - points[1], space[0] - space[1]);
                let step = params[1];
                for k in [1.0, 2.0] {
                    let f = reach * k / 2.0 / step;
                    uv.insert(0, points[0] + du * f);
                    xyz.insert(0, space[0] + dx * f);
                    at.insert(0, -reach * k / 2.0);
                }
            }
            if extend.1 {
                let (du, dx) = (points[n - 1] - points[n - 2], space[n - 1] - space[n - 2]);
                let step = params[n - 1] - params[n - 2];
                for k in [1.0, 2.0] {
                    let f = reach * k / 2.0 / step;
                    uv.push(points[n - 1] + du * f);
                    xyz.push(space[n - 1] + dx * f);
                    at.push(params[n - 1] + reach * k / 2.0);
                }
            }
            let pc = crate::geometry::fit::interpolate_with_parameters(&uv, &at)?;
            let c = crate::geometry::fit::interpolate_with_parameters(&xyz, &at)?;
            (pc, c, at[0], *at.last().unwrap())
        };
        let last_exact = if closed { last } else { *params.last().unwrap() };
        // Midway between each two points: the fits against the zero set there.
        let n = points.len();
        let spans = if closed { n } else { n - 1 };
        let mut insert: Vec<(usize, Point2)> = Vec::new();
        for i in 0..spans {
            let (a, b) = (params[i], if i + 1 < n { params[i + 1] } else { last_exact });
            for f in [0.25, 0.5, 0.75] {
                let s = a + (b - a) * f;
                let uv = pcurve.point(s);
                let Some(q) = project(gauge, uv, h) else {
                    return Err(failed("a point that does not settle on the curve"));
                };
                let exact = gauge.surface.point(q.x, q.y);
                let error = exact.distance(curve.point(s)).max(exact.distance(gauge.surface.point(uv.x, uv.y)));
                if error > tolerance {
                    if f == 0.5 || !insert.iter().any(|&(k, _)| k == i) {
                        let q = project(gauge, pcurve.point(0.5 * (a + b)), h).unwrap_or(q);
                        insert.push((i, q));
                    }
                    break;
                }
            }
        }
        if insert.is_empty() {
            return Ok(Fitted { pcurve, curve, first, last });
        }
        insert.dedup_by_key(|x| x.0);
        let mut next = Vec::with_capacity(n + insert.len());
        let mut k = 0;
        for (i, &p) in points.iter().enumerate() {
            next.push(p);
            while k < insert.len() && insert[k].0 == i {
                next.push(insert[k].1);
                k += 1;
            }
        }
        points = next;
        if points.len() > 20000 {
            break;
        }
    }
    Err(failed("a curve the fits do not close on"))
}

/// Whether `q` lies within the closed polylines (even-odd).
pub(super) fn within(polygons: &[Vec<Point2>], q: Point2) -> bool {
    let mut inside = false;
    for poly in polygons {
        for i in 0..poly.len() {
            let (a, b) = (poly[i], poly[(i + 1) % poly.len()]);
            if (a.y > q.y) != (b.y > q.y) && q.x < a.x + (q.y - a.y) / (b.y - a.y) * (b.x - a.x) {
                inside = !inside;
            }
        }
    }
    inside
}

/// The face's loops as polylines in its parameters.
pub(super) fn loop_polygons(model: &Model, face: Shape) -> Option<Vec<Vec<Point2>>> {
    let loops = crate::heal::face_loops(model, face)?;
    Some(loops.iter().map(|l| l.iter().flat_map(|u| (0..32).map(move |i| u.at(i as f64 / 32.0))).collect()).collect())
}

/// Distance from `q` to the polylines, in the parameters.
fn distance_to(polygons: &[Vec<Point2>], q: Point2) -> f64 {
    let mut best = f64::INFINITY;
    for poly in polygons {
        for i in 0..poly.len() {
            let (a, b) = (poly[i], poly[(i + 1) % poly.len()]);
            let d = b - a;
            let l = d.x * d.x + d.y * d.y;
            let s = if l > 0.0 { (((q.x - a.x) * d.x + (q.y - a.y) * d.y) / l).clamp(0.0, 1.0) } else { 0.0 };
            best = best.min(q.distance(Point2::new(a.x + d.x * s, a.y + d.y * s)));
        }
    }
    best
}

/// The cuts dividing `face` along its characteristic curves under `steady`, each fitted within
/// `tolerance` in space; none where the motion runs along the face everywhere (it slides
/// within its own surface) or nowhere.
pub(super) fn cuts(model: &Model, face: Shape, steady: &Steady, tolerance: f64) -> Result<Vec<Cut>> {
    let Some(surface) = model.surface(face) else { return Ok(Vec::new()) };
    let gauge = Gauge { surface, steady };
    let Some(polygons) = loop_polygons(model, face) else {
        return Err(Error::Failed("a face without pcurves".into()));
    };
    let bounds = model.face_uv_bounds(face);
    let (b, held) = search_box(surface, bounds);
    let (nu, nv) = grid_size(&gauge, b);
    let h = (1e-7 * (b.1 - b.0), 1e-7 * (b.3 - b.2));
    let scale = (bounds.1 - bounds.0).hypot(bounds.3 - bounds.2);
    let mut out = Vec::new();
    for chain in chains(&gauge, b, nu, nv) {
        // A curve along the face's own boundary, or off it, divides nothing.
        if !chain.points.iter().any(|&q| within(&polygons, q) && distance_to(&polygons, q) > 1e-7 * scale) {
            continue;
        }
        let inside = chain.points.iter().filter(|&&q| within(&polygons, q)).count();
        let mut points = chain.points;
        let mut closed = chain.closed;
        if closed && inside < points.len() {
            // Closed but crossing out of the face: opened at a point outside it, ending there.
            let k = points.iter().position(|&q| !within(&polygons, q) && distance_to(&polygons, q) > 1e-6 * scale).unwrap_or(0);
            points.rotate_left(k);
            points.push(points[0]);
            closed = false;
        }
        let extend = if closed { (false, false) } else { (on_held(b, held, points[0]), on_held(b, held, *points.last().unwrap())) };
        let fitted = fit(&gauge, points, closed, extend, h, tolerance)?;
        out.push(Cut { pcurve: Curve2::BSpline(fitted.pcurve), curve: Curve::BSpline(fitted.curve), first: fitted.first, last: fitted.last });
    }
    Ok(out)
}

/// Whether `q` lies on a side of the box the surface ends at (where a curve could not be
/// followed past the face, so is carried on straight).
fn on_held(b: (f64, f64, f64, f64), held: [bool; 4], q: Point2) -> bool {
    let (u0, u1, v0, v1) = b;
    let eps = 1e-9 * ((u1 - u0) + (v1 - v0));
    (held[0] && (q.x - u0).abs() <= eps) || (held[1] && (q.x - u1).abs() <= eps) || (held[2] && (q.y - v0).abs() <= eps) || (held[3] && (q.y - v1).abs() <= eps)
}
