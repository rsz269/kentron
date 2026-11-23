//! Sweeps under a general rigid motion, given by where it puts the body at each instant.
//!
//! Seen from the body, its points move at `V(t, p) = ω(t) × p + v(t)`: the field changes with
//! `t`, and so may where it runs along the faces. Over a stretch of the motion where it does
//! not (each face's `V·n` keeps its sign everywhere, or vanishes along curves that stay where
//! they are, its values only scaled), `W` and its boundary are fixed on the body as under a
//! steady motion, and the swept solid of the stretch is built the same way, the boundary's
//! curves carried by the motion (B-spline surfaces through their poles' paths). The motion is
//! cut where a face turns over (where its `V·n`, scaled, passes through zero), and the whole
//! sweep is the union of its stretches'. A motion along which a curve where it grazes a face
//! moves over the face is refused: its envelope is not a product.

use std::sync::Arc;

use crate::error::{Error, Result};
use crate::geometry::ParametricSurface;
use crate::math::{Point, Point2, Transform, Transformable, Vector};
use crate::topology::{Model, NodeKind, Shape};

use super::{contour, Carriage, Steady};

/// One stretch of a general motion, over which `W` stays where it is on the body.
pub(in crate::sweep) struct Stage {
    place: Arc<dyn Fn(f64) -> Transform + Send + Sync>,
    from: f64,
    to: f64,
    /// Back from where the motion puts the body at `from`.
    start_inverse: Transform,
    /// The body's twist midway, as the body stands at the start of the whole motion, against
    /// the stretch's fraction.
    omega: Vector,
    linear: Vector,
}

impl Stage {
    /// Where the stretch puts the body as it stands at its start, at fraction `f` of it.
    pub(super) fn at(&self, f: f64) -> Transform {
        (self.place)(self.from + (self.to - self.from) * f) * self.start_inverse
    }

    /// The velocity midway of the body's point `p` (as the body stands at the start of the
    /// whole motion): the
    /// characteristic curves are where it runs along the faces, all over the stretch.
    pub(super) fn velocity(&self, p: Point) -> Vector {
        self.omega.cross(p - Point::ORIGIN) + self.linear
    }
}

impl std::fmt::Debug for Stage {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Stage({}..{})", self.from, self.to)
    }
}

/// The body's twist at `t` as it stands where `frame` puts it (`frame_inverse` taking it
/// back): `(ω, v)`, its points `p` moving at `ω × p + v` (per unit of `t`), from central
/// differences of the places.
fn twist(place: &(dyn Fn(f64) -> Transform + Send + Sync), t: f64, frame_inverse: &Transform) -> (Vector, Vector) {
    let h = 1e-5;
    let n = |s: f64| place(s) * *frame_inverse;
    let back = n(t).inverse();
    // (Second-order differences: central inside, one-sided at the ends, where the places
    // stop.)
    let (offsets, weights): ([f64; 3], [f64; 3]) = if t - h < 0.0 {
        ([0.0, h, 2.0 * h], [-1.5, 2.0, -0.5])
    } else if t + h > 1.0 {
        ([0.0, -h, -2.0 * h], [1.5, -2.0, 0.5])
    } else {
        ([-h, 0.0, h], [-0.5, 0.0, 0.5])
    };
    let places: Vec<Transform> = offsets.iter().map(|&d| n(t + d)).collect();
    // Each probe's velocity, in the body's own axes.
    let w = |p: Point| {
        let moved: Vec<Point> = places.iter().map(|m| p.transformed(m)).collect();
        let v = (moved[0] - Point::ORIGIN) * weights[0] + (moved[1] - Point::ORIGIN) * weights[1] + (moved[2] - Point::ORIGIN) * weights[2];
        v.transformed(&back) * (1.0 / h)
    };
    let origin = w(Point::ORIGIN);
    let (ex, ey, ez) = (w(Point::new(1.0, 0.0, 0.0)) - origin, w(Point::new(0.0, 1.0, 0.0)) - origin, w(Point::new(0.0, 0.0, 1.0)) - origin);
    // (The columns of the skew matrix ω×: ω = (A₃₂, A₁₃, A₂₁), averaged with the transposes.)
    let omega = Vector::new(0.5 * (ey.z - ez.y), 0.5 * (ez.x - ex.z), 0.5 * (ex.y - ey.x));
    (omega, origin)
}

/// Points inside a face with its outward normals there (at most a grid's worth).
fn face_samples(model: &Model, face: Shape) -> Result<Vec<(Point, Vector)>> {
    let surface = model.surface(face).ok_or_else(|| Error::InvalidInput("a face without a surface".into()))?;
    let polygons = contour::loop_polygons(model, face).ok_or_else(|| Error::Failed("a face without pcurves".into()))?;
    let sign = if face.orientation == crate::topology::Orientation::Reversed { -1.0 } else { 1.0 };
    let (u0, u1, v0, v1) = model.face_uv_bounds(face);
    const N: usize = 16;
    let mut out = Vec::new();
    for i in 0..N {
        for j in 0..N {
            let q = Point2::new(u0 + (u1 - u0) * (i as f64 + 0.5) / N as f64, v0 + (v1 - v0) * (j as f64 + 0.5) / N as f64);
            if !contour::within(&polygons, q) {
                continue;
            }
            let (p, su, sv) = surface.d1(q.x, q.y);
            let n = su.cross(sv);
            if n.length() > 0.0 {
                out.push((p, n * (sign / n.length())));
            }
        }
    }
    Ok(out)
}

/// The swept solid of `tool` carried by `carriage`.
pub(super) fn sweep(model: &mut Model, tool: Shape, carriage: &Carriage) -> Result<Shape> {
    let place = carriage.0.clone();
    let start = place(0.0);
    let placed = crate::modify::transform(model, tool, &start, Default::default())?;
    let start_inverse = start.inverse();
    let bounds = crate::properties::bounding_box(model, placed);
    let corners: Vec<Point> = (0..8).map(|i| Point::new(if i & 1 == 0 { bounds.min.x } else { bounds.max.x }, if i & 2 == 0 { bounds.min.y } else { bounds.max.y }, if i & 4 == 0 { bounds.min.z } else { bounds.max.z })).collect();
    let faces = model.unique_subshapes(placed, NodeKind::Face);
    let samples: Vec<Vec<(Point, Vector)>> = faces.iter().map(|&f| face_samples(model, f)).collect::<Result<_>>()?;
    // The face's V·n at each sample, at `t` (with the speeds' scale then).
    let field = |t: f64| twist(&*place, t, &start_inverse);
    let values = |t: f64, k: usize| -> (Vec<f64>, f64) {
        let (omega, linear) = field(t);
        let speed = corners.iter().map(|&p| (omega.cross(p - Point::ORIGIN) + linear).length()).fold(0.0, f64::max);
        (samples[k].iter().map(|&(p, n)| n.dot(omega.cross(p - Point::ORIGIN) + linear)).collect(), speed)
    };
    // Each face's state at each sampled instant: still (0), carried one way or the other
    // (±1), the curves where it is grazed fixed.
    const TIMES: usize = 128;
    let times: Vec<f64> = (0..=TIMES).map(|j| j as f64 / TIMES as f64).collect();
    let mut cuts: Vec<f64> = Vec::new();
    let mut still_faces: Vec<usize> = Vec::new();
    for k in 0..faces.len() {
        if samples[k].is_empty() {
            continue;
        }
        let rows: Vec<(Vec<f64>, f64)> = times.iter().map(|&t| values(t, k)).collect();
        let still = |g: &[f64], speed: f64| g.iter().all(|x| x.abs() <= 1e-9 * speed);
        let mixed = |g: &[f64], speed: f64| g.iter().any(|&x| x > 1e-9 * speed) && g.iter().any(|&x| x < -1e-9 * speed);
        // (Where it is grazed somewhere, the instant it is carried most: the reference its
        // values must be a multiple of all along.)
        let reference = rows.iter().filter(|(g, s)| mixed(g, *s)).max_by(|a, b| norm(&a.0).total_cmp(&norm(&b.0))).map(|(g, _)| g.clone());
        let mut state = Vec::with_capacity(rows.len());
        for (g, speed) in &rows {
            if still(g, *speed) {
                state.push(0);
                continue;
            }
            match &reference {
                Some(r) => {
                    let lambda = dot(g, r) / dot(r, r);
                    let off = g.iter().zip(r).map(|(x, y)| (x - lambda * y).powi(2)).sum::<f64>().sqrt();
                    // (A multiple within the precision of the differences the twist comes from,
                    // against the speeds.)
                    if off > 1e-6 * speed * (g.len() as f64).sqrt() {
                        return Err(Error::InvalidInput("a motion under which where it grazes the tool moves over it: its envelope is no product of a curve and the motion".into()));
                    }
                    state.push(if lambda > 0.0 { 1 } else { -1 });
                }
                None => state.push(if g.iter().any(|&x| x > 0.0) { 1 } else { -1 }),
            }
        }
        if state.iter().all(|&x| x == 0) {
            still_faces.push(k);
        }
        // Where it turns over: the instant its value at its most carried sample passes zero.
        let probe = {
            let r = reference.clone().unwrap_or_else(|| rows.iter().max_by(|a, b| norm(&a.0).total_cmp(&norm(&b.0))).map(|(g, _)| g.clone()).unwrap_or_default());
            (0..r.len()).max_by(|&a, &b| r[a].abs().total_cmp(&r[b].abs())).unwrap_or(0)
        };
        for j in 1..times.len() {
            if state[j] != state[j - 1] && state[j] != 0 && state[j - 1] != 0 {
                let f = |x: f64| values(times[j - 1] + (times[j] - times[j - 1]) * x, k).0[probe];
                let x = contour::root(&f, f(0.0), f(1.0));
                cuts.push(times[j - 1] + (times[j] - times[j - 1]) * x);
            } else if state[j] != state[j - 1] && (state[j] == 0) != (state[j - 1] == 0) {
                // (Still only at an instant between: there it turns over too, or stops.)
                let at = if state[j] == 0 { times[j] } else { times[j - 1] };
                cuts.push(at);
            }
        }
    }
    // A face the motion runs along throughout, beside one it carries, while turning: the
    // face's edges sweep regions of its surface overlapping one another as the motion turns,
    // which the product does not resolve.
    let turning = times.iter().any(|&t| field(t).0.length() > 1e-9 * field(t).1.length().max(1.0));
    if turning {
        for &k in &still_faces {
            let edges = model.unique_subshapes(faces[k], NodeKind::Edge);
            let beside_carried = (0..faces.len()).filter(|j| !still_faces.contains(j)).any(|j| model.unique_subshapes(faces[j], NodeKind::Edge).iter().any(|e| edges.iter().any(|f| f.node == e.node)));
            if beside_carried {
                return Err(Error::InvalidInput("a turning motion running along a face of the tool beside faces it carries: the face's edges sweep regions of its surface over one another".into()));
            }
        }
    }
    cuts.sort_by(f64::total_cmp);
    cuts.dedup_by(|a, b| (*a - *b).abs() <= 1e-9);
    let mut bounds_t = vec![0.0];
    bounds_t.extend(cuts.into_iter().filter(|&t| t > 1e-9 && t < 1.0 - 1e-9));
    bounds_t.push(1.0);
    // The stages, each seen from the body where it stands at the start, carrying it from
    // where the stage starts.
    let mut stages = Vec::new();
    let mut junctions = vec![Transform::IDENTITY];
    for w in bounds_t.windows(2) {
        let (from, to) = (w[0], w[1]);
        let (omega, linear) = stage_twist(&*place, from, to, &start_inverse);
        stages.push(Steady::Carried(Arc::new(Stage { place: place.clone(), from, to, start_inverse: place(from).inverse(), omega, linear })));
        junctions.push(place(to) * start_inverse);
    }
    super::product::sweep_in_stages(model, placed, &stages, junctions)
}

/// The twist midway along the stretch, on the body as it stands at the start of the motion:
/// its points `p` move (per unit of the stretch's fraction) at `ω × p + v` midway, in the body's
/// own axes then.
fn stage_twist(place: &(dyn Fn(f64) -> Transform + Send + Sync), from: f64, to: f64, start_inverse: &Transform) -> (Vector, Vector) {
    let (omega, linear) = twist(place, 0.5 * (from + to), start_inverse);
    (omega * (to - from), linear * (to - from))
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

fn norm(a: &[f64]) -> f64 {
    dot(a, a).sqrt()
}
