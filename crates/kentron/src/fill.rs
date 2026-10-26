//! Filling surfaces: a face spanning a closed loop of edges, through given points within it.
//!
//! The surface is the plate of least bending energy over the loop's mean plane: its height
//! above the plane interpolates the loop (sampled densely along each edge) and the points,
//! a polyharmonic spline (for an energy of order `m`, the kernel `r^(2m-2) ln r` with
//! polynomials of degree `m - 1`; `m = 2` the thin plate). It is fitted by a bicubic
//! B-spline whose parameters are the plane's own coordinates, so that each edge's curve on
//! it is the edge projected onto the plane. The loop must project onto the plane as a simple
//! loop (the face a graph over it).

use crate::error::{Error, Result};
use crate::geometry::curve2::Curve2;
use crate::geometry::fit;
use crate::geometry::surface::{BSplineSurface, KnotVector, Surface};
use crate::geometry::{ParametricCurve, ParametricSurface};
use crate::math::{Point, Point2, Vector};
use crate::topology::{CoEdgeId, Model, NodeKind, Orientation, Shape};

/// How a filling surface is made.
#[derive(Debug, Clone)]
pub struct FillingOptions {
    /// The order of the bending energy minimized (2: the thin plate; OCCT's `Degree`).
    pub order: usize,
    /// Points within the loop the surface passes through.
    pub points: Vec<Point>,
    /// How many points along an edge the surface passes through on average (OCCT's
    /// `NbPtsOnCur`): each edge takes its share by its length, between them the surface may
    /// stray from it, the face's tolerance covering that.
    pub points_per_edge: usize,
}

impl Default for FillingOptions {
    fn default() -> Self {
        FillingOptions { order: 3, points: Vec::new(), points_per_edge: 15 }
    }
}

/// Grid lines each way for the B-spline fitting the plate.
const GRID: usize = 41;
/// The most points along the loop the plate passes through.
const BUDGET: usize = 480;

/// A face spanning the closed loop of `edges` (in any order and direction), passing through
/// `options.points`.
#[doc(alias = "BRepOffsetAPI_MakeFilling")]
#[doc(alias = "GeomPlate_BuildPlateSurface")]
pub fn n_sided(model: &mut Model, edges: &[Shape], options: &FillingOptions) -> Result<Shape> {
    let loop_edges = ordered(model, edges)?;
    let wire = crate::make::wire_within(model, &loop_edges, 1e-5)?;
    let wire_id = crate::topology::WireId::try_from(wire.node).expect("a wire");
    if !model.is_closed_wire(wire_id) {
        return Err(Error::InvalidInput("the edges of a filling do not close a loop".into()));
    }
    // The loop's edges in order, each with its curve and the way the loop runs along it.
    let coedges: Vec<CoEdgeId> = model.explore(wire, NodeKind::CoEdge).map(|c| CoEdgeId::try_from(c.node).expect("a coedge")).collect();
    let mut runs = Vec::with_capacity(coedges.len());
    for &c in &coedges {
        let (edge, o) = (model[c].edge, model[c].orientation);
        let (curve, first, last) = model.edge_curve(Shape::from(edge)).ok_or_else(|| Error::InvalidInput("a degenerate edge in a filling's loop".into()))?;
        runs.push((edge, o, curve, first, last));
    }
    face_on_loop(model, &runs, options)
}

/// A filling face bounded by `runs` (edges, each with the way the loop runs along it, its
/// curve and range), in order around a closed loop, through `options.points`: on the edges
/// themselves (as a solid's cap shares its edges with the faces beside it), running the way
/// given (the face reversed where its surface faces the other way).
pub(crate) fn face_on_loop(model: &mut Model, runs: &[(crate::topology::EdgeId, Orientation, crate::geometry::curve::Curve, f64, f64)], options: &FillingOptions) -> Result<Shape> {
    let order = options.order.clamp(2, 4);
    // The loop as a polygon (20 points along each edge, as OCCT's average plane takes them),
    // and its mean plane (Newell's normal through the centroid).
    let per_edge = 20;
    let mut polygon = Vec::new();
    for (_, o, curve, first, last) in runs {
        let mut ts: Vec<f64> = (0..per_edge).map(|k| first + (last - first) * k as f64 / per_edge as f64).collect();
        if *o == Orientation::Reversed {
            ts = ts.iter().map(|t| first + last - t).collect();
        }
        polygon.extend(ts.into_iter().map(|t| curve.point(t)));
    }
    let n = polygon.len() as f64;
    let center = Point::ORIGIN + polygon.iter().fold(Vector::new(0.0, 0.0, 0.0), |s, p| s + (*p - Point::ORIGIN)) * (1.0 / n);
    let mut normal = Vector::new(0.0, 0.0, 0.0);
    for (i, p) in polygon.iter().enumerate() {
        let q = polygon[(i + 1) % polygon.len()];
        normal = normal + (*p - center).cross(q - center);
    }
    if normal.length() <= 1e-12 * polygon.iter().map(|p| p.distance(center)).fold(0.0, f64::max).powi(2) {
        return Err(Error::InvalidInput("a filling's loop encloses no area".into()));
    }
    // (The plane fitting all the constraints best, as OCCT's average plane: its normal the
    // least-spread direction of the loop's samples and the points; turned with the loop.)
    let normal = {
        let all: Vec<Point> = polygon.iter().chain(&options.points).copied().collect();
        let mean = Point::ORIGIN + all.iter().fold(Vector::new(0.0, 0.0, 0.0), |s, p| s + (*p - Point::ORIGIN)) * (1.0 / all.len() as f64);
        let mut c = [[0.0f64; 3]; 3];
        for p in &all {
            let d = *p - mean;
            let d = [d.x, d.y, d.z];
            for i in 0..3 {
                for j in 0..3 {
                    c[i][j] += d[i] * d[j];
                }
            }
        }
        let least = least_spread(c);
        let least = Vector::new(least[0], least[1], least[2]);
        if least.dot(normal) < 0.0 { -least } else { least }
    };
    // (In the plane, the axes the constraints spread most along and least, as OCCT's average
    // plane takes them: the plate is scaled along each by the constraints' extent.)
    let e1 = {
        let a = if normal.x.abs() < 0.9 { Vector::new(1.0, 0.0, 0.0) } else { Vector::new(0.0, 1.0, 0.0) };
        let a = (a - normal * a.dot(normal)).normalized();
        let b = normal.cross(a);
        let (mut cxx, mut cyy, mut cxy) = (0.0, 0.0, 0.0);
        let all: Vec<Point> = polygon.iter().chain(&options.points).copied().collect();
        let mean = Point::ORIGIN + all.iter().fold(Vector::new(0.0, 0.0, 0.0), |s, p| s + (*p - Point::ORIGIN)) * (1.0 / all.len() as f64);
        for p in &all {
            let (x, y) = ((*p - mean).dot(a), (*p - mean).dot(b));
            (cxx, cyy, cxy) = (cxx + x * x, cyy + y * y, cxy + x * y);
        }
        let angle = 0.5 * (2.0 * cxy).atan2(cxx - cyy);
        a * angle.cos() + b * angle.sin()
    };
    let e2 = normal.cross(e1);
    let flat = |p: Point| Point2::new((p - center).dot(e1), (p - center).dot(e2));
    let height = |p: Point| (p - center).dot(normal);
    // The constraints: points along the edges and the points given, each once. Along the
    // edges as OCCT's plate takes them: each edge its share of 1.2 times as many points as
    // asked for in all, by its length (3 at least), closer together towards its ends (by
    // the cosine of even steps along its length on the plane).
    // (Fewer in all with many edges: the plate's system grows with the cube of the points.)
    let size = polygon.iter().map(|p| p.distance(center)).fold(0.0, f64::max);
    let lengths: Vec<Vec<(f64, f64)>> = runs
        .iter()
        .map(|(_, _, curve, first, last)| {
            let k = 64;
            let mut acc = vec![(*first, 0.0)];
            for i in 1..=k {
                let (t0, t1) = (first + (last - first) * (i - 1) as f64 / k as f64, first + (last - first) * i as f64 / k as f64);
                let (a, b) = (flat(curve.point(t0)), flat(curve.point(t1)));
                acc.push((t1, acc[i - 1].1 + (b.x - a.x).hypot(b.y - a.y)));
            }
            acc
        })
        .collect();
    let total: f64 = lengths.iter().map(|l| l[l.len() - 1].1).sum();
    let wanted = (options.points_per_edge * runs.len()) as f64;
    let thin = (BUDGET as f64 / (1.2 * wanted).max(1.0)).min(1.0);
    let mut along = Vec::new();
    for ((_, _, curve, _, _), l) in runs.iter().zip(&lengths) {
        let length = l[l.len() - 1].1;
        let count = ((1.2 * wanted * thin * length / total.max(1e-300)) as usize).max(3);
        for j in 0..count {
            let s = length * (1.0 - (j as f64 * std::f64::consts::PI / (count - 1) as f64).cos()) / 2.0;
            // (The parameter at that length, between the samples about it.)
            let i = l.partition_point(|&(_, a)| a < s).clamp(1, l.len() - 1);
            let ((t0, a0), (t1, a1)) = (l[i - 1], l[i]);
            let t = if a1 > a0 { t0 + (t1 - t0) * (s - a0) / (a1 - a0) } else { t1 };
            along.push(curve.point(t));
        }
    }
    let mut sites: Vec<(Point2, f64)> = Vec::new();
    for p in along.iter().chain(&options.points) {
        let (q, h) = (flat(*p), height(*p));
        if !sites.iter().any(|(s, _)| (s.x - q.x).hypot(s.y - q.y) <= 1e-9 * size) {
            sites.push((q, h));
        }
    }
    let plate = Plate::through(&sites, order, size)?;
    let at = |q: Point2| center + e1 * q.x + e2 * q.y + normal * plate.at(q);
    // The surface fitted as a bicubic B-spline over the loop's box in the plane (a little
    // larger), its parameters the plane's coordinates.
    let (mut lo, mut hi) = (Point2::new(f64::INFINITY, f64::INFINITY), Point2::new(f64::NEG_INFINITY, f64::NEG_INFINITY));
    for (s, _) in &sites {
        lo = Point2::new(lo.x.min(s.x), lo.y.min(s.y));
        hi = Point2::new(hi.x.max(s.x), hi.y.max(s.y));
    }
    let margin = 0.02 * (hi.x - lo.x).max(hi.y - lo.y);
    let (lo, hi) = (Point2::new(lo.x - margin, lo.y - margin), Point2::new(hi.x + margin, hi.y + margin));
    let us: Vec<f64> = (0..GRID).map(|i| lo.x + (hi.x - lo.x) * i as f64 / (GRID - 1) as f64).collect();
    let vs: Vec<f64> = (0..GRID).map(|j| lo.y + (hi.y - lo.y) * j as f64 / (GRID - 1) as f64).collect();
    // Points interpolated along u on each grid line, then the poles along v.
    let rows: Vec<crate::geometry::curve::BSplineCurve> = vs
        .iter()
        .map(|&v| fit::interpolate_with_parameters(&us.iter().map(|&u| at(Point2::new(u, v))).collect::<Vec<_>>(), &us))
        .collect::<Result<_>>()?;
    let nu = rows[0].poles().len();
    let columns: Vec<crate::geometry::curve::BSplineCurve> = (0..nu)
        .map(|i| fit::interpolate_with_parameters(&rows.iter().map(|r| r.poles()[i]).collect::<Vec<_>>(), &vs))
        .collect::<Result<_>>()?;
    let poles: Vec<Vec<Point>> = columns.iter().map(|c| c.poles().to_vec()).collect();
    let knots = |c: &crate::geometry::curve::BSplineCurve| KnotVector::new(c.knots().to_vec(), c.multiplicities().to_vec(), c.degree());
    let surface = Surface::BSpline(BSplineSurface::try_new(poles, knots(&rows[0]), knots(&columns[0]))?);
    // The face: each edge's curve on the surface its projection onto the plane; the loop
    // counter-clockwise about the normal (else reversed), and its edges' tolerances raised to
    // how far they stray from the surface.
    let area: f64 = (0..polygon.len()).map(|i| {
        let (a, b) = (flat(polygon[i]), flat(polygon[(i + 1) % polygon.len()]));
        a.x * b.y - b.x * a.y
    }).sum();
    let mut uses = Vec::with_capacity(runs.len());
    let mut worst: f64 = 0.0;
    for (edge, o, curve, first, last) in runs {
        let k = 96;
        let ts: Vec<f64> = (0..=k).map(|i| first + (last - first) * i as f64 / k as f64).collect();
        let on_plane: Vec<Point2> = ts.iter().map(|&t| flat(curve.point(t))).collect();
        for (t, q) in ts.iter().zip(&on_plane) {
            worst = worst.max(surface.point(q.x, q.y).distance(curve.point(*t)));
        }
        let pcurve = Curve2::BSpline(fit::interpolate_with_parameters(&on_plane, &ts)?);
        uses.push((*edge, *o, pcurve, *first, *last));
    }
    let flipped = area < 0.0;
    if flipped {
        uses.reverse();
        for u in &mut uses {
            u.1 = u.1.reversed();
        }
    }
    let tolerance = (1.5 * worst).max(crate::tolerance::LINEAR);
    let face = model.new_empty(NodeKind::Face);
    model.set_face_surface(face, surface, tolerance)?;
    let mut coedges = Vec::with_capacity(uses.len());
    for (edge, o, pcurve, first, last) in uses {
        if model[edge].tolerance < tolerance {
            model.set_edge_tolerance(edge, tolerance)?;
        }
        let c = model.add_coedge(edge, o)?;
        model.set_pcurve_with_range(c, pcurve, first, last)?;
        coedges.push(c);
    }
    let w = model.add_wire(&coedges)?;
    model.add(face, Shape::from(w))?;
    for v in model.unique_subshapes(face, NodeKind::Vertex) {
        let id = crate::topology::VertexId::try_from(v.node).expect("a vertex");
        if model[id].tolerance < tolerance {
            model.set_vertex_tolerance(id, tolerance)?;
        }
    }
    Ok(if flipped { face.reversed() } else { face })
}

/// Copies of `edges` (each placed where it is: the face's own, joined at their ends without
/// changing the edges given, which other faces may share) in order end to end, by their end
/// points.
fn ordered(model: &mut Model, edges: &[Shape]) -> Result<Vec<Shape>> {
    let mut placed = Vec::with_capacity(edges.len());
    for &e in edges {
        let given: Vec<Shape> = if e.kind() == NodeKind::Wire { model.explore(e, NodeKind::Edge).collect() } else { vec![e] };
        for x in given {
            placed.push(crate::modify::copy(model, x, false)?);
        }
    }
    let ends = |model: &Model, e: Shape| {
        let (a, b) = model.edge_vertices(e);
        (model.point(a), model.point(b))
    };
    let size = placed.iter().map(|&e| { let (a, b) = ends(model, e); a.distance(b) }).fold(0.0, f64::max);
    // (Ends a little apart, as edges of sections crossing a seam may leave them, still meet.)
    let near = (1e-6 * size).max(1e-5);
    let mut out = vec![placed.remove(0)];
    let (first, mut tip) = ends(model, out[0]);
    while !placed.is_empty() {
        let next = placed.iter().position(|&e| { let (a, b) = ends(model, e); a.distance(tip) <= near || b.distance(tip) <= near });
        let Some(i) = next else {
            // (The first edge the other way round, its start the loop's tip.)
            if out.len() == 1 && placed.iter().any(|&e| { let (a, b) = ends(model, e); a.distance(first) <= near || b.distance(first) <= near }) {
                tip = first;
                continue;
            }
            return Err(Error::InvalidInput("the edges of a filling do not join end to end".into()));
        };
        let e = placed.remove(i);
        let (a, b) = ends(model, e);
        tip = if a.distance(tip) <= near { b } else { a };
        out.push(e);
    }
    Ok(out)
}

/// A polyharmonic spline over the plane: heights interpolating the sites, of least energy
/// of its order.
struct Plate {
    sites: Vec<Point2>,
    weights: Vec<f64>,
    /// The polynomial's coefficients, over monomials `x^a y^b` with `a + b < order`.
    polynomial: Vec<f64>,
    order: usize,
    /// Coordinates scaled by these, each way (the sites' box within about the unit square:
    /// the plate's energy measured so, as OCCT's is).
    scale: (f64, f64),
}

impl Plate {
    fn through(sites: &[(Point2, f64)], order: usize, size: f64) -> Result<Plate> {
        // (Each way by half the sites' extent that way.)
        let half = |f: &dyn Fn(&Point2) -> f64| {
            let (lo, hi) = sites.iter().map(|(p, _)| f(p)).fold((f64::INFINITY, f64::NEG_INFINITY), |(l, h), x| (l.min(x), h.max(x)));
            (0.5 * (hi - lo)).max(1e-9 * size).max(1e-300)
        };
        let scale = (1.0 / half(&|p| p.x), 1.0 / half(&|p| p.y));
        let at: Vec<Point2> = sites.iter().map(|(p, _)| Point2::new(p.x * scale.0, p.y * scale.1)).collect();
        let n = at.len();
        let monomials = Self::monomials(order);
        let m = monomials.len();
        if n < m + 1 {
            return Err(Error::InvalidInput("too few constraints for a filling".into()));
        }
        // [K P; Pᵀ 0] [w; c] = [h; 0]
        let size = n + m;
        let mut a = vec![vec![0.0; size]; size];
        let mut b = vec![0.0; size];
        for i in 0..n {
            for j in 0..n {
                a[i][j] = Self::kernel(at[i], at[j], order);
            }
            for (k, &(px, py)) in monomials.iter().enumerate() {
                let v = at[i].x.powi(px as i32) * at[i].y.powi(py as i32);
                a[i][n + k] = v;
                a[n + k][i] = v;
            }
            b[i] = sites[i].1;
        }
        let x = solve(a, b).ok_or_else(|| Error::Failed("a filling's plate could not be solved (constraints all on a curve of low degree?)".into()))?;
        Ok(Plate { sites: at, weights: x[..n].to_vec(), polynomial: x[n..].to_vec(), order, scale })
    }

    fn monomials(order: usize) -> Vec<(usize, usize)> {
        (0..order).flat_map(|d| (0..=d).map(move |b| (d - b, b))).collect()
    }

    fn kernel(p: Point2, q: Point2, order: usize) -> f64 {
        let r2 = (p.x - q.x).powi(2) + (p.y - q.y).powi(2);
        if r2 <= 0.0 {
            return 0.0;
        }
        // r^(2m-2) ln r = r2^(m-1) ln(r2) / 2
        r2.powi(order as i32 - 1) * 0.5 * r2.ln()
    }

    fn at(&self, p: Point2) -> f64 {
        let p = Point2::new(p.x * self.scale.0, p.y * self.scale.1);
        let mut h: f64 = self.sites.iter().zip(&self.weights).map(|(s, w)| w * Self::kernel(p, *s, self.order)).sum();
        for (k, &(px, py)) in Self::monomials(self.order).iter().enumerate() {
            h += self.polynomial[k] * p.x.powi(px as i32) * p.y.powi(py as i32);
        }
        h
    }
}

/// The unit eigenvector of the symmetric 3×3 matrix `c` with the least eigenvalue (Jacobi
/// rotations).
fn least_spread(mut c: [[f64; 3]; 3]) -> [f64; 3] {
    let mut v = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    for _ in 0..50 {
        let (mut p, mut q, mut big) = (0, 1, 0.0);
        for i in 0..3 {
            for j in i + 1..3 {
                if c[i][j].abs() > big {
                    (p, q, big) = (i, j, c[i][j].abs());
                }
            }
        }
        if big <= 1e-15 * (c[0][0].abs() + c[1][1].abs() + c[2][2].abs()) {
            break;
        }
        let theta = 0.5 * (2.0 * c[p][q]).atan2(c[q][q] - c[p][p]);
        let (sn, cs) = theta.sin_cos();
        let mut r = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        r[p][p] = cs;
        r[q][q] = cs;
        r[p][q] = sn;
        r[q][p] = -sn;
        // c = rᵀ c r, v = v r
        let mul = |a: [[f64; 3]; 3], b: [[f64; 3]; 3]| {
            let mut m = [[0.0; 3]; 3];
            for i in 0..3 {
                for j in 0..3 {
                    m[i][j] = (0..3).map(|k| a[i][k] * b[k][j]).sum();
                }
            }
            m
        };
        let rt = [[r[0][0], r[1][0], r[2][0]], [r[0][1], r[1][1], r[2][1]], [r[0][2], r[1][2], r[2][2]]];
        c = mul(mul(rt, c), r);
        v = mul(v, r);
    }
    let k = (0..3).min_by(|&i, &j| c[i][i].total_cmp(&c[j][j])).expect("three");
    [v[0][k], v[1][k], v[2][k]]
}

/// The solution of the dense system `a x = b` (LU with partial pivoting), if not singular.
fn solve(mut a: Vec<Vec<f64>>, mut b: Vec<f64>) -> Option<Vec<f64>> {
    let n = b.len();
    let largest = a.iter().flatten().fold(0.0f64, |m, v| m.max(v.abs()));
    for k in 0..n {
        let p = (k..n).max_by(|&i, &j| a[i][k].abs().total_cmp(&a[j][k].abs()))?;
        if a[p][k].abs() <= 1e-300 * largest.max(1.0) || !a[p][k].is_finite() {
            return None;
        }
        a.swap(k, p);
        b.swap(k, p);
        for i in k + 1..n {
            let f = a[i][k] / a[k][k];
            if f != 0.0 {
                for j in k..n {
                    a[i][j] -= f * a[k][j];
                }
                b[i] -= f * b[k];
            }
        }
    }
    let mut x = vec![0.0; n];
    for k in (0..n).rev() {
        let s: f64 = (k + 1..n).map(|j| a[k][j] * x[j]).sum();
        x[k] = (b[k] - s) / a[k][k];
    }
    Some(x)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn square(model: &mut Model, lift: f64) -> Vec<Shape> {
        let corners = [Point::new(0.0, 0.0, 0.0), Point::new(10.0, 0.0, lift), Point::new(10.0, 10.0, 0.0), Point::new(0.0, 10.0, lift)];
        // (Given out of order, one the other way round.)
        [(0, 1), (2, 3), (2, 1), (3, 0)].iter().map(|&(i, j)| crate::make::edge_between(model, corners[i], corners[j]).unwrap()).collect()
    }

    /// A 20 x 10 rectangle's filling through a point 3 above its middle: its plate measured
    /// along the rectangle's own axes, each scaled by its half width, as OCCT's is (OCCT's
    /// heights at points within it, and its area 221.3898).
    #[test]
    fn a_long_loop_fills_as_occt_measures_its_plate() {
        let mut model = Model::new();
        let corners = [Point::new(0.0, 0.0, 0.0), Point::new(20.0, 0.0, 0.0), Point::new(20.0, 10.0, 0.0), Point::new(0.0, 10.0, 0.0)];
        let edges: Vec<Shape> = (0..4).map(|i| crate::make::edge_between(&mut model, corners[i], corners[(i + 1) % 4]).unwrap()).collect();
        let options = FillingOptions { points: vec![Point::new(10.0, 5.0, 3.0)], ..Default::default() };
        let face = n_sided(&mut model, &edges, &options).unwrap();
        let surface = model.surface(face).unwrap().clone();
        for (x, y, z) in [(10.0, 5.0, 3.0), (5.0, 5.0, 1.8935), (15.0, 2.5, 1.2191), (3.0, 7.0, 0.8482)] {
            let (_, _, d) = crate::geometry::extrema::nearest_on_surface(Point::new(x, y, z), &surface, None);
            assert!(d < 2e-4, "{x} {y}: {d}");
        }
        let area = crate::properties::area(&model, face).mass;
        assert!((area - 221.3898).abs() < 1e-3, "{area}");
    }

    #[test]
    fn a_flat_loop_fills_with_its_plane() {
        let mut model = Model::new();
        let edges = square(&mut model, 0.0);
        let face = n_sided(&mut model, &edges, &FillingOptions::default()).unwrap();
        assert!(crate::validity::check(&model, face).is_valid());
        let area = crate::properties::area(&model, face).mass;
        assert!((area - 100.0).abs() < 1e-6, "{area}");
    }

    #[test]
    fn fillings_pass_through_their_points_and_edges() {
        let mut model = Model::new();
        let edges = square(&mut model, 2.0);
        let top = Point::new(5.0, 5.0, 4.0);
        let face = n_sided(&mut model, &edges, &FillingOptions { points: vec![top], ..Default::default() }).unwrap();
        assert!(crate::validity::check(&model, face).is_valid());
        let surface = model.surface(face).unwrap().clone();
        let (u, v, d) = crate::geometry::extrema::nearest_on_surface(top, &surface, None);
        assert!(d < 1e-4, "{d} at {u} {v}");
        // The edges lie on it, within the face's tolerance.
        assert!(model.tolerance(face) < 1e-4, "{}", model.tolerance(face));
    }
}
