//! The operands' topology as the boolean operation sees it: vertices merged across both
//! operands, edges with their split points, faces with their boundaries in parameter space.

use crate::hashing::FastMap as HashMap;

use crate::error::{Error, Result};
use crate::geometry::curve::{Circle, Curve};
use crate::geometry::curve2::Curve2;
use crate::geometry::surface::{RevolutionSurface, Surface, Torus};
use crate::geometry::{ParametricCurve, ParametricSurface};
use crate::math::{BoundingBox, Bvh, Direction, Frame, Point, Transform, Transformable, Vector};
use crate::topology::{CoEdgeId, EdgeId, FaceId, Model, NodeId, NodeKind, Orientation, Shape, VertexId};

use super::geom::{self, Uv, on_tol};

/// Vertices, merged when closer than [`on_tol`] (the distance at which a point counts as
/// lying on a curve, so that the same point is never two vertices).
#[derive(Default)]
pub(super) struct Vertices {
    pub points: Vec<Point>,
    grid: crate::hashing::FastMap<[i64; 3], Vec<usize>>,
}

const CELL: f64 = 1.0e-5;

impl Vertices {
    fn cell(p: Point) -> [i64; 3] {
        [(p.x / CELL).floor() as i64, (p.y / CELL).floor() as i64, (p.z / CELL).floor() as i64]
    }

    /// The vertex at `p`: an existing one within [`on_tol`], or a new one.
    pub fn add(&mut self, p: Point) -> usize {
        if let Some(i) = self.find(p, on_tol()) {
            return i;
        }
        let i = self.points.len();
        self.points.push(p);
        self.grid.entry(Self::cell(p)).or_default().push(i);
        i
    }

    /// The vertex closest to `p` within `tolerance` (at most a grid cell).
    pub fn find(&self, p: Point, tolerance: f64) -> Option<usize> {
        let c = Self::cell(p);
        // (Only the cells within `tolerance` of `p`: mostly its own.)
        let (lo, hi) = (Self::cell(p - Vector::new(tolerance, tolerance, tolerance)), Self::cell(p + Vector::new(tolerance, tolerance, tolerance)));
        let near = |k: usize| (lo[k].max(c[k] - 1) - c[k])..=(hi[k].min(c[k] + 1) - c[k]);
        let mut best: Option<(usize, f64)> = None;
        for dx in near(0) {
            for dy in near(1) {
                for dz in near(2) {
                    let Some(ids) = self.grid.get(&[c[0] + dx, c[1] + dy, c[2] + dz]) else { continue };
                    for &i in ids {
                        let d = self.points[i].distance(p);
                        if d <= tolerance && best.is_none_or(|(_, bd)| d < bd) {
                            best = Some((i, d));
                        }
                    }
                }
            }
        }
        best.map(|(i, _)| i)
    }
}

/// What a face's samples give: its boundary in (u, v), parameter ranges and boxes (see
/// [`Face`]).
#[derive(Clone)]
pub(super) struct FaceSample {
    pub polylines: Vec<Vec<Uv>>,
    pub u_range: (f64, f64),
    pub v_range: (f64, f64),
    pub bbox: (Point, Point),
    pub extent: (Point, Point),
}

/// The box of an edge's samples, kept in the model for the model's edge.
struct EdgeBox((Point, Point));

/// The box of the model's (unplaced) edge `id`'s samples, kept in the model with it.
/// Makes sure the model keeps an edge's box ([`model_edge_box`]), fetching its curve only
/// if it does not yet.
/// (Edges with curves only: a degenerate edge's box is its vertex, as the graph merges it.)
pub(super) fn warm_edge_box(model: &Model, id: EdgeId) {
    let edge = Shape::from(id);
    if !model.has_edge_curve(edge) {
        return;
    }
    model.derived(id.into(), model.edge_fingerprint(id), || {
        let (curve, first, last) = model.edge_curve(edge).expect("an edge with a curve");
        let start = model.point(model.edge_vertices(edge).0);
        EdgeBox(Edge::sampled_box(Some(&curve), first, last, start))
    });
}

pub(super) fn model_edge_box(model: &Model, id: EdgeId, curve: Option<&Curve>, first: f64, last: f64, start: Point) -> (Point, Point) {
    model.derived(id.into(), model.edge_fingerprint(id), || EdgeBox(Edge::sampled_box(curve, first, last, start))).0
}

/// An edge: of an operand, or a section edge where faces of the two operands meet.
pub(super) struct Edge {
    /// `None` for a degenerate edge (collapsed to its vertex, e.g. at a sphere's pole).
    pub curve: Option<Curve>,
    pub first: f64,
    pub last: f64,
    pub start: usize,
    pub end: usize,
    /// Split points: (parameter, vertex), the ends included once finalized.
    pub paves: Vec<(f64, usize)>,
    pub bbox: (Point, Point),
    /// The operand it comes from; `None` for section edges.
    pub operand: Option<usize>,
    /// Not on a face of its operand: of a wire or an edge argument.
    pub free: bool,
    /// For a free edge, the wire argument it belongs to (numbered across both operands).
    pub container: Option<usize>,
    /// For a free edge, whether its wire uses it reversed.
    pub reversed: bool,
    /// The operand's edge it is (unplaced), which the result may reuse whole.
    pub origin: Option<EdgeId>,
    /// The operand's edge it is, however placed.
    pub source: Option<EdgeId>,
    /// How far it may stray from the surfaces of its faces: its tolerance in the model (0
    /// for section edges).
    pub tolerance: f64,
    /// Its samples and chords for finding crossings (made on first use; the curve and its
    /// range do not change).
    pub chords: std::sync::OnceLock<Chords>,
}

/// An edge's sample parameters and points, the boxes of its chords widened by half their
/// length, and the boxes holding each run of [`Chords::RUN`] of those.
pub(super) struct Chords {
    pub params: Vec<f64>,
    pub points: Vec<Point>,
    pub boxes: Vec<([f64; 3], [f64; 3])>,
    pub runs: Vec<([f64; 3], [f64; 3])>,
}

impl Chords {
    pub const RUN: usize = 8;

    pub fn of(curve: &Curve, first: f64, last: f64) -> Chords {
        let params = geom::sample_parameters(Some(curve), first, last);
        let points = curve.points_along(&params);
        let chord_box = |a: Point, b: Point| {
            let r = 0.5 * a.distance(b);
            ([a.x.min(b.x) - r, a.y.min(b.y) - r, a.z.min(b.z) - r], [a.x.max(b.x) + r, a.y.max(b.y) + r, a.z.max(b.z) + r])
        };
        let boxes: Vec<([f64; 3], [f64; 3])> = points.windows(2).map(|w| chord_box(w[0], w[1])).collect();
        let runs = boxes
            .chunks(Self::RUN)
            .map(|run| {
                run.iter().fold(([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]), |(lo, hi), (l, h)| {
                    (std::array::from_fn(|k| lo[k].min(l[k])), std::array::from_fn(|k| hi[k].max(h[k])))
                })
            })
            .collect();
        Chords { params, points, boxes, runs }
    }
}

impl Edge {
    pub fn new(curve: Option<Curve>, first: f64, last: f64, start: usize, end: usize, vertices: &Vertices, operand: Option<usize>) -> Edge {
        let bbox = Edge::sampled_box(curve.as_ref(), first, last, vertices.points[start]);
        Edge::with_box(curve, first, last, start, end, bbox, operand)
    }

    pub fn with_box(curve: Option<Curve>, first: f64, last: f64, start: usize, end: usize, bbox: (Point, Point), operand: Option<usize>) -> Edge {
        Edge { curve, first, last, start, end, paves: Vec::new(), bbox, operand, free: false, container: None, reversed: false, origin: None, source: None, tolerance: 0.0, chords: Default::default() }
    }

    /// The box of the curve's samples (the start point for a degenerate edge).
    pub fn sampled_box(curve: Option<&Curve>, first: f64, last: f64, start: Point) -> (Point, Point) {
        let Some(c) = curve else { return geom::bounding_box(vec![start], geom::LOOSE_ON_TOL) };
        let points = c.points_along(&geom::sample_parameters(Some(c), first, last));
        // Chords may miss a little of the curve; the margin covers it: a hundredth of the
        // samples' extent (of the box's, as a closed curve's ends are no measure of it).
        // (Kept in the model: wide enough for the looser tolerances too.)
        let (lo, hi) = geom::bounding_box(points.iter().copied(), 0.0);
        geom::bounding_box(points, geom::LOOSE_ON_TOL + 0.01 * lo.distance(hi))
    }

    pub fn point(&self, t: f64, vertices: &Vertices) -> Point {
        match &self.curve {
            Some(c) => c.point(t),
            None => vertices.points[self.start],
        }
    }
}

/// An edge of a face's boundary, with its curve in the face's parameter space (same
/// parametrization as the edge).
pub(super) struct BoundaryUse {
    pub edge: usize,
    pub pcurve: Curve2,
    /// False for an edge inside the face (internal or external), bounding no region of it.
    pub bounds: bool,
    /// Whether the face runs along the edge (its region on the left of the pcurve as it
    /// goes), else against it.
    pub along: bool,
}

/// A face of an operand.
pub(super) struct Face {
    pub operand: usize,
    /// On the boundary of a solid of its operand (else a face or shell argument).
    pub solid: bool,
    /// For a face not on a solid, the shell argument it belongs to (numbered across both
    /// operands).
    pub container: Option<usize>,
    /// The surface the operation computes with.
    pub surface: Surface,
    /// The surface result faces get: the face's own where `surface` stands in for it (a
    /// torus is computed as a surface of revolution), else `surface`.
    pub output: Surface,
    /// +1 if `∂u × ∂v` points out of the operand, -1 if into it.
    pub outward: f64,
    pub boundary: Vec<BoundaryUse>,
    /// The boundary sampled in (u, v), for locating points.
    pub polylines: Vec<Vec<Uv>>,
    pub u_period: Option<f64>,
    /// The period of `v` (a surface of revolution of a circle), if periodic.
    pub v_period: Option<f64>,
    pub u_range: (f64, f64),
    pub v_range: (f64, f64),
    /// The surface's bounds in u and in v where it closes round across them without a
    /// period (a clamped B-spline swept round a tube): a point on that seam has either.
    pub closed_bounds: [Option<(f64, f64)>; 2],
    pub bbox: (Point, Point),
    /// A box holding the face, nearly tight: its grid's points, widened on a curved surface
    /// by more than the surface bulges between them.
    pub extent: (Point, Point),
    /// Vertices where edges of the other operand meet the face.
    pub vertices_on: Vec<usize>,
    /// Section edges lying on the face.
    pub sections: Vec<usize>,
    /// Edges of the other operand lying on the face's surface.
    pub on_edges: Vec<usize>,
    /// The operand's face it is (unplaced, with its own surface), which the result may reuse
    /// where it keeps the face whole.
    pub origin: Option<FaceId>,
    /// The operand's face whose regions (when untouched) it may keep with it, however placed.
    pub regions_key: Option<FaceId>,
    /// Its place among the operands' faces (far ones included).
    pub seq: usize,
    /// How far its edges may stray from its surface (their tolerances, read from a file).
    pub drift: f64,
    /// Points of a grid over the face's traced range and their largest spacing, as section
    /// tracing samples it to find the part near another face (made on first use).
    pub samples: std::sync::OnceLock<(Vec<Vec<Point>>, f64)>,
    /// Its boundary's polylines indexed for telling points inside (made on first use).
    pub inside_index: std::sync::OnceLock<geom::InsidePolylines>,
    /// Boxes of the cells of a grid over its surface's range, each grown by the grid's
    /// spacing, indexed: what passes none is clear of the surface there (made on first use,
    /// for free-form faces).
    pub cells: std::sync::OnceLock<crate::math::Bvh>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum Location {
    Inside,
    Boundary,
    Outside,
}

impl Face {
    /// The parameters of a point of the surface, with the periodic `u` brought into the
    /// face's range where possible.
    pub fn uv(&self, p: Point) -> Uv {
        let [u, v] = geom::uv_of(&self.surface, p);
        [self.shift_u(u), self.shift_v(v)]
    }

    /// `u` moved by whole periods into (or next to) the face's u range.
    pub fn shift_u(&self, u: f64) -> f64 {
        shift(u, self.u_period, self.u_range)
    }

    /// `v` moved by whole periods into (or next to) the face's v range.
    pub fn shift_v(&self, v: f64) -> f64 {
        shift(v, self.v_period, self.v_range)
    }

    /// `uv` moved by whole periods to be nearest `near` (or, on a seam of a surface closed
    /// without a period, to the bound nearer `near`).
    pub fn unwrap(&self, mut uv: Uv, near: Uv) -> Uv {
        for (k, period) in [(0, self.u_period), (1, self.v_period)] {
            if let Some(period) = period {
                uv[k] += ((near[k] - uv[k]) / period).round() * period;
            } else if let Some((lo, hi)) = self.closed_bounds[k] {
                let eps = 1e-9 * (hi - lo).abs();
                if (uv[k] - lo).abs() <= eps && (near[k] - hi).abs() < (near[k] - lo).abs() {
                    uv[k] = hi;
                } else if (uv[k] - hi).abs() <= eps && (near[k] - lo).abs() < (near[k] - hi).abs() {
                    uv[k] = lo;
                }
            }
        }
        uv
    }

    /// The copies of `uv` (shifted by whole periods) within the face's parameter ranges.
    pub fn representatives(&self, uv: Uv) -> impl Iterator<Item = Uv> {
        representatives(uv, [self.u_period, self.v_period], self.u_range, self.v_range)
    }

    /// Whether `uv` lies inside the face's boundary polylines ([`geom::inside_polylines`]; a
    /// long boundary indexed, as a face may be asked about many points).
    pub fn inside_polylines(&self, uv: Uv) -> bool {
        if self.polylines.iter().map(Vec::len).sum::<usize>() < 64 {
            return geom::inside_polylines(&self.polylines, uv);
        }
        self.inside_index.get_or_init(|| geom::InsidePolylines::new(&self.polylines)).contains(uv)
    }

    /// Where `p` (on the face's surface) lies with respect to the face: on its boundary
    /// within [`on_tol`], else inside or outside.
    pub fn locate(&self, p: Point, edges: &[Edge], vertices: &Vertices) -> Location {
        if self.near_boundary(p, on_tol(), edges, vertices) {
            return Location::Boundary;
        }
        // A face without a boundary is its whole surface (an unbounded plane).
        if self.boundary.is_empty() {
            return Location::Inside;
        }
        // Off the surface by more than sections stray from it (fitted to a thousandth of
        // their size): beyond its edge.
        let reach = 1e-3 * self.extent.0.distance(self.extent.1) + on_tol();
        let Some(uv) = geom::uv_within(&self.surface, p, reach) else { return Location::Outside };
        let inside_at = |uv: Uv| self.representatives(uv).into_iter().any(|uv| self.inside_polylines(uv));
        if inside_at(uv) {
            return Location::Inside;
        }
        // On a seam (both its copies in the polylines, neither side counting it inside):
        // inside if a step either way across it is (a true boundary has outside on one side).
        for (k, period) in [(0, self.u_period), (1, self.v_period)] {
            let Some(period) = period else { continue };
            let step = 1e-7 * period;
            let nudged = |d: f64| {
                let mut q = uv;
                q[k] += d;
                q
            };
            if inside_at(nudged(step)) && inside_at(nudged(-step)) {
                return Location::Inside;
            }
        }
        Location::Outside
    }

    /// Whether `p` (on the face's surface) lies on the face's side of its nearest boundary
    /// edge, where that is within `near` (and not a seam, nor nearest at an end of it): by the
    /// edge's own curve, where the boundary's polylines, their chords cutting across its
    /// bends, cannot tell a point that close.
    pub fn inside_by_nearest_edge(&self, p: Point, near: f64, edges: &[Edge]) -> Option<bool> {
        let mut nearest: Option<(f64, usize, f64)> = None;
        for (k, b) in self.boundary.iter().enumerate().filter(|(_, b)| b.bounds) {
            let e = &edges[b.edge];
            let Some(curve) = &e.curve else { continue };
            let (lo, hi) = e.bbox;
            let outside = |x: f64, a: f64, b: f64| x < a - near || x > b + near;
            if outside(p.x, lo.x, hi.x) || outside(p.y, lo.y, hi.y) || outside(p.z, lo.z, hi.z) {
                continue;
            }
            let (t, d) = geom::project(curve, e.first, e.last, p);
            if d <= near && nearest.is_none_or(|(best, _, _)| d < best) {
                nearest = Some((d, k, t));
            }
        }
        let (_, k, t) = nearest?;
        let b = &self.boundary[k];
        let e = &edges[b.edge];
        let span = e.last - e.first;
        if self.boundary.iter().filter(|c| c.edge == b.edge).count() > 1 || (t - e.first).abs() <= 1e-6 * span.abs() || (e.last - t).abs() <= 1e-6 * span.abs() {
            return None;
        }
        let (q, d) = b.pcurve.d1(t);
        let foot = [q.x, q.y];
        let uv = geom::uv_of(&self.surface, p);
        // (The copy of the point's parameters nearest the foot's, across periods.)
        let near_copy = |x: f64, f: f64, period: Option<f64>| period.map_or(x, |t| x + ((f - x) / t).round() * t);
        let uv = [near_copy(uv[0], foot[0], self.u_period), near_copy(uv[1], foot[1], self.v_period)];
        let cross = d.x * (uv[1] - foot[1]) - d.y * (uv[0] - foot[0]);
        if cross == 0.0 {
            return None;
        }
        Some((cross > 0.0) == b.along)
    }

    /// Whether `p` lies within `near` of one of the face's edges.
    pub fn near_boundary(&self, p: Point, near: f64, edges: &[Edge], vertices: &Vertices) -> bool {
        self.boundary.iter().filter(|b| b.bounds).any(|b| {
            let e = &edges[b.edge];
            let (lo, hi) = e.bbox;
            let outside = |x: f64, a: f64, b: f64| x < a - near || x > b + near;
            if outside(p.x, lo.x, hi.x) || outside(p.y, lo.y, hi.y) || outside(p.z, lo.z, hi.z) {
                return false;
            }
            let d = match &e.curve {
                Some(c) => geom::project(c, e.first, e.last, p).1,
                None => vertices.points[e.start].distance(p),
            };
            d <= near
        })
    }

    /// The outward unit normal at a point of the face.
    pub fn outward_normal(&self, p: Point) -> crate::math::Vector {
        geom::normal_at_point(&self.surface, p) * self.outward
    }
}

/// The copies of `(u, v)` by whole periods that fall within a face's parameter ranges
/// (a face winding round more than once, as a band round a cylinder, holds several).
/// (At most four each way, in a fixed array: asked for every point located.)
pub(super) fn representatives([u, v]: Uv, periods: [Option<f64>; 2], u_range: (f64, f64), v_range: (f64, f64)) -> impl Iterator<Item = Uv> {
    let copies = |x: f64, period: Option<f64>, (lo, hi): (f64, f64)| -> ([f64; 4], usize) {
        match period {
            Some(p) => {
                let mut out = ([0.0; 4], 0);
                // (From the lowest copy in the range up.)
                let mut y = x + ((lo - 1e-9 - x) / p).ceil() * p;
                while y <= hi + 1e-9 && out.1 < 4 {
                    out.0[out.1] = y;
                    out.1 += 1;
                    y += p;
                }
                out
            }
            None => ([x, 0.0, 0.0, 0.0], 1),
        }
    };
    let (us, nu) = copies(u, periods[0], u_range);
    let (vs, nv) = copies(v, periods[1], v_range);
    us.into_iter().take(nu).flat_map(move |u| vs.into_iter().take(nv).map(move |v| [u, v]))
}

/// The surface of revolution of the torus's meridian circle (at u = 0), parametrized like
/// the torus: v runs around the tube from the outer equator towards the frame's z, and u
/// turns x towards y (about the frame's z, or against it where the frame is left-handed).
pub(super) fn torus_as_revolution(t: &Torus) -> Surface {
    let f = t.frame();
    let (x, y, z) = (Vector::from(f.x()), Vector::from(f.y()), Vector::from(f.z()));
    let center = f.origin() + x * t.major_radius();
    let meridian = Circle::new(Frame::new(center, Direction::from_vector(x.cross(z)), f.x()), t.minor_radius());
    let axis = crate::math::Axis::new(f.origin(), Direction::from_vector(x.cross(y)));
    Surface::Revolution(RevolutionSurface::new(Curve::Circle(meridian), axis))
}

/// The plane of a B-spline surface whose poles all lie in one, oriented like the
/// surface's `∂u × ∂v`.
pub(super) fn planar_bspline(surface: &Surface) -> Option<Surface> {
    let Surface::BSpline(b) = surface else { return None };
    let poles: Vec<Point> = b.poles().iter().flatten().copied().collect();
    let (u1, u2, v1, v2) = b.bounds();
    let (p, su, sv) = b.d1(0.5 * (u1 + u2), 0.5 * (v1 + v2));
    let normal = crate::math::Direction::try_from_vector(su.cross(sv))?;
    let plane = crate::geometry::surface::Plane::new(p, normal);
    poles.iter().all(|&q| plane.distance(q) <= crate::tolerance::LINEAR).then_some(Surface::Plane(plane))
}

/// `x` moved by whole periods into (or next to) `range`.
fn shift(x: f64, period: Option<f64>, (lo, hi): (f64, f64)) -> f64 {
    match period {
        Some(period) => x + ((0.5 * (lo + hi) - x) / period).round() * period,
        None => x,
    }
}

/// The operands' vertices, edges and faces (two operands, or more for a general fuse).
pub(super) struct Graph {
    pub vertices: Vertices,
    pub edges: Vec<Edge>,
    pub faces: Vec<Face>,
    /// Per operand: whether its shell is open (a half-space bounded by its faces).
    pub open: Vec<bool>,
    /// Per operand: its tight bounding box.
    pub tight_boxes: Vec<Option<(Point, Point)>>,
    /// Per operand: its dimension (3 with solids, 2 with faces, 1 with edges, 0 with
    /// vertices only), as OCCT's boolean operations define it.
    pub dimension: Vec<u8>,
    /// Vertex arguments: (operand, vertex).
    pub free_vertices: Vec<(usize, usize)>,
    /// The operands' vertices (unplaced) the graph's vertices are, which the result reuses.
    pub vertex_origin: HashMap<usize, VertexId>,
    /// Faces left out, far from the other operand.
    pub far: Vec<super::far::FarFace>,
    /// Per operand: a box holding each of its faces, in the order `Model::explore` lists
    /// them (graph faces and far faces alike).
    pub face_bounds: Vec<Vec<BoundingBox>>,
    /// The number of the operands' faces met so far.
    seq: usize,
    /// The number of wire and shell arguments seen.
    containers: usize,
    /// Trees over the operands' faces and edges, built once both are added.
    index: std::sync::OnceLock<Index>,
    /// Per operand: its faces' indices, ascending (found once).
    operand_faces: Vec<std::sync::OnceLock<Vec<usize>>>,
    /// Per operand: those of its faces on a solid, ascending (found once).
    operand_solid_faces: Vec<std::sync::OnceLock<Vec<usize>>>,
}

impl Default for Graph {
    fn default() -> Graph {
        Graph::with_operands(2)
    }
}

impl Graph {
    /// An empty graph for `n` operands.
    pub fn with_operands(n: usize) -> Graph {
        Graph {
            vertices: Vertices::default(),
            edges: Vec::new(),
            faces: Vec::new(),
            open: vec![false; n],
            tight_boxes: vec![None; n],
            dimension: vec![0; n],
            free_vertices: Vec::new(),
            vertex_origin: HashMap::default(),
            far: Vec::new(),
            face_bounds: vec![Vec::new(); n],
            seq: 0,
            containers: 0,
            index: std::sync::OnceLock::new(),
            operand_faces: (0..n).map(|_| std::sync::OnceLock::new()).collect(),
            operand_solid_faces: (0..n).map(|_| std::sync::OnceLock::new()).collect(),
        }
    }

    /// The number of operands.
    pub fn operands(&self) -> usize {
        self.open.len()
    }
}

/// Bounding volume hierarchies over each operand's faces and edges (those of the operands,
/// not section edges), for finding candidates without testing every face or edge; with more
/// than two operands, over all of them too.
pub(super) struct Index {
    faces: Vec<Bvh>,
    face_ids: Vec<Vec<usize>>,
    edges: Vec<Bvh>,
    edge_ids: Vec<Vec<usize>>,
    far: Vec<Bvh>,
    far_ids: Vec<Vec<usize>>,
    /// (More than two operands: all faces and edges, with each one's operand.)
    all_faces: Option<(Bvh, Vec<usize>)>,
    all_edges: Option<(Bvh, Vec<usize>)>,
    face_operand: Vec<usize>,
    edge_operand: Vec<Option<usize>>,
}

/// True if `shape` is a solid, or a compound of nothing but solids (so all its edges and
/// vertices lie on faces).
pub(super) fn made_of_solids(model: &Model, shape: Shape) -> bool {
    match shape.node.kind {
        NodeKind::Solid => true,
        NodeKind::Compound | NodeKind::CompSolid => {
            let mut children = model.children(shape).peekable();
            children.peek().is_some() && children.all(|c| made_of_solids(model, c))
        }
        _ => false,
    }
}

pub(super) fn to_box(b: &(Point, Point)) -> BoundingBox {
    BoundingBox::new(b.0, b.1)
}

impl Index {
    fn ids(ids: &[usize], mut found: Vec<usize>) -> Vec<usize> {
        for i in &mut found {
            *i = ids[*i];
        }
        found.sort_unstable();
        found
    }

    /// Faces of operand `k` whose box overlaps `b`, in increasing order.
    pub fn faces_overlapping(&self, k: usize, b: &(Point, Point)) -> Vec<usize> {
        Self::ids(&self.face_ids[k], self.faces[k].overlapping(&to_box(b)))
    }

    /// Faces of operand `k` whose box contains `p`, in increasing order.
    pub fn faces_containing(&self, k: usize, p: Point) -> Vec<usize> {
        Self::ids(&self.face_ids[k], self.faces[k].containing(p))
    }

    /// Faces of operand `k` whose box the ray crosses, in increasing order.
    pub fn faces_along_ray(&self, k: usize, origin: Point, direction: Vector) -> Vec<usize> {
        Self::ids(&self.face_ids[k], self.faces[k].along_ray(origin, direction))
    }

    /// Far faces of operand `k` whose box the ray crosses, in increasing order.
    pub fn far_along_ray(&self, k: usize, origin: Point, direction: Vector) -> Vec<usize> {
        Self::ids(&self.far_ids[k], self.far[k].along_ray(origin, direction))
    }

    /// Edges of operand `k` whose box contains `p`, in increasing order.
    pub fn edges_containing(&self, k: usize, p: Point) -> Vec<usize> {
        Self::ids(&self.edge_ids[k], self.edges[k].containing(p))
    }

    /// Faces of operands other than `k` whose box overlaps `b`, in increasing order.
    pub fn faces_overlapping_others(&self, k: usize, b: &(Point, Point)) -> Vec<usize> {
        match &self.all_faces {
            None => self.faces_overlapping(1 - k, b),
            Some((tree, ids)) => {
                let mut found = tree.overlapping(&to_box(b));
                found.iter_mut().for_each(|i| *i = ids[*i]);
                found.retain(|&f| self.face_operand[f] != k);
                found.sort_unstable();
                found
            }
        }
    }

    /// Faces of operands other than `k` whose box contains `p`, in increasing order.
    pub fn faces_containing_others(&self, k: usize, p: Point) -> Vec<usize> {
        match &self.all_faces {
            None => self.faces_containing(1 - k, p),
            Some((tree, ids)) => {
                let mut found = tree.containing(p);
                found.iter_mut().for_each(|i| *i = ids[*i]);
                found.retain(|&f| self.face_operand[f] != k);
                found.sort_unstable();
                found
            }
        }
    }

    /// Edges of operands other than `k` whose box contains `p`, in increasing order.
    pub fn edges_containing_others(&self, k: usize, p: Point) -> Vec<usize> {
        match &self.all_edges {
            None => self.edges_containing(1 - k, p),
            Some((tree, ids)) => {
                let mut found = tree.containing(p);
                found.iter_mut().for_each(|i| *i = ids[*i]);
                found.retain(|&e| self.edge_operand[e] != Some(k));
                found.sort_unstable();
                found
            }
        }
    }

    /// Pairs of faces of different operands with overlapping boxes (the one of the lower
    /// operand first), in increasing order.
    pub fn face_pairs(&self) -> Vec<(usize, usize)> {
        let mut pairs: Vec<(usize, usize)> = match &self.all_faces {
            None => self.faces[0].overlapping_pairs(&self.faces[1]).into_iter().map(|(a, b)| (self.face_ids[0][a], self.face_ids[1][b])).collect(),
            Some((tree, ids)) => {
                // (The tree against itself may give a pair either way round, or both.)
                let mut found: Vec<(usize, usize)> = tree
                    .overlapping_pairs(tree)
                    .into_iter()
                    .map(|(a, b)| (ids[a], ids[b]))
                    .filter(|&(a, b)| self.face_operand[a] != self.face_operand[b])
                    .map(|(a, b)| if self.face_operand[a] < self.face_operand[b] { (a, b) } else { (b, a) })
                    .collect();
                found.sort_unstable();
                found.dedup();
                found
            }
        };
        pairs.sort_unstable();
        pairs
    }

    /// Pairs of edges of different operands with overlapping boxes (the one of the lower
    /// operand first), in increasing order.
    pub fn edge_pairs(&self) -> Vec<(usize, usize)> {
        let mut pairs: Vec<(usize, usize)> = match &self.all_edges {
            None => self.edges[0].overlapping_pairs(&self.edges[1]).into_iter().map(|(a, b)| (self.edge_ids[0][a], self.edge_ids[1][b])).collect(),
            Some((tree, ids)) => {
                let mut found: Vec<(usize, usize)> = tree
                    .overlapping_pairs(tree)
                    .into_iter()
                    .map(|(a, b)| (ids[a], ids[b]))
                    .filter_map(|(a, b)| match (self.edge_operand[a], self.edge_operand[b]) {
                        (Some(x), Some(y)) if x < y => Some((a, b)),
                        (Some(x), Some(y)) if x > y => Some((b, a)),
                        _ => None,
                    })
                    .collect();
                found.sort_unstable();
                found.dedup();
                found
            }
        };
        pairs.sort_unstable();
        pairs
    }
}

impl Graph {
    /// The indices of `operand`'s faces, ascending.
    pub fn faces_of(&self, operand: usize) -> &[usize] {
        self.operand_faces[operand].get_or_init(|| (0..self.faces.len()).filter(|&i| self.faces[i].operand == operand).collect())
    }

    /// The indices of `operand`'s faces on a solid of it, ascending.
    pub fn solid_faces_of(&self, operand: usize) -> &[usize] {
        self.operand_solid_faces[operand].get_or_init(|| self.faces_of(operand).iter().copied().filter(|&i| self.faces[i].solid).collect())
    }
}

impl Graph {
    /// The trees over the operands' faces and edges (built on first use, once both operands
    /// are added).
    pub fn index(&self) -> &Index {
        self.index.get_or_init(|| {
            let n = self.operands();
            // Each operand's faces, edges and far faces, in increasing order.
            let mut face_ids: Vec<Vec<usize>> = vec![Vec::new(); n];
            for (i, f) in self.faces.iter().enumerate() {
                face_ids[f.operand].push(i);
            }
            let mut edge_ids: Vec<Vec<usize>> = vec![Vec::new(); n];
            for (i, e) in self.edges.iter().enumerate() {
                if let Some(k) = e.operand {
                    edge_ids[k].push(i);
                }
            }
            let mut far_ids: Vec<Vec<usize>> = vec![Vec::new(); n];
            for (i, f) in self.far.iter().enumerate() {
                far_ids[f.operand].push(i);
            }
            let tree = |ids: &[usize], bbox: &dyn Fn(usize) -> (Point, Point)| Bvh::new(&ids.iter().map(|&i| to_box(&bbox(i))).collect::<Vec<_>>());
            let face_box = |i: usize| self.faces[i].bbox;
            let edge_box = |i: usize| self.edges[i].bbox;
            let far_box = |i: usize| self.far[i].placed(self.far[i].data.sample.bbox);
            let (all_faces, all_edges) = if n > 2 {
                let faces: Vec<usize> = (0..self.faces.len()).collect();
                let edges: Vec<usize> = (0..self.edges.len()).filter(|&i| self.edges[i].operand.is_some()).collect();
                (Some((tree(&faces, &face_box), faces)), Some((tree(&edges, &edge_box), edges)))
            } else {
                (None, None)
            };
            Index {
                faces: face_ids.iter().map(|ids| tree(ids, &face_box)).collect(),
                edges: edge_ids.iter().map(|ids| tree(ids, &edge_box)).collect(),
                far: far_ids.iter().map(|ids| tree(ids, &far_box)).collect(),
                face_ids,
                edge_ids,
                far_ids,
                all_faces,
                all_edges,
                face_operand: self.faces.iter().map(|f| f.operand).collect(),
                edge_operand: self.edges.iter().map(|e| e.operand).collect(),
            }
        })
    }
}

/// Shapes met while exploring: their nodes and placements.
#[derive(Default)]
struct Seen(crate::hashing::FastMap<NodeId, Vec<Transform>>);

impl Seen {
    fn has(&self, s: &Shape) -> bool {
        self.0.get(&s.node).is_some_and(|ts| ts.contains(&s.placement))
    }

    fn add(&mut self, s: &Shape) {
        let ts = self.0.entry(s.node).or_default();
        if !ts.contains(&s.placement) {
            ts.push(s.placement);
        }
    }
}

impl Graph {
    /// Adds an operand (`operand` 0 for the object, 1 for the tool): the faces of its solids,
    /// its other faces (grouped by shell), its edges not on faces (grouped by wire) and its
    /// vertices not on edges.
    /// Adds an operand's faces and edges; faces in `far` (unplaced) are left out, as far
    /// faces.
    pub fn add_operand(&mut self, model: &Model, shape: Shape, operand: usize, far: &super::far::FarData) -> Result<()> {
        // Edge uses by the operand's faces on solids (graph edges by the model's edge where
        // they are one), to tell an open shell.
        #[derive(PartialEq, Eq, Hash)]
        enum Use {
            Model(EdgeId),
            Graph(usize),
        }
        let mut far_uses: Vec<Use> = Vec::new();
        let mut on_solids = Seen::default();
        let mut shell_of: HashMap<NodeId, Vec<(Transform, usize)>> = HashMap::default();
        for solid in model.explore(shape, NodeKind::Solid) {
            for face in model.explore(solid, NodeKind::Face) {
                on_solids.add(&face);
            }
        }
        for shell in model.explore(shape, NodeKind::Shell) {
            let faces: Vec<Shape> = model.explore(shell, NodeKind::Face).collect();
            if faces.iter().all(|f| on_solids.has(f)) {
                continue;
            }
            let container = self.containers;
            self.containers += 1;
            for face in faces {
                shell_of.entry(face.node).or_default().push((face.placement, container));
            }
        }
        let mut known: HashMap<NodeId, Vec<(Transform, usize)>> = HashMap::default();
        let mut on_faces = Seen::default();
        for face in model.explore(shape, NodeKind::Face) {
            // (A far face is unplaced, or placed where the operation only counts it.)
            if let (Some(data), Ok(id)) = (far.get(&face.node).cloned(), FaceId::try_from(face.node)) {
                let outward = if face.orientation == Orientation::Reversed { -1.0 } else { 1.0 };
                far_uses.extend(data.edges.iter().filter(|e| e.curve.is_some()).map(|e| Use::Model(e.edge)));
                let far_face = super::far::FarFace { operand, face: id, placement: face.placement, outward, data, seq: self.seq };
                let extent = far_face.placed(far_face.data.sample.extent);
                self.face_bounds[operand].push(BoundingBox::new(extent.0, extent.1));
                self.far.push(far_face);
                self.seq += 1;
                continue;
            }
            for edge in model.explore(face, NodeKind::Edge) {
                on_faces.add(&edge);
            }
            let solid = on_solids.has(&face);
            let container = shell_of.get(&face.node).and_then(|l| l.iter().find(|(t, _)| *t == face.placement)).map(|&(_, c)| c);
            let original = model
                .surface(face)
                .ok_or_else(|| Error::InvalidInput("a face without a surface".into()))?
                .transformed(&face.placement);
            // A trimmed surface is its basis (with the same parameters; the face's boundary
            // trims it), an offset surface its equivalent where it has one.
            let original = untrimmed(original);
            // A planar B-spline face (all poles in a plane) is a plane, with its pcurves
            // recomputed from its edges.
            let recognized = planar_bspline(&original);
            let output = recognized.clone().unwrap_or(original);
            // A torus is computed as the surface of revolution of its meridian circle, with
            // the same parametrization.
            let surface = match &output {
                Surface::Torus(t) => torus_as_revolution(t),
                s => s.clone(),
            };
            if !matches!(surface, Surface::Plane(_) | Surface::Cylinder(_) | Surface::Cone(_) | Surface::Sphere(_) | Surface::Revolution(_) | Surface::Extrusion(_))
                && !geom::freeform(&surface)
            {
                geom::unsupported(&surface);
            }
            let outward = if face.orientation == Orientation::Reversed { -1.0 } else { 1.0 };
            let mut boundary = Vec::new();
            let mut drift: f64 = 0.0;
            for wire in model.children(Shape { orientation: Orientation::Forward, ..face }) {
                for coedge in model.children(wire) {
                    let coedge_id = CoEdgeId::try_from(coedge.node).expect("a wire holds coedges");
                    let edge_id: EdgeId = model[coedge_id].edge;
                    // (An edge without a pcurve on the face: its 3D curve projected onto the
                    // surface.)
                    let (mut pcurve, first, last) = match model.effective_pcurve(coedge_id) {
                        Some(found) => found,
                        None => {
                            let edge_shape = Shape { placement: coedge.placement, ..Shape::from(edge_id) };
                            let (curve, a, b) = model.edge_curve(edge_shape).ok_or_else(|| Error::InvalidInput("an edge of a face without a curve or a pcurve".into()))?;
                            let curve = curve.transformed(&face.placement);
                            let projected = crate::geometry::project::analytic_pcurve(&curve, &surface)
                                .or_else(|| crate::geometry::project::pcurve_on_range_within(&curve, &surface, a, b, &[1e-9, 1e-7, 1e-5]).map(Curve2::BSpline))
                                .ok_or_else(|| Error::InvalidInput("no pcurve for an edge of a face".into()))?;
                            (projected, a, b)
                        }
                    };
                    if recognized.is_some() {
                        let edge_shape = Shape { placement: coedge.placement, ..Shape::from(edge_id) };
                        let (curve, a, b) = model.edge_curve(edge_shape).expect("an edge of a planar face has a curve");
                        pcurve = crate::geometry::project::analytic_pcurve(&curve, &surface)
                            .or_else(|| crate::geometry::project::pcurve_on_range(&curve, &surface, a, b, 1e-9).map(Curve2::BSpline))
                            .ok_or_else(|| Error::InvalidInput("no pcurve on a planar B-spline face".into()))?;
                    }
                    let placement = coedge.placement;
                    let entry = known.entry(edge_id.into()).or_default();
                    let index = match entry.iter().find(|(t, _)| *t == placement) {
                        Some(&(_, i)) => i,
                        None => {
                            let i = self.add_model_edge(model, edge_id, placement, operand, (first, last));
                            entry.push((placement, i));
                            i
                        }
                    };
                    let bounds = !matches!(model[coedge_id].orientation, Orientation::Internal | Orientation::External);
                    let along = model[coedge_id].orientation.compose(wire.orientation) == Orientation::Forward;
                    boundary.push(BoundaryUse { edge: index, pcurve, bounds, along });
                    drift = drift.max(model.tolerance(Shape::from(edge_id)));
                }
            }
            let origin = FaceId::try_from(face.node).ok().filter(|_| recognized.is_none() && face.placement == Transform::IDENTITY);
            // (Samples kept with the model's face also where it is only moved.)
            let by = face.placement.translation_part();
            let moved_only = Transform::translation(by) == face.placement;
            let cached = FaceId::try_from(face.node).ok().filter(|_| recognized.is_none() && moved_only);
            self.add_face(operand, surface, output, outward, boundary, cached.map(|id| (model, id, by)));
            let added = self.faces.last_mut().expect("face");
            added.drift = drift;
            (added.solid, added.container) = (solid, if solid { None } else { container });
            added.origin = origin;
            // (Its regions, in its own parameters, are the same wherever it is placed.)
            added.regions_key = FaceId::try_from(face.node).ok().filter(|_| recognized.is_none());
            added.seq = self.seq;
            self.seq += 1;
            let extent = BoundingBox::new(added.extent.0, added.extent.1);
            self.face_bounds[operand].push(extent);
        }
        // Edges not on faces (of wires, then alone) and vertices not on edges: none in
        // solids.
        if !made_of_solids(model, shape) {
            // Edges not on faces: of wires, then alone.
            let mut in_wires = Seen::default();
            for wire in model.explore(shape, NodeKind::Wire) {
                let edges: Vec<Shape> = model.explore(wire, NodeKind::Edge).filter(|e| !on_faces.has(e)).collect();
                if edges.is_empty() {
                    continue;
                }
                let container = self.containers;
                self.containers += 1;
                for edge in edges {
                    in_wires.add(&edge);
                    self.add_free_edge(model, edge, operand, Some(container), &mut known);
                }
            }
            for edge in model.explore(shape, NodeKind::Edge) {
                if !on_faces.has(&edge) && !in_wires.has(&edge) {
                    self.add_free_edge(model, edge, operand, None, &mut known);
                }
            }
            let mut on_edges = Seen::default();
            for edge in model.explore(shape, NodeKind::Edge) {
                for v in model.explore(edge, NodeKind::Vertex) {
                    on_edges.add(&v);
                }
            }
            for v in model.explore(shape, NodeKind::Vertex) {
                if !on_edges.has(&v) {
                    on_edges.add(&v);
                    let index = self.vertices.add(model.point(v));
                    self.free_vertices.push((operand, index));
                }
            }
        }
        let has = |k| model.explore(shape, k).next().is_some();
        self.dimension[operand] = if has(NodeKind::Solid) {
            3
        } else if has(NodeKind::Face) {
            2
        } else if has(NodeKind::Edge) {
            1
        } else if has(NodeKind::Vertex) {
            0
        } else {
            return Err(Error::InvalidInput("a boolean operation on an empty shape".into()));
        };
        // An edge used by one face only (not a seam) leaves a solid's shell open: the
        // operand is a half-space, bounded by its faces' surfaces.
        let mut uses: HashMap<Use, usize> = HashMap::default();
        for face in self.faces.iter().filter(|f| f.operand == operand && f.solid) {
            for b in &face.boundary {
                let e = &self.edges[b.edge];
                if e.curve.is_some() {
                    *uses.entry(e.source.map_or(Use::Graph(b.edge), Use::Model)).or_default() += 1;
                }
            }
        }
        for u in far_uses {
            *uses.entry(u).or_default() += 1;
        }
        self.open[operand] = uses.values().any(|&n| n < 2);
        Ok(())
    }

    /// An edge of a wire or edge argument.
    fn add_free_edge(&mut self, model: &Model, edge: Shape, operand: usize, container: Option<usize>, known: &mut HashMap<NodeId, Vec<(Transform, usize)>>) {
        let id = EdgeId::try_from(edge.node).expect("an edge");
        if model.edge_curve(edge).is_none() {
            return; // degenerate: nothing of it remains in a result
        }
        let entry = known.entry(edge.node).or_default();
        if entry.iter().any(|(t, _)| *t == edge.placement) {
            return;
        }
        let i = self.add_model_edge(model, id, edge.placement, operand, (0.0, 0.0));
        entry.push((edge.placement, i));
        self.edges[i].free = true;
        self.edges[i].container = container;
        self.edges[i].reversed = edge.orientation == Orientation::Reversed;
    }

    fn add_model_edge(&mut self, model: &Model, id: EdgeId, placement: Transform, operand: usize, range: (f64, f64)) -> usize {
        let shape = Shape { placement, ..Shape::from(id) };
        let (start_vertex, end_vertex) = model.edge_vertices(shape);
        let (start, end) = (self.vertices.add(model.point(start_vertex)), self.vertices.add(model.point(end_vertex)));
        let (curve, first, last) = match model.edge_curve(shape) {
            Some((c, a, b)) => (Some(c), a, b),
            None => (None, range.0, range.1),
        };
        let start_point = self.vertices.points[start];
        // (Kept with the model's edge, unplaced, also where it is only moved.)
        let by = placement.translation_part();
        let bbox = if placement == Transform::IDENTITY {
            model_edge_box(model, id, curve.as_ref(), first, last, start_point)
        } else if Transform::translation(by) == placement {
            let unplaced = |model: &Model| model.edge_curve(Shape::from(id)).map(|(c, _, _)| c);
            let (lo, hi) = model.derived(id.into(), model.edge_fingerprint(id), || EdgeBox(Edge::sampled_box(unplaced(model).as_ref(), first, last, start_point - by))).0;
            (lo + by, hi + by)
        } else {
            Edge::sampled_box(curve.as_ref(), first, last, start_point)
        };
        self.edges.push(Edge::with_box(curve, first, last, start, end, bbox, Some(operand)));
        self.edges.last_mut().expect("edge").source = Some(id);
        self.edges.last_mut().expect("edge").tolerance = model.tolerance(Shape::from(id));
        if placement == Transform::IDENTITY {
            self.edges.last_mut().expect("edge").origin = Some(id);
            for (v, shape) in [(start, start_vertex), (end, end_vertex)] {
                if let (Ok(vertex), true) = (VertexId::try_from(shape.node), shape.placement == Transform::IDENTITY) {
                    self.vertex_origin.entry(v).or_insert(vertex);
                }
            }
        }
        self.edges.len() - 1
    }

    /// Adds a face; with `cached` (the model's face, moved by a translation), its samples
    /// (boundary polylines, ranges, boxes) are kept in the model for the model's face, as
    /// unplaced, while that face is unchanged, then moved.
    fn add_face(&mut self, operand: usize, surface: Surface, output: Surface, outward: f64, boundary: Vec<BoundaryUse>, cached: Option<(&Model, FaceId, Vector)>) {
        let sample = match cached {
            Some((model, id, by)) => {
                let kept = model.derived(id.into(), model.face_fingerprint(id), || match unplaced_surface(model, id) {
                    Some(local) => self.sample_face(&local, &boundary),
                    None => self.sample_face(&surface.transformed(&Transform::translation(-by)), &boundary),
                });
                let moved = |(lo, hi): (Point, Point)| if by == Vector::ZERO { (lo, hi) } else { (lo + by, hi + by) };
                FaceSample { polylines: kept.polylines.clone(), u_range: kept.u_range, v_range: kept.v_range, bbox: moved(kept.bbox), extent: moved(kept.extent) }
            }
            None => self.sample_face(&surface, &boundary),
        };
        let FaceSample { polylines, u_range, v_range, bbox, extent } = sample;
        let (u_period, v_period) = (surface.u_period(), surface.v_period());
        let closed_bounds = {
            use crate::geometry::ParametricSurface;
            let (u0, u1, v0, v1) = surface.bounds();
            let finite = |a: f64, b: f64| a.is_finite() && b.is_finite();
            [
                (u_period.is_none() && finite(u0, u1) && surface.is_u_closed()).then_some((u0, u1)),
                (v_period.is_none() && finite(v0, v1) && surface.is_v_closed()).then_some((v0, v1)),
            ]
        };
        self.faces.push(Face {
            operand,
            solid: true,
            container: None,
            surface,
            output,
            outward,
            boundary,
            polylines,
            u_period,
            v_period,
            u_range,
            v_range,
            closed_bounds,
            bbox,
            extent,
            vertices_on: Vec::new(),
            sections: Vec::new(),
            on_edges: Vec::new(),
            origin: None,
            regions_key: None,
            seq: 0,
            drift: 0.0,
            samples: Default::default(),
            inside_index: Default::default(),
            cells: Default::default(),
        });
    }

    /// A face's boundary sampled in (u, v), its parameter ranges and boxes.
    fn sample_face(&self, surface: &Surface, boundary: &[BoundaryUse]) -> FaceSample {
        let ranged: Vec<(&Curve2, f64, f64)> = boundary.iter().filter(|b| b.bounds).map(|b| (&b.pcurve, self.edges[b.edge].first, self.edges[b.edge].last)).collect();
        sample_face(surface, &ranged)
    }
}

/// A face's boundary sampled in (u, v) (each pcurve over its range), its parameter ranges
/// and boxes.
pub(super) fn sample_face(surface: &Surface, boundary: &[(&Curve2, f64, f64)]) -> FaceSample {
    {
        let polylines: Vec<Vec<Uv>> = boundary.iter().map(|&(pcurve, first, last)| geom::sample_pcurve_closely(pcurve, first, last)).collect();
        let us = polylines.iter().flatten().map(|p| p[0]);
        let u_range = us.clone().fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), u| (lo.min(u), hi.max(u)));
        let vs = polylines.iter().flatten().map(|p| p[1]);
        let v_range = vs.fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), v| (lo.min(v), hi.max(v)));
        if polylines.is_empty() {
            // Without a boundary: the whole (unbounded) surface.
            let far = crate::tolerance::INFINITE;
            let bbox = (Point::new(-far, -far, -far), Point::new(far, far, far));
            return FaceSample { polylines, u_range: (-far, far), v_range: (-far, far), bbox, extent: bbox };
        }
        // The box of a grid over the face's parameter rectangle (a curved face may bulge
        // beyond its boundary).
        let grid = |n: (usize, usize)| -> Vec<Point> {
            let us: Vec<f64> = (0..=n.0).map(|i| u_range.0 + (u_range.1 - u_range.0) * i as f64 / n.0 as f64).collect();
            let vs: Vec<f64> = (0..=n.1).map(|j| v_range.0 + (v_range.1 - v_range.0) * j as f64 / n.1 as f64).collect();
            // (Evaluated row by row along u, as a B-spline batches them; kept u-major.)
            let rows = surface.points_on_grid(&us, &vs);
            (0..=n.0).flat_map(|i| (0..=n.1).map(move |j| (i, j))).map(|(i, j)| rows[j][i]).collect()
        };
        let mut n = (12, 12);
        let mut points = grid(n);
        // (On a B-spline turning far between those points, as a thread's flank winding many
        // times round, whose points a turn apart are alike: the box missing much of its
        // poles' box, which holds it. A few points in each of its spans across the face then.)
        if let Surface::BSpline(b) = surface {
            let (lo, hi) = geom::bounding_box(points.iter().copied(), 0.0);
            let (plo, phi) = geom::bounding_box(b.poles().iter().flatten().copied(), 0.0);
            let reach = 0.25 * plo.distance(phi);
            if plo.x < lo.x - reach || plo.y < lo.y - reach || plo.z < lo.z - reach || phi.x > hi.x + reach || phi.y > hi.y + reach || phi.z > hi.z + reach {
                let spans = |knots: &[f64], (lo, hi): (f64, f64)| (2 * knots.iter().filter(|&&k| k > lo && k < hi).count()).clamp(12, 512);
                n = (spans(b.u_knots(), u_range), spans(b.v_knots(), v_range));
                points = grid(n);
            }
        }
        let diagonal = geom::bounding_box(points.iter().copied(), 0.0);
        // (Kept in the model: wide enough for the looser tolerances too.)
        let margin = 0.05 * diagonal.0.distance(diagonal.1) + geom::LOOSE_ON_TOL;
        // Between grid points a curved surface strays from them by a small part of the cell
        // (an eighth of the chord times the angle turned): a quarter of the longest cell
        // side covers it.
        let bulge = if matches!(surface, Surface::Plane(_)) {
            0.0
        } else {
            let at = |i: usize, j: usize| points[i * (n.1 + 1) + j];
            (0..n.0).flat_map(|i| (0..n.1).map(move |j| (i, j))).map(|(i, j)| at(i, j).distance(at(i + 1, j)).max(at(i, j).distance(at(i, j + 1)))).fold(0.0, f64::max) * 0.25
        };
        let extent = geom::bounding_box(points.iter().copied(), bulge + geom::LOOSE_ON_TOL);
        let bbox = geom::bounding_box(points, margin);
        FaceSample { polylines, u_range, v_range, bbox, extent }
    }
}

/// The surface the graph computes with for the model's face `id`, unplaced (as
/// [`Graph::add_operand`] makes it of the placed one).
fn unplaced_surface(model: &Model, id: FaceId) -> Option<Surface> {
    Some(match untrimmed(model.surface(Shape::from(id))?.clone()) {
        Surface::Torus(t) => torus_as_revolution(&t),
        s => s,
    })
}

/// `surface` without the trimming of a trimmed surface (its basis has the same
/// parameters), and an offset surface as its equivalent surface where it has one.
fn untrimmed(surface: Surface) -> Surface {
    match surface {
        Surface::Trimmed(t) => untrimmed(t.basis().clone()),
        Surface::Offset(o) => match o.equivalent_surface() {
            Some(e) => untrimmed(e),
            None => Surface::Offset(o),
        },
        s => s,
    }
}
