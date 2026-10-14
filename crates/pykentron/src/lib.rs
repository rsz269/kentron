//! The native part of `pykentron`: kentron's modeling operations on shape handles, for the
//! Python package (which provides the CadQuery-compatible API on top).
//!
//! Every shape lives in one model shared by the Python process (kentron's arena); a
//! [`Handle`] names a shape in it, with its orientation and placement. What no handle
//! reaches any more is collected from time to time ([`with_model`]). Points, vectors and
//! placements cross the boundary as tuples (placements as the three rows of a 3×4 matrix).
//! kentron's errors become `ValueError`s, and operations it does not implement yet (its
//! `todo!()`s) `NotImplementedError`s.

use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Mutex, MutexGuard, OnceLock};

use pyo3::exceptions::{PyNotImplementedError, PyRuntimeError, PyValueError};
use pyo3::prelude::*;

use kentron::geometry::curve::{Circle, Curve};
use kentron::geometry::surface::Surface;
use kentron::geometry::{ParametricCurve, ParametricSurface};
use kentron::math::{Axis, Direction, Frame, Point, Transform, Vector};
use kentron::hashing::FastMap;
use kentron::topology::{CollectOptions, Model, NodeId, NodeKind, Orientation, Shape, WireId};
use kentron::{blend, boolean, loft, make, modify, offset, primitives, properties, sweep, validity};

mod cadquery;
mod profile;

type Xyz = (f64, f64, f64);
type Rows = [[f64; 4]; 3];

fn model() -> &'static Mutex<Model> {
    static MODEL: OnceLock<Mutex<Model>> = OnceLock::new();
    MODEL.get_or_init(|| Mutex::new(Model::new()))
}

/// How many handles name each node: what the model keeps when it collects garbage.
fn handled() -> MutexGuard<'static, FastMap<NodeId, u32>> {
    static HANDLED: OnceLock<Mutex<FastMap<NodeId, u32>>> = OnceLock::new();
    HANDLED.get_or_init(Default::default).lock().unwrap_or_else(|poisoned| poisoned.into_inner())
}

/// Nodes the model may hold before it first collects garbage.
const FIRST_COLLECTION: usize = 50_000;

/// Collects what no handle reaches once the model has twice the nodes it kept at the last
/// collection (so collecting costs a bounded share of the work that made the nodes), as the
/// [profile](profile) says.
fn collect_if_grown(model: &mut Model) {
    static NEXT: AtomicUsize = AtomicUsize::new(FIRST_COLLECTION);
    let nodes = |model: &Model| NodeKind::ALL.iter().map(|&k| model.count(k)).sum::<usize>();
    if nodes(model) < NEXT.load(Ordering::Relaxed) {
        return;
    }
    let roots: Vec<NodeId> = handled().keys().copied().collect();
    model.collect_garbage_with(roots, CollectOptions { forget_derived: profile::current().forgets_derived() });
    NEXT.store((2 * nodes(model)).max(FIRST_COLLECTION), Ordering::Relaxed);
}

/// Runs `f` on the shared model, turning kentron's errors and unimplemented operations into
/// Python exceptions.
///
/// Garbage is collected first, if it is time: every shape `f` is given comes from a live
/// handle, and so is kept. (A binding must not carry a shape, rather than a handle, from
/// one call to the next.)
fn with_model<R>(f: impl FnOnce(&mut Model) -> kentron::Result<R>) -> PyResult<R> {
    let mut guard = model().lock().unwrap_or_else(|poisoned| poisoned.into_inner());
    let outcome = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        collect_if_grown(&mut guard);
        f(&mut guard)
    }));
    match outcome {
        Ok(Ok(r)) => Ok(r),
        Ok(Err(e)) => Err(PyValueError::new_err(e.to_string())),
        Err(panic) => {
            let message = panic
                .downcast_ref::<String>()
                .cloned()
                .or_else(|| panic.downcast_ref::<&str>().map(|s| s.to_string()))
                .unwrap_or_else(|| "kentron panicked".into());
            if message.starts_with("not yet implemented") || message.starts_with("not implemented") {
                Err(PyNotImplementedError::new_err(message))
            } else {
                Err(PyRuntimeError::new_err(message))
            }
        }
    }
}

fn point(p: Xyz) -> Point {
    Point::new(p.0, p.1, p.2)
}

fn vector(v: Xyz) -> Vector {
    Vector::new(v.0, v.1, v.2)
}

fn xyz(p: Point) -> Xyz {
    (p.x, p.y, p.z)
}

fn direction(v: Xyz) -> PyResult<Direction> {
    Direction::try_from_vector(vector(v)).ok_or_else(|| PyValueError::new_err("a zero direction"))
}

/// The frame at `origin` with Z along `z` (and X along `x` if given, else chosen like
/// OCCT's `gp_Ax2(P, V)`).
fn frame(origin: Xyz, z: Xyz, x: Option<Xyz>) -> PyResult<Frame> {
    let z = direction(z)?;
    Ok(match x {
        Some(x) => Frame::new(point(origin), z, direction(x)?),
        None => Frame::from_z(point(origin), z),
    })
}

fn kind_name(kind: NodeKind) -> &'static str {
    match kind {
        NodeKind::Vertex => "Vertex",
        NodeKind::Edge => "Edge",
        NodeKind::Wire => "Wire",
        NodeKind::Face => "Face",
        NodeKind::Shell => "Shell",
        NodeKind::Solid => "Solid",
        NodeKind::CompSolid => "CompSolid",
        NodeKind::Compound => "Compound",
        _ => "Other",
    }
}

fn kind_of(name: &str) -> PyResult<NodeKind> {
    Ok(match name {
        "Vertex" => NodeKind::Vertex,
        "Edge" => NodeKind::Edge,
        "Wire" => NodeKind::Wire,
        "Face" => NodeKind::Face,
        "Shell" => NodeKind::Shell,
        "Solid" => NodeKind::Solid,
        "CompSolid" => NodeKind::CompSolid,
        "Compound" => NodeKind::Compound,
        other => return Err(PyValueError::new_err(format!("unknown shape kind {other}"))),
    })
}

/// A shape of the shared model: a node with its orientation and placement. The node, and
/// what it is made of, stay in the model while a handle names it.
#[pyclass(frozen, module = "pykentron._kentron")]
struct Handle {
    shape: Shape,
}

fn handle(shape: Shape) -> Handle {
    *handled().entry(shape.node).or_default() += 1;
    Handle { shape }
}

impl Clone for Handle {
    fn clone(&self) -> Self {
        handle(self.shape)
    }
}

impl Drop for Handle {
    fn drop(&mut self) {
        let mut handled = handled();
        if let Some(count) = handled.get_mut(&self.shape.node) {
            *count -= 1;
            if *count == 0 {
                handled.remove(&self.shape.node);
            }
        }
    }
}

#[pymethods]
impl Handle {
    /// The kind: "Vertex", "Edge", "Wire", "Face", "Shell", "Solid", "CompSolid" or
    /// "Compound".
    fn kind(&self) -> &'static str {
        kind_name(self.shape.kind())
    }

    /// "Forward", "Reversed", "Internal" or "External".
    fn orientation(&self) -> &'static str {
        match self.shape.orientation {
            Orientation::Forward => "Forward",
            Orientation::Reversed => "Reversed",
            Orientation::Internal => "Internal",
            Orientation::External => "External",
        }
    }

    fn reversed(&self) -> Handle {
        handle(self.shape.reversed())
    }

    fn with_orientation(&self, orientation: &str) -> PyResult<Handle> {
        let orientation = match orientation {
            "Forward" => Orientation::Forward,
            "Reversed" => Orientation::Reversed,
            "Internal" => Orientation::Internal,
            "External" => Orientation::External,
            other => return Err(PyValueError::new_err(format!("unknown orientation {other}"))),
        };
        Ok(handle(Shape { orientation, ..self.shape }))
    }

    /// The placement, as the rows of a 3×4 matrix.
    fn placement(&self) -> Rows {
        self.shape.placement.rows()
    }

    /// The shape placed by `rows` (a rigid motion or uniform scale) after its own placement.
    fn moved(&self, rows: Rows) -> PyResult<Handle> {
        let t = Transform::from_rows(rows).ok_or_else(|| PyValueError::new_err("a singular placement"))?;
        Ok(handle(self.shape.moved(&t)))
    }

    /// The shape with the placement `rows` instead of its own.
    fn located(&self, rows: Rows) -> PyResult<Handle> {
        let t = Transform::from_rows(rows).ok_or_else(|| PyValueError::new_err("a singular placement"))?;
        Ok(handle(self.shape.with_placement(t)))
    }

    /// Same node and placement (whatever the orientation).
    fn is_same(&self, other: &Handle) -> bool {
        self.shape.is_same(&other.shape)
    }

    /// Same node, placement and orientation.
    fn is_equal(&self, other: &Handle) -> bool {
        self.shape == other.shape
    }

    /// A hash of the node (equal for the same node, whatever its placement).
    fn node_hash(&self) -> u64 {
        let n = self.shape.node;
        (n.index as u64) << 4 | n.kind as u64
    }
}

// --- construction --------------------------------------------------------------------------

/// A box of sizes `dx × dy × dz` with a corner at the origin of the frame.
#[pyfunction]
#[pyo3(signature = (dx, dy, dz, origin, z=(0.0, 0.0, 1.0), x=None))]
fn make_box(dx: f64, dy: f64, dz: f64, origin: Xyz, z: Xyz, x: Option<Xyz>) -> PyResult<Handle> {
    let f = frame(origin, z, x)?;
    with_model(|m| primitives::Cuboid::new(dx, dy, dz).with_frame(f).build(m)).map(handle)
}

/// A cylinder (a sector of `angle` degrees) on the frame's Z axis.
#[pyfunction]
#[pyo3(signature = (radius, height, origin, z, angle=360.0))]
fn make_cylinder(radius: f64, height: f64, origin: Xyz, z: Xyz, angle: f64) -> PyResult<Handle> {
    let f = frame(origin, z, None)?;
    with_model(|m| {
        let mut c = primitives::Cylinder::new(radius, height).with_frame(f);
        if angle < 360.0 {
            c = c.with_angle(angle.to_radians());
        }
        c.build(m)
    })
    .map(handle)
}

/// A cone (or frustum) on the frame's Z axis.
#[pyfunction]
#[pyo3(signature = (radius1, radius2, height, origin, z, angle=360.0))]
fn make_cone(radius1: f64, radius2: f64, height: f64, origin: Xyz, z: Xyz, angle: f64) -> PyResult<Handle> {
    let f = frame(origin, z, None)?;
    with_model(|m| {
        let mut c = primitives::Cone::new(radius1, radius2, height).with_frame(f);
        if angle < 360.0 {
            c = c.with_angle(angle.to_radians());
        }
        c.build(m)
    })
    .map(handle)
}

/// A sphere (between latitudes, a sector of `angle` degrees) centred on the frame's origin.
#[pyfunction]
#[pyo3(signature = (radius, origin, z, latitude1=-90.0, latitude2=90.0, angle=360.0))]
fn make_sphere(radius: f64, origin: Xyz, z: Xyz, latitude1: f64, latitude2: f64, angle: f64) -> PyResult<Handle> {
    let f = frame(origin, z, None)?;
    with_model(|m| {
        let mut s = primitives::Sphere::new(radius).with_frame(f);
        if latitude1 > -90.0 || latitude2 < 90.0 {
            s = s.with_latitudes(latitude1.to_radians(), latitude2.to_radians());
        }
        if angle < 360.0 {
            s = s.with_angle(angle.to_radians());
        }
        s.build(m)
    })
    .map(handle)
}

/// A torus around the frame's Z axis.
#[pyfunction]
#[pyo3(signature = (major, minor, origin, z))]
fn make_torus(major: f64, minor: f64, origin: Xyz, z: Xyz) -> PyResult<Handle> {
    let f = frame(origin, z, None)?;
    with_model(|m| primitives::Torus::new(major, minor).with_frame(f).build(m)).map(handle)
}

#[pyfunction]
fn make_vertex(p: Xyz) -> PyResult<Handle> {
    with_model(|m| Ok(make::vertex(m, point(p)))).map(handle)
}

/// A straight edge.
#[pyfunction]
fn make_line(a: Xyz, b: Xyz) -> PyResult<Handle> {
    with_model(|m| make::edge_between(m, point(a), point(b))).map(handle)
}

/// A circular edge (an arc from `angle1` to `angle2` degrees, measured from the frame's X).
#[pyfunction]
#[pyo3(signature = (radius, centre, normal, angle1=0.0, angle2=360.0, x=None))]
fn make_circle(radius: f64, centre: Xyz, normal: Xyz, angle1: f64, angle2: f64, x: Option<Xyz>) -> PyResult<Handle> {
    let f = frame(centre, normal, x)?;
    with_model(|m| {
        let curve = Curve::Circle(Circle::new(f, radius));
        if angle2 - angle1 >= 360.0 - 1e-9 {
            make::edge_from_curve(m, curve)
        } else {
            make::edge_from_curve_range(m, curve, angle1.to_radians(), angle2.to_radians())
        }
    })
    .map(handle)
}

/// The arc from `a` through `b` to `c`.
#[pyfunction]
fn make_three_point_arc(a: Xyz, b: Xyz, c: Xyz) -> PyResult<Handle> {
    let (a, b, c) = (point(a), point(b), point(c));
    // The circumcentre, in the plane of the points.
    let (ab, ac) = (b - a, c - a);
    let n = ab.cross(ac);
    if n.length() <= 1e-12 * ab.length() * ac.length() {
        return Err(PyValueError::new_err("an arc through collinear points"));
    }
    let centre = a + (n.cross(ab) * ac.dot(ac) + ac.cross(n) * ab.dot(ab)) * (1.0 / (2.0 * n.dot(n)));
    let radius = centre.distance(a);
    let f = Frame::new(centre, direction((n.x, n.y, n.z))?, direction(((a - centre).x, (a - centre).y, (a - centre).z))?);
    let angle = |p: Point| {
        let d = p - centre;
        d.dot(Vector::from(f.y())).atan2(d.dot(Vector::from(f.x()))).rem_euclid(std::f64::consts::TAU)
    };
    // Counter-clockwise about the normal from a (angle 0), b before c.
    let end = angle(c);
    with_model(|m| make::edge_from_curve_range(m, Curve::Circle(Circle::new(f, radius)), 0.0, end)).map(handle)
}

/// A wire through `edges`, in order.
#[pyfunction]
#[pyo3(signature = (edges, tolerance=0.0))]
fn make_wire(edges: Vec<Handle>, tolerance: f64) -> PyResult<Handle> {
    let edges: Vec<Shape> = edges.iter().map(|h| h.shape).collect();
    with_model(|m| make::wire_within(m, &edges, tolerance)).map(handle)
}

/// Whether an edge or wire is closed topologically: an edge starting and ending at one
/// vertex, a wire whose last edge returns to its first vertex (OCCT's `BRep_Tool::IsClosed`).
#[pyfunction]
fn is_closed(shape: Handle) -> PyResult<bool> {
    with_model(|m| {
        Ok(match shape.shape.kind() {
            NodeKind::Wire => WireId::try_from(shape.shape.node).is_ok_and(|w| m.is_closed_wire(w)),
            NodeKind::Edge => {
                let (a, b) = m.edge_vertices(shape.shape);
                a.node == b.node
            }
            _ => false,
        })
    })
}

/// The polygon through `points` (closed back to the first if `close`).
#[pyfunction]
#[pyo3(signature = (points, close=false))]
fn make_polygon(points: Vec<Xyz>, close: bool) -> PyResult<Handle> {
    let mut pts: Vec<Point> = points.into_iter().map(point).collect();
    pts.dedup_by(|a, b| a.distance(*b) <= kentron::tolerance::LINEAR);
    if close && pts.len() > 1 && pts[0].distance(pts[pts.len() - 1]) <= kentron::tolerance::LINEAR {
        pts.pop();
    }
    with_model(|m| {
        let n = if close { pts.len() } else { pts.len() - 1 };
        let edges = (0..n).map(|i| make::edge_between(m, pts[i], pts[(i + 1) % pts.len()])).collect::<kentron::Result<Vec<_>>>()?;
        make::wire(m, &edges)
    })
    .map(handle)
}

/// The planar face inside `outer`, with holes inside `inners`.
#[pyfunction]
#[pyo3(signature = (outer, inners=Vec::new()))]
fn make_face(outer: Handle, inners: Vec<Handle>) -> PyResult<Handle> {
    with_model(|m| {
        // Fresh wires through the given ones' edges (a wire bounds one face only).
        let rewire = |m: &mut Model, w: Shape| {
            let edges: Vec<Shape> = m.explore(w, NodeKind::Edge).collect();
            make::wire(m, &edges)
        };
        let outer = rewire(m, outer.shape)?;
        // An open outer wire makes an invalid face, as in OCCT (the caller checks validity).
        let open = WireId::try_from(outer.node).is_ok_and(|w| !m.is_closed_wire(w));
        if open && inners.is_empty() {
            return make::open_face(m, outer);
        }
        let face = make::planar_face(m, outer)?;
        let normal = face_normal(m, face);
        for inner in inners {
            let wire = rewire(m, inner.shape)?;
            // Holes run clockwise about the face's normal.
            let probe = make::planar_face(m, wire)?;
            let along = face_normal(m, probe).dot(normal) > 0.0;
            let wire = rewire(m, inner.shape)?;
            m.add(face, if along { wire.reversed() } else { wire })?;
        }
        Ok(face)
    })
    .map(handle)
}

/// The normal of a planar face as it is used.
fn face_normal(m: &Model, face: Shape) -> Vector {
    match m.surface(face) {
        Some(Surface::Plane(p)) => {
            let n = Vector::from(p.normal());
            if face.orientation == Orientation::Reversed { -n } else { n }
        }
        _ => Vector::ZERO,
    }
}

/// `profile` (a face or wire) swept along `v`.
#[pyfunction]
fn extrude(profile: Handle, v: Xyz) -> PyResult<Handle> {
    with_model(|m| sweep::extrude(m, profile.shape, vector(v)).map(|s| s.shape)).map(handle)
}

/// `profile` revolved by `angle` degrees about the axis through `origin` along `axis`.
#[pyfunction]
fn revolve(profile: Handle, origin: Xyz, axis: Xyz, angle: f64) -> PyResult<Handle> {
    let axis = Axis::new(point(origin), direction(axis)?);
    with_model(|m| sweep::revolve(m, profile.shape, axis, angle.to_radians()).map(|s| s.shape)).map(handle)
}

/// The loft through `sections` (wires, or vertices at the ends).
#[pyfunction]
#[pyo3(signature = (sections, solid=true, ruled=false))]
fn make_loft(sections: Vec<Handle>, solid: bool, ruled: bool) -> PyResult<Handle> {
    let sections: Vec<Shape> = sections.iter().map(|h| h.shape).collect();
    with_model(|m| loft::loft(m, &sections, &loft::LoftOptions { solid, ruled, ..Default::default() }).map(|l| l.shape)).map(handle)
}

/// A planar wire offset in its plane ("arc" or "intersection" joins).
#[pyfunction]
#[pyo3(signature = (wire, distance, join="arc"))]
fn offset_wire(wire: Handle, distance: f64, join: &str) -> PyResult<Handle> {
    let join = join_of(join)?;
    with_model(|m| offset::planar_wire(m, wire.shape, distance, join)).map(handle)
}

/// A solid or shell offset along its faces' normals.
#[pyfunction]
#[pyo3(signature = (shape, distance, join="arc"))]
fn offset_shape(shape: Handle, distance: f64, join: &str) -> PyResult<Handle> {
    let join = join_of(join)?;
    with_model(|m| offset::shape(m, shape.shape, distance, join)).map(handle)
}

/// The smallest distance between two shapes, and the closest points.
#[pyfunction]
fn distance(a: Handle, b: Handle) -> PyResult<(f64, Xyz, Xyz)> {
    let d = with_model(|m| kentron::distance::distance(m, a.shape, b.shape))?;
    Ok((d.value, (d.on_first.x, d.on_first.y, d.on_first.z), (d.on_second.x, d.on_second.y, d.on_second.z)))
}

/// Where a point lies with respect to a shape's solids: "in", "out" or "on".
#[pyfunction]
fn classify_point(shape: Handle, point: Xyz) -> PyResult<&'static str> {
    let p = kentron::math::Point::new(point.0, point.1, point.2);
    let at = with_model(|m| kentron::boolean::classify_point(m, shape.shape, p))?;
    Ok(match at {
        kentron::boolean::PointLocation::Inside => "in",
        kentron::boolean::PointLocation::Outside => "out",
        kentron::boolean::PointLocation::Boundary => "on",
    })
}

/// The shape's faces as one triangle mesh: points and index triples.
#[pyfunction]
#[pyo3(signature = (shape, tolerance, angular=0.5))]
fn triangulate(shape: Handle, tolerance: f64, angular: f64) -> PyResult<(Vec<Xyz>, Vec<(u32, u32, u32)>)> {
    let options = kentron::meshing::MeshOptions { deflection: tolerance, angle: angular };
    let (nodes, triangles) = with_model(|m| kentron::meshing::triangulate(m, shape.shape, options))?;
    Ok((nodes.iter().map(|p| (p.x, p.y, p.z)).collect(), triangles.iter().map(|t| (t[0], t[1], t[2])).collect()))
}

/// A shape's mesh kept on the Rust side (for analyses over it, without a Python object per
/// node).
#[pyclass(frozen, module = "pykentron._kentron")]
struct Mesh {
    nodes: Vec<kentron::math::Point>,
    triangles: Vec<[u32; 3]>,
}

#[pymethods]
impl Mesh {
    /// How many nodes and triangles.
    fn counts(&self) -> (usize, usize) {
        (self.nodes.len(), self.triangles.len())
    }

    /// The nodes' coordinates (f64, x y z per node) and the triangles' node indices (u32,
    /// three per triangle), as native-endian bytes.
    fn arrays<'py>(&self, py: Python<'py>) -> (Bound<'py, pyo3::types::PyBytes>, Bound<'py, pyo3::types::PyBytes>) {
        let nodes: Vec<u8> = self.nodes.iter().flat_map(|p| [p.x, p.y, p.z]).flat_map(f64::to_ne_bytes).collect();
        let triangles: Vec<u8> = self.triangles.iter().flatten().flat_map(|i| i.to_ne_bytes()).collect();
        (pyo3::types::PyBytes::new(py, &nodes), pyo3::types::PyBytes::new(py, &triangles))
    }
}

/// The shape's mesh, as `triangulate` makes it, kept on the Rust side.
#[pyfunction]
fn mesh(shape: Handle, tolerance: f64, angular: f64) -> PyResult<Mesh> {
    let options = kentron::meshing::MeshOptions { deflection: tolerance, angle: angular };
    let (nodes, triangles) = with_model(|m| kentron::meshing::triangulate(m, shape.shape, options))?;
    Ok(Mesh { nodes, triangles })
}

/// Which points of the grid `xs` x `ys` x `zs` lie inside the mesh (by its crossings along
/// +Z): one byte (0 or 1) per point, x-major.
#[pyfunction]
fn voxelize<'py>(py: Python<'py>, mesh: &Mesh, xs: Vec<f64>, ys: Vec<f64>, zs: Vec<f64>) -> Bound<'py, pyo3::types::PyBytes> {
    let inside = kentron::mesh::voxelize(&mesh.nodes, &mesh.triangles, &xs, &ys, &zs);
    let bytes: Vec<u8> = inside.into_iter().map(u8::from).collect();
    pyo3::types::PyBytes::new(py, &bytes)
}

/// A voxel grid's morphology: `operation` ("erode", "dilate", "open" or "close") of the grid
/// `cells` (one byte per cell, x-major, of sizes `dims`) by the structuring element
/// `element` (likewise, of odd sizes `element_dims`, centred), cells outside empty; one
/// byte (0 or 1) per cell.
#[pyfunction]
fn morphology<'py>(py: Python<'py>, operation: &str, cells: &[u8], dims: [usize; 3], element: &[u8], element_dims: [usize; 3]) -> PyResult<Bound<'py, pyo3::types::PyByteArray>> {
    use kentron::voxels;
    if cells.len() != dims.iter().product::<usize>() {
        return Err(PyValueError::new_err("the cells do not fill the grid's sizes"));
    }
    let to_bool = |bytes: &[u8]| bytes.iter().map(|&b| b != 0).collect::<Vec<bool>>();
    let shape = voxels::Element::new(&to_bool(element), element_dims).ok_or_else(|| PyValueError::new_err("a structuring element of odd sizes, filled"))?;
    let grid = to_bool(cells);
    let out = match operation {
        "erode" => voxels::erode(&grid, dims, &shape),
        "dilate" => voxels::dilate(&grid, dims, &shape),
        "open" => voxels::open(&grid, dims, &shape),
        "close" => voxels::close(&grid, dims, &shape),
        other => return Err(PyValueError::new_err(format!("no morphology operation {other:?}"))),
    };
    pyo3::types::PyByteArray::new_with(py, out.len(), |bytes| {
        for (b, &v) in bytes.iter_mut().zip(&out) {
            *b = u8::from(v);
        }
        Ok(())
    })
}

/// A voxel grid's connected parts (cells sharing a face), in the order of their first cells:
/// their numbers of cells and their centres (mean cell indices).
#[pyfunction]
fn components(cells: &[u8], dims: [usize; 3]) -> PyResult<(Vec<usize>, Vec<Xyz>)> {
    if cells.len() != dims.iter().product::<usize>() {
        return Err(PyValueError::new_err("the cells do not fill the grid's sizes"));
    }
    let grid: Vec<bool> = cells.iter().map(|&b| b != 0).collect();
    let found = kentron::voxels::components(&grid, dims);
    Ok(found.into_iter().map(|(n, c)| (n, (c[0], c[1], c[2]))).unzip())
}

/// Writes the shape's mesh to an STL file.
#[pyfunction]
#[pyo3(signature = (shape, path, tolerance, angular=0.5, ascii=false))]
fn write_stl(shape: Handle, path: &str, tolerance: f64, angular: f64, ascii: bool) -> PyResult<()> {
    let options = kentron::meshing::MeshOptions { deflection: tolerance, angle: angular };
    let (nodes, triangles) = with_model(|m| kentron::meshing::triangulate(m, shape.shape, options))?;
    let mut f = std::fs::File::create(path).map_err(|e| PyValueError::new_err(format!("{path}: {e}")))?;
    kentron::io::stl::write(&nodes, &triangles, &mut f, ascii).map_err(|e| PyValueError::new_err(e.to_string()))
}

/// The shape with faces on the same surfaces, and edges on the same curves, merged.
#[pyfunction]
fn unify_same_domain(shape: Handle) -> PyResult<Handle> {
    with_model(|m| kentron::unify::unify_same_domain(m, shape.shape)).map(handle)
}

/// A solid hollowed out to walls of `thickness`, open at `faces`.
#[pyfunction]
#[pyo3(signature = (shape, faces, thickness, join="arc"))]
fn thick_solid(shape: Handle, faces: Vec<Handle>, thickness: f64, join: &str) -> PyResult<Handle> {
    let join = join_of(join)?;
    let faces: Vec<Shape> = faces.iter().map(|h| h.shape).collect();
    with_model(|m| offset::thick_solid(m, shape.shape, &faces, thickness, join)).map(handle)
}

fn join_of(name: &str) -> PyResult<offset::Join> {
    match name {
        "arc" => Ok(offset::Join::Arc),
        "intersection" => Ok(offset::Join::Intersection),
        other => Err(PyValueError::new_err(format!("unknown join {other}"))),
    }
}

/// `shape` with `edges` rounded by `radius`.
#[pyfunction]
fn fillet(shape: Handle, edges: Vec<Handle>, radius: f64) -> PyResult<Handle> {
    let edges: Vec<Shape> = edges.iter().map(|h| h.shape).collect();
    with_model(|m| blend::fillet(m, shape.shape, &edges, radius)).map(handle)
}

/// A compound of `shapes`.
#[pyfunction]
fn make_compound(shapes: Vec<Handle>) -> PyResult<Handle> {
    with_model(|m| {
        let c = m.new_empty(NodeKind::Compound);
        for s in shapes {
            m.add(c, s.shape)?;
        }
        Ok(c)
    })
    .map(handle)
}

/// A copy of `shape` with its geometry moved by `rows` (any similarity, reflections
/// included), its own placement applied too.
#[pyfunction]
fn transformed(shape: Handle, rows: Rows) -> PyResult<Handle> {
    let t = Transform::from_rows(rows).ok_or_else(|| PyValueError::new_err("a singular transformation"))?;
    with_model(|m| {
        let placement = shape.shape.placement;
        let unplaced = Shape { placement: Transform::IDENTITY, ..shape.shape };
        modify::transform(m, unplaced, &(t * placement), Default::default())
    })
    .map(handle)
}

// --- booleans ------------------------------------------------------------------------------

#[pyfunction]
fn union(a: Handle, b: Handle) -> PyResult<Handle> {
    with_model(|m| boolean::union(m, a.shape, b.shape)).map(handle)
}

#[pyfunction]
fn difference(a: Handle, b: Handle) -> PyResult<Handle> {
    with_model(|m| boolean::difference(m, a.shape, b.shape)).map(handle)
}

#[pyfunction]
fn intersection(a: Handle, b: Handle) -> PyResult<Handle> {
    with_model(|m| boolean::intersection(m, a.shape, b.shape)).map(handle)
}

#[pyfunction]
fn section(a: Handle, b: Handle) -> PyResult<Handle> {
    with_model(|m| boolean::section(m, a.shape, b.shape)).map(handle)
}

// --- queries -------------------------------------------------------------------------------

/// (volume, centre of mass).
#[pyfunction]
fn volume(shape: Handle) -> PyResult<(f64, Xyz)> {
    with_model(|m| {
        // Solids from their faces' terms, which the model keeps.
        if let Some((volume, center)) = properties::solid_volume(m, shape.shape) {
            return Ok((volume, xyz(center)));
        }
        let p = properties::volume(m, shape.shape);
        Ok((p.mass, xyz(p.center_of_mass)))
    })
}

/// (area, centre of mass).
#[pyfunction]
fn area(shape: Handle) -> PyResult<(f64, Xyz)> {
    with_model(|m| {
        let p = properties::area(m, shape.shape);
        Ok((p.mass, xyz(p.center_of_mass)))
    })
}

/// (length, centre of mass).
#[pyfunction]
fn length(shape: Handle) -> PyResult<(f64, Xyz)> {
    with_model(|m| {
        let p = properties::length(m, shape.shape);
        Ok((p.mass, xyz(p.center_of_mass)))
    })
}

/// (min corner, max corner); with `use_meshes`, faces carrying meshes count by them.
#[pyfunction]
#[pyo3(signature = (shape, use_meshes=false))]
fn bounding_box(shape: Handle, use_meshes: bool) -> PyResult<(Xyz, Xyz)> {
    with_model(|m| {
        let b = if use_meshes { properties::bounding_box_using_meshes(m, shape.shape) } else { properties::tight_bounding_box(m, shape.shape) };
        Ok((xyz(b.min), xyz(b.max)))
    })
}

/// No validity errors.
#[pyfunction]
fn is_valid(shape: Handle) -> PyResult<bool> {
    with_model(|m| Ok(validity::check(m, shape.shape).errors().next().is_none()))
}

/// The validity errors' descriptions.
#[pyfunction]
fn problems(shape: Handle) -> PyResult<Vec<String>> {
    with_model(|m| Ok(validity::check(m, shape.shape).errors().map(|p| p.description.clone()).collect()))
}

/// The distinct sub-shapes of a kind, in order of first occurrence (edges collapsed to a
/// point left out if `skip_degenerate`).
#[pyfunction]
#[pyo3(signature = (shape, kind, skip_degenerate=false))]
fn subshapes(shape: Handle, kind: &str, skip_degenerate: bool) -> PyResult<Vec<Handle>> {
    let kind = kind_of(kind)?;
    with_model(|m| {
        let mut found = m.unique_subshapes(shape.shape, kind);
        if skip_degenerate && kind == NodeKind::Edge {
            found.retain(|&e| !m.is_degenerate(e));
        }
        Ok(found.into_iter().map(handle).collect())
    })
}

/// The direct children.
#[pyfunction]
fn children(shape: Handle) -> PyResult<Vec<Handle>> {
    with_model(|m| Ok(m.children(shape.shape).map(handle).collect()))
}

#[pyfunction]
fn vertex_point(vertex: Handle) -> PyResult<Xyz> {
    with_model(|m| Ok(xyz(m.point(vertex.shape))))
}

#[pyfunction]
fn tolerance(shape: Handle) -> PyResult<f64> {
    with_model(|m| Ok(m.tolerance(shape.shape)))
}

/// The kind of a face's surface: "PLANE", "CYLINDER", "CONE", "SPHERE", "TORUS",
/// "BEZIER", "BSPLINE", "REVOLUTION", "EXTRUSION", "OFFSET" or "OTHER" (CadQuery's names).
#[pyfunction]
fn surface_type(face: Handle) -> PyResult<&'static str> {
    with_model(|m| {
        fn name(s: &Surface) -> &'static str {
            match s {
                Surface::Plane(_) => "PLANE",
                Surface::Cylinder(_) => "CYLINDER",
                Surface::Cone(_) => "CONE",
                Surface::Sphere(_) => "SPHERE",
                Surface::Torus(_) => "TORUS",
                Surface::Bezier(_) => "BEZIER",
                Surface::BSpline(_) => "BSPLINE",
                Surface::Revolution(_) => "REVOLUTION",
                Surface::Extrusion(_) => "EXTRUSION",
                Surface::Offset(_) => "OFFSET",
                Surface::Trimmed(t) => name(t.basis()),
            }
        }
        Ok(m.surface(face.shape).map_or("OTHER", name))
    })
}

/// The kind of an edge's curve: "LINE", "CIRCLE", "ELLIPSE", "HYPERBOLA", "PARABOLA",
/// "BEZIER", "BSPLINE", "OFFSET" or "OTHER".
#[pyfunction]
fn curve_type(edge: Handle) -> PyResult<&'static str> {
    with_model(|m| {
        fn name(c: &Curve) -> &'static str {
            match c {
                Curve::Line(_) => "LINE",
                Curve::Circle(_) => "CIRCLE",
                Curve::Ellipse(_) => "ELLIPSE",
                Curve::Hyperbola(_) => "HYPERBOLA",
                Curve::Parabola(_) => "PARABOLA",
                Curve::Bezier(_) => "BEZIER",
                Curve::BSpline(_) => "BSPLINE",
                Curve::Offset(_) => "OFFSET",
                Curve::Trimmed(t) => name(t.basis()),
            }
        }
        Ok(m.edge_curve(edge.shape).map_or("OTHER", |(c, _, _)| name(&c)))
    })
}

/// Points along an edge (as used: reversed edges run backwards), at fractions `ts` of its
/// range.
#[pyfunction]
fn edge_points(edge: Handle, ts: Vec<f64>) -> PyResult<Vec<Xyz>> {
    with_model(|m| {
        let Some((c, a, b)) = m.edge_curve(edge.shape) else {
            let (v, _) = m.edge_vertices(edge.shape);
            return Ok(vec![xyz(m.point(v)); ts.len()]);
        };
        let reversed = edge.shape.orientation == Orientation::Reversed;
        Ok(ts
            .iter()
            .map(|&t| {
                let s = if reversed { 1.0 - t } else { t };
                xyz(c.point(a + (b - a) * s))
            })
            .collect())
    })
}

/// The unit tangent of an edge (as used) at fraction `t` of its range.
#[pyfunction]
fn edge_tangent(edge: Handle, t: f64) -> PyResult<Xyz> {
    with_model(|m| {
        let (c, a, b) = m.edge_curve(edge.shape).ok_or_else(|| kentron::Error::InvalidInput("a degenerate edge".into()))?;
        let reversed = edge.shape.orientation == Orientation::Reversed;
        let s = if reversed { 1.0 - t } else { t };
        let d = c.d1(a + (b - a) * s).1.normalized();
        Ok(xyz(Point::ORIGIN + if reversed { -d } else { d }))
    })
}

/// (radius, centre) of a circular edge.
#[pyfunction]
fn edge_circle(edge: Handle) -> PyResult<Option<(f64, Xyz)>> {
    with_model(|m| {
        Ok(match m.edge_curve(edge.shape) {
            Some((Curve::Circle(c), _, _)) => Some((c.radius(), xyz(c.center()))),
            _ => None,
        })
    })
}

/// The face's outward unit normal (as used) at the surface point nearest `p`.
#[pyfunction]
fn face_normal_at(face: Handle, p: Xyz) -> PyResult<Xyz> {
    with_model(|m| {
        let s = m.surface(face.shape).ok_or_else(|| kentron::Error::InvalidInput("a face without a surface".into()))?.clone();
        let s = s.transformed_by(&face.shape.placement);
        let (u, v) = nearest_parameters(&s, point(p));
        let (_, du, dv) = s.d1(u, v);
        let n = du.cross(dv).normalized();
        let n = if face.shape.orientation == Orientation::Reversed { -n } else { n };
        Ok(xyz(Point::ORIGIN + n))
    })
}

/// The face's centre of its parameter range (for a normal "at the centre").
#[pyfunction]
fn face_uv_center_point(face: Handle) -> PyResult<Xyz> {
    with_model(|m| {
        let s = m.surface(face.shape).ok_or_else(|| kentron::Error::InvalidInput("a face without a surface".into()))?.clone();
        let s = s.transformed_by(&face.shape.placement);
        let (u0, u1, v0, v1) = m.face_uv_bounds(face.shape);
        Ok(xyz(s.point(0.5 * (u0 + u1), 0.5 * (v0 + v1))))
    })
}

/// The parameters of the point of `s` nearest `p`: closed-form where the surface has them,
/// refined by Gauss-Newton steps.
fn nearest_parameters(s: &Surface, p: Point) -> (f64, f64) {
    let (mut u, mut v) = s.parameters(p).unwrap_or_else(|| {
        let (u0, u1, v0, v1) = s.bounds();
        let clamp = |a: f64, b: f64| if a.is_finite() && b.is_finite() { 0.5 * (a + b) } else { 0.0 };
        (clamp(u0, u1), clamp(v0, v1))
    });
    for _ in 0..30 {
        let (q, du, dv) = s.d1(u, v);
        let r = q - p;
        let (a, b, c) = (du.dot(du), du.dot(dv), dv.dot(dv));
        let det = a * c - b * b;
        if det.abs() <= 1e-300 {
            break;
        }
        let (g1, g2) = (r.dot(du), r.dot(dv));
        let (x, y) = ((-g1 * c + g2 * b) / det, (-g2 * a + g1 * b) / det);
        u += x;
        v += y;
        if x.abs() + y.abs() <= 1e-14 {
            break;
        }
    }
    (u, v)
}

trait Placed {
    fn transformed_by(&self, t: &Transform) -> Self;
}

impl Placed for Surface {
    fn transformed_by(&self, t: &Transform) -> Surface {
        use kentron::math::Transformable;
        if *t == Transform::IDENTITY { self.clone() } else { self.transformed(t) }
    }
}

// --- files ---------------------------------------------------------------------------------

/// Reads a BRep file.
#[pyfunction]
fn read_brep(path: &str) -> PyResult<Handle> {
    let mut f = std::fs::File::open(path).map_err(|e| PyValueError::new_err(format!("{path}: {e}")))?;
    with_model(|m| kentron::io::brep::read(m, &mut f)).map(handle)
}

/// Writes a BRep file.
#[pyfunction]
fn write_brep(shape: Handle, path: &str) -> PyResult<()> {
    let mut f = std::fs::File::create(path).map_err(|e| PyValueError::new_err(format!("{path}: {e}")))?;
    with_model(|m| kentron::io::brep::write(m, shape.shape, &mut f, &Default::default()))
}

/// True for a path ending in `.binstep` (STEP in binary; read either way by its content).
fn is_binstep_path(path: &str) -> bool {
    std::path::Path::new(path).extension().is_some_and(|e| e.eq_ignore_ascii_case("binstep"))
}

/// Reads a STEP file's shapes as one (a compound of its products' shapes, placed).
#[pyfunction]
fn read_step(path: &str) -> PyResult<Handle> {
    let mut f = std::fs::File::open(path).map_err(|e| PyValueError::new_err(format!("{path}: {e}")))?;
    with_model(|m| kentron::io::step::read(m, &mut f)).map(handle)
}

/// A STEP file's products: per root, (name, shape or None, placement rows, colour or None,
/// children).
#[pyfunction]
fn read_step_assembly(py: Python<'_>, path: &str) -> PyResult<Vec<Py<PyAny>>> {
    let mut f = std::fs::File::open(path).map_err(|e| PyValueError::new_err(format!("{path}: {e}")))?;
    let roots = with_model(|m| kentron::io::step::read_assembly(m, &mut f))?;
    fn convert(py: Python<'_>, node: &kentron::io::step::Node) -> PyResult<Py<PyAny>> {
        let children = node.children.iter().map(|c| convert(py, c)).collect::<PyResult<Vec<_>>>()?;
        let shape = node.shape.map(handle);
        let color = node.color.map(|c| (c[0], c[1], c[2]));
        Ok((node.name.clone(), shape, node.placement.rows(), color, children).into_pyobject(py)?.into_any().unbind())
    }
    roots.iter().map(|n| convert(py, n)).collect()
}

/// A STEP length unit by CadQuery's name for it ("MM", "CM", "M", "KM", "INCH", "FT", "MI",
/// "UM", "NM").
fn length_unit(name: &str) -> PyResult<kentron::io::step::LengthUnit> {
    use kentron::io::step::LengthUnit::*;
    Ok(match name.to_ascii_uppercase().as_str() {
        "MM" => Millimetre,
        "CM" => Centimetre,
        "M" => Metre,
        "KM" => Kilometre,
        "UM" => Micrometre,
        "NM" => Nanometre,
        "INCH" | "IN" => Inch,
        "FT" => Foot,
        "MI" => Mile,
        other => return Err(PyValueError::new_err(format!("unknown length unit {other}"))),
    })
}

/// Writes one shape to a STEP file as a product named `name`, its lengths in `unit` (the
/// file declaring it).
/// `shape` (a wire or an edge: a wire; a face: a face) projected onto the face `target`
/// along `direction`.
#[pyfunction]
fn project_along(shape: Handle, target: Handle, direction: Xyz) -> PyResult<Handle> {
    let d = vector(direction);
    with_model(|m| match shape.shape.kind() {
        NodeKind::Face => kentron::project::face_along(m, shape.shape, target.shape, d),
        NodeKind::Edge => {
            let w = make::wire(m, &[shape.shape])?;
            kentron::project::wire_along(m, w, target.shape, d)
        }
        _ => kentron::project::wire_along(m, shape.shape, target.shape, d),
    })
    .map(handle)
}

/// The faces of `text` set in the font file at `font_path`, `size` high (its em), aligned
/// across (`left`, `center`, `right`) and up and down (`top`, `center`, `bottom`) as OCCT's
/// text builder aligns: a compound of faces in the XY plane.
#[pyfunction]
fn text(text: &str, size: f64, font_path: &str, halign: &str, valign: &str) -> PyResult<Handle> {
    use kentron::text::{HorizontalAlignment as H, VerticalAlignment as V};
    let font = std::fs::read(font_path).map_err(|e| PyValueError::new_err(format!("{font_path}: {e}")))?;
    let h = match halign {
        "left" => H::Left,
        "center" => H::Center,
        "right" => H::Right,
        other => return Err(PyValueError::new_err(format!("unknown horizontal alignment {other}"))),
    };
    let v = match valign {
        "top" => V::Top,
        "center" => V::Center,
        "bottom" => V::Bottom,
        other => return Err(PyValueError::new_err(format!("unknown vertical alignment {other}"))),
    };
    with_model(|m| kentron::text::text(m, &font, text, size, h, v)).map(handle)
}

#[pyfunction]
#[pyo3(signature = (shape, path, name, color=None, unit="MM"))]
fn write_step(shape: Handle, path: &str, name: &str, color: Option<Xyz>, unit: &str) -> PyResult<()> {
    let unit = length_unit(unit)?;
    let mut f = std::fs::File::create(path).map_err(|e| PyValueError::new_err(format!("{path}: {e}")))?;
    let color = color.map(|c| [c.0, c.1, c.2]);
    if is_binstep_path(path) {
        if unit != kentron::io::step::LengthUnit::Millimetre {
            return Err(PyNotImplementedError::new_err("binary STEP in units other than millimetres"));
        }
        return with_model(|m| kentron::io::step::write_binstep(m, shape.shape, name, color, &mut f));
    }
    with_model(|m| kentron::io::step::write_in(m, shape.shape, name, color, unit, &mut f))
}

/// Writes products to a STEP file: nested `(name, shape or None, placement rows, rgb or
/// None, children)` tuples, as `read_step_assembly` gives them.
#[pyfunction]
#[pyo3(signature = (path, roots, unit="MM"))]
fn write_step_assembly(path: &str, roots: Vec<Bound<'_, PyAny>>, unit: &str) -> PyResult<()> {
    let unit = length_unit(unit)?;
    fn convert(node: &Bound<'_, PyAny>) -> PyResult<kentron::io::step::Node> {
        let (name, shape, rows, color, children): (String, Option<Handle>, [[f64; 4]; 3], Option<Xyz>, Vec<Bound<'_, PyAny>>) = node.extract()?;
        let placement = Transform::from_rows(rows).ok_or_else(|| PyValueError::new_err("a singular placement"))?;
        Ok(kentron::io::step::Node {
            name,
            shape: shape.map(|h| h.shape),
            color: color.map(|c| [c.0, c.1, c.2]),
            placement,
            children: children.iter().map(convert).collect::<PyResult<Vec<_>>>()?,
        })
    }
    let roots = roots.iter().map(convert).collect::<PyResult<Vec<_>>>()?;
    let mut f = std::fs::File::create(path).map_err(|e| PyValueError::new_err(format!("{path}: {e}")))?;
    if is_binstep_path(path) {
        if unit != kentron::io::step::LengthUnit::Millimetre {
            return Err(PyNotImplementedError::new_err("binary STEP in units other than millimetres"));
        }
        return with_model(|m| kentron::io::step::write_assembly_binstep(m, &roots, &mut f));
    }
    with_model(|m| kentron::io::step::write_assembly_in(m, &roots, unit, &mut f))
}

#[pymodule]
fn _kentron(m: &Bound<'_, PyModule>) -> PyResult<()> {
    profile::check_environment()?;
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    m.add_function(wrap_pyfunction!(profile::profile, m)?)?;
    m.add_function(wrap_pyfunction!(profile::set_profile, m)?)?;
    m.add_class::<Handle>()?;
    m.add_class::<Mesh>()?;
    m.add_function(wrap_pyfunction!(is_closed, m)?)?;
    m.add_function(wrap_pyfunction!(make_box, m)?)?;
    m.add_function(wrap_pyfunction!(make_cylinder, m)?)?;
    m.add_function(wrap_pyfunction!(make_cone, m)?)?;
    m.add_function(wrap_pyfunction!(make_sphere, m)?)?;
    m.add_function(wrap_pyfunction!(make_torus, m)?)?;
    m.add_function(wrap_pyfunction!(make_vertex, m)?)?;
    m.add_function(wrap_pyfunction!(make_line, m)?)?;
    m.add_function(wrap_pyfunction!(make_circle, m)?)?;
    m.add_function(wrap_pyfunction!(make_three_point_arc, m)?)?;
    m.add_function(wrap_pyfunction!(make_wire, m)?)?;
    m.add_function(wrap_pyfunction!(make_polygon, m)?)?;
    m.add_function(wrap_pyfunction!(make_face, m)?)?;
    m.add_function(wrap_pyfunction!(extrude, m)?)?;
    m.add_function(wrap_pyfunction!(revolve, m)?)?;
    m.add_function(wrap_pyfunction!(make_loft, m)?)?;
    m.add_function(wrap_pyfunction!(offset_wire, m)?)?;
    m.add_function(wrap_pyfunction!(offset_shape, m)?)?;
    m.add_function(wrap_pyfunction!(thick_solid, m)?)?;
    m.add_function(wrap_pyfunction!(unify_same_domain, m)?)?;
    m.add_function(wrap_pyfunction!(triangulate, m)?)?;
    m.add_function(wrap_pyfunction!(mesh, m)?)?;
    m.add_function(wrap_pyfunction!(voxelize, m)?)?;
    m.add_function(wrap_pyfunction!(morphology, m)?)?;
    m.add_function(wrap_pyfunction!(components, m)?)?;
    m.add_function(wrap_pyfunction!(text, m)?)?;
    m.add_function(wrap_pyfunction!(project_along, m)?)?;
    m.add_function(wrap_pyfunction!(write_step, m)?)?;
    m.add_function(wrap_pyfunction!(write_step_assembly, m)?)?;
    m.add_function(wrap_pyfunction!(distance, m)?)?;
    m.add_function(wrap_pyfunction!(classify_point, m)?)?;
    m.add_function(wrap_pyfunction!(write_stl, m)?)?;
    m.add_function(wrap_pyfunction!(fillet, m)?)?;
    m.add_function(wrap_pyfunction!(make_compound, m)?)?;
    m.add_function(wrap_pyfunction!(transformed, m)?)?;
    m.add_function(wrap_pyfunction!(union, m)?)?;
    m.add_function(wrap_pyfunction!(difference, m)?)?;
    m.add_function(wrap_pyfunction!(intersection, m)?)?;
    m.add_function(wrap_pyfunction!(section, m)?)?;
    m.add_function(wrap_pyfunction!(volume, m)?)?;
    m.add_function(wrap_pyfunction!(area, m)?)?;
    m.add_function(wrap_pyfunction!(length, m)?)?;
    m.add_function(wrap_pyfunction!(bounding_box, m)?)?;
    m.add_function(wrap_pyfunction!(is_valid, m)?)?;
    m.add_function(wrap_pyfunction!(problems, m)?)?;
    m.add_function(wrap_pyfunction!(subshapes, m)?)?;
    m.add_function(wrap_pyfunction!(children, m)?)?;
    m.add_function(wrap_pyfunction!(vertex_point, m)?)?;
    m.add_function(wrap_pyfunction!(tolerance, m)?)?;
    m.add_function(wrap_pyfunction!(surface_type, m)?)?;
    m.add_function(wrap_pyfunction!(curve_type, m)?)?;
    m.add_function(wrap_pyfunction!(edge_points, m)?)?;
    m.add_function(wrap_pyfunction!(edge_tangent, m)?)?;
    m.add_function(wrap_pyfunction!(edge_circle, m)?)?;
    m.add_function(wrap_pyfunction!(face_normal_at, m)?)?;
    m.add_function(wrap_pyfunction!(face_uv_center_point, m)?)?;
    m.add_function(wrap_pyfunction!(read_brep, m)?)?;
    m.add_function(wrap_pyfunction!(write_brep, m)?)?;
    m.add_function(wrap_pyfunction!(read_step, m)?)?;
    m.add_function(wrap_pyfunction!(read_step_assembly, m)?)?;
    cadquery::register(m)?;
    Ok(())
}
