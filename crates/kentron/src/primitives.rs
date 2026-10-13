//! Primitive solids: box, wedge, cylinder, cone, sphere, torus.
//!
//! Each primitive is described by a small builder value and created in a model with
//! `build`:
//!
//! ```no_run
//! # use kentron::{primitives, topology::Model};
//! let mut model = Model::new();
//! let half_cylinder = primitives::Cylinder::new(5.0, 10.0)
//!     .with_angle(std::f64::consts::PI)
//!     .build(&mut model)?;
//! # Ok::<(), kentron::Error>(())
//! ```
//!
//! Primitives are positioned by a [`Frame`] (default: [`Frame::WORLD`]). Rotational
//! primitives revolve about the frame's Z axis, starting from its X axis; `with_angle`
//! limits the revolution, which adds two planar side faces.

use std::f64::consts::{FRAC_PI_2, TAU};

use std::collections::HashMap;

use crate::error::{Error, Result};
use crate::geometry::curve::{Circle, Curve, Line};
use crate::geometry::curve2::{Curve2, Line2};
use crate::geometry::project;
use crate::geometry::surface::{
    Cone as ConeSurface, Cylinder as CylinderSurface, Plane, Sphere as SphereSurface, Surface, Torus as TorusSurface,
};
use crate::geometry::ParametricSurface;
use crate::math::{Direction, Direction2, Frame, Point, Point2, Vector};
use crate::tolerance;
use crate::topology::{EdgeId, Model, NodeKind, Orientation, Shape, VertexId};

/// A rectangular box: `[0, dx] x [0, dy] x [0, dz]` in its frame.
#[doc(alias = "BRepPrimAPI_MakeBox")]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Cuboid {
    frame: Frame,
    dx: f64,
    dy: f64,
    dz: f64,
}

impl Cuboid {
    pub fn new(dx: f64, dy: f64, dz: f64) -> Self {
        Self { frame: Frame::WORLD, dx, dy, dz }
    }

    /// The axis-aligned box with opposite corners `p1` and `p2`.
    pub fn from_corners(p1: Point, p2: Point) -> Self {
        let origin = Point::new(p1.x.min(p2.x), p1.y.min(p2.y), p1.z.min(p2.z));
        let frame = Frame::new(origin, Direction::Z, Direction::X);
        Self { frame, dx: (p2.x - p1.x).abs(), dy: (p2.y - p1.y).abs(), dz: (p2.z - p1.z).abs() }
    }

    pub fn with_frame(self, frame: Frame) -> Self {
        Self { frame, ..self }
    }

    /// Fails unless all dimensions are positive.
    /// Builds the shape in `model`; its nodes start unedited (revision 0).
    pub fn build(&self, model: &mut Model) -> Result<Shape> {
        let shape = self.build_fresh(model)?;
        model.settle(shape);
        Ok(shape)
    }

    fn build_fresh(&self, model: &mut Model) -> Result<Shape> {
        if !(self.dx > 0.0 && self.dy > 0.0 && self.dz > 0.0) {
            return Err(Error::InvalidInput(format!(
                "box dimensions must be positive, got {} x {} x {}",
                self.dx, self.dy, self.dz
            )));
        }
        let f = &self.frame;
        let (x, y, z) = (Vector::from(f.x()), Vector::from(f.y()), Vector::from(f.z()));
        // Corner `n` has offsets (n & 1, n >> 1 & 1, n >> 2 & 1) along (x, y, z).
        let corners: [Point; 8] = std::array::from_fn(|n| {
            let (i, j, k) = ((n & 1) as f64, ((n >> 1) & 1) as f64, ((n >> 2) & 1) as f64);
            f.origin() + x * (self.dx * i) + y * (self.dy * j) + z * (self.dz * k)
        });
        hexahedron(model, f, corners)
    }
}

/// The solid with the eight corners of a box (corner `n` at the low or high end of each
/// axis by the bits of `n`: x, y, z), each face planar; corners may coincide (a wedge's top
/// collapsing to an edge), dropping the edges between them and faces left without area.
///
/// The faces are stored in OCCT's order (XMin, XMax, YMin, YMax, ZMin, ZMax) and, for
/// each, the faces its edges are shared with in the order of OCCT's wires
/// (BRepPrim_GWedge), so that exploring gives the same order. Each edge runs from its low
/// corner to its high one; the wire follows its first two edges forward and the last two
/// backward, which turns clockwise about the axis: it bounds the Max faces reversed.
fn hexahedron(model: &mut Model, frame: &Frame, corners: [Point; 8]) -> Result<Shape> {
    const NEIGHBOURS: [[(usize, usize); 4]; 3] = [
        [(1, 0), (2, 1), (1, 1), (2, 0)],
        [(2, 0), (0, 1), (2, 1), (0, 0)],
        [(0, 0), (1, 1), (0, 1), (1, 0)],
    ];
    let (x, y, z) = (Vector::from(frame.x()), Vector::from(frame.y()), Vector::from(frame.z()));
    let right_handed = x.cross(y).dot(z) > 0.0;
    let axes = [x, y, z];

    let mut b = Builder::new(model);
    // One vertex per distinct point.
    let mut points: Vec<(Point, VertexId)> = Vec::new();
    let mut vertex_of = [0usize; 8];
    for (n, &p) in corners.iter().enumerate() {
        vertex_of[n] = match points.iter().position(|(q, _)| q.distance(p) <= tolerance::LINEAR) {
            Some(i) => i,
            None => {
                points.push((p, b.vertex(p)));
                points.len() - 1
            }
        };
    }
    // Edges by their vertices, with the vertex they start from.
    let mut edges: HashMap<(usize, usize), (EdgeId, usize)> = HashMap::new();
    for (axis, side) in (0..3).flat_map(|a| [(a, 0), (a, 1)]) {
        let mut uses = Vec::with_capacity(4);
        for (k, &(other, other_side)) in NEIGHBOURS[axis].iter().enumerate() {
            // The edge's corners: on both faces, low and high along the third axis.
            let third = 3 - axis - other;
            let base = (side << axis) | (other_side << other);
            let (from, to) = (vertex_of[base], vertex_of[base | (1 << third)]);
            if from == to {
                continue;
            }
            let (edge, start) = match edges.get(&(from.min(to), from.max(to))) {
                Some(&e) => e,
                None => {
                    let edge = b.segment(points[from].1, points[to].1)?;
                    edges.insert((from.min(to), from.max(to)), (edge, from));
                    (edge, from)
                }
            };
            let forward = (start == from) == (k < 2);
            if !uses.iter().any(|u: &Use| u.edge == edge) {
                uses.push(if forward { fwd(edge) } else { rev(edge) });
            }
        }
        if uses.len() < 3 {
            continue;
        }
        // Clockwise about the axis: right for a Min face (outward normal against the axis)
        // of a right-handed frame.
        let reversed = (side == 1) == right_handed;
        // The face's plane, through its corners, its normal outward (along the axis's
        // direction for a Max face).
        let ring: Vec<Point> = (0..4)
            .map(|i| {
                let (a, b) = ((axis + 1) % 3, (axis + 2) % 3);
                let bits = [(0, 0), (1, 0), (1, 1), (0, 1)][i];
                corners[(side << axis) | (bits.0 << a) | (bits.1 << b)]
            })
            .collect();
        let mut normal = Vector::new(0.0, 0.0, 0.0);
        for i in 0..4 {
            normal = normal + (ring[i] - ring[0]).cross(ring[(i + 1) % 4] - ring[0]);
        }
        let toward = axes[axis] * if side == 1 { 1.0 } else { -1.0 };
        let outward = if normal.dot(toward) < 0.0 { -normal } else { normal };
        let along = axes[(axis + 1) % 3];
        let in_plane = along - outward.normalized() * along.dot(outward.normalized());
        b.face_with_wires(plane(ring[0], outward, in_plane), vec![(uses, reversed)])?;
    }
    b.solid()
}

/// A box whose top face (at `y = dy`) is the rectangle `[x_min, x_max] x [z_min, z_max]`;
/// the bottom face is `[0, dx] x [0, dz]` at `y = 0`. A degenerate top (zero width)
/// collapses to an edge, leaving five faces.
#[doc(alias = "BRepPrimAPI_MakeWedge")]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Wedge {
    frame: Frame,
    dx: f64,
    dy: f64,
    dz: f64,
    x_min: f64,
    z_min: f64,
    x_max: f64,
    z_max: f64,
}

impl Wedge {
    /// The top face spans `[0, top_dx] x [0, dz]`.
    pub fn new(dx: f64, dy: f64, dz: f64, top_dx: f64) -> Self {
        Self::with_top(dx, dy, dz, 0.0, 0.0, top_dx, dz)
    }

    pub fn with_top(
        dx: f64,
        dy: f64,
        dz: f64,
        x_min: f64,
        z_min: f64,
        x_max: f64,
        z_max: f64,
    ) -> Self {
        Self { frame: Frame::WORLD, dx, dy, dz, x_min, z_min, x_max, z_max }
    }

    pub fn with_frame(self, frame: Frame) -> Self {
        Self { frame, ..self }
    }

    /// Fails unless `dx`, `dy` and `dz` are positive and the top is a rectangle, a segment
    /// or a point (`x_min <= x_max`, `z_min <= z_max`).
    /// Builds the shape in `model`; its nodes start unedited (revision 0).
    pub fn build(&self, model: &mut Model) -> Result<Shape> {
        let shape = self.build_fresh(model)?;
        model.settle(shape);
        Ok(shape)
    }

    fn build_fresh(&self, model: &mut Model) -> Result<Shape> {
        if !(self.dx > 0.0 && self.dy > 0.0 && self.dz > 0.0 && self.x_min <= self.x_max && self.z_min <= self.z_max) {
            return Err(Error::InvalidInput(format!("invalid wedge {self:?}")));
        }
        let f = &self.frame;
        let (x, y, z) = (Vector::from(f.x()), Vector::from(f.y()), Vector::from(f.z()));
        // As for a box: corner `n` at the low or high end of x, y, z by its bits; the top
        // (high y) corners span the top rectangle.
        let corners: [Point; 8] = std::array::from_fn(|n| {
            let (i, j, k) = (n & 1, (n >> 1) & 1, (n >> 2) & 1);
            if j == 0 {
                f.origin() + x * (self.dx * i as f64) + z * (self.dz * k as f64)
            } else {
                let xs = if i == 0 { self.x_min } else { self.x_max };
                let zs = if k == 0 { self.z_min } else { self.z_max };
                f.origin() + x * xs + y * self.dy + z * zs
            }
        });
        hexahedron(model, f, corners)
    }
}

/// A cylinder of the given radius along the frame's Z axis, from `z = 0` to `z = height`.
#[doc(alias = "BRepPrimAPI_MakeCylinder")]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Cylinder {
    frame: Frame,
    radius: f64,
    height: f64,
    angle: f64,
}

impl Cylinder {
    pub fn new(radius: f64, height: f64) -> Self {
        Self { frame: Frame::WORLD, radius, height, angle: TAU }
    }

    pub fn with_frame(self, frame: Frame) -> Self {
        Self { frame, ..self }
    }

    /// Revolves only through `angle` (radians, in (0, 2*PI]).
    pub fn with_angle(self, angle: f64) -> Self {
        Self { angle, ..self }
    }

    /// Fails unless the radius and height are positive and the angle is in (0, 2π].
    /// Builds the shape in `model`; its nodes start unedited (revision 0).
    pub fn build(&self, model: &mut Model) -> Result<Shape> {
        let shape = self.build_fresh(model)?;
        model.settle(shape);
        Ok(shape)
    }

    fn build_fresh(&self, model: &mut Model) -> Result<Shape> {
        if !(self.radius > 0.0 && self.height > 0.0) {
            return Err(Error::InvalidInput(format!(
                "cylinder radius and height must be positive, got {} and {}",
                self.radius, self.height
            )));
        }
        let (r, h) = (self.radius, self.height);
        Revolved {
            frame: self.frame,
            surface: Surface::Cylinder(CylinderSurface::new(self.frame, r)),
            v_range: (0.0, h),
            ends: [(r, 0.0), (r, h)],
            angle: self.angle,
        }
        .build(model)
    }
}

/// A (possibly truncated) cone along the frame's Z axis: radius `bottom_radius` at `z = 0`,
/// `top_radius` at `z = height`. A zero radius makes an apex.
#[doc(alias = "BRepPrimAPI_MakeCone")]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Cone {
    frame: Frame,
    bottom_radius: f64,
    top_radius: f64,
    height: f64,
    angle: f64,
}

impl Cone {
    pub fn new(bottom_radius: f64, top_radius: f64, height: f64) -> Self {
        Self { frame: Frame::WORLD, bottom_radius, top_radius, height, angle: TAU }
    }

    pub fn with_frame(self, frame: Frame) -> Self {
        Self { frame, ..self }
    }

    /// Revolves only through `angle` (radians, in (0, 2*PI]).
    pub fn with_angle(self, angle: f64) -> Self {
        Self { angle, ..self }
    }

    /// Fails unless the height is positive, the radii are non-negative, different, and
    /// not both zero, and the angle is in (0, 2π]. A zero radius makes an apex.
    /// Builds the shape in `model`; its nodes start unedited (revision 0).
    pub fn build(&self, model: &mut Model) -> Result<Shape> {
        let shape = self.build_fresh(model)?;
        model.settle(shape);
        Ok(shape)
    }

    fn build_fresh(&self, model: &mut Model) -> Result<Shape> {
        let (r1, r2, h) = (self.bottom_radius, self.top_radius, self.height);
        if !(h > 0.0 && r1 >= 0.0 && r2 >= 0.0) || (r1 - r2).abs() <= tolerance::LINEAR {
            return Err(Error::InvalidInput(format!(
                "a cone needs a positive height and two different non-negative radii, got {r1}, {r2} and {h}"
            )));
        }
        // The side's generating lines, from the bottom circle to the top one.
        let semi_angle = ((r2 - r1) / h).atan();
        let slant = h / semi_angle.cos();
        Revolved {
            frame: self.frame,
            surface: Surface::Cone(ConeSurface::new(self.frame, semi_angle, r1)),
            v_range: (0.0, slant),
            ends: [(r1, 0.0), (r2, h)],
            angle: self.angle,
        }
        .build(model)
    }
}

/// A sphere centered at the frame's origin, optionally limited in latitude (from the XY
/// plane, in [-PI/2, PI/2]) and longitude (`with_angle`, from the X axis).
#[doc(alias = "BRepPrimAPI_MakeSphere")]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Sphere {
    frame: Frame,
    radius: f64,
    min_latitude: f64,
    max_latitude: f64,
    angle: f64,
}

impl Sphere {
    pub fn new(radius: f64) -> Self {
        Self {
            frame: Frame::WORLD,
            radius,
            min_latitude: -FRAC_PI_2,
            max_latitude: FRAC_PI_2,
            angle: TAU,
        }
    }

    pub fn with_frame(self, frame: Frame) -> Self {
        Self { frame, ..self }
    }

    /// Keeps only the zone between the two latitudes.
    pub fn with_latitudes(self, min_latitude: f64, max_latitude: f64) -> Self {
        Self { min_latitude, max_latitude, ..self }
    }

    /// Revolves only through `angle` (radians, in (0, 2*PI]).
    pub fn with_angle(self, angle: f64) -> Self {
        Self { angle, ..self }
    }

    /// Fails unless the radius is positive, the latitudes are increasing within
    /// [-π/2, π/2], and the angle is in (0, 2π].
    /// Builds the shape in `model`; its nodes start unedited (revision 0).
    pub fn build(&self, model: &mut Model) -> Result<Shape> {
        let shape = self.build_fresh(model)?;
        model.settle(shape);
        Ok(shape)
    }

    fn build_fresh(&self, model: &mut Model) -> Result<Shape> {
        let (r, v1, v2) = (self.radius, self.min_latitude, self.max_latitude);
        if !(r > 0.0) {
            return Err(Error::InvalidInput(format!("sphere radius must be positive, got {r}")));
        }
        let slack = tolerance::PARAMETRIC;
        if !(v1 < v2 && v1 >= -FRAC_PI_2 - slack && v2 <= FRAC_PI_2 + slack) {
            return Err(Error::InvalidInput(format!("invalid sphere latitudes [{v1}, {v2}]")));
        }
        let (v1, v2) = (v1.max(-FRAC_PI_2), v2.min(FRAC_PI_2));
        // At a pole the parallel collapses onto the axis.
        let parallel = |v: f64| {
            let radius = if v.abs() >= FRAC_PI_2 { 0.0 } else { r * v.cos() };
            (radius, r * v.sin())
        };
        Revolved {
            frame: self.frame,
            surface: Surface::Sphere(SphereSurface::new(self.frame, r)),
            v_range: (v1, v2),
            ends: [parallel(v1), parallel(v2)],
            angle: self.angle,
        }
        .build(model)
    }
}

/// A torus about the frame's Z axis (see [`crate::geometry::surface::Torus`] for its
/// parametrization). `with_angle` limits the revolution. `with_tube_angles` replaces the
/// tube cross-section by the region bounded by the tube arc between the two angles, the
/// horizontal segments from the arc's ends toward the axis, and the axis. When both ends
/// are at the same height the segments overlap and the region is the one between the arc
/// and its chord.
#[doc(alias = "BRepPrimAPI_MakeTorus")]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Torus {
    frame: Frame,
    major_radius: f64,
    minor_radius: f64,
    tube_min: f64,
    tube_max: f64,
    angle: f64,
}

impl Torus {
    pub fn new(major_radius: f64, minor_radius: f64) -> Self {
        Self { frame: Frame::WORLD, major_radius, minor_radius, tube_min: 0.0, tube_max: TAU, angle: TAU }
    }

    pub fn with_frame(self, frame: Frame) -> Self {
        Self { frame, ..self }
    }

    /// Keeps the part of the tube cross-section between the angles `min` and `max`
    /// (radians, measured from the frame's XY plane outward; PI/2 is the top of the tube).
    pub fn with_tube_angles(self, min: f64, max: f64) -> Self {
        Self { tube_min: min, tube_max: max, ..self }
    }

    /// Revolves only through `angle` (radians, in (0, 2*PI]).
    pub fn with_angle(self, angle: f64) -> Self {
        Self { angle, ..self }
    }

    /// Fails unless `0 < minor_radius < major_radius` and the angle is in (0, 2π].
    /// Builds the shape in `model`; its nodes start unedited (revision 0).
    pub fn build(&self, model: &mut Model) -> Result<Shape> {
        let shape = self.build_fresh(model)?;
        model.settle(shape);
        Ok(shape)
    }

    fn build_fresh(&self, model: &mut Model) -> Result<Shape> {
        let (big, small) = (self.major_radius, self.minor_radius);
        let full_tube = self.tube_max - self.tube_min >= TAU - tolerance::PARAMETRIC;
        if !(small > 0.0 && (small < big || (small <= big && !full_tube))) {
            return Err(Error::InvalidInput(format!("invalid torus radii {big}, {small}")));
        }
        if !full_tube {
            // Like a sphere's zone: the tube's arc, closed by segments to the axis.
            let (v1, v2) = (self.tube_min, self.tube_max);
            if !(v1 < v2) {
                return Err(Error::InvalidInput(format!("invalid torus tube angles [{v1}, {v2}]")));
            }
            let end = |v: f64| ((big + small * v.cos()).max(0.0), small * v.sin());
            return Revolved {
                frame: self.frame,
                surface: Surface::Torus(TorusSurface::new(self.frame, big, small)),
                v_range: (v1, v2),
                ends: [end(v1), end(v2)],
                angle: self.angle,
            }
            .build(model);
        }
        // The disk of the tube's cross-section in the XZ half-plane, revolved.
        let f = self.frame;
        let center = f.origin() + Vector::from(f.x()) * big;
        let circle = Circle::new(Frame::new(center, f.y(), f.x()), small);
        let edge = crate::make::edge_from_curve(model, Curve::Circle(circle))?;
        let wire = crate::make::wire(model, &[edge])?;
        let disk = crate::make::planar_face(model, wire)?;
        Ok(crate::sweep::revolve(model, disk, f.z_axis(), self.angle)?.shape)
    }
}

/// A solid of revolution about the Z axis of `frame`: the side face covers
/// `[0, angle] x v_range` of `surface`, a surface of revolution about that axis whose
/// u-isos are the profile (lines for cylinders and cones, meridians for spheres). Each end
/// of the profile is a circle of the given radius at the given height (closed by a planar
/// cap) or, for a zero radius, a point on the axis (a degenerate edge of the side). A
/// partial revolution adds two planar sides through the axis.
///
/// The surface must be oriented so that its natural normal points out of the solid, and
/// heights must increase with `v`.
struct Revolved {
    frame: Frame,
    surface: Surface,
    v_range: (f64, f64),
    /// (radius, height) at the start and end of `v_range`.
    ends: [(f64, f64); 2],
    angle: f64,
}

impl Revolved {
    fn build(&self, model: &mut Model) -> Result<Shape> {
        check_angle(self.angle)?;
        let angle = self.angle.min(TAU);
        let full = angle >= TAU - tolerance::PARAMETRIC;
        let f = self.frame;
        let (x, y, z) = (Vector::from(f.x()), Vector::from(f.y()), Vector::from(f.z()));
        let (v_lo, v_hi) = self.v_range;
        let on_axis = |height: f64| f.origin() + z * height;
        let is_point = |end: usize| self.ends[end].0 <= tolerance::LINEAR;
        let surface = &self.surface;
        // The meridian with the ends' segments to the axis and the axis bound the profile.
        // Its normal (∂u × ∂v) points outward where that loop turns counter-clockwise in
        // (radius, height); otherwise (an arc of a torus turning towards the axis, say) the
        // material lies on the other side: the side and cap faces are used reversed, and
        // the sides' boundaries turn the other way.
        let meridian: Vec<(f64, f64)> = (0..=64)
            .map(|i| {
                let p = surface.point(0.0, v_lo + (v_hi - v_lo) * i as f64 / 64.0) - f.origin();
                (p.dot(x), p.dot(z))
            })
            .chain([(0.0, self.ends[1].1), (0.0, self.ends[0].1)])
            .collect();
        let twice_area: f64 = (0..meridian.len())
            .map(|i| {
                let (a, b) = (meridian[i], meridian[(i + 1) % meridian.len()]);
                a.0 * b.1 - b.0 * a.1
            })
            .sum();
        let flip = twice_area < 0.0;
        if !full && !is_point(0) && !is_point(1) && (self.ends[0].1 - self.ends[1].1).abs() <= tolerance::LINEAR {
            todo!("solids of revolution whose profile's ends are at the same height");
        }

        let mut b = Builder::new(model);
        // Per end: the axis vertex, the side vertices at u = 0 and u = angle (the axis
        // vertex where the end is a point), and the parallel or degenerate edge.
        let mut axis_vertex = [VertexId::INVALID; 2];
        let mut side_vertices = [[VertexId::INVALID; 2]; 2];
        let mut parallel = [EdgeId::INVALID; 2];
        for (end, &v) in [v_lo, v_hi].iter().enumerate() {
            let (radius, height) = self.ends[end];
            if is_point(end) {
                let apex = b.vertex(on_axis(height));
                axis_vertex[end] = apex;
                side_vertices[end] = [apex, apex];
                parallel[end] = b.degenerate_edge(apex)?;
                continue;
            }
            if !full {
                axis_vertex[end] = b.vertex(on_axis(height));
            }
            let start = b.vertex(surface.point(0.0, v));
            let stop = if full { start } else { b.vertex(surface.point(angle, v)) };
            side_vertices[end] = [start, stop];
            let circle = Circle::new(Frame::new(on_axis(height), f.z(), f.x()), radius);
            parallel[end] = b.edge(start, stop, Curve::Circle(circle), 0.0, angle)?;
        }
        let generator = |b: &mut Builder, u: f64, i: usize| {
            b.edge(side_vertices[0][i], side_vertices[1][i], surface.u_iso(u), v_lo, v_hi)
        };
        let generator0 = generator(&mut b, 0.0, 0)?;
        let generator1 = if full { generator0 } else { generator(&mut b, angle, 1)? };

        // Faces and wires are stored in OCCT's order (BRepPrim_OneAxis: lateral, top,
        // bottom, start and end faces), so that exploring gives the same order. OCCT stores
        // each wire clockwise about the outward normal; here the same edges in the same order
        // run clockwise, the wire being used reversed.
        //
        // The side face: top, end generator, bottom, start generator.
        b.oriented_face(
            surface.clone(),
            vec![(
                vec![
                    fwd(parallel[1]).on_line((0.0, v_hi), (1.0, 0.0), 0.0, angle),
                    rev(generator1).on_line((angle, 0.0), (0.0, 1.0), v_lo, v_hi),
                    rev(parallel[0]).on_line((0.0, v_lo), (1.0, 0.0), 0.0, angle),
                    fwd(generator0).on_line((0.0, 0.0), (0.0, 1.0), v_lo, v_hi),
                ],
                true,
            )],
            flip,
        )?;

        let [bottom_center, top_center] = [on_axis(self.ends[0].1), on_axis(self.ends[1].1)];
        if full {
            if !is_point(1) {
                b.oriented_face(plane(top_center, z, x), vec![(vec![fwd(parallel[1])], false)], flip)?;
            }
            if !is_point(0) {
                b.oriented_face(plane(bottom_center, -z, x), vec![(vec![rev(parallel[0])], false)], flip)?;
            }
            return b.solid();
        }

        // Sector caps, then the two sides through the axis.
        let axis = b.segment(axis_vertex[0], axis_vertex[1])?;
        let mut radials = [[None; 2]; 2];
        for end in 0..2 {
            if !is_point(end) {
                for i in 0..2 {
                    radials[end][i] = Some(b.segment(axis_vertex[end], side_vertices[end][i])?);
                }
            }
        }
        // Top: the parallel, the start radial, the end radial; bottom: the parallel, the
        // end radial, the start radial.
        if let [Some(top0), Some(top1)] = radials[1] {
            b.oriented_face(plane(top_center, z, x), vec![(vec![rev(parallel[1]), rev(top0), fwd(top1)], true)], flip)?;
        }
        if let [Some(bottom0), Some(bottom1)] = radials[0] {
            b.oriented_face(plane(bottom_center, -z, x), vec![(vec![fwd(parallel[0]), rev(bottom1), fwd(bottom0)], true)], flip)?;
        }
        let (sin, cos) = angle.sin_cos();
        let radial = x * cos + y * sin;
        let tangent = y * cos - x * sin;
        // At u = 0 (outward normal -y): the bottom radial, the axis, the top radial, the
        // generator.
        let mut side0 = Vec::new();
        side0.extend(radials[0][0].map(rev));
        side0.push(fwd(axis));
        side0.extend(radials[1][0].map(fwd));
        side0.push(rev(generator0));
        b.face_with_wires(plane(bottom_center, -y, x), vec![(side0, !flip)])?;
        // At u = angle (outward normal along the direction of revolution): the top radial,
        // the axis, the bottom radial, the generator.
        let mut side1 = Vec::new();
        side1.extend(radials[1][1].map(rev));
        side1.push(rev(axis));
        side1.extend(radials[0][1].map(fwd));
        side1.push(fwd(generator1));
        b.face_with_wires(plane(bottom_center, tangent, radial), vec![(side1, !flip)])?;
        b.solid()
    }
}

/// Builds the topology of a primitive: vertices, edges with their curves, faces with
/// their pcurves (computed for planar faces), and the closed solid.
struct Builder<'m> {
    model: &'m mut Model,
    tolerance: f64,
    shell: Shape,
}

/// An edge as used by a face boundary: its direction of travel and, unless it can be
/// derived (straight or circular edges on planes), its pcurve with range.
struct Use {
    edge: EdgeId,
    orientation: Orientation,
    pcurve: Option<(Curve2, f64, f64)>,
}

fn fwd(edge: EdgeId) -> Use {
    Use { edge, orientation: Orientation::Forward, pcurve: None }
}

fn rev(edge: EdgeId) -> Use {
    Use { edge, orientation: Orientation::Reversed, pcurve: None }
}

impl Use {
    /// The pcurve `origin + t·direction` over `[first, last]`.
    fn on_line(mut self, origin: (f64, f64), direction: (f64, f64), first: f64, last: f64) -> Use {
        let line = Line2::new(Point2::new(origin.0, origin.1), Direction2::new(direction.0, direction.1));
        self.pcurve = Some((Curve2::Line(line), first, last));
        self
    }
}

impl<'m> Builder<'m> {
    fn new(model: &'m mut Model) -> Self {
        let shell = model.new_empty(NodeKind::Shell);
        Self { model, tolerance: tolerance::LINEAR, shell }
    }

    fn vertex(&mut self, p: Point) -> VertexId {
        self.model.add_vertex(p, self.tolerance)
    }

    /// An edge from `start` to `end` along `curve` over `[first, last]`.
    fn edge(&mut self, start: VertexId, end: VertexId, curve: Curve, first: f64, last: f64) -> Result<EdgeId> {
        let edge = self.model.add_edge(start, end, self.tolerance)?;
        self.model.set_edge_curve(edge.into(), curve, self.tolerance)?;
        self.model.set_edge_range(edge.into(), first, last)?;
        Ok(edge)
    }

    /// The straight edge between two vertices.
    fn segment(&mut self, start: VertexId, end: VertexId) -> Result<EdgeId> {
        let (a, b) = (self.point(start), self.point(end));
        let line = Line::new(a, Direction::from_vector(b - a));
        self.edge(start, end, Curve::Line(line), 0.0, a.distance(b))
    }

    /// An edge collapsed to `vertex` (no 3D curve), e.g. at a sphere's pole.
    fn degenerate_edge(&mut self, vertex: VertexId) -> Result<EdgeId> {
        self.model.add_edge(vertex, vertex, self.tolerance)
    }

    fn point(&self, vertex: VertexId) -> Point {
        self.model.get(vertex).expect("vertex").point
    }

    /// A face on `surface` bounded by `loops` (outer first), added to the shell.
    fn face(&mut self, surface: Surface, loops: Vec<Vec<Use>>) -> Result<Shape> {
        self.face_with_wires(surface, loops.into_iter().map(|uses| (uses, false)).collect())
    }

    /// A face bounded by wires, each used reversed if its flag is set.
    fn face_with_wires(&mut self, surface: Surface, loops: Vec<(Vec<Use>, bool)>) -> Result<Shape> {
        self.oriented_face(surface, loops, false)
    }

    /// A face bounded by wires (each used reversed if its flag is set), itself used
    /// reversed in the shell if `reversed`.
    fn oriented_face(&mut self, surface: Surface, loops: Vec<(Vec<Use>, bool)>, reversed: bool) -> Result<Shape> {
        let face = self.model.new_empty(NodeKind::Face);
        self.model.set_face_surface(face, surface.clone(), self.tolerance)?;
        for (uses, reversed) in loops {
            let mut coedges = Vec::with_capacity(uses.len());
            for u in uses {
                let coedge = self.model.add_coedge(u.edge, u.orientation)?;
                let (curve, first, last) = match u.pcurve {
                    Some(pcurve) => pcurve,
                    None => {
                        let slot = self.model.get(self.model.get(u.edge).expect("edge").curve).expect("edge curve");
                        let curve = project::analytic_pcurve(&slot.curve, &surface)
                            .expect("edge lies on the face's plane");
                        (curve, slot.first, slot.last)
                    }
                };
                self.model.set_pcurve_with_range(coedge, curve, first, last)?;
                coedges.push(coedge);
            }
            let wire = Shape::from(self.model.add_wire(&coedges)?);
            let wire = if reversed { Shape { orientation: Orientation::Reversed, ..wire } } else { wire };
            self.model.add(face, wire)?;
        }
        let face = if reversed { Shape { orientation: Orientation::Reversed, ..face } } else { face };
        self.model.add(self.shell, face)?;
        Ok(face)
    }

    fn solid(self) -> Result<Shape> {
        let solid = self.model.new_empty(NodeKind::Solid);
        self.model.add(solid, self.shell)?;
        Ok(solid)
    }
}

/// The solid bounded by `face` alone, on the side of it where `reference` lies (an
/// unbounded solid, as a tool for booleans).
#[doc(alias = "BRepPrimAPI_MakeHalfSpace")]
pub fn half_space(model: &mut Model, face: Shape, reference: Point) -> Result<Shape> {
    use crate::geometry::ParametricSurface;
    use crate::math::Transformable;
    if face.kind() != NodeKind::Face {
        return Err(Error::InvalidInput("a half-space is bounded by a face".into()));
    }
    // A placed face: a copy moved where it is placed (a shell holds faces as they are).
    if face.placement != crate::math::Transform::IDENTITY {
        let unplaced = Shape { placement: crate::math::Transform::IDENTITY, ..face };
        let moved = crate::modify::transform(model, unplaced, &face.placement, Default::default())?;
        return half_space(model, moved, reference);
    }
    let surface = model
        .surface(face)
        .ok_or_else(|| Error::InvalidInput("a half-space needs a face with a surface".into()))?
        .transformed(&face.placement);
    // The face's normal, as used, where the face is nearest the reference point (a curved
    // face can turn round: its middle can face the other way); at the middle of its
    // parameter range if that is off the face.
    let (u1, u2, v1, v2) = model.face_uv_bounds(face);
    let (u, v, _) = crate::geometry::extrema::nearest_on_surface(reference, &surface, None);
    let on_face = (u1..=u2).contains(&u) && (v1..=v2).contains(&v);
    let (u, v) = if on_face { (u, v) } else { (0.5 * (u1 + u2), 0.5 * (v1 + v2)) };
    let (p, su, sv) = surface.d1(u, v);
    let normal = su.cross(sv) * if face.orientation == Orientation::Reversed { -1.0 } else { 1.0 };
    let side = (reference - p).dot(normal);
    if side.abs() <= tolerance::LINEAR {
        return Err(Error::InvalidInput("the reference point lies on the half-space's boundary".into()));
    }
    // The face's normal points out of the material, away from the reference point.
    let shell = model.new_empty(NodeKind::Shell);
    model.add(shell, if side > 0.0 { face.reversed() } else { face })?;
    let solid = model.new_empty(NodeKind::Solid);
    model.add(solid, shell)?;
    Ok(solid)
}

/// The plane through `origin` with normal `normal` and X axis `x` (the parametrization of
/// a primitive's planar face).
fn plane(origin: Point, normal: Vector, x: Vector) -> Surface {
    Surface::Plane(Plane::from_frame(Frame::new(origin, Direction::from_vector(normal), Direction::from_vector(x))))
}

fn check_angle(angle: f64) -> Result<()> {
    if angle > 0.0 && angle <= TAU + tolerance::PARAMETRIC {
        Ok(())
    } else {
        Err(Error::InvalidInput(format!("revolution angle must be in (0, 2π], got {angle}")))
    }
}

#[cfg(test)]
mod tests {
    use std::f64::consts::PI;

    use super::*;

    #[test]
    fn tori_with_partial_tubes() {
        // Pappus: 2π ∫ρ dA over the profile. Tube [0, π/2] (R 20, r 5): the rectangle
        // [0, 20] x [0, 5] and the quarter disk beyond it; tube [π/2, π] (R 20, r 10): the
        // rectangle [0, 20] x [0, 10] less the quarter disk turned towards the axis.
        let quarter = |r: f64, centroid: f64| PI * r * r / 4.0 * centroid;
        let cases = [
            ((20.0, 5.0, 0.0, PI / 2.0), 2.0 * PI * (200.0 * 5.0 + quarter(5.0, 20.0 + 20.0 / (3.0 * PI)))),
            ((20.0, 10.0, PI / 2.0, PI), 2.0 * PI * (200.0 * 10.0 - quarter(10.0, 20.0 - 40.0 / (3.0 * PI)))),
        ];
        for ((big, small, v1, v2), volume) in cases {
            for angle in [2.0 * PI, PI / 2.0] {
                let mut model = Model::new();
                let torus = Torus::new(big, small).with_tube_angles(v1, v2).with_angle(angle).build(&mut model).unwrap();
                let report = crate::validity::check(&model, torus);
                assert!(report.problems.is_empty(), "{report:?}");
                let v = crate::properties::volume(&model, torus).mass;
                let expected = volume * angle / (2.0 * PI);
                assert!((v - expected).abs() < 1e-9 * expected, "tube [{v1}, {v2}], angle {angle}: {v} != {expected}");
            }
        }
    }
}
