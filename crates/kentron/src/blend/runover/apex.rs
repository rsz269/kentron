//! A concave fillet rolling on past its edge's end, where the edge ends at the apex of a
//! horn torus: the corner an earlier fillet's ball swings round in (two fillets of one
//! radius r over a concave edge square to the plane they share, a pocket's rim at its
//! corners, OCCT's `ChFiKPart_MakeRotule` torus). The edge's two faces, planar walls, go on
//! past the apex into the rim's two cylinders (tangent to them along lines through it), the
//! torus between those.
//!
//! OCCT has no counterpart: at such a vertex it mitres the rim (see `toric_corners`), and
//! the concave edge then goes on, tangentially, along the mitre's crease, which its fillet
//! runs on along (a stripe whose spine goes on through the vertex). Its own way past a
//! spine's end, `ChFi3d_FilBuilder::ExtentOneCorner` running the guide on past the vertex
//! and the walking crossing onto faces tangent across the arcs it reaches, would carry the
//! ball on here too, but never meets this vertex, and leaves its end to the corner's
//! intersection with the face it runs into, which here is no face but a ball sitting in the
//! torus. kentron's rim swings round in the torus (where OCCT's mitre leaves out material a
//! rolling ball keeps), and the crease is gone: what the ball of the concave fillet's
//! radius R does there is computed from its definition, the union of the balls in the air
//! not reaching into the material, and checked against a voxel closing of the solid.
//!
//! Past the apex the ball touches the two cylinders: its centre lies at R + r from both
//! axes (lines at r from the walls, in the plane square to the edge through the apex), on
//! the ellipse where the two cylinders of radius R + r about them cross, in the plane
//! bisecting the walls. It rises until it lies on the edge's line, at `√(R² + 2Rr)` past the
//! apex: there it touches the torus all along a circle, its contacts on the cylinders at
//! their ends, and can go no further. The fillet is the envelope of the balls along the
//! ellipse, from the edge's own fillet at the apex up to that last ball, and the last ball's
//! sphere between its section there and the circle on the torus.
//!
//! It is built as the run-over fillets are (see the [module](super) documentation): a wedge
//! from within the balls out through the contacts into the material is cut away (in each
//! section, square to the ellipse, from within the ball through each contact to the
//! cylinder's axis; past the last section, the cone from the last ball's centre through its
//! circle on the torus, turned about the edge's line), the faces of the solid it leaves are
//! sewn to the balls' arcs lofted and the last ball's sphere.

use crate::error::{Error, Result};
use crate::geometry::curve::{Circle, Curve};
use crate::geometry::surface::{Sphere, Surface};
use crate::math::{Axis, Direction, Frame, Point, Transform, Transformable, Vector};
use crate::topology::{Model, NodeKind, Shape};
use crate::{boolean, loft, make, sew, unify};

/// Where a concave fillet's edge ends at the apex of a horn torus (see the
/// [module](self) documentation).
#[derive(Clone, Debug)]
pub(in crate::blend) struct Apex {
    /// The edge's end, the torus's apex.
    pub vertex: Point,
    /// Along the edge, out of it past its end (the torus's axis).
    pub out: Vector,
    /// Each wall's normal into the material.
    pub inward: [Vector; 2],
    /// Along each wall from the edge, square to it (the way the rim's cylinder runs).
    pub along: [Vector; 2],
    /// The rim's radius (the torus's and the cylinders').
    pub rim: f64,
    /// The concave fillet's radius.
    pub radius: f64,
    /// The rim's cylinders on each wall's side and the torus, in the solid before the
    /// fillet.
    pub faces: [Shape; 3],
}

/// A section of the fillet past the apex: the ball's centre, the feet on the cylinders'
/// axes and the directions to them.
#[derive(Clone, Copy)]
struct Section {
    center: Point,
    feet: [Point; 2],
    toward: [Vector; 2],
}

/// The sections along the ellipse, evenly in its angle (from the apex's plane, the edge's
/// own fillet's end, to the last ball on the edge's line).
const SECTIONS: usize = 32;

impl Apex {
    /// The ball's section at angle `psi` of the ellipse (0 at the apex's plane).
    fn section(&self, psi: f64) -> Section {
        let (r, big) = (self.rim, self.radius + self.rim);
        let [d0, d1] = self.inward;
        // (Into the air, bisecting the walls, and the cosine of half the angle between their
        // normals.)
        let m = -(d0 + d1).normalized();
        let cos_half = (0.5 * d0.angle(d1)).cos();
        let s = (big * psi.cos() - r) / cos_half;
        let center = self.vertex + self.out * (big * psi.sin()) + m * s;
        let feet = [0, 1].map(|i| {
            let a = self.vertex + self.inward[i] * r;
            a + self.along[i] * (center - a).dot(self.along[i])
        });
        Section { center, feet, toward: feet.map(|f| (f - center) * (1.0 / big)) }
    }

    /// The ellipse's angle where the ball reaches the edge's line.
    fn last_angle(&self) -> f64 {
        (self.rim / (self.radius + self.rim)).acos()
    }
}

/// `shape` (the solid filleted, the edge's own fillet ending square at the apex) with the
/// fillet rolled on past `apex` (see the [module](self) documentation).
pub(in crate::blend) fn roll_past(model: &mut Model, shape: Shape, apex: &Apex) -> Result<Shape> {
    let radius = apex.radius;
    let last = apex.last_angle();
    let sections: Vec<Section> = (0..=SECTIONS).map(|j| apex.section(last * j as f64 / SECTIONS as f64)).collect();
    // (The ball's contacts on the cylinders lie on their faces, and on the torus on its: the
    // rim is long enough for the ball to roll on it.)
    let on = |model: &Model, face: Shape, p: Point| matches!(boolean::PointClassifier::new(model, face).map(|c| c.locate(p)), Ok(boolean::PointLocation::Boundary));
    for s in &sections[1..SECTIONS] {
        for i in 0..2 {
            if (s.feet[i] - apex.vertex - apex.inward[i] * apex.rim).dot(apex.along[i]) < 0.0 || !on(model, apex.faces[i], s.center + s.toward[i] * radius) {
                return Err(Error::InvalidInput("a fillet rolling on past a rim's torus off the rim's faces".into()));
            }
        }
    }
    let q = sections[SECTIONS].center;
    let turn = Vector::from(Direction::from_vector(apex.inward[0].cross(apex.inward[1])));
    let alpha = apex.inward[0].angle(apex.inward[1]);
    let axis = Axis::new(apex.vertex, Direction::from_vector(turn));
    for k in 1..8 {
        let a = (apex.vertex + apex.inward[0] * apex.rim).transformed(&Transform::rotation(axis, alpha * k as f64 / 8.0));
        if !on(model, apex.faces[2], q + (a - q) * (radius / (radius + apex.rim))) {
            return Err(Error::InvalidInput("a fillet rolling on past a rim's torus off the torus".into()));
        }
    }

    // The wedge through the sections, from within each ball through its contacts to the
    // axes' feet, closed between those.
    let inner = 0.75 * radius;
    let mut wedge_wires = Vec::with_capacity(sections.len());
    let mut arc_wires = Vec::with_capacity(sections.len());
    for s in &sections {
        let points = [s.center + s.toward[0] * inner, s.feet[0], s.feet[1], s.center + s.toward[1] * inner];
        wedge_wires.push(super::polygon(model, &points)?);
        let (a, b) = (s.toward[0], s.toward[1]);
        let circle = Circle::new(Frame::new(s.center, Direction::from_vector(a.cross(b)), Direction::from_vector(a)), radius);
        let edge = make::edge_from_curve_range(model, Curve::Circle(circle), 0.0, a.angle(b))?;
        arc_wires.push(make::wire(model, &[edge])?);
    }
    let spacing = loft::Spacing::Even;
    let wedge = loft::loft_breaking(model, &wedge_wires, &loft::LoftOptions { solid: true, spacing, ..Default::default() }, &[]).map_err(|e| Error::Failed(format!("lofting the wedge past a rim's torus: {e}")))?;
    let wedge = if super::super::volume_of(model, wedge) < 0.0 { wedge.reversed() } else { wedge };
    // Past the last section: the cone from the last ball's centre through its circle on the
    // torus, and on in to the axis and down to the apex, turned about the edge's line from
    // one wall's normal to the other's.
    let g0 = (apex.vertex + apex.inward[0] * apex.rim - q) * (1.0 / (radius + apex.rim));
    let up = q - apex.vertex;
    let start = q + g0 * inner;
    let on_axis = apex.vertex + up * ((start - apex.vertex).dot(up) / up.dot(up));
    let profile = super::polygon(model, &[on_axis, start, apex.vertex + apex.inward[0] * apex.rim, apex.vertex])?;
    let profile = make::planar_face(model, profile)?;
    let cap_wedge = crate::sweep::revolve(model, profile, axis, alpha)?.shape;
    let cut = boolean::difference(model, shape, wedge).map_err(|e| Error::Failed(format!("cutting the wedge past a rim's torus: {e}")))?;
    let cut = boolean::difference(model, cut, cap_wedge).map_err(|e| Error::Failed(format!("cutting the wedge past a rim's torus: {e}")))?;
    // The faces the wedges leave: not those on their boundaries.
    let (in_wedge, in_cap) = (boolean::PointClassifier::new(model, wedge)?, boolean::PointClassifier::new(model, cap_wedge)?);
    let mut faces = Vec::new();
    for f in model.unique_subshapes(cut, NodeKind::Face) {
        let Some(p) = crate::sweep::interior_point(model, f) else {
            return Err(Error::Failed("a face cut past a rim's torus without a point inside".into()));
        };
        if in_wedge.locate(p) != boolean::PointLocation::Boundary && in_cap.locate(p) != boolean::PointLocation::Boundary {
            faces.push(f);
        }
    }
    // The balls' arcs, and the last ball's sphere from its arc to its circle on the torus.
    let arcs = loft::loft_breaking(model, &arc_wires, &loft::LoftOptions { solid: false, spacing, ..Default::default() }, &[]).map_err(|e| Error::Failed(format!("lofting the fillet past a rim's torus: {e}")))?;
    faces.extend(model.unique_subshapes(arcs, NodeKind::Face));
    let last_section = &sections[SECTIONS];
    let cos_phi = up.length() / (radius + apex.rim);
    let small = Circle::new(Frame::new(q - up.normalized() * (radius * cos_phi), Direction::from_vector(turn), Direction::from_vector(apex.inward[0])), radius * (1.0 - cos_phi * cos_phi).sqrt());
    let small = make::edge_from_curve_range(model, Curve::Circle(small), 0.0, alpha)?;
    let (a, b) = (last_section.toward[1], last_section.toward[0]);
    let great = Circle::new(Frame::new(q, Direction::from_vector(a.cross(b)), Direction::from_vector(a)), radius);
    let great = make::edge_from_curve_range(model, Curve::Circle(great), 0.0, a.angle(b))?;
    let rim = make::wire_within(model, &[small, great], 1e-7 * (1.0 + radius))?;
    let sphere = Surface::Sphere(Sphere::new(Frame::new(q, Direction::from_vector(turn), Direction::from_vector(apex.inward[0])), radius));
    let cap = make::face_bounded_on(model, sphere, &[rim])?;
    // (The lens between the two arcs, on one side of the great circle, not the rest of the
    // sphere.)
    if crate::properties::area(model, cap).mass > 2.0 * std::f64::consts::PI * radius * radius {
        return Err(Error::Failed("the last ball's face past a rim's torus came out the rest of its sphere".into()));
    }
    faces.push(cap);
    let sewn = sew::sew(model, &faces, &sew::SewOptions { tolerance: 1e-5 * radius.max(1.0), ..Default::default() })?;
    let shells = model.unique_subshapes(sewn, NodeKind::Shell);
    let [shell] = shells[..] else {
        return Err(Error::Failed(format!("a fillet past a rim's torus whose pieces do not close ({} shells)", shells.len())));
    };
    let solid = super::solid_of(model, shell)?;
    // (Each edge run along once each way; a degenerate one, as the apex of another of the
    // rim's tori, is used by one face.)
    let mut uses: crate::hashing::FastMap<crate::topology::NodeId, (usize, usize)> = Default::default();
    for e in model.explore(shell, NodeKind::Edge) {
        if !model.is_degenerate(e) {
            let used = uses.entry(e.node).or_default();
            if e.orientation == crate::topology::Orientation::Reversed { used.1 += 1 } else { used.0 += 1 }
        }
    }
    if uses.values().any(|&used| used != (1, 1)) {
        return Err(Error::Failed("a fillet past a rim's torus whose pieces leave the shell open".into()));
    }
    let result = unify::unify_same_domain(model, solid)?;
    super::tolerant_of_pcurves(model, result, 1e-5 * (1.0 + radius))?;
    // (Checked: valid, and what it adds no more than the wedges hold, nor less than nothing.)
    let (before, after) = (super::super::volume_of(model, shape), super::super::volume_of(model, result));
    let room = super::super::volume_of(model, wedge) + super::super::volume_of(model, cap_wedge);
    let slack = 1e-6 * before.abs().max(1.0);
    if !crate::validity::check(model, result).is_valid() || after - before < -slack || after - before > room + slack {
        return Err(Error::Failed("a fillet past a rim's torus that does not come out right".into()));
    }
    Ok(result)
}
