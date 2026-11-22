//! Rounds named by faces rather than by an edge: a ball rolling between two faces that need
//! not meet, and a full round taking the place of a face between two others.
//!
//! OCCT's rolling ball between faces is `BiTgte_Blend` (DRAW's `rollingball`), built on the
//! offset algorithm's machinery (`BRepOffset_Offset`, `BRepOffset_Inter3d`/`Inter2d`,
//! `BRepOffset_MakeLoops`), none of which kentron has. These are instead built from the
//! run-over fillet's: the ball rolled section by section between the two faces (as bounded),
//! the wedge beyond the chord between its contacts cut away and the lens between the chord
//! and its arc kept within the solid (see the [module](super) documentation). For a full round
//! OCCT has no call at all; the round is the circle touching all three faces, by its
//! definition.
//!
//! For now only planar faces: the line the ball runs along is then where their planes meet,
//! and (for a full round) the circle touching three planes along one direction is the same
//! all along it.

use super::*;

/// `face` as `shape` holds it: with the orientation it has in its shell.
fn held(model: &Model, shape: Shape, face: Shape) -> Result<Shape> {
    model
        .explore(shape, NodeKind::Face)
        .find(|f| f.node == face.node)
        .ok_or_else(|| Error::InvalidInput("the face is not one of the solid's".into()))
}

/// The points of `face`'s vertices.
fn corners(model: &Model, face: Shape) -> Vec<Point> {
    model.unique_subshapes(face, NodeKind::Vertex).iter().map(|&v| model.point(v)).collect()
}

/// A planar face's point (its first vertex) and its normal out of the solid.
fn plane_of(model: &Model, face: Shape) -> Result<(Point, Vector)> {
    let planar = model.surface(face).is_some_and(|s| matches!(s.transformed(&face.placement), Surface::Plane(_)));
    let at = corners(model, face).first().copied();
    match (planar, at) {
        (true, Some(p)) => outward_normal(model, face, p).map(|n| (p, n)).ok_or_else(|| Error::Failed("a planar face without a normal".into())),
        _ => Err(Error::Failed("rounds named by faces other than planes are not supported yet".into())),
    }
}

/// The range of `points` along `d`.
fn range(points: &[Point], d: Vector) -> (f64, f64) {
    points.iter().map(|p| (*p - Point::ORIGIN).dot(d)).fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), t| (lo.min(t), hi.max(t)))
}

/// The point `x` with `n[i]·x = h[i]` for the three (independent) normals `n`.
fn meet(n: [Vector; 3], h: [f64; 3]) -> Option<Point> {
    let det = n[0].dot(n[1].cross(n[2]));
    if !(det.abs() > 1e-12) {
        return None;
    }
    let v = (n[1].cross(n[2]) * h[0] + n[2].cross(n[0]) * h[1] + n[0].cross(n[1]) * h[2]) * (1.0 / det);
    Some(Point::ORIGIN + v)
}

/// Whether the two faces share an edge.
fn adjacent(model: &Model, a: Shape, b: Shape) -> bool {
    let edges = model.unique_subshapes(a, NodeKind::Edge);
    model.unique_subshapes(b, NodeKind::Edge).iter().any(|e| edges.iter().any(|f| f.node == e.node))
}

/// `shape` rounded between `face_a` and `face_b` by a ball of `radius` rolling touching both:
/// what lies beyond its path between the contacts, faces between them included, is cut away.
pub(in crate::blend) fn fillet_faces(model: &mut Model, shape: Shape, face_a: Shape, face_b: Shape, radius: f64) -> Result<Shape> {
    if !(radius > 0.0 && radius.is_finite()) {
        return Err(Error::InvalidInput(format!("the radius must be positive, got {radius}")));
    }
    if face_a.node == face_b.node {
        return Err(Error::InvalidInput("the two faces must differ".into()));
    }
    let faces = [held(model, shape, face_a)?, held(model, shape, face_b)?];
    let [(pa, na), (pb, nb)] = [plane_of(model, faces[0])?, plane_of(model, faces[1])?];
    // The ball runs along where the planes meet.
    let cross = na.cross(nb);
    if cross.length() <= 1e-9 {
        return Err(Error::Failed("a round between parallel faces is not supported yet".into()));
    }
    let d = cross.normalized();
    let origin = meet([na, nb, d], [(pa - Point::ORIGIN).dot(na), (pb - Point::ORIGIN).dot(nb), 0.0]).ok_or_else(|| Error::Failed("the faces' planes do not meet".into()))?;
    // Along as much of it as both faces lie beside.
    let (ra, rb) = (range(&corners(model, faces[0]), d), range(&corners(model, faces[1]), d));
    let (t0, t1) = (ra.0.max(rb.0), ra.1.min(rb.1));
    if !(t1 - t0 > 1e-7 * (1.0 + radius)) {
        return Err(Error::InvalidInput("the two faces lie side by side nowhere".into()));
    }
    let line = Curve::Line(crate::geometry::curve::Line::new(origin, Direction::from_vector(d)));
    let spine = Spine::new(vec![(line, t0, t1)], false);
    // The direction of the line in each face's boundary, were it an edge of both: +1 where the
    // face lies to its left (seen from outside).
    let along = |model: &Model, face: Shape, n: Vector| -> Result<f64> {
        let points = corners(model, face);
        let middle = points.iter().fold(Vector::new(0.0, 0.0, 0.0), |s, p| s + (*p - Point::ORIGIN)) * (1.0 / points.len() as f64);
        let side = n.cross(d).dot(Point::ORIGIN + middle - origin);
        if side.abs() <= 1e-9 * (1.0 + radius) {
            return Err(Error::InvalidInput("a face lies across the line its plane meets the other's on".into()));
        }
        Ok(side.signum())
    };
    let alongs = [along(model, faces[0], na)?, along(model, faces[1], nb)?];
    let side = |model: &Model, f: Shape| Bounded::of(model, f).map(|b| Side(vec![b])).ok_or_else(|| Error::Failed("a face without a boundary".into()));
    let side_set = [side(model, faces[0])?, side(model, faces[1])?];
    let rolled = roll_on(model, &spine, side_set, [(faces[0], alongs[0]), (faces[1], alongs[1])], &|_| radius, 48, Switches::Bisected)
        .map_err(|_| Error::InvalidInput("the ball does not fit between the two faces".into()))?;
    if !rolled.convex {
        return Err(Error::Failed("a round adding material between two faces is not supported yet".into()));
    }
    // (The ball on the two faces themselves, not over their rims.)
    if rolled.supports.iter().any(|s| s.iter().any(|&on| on != (0, None))) {
        return Err(Error::InvalidInput("the ball does not fit between the two faces".into()));
    }
    by_booleans(model, shape, rolled, radius, &Drawing { wedge: None, straight: true })
}

/// `shape` with `middle` replaced by a round tangent to it and to `side_a` and `side_b`: the
/// circle touching all three faces, kept from one side's contact over the middle to the other's.
pub(in crate::blend) fn full_round(model: &mut Model, shape: Shape, side_a: Shape, middle: Shape, side_b: Shape) -> Result<Shape> {
    if side_a.node == middle.node || side_b.node == middle.node || side_a.node == side_b.node {
        return Err(Error::InvalidInput("a full round needs three different faces".into()));
    }
    let [fa, fm, fb] = [held(model, shape, side_a)?, held(model, shape, middle)?, held(model, shape, side_b)?];
    if !adjacent(model, fa, fm) || !adjacent(model, fb, fm) {
        return Err(Error::InvalidInput("the middle face must share an edge with each side face".into()));
    }
    let [(pa, na), (pm, nm), (pb, nb)] = [plane_of(model, fa)?, plane_of(model, fm)?, plane_of(model, fb)?];
    // Along the rib, square to all three normals.
    let cross = na.cross(nm);
    if cross.length() <= 1e-9 {
        return Err(Error::Failed("a full round with a side tangent to the middle face".into()));
    }
    let d = cross.normalized();
    if nb.dot(d).abs() > 1e-9 {
        return Err(Error::Failed("a full round whose circle changes along the rib is not supported yet".into()));
    }
    let h = |p: Point, n: Vector| (p - Point::ORIGIN).dot(n);
    // The centre `c` and radius `r` of the circle touching the three planes from within (in the
    // plane square to `d` through the origin): n·c + r = n·p for each, in coordinates
    // `c = x·e1 + y·e2`.
    let e1 = nm;
    let e2 = d.cross(e1);
    let rows = [na, nm, nb].map(|n| Vector::new(n.dot(e1), n.dot(e2), 1.0));
    let solved = meet(rows, [h(pa, na), h(pm, nm), h(pb, nb)]).ok_or_else(|| Error::InvalidInput("no circle touches the three faces".into()))?;
    let (x, y, r) = (solved.x, solved.y, solved.z);
    if !(r > 1e-9 * (1.0 + x.abs() + y.abs())) {
        return Err(Error::InvalidInput("no circle touches the three faces from within the solid".into()));
    }
    let center0 = Point::ORIGIN + e1 * x + e2 * y;
    // Along the middle face, as far as the sides go too.
    let ranges = [range(&corners(model, fm), d), range(&corners(model, fa), d), range(&corners(model, fb), d)];
    let (t0, t1) = (ranges.iter().map(|r| r.0).fold(f64::NEG_INFINITY, f64::max), ranges.iter().map(|r| r.1).fold(f64::INFINITY, f64::min));
    if !(t1 - t0 > 1e-7 * (1.0 + r)) {
        return Err(Error::InvalidInput("the three faces lie side by side nowhere".into()));
    }
    let section = |t: f64| {
        let center = center0 + d * t;
        Section { t, at: center + nm * r, center, contacts: [center + na * r, center + nb * r], uv: [None, None], radius: r }
    };
    const COUNT: usize = 4;
    let sections: Vec<Section> = (0..=COUNT).map(|i| section(t0 + (t1 - t0) * i as f64 / COUNT as f64)).collect();
    // The circle touching each face on it, not past its rim (a side shorter than the radius).
    let bounded = [fa, fm, fb].map(|f| Bounded::of(model, f));
    let [Some(ba), Some(bm), Some(bb)] = bounded else {
        return Err(Error::Failed("a face without a boundary".into()));
    };
    let near = 1e-7 * (1.0 + r);
    for s in &sections {
        for (b, p) in [(&ba, s.contacts[0]), (&bm, s.at), (&bb, s.contacts[1])] {
            if b.nearest(p).0 > near {
                return Err(Error::InvalidInput("the circle touching the three faces touches one off it".into()));
            }
        }
    }
    // The wedge: from each contact out of the solid, along its side past the middle face, and
    // over it to the other side.
    let margin = 0.5 * r;
    let wedge = move |s: &Section| -> Vec<Point> {
        let level = (s.at - Point::ORIGIN).dot(d);
        let corner = |p: Point, n: Vector| meet([n, nm, d], [h(p, n), h(pm, nm), level]);
        let (Some(qa), Some(qb)) = (corner(pa, na), corner(pb, nb)) else { return Vec::new() };
        let [ca, cb] = s.contacts;
        vec![ca, ca + na * margin, qa + (na + nm) * margin, qb + (nb + nm) * margin, cb + nb * margin, cb]
    };
    let rolled = Rolled {
        supports: vec![[(0, None); 2]; sections.len()],
        sections,
        breaks: Vec::new(),
        side_set: [Side(vec![ba]), Side(vec![bb])],
        convex: true,
        closed: false,
        tangents: [d, d],
    };
    by_booleans(model, shape, rolled, r, &Drawing { wedge: Some(&wedge), straight: true })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::primitives::Cuboid;
    use std::f64::consts::PI;

    /// The one face of `shape` through `p`.
    fn face_at(model: &mut Model, shape: Shape, p: [f64; 3]) -> Shape {
        let vertex = make::vertex(model, Point::new(p[0], p[1], p[2]));
        let found: Vec<Shape> = model
            .unique_subshapes(shape, NodeKind::Face)
            .into_iter()
            .filter(|&f| crate::distance::distance(model, vertex, f).unwrap().value <= 1e-9)
            .collect();
        assert_eq!(found.len(), 1, "faces through {p:?}");
        found[0]
    }

    fn volume(model: &Model, s: Shape) -> f64 {
        crate::properties::volume(model, s).mass
    }

    fn area(model: &Model, s: Shape) -> f64 {
        crate::properties::area(model, s).mass
    }

    fn normal(model: &mut Model, s: Shape, p: [f64; 3]) -> Option<Direction> {
        crate::properties::normals(model, s, &[Point::new(p[0], p[1], p[2])]).unwrap()[0]
    }

    fn assert_normal(model: &mut Model, s: Shape, p: [f64; 3], n: Vector) {
        let got = normal(model, s, p).unwrap_or_else(|| panic!("no normal at {p:?}"));
        assert!((Vector::from(got) - n.normalized()).length() < 1e-9, "{got:?} at {p:?}");
    }

    fn the_box(model: &mut Model) -> Shape {
        Cuboid::new(10.0, 20.0, 30.0).build(model).unwrap()
    }

    /// Between the bottom and the side x = 0 of a 10 × 20 × 30 box (which share the edge along
    /// y), a ball of radius 2 takes off the corner less the quarter disc, (4 − π)·20, and the
    /// strips of width 2 on both faces for the quarter cylinder; both ends lose the corner too.
    fn assert_rounded_at_the_corner(model: &mut Model, rounded: Shape) {
        let r: f64 = 2.0;
        let corner = r * r - PI * r * r / 4.0;
        assert!((volume(model, rounded) - (6000.0 - corner * 20.0)).abs() < 1e-6, "{}", volume(model, rounded));
        let area_before = 2.0 * (10.0 * 20.0 + 10.0 * 30.0 + 20.0 * 30.0);
        let want = area_before - 2.0 * r * 20.0 + PI * r / 2.0 * 20.0 - 2.0 * corner;
        assert!((area(model, rounded) - want).abs() < 1e-6, "{} against {want}", area(model, rounded));
        assert!(crate::validity::check(model, rounded).is_valid());
        // Smooth where the round meets each face, the normal that face's; halfway round, the
        // normal from the axis.
        assert_normal(model, rounded, [2.0, 10.0, 0.0], Vector::new(0.0, 0.0, -1.0));
        assert_normal(model, rounded, [0.0, 10.0, 2.0], Vector::new(-1.0, 0.0, 0.0));
        let h = 2.0 - 2.0 * (PI / 4.0).cos();
        assert_normal(model, rounded, [h, 7.0, h], Vector::new(-1.0, 0.0, -1.0));
    }

    #[test]
    fn faces_sharing_an_edge_are_rounded_as_the_edge_is() {
        let mut model = Model::new();
        let b = the_box(&mut model);
        let (bottom, side) = (face_at(&mut model, b, [5.0, 10.0, 0.0]), face_at(&mut model, b, [0.0, 10.0, 15.0]));
        let rounded = fillet_faces(&mut model, b, bottom, side, 2.0).unwrap();
        assert_rounded_at_the_corner(&mut model, rounded);
    }

    /// A chamfer 0.5 back from the corner lies within the corner the round of radius 2 takes
    /// off (the chamfer passes 0.5·sin 45° from the corner, the round 2·(√2 − 1)): rounding
    /// between the faces either side of it leaves what the sharp corner's round does.
    #[test]
    fn a_chamfer_narrower_than_the_round_is_consumed() {
        let mut model = Model::new();
        let b = the_box(&mut model);
        let edge = model
            .unique_subshapes(b, NodeKind::Edge)
            .into_iter()
            .find(|&e| {
                let (c, t0, t1) = model.edge_curve(e).unwrap();
                c.point(0.5 * (t0 + t1)).distance(Point::new(0.0, 10.0, 0.0)) < 1e-9
            })
            .unwrap();
        let chamfered = crate::blend::chamfer(&mut model, b, &[edge], 0.5).unwrap();
        assert_eq!(model.unique_subshapes(chamfered, NodeKind::Face).len(), 7);
        let (bottom, side) = (face_at(&mut model, chamfered, [5.0, 10.0, 0.0]), face_at(&mut model, chamfered, [0.0, 10.0, 15.0]));
        let rounded = fillet_faces(&mut model, chamfered, bottom, side, 2.0).unwrap();
        assert_eq!(model.unique_subshapes(rounded, NodeKind::Face).len(), 7);
        assert_rounded_at_the_corner(&mut model, rounded);
    }

    #[test]
    fn a_round_between_faces_is_refused_where_it_cannot_be() {
        let mut model = Model::new();
        let b = the_box(&mut model);
        let (bottom, side, top) = (face_at(&mut model, b, [5.0, 10.0, 0.0]), face_at(&mut model, b, [0.0, 10.0, 15.0]), face_at(&mut model, b, [5.0, 10.0, 30.0]));
        for radius in [0.0, -1.0, f64::NAN] {
            assert!(matches!(fillet_faces(&mut model, b, bottom, side, radius), Err(Error::InvalidInput(_))));
        }
        assert!(matches!(fillet_faces(&mut model, b, bottom, bottom, 1.0), Err(Error::InvalidInput(_))));
        // A ball wider than the bottom face.
        assert!(matches!(fillet_faces(&mut model, b, bottom, side, 40.0), Err(Error::InvalidInput(_))));
        // Parallel faces: no one line for the ball.
        assert!(fillet_faces(&mut model, b, bottom, top, 1.0).is_err());
    }

    /// The face y = 0 of a 10 × 20 × 30 box between its sides x = 0 and x = 10, 10 apart: a
    /// half cylinder of radius 5 about the line x = 5, y = 5, taking off two corners of
    /// 25 − 25π/4 each along the 30 of it.
    #[test]
    fn a_full_round_between_parallel_sides_is_a_half_cylinder() {
        let mut model = Model::new();
        let b = the_box(&mut model);
        let (a, m, c) = (face_at(&mut model, b, [0.0, 10.0, 15.0]), face_at(&mut model, b, [5.0, 0.0, 15.0]), face_at(&mut model, b, [10.0, 10.0, 15.0]));
        let rounded = full_round(&mut model, b, a, m, c).unwrap();
        let corners = 2.0 * (25.0 - 25.0 * PI / 4.0);
        assert!((volume(&model, rounded) - (6000.0 - corners * 30.0)).abs() < 1e-6, "{}", volume(&model, rounded));
        let area_before = 2.0 * (10.0 * 20.0 + 10.0 * 30.0 + 20.0 * 30.0);
        let want = area_before - (10.0 + 2.0 * 5.0) * 30.0 + 5.0 * PI * 30.0 - 2.0 * corners;
        assert!((area(&model, rounded) - want).abs() < 1e-6, "{} against {want}", area(&model, rounded));
        assert!(crate::validity::check(&model, rounded).is_valid());
        // Smooth where it meets the sides; where the middle face was, the cylinder's lowest
        // line.
        assert_normal(&mut model, rounded, [0.0, 5.0, 15.0], Vector::new(-1.0, 0.0, 0.0));
        assert_normal(&mut model, rounded, [10.0, 5.0, 15.0], Vector::new(1.0, 0.0, 0.0));
        assert_normal(&mut model, rounded, [5.0, 0.0, 15.0], Vector::new(0.0, -1.0, 0.0));
    }

    #[test]
    fn a_full_round_is_refused_where_the_faces_do_not_make_one() {
        let mut model = Model::new();
        let b = the_box(&mut model);
        let (a, m, c) = (face_at(&mut model, b, [0.0, 10.0, 15.0]), face_at(&mut model, b, [5.0, 0.0, 15.0]), face_at(&mut model, b, [10.0, 10.0, 15.0]));
        assert!(matches!(full_round(&mut model, b, a, a, c), Err(Error::InvalidInput(_))));
        // The middle face sharing no edge with a side.
        assert!(matches!(full_round(&mut model, b, m, a, c), Err(Error::InvalidInput(_))));
        // Sides 10 apart and only 4 high: the circle would touch them off their faces.
        let low = Cuboid::new(10.0, 4.0, 30.0).build(&mut model).unwrap();
        let (a, m, c) = (face_at(&mut model, low, [0.0, 2.0, 15.0]), face_at(&mut model, low, [5.0, 0.0, 15.0]), face_at(&mut model, low, [10.0, 2.0, 15.0]));
        assert!(matches!(full_round(&mut model, low, a, m, c), Err(Error::InvalidInput(_))));
    }
}
