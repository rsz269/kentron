//! Drafted extrusions: a planar face extruded with its sides leaning in (or out) by an angle,
//! its rounded corners narrowing to points where the draft closes them.

use std::f64::consts::FRAC_PI_2;

use crate::error::{Error, Result};
use crate::geometry::curve::{Circle, Curve};
use crate::geometry::ParametricCurve;
use crate::make;
use crate::math::{Axis, Direction, Point, Transform, Transformable, Vector};
use crate::properties;
use crate::sew::{self, SewOptions};
use crate::tolerance::{self, ANGULAR};
use crate::topology::{FaceId, Model, NodeKind, Orientation, Shape, VertexId};

/// What an edge of the profile is.
enum Kind {
    Line,
    /// A circular arc about `center` (its circle and parameter range, as the edge has them),
    /// turning by `turn` about the extrusion's direction from start to end; at the top its
    /// radius is `top_radius`, at most zero where the draft closes it.
    Arc { circle: Circle, first: f64, last: f64, center: Point, radius: f64, top_radius: f64, turn: f64 },
}

/// An edge of the profile in the order of its outline, with its ends and tangents that way.
struct Piece {
    start: Point,
    end: Point,
    tangent_start: Vector,
    tangent_end: Vector,
    kind: Kind,
}

impl Piece {
    fn vanishes(&self) -> bool {
        matches!(self.kind, Kind::Arc { top_radius, .. } if top_radius <= tolerance::LINEAR)
    }
}

/// The planar face `profile` swept along `v` (square to its plane) with its sides leaning
/// in by `angle` (out for a negative one), as OCCT's drafted prism: straight edges sweep
/// drafted planes and arcs cones. An arc the draft closes narrows to a point, the cone's
/// apex, above which the planes on either side meet along a corner line; the edges around
/// such an arc must be straight. Only faces without holes bounded by lines and arcs.
#[doc(alias = "LocOpe_DPrism")]
#[doc(alias = "BRepOffsetAPI_DraftAngle")]
pub fn extrude_drafted(model: &mut Model, profile: Shape, v: Vector, angle: f64) -> Result<Shape> {
    if profile.placement != Transform::IDENTITY {
        // Sweep the face where it is defined, then place the result like the profile.
        let local = v.transformed(&profile.placement.inverse());
        let solid = extrude_drafted(model, Shape { placement: Transform::IDENTITY, ..profile }, local, angle)?;
        return Ok(solid.moved(&profile.placement));
    }
    let face = FaceId::try_from(profile.node).map_err(|n| Error::InvalidInput(format!("{n:?} is not a face")))?;
    let height = v.length();
    let z = Direction::try_from_vector(v)
        .filter(|_| height > tolerance::LINEAR && height.is_finite())
        .ok_or_else(|| Error::InvalidInput(format!("extrusion vector must be non-zero, got {v:?}")))?;
    if !(angle.abs() < FRAC_PI_2 - ANGULAR) {
        return Err(Error::InvalidInput(format!("draft angle must be less than a right angle, got {angle}")));
    }
    if model.children(profile).count() != 1 {
        return Err(Error::InvalidInput("a drafted extrusion of a face with holes".into()));
    }
    let wire = model.outer_wire(face).ok_or_else(|| Error::InvalidInput("a face without a boundary".into()))?;
    let z_vector = Vector::from(z);
    let d = height * angle.tan();
    let point = |model: &Model, s: Shape| model[VertexId::from_index(s.node.index)].point;

    // The outline's pieces, and its sense about `z`.
    let edges = model.ordered_edges(Shape::from(wire));
    let mut pieces = Vec::with_capacity(edges.len());
    let mut outline = Vec::new();
    for &edge in &edges {
        let (a, b) = model.edge_vertices(edge);
        let (start, end) = (point(model, a), point(model, b));
        let (curve, first, last) =
            model.edge_curve(edge).ok_or_else(|| Error::InvalidInput("a profile edge without a curve".into()))?;
        let reversed = edge.orientation == Orientation::Reversed;
        let tangent_at = |t: f64| {
            let d1 = curve.d1(t).1.normalized();
            if reversed { -d1 } else { d1 }
        };
        let (tangent_start, tangent_end) =
            if reversed { (tangent_at(last), tangent_at(first)) } else { (tangent_at(first), tangent_at(last)) };
        let kind = match curve {
            Curve::Line(_) => Kind::Line,
            Curve::Circle(circle) if circle.frame().z().cross(z).length() <= ANGULAR => {
                let center = circle.frame().origin();
                let radius = circle.radius();
                let sense = z_vector.dot((start - center).cross(tangent_start)).signum();
                outline.push(curve.point(0.5 * (first + last)));
                Kind::Arc { circle, first, last, center, radius, top_radius: radius, turn: sense * (last - first) }
            }
            _ => return Err(Error::InvalidInput("a drafted extrusion of edges other than lines and arcs square to it".into())),
        };
        outline.push(start);
        pieces.push(Piece { start, end, tangent_start, tangent_end, kind });
    }
    let n = pieces.len();
    if n < 2 || pieces[n - 1].end.distance(pieces[0].start) > tolerance::LINEAR {
        return Err(Error::InvalidInput("a drafted extrusion of an open outline".into()));
    }
    let origin = pieces[0].start;
    if outline.iter().any(|&p| (p - origin).dot(z_vector).abs() > tolerance::LINEAR) {
        return Err(Error::InvalidInput("the extrusion is not square to the profile".into()));
    }
    let area: f64 = (0..outline.len())
        .map(|i| (outline[i] - origin).cross(outline[(i + 1) % outline.len()] - origin).dot(z_vector))
        .sum();
    // Inward from the outline, square to its tangent `t`.
    let sense = area.signum();
    let inward = |t: Vector| z_vector.cross(t) * sense;
    for piece in &mut pieces {
        let start = piece.start;
        let leans = inward(piece.tangent_start);
        if let Kind::Arc { center, radius, top_radius, .. } = &mut piece.kind {
            *top_radius = if leans.dot(*center - start) > 0.0 { *radius - d } else { *radius + d };
        }
    }

    // The top of each joint (piece i's start), in the profile's plane.
    let offset_line = |piece: &Piece| (piece.start + inward(piece.tangent_start) * d, piece.tangent_start);
    let meet = |(p1, t1): (Point, Vector), (p2, t2): (Point, Vector)| -> Result<Point> {
        let across = t1.cross(t2).dot(z_vector);
        if across.abs() <= ANGULAR {
            return Err(Error::InvalidInput("a drafted extrusion's corner between edges in line".into()));
        }
        Ok(p1 + t1 * ((p2 - p1).cross(t2).dot(z_vector) / across))
    };
    let line_beside = |i: usize| -> Result<&Piece> {
        let piece = &pieces[i % n];
        match piece.kind {
            Kind::Line => Ok(piece),
            _ => Err(Error::InvalidInput("an arc a draft closes between edges not straight".into())),
        }
    };
    let mut tops = Vec::with_capacity(n);
    for i in 0..n {
        let (before, after) = (&pieces[(i + n - 1) % n], &pieces[i]);
        let top = if after.vanishes() {
            meet(offset_line(line_beside(i + n - 1)?), offset_line(line_beside(i + 1)?))?
        } else if before.vanishes() {
            meet(offset_line(line_beside(i + n - 2)?), offset_line(line_beside(i)?))?
        } else if before.tangent_end.dot(after.tangent_start) >= 1.0 - ANGULAR {
            after.start + inward(after.tangent_start) * d
        } else if matches!((&before.kind, &after.kind), (Kind::Line, Kind::Line)) {
            meet(offset_line(before), offset_line(after))?
        } else {
            return Err(Error::InvalidInput("a drafted extrusion's sharp corner at an arc".into()));
        };
        tops.push(top + v);
    }
    let apex = |piece: &Piece| match piece.kind {
        Kind::Arc { center, radius, .. } => center + z_vector * (radius * height / d.abs()),
        Kind::Line => unreachable!(),
    };

    // The faces: the profile, the sides and the top.
    let mut faces = vec![profile];
    let mut top_edges = Vec::new();
    for i in 0..n {
        let piece = &pieces[i];
        let (top_start, top_end) = (tops[i], tops[(i + 1) % n]);
        match piece.kind {
            Kind::Line => {
                if (top_end - top_start).dot(piece.tangent_start) <= tolerance::LINEAR {
                    return Err(Error::InvalidInput("the draft closes a straight edge of the profile".into()));
                }
                let mut corners = vec![piece.start, piece.end];
                if pieces[(i + 1) % n].vanishes() {
                    corners.push(apex(&pieces[(i + 1) % n]));
                }
                corners.extend([top_end, top_start]);
                if pieces[(i + n - 1) % n].vanishes() {
                    corners.push(apex(&pieces[(i + n - 1) % n]));
                }
                corners.dedup_by(|a, b| a.distance(*b) <= tolerance::LINEAR);
                let vertices: Vec<Shape> = corners.iter().map(|&p| make::vertex(model, p)).collect();
                let sides = (0..vertices.len())
                    .map(|k| make::edge(model, vertices[k], vertices[(k + 1) % vertices.len()]))
                    .collect::<Result<Vec<_>>>()?;
                let wire = make::wire(model, &sides)?;
                faces.push(make::planar_face(model, wire)?);
                top_edges.push(make::edge_between(model, top_start, top_end)?);
            }
            Kind::Arc { ref circle, first, last, center, top_radius, turn, .. } => {
                let vanishes = piece.vanishes();
                let generator_end = if vanishes { apex(piece) } else { top_start };
                let generator = make::edge_between(model, piece.start, generator_end)?;
                let cone = super::revolve(model, generator, Axis::new(center, z), turn)?.shape;
                faces.extend(model.explore(cone, NodeKind::Face).collect::<Vec<_>>());
                if !vanishes {
                    let lifted = circle.transformed(&Transform::translation(v));
                    let top = Circle::new(*lifted.frame(), top_radius);
                    top_edges.push(make::edge_from_curve_range(model, Curve::Circle(top), first, last)?);
                }
            }
        }
    }
    let top_wire = make::wire(model, &top_edges)?;
    faces.push(make::planar_face(model, top_wire)?);

    let size = pieces.iter().map(|p| p.start.distance(origin)).fold(height, f64::max);
    let sewn = sew::sew(model, &faces, &SewOptions { tolerance: 1e-9 * size + tolerance::LINEAR, ..Default::default() })?;
    let shells = model.unique_subshapes(sewn, NodeKind::Shell);
    let [shell] = shells[..] else {
        return Err(Error::Failed(format!("a drafted extrusion's faces sewn into {} shells, not one", shells.len())));
    };
    let solid = model.new_empty(NodeKind::Solid);
    model.add(solid, shell)?;
    if properties::volume(model, solid).mass < 0.0 {
        let turned = model.new_empty(NodeKind::Solid);
        model.add(turned, shell.reversed())?;
        return Ok(turned);
    }
    Ok(solid)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::Frame;
    use crate::validity;

    /// A rounded square of side `side` with corners of `radius`, about the origin in z = 0.
    fn rounded_square(model: &mut Model, side: f64, radius: f64) -> Shape {
        let h = side / 2.0;
        let mut edges = Vec::new();
        for k in 0..4 {
            let turn = Transform::rotation(Axis::new(Point::ORIGIN, Direction::new(0.0, 0.0, 1.0)), k as f64 * FRAC_PI_2);
            // The side along +x at y = -h, then the corner round to the side along +y.
            let a = Point::new(-h + radius, -h, 0.0).transformed(&turn);
            let b = Point::new(h - radius, -h, 0.0).transformed(&turn);
            edges.push(make::edge_between(model, a, b).unwrap());
            let center = Point::new(h - radius, -h + radius, 0.0);
            let frame = Frame::new(center, Direction::new(0.0, 0.0, 1.0), Direction::new(0.0, -1.0, 0.0)).transformed(&turn);
            let corner = Curve::Circle(Circle::new(frame, radius));
            edges.push(make::edge_from_curve_range(model, corner, 0.0, FRAC_PI_2).unwrap());
        }
        let wire = make::wire(model, &edges).unwrap();
        make::planar_face(model, wire).unwrap()
    }

    #[test]
    fn a_rounded_square_drafted_until_its_corners_close() {
        // As OCCT's drafted prism: 10 faces, the corners' cones ending in apexes below the
        // top.
        let mut model = Model::new();
        let square = rounded_square(&mut model, 41.5, 0.605);
        let solid = extrude_drafted(&mut model, square, Vector::new(0.0, 0.0, 1.13), 45f64.to_radians()).unwrap();
        let report = validity::check(&model, solid);
        assert!(report.is_valid(), "{report:?}");
        assert_eq!(model.explore(solid, NodeKind::Face).count(), 10);
        let volume = properties::volume(&model, solid).mass;
        // The square frustum less its rounded corners: (4 − π)(r − t)² across at height t,
        // until they close at t = r.
        let (s, r, h) = (41.5f64, 0.605f64, 1.13f64);
        let frustum = h / 3.0 * (s * s + s * (s - 2.0 * h) + (s - 2.0 * h).powi(2));
        let expected = frustum - (4.0 - std::f64::consts::PI) * r.powi(3) / 3.0;
        assert!((volume - expected).abs() < 1e-6 * expected, "{volume} vs {expected}");
    }

    #[test]
    fn a_rounded_square_drafted_a_little_keeps_its_corners_round() {
        let mut model = Model::new();
        let square = rounded_square(&mut model, 10.0, 2.0);
        let solid = extrude_drafted(&mut model, square, Vector::new(0.0, 0.0, 1.0), 0.5).unwrap();
        let report = validity::check(&model, solid);
        assert!(report.is_valid(), "{report:?}");
        assert_eq!(model.explore(solid, NodeKind::Face).count(), 10);
        // Each slice is the rounded square offset in: area s² − (4 − π) r².
        let (d, s, r) = (0.5f64.tan(), 10.0f64, 2.0f64);
        let area = |t: f64| (s - 2.0 * d * t).powi(2) - (4.0 - std::f64::consts::PI) * (r - d * t).powi(2);
        let expected = (area(0.0) + 4.0 * area(0.5) + area(1.0)) / 6.0;
        let volume = properties::volume(&model, solid).mass;
        assert!((volume - expected).abs() < 1e-6 * expected, "{volume} vs {expected}");
    }
}
