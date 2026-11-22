//! Closed-form checks of [`super::setback::fillet_setback`].

use std::f64::consts::PI;

use crate::boolean::{classify_point, PointLocation};
use crate::error::Error;
use crate::math::Point;
use crate::topology::{Model, NodeKind, Shape};

use super::fillet_setback;

fn edges_of(model: &Model, shape: Shape) -> Vec<Shape> {
    model.unique_subshapes(shape, NodeKind::Edge)
}

/// A box of sides `w` with every edge set back at every corner is the box shrunk by `r` with
/// its corners chamfered back `setback - r` along its edges, swelled by `r`: a convex body
/// whose volume Steiner's formula gives.
#[test]
fn whole_box_is_the_chamfered_core_swelled() {
    let (w, r, s) = ([10.0, 12.0, 14.0], 1.0, 2.5);
    let mut model = Model::new();
    let block = crate::primitives::Cuboid::from_corners(Point::new(0.0, 0.0, 0.0), Point::new(w[0], w[1], w[2])).build(&mut model).unwrap();
    let edges = edges_of(&model, block);
    assert_eq!(edges.len(), 12);
    let out = fillet_setback(&mut model, block, &edges, r, s).unwrap();
    let core = w.map(|x| x - 2.0 * r);
    let d = s - r;
    let volume = core[0] * core[1] * core[2] - 8.0 * d.powi(3) / 6.0;
    let area = 2.0 * (core[0] * core[1] + core[1] * core[2] + core[2] * core[0]) - 12.0 * d * d + 4.0 * 3f64.sqrt() * d * d;
    // Each edge's length by the angle its faces' normals turn through.
    let turns = 4.0 * core.iter().map(|x| x - 2.0 * d).sum::<f64>() * PI / 2.0 + 24.0 * d * 2f64.sqrt() * (1.0 / 3f64.sqrt()).acos();
    let want = volume + area * r + turns * r * r / 2.0 + 4.0 * PI * r.powi(3) / 3.0;
    let got = super::volume_of(&model, out);
    assert!((got - want).abs() <= 1e-4 * want, "volume {got} against {want}");
}

/// One corner of a box: where its three fillets' contact lines would cross on a face, the
/// result holds nothing; the flat holds points just beyond the chamfer's trace; and past the
/// setback each edge is the ordinary fillet.
#[test]
fn corner_flat_stops_short() {
    let (w, r, s) = (10.0, 1.0, 3.0);
    let mut model = Model::new();
    let block = crate::primitives::Cuboid::from_corners(Point::new(0.0, 0.0, 0.0), Point::new(w, w, w)).build(&mut model).unwrap();
    let at = Point::new(w, w, w);
    let edges: Vec<Shape> = edges_of(&model, block)
        .into_iter()
        .filter(|&e| {
            let (a, b) = model.edge_vertices(e);
            model.point(a).distance(at) < 1e-9 || model.point(b).distance(at) < 1e-9
        })
        .collect();
    assert_eq!(edges.len(), 3);
    let out = fillet_setback(&mut model, block, &edges, r, s).unwrap();
    let inside = |model: &Model, p: Point| matches!(classify_point(model, out, p).unwrap(), PointLocation::Inside);
    // (Depths in from the corner.)
    let at = |a: f64, b: f64, c: f64| Point::new(w - a, w - b, w - c);
    let e = 1e-4;
    assert!(!inside(&model, at(r, r, e)));
    // The flat on the face square to z reaches the line a + b = s + r.
    assert!(inside(&model, at(0.5 * (s + r) + 0.05, 0.5 * (s + r) + 0.05, e)));
    assert!(!inside(&model, at(0.5 * (s + r) - 0.05, 0.5 * (s + r) - 0.05, e)));
    // Past the setback, the round of radius r about its axis r in from both faces.
    let k = 1.0 - 0.5f64.sqrt();
    assert!(inside(&model, at(s + 0.5, r * k + 0.01, r * k + 0.01)));
    assert!(!inside(&model, at(s + 0.5, r * k - 0.01, r * k - 0.01)));
    // The flat of the corner: the core's cut a + b + c = s + 2r swelled by r.
    let m = (s + 2.0 * r - 3f64.sqrt() * r) / 3.0;
    assert!(inside(&model, at(m + 0.01, m + 0.01, m + 0.01)));
    assert!(!inside(&model, at(m - 0.01, m - 0.01, m - 0.01)));
}

#[test]
fn setback_must_pass_the_radius() {
    let mut model = Model::new();
    let block = crate::primitives::Cuboid::from_corners(Point::new(0.0, 0.0, 0.0), Point::new(10.0, 10.0, 10.0)).build(&mut model).unwrap();
    let edges = edges_of(&model, block);
    assert!(matches!(fillet_setback(&mut model, block, &edges, 1.0, 1.0), Err(Error::InvalidInput(_))));
    assert!(matches!(fillet_setback(&mut model, block, &edges, 0.0, 2.0), Err(Error::InvalidInput(_))));
    // Two stops of 6 + (2 - √3) on a side of 10 do not fit.
    assert!(matches!(fillet_setback(&mut model, block, &edges, 1.0, 6.0), Err(Error::InvalidInput(_))));
}

