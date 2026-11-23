//! Guide preparation against OCCT 8.1's (`sweep/tests/guide_occt.txt`): reduced arc
//! lengths, closed guides opened at an origin, open guides turned to run with the spine.

use super::*;
use crate::geometry::abscissa;
use crate::geometry::adaptor::CurveAdaptor;
use crate::topology::comp_curve::tests::{compare, occt_output, read_wire, section};
use crate::topology::comp_curve::CompCurve;
use crate::topology::Model;

fn length(c: &GeomAdaptor) -> f64 {
    abscissa::length(c, None)
}

/// The edges of the wire dumped in section `title`: orientation flag, range, start, end.
fn dumped_edges(title: &str) -> Vec<(bool, f64, f64, Point, Point)> {
    let sections = occt_output();
    edges_of(&section(&sections, title).1)
}

/// The edges dumped in `lines`.
fn edges_of(lines: &[String]) -> Vec<(bool, f64, f64, Point, Point)> {
    let nums = |s: &str| -> Vec<f64> { s.split_whitespace().filter_map(|w| w.parse().ok()).collect() };
    let mut out = Vec::new();
    let mut i = 0;
    while i < lines.len() {
        if let Some(rest) = lines[i].strip_prefix("edge or ") {
            let head = nums(rest);
            let (s, e) = (nums(lines[i + 1].strip_prefix("start ").unwrap()), nums(lines[i + 2].strip_prefix("end ").unwrap()));
            out.push((head[0] == 1.0, head[1], head[2], Point::new(s[0], s[1], s[2]), Point::new(e[0], e[1], e[2])));
            i += 3;
        } else {
            i += 1;
        }
    }
    out
}

/// `wire` against the dumped wire `title`: each edge's direction, range and ends.
#[track_caller]
fn same_wire(wire: &EdgeWire, title: &str) {
    let expected = dumped_edges(title);
    assert_eq!(wire.edges.len(), expected.len(), "{title}: edge count");
    for (k, (edge, (reversed, first, last, start, end))) in wire.edges.iter().zip(expected).enumerate() {
        assert_eq!(edge.reversed, reversed, "{title} edge {k}: direction");
        assert!((edge.first - first).abs() < 1e-12 && (edge.last - last).abs() < 1e-12, "{title} edge {k}: range {} {} vs {first} {last}", edge.first, edge.last);
        let (t0, t1) = edge.run();
        assert!(edge.point(t0).unwrap().distance(start) < 1e-12 && edge.point(t1).unwrap().distance(end) < 1e-12, "{title} edge {k}: ends");
        assert!(edge.start.unwrap().point.distance(start) < 1e-7 && edge.end.unwrap().point.distance(end) < 1e-7, "{title} edge {k}: vertices");
    }
}

/// The reduced arc lengths of wire A as OCCT's, to 1e-6; its whole length (the first
/// entry) as an independent fine Gauss integration of the speed gives it, to 1e-9 relative
/// (OCCT's, 37.447968511, is 7e-8 short); of a wire of no length, zeros and a final 1.
#[test]
fn reduced_abscissae_as_occt() {
    let mut model = Model::new();
    let a = read_wire(&mut model, include_str!("tests/guide_a.brep"));
    let wire = EdgeWire::from_model(&model, a, &mut ShapeKeys::new()).unwrap();
    let acr = compute_acr(&wire, &length);
    let sections = occt_output();
    let line = sections.iter().flat_map(|(_, l)| l.iter()).find(|l| l.starts_with("acr ")).unwrap();
    let expected: Vec<f64> = line.split_whitespace().skip(1).map(|w| w.parse().unwrap()).collect();
    assert_eq!(acr.len(), expected.len());
    for (a, b) in acr.iter().zip(&expected).skip(1) {
        assert!((a - b).abs() < 1e-6, "{acr:?} vs {expected:?}");
    }
    // (Gauss rules on 2000 pieces of each edge's range.)
    let total: f64 = wire
        .edges
        .iter()
        .map(|e| {
            let c = e.curve.as_ref().unwrap();
            let pieces = 2000;
            (0..pieces)
                .map(|k| {
                    let (p, q) = (e.first + (e.last - e.first) * k as f64 / pieces as f64, e.first + (e.last - e.first) * (k + 1) as f64 / pieces as f64);
                    let (m, h) = (0.5 * (p + q), 0.5 * (q - p));
                    crate::geometry::measure::GAUSS.iter().map(|&(x, w)| w * h * c.d1(m + h * x).1.length()).sum::<f64>()
                })
                .sum::<f64>()
        })
        .sum();
    assert!((acr[0] - total).abs() < 1e-9 * total, "{} vs {total}", acr[0]);
    assert!((acr[0] - expected[0]).abs() < 1e-7 * total);
    assert_eq!(compute_acr(&wire, &|_| 0.0), vec![0.0, 0.0, 0.0, 1.0]);
}

/// Closed wire B (a line, an arc, a reversed line, a line) opened at its point nearest
/// points: inside the first line both ways (run backwards, the wire begins at the far end of
/// the cut edge), on the arc both ways, at a vertex exactly, and next to a vertex within the
/// tolerance (the edge's other vertex then becomes the origin); as OCCT opens it, wires and
/// composite curves.
#[test]
fn search_origin_as_occt() {
    let mut model = Model::new();
    let b = read_wire(&mut model, include_str!("tests/guide_b.brep"));
    same_wire(&EdgeWire::from_model(&model, b, &mut ShapeKeys::new()).unwrap(), "wire B");
    let arc_point = Point::new(14.931969619160718, 5.8219949365267869, 0.0);
    let outward = (arc_point - Point::new(10.0, 5.0, 0.0)) * 0.2;
    let arc_tangent = Vector::new(-outward.y, outward.x, 0.0);
    let cases = [
        (0, Point::new(4.0, -1.0, 0.0), Vector::new(-1.0, 0.0, 0.0)),
        (1, Point::new(4.0, -1.0, 0.0), Vector::new(1.0, 0.2, 0.0)),
        (2, arc_point + outward, arc_tangent),
        (3, arc_point + outward, -arc_tangent),
        (5, Point::new(10.0, 10.0, 0.0), Vector::new(-1.0, 0.0, 0.0)),
        (6, Point::new(9.5, 10.3, 0.0), Vector::new(0.0, 1.0, 0.0)),
        (6, Point::new(0.0, 10.0, 0.0), Vector::new(0.0, 1.0, 0.0)),
    ];
    for (k, p, dir) in cases {
        let wire = search_origin(&model, b, p, dir, 1.0).unwrap();
        same_wire(&wire, &format!("wire B origin {k}"));
        compare(&CompCurve::from_wire(wire).unwrap(), &format!("B origin {k}"), true);
    }
    let open = read_wire(&mut model, include_str!("tests/guide_a.brep"));
    assert!(matches!(search_origin(&model, open, Point::ORIGIN, Vector::new(1.0, 0.0, 0.0), 1.0), Err(Error::InvalidInput(_))));
}

/// Open guides along open spines, turned to run the spine's way as OCCT turns them (a
/// straight spine, a straight guide, both bent): the guide's start and end.
#[test]
fn compatible_wires_as_occt() {
    let mut model = Model::new();
    let s = [
        read_wire(&mut model, include_str!("tests/guide_s1.brep")),
        read_wire(&mut model, include_str!("tests/guide_s2.brep")),
        read_wire(&mut model, include_str!("tests/guide_s3.brep")),
        read_wire(&mut model, include_str!("tests/guide_s4.brep")),
    ];
    let g = [
        read_wire(&mut model, include_str!("tests/guide_g1.brep")),
        read_wire(&mut model, include_str!("tests/guide_g2.brep")),
        read_wire(&mut model, include_str!("tests/guide_g3.brep")),
        read_wire(&mut model, include_str!("tests/guide_g4.brep")),
        read_wire(&mut model, include_str!("tests/guide_g5.brep")),
    ];
    for (k, (spine, guide)) in [(s[0], g[0]), (s[1], g[1]), (s[2], g[2]), (s[3], g[3]), (s[2], g[4])].into_iter().enumerate() {
        let wire = compatible_guide(&model, spine, guide, &length).unwrap();
        // (OCCT's dump: the spine's wire, then the guide's; its edges cut to match.)
        let sections = occt_output();
        let expected = sections.iter().filter(|(t, _)| t.starts_with(&format!("wire {} ", k + 1))).nth(1).map(|(_, l)| edges_of(l)).unwrap();
        let start = wire.edges[0].start.unwrap().point;
        let end = wire.edges.last().unwrap().end.unwrap().point;
        assert!(start.distance(expected[0].3) < 1e-7 && end.distance(expected.last().unwrap().4) < 1e-7, "compat {}: {start:?} {end:?}", k + 1);
        let curve = CompCurve::from_wire(wire).unwrap();
        assert!(curve.value(curve.first_parameter()).distance(expected[0].3) < 1e-7);
    }
}

/// A closed circular guide around a straight spine, queried from opposite directions: the
/// forward wire starts at the point nearest, the backward one at the far end of the cut
/// edge; an open guide drawn the other way round is turned.
#[test]
fn closed_and_reversed_guides() {
    use crate::geometry::curve::{Circle, Curve};
    use crate::math::{Direction, Frame};
    let mut model = Model::new();
    let circle = Curve::Circle(Circle::new(Frame::from_z(Point::ORIGIN, Direction::Z), 3.0));
    let ring = crate::make::edge_from_curve(&mut model, circle).unwrap();
    let ring = crate::make::wire(&mut model, &[ring]).unwrap();
    let p = Point::new(0.0, 5.0, 0.0);
    let forward = search_origin(&model, ring, p, Vector::new(-1.0, 0.0, 0.0), 1e-2).unwrap();
    let first = &forward.edges[0];
    assert!(first.point(first.run().0).unwrap().distance(Point::new(0.0, 3.0, 0.0)) < 1e-9);
    let backward = search_origin(&model, ring, p, Vector::new(1.0, 0.0, 0.0), 1e-2).unwrap();
    let first = &backward.edges[0];
    // (The cut edge runs from the nearest point to the circle's seam at (3, 0, 0).)
    assert!(first.point(first.run().0).unwrap().distance(Point::new(3.0, 0.0, 0.0)) < 1e-9);
    let spine = crate::make::edge_between(&mut model, Point::ORIGIN, Point::new(10.0, 0.0, 0.0)).unwrap();
    let guide = crate::make::edge_between(&mut model, Point::new(10.0, 2.0, 1.0), Point::new(0.0, 2.0, 1.0)).unwrap();
    let wire = compatible_guide(&model, spine, guide, &length).unwrap();
    assert!(wire.edges[0].start.unwrap().point.distance(Point::new(0.0, 2.0, 1.0)) < 1e-12);
}
