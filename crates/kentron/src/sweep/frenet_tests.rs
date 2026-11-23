//! The Frenet frame against the textbook formulas, the default normal of straight pieces,
//! one-sided limits at an inflection and at a C1 knot.

use super::*;
use crate::geometry::curve::{BSplineCurve, Circle, Line};
use crate::math::{Point, Transformable};

fn p(x: f64, y: f64, z: f64) -> Point {
    Point::new(x, y, z)
}

fn planar_s() -> Curve {
    Curve::BSpline(BSplineCurve::new(vec![p(0.0, 0.0, 0.0), p(1.0, 2.0, 0.0), p(2.0, 2.0, 0.0), p(3.0, -2.0, 0.0), p(4.0, -2.0, 0.0), p(5.0, 0.0, 0.0)], vec![0.0, 1.0, 2.0, 3.0], vec![4, 1, 1, 4], 3))
}

/// The textbook frame from the curve's own derivatives at `u` (forward travel).
fn textbook(curve: &Curve, u: f64) -> (Vector, Vector, Vector) {
    let (_, d1, d2) = curve.d2(u);
    let t = d1.normalized();
    let b = d1.cross(d2).normalized();
    (t, b.cross(t), b)
}

/// Along a helix (a fitted B-spline) the normal points at the axis and the binormal is
/// `(c sin t, −c cos t, r)/√(r²+c²)`, to the fit's accuracy; everywhere the frame is the
/// textbook one to rounding.
#[test]
fn along_a_helix_the_normal_points_at_the_axis() {
    let (r, c) = (2.0, 0.5);
    let curve = super::super::corrected::tests::helix(r, c, 2.0, 400);
    let Curve::BSpline(h) = &curve else { unreachable!() };
    let (first, last) = (h.first_parameter(), h.last_parameter());
    let law = FrenetLaw::new(&curve, first, last);
    for i in 0..=50 {
        let f = i as f64 / 50.0;
        let u = first + (last - first) * f;
        let (t, n) = law.at_fraction(f);
        let (et, en, _) = textbook(&curve, u);
        assert!((t - et).length() < 1e-12 && (n - en).length() < 1e-10, "at {f}");
        let q = curve.point(u);
        let angle = q.y.atan2(q.x);
        let towards_axis = Vector::new(-angle.cos(), -angle.sin(), 0.0);
        assert!((n - towards_axis).length() < 1e-4, "at {f}: {n:?}");
        let b = Vector::new(c * angle.sin(), -c * angle.cos(), r) * (1.0 / (r * r + c * c).sqrt());
        assert!((t.cross(n) - b).length() < 1e-4, "at {f}");
    }
}

/// Along a circle the normal points at the centre and the binormal is the axis in the sense
/// of travel; run backwards, tangent and binormal flip and the normal stays.
#[test]
fn along_a_circle_the_normal_points_at_the_centre() {
    let frame = crate::math::Frame::new(p(1.0, -2.0, 0.5), Direction::from_vector(Vector::new(0.3, -0.4, 1.0)), Direction::from_vector(Vector::new(1.0, 0.0, -0.3)));
    let circle = Curve::Circle(Circle::new(frame, 3.0));
    let axis = Vector::from(frame.z());
    let forward = FrenetLaw::new(&circle, 0.3, 4.0);
    let backward = FrenetLaw::new(&circle, 4.0, 0.3);
    for i in 0..=10 {
        let f = i as f64 / 10.0;
        let u = 0.3 + 3.7 * f;
        let (t, n) = forward.at_fraction(f);
        let centre = (frame.origin() - circle.point(u)).normalized();
        assert!((n - centre).length() < 1e-12 && (t.cross(n) - axis).length() < 1e-12);
        let (tb, nb) = backward.at_fraction(1.0 - f);
        assert!((tb + t).length() < 1e-12 && (nb - n).length() < 1e-12 && (tb.cross(nb) + axis).length() < 1e-12);
    }
}

/// Along a line in each octant the normal is the default X axis of the tangent, exactly;
/// along +Z it is +X, along −Y it is −Z.
#[test]
fn along_a_line_the_normal_is_the_default_one() {
    for sx in [-1.0, 1.0] {
        for sy in [-1.0, 1.0] {
            for sz in [-1.0, 1.0] {
                let d = Direction::from_vector(Vector::new(0.3 * sx, 0.5 * sy, 0.8 * sz));
                let line = Curve::Line(Line::new(p(1.0, 2.0, 3.0), d));
                let (t, n) = FrenetLaw::new(&line, 0.0, 4.0).at_fraction(0.4);
                assert_eq!(t, Vector::from(d));
                assert_eq!(n, Vector::from(crate::math::Frame::from_z(Point::ORIGIN, d).x()));
            }
        }
    }
    let along = |d: Direction| FrenetLaw::new(&Curve::Line(Line::new(Point::ORIGIN, d)), 0.0, 1.0).at_fraction(0.5).1;
    assert_eq!(along(Direction::Z), Vector::new(1.0, 0.0, 0.0));
    assert_eq!(along(Direction::NEG_Y), Vector::new(0.0, 0.0, -1.0));
}

/// The inflection of a planar S (found here by bisection on the cross product's normal
/// component) is a singular point: at it the frame is the limit from before, just after it
/// the binormal is the opposite one; and the S's frame flips by a half turn there.
#[test]
fn at_an_inflection_the_frame_is_the_limit_from_before() {
    let curve = planar_s();
    let z = |u: f64| {
        let (_, d1, d2) = curve.d2(u);
        d1.cross(d2).z
    };
    let (mut lo, mut hi) = (1.0, 2.0);
    assert!(z(lo) * z(hi) < 0.0);
    for _ in 0..200 {
        let m = 0.5 * (lo + hi);
        if m <= lo || m >= hi {
            break;
        }
        if z(m) * z(lo) > 0.0 { lo = m } else { hi = m }
    }
    let s = 0.5 * (lo + hi);
    let frenet = Frenet::new(GeomAdaptor::new(&curve, 0.0, 3.0));
    assert!(frenet.intervals().iter().any(|&b| (b - s).abs() < 1e-9), "{:?} vs {s}", frenet.intervals());
    let law = FrenetLaw::new(&curve, 0.0, 3.0);
    let at = |u: f64| law.at_fraction(u / 3.0);
    let binormal = |(t, n): (Vector, Vector)| t.cross(n);
    let before = binormal(at(s - 1e-3));
    let after = binormal(at(s + 1e-3));
    assert!((before + after).length() < 1e-12);
    let exact = frenet.intervals().iter().copied().find(|&b| (b - s).abs() < 1e-9).unwrap();
    assert!((binormal(law.at_fraction(exact / 3.0)) - before).length() < 1e-12);
    assert!((binormal(at(exact + 1e-7)) - after).length() < 1e-12);
    // (Seen from the interval after it, the limit from after.)
    let i = frenet.interval_from_right(exact);
    assert!((frenet.frame_in(exact, i).b - after).length() < 1e-12);
}

/// At a C1 knot the frame from the left is the left-hand curvature's; the whole-curve frame
/// takes the right-hand side.
#[test]
fn at_a_c1_knot_the_frame_is_the_left_one() {
    let curve = Curve::BSpline(BSplineCurve::new(vec![p(0.0, 0.0, 0.0), p(1.0, 2.0, 0.0), p(2.0, 2.0, 1.0), p(3.0, 0.0, 2.0), p(4.0, -1.0, 1.0), p(5.0, 1.0, 0.0), p(6.0, 2.0, 2.0)], vec![0.0, 1.0, 2.0, 3.0], vec![4, 2, 1, 4], 3));
    let side = |a: f64, b: f64| {
        let (_, d1, d2) = GeomAdaptor::new(&curve, a, b).d2(1.0);
        d1.cross(d2).normalized().cross(d1.normalized())
    };
    let (left, right) = (side(0.0, 1.0), side(1.0, 3.0));
    assert!((left - right).length() > 1e-2);
    let (_, n) = FrenetLaw::new(&curve, 0.0, 3.0).at_fraction(1.0 / 3.0);
    assert!((n - left).length() < 1e-12);
    let frenet = Frenet::new(GeomAdaptor::new(&curve, 0.0, 3.0));
    assert!((frenet.frame(1.0).n - right).length() < 1e-12);
    // (Run backwards, the knot is met from the other side first.)
    let (_, n) = FrenetLaw::new(&curve, 3.0, 0.0).at_fraction(2.0 / 3.0);
    assert!((n - right).length() < 1e-12);
}

/// A straight stretch of a B-spline (a line's poles) takes the default frame, and a tilted
/// copy of the S keeps the textbook frame turned with it.
#[test]
fn straight_stretches_take_the_default_frame() {
    let straight = Curve::BSpline(BSplineCurve::new(vec![p(0.0, 0.0, 0.0), p(1.0, 1.0, 1.0), p(2.0, 2.0, 2.0), p(4.0, 4.0, 4.0)], vec![0.0, 1.0], vec![4, 4], 3));
    let (t, n) = FrenetLaw::new(&straight, 0.0, 1.0).at_fraction(0.3);
    let d = Direction::from_vector(t);
    assert!((t - Vector::new(1.0, 1.0, 1.0).normalized()).length() < 1e-15);
    assert!((n - Vector::from(crate::math::Frame::from_z(Point::ORIGIN, d).x())).length() < 1e-15);
    let tilt = crate::math::Transform::rotation(crate::math::Axis::new(Point::ORIGIN, Direction::X), 0.4);
    let curve = planar_s().transformed(&tilt);
    let law = FrenetLaw::new(&curve, 0.0, 3.0);
    for u in [0.1, 0.7, 2.2, 2.9] {
        let (t, n) = law.at_fraction(u / 3.0);
        let (et, en, _) = textbook(&curve, u);
        assert!((t - et).length() < 1e-12 && (n - en).length() < 1e-12);
    }
}
