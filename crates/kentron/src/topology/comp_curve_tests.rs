//! Checked against OCCT 8.1: the wires are `sweep/tests/guide_*.brep`, and
//! `sweep/tests/guide_occt.txt` is what OCCT 8.1 prints for them (`%.17g`), in the format
//! [`dump_comp`] and [`dump_wire`] write.

use super::*;

/// `GeomAbs_Shape` as OCCT numbers it.
fn shape_number(c: Continuity) -> usize {
    match c {
        Continuity::C0 => 0,
        Continuity::G1 => 1,
        Continuity::C1 => 2,
        Continuity::G2 => 3,
        Continuity::C2 => 4,
        Continuity::C3 => 5,
        Continuity::CN => 6,
    }
}

fn orientation_number(o: Orientation) -> usize {
    match o {
        Orientation::Forward => 0,
        Orientation::Reversed => 1,
        Orientation::Internal => 2,
        Orientation::External => 3,
    }
}

fn xyz(name: &str, x: f64, y: f64, z: f64) -> String {
    format!("{name} {x:e} {y:e} {z:e}\n")
}

/// A composite curve: its bounds, continuity, resolution, intervals and derivatives at `us`.
pub(crate) fn dump_comp(tag: &str, c: &CompCurve, us: &[f64]) -> String {
    let mut out = format!(
        "== {tag} first {:e} last {:e} closed {} cont {} res {:e}\n",
        c.first_parameter(),
        c.last_parameter(),
        c.is_closed() as usize,
        shape_number(c.continuity()),
        c.resolution(1e-3)
    );
    for s in [Continuity::C0, Continuity::C1, Continuity::C2, Continuity::C3, Continuity::CN] {
        out += &format!("intervals {}:", shape_number(s));
        for t in c.intervals(s) {
            out += &format!(" {t:e}");
        }
        out += "\n";
    }
    for &u in us {
        let (p, d1, d2, d3) = c.d3(u);
        let d4 = c.dn(u, 4);
        out += &format!("u {u:e}\n");
        out += &xyz("p", p.x, p.y, p.z);
        out += &xyz("d1", d1.x, d1.y, d1.z);
        out += &xyz("d2", d2.x, d2.y, d2.z);
        out += &xyz("d3", d3.x, d3.y, d3.z);
        out += &xyz("dn4", d4.x, d4.y, d4.z);
    }
    out
}

/// A wire's edges in explorer order: orientation, range, and the points they run between.
pub(crate) fn dump_wire(tag: &str, w: &EdgeWire) -> String {
    let mut out = format!("== wire {tag} closed {}\n", w.closed as usize);
    for (e, _) in w.explored_edges() {
        out += &format!("edge or {} range {:e} {:e}\n", orientation_number(e.orientation), e.first, e.last);
        let (a, b) = e.ends(true);
        let (a, b) = (a.unwrap().point, b.unwrap().point);
        out += &xyz("start", a.x, a.y, a.z);
        out += &xyz("end", b.x, b.y, b.z);
    }
    out
}

/// OCCT's output for the fixtures (comments left out).
pub(crate) fn occt_output() -> Vec<&'static str> {
    include_str!("../sweep/tests/guide_occt.txt").lines().filter(|l| !l.starts_with('#')).collect()
}

/// Compares the lines of `actual` with those of `expected` (from the first): words alike,
/// numbers within `tolerance` of each other relative to their size (at least 1). Returns
/// the largest difference (so relative) met. The resolutions (after `res`) are not
/// compared: `GeomAdaptor::resolution` is kentron's own for B-splines and circles.
pub(crate) fn compare(expected: &[&str], actual: &str, tolerance: f64) -> f64 {
    let actual: Vec<&str> = actual.lines().collect();
    assert!(expected.len() >= actual.len(), "{} lines expected, {} written", expected.len(), actual.len());
    let mut worst = 0.0f64;
    for (n, (e, a)) in expected.iter().zip(&actual).enumerate() {
        let (ew, aw): (Vec<&str>, Vec<&str>) = (e.split_whitespace().collect(), a.split_whitespace().collect());
        assert_eq!(ew.len(), aw.len(), "line {n}: expected `{e}`, got `{a}`");
        for (k, (x, y)) in ew.iter().zip(&aw).enumerate() {
            if k > 0 && ew[k - 1] == "res" {
                continue;
            }
            match (x.parse::<f64>(), y.parse::<f64>()) {
                (Ok(x), Ok(y)) => {
                    let d = (x - y).abs() / x.abs().max(y.abs()).max(1.0);
                    worst = worst.max(d);
                    assert!(d <= tolerance, "line {n}: expected `{e}`, got `{a}` ({d:e})");
                }
                _ => assert_eq!(x, y, "line {n}: expected `{e}`, got `{a}`"),
            }
        }
    }
    worst
}

pub(crate) fn read_wire(model: &mut Model, text: &str) -> Shape {
    crate::io::brep::read(model, &mut text.as_bytes()).unwrap()
}

/// The length of an edge's curve over its range (kentron's own measure, standing for
/// `GCPnts_AbscissaPoint::Length` in these tests).
pub(crate) fn length(c: &GeomAdaptor) -> f64 {
    crate::geometry::measure::length_between(c.curve(), c.first_parameter(), c.last_parameter())
}

/// Wire A (a line, a reversed arc and a C1 B-spline, placed), its trimmed curves and
/// wire B (closed: line, arc, reversed line, line) against OCCT's composite curves.
#[test]
fn composite_curves_as_occt() {
    let expected = occt_output();
    let mut model = Model::new();
    let a = read_wire(&mut model, include_str!("../sweep/tests/guide_a.brep"));
    let b = read_wire(&mut model, include_str!("../sweep/tests/guide_b.brep"));
    let wire_a = EdgeWire::from_model(&model, a, &mut ShapeKeys::new()).unwrap();
    let c = CompCurve::from_wire(wire_a.clone()).unwrap();
    let mut out = dump_wire("A", &wire_a);
    out += &dump_comp("A", &c, &[0., 0.3, 1., 1.5, 2., 2.25, 2.7, 3.]);
    out += &dump_comp("A trimmed", &c.trim(0.5, 2.5, 1e-9), &[0.5, 1., 1.75, 2., 2.5]);
    out += &dump_comp("A trimmed2", &c.trim(1.2, 1.6, 0.), &[1.2, 1.4, 1.6]);
    let start = out.lines().count();
    // (Knots at curvilinear abscissae: OCCT's lengths against kentron's own.)
    let ac = CompCurve::by_abscissa(wire_a, &length).unwrap();
    let ac_dump = dump_comp("A ac", &ac, &[0., 5., 10., 17., ac.last_parameter()]);
    let worst = compare(&expected, &out, 1e-12);
    // (Kentron's own lengths differ from GCPnts_AbscissaPoint's by up to 1e-7.)
    let worst_ac = compare(&expected[start..], &ac_dump, 1e-6);
    let start = start + ac_dump.lines().count() + 1;
    let wire_b = EdgeWire::from_model(&model, b, &mut ShapeKeys::new()).unwrap();
    let mut out = dump_wire("B", &wire_b);
    out += &dump_comp("B", &CompCurve::from_wire(wire_b).unwrap(), &[0., 0.5, 1., 1.25, 2., 2.5, 3.5, 4.]);
    let worst_b = compare(&expected[start..], &out, 1e-12);
    eprintln!("largest relative deviation from OCCT: {worst:e} (A), {worst_ac:e} (A by abscissa), {worst_b:e} (B)");
}

/// The edge carrying a parameter, near the knots, and `PTol` choosing the edge there.
#[test]
fn edges_at_parameters() {
    let mut model = Model::new();
    let a = read_wire(&mut model, include_str!("../sweep/tests/guide_a.brep"));
    let c = CompCurve::new(&model, a).unwrap();
    let (e, u) = c.edge(0.5);
    assert_eq!((e.first, e.last, u), (0.0, 10.0, 5.0));
    // The reversed arc runs from its last parameter.
    let (e, u) = c.edge(1.25);
    assert_eq!(e.orientation, Orientation::Reversed);
    assert!((u - 0.75 * e.last).abs() < 1e-15);
    // At a knot, the edge after it (with no tolerance).
    let (e, _) = c.edge(1.0);
    assert_eq!(e.orientation, Orientation::Reversed);
    assert_eq!(c.curve_type(), None);
    assert!(matches!(c.basis_curve(), Some(Curve::Line(_))));
}
