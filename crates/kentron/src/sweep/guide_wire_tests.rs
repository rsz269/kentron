//! Checked against OCCT 8.1 on the wires of `tests/guide_*.brep`: `tests/guide_occt.txt`
//! is what OCCT 8.1 prints (see `topology::comp_curve`'s tests for the format).

use super::*;
use crate::topology::comp_curve::tests::{compare, dump_comp, dump_wire, length, occt_output, read_wire};
use crate::topology::comp_curve::CompCurve;

/// OCCT's lines from the one starting with `first` on.
fn from<'a>(expected: &'a [&'a str], first: &str) -> &'a [&'a str] {
    let at = expected.iter().position(|l| l.starts_with(first)).unwrap_or_else(|| panic!("no `{first}` in OCCT's output"));
    &expected[at..]
}

#[test]
fn reduced_abscissae_as_occt() {
    let expected = occt_output();
    let mut model = Model::new();
    let a = read_wire(&mut model, include_str!("tests/guide_a.brep"));
    let wire = EdgeWire::from_model(&model, a, &mut ShapeKeys::new()).unwrap();
    let acr = compute_acr(&wire, &length);
    let out = format!("acr {:e} {:e} {:e} {:e}\n", acr[0], acr[1], acr[2], acr[3]);
    // (Kentron's own lengths, standing for GCPnts_AbscissaPoint's, differ by up to 1e-7.)
    let worst = compare(from(&expected, "acr "), &out, 1e-6);
    eprintln!("ComputeACR: largest relative deviation from OCCT {worst:e}");
    // A wire reduced to a point: 1 at its end.
    let point = compute_acr(&wire, &|_| 0.0);
    assert_eq!(point, vec![0.0, 0.0, 0.0, 1.0]);
}

/// `BRepFill::SearchOrigin` on the closed wire B: cut at a new vertex (on a line, on the
/// arc), started at a vertex (within the tolerance, or exactly there), either way round.
#[test]
fn search_origin_as_occt() {
    let expected = occt_output();
    let mut model = Model::new();
    let b = read_wire(&mut model, include_str!("tests/guide_b.brep"));
    let cases = [
        (Point::new(4., -1., 0.), Vector::new(-1., 0., 0.), 1e-5),
        (Point::new(4., -1., 0.), Vector::new(1., 0.2, 0.), 1e-5),
        (Point::new(16., 6., 1.), Vector::new(0., 1., 0.), 1e-5),
        (Point::new(16., 6., 1.), Vector::new(0., -1., 0.), 1e-5),
        (Point::new(5., 10.5, 0.), Vector::new(1., 0., 0.), 1e-5),
        (Point::new(0.2, 10.1, 0.), Vector::new(0., -1., 0.), 1.),
        (Point::new(0.0, 10., 0.), Vector::new(0., 1., 0.), 1e-5),
    ];
    let mut worst = 0.0f64;
    for (k, (p, dir, tolerance)) in cases.into_iter().enumerate() {
        let tag = format!("B origin {k}");
        let wire = search_origin(&model, b, p, dir, tolerance).unwrap();
        let mut out = dump_wire(&tag, &wire);
        let c = CompCurve::from_wire(wire).unwrap();
        out += &dump_comp(&tag, &c, &[0., 0.5, 1., 2.5, c.last_parameter()]);
        worst = worst.max(compare(from(&expected, &format!("== wire {tag} ")), &out, 1e-12));
    }
    eprintln!("SearchOrigin: largest relative deviation from OCCT {worst:e}");
    // An open wire is refused.
    let a = read_wire(&mut model, include_str!("tests/guide_a.brep"));
    assert!(search_origin(&model, a, Point::ORIGIN, Vector::new(1., 0., 0.), 1e-5).is_err());
}

/// `BRepFill_CompatibleWires` on an open spine and guide running the other way with a
/// different number of edges: (1) a line and an arc against three lines (one reversed);
/// (2) three lines and a B-spline against a C1 B-spline and a line, cut on reversed edges;
/// (3) two straight wires, turned by the distances between their ends; (4) two cornered
/// wires running the same way, left as they are; (5) a straight spine and an arc.
#[test]
fn compatible_wires_as_occt() {
    let expected = occt_output();
    let mut model = Model::new();
    let mut worst = 0.0f64;
    for (tag, spine, guide) in [
        ("1", include_str!("tests/guide_s1.brep"), include_str!("tests/guide_g1.brep")),
        ("2", include_str!("tests/guide_s2.brep"), include_str!("tests/guide_g2.brep")),
        ("3", include_str!("tests/guide_s3.brep"), include_str!("tests/guide_g3.brep")),
        ("4", include_str!("tests/guide_s4.brep"), include_str!("tests/guide_g4.brep")),
        ("5", include_str!("tests/guide_s3.brep"), include_str!("tests/guide_g5.brep")),
    ] {
        let spine = read_wire(&mut model, spine);
        let guide = read_wire(&mut model, guide);
        let (s, g) = compatible_wires(&model, spine, guide, &length).unwrap();
        let mut out = dump_wire(tag, &s);
        out += &dump_wire(tag, &g);
        let c = CompCurve::from_wire(g.clone()).unwrap();
        out += &dump_comp(tag, &c, &[0., 0.5, 1., 1.5, 2.5, c.last_parameter()]);
        // (The cuts are at abscissae from kentron's own lengths.)
        worst = worst.max(compare(&from(&expected, &format!("== compat {tag} done 1"))[1..], &out, 1e-6));
        // The guide alone, as the pipe takes it.
        let alone = compatible_guide(&model, spine, guide, &length).unwrap();
        assert_eq!(dump_wire(tag, &alone), dump_wire(tag, &g));
    }
    eprintln!("CompatibleWires: largest relative deviation from OCCT {worst:e}");
    // A closed wire with an open one: not the same topology.
    let open = read_wire(&mut model, include_str!("tests/guide_s1.brep"));
    let closed = read_wire(&mut model, include_str!("tests/guide_b.brep"));
    assert!(compatible_guide(&model, open, closed, &length).is_err());
}

/// The continuity where a wire's edges meet (`WireContinuity`), and `MakeWire` turning the
/// pieces it is given to run on from the wire's end.
#[test]
fn wire_continuity_and_make_wire() {
    let mut model = Model::new();
    let mut keys = ShapeKeys::new();
    let a = read_wire(&mut model, include_str!("tests/guide_a.brep"));
    let b = read_wire(&mut model, include_str!("tests/guide_b.brep"));
    let wa = EdgeWire::from_model(&model, a, &mut keys).unwrap();
    let wb = EdgeWire::from_model(&model, b, &mut keys).unwrap();
    // A: a corner between the line and the arc; B: corners too.
    assert_eq!(wire_continuity(&wa).unwrap(), Continuity::C0);
    assert_eq!(wire_continuity(&wb).unwrap(), Continuity::C0);
    // The pieces of A's first edge, cut at a quarter, and its others: one closed wire again.
    let e = &wa.edges[0];
    let pieces = trim_edge(e, &[0.25], 0.0, 1.0, true, &mut keys).unwrap();
    assert_eq!(pieces.len(), 2);
    assert!((pieces[0].last - 2.5).abs() < 1e-15 && pieces[1].first == pieces[0].last);
    let mut made = MakeWire::default();
    for piece in pieces.iter().rev() {
        made.add(piece, &mut keys).unwrap();
    }
    let w = made.wire();
    assert!(!w.closed);
    assert_eq!(w.explore().len(), 2);
    // An edge away from the wire is refused.
    let far = trim_edge(&wa.edges[2], &[0.5], 0.0, 1.0, true, &mut keys).unwrap();
    assert!(made_with(&pieces[..1], &far[1]).is_err());
}

fn made_with(edges: &[WireEdge], more: &WireEdge) -> Result<()> {
    let mut keys = ShapeKeys::new();
    let mut made = MakeWire::default();
    for e in edges {
        made.add(e, &mut keys)?;
    }
    made.add(more, &mut keys)
}
