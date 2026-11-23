//! Composite curves of wires against OCCT 8.1's view of the same wires (bounds, closure,
//! interval lists, derivatives at given parameters, from `sweep/tests/guide_occt.txt`), and
//! the composite's own properties.

use super::*;
use crate::math::Point;

/// The sections of the OCCT dump, by title (the text after `== `).
pub(crate) fn occt_output() -> Vec<(String, Vec<String>)> {
    let mut out: Vec<(String, Vec<String>)> = Vec::new();
    for line in include_str!("../sweep/tests/guide_occt.txt").lines() {
        if let Some(title) = line.strip_prefix("== ") {
            out.push((title.to_string(), Vec::new()));
        } else if let Some(last) = out.last_mut() {
            last.1.push(line.to_string());
        }
    }
    out
}

/// The section titled `title` (up to its first field after the name).
pub(crate) fn section<'a>(sections: &'a [(String, Vec<String>)], title: &str) -> &'a (String, Vec<String>) {
    sections.iter().find(|(t, _)| t == title || t.starts_with(&format!("{title} "))).unwrap_or_else(|| panic!("no section {title}"))
}

fn numbers(text: &str) -> Vec<f64> {
    text.split_whitespace().filter_map(|w| w.parse().ok()).collect()
}

/// The value of `key` in a title like `A first 0 last 3 closed 0`.
fn field(title: &str, key: &str) -> f64 {
    let words: Vec<&str> = title.split_whitespace().collect();
    let i = words.iter().position(|&w| w == key).unwrap_or_else(|| panic!("no {key} in {title}"));
    words[i + 1].parse().unwrap()
}

/// The first wire (or edge) in a BREP text.
pub(crate) fn read_wire(model: &mut Model, text: &str) -> Shape {
    let shape = crate::io::brep::read(model, &mut text.as_bytes()).unwrap();
    model.explore(shape, NodeKind::Wire).next().unwrap_or(shape)
}

#[track_caller]
fn close(a: f64, b: f64, what: &str) {
    assert!((a - b).abs() <= 1e-12 * b.abs().max(1.0), "{what}: {a} vs OCCT's {b}");
}

/// Checks `curve` against the section titled `title`: bounds, closure, the interval lists
/// (where `intervals`) and every dumped derivative.
pub(crate) fn compare(curve: &CompCurve, title: &str, intervals: bool) {
    let sections = occt_output();
    let (head, lines) = section(&sections, title);
    close(curve.first_parameter(), field(head, "first"), "first");
    close(curve.last_parameter(), field(head, "last"), "last");
    assert_eq!(curve.is_closed(), field(head, "closed") == 1.0, "{title}: closed");
    let continuity = |c: u32| match c {
        0 => Continuity::C0,
        2 => Continuity::C1,
        4 => Continuity::C2,
        5 => Continuity::C3,
        _ => Continuity::CN,
    };
    let mut u = f64::NAN;
    for line in lines {
        let (key, rest) = line.split_once(' ').unwrap_or((line, ""));
        match key {
            "intervals" if intervals => {
                let (c, list) = rest.split_once(": ").unwrap();
                let expected = numbers(list);
                let ours = curve.intervals(continuity(c.parse().unwrap()));
                assert_eq!(ours.len(), expected.len(), "{title} intervals {c}: {ours:?} vs {expected:?}");
                for (a, b) in ours.iter().zip(&expected) {
                    close(*a, *b, &format!("{title} intervals {c}"));
                }
            }
            "u" => u = numbers(rest)[0],
            "p" | "d1" | "d2" | "d3" | "dn4" => {
                let e = numbers(rest);
                let v = match key {
                    "p" => curve.value(u) - Point::ORIGIN,
                    "d1" => curve.d1(u).1,
                    "d2" => curve.d2(u).2,
                    "d3" => curve.d3(u).3,
                    _ => curve.dn(u, 4),
                };
                for (a, b) in [v.x, v.y, v.z].iter().zip(&e) {
                    close(*a, *b, &format!("{title} {key} at {u}"));
                }
            }
            _ => {}
        }
    }
}

/// Wire A (a line, a reversed arc and a C1 B-spline, placed), two trims of it, and the
/// closed wire B, as OCCT reads them: to 1e-12. (OCCT lists the second trim's intervals
/// over the untrimmed range; here they span the trim, as everywhere else.)
#[test]
fn composite_curves_as_occt() {
    let mut model = Model::new();
    let a = read_wire(&mut model, include_str!("../sweep/tests/guide_a.brep"));
    let b = read_wire(&mut model, include_str!("../sweep/tests/guide_b.brep"));
    let curve = CompCurve::new(&model, a).unwrap();
    compare(&curve, "A", true);
    compare(&curve.trim(0.5, 2.5, 1e-9), "A trimmed", true);
    compare(&curve.trim(1.2, 1.6, 1e-9), "A trimmed2", false);
    compare(&CompCurve::new(&model, b).unwrap(), "B", true);
}

/// The edge holding a parameter: inside an edge, on the reversed arc (its parameter running
/// down), at a knot (the edge after it); not an elementary curve; the basis is the first
/// edge's (a line).
#[test]
fn edges_at_parameters() {
    let mut model = Model::new();
    let a = read_wire(&mut model, include_str!("../sweep/tests/guide_a.brep"));
    let curve = CompCurve::new(&model, a).unwrap();
    let (edge, x) = curve.edge(0.3);
    assert!(matches!(edge.curve(), Curve::Line(_)) && (x - 3.0).abs() < 1e-12);
    let (edge, x) = curve.edge(1.25);
    assert!(matches!(edge.curve(), Curve::Circle(_)));
    let range = edge.last_parameter() - edge.first_parameter();
    assert!((x - (edge.last_parameter() - 0.25 * range)).abs() < 1e-12, "{x}");
    let (edge, _) = curve.edge(2.0);
    assert!(matches!(edge.curve(), Curve::BSpline(_)));
    assert!(curve.curve_type().is_none());
    assert!(matches!(curve.basis_curve(), Some(Curve::Line(_))));
}

/// A closed square of four lines, one of them reversed in the wire: points continuous at the
/// knots, derivatives along the wire, a trim agreeing with the whole inside it.
#[test]
fn a_square_runs_along_its_wire() {
    let mut model = Model::new();
    let corners = [Point::new(0.0, 0.0, 0.0), Point::new(2.0, 0.0, 0.0), Point::new(2.0, 3.0, 0.0), Point::new(0.0, 3.0, 0.0)];
    let v: Vec<Shape> = corners.iter().map(|&p| crate::make::vertex(&mut model, p)).collect();
    let mut edges: Vec<Shape> = (0..4).map(|i| crate::make::edge(&mut model, v[i], v[(i + 1) % 4]).unwrap()).collect();
    // (The third edge made the other way, and used reversed.)
    edges[2] = crate::make::edge(&mut model, v[3], v[2]).unwrap().reversed();
    let wire = crate::make::wire(&mut model, &edges).unwrap();
    let curve = CompCurve::new(&model, wire).unwrap();
    assert!(curve.is_closed());
    assert_eq!((curve.first_parameter(), curve.last_parameter()), (0.0, 4.0));
    let start = curve.value(0.0);
    for k in 1..4 {
        let w = k as f64;
        assert!(curve.value(w - 1e-12).distance(curve.value(w)) < 1e-9);
    }
    assert!(curve.value(4.0).distance(start) < 1e-12);
    for i in 0..40 {
        let w = 0.05 + 0.1 * i as f64;
        let (p, d) = curve.d1(w);
        let ahead = curve.value(w + 1e-6);
        assert!(((ahead - p) * 1e6 - d).length() < 1e-5, "at {w}");
    }
    let trimmed = curve.trim(0.7, 3.2, 0.0);
    for i in 0..=10 {
        let w = 0.7 + 0.25 * i as f64;
        assert!(trimmed.value(w).distance(curve.value(w)) < 1e-15);
    }
}
