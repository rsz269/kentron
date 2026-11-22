use super::*;
use crate::geometry::surface::Torus;
use crate::math::{Direction, Frame};

fn line(o: [f64; 3], d: [f64; 3]) -> LineCurve {
    LineCurve { origin: Point::new(o[0], o[1], o[2]), direction: Vector::new(d[0], d[1], d[2]) }
}

#[test]
fn line_meets_a_tilted_cylinder() {
    // A cylinder of radius 2 about the axis through (1, 0, 0) along (0, 1, 1)/√2: the line
    // along x through the origin meets it where its distance from the axis is 2, at
    // x = 1 ± 2 (the line is square to the axis).
    let frame = Frame::new(Point::new(1.0, 0.0, 0.0), Direction::new(0.0, 1.0, 1.0), Direction::new(1.0, 0.0, 0.0));
    let quad = cylinder_quadric(&Ax3::from_frame(&frame), 2.0);
    let ana = line_quadric([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], &quad);
    assert!(ana.done && !ana.inquadric);
    let mut ws: Vec<f64> = ana.points.iter().map(|p| p.1).collect();
    ws.sort_by(f64::total_cmp);
    assert_eq!(ws.len(), 2);
    assert!((ws[0] + 1.0).abs() < 1e-12 && (ws[1] - 3.0).abs() < 1e-12, "{ws:?}");
    // Along the axis, half a unit from it: a constant polynomial, no root.
    let along = line_quadric([1.5, 0.0, 0.0], [0.0, 0.5f64.sqrt(), 0.5f64.sqrt()], &quad);
    assert!(along.points.is_empty());
}

#[test]
fn line_meets_a_torus_four_times() {
    let torus = Torus::new(Frame::WORLD, 3.0, 1.0);
    let mut ws: Vec<f64> = line_torus(&line([-10.0, 0.0, 0.0], [1.0, 0.0, 0.0]), &torus).iter().map(|h| h.0).collect();
    ws.sort_by(f64::total_cmp);
    let want = [6.0, 8.0, 12.0, 14.0];
    assert_eq!(ws.len(), 4, "{ws:?}");
    for (w, e) in ws.iter().zip(want) {
        assert!((w - e).abs() < 1e-9, "{ws:?}");
    }
}

#[test]
fn polyhedron_triangles_cover_the_grid() {
    // OCCT's numbering: each cell of the grid two triangles, every point used, every
    // triangle's three points distinct and within one cell.
    let surface = Surface::Plane(crate::geometry::surface::Plane::new(Point::ORIGIN, Direction::Z));
    let us: Vec<f64> = (0..=4).map(|i| i as f64).collect();
    let vs: Vec<f64> = (0..=3).map(|j| j as f64).collect();
    let poly = Polyhedron::new(&surface, &us, &vs);
    assert_eq!(poly.nb_triangles(), 24);
    let mut used = vec![false; poly.pnts.len()];
    let mut area = 0.0;
    for t in 1..=poly.nb_triangles() {
        let [a, b, c] = poly.triangle(t);
        assert!(a != b && b != c && a != c);
        for i in [a, b, c] {
            used[i] = true;
        }
        let (pa, pb, pc) = (poly.pnts[a], poly.pnts[b], poly.pnts[c]);
        area += 0.5 * (pb - pa).cross(pc - pa).length();
    }
    assert!(used[1..].iter().all(|&u| u));
    assert!((area - 12.0).abs() < 1e-12);
    // A flat grid's deflection is the floor, 1e-4.
    assert_eq!(poly.deflection, 0.0001);
}

#[test]
fn line_box_clips_to_the_box() {
    let b = Bnd { lo: Point::new(0.0, 0.0, 0.0), hi: Point::new(1.0, 2.0, 3.0) };
    let (a, z) = lin_box(Point::new(-1.0, 1.0, 1.0), Vector::new(1.0, 0.0, 0.0), &b).unwrap();
    assert!((a - 1.0).abs() < 1e-15 && (z - 2.0).abs() < 1e-15);
    assert!(lin_box(Point::new(-1.0, 5.0, 1.0), Vector::new(1.0, 0.0, 0.0), &b).is_none());
}
