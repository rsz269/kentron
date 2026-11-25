//! Fillets along loops whose rounds are tighter than the ball (see [`super::tight`]).
//!
//! The references are independent of the fillet: every ball of the radius inside the prism
//! touching its top and sides (the union of balls along the top loop shrunk by the radius) must
//! be kept; outside the stretches where the ball changes, each straight edge's fillet is the
//! ordinary one, a quarter of the radius's disc from each side's square taken away; the round's
//! fillet is the sphere as wide as the round, on its axis; and from above, the top, the fillets
//! and the corners form one smooth surface (the scan of the conformance suite's crease check:
//! normals of neighbouring points apart by more than 0.1 rad, halved down to a thousandth of
//! the round's radius). The square's, the triangle's and the narrow rectangle's results also
//! hold every voxel of the voxel opening of their prisms (`tools/brute.py holds`, voxels of
//! 0.025 to 0.04): none of its material missing, none outside the prism.

use std::f64::consts::PI;

use crate::boolean::{PointClassifier, PointLocation};
use crate::geometry::curve::{Circle, Curve};
use crate::make;
use crate::math::{Direction, Frame, Point, Vector};
use crate::sweep;
use crate::topology::{Model, NodeKind, Shape};
use crate::{properties, raycast, validity};

use super::fillet;

/// A prism `h` high on the convex polygon `pts` (counter-clockwise), its corner `i` rounded by
/// `rho[i]`, and the edge of its top loop at the middle of the first side.
fn rounded_prism(model: &mut Model, pts: &[[f64; 2]], rho: &[f64], h: f64) -> (Shape, Shape) {
    let n = pts.len();
    // Each corner's round: its centre, its ends and its start angle and sweep.
    let rounds: Vec<([f64; 2], [f64; 2], [f64; 2], f64, f64)> = (0..n)
        .map(|i| {
            let (p, c, q) = (pts[(i + n - 1) % n], pts[i], pts[(i + 1) % n]);
            let unit = |a: [f64; 2]| {
                let l = (a[0] - c[0]).hypot(a[1] - c[1]);
                [(a[0] - c[0]) / l, (a[1] - c[1]) / l]
            };
            let (u, v) = (unit(p), unit(q));
            let half = (u[0] * v[0] + u[1] * v[1]).clamp(-1.0, 1.0).acos() / 2.0;
            let back = rho[i] / half.tan();
            let b = [u[0] + v[0], u[1] + v[1]];
            let lb = b[0].hypot(b[1]);
            let centre = [c[0] + b[0] / lb * rho[i] / half.sin(), c[1] + b[1] / lb * rho[i] / half.sin()];
            let (a0, a1) = ([c[0] + u[0] * back, c[1] + u[1] * back], [c[0] + v[0] * back, c[1] + v[1] * back]);
            (centre, a0, a1, (a0[1] - centre[1]).atan2(a0[0] - centre[0]), PI - 2.0 * half)
        })
        .collect();
    let mut edges = Vec::new();
    for i in 0..n {
        let (centre, _, a1, start, sweep) = rounds[i];
        let circle = Circle::new(Frame::new(Point::new(centre[0], centre[1], 0.0), Direction::Z, Direction::X), rho[i]);
        edges.push(make::edge_from_curve_range(model, Curve::Circle(circle), start, start + sweep).unwrap());
        let next = rounds[(i + 1) % n].1;
        edges.push(make::edge_between(model, Point::new(a1[0], a1[1], 0.0), Point::new(next[0], next[1], 0.0)).unwrap());
    }
    let wire = make::wire(model, &edges).unwrap();
    let face = make::planar_face(model, wire).unwrap();
    let prism = sweep::extrude(model, face, Vector::new(0.0, 0.0, h)).unwrap().shape;
    let (a, b) = (rounds[0].2, rounds[1].1);
    let m = Point::new((a[0] + b[0]) / 2.0, (a[1] + b[1]) / 2.0, h);
    let edge = model
        .unique_subshapes(prism, NodeKind::Edge)
        .into_iter()
        .find(|&e| {
            let (a, b) = model.edge_vertices(e);
            let (a, b) = (model.point(a), model.point(b));
            Point::new((a.x + b.x) / 2.0, (a.y + b.y) / 2.0, (a.z + b.z) / 2.0).distance(m) < 1e-9
        })
        .unwrap();
    (prism, edge)
}

/// The polygon `pts` moved in by `d` along each side's inward normal (its corners where the
/// moved sides meet).
fn inset(pts: &[[f64; 2]], d: f64) -> Vec<[f64; 2]> {
    let n = pts.len();
    let side = |i: usize| {
        let (p, q) = (pts[i], pts[(i + 1) % n]);
        let l = (q[0] - p[0]).hypot(q[1] - p[1]);
        let (t, m) = ([(q[0] - p[0]) / l, (q[1] - p[1]) / l], [-(q[1] - p[1]) / l, (q[0] - p[0]) / l]);
        ([p[0] + m[0] * d, p[1] + m[1] * d], t)
    };
    (0..n)
        .map(|i| {
            let ((p, t), (q, s)) = (side((i + n - 1) % n), side(i));
            // p + a t = q + b s.
            let det = t[0] * (-s[1]) - t[1] * (-s[0]);
            let a = ((q[0] - p[0]) * (-s[1]) - (q[1] - p[1]) * (-s[0])) / det;
            [p[0] + a * t[0], p[1] + a * t[1]]
        })
        .collect()
}

fn area(pts: &[[f64; 2]]) -> f64 {
    let n = pts.len();
    0.5 * (0..n).map(|i| pts[i][0] * pts[(i + 1) % n][1] - pts[(i + 1) % n][0] * pts[i][1]).sum::<f64>()
}

fn perimeter(pts: &[[f64; 2]]) -> f64 {
    let n = pts.len();
    (0..n).map(|i| (pts[(i + 1) % n][0] - pts[i][0]).hypot(pts[(i + 1) % n][1] - pts[i][1])).sum()
}

/// Creases (and points without a normal, and points in all) seen from above over a disc of
/// radius `half` about `(x, y)`: a 24 × 24 grid of rays down, neighbours whose normals are more
/// than 0.1 rad apart halved until they agree or lie within `fine` of each other (a crease).
fn creases(model: &mut Model, shape: Shape, x: f64, y: f64, top: f64, half: f64, fine: f64) -> (usize, usize, usize) {
    let caster = raycast::Raycaster::new(model, shape).unwrap();
    let hit = |model: &mut Model, p: [f64; 2]| -> Option<Option<Vector>> {
        let t = caster.cast(Point::new(p[0], p[1], top + 1.0), Vector::new(0.0, 0.0, -1.0)).unwrap();
        let t = t.into_iter().fold(f64::INFINITY, f64::min);
        t.is_finite().then(|| {
            let q = Point::new(p[0], p[1], top + 1.0 - t);
            properties::normals(model, shape, &[q]).unwrap()[0].map(Vector::from)
        })
    };
    let n = 24;
    let mut grid: Vec<Vec<Option<([f64; 2], Option<Vector>)>>> = vec![vec![None; n]; n];
    let (mut points, mut missing) = (0, 0);
    for (i, row) in grid.iter_mut().enumerate() {
        for (j, cell) in row.iter_mut().enumerate() {
            let (u, v) = ((i as f64 + 0.5) / n as f64 * 2.0 - 1.0, (j as f64 + 0.5) / n as f64 * 2.0 - 1.0);
            if u * u + v * v > 1.0 {
                continue;
            }
            let p = [x + u * half, y + v * half];
            if let Some(normal) = hit(model, p) {
                points += 1;
                missing += usize::from(normal.is_none());
                *cell = Some((p, normal));
            }
        }
    }
    let mut found = 0;
    for i in 0..n {
        for j in 0..n {
            for (di, dj) in [(1, 0), (0, 1)] {
                let next = grid.get(i + di).and_then(|row| row.get(j + dj)).copied().flatten();
                let (Some((mut p, Some(mut a))), Some((mut q, Some(mut b)))) = (grid[i][j], next) else { continue };
                // Halved while apart.
                while a.angle(b) > 0.1 {
                    if (p[0] - q[0]).hypot(p[1] - q[1]) < fine {
                        found += 1;
                        break;
                    }
                    let m = [(p[0] + q[0]) / 2.0, (p[1] + q[1]) / 2.0];
                    let Some(Some(c)) = hit(model, m) else { break };
                    if a.angle(c) > b.angle(c) {
                        (q, b) = (m, c);
                    } else {
                        (p, a) = (m, c);
                    }
                }
            }
        }
    }
    (found, missing, points)
}

/// `pts` rounded by `rho` at each corner, its top loop filleted by `r` (larger than some of the
/// rounds): the fillet holds every ball of radius `r` inside the prism, is the ordinary one
/// along the straight edges outside the stretches the ball changes in, the sphere as wide as
/// the round round each tighter round, and smooth, from above, round every corner.
fn check_tight_loop(pts: &[[f64; 2]], rho: &[f64], r: f64, h: f64) {
    let n = pts.len();
    let mut model = Model::new();
    let (prism, edge) = rounded_prism(&mut model, pts, rho, h);
    let out = fillet(&mut model, prism, &[edge], r).unwrap();
    assert!(validity::check(&model, out).is_valid());
    let classifier = PointClassifier::new(&model, out).unwrap();
    let inside = |p: Point| classifier.locate(p) == PointLocation::Inside;
    // The turn at each corner.
    let turn: Vec<f64> = (0..n)
        .map(|i| {
            let (p, c, q) = (pts[(i + n - 1) % n], pts[i], pts[(i + 1) % n]);
            let (a, b) = ((c[1] - p[1]).atan2(c[0] - p[0]), (q[1] - c[1]).atan2(q[0] - c[0]));
            (b - a).rem_euclid(2.0 * PI)
        })
        .collect();
    // The top loop moved in by r: the polygon moved in, sharp at the rounds tighter than r,
    // rounded by rho - r about the same centre at the others.
    let core = inset(pts, r);
    let centre_of = |i: usize| {
        let (p, c, q) = (pts[(i + n - 1) % n], pts[i], pts[(i + 1) % n]);
        let (u, v) = ([p[0] - c[0], p[1] - c[1]], [q[0] - c[0], q[1] - c[1]]);
        let (lu, lv) = (u[0].hypot(u[1]), v[0].hypot(v[1]));
        let b = [u[0] / lu + v[0] / lv, u[1] / lu + v[1] / lv];
        let lb = b[0].hypot(b[1]);
        let s = rho[i] / (0.5 * (PI - turn[i])).sin();
        [c[0] + b[0] / lb * s, c[1] + b[1] / lb * s]
    };
    let mut ring: Vec<[f64; 2]> = Vec::new();
    for i in 0..n {
        if rho[i] < r {
            ring.push(core[i]);
            continue;
        }
        // From square to the side before to square to the side after.
        let (p, c) = (pts[(i + n - 1) % n], pts[i]);
        let a = (c[1] - p[1]).atan2(c[0] - p[0]) - PI / 2.0;
        let o = centre_of(i);
        for k in 0..=16 {
            let t = a + turn[i] * k as f64 / 16.0;
            ring.push([o[0] + (rho[i] - r) * t.cos(), o[1] + (rho[i] - r) * t.sin()]);
        }
    }
    // Every ball of radius r inside: about the loop moved in, at h - r.
    let z = h - r;
    for i in 0..ring.len() {
        let (c, d) = (ring[i], ring[(i + 1) % ring.len()]);
        for k in 0..=8 {
            let t = k as f64 / 8.0;
            let o = Point::new(c[0] + (d[0] - c[0]) * t, c[1] + (d[1] - c[1]) * t, z);
            for (a, e) in [(0.0f64, 0.0f64), (1.0, 0.3), (2.0, 0.6), (3.0, 1.0), (4.0, 1.3), (5.0, 0.2)] {
                let w = Vector::new(a.cos() * e.cos(), a.sin() * e.cos(), e.sin());
                assert!(inside(o + w * (0.98 * r)), "ball about {o:?} not kept");
            }
        }
    }
    for i in 0..n {
        let (p, c, q) = (pts[(i + n - 1) % n], pts[i], pts[(i + 1) % n]);
        let (u, v) = ([p[0] - c[0], p[1] - c[1]], [q[0] - c[0], q[1] - c[1]]);
        let (lu, lv) = (u[0].hypot(u[1]), v[0].hypot(v[1]));
        let (u, v) = ([u[0] / lu, u[1] / lu], [v[0] / lv, v[1] / lv]);
        let half = 0.5 * (PI - turn[i]);
        if rho[i] < r {
            // The round's sphere: as wide as it, centred on its axis rho under the top.
            let b = [u[0] + v[0], u[1] + v[1]];
            let lb = b[0].hypot(b[1]);
            let (a, away) = ([c[0] + b[0] / lb * rho[i] / half.sin(), c[1] + b[1] / lb * rho[i] / half.sin()], [-b[0] / lb, -b[1] / lb]);
            let centre = Point::new(a[0], a[1], h - rho[i]);
            for e in [0.3f64, 0.8, 1.2] {
                let w = Vector::new(away[0] * e.cos(), away[1] * e.cos(), e.sin());
                assert!(inside(centre + w * (0.99 * rho[i])), "short of the round's sphere at {e}");
                assert!(!inside(centre + w * (1.01 * rho[i])), "past the round's sphere at {e}");
            }
            // Smooth round the corner, from above.
            let (found, missing, points) = creases(&mut model, out, c[0], c[1], h, 2.0 * r, rho[i] / 1000.0);
            assert_eq!(found, 0, "creases at corner {i}");
            assert!(missing * 50 <= points, "{missing} of {points} points without a normal at corner {i}");
        }
        // The ordinary fillet halfway along the side on to the next corner, clear of the
        // stretches the ball changes in: its section's corner square less the quarter disc.
        let m = [(c[0] + q[0]) / 2.0, (c[1] + q[1]) / 2.0];
        let inward = [-v[1], v[0]];
        let at = |s: f64| Point::new(m[0] + inward[0] * s, m[1] + inward[1] * s, h - s);
        let k = 1.0 - 1.0 / 2f64.sqrt();
        assert!(!inside(at(0.98 * k * r)), "the ordinary fillet not taken away");
        assert!(inside(at(1.02 * k * r)), "more than the ordinary fillet taken away");
    }
    // Its volume: at least the union of the balls of radius r inside (Steiner's formula over the
    // top slab), at most the prism less the ordinary fillet along each straight edge outside
    // the stretches the ball changes in next to each tighter round (twice the wedged ball's
    // growth, or twice the growth, the more).
    let rounded_area = area(pts) - (0..n).map(|i| rho[i] * rho[i] * ((0.5 * turn[i]).tan() - 0.5 * turn[i])).sum::<f64>();
    let wide = |i: usize| (rho[i] - r).max(0.0);
    let core_area = area(&core) - (0..n).map(|i| wide(i) * wide(i) * ((0.5 * turn[i]).tan() - 0.5 * turn[i])).sum::<f64>();
    let core_length = perimeter(&core) - (0..n).map(|i| wide(i) * 2.0 * ((0.5 * turn[i]).tan() - 0.5 * turn[i])).sum::<f64>();
    let least = rounded_area * z + core_area * r + PI * r * r / 4.0 * core_length + 2.0 * PI / 3.0 * r.powi(3);
    let change = |i: usize| if rho[i] < r { 2.0 * (r - rho[i]) * (0.5 * turn[i]).tan().max(1.0) } else { 0.0 };
    let straight: f64 = (0..n)
        .map(|i| {
            let j = (i + 1) % n;
            let len = (pts[j][0] - pts[i][0]).hypot(pts[j][1] - pts[i][1]) - rho[i] * (0.5 * turn[i]).tan() - rho[j] * (0.5 * turn[j]).tan();
            (len - change(i) - change(j)).max(0.0)
        })
        .sum();
    let most = rounded_area * h - (1.0 - PI / 4.0) * r * r * straight;
    let v = properties::volume(&model, out).mass;
    assert!(v >= least && v <= most, "volume {v} not within [{least}, {most}]");
}

/// A square's rounds (turns of a quarter) tighter than the ball.
#[test]
fn square_top_loop_rounds_tighter_than_the_ball() {
    let w = 7.6;
    check_tight_loop(&[[0.0, 0.0], [w, 0.0], [w, w], [0.0, w]], &[0.412; 4], 0.694, 3.34);
}

/// A hexagon's (turns of a sixth).
#[test]
fn hexagon_top_loop_rounds_tighter_than_the_ball() {
    let pts: Vec<[f64; 2]> = (0..6).map(|k| [5.0 * (k as f64 * PI / 3.0).cos(), 5.0 * (k as f64 * PI / 3.0).sin()]).collect();
    check_tight_loop(&pts, &[0.5; 6], 1.3, 6.0);
}

/// A triangle's (turns of a third): the wedged ball's growth longer than the growth.
#[test]
fn triangle_top_loop_rounds_tighter_than_the_ball() {
    let pts: Vec<[f64; 2]> = (0..3).map(|k| [8.0 * (PI / 2.0 + k as f64 * 2.0 * PI / 3.0).cos(), 8.0 * (PI / 2.0 + k as f64 * 2.0 * PI / 3.0).sin()]).collect();
    check_tight_loop(&pts, &[0.5; 3], 1.2, 5.0);
}

/// A rectangle with two rounds tighter than the ball and two wider, the ball changing only next
/// to the tighter ones; and a narrow one whose short sides are too short for both changes, both
/// shortened alike.
#[test]
fn top_loop_with_some_rounds_tighter_and_short_sides() {
    check_tight_loop(&[[0.0, 0.0], [9.0, 0.0], [9.0, 6.0], [0.0, 6.0]], &[0.3, 2.0, 0.3, 2.0], 1.0, 4.0);
    check_tight_loop(&[[0.0, 0.0], [3.0, 0.0], [3.0, 8.0], [0.0, 8.0]], &[0.4; 4], 1.0, 4.0);
}
