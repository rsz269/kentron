use std::f64::consts::PI;

use super::*;
use crate::math::{Direction, Frame, Transform};
use crate::primitives::{Cuboid, Cylinder, Sphere, Torus};
use crate::{properties, validity};

fn measure(model: &Model, shape: Shape) -> (f64, f64) {
    (properties::volume(model, shape).mass, properties::area(model, shape).mass)
}

fn assert_close(what: &str, value: f64, expected: f64, relative: f64) {
    assert!((value - expected).abs() <= relative * expected.abs(), "{what} {value} against {expected} ({:e} off)", (value - expected) / expected);
}

fn check(model: &Model, shape: Shape) {
    let report = validity::check(model, shape);
    assert!(report.is_valid(), "{report:?}");
}

/// A ball slid is a capsule: a cylinder between two half balls.
#[test]
fn a_ball_slid_is_a_capsule() {
    let mut model = Model::new();
    let ball = Sphere::new(1.5).with_frame(Frame::from_z(Point::new(1.0, 2.0, 3.0), Direction::Z)).build(&mut model).unwrap();
    let by = Vector::new(3.0, -2.0, 4.0);
    let swept = solid(&mut model, ball, &RigidMotion::Slide { by }).unwrap();
    let (volume, area) = measure(&model, swept);
    let (r, l) = (1.5f64, by.length());
    assert_close("volume", volume, 4.0 / 3.0 * PI * r.powi(3) + PI * r * r * l, 1e-8);
    assert_close("area", area, 4.0 * PI * r * r + 2.0 * PI * r * l, 1e-8);
    check(&model, swept);
}

/// A box slid aslant fills the box and, for each axis, its face square to it times the
/// slide's reach along it.
#[test]
fn a_box_slid_aslant() {
    let mut model = Model::new();
    let (a, b, c) = (2.0, 3.0, 1.5);
    let block = Cuboid::new(a, b, c).build(&mut model).unwrap();
    let by = Vector::new(1.2, -0.7, 2.5);
    let swept = solid(&mut model, block, &RigidMotion::Slide { by }).unwrap();
    let (volume, _) = measure(&model, swept);
    assert_close("volume", volume, a * b * c + by.x.abs() * b * c + by.y.abs() * a * c + by.z.abs() * a * b, 1e-10);
    check(&model, swept);
}

/// A torus slid along its axis: itself and, for the reach, the ring between its outer and
/// inner equators.
#[test]
fn a_torus_slid_along_its_axis() {
    let mut model = Model::new();
    let (big, small) = (3.0, 1.0);
    let axis = Direction::new(0.3, -0.2, 0.9);
    let ring = Torus::new(big, small).with_frame(Frame::from_z(Point::new(0.5, 0.0, -1.0), axis)).build(&mut model).unwrap();
    let length = 5.0;
    let swept = solid(&mut model, ring, &RigidMotion::Slide { by: Vector::from(axis) * length }).unwrap();
    let (volume, area) = measure(&model, swept);
    assert_close("volume", volume, 2.0 * PI * PI * big * small * small + length * 4.0 * PI * big * small, 1e-8);
    assert_close("area", area, 4.0 * PI * PI * big * small + 4.0 * PI * big * length, 1e-8);
    check(&model, swept);
}

/// A drum turned about a line beside it and parallel to its axis sweeps the annular sector
/// between its nearest and farthest circles, and itself at each end.
#[test]
fn a_drum_turned_about_a_parallel_line() {
    let mut model = Model::new();
    let (r, h, big) = (1.0, 2.0, 4.0);
    let drum = Cylinder::new(r, h).with_frame(Frame::from_z(Point::new(big, 0.0, 0.0), Direction::Z)).build(&mut model).unwrap();
    let angle = 1.3;
    let swept = solid(&mut model, drum, &RigidMotion::Turn { axis: Axis::new(Point::ORIGIN, Direction::Z), from: -0.4, to: angle - 0.4 }).unwrap();
    let (volume, area) = measure(&model, swept);
    let section = PI * r * r + 2.0 * big * r * angle;
    assert_close("volume", volume, section * h, 1e-8);
    assert_close("area", area, 2.0 * section + h * (2.0 * PI * r + 2.0 * big * angle), 1e-8);
    check(&model, swept);
}

/// A ball screwed sweeps a tube round the helix its centre runs, closed by half balls:
/// `πr²·L + 4/3·πr³` (its turns apart, its curvature gentler than the ball).
#[test]
fn a_ball_screwed_is_a_helical_tube() {
    let mut model = Model::new();
    let (r, big, lead) = (0.8, 4.0, 3.0);
    let ball = Sphere::new(r).with_frame(Frame::from_z(Point::new(big, 0.0, 0.0), Direction::new(0.2, 0.3, 1.0))).build(&mut model).unwrap();
    let (from, to) = (-0.5, 4.0);
    let swept = solid(&mut model, ball, &RigidMotion::Screw { axis: Axis::new(Point::ORIGIN, Direction::Z), lead, from, to }).unwrap();
    let (volume, area) = measure(&model, swept);
    let length = (to - from) * (big * big + (lead / TAU).powi(2)).sqrt();
    assert_close("volume", volume, PI * r * r * length + 4.0 / 3.0 * PI * r.powi(3), 1e-7);
    assert_close("area", area, 2.0 * PI * r * length + 4.0 * PI * r * r, 1e-7);
    check(&model, swept);
}

/// A ball carried straight, ever faster: still a capsule (the motion given by its places).
#[test]
fn a_ball_carried_straight_unevenly_is_a_capsule() {
    let mut model = Model::new();
    let ball = Sphere::new(1.0).build(&mut model).unwrap();
    let by = Vector::new(2.0, 1.0, -2.0);
    let carried = RigidMotion::Carried(Carriage(std::sync::Arc::new(move |t: f64| Transform::translation(by * (t * t)))));
    let swept = solid(&mut model, ball, &carried).unwrap();
    let (volume, area) = measure(&model, swept);
    assert_close("volume", volume, 4.0 / 3.0 * PI + PI * 3.0, 1e-8);
    assert_close("area", area, 4.0 * PI + 2.0 * PI * 3.0, 1e-8);
    check(&model, swept);
}

/// A bar carried by a screw given by its places sweeps what the screw sweeps it through.
#[test]
fn a_screw_given_by_its_places_sweeps_as_the_screw() {
    let axis = Axis::new(Point::ORIGIN, Direction::Z);
    let (lead, from, to) = (5.0, -0.3, 2.9);
    let volumes: Vec<f64> = (0..2)
        .map(|k| {
            let mut model = Model::new();
            let bar = Cuboid::new(3.0, 1.0, 1.5).with_frame(Frame::from_z(Point::new(4.0, -0.5, 0.0), Direction::Z)).build(&mut model).unwrap();
            let screw = RigidMotion::Screw { axis, lead, from, to };
            let motion = if k == 0 { screw } else { RigidMotion::Carried(Carriage(std::sync::Arc::new(move |t: f64| screw.at(t)))) };
            let swept = solid(&mut model, bar, &motion).unwrap_or_else(|e| panic!("motion {k}: {e:?}"));
            check(&model, swept);
            measure(&model, swept).0
        })
        .collect();
    assert_close("volume", volumes[1], volumes[0], 1e-8);
}

/// A box carried round a circle, turning with it, runs along its top and bottom beside faces
/// it carries out: refused; as is a ring whose contact moves over it (tumbling about one line
/// while carried along another).
#[test]
fn motions_beyond_the_product_are_refused() {
    let mut model = Model::new();
    let key = Cuboid::new(3.0, 1.6, 4.0).build(&mut model).unwrap();
    let round = RigidMotion::Carried(Carriage(std::sync::Arc::new(|t: f64| Transform::rotation(Axis::new(Point::new(0.0, 10.0, 0.0), Direction::Z), t))));
    assert!(solid(&mut model, key, &round).is_err());
    let torus = Torus::new(3.0, 1.0).build(&mut model).unwrap();
    let tumble = RigidMotion::Carried(Carriage(std::sync::Arc::new(|t: f64| {
        Transform::translation(Vector::new(0.0, 10.0 * t, 0.0)) * Transform::rotation(Axis::new(Point::ORIGIN, Direction::X), 1.2 * t)
    })));
    assert!(solid(&mut model, torus, &tumble).is_err());
}

/// A torus leaning 35° pushed straight down 8, against the volume its slide covers summed
/// column by column (each vertical line meets the torus in intervals, each reaching 8 lower),
/// on a grid fine enough for 1e-3.
#[test]
fn a_leaning_torus_pushed_down_against_its_columns() {
    let mut model = Model::new();
    let lean = 35f64.to_radians();
    let a = Vector::new(lean.sin(), 0.0, lean.cos());
    let ring = Torus::new(3.0, 1.0).with_frame(Frame::from_z(Point::new(0.0, 0.0, 4.0), Direction::from_vector(a))).build(&mut model).unwrap();
    let swept = solid(&mut model, ring, &RigidMotion::Slide { by: Vector::new(0.0, 0.0, -8.0) }).unwrap();
    check(&model, swept);
    let (volume, _) = measure(&model, swept);
    // (The torus as (|q⊥| − 3)² + (q·a)² ≤ 1, q from its centre.)
    let inside = |p: Point| {
        let q = p - Point::new(0.0, 0.0, 4.0);
        let along = q.dot(a);
        let across = (q - a * along).length();
        (across - 3.0).powi(2) + along * along - 1.0
    };
    let n = 300;
    let h = 9.0 / n as f64;
    let mut total = 0.0;
    for i in 0..n {
        for j in 0..n {
            let (x, y) = (-4.5 + h * (i as f64 + 0.5), -4.5 + h * (j as f64 + 0.5));
            let f = |z: f64| inside(Point::new(x, y, z));
            // The column's crossings, by sign changes on fine steps, each narrowed by bisection.
            let steps = 300;
            let z = |k: usize| -0.5 + 9.0 * k as f64 / steps as f64;
            let mut crossings = Vec::new();
            for k in 0..steps {
                let f0 = f(z(k));
                if (f0 < 0.0) != (f(z(k + 1)) < 0.0) {
                    let (mut l, mut r) = (z(k), z(k + 1));
                    for _ in 0..50 {
                        let m = 0.5 * (l + r);
                        if (f(m) < 0.0) == (f0 < 0.0) { l = m } else { r = m }
                    }
                    crossings.push(0.5 * (l + r));
                }
            }
            let mut covered: Vec<(f64, f64)> = crossings.chunks(2).filter(|c| c.len() == 2).map(|c| (c[0] - 8.0, c[1])).collect();
            covered.sort_by(|p, q| p.0.total_cmp(&q.0));
            let mut end = f64::NEG_INFINITY;
            for (lo, hi) in covered {
                let start = lo.max(end);
                if hi > start {
                    total += hi - start;
                }
                end = end.max(hi);
            }
        }
    }
    assert_close("volume", volume, total * h * h, 1e-3);
}

/// A torus leaning 65° pushed straight down: past 45° its characteristic curves (where its
/// surface runs upright) turn past the slide's direction, and the sheets they sweep fold over
/// one another. Refused as soon as the fold is found on the curves, before the product is made
/// and split where it crosses itself (which found no piece, after 1 s here and 37 s for a
/// spline bar screwed through a rod).
#[test]
fn a_fold_is_refused_before_the_product_is_made() {
    let mut model = Model::new();
    let lean = 65f64.to_radians();
    let ring = Torus::new(3.0, 1.0).with_frame(Frame::from_z(Point::new(0.0, 0.0, 4.0), Direction::new(lean.sin(), 0.0, lean.cos()))).build(&mut model).unwrap();
    let refused = solid(&mut model, ring, &RigidMotion::Slide { by: Vector::new(0.0, 0.0, -8.0) });
    let Err(Error::Failed(why)) = refused else { panic!("a fold refused, not {refused:?}") };
    assert!(why.contains("folding over itself"), "{why}");
}

/// A torus slid square to its axis, far enough that its hole is passed over: it grazes the
/// motion along its two rims and two of its tube's circles, which cross one another, and its
/// near side and then its far side pass over its hole (the product runs into itself there).
/// The swept solid is, along each line in the motion's direction, the torus's span there and
/// the slide: `π(2rR² + πRr² + 4r³/3)` (the torus with its hole filled across the motion)
/// and the slide times the torus's shadow `4Rr + πr²`.
#[test]
fn a_torus_slid_square_to_its_axis_past_its_hole() {
    let mut model = Model::new();
    let (big, small, length) = (6.0, 2.0, 40.0);
    let ring = Torus::new(big, small).with_frame(Frame::from_z(Point::new(11.5, 0.0, -9.0), Direction::Y)).build(&mut model).unwrap();
    let swept = solid(&mut model, ring, &RigidMotion::Slide { by: Vector::new(0.0, 0.0, length) }).unwrap();
    let (volume, _) = measure(&model, swept);
    let filled = PI * (2.0 * small * big * big + PI * big * small * small + 4.0 * small.powi(3) / 3.0);
    assert_close("volume", volume, filled + length * (4.0 * big * small + PI * small * small), 1e-8);
    check(&model, swept);
}

/// A square bar screwed about a line it points at, rising faster as it turns (no screw
/// carries it), comes back over where it started: the leading half of its far end, turned
/// part of a turn, lands inside its start before it has risen past it. The product winds
/// twice round what is covered twice, and is split there. Against the volume summed column by
/// column: each line parallel to the axis meets the bar, turned by θ, for θ in intervals (where
/// its turned point lies in the bar's section), each holding the heights the bar spans risen
/// from the interval's start to its end; the line holds their union (tools/brute.py `turned`
/// the same, finer, 8841.952 ± 0.001). The area against its pieces, found apart: the start's
/// underside and the end's top (2 × 128); the sheets the bar's level edges sweep, whole
/// (∫∫ √(h′² + x²) dx dθ: 2 × 1523.9247 + 2 × 800.2243); the cylinders r = 4 (60π) and
/// r = √208 (its corners' ribbons' union, 1135.0351); the sides of the bar at the start and
/// the end and the sheets its far end's edges sweep, where no place covers them (by the
/// columns, 75.0907 and 54.0670).
#[test]
fn a_bar_screwed_at_a_growing_lead_over_its_own_start() {
    use std::f64::consts::TAU;
    let mut model = Model::new();
    let bar = Cuboid::new(8.0, 16.0, 2.0).with_frame(Frame::from_z(Point::new(4.0, -8.0, -1.0), Direction::Z)).build(&mut model).unwrap();
    let (from, to) = ((-90f64).to_radians(), 1260f64.to_radians());
    let rise = |theta: f64| {
        let s = theta / TAU;
        6.0 * s + 0.3 * s * s
    };
    let place = move |t: f64| {
        let theta = from + (to - from) * t;
        Transform::translation(Vector::new(0.0, 0.0, rise(theta))) * Transform::rotation(Axis::new(Point::ORIGIN, Direction::Z), theta)
    };
    let swept = solid(&mut model, bar, &RigidMotion::Carried(Carriage(std::sync::Arc::new(place)))).unwrap();
    check(&model, swept);
    let (volume, area) = measure(&model, swept);
    // The angles at which the circle of radius r lies in the section 4 ≤ x ≤ 12, |y| ≤ 8:
    // between where it crosses the section's sides.
    let section = [(4.0, -8.0), (12.0, -8.0), (12.0, 8.0), (4.0, 8.0)];
    let inside = |x: f64, y: f64| (4.0..=12.0).contains(&x) && (-8.0..=8.0).contains(&y);
    let arcs = |r: f64| -> Vec<(f64, f64)> {
        let mut cuts = vec![-PI, PI];
        for k in 0..4 {
            let ((ax, ay), (bx, by)) = (section[k], section[(k + 1) % 4]);
            let (dx, dy) = (bx - ax, by - ay);
            let (a, b, c) = (dx * dx + dy * dy, 2.0 * (ax * dx + ay * dy), ax * ax + ay * ay - r * r);
            let disc = b * b - 4.0 * a * c;
            if disc >= 0.0 {
                for s in [(-b - disc.sqrt()) / (2.0 * a), (-b + disc.sqrt()) / (2.0 * a)] {
                    if (0.0..=1.0).contains(&s) {
                        cuts.push((ay + s * dy).atan2(ax + s * dx));
                    }
                }
            }
        }
        cuts.sort_by(f64::total_cmp);
        cuts.windows(2).filter(|w| w[1] > w[0] && inside(r * (0.5 * (w[0] + w[1])).cos(), r * (0.5 * (w[0] + w[1])).sin())).map(|w| (w[0], w[1])).collect()
    };
    // Across r, Gauss–Legendre between the radii where the arcs change (r = 4, √80, 12,
    // √208), on r = a + (b − a)(1 − cos πs)/2, as the arcs' ends move as the square root of
    // the distance from them; around, the midpoint rule.
    let (nodes, weights) = (
        [-0.9815606342467192, -0.9041172563704749, -0.7699026741943047, -0.5873179542866175, -0.3678314989981802, -0.1252334085114689, 0.1252334085114689, 0.3678314989981802, 0.5873179542866175, 0.7699026741943047, 0.9041172563704749, 0.9815606342467192],
        [0.0471753363865118, 0.1069393259953184, 0.1600783285433462, 0.2031674267230659, 0.2334925365383548, 0.2491470458134028, 0.2491470458134028, 0.2334925365383548, 0.2031674267230659, 0.1600783285433462, 0.1069393259953184, 0.0471753363865118],
    );
    let breaks = [4.0, 80f64.sqrt(), 12.0, 208f64.sqrt()];
    let around = 20000;
    let mut total = 0.0;
    for w in breaks.windows(2) {
        for sub in 0..4 {
            // (Each stretch in four, twelve points each.)
            let (a, b) = (w[0], w[1]);
            for (x, wx) in nodes.iter().zip(weights) {
                let s = (sub as f64 + 0.5 * (x + 1.0)) / 4.0;
                let r = a + (b - a) * 0.5 * (1.0 - (PI * s).cos());
                let dr = (b - a) * 0.5 * PI * (PI * s).sin() * 0.5 * wx / 4.0;
                let arcs = arcs(r);
                let mut length = 0.0;
                for i in 0..around {
                    let phi = TAU * (i as f64 + 0.5) / around as f64;
                    let mut held: Vec<(f64, f64)> = Vec::new();
                    for &(lo, hi) in &arcs {
                        for k in -2..6 {
                            let turn = TAU * k as f64;
                            let (ta, tb) = ((phi - hi + turn).clamp(from, to), (phi - lo + turn).clamp(from, to));
                            if tb > ta {
                                held.push((rise(ta) - 1.0, rise(tb) + 1.0));
                            }
                        }
                    }
                    held.sort_by(|p, q| p.0.total_cmp(&q.0));
                    let mut end = f64::NEG_INFINITY;
                    for (lo, hi) in held {
                        let start = lo.max(end);
                        if hi > start {
                            length += hi - start;
                        }
                        end = end.max(hi);
                    }
                }
                total += dr * r * length * TAU / around as f64;
            }
        }
    }
    // (Unsplit, the product's is 9635.86: what is covered twice counted twice.)
    assert_close("volume", volume, total, 1e-6);
    assert_close("area", area, 6356.9864, 1e-6);
}

/// A drum turned partway about a parallel axis it does not meet, then cut by a box holding
/// part of its path: the faces its rims sweep are planes and its sides' sweeps cylinders (not
/// surfaces of revolution, a rim's singular where it runs along the turn), and the boolean
/// takes the sector of the annulus inside the box and the half drum at the end.
#[test]
fn a_drum_turned_partway_is_cut_by_a_box() {
    use crate::geometry::surface::Surface;
    use crate::topology::NodeKind;
    let mut model = Model::new();
    let drum = Cylinder::new(1.0, 1.0).with_frame(Frame::from_z(Point::new(10.0, 0.0, -0.5), Direction::Z)).build(&mut model).unwrap();
    let turn = RigidMotion::Turn { axis: crate::math::Axis::Z, from: (-140f64).to_radians(), to: 30f64.to_radians() };
    let swept = solid(&mut model, drum, &turn).unwrap();
    for face in model.unique_subshapes(swept, NodeKind::Face) {
        assert!(matches!(model.surface(face), Some(Surface::Plane(_) | Surface::Cylinder(_))), "{:?}", model.surface(face));
    }
    let sector = |degrees: f64| PI * (11.0f64.powi(2) - 9.0f64.powi(2)) * degrees / 360.0;
    let (volume, _) = measure(&model, swept);
    assert_close("swept volume", volume, sector(170.0) + PI, 1e-9);
    check(&model, swept);
    let block = Cuboid::new(12.0, 12.0, 2.0).with_frame(Frame::from_z(Point::new(0.0, 0.0, -1.0), Direction::Z)).build(&mut model).unwrap();
    let inside = sector(30.0) + 0.5 * PI;
    let common = crate::boolean::intersection(&mut model, block, swept).unwrap();
    assert_close("common", measure(&model, common).0, inside, 1e-9);
    let cut = crate::boolean::difference(&mut model, block, swept).unwrap();
    assert_close("cut", measure(&model, cut).0, 288.0 - inside, 1e-9);
    check(&model, cut);
}

/// The sector of the ring about the z axis between radii 8 and 12, |z| ≤ 1, polar angle 0°
/// to 60°.
fn ring_sector(model: &mut Model) -> Shape {
    let q = [Point::new(8.0, 0.0, -1.0), Point::new(12.0, 0.0, -1.0), Point::new(12.0, 0.0, 1.0), Point::new(8.0, 0.0, 1.0)];
    let edges: Vec<Shape> = (0..4).map(|i| crate::make::edge_between(model, q[i], q[(i + 1) % 4]).unwrap()).collect();
    let wire = crate::make::wire(model, &edges).unwrap();
    let face = crate::make::planar_face(model, wire).unwrap();
    crate::sweep::revolve(model, face, crate::math::Axis::Z, PI / 3.0).unwrap().shape
}

/// A pin (radius 0.5 along x from 6 to 11) turned about z through the whole of a ring
/// sector: the rim of its far end turns on the sphere about the origin through it, and both
/// halves of the rim sweep that sphere (each its own points: as surfaces of revolution of the
/// circle they would sweep the same points twice, the product not split there). In each
/// meridian half-plane the cut is r from 8 out to √(121.25 − z²), |z| ≤ 1/2: the sector less
/// π/3 · (57 + 1/6) / 2.
#[test]
fn a_pin_turned_through_a_ring_sector() {
    let mut model = Model::new();
    let pin = Cylinder::new(0.5, 5.0).with_frame(Frame::from_z(Point::new(6.0, 0.0, 0.0), Direction::X)).build(&mut model).unwrap();
    let turn = RigidMotion::Turn { axis: crate::math::Axis::Z, from: (-140f64).to_radians(), to: 140f64.to_radians() };
    let swept = solid(&mut model, pin, &turn).unwrap();
    check(&model, swept);
    let block = ring_sector(&mut model);
    let cut = crate::boolean::difference(&mut model, block, swept).unwrap();
    check(&model, cut);
    let sector = PI * (144.0 - 64.0) / 6.0 * 2.0;
    assert_close("cut", measure(&model, cut).0, sector - PI / 3.0 * (57.0 + 1.0 / 6.0) / 2.0, 1e-9);
}

/// A bar (10 ≤ x ≤ 11, |y| ≤ 1/2, |z| ≤ 0.6) turned about z from −140° to 30°, partway into a
/// ring sector: its two far corners sweep one cylinder, the product running into itself
/// between it and the bar's far face (an overlap a few degrees wide at each point, found by
/// stepping along the turn as finely as the body is sampled) and resolved. The reference: in
/// the sector, at radius r, the turn covers 30° and the bar's reach past it,
/// min(asin(1/2r), acos(10/r)), over |z| ≤ 0.6 (by quadrature).
#[test]
fn a_bar_turned_partway_into_a_ring_sector() {
    let mut model = Model::new();
    let bar = Cuboid::new(1.0, 1.0, 1.2).with_frame(Frame::from_z(Point::new(10.0, -0.5, -0.6), Direction::Z)).build(&mut model).unwrap();
    let turn = RigidMotion::Turn { axis: crate::math::Axis::Z, from: (-140f64).to_radians(), to: 30f64.to_radians() };
    let swept = solid(&mut model, bar, &turn).unwrap();
    check(&model, swept);
    let block = ring_sector(&mut model);
    let cut = crate::boolean::difference(&mut model, block, swept).unwrap();
    check(&model, cut);
    assert_close("cut", measure(&model, cut).0, 76.49537612378384, 1e-9);
}

/// A box (2.5 ≤ x ≤ 3.5, |y| ≤ 1/2, 2 ≤ z ≤ 3) screwed about z through two turns of lead 4:
/// its outer upright edges slide along themselves, each sweeping (exactly) the cylinder of
/// radius √12.5, and the product, running into itself between it and the box's outer face,
/// is resolved. The volume against tools/brute.py's `turned` (48.48010, its quadrature
/// converged to 1e-6).
#[test]
fn a_box_screwed_two_turns() {
    use crate::geometry::surface::Surface;
    use crate::topology::NodeKind;
    let mut model = Model::new();
    let block = Cuboid::new(1.0, 1.0, 1.0).with_frame(Frame::from_z(Point::new(2.5, -0.5, 2.0), Direction::Z)).build(&mut model).unwrap();
    let screw = RigidMotion::Screw { axis: crate::math::Axis::Z, lead: 4.0, from: 0.0, to: 4.0 * PI };
    let swept = solid(&mut model, block, &screw).unwrap();
    check(&model, swept);
    let corners = model.unique_subshapes(swept, NodeKind::Face).into_iter().filter(|&f| matches!(model.surface(f), Some(Surface::Cylinder(c)) if (c.radius() - 12.5f64.sqrt()).abs() < 1e-12)).count();
    assert!(corners > 0, "no cylinder swept by the outer edges");
    assert_close("volume", measure(&model, swept).0, 48.48010, 1e-6);
}
