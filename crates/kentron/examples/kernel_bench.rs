//! Times the kernel's heavily used parts, one thread, each workload the least CPU time of
//! several runs, with a checksum of its results (to see whether a build changes them: fused
//! multiply-adds round differently).
//!
//! `cargo run --release -p kentron --no-default-features --example kernel_bench [runs] [names...]`

use std::hint::black_box;

use kentron::geometry::curve::{BSplineCurve, Curve};
use kentron::geometry::extrema;
use kentron::geometry::measure;
use kentron::geometry::surface::{BSplineSurface, KnotVector, Surface};
use kentron::geometry::{ParametricCurve, ParametricSurface};
use kentron::math::{BoundingBox, Bvh, Direction, Frame, Point, Vector};
use kentron::primitives::{Cuboid, Cylinder, Sphere, Torus};
use kentron::topology::{Model, NodeKind, Shape};
use kentron::{blend, boolean, meshing, properties};

/// This thread's CPU time in seconds.
fn cpu_time() -> f64 {
    #[repr(C)]
    struct Timespec {
        seconds: i64,
        nanoseconds: i64,
    }
    extern "C" {
        fn clock_gettime(clock: i32, time: *mut Timespec) -> i32;
    }
    const CLOCK_THREAD_CPUTIME_ID: i32 = 3;
    let mut t = Timespec { seconds: 0, nanoseconds: 0 };
    unsafe { clock_gettime(CLOCK_THREAD_CPUTIME_ID, &mut t) };
    t.seconds as f64 + 1e-9 * t.nanoseconds as f64
}

/// A small deterministic generator (xorshift), for reproducible inputs.
struct Random(u64);

impl Random {
    fn next(&mut self) -> f64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        (self.0 >> 11) as f64 / (1u64 << 53) as f64
    }
}

/// Uniform knots `0..=spans` with clamped ends.
fn knots(spans: usize, degree: usize) -> (Vec<f64>, Vec<usize>) {
    let knots: Vec<f64> = (0..=spans).map(|i| i as f64).collect();
    let mut multiplicities = vec![1; spans + 1];
    multiplicities[0] = degree + 1;
    multiplicities[spans] = degree + 1;
    (knots, multiplicities)
}

/// A wavy degree-`degree` curve through `count` poles, rational if `rational`.
fn wavy_curve(count: usize, degree: usize, rational: bool) -> BSplineCurve {
    let poles: Vec<Point> = (0..count).map(|i| {
        let t = i as f64 * 0.3;
        Point::new(t * 3.0, (2.0 * t).sin() * 5.0, (1.3 * t).cos() * 4.0)
    }).collect();
    let (k, m) = knots(count - degree, degree);
    let curve = BSplineCurve::new(poles, k, m, degree);
    if rational { curve.with_weights((0..count).map(|i| 1.0 + 0.4 * (i as f64 * 0.7).sin().abs()).collect()) } else { curve }
}

/// A wavy surface over `count` × `count` poles.
fn wavy_surface(count: usize, degree: usize) -> BSplineSurface {
    let poles: Vec<Vec<Point>> = (0..count)
        .map(|i| (0..count).map(|j| {
            let (x, y) = (i as f64, j as f64);
            Point::new(x * 2.0, y * 2.0, (0.7 * x).sin() * 3.0 + (0.5 * y).cos() * 2.0)
        }).collect())
        .collect();
    let (k, m) = knots(count - degree, degree);
    BSplineSurface::new(poles, KnotVector::new(k.clone(), m.clone(), degree), KnotVector::new(k, m, degree))
}

fn mix(sum: &mut f64, p: Point) {
    *sum += p.x + 1.7 * p.y + 2.3 * p.z;
}

fn mix_vector(sum: &mut f64, v: Vector) {
    *sum += v.x + 1.7 * v.y + 2.3 * v.z;
}

/// Curve point, first and second derivatives along the whole range.
fn curve_eval(rational: bool, degree: usize) -> f64 {
    let curve = wavy_curve(200, degree, rational);
    let (a, b) = (curve.first_parameter(), curve.last_parameter());
    let mut sum = 0.0;
    let count = 400_000;
    for i in 0..count {
        let t = a + (b - a) * (i as f64 + 0.5) / count as f64;
        let (p, d1, d2) = black_box(&curve).d2(t);
        mix(&mut sum, p);
        mix_vector(&mut sum, d1);
        mix_vector(&mut sum, d2);
    }
    sum
}

/// Surface points and second derivatives on a grid.
fn surface_eval(degree: usize) -> f64 {
    let surface = wavy_surface(24, degree);
    let (u0, u1, v0, v1) = surface.bounds();
    let mut sum = 0.0;
    let n = 400;
    for i in 0..n {
        for j in 0..n {
            let (u, v) = (u0 + (u1 - u0) * (i as f64 + 0.5) / n as f64, v0 + (v1 - v0) * (j as f64 + 0.5) / n as f64);
            let d = black_box(&surface).d2(u, v);
            mix(&mut sum, d.point);
            mix_vector(&mut sum, d.du);
            mix_vector(&mut sum, d.dv);
            mix_vector(&mut sum, d.duv);
        }
    }
    sum
}

/// The `k` grid points nearest each of many points (the start of a projection): distances
/// then a full sort (`way` 0), a selection of the `k` (1), or the distances four at a time
/// with AVX2 from coordinates kept apart, then the selection (2).
fn nearest_grid(way: usize) -> f64 {
    let surface = wavy_surface(16, 3);
    let (u0, u1, v0, v1) = surface.bounds();
    let n = 64;
    let mut grid: Vec<(f64, f64, Point)> = Vec::new();
    for i in 0..=n {
        for j in 0..=n {
            let (u, v) = (u0 + (u1 - u0) * i as f64 / n as f64, v0 + (v1 - v0) * j as f64 / n as f64);
            grid.push((u, v, surface.point(u, v)));
        }
    }
    let (xs, ys, zs): (Vec<f64>, Vec<f64>, Vec<f64>) = (grid.iter().map(|g| g.2.x).collect(), grid.iter().map(|g| g.2.y).collect(), grid.iter().map(|g| g.2.z).collect());
    let mut random = Random(0x9e37_79b9_7f4a_7c15);
    let mut sum = 0.0;
    let k = 8;
    let mut d: Vec<(f64, usize)> = vec![(0.0, 0); grid.len()];
    for _ in 0..4_000 {
        let p = Point::new(30.0 * random.next(), 30.0 * random.next(), 8.0 * random.next() - 4.0);
        match way {
            0 | 1 => {
                for (i, g) in grid.iter().enumerate() {
                    d[i] = (g.2.distance(p), i);
                }
            }
            3 | 4 => {}
            _ => {
                #[cfg(target_arch = "x86_64")]
                unsafe {
                    squared_distances(&xs, &ys, &zs, p, &mut d)
                }
            }
        }
        if way == 4 {
            // The same filter one at a time, from the points as they are.
            let mut best = [(f64::INFINITY, usize::MAX); 8];
            for (i, g) in grid.iter().enumerate() {
                let s = (g.2 - p).length_squared();
                if s < best[7].0 {
                    let mut at = 7;
                    while at > 0 && best[at - 1].0 > s {
                        best[at] = best[at - 1];
                        at -= 1;
                    }
                    best[at] = (s, i);
                }
            }
            for &(_, i) in &best {
                sum += grid[i].0 + grid[i].1;
            }
            continue;
        }
        if way == 3 {
            #[cfg(target_arch = "x86_64")]
            let best = unsafe { nearest_filtered::<8>(&xs, &ys, &zs, p) };
            for &(_, i) in &best {
                sum += grid[i].0 + grid[i].1;
            }
            continue;
        }
        if way == 0 {
            d.sort_by(|a, b| a.0.total_cmp(&b.0));
        } else {
            d.select_nth_unstable_by(k, |a, b| a.0.total_cmp(&b.0));
            d[..k].sort_by(|a, b| a.0.total_cmp(&b.0));
        }
        for &(_, i) in &d[..k] {
            sum += grid[i].0 + grid[i].1;
        }
    }
    sum
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2,fma")]
unsafe fn squared_distances(xs: &[f64], ys: &[f64], zs: &[f64], p: Point, out: &mut [(f64, usize)]) {
    use std::arch::x86_64::*;
    let (px, py, pz) = (_mm256_set1_pd(p.x), _mm256_set1_pd(p.y), _mm256_set1_pd(p.z));
    let mut i = 0;
    let mut lanes = [0.0f64; 4];
    while i + 4 <= xs.len() {
        let dx = _mm256_sub_pd(_mm256_loadu_pd(xs.as_ptr().add(i)), px);
        let dy = _mm256_sub_pd(_mm256_loadu_pd(ys.as_ptr().add(i)), py);
        let dz = _mm256_sub_pd(_mm256_loadu_pd(zs.as_ptr().add(i)), pz);
        let s = _mm256_fmadd_pd(dz, dz, _mm256_fmadd_pd(dy, dy, _mm256_mul_pd(dx, dx)));
        _mm256_storeu_pd(lanes.as_mut_ptr(), s);
        for l in 0..4 {
            out[i + l] = (lanes[l], i + l);
        }
        i += 4;
    }
    for j in i..xs.len() {
        let (dx, dy, dz) = (xs[j] - p.x, ys[j] - p.y, zs[j] - p.z);
        out[j] = (dx * dx + dy * dy + dz * dz, j);
    }
}

/// The `K` nearest by squared distance, sorted: four at a time against the worst kept, the
/// few lanes that beat it (a mask's set bits, found by counting trailing zeros) inserted.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2,fma,bmi1")]
unsafe fn nearest_filtered<const K: usize>(xs: &[f64], ys: &[f64], zs: &[f64], p: Point) -> [(f64, usize); K] {
    use std::arch::x86_64::*;
    let mut best = [(f64::INFINITY, usize::MAX); K];
    let insert = |best: &mut [(f64, usize); K], s: f64, i: usize| {
        let mut at = K - 1;
        while at > 0 && best[at - 1].0 > s {
            best[at] = best[at - 1];
            at -= 1;
        }
        best[at] = (s, i);
    };
    let (px, py, pz) = (_mm256_set1_pd(p.x), _mm256_set1_pd(p.y), _mm256_set1_pd(p.z));
    let mut lanes = [0.0f64; 4];
    let mut i = 0;
    while i + 4 <= xs.len() {
        let dx = _mm256_sub_pd(_mm256_loadu_pd(xs.as_ptr().add(i)), px);
        let dy = _mm256_sub_pd(_mm256_loadu_pd(ys.as_ptr().add(i)), py);
        let dz = _mm256_sub_pd(_mm256_loadu_pd(zs.as_ptr().add(i)), pz);
        let s = _mm256_fmadd_pd(dz, dz, _mm256_fmadd_pd(dy, dy, _mm256_mul_pd(dx, dx)));
        let mut mask = _mm256_movemask_pd(_mm256_cmp_pd(s, _mm256_set1_pd(best[K - 1].0), _CMP_LT_OQ)) as u32;
        if mask != 0 {
            _mm256_storeu_pd(lanes.as_mut_ptr(), s);
            while mask != 0 {
                let l = mask.trailing_zeros() as usize;
                if lanes[l] < best[K - 1].0 {
                    insert(&mut best, lanes[l], i + l);
                }
                mask &= mask - 1;
            }
        }
        i += 4;
    }
    for j in i..xs.len() {
        let (dx, dy, dz) = (xs[j] - p.x, ys[j] - p.y, zs[j] - p.z);
        let s = dx * dx + dy * dy + dz * dz;
        if s < best[K - 1].0 {
            insert(&mut best, s, j);
        }
    }
    best
}

/// Nearest points on a free-form surface of points scattered above it.
fn surface_projection() -> f64 {
    let surface = Surface::BSpline(wavy_surface(16, 3));
    let mut random = Random(0x9e37_79b9_7f4a_7c15);
    let mut sum = 0.0;
    for _ in 0..4_000 {
        let p = Point::new(30.0 * random.next(), 30.0 * random.next(), 8.0 * random.next() - 4.0);
        let (u, v, d) = extrema::nearest_on_surface(p, &surface, None);
        sum += u + v + d;
    }
    sum
}

/// Nearest points on a free-form curve.
fn curve_projection() -> f64 {
    let curve = Curve::BSpline(wavy_curve(60, 3, false));
    let (a, b) = (curve.first_parameter(), curve.last_parameter());
    let mut random = Random(0x2545_f491_4f6c_dd1d);
    let mut sum = 0.0;
    for _ in 0..20_000 {
        let p = Point::new(55.0 * random.next(), 12.0 * random.next() - 6.0, 10.0 * random.next() - 5.0);
        let found = extrema::nearest_on_curve(p, &curve, a, b);
        sum += found.parameter + found.distance;
    }
    sum
}

/// Arc lengths of stretches of a rational curve.
fn arc_length() -> f64 {
    let curve = Curve::BSpline(wavy_curve(80, 4, true));
    let (a, b) = (curve.first_parameter(), curve.last_parameter());
    let mut random = Random(0x1234_5678_9abc_def1);
    let mut sum = 0.0;
    for _ in 0..500 {
        let (x, y) = (a + (b - a) * random.next(), a + (b - a) * random.next());
        sum += measure::length_between(&curve, x.min(y), x.max(y));
    }
    sum
}

/// Building a tree of many boxes and finding those overlapping others.
fn bvh_queries() -> f64 {
    let mut random = Random(0xdead_beef_cafe_f00d);
    let boxes: Vec<BoundingBox> = (0..100_000)
        .map(|_| {
            let c = Point::new(100.0 * random.next(), 100.0 * random.next(), 100.0 * random.next());
            let h = Vector::new(random.next(), random.next(), random.next());
            BoundingBox::new(c - h, c + h)
        })
        .collect();
    let tree = Bvh::new(&boxes);
    let mut sum = 0.0;
    for b in boxes.iter().step_by(2) {
        let found = tree.overlapping(b);
        sum += found.len() as f64 + found.iter().sum::<usize>() as f64 * 1e-6;
    }
    for b in boxes.iter().step_by(20) {
        let p = b.center();
        if let Some((i, d)) = tree.nearest(p + Vector::new(0.5, 0.5, 0.5), |i| boxes[i].center().distance(p)) {
            sum += i as f64 * 1e-6 + d;
        }
    }
    sum
}

/// A plate cut by a grid of pins (analytic faces).
fn boolean_analytic() -> f64 {
    let mut model = Model::new();
    let plate = Cuboid::from_corners(Point::new(0.0, 0.0, 0.0), Point::new(100.0, 100.0, 2.0)).build(&mut model).unwrap();
    let pins = model.new_empty(NodeKind::Compound);
    let k = 10;
    let pitch = 100.0 / k as f64;
    for i in 0..k {
        for j in 0..k {
            let center = Point::new(pitch * (i as f64 + 0.5), pitch * (j as f64 + 0.5), -1.0);
            let pin = Cylinder::new(0.3 * pitch, 4.0).with_frame(Frame::from_z(center, Direction::Z)).build(&mut model).unwrap();
            model.add(pins, pin).unwrap();
        }
    }
    let result = boolean::difference(&mut model, plate, pins).unwrap();
    properties::volume(&model, result).mass
}

/// A closed solid lofted through wavy closed B-spline sections around `center`.
fn blob(model: &mut Model, center: Point, radius: f64, phase: f64) -> Shape {
    let sections: Vec<Shape> = (0..5)
        .map(|k| {
            let z = center.z + 6.0 * k as f64;
            let count = 12;
            let poles: Vec<Point> = (0..count)
                .map(|i| {
                    let a = std::f64::consts::TAU * i as f64 / count as f64;
                    let r = radius * (1.0 + 0.15 * (3.0 * a + phase + k as f64).sin()) * (1.0 - 0.02 * (k as f64 - 2.0).powi(2));
                    Point::new(center.x + r * a.cos(), center.y + r * a.sin(), z)
                })
                .collect();
            let curve = BSplineCurve::new_periodic(poles, (0..=count).map(|i| i as f64).collect(), vec![1; count + 1], 3);
            let edge = kentron::make::edge_from_curve(model, Curve::BSpline(curve)).unwrap();
            kentron::make::wire(model, &[edge]).unwrap()
        })
        .collect();
    kentron::loft::loft(model, &sections, &kentron::loft::LoftOptions { solid: true, ..Default::default() }).unwrap().shape
}

/// Free-form solids (lofted through B-spline sections) fused and cut.
fn boolean_freeform() -> f64 {
    let mut model = Model::new();
    let a = blob(&mut model, Point::new(0.0, 0.0, 0.0), 10.0, 0.0);
    let b = blob(&mut model, Point::new(7.0, 3.0, 4.0), 8.0, 1.0);
    let fused = boolean::union(&mut model, a, b).unwrap();
    let cut = boolean::difference(&mut model, a, b).unwrap();
    properties::volume(&model, fused).mass + properties::volume(&model, cut).mass
}

/// Fillets: every edge of a box with a hole.
fn fillet() -> f64 {
    let mut model = Model::new();
    let block = Cuboid::new(40.0, 30.0, 20.0).build(&mut model).unwrap();
    let pin = Cylinder::new(6.0, 40.0).with_frame(Frame::from_z(Point::new(20.0, 15.0, -10.0), Direction::Z)).build(&mut model).unwrap();
    let holed = boolean::difference(&mut model, block, pin).unwrap();
    // (Not the hole's seam, where the cylinder meets itself.)
    let edges: Vec<Shape> = model
        .explore(holed, NodeKind::Edge)
        .filter(|&e| {
            let (curve, first, last) = model.edge_curve(e).unwrap();
            let p = curve.point(0.5 * (first + last));
            !matches!(curve, Curve::Line(_)) || (p.x - 20.0).hypot(p.y - 15.0) > 6.5
        })
        .collect();
    let rounded = blend::fillet(&mut model, holed, &edges, 2.0).unwrap();
    properties::volume(&model, rounded).mass
}

/// Meshing curved faces finely.
fn mesh() -> f64 {
    let mut model = Model::new();
    let torus = Torus::new(10.0, 3.0).build(&mut model).unwrap();
    let sphere = Sphere::new(8.0).build(&mut model).unwrap();
    let lofted = blob(&mut model, Point::new(0.0, 0.0, 0.0), 10.0, 0.0);
    let mut sum = 0.0;
    for shape in [torus, sphere, lofted] {
        let (points, triangles) = meshing::triangulate(&mut model, shape, meshing::MeshOptions { deflection: 2e-3, angle: 0.3 }).unwrap();
        sum += triangles.len() as f64;
        for p in points.iter().step_by(7) {
            mix(&mut sum, *p);
        }
    }
    sum
}

type Workload = (&'static str, fn() -> f64);

fn main() {
    let workloads: Vec<Workload> = vec![
        ("curve_eval_cubic", || curve_eval(false, 3)),
        ("curve_eval_rational5", || curve_eval(true, 5)),
        ("surface_eval_bicubic", || surface_eval(3)),
        ("surface_eval_deg5", || surface_eval(5)),
        ("nearest_grid_sort", || nearest_grid(0)),
        ("nearest_grid_select", || nearest_grid(1)),
        ("nearest_grid_avx", || nearest_grid(2)),
        ("nearest_grid_filter", || nearest_grid(3)),
        ("nearest_grid_filter_scalar", || nearest_grid(4)),
        ("surface_projection", surface_projection),
        ("curve_projection", curve_projection),
        ("arc_length", arc_length),
        ("bvh_queries", bvh_queries),
        ("boolean_analytic", boolean_analytic),
        ("boolean_freeform", boolean_freeform),
        ("fillet", fillet),
        ("mesh", mesh),
    ];
    let mut args = std::env::args().skip(1);
    let runs: usize = args.next().map(|a| a.parse().expect("a run count")).unwrap_or(5);
    let names: Vec<String> = args.collect();
    for (name, work) in workloads {
        if !names.is_empty() && !names.iter().any(|n| name.contains(n.as_str())) {
            continue;
        }
        let mut best = f64::INFINITY;
        let mut checksum = 0.0;
        for _ in 0..runs {
            let start = cpu_time();
            checksum = black_box(work());
            best = best.min(cpu_time() - start);
        }
        println!("{name:<22} {:>10.2} ms  {checksum:.12e}", best * 1e3);
    }
}
