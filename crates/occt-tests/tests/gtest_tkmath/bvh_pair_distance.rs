//! BVH_PairDistance_Test.cxx (BVH_PairDistance: `math::Bvh::closest_pair` over triangle meshes)
//!
//! OCCT meshes the shapes with BRepMesh; here the meshes are built directly: a box's faces as
//! two triangles each, and a sphere as latitude and longitude bands with a pole facing the
//! other sphere (where the meshes, inscribed in the spheres, come closest, at the analytic
//! distance). The distance between two triangles is OCCT's test's: from corners and points
//! along the edges of each to the other.

use kentron::math::{closest_point_on_triangle, BoundingBox, Bvh, Direction, Frame, Point, Vector};

type Triangle = [Point; 3];

fn point_triangle_distance(p: Point, t: &Triangle) -> f64 {
    closest_point_on_triangle(p, t[0], t[1], t[2]).distance(p)
}

/// As OCCT's `triangleTriangleSqDistance` (not squared): corners, then 100 points per edge.
fn triangle_distance(a: &Triangle, b: &Triangle) -> f64 {
    let mut best = f64::INFINITY;
    for (t, other) in [(a, b), (b, a)] {
        for &p in t {
            best = best.min(point_triangle_distance(p, other));
        }
        for k in 0..3 {
            let (p, q) = (t[k], t[(k + 1) % 3]);
            for i in 1..100 {
                best = best.min(point_triangle_distance(p + (q - p) * (i as f64 / 100.0), other));
            }
        }
    }
    best
}

/// A box's faces, two triangles each.
fn box_mesh(corner: Point, size: f64) -> Vec<Triangle> {
    let p = |i: usize, j: usize, k: usize| corner + Vector::new(i as f64, j as f64, k as f64) * size;
    let quads = [
        [p(0, 0, 0), p(1, 0, 0), p(1, 1, 0), p(0, 1, 0)],
        [p(0, 0, 1), p(1, 0, 1), p(1, 1, 1), p(0, 1, 1)],
        [p(0, 0, 0), p(1, 0, 0), p(1, 0, 1), p(0, 0, 1)],
        [p(0, 1, 0), p(1, 1, 0), p(1, 1, 1), p(0, 1, 1)],
        [p(0, 0, 0), p(0, 1, 0), p(0, 1, 1), p(0, 0, 1)],
        [p(1, 0, 0), p(1, 1, 0), p(1, 1, 1), p(1, 0, 1)],
    ];
    quads.iter().flat_map(|q| [[q[0], q[1], q[2]], [q[0], q[2], q[3]]]).collect()
}

/// A sphere in bands, its north pole along `pole`.
fn sphere_mesh(center: Point, radius: f64, pole: Direction, bands: usize) -> Vec<Triangle> {
    let frame = Frame::from_z(center, pole);
    let (x, y, z) = (Vector::from(frame.x()), Vector::from(frame.y()), Vector::from(frame.z()));
    let at = |i: usize, j: usize| {
        let (theta, phi) = (std::f64::consts::PI * i as f64 / bands as f64, std::f64::consts::TAU * j as f64 / (2 * bands) as f64);
        center + (z * theta.cos() + (x * phi.cos() + y * phi.sin()) * theta.sin()) * radius
    };
    let mut out = Vec::new();
    for i in 0..bands {
        for j in 0..2 * bands {
            let (a, b, c, d) = (at(i, j), at(i, j + 1), at(i + 1, j + 1), at(i + 1, j));
            if i > 0 {
                out.push([a, b, c]);
            }
            if i + 1 < bands {
                out.push([a, c, d]);
            }
        }
    }
    out
}

fn tree(mesh: &[Triangle]) -> Bvh {
    Bvh::new(&mesh.iter().map(|t| BoundingBox::from_points(*t)).collect::<Vec<_>>())
}

fn mesh_distance(a: &[Triangle], b: &[Triangle]) -> f64 {
    let (ta, tb) = (tree(a), tree(b));
    let (_, _, d) = ta.closest_pair(&tb, |i, j| triangle_distance(&a[i], &b[j])).expect("a closest pair");
    d
}

/// Migrated from tests/lowalgos/bvh/bug30655_2: two 10 mm boxes 10 mm apart.
#[test]
fn bug30655_box_distance() {
    let d = mesh_distance(&box_mesh(Point::ORIGIN, 10.0), &box_mesh(Point::new(20.0, 0.0, 0.0), 10.0));
    assert!((d - 10.0).abs() <= 1e-2, "{d}");
}

/// Migrated from tests/lowalgos/bvh/bug30655_3: spheres of radii 100 and 110, centres
/// (0, 0, 0) and (150, 150, 150).
#[test]
fn bug30655_sphere_distance() {
    let towards = Direction::new(1.0, 1.0, 1.0);
    let a = sphere_mesh(Point::ORIGIN, 100.0, towards, 48);
    let b = sphere_mesh(Point::new(150.0, 150.0, 150.0), 110.0, -towards, 48);
    let analytic = (3.0 * 150.0 * 150.0f64).sqrt() - 210.0;
    let d = mesh_distance(&a, &b);
    assert!((d - analytic).abs() <= 1e-2, "{d} vs {analytic}");
}
