//! Linear algebra and positioning: points, vectors, directions, axes, frames, matrices and
//! transformations, in 3D and (with a `2` suffix) in 2D.

mod axis;
mod bounding_box;
mod bvh;
mod direction;
mod frame;
mod matrix;
mod point;
mod quaternion;
mod transform;
mod vector;

pub use axis::{Axis, Axis2};
pub use bounding_box::{BoundingBox, BoundingBox2};
pub use bvh::{Bvh, BvhNode, BvhOptions, PairSelector, Selector};
pub use direction::{Direction, Direction2};
pub use frame::{Frame, Frame2};
pub use matrix::Matrix3;
pub use point::{Point, Point2};
pub use quaternion::{EulerSequence, Quaternion};
pub use transform::{
    AffineTransform, AffineTransformable, Transform, Transform2, Transformable, Transformable2,
};
pub use vector::{Vector, Vector2};

/// The point of the triangle `abc` closest to `p` (on its plane inside it, else on an edge or
/// at a corner).
#[doc(alias = "BVH_Tools::PointTriangleProjection")]
pub fn closest_point_on_triangle(p: Point, a: Point, b: Point, c: Point) -> Point {
    // Ericson, Real-Time Collision Detection, 5.1.5: by the Voronoi regions of the corners
    // and edges.
    let (ab, ac, ap) = (b - a, c - a, p - a);
    let (d1, d2) = (ab.dot(ap), ac.dot(ap));
    if d1 <= 0.0 && d2 <= 0.0 {
        return a;
    }
    let bp = p - b;
    let (d3, d4) = (ab.dot(bp), ac.dot(bp));
    if d3 >= 0.0 && d4 <= d3 {
        return b;
    }
    let vc = d1 * d4 - d3 * d2;
    if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 {
        return a + ab * (d1 / (d1 - d3));
    }
    let cp = p - c;
    let (d5, d6) = (ab.dot(cp), ac.dot(cp));
    if d6 >= 0.0 && d5 <= d6 {
        return c;
    }
    let vb = d5 * d2 - d1 * d6;
    if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 {
        return a + ac * (d2 / (d2 - d6));
    }
    let va = d3 * d6 - d5 * d4;
    if va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0 {
        return b + (c - b) * ((d4 - d3) / ((d4 - d3) + (d5 - d6)));
    }
    let denom = 1.0 / (va + vb + vc);
    a + ab * (vb * denom) + ac * (vc * denom)
}

/// The area of the triangle `abc` (half the length of the cross product of two sides).
#[doc(alias = "bg_area_of_triangle")]
pub fn triangle_area(a: Point, b: Point, c: Point) -> f64 {
    0.5 * (b - a).cross(c - a).length()
}

/// The point of the segment `ab` closest to `p` (`a` for a segment of no length).
#[doc(alias = "bg_dist_pt3_lseg3")]
pub fn closest_point_on_segment(p: Point, a: Point, b: Point) -> Point {
    let ab = b - a;
    let length_squared = ab.dot(ab);
    if length_squared <= 0.0 {
        return a;
    }
    a + ab * ((p - a).dot(ab) / length_squared).clamp(0.0, 1.0)
}

/// The parameters `(s, t)` in `[0, 1]` of the points `p0 + s·(p1 - p0)` and
/// `q0 + t·(q1 - q0)` of two segments closest to each other (for parallel segments, one
/// such pair).
#[doc(alias = "bg_distsq_lseg3_lseg3")]
pub fn closest_points_on_segments(p0: Point, p1: Point, q0: Point, q1: Point) -> (f64, f64) {
    // Ericson, Real-Time Collision Detection, 5.1.9.
    let (d1, d2, r) = (p1 - p0, q1 - q0, p0 - q0);
    let (a, e, f) = (d1.dot(d1), d2.dot(d2), d2.dot(r));
    if a <= f64::MIN_POSITIVE && e <= f64::MIN_POSITIVE {
        return (0.0, 0.0);
    }
    if a <= f64::MIN_POSITIVE {
        return (0.0, (f / e).clamp(0.0, 1.0));
    }
    let c = d1.dot(r);
    if e <= f64::MIN_POSITIVE {
        return ((-c / a).clamp(0.0, 1.0), 0.0);
    }
    let b = d1.dot(d2);
    let denom = a * e - b * b;
    let mut s = if denom > 1e-15 * a * e { ((b * f - c * e) / denom).clamp(0.0, 1.0) } else { 0.0 };
    let mut t = (b * s + f) / e;
    if t < 0.0 {
        t = 0.0;
        s = (-c / a).clamp(0.0, 1.0);
    } else if t > 1.0 {
        t = 1.0;
        s = ((b - c) / a).clamp(0.0, 1.0);
    }
    (s, t)
}

/// The distance between two (filled) triangles: zero where they cross or touch.
pub fn triangle_distance(t: [Point; 3], u: [Point; 3]) -> f64 {
    // Where they meet, an edge of one passes through the other; else the nearest points
    // are a corner and the other triangle, or a point of an edge of each.
    let pierces = |edge: (Point, Point), tri: [Point; 3]| {
        let n = (tri[1] - tri[0]).cross(tri[2] - tri[0]);
        let (da, db) = (n.dot(edge.0 - tri[0]), n.dot(edge.1 - tri[0]));
        if n.length() <= f64::MIN_POSITIVE || da * db > 0.0 || da == db {
            return false;
        }
        let p = edge.0 + (edge.1 - edge.0) * (da / (da - db));
        closest_point_on_triangle(p, tri[0], tri[1], tri[2]).distance(p) <= 1e-12 * (1.0 + (p - Point::ORIGIN).length())
    };
    let edges = |tri: [Point; 3]| [(tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])];
    if edges(t).iter().any(|&e| pierces(e, u)) || edges(u).iter().any(|&e| pierces(e, t)) {
        return 0.0;
    }
    let mut d = f64::INFINITY;
    for (corners, tri) in [(t, u), (u, t)] {
        for p in corners {
            d = d.min(closest_point_on_triangle(p, tri[0], tri[1], tri[2]).distance(p));
        }
    }
    for (p0, p1) in edges(t) {
        for (q0, q1) in edges(u) {
            let (s, w) = closest_points_on_segments(p0, p1, q0, q1);
            d = d.min((p0 + (p1 - p0) * s).distance(q0 + (q1 - q0) * w));
        }
    }
    d
}

/// Whether two (filled) triangles cross or come within `tolerance` of each other.
#[doc(alias = "bg_tri_tri_isect")]
pub fn triangles_intersect(t: [Point; 3], u: [Point; 3], tolerance: f64) -> bool {
    triangle_distance(t, u) <= tolerance
}
