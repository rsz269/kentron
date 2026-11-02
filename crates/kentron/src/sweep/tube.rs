//! Tubes between two curves: the surface of a pipe of constant radius along a path, cut
//! lengthwise along two curves lying on it.

use crate::error::{Error, Result};
use crate::geometry::curve::Curve;
use crate::geometry::surface::{BSplineSurface, KnotVector, Surface};
use crate::geometry::ParametricCurve;
use crate::math::{Point, Vector};

/// Sections along the path, and points along each section's arc.
const ALONG: usize = 128;
const ACROSS: usize = 32;

/// The surface of a pipe of `radius` along `path` (over `[first, last]`) between `side1` and
/// `side2`, curves on it: at each point of the path, the arc of the circle of `radius` round
/// it, square to the path, from where `side1` crosses that plane to where `side2` does (the
/// shorter way round). Its first parameter runs along the arcs (0 to 1), its second along
/// the path (over `[first, last]`).
#[doc(alias = "GeomFill_Pipe")]
pub fn tube_between(path: &Curve, first: f64, last: f64, side1: &Curve, side2: &Curve, radius: f64) -> Result<Surface> {
    if !(radius > 0.0) || !(last > first) {
        return Err(Error::InvalidInput("a tube needs a radius and a path".into()));
    }
    // Where a side crosses the plane square to the path at `p` (tangent `t`): from the
    // previous crossing (or, at first, the side's point nearest `p`), Newton's method on the
    // distance from the plane.
    let crossing = |side: &Curve, p: Point, t: Vector, from: f64| -> Result<f64> {
        let mut s = from;
        for _ in 0..50 {
            let (q, d) = side.d1(s);
            let (g, dg) = ((q - p).dot(t), d.dot(t));
            if g.abs() <= 1e-12 * (1.0 + radius) {
                return Ok(s);
            }
            if dg.abs() < 1e-300 {
                break;
            }
            s -= g / dg;
        }
        let (q, _) = side.d1(s);
        if (q - p).dot(t).abs() <= 1e-9 * (1.0 + radius) {
            Ok(s)
        } else {
            Err(Error::Failed("a tube's side does not cross the plane square to its path".into()))
        }
    };
    let start = |side: &Curve, p: Point| -> f64 {
        let (lo, hi) = (side.first_parameter(), side.last_parameter());
        // (An unbounded side looked for near the path.)
        let (lo, hi) = if lo.is_finite() && hi.is_finite() { (lo, hi) } else { (-1e3 * (1.0 + radius), 1e3 * (1.0 + radius)) };
        crate::geometry::extrema::nearest_on_curve(p, side, lo, hi).parameter
    };
    let params: Vec<f64> = (0..=ALONG).map(|i| first + (last - first) * i as f64 / ALONG as f64).collect();
    let (mut s1, mut s2) = (None, None);
    // rows[j][i]: point j along the arc at section i.
    let mut rows: Vec<Vec<Point>> = vec![Vec::with_capacity(params.len()); ACROSS + 1];
    for &v in &params {
        let (p, d) = path.d1(v);
        let t = d * (1.0 / d.length().max(1e-300));
        let a = crossing(side1, p, t, s1.unwrap_or_else(|| start(side1, p)))?;
        let b = crossing(side2, p, t, s2.unwrap_or_else(|| start(side2, p)))?;
        (s1, s2) = (Some(a), Some(b));
        let (q1, q2) = (side1.point(a) - p, side2.point(b) - p);
        let flat = |q: Vector| q - t * q.dot(t);
        let (q1, q2) = (flat(q1), flat(q2));
        if q1.length() < 1e-12 || q2.length() < 1e-12 {
            return Err(Error::Failed("a tube's side meets its path".into()));
        }
        let e1 = q1 * (1.0 / q1.length());
        let e2 = {
            let w = t.cross(e1);
            if w.dot(q2) < 0.0 { -w } else { w }
        };
        let angle = q2.dot(e2).atan2(q2.dot(e1));
        for (j, row) in rows.iter_mut().enumerate() {
            let phi = angle * j as f64 / ACROSS as f64;
            row.push(p + (e1 * phi.cos() + e2 * phi.sin()) * radius);
        }
    }
    // Along the path: the rows fitted alike; across, the poles of each column fitted alike.
    let along = crate::geometry::fit::approximate_alike(&rows, &params, 3, 1e-8)?;
    let across_params: Vec<f64> = (0..=ACROSS).map(|j| j as f64 / ACROSS as f64).collect();
    let pole_count = along[0].poles().len();
    let columns: Vec<Vec<Point>> = (0..pole_count).map(|k| along.iter().map(|c| c.poles()[k]).collect()).collect();
    let across = crate::geometry::fit::approximate_alike(&columns, &across_params, 3, 1e-8)?;
    // (poles[i][k]: across index i, along index k.)
    let poles: Vec<Vec<Point>> = (0..across[0].poles().len()).map(|i| across.iter().map(|c| c.poles()[i]).collect()).collect();
    let u = KnotVector::new(across[0].knots().to_vec(), across[0].multiplicities().to_vec(), across[0].degree());
    let v = KnotVector::new(along[0].knots().to_vec(), along[0].multiplicities().to_vec(), along[0].degree());
    Ok(Surface::BSpline(BSplineSurface::try_new(poles, u, v)?))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::curve::Line;
    use crate::math::Direction;

    #[test]
    fn a_quarter_of_a_cylinder_between_two_lines() {
        // Along the z axis, radius 2, between the lines at x = 2 and y = 2: a quarter of the
        // cylinder, its area a quarter of its side's.
        let path = Curve::Line(Line::new(Point::ORIGIN, Direction::Z));
        let side1 = Curve::Line(Line::new(Point::new(2.0, 0.0, 0.0), Direction::Z));
        let side2 = Curve::Line(Line::new(Point::new(0.0, 2.0, 0.0), Direction::Z));
        let surface = tube_between(&path, 0.0, 10.0, &side1, &side2, 2.0).unwrap();
        use crate::geometry::ParametricSurface;
        let p = surface.point(0.5, 5.0);
        assert!((p - Point::new(2f64.sqrt(), 2f64.sqrt(), 5.0)).length() < 1e-6, "{p:?}");
        let mut model = crate::topology::Model::new();
        let face = crate::make::face_on_surface(&mut model, surface).unwrap();
        let area = crate::properties::area(&model, face).mass;
        assert!((area - 0.25 * 2.0 * std::f64::consts::PI * 2.0 * 10.0).abs() < 1e-4, "{area}");
    }
}
