//! Extending and reshaping bounded curves and surfaces.

use crate::math::{Point, Vector};

use super::convert;
use super::curve::{BSplineCurve, Curve};
use super::knots::Spline;
use super::surface::{BSplineSurface, Surface};
use super::{ParametricCurve, ParametricSurface};

/// A boundary of a surface's parameter rectangle.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum SurfaceSide {
    UMin,
    UMax,
    VMin,
    VMax,
}

/// Extends a bounded curve to reach `point`, at its end (`at_end`) or start, joining with
/// continuity C`continuity` (0 to 3). The other end is unchanged. `None` if the curve is
/// unbounded.
#[doc(alias = "GeomLib::ExtendCurveToPoint")]
pub fn curve_to_point(
    curve: &Curve,
    point: Point,
    continuity: u32,
    at_end: bool,
) -> Option<BSplineCurve> {
    let k = continuity.min(3) as usize;
    let mut bspline = clamped_bspline(curve)?;
    let mut spline = bspline.spline();
    spline.elevate((k + 1).saturating_sub(spline.degree));
    // The piece is parametrized at the curve's speed at the joined end.
    let t = if at_end { bspline.last_parameter() } else { bspline.first_parameter() };
    let (q, d) = bspline.d1(t);
    let speed = d.length();
    let gap = q.distance(point);
    if gap == 0.0 {
        return Some(bspline);
    }
    let delta = if speed > 0.0 { gap / speed } else { 1.0 };
    let w = if at_end { spline.poles[spline.poles.len() - 1][3] } else { spline.poles[0][3] };
    let target = vec![w * point.x, w * point.y, w * point.z, w];
    if at_end {
        spline.append_piece(delta, k, target);
    } else {
        spline.prepend_piece(delta, k, target);
    }
    bspline.set_spline(spline);
    Some(bspline)
}

/// Extends a bounded surface across `side` by about `length` (measured along the
/// surface), joining with continuity C`continuity` (0 to 3). The original surface is
/// reproduced exactly over its old domain. `None` if the surface is unbounded.
#[doc(alias = "GeomLib::ExtendSurfByLength")]
pub fn surface_by_length(
    surface: &Surface,
    length: f64,
    continuity: u32,
    side: SurfaceSide,
) -> Option<BSplineSurface> {
    let (u0, u1, v0, v1) = surface.bounds();
    if ![u0, u1, v0, v1].iter().all(|x| x.is_finite() && x.abs() < 1e50) {
        return None;
    }
    let mut bspline = convert::surface_to_bspline(surface)?;
    bspline.set_u_not_periodic();
    bspline.set_v_not_periodic();
    let k = continuity.min(3) as usize;
    let at_end = matches!(side, SurfaceSide::UMax | SurfaceSide::VMax);
    let extend = |s: &mut BSplineSurface| {
        s.edit_u(|spline| {
            let (lo, hi) = spline.domain();
            spline.segment(lo, hi);
            spline.elevate((k + 1).saturating_sub(spline.degree));
            if at_end {
                extend_rows(spline, length, k);
            } else {
                *spline = spline.reversed();
                extend_rows(spline, length, k);
                *spline = spline.reversed();
                // Reversing back moved the old part on by the piece's length.
                let shift = spline.domain().1 - hi;
                spline.flat.iter_mut().for_each(|x| *x -= shift);
            }
        });
    };
    match side {
        SurfaceSide::UMin | SurfaceSide::UMax => extend(&mut bspline),
        SurfaceSide::VMin | SurfaceSide::VMax => bspline.along_v(|s| extend(s)),
    }
    Some(bspline)
}

/// Appends to a spline whose rows are homogeneous points (one per pole curve) a piece
/// carrying each pole curve about `length` further along its end tangent.
fn extend_rows(spline: &mut Spline, length: f64, k: usize) {
    let derivatives = spline.end_derivatives(1);
    let (end, slope) = (&derivatives[0], &derivatives[1]);
    // Per pole curve: the end point and the unit tangent there.
    let blocks: Vec<(Vector, f64, Vector)> = end
        .chunks(4)
        .zip(slope.chunks(4))
        .map(|(a, d)| {
            let w = a[3];
            let c = Vector::new(a[0] / w, a[1] / w, a[2] / w);
            let dc = (Vector::new(d[0], d[1], d[2]) - c * d[3]) * (1.0 / w);
            (c, w, dc)
        })
        .collect();
    let mean_speed = blocks.iter().map(|b| b.2.length()).sum::<f64>() / blocks.len() as f64;
    let delta = if mean_speed > 0.0 { length / mean_speed } else { length };
    let last: Vec<f64> = blocks
        .iter()
        .flat_map(|&(c, w, dc)| {
            let s = dc.length();
            let t = if s > 0.0 { c + dc * (length / s) } else { c };
            [w * t.x, w * t.y, w * t.z, w]
        })
        .collect();
    spline.append_piece(delta, k, last);
}

/// Deforms a bounded curve so that it starts at `first` with tangent direction
/// `first_tangent` and ends at `last` with tangent direction `last_tangent`.
#[doc(alias = "GeomLib::AdjustExtremity")]
pub fn adjust_ends(
    curve: &Curve,
    first: Point,
    last: Point,
    first_tangent: Vector,
    last_tangent: Vector,
) -> Option<BSplineCurve> {
    let mut bspline = clamped_bspline(curve)?;
    let mut spline = bspline.spline();
    // Enough poles that the two at each end are distinct: at least cubic.
    spline.elevate(3usize.saturating_sub(spline.degree));
    let n = spline.poles.len();
    let point = |row: &[f64]| Point::new(row[0] / row[3], row[1] / row[3], row[2] / row[3]);
    let chord = first.distance(last);
    let reach = |a: Point, b: Point| {
        let d = a.distance(b);
        if d > 0.0 { d } else { chord / 3.0 }
    };
    let (start_reach, end_reach) = (reach(point(&spline.poles[0]), point(&spline.poles[1])), reach(point(&spline.poles[n - 1]), point(&spline.poles[n - 2])));
    let unit = |v: Vector| if v.length() > 0.0 { v * (1.0 / v.length()) } else { v };
    let set = |row: &mut Vec<f64>, p: Point| {
        let w = row[3];
        *row = vec![w * p.x, w * p.y, w * p.z, w];
    };
    set(&mut spline.poles[0], first);
    set(&mut spline.poles[1], first + unit(first_tangent) * start_reach);
    set(&mut spline.poles[n - 1], last);
    set(&mut spline.poles[n - 2], last - unit(last_tangent) * end_reach);
    bspline.set_spline(spline);
    Some(bspline)
}

/// A bounded curve as a B-spline clamped at its ends.
fn clamped_bspline(curve: &Curve) -> Option<BSplineCurve> {
    let (a, b) = (curve.first_parameter(), curve.last_parameter());
    if !(a.is_finite() && b.is_finite() && a.abs() < 1e50 && b.abs() < 1e50) {
        return None;
    }
    let mut bspline = convert::curve_to_bspline(curve)?;
    bspline.set_not_periodic();
    Some(bspline.clamped())
}
