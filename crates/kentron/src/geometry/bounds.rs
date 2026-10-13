//! Axis-aligned bounding boxes of curves and surfaces.
//!
//! The plain functions are fast and may be somewhat loose for free-form geometry (e.g.
//! using control polygons); the `_tight` variants search for the true extent. Both are
//! exact for lines, conics and elementary surfaces, and both enlarge the box by
//! `tolerance`. Parameters at or beyond ±[`INFINITE`](crate::tolerance::INFINITE) are
//! infinite, and give infinite ("open") sides where the geometry is unbounded.
//!
//! Each coordinate's extremes are found where they are: at the ends of the range, at the
//! curve's breaks, and where the coordinate's derivative vanishes (refined from samples),
//! so the boxes are exact wherever that search finds every critical point, which it does
//! for lines, conics, elementary surfaces and smooth free-form pieces alike. An infinite
//! end opens the box on the sides the geometry runs off to there.

use crate::math::{BoundingBox, BoundingBox2, Point, Point2};
use crate::tolerance::INFINITE;

use super::control::ControlPoint;
use super::curve::Curve;
use super::curve2::Curve2;
use super::parametric::{ParametricCurve, ParametricSurface};
use super::surface::Surface;

/// Per coordinate, the smallest and largest value (±infinity for open sides).
type Extent = [[f64; 2]; 3];

const NONE: Extent = [[f64::INFINITY, f64::NEG_INFINITY]; 3];

fn add(extent: &mut Extent, c: [f64; 3]) {
    for k in 0..3 {
        if c[k].is_finite() {
            extent[k][0] = extent[k][0].min(c[k]);
            extent[k][1] = extent[k][1].max(c[k]);
        }
    }
}

fn is_infinite(t: f64) -> bool {
    t.abs() >= 0.5 * INFINITE
}

/// The extent of a curve over `[t1, t2]` (in either order).
fn curve_extent<C>(curve: &C, t1: f64, t2: f64) -> Extent
where
    C: ParametricCurve,
    C::Point: ControlPoint<Vector = C::Vector>,
{
    let (lo, hi) = (t1.min(t2), t1.max(t2));
    let mut extent = NONE;
    let (lo_open, hi_open) = (is_infinite(lo), is_infinite(hi));
    // The open sides at infinite ends: where the derivative points far out.
    let derivative = |t: f64| C::Point::vector_coords(curve.d1(t).1);
    let far = |sign: f64| -> Option<[f64; 3]> {
        [1e6, 1e4, 1e3, 1e2, 30.0].into_iter().map(|t| derivative(sign * t)).find(|d| d.iter().all(|x| x.is_finite()))
    };
    let mut open = [[false; 2]; 3];
    for (is_open, sign) in [(lo_open, -1.0), (hi_open, 1.0)] {
        if !is_open {
            continue;
        }
        let Some(d) = far(sign) else {
            open = [[true; 2]; 3];
            continue;
        };
        for k in 0..3 {
            // Going out to +∞ the coordinate follows the derivative, going to -∞ against it.
            let direction = d[k] * sign;
            if direction > 1e-12 {
                open[k][1] = true;
            } else if direction < -1e-12 {
                open[k][0] = true;
            }
        }
    }
    // The finite part: the range, an infinite end replaced by a window.
    let (a, b) = match (lo_open, hi_open) {
        (false, false) => (lo, hi),
        (true, false) => (hi - 1e3, hi),
        (false, true) => (lo, lo + 1e3),
        (true, true) => (-1e3, 1e3),
    };
    for t in critical_parameters(curve, a, b) {
        add(&mut extent, curve.point(t).coords());
    }
    for k in 0..3 {
        if open[k][0] {
            extent[k][0] = f64::NEG_INFINITY;
        }
        if open[k][1] {
            extent[k][1] = f64::INFINITY;
        }
    }
    extent
}

/// The ends of `[a, b]`, the breaks inside it and every parameter where a coordinate's
/// derivative changes sign (found between samples, refined by halving).
fn critical_parameters<C>(curve: &C, a: f64, b: f64) -> Vec<f64>
where
    C: ParametricCurve,
    C::Point: ControlPoint<Vector = C::Vector>,
{
    let n = match curve.period().filter(|_| curve.is_periodic()) {
        Some(p) => ((b - a) / p * 64.0).ceil().clamp(16.0, 4096.0) as usize,
        None => 96,
    };
    let mut extra = Vec::new();
    let period = curve.period().filter(|_| curve.is_periodic());
    for k in curve.breaks() {
        match period {
            Some(p) => {
                let mut t = k + ((a - k) / p).ceil() * p;
                while t < b {
                    extra.push(t);
                    t += p;
                }
            }
            None if a < k && k < b => extra.push(k),
            None => {}
        }
    }
    critical_1d(&|t| C::Point::vector_coords(curve.d1(t).1), a, b, n, &extra)
}

/// The ends of `[a, b]`, the `extra` parameters and where a coordinate of `derivative`
/// changes sign between `n` even steps (refined by halving).
fn critical_1d(derivative: &dyn Fn(f64) -> [f64; 3], a: f64, b: f64, n: usize, extra: &[f64]) -> Vec<f64> {
    let mut grid: Vec<f64> = (0..=n).map(|i| if i == n { b } else { a + (b - a) * i as f64 / n as f64 }).collect();
    grid.extend_from_slice(extra);
    grid.sort_by(f64::total_cmp);
    grid.dedup();
    let mut out = grid.clone();
    let values: Vec<[f64; 3]> = grid.iter().map(|&t| derivative(t)).collect();
    for i in 0..grid.len().saturating_sub(1) {
        for k in 0..3 {
            let (f0, f1) = (values[i][k], values[i + 1][k]);
            if f0 == 0.0 || f0 * f1 >= 0.0 {
                continue;
            }
            let (mut s, mut t, mut fs) = (grid[i], grid[i + 1], f0);
            for _ in 0..80 {
                let m = 0.5 * (s + t);
                let fm = derivative(m)[k];
                if fm == 0.0 {
                    (s, t) = (m, m);
                    break;
                }
                if (fm < 0.0) == (fs < 0.0) {
                    (s, fs) = (m, fm);
                } else {
                    t = m;
                }
            }
            out.push(0.5 * (s + t));
        }
    }
    out
}

fn to_box(e: Extent, tolerance: f64) -> BoundingBox {
    if e.iter().any(|r| r[0] > r[1]) {
        return BoundingBox::EMPTY;
    }
    BoundingBox::new(Point::new(e[0][0], e[1][0], e[2][0]), Point::new(e[0][1], e[1][1], e[2][1])).enlarged(tolerance)
}

fn to_box2(e: Extent, tolerance: f64) -> BoundingBox2 {
    if e[..2].iter().any(|r| r[0] > r[1]) {
        return BoundingBox2::EMPTY;
    }
    BoundingBox2::new(Point2::new(e[0][0], e[1][0]), Point2::new(e[0][1], e[1][1])).enlarged(tolerance)
}

/// The box of a curve over its whole domain.
#[doc(alias = "BndLib_Add3dCurve::Add")]
#[doc(alias = "GeomBndLib_Curve::Add")]
#[doc(alias = "BndLib::Add")]
pub fn curve(curve: &Curve, tolerance: f64) -> BoundingBox {
    curve_between(curve, curve.first_parameter(), curve.last_parameter(), tolerance)
}

/// The box of the part of a curve between parameters `t1` and `t2`.
#[doc(alias = "BndLib_Add3dCurve::Add")]
pub fn curve_between(curve: &Curve, t1: f64, t2: f64, tolerance: f64) -> BoundingBox {
    to_box(curve_extent(curve, t1, t2), tolerance)
}

/// [`curve`], as tight as possible.
#[doc(alias = "BndLib_Add3dCurve::AddOptimal")]
#[doc(alias = "GeomBndLib_Curve::AddOptimal")]
pub fn curve_tight(curve: &Curve, tolerance: f64) -> BoundingBox {
    self::curve(curve, tolerance)
}

/// [`curve_between`], as tight as possible.
pub fn curve_between_tight(curve: &Curve, t1: f64, t2: f64, tolerance: f64) -> BoundingBox {
    curve_between(curve, t1, t2, tolerance)
}

/// The box of a 2D curve over its whole domain.
#[doc(alias = "BndLib_Add2dCurve::Add")]
#[doc(alias = "GeomBndLib_Curve2d::Add")]
pub fn curve2(curve: &Curve2, tolerance: f64) -> BoundingBox2 {
    curve2_between(curve, curve.first_parameter(), curve.last_parameter(), tolerance)
}

pub fn curve2_between(curve: &Curve2, t1: f64, t2: f64, tolerance: f64) -> BoundingBox2 {
    to_box2(curve_extent(curve, t1, t2), tolerance)
}

#[doc(alias = "BndLib_Add2dCurve::AddOptimal")]
pub fn curve2_tight(curve: &Curve2, tolerance: f64) -> BoundingBox2 {
    curve2(curve, tolerance)
}

pub fn curve2_between_tight(curve: &Curve2, t1: f64, t2: f64, tolerance: f64) -> BoundingBox2 {
    curve2_between(curve, t1, t2, tolerance)
}

/// The box of a surface over its whole domain.
#[doc(alias = "BndLib_AddSurface::Add")]
#[doc(alias = "GeomBndLib_Surface::Add")]
pub fn surface(surface: &Surface, tolerance: f64) -> BoundingBox {
    let (u1, u2, v1, v2) = surface.bounds();
    surface_within(surface, u1, u2, v1, v2, tolerance)
}

/// The box of the part of a surface over `[u1, u2] x [v1, v2]`.
#[doc(alias = "BndLib_AddSurface::Add")]
pub fn surface_within(surface: &Surface, u1: f64, u2: f64, v1: f64, v2: f64, tolerance: f64) -> BoundingBox {
    to_box(surface_extent(surface, (u1.min(u2), u1.max(u2)), (v1.min(v2), v1.max(v2))), tolerance)
}

/// [`surface`], as tight as possible.
#[doc(alias = "BndLib_AddSurface::AddOptimal")]
#[doc(alias = "GeomBndLib_Surface::AddOptimal")]
pub fn surface_tight(surface: &Surface, tolerance: f64) -> BoundingBox {
    self::surface(surface, tolerance)
}

/// [`surface_within`], as tight as possible.
pub fn surface_within_tight(surface: &Surface, u1: f64, u2: f64, v1: f64, v2: f64, tolerance: f64) -> BoundingBox {
    surface_within(surface, u1, u2, v1, v2, tolerance)
}

/// The extent of a surface over a parameter rectangle: its boundary curves, the interior
/// critical points of each coordinate, and the open sides of infinite ranges.
fn surface_extent(surface: &Surface, (u1, u2): (f64, f64), (v1, v2): (f64, f64)) -> Extent {
    let window = |a: f64, b: f64| match (is_infinite(a), is_infinite(b)) {
        (false, false) => (a, b),
        (true, false) => (b - 1e3, b),
        (false, true) => (a, a + 1e3),
        (true, true) => (-1e3, 1e3),
    };
    let ((a, b), (c, d)) = (window(u1, u2), window(v1, v2));
    let mut extent = NONE;
    // The boundary: along the four sides.
    let vec = |v: crate::math::Vector| [v.x, v.y, v.z];
    for v in [c, d] {
        for u in critical_1d(&|u| vec(surface.d1(u, v).1), a, b, 128, &[]) {
            add(&mut extent, surface.point(u, v).coords());
        }
    }
    for u in [a, b] {
        for v in critical_1d(&|v| vec(surface.d1(u, v).2), c, d, 128, &[]) {
            add(&mut extent, surface.point(u, v).coords());
        }
    }
    // Inside: where a coordinate's gradient vanishes, by Newton steps from a grid.
    const N: usize = 16;
    for i in 1..N {
        for j in 1..N {
            let (u0, v0) = (a + (b - a) * i as f64 / N as f64, c + (d - c) * j as f64 / N as f64);
            add(&mut extent, surface.point(u0, v0).coords());
            for k in 0..3 {
                if let Some((u, v)) = critical_point(surface, k, (u0, v0), (a, b, c, d)) {
                    add(&mut extent, surface.point(u, v).coords());
                }
            }
        }
    }
    // Open sides: along an infinite parameter, the coordinates the surface runs out on.
    for (open, along_u, sign) in [(is_infinite(u1), true, -1.0), (is_infinite(u2), true, 1.0), (is_infinite(v1), false, -1.0), (is_infinite(v2), false, 1.0)] {
        if !open {
            continue;
        }
        for i in 0..=8 {
            let s = i as f64 / 8.0;
            let (u, v) = if along_u { (sign * 1e3, c + (d - c) * s) } else { (a + (b - a) * s, sign * 1e3) };
            let (_, du, dv) = surface.d1(u, v);
            let d = if along_u { du } else { dv };
            let d = [d.x, d.y, d.z];
            for k in 0..3 {
                let direction = d[k] * sign;
                if direction > 1e-12 {
                    extent[k][1] = f64::INFINITY;
                } else if direction < -1e-12 {
                    extent[k][0] = f64::NEG_INFINITY;
                }
            }
        }
    }
    extent
}

/// A point inside the rectangle where coordinate `k` of the surface has zero gradient,
/// by Newton steps from `start`.
fn critical_point(surface: &Surface, k: usize, start: (f64, f64), (a, b, c, d): (f64, f64, f64, f64)) -> Option<(f64, f64)> {
    let coord = |v: crate::math::Vector| [v.x, v.y, v.z][k];
    let (mut u, mut v) = start;
    for _ in 0..30 {
        let s = surface.d2(u, v);
        let (gu, gv) = (coord(s.du), coord(s.dv));
        let (huu, huv, hvv) = (coord(s.duu), coord(s.duv), coord(s.dvv));
        let det = huu * hvv - huv * huv;
        if det.abs() <= 1e-300 {
            return None;
        }
        let (du, dv) = ((-gu * hvv + gv * huv) / det, (-gv * huu + gu * huv) / det);
        (u, v) = (u + du, v + dv);
        if !(a <= u && u <= b && c <= v && v <= d) {
            return None;
        }
        if du.abs() <= 1e-14 * (1.0 + u.abs()) && dv.abs() <= 1e-14 * (1.0 + v.abs()) {
            return Some((u, v));
        }
    }
    let s = surface.d1(u, v);
    (coord(s.1).abs() + coord(s.2).abs() <= 1e-9).then_some((u, v))
}
