//! Parametric resolutions: a parameter step that moves a point of a curve or surface by about
//! a given distance, from bounds of the first derivatives.
//!
//! These are the estimates that the surface analysis and line–face intersection size their
//! steps and solver tolerances with: for B-splines and Béziers the degree times the largest
//! pole difference (1-norm) over its knot span (a bound of the derivative), for rational ones
//! a bound through the weighted pole differences in a window of poles around each span; for
//! the elementary curves and surfaces closed forms (an arc's angle for a chord). They differ
//! from the curves' own `ParametricCurve::resolution`, which estimates differently.

use crate::geometry::bezier::{BezierCurve, BezierSurface};
use crate::geometry::bspline::{BSplineCurve, BSplineSurface};
use crate::geometry::curve::Curve;
use crate::geometry::surface::Surface;
use crate::geometry::ControlPoint;
use std::f64::consts::TAU;

const CONFUSION: f64 = crate::tolerance::LINEAR;

/// The parametric counterpart of a distance, for curves and surfaces without a closed form.
fn parametric(r3d: f64) -> f64 {
    r3d / 100.0
}

/// The gap from `x ≥ 0` to the next double above it.
fn ulp(x: f64) -> f64 {
    x.next_up() - x
}

fn norm1(c: [f64; 3]) -> f64 {
    c[0].abs() + c[1].abs() + c[2].abs()
}

fn minus(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// `a·s − b·t`, by coordinates.
fn combination(a: [f64; 3], s: f64, b: [f64; 3], t: f64) -> [f64; 3] {
    [a[0] * s - b[0] * t, a[1] * s - b[1] * t, a[2] * s - b[2] * t]
}

/// `1/d`, or `1/DBL_MIN` for a vanishing `d`.
fn inverse(d: f64) -> f64 {
    if d > f64::MIN_POSITIVE { 1.0 / d } else { 1.0 / f64::MIN_POSITIVE }
}

/// The knots repeated by multiplicity; for a periodic vector (whose last knot is the first one
/// a period on) the period's knots extended cyclically so that the sequence has
/// `Σ mults + 2(degree + 1 − mults[0])` entries, the last copy of the first knot at index
/// `degree`.
fn flat_knots(knots: &[f64], mults: &[usize], degree: usize, periodic: bool) -> Vec<f64> {
    if !periodic {
        return knots.iter().zip(mults).flat_map(|(&k, &m)| std::iter::repeat_n(k, m)).collect();
    }
    let nk = knots.len();
    let period = knots[nk - 1] - knots[0];
    let one: Vec<f64> = knots[..nk - 1].iter().zip(&mults[..nk - 1]).flat_map(|(&k, &m)| std::iter::repeat_n(k, m)).collect();
    let q = one.len() as i64;
    let shift = degree as i64 + 1 - mults[0] as i64;
    let total: usize = mults.iter().sum::<usize>() + 2 * shift.max(0) as usize;
    (0..total as i64)
        .map(|i| {
            let j = i - shift;
            let (r, m) = (j.div_euclid(q), j.rem_euclid(q));
            one[m as usize] + r as f64 * period
        })
        .collect()
}

/// The pole count of a periodic B-spline made non-periodic (its poles then repeat cyclically).
fn unperiodic_pole_count(mults: &[usize], degree: usize) -> usize {
    let nk = mults.len();
    let order = degree + 1;
    // From the far end inwards, until the multiplicities reach the order.
    let mut s = mults[0];
    let mut j = nk as i64 - 2;
    while s < order && j >= 0 {
        s += mults[j as usize];
        j -= 1;
    }
    let a = s.min(order) - mults[0];
    let mut s = mults[nk - 1];
    let mut j = 1;
    while s < order && j <= nk - 1 {
        s += mults[j];
        j += 1;
    }
    let b = s.min(order) - mults[nk - 1];
    (mults.iter().sum::<usize>() + a + b).saturating_sub(order)
}

/// The resolution of a curve given by poles `poles[j mod c]`, weights (when rational), flat
/// knots and degree.
fn curve_resolution_from(poles: &[[f64; 3]], weights: Option<&[f64]>, flat: &[f64], degree: usize, tolerance_3d: f64) -> f64 {
    let c = poles.len();
    let p = degree;
    let n_flat = flat.len().saturating_sub(p + 1);
    let rational = weights.is_some_and(|w| w.windows(2).any(|x| (x[0] - x[1]).abs() > f64::MIN_POSITIVE));
    let mut d: f64 = 0.0;
    if c == 0 {
        return tolerance_3d * inverse(d);
    }
    if rational {
        let w = weights.expect("rational");
        let wj = |j: usize| w[j % w.len()];
        let w_min = (0..c).map(wj).fold(f64::INFINITY, f64::min);
        let mut most: f64 = 0.0;
        for i in 1..n_flat {
            let h = flat[i + p] - flat[i];
            let (qi, qh) = (poles[i % c], poles[(i - 1) % c]);
            let (wi, wh) = (wj(i % c), wj((i - 1) % c));
            let lo = i.saturating_sub(p + 1);
            let hi = (i + 2 * p + 1).min(n_flat);
            for j in lo..hi {
                let qj = poles[j % c];
                let v = combination(minus(qj, qi), wi, minus(qj, qh), wh);
                most = most.max(norm1(v) / h);
            }
        }
        d = p as f64 / w_min * most;
    } else {
        for i in 1..n_flat {
            let h = flat[i + p] - flat[i];
            d = d.max(norm1(minus(poles[i % c], poles[(i - 1) % c])) / h);
        }
        d *= p as f64;
    }
    tolerance_3d * inverse(d)
}

/// The parameter step of a B-spline curve for a 3D distance.
pub(crate) fn bspline_curve_resolution<P: ControlPoint>(curve: &BSplineCurve<P>, tolerance_3d: f64) -> f64 {
    let p = curve.degree();
    let periodic = crate::geometry::ParametricCurve::is_periodic(curve);
    let flat = flat_knots(curve.knots(), curve.multiplicities(), p, periodic);
    let stored: Vec<[f64; 3]> = curve.poles().iter().map(|q| q.coords()).collect();
    let n = stored.len();
    let (poles, weights): (Vec<[f64; 3]>, Option<Vec<f64>>) = if periodic {
        let c = unperiodic_pole_count(curve.multiplicities(), p);
        ((0..c).map(|j| stored[j % n]).collect(), curve.weights().map(|w| (0..c).map(|j| w[j % n]).collect()))
    } else {
        (stored, curve.weights().map(|w| w.to_vec()))
    };
    let rational = curve.weights().is_some_and(|w| w.windows(2).any(|x| (x[0] - x[1]).abs() > f64::MIN_POSITIVE));
    curve_resolution_from(&poles, if rational { weights.as_deref() } else { None }, &flat, p, tolerance_3d)
}

/// The parameter step of a Bézier curve for a 3D distance.
pub(crate) fn bezier_curve_resolution<P: ControlPoint>(curve: &BezierCurve<P>, tolerance_3d: f64) -> f64 {
    let p = curve.poles().len() - 1;
    let flat: Vec<f64> = std::iter::repeat_n(0.0, p + 1).chain(std::iter::repeat_n(1.0, p + 1)).collect();
    let poles: Vec<[f64; 3]> = curve.poles().iter().map(|q| q.coords()).collect();
    curve_resolution_from(&poles, curve.weights(), &flat, p, tolerance_3d)
}

/// The data of a pole net for the surface resolutions.
struct Net<'a> {
    poles: &'a [Vec<crate::math::Point>],
    weights: Option<&'a [Vec<f64>]>,
    u_flat: Vec<f64>,
    v_flat: Vec<f64>,
    u_degree: usize,
    v_degree: usize,
}

fn net_resolution(net: &Net, tolerance_3d: f64) -> (f64, f64) {
    let nu = net.poles.len();
    let nv = net.poles.first().map_or(0, |r| r.len());
    if nu == 0 || nv == 0 {
        return (0.0, 0.0);
    }
    let (pu, pv) = (net.u_degree, net.v_degree);
    let n_u = net.u_flat.len().saturating_sub(pu + 1);
    let n_v = net.v_flat.len().saturating_sub(pv + 1);
    let pole = |i: usize, j: usize| {
        let q = net.poles[i % nu][j % nv];
        [q.x, q.y, q.z]
    };
    // (The flags cross: weights varying along the first index make the v direction rational.)
    let (mut u_rational, mut v_rational) = (false, false);
    if let Some(w) = net.weights {
        for i in 0..nu {
            for j in 0..nv {
                let here = w[i][j];
                if i + 1 < nu && (here - w[i + 1][j]).abs() > ulp(here.abs()) {
                    v_rational = true;
                }
                if j + 1 < nv && (here - w[i][j + 1]).abs() > ulp(here.abs()) {
                    u_rational = true;
                }
            }
        }
    }
    let weight = |i: usize, j: usize| net.weights.map_or(1.0, |w| w[i % nu][j % nv]);
    let w_min = if u_rational || v_rational { net.weights.expect("weights").iter().flatten().copied().fold(f64::INFINITY, f64::min) } else { 1.0 };
    let window = |i: usize, p: usize, n: usize| (i.saturating_sub(p + 1), (i + 2 * p + 1).min(n.saturating_sub(1)));

    let mut d_u: f64 = 0.0;
    if u_rational {
        for i in 1..n_u {
            let h = net.u_flat[i + pu] - net.u_flat[i];
            for j in 0..n_v {
                let (pi, ph) = (pole(i, j), pole(i - 1, j));
                let (wi, wh) = (weight(i, j), weight(i - 1, j));
                let (a0, a1) = window(i, pu, n_u);
                let (b0, b1) = window(j, pv, n_v);
                for a in a0..=a1 {
                    for b in b0..=b1 {
                        let q = pole(a, b);
                        d_u = d_u.max(norm1(combination(minus(q, pi), wi, minus(q, ph), wh)) / h);
                    }
                }
            }
        }
        d_u *= pu as f64 / w_min;
    } else {
        for i in 1..n_u {
            let h = net.u_flat[i + pu] - net.u_flat[i];
            for j in 0..n_v {
                d_u = d_u.max(norm1(minus(pole(i, j), pole(i - 1, j))) / h);
            }
        }
        d_u *= pu as f64;
    }
    let mut d_v: f64 = 0.0;
    if v_rational {
        for i in 1..n_v {
            let h = net.v_flat[i + pv] - net.v_flat[i];
            for j in 0..n_u {
                let (pi, ph) = (pole(j, i), pole(j, i - 1));
                let (wi, wh) = (weight(j, i), weight(j, i - 1));
                let (a0, a1) = window(j, pu, n_u);
                let (b0, b1) = window(i, pv, n_v);
                for a in a0..=a1 {
                    for b in b0..=b1 {
                        // (The v-window index as the row, the u-window index as the column.)
                        let q = pole(b % nu, a % nv);
                        d_v = d_v.max(norm1(combination(minus(q, pi), wi, minus(q, ph), wh)) / h);
                    }
                }
            }
        }
        d_v *= pv as f64 / w_min;
    } else {
        for i in 1..n_v {
            let h = net.v_flat[i + pv] - net.v_flat[i];
            for j in 0..n_u {
                d_v = d_v.max(norm1(minus(pole(j, i), pole(j, i - 1))) / h);
            }
        }
        d_v *= pv as f64;
    }
    d_u *= std::f64::consts::SQRT_2;
    d_v *= std::f64::consts::SQRT_2;
    if d_u != 0.0 && d_v != 0.0 {
        (tolerance_3d * (1.0 / d_u), tolerance_3d * (1.0 / d_v))
    } else {
        (0.0, 0.0)
    }
}

/// The parameter steps `(u, v)` of a B-spline surface for a 3D distance.
pub(crate) fn bspline_surface_resolution(surface: &BSplineSurface, tolerance_3d: f64) -> (f64, f64) {
    let (u, v) = (surface.u_knot_vector(), surface.v_knot_vector());
    let net = Net {
        poles: surface.poles(),
        weights: surface.weights(),
        u_flat: flat_knots(&u.knots, &u.multiplicities, u.degree, u.periodic),
        v_flat: flat_knots(&v.knots, &v.multiplicities, v.degree, v.periodic),
        u_degree: u.degree,
        v_degree: v.degree,
    };
    net_resolution(&net, tolerance_3d)
}

/// The parameter steps `(u, v)` of a Bézier surface for a 3D distance.
pub(crate) fn bezier_surface_resolution(surface: &BezierSurface, tolerance_3d: f64) -> (f64, f64) {
    let bezier_flat = |p: usize| -> Vec<f64> { std::iter::repeat_n(0.0, p + 1).chain(std::iter::repeat_n(1.0, p + 1)).collect() };
    let pu = surface.poles().len() - 1;
    let pv = surface.poles()[0].len() - 1;
    let net = Net { poles: surface.poles(), weights: surface.weights(), u_flat: bezier_flat(pu), v_flat: bezier_flat(pv), u_degree: pu, v_degree: pv };
    net_resolution(&net, tolerance_3d)
}

/// The parameter step of a curve for a 3D distance `r3d`.
pub(crate) fn curve_resolution(curve: &Curve, r3d: f64) -> f64 {
    match curve {
        Curve::Line(_) => r3d,
        Curve::Circle(c) => {
            let r = c.radius();
            if r > r3d / 2.0 { 2.0 * (r3d / (2.0 * r)).asin() } else { TAU }
        }
        Curve::Ellipse(e) => r3d / e.major_radius(),
        Curve::Hyperbola(_) | Curve::Parabola(_) | Curve::Offset(_) => parametric(r3d),
        Curve::Trimmed(t) => curve_resolution(t.basis(), r3d),
        Curve::Bezier(b) => bezier_curve_resolution(b, r3d),
        Curve::BSpline(b) => bspline_curve_resolution(b, r3d),
    }
}

/// `2 asin ρ`, or a whole turn beyond `ρ = 1`.
fn angle(rho: f64) -> f64 {
    if rho <= 1.0 { 2.0 * rho.asin() } else { TAU }
}

/// The angular step of a circle of radius `r` for a chord `r3d` (0 for a vanishing radius).
fn turn(r3d: f64, r: f64) -> f64 {
    if r > CONFUSION { angle(r3d / (2.0 * r)) } else { 0.0 }
}

/// The parameter step along u of `surface` for a 3D distance `r3d`; `bounds` `(u1, u2, v1, v2)`
/// are the range of the surface in use (a cone's v range sets its widest circle).
pub(crate) fn u_resolution(surface: &Surface, bounds: (f64, f64, f64, f64), r3d: f64) -> f64 {
    match surface {
        Surface::Plane(_) => r3d,
        Surface::Cylinder(c) => turn(r3d, c.radius()),
        Surface::Sphere(s) => turn(r3d, s.radius()),
        Surface::Torus(t) => turn(r3d, t.major_radius() + t.minor_radius()),
        Surface::Cone(c) => {
            let (_, _, v1, v2) = bounds;
            if v2 - v1 > 1e10 {
                return parametric(r3d);
            }
            let sin = c.semi_angle_sin_cos().0;
            let r0 = c.radius();
            let rho = (r0 + v2 * sin).abs().max((r0 + v1 * sin).abs());
            if rho > CONFUSION { r3d / rho } else { 0.0 }
        }
        Surface::Revolution(_) => parametric(r3d),
        Surface::Extrusion(e) => curve_resolution(e.basis(), r3d),
        Surface::Offset(o) => u_resolution(o.basis(), bounds, r3d),
        Surface::Trimmed(t) => u_resolution(t.basis(), bounds, r3d),
        Surface::Bezier(b) => bezier_surface_resolution(b, r3d).0,
        Surface::BSpline(b) => bspline_surface_resolution(b, r3d).0,
    }
}

/// The parameter step along v of `surface` for a 3D distance `r3d` (see [`u_resolution`]).
pub(crate) fn v_resolution(surface: &Surface, bounds: (f64, f64, f64, f64), r3d: f64) -> f64 {
    match surface {
        Surface::Plane(_) | Surface::Cylinder(_) | Surface::Cone(_) | Surface::Extrusion(_) => r3d,
        Surface::Sphere(s) => turn(r3d, s.radius()),
        Surface::Torus(t) => turn(r3d, t.minor_radius()),
        Surface::Revolution(r) => curve_resolution(r.basis(), r3d),
        Surface::Offset(o) => v_resolution(o.basis(), bounds, r3d),
        Surface::Trimmed(t) => v_resolution(t.basis(), bounds, r3d),
        Surface::Bezier(b) => bezier_surface_resolution(b, r3d).1,
        Surface::BSpline(b) => bspline_surface_resolution(b, r3d).1,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::bspline::KnotVector;
    use crate::math::{Direction, Frame, Point};

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() <= 1e-13 * a.abs().max(b.abs())
    }

    #[test]
    fn circles_and_lines() {
        let frame = Frame::new(Point::ORIGIN, Direction::new(0.0, 0.0, 1.0), Direction::new(1.0, 0.0, 0.0));
        let small = Curve::Circle(crate::geometry::curve::Circle::new(frame, 0.01));
        assert_eq!(curve_resolution(&small, 0.7), TAU);
        let big = Curve::Circle(crate::geometry::curve::Circle::new(frame, 10.0));
        assert!(close(curve_resolution(&big, 0.01), 2.0 * (0.01f64 / 20.0).asin()));
        let line = Curve::Line(crate::geometry::curve::Line::new(Point::ORIGIN, Direction::new(1.0, 0.0, 0.0)));
        assert_eq!(curve_resolution(&line, 0.3), 0.3);
    }

    #[test]
    fn straight_degree_one_bspline() {
        // Uniform knots of span h, poles a step d apart: r3d·h / (p·‖ΔP‖₁).
        let poles: Vec<Point> = (0..5).map(|i| Point::new(2.0 * i as f64, i as f64, 0.0)).collect();
        let knots: Vec<f64> = (0..5).map(|i| 0.5 * i as f64).collect();
        let b = BSplineCurve::new(poles, knots, vec![2, 1, 1, 1, 2], 1);
        assert!(close(bspline_curve_resolution(&b, 0.01), 0.01 * 0.5 / 3.0));
    }

    #[test]
    fn equal_weights_are_not_rational() {
        let poles = vec![Point::new(0.0, 0.0, 0.0), Point::new(1.0, 2.0, 0.0), Point::new(3.0, 2.0, 1.0), Point::new(4.0, 0.0, 0.0)];
        let plain = BezierCurve::new(poles.clone());
        let weighted = BezierCurve::rational(poles, vec![2.0; 4]);
        assert_eq!(bezier_curve_resolution(&plain, 0.1), bezier_curve_resolution(&weighted, 0.1));
        // Degree 3 over [0, 1]: largest pole step 1-norm 4 (from (3,2,1) to (4,0,0)).
        assert!(close(bezier_curve_resolution(&plain, 0.1), 0.1 / (3.0 * 4.0)));
    }

    #[test]
    fn cone_radius_changing_sign() {
        let frame = Frame::new(Point::ORIGIN, Direction::new(0.0, 0.0, 1.0), Direction::new(1.0, 0.0, 0.0));
        let cone = crate::geometry::surface::Cone::new(frame, 0.5, 1.0);
        let sin = 0.5f64.sin();
        let s = Surface::Cone(cone);
        let r = u_resolution(&s, (0.0, TAU, -5.0, 2.0), 0.01);
        let rho = (1.0 + 2.0 * sin).abs().max((1.0 - 5.0 * sin).abs());
        assert!(close(r, 0.01 / rho));
        assert_eq!(u_resolution(&s, (0.0, TAU, -2e100, 2e100), 0.01), 0.0001);
        assert_eq!(v_resolution(&s, (0.0, TAU, -5.0, 2.0), 0.01), 0.01);
    }

    #[test]
    fn bilinear_patch() {
        // A flat 2×2 patch of unit side: Du = Dv = √2·1·1.
        let poles = vec![vec![Point::new(0.0, 0.0, 0.0), Point::new(0.0, 1.0, 0.0)], vec![Point::new(1.0, 0.0, 0.0), Point::new(1.0, 1.0, 0.0)]];
        let b = BSplineSurface::new(poles.clone(), KnotVector::new(vec![0.0, 1.0], vec![2, 2], 1), KnotVector::new(vec![0.0, 1.0], vec![2, 2], 1));
        let (u, v) = bspline_surface_resolution(&b, 0.1);
        assert!(close(u, 0.1 / 2f64.sqrt()) && close(v, 0.1 / 2f64.sqrt()), "{u} {v}");
        let z = BezierSurface::new(poles);
        assert_eq!(bezier_surface_resolution(&z, 0.1), (u, v));
        // Degenerate (all poles equal): nothing.
        let same = vec![vec![Point::ORIGIN; 2]; 2];
        assert_eq!(bezier_surface_resolution(&BezierSurface::new(same), 0.1), (0.0, 0.0));
    }

    #[test]
    fn periodic_flat_knots_and_poles() {
        let flat = flat_knots(&[0.0, 1.0, 2.0, 3.0, 4.0], &[1, 1, 1, 1, 1], 3, true);
        assert_eq!(flat.len(), 5 + 6);
        assert_eq!(flat[3], 0.0);
        assert_eq!(flat[0], -3.0);
        assert_eq!(unperiodic_pole_count(&[1, 1, 1, 1, 1], 3), 7);
    }
}
