//! OCCT's parametric resolutions, ported: the parameter steps that move a point of a curve
//! or surface by at most a given 3D distance, as `GeomAdaptor_Curve::Resolution` and
//! `GeomAdaptor_Surface::UResolution` / `VResolution` give them, with the B-spline and
//! Bézier estimates under them (`BSplCLib::Resolution`, `BSplSLib::Resolution`, through
//! `Geom_BSplineCurve`, `Geom_BezierCurve`, `Geom_BSplineSurface` and `Geom_BezierSurface`).
//!
//! These are not kentron's own [`ParametricCurve::resolution`](super::ParametricCurve) and
//! [`ParametricSurface::resolution`](super::ParametricSurface), which estimate differently
//! (a Euclidean pole step, an ellipse taken as a circle of its major radius, ...): code
//! that ports OCCT and must take OCCT's steps uses these.
//!
//! As OCCT's do, the B-spline and Bézier ones take the weights as rational only when they
//! differ (`Geom_BSplineCurve`'s and `Geom_BezierCurve`'s `Rational`, by more than
//! `gp::Resolution()`; `Geom_BSplineSurface`'s and `Geom_BezierSurface`'s, by more than
//! `Epsilon`), where kentron keeps any weights it is given.

use std::f64::consts::{PI, SQRT_2};

use super::ControlPoint;
use super::adaptor::epsilon;
use super::bezier::{BezierCurve, BezierSurface};
use super::bspline::{BSplineCurve, BSplineSurface, periodic_flatten};
use super::curve::Curve;
use super::surface::{Cone, Surface};
use crate::math::Point;

/// `Precision::Confusion()`.
const CONFUSION: f64 = 1.0e-7;
/// `RealSmall()` (`DBL_MIN`), which is also `gp::Resolution()`.
const REAL_SMALL: f64 = f64::MIN_POSITIVE;

/// `Precision::Parametric(P)`.
fn parametric(p: f64) -> f64 {
    p * 0.01
}

/// `BSplCLib::KnotSequence(Knots, Mults, Degree, Periodic)`: the flat knots, of
/// `BSplCLib::KnotSequenceLength` (for a periodic sequence, `Degree + 1 - Mults(first)`
/// knots of the period before and as many of the period after added).
#[doc(alias = "BSplCLib::KnotSequence")]
fn knot_sequence(knots: &[f64], multiplicities: &[usize], degree: usize, periodic: bool) -> Vec<f64> {
    if !periodic {
        return knots.iter().zip(multiplicities).flat_map(|(&k, &m)| std::iter::repeat(k).take(m)).collect();
    }
    // (kentron's periodic flat knots are OCCT's, `pkn[j] - period` and `pkn[j] + period`,
    // from the same position, and run on for `Mults(first) - 1` more knots; OCCT's also
    // read past the knots where too few follow the first, which these extend periodically.)
    let length = multiplicities.iter().sum::<usize>() + 2 * (degree + 1 - multiplicities[0]);
    let mut flat = periodic_flatten(knots, multiplicities, degree);
    flat.truncate(length);
    flat
}

/// `BSplCLib::Resolution(Poles, ArrayDimension, NumPoles, Weights, FlatKnots, Degree,
/// Tolerance3D, UTolerance)`: `tolerance_3d` over a bound of the first derivative,
/// `Degree` times the largest (L1) pole step over its knot span, for a rational curve the
/// weighted steps to every pole near over the smallest weight. The poles (`NumPoles` of
/// them) are taken cyclically along the flat knots; a 2D curve's poles have a zero third
/// coordinate, which adds nothing to the sums (OCCT's dimension-2 branch).
#[doc(alias = "BSplCLib::Resolution")]
fn bsplclib_resolution(poles: &[[f64; 3]], weights: Option<&[f64]>, flat_knots: &[f64], degree: usize, tolerance_3d: f64) -> f64 {
    let num_poles_arg = poles.len();
    let deg1 = degree + 1;
    let deg2 = (degree << 1) + 1;
    let num_poles = flat_knots.len() - deg1;
    let mut max_derivative = 0.0e0;
    match weights {
        Some(wg) => {
            let mut min_weights = wg[0];
            for &w in &wg[1..num_poles_arg] {
                if w < min_weights {
                    min_weights = w;
                }
            }
            for ii in 1..num_poles {
                let ii_index = ii % num_poles_arg;
                let ii_minus = (ii - 1) % num_poles_arg;
                let pa_ii = poles[ii_index];
                let pa_mi = poles[ii_minus];
                let wg_ii_index = wg[ii_index];
                let wg_ii_minus = wg[ii_minus];
                let mut inverse = flat_knots[ii + degree] - flat_knots[ii];
                inverse = 1.0e0 / inverse;
                let lower = ii.saturating_sub(deg1);
                let upper = (deg2 + ii).min(num_poles);
                for jj in lower..upper {
                    let pa_jj = poles[jj % num_poles_arg];
                    let mut value = 0.0e0;
                    for k in 0..3 {
                        let factor = ((pa_jj[k] - pa_ii[k]) * wg_ii_index) - ((pa_jj[k] - pa_mi[k]) * wg_ii_minus);
                        value += factor.abs();
                    }
                    value *= inverse;
                    if max_derivative < value {
                        max_derivative = value;
                    }
                }
            }
            max_derivative /= min_weights;
        }
        None => {
            for ii in 1..num_poles {
                let pa_ii = poles[ii % num_poles_arg];
                let pa_mi = poles[(ii - 1) % num_poles_arg];
                let mut inverse = flat_knots[ii + degree] - flat_knots[ii];
                inverse = 1.0e0 / inverse;
                let mut value = 0.0e0;
                for k in 0..3 {
                    value += (pa_ii[k] - pa_mi[k]).abs();
                }
                value *= inverse;
                if max_derivative < value {
                    max_derivative = value;
                }
            }
        }
    }
    max_derivative *= degree as f64;
    if max_derivative > REAL_SMALL { tolerance_3d / max_derivative } else { tolerance_3d / REAL_SMALL }
}

/// `Rational(Weights)` of `Geom_BSplineCurve` and `Geom_BezierCurve`: some two weights in a
/// row differ by more than `gp::Resolution()`.
fn curve_weights_rational(weights: Option<&[f64]>) -> Option<&[f64]> {
    weights.filter(|w| w.windows(2).any(|p| (p[0] - p[1]).abs() > REAL_SMALL))
}

/// `Geom_BSplineCurve::Resolution(Tolerance3D, UTolerance)`: of a periodic curve on its
/// poles unperiodized (`BSplCLib::PrepareUnperiodize`'s count, taken cyclically).
#[doc(alias = "Geom_BSplineCurve::Resolution")]
pub(crate) fn bspline_curve_resolution<P: ControlPoint>(curve: &BSplineCurve<P>, tolerance_3d: f64) -> f64 {
    use super::ParametricCurve;
    let degree = curve.degree();
    let weights = curve_weights_rational(curve.weights());
    let periodic = curve.is_periodic();
    let flat = knot_sequence(curve.knots(), curve.multiplicities(), degree, periodic);
    let coords: Vec<[f64; 3]> = curve.poles().iter().map(|p| p.coords()).collect();
    let max_deriv_inv = if periodic {
        let count = prepare_unperiodize_pole_count(degree, curve.multiplicities());
        let n = coords.len();
        let new_poles: Vec<[f64; 3]> = (0..count).map(|ii| coords[ii % n]).collect();
        let new_weights: Option<Vec<f64>> = weights.map(|w| (0..count).map(|ii| w[ii % n]).collect());
        bsplclib_resolution(&new_poles, new_weights.as_deref(), &flat, degree, 1.0)
    } else {
        bsplclib_resolution(&coords, weights, &flat, degree, 1.0)
    };
    tolerance_3d * max_deriv_inv
}

/// `BSplCLib::PrepareUnperiodize(Degree, Mults, NbKnots, NbPoles)`'s `NbPoles`: the poles
/// of the periodic curve made non-periodic, with knots added at both ends until the end
/// multiplicities reach `Degree + 1`.
#[doc(alias = "BSplCLib::PrepareUnperiodize")]
fn prepare_unperiodize_pole_count(degree: usize, multiplicities: &[usize]) -> usize {
    let mults: Vec<i64> = multiplicities.iter().map(|&m| m as i64).collect();
    let (lower, upper) = (0usize, mults.len() - 1);
    let degree = degree as i64;
    let mut nb_poles = -degree - 1;
    for &m in &mults {
        nb_poles += m;
    }
    // Add knots at the beginning of the curve to raise multiplicities to Degree + 1.
    // (Where too few knots precede, OCCT reads before its array: the count stops there.)
    let mut sigma = mults[lower];
    let mut k = upper as i64 - 1;
    while sigma < degree + 1 && k >= lower as i64 {
        sigma += mults[k as usize];
        nb_poles += mults[k as usize];
        k -= 1;
    }
    // We must add exactly until Degree + 1: suppress the excess.
    if sigma > degree + 1 {
        nb_poles -= sigma - degree - 1;
    }
    // Add knots at the end of the curve to raise multiplicities to Degree + 1.
    sigma = mults[upper];
    k = lower as i64 + 1;
    while sigma < degree + 1 && k <= upper as i64 {
        sigma += mults[k as usize];
        nb_poles += mults[k as usize];
        k += 1;
    }
    if sigma > degree + 1 {
        nb_poles -= sigma - degree - 1;
    }
    nb_poles as usize
}

/// The flat knots of a Bézier of degree `degree`: `[0; degree + 1]` then `[1; degree + 1]`
/// (`Geom_BezierCurve::KnotSequence`).
fn bezier_flat_knots(degree: usize) -> Vec<f64> {
    let mut flat = vec![0.0; degree + 1];
    flat.extend(std::iter::repeat(1.0).take(degree + 1));
    flat
}

/// `Geom_BezierCurve::Resolution(Tolerance3D, UTolerance)`.
#[doc(alias = "Geom_BezierCurve::Resolution")]
pub(crate) fn bezier_curve_resolution<P: ControlPoint>(curve: &BezierCurve<P>, tolerance_3d: f64) -> f64 {
    let degree = curve.poles().len() - 1;
    let coords: Vec<[f64; 3]> = curve.poles().iter().map(|p| p.coords()).collect();
    let weights = curve_weights_rational(curve.weights());
    let max_deriv_inv = bsplclib_resolution(&coords, weights, &bezier_flat_knots(degree), degree, 1.0);
    tolerance_3d * max_deriv_inv
}

/// `Rational(Weights, Urational, Vrational)` of `Geom_BSplineSurface` (and the same test
/// in `Geom_BezierSurface`): `(URational, VRational)`. As in OCCT, `VRational` is set by
/// weights differing along the first index (u, `poles[i]`), `URational` by weights
/// differing along the second.
fn surface_weights_rational(weights: &[Vec<f64>]) -> (bool, bool) {
    let (rows, cols) = (weights.len(), weights[0].len());
    let mut v_rational = false;
    let mut j = 0;
    while !v_rational && j < cols {
        let mut i = 0;
        while !v_rational && i + 1 < rows {
            v_rational = (weights[i][j] - weights[i + 1][j]).abs() > epsilon(weights[i][j].abs());
            i += 1;
        }
        j += 1;
    }
    let mut u_rational = false;
    let mut i = 0;
    while !u_rational && i < rows {
        let mut j = 0;
        while !u_rational && j + 1 < cols {
            u_rational = (weights[i][j] - weights[i][j + 1]).abs() > epsilon(weights[i][j].abs());
            j += 1;
        }
        i += 1;
    }
    (u_rational, v_rational)
}

/// `BSplSLib::Resolution(Poles, Weights, UKnots, VKnots, UMults, VMults, UDegree, VDegree,
/// URational, VRational, UPeriodic, VPeriodic, Tolerance3D, UTolerance, VTolerance)`, on
/// the flat knots: `(UTolerance, VTolerance)`, both 0 when either derivative bound is.
/// (OCCT's rational v branch takes the neighbouring poles' rows from the v range and their
/// columns from the u range, each index cycled by the other direction's pole count; so
/// does this.)
#[doc(alias = "BSplSLib::Resolution")]
#[allow(clippy::too_many_arguments)]
fn bsplslib_resolution(
    poles: &[Vec<Point>],
    weights: Option<&[Vec<f64>]>,
    flat_knots_in_u: &[f64],
    flat_knots_in_v: &[f64],
    u_degree: usize,
    v_degree: usize,
    u_rational: bool,
    v_rational: bool,
    tolerance_3d: f64,
) -> (f64, f64) {
    // (1-based indices as OCCT's, `at` and `w` taking them.)
    let at = |i: i64, j: i64| -> [f64; 3] { poles[(i - 1) as usize][(j - 1) as usize].coords() };
    let w = |i: i64, j: i64| -> f64 { weights.expect("weights of a rational surface")[(i - 1) as usize][(j - 1) as usize] };
    let flat_u = |i: i64| flat_knots_in_u[(i - 1) as usize];
    let flat_v = |i: i64| flat_knots_in_v[(i - 1) as usize];
    let mut max_derivative = [0.0e0f64; 2];
    let mut min_weights = 0.0;
    if u_rational || v_rational {
        let wg = weights.expect("weights of a rational surface");
        min_weights = wg[0][0];
        for &m in wg.iter().flatten().skip(1) {
            if min_weights > m {
                min_weights = m;
            }
        }
    }
    let (u_degree, v_degree) = (u_degree as i64, v_degree as i64);
    let ud1 = u_degree + 1;
    let vd1 = v_degree + 1;
    let num_poles = [flat_knots_in_u.len() as i64 - ud1, flat_knots_in_v.len() as i64 - vd1];
    let poles_length = [poles.len() as i64, poles[0].len() as i64];
    if u_rational {
        let ud2 = u_degree << 1;
        let vd2 = v_degree << 1;
        for ii in 2..=num_poles[0] {
            let ii_index = (ii - 1) % poles_length[0] + 1;
            let ii_minus = (ii - 2) % poles_length[0] + 1;
            let mut inverse = flat_u(ii + u_degree) - flat_u(ii);
            inverse = 1.0e0 / inverse;
            let lower0 = (ii - ud1).max(1);
            let upper0 = (ii + ud2 + 1).min(num_poles[0]);
            for jj in 1..=num_poles[1] {
                let jj_index = (jj - 1) % poles_length[1] + 1;
                let lower1 = (jj - vd1).max(1);
                let upper1 = (jj + vd2 + 1).min(num_poles[1]);
                let (pij, wij) = (at(ii_index, jj_index), w(ii_index, jj_index));
                let (pmj, wmj) = (at(ii_minus, jj_index), w(ii_minus, jj_index));
                for pp in lower0..=upper0 {
                    let pp_index = (pp - 1) % poles_length[0] + 1;
                    for qq in lower1..=upper1 {
                        let mut value = 0.0e0;
                        let qq_index = (qq - 1) % poles_length[1] + 1;
                        let ppq = at(pp_index, qq_index);
                        for k in 0..3 {
                            let mut factor = (ppq[k] - pij[k]) * wij;
                            factor -= (ppq[k] - pmj[k]) * wmj;
                            value += factor.abs();
                        }
                        value *= inverse;
                        if max_derivative[0] < value {
                            max_derivative[0] = value;
                        }
                    }
                }
            }
        }
        max_derivative[0] /= min_weights;
    } else {
        for ii in 2..=num_poles[0] {
            let ii_index = (ii - 1) % poles_length[0] + 1;
            let ii_minus = (ii - 2) % poles_length[0] + 1;
            let mut inverse = flat_u(ii + u_degree) - flat_u(ii);
            inverse = 1.0e0 / inverse;
            for jj in 1..=num_poles[1] {
                let jj_index = (jj - 1) % poles_length[1] + 1;
                let mut value = 0.0e0;
                let (pij, pmj) = (at(ii_index, jj_index), at(ii_minus, jj_index));
                for k in 0..3 {
                    value += (pij[k] - pmj[k]).abs();
                }
                value *= inverse;
                if max_derivative[0] < value {
                    max_derivative[0] = value;
                }
            }
        }
    }
    max_derivative[0] *= u_degree as f64;
    if v_rational {
        let ud2 = u_degree << 1;
        let vd2 = v_degree << 1;
        for ii in 2..=num_poles[1] {
            let ii_index = (ii - 1) % poles_length[1] + 1;
            let ii_minus = (ii - 2) % poles_length[1] + 1;
            let mut inverse = flat_v(ii + v_degree) - flat_v(ii);
            inverse = 1.0e0 / inverse;
            let lower0 = (ii - vd1).max(1);
            let upper0 = (ii + vd2 + 1).min(num_poles[1]);
            for jj in 1..=num_poles[0] {
                let jj_index = (jj - 1) % poles_length[0] + 1;
                let lower1 = (jj - ud1).max(1);
                let upper1 = (jj + ud2 + 1).min(num_poles[0]);
                let (pji, wji) = (at(jj_index, ii_index), w(jj_index, ii_index));
                let (pjm, wjm) = (at(jj_index, ii_minus), w(jj_index, ii_minus));
                for pp in lower1..=upper1 {
                    let pp_index = (pp - 1) % poles_length[1] + 1;
                    for qq in lower0..=upper0 {
                        let mut value = 0.0e0;
                        let qq_index = (qq - 1) % poles_length[0] + 1;
                        let pqp = at(qq_index, pp_index);
                        for k in 0..3 {
                            let mut factor = (pqp[k] - pji[k]) * wji;
                            factor -= (pqp[k] - pjm[k]) * wjm;
                            value += factor.abs();
                        }
                        value *= inverse;
                        if max_derivative[1] < value {
                            max_derivative[1] = value;
                        }
                    }
                }
            }
        }
        max_derivative[1] /= min_weights;
    } else {
        for ii in 2..=num_poles[1] {
            let ii_index = (ii - 1) % poles_length[1] + 1;
            let ii_minus = (ii - 2) % poles_length[1] + 1;
            let mut inverse = flat_v(ii + v_degree) - flat_v(ii);
            inverse = 1.0e0 / inverse;
            for jj in 1..=num_poles[0] {
                let jj_index = (jj - 1) % poles_length[0] + 1;
                let mut value = 0.0e0;
                let (pji, pjm) = (at(jj_index, ii_index), at(jj_index, ii_minus));
                for k in 0..3 {
                    value += (pji[k] - pjm[k]).abs();
                }
                value *= inverse;
                if max_derivative[1] < value {
                    max_derivative[1] = value;
                }
            }
        }
    }
    max_derivative[1] *= v_degree as f64;
    max_derivative[0] *= SQRT_2;
    max_derivative[1] *= SQRT_2;
    if max_derivative[0] != 0.0 && max_derivative[1] != 0.0 {
        (tolerance_3d / max_derivative[0], tolerance_3d / max_derivative[1])
    } else {
        (0.0, 0.0)
    }
}

/// `Geom_BSplineSurface::Resolution(Tolerance3D, UTolerance, VTolerance)`:
/// `(UTolerance, VTolerance)`.
#[doc(alias = "Geom_BSplineSurface::Resolution")]
pub(crate) fn bspline_surface_resolution(surface: &BSplineSurface, tolerance_3d: f64) -> (f64, f64) {
    let (u, v) = (surface.u_knot_vector(), surface.v_knot_vector());
    let (u_rational, v_rational) = surface.weights().map_or((false, false), surface_weights_rational);
    let weights = if u_rational || v_rational { surface.weights() } else { None };
    let (u_inv, v_inv) = bsplslib_resolution(
        surface.poles(),
        weights,
        &knot_sequence(&u.knots, &u.multiplicities, u.degree, u.periodic),
        &knot_sequence(&v.knots, &v.multiplicities, v.degree, v.periodic),
        u.degree,
        v.degree,
        u_rational,
        v_rational,
        1.0,
    );
    (tolerance_3d * u_inv, tolerance_3d * v_inv)
}

/// `Geom_BezierSurface::Resolution(Tolerance3D, UTolerance, VTolerance)`:
/// `(UTolerance, VTolerance)`.
#[doc(alias = "Geom_BezierSurface::Resolution")]
pub(crate) fn bezier_surface_resolution(surface: &BezierSurface, tolerance_3d: f64) -> (f64, f64) {
    let (u_degree, v_degree) = (surface.poles().len() - 1, surface.poles()[0].len() - 1);
    let (u_rational, v_rational) = surface.weights().map_or((false, false), surface_weights_rational);
    let weights = if u_rational || v_rational { surface.weights() } else { None };
    // (OCCT passes `UKnots()` for both directions: both are `[0, 1]`.)
    let (u_inv, v_inv) = bsplslib_resolution(
        surface.poles(),
        weights,
        &bezier_flat_knots(u_degree),
        &bezier_flat_knots(v_degree),
        u_degree,
        v_degree,
        u_rational,
        v_rational,
        1.0,
    );
    (tolerance_3d * u_inv, tolerance_3d * v_inv)
}

/// `GeomAdaptor_Curve::Resolution(R3D)`: a trimmed curve's is its basis', an offset
/// curve's, a hyperbola's and a parabola's `Precision::Parametric(R3D)`.
#[doc(alias = "GeomAdaptor_Curve::Resolution")]
pub(crate) fn curve_resolution(curve: &Curve, r3d: f64) -> f64 {
    match curve {
        Curve::Trimmed(t) => curve_resolution(t.basis(), r3d),
        Curve::Line(_) => r3d,
        Curve::Circle(c) => {
            let r = c.radius();
            if r > r3d / 2.0 { 2.0 * (r3d / (2.0 * r)).asin() } else { 2.0 * PI }
        }
        Curve::Ellipse(e) => r3d / e.major_radius(),
        Curve::Bezier(b) => bezier_curve_resolution(b, r3d),
        Curve::BSpline(b) => bspline_curve_resolution(b, r3d),
        Curve::Hyperbola(_) | Curve::Parabola(_) | Curve::Offset(_) => parametric(r3d),
    }
}

/// The radius of the cone's circle at `v` (`Geom_ConicalSurface::VIso(V)`'s, through
/// `ElSLib::ConeVIso`: made positive).
fn cone_v_iso_radius(cone: &Cone, v: f64) -> f64 {
    let r = cone.radius() + v * cone.semi_angle().sin();
    if r < 0.0 { -r } else { r }
}

/// `2 * asin(Res)` (`2π` from `Res > 1`), the elementary surfaces' tail of `UResolution`
/// and `VResolution`.
fn angular(res: f64) -> f64 {
    if res <= 1.0 { 2.0 * res.asin() } else { 2.0 * PI }
}

/// `GeomAdaptor_Surface::UResolution(R3d)` of an adaptor on `surface` over `bounds`
/// (`(u1, u2, v1, v2)`, the adaptor's, normally the surface's own; only a cone's v bounds
/// matter). A trimmed surface's is its basis' over the same bounds (the adaptor loads the
/// basis), an offset surface's its basis', a surface of revolution's
/// `Precision::Parametric(R3d)`.
#[doc(alias = "GeomAdaptor_Surface::UResolution")]
pub(crate) fn u_resolution(surface: &Surface, bounds: (f64, f64, f64, f64), r3d: f64) -> f64 {
    let mut res = 0.0;
    match surface {
        Surface::Trimmed(t) => return u_resolution(t.basis(), bounds, r3d),
        Surface::Extrusion(e) => return curve_resolution(e.basis(), r3d),
        Surface::Torus(t) => {
            let r = t.major_radius() + t.minor_radius();
            if r > CONFUSION {
                res = r3d / (2.0 * r);
            }
        }
        Surface::Sphere(s) => {
            let r = s.radius();
            if r > CONFUSION {
                res = r3d / (2.0 * r);
            }
        }
        Surface::Cylinder(c) => {
            let r = c.radius();
            if r > CONFUSION {
                res = r3d / (2.0 * r);
            }
        }
        Surface::Cone(c) => {
            let (_, _, v_first, v_last) = bounds;
            if v_last - v_first > 1.0e10 {
                // Not truly bounded: unknown resolution.
                return parametric(r3d);
            }
            let rayon1 = cone_v_iso_radius(c, v_last);
            let rayon2 = cone_v_iso_radius(c, v_first);
            let r = if rayon1 > rayon2 { rayon1 } else { rayon2 };
            return if r > CONFUSION { r3d / r } else { 0.0 };
        }
        Surface::Plane(_) => return r3d,
        Surface::Bezier(b) => return bezier_surface_resolution(b, r3d).0,
        Surface::BSpline(b) => return bspline_surface_resolution(b, r3d).0,
        Surface::Offset(o) => return u_resolution(o.basis(), bounds, r3d),
        Surface::Revolution(_) => return parametric(r3d),
    }
    angular(res)
}

/// `GeomAdaptor_Surface::VResolution(R3d)` of an adaptor on `surface` over `bounds` (as
/// [`u_resolution`]'s). A surface of revolution's is its meridian's, an extrusion's,
/// cylinder's, cone's and plane's `R3d`.
#[doc(alias = "GeomAdaptor_Surface::VResolution")]
pub(crate) fn v_resolution(surface: &Surface, bounds: (f64, f64, f64, f64), r3d: f64) -> f64 {
    let mut res = 0.0;
    match surface {
        Surface::Trimmed(t) => return v_resolution(t.basis(), bounds, r3d),
        Surface::Revolution(r) => return curve_resolution(r.basis(), r3d),
        Surface::Torus(t) => {
            let r = t.minor_radius();
            if r > CONFUSION {
                res = r3d / (2.0 * r);
            }
        }
        Surface::Sphere(s) => {
            let r = s.radius();
            if r > CONFUSION {
                res = r3d / (2.0 * r);
            }
        }
        Surface::Extrusion(_) | Surface::Cylinder(_) | Surface::Cone(_) | Surface::Plane(_) => return r3d,
        Surface::Bezier(b) => return bezier_surface_resolution(b, r3d).1,
        Surface::BSpline(b) => return bspline_surface_resolution(b, r3d).1,
        Surface::Offset(o) => return v_resolution(o.basis(), bounds, r3d),
    }
    angular(res)
}

#[cfg(test)]
mod tests {
    use std::f64::consts::PI;

    use super::*;
    use crate::geometry::curve::{Circle, Ellipse, Hyperbola, Line, OffsetCurve, TrimmedCurve};
    use crate::geometry::bspline::KnotVector;
    use crate::geometry::surface::{Cylinder, ExtrusionSurface, OffsetSurface, Plane, RevolutionSurface, Sphere, Torus, TrimmedSurface};
    use crate::math::{Axis, Direction, Frame};

    /// The 3D tolerances each case is asked about.
    const R3D: [f64; 4] = [1e-7, 0.01, 0.7, 10.0];

    // OCCT 8.1's `GeomAdaptor_Curve::Resolution` and `GeomAdaptor_Surface::UResolution` /
    // `VResolution` of the same curves and surfaces (a throwaway driver printing them).
    #[allow(clippy::excessive_precision)]
    const CURVES: &[(&str, [f64; 4])] = &[
        ("bspline", [7.1428571428571438e-09, 0.00071428571428571439, 0.050000000000000003, 0.71428571428571441]),
        ("bspline_rational", [8.6206896551724135e-10, 8.6206896551724145e-05, 0.0060344827586206896, 0.086206896551724144]),
        ("bspline_unit_weights", [7.1428571428571438e-09, 0.00071428571428571439, 0.050000000000000003, 0.71428571428571441]),
        ("bspline_periodic", [1.3333333333333332e-08, 0.0013333333333333333, 0.093333333333333324, 1.3333333333333333]),
        ("bspline_periodic_rational", [5.0000000000000001e-09, 0.00050000000000000001, 0.034999999999999996, 0.5]),
        ("bspline_periodic_mult2", [2.8985507246376803e-09, 0.00028985507246376805, 0.020289855072463763, 0.28985507246376807]),
        ("bezier", [3.333333333333333e-09, 0.00033333333333333332, 0.023333333333333331, 0.33333333333333331]),
        ("bezier_rational", [1.461988304093567e-09, 0.00014619883040935673, 0.01023391812865497, 0.14619883040935672]),
        ("circle", [5.0000000000000004e-08, 0.0050000052083479822, 0.35181153632743256, 6.2831853071795862]),
        ("circle_small", [2.0000000000000034e-07, 0.020000333348334228, 1.5507949932215059, 6.2831853071795862]),
        ("ellipse", [2e-08, 0.002, 0.13999999999999999, 2.0]),
        ("line", [9.9999999999999995e-08, 0.01, 0.69999999999999996, 10.0]),
        ("hyperbola", [1.0000000000000001e-09, 0.0001, 0.0069999999999999993, 0.10000000000000001]),
        ("trimmed_bspline", [7.1428571428571438e-09, 0.00071428571428571439, 0.050000000000000003, 0.71428571428571441]),
        ("offset", [1.0000000000000001e-09, 0.0001, 0.0069999999999999993, 0.10000000000000001]),
    ];
    const SURFACES: &[(&str, [(f64, f64); 4])] = &[
        ("bspline", [(1.2740662724081939e-08, 4.4194173824159215e-09), (0.0012740662724081938, 0.00044194173824159215), (0.08918463906857356, 0.030935921676911449), (1.2740662724081939, 0.44194173824159216)]),
        ("bspline_rational", [(2.9648083068618343e-09, 1.931985740946851e-09), (0.00029648083068618342, 0.0001931985740946851), (0.020753658148032837, 0.013523900186627956), (0.29648083068618342, 0.1931985740946851)]),
        ("bspline_rational_rows", [(1.2740662724081939e-08, 1.4731391274719735e-09), (0.0012740662724081938, 0.00014731391274719734), (0.08918463906857356, 0.010311973892303814), (1.2740662724081939, 0.14731391274719735)]),
        ("bspline_rational_cols", [(5.0962650896327724e-09, 4.4194173824159215e-09), (0.00050962650896327723, 0.00044194173824159215), (0.035673855627429404, 0.030935921676911449), (0.50962650896327721, 0.44194173824159216)]),
        ("bspline_periodic", [(9.3040365945598352e-09, 1.0398629135096284e-08), (0.00093040365945598364, 0.0010398629135096285), (0.065128256161918849, 0.072790403945673987), (0.93040365945598369, 1.0398629135096285)]),
        ("bspline_periodic_rational", [(5.7723002545840605e-09, 4.1594516540385136e-09), (0.00057723002545840608, 0.0004159451654038514), (0.040406101782088429, 0.029116161578269593), (0.57723002545840607, 0.4159451654038514)]),
        ("bspline_periodic_v_mult2", [(6.0957481136771307e-09, 1.0545962433803836e-09), (0.00060957481136771313, 0.00010545962433803838), (0.042670236795739912, 0.007382173703662686), (0.60957481136771308, 0.10545962433803838)]),
        ("bezier", [(1.1785113019775788e-08, 1.0878565864408422e-08), (0.0011785113019775787, 0.0010878565864408422), (0.082495791138430516, 0.076149961050858958), (1.1785113019775788, 1.0878565864408423)]),
        ("bezier_rational", [(4.6673714929805104e-09, 5.8438576957565901e-09), (0.00046673714929805108, 0.00058438576957565904), (0.032671600450863571, 0.040907003870296131), (0.46673714929805105, 0.58438576957565913)]),
        ("plane", [(9.9999999999999995e-08, 9.9999999999999995e-08), (0.01, 0.01), (0.69999999999999996, 0.69999999999999996), (10.0, 10.0)]),
        ("cylinder", [(3.3333333333333334e-08, 9.9999999999999995e-08), (0.0033333348765451393, 0.01), (0.23386592292475683, 0.69999999999999996), (6.2831853071795862, 10.0)]),
        ("cone_infinite", [(1.0000000000000001e-09, 9.9999999999999995e-08), (0.0001, 0.01), (0.0069999999999999993, 0.69999999999999996), (0.10000000000000001, 10.0)]),
        ("cone_bounded", [(3.1425977355019886e-08, 9.9999999999999995e-08), (0.0031425977355019888, 0.01), (0.21998184148513919, 0.69999999999999996), (3.142597735501989, 10.0)]),
        ("cone_bounded_negative", [(4.3563110318006059e-08, 9.9999999999999995e-08), (0.0043563110318006063, 0.01), (0.30494177222604241, 0.69999999999999996), (4.3563110318006064, 10.0)]),
        ("cone_trimmed", [(3.4643305014170685e-08, 9.9999999999999995e-08), (0.0034643305014170687, 0.01), (0.2425031350991948, 0.69999999999999996), (3.4643305014170687, 10.0)]),
        ("sphere", [(6.6666666666666681e-08, 6.6666666666666681e-08), (0.0066666790124074084, 0.0066666790124074084), (0.47100847344159952, 0.47100847344159952), (6.2831853071795862, 6.2831853071795862)]),
        ("torus", [(1.6666666666666667e-08, 1.0000000000000004e-07), (0.0016666668595679615, 0.010000041667135424), (0.11673293331115564, 0.71514220729102052), (1.9702215666754914, 6.2831853071795862)]),
        ("torus_small", [(3.3333333333333334e-08, 1.0000000000000004e-07), (0.0033333348765451393, 0.010000041667135424), (0.23386592292475683, 0.71514220729102052), (6.2831853071795862, 6.2831853071795862)]),
        ("revolution_bspline", [(1.0000000000000001e-09, 7.1428571428571438e-09), (0.0001, 0.00071428571428571439), (0.0069999999999999993, 0.050000000000000003), (0.10000000000000001, 0.71428571428571441)]),
        ("revolution_circle", [(1.0000000000000001e-09, 1.0000000000000004e-07), (0.0001, 0.010000041667135424), (0.0069999999999999993, 0.71514220729102052), (0.10000000000000001, 6.2831853071795862)]),
        ("extrusion_ellipse", [(2e-08, 9.9999999999999995e-08), (0.002, 0.01), (0.13999999999999999, 0.69999999999999996), (2.0, 10.0)]),
        ("extrusion_bspline", [(7.1428571428571438e-09, 9.9999999999999995e-08), (0.00071428571428571439, 0.01), (0.050000000000000003, 0.69999999999999996), (0.71428571428571441, 10.0)]),
        ("offset_bspline", [(2.9648083068618343e-09, 1.931985740946851e-09), (0.00029648083068618342, 0.0001931985740946851), (0.020753658148032837, 0.013523900186627956), (0.29648083068618342, 0.1931985740946851)]),
        ("offset_cylinder", [(3.3333333333333334e-08, 9.9999999999999995e-08), (0.0033333348765451393, 0.01), (0.23386592292475683, 0.69999999999999996), (6.2831853071795862, 10.0)]),
    ];

    fn bspline_curve(poles: &[(f64, f64, f64)], weights: Option<Vec<f64>>, knots: Vec<f64>, multiplicities: Vec<usize>, periodic: bool) -> BSplineCurve {
        let poles: Vec<Point> = poles.iter().map(|&(x, y, z)| Point::new(x, y, z)).collect();
        let curve = if periodic { BSplineCurve::new_periodic(poles, knots, multiplicities, 3) } else { BSplineCurve::new(poles, knots, multiplicities, 3) };
        match weights {
            Some(w) => curve.with_weights(w),
            None => curve,
        }
    }

    const P6: [(f64, f64, f64); 6] = [(0.0, 0.0, 0.0), (1.0, 2.0, 0.0), (3.0, 3.0, 1.0), (5.0, 1.0, 2.0), (7.0, 4.0, 0.0), (9.0, 0.0, -1.0)];
    const P5: [(f64, f64, f64); 5] = [(0.0, 0.0, 0.0), (4.0, 1.0, 0.0), (5.0, 5.0, 1.0), (1.0, 6.0, 2.0), (-2.0, 3.0, 0.5)];

    fn c1() -> BSplineCurve {
        bspline_curve(&P6, None, vec![0.0, 1.0, 2.5, 4.0], vec![4, 1, 1, 4], false)
    }

    fn frame() -> Frame {
        Frame::new(Point::new(1.0, 2.0, 3.0), Direction::new(0.0, 0.6, 0.8), Direction::new(1.0, 0.0, 0.0))
    }

    fn curve(name: &str) -> Curve {
        let p5 = || P5.iter().map(|&(x, y, z)| Point::new(x, y, z)).collect::<Vec<_>>();
        match name {
            "bspline" => Curve::BSpline(c1()),
            "bspline_rational" => Curve::BSpline(bspline_curve(&P6, Some(vec![1.0, 2.0, 0.5, 1.5, 1.0, 3.0]), vec![0.0, 1.0, 2.5, 4.0], vec![4, 1, 1, 4], false)),
            "bspline_unit_weights" => Curve::BSpline(bspline_curve(&P6, Some(vec![2.0; 6]), vec![0.0, 1.0, 2.5, 4.0], vec![4, 1, 1, 4], false)),
            "bspline_periodic" => Curve::BSpline(bspline_curve(&P5, None, vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0], vec![1; 6], true)),
            "bspline_periodic_rational" => Curve::BSpline(bspline_curve(&P5, Some(vec![1.0, 2.0, 1.0, 0.7, 1.3]), vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0], vec![1; 6], true)),
            "bspline_periodic_mult2" => Curve::BSpline(bspline_curve(&P5, Some(vec![1.0, 2.0, 1.0, 0.7, 1.3]), vec![0.0, 1.0, 2.0, 3.0], vec![2, 1, 2, 2], true)),
            "bezier" => Curve::Bezier(BezierCurve::new(p5())),
            "bezier_rational" => Curve::Bezier(BezierCurve::rational(p5(), (0..5).map(|i| 1.0 + 0.3 * i as f64).collect())),
            "circle" => Curve::Circle(Circle::new(frame(), 2.0)),
            "circle_small" => Curve::Circle(Circle::new(frame(), 0.5)),
            "ellipse" => Curve::Ellipse(Ellipse::new(frame(), 5.0, 2.0)),
            "line" => Curve::Line(Line::new(Point::new(1.0, 2.0, 3.0), Direction::new(0.0, 0.6, 0.8))),
            "hyperbola" => Curve::Hyperbola(Hyperbola::new(frame(), 3.0, 1.0)),
            "trimmed_bspline" => Curve::Trimmed(TrimmedCurve::new(Curve::BSpline(c1()), 0.5, 3.0)),
            "offset" => Curve::Offset(OffsetCurve::new(Curve::Circle(Circle::new(frame(), 2.0)), 0.5, Direction::new(0.0, 0.6, 0.8))),
            _ => panic!("no curve {name}"),
        }
    }

    /// The test grid of poles (as the OCCT driver's).
    fn grid(i: usize, j: usize) -> Point {
        let (di, dj) = (i as f64, j as f64);
        Point::new(di + 0.3 * dj * dj, 1.5 * dj + 0.2 * di * dj, 0.1 * di * di - 0.4 * dj + 0.05 * di * dj)
    }

    /// Weights: 0 none, 1 varying both ways, 2 along u (the first index), 3 along v.
    fn grid_weight(mode: u32, i: usize, j: usize) -> f64 {
        match mode {
            1 => 1.0 + 0.25 * ((i + 2 * j) % 3) as f64,
            2 => 1.0 + 0.5 * i as f64,
            3 => 1.0 + 0.5 * j as f64,
            _ => 1.0,
        }
    }

    fn bspline_surface(nu: usize, nv: usize, u: KnotVector, v: KnotVector, mode: u32) -> BSplineSurface {
        let poles = (0..nu).map(|i| (0..nv).map(|j| grid(i, j)).collect()).collect();
        let surface = BSplineSurface::new(poles, u, v);
        if mode == 0 { surface } else { surface.with_weights((0..nu).map(|i| (0..nv).map(|j| grid_weight(mode, i, j)).collect()).collect()) }
    }

    fn plain_bspline_surface(mode: u32) -> BSplineSurface {
        bspline_surface(5, 4, KnotVector::new(vec![0.0, 1.0, 2.5], vec![4, 1, 4], 3), KnotVector::new(vec![0.0, 0.5, 1.0], vec![3, 1, 3], 2), mode)
    }

    /// The surface and the adaptor's bounds.
    fn surface(name: &str) -> (Surface, (f64, f64, f64, f64)) {
        let cone = || Cone::new(frame(), 0.3, 2.0);
        let surface = match name {
            "bspline" => Surface::BSpline(plain_bspline_surface(0)),
            "bspline_rational" => Surface::BSpline(plain_bspline_surface(1)),
            "bspline_rational_rows" => Surface::BSpline(plain_bspline_surface(2)),
            "bspline_rational_cols" => Surface::BSpline(plain_bspline_surface(3)),
            "bspline_periodic" | "bspline_periodic_rational" => {
                let mode = if name == "bspline_periodic" { 0 } else { 1 };
                Surface::BSpline(bspline_surface(5, 3, KnotVector::periodic(vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0], vec![1; 6], 3), KnotVector::new(vec![0.0, 1.0], vec![3, 3], 2), mode))
            }
            "bspline_periodic_v_mult2" => Surface::BSpline(bspline_surface(4, 5, KnotVector::new(vec![0.0, 1.5, 3.0], vec![3, 1, 3], 2), KnotVector::periodic(vec![0.0, 1.0, 2.0, 3.0], vec![2, 1, 2, 2], 3), 1)),
            "bezier" | "bezier_rational" => {
                let poles: Vec<Vec<Point>> = (0..4).map(|i| (0..3).map(|j| grid(i, j)).collect()).collect();
                if name == "bezier" {
                    Surface::Bezier(BezierSurface::new(poles))
                } else {
                    Surface::Bezier(BezierSurface::rational(poles, (0..4).map(|i| (0..3).map(|j| grid_weight(1, i, j)).collect()).collect()))
                }
            }
            "plane" => Surface::Plane(Plane::from_frame(frame())),
            "cylinder" => Surface::Cylinder(Cylinder::new(frame(), 3.0)),
            "cone_infinite" => Surface::Cone(cone()),
            "cone_bounded" => return (Surface::Cone(cone()), (0.0, 2.0 * PI, -1.0, 4.0)),
            "cone_bounded_negative" => return (Surface::Cone(cone()), (0.0, 2.0 * PI, -10.0, 1.0)),
            "cone_trimmed" => Surface::Trimmed(TrimmedSurface::new(Surface::Cone(cone()), 0.0, PI, 0.5, 3.0)),
            "sphere" => Surface::Sphere(Sphere::new(frame(), 1.5)),
            "torus" => Surface::Torus(Torus::new(frame(), 5.0, 1.0)),
            "torus_small" => Surface::Torus(Torus::new(frame(), 2.0, 1.0)),
            "revolution_bspline" => Surface::Revolution(RevolutionSurface::new(Curve::BSpline(c1()), Axis::new(Point::ORIGIN, Direction::new(0.0, 0.0, 1.0)))),
            "revolution_circle" => {
                let meridian = Circle::new(Frame::new(Point::new(5.0, 0.0, 0.0), Direction::new(0.0, 1.0, 0.0), Direction::new(1.0, 0.0, 0.0)), 1.0);
                Surface::Revolution(RevolutionSurface::new(Curve::Circle(meridian), Axis::new(Point::ORIGIN, Direction::new(0.0, 0.0, 1.0))))
            }
            "extrusion_ellipse" => Surface::Extrusion(ExtrusionSurface::new(Curve::Ellipse(Ellipse::new(frame(), 5.0, 2.0)), Direction::new(0.0, 0.6, 0.8))),
            "extrusion_bspline" => Surface::Extrusion(ExtrusionSurface::new(Curve::BSpline(c1()), Direction::new(0.0, 0.6, 0.8))),
            "offset_bspline" => Surface::Offset(OffsetSurface::new(Surface::BSpline(plain_bspline_surface(1)), 0.5)),
            "offset_cylinder" => Surface::Offset(OffsetSurface::new(Surface::Cylinder(Cylinder::new(frame(), 3.0)), 0.5)),
            _ => panic!("no surface {name}"),
        };
        use crate::geometry::ParametricSurface;
        let bounds = surface.bounds();
        (surface, bounds)
    }

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() <= 1e-12 * b.abs()
    }

    #[test]
    fn curve_resolutions_match_occt() {
        for (name, expected) in CURVES {
            let c = curve(name);
            for (&r3d, &want) in R3D.iter().zip(expected) {
                let got = curve_resolution(&c, r3d);
                assert!(close(got, want), "{name} at {r3d}: {got:e}, OCCT {want:e}");
            }
        }
    }

    #[test]
    fn surface_resolutions_match_occt() {
        for (name, expected) in SURFACES {
            let (s, bounds) = surface(name);
            for (&r3d, &(want_u, want_v)) in R3D.iter().zip(expected) {
                let (u, v) = (u_resolution(&s, bounds, r3d), v_resolution(&s, bounds, r3d));
                assert!(close(u, want_u), "{name} u at {r3d}: {u:e}, OCCT {want_u:e}");
                assert!(close(v, want_v), "{name} v at {r3d}: {v:e}, OCCT {want_v:e}");
            }
        }
    }

    /// OCCT takes a surface's weights as rational in v when they differ along u (the first
    /// index), and the other way round.
    #[test]
    fn surface_rational_flags_as_occt() {
        let rows = vec![vec![1.0, 1.0], vec![2.0, 2.0]];
        assert_eq!(surface_weights_rational(&rows), (false, true));
        let cols = vec![vec![1.0, 2.0], vec![1.0, 2.0]];
        assert_eq!(surface_weights_rational(&cols), (true, false));
    }
}
