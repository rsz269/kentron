//! Local differential properties of curves and surfaces: tangents, curvature, normals.

use crate::math::{Direction, Direction2, Point2, Vector2};

use super::curve2::Curve2;
use super::parametric::{ParametricCurve, ParametricSurface};
use super::surface::Surface;

/// Where `f` changes sign over the curve's domain (a window round 0 for unbounded
/// curves), each with the sampling step; none if `f` vanishes everywhere.
fn roots(curve: &Curve2, f: &dyn Fn(f64) -> f64) -> Vec<(f64, f64)> {
    let (mut a, mut b) = (curve.first_parameter(), curve.last_parameter());
    if !(a.is_finite() && a.abs() < 1e50) {
        a = -60.0;
    }
    if !(b.is_finite() && b.abs() < 1e50) {
        b = 60.0;
    }
    let period = curve.period().filter(|_| curve.is_periodic());
    const N: usize = 720;
    let h = (b - a) / N as f64;
    // A periodic curve is sampled half a step off, so that roots at its origin fall
    // between samples.
    let shift = if period.is_some() { -0.5 * h } else { 0.0 };
    let ts: Vec<f64> = (0..=N).map(|i| a + shift + h * i as f64).collect();
    let fs: Vec<f64> = ts.iter().map(|&t| f(t)).collect();
    let scale = fs.iter().fold(0.0f64, |m, v| m.max(v.abs()));
    if !(scale > 1e-12) {
        return Vec::new();
    }
    let mut out: Vec<(f64, f64)> = Vec::new();
    for i in 0..=N {
        // A root right at a sample (where the sign changes across it).
        if fs[i] == 0.0 && i > 0 && i < N && fs[i - 1] * fs[i + 1] < 0.0 {
            out.push((ts[i], h));
            continue;
        }
        if i == N {
            break;
        }
        let (f0, f1) = (fs[i], fs[i + 1]);
        if !(f0 * f1 < 0.0) {
            continue;
        }
        let (mut s, mut t, mut fs0) = (ts[i], ts[i + 1], f0);
        for _ in 0..100 {
            let m = 0.5 * (s + t);
            let fm = f(m);
            if fm == 0.0 {
                (s, t) = (m, m);
                break;
            }
            if (fm < 0.0) == (fs0 < 0.0) {
                (s, fs0) = (m, fm);
            } else {
                t = m;
            }
        }
        let mut root = 0.5 * (s + t);
        if let Some(p) = period {
            root = a + (root - a).rem_euclid(p);
            if (root - (a + p)).abs() <= 1e-12 * p {
                root = a;
            }
        }
        // Exactly at a domain end, a root is at the end.
        if root.abs() <= 1e-12 {
            root = 0.0;
        }
        out.push((root, h));
    }
    out.sort_by(|x, y| x.0.total_cmp(&y.0));
    out.dedup_by(|x, y| (x.0 - y.0).abs() <= 1e-9);
    out
}

/// Local properties of a 2D curve at one parameter. Quantities that are undefined there
/// (a tangent where the derivatives vanish, a normal or center of curvature where the
/// curvature is zero) are `None`.
#[doc(alias = "GeomLProp_CLProps2d")]
#[doc(alias = "Geom2dLProp_CLProps2d")]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CurvePoint2 {
    pub parameter: f64,
    pub point: Point2,
    pub d1: Vector2,
    pub d2: Vector2,
    /// The unit tangent, from the first non-vanishing derivative.
    pub tangent: Option<Direction2>,
    /// Unsigned curvature (0 if the tangent is undefined).
    pub curvature: f64,
    /// The principal normal, pointing towards the center of curvature.
    pub normal: Option<Direction2>,
    pub center_of_curvature: Option<Point2>,
}

/// The local properties of `curve` at `t`; derivatives shorter than `tolerance` count as
/// vanishing.
#[doc(alias = "GeomLProp_CLProps2d")]
pub fn curve_point2(curve: &Curve2, t: f64, tolerance: f64) -> CurvePoint2 {
    let (point, d1, d2, d3) = curve.d3(t);
    // The tangent from the first derivative that does not vanish.
    let tangent = [d1, d2, d3].into_iter().find(|d| d.length() > tolerance).and_then(|d| Direction2::try_new(d.x, d.y));
    let speed = d1.length();
    let cross = d1.x * d2.y - d1.y * d2.x;
    let curvature = if speed > tolerance { (cross / (speed * speed * speed)).abs() } else { 0.0 };
    // The normal turns from the tangent towards the side the curve bends to.
    let normal = match tangent {
        Some(t) if curvature > tolerance.min(1e-12) && cross != 0.0 => {
            Some(if cross > 0.0 { Direction2::new(-t.y(), t.x()) } else { Direction2::new(t.y(), -t.x()) })
        }
        _ => None,
    };
    let center_of_curvature = normal.map(|n| Point2::new(point.x + n.x() / curvature, point.y + n.y() / curvature));
    CurvePoint2 { parameter: t, point, d1, d2, tangent, curvature, normal, center_of_curvature }
}

/// The kind of a [`CurvatureFeature`].
#[doc(alias = "LProp_CIType")]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum CurvatureFeatureKind {
    /// The curvature changes sign.
    Inflection,
    /// A local minimum of the radius of curvature (a curvature maximum, e.g. the apex of
    /// a parabola).
    #[doc(alias = "LProp_MinCur")]
    MinRadius,
    /// A local maximum of the radius of curvature (a curvature minimum).
    #[doc(alias = "LProp_MaxCur")]
    MaxRadius,
}

/// A point of interest for a curve's curvature.
#[doc(alias = "LProp_CurAndInf")]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CurvatureFeature {
    pub parameter: f64,
    pub kind: CurvatureFeatureKind,
}

/// The inflections and curvature extrema of a 2D curve over its domain, sorted by
/// parameter. Curves of constant curvature (lines, circles) have none.
#[doc(alias = "GeomLProp_CurAndInf2d::Perform")]
#[doc(alias = "Geom2dLProp_CurAndInf2d::Perform")]
pub fn curvature_features2(curve: &Curve2) -> Vec<CurvatureFeature> {
    let mut features: Vec<CurvatureFeature> = inflections2(curve)
        .into_iter()
        .map(|parameter| CurvatureFeature { parameter, kind: CurvatureFeatureKind::Inflection })
        .collect();
    features.extend(curvature_extrema2(curve));
    features.sort_by(|a, b| a.parameter.total_cmp(&b.parameter));
    features
}

/// The parameters of the inflections of a 2D curve, increasing.
#[doc(alias = "GeomLProp_CurAndInf2d::PerformInf")]
pub fn inflections2(curve: &Curve2) -> Vec<f64> {
    // Conics bend one way throughout.
    if let Curve2::Line(_) | Curve2::Circle(_) | Curve2::Ellipse(_) | Curve2::Parabola(_) | Curve2::Hyperbola(_) = curve {
        return Vec::new();
    }
    // Where d1 × d2 changes sign.
    let cross = |t: f64| {
        let (_, d1, d2) = curve.d2(t);
        d1.x * d2.y - d1.y * d2.x
    };
    roots(curve, &cross).into_iter().map(|(t, _)| t).collect()
}

/// The curvature extrema of a 2D curve, sorted by parameter.
#[doc(alias = "GeomLProp_CurAndInf2d::PerformCurExt")]
pub fn curvature_extrema2(curve: &Curve2) -> Vec<CurvatureFeature> {
    // A parabola's or hyperbola's only one is its vertex (where the curvature underflows
    // far along the branches, sampling would find noise).
    if let Curve2::Parabola(_) | Curve2::Hyperbola(_) = curve {
        return vec![CurvatureFeature { parameter: 0.0, kind: CurvatureFeatureKind::MinRadius }];
    }
    // Where the signed curvature's derivative changes sign; none if it is constant.
    let signed = |t: f64| {
        let (_, d1, d2) = curve.d2(t);
        let s = d1.length();
        (d1.x * d2.y - d1.y * d2.x) / (s * s * s)
    };
    let slope = |t: f64| {
        let (_, d1, d2, d3) = curve.d3(t);
        let s2 = d1.dot(d1);
        let cross12 = d1.x * d2.y - d1.y * d2.x;
        let cross13 = d1.x * d3.y - d1.y * d3.x;
        (cross13 * s2 - 3.0 * cross12 * d1.dot(d2)) / (s2 * s2 * s2.sqrt())
    };
    roots(curve, &slope)
        .into_iter()
        .map(|(t, h)| {
            let here = signed(t).abs();
            let kind = if here >= signed(t - h).abs() && here >= signed(t + h).abs() {
                CurvatureFeatureKind::MinRadius
            } else {
                CurvatureFeatureKind::MaxRadius
            };
            CurvatureFeature { parameter: t, kind }
        })
        .collect()
}

/// A surface normal, possibly estimated.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SurfaceNormal {
    /// The normal from the first derivatives.
    Regular(Direction),
    /// The first derivatives are parallel or vanish there (e.g. at the pole of a sphere);
    /// the normal is the limit from higher derivatives or nearby points.
    Singular(Direction),
}

impl SurfaceNormal {
    pub fn direction(&self) -> Direction {
        match *self {
            Self::Regular(d) | Self::Singular(d) => d,
        }
    }
}

/// The normal of `surface` at `(u, v)`, estimated at singular points; `None` if it cannot
/// be determined. Derivatives shorter than `tolerance` count as vanishing.
#[doc(alias = "GeomLib::NormEstim")]
pub fn surface_normal(surface: &Surface, u: f64, v: f64, tolerance: f64) -> Option<SurfaceNormal> {
    let (_, du, dv) = surface.d1(u, v);
    let n = du.cross(dv);
    if du.length() > tolerance && dv.length() > tolerance && n.length() > tolerance * du.length().max(dv.length()) {
        return Direction::try_from_vector(n).map(SurfaceNormal::Regular);
    }
    // Singular: the normal of nearby points, stepping into the domain.
    let (u0, u1, v0, v1) = surface.bounds();
    let inward = |x: f64, lo: f64, hi: f64, h: f64| if x + h <= hi || !hi.is_finite() { x + h } else if x - h >= lo { x - h } else { x };
    for h in [1e-8, 1e-6, 1e-4, 1e-3] {
        let (un, vn) = (inward(u, u0, u1, h), inward(v, v0, v1, h));
        for (uu, vv) in [(u, vn), (un, v), (un, vn)] {
            let (_, a, b) = surface.d1(uu, vv);
            let m = a.cross(b);
            if m.length() > tolerance * a.length().max(b.length()) && m.length() > 0.0 {
                return Direction::try_from_vector(m).map(SurfaceNormal::Singular);
            }
        }
    }
    None
}
