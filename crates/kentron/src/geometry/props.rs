//! Local differential properties of curves and surfaces: tangents, curvature, normals.

use crate::math::{Direction, Direction2, Point2, Vector, Vector2};

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

/// `GeomLib::NormEstim`'s step into the parameters off a singular point.
const NORMAL_PARAMETER_STEP: f64 = 1.0e-5;
/// `GeomLib::NormEstim`'s least squared length of a normal found off a singular point.
const NORMAL_SQUARED_MAGNITUDE_TOL: f64 = 1.0e-16;
/// `GeomLib::NormEstim`'s angle within which the normals either side of a cone's apex are
/// opposite.
const CONE_NORMAL_ANGULAR_TOLERANCE: f64 = 1.0e-4;

/// `CSLib::Normal` from the second derivatives, at a point where the first ones are
/// parallel or vanish: the normal as the limit of the normals nearby (along the derivatives
/// of `D1U ^ D1V`), unless they turn round the point (`CSLib_InfinityOfSolutions`) or
/// vanish too.
#[doc(alias = "CSLib::Normal")]
fn limit_normal(d: &super::SurfaceD2, sin_tol: f64) -> Option<Direction> {
    let d1nu = d.duu.cross(d.dv) + d.du.cross(d.duv);
    let d1nv = d.duv.cross(d.dv) + d.du.cross(d.dvv);
    let (lu, lv) = (d1nu.length_squared(), d1nv.length_squared());
    let eps = f64::EPSILON;
    if lu <= eps && lv <= eps {
        return None;
    }
    if lu < eps {
        return Direction::try_from_vector(d1nv);
    }
    if lv < eps {
        return Direction::try_from_vector(d1nu);
    }
    if lv / lu <= eps || lu / lv <= eps {
        return None;
    }
    let sin2 = d1nu.cross(d1nv).length_squared() / (lu * lv);
    if sin2 < sin_tol * sin_tol {
        return Direction::try_from_vector(d1nu);
    }
    None
}

/// The unit normal of `surface` at `(u, v)` (along `D1U ^ D1V`), as OCCT's
/// `GeomLib::NormEstim` finds it for a mesh's nodes: from the first derivatives where
/// neither is shorter than `tolerance`; where one is (as at a sphere's pole), the limit of
/// the normals nearby, turned the way the normal a step into the parameters points. `None`
/// where the surface has no normal there: the first derivatives' cross product shorter than
/// `tolerance`, or a point the normals turn round (a cone's apex: the normals a step either
/// way along v opposite).
#[doc(alias = "GeomLib::NormEstim")]
pub fn estimated_normal(surface: &Surface, u: f64, v: f64, tolerance: f64) -> Option<Direction> {
    let tol2 = tolerance * tolerance;
    let (_, du, dv) = surface.d1(u, v);
    let (mdu, mdv) = (du.length_squared(), dv.length_squared());
    if mdu >= tol2 && mdv >= tol2 {
        let n = du.cross(dv);
        if n.length_squared() < tol2 {
            return None;
        }
        return Direction::try_from_vector(n);
    }
    let mut normal = limit_normal(&surface.d2(u, v), tolerance)?;
    let (umin, umax, vmin, vmax) = surface.bounds();
    let step = NORMAL_PARAMETER_STEP;
    let normal_at = |u: f64, v: f64| {
        let (_, a, b) = surface.d1(u, v);
        (a, b)
    };
    let long_enough = |a: Vector, b: Vector| a.length_squared() > NORMAL_SQUARED_MAGNITUDE_TOL && b.length_squared() > NORMAL_SQUARED_MAGNITUDE_TOL;
    let mut sign = -1.0;
    // A cone's apex: the normals a step either way along v opposite.
    if v > vmin + step && v < vmax - step {
        let (a, b) = normal_at(u, v - sign * step);
        if long_enough(a, b) {
            let before = Direction::try_from_vector(a.cross(b));
            let (a, b) = normal_at(u, v + sign * step);
            if long_enough(a, b) {
                let after = Direction::try_from_vector(a.cross(b));
                if let (Some(n1), Some(n2)) = (before, after) {
                    if n1.is_opposite(n2, CONE_NORMAL_ANGULAR_TOLERANCE) {
                        return None;
                    }
                }
            }
        }
    }
    // Along v.
    if mdu < tol2 && mdv >= tol2 {
        if vmax - v > v - vmin {
            sign = 1.0;
        }
        let (a, b) = normal_at(u, v + sign * step);
        let mut n = a.cross(b);
        if n.length_squared() < NORMAL_SQUARED_MAGNITUDE_TOL {
            let sign1 = if umax - u > u - umin { 1.0 } else { -1.0 };
            let (a, b) = normal_at(u + sign1 * step, v + sign * step);
            n = a.cross(b);
        }
        if n.length_squared() >= NORMAL_SQUARED_MAGNITUDE_TOL && n.dot(Vector::from(normal)) < 0.0 {
            normal = -normal;
        }
    }
    // Along u.
    if mdv < tol2 && mdu >= tol2 {
        if umax - u > u - umin {
            sign = 1.0;
        }
        let (a, b) = normal_at(u + sign * step, v);
        let mut n = a.cross(b);
        if n.length_squared() < NORMAL_SQUARED_MAGNITUDE_TOL {
            let sign1 = if vmax - v > v - vmin { 1.0 } else { -1.0 };
            let (a, b) = normal_at(u + sign * step, v + sign1 * step);
            n = a.cross(b);
        }
        if n.length_squared() >= NORMAL_SQUARED_MAGNITUDE_TOL && n.dot(Vector::from(normal)) < 0.0 {
            normal = -normal;
        }
    }
    Some(normal)
}

/// The principal curvatures of a surface at a point, with their directions, as
/// `GeomLProp_SLProps` gives them: signed along the surface's normal `D1U ^ D1V` (positive
/// where the surface bends towards it).
#[doc(alias = "GeomLProp_SLProps::MaxCurvature")]
#[doc(alias = "GeomLProp_SLProps::MinCurvature")]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SurfaceCurvatures {
    pub max: f64,
    pub min: f64,
    pub max_direction: Direction,
    pub min_direction: Direction,
    #[doc(alias = "GeomLProp_SLProps::MeanCurvature")]
    pub mean: f64,
    #[doc(alias = "GeomLProp_SLProps::GaussianCurvature")]
    pub gaussian: f64,
}

/// Local properties of a surface at one parameter pair, as `GeomLProp_SLProps` (to the
/// second derivatives) gives them. Quantities undefined there are `None`.
#[doc(alias = "GeomLProp_SLProps")]
#[doc(alias = "BRepLProp_SLProps")]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SurfacePoint {
    pub u: f64,
    pub v: f64,
    pub derivatives: super::SurfaceD2,
    /// The unit normal along `D1U ^ D1V`: none where a first derivative vanishes or the two
    /// are parallel (`IsNormalDefined`).
    pub normal: Option<Direction>,
    /// None where the normal, or a tangent along u or v (from the first of the first and
    /// second derivatives not vanishing), is undefined, or the shape operator has no two
    /// eigenvalues (`IsCurvatureDefined`).
    pub curvatures: Option<SurfaceCurvatures>,
}

/// The local properties of `surface` at `(u, v)`; derivatives (and the sine between the
/// first two) below `tolerance` count as vanishing (`GeomLProp_SLProps(S, u, v, 2, tol)`).
#[doc(alias = "GeomLProp_SLProps")]
#[doc(alias = "BRepLProp_SLProps")]
pub fn surface_point(surface: &Surface, u: f64, v: f64, tolerance: f64) -> SurfacePoint {
    let d = surface.d2(u, v);
    let normal = sl_normal(d.du, d.dv, tolerance);
    // (`IsTangentUDefined`, `IsTangentVDefined`: the first or the second derivative longer
    // than the tolerance.)
    let tol2 = tolerance * tolerance;
    let tangent_u = d.du.length_squared() > tol2 || d.duu.length_squared() > tol2;
    let tangent_v = d.dv.length_squared() > tol2 || d.dvv.length_squared() > tol2;
    let curvatures = match normal {
        Some(n) if tangent_u && tangent_v => sl_curvatures(&d, n),
        _ => None,
    };
    SurfacePoint { u, v, derivatives: d, normal, curvatures }
}

/// `CSLib::Normal(D1U, D1V, SinTol)`: none where a derivative vanishes or the sine between
/// them is below `sin_tol`.
#[doc(alias = "CSLib::Normal")]
fn sl_normal(du: Vector, dv: Vector, sin_tol: f64) -> Option<Direction> {
    let (mu, mv) = (du.length_squared(), dv.length_squared());
    if mu <= crate::tolerance::RESOLUTION || mv <= crate::tolerance::RESOLUTION {
        return None;
    }
    let n = du.cross(dv);
    if n.length_squared() / (mu * mv) < sin_tol * sin_tol {
        return None;
    }
    Direction::try_from_vector(n)
}

/// `LProp_SurfaceUtils::ComputeSurfCurvatures`: the shape operator's eigenvalues and
/// eigenvectors from the fundamental forms.
fn sl_curvatures(d: &super::SurfaceD2, normal: Direction) -> Option<SurfaceCurvatures> {
    let n = Vector::from(normal);
    let (e, f, g) = (d.du.length_squared(), d.du.dot(d.dv), d.dv.length_squared());
    let (l, m, nn) = (n.dot(d.duu), n.dot(d.duv), n.dot(d.dvv));
    let (a0, b0, c0) = (e * m - f * l, e * nn - g * l, f * nn - g * m);
    let largest = a0.abs().max(b0.abs()).max(c0.abs());
    let eps = f64::EPSILON;
    if largest < eps {
        // An umbilic point.
        if g < eps {
            return None;
        }
        let k = nn / g;
        return Some(SurfaceCurvatures {
            max: k,
            min: k,
            max_direction: Direction::try_from_vector(d.du.cross(n))?,
            min_direction: Direction::try_from_vector(d.du)?,
            mean: k,
            gaussian: k * k,
        });
    }
    let (a, b, c) = (a0 / largest, b0 / largest, c0 / largest);
    let two_roots = |a: f64, b: f64, c: f64| {
        let r = super::polynomial_roots::PolynomialRoots::quadratic(a, b, c);
        (r.done && !r.infinite && r.roots.len() == 2).then(|| (r.roots[0], r.roots[1]))
    };
    let (k1, k2, v1, v2) = if a.abs() > eps {
        let (r1, r2) = two_roots(a, b, c)?;
        let k = |r: f64| ((l * r + 2.0 * m) * r + nn) / ((e * r + 2.0 * f) * r + g);
        (k(r1), k(r2), d.du * r1 + d.dv, d.du * r2 + d.dv)
    } else if c.abs() > eps {
        let (r1, r2) = two_roots(c, b, a)?;
        let k = |r: f64| ((nn * r + 2.0 * m) * r + l) / ((g * r + 2.0 * f) * r + e);
        (k(r1), k(r2), d.du + d.dv * r1, d.du + d.dv * r2)
    } else {
        (l / e, nn / g, d.du, d.dv)
    };
    let ((min, min_v), (max, max_v)) = if k1 < k2 { ((k1, v1), (k2, v2)) } else { ((k2, v2), (k1, v1)) };
    let eg_ff = e * g - f * f;
    Some(SurfaceCurvatures {
        max,
        min,
        max_direction: Direction::try_from_vector(max_v)?,
        min_direction: Direction::try_from_vector(min_v)?,
        mean: (nn * e - 2.0 * m * f + l * g) / (2.0 * eg_ff),
        gaussian: (l * nn - m * m) / eg_ff,
    })
}
