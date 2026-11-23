//! Real roots of polynomials of degree up to four in closed form, and of the trigonometric
//! equation `A cos²x + 2B cos x sin x + C cos x + D sin x + E = 0` on an interval.
//!
//! Quadratics by the cancellation-free formula, cubics by Cardano's formula (one real root) or
//! Viète's trigonometric form (three), quartics by Ferrari's factorisation into two quadratics
//! through the largest root of its resolvent cubic. The variable is scaled by a power of two
//! so that the monic coefficients are of order one, discriminants within their own rounding
//! error count as zero (a double root, not none), and every root is polished by Newton steps on
//! the original polynomial that are kept only while they lower its value.

use super::adaptor::epsilon;
use std::f64::consts::{PI, TAU};

/// The real roots of a polynomial.
#[derive(Clone, Debug, Default)]
pub(crate) struct PolynomialRoots {
    /// False when the coefficients overflow the computation.
    pub done: bool,
    /// The polynomial is identically zero.
    pub infinite: bool,
    /// The real roots, a multiple root repeated by its multiplicity.
    pub roots: Vec<f64>,
}

/// A coefficient this small is no coefficient.
const NEGLIGIBLE: f64 = 1e-30;

/// The unit roundoff.
const U: f64 = f64::EPSILON / 2.0;

/// `p(x)` and `p'(x)` by Horner's rule, coefficients from the leading one.
fn horner(coefficients: &[f64], x: f64) -> (f64, f64) {
    let (mut p, mut d) = (0.0, 0.0);
    for &c in coefficients {
        d = d * x + p;
        p = p * x + c;
    }
    (p, d)
}

/// Newton steps on `coefficients` from `x`, kept while `|p|` drops.
fn polish(coefficients: &[f64], mut x: f64) -> f64 {
    if !x.is_finite() {
        return x;
    }
    let (mut p, mut d) = horner(coefficients, x);
    for _ in 0..6 {
        if p == 0.0 || d == 0.0 || !d.is_finite() {
            break;
        }
        let next = x - p / d;
        if next == x || !next.is_finite() {
            break;
        }
        let (pn, dn) = horner(coefficients, next);
        if pn.abs() >= p.abs() {
            break;
        }
        (x, p, d) = (next, pn, dn);
    }
    x
}

/// The power of two nearest to `v > 0`.
fn power_of_two(v: f64) -> f64 {
    if !(v > 0.0) || !v.is_finite() {
        return 1.0;
    }
    2f64.powi(v.log2().round().clamp(-1000.0, 1000.0) as i32)
}

/// The power of two scaling the roots of the monic polynomial `xⁿ + c₁xⁿ⁻¹ + … + cₙ` to order
/// one: near `max |c_i|^{1/i}`.
fn root_scale(monic: &[f64]) -> f64 {
    let size = monic.iter().enumerate().map(|(i, c)| c.abs().powf(1.0 / (i + 1) as f64)).fold(0.0, f64::max);
    power_of_two(size)
}

impl PolynomialRoots {
    fn failed() -> Self {
        PolynomialRoots { done: false, infinite: false, roots: Vec::new() }
    }

    fn found(roots: Vec<f64>) -> Self {
        PolynomialRoots { done: true, infinite: false, roots }
    }

    /// The roots of `a x + b`.
    pub(crate) fn linear(a: f64, b: f64) -> Self {
        if a.abs() <= NEGLIGIBLE {
            if b.abs() <= NEGLIGIBLE {
                return PolynomialRoots { done: true, infinite: true, roots: Vec::new() };
            }
            return Self::found(Vec::new());
        }
        let x = -b / a;
        if !x.is_finite() {
            return Self::failed();
        }
        Self::found(vec![x])
    }

    /// The roots of `a x² + b x + c`: two distinct ones larger in magnitude first, a double one
    /// twice.
    pub(crate) fn quadratic(a: f64, b: f64, c: f64) -> Self {
        if a.abs() <= NEGLIGIBLE {
            return Self::linear(b, c);
        }
        let (p, q) = (b / a, c / a);
        if !p.is_finite() || !q.is_finite() {
            return Self::failed();
        }
        let Some(roots) = monic_quadratic(p, q) else { return Self::failed() };
        let coefficients = [a, b, c];
        Self::found(roots.into_iter().map(|x| polish(&coefficients, x)).collect())
    }

    /// The roots of `a x³ + b x² + c x + d` (see the module notes for their order).
    pub(crate) fn cubic(a: f64, b: f64, c: f64, d: f64) -> Self {
        if a.abs() <= NEGLIGIBLE {
            return Self::quadratic(b, c, d);
        }
        let monic = [b / a, c / a, d / a];
        if monic.iter().any(|v| !v.is_finite()) {
            return Self::failed();
        }
        let Some(roots) = monic_cubic(monic[0], monic[1], monic[2]) else { return Self::failed() };
        let coefficients = [a, b, c, d];
        Self::found(roots.into_iter().map(|x| polish(&coefficients, x)).collect())
    }

    /// The roots of `a x⁴ + b x³ + c x² + d x + e`: those of the first of Ferrari's quadratic
    /// factors, then the second's. A leading coefficient within a few hundred ulps of the
    /// largest other one (and not the scale of the whole polynomial) is dropped.
    pub(crate) fn quartic(a: f64, b: f64, c: f64, d: f64, e: f64) -> Self {
        if a.abs() <= NEGLIGIBLE {
            return Self::cubic(b, c, d, e);
        }
        let others = [b, c, d, e];
        let largest = others.iter().fold(0.0f64, |m, v| m.max(v.abs()));
        let tiny = 256.0 * f64::EPSILON * largest;
        if a.abs() <= tiny && !others.iter().any(|&v| v != 0.0 && v.abs() <= tiny) {
            return Self::cubic(b, c, d, e);
        }
        let monic = [b / a, c / a, d / a, e / a];
        if monic.iter().any(|v| !v.is_finite()) {
            return Self::failed();
        }
        let s = root_scale(&monic);
        let (a1, b1, c1, d1) = (monic[0] / s, monic[1] / (s * s), monic[2] / (s * s * s), monic[3] / (s * s * s * s));
        let Some(factors) = ferrari(a1, b1, c1, d1) else { return Self::failed() };
        let coefficients = [a, b, c, d, e];
        let mut roots = Vec::with_capacity(4);
        for (p, q) in factors {
            let Some(pair) = monic_quadratic(p, q) else { return Self::failed() };
            roots.extend(pair.into_iter().map(|y| polish(&coefficients, y * s)));
        }
        Self::found(roots)
    }
}

/// The roots of `x² + p x + q` in the order of [`PolynomialRoots::quadratic`] (unpolished);
/// `None` on overflow.
fn monic_quadratic(p: f64, q: f64) -> Option<Vec<f64>> {
    // Scaled so that p and q are of order one.
    let s = root_scale(&[p, q]);
    let (ps, qs) = (p / s, q / (s * s));
    let half = 0.5 * ps;
    let disc = half * half - qs;
    // (The rounding error of the discriminant.)
    let error = 4.0 * U * (half * half + qs.abs());
    if disc.abs() <= error {
        let x = -half * s;
        return x.is_finite().then(|| vec![x, x]);
    }
    if disc < 0.0 {
        return Some(Vec::new());
    }
    let root = disc.sqrt();
    // The root larger in magnitude first (−p/2 ∓ √disc with the sign of p), then by the product.
    let first = if ps > 0.0 { -half - root } else { -half + root };
    let second = qs / first;
    let (x1, x2) = (first * s, second * s);
    (x1.is_finite() && x2.is_finite()).then(|| vec![x1, x2])
}

/// The roots of `x³ + β x² + γ x + δ` (unpolished): three distinct ones the extreme one
/// opposite to the sign of `β` first, then the middle one (but `√−γ, −√−γ, 0` for
/// `x³ + γx`); a double root twice before the simple one; a single or triple root.
fn monic_cubic(beta: f64, gamma: f64, delta: f64) -> Option<Vec<f64>> {
    if beta == 0.0 && delta == 0.0 {
        if gamma < 0.0 {
            let r = (-gamma).sqrt();
            return Some(vec![r, -r, 0.0]);
        }
        if gamma == 0.0 {
            return Some(vec![0.0; 3]);
        }
    }
    let s = root_scale(&[beta, gamma, delta]);
    let (b, c, d) = (beta / s, gamma / (s * s), delta / (s * s * s));
    // Depressed: x = y − b/3, y³ + P y + Q.
    let shift = b / 3.0;
    let p = c - b * shift;
    let q = (2.0 * shift * shift - c) * shift + d;
    let p_error = 2.0 * U * (c.abs() + (b * shift).abs());
    let q_error = 4.0 * U * ((2.0 * shift * shift * shift).abs() + (c * shift).abs() + d.abs());
    let (h, t) = (0.5 * q, p / 3.0);
    let disc = h * h + t * t * t;
    let disc_error = 4.0 * (h.abs() * q_error + t * t * p_error + U * (h * h + (t * t * t).abs())) ;
    let out = |ys: &[f64]| -> Option<Vec<f64>> {
        let xs: Vec<f64> = ys.iter().map(|y| (y - shift) * s).collect();
        xs.iter().all(|x| x.is_finite()).then_some(xs)
    };
    if disc.abs() <= disc_error {
        // A double root (or triple, P and Q both vanishing).
        if p.abs() <= p_error && q.abs() <= q_error {
            return out(&[0.0, 0.0, 0.0]);
        }
        let a = (-h).cbrt();
        return out(&[-a, -a, 2.0 * a]);
    }
    if disc > 0.0 {
        // One real root (Cardano, the cube root of the larger term).
        let a = -(h.abs() + disc.sqrt()).cbrt().copysign(h);
        let y = if a != 0.0 { a - t / a } else { 0.0 };
        return out(&[y]);
    }
    // Three real roots (Viète): y_k = 2√(−t) cos(θ − 2πk/3), largest, middle, smallest.
    let r = (-t).sqrt();
    let cos3 = (-h / (r * r * r)).clamp(-1.0, 1.0);
    let theta = cos3.acos() / 3.0;
    let largest = 2.0 * r * theta.cos();
    let middle = 2.0 * r * (theta - TAU / 3.0).cos();
    let smallest = 2.0 * r * (theta - 2.0 * TAU / 3.0).cos();
    if beta >= 0.0 { out(&[smallest, middle, largest]) } else { out(&[largest, middle, smallest]) }
}

/// The two monic quadratic factors `(p, q)` of `x⁴ + a x³ + b x² + c x + d` (coefficients of
/// order one), the first with `+P`, from the largest root of the resolvent cubic.
fn ferrari(a: f64, b: f64, c: f64, d: f64) -> Option<[(f64, f64); 2]> {
    let resolvent = [1.0, -b, a * c - 4.0 * d, d * (4.0 * b - a * a) - c * c];
    let ys = monic_cubic(resolvent[1], resolvent[2], resolvent[3])?;
    let y = polish(&resolvent, ys.iter().copied().fold(f64::NEG_INFINITY, f64::max));
    let p2 = (0.25 * a * a - b + y).max(0.0);
    let q2 = (0.25 * y * y - d).max(0.0);
    // 2PQ' = a y/2 − c: the smaller of P, Q' from the larger, which is better known.
    let cross = 0.5 * a * y - c;
    let sign = if cross >= 0.0 { 1.0 } else { -1.0 };
    let (pp, qq) = if p2 >= q2 {
        let pp = p2.sqrt();
        (pp, if pp > 0.0 { cross / (2.0 * pp) } else { sign * q2.sqrt() })
    } else {
        let qq = sign * q2.sqrt();
        (if qq != 0.0 { cross / (2.0 * qq) } else { p2.sqrt() }, qq)
    };
    let factors = [(0.5 * a + pp, 0.5 * y + qq), (0.5 * a - pp, 0.5 * y - qq)];
    factors.iter().all(|(p, q)| p.is_finite() && q.is_finite()).then_some(factors)
}

/// The roots of the trigonometric equation.
#[derive(Clone, Debug, Default)]
pub(crate) struct TrigonometricRoots {
    pub done: bool,
    /// Every x is a root (all coefficients negligible).
    pub infinite: bool,
    pub roots: Vec<f64>,
}

/// A trigonometric coefficient this small is none.
const TRIG_NEGLIGIBLE: f64 = 1e-12;

/// The roots of `A cos²x + 2B cos x sin x + C cos x + D sin x + E = 0` in the window
/// `[inf, inf + min(sup − inf, 2π)]` (`[0, 2π]` with both bounds infinite, the period before
/// `sup` with only `inf` infinite).
pub(crate) fn trigonometric_roots(a: f64, b: f64, c: f64, d: f64, e: f64, inf_bound: f64, sup_bound: f64) -> TrigonometricRoots {
    let negligible = |v: f64| v.abs() <= TRIG_NEGLIGIBLE;
    let (w0, w1) = match (inf_bound <= -f64::MAX, sup_bound >= f64::MAX) {
        (true, true) => (0.0, TAU),
        (false, true) => (inf_bound, inf_bound + TAU),
        (true, false) => (sup_bound - TAU, sup_bound),
        (false, false) => (inf_bound, inf_bound + (sup_bound - inf_bound).min(TAU)),
    };
    let width = w1 - w0;
    let tol = epsilon(width);
    let mut out = TrigonometricRoots { done: true, infinite: false, roots: Vec::new() };
    if [a, b, c, d, e].iter().all(|&v| negligible(v)) {
        out.infinite = true;
        return out;
    }
    // A candidate shifted by whole turns to the first place at or after the window's start.
    let shifted = |x: f64, slack: f64| -> Option<f64> {
        if x >= w0 - slack && x <= w1 + slack {
            return Some(x);
        }
        let k = ((w0 - slack - x) / TAU).ceil();
        let y = x + k * TAU;
        (y >= w0 - slack && y <= w1 + slack).then_some(y)
    };
    let f = |x: f64| {
        let (s, co) = x.sin_cos();
        a * co * co + 2.0 * b * co * s + c * co + d * s + e
    };
    let df = |x: f64| {
        let (s, co) = x.sin_cos();
        -2.0 * a * co * s + 2.0 * b * (co * co - s * s) - c * s + d * co
    };
    if negligible(a) && negligible(b) {
        if negligible(c) && negligible(d) {
            return out;
        }
        if negligible(c) {
            let s = -e / d;
            if s.abs() > 1.0 {
                return out;
            }
            let x = s.asin();
            let first = if x < 0.0 { TAU + x } else { x };
            out.roots.extend([first, PI - x].into_iter().filter_map(|x| shifted(x, tol)));
            return out;
        }
        if negligible(d) {
            let s = -e / c;
            if s.abs() > 1.0 {
                return out;
            }
            let x = s.acos();
            out.roots.extend([x, TAU - x].into_iter().filter_map(|x| shifted(x, tol)));
            return out;
        }
        let q = PolynomialRoots::quadratic(e - c, 2.0 * d, e + c);
        return generic(q, negligible(e - c), w0, w1, tol, &f, &df, &shifted);
    }
    if negligible(a) && negligible(e) {
        let clamp = |x: f64| x.clamp(inf_bound, sup_bound);
        let slack = crate::tolerance::PARAMETRIC;
        if negligible(c) {
            let mut xs = vec![0.0, PI];
            let ratio = d / (2.0 * b);
            if ratio.abs() <= 1.0 + 1e-9 {
                let x = (-ratio).clamp(-1.0, 1.0).acos();
                xs.extend([x, TAU - x]);
            }
            out.roots.extend(xs.into_iter().filter_map(|x| shifted(x, slack)).map(clamp));
            return out;
        }
        if negligible(d) {
            let mut xs = vec![PI / 2.0, 1.5 * PI];
            let ratio = c / (2.0 * b);
            if ratio.abs() <= 1.0 + 1e-9 {
                let x = (-ratio).clamp(-1.0, 1.0).asin();
                xs.extend([x, PI - x]);
            }
            out.roots.extend(xs.into_iter().filter_map(|x| shifted(x, slack)).map(clamp));
            return out;
        }
    }
    // t = tan(x/2): cos x = (1 − t²)/(1 + t²), sin x = 2t/(1 + t²).
    let q = PolynomialRoots::quartic(a - c + e, 2.0 * d - 4.0 * b, 2.0 * (e - a), 4.0 * b + 2.0 * d, a + c + e);
    generic(q, negligible(a - c + e), w0, w1, tol, &f, &df, &shifted)
}

/// The roots `x = 2 atan t` of the half-angle polynomial's roots `t`, with `x = π` where the
/// leading coefficient vanishes (`t = ∞`), refined, in the window, ascending.
#[allow(clippy::too_many_arguments)]
fn generic(q: PolynomialRoots, at_pi: bool, w0: f64, w1: f64, tol: f64, f: &dyn Fn(f64) -> f64, df: &dyn Fn(f64) -> f64, shifted: &dyn Fn(f64, f64) -> Option<f64>) -> TrigonometricRoots {
    let mut out = TrigonometricRoots { done: q.done, infinite: false, roots: Vec::new() };
    if !q.done {
        return out;
    }
    let mut xs: Vec<f64> = q.roots.iter().map(|t| 2.0 * t.atan()).collect();
    if at_pi {
        xs.push(PI);
    }
    let width = w1 - w0;
    let from_polynomial = xs.len() - usize::from(at_pi);
    for (i, x) in xs.into_iter().enumerate() {
        // Newton on f, kept while |f| drops and the root stays put.
        let (mut x, mut v) = (x, f(x));
        for _ in 0..6 {
            let dv = df(x);
            if v == 0.0 || dv == 0.0 {
                break;
            }
            let next = x - v / dv;
            let vn = f(next);
            if vn.abs() >= v.abs() || (next - x).abs() > 1e-3 * width.max(1e-3) {
                break;
            }
            (x, v) = (next, vn);
        }
        let Some(y) = shifted(x, tol) else { continue };
        // (π once, though a polynomial root lies there too.)
        if i >= from_polynomial && out.roots.iter().any(|&r| (r - y).abs() <= 1e-9 * (1.0 + y.abs())) {
            continue;
        }
        out.roots.push(y);
    }
    out.roots.sort_by(f64::total_cmp);
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn assert_roots(got: &[f64], expected: &[f64], tol: f64) {
        assert_eq!(got.len(), expected.len(), "{got:?} vs {expected:?}");
        for (g, e) in got.iter().zip(expected) {
            assert!((g - e).abs() <= tol * (1.0 + e.abs()), "{got:?} vs {expected:?}");
        }
    }

    #[test]
    fn orders_of_the_closed_forms() {
        assert_roots(&PolynomialRoots::quartic(1.0, -10.0, 35.0, -50.0, 24.0).roots, &[2.0, 1.0, 4.0, 3.0], 1e-14);
        assert_eq!(PolynomialRoots::quartic(1.0, 0.0, -2.0, 0.0, 1.0).roots, [-1.0, -1.0, 1.0, 1.0]);
        assert_roots(&PolynomialRoots::cubic(1.0, -6.0, 11.0, -6.0).roots, &[3.0, 2.0, 1.0], 1e-15);
        assert_roots(&PolynomialRoots::cubic(1.0, 6.0, 11.0, 6.0).roots, &[-3.0, -2.0, -1.0], 1e-15);
        assert_roots(&PolynomialRoots::cubic(2.0, -3.0, 0.0, 1.0).roots, &[1.0, 1.0, -0.5], 1e-15);
        assert_eq!(PolynomialRoots::cubic(1.0, 0.0, -4.0, 0.0).roots, [2.0, -2.0, 0.0]);
        assert_eq!(PolynomialRoots::cubic(1.0, 0.0, 0.0, -8.0).roots, [2.0]);
        assert_roots(&PolynomialRoots::quadratic(1.0, -3.0, 2.0).roots, &[2.0, 1.0], 0.0);
        assert_roots(&PolynomialRoots::quadratic(1.0, 3.0, 2.0).roots, &[-2.0, -1.0], 0.0);
        assert_eq!(PolynomialRoots::quadratic(1.0, 2.0, 1.0).roots, [-1.0, -1.0]);
        assert!(PolynomialRoots::quadratic(1.0, 0.0, 1.0).roots.is_empty());
        assert_roots(&PolynomialRoots::cubic(1.0, -3.0, 3.0, -1.0).roots, &[1.0, 1.0, 1.0], 1e-15);
    }

    #[test]
    fn well_separated_and_clustered_roots() {
        let quartic = |r: [f64; 4]| {
            let (s1, s2, s3, s4) = (
                r[0] + r[1] + r[2] + r[3],
                r[0] * r[1] + r[0] * r[2] + r[0] * r[3] + r[1] * r[2] + r[1] * r[3] + r[2] * r[3],
                r[0] * r[1] * r[2] + r[0] * r[1] * r[3] + r[0] * r[2] * r[3] + r[1] * r[2] * r[3],
                r[0] * r[1] * r[2] * r[3],
            );
            PolynomialRoots::quartic(1.0, -s1, s2, -s3, s4)
        };
        for r in [[-3.7, 0.25, 1.9, 12.5], [-0.01, 0.3, 0.31, -50.0], [1e-3, 2e-3, 7.0, -7.0]] {
            let mut got = quartic(r).roots;
            got.sort_by(f64::total_cmp);
            let mut e = r.to_vec();
            e.sort_by(f64::total_cmp);
            for (g, x) in got.iter().zip(&e) {
                assert!((g - x).abs() <= 1e-13 * x.abs().max(1.0), "{got:?} {e:?}");
            }
        }
        let mut got = quartic([1.0, 1.001, 1.002, 1.003]).roots;
        got.sort_by(f64::total_cmp);
        // (Conditioning: coefficient rounding ~1e-15 over |p′| ~6e-9 moves the roots ~2e-7.)
        assert_roots(&got, &[1.0, 1.001, 1.002, 1.003], 1e-6);
    }

    #[test]
    fn scaling_and_reduction() {
        let mut r = PolynomialRoots::quartic(1e-20, 0.0, 0.0, 0.0, -1e-20).roots;
        r.sort_by(f64::total_cmp);
        assert_roots(&r, &[-1.0, 1.0], 1e-15);
        let small = PolynomialRoots::quadratic(1.0, -1e8, 1.0).roots;
        assert_eq!(small.len(), 2);
        assert!((small[1] - 1e-8).abs() <= 1e-23, "{small:?}");
        assert_roots(&PolynomialRoots::quartic(0.0, 1.0, -6.0, 11.0, -6.0).roots, &[3.0, 2.0, 1.0], 1e-15);
        let l = PolynomialRoots::linear(0.0, 0.0);
        assert!(l.done && l.infinite && l.roots.is_empty());
        assert!(PolynomialRoots::linear(0.0, 1.0).roots.is_empty());
        assert_eq!(PolynomialRoots::linear(2.0, -1.0).roots, [0.5]);
        // Tangency survives rounding.
        assert_eq!(PolynomialRoots::quadratic(1.0, -2.0, 1.0 + 1e-17).roots.len(), 2);
        assert_eq!(PolynomialRoots::quadratic(1.0, -0.2, 0.01).roots.len(), 2);
        assert!(PolynomialRoots::quadratic(1.0, -2.0, 1.0 + 1e-12).roots.is_empty());
        // Huge roots of a quartic whose leading coefficient is rounding noise: the cubic's.
        assert_eq!(PolynomialRoots::quartic(1e-17, 1.0, -6.0, 11.0, -6.0).roots.len(), 3);
    }

    #[test]
    fn trigonometric_forms() {
        let r = trigonometric_roots(0.0, 0.0, 0.0, 2.0, 1.0, 0.0, TAU).roots;
        assert_roots(&r, &[11.0 * PI / 6.0, 7.0 * PI / 6.0], 1e-14);
        let r = trigonometric_roots(0.0, 0.5, 0.0, 0.3, 0.0, 0.0, TAU).roots;
        assert_roots(&r, &[0.0, PI, (-0.3f64).acos(), TAU - (-0.3f64).acos()], 1e-14);
        let r = trigonometric_roots(0.0, 0.0, 1.0, 0.0, -1.0, 0.0, TAU).roots;
        assert_roots(&r, &[0.0, TAU], 1e-14);
        // cos²x + cos x: π/2, π, 3π/2 (π from the vanishing leading coefficient).
        let r = trigonometric_roots(1.0, 0.0, 1.0, 0.0, 0.0, 0.0, TAU).roots;
        assert_roots(&r, &[PI / 2.0, PI, 1.5 * PI], 1e-14);
        // A window not starting in [0, 2π).
        let r = trigonometric_roots(0.0, 0.0, 0.0, 2.0, 1.0, -1.0, 3.0).roots;
        assert_roots(&r, &[-PI / 6.0], 1e-14);
        let r = trigonometric_roots(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, TAU);
        assert!(r.infinite && r.roots.is_empty());
        // Generic: cos²x + 2·0.3 cos x sin x + 0.1 cos x − 0.4 sin x − 0.2.
        let (a, b, c, d, e) = (1.0, 0.3, 0.1, -0.4, -0.2);
        let r = trigonometric_roots(a, b, c, d, e, 0.0, TAU).roots;
        assert!(!r.is_empty());
        assert!(r.windows(2).all(|w| w[0] <= w[1]));
        for x in &r {
            let (s, co) = x.sin_cos();
            assert!((a * co * co + 2.0 * b * co * s + c * co + d * s + e).abs() < 1e-14);
        }
        // Every root found: by sampling.
        let f = |x: f64| {
            let (s, co) = x.sin_cos();
            a * co * co + 2.0 * b * co * s + c * co + d * s + e
        };
        let changes = (0..1000).filter(|&i| f(TAU * i as f64 / 1000.0) * f(TAU * (i + 1) as f64 / 1000.0) < 0.0).count();
        assert_eq!(r.len(), changes);
    }
}
