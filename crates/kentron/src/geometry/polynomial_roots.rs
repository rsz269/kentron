//! OCCT's closed-form root solvers, ported line by line: `math_DirectPolynomialRoots` (the
//! real roots of polynomials up to degree 4: Ferrari, Cardano, the stable quadratic formula,
//! each root polished by Newton steps) and `math_TrigonometricFunctionRoots` (the roots of
//! `a cos²x + 2b cos x sin x + c cos x + d sin x + e` over an interval, through the tangent
//! of the half angle, polished by `math_NewtonFunctionRoot`, ported here too).

use std::f64::consts::PI;

use super::adaptor::epsilon;
use super::function_roots::FunctionWithDerivative;

/// `ZERO_THRESHOLD`: coefficients at most this count as zero.
const ZERO_THRESHOLD: f64 = 1.0e-30;
/// `MACHINE_EPSILON` (`RealEpsilon()`).
const MACHINE_EPSILON: f64 = f64::EPSILON;
/// `FLOATING_RADIX`.
const FLOATING_RADIX: f64 = 2.0;
/// `MAX_NEWTON_ITERATIONS`.
const MAX_NEWTON_ITERATIONS: i32 = 10;
/// `OVERFLOW_LIMIT`.
const OVERFLOW_LIMIT: f64 = 1.0e+80;
/// `RealSmall()`.
const REAL_SMALL: f64 = f64::MIN_POSITIVE;
/// `Precision::PConfusion()`.
const P_CONFUSION: f64 = 1.0e-9;
/// `RealLast()`.
const REAL_LAST: f64 = f64::MAX;
/// `RealFirst()`.
const REAL_FIRST: f64 = -f64::MAX;

/// `INV_LOG_RADIX`: 1 / ln 2, as OCCT computes it.
fn inv_log_radix() -> f64 {
    1.0 / 2.0f64.ln()
}

/// The outcome of `math_DirectPolynomialRoots`.
#[doc(alias = "math_DirectPolynomialRoots")]
#[derive(Clone, Debug, Default)]
pub(crate) struct PolynomialRoots {
    /// `IsDone()`.
    pub done: bool,
    /// `InfiniteRoots()`: the polynomial is zero.
    pub infinite: bool,
    /// `Value(1..=NbSolutions())`, in OCCT's order (double roots repeated).
    pub roots: Vec<f64>,
}

/// `EvaluatePolynomial`: Horner's scheme, highest coefficient first.
fn evaluate_polynomial(poly: &[f64], x: f64) -> f64 {
    let mut result = poly[0];
    for &c in &poly[1..] {
        result = result * x + c;
    }
    result
}

/// `EvaluatePolynomialWithDerivative`.
fn evaluate_polynomial_with_derivative(poly: &[f64], x: f64) -> (f64, f64) {
    let mut value = poly[0] * x + poly[1];
    let mut derivative = poly[0];
    for &c in &poly[2..] {
        derivative = derivative * x + value;
        value = value * x + c;
    }
    (value, derivative)
}

/// `RefineRoot`: Newton steps from `guess`, kept only if the last value computed is no
/// larger than the value at the guess.
fn refine_root(poly: &[f64], guess: f64) -> f64 {
    let mut value = 0.0;
    let mut solution = guess;
    let initial_value = evaluate_polynomial(poly, guess);
    for _ in 1..MAX_NEWTON_ITERATIONS {
        let derivative;
        (value, derivative) = evaluate_polynomial_with_derivative(poly, solution);
        if derivative.abs() <= ZERO_THRESHOLD {
            break;
        }
        let delta = -value / derivative;
        if delta.abs() <= MACHINE_EPSILON * solution.abs() {
            break;
        }
        solution += delta;
    }
    // Return the improved solution only if it is better.
    if value.abs() <= initial_value.abs() { solution } else { guess }
}

/// `ScaleAndRefineAllRoots`.
fn scale_and_refine_all_roots(roots: &mut [f64], scale_factor: f64, poly: &[f64]) {
    for root in roots.iter_mut() {
        *root *= scale_factor;
        *root = refine_root(poly, *root);
    }
}

/// `ComputeBaseExponent`: the base-2 exponent used to scale the coefficients.
fn compute_base_exponent(value: f64) -> i32 {
    if value > 1.0 {
        (value.ln() * inv_log_radix()) as i32
    } else if value < -1.0 {
        (-(-value).ln() * inv_log_radix()) as i32
    } else {
        0
    }
}

/// `ComputeSpecialDiscriminant`: the cubic's discriminant for `P < 0`, computed stably.
fn compute_special_discriminant(beta: f64, gamma: f64, del: f64, a1: f64) -> f64 {
    let sigma = beta * gamma / 3.0 - 2.0 * beta * beta * beta / 27.0;
    let psi = gamma * gamma * (4.0 * gamma - beta * beta) / 27.0;
    let d1 = if sigma >= 0.0 { sigma + 2.0 * (-a1).sqrt() } else { sigma - 2.0 * (-a1).sqrt() };
    let d2 = psi / d1;
    if (del - d1).abs() >= 18.0 * MACHINE_EPSILON * (del.abs() + d1.abs()) && (del - d2).abs() >= 24.0 * MACHINE_EPSILON * (del.abs() + d2.abs()) {
        return (del - d1) * (del - d2) / 4.0;
    }
    0.0
}

/// `SolveCubicThreeRealRoots`: the trigonometric method.
fn solve_cubic_three_real_roots(beta: f64, gamma: f64, del: f64, p: f64, q: f64, discr: f64, roots: &mut [f64; 4]) {
    if beta == 0.0 && q == 0.0 {
        // The special case x^3 + Px = 0.
        roots[0] = (-p).sqrt();
        roots[1] = -roots[0];
        roots[2] = 0.0;
    } else {
        let sb = if beta >= 0.0 { 1.0 } else { -1.0 };
        let omega = (0.5 * q / (-discr).sqrt()).atan();
        let sp3 = (-p / 3.0).sqrt();
        let y1 = -2.0 * sb * sp3 * (PI / 6.0 - sb * omega / 3.0).cos();
        roots[0] = -beta / 3.0 + y1;
        if beta * q <= 0.0 {
            roots[1] = -beta / 3.0 + 2.0 * sp3 * (omega / 3.0).sin();
        } else {
            // An alternative formula, more accurate.
            let dbg = del - beta * gamma;
            let sdbg = if dbg >= 0.0 { 1.0 } else { -1.0 };
            let den1 = 8.0 * beta * beta / 9.0 - 4.0 * beta * y1 / 3.0 - 2.0 * q / y1;
            let den2 = 2.0 * y1 * y1 - q / y1;
            roots[1] = dbg / den1 + sdbg * (-27.0 * discr).sqrt() / den2;
        }
        // Vieta's formula for the third root.
        roots[2] = -del / (roots[0] * roots[1]);
    }
}

/// `SolveCubicOneRealRoot`: Cardano's formula.
fn solve_cubic_one_real_root(beta: f64, del: f64, p: f64, q: f64, discr: f64, roots: &mut [f64; 4]) {
    let mut u = discr.sqrt() + (q / 2.0).abs();
    u = if u >= 0.0 { u.powf(1.0 / 3.0) } else { -(u.abs().powf(1.0 / 3.0)) };
    let h = if p >= 0.0 { u * u + p / 3.0 + (p / u) * (p / u) / 9.0 } else { u * q.abs() / (u * u - p / 3.0) };
    if beta * q >= 0.0 {
        if h.abs() <= REAL_SMALL && q.abs() <= REAL_SMALL {
            roots[0] = -beta / 3.0 - u + p / (3.0 * u);
        } else {
            roots[0] = -beta / 3.0 - q / h;
        }
    } else {
        roots[0] = -del / (beta * beta / 9.0 + h - beta * q / (3.0 * h));
    }
}

/// `SolveCubicMultipleRoots`: a zero discriminant.
fn solve_cubic_multiple_roots(beta: f64, gamma: f64, del: f64, p: f64, q: f64, roots: &mut [f64; 4]) -> usize {
    let sq = if q >= 0.0 { 1.0 } else { -1.0 };
    let sp3 = (-p / 3.0).sqrt();
    if beta * q <= 0.0 {
        roots[0] = -beta / 3.0 + sq * sp3;
        roots[1] = roots[0];
        if beta * q == 0.0 {
            roots[2] = -beta / 3.0 - 2.0 * sq * sp3;
        } else {
            roots[2] = -del / (roots[0] * roots[1]);
        }
    } else {
        roots[0] = -gamma / (beta + 3.0 * sq * sp3);
        roots[1] = roots[0];
        roots[2] = -beta / 3.0 - 2.0 * sq * sp3;
    }
    3
}

/// `ShouldReduceDegreeQuartic`: whether the leading coefficient counts as zero against the
/// others (jgv, 22.01.09, for numerical stability).
fn should_reduce_degree_quartic(a: f64, b: f64, c: f64, d: f64, e: f64) -> bool {
    if a.abs() <= ZERO_THRESHOLD {
        return true;
    }
    let mut max_coeff = ZERO_THRESHOLD;
    max_coeff = max_coeff.max(b.abs());
    max_coeff = max_coeff.max(c.abs());
    max_coeff = max_coeff.max(d.abs());
    max_coeff = max_coeff.max(e.abs());
    if max_coeff > ZERO_THRESHOLD {
        max_coeff = epsilon(100.0 * max_coeff);
    }
    if a.abs() <= max_coeff {
        let max_coeff1000 = 1000.0 * max_coeff;
        let mut with_a = false;
        for x in [b, c, d, e] {
            if x.abs() > ZERO_THRESHOLD && x.abs() <= max_coeff1000 {
                with_a = true;
            }
        }
        return !with_a;
    }
    false
}

impl PolynomialRoots {
    fn new() -> Self {
        PolynomialRoots { done: true, infinite: false, roots: Vec::new() }
    }

    /// `math_DirectPolynomialRoots(A, B, C, D, E)`: `A x⁴ + B x³ + C x² + D x + E`.
    pub(crate) fn quartic(a: f64, b: f64, c: f64, d: f64, e: f64) -> Self {
        let mut out = Self::new();
        out.solve_quartic(a, b, c, d, e);
        out
    }

    /// `math_DirectPolynomialRoots(A, B, C, D)`: `A x³ + B x² + C x + D`.
    pub(crate) fn cubic(a: f64, b: f64, c: f64, d: f64) -> Self {
        let mut out = Self::new();
        out.solve_cubic(a, b, c, d);
        out
    }

    /// `math_DirectPolynomialRoots(A, B, C)`: `A x² + B x + C`.
    pub(crate) fn quadratic(a: f64, b: f64, c: f64) -> Self {
        let mut out = Self::new();
        out.solve_quadratic(a, b, c);
        out
    }

    /// `math_DirectPolynomialRoots(A, B)`: `A x + B`.
    pub(crate) fn linear(a: f64, b: f64) -> Self {
        let mut out = Self::new();
        out.solve_linear(a, b);
        out
    }

    fn solve_quartic(&mut self, a0: f64, b0: f64, c0: f64, d0: f64, e0: f64) {
        // Check for degree reduction.
        if should_reduce_degree_quartic(a0, b0, c0, d0, e0) {
            self.solve_cubic(b0, c0, d0, e0);
            return;
        }
        // Normalize the coefficients.
        let a = b0 / a0;
        let b = c0 / a0;
        let c = d0 / a0;
        let d = e0 / a0;
        // Scale them against overflow and underflow (`ScaleQuartic(A, B, C, D, D)`).
        let exp = compute_base_exponent(d) / 4;
        let scale_factor = FLOATING_RADIX.powi(exp);
        let scale_factor2 = scale_factor * scale_factor;
        let sa = a / scale_factor;
        let sb = b / scale_factor2;
        let sc = c / (scale_factor2 * scale_factor);
        let sd = d / (scale_factor2 * scale_factor2);
        // Solve Ferrari's resolvent cubic: Y^3 + R3 Y^2 + S3 Y + T3 = 0.
        let r3 = -sb;
        let s3 = sa * sc - 4.0 * sd;
        let t3 = sd * (4.0 * sb - sa * sa) - sc * sc;
        let cubic = PolynomialRoots::cubic(1.0, r3, s3, t3);
        if !cubic.done {
            self.done = false;
            return;
        }
        // The largest root, for numerical stability.
        let mut y0 = cubic.roots[0];
        for &y in &cubic.roots[1..] {
            if y > y0 {
                y0 = y;
            }
        }
        // Factor into two quadratics (`FactorQuarticViaFerrari`).
        let discr = sa * y0 * 0.5 - sc;
        let sdiscr = if discr >= 0.0 { 1.0 } else { -1.0 };
        let mut p0 = sa * sa * 0.25 - sb + y0;
        p0 = if p0 < 0.0 { 0.0 } else { p0.sqrt() };
        let mut q0 = y0 * y0 * 0.25 - sd;
        // Q0 very close to zero is taken as zero.
        if q0.abs() < 10.0 * MACHINE_EPSILON {
            q0 = 0.0;
        } else {
            q0 = if q0 < 0.0 { 0.0 } else { q0.sqrt() };
        }
        let ademi = sa * 0.5;
        let ydemi = y0 * 0.5;
        let sdiscr_q0 = sdiscr * q0;
        let mut p1 = ademi + p0;
        let mut q1 = ydemi + sdiscr_q0;
        let mut p2 = ademi - p0;
        let mut q2 = ydemi - sdiscr_q0;
        // Clean up the near-zero coefficients.
        let eps = 100.0 * MACHINE_EPSILON;
        if p1.abs() <= eps {
            p1 = 0.0;
        }
        if p2.abs() <= eps {
            p2 = 0.0;
        }
        if q1.abs() <= eps {
            q1 = 0.0;
        }
        if q2.abs() <= eps {
            q2 = 0.0;
        }
        let quadratic1 = PolynomialRoots::quadratic(1.0, p1, q1);
        if !quadratic1.done {
            self.done = false;
            return;
        }
        let quadratic2 = PolynomialRoots::quadratic(1.0, p2, q2);
        if !quadratic2.done {
            self.done = false;
            return;
        }
        // Collect all roots, then scale back and refine them.
        self.roots = quadratic1.roots.iter().chain(quadratic2.roots.iter()).copied().collect();
        scale_and_refine_all_roots(&mut self.roots, scale_factor, &[a0, b0, c0, d0, e0]);
    }

    fn solve_cubic(&mut self, a0: f64, b0: f64, c0: f64, d0: f64) {
        // Check for degree reduction.
        if a0.abs() <= ZERO_THRESHOLD {
            self.solve_quadratic(b0, c0, d0);
            return;
        }
        // Normalize the coefficients.
        let beta = b0 / a0;
        let gamma = c0 / a0;
        let del = d0 / a0;
        // Scale against overflow and underflow (`ScaleCubic(Beta, Gamma, Del, Del)`).
        let exp = compute_base_exponent(del) / 3;
        let scale_factor = FLOATING_RADIX.powi(exp);
        let scale_factor2 = scale_factor * scale_factor;
        let sa = beta / scale_factor;
        let sb = gamma / scale_factor2;
        let sc = del / (scale_factor2 * scale_factor);
        // To the depressed cubic t^3 + Pt + Q = 0.
        let p1 = sb;
        let p2 = -(sa * sa) / 3.0;
        let mut p = p1 + p2;
        let ep = 5.0 * MACHINE_EPSILON * (p1.abs() + p2.abs());
        if p.abs() <= ep {
            p = 0.0;
        }
        let q1 = sc;
        let q2 = -sa * sb / 3.0;
        let q3 = 2.0 * (sa * sa * sa) / 27.0;
        let mut q = q1 + q2 + q3;
        let eq = 10.0 * MACHINE_EPSILON * (q1.abs() + q2.abs() + q3.abs());
        if q.abs() <= eq {
            q = 0.0;
        }
        // Check for overflow.
        if p.abs() > OVERFLOW_LIMIT {
            self.done = false;
            return;
        }
        // The discriminant.
        let a1 = (p * p * p) / 27.0;
        let a2 = (q * q) / 4.0;
        let mut discr = a1 + a2;
        // Special handling for P < 0.
        if p < 0.0 {
            discr = compute_special_discriminant(sa, sb, sc, a1);
        }
        let mut roots = [0.0; 4];
        let count = if discr < 0.0 {
            // Three distinct real roots.
            solve_cubic_three_real_roots(sa, sb, sc, p, q, discr, &mut roots);
            3
        } else if discr > 0.0 {
            // One real root.
            solve_cubic_one_real_root(sa, sc, p, q, discr, &mut roots);
            1
        } else {
            // Multiple roots.
            solve_cubic_multiple_roots(sa, sb, sc, p, q, &mut roots)
        };
        self.roots = roots[..count].to_vec();
        scale_and_refine_all_roots(&mut self.roots, scale_factor, &[a0, b0, c0, d0]);
    }

    fn solve_quadratic(&mut self, a0: f64, b0: f64, c0: f64) {
        // Check for degree reduction.
        if a0.abs() <= ZERO_THRESHOLD {
            self.solve_linear(b0, c0);
            return;
        }
        // The normalized x^2 + P x + Q = 0.
        let p = b0 / a0;
        let q = c0 / a0;
        // The discriminant, with its error bound.
        let eps_d = 3.0 * MACHINE_EPSILON * (p * p + (4.0 * q).abs());
        let mut discrim = p * p - 4.0 * q;
        if discrim.abs() <= eps_d {
            discrim = 0.0;
        }
        if discrim < 0.0 {
            // No real roots.
            self.roots.clear();
        } else if discrim == 0.0 {
            // A double root.
            let r = refine_root(&[1.0, p, q], -0.5 * p);
            self.roots = vec![r, r];
        } else {
            // Two distinct real roots, by the numerically stable formula.
            let mut r0 = if p > 0.0 { -(p + discrim.sqrt()) / 2.0 } else { -(p - discrim.sqrt()) / 2.0 };
            r0 = refine_root(&[1.0, p, q], r0);
            let mut r1 = q / r0;
            r1 = refine_root(&[1.0, p, q], r1);
            self.roots = vec![r0, r1];
        }
    }

    fn solve_linear(&mut self, a: f64, b: f64) {
        if a.abs() <= ZERO_THRESHOLD {
            if b.abs() <= ZERO_THRESHOLD {
                // 0 = 0: infinitely many solutions.
                self.infinite = true;
                return;
            }
            // 0 x + B = 0: no solution.
            self.roots.clear();
            return;
        }
        // The normal case: a single solution.
        self.roots = vec![-b / a];
    }
}

/// `math_TrigonometricEquationFunction`: `a cos²x + 2b cos x sin x + c cos x + d sin x + e`.
#[doc(alias = "math_TrigonometricEquationFunction")]
struct TrigonometricEquation {
    a: f64,
    b: f64,
    c: f64,
    d: f64,
    e: f64,
}

impl FunctionWithDerivative for TrigonometricEquation {
    fn value(&mut self, x: f64, f: &mut f64) -> bool {
        let (cn, sn) = (x.cos(), x.sin());
        *f = cn * (self.a * cn + (self.b + self.b) * sn + self.c) + self.d * sn + self.e;
        true
    }

    fn derivative(&mut self, x: f64, d: &mut f64) -> bool {
        let (cn, sn) = (x.cos(), x.sin());
        *d = -self.a * cn * sn + self.b * (cn * cn - sn * sn);
        *d += *d;
        *d += -self.c * sn + self.d * cn;
        true
    }

    fn values(&mut self, x: f64, f: &mut f64, d: &mut f64) -> bool {
        let (cn, sn) = (x.cos(), x.sin());
        let aacn = self.a * cn;
        let bbsn = self.b * sn;
        *f = aacn * cn + bbsn * (cn + cn) + self.c * cn + self.d * sn + self.e;
        *d = -aacn * sn + self.b * (cn * cn - sn * sn);
        *d += *d;
        *d += -self.c * sn + self.d * cn;
        true
    }

    fn get_state_number(&mut self) -> i32 {
        0
    }
}

/// The outcome of `math_TrigonometricFunctionRoots`.
#[doc(alias = "math_TrigonometricFunctionRoots")]
#[derive(Clone, Debug, Default)]
pub(crate) struct TrigonometricRoots {
    /// `IsDone()`.
    pub done: bool,
    /// `InfiniteRoots()`: the equation is zero.
    pub infinite: bool,
    /// `Value(1..=NbSolutions())`.
    pub roots: Vec<f64>,
}

/// `math_TrigonometricFunctionRoots(A, B, C, D, E, InfBound, SupBound)`: the roots of
/// `A cos²x + 2B cos x sin x + C cos x + D sin x + E` in `[inf_bound, sup_bound]` (an
/// interval of at most a period is searched).
#[doc(alias = "math_TrigonometricFunctionRoots")]
// (OCCT's 1-based index loops, kept as written.)
#[allow(clippy::too_many_arguments, clippy::needless_range_loop, clippy::manual_memcpy, clippy::mut_range_bound)]
pub(crate) fn trigonometric_roots(a: f64, b: f64, c: f64, d: f64, e: f64, inf_bound: f64, sup_bound: f64) -> TrigonometricRoots {
    // `Sol(1..4)`, 1-based (with room for OCCT's shifts past the fourth).
    let mut sol = [0.0f64; 7];
    let mut nb_sol: usize = 0;
    let mut zer = [0.0f64; 5];
    let mut nzer: usize;
    let nit = 10;
    let tol1 = 1.0e-15;
    let eps = 1.5e-12;
    let depi = PI + PI;
    let mut out = TrigonometricRoots { done: true, infinite: false, roots: Vec::new() };
    let finish = |out: &mut TrigonometricRoots, sol: &[f64; 7], nb_sol: usize| {
        out.roots = sol[1..=nb_sol].to_vec();
    };

    let (my_borne_inf, delta, modulo) = if inf_bound <= -f64::MAX && sup_bound >= f64::MAX {
        (0.0, depi, 0.0)
    } else if sup_bound >= f64::MAX {
        (inf_bound, depi, inf_bound / depi)
    } else if inf_bound <= -f64::MAX {
        let m = sup_bound - depi;
        (m, depi, m / depi)
    } else {
        let mut delta = sup_bound - inf_bound;
        if sup_bound - inf_bound > depi {
            delta = depi;
        }
        (inf_bound, delta, inf_bound / depi)
    };

    if a.abs() <= eps && b.abs() <= eps {
        if c.abs() <= eps {
            if d.abs() <= eps {
                if e.abs() <= eps {
                    // Infinitely many solutions.
                    out.infinite = true;
                    return out;
                }
                return out;
            }
            // An equation d sin(x) + e = 0.
            let aa = -e / d;
            if aa.abs() > 1.0 {
                return out;
            }
            zer[1] = aa.asin();
            zer[2] = PI - zer[1];
            nzer = 2;
            for i in 1..=nzer {
                if zer[i] <= -eps {
                    zer[i] = depi - zer[i].abs();
                }
                // The solutions brought between InfBound and SupBound.
                zer[i] += modulo.trunc() * depi;
                let x = zer[i] - my_borne_inf;
                if x > -epsilon(delta) && x < delta + epsilon(delta) {
                    nb_sol += 1;
                    sol[nb_sol] = zer[i];
                }
            }
            finish(&mut out, &sol, nb_sol);
            return out;
        } else if d.abs() <= eps {
            // An equation of the first degree c cos(x) + e = 0.
            let aa = -e / c;
            if aa.abs() > 1.0 {
                return out;
            }
            zer[1] = aa.acos();
            zer[2] = -zer[1];
            nzer = 2;
            for i in 1..=nzer {
                if zer[i] <= -eps {
                    zer[i] = depi - zer[i].abs();
                }
                zer[i] += modulo.trunc() * 2.0 * PI;
                let x = zer[i] - my_borne_inf;
                if x >= -epsilon(delta) && x <= delta + epsilon(delta) {
                    nb_sol += 1;
                    sol[nb_sol] = zer[i];
                }
            }
            finish(&mut out, &sol, nb_sol);
            return out;
        } else {
            // An equation of the second degree in the tangent of the half angle.
            let aa = e - c;
            let bb = 2.0 * d;
            let cc = e + c;
            let resol = PolynomialRoots::quadratic(aa, bb, cc);
            if !resol.done {
                out.done = false;
                return out;
            } else if !resol.infinite {
                nzer = resol.roots.len();
                for i in 1..=nzer {
                    zer[i] = resol.roots[i - 1];
                }
            } else {
                out.infinite = true;
                return out;
            }
        }
    } else {
        // Two additional analytical cases.
        if a.abs() <= eps && e.abs() <= eps {
            if c.abs() <= eps {
                // 2 B sin cos + D sin = 0.
                nzer = 2;
                zer[1] = 0.0;
                zer[2] = PI;
                let aa = -d / (b * 2.0);
                if aa.abs() <= 1.0 + P_CONFUSION {
                    nzer = 4;
                    if aa >= 1.0 {
                        zer[3] = 0.0;
                        zer[4] = 0.0;
                    } else if aa <= -1.0 {
                        zer[3] = PI;
                        zer[4] = PI;
                    } else {
                        zer[3] = aa.acos();
                        zer[4] = depi - zer[3];
                    }
                }
                special_case_solutions(&mut zer, nzer, my_borne_inf, eps, depi, modulo, delta, inf_bound, sup_bound, &mut sol, &mut nb_sol);
                finish(&mut out, &sol, nb_sol);
                return out;
            }
            if d.abs() <= eps {
                // 2 B sin cos + C cos = 0.
                nzer = 2;
                zer[1] = PI / 2.0;
                zer[2] = PI * 3.0 / 2.0;
                let aa = -c / (b * 2.0);
                if aa.abs() <= 1.0 + P_CONFUSION {
                    nzer = 4;
                    if aa >= 1.0 {
                        zer[3] = PI / 2.0;
                        zer[4] = PI / 2.0;
                    } else if aa <= -1.0 {
                        zer[3] = PI * 3.0 / 2.0;
                        zer[4] = PI * 3.0 / 2.0;
                    } else {
                        zer[3] = aa.asin();
                        zer[4] = PI - zer[3];
                    }
                }
                special_case_solutions(&mut zer, nzer, my_borne_inf, eps, depi, modulo, delta, inf_bound, sup_bound, &mut sol, &mut nb_sol);
                finish(&mut out, &sol, nb_sol);
                return out;
            }
        }

        // An equation of the fourth degree.
        let mut ko = [0.0f64; 6];
        ko[1] = a - c + e;
        ko[2] = 2.0 * d - 4.0 * b;
        ko[3] = 2.0 * e - 2.0 * a;
        ko[4] = 4.0 * b + 2.0 * d;
        ko[5] = a + c + e;
        loop {
            let mut bko = false;
            let resol4 = PolynomialRoots::quartic(ko[1], ko[2], ko[3], ko[4], ko[5]);
            if !resol4.done {
                out.done = false;
                return out;
            } else if !resol4.infinite {
                nzer = resol4.roots.len();
                for i in 1..=nzer {
                    zer[i] = resol4.roots[i - 1];
                }
            } else {
                out.infinite = true;
                return out;
            }
            // Sort the roots.
            loop {
                let mut triok = true;
                for i in 1..nzer {
                    if zer[i] > zer[i + 1] {
                        zer.swap(i, i + 1);
                        triok = false;
                    }
                }
                if triok {
                    break;
                }
            }
            for i in 1..nzer {
                if (zer[i + 1] - zer[i]).abs() < eps {
                    // A double root, or a numerical error?
                    let qw = zer[i + 1];
                    let va = ko[4] + qw * (2.0 * ko[3] + qw * (3.0 * ko[2] + qw * (4.0 * ko[1])));
                    if va.abs() > eps {
                        bko = true;
                        break;
                    }
                }
            }
            if bko {
                // A small coefficient: divide.
                for k in &mut ko[1..] {
                    *k *= 0.0001;
                }
            } else {
                break;
            }
        }
    }

    // Check the solutions against the bounds.
    let supm_infs100 = (sup_bound - inf_bound) * 0.01;
    nb_sol = 0;
    for &z in &zer[1..=nzer] {
        let mut teta = z.atan();
        teta += teta;
        if z <= -eps {
            teta = depi - teta.abs();
        }
        teta += modulo.trunc() * depi;
        if teta - my_borne_inf < 0.0 {
            teta += depi;
        }
        let x = teta - my_borne_inf;
        if x >= -epsilon(delta) && x <= delta + epsilon(delta) {
            // Newton's method (OCC541(apo)).
            let mut teta_newton = teta;
            let mut my_f = TrigonometricEquation { a, b, c, d, e };
            let resol = newton_function_root(&mut my_f, teta, tol1, eps, nit);
            if resol.done {
                teta_newton = resol.root;
            }
            // (lbr, 7 March 97: Newton converging very far from the starting point is
            // ignored.)
            let delta_newton = teta_newton - teta;
            if !(delta_newton > supm_infs100 || delta_newton < -supm_infs100) {
                teta = teta_newton;
            }
            // Insert in increasing order.
            let mut flag4 = false;
            for k in 1..=nb_sol {
                if teta < sol[k] {
                    for l in k..=nb_sol {
                        let j = nb_sol - l + k;
                        sol[j + 1] = sol[j];
                    }
                    sol[k] = teta;
                    nb_sol += 1;
                    flag4 = true;
                    break;
                }
            }
            if !flag4 {
                nb_sol += 1;
                sol[nb_sol] = teta;
            }
        }
    }
    // The special case of PI.
    if nb_sol < 4 {
        let start_index = nb_sol + 1;
        for sol_it in start_index..=4 {
            let teta = PI + modulo.trunc() * 2.0 * PI;
            let x = teta - my_borne_inf;
            if x >= -epsilon(delta) && x <= delta + epsilon(delta) && (a - c + e).abs() <= eps {
                let mut flag4 = false;
                let mut j = 0;
                for k in 1..=nb_sol {
                    j = k;
                    if teta < sol[k] {
                        flag4 = true;
                        break;
                    }
                    if sol_it == start_index && (teta - sol[k]).abs() <= eps {
                        finish(&mut out, &sol, nb_sol);
                        return out;
                    }
                }
                if !flag4 {
                    nb_sol += 1;
                    sol[nb_sol] = teta;
                } else {
                    for k in j..=nb_sol {
                        let i = nb_sol - k + j;
                        sol[i + 1] = sol[i];
                    }
                    sol[j] = teta;
                    nb_sol += 1;
                }
            }
        }
    }
    finish(&mut out, &sol, nb_sol);
    out
}

/// The solution loop shared by the two `2 B sin cos + ...` cases: each candidate brought
/// into the bounds by whole periods, kept within `PConfusion` of them and clamped.
#[allow(clippy::too_many_arguments, clippy::needless_range_loop)]
fn special_case_solutions(zer: &mut [f64; 5], nzer: usize, my_borne_inf: f64, eps: f64, depi: f64, modulo: f64, delta: f64, inf_bound: f64, sup_bound: f64, sol: &mut [f64; 7], nb_sol: &mut usize) {
    *nb_sol = 0;
    for i in 1..=nzer {
        if zer[i] <= my_borne_inf - eps {
            zer[i] += depi;
        }
        // The solutions brought between InfBound and SupBound.
        zer[i] += modulo.trunc() * 2.0 * PI;
        let x = zer[i] - my_borne_inf;
        if x >= -P_CONFUSION && x <= delta + P_CONFUSION {
            if zer[i] < inf_bound {
                zer[i] = inf_bound;
            }
            if zer[i] > sup_bound {
                zer[i] = sup_bound;
            }
            *nb_sol += 1;
            sol[*nb_sol] = zer[i];
        }
    }
}

/// The outcome of `math_NewtonFunctionRoot`.
#[doc(alias = "math_NewtonFunctionRoot")]
#[derive(Clone, Copy, Debug)]
pub(crate) struct NewtonFunctionRoot {
    /// `IsDone()`.
    pub done: bool,
    /// `Root()`: the best point met (smallest |F|).
    pub root: f64,
    /// `Derivative()`.
    pub derivative: f64,
    /// `Value()`.
    pub value: f64,
    /// `NbIterations()`.
    pub iterations: i32,
}

/// `math_NewtonFunctionRoot(F, Guess, EpsX, EpsF, A, B, NbIterations)`: Newton's method from
/// `guess`, its iterates clamped into `[a, b]`, until the step is within `eps_x` and the
/// value within `eps_f`.
#[doc(alias = "math_NewtonFunctionRoot")]
#[allow(clippy::too_many_arguments)]
pub(crate) fn newton_function_root_bounded(function: &mut dyn FunctionWithDerivative, guess: f64, eps_x: f64, eps_f: f64, a: f64, b: f64, iteration_limit: i32) -> NewtonFunctionRoot {
    newton_perform(function, guess, eps_x, eps_f, a, b, iteration_limit)
}

/// `math_NewtonFunctionRoot(F, Guess, EpsX, EpsF, NbIterations)`: between `RealFirst()`
/// and `RealLast()`.
#[doc(alias = "math_NewtonFunctionRoot")]
pub(crate) fn newton_function_root(function: &mut dyn FunctionWithDerivative, guess: f64, eps_x: f64, eps_f: f64, iteration_limit: i32) -> NewtonFunctionRoot {
    newton_perform(function, guess, eps_x, eps_f, REAL_FIRST, REAL_LAST, iteration_limit)
}

/// `math_NewtonFunctionRoot::Perform(F, Guess)`.
fn newton_perform(function: &mut dyn FunctionWithDerivative, guess: f64, eps_x: f64, eps_f: f64, binf: f64, bsup: f64, iteration_limit: i32) -> NewtonFunctionRoot {
    let mut x;
    let mut dfx = 0.0;
    // (lbr, 12 Nov 97: the best estimate was not saved, and a solution worse than the guess
    // could come back.)
    let mut best_x = REAL_LAST;
    let mut best_fx = REAL_LAST;
    let (aa, bb) = if binf < bsup { (binf, bsup) } else { (bsup, binf) };
    let mut dx = REAL_LAST;
    let mut fx = REAL_LAST;
    x = guess;
    let mut it = 1;
    while it <= iteration_limit && (dx.abs() > eps_x || fx.abs() > eps_f) {
        let ok = function.values(x, &mut fx, &mut dfx);
        let mut abs_fx = fx;
        if abs_fx < 0.0 {
            abs_fx = -abs_fx;
        }
        if abs_fx < best_fx {
            best_fx = abs_fx;
            best_x = x;
        }
        if ok {
            if dfx == 0.0 {
                it = iteration_limit + 1;
            } else {
                dx = fx / dfx;
                x -= dx;
                // Limit the variations of X.
                if x <= aa {
                    x = aa;
                }
                if x >= bb {
                    x = bb;
                }
                it += 1;
            }
        } else {
            it = iteration_limit + 1;
        }
    }
    x = best_x;
    // (The `Done = false` OCCT sets inside the loop is overwritten here.)
    let done = it <= iteration_limit;
    NewtonFunctionRoot { done, root: x, derivative: dfx, value: fx, iterations: it }
}
