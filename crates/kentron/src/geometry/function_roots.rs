//! OCCT's `math_FunctionRoots` (the `MATH_FUNCTIONROOTS_NEWCODE` path, the one OCCT compiles;
//! the `MATH_FUNCTIONROOTS_OLDCODE` path and its `math_DirectPolynomialRoots` are compiled out)
//! and the `math_BracketedRoot` it calls, ported line by line.
//!
//! The functions follow OCCT's `math_Function` protocol: `value(x, &mut f)` returns false when
//! it cannot evaluate, and then leaves in `f` whatever it left there (often the caller's
//! previous value). The roots finder reads such values anyway, as OCCT does; where OCCT reads
//! an uninitialized local this way, the port starts the local at 0.

/// `ITMAX`: the iteration limit of `Solve`.
const ITERATION_LIMIT: usize = 100;
/// `EPSEPS`: the relative part of `Solve`'s convergence tolerance.
const RELATIVE_EPSILON: f64 = 2e-14;

/// A function of one variable.
#[doc(alias = "math_Function")]
pub(crate) trait Function {
    /// `Value(X, F)`: false when the value cannot be computed (`value` is then left as the
    /// function leaves it).
    fn value(&mut self, x: f64, value: &mut f64) -> bool;
}

/// A function of one variable with its derivative.
#[doc(alias = "math_FunctionWithDerivative")]
pub(crate) trait FunctionWithDerivative {
    /// `Value(X, F)`.
    fn value(&mut self, x: f64, value: &mut f64) -> bool;

    /// `Derivative(X, D)`.
    fn derivative(&mut self, x: f64, derivative: &mut f64) -> bool;

    /// `Values(X, F, D)`.
    fn values(&mut self, x: f64, value: &mut f64, derivative: &mut f64) -> bool;

    /// `GetStateNumber()`: called on each root found, right after the root's `value`, so the
    /// function can save its state there.
    fn get_state_number(&mut self) -> i32;
}

/// The outcome of `math_BracketedRoot`.
#[doc(alias = "math_BracketedRoot")]
#[derive(Clone, Copy, Debug)]
pub(crate) struct BracketedRoot {
    /// `IsDone()`.
    pub done: bool,
    /// `Root()`.
    pub root: f64,
    /// `Value()`: the function's value at the root.
    pub value: f64,
    /// `NbIterations()`.
    pub iterations: i32,
}

/// `math_BracketedRoot(F, Bound1, Bound2, Tolerance, NbIterations, ZEPS)`: Brent's method for a
/// root of `function` between two bounds where it takes values of opposite signs.
#[doc(alias = "math_BracketedRoot")]
pub(crate) fn bracketed_root(function: &mut dyn Function, bound1: f64, bound2: f64, tolerance: f64, iteration_limit: i32, zeps: f64) -> BracketedRoot {
    // `Fa`, `c`, `d`, `e`, `TheError`; OCCT leaves `Fa` and `TheError` uninitialized until
    // the function sets them.
    let mut fa = 0.0;
    let mut c = 0.0;
    let mut d = 0.0;
    let mut e = 0.0;
    let mut a = bound1;
    let mut the_root = bound2;
    let mut the_error = 0.0;
    function.value(a, &mut fa);
    function.value(the_root, &mut the_error);
    if fa * the_error > 0.0 {
        return BracketedRoot { done: false, root: the_root, value: the_error, iterations: 0 };
    }
    let mut fc = the_error;
    let mut iteration = 1;
    while iteration <= iteration_limit {
        if the_error * fc > 0.0 {
            // Rename a, TheRoot, c and adjust the bounding interval d.
            c = a;
            fc = fa;
            d = the_root - a;
            e = d;
        }
        if fc.abs() < fa.abs() {
            a = the_root;
            the_root = c;
            c = a;
            fa = the_error;
            the_error = fc;
            fc = fa;
        }
        // Convergence check.
        let tol1 = 2.0 * zeps * the_root.abs() + 0.5 * tolerance;
        let xm = 0.5 * (c - the_root);
        if xm.abs() <= tol1 || the_error == 0.0 {
            return BracketedRoot { done: true, root: the_root, value: the_error, iterations: iteration };
        }
        if e.abs() >= tol1 && fa.abs() > the_error.abs() {
            // Attempt inverse quadratic interpolation.
            let s = the_error / fa;
            let mut p;
            let mut q;
            if a == c {
                p = 2.0 * xm * s;
                q = 1.0 - s;
            } else {
                q = fa / fc;
                let r = the_error / fc;
                p = s * (2.0 * xm * q * (q - r) - (the_root - a) * (r - 1.0));
                q = (q - 1.0) * (r - 1.0) * (s - 1.0);
            }
            // Check whether in bounds.
            if p > 0.0 {
                q = -q;
            }
            p = p.abs();
            let min1 = 3.0 * xm * q - (tol1 * q).abs();
            let min2 = (e * q).abs();
            if 2.0 * p < (if min1 < min2 { min1 } else { min2 }) {
                // Accept the interpolation.
                e = d;
                d = p / q;
            } else {
                // The interpolation failed, use bisection.
                d = xm;
                e = d;
            }
        } else {
            // Bounds decreasing too slowly, use bisection.
            d = xm;
            e = d;
        }
        // Move the last best guess to a.
        a = the_root;
        fa = the_error;
        if d.abs() > tol1 {
            the_root += d;
        } else {
            the_root += if xm > 0.0 { tol1.abs() } else { -tol1.abs() };
        }
        function.value(the_root, &mut the_error);
        iteration += 1;
    }
    BracketedRoot { done: false, root: the_root, value: the_error, iterations: iteration }
}

/// The outcome of `math_FunctionRoots`.
#[doc(alias = "math_FunctionRoots")]
#[derive(Clone, Debug, Default)]
pub(crate) struct FunctionRoots {
    /// `IsDone()`.
    pub done: bool,
    /// `IsAllNull()` (false where `done` is false: OCCT leaves it unset then).
    pub all_null: bool,
    /// `Value(1..=NbSolutions())`, in increasing order.
    pub solutions: Vec<f64>,
    /// `StateNumber(1..=NbSolutions())`.
    pub state_numbers: Vec<i32>,
}

/// `DerivFunction` (local to math_FunctionRoots.cxx): the derivative of a function, as a
/// function.
struct DerivativeFunction<'a> {
    function: &'a mut dyn FunctionWithDerivative,
}

impl Function for DerivativeFunction<'_> {
    fn value(&mut self, x: f64, value: &mut f64) -> bool {
        self.function.derivative(x, value)
    }
}

/// `AppendRoot` (static in math_FunctionRoots.cxx): inserts `x` among the sorted roots unless
/// one is within `dx` of it, saving the function's state at it.
fn append_root(solutions: &mut Vec<f64>, state_numbers: &mut Vec<i32>, x: f64, function: &mut dyn FunctionWithDerivative, _k: f64, dx: f64) {
    let n = solutions.len();
    let mut t = 0.0;
    if n == 0 {
        solutions.push(x);
        function.value(x, &mut t);
        state_numbers.push(function.get_state_number());
    } else {
        let mut i = 1;
        let mut place = n + 1;
        while i <= n {
            t = solutions[i - 1];
            if t >= x {
                place = i;
                i = n;
            }
            if (x - t).abs() <= dx {
                place = 0;
                i = n;
            }
            i += 1;
        }
        if place > n {
            solutions.push(x);
            function.value(x, &mut t);
            state_numbers.push(function.get_state_number());
        } else if place > 0 {
            solutions.insert(place - 1, x);
            function.value(x, &mut t);
            state_numbers.insert(place - 1, function.get_state_number());
        }
    }
}

/// `Solve` (static in math_FunctionRoots.cxx): Brent's method for the root of `F - K` between
/// `x1` and `x2`, polished by up to six Newton steps, then appended with [`append_root`].
#[allow(clippy::too_many_arguments)]
fn solve(function: &mut dyn FunctionWithDerivative, k: f64, x1: f64, y1: f64, x2: f64, y2: f64, tol: f64, dx: f64, solutions: &mut Vec<f64>, state_numbers: &mut Vec<i32>) {
    let half_tolerance = 0.5 * tol;
    let mut d = 0.0;
    let mut e = 0.0;
    let mut a = x1;
    let mut b = x2;
    let mut c = x2;
    let mut fa = y1;
    let mut fb = y2;
    let mut fc = y2;
    for _iteration in 1..=ITERATION_LIMIT {
        if (fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0) {
            c = a;
            fc = fa;
            d = b - a;
            e = d;
        }
        if fc.abs() < fb.abs() {
            a = b;
            b = c;
            c = a;
            fa = fb;
            fb = fc;
            fc = fa;
        }
        let tol1 = RELATIVE_EPSILON * b.abs() + half_tolerance;
        let xm = 0.5 * (c - b);
        if xm.abs() < tol1 || fb == 0.0 {
            // Try Newton iterations.
            let mut xp = b;
            let mut yp = 0.0;
            let mut dp = 0.0;
            let mut newton_steps_left = 5;
            loop {
                let mut ok = function.values(xp, &mut yp, &mut dp);
                if ok {
                    ok = false;
                    if dp.abs() > 1e-10 {
                        xp -= (yp - k) / dp;
                    }
                    if xp <= x2 && xp >= x1 {
                        function.value(xp, &mut yp);
                        yp -= k;
                        if yp.abs() < fb.abs() {
                            b = xp;
                            fb = yp;
                            ok = true;
                        }
                    }
                }
                // `while (Ok && --itern >= 0)`
                if !ok {
                    break;
                }
                newton_steps_left -= 1;
                if newton_steps_left < 0 {
                    break;
                }
            }
            append_root(solutions, state_numbers, b, function, k, dx);
            return;
        }
        if e.abs() >= tol1 && fa.abs() > fb.abs() {
            let s = fb / fa;
            let mut p;
            let mut q;
            if a == c {
                p = xm * s;
                p += p;
                q = 1.0 - s;
            } else {
                q = fa / fc;
                let r = fb / fc;
                p = s * ((xm + xm) * q * (q - r) - (b - a) * (r - 1.0));
                q = (q - 1.0) * (r - 1.0) * (s - 1.0);
            }
            if p > 0.0 {
                q = -q;
            }
            p = p.abs();
            let min1 = 3.0 * xm * q - (tol1 * q).abs();
            let min2 = (e * q).abs();
            if (p + p) < (if min1 < min2 { min1 } else { min2 }) {
                e = d;
                d = p / q;
            } else {
                d = xm;
                e = d;
            }
        } else {
            d = xm;
            e = d;
        }
        a = b;
        fa = fb;
        if d.abs() > tol1 {
            b += d;
        } else if xm >= 0.0 {
            b += tol1.abs();
        } else {
            b += -tol1.abs();
        }
        function.value(b, &mut fb);
        fb -= k;
    }
    // No convergence: OCCT appends nothing.
}

/// `math_FunctionRoots(F, A, B, NbSample, EpsX, EpsF, EpsNull, K)`: the roots of `F - K` on
/// `[A, B]`, found by sampling `max(2 NbSample, 20)` intervals, solving every sign change,
/// keeping sampled zeros and end values within `EpsF`, and searching the neighbourhood of
/// sampled minima of `|F - K|` for tangential roots.
#[doc(alias = "math_FunctionRoots")]
#[allow(clippy::too_many_arguments)]
pub(crate) fn function_roots(function: &mut dyn FunctionWithDerivative, a: f64, b: f64, sample_count: i32, eps_x_given: f64, eps_f: f64, eps_null: f64, k: f64) -> FunctionRoots {
    let mut solutions = Vec::new();
    let mut state_numbers = Vec::new();
    let mut x0 = a;
    let mut xn = b;
    let mut n = sample_count;
    //-- Low level checks.
    if b < a {
        x0 = b;
        xn = a;
    }
    n *= 2;
    if n < 20 {
        n = 20;
    }
    //-- Is EpsX too small (U + Nn*EpsX == U)?
    let mut eps_x = eps_x_given;
    let delta_u = x0.abs() + xn.abs();
    let n_eps_x = 0.0000000001 * delta_u;
    if eps_x < n_eps_x {
        eps_x = n_eps_x;
    }

    //-- Sample F - K at X0, X1, ..., XN.
    let count = n as usize;
    let mut x = x0;
    let dx = (xn - x0) / n as f64;
    let mut sampled = vec![0.0; count + 1];
    let mut valid: i64 = -1;
    let mut aux = 0.0;
    for _ in 0..=count {
        if x > xn {
            x = xn;
        }
        if function.value(x, &mut aux) {
            valid += 1;
            sampled[valid as usize] = aux - k;
        }
        x += dx;
    }
    if valid < n as i64 {
        return FunctionRoots { done: false, all_null: false, solutions, state_numbers };
    }

    //-- Is the whole function null?
    let mut all_null = true;
    let mut i = 0;
    while all_null && i <= count {
        if sampled[i] > eps_null || sampled[i] < -eps_null {
            all_null = false;
        }
        i += 1;
    }
    if !all_null {
        //-- Some samples are out of tolerance: solve the STRICT sign changes.
        let tol = eps_x;
        x = x0;
        for i in 0..count {
            let ip1 = i + 1;
            let mut x2 = x + dx;
            if x2 > xn {
                x2 = xn;
            }
            if sampled[i] < 0.0 {
                if sampled[ip1] > 0.0 {
                    solve(function, k, x, sampled[i], x2, sampled[ip1], tol, n_eps_x, &mut solutions, &mut state_numbers);
                }
            } else if sampled[ip1] < 0.0 {
                solve(function, k, x, sampled[i], x2, sampled[ip1], tol, n_eps_x, &mut solutions, &mut state_numbers);
            }
            x += dx;
        }
        //-- Samples where F - K vanishes exactly: with two points u0, u1 around Xi, solve if
        //-- (F(u0)-K)*(F(u1)-K) < 0, else insert Xi itself unless both are zero.
        for (i, &sample) in sampled.iter().enumerate() {
            if sample == 0.0 {
                x = x0 + i as f64 * dx;
                if x > xn {
                    x = xn;
                }
                let mut u0 = dx * 0.5;
                let mut u1 = x + u0;
                u0 += x;
                if u0 < x0 {
                    u0 = x0;
                }
                if u0 > xn {
                    u0 = xn;
                }
                if u1 < x0 {
                    u1 = x0;
                }
                if u1 > xn {
                    u1 = xn;
                }
                let mut y0 = 0.0;
                let mut y1 = 0.0;
                function.value(u0, &mut y0);
                y0 -= k;
                function.value(u1, &mut y1);
                y1 -= k;
                if y0 * y1 < 0.0 {
                    solve(function, k, u0, y0, u1, y1, tol, n_eps_x, &mut solutions, &mut state_numbers);
                } else if y0 != 0.0 || y1 != 0.0 {
                    append_root(&mut solutions, &mut state_numbers, x, function, k, n_eps_x);
                }
            }
        }
        //-- The end points.
        if sampled[0] <= eps_f && sampled[0] >= -eps_f {
            append_root(&mut solutions, &mut state_numbers, x0, function, k, n_eps_x);
        }
        if sampled[count] <= eps_f && sampled[count] >= -eps_f {
            append_root(&mut solutions, &mut state_numbers, xn, function, k, n_eps_x);
        }

        //-- Sampled minima with f(x) > 0 and maxima with f(x) < 0: search their neighbourhood
        //-- more finely.
        // `ym` and `dym` live across iterations in OCCT (declared before the loop).
        let mut ym = 0.0;
        let mut dym = 0.0;
        let major_dx = 5.0 * dx;
        let mut im1 = 0;
        let mut ip1 = 2;
        let mut xm = x0 + dx;
        let mut i = 1;
        while i < count {
            let mut rediscretize = false;
            if xm > xn {
                xm = xn;
            }
            if sampled[i] > 0.0 {
                if sampled[im1] > sampled[i] && sampled[ip1] > sampled[i] {
                    //-- Can it cross the axis? Estimate from Xim1.
                    let mut xm1 = xm - dx;
                    if xm1 < x0 {
                        xm1 = x0;
                    }
                    function.values(xm1, &mut ym, &mut dym);
                    ym -= k;
                    if dym.abs() > 1e-10 {
                        // Normally dym < 0.
                        let t = ym / dym;
                        if t < major_dx && t > -major_dx {
                            rediscretize = true;
                        }
                    }
                    //-- Estimate from Xip1.
                    if !rediscretize {
                        let mut xp1 = xm + dx;
                        if xp1 > xn {
                            xp1 = xn;
                        }
                        function.values(xp1, &mut ym, &mut dym);
                        ym -= k;
                        if dym.abs() > 1e-10 {
                            // Normally dym > 0.
                            let t = ym / dym;
                            if t < major_dx && t > -major_dx {
                                rediscretize = true;
                            }
                        }
                    }
                }
            } else if sampled[i] < 0.0 && sampled[im1] < sampled[i] && sampled[ip1] < sampled[i] {
                //-- Can it cross the axis? Estimate from Xim1.
                let mut xm1 = xm - dx;
                if xm1 < x0 {
                    xm1 = x0;
                }
                function.values(xm1, &mut ym, &mut dym);
                ym -= k;
                if dym.abs() > 1e-10 {
                    // Normally dym > 0.
                    let t = ym / dym;
                    if t < major_dx && t > -major_dx {
                        rediscretize = true;
                    }
                }
                //-- OCCT's second estimate is again from Xim1 (not Xip1), kept as is.
                if !rediscretize {
                    let mut xm1 = xm - dx;
                    if xm1 < x0 {
                        xm1 = x0;
                    }
                    function.values(xm1, &mut ym, &mut dym);
                    ym -= k;
                    if dym.abs() > 1e-10 {
                        // Normally dym < 0.
                        let t = ym / dym;
                        if t < major_dx && t > -major_dx {
                            rediscretize = true;
                        }
                    }
                }
            }
            if rediscretize {
                search_near_extremum(function, k, &sampled, im1, ip1, xm, dx, x0, xn, eps_x_given, eps_f, tol, n_eps_x, &mut solutions, &mut state_numbers);
            }
            xm += dx;
            i += 1;
            im1 += 1;
            ip1 += 1;
        }
    }
    FunctionRoots { done: true, all_null, solutions, state_numbers }
}

/// The `if (Rediscr)` block of `math_FunctionRoots`: around the sample `xm` (a positive
/// minimum or negative maximum of the samples), the zero of the derivative by
/// `math_BracketedRoot` and the extremum by golden-section search (solving any sign change
/// it meets), the better of the two appended when within `EpsF`.
#[allow(clippy::too_many_arguments)]
fn search_near_extremum(
    function: &mut dyn FunctionWithDerivative,
    k: f64,
    sampled: &[f64],
    im1: usize,
    ip1: usize,
    xm: f64,
    dx: f64,
    lower: f64,
    upper: f64,
    eps_x_given: f64,
    eps_f: f64,
    tol: f64,
    n_eps_x: f64,
    solutions: &mut Vec<f64>,
    state_numbers: &mut Vec<i32>,
) {
    let mut x0 = xm - dx;
    let mut x3 = xm + dx;
    if x0 < lower {
        x0 = lower;
    }
    if x3 > upper {
        x3 = upper;
    }
    let mut solution1 = 0.0;
    let mut solution2 = 0.0;
    let mut value1 = 0.0;
    let mut value2 = 0.0;
    let mut derivative1 = 0.0;
    let mut derivative2 = 0.0;
    let mut is_solution1 = false;
    let mut is_solution2 = false;
    //-- The minimum of |F| between x0 and x3, as the zero of the derivative.
    let root = bracketed_root(&mut DerivativeFunction { function: &mut *function }, x0, x3, eps_x_given, 100, 1.0e-12);
    if root.done {
        solution1 = root.root;
        function.value(solution1, &mut value1);
        value1 = value1.abs();
        if value1 < eps_f {
            is_solution1 = true;
            derivative1 = root.value;
        }
    }

    //-- An extremum between x0 and x3, with x0 < x1 < x2 < x3, |f(x0)| > |f(x1)| and
    //-- |f(x3)| > |f(x2)|. On entry a = xm - dx, b = xm, c = xm + dx.
    let mut x1;
    let mut x2;
    let golden = 0.61803399;
    let complement = 1.0 - golden;
    let tol_cr = n_eps_x * 10.0;
    let mut f0 = sampled[im1];
    let mut f3 = sampled[ip1];
    let searching_minimum = f0 > 0.0;
    if (x3 - xm).abs() > (x0 - xm).abs() {
        x1 = xm;
        x2 = xm + complement * (x3 - xm);
    } else {
        x2 = xm;
        x1 = xm - complement * (xm - x0);
    }
    let mut f1 = 0.0;
    let mut f2 = 0.0;
    function.value(x1, &mut f1);
    f1 -= k;
    function.value(x2, &mut f2);
    f2 -= k;
    let tol_x = 0.001 * n_eps_x;
    while (x3 - x0).abs() > tol_cr * (x1.abs() + x2.abs()) && (x1 - x2).abs() > tol_x {
        let toward_x3 = if searching_minimum { f2 < f1 } else { f2 > f1 };
        if toward_x3 {
            x0 = x1;
            x1 = x2;
            x2 = golden * x1 + complement * x3;
            f0 = f1;
            f1 = f2;
            function.value(x2, &mut f2);
            f2 -= k;
        } else {
            x3 = x2;
            x2 = x1;
            x1 = golden * x2 + complement * x0;
            f3 = f2;
            f2 = f1;
            function.value(x1, &mut f1);
            f1 -= k;
        }
        //-- Not only extrema: check whether a root is crossed.
        if f1 * f0 < 0.0 {
            solve(function, k, x0, f0, x1, f1, tol, n_eps_x, solutions, state_numbers);
        }
        if f2 * f3 < 0.0 {
            solve(function, k, x2, f2, x3, f3, tol, n_eps_x, solutions, state_numbers);
        }
    }
    if (searching_minimum && f1 < f2) || (!searching_minimum && f1 > f2) {
        //-- x1, f(x1) minimum.
        if f1.abs() < eps_f {
            is_solution2 = true;
            solution2 = x1;
            value2 = f1.abs();
        }
    } else {
        //-- x2, f(x2) minimum.
        if f2.abs() < eps_f {
            is_solution2 = true;
            solution2 = x2;
            value2 = f2.abs();
        }
    }
    // The better of solution1 and solution2.
    if is_solution1 && is_solution2 {
        if value2 - value1 > eps_f {
            append_root(solutions, state_numbers, solution1, function, k, n_eps_x);
        } else if value1 - value2 > eps_f {
            append_root(solutions, state_numbers, solution2, function, k, n_eps_x);
        } else {
            derivative1 = derivative1.abs();
            function.derivative(solution2, &mut derivative2);
            derivative2 = derivative2.abs();
            if derivative1 < derivative2 {
                append_root(solutions, state_numbers, solution1, function, k, n_eps_x);
            } else {
                append_root(solutions, state_numbers, solution2, function, k, n_eps_x);
            }
        }
    } else if is_solution1 {
        append_root(solutions, state_numbers, solution1, function, k, n_eps_x);
    } else if is_solution2 {
        append_root(solutions, state_numbers, solution2, function, k, n_eps_x);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A plain function with its derivative, for the roots finder.
    struct Plain<F: Fn(f64) -> (f64, f64)> {
        f: F,
        states: i32,
    }

    impl<F: Fn(f64) -> (f64, f64)> FunctionWithDerivative for Plain<F> {
        fn value(&mut self, x: f64, value: &mut f64) -> bool {
            *value = (self.f)(x).0;
            true
        }
        fn derivative(&mut self, x: f64, derivative: &mut f64) -> bool {
            *derivative = (self.f)(x).1;
            true
        }
        fn values(&mut self, x: f64, value: &mut f64, derivative: &mut f64) -> bool {
            (*value, *derivative) = (self.f)(x);
            true
        }
        fn get_state_number(&mut self) -> i32 {
            self.states += 1;
            self.states
        }
    }

    #[test]
    fn roots_of_a_cubic_in_order() {
        let mut f = Plain { f: |x: f64| ((x - 0.3) * (x + 0.5) * (x - 1.2), 3.0 * x * x - 2.0 * x - 0.39), states: 0 };
        let roots = function_roots(&mut f, -1.0, 2.0, 10, 1e-9, 1e-10, 1e-10, 0.0);
        assert!(roots.done && !roots.all_null);
        assert_eq!(roots.solutions.len(), 3);
        for (found, expected) in roots.solutions.iter().zip([-0.5, 0.3, 1.2]) {
            assert!((found - expected).abs() < 1e-12, "{found} vs {expected}");
        }
        assert_eq!(roots.state_numbers.len(), 3);
    }

    #[test]
    fn tangential_root_found_by_rediscretization() {
        // (x - 0.31)^2 touches zero between samples without changing sign.
        let mut f = Plain { f: |x: f64| ((x - 0.31) * (x - 0.31), 2.0 * (x - 0.31)), states: 0 };
        let roots = function_roots(&mut f, 0.0, 1.0, 10, 1e-9, 1e-10, 1e-10, 0.0);
        assert!(roots.done);
        assert_eq!(roots.solutions.len(), 1);
        assert!((roots.solutions[0] - 0.31).abs() < 1e-5, "{:?}", roots.solutions);
    }

    #[test]
    fn brent_on_a_bracket() {
        struct Shifted;
        impl Function for Shifted {
            fn value(&mut self, x: f64, value: &mut f64) -> bool {
                *value = x * x - 2.0;
                true
            }
        }
        let root = bracketed_root(&mut Shifted, 0.0, 2.0, 1e-12, 100, 1e-12);
        assert!(root.done);
        assert!((root.root - 2f64.sqrt()).abs() < 1e-11);
    }
}
