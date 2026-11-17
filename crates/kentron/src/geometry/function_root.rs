//! OCCT's `math_FunctionRoot` (the root of a function of one variable near a guess, between
//! bounds) and the solver it hands the work to, `math_FunctionSetRoot`, ported for the one
//! case `math_FunctionRoot` uses: one variable, one equation.
//!
//! With a single unknown every vector of `math_FunctionSetRoot` is one number and every
//! matrix one entry, so the port keeps scalars; the arithmetic is OCCT's (`Norm2()` of a
//! one-vector is its square, `Norm()` the square root of that). What the one-variable case
//! cannot reach is left out: the constrained sub-problem of `SearchDirection` for a strict
//! subset of the unknowns at a bound (with one unknown it is either free or at a bound), the
//! non-square solvers (`math_GaussLeastSquare`, `math_SVD` of a wide matrix), and the
//! `theStopOnDivergent` mode (`math_FunctionRoot` calls `Perform` with its default, false).

use super::adaptor::epsilon;
use super::brent_minimum::BrentMinimum;
use super::function_roots::FunctionWithDerivative;

/// `Precision::IsInfinite`: beyond half of `Precision::Infinite()` (2e100).
fn is_infinite(x: f64) -> bool {
    x.abs() >= 0.5 * 2e100
}

/// `Precision::SquareConfusion()`.
const SQUARE_CONFUSION: f64 = 1e-14;

/// The outcome of `math_FunctionRoot`.
#[doc(alias = "math_FunctionRoot")]
#[derive(Clone, Copy, Debug)]
pub(crate) struct FunctionRoot {
    /// `IsDone()`.
    pub done: bool,
    /// `Root()`.
    pub root: f64,
    /// `Derivative()`: the derivative at the last point the solver evaluated.
    pub derivative: f64,
    /// `NbIterations()`.
    pub iterations: usize,
}

/// `math_FunctionRoot(F, Guess, Tolerance, A, B, NbIterations)`: the root of `function`
/// starting from `guess`, kept inside `[a, b]`, to within `tolerance` on the variable.
/// (`math_FunctionRoot` then evaluates the function once more at the root for `Value()`;
/// nothing here reads that value, so the port skips the evaluation.)
#[doc(alias = "math_FunctionRoot")]
pub(crate) fn function_root(function: &mut dyn FunctionWithDerivative, guess: f64, tolerance: f64, a: f64, b: f64, iteration_limit: usize) -> FunctionRoot {
    let mut solver = FunctionSetRoot::new(tolerance, iteration_limit);
    solver.perform(function, guess, a, b);
    if solver.done {
        function.get_state_number();
    }
    FunctionRoot { done: solver.done, root: solver.sol, derivative: solver.df, iterations: solver.kount }
}

/// `MyDirFunction`: the function along the line `p0 + x * dir`, as half its square.
#[derive(Clone, Copy)]
struct DirFunction {
    p0: f64,
    dir: f64,
}

impl DirFunction {
    /// `MyDirFunction::Value(x, fval)`.
    fn value(&self, function: &mut dyn FunctionWithDerivative, x: f64) -> Option<f64> {
        let p = self.dir * x + self.p0;
        let mut fv = 0.0;
        if !function.value(p, &mut fv) {
            return None;
        }
        // (Precision::HalfInfinite later, says OCCT.)
        if fv <= -1.0e100 || fv >= 1.0e100 {
            return None;
        }
        Some(0.5 * (fv * fv))
    }
}

/// `MyDirFunction::Value(Sol, FF, DF, GH, F2, Gnr1)`: the function, its derivative, the
/// gradient of half its square and that gradient's square norm at `sol`.
struct Evaluation {
    ff: f64,
    df: f64,
    gh: f64,
    f2: f64,
    gnr1: f64,
}

fn evaluate(function: &mut dyn FunctionWithDerivative, sol: f64) -> Option<Evaluation> {
    let (mut ff, mut df) = (0.0, 0.0);
    if !function.values(sol, &mut ff, &mut df) {
        return None;
    }
    if ff < 0.0 {
        if ff <= -1.0e100 {
            return None;
        }
    } else if ff >= 1.0e100 {
        return None;
    }
    let f2 = 0.5 * (ff * ff);
    // GH = DF^T FF.
    let gh = df * ff;
    if is_infinite(gh) {
        return None;
    }
    Some(Evaluation { ff, df, gh, f2, gnr1: gh * gh })
}

/// `math_SVD(A).Solve(B, X)` for a 1x1 matrix, through `SVD_Decompose` and `SVD_Solve` of
/// math_Recipes as they run on one entry.
fn svd_solve(a: f64, b: f64) -> f64 {
    // SVD_Decompose: the Householder reduction of the single column.
    let mut u = a;
    let mut g = 0.0;
    let scale = if u > 0.0 { u } else { -u };
    if scale != 0.0 {
        u /= scale;
        let s = u * u;
        let f = u;
        // SIGN(a, b) is (b >= 0 ? |a| : -|a|).
        g = -(if f >= 0.0 { s.sqrt().abs() } else { -s.sqrt().abs() });
        u = f - g;
        u *= scale;
    }
    let mut w = scale * g;
    let mut v = 1.0;
    // The accumulation of the left-hand transformation.
    let g = w;
    if g != 0.0 {
        let g = 1.0 / g;
        u *= g;
    } else {
        u = 0.0;
    }
    u += 1.0;
    // The diagonalization stops at once (no super-diagonal); a negative singular value is
    // made positive with the sign moved into V.
    if w < 0.0 {
        w = -w;
        v = -v;
    }
    // math_SVD::Solve zeroes singular values below 1e-6 times the largest (never the case
    // for a single one), then SVD_Solve.
    let wmin = 1.0e-6 * w;
    if w < wmin {
        w = 0.0;
    }
    let mut s = 0.0;
    if w != 0.0 {
        s += u * b;
        s /= w;
    }
    let mut x = 0.0;
    x += v * s;
    x
}

/// `SearchDirection(DF, GH, FF, ChangeDirection, InvLengthMax, Direction, Dy)` (the
/// unconstrained one): Newton's step, shortened to a quarter of the domain, or the steepest
/// descent where Newton's does not descend. `change_direction` is taken by value, as in
/// OCCT. Returns `(Direction, Dy)`.
fn search_direction(df: f64, gh: f64, ff: f64, mut change_direction: bool, inv_length_max: f64, previous_direction: f64) -> (f64, f64) {
    let eps = 1.0e-32;
    let mut direction = previous_direction;
    if !change_direction {
        direction = -ff;
        // math_Gauss(DF, 1.e-9): LU_Decompose finds a 1x1 matrix singular when its entry is
        // not larger than the pivot threshold in size (NaN included).
        if df.abs() > 1.0e-9 {
            // LU_Solve's back-substitution.
            direction /= df;
        } else {
            // "We have to forget singular directions": math_SVD (which always succeeds on
            // one entry).
            direction = svd_solve(df, -ff);
        }
    }
    // Directions too long are cut to a quarter of the domain.
    let ratio = (direction * inv_length_max).abs();
    if ratio > 1.0 {
        direction /= ratio;
    }
    let mut dy = direction * gh;
    if dy >= -eps {
        // Newton does not descend: take the gradient.
        change_direction = true;
    }
    if change_direction {
        direction = -gh;
        dy = -(gh * gh);
    }
    (direction, dy)
}

/// The constrained `SearchDirection(DF, GH, FF, Constraints, X, ...)`: with the one unknown
/// at a bound there is nothing left to do.
fn search_direction_constrained(df: f64, gh: f64, ff: f64, constrained: bool, change_direction: bool, inv_length_max: f64, previous_direction: f64) -> (f64, f64) {
    if !constrained {
        search_direction(df, gh, ff, change_direction, inv_length_max, previous_direction)
    } else {
        (0.0, 0.0)
    }
}

/// `MinimizeDirection(P0, P1, P2, F1, Delta, Tol, F)`: a Brent search along `p1 - p0`
/// through `p1` for a value below `f1`; on success `delta` becomes the step from `p1`.
fn minimize_direction_brent(function: &mut dyn FunctionWithDerivative, p0: f64, p1: f64, p2: f64, f1: f64, delta: &mut f64, tol: f64) -> bool {
    // (1) A one-dimensional parametric tolerance.
    let mut tol1d: f64 = 2.1;
    let eps = 1.0e-16;
    let mut invnorme = delta.abs();
    if invnorme > eps {
        tol1d = tol1d.min(tol / invnorme);
    }
    if tol1d > 1.9 {
        return false;
    }
    tol1d /= 3.0;
    *delta = p1 - p0;
    invnorme = (*delta * *delta).sqrt();
    if invnorme <= eps {
        return false;
    }
    invnorme = 1.0 / invnorme;
    let line = DirFunction { p0: p1, dir: *delta };
    // (2) Minimize.
    let ax = -1.0;
    let bx = 0.0;
    let cx = ((p2 - p1) * (p2 - p1)).sqrt() * invnorme;
    if cx < 1.0e-2 {
        return false;
    }
    let mut solver = BrentMinimum::new(tol1d, 100, tol1d);
    solver.perform(&mut |x| line.value(function, x), ax, bx, cx);
    if solver.is_done() {
        let tsol = solver.location();
        if solver.minimum() < f1 {
            *delta *= tsol;
            return true;
        }
    }
    false
}

/// `MinimizeDirection(P, Dir, PValue, PDirValue, Gradient, DGradient, Tol, F)`: a quadratic
/// interpolation along `dir` from `p`, then a Brent search if it gained too little; on
/// success `dir` is scaled to the step taken.
#[allow(clippy::too_many_arguments)]
fn minimize_direction_interpolated(function: &mut dyn FunctionWithDerivative, p: f64, dir: &mut f64, p_value: f64, p_dir_value: f64, gradient: f64, d_gradient: f64, tol: f64) -> bool {
    if is_infinite(p_value) || is_infinite(p_dir_value) {
        return false;
    }
    // (0) A one-dimensional parametric tolerance.
    let mut good = false;
    let eps = 1.0e-20;
    let mut tol1d: f64 = 1.1;
    let mut result = p_value;
    let absdir = dir.abs();
    if absdir > eps {
        tol1d = tol1d.min(tol / absdir);
    }
    if tol1d > 0.9 {
        return false;
    }

    // (1) A first quadratic interpolation.
    let (mut ax, mut bx, mut cx);
    let mut tsol;
    let df1 = gradient * *dir;
    let df2 = d_gradient * *dir;
    if df1 < -eps && df2 > eps {
        // A basin.
        tsol = -df1 / (df2 - df1);
    } else {
        cx = p_value;
        bx = df1;
        ax = p_dir_value - (bx + cx);
        if ax.abs() <= eps {
            // Linear.
            tsol = if bx.abs() >= eps { -cx / bx } else { 0.0 };
        } else {
            // Quadratic.
            let mut delta = bx * bx - 4.0 * ax * cx;
            if delta > 1.0e-9 {
                // Roots: take the one nearest 0.
                delta = delta.sqrt();
                tsol = -(bx + delta);
                let tsolbis = delta - bx;
                if tsolbis.abs() < tsol.abs() {
                    tsol = tsolbis;
                }
                tsol /= 2.0 * ax;
            } else {
                // No root or hardly any: take the extremum.
                tsol = -(0.5 * bx) / ax;
            }
        }
    }
    if tsol.abs() >= 1.0 {
        return false;
    }

    let line = DirFunction { p0: p, dir: *dir };
    // (OCCT ignores whether this evaluation succeeds; its value is then uninitialized. The
    // port reads it as infinite, which keeps the point from being taken.)
    let fsol = line.value(function, tsol).unwrap_or(f64::INFINITY);
    if fsol < p_value {
        good = true;
        result = fsol;
    }

    // (2) Too little progress: a proper search from the same start.
    if fsol > 0.2 * p_value && tol1d < 0.5 {
        if tsol < 0.0 {
            ax = tsol;
            bx = 0.0;
            cx = 1.0;
        } else {
            ax = 0.0;
            bx = tsol;
            cx = 1.0;
        }
        let mut solver = BrentMinimum::new(tol1d, 100, tol1d);
        solver.perform(&mut |x| line.value(function, x), ax, bx, cx);
        if solver.is_done() && solver.minimum() <= result {
            tsol = solver.location();
            good = true;
            result = solver.minimum();
            // The objective changes too fast: search both halves again.
            if gradient * gradient > 1.0 / SQUARE_CONFUSION && tsol > ax && tsol < cx {
                solver.perform(&mut |x| line.value(function, x), ax, (ax + tsol) / 2.0, tsol);
                if solver.is_done() && solver.minimum() <= result {
                    tsol = solver.location();
                    good = true;
                    result = solver.minimum();
                }
                solver.perform(&mut |x| line.value(function, x), tsol, (cx + tsol) / 2.0, cx);
                if solver.is_done() && solver.minimum() <= result {
                    tsol = solver.location();
                    good = true;
                }
            }
        }
    }
    if good {
        *dir *= tsol;
    }
    good
}

/// `math_FunctionSetRoot` for one variable and one equation.
#[doc(alias = "math_FunctionSetRoot")]
struct FunctionSetRoot {
    /// `Delta`.
    delta: f64,
    /// `Sol`.
    sol: f64,
    /// `DF`.
    df: f64,
    /// `Tol`.
    tol: f64,
    done: bool,
    /// `Kount`: the iteration reached.
    kount: usize,
    /// `Itermax`.
    iteration_limit: usize,
    sol_save: f64,
    gh: f64,
    dh: f64,
    dh_save: f64,
    ff: f64,
    previous_solution: f64,
    /// `Save(0..Itermax)`: the half square of the function after each iteration.
    save: Vec<f64>,
}

/// `Bounds(InfBound, SupBound, Tol, Sol, SolSave, Constraints, Delta, isNewSol)`: brings a
/// step that left `[inf, sup]` back onto the bound. Returns `(Out, constrained, isNewSol)`.
fn bounds(inf: f64, sup: f64, tol: f64, sol: &mut f64, sol_save: f64, delta: &mut f64) -> (bool, bool, bool) {
    let mut out = false;
    let mut monratio: f64 = 1.0;
    let mut is_new_sol = true;
    let mut constrained = false;
    // The ratio to bring the step back by.
    *delta = *sol - sol_save;
    if inf == sup {
        constrained = true;
        // ("Ok, but this should be avoided.")
        out = true;
    } else if *sol < inf {
        constrained = true;
        out = true;
        // (Delta is negative.) Avoid zero ratios for nothing.
        if -*delta > tol {
            monratio = monratio.min((inf - sol_save) / *delta);
        }
    } else if *sol > sup {
        constrained = true;
        out = true;
        if *delta > tol {
            monratio = monratio.min((sup - sol_save) / *delta);
        }
    }
    if out {
        // Truncate, and clamp once more against numerical trouble.
        if monratio == 0.0 {
            is_new_sol = false;
            *sol = sol_save;
            *delta = 0.0;
        } else {
            *delta *= monratio;
            *sol = sol_save + *delta;
            if *sol < inf {
                *sol = inf;
                *delta = *sol - sol_save;
            } else if *sol > sup {
                *sol = sup;
                *delta = *sol - sol_save;
            }
        }
    }
    (out, constrained, is_new_sol)
}

impl FunctionSetRoot {
    /// `math_FunctionSetRoot(F, Tolerance, NbIterations)`.
    fn new(tolerance: f64, iteration_limit: usize) -> Self {
        FunctionSetRoot {
            delta: 0.0,
            sol: 0.0,
            df: 0.0,
            tol: tolerance,
            done: false,
            kount: 0,
            iteration_limit,
            sol_save: 0.0,
            gh: 0.0,
            dh: 0.0,
            dh_save: 0.0,
            ff: 0.0,
            previous_solution: 0.0,
            save: vec![0.0; iteration_limit + 1],
        }
    }

    /// `IsSolutionReached`: the last step within the tolerance.
    fn is_solution_reached(&self) -> bool {
        self.delta.abs() <= self.tol
    }

    /// Stores an evaluation's results (`FF`, `DF`, `GH`), returning `(F2, Gnr1)`.
    fn take(&mut self, e: Evaluation) -> (f64, f64) {
        self.ff = e.ff;
        self.df = e.df;
        self.gh = e.gh;
        (e.f2, e.gnr1)
    }

    /// `Perform(F, StartingPoint, InfBound, SupBound)`.
    #[doc(alias = "math_FunctionSetRoot::Perform")]
    // (Some assignments are overwritten before being read, as in OCCT.)
    #[allow(unused_assignments)]
    fn perform(&mut self, f: &mut dyn FunctionWithDerivative, starting_point: f64, inf_bound: f64, sup_bound: f64) {
        let mut change_direction = false;
        let mut sort;
        let mut is_new_sol;
        let mut good;
        let mut verif;
        let mut stop;
        const EPS_SQRT: f64 = 1.0e-16;
        const EPS: f64 = 1.0e-32;
        const EPS2: f64 = 1.0e-64;
        const PROGRES: f64 = 0.005;
        let (mut f2, mut gnr1);
        let mut previous_minimum;
        let mut dy;
        let mut old_f;
        let mut ambda;
        let mut ambda2;
        let mut oldgr;
        // To keep steps to a quarter of the domain.
        let a_sup_bound = sup_bound.min(2e100);
        let an_inf_bound = inf_bound.max(-2e100);
        let inv_length_max = 1.0 / ((a_sup_bound - an_inf_bound) / 4.0).max(1.0e-9);
        // `aConstraints`: whether the unknown is on a bound.
        let mut constrained;
        let mut descente_iter;

        self.done = false;
        self.sol = starting_point;
        self.kount = 0;

        // Start inside the bounds.
        if self.sol <= inf_bound {
            self.sol = inf_bound;
        } else if self.sol > sup_bound {
            self.sol = sup_bound;
        }

        // The first value of F and its gradient.
        match evaluate(f, self.sol) {
            Some(e) => (f2, gnr1) = self.take(e),
            None => {
                self.done = false;
                return;
            }
        }
        ambda2 = gnr1;
        // Save(0) is used only by the accelerating test at the end of the loop; on a
        // solution already, it must not force a second iteration.
        self.save[0] = f2.max(EPS_SQRT);
        let a_tol_func = epsilon(f2);

        if f2 <= EPS || gnr1 <= EPS2 {
            self.done = true;
            return;
        }

        self.kount = 1;
        while self.kount <= self.iteration_limit {
            let kount = self.kount;
            previous_minimum = f2;
            oldgr = gnr1;
            self.previous_solution = self.sol;
            self.sol_save = self.sol;

            (self.dh, dy) = search_direction(self.df, self.gh, self.ff, change_direction, inv_length_max, self.dh);
            if dy.abs() <= EPS {
                self.done = true;
                // (OCCT updates F here before GetStateNumber.)
                return;
            }
            if change_direction {
                ambda = ambda2 / dy.abs().sqrt();
                if ambda > 1.0 {
                    ambda = 1.0;
                }
            } else {
                ambda = 1.0;
                ambda2 = 0.5 * ambda / (self.dh * self.dh).sqrt();
            }
            self.sol += ambda * self.dh;
            (sort, constrained, is_new_sol) = bounds(inf_bound, sup_bound, self.tol, &mut self.sol, self.sol_save, &mut self.delta);

            self.dh_save = self.gh;
            if is_new_sol {
                match evaluate(f, self.sol) {
                    Some(e) => (f2, gnr1) = self.take(e),
                    None => {
                        self.done = false;
                        return;
                    }
                }
            }

            if f2 <= EPS || gnr1 <= EPS2 {
                self.done = true;
                return;
            }

            if sort || f2 / previous_minimum > PROGRES {
                dy = self.gh * self.dh;
                old_f = previous_minimum;
                stop = false;
                good = false;
                descente_iter = 0;
                let mut sortbis;

                // Standard processing, away from the bounds.
                if !sort {
                    // Not out: try to go further.
                    while f2 / previous_minimum > PROGRES && !stop {
                        if f2 < old_f && dy < 0.0 {
                            // Try to progress in this direction.
                            descente_iter += 1;
                            self.sol_save = self.sol;
                            old_f = f2;
                            self.sol += ambda * self.dh;
                            (stop, constrained, is_new_sol) = bounds(inf_bound, sup_bound, self.tol, &mut self.sol, self.sol_save, &mut self.delta);
                            ambda *= 1.7;
                        } else {
                            if f2 >= old_f || f2 >= previous_minimum {
                                good = false;
                                if descente_iter == 0 {
                                    // The first step fails: interpolate, and minimize if need be.
                                    descente_iter += 1;
                                    good = minimize_direction_interpolated(f, self.sol_save, &mut self.delta, old_f, f2, self.dh_save, self.gh, self.tol);
                                } else if change_direction || descente_iter > 1 || old_f > previous_minimum {
                                    // The progress was useful: minimize.
                                    descente_iter += 1;
                                    good = minimize_direction_brent(f, self.previous_solution, self.sol_save, self.sol, old_f, &mut self.delta, self.tol);
                                }
                                if !good {
                                    self.sol = self.sol_save;
                                    f2 = old_f;
                                } else {
                                    self.sol = self.sol_save + self.delta;
                                    (_, constrained, is_new_sol) = bounds(inf_bound, sup_bound, self.tol, &mut self.sol, self.sol_save, &mut self.delta);
                                }
                                // The point on the bound is rejected.
                                sort = false;
                            }
                            // And leave in any case.
                            stop = true;
                        }
                        self.dh_save = self.gh;
                        if is_new_sol {
                            match evaluate(f, self.sol) {
                                Some(e) => (f2, gnr1) = self.take(e),
                                None => {
                                    self.done = false;
                                    return;
                                }
                            }
                        }
                        dy = self.gh * self.dh;
                        if dy.abs() <= EPS {
                            if f2 > old_f {
                                self.sol = self.sol_save;
                            }
                            self.done = true;
                            return;
                        }
                        if descente_iter >= 100 {
                            stop = true;
                        }
                    }
                }

                // Processing on the bounds.
                if sort {
                    // Not to progress along the bound.
                    stop = f2 > 1.001 * old_f;
                    sortbis = sort;
                    descente_iter = 0;
                    while sortbis && (f2 < old_f || descente_iter == 0) && !stop {
                        descente_iter += 1;
                        // Try to progress along the bound.
                        self.sol_save = self.sol;
                        old_f = f2;
                        (self.dh, dy) = search_direction_constrained(self.df, self.gh, self.ff, constrained, change_direction, inv_length_max, self.dh);
                        if dy < -EPS {
                            // (Avoids useless computations and division by zero.)
                            if change_direction {
                                ambda = ambda2 / (-dy).sqrt();
                                if ambda > 1.0 {
                                    ambda = 1.0;
                                }
                            } else {
                                ambda = 1.0;
                                ambda2 = 0.5 * ambda / (self.dh * self.dh).sqrt();
                            }
                            self.sol += ambda * self.dh;
                            (sortbis, constrained, is_new_sol) = bounds(inf_bound, sup_bound, self.tol, &mut self.sol, self.sol_save, &mut self.delta);

                            self.dh_save = self.gh;
                            if is_new_sol {
                                match evaluate(f, self.sol) {
                                    Some(e) => (f2, gnr1) = self.take(e),
                                    None => {
                                        self.done = false;
                                        return;
                                    }
                                }
                            }
                            ambda2 = gnr1;
                        } else {
                            stop = true;
                        }

                        while f2 / previous_minimum > PROGRES && f2 < old_f && !stop {
                            descente_iter += 1;
                            if f2 < old_f && dy < 0.0 {
                                // Try to progress in this direction.
                                self.sol_save = self.sol;
                                old_f = f2;
                                self.sol += ambda * self.dh;
                                (sortbis, constrained, is_new_sol) = bounds(inf_bound, sup_bound, self.tol, &mut self.sol, self.sol_save, &mut self.delta);
                            }
                            self.dh_save = self.gh;
                            if is_new_sol {
                                match evaluate(f, self.sol) {
                                    Some(e) => (f2, gnr1) = self.take(e),
                                    None => {
                                        self.done = false;
                                        return;
                                    }
                                }
                            }
                            ambda2 = gnr1;
                            dy = self.gh * self.dh;
                            stop = dy >= 0.0 || descente_iter >= 10 || sortbis;
                        }
                        stop = dy >= 0.0 || descente_iter >= 10;
                    }
                    if (f2 / previous_minimum > PROGRES && f2 >= old_f) || f2 >= previous_minimum {
                        // Minimize by Brent.
                        descente_iter += 1;
                        good = minimize_direction_interpolated(f, self.sol_save, &mut self.delta, old_f, f2, self.dh_save, self.gh, self.tol);
                        if !good {
                            self.sol = self.sol_save;
                        } else {
                            self.sol = self.sol_save + self.delta;
                            (_, constrained, is_new_sol) = bounds(inf_bound, sup_bound, self.tol, &mut self.sol, self.sol_save, &mut self.delta);
                            if is_new_sol {
                                match evaluate(f, self.sol) {
                                    Some(e) => (f2, gnr1) = self.take(e),
                                    None => {
                                        self.done = false;
                                        return;
                                    }
                                }
                            }
                        }
                    }
                }
            }

            // The stopping tests.
            self.save[kount] = f2;
            // Is this the solution?
            if change_direction {
                // Gradient: avoid looping.
                verif = true;
            } else {
                verif = false;
                if kount > 1 {
                    // Speeds up the nearly quadratic cases.
                    if self.save[kount - 1] < 1.0e-4 * self.save[kount - 2] {
                        verif = true;
                    }
                } else {
                    // For the cases already solved.
                    verif = f2 < 1.0e-6 * self.save[0];
                }
            }
            if verif {
                self.delta = self.previous_solution - self.sol;
                if self.is_solution_reached() {
                    if previous_minimum < f2 {
                        self.sol = self.sol_save;
                    }
                    self.done = true;
                    return;
                }
            }

            // The progress: the current minimum against the previous one.
            if f2 - previous_minimum <= a_tol_func {
                if kount > 5 {
                    // Is the history good?
                    if f2 >= 0.95 * self.save[kount - 5] {
                        if !change_direction {
                            change_direction = true;
                        } else {
                            // Less than 5% gained: leave.
                            self.done = true;
                            return;
                        }
                    } else {
                        // It is: start again.
                        change_direction = false;
                    }
                } else {
                    // No history: go on.
                    change_direction = false;
                }
                // If Newton does not decrease the gradient enough, try the gradient method,
                // unless f decreases (with Newton the gradient of f can grow while f
                // decreases: Newton must then be kept).
                if gnr1 > 0.9 * oldgr && f2 > 0.5 * previous_minimum {
                    change_direction = true;
                }
                // Not changing strategy: check, if not done already.
                if !change_direction && !verif {
                    self.delta = self.previous_solution - self.sol;
                    if self.is_solution_reached() {
                        self.done = true;
                        return;
                    }
                }
            } else {
                // A regression.
                if !change_direction {
                    // Switch to the gradient.
                    change_direction = true;
                    self.sol = self.previous_solution;
                    match evaluate(f, self.sol) {
                        Some(e) => (f2, gnr1) = self.take(e),
                        None => {
                            self.done = false;
                            return;
                        }
                    }
                } else {
                    // No way out left.
                    return;
                }
            }
            self.kount += 1;
        }
    }
}
