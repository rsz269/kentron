//! Sequential least squares programming (SLSQP): the minimum of an objective with a gradient,
//! within lower and upper bounds on each variable, taken step by step exactly as nlopt 2.11.0's
//! `NLOPT_LD_SLSQP` takes it, so that where the minimum is not unique (a sketch whose
//! constraints leave it free to move) it ends where CadQuery's solver, which calls nlopt, ends.
//!
//! The method is Dieter Kraft's (ACM TOMS algorithm 733, 1994), in nlopt's C translation
//! (slsqp.c, by S. G. Johnson via f2c). Each iteration solves a quadratic model of the
//! objective, with a BFGS approximation of its Hessian kept as LDLᵀ factors, for a search
//! direction within the bounds, a linear least squares problem with inequalities (Lawson &
//! Hanson's LSEI, LSI, LDP and NNLS); then searches along the direction for a step that lowers
//! the objective enough (Armijo-like, by safeguarded quadratic interpolation), and updates the
//! factors by the damped BFGS formula. Only what CadQuery asks of nlopt is kept: no nonlinear
//! constraints (nlopt's m = p = 0), bounds that may be infinite, nlopt's own stopping tests
//! (stop.c) in place of Kraft's accuracy test (nlopt passes it acc = 0).
//!
//! To reproduce nlopt's iterates to the last bit, the port is a line-by-line translation: the
//! routines keep Kraft's names, arguments, arithmetic order and one-based indexing, and all their
//! real arrays live in one `memory` vector, addressed by index where the C code passes pointers
//! (the C code passes overlapping pieces of one work array to the same routine, which Rust
//! slices cannot express). The reverse-communication state machine of slsqpb (its `mode`) is
//! kept too, with the driver (slsqp.c: nlopt_slsqp) calling it until it reports an end.

/// When to stop, as nlopt's stopping settings (stop.c). CadQuery sets `ftol_rel` = `ftol_abs` = 0,
/// `xtol_rel` = 1e-9, `xtol_abs` = 1e-12 for every variable, and no evaluation limit.
#[derive(Debug, Clone)]
pub(crate) struct Stopping {
    /// Stop once an iteration changes the objective by less than `ftol_rel` times its size
    /// (the mean of the old and new values' magnitudes) or less than `ftol_abs`.
    pub ftol_rel: f64,
    pub ftol_abs: f64,
    /// Stop once an iteration moves the variables by less than `xtol_rel` times their size (both
    /// as L1 norms), or each variable by less than its `xtol_abs` (empty: no such test).
    pub xtol_rel: f64,
    pub xtol_abs: Vec<f64>,
    /// Stop after this many evaluations (0: no limit).
    pub max_evaluations: usize,
}

/// Why the minimization stopped: nlopt's result codes (nlopt.h: nlopt_result), those that
/// SLSQP without nonlinear constraints can give.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Result {
    Success = 1,
    FtolReached = 3,
    XtolReached = 4,
    MaxevalReached = 5,
    Failure = -1,
    /// A bound below its partner, or a start outside the bounds.
    InvalidArgs = -2,
    RoundoffLimited = -4,
}

#[derive(Debug, Clone)]
pub(crate) struct Outcome {
    /// The best point evaluated (not necessarily the last), and the objective there.
    pub x: Vec<f64>,
    pub minimum: f64,
    pub evaluations: usize,
    pub result: Result,
}

/// The minimum of `objective` from `x0`, each variable within `lower[i]..=upper[i]` (which may be
/// infinite), as nlopt's SLSQP finds it (optimize.c: nlopt_optimize, slsqp.c: nlopt_slsqp).
/// `objective(x, gradient)` returns f(x), filling `gradient` (of length n) when it is `Some`.
///
/// Like nlopt, it refuses (`InvalidArgs`, without evaluating) a start outside the bounds rather
/// than moving it within them; with no variables it evaluates the objective once.
pub(crate) fn minimize(x0: &[f64], lower: &[f64], upper: &[f64], mut objective: impl FnMut(&[f64], Option<&mut [f64]>) -> f64, stop: &Stopping) -> Outcome {
    let n = x0.len();
    assert!(lower.len() == n && upper.len() == n && (stop.xtol_abs.is_empty() || stop.xtol_abs.len() == n));
    // (optimize.c: nlopt_optimize_)
    if n == 0 {
        let minimum = objective(x0, None);
        return Outcome { x: Vec::new(), minimum, evaluations: 1, result: Result::Success };
    }
    for i in 0..n {
        if lower[i] > upper[i] || x0[i] < lower[i] || x0[i] > upper[i] {
            return Outcome { x: x0.to_vec(), minimum: f64::INFINITY, evaluations: 0, result: Result::InvalidArgs };
        }
    }
    nlopt_slsqp(x0, lower, upper, &mut objective, stop)
}

// ---------------------------------------------------------------------------------------------
// The driver (slsqp.c: nlopt_slsqp), for m = p = 0: the constraint values and normals are empty,
// every point is feasible, and slsqp is told m = meq = 0, la = 1.

fn nlopt_slsqp(x: &[f64], lb: &[f64], ub: &[f64], f: &mut dyn FnMut(&[f64], Option<&mut [f64]>) -> f64, stop: &Stopping) -> Outcome {
    let n = x.len();
    let ni = n as i32;
    let (mpi, pi) = (0_i32, 0_i32);
    let mpi1 = 1_i32;
    let (len_w, len_jw) = length_work(mpi, pi, ni);
    // The C code's work arrays, in its order (cgrad, c, grad, xcur, xprev, w), with the bounds
    // (the caller's arrays in C) placed before w so that w is last.
    let cgrad = 0;
    let c = cgrad + mpi1 as usize * (n + 1);
    let grad = c + mpi as usize;
    let xcur = grad + n + 1;
    let xprev = xcur + n;
    let xl = xprev + n;
    let xu = xl + n;
    let w = xu + n;
    let mut memory = vec![0.0; w + len_w as usize];
    let mut jw = vec![0_i32; len_jw as usize];
    memory[xcur..xcur + n].copy_from_slice(x);
    memory[xprev..xprev + n].copy_from_slice(x);
    memory[xl..xl + n].copy_from_slice(lb);
    memory[xu..xu + n].copy_from_slice(ub);
    let mut state = SlsqpbState::default();
    let mut best = x.to_vec();
    let mut minf = f64::INFINITY;
    let mut fcur = f64::INFINITY;
    let mut fprev = f64::INFINITY;
    let mut feasible = false;
    let mut infeasibility = f64::INFINITY;
    let feasible_cur = true;
    let infeasibility_cur = 0.0;
    let mut mode = 0_i32;
    let mut prev_mode = 0_i32;
    let mut acc = 0.0; // (nlopt does its own convergence tests)
    let mut iter = 0_i32; // (no iteration limit: nlopt counts evaluations instead)
    let mut evaluations = 0_usize;
    let mut ret = Result::Success;
    let mut want_grad = true;
    let mut first = true; // (the C code evaluates before calling slsqp the first time)
    'done: loop {
        let evaluate = if first {
            first = false;
            true
        } else {
            slsqp(&mut memory, &mut jw, mpi, pi, mpi1, ni, xcur, xl, xu, fcur, c, grad, cgrad, &mut acc, &mut iter, &mut mode, w, len_w, 0, len_jw, &mut state);
            match mode {
                // Objective and gradient (unless a mode −2 evaluation just gave them here).
                -1 => !(prev_mode == -2 && !want_grad),
                -2 | 1 => true,
                0 => break 'done,
                8 => {
                    // A positive directional derivative in the line search: nlopt's relaxed test
                    // (ten times the tolerances, except xtol_abs) against slsqpb's last iterate.
                    ret = Result::RoundoffLimited;
                    if feasible_cur {
                        let relaxed = Stopping { ftol_rel: stop.ftol_rel * 10.0, ftol_abs: stop.ftol_abs * 10.0, xtol_rel: stop.xtol_rel * 10.0, ..stop.clone() };
                        let x0 = state.x0.expect("slsqp ran");
                        if stop_ftol(&relaxed, fcur, state.f0) {
                            ret = Result::FtolReached;
                        } else if stop_x(&relaxed, &memory[xcur..xcur + n], &memory[x0..x0 + n]) {
                            ret = Result::XtolReached;
                        }
                    }
                    break 'done;
                }
                5 | 6 | 7 => {
                    ret = Result::RoundoffLimited;
                    break 'done;
                }
                4 | 3 | 9 => {
                    ret = Result::Failure;
                    break 'done;
                }
                _ => {
                    ret = Result::InvalidArgs;
                    break 'done;
                }
            }
        };
        if evaluate {
            if mode == -1 || mode == -2 || mode == 0 {
                want_grad = true;
            }
            let (head, tail) = memory.split_at_mut(xcur);
            let point = &tail[..n];
            fcur = if want_grad { f(point, Some(&mut head[grad..grad + n])) } else { f(point, None) };
            evaluations += 1;
            if fcur.is_finite() {
                want_grad = false;
            }
        }
        prev_mode = mode;
        // The best point so far.
        if fcur.is_finite() && ((fcur < minf && (feasible_cur || !feasible)) || (!feasible && infeasibility_cur < infeasibility)) {
            minf = fcur;
            feasible = feasible_cur;
            infeasibility = infeasibility_cur;
            best.copy_from_slice(&memory[xcur..xcur + n]);
        }
        // Mode −1 ends a line search: the only time to test for convergence (as Kraft's code).
        if mode == -1 {
            if !fprev.is_infinite() && feasible_cur {
                if stop_ftol(stop, fcur, fprev) {
                    ret = Result::FtolReached;
                } else if stop_x(stop, &memory[xcur..xcur + n], &memory[xprev..xprev + n]) {
                    ret = Result::XtolReached;
                }
            }
            fprev = fcur;
            memory.copy_within(xcur..xcur + n, xprev);
        }
        if stop.max_evaluations > 0 && evaluations >= stop.max_evaluations {
            ret = Result::MaxevalReached;
        }
        if ret != Result::Success {
            break;
        }
    }
    if minf.is_infinite() {
        // No finite value anywhere: the last point evaluated (or the one before, if its value is
        // infinite).
        if fcur.is_infinite() {
            minf = fprev;
            best.copy_from_slice(&memory[xprev..xprev + n]);
        } else {
            minf = fcur;
            best.copy_from_slice(&memory[xcur..xcur + n]);
        }
    }
    Outcome { x: best, minimum: minf, evaluations, result: ret }
}

/// Whether `vnew` is within the tolerances of `vold` (stop.c: relstop).
fn relstop(vold: f64, vnew: f64, reltol: f64, abstol: f64) -> bool {
    if vold.is_infinite() {
        return false;
    }
    (vnew - vold).abs() < abstol || (vnew - vold).abs() < reltol * (vnew.abs() + vold.abs()) * 0.5 || (reltol > 0.0 && vnew == vold)
}

/// (stop.c: nlopt_stop_ftol)
fn stop_ftol(s: &Stopping, f: f64, oldf: f64) -> bool {
    relstop(oldf, f, s.ftol_rel, s.ftol_abs)
}

/// (stop.c: nlopt_stop_x, without weights)
fn stop_x(s: &Stopping, x: &[f64], oldx: &[f64]) -> bool {
    let mut diff_norm = 0.0;
    for i in 0..x.len() {
        diff_norm += (x[i] - oldx[i]).abs();
    }
    let mut vector_norm = 0.0;
    for i in 0..x.len() {
        vector_norm += x[i].abs();
    }
    if diff_norm < s.xtol_rel * vector_norm {
        return true;
    }
    if s.xtol_abs.is_empty() {
        return false;
    }
    for i in 0..x.len() {
        if (x[i] - oldx[i]).abs() >= s.xtol_abs[i] {
            return false;
        }
    }
    true
}

/// The lengths of the real and integer work arrays (slsqp.c: length_work).
fn length_work(m: i32, meq: i32, n: i32) -> (i32, i32) {
    let n1 = n + 1;
    let mineq = m - meq + n1 + n1;
    let len_w = (3 * n1 + m) * (n1 + 1) + (n1 - meq + 1) * (mineq + 2) + 2 * mineq + (n1 + mineq) * (n1 - meq) + 2 * meq + n1 + (n + 1) * n / 2 + 2 * m + 3 * n + 3 * n1 + 1;
    (len_w, mineq)
}

// ---------------------------------------------------------------------------------------------
// Indexing. A C pointer p into the memory becomes an index; the f2c code's one-based access
// p[k] after `p -= offset` becomes memory[at(base, k)] with base = p − offset.

#[inline]
fn at(base: isize, k: i32) -> usize {
    (base + k as isize) as usize
}

#[inline]
fn base(pointer: usize, offset: i32) -> isize {
    pointer as isize - offset as isize
}

/// (slsqp.c: MIN2, MAX2, with their behaviour on ties and NaNs)
#[inline]
fn min2<T: PartialOrd>(a: T, b: T) -> T {
    if a <= b { a } else { b }
}
#[inline]
fn max2<T: PartialOrd>(a: T, b: T) -> T {
    if a >= b { a } else { b }
}

// BLAS level 1, as slsqp.c has them (zero-based, with increments).

/// y ← x (slsqp.c: dcopy___; an increment of 0 for x copies its first element everywhere).
fn dcopy(n: i32, mem: &mut [f64], dx: usize, incx: i32, dy: usize, incy: i32) {
    if n <= 0 {
        return;
    }
    if incx == 0 {
        let x = mem[dx];
        for i in 0..n as usize {
            mem[dy + i * incy as usize] = x;
        }
    } else {
        for i in 0..n as usize {
            mem[dy + i * incy as usize] = mem[dx + i * incx as usize];
        }
    }
}

/// y ← y + a x (slsqp.c: daxpy_sl__)
fn daxpy(n: i32, da: f64, mem: &mut [f64], dx: usize, incx: i32, dy: usize, incy: i32) {
    if n <= 0 || da == 0.0 {
        return;
    }
    for i in 0..n as usize {
        mem[dy + i * incy as usize] += da * mem[dx + i * incx as usize];
    }
}

/// x · y (slsqp.c: ddot_sl__)
fn ddot(n: i32, mem: &[f64], dx: usize, incx: i32, dy: usize, incy: i32) -> f64 {
    let mut sum = 0.0;
    if n <= 0 {
        return 0.0;
    }
    for i in 0..n as usize {
        sum += mem[dx + i * incx as usize] * mem[dy + i * incy as usize];
    }
    sum
}

/// |x|, scaled by the largest magnitude (slsqp.c: dnrm2___)
fn dnrm2(n: i32, mem: &[f64], dx: usize, incx: i32) -> f64 {
    let mut xmax = 0.0_f64;
    let mut sum = 0.0;
    for i in 0..n.max(0) as usize {
        let xabs = mem[dx + incx as usize * i].abs();
        if xmax < xabs {
            xmax = xabs;
        }
    }
    if xmax == 0.0 {
        return 0.0;
    }
    let scale = 1.0 / xmax;
    for i in 0..n.max(0) as usize {
        let xs = scale * mem[dx + incx as usize * i];
        sum += xs * xs;
    }
    xmax * sum.sqrt()
}

/// Applies a Givens rotation (slsqp.c: dsrot_).
fn dsrot(n: i32, mem: &mut [f64], dx: usize, incx: i32, dy: usize, incy: i32, c: f64, s: f64) {
    for i in 0..n.max(0) as usize {
        let (xi, yi) = (dx + incx as usize * i, dy + incy as usize * i);
        let (x, y) = (mem[xi], mem[yi]);
        mem[xi] = c * x + s * y;
        mem[yi] = c * y - s * x;
    }
}

/// Constructs a Givens rotation, returning (c, s) (slsqp.c: dsrotg_).
fn dsrotg(da: &mut f64, db: &mut f64) -> (f64, f64) {
    let (absa, absb) = (da.abs(), db.abs());
    let (roe, scale) = if absa > absb { (*da, absa) } else { (*db, absb) };
    if scale != 0.0 {
        let iscale = 1.0 / scale;
        let (tmpa, tmpb) = (*da * iscale, *db * iscale);
        let r = (if roe < 0.0 { -scale } else { scale }) * ((tmpa * tmpa) + (tmpb * tmpb)).sqrt();
        let c = *da / r;
        let s = *db / r;
        *da = r;
        if c != 0.0 && c.abs() <= s {
            *db = 1.0 / c;
        } else {
            *db = s;
        }
        (c, s)
    } else {
        *da = 0.0;
        *db = 0.0;
        (1.0, 0.0)
    }
}

/// x ← a x (slsqp.c: dscal_sl__)
fn dscal(n: i32, da: f64, mem: &mut [f64], dx: usize, incx: i32) {
    for i in 0..n.max(0) as usize {
        mem[dx + i * incx as usize] *= da;
    }
}

// ---------------------------------------------------------------------------------------------
// Lawson & Hanson's least squares routines, with Kraft's revisions (slsqp.c).

/// Constructs (mode 1) or applies (mode 2) the Householder transformation I + u uᵀ / b that zeros
/// the elements `l1..=m` of the pivot vector u (stride `iue`) against element `lpivot`, applying it
/// to `ncv` vectors of c (element stride `ice`, vector stride `icv`). (slsqp.c: h12_)
#[allow(clippy::too_many_arguments)]
fn h12(mode: i32, lpivot: i32, l1: i32, m: i32, mem: &mut [f64], u: usize, iue: i32, up: &mut f64, c: usize, ice: i32, icv: i32, ncv: i32) {
    let one = 1.0;
    let u_dim1 = iue;
    let u_b = base(u, 1 + u_dim1);
    let c_b = base(c, 1);
    macro_rules! u {
        ($j:expr) => {
            mem[at(u_b, ($j) * u_dim1 + 1)]
        };
    }
    macro_rules! c {
        ($i:expr) => {
            mem[at(c_b, $i)]
        };
    }
    if 0 >= lpivot || lpivot >= l1 || l1 > m {
        return;
    }
    let mut cl = u!(lpivot).abs();
    if mode != 2 {
        // Construct the transformation.
        for j in l1..=m {
            let sm = u!(j).abs();
            cl = max2(sm, cl);
        }
        if cl <= 0.0 {
            return;
        }
        let clinv = one / cl;
        let d = u!(lpivot) * clinv;
        let mut sm = d * d;
        for j in l1..=m {
            let d = u!(j) * clinv;
            sm += d * d;
        }
        cl *= sm.sqrt();
        if u!(lpivot) > 0.0 {
            cl = -cl;
        }
        *up = u!(lpivot) - cl;
        u!(lpivot) = cl;
    } else if cl <= 0.0 {
        return;
    }
    // Apply the transformation I + u uᵀ / b to c.
    if ncv <= 0 {
        return;
    }
    let mut b = *up * u!(lpivot);
    if b >= 0.0 {
        return;
    }
    b = one / b;
    let mut i2 = 1 - icv + ice * (lpivot - 1);
    let incr = ice * (l1 - lpivot);
    for _ in 1..=ncv {
        i2 += icv;
        let mut i3 = i2 + incr;
        let mut i4 = i3;
        let mut sm = c!(i2) * *up;
        for i in l1..=m {
            sm += c!(i3) * u!(i);
            i3 += ice;
        }
        if sm == 0.0 {
            continue;
        }
        sm *= b;
        c!(i2) += sm * *up;
        for i in l1..=m {
            c!(i4) += sm * u!(i);
            i4 += ice;
        }
    }
}

/// Nonnegative least squares: x ≥ 0 minimizing |A x − b| for the m × n matrix A (leading
/// dimension `mda`); A and b are overwritten by Q A and Q b, w gets the dual solution, z is
/// working space, `indx` (in `jw`) the sets P and Z. Mode 1: solved; 2: bad dimensions; 3: more
/// than 3n iterations. (slsqp.c: nnls_)
#[allow(clippy::too_many_arguments)]
fn nnls(mem: &mut [f64], jw: &mut [i32], a: usize, mda: i32, m: i32, n: i32, b: usize, x: usize, rnorm: &mut f64, w: usize, z: usize, indx: usize, mode: &mut i32) {
    let one = 1.0;
    let factor = 0.01;
    let z_b = base(z, 1);
    let b_b = base(b, 1);
    let indx_b = indx as isize - 1;
    let w_b = base(w, 1);
    let x_b = base(x, 1);
    let a_dim1 = mda;
    let a_b = base(a, 1 + a_dim1);
    macro_rules! z {
        ($i:expr) => {
            mem[at(z_b, $i)]
        };
    }
    macro_rules! w {
        ($i:expr) => {
            mem[at(w_b, $i)]
        };
    }
    macro_rules! x {
        ($i:expr) => {
            mem[at(x_b, $i)]
        };
    }
    macro_rules! a {
        ($i:expr) => {
            mem[at(a_b, $i)]
        };
    }
    macro_rules! indx {
        ($i:expr) => {
            jw[at(indx_b, $i)]
        };
    }
    let p = |b: isize, k: i32| at(b, k); // (the address &p[k])

    *mode = 2;
    if m <= 0 || n <= 0 {
        return;
    }
    *mode = 1;
    let mut iter = 0;
    let itmax = n * 3;
    // Step one (initialize).
    for i in 1..=n {
        indx!(i) = i;
    }
    let mut iz1 = 1;
    let iz2 = n;
    let mut nsetp = 0;
    let mut npp1 = 1;
    x!(1) = 0.0;
    dcopy(n, mem, p(x_b, 1), 0, p(x_b, 1), 1);
    let mut izmax = 0;
    let mut up = 0.0;
    let mut jj = 0;
    'l280: {
        // Step two (compute dual variables): entry of loop A.
        'l110: loop {
            if iz1 > iz2 || nsetp >= m {
                break 'l280;
            }
            for iz in iz1..=iz2 {
                let j = indx!(iz);
                w!(j) = ddot(m - nsetp, mem, p(a_b, npp1 + j * a_dim1), 1, p(b_b, npp1), 1);
            }
            // Step three (test dual variables).
            let (iz, j) = 'l160: loop {
                let mut wmax = 0.0;
                for iz in iz1..=iz2 {
                    let j = indx!(iz);
                    if !(w!(j) <= wmax) {
                        wmax = w!(j);
                        izmax = iz;
                    }
                }
                // Exit of loop A.
                if wmax <= 0.0 {
                    break 'l280;
                }
                let iz = izmax;
                let j = indx!(iz);
                // Step four (test index j for linear dependency).
                let asave = a!(npp1 + j * a_dim1);
                h12(1, npp1, npp1 + 1, m, mem, p(a_b, j * a_dim1 + 1), 1, &mut up, p(z_b, 1), 1, 1, 0);
                let unorm = dnrm2(nsetp, mem, p(a_b, j * a_dim1 + 1), 1);
                let t = factor * a!(npp1 + j * a_dim1).abs();
                let d = unorm + t;
                if !(d - unorm <= 0.0) {
                    dcopy(m, mem, p(b_b, 1), 1, p(z_b, 1), 1);
                    h12(2, npp1, npp1 + 1, m, mem, p(a_b, j * a_dim1 + 1), 1, &mut up, p(z_b, 1), 1, 1, 1);
                    if z!(npp1) / a!(npp1 + j * a_dim1) > 0.0 {
                        break 'l160 (iz, j);
                    }
                }
                a!(npp1 + j * a_dim1) = asave;
                w!(j) = 0.0;
            };
            // Step five (add column).
            dcopy(m, mem, p(z_b, 1), 1, p(b_b, 1), 1);
            indx!(iz) = indx!(iz1);
            indx!(iz1) = j;
            iz1 += 1;
            nsetp = npp1;
            npp1 += 1;
            for jz in iz1..=iz2 {
                jj = indx!(jz);
                h12(2, nsetp, npp1, m, mem, p(a_b, j * a_dim1 + 1), 1, &mut up, p(a_b, jj * a_dim1 + 1), 1, mda, 1);
            }
            let k = min2(npp1, mda);
            w!(j) = 0.0;
            dcopy(m - nsetp, mem, p(w_b, j), 0, p(a_b, k + j * a_dim1), 1);
            // Step six (solve least squares subproblem): entry of loop B.
            loop {
                for ip in (1..=nsetp).rev() {
                    if ip != nsetp {
                        let d = -z!(ip + 1);
                        daxpy(ip, d, mem, p(a_b, jj * a_dim1 + 1), 1, p(z_b, 1), 1);
                    }
                    jj = indx!(ip);
                    z!(ip) /= a!(ip + jj * a_dim1);
                }
                iter += 1;
                if iter > itmax {
                    *mode = 3;
                    break 'l280;
                }
                // Steps seven to ten (step length).
                let mut alpha = one;
                jj = 0;
                for ip in 1..=nsetp {
                    if z!(ip) > 0.0 {
                        continue;
                    }
                    let l = indx!(ip);
                    let t = -x!(l) / (z!(ip) - x!(l));
                    if alpha < t {
                        continue;
                    }
                    alpha = t;
                    jj = ip;
                }
                for ip in 1..=nsetp {
                    let l = indx!(ip);
                    x!(l) = (one - alpha) * x!(l) + alpha * z!(ip);
                }
                // Exit of loop B.
                if jj == 0 {
                    continue 'l110;
                }
                // Step eleven (delete column).
                let mut i = indx!(jj);
                loop {
                    x!(i) = 0.0;
                    jj += 1;
                    for j in jj..=nsetp {
                        let ii = indx!(j);
                        indx!(j - 1) = ii;
                        let (pa, pb) = (p(a_b, j - 1 + ii * a_dim1), p(a_b, j + ii * a_dim1));
                        let (mut da, mut db) = (mem[pa], mem[pb]);
                        let (c, s) = dsrotg(&mut da, &mut db);
                        mem[pa] = da;
                        mem[pb] = db;
                        let t = a!(j - 1 + ii * a_dim1);
                        dsrot(n, mem, p(a_b, j - 1 + a_dim1), mda, p(a_b, j + a_dim1), mda, c, s);
                        a!(j - 1 + ii * a_dim1) = t;
                        a!(j + ii * a_dim1) = 0.0;
                        dsrot(1, mem, p(b_b, j - 1), 1, p(b_b, j), 1, c, s);
                    }
                    npp1 = nsetp;
                    nsetp -= 1;
                    iz1 -= 1;
                    indx!(iz1) = i;
                    if nsetp <= 0 {
                        *mode = 3;
                        break 'l280;
                    }
                    let mut again = false;
                    for j in 1..=nsetp {
                        jj = j;
                        i = indx!(jj);
                        if x!(i) <= 0.0 {
                            again = true;
                            break;
                        }
                    }
                    if !again {
                        break;
                    }
                }
                dcopy(m, mem, p(b_b, 1), 1, p(z_b, 1), 1);
            }
        }
    }
    // Step twelve (solution).
    let k = min2(npp1, m);
    *rnorm = dnrm2(m - nsetp, mem, p(b_b, k), 1);
    if npp1 > m {
        w!(1) = 0.0;
        dcopy(n, mem, p(w_b, 1), 0, p(w_b, 1), 1);
    }
}

/// Least distance programming: x minimizing ½ xᵀx subject to G x ≥ h (G m × n, leading
/// dimension `mg`), by NNLS on the dual. w gets the Lagrange multipliers. Mode 1: solved; 2: bad
/// dimensions; 3: NNLS iteration limit; 4: incompatible constraints. (slsqp.c: ldp_)
#[allow(clippy::too_many_arguments)]
fn ldp(mem: &mut [f64], jw: &mut [i32], g: usize, mg: i32, m: i32, n: i32, h: usize, x: usize, xnorm: &mut f64, w: usize, indx: usize, mode: &mut i32) {
    let one = 1.0;
    let h_b = base(h, 1);
    let x_b = base(x, 1);
    let g_dim1 = mg;
    let g_b = base(g, 1 + g_dim1);
    let w_b = base(w, 1);
    let p = |b: isize, k: i32| at(b, k);
    macro_rules! w {
        ($i:expr) => {
            mem[at(w_b, $i)]
        };
    }
    macro_rules! x {
        ($i:expr) => {
            mem[at(x_b, $i)]
        };
    }
    *mode = 2;
    if n <= 0 {
        return;
    }
    // State the dual problem.
    *mode = 1;
    x!(1) = 0.0;
    dcopy(n, mem, p(x_b, 1), 0, p(x_b, 1), 1);
    *xnorm = 0.0;
    if m == 0 {
        return;
    }
    let mut iw = 0;
    for j in 1..=m {
        for i in 1..=n {
            iw += 1;
            w!(iw) = mem[at(g_b, j + i * g_dim1)];
        }
        iw += 1;
        w!(iw) = mem[at(h_b, j)];
    }
    let if_ = iw + 1;
    for _ in 1..=n {
        iw += 1;
        w!(iw) = 0.0;
    }
    w!(iw + 1) = one;
    let n1 = n + 1;
    let iz = iw + 2;
    let iy = iz + n1;
    let iwdual = iy + m;
    // Solve the dual problem.
    let mut rnorm = 0.0;
    nnls(mem, jw, p(w_b, 1), n1, n1, m, p(w_b, if_), p(w_b, iy), &mut rnorm, p(w_b, iwdual), p(w_b, iz), indx, mode);
    if *mode != 1 {
        return;
    }
    *mode = 4;
    if rnorm <= 0.0 {
        return;
    }
    // Compute the solution of the primal problem.
    let mut fac = one - ddot(m, mem, p(h_b, 1), 1, p(w_b, iy), 1);
    let d = one + fac;
    if d - one <= 0.0 {
        return;
    }
    *mode = 1;
    fac = one / fac;
    for j in 1..=n {
        x!(j) = fac * ddot(m, mem, p(g_b, j * g_dim1 + 1), 1, p(w_b, iy), 1);
    }
    *xnorm = dnrm2(n, mem, p(x_b, 1), 1);
    // Compute the Lagrange multipliers of the primal problem.
    w!(1) = 0.0;
    dcopy(m, mem, p(w_b, 1), 0, p(w_b, 1), 1);
    daxpy(m, fac, mem, p(w_b, iy), 1, p(w_b, 1), 1);
}

/// Least squares with inequalities: x minimizing |E x − f| subject to G x ≥ h, by QR of E and
/// LDP (Lawson & Hanson 23.5). w gets the multipliers. Mode 1: solved; 2–4 as LDP; 5: E not of
/// full rank. (slsqp.c: lsi_)
#[allow(clippy::too_many_arguments)]
fn lsi(mem: &mut [f64], jw: &mut [i32], e: usize, f: usize, g: usize, h: usize, le: i32, me: i32, lg: i32, mg: i32, n: i32, x: usize, xnorm: &mut f64, w: usize, jw_p: usize, mode: &mut i32) {
    let epmach = 2.22e-16;
    let one = 1.0;
    let f_b = base(f, 1);
    let h_b = base(h, 1);
    let x_b = base(x, 1);
    let g_dim1 = lg;
    let g_b = base(g, 1 + g_dim1);
    let e_dim1 = le;
    let e_b = base(e, 1 + e_dim1);
    let w_b = base(w, 1);
    let p = |b: isize, k: i32| at(b, k);
    macro_rules! e {
        ($i:expr) => {
            mem[at(e_b, $i)]
        };
    }
    macro_rules! g {
        ($i:expr) => {
            mem[at(g_b, $i)]
        };
    }
    macro_rules! x {
        ($i:expr) => {
            mem[at(x_b, $i)]
        };
    }
    let mut t = 0.0;
    // QR factors of E, applied to f.
    for i in 1..=n {
        let j = min2(i + 1, n);
        h12(1, i, i + 1, me, mem, p(e_b, i * e_dim1 + 1), 1, &mut t, p(e_b, j * e_dim1 + 1), 1, le, n - i);
        h12(2, i, i + 1, me, mem, p(e_b, i * e_dim1 + 1), 1, &mut t, p(f_b, 1), 1, 1, 1);
    }
    // Transform G and h to get a least distance problem.
    *mode = 5;
    for i in 1..=mg {
        for j in 1..=n {
            if e!(j + j * e_dim1).abs() < epmach {
                return;
            }
            let d = ddot(j - 1, mem, p(g_b, i + g_dim1), lg, p(e_b, j * e_dim1 + 1), 1);
            g!(i + j * g_dim1) = (g!(i + j * g_dim1) - d) / e!(j + j * e_dim1);
        }
        let d = ddot(n, mem, p(g_b, i + g_dim1), lg, p(f_b, 1), 1);
        mem[at(h_b, i)] -= d;
    }
    // Solve the least distance problem.
    ldp(mem, jw, g, lg, mg, n, p(h_b, 1), p(x_b, 1), xnorm, p(w_b, 1), jw_p, mode);
    if *mode != 1 {
        return;
    }
    // The solution of the original problem.
    daxpy(n, one, mem, p(f_b, 1), 1, p(x_b, 1), 1);
    for i in (1..=n).rev() {
        let j = min2(i + 1, n);
        let d = ddot(n - i, mem, p(e_b, i + j * e_dim1), le, p(x_b, j), 1);
        x!(i) = (x!(i) - d) / e!(i + i * e_dim1);
    }
    let j = min2(n + 1, me);
    let t = dnrm2(me - n, mem, p(f_b, j), 1);
    *xnorm = (*xnorm * *xnorm + t * t).sqrt();
}

/// Rank-deficient least squares A x = b by Householder triangularization with column pivoting
/// (A m × n, leading dimension `mda`; b, leading dimension `mdb`, overwritten by the solution;
/// `tau` the tolerance for the pseudorank `krank`). h and g are working space, `ip` (in `jw`) the
/// column interchanges. Only `nb` = 1 is used. (slsqp.c: hfti_)
#[allow(clippy::too_many_arguments)]
fn hfti(mem: &mut [f64], jw: &mut [i32], a: usize, mda: i32, m: i32, n: i32, b: usize, mdb: i32, nb: i32, tau: f64, krank: &mut i32, rnorm: &mut f64, h: usize, g: usize, ip: usize) {
    assert!(nb == 1);
    let factor = 0.001;
    let ip_b = ip as isize - 1;
    let g_b = base(g, 1);
    let h_b = base(h, 1);
    let a_dim1 = mda;
    let a_b = base(a, 1 + a_dim1);
    let b_dim1 = mdb;
    let b_offset = 1 + b_dim1;
    let b_b = base(b, b_offset);
    let p = |b: isize, k: i32| at(b, k);
    macro_rules! a {
        ($i:expr) => {
            mem[at(a_b, $i)]
        };
    }
    macro_rules! b {
        ($i:expr) => {
            mem[at(b_b, $i)]
        };
    }
    macro_rules! h {
        ($i:expr) => {
            mem[at(h_b, $i)]
        };
    }
    macro_rules! ip {
        ($i:expr) => {
            jw[at(ip_b, $i)]
        };
    }
    let mut k = 0;
    let ldiag = min2(m, n);
    'l270: {
        if ldiag <= 0 {
            break 'l270;
        }
        let mut hmax = 0.0;
        for j in 1..=ldiag {
            let mut lmax;
            let mut recompute = true;
            if j != 1 {
                // Update the squared column norms; recompute them if they have lost accuracy.
                lmax = j;
                for l in j..=n {
                    let d = a!(j - 1 + l * a_dim1);
                    h!(l) -= d * d;
                    if h!(l) > h!(lmax) {
                        lmax = l;
                    }
                }
                let d = hmax + factor * h!(lmax);
                if d - hmax > 0.0 {
                    recompute = false;
                }
            } else {
                lmax = j;
            }
            if recompute {
                lmax = j;
                for l in j..=n {
                    h!(l) = 0.0;
                    for i in j..=m {
                        let d = a!(i + l * a_dim1);
                        h!(l) += d * d;
                    }
                    if h!(l) > h!(lmax) {
                        lmax = l;
                    }
                }
                hmax = h!(lmax);
            }
            // Column interchanges if needed.
            ip!(j) = lmax;
            if ip!(j) != j {
                for i in 1..=m {
                    let tmp = a!(i + j * a_dim1);
                    a!(i + j * a_dim1) = a!(i + lmax * a_dim1);
                    a!(i + lmax * a_dim1) = tmp;
                }
                h!(lmax) = h!(j);
            }
            // The j-th transformation, applied to A and b.
            let i = min2(j + 1, n);
            let mut up = h!(j);
            h12(1, j, j + 1, m, mem, p(a_b, j * a_dim1 + 1), 1, &mut up, p(a_b, i * a_dim1 + 1), 1, mda, n - j);
            h12(2, j, j + 1, m, mem, p(a_b, j * a_dim1 + 1), 1, &mut up, p(b_b, b_offset), 1, mdb, nb);
            h!(j) = up;
        }
        // Determine the pseudorank.
        k = ldiag;
        for j in 1..=ldiag {
            if a!(j + j * a_dim1).abs() <= tau {
                k = j - 1;
                break;
            }
        }
        let kp1 = k + 1;
        // Norm of the residuals.
        for jb in 1..=nb {
            *rnorm = dnrm2(m - k, mem, p(b_b, kp1 + jb * b_dim1), 1);
        }
        if k <= 0 {
            for jb in 1..=nb {
                for i in 1..=n {
                    b!(i + jb * b_dim1) = 0.0;
                }
            }
            break 'l270;
        }
        if k != n {
            // Householder decomposition of the first k rows.
            for i in (1..=k).rev() {
                let mut up = mem[at(g_b, i)];
                h12(1, i, kp1, n, mem, p(a_b, i + a_dim1), mda, &mut up, p(a_b, 1 + a_dim1), mda, 1, i - 1);
                mem[at(g_b, i)] = up;
            }
        }
        for jb in 1..=nb {
            // Solve the k × k triangular system.
            for i in (1..=k).rev() {
                let j = min2(i + 1, n);
                let d = ddot(k - i, mem, p(a_b, i + j * a_dim1), mda, p(b_b, j + jb * b_dim1), 1);
                b!(i + jb * b_dim1) = (b!(i + jb * b_dim1) - d) / a!(i + i * a_dim1);
            }
            // Complete the solution vector.
            if k != n {
                for j in kp1..=n {
                    b!(j + jb * b_dim1) = 0.0;
                }
                for i in 1..=k {
                    let mut up = mem[at(g_b, i)];
                    h12(2, i, kp1, n, mem, p(a_b, i + a_dim1), mda, &mut up, p(b_b, jb * b_dim1 + 1), 1, mdb, 1);
                    mem[at(g_b, i)] = up;
                }
            }
            // Reorder the solution by the column interchanges.
            for j in (1..=ldiag).rev() {
                if ip!(j) == j {
                    continue;
                }
                let l = ip!(j);
                let tmp = b!(l + jb * b_dim1);
                b!(l + jb * b_dim1) = b!(j + jb * b_dim1);
                b!(j + jb * b_dim1) = tmp;
            }
        }
    }
    *krank = k;
}

/// Least squares with equalities and inequalities: x minimizing |E x − f| subject to C x = d,
/// G x ≥ h, by QR of C and a basis of its null space (Lawson & Hanson 23.6). w gets the
/// multipliers in its first mc + mg elements. Mode 1: solved; 2–5 as LSI; 6: C not of full
/// rank; 7: rank defect in HFTI. (slsqp.c: lsei_)
#[allow(clippy::too_many_arguments)]
fn lsei(mem: &mut [f64], jw: &mut [i32], c: usize, d: usize, e: usize, f: usize, g: usize, h: usize, lc: i32, mc: i32, le: i32, me: i32, lg: i32, mg: i32, n: i32, x: usize, xnrm: &mut f64, w: usize, jw_p: usize, mode: &mut i32) {
    let epmach = 2.22e-16;
    let d_b = base(d, 1);
    let f_b = base(f, 1);
    let h_b = base(h, 1);
    let x_b = base(x, 1);
    let g_dim1 = lg;
    let g_b = base(g, 1 + g_dim1);
    let e_dim1 = le;
    let e_b = base(e, 1 + e_dim1);
    let c_dim1 = lc;
    let c_b = base(c, 1 + c_dim1);
    let w_b = base(w, 1);
    let p = |b: isize, k: i32| at(b, k);
    macro_rules! c {
        ($i:expr) => {
            mem[at(c_b, $i)]
        };
    }
    macro_rules! d {
        ($i:expr) => {
            mem[at(d_b, $i)]
        };
    }
    macro_rules! f {
        ($i:expr) => {
            mem[at(f_b, $i)]
        };
    }
    macro_rules! w {
        ($i:expr) => {
            mem[at(w_b, $i)]
        };
    }
    macro_rules! x {
        ($i:expr) => {
            mem[at(x_b, $i)]
        };
    }
    *mode = 2;
    if mc > n {
        return;
    }
    let l = n - mc;
    let mc1 = mc + 1;
    let iw = (l + 1) * (mg + 2) + (mg << 1) + mc;
    let ie = iw + mc + 1;
    let if_ = ie + me * l;
    let ig = if_ + me;
    // Triangularize C and apply the factors to E and G.
    for i in 1..=mc {
        let j = min2(i + 1, lc);
        let mut up = w!(iw + i);
        h12(1, i, i + 1, n, mem, p(c_b, i + c_dim1), lc, &mut up, p(c_b, j + c_dim1), lc, 1, mc - i);
        h12(2, i, i + 1, n, mem, p(c_b, i + c_dim1), lc, &mut up, e, le, 1, me);
        h12(2, i, i + 1, n, mem, p(c_b, i + c_dim1), lc, &mut up, g, lg, 1, mg);
        w!(iw + i) = up;
    }
    // Solve C x = d and modify f.
    *mode = 6;
    for i in 1..=mc {
        if c!(i + i * c_dim1).abs() < epmach {
            return;
        }
        let s = ddot(i - 1, mem, p(c_b, i + c_dim1), lc, p(x_b, 1), 1);
        x!(i) = (d!(i) - s) / c!(i + i * c_dim1);
    }
    *mode = 1;
    w!(mc1) = 0.0;
    dcopy(mg, mem, p(w_b, mc1), 0, p(w_b, mc1), 1); // (nlopt's 2010 fix: mg, not mg − mc)
    'l50: {
        if mc == n {
            break 'l50;
        }
        for i in 1..=me {
            let s = ddot(mc, mem, p(e_b, i + e_dim1), le, p(x_b, 1), 1);
            w!(if_ - 1 + i) = f!(i) - s;
        }
        // Store the transformed E and G.
        for i in 1..=me {
            dcopy(l, mem, p(e_b, i + mc1 * e_dim1), le, p(w_b, ie - 1 + i), me);
        }
        for i in 1..=mg {
            dcopy(l, mem, p(g_b, i + mc1 * g_dim1), lg, p(w_b, ig - 1 + i), mg);
        }
        if mg <= 0 {
            // Solve the least squares problem without inequality constraints.
            *mode = 7;
            let k = max2(le, n);
            let t = epmach.sqrt();
            let mut krank = 0;
            hfti(mem, jw, p(w_b, ie), me, me, l, p(w_b, if_), k, 1, t, &mut krank, xnrm, p(w_b, 1), p(w_b, l + 1), jw_p);
            dcopy(l, mem, p(w_b, if_), 1, p(x_b, mc1), 1);
            if krank != l {
                return;
            }
            *mode = 1;
            break 'l50;
        }
        // Modify h and solve the inequality constrained least squares problem.
        for i in 1..=mg {
            let s = ddot(mc, mem, p(g_b, i + g_dim1), lg, p(x_b, 1), 1);
            mem[at(h_b, i)] -= s;
        }
        lsi(mem, jw, p(w_b, ie), p(w_b, if_), p(w_b, ig), p(h_b, 1), me, me, mg, mg, l, p(x_b, mc1), xnrm, p(w_b, mc1), jw_p, mode);
        if mc == 0 {
            return;
        }
        let t = dnrm2(mc, mem, p(x_b, 1), 1);
        *xnrm = (*xnrm * *xnrm + t * t).sqrt();
        if *mode != 1 {
            return;
        }
    }
    // The solution of the original problem, and the Lagrange multipliers.
    for i in 1..=me {
        let s = ddot(n, mem, p(e_b, i + e_dim1), le, p(x_b, 1), 1);
        f!(i) = s - f!(i);
    }
    for i in 1..=mc {
        let s = ddot(me, mem, p(e_b, i * e_dim1 + 1), 1, p(f_b, 1), 1) - ddot(mg, mem, p(g_b, i * g_dim1 + 1), 1, p(w_b, mc1), 1);
        d!(i) = s;
    }
    for i in (1..=mc).rev() {
        let mut up = w!(iw + i);
        h12(2, i, i + 1, n, mem, p(c_b, i + c_dim1), lc, &mut up, p(x_b, 1), 1, 1, 1);
        w!(iw + i) = up;
    }
    for i in (1..=mc).rev() {
        let j = min2(i + 1, lc);
        let s = ddot(mc - i, mem, p(c_b, j + i * c_dim1), 1, p(w_b, j), 1);
        w!(i) = (d!(i) - s) / c!(i + i * c_dim1);
    }
}

/// The quadratic subproblem as a least squares problem: x minimizing |E x − f| with E = D^½ Lᵀ
/// and f = −D^−½ L⁻¹ g, from the LDLᵀ factors l (n(n+1)/2, columnwise, D on the diagonal) and the
/// gradient g, subject to A x − b = 0 (the first meq rows), A x − b ≥ 0 (the other m − meq) and
/// xl ≤ x ≤ xu (infinite bounds left out, nlopt's 2010 change). y gets the multipliers. When
/// `nl` is not n(n+1)/2 + 1 the problem is the augmented one of an inconsistent linearization,
/// whose last variable has its diagonal in l[nl]. Mode as LSEI. (slsqp.c: lsq_)
#[allow(clippy::too_many_arguments)]
fn lsq(mem: &mut [f64], jw: &mut [i32], m: i32, meq: i32, n: i32, nl: i32, la: i32, l: usize, g: usize, a: usize, b: usize, xl: usize, xu: usize, x: usize, y: usize, w: usize, jw_p: usize, mode: &mut i32) {
    let one = 1.0;
    let y_b = base(y, 1);
    let x_b = base(x, 1);
    let xu_b = base(xu, 1);
    let xl_b = base(xl, 1);
    let g_b = base(g, 1);
    let l_b = base(l, 1);
    let b_b = base(b, 1);
    let a_dim1 = la;
    let a_b = base(a, 1 + a_dim1);
    let w_b = base(w, 1);
    let p = |b: isize, k: i32| at(b, k);
    macro_rules! w {
        ($i:expr) => {
            mem[at(w_b, $i)]
        };
    }
    macro_rules! l {
        ($i:expr) => {
            mem[at(l_b, $i)]
        };
    }
    macro_rules! xl {
        ($i:expr) => {
            mem[at(xl_b, $i)]
        };
    }
    macro_rules! xu {
        ($i:expr) => {
            mem[at(xu_b, $i)]
        };
    }
    macro_rules! x {
        ($i:expr) => {
            mem[at(x_b, $i)]
        };
    }
    let n1 = n + 1;
    let mineq = m - meq;
    let m1 = mineq + n + n;
    // Whether to solve the problem with an inconsistent linearization (n2 = 1) or not (n2 = 0).
    let mut n2 = n1 * n / 2 + 1;
    n2 = if n2 == nl { 0 } else { 1 };
    let n3 = n - n2;
    // Recover the matrix E and the vector f from l and g.
    let mut i2 = 1;
    let mut i3 = 1;
    let mut i4 = 1;
    let ie = 1;
    let if_ = n * n + 1;
    for i in 1..=n3 {
        let i1 = n1 - i;
        let diag = l!(i2).sqrt();
        w!(i3) = 0.0;
        dcopy(i1, mem, p(w_b, i3), 0, p(w_b, i3), 1);
        dcopy(i1 - n2, mem, p(l_b, i2), 1, p(w_b, i3), n);
        dscal(i1 - n2, diag, mem, p(w_b, i3), n);
        w!(i3) = diag;
        let s = ddot(i - 1, mem, p(w_b, i4), 1, p(w_b, if_), 1);
        w!(if_ - 1 + i) = (mem[at(g_b, i)] - s) / diag;
        i2 = i2 + i1 - n2;
        i3 += n1;
        i4 += n;
    }
    if n2 == 1 {
        w!(i3) = l!(nl);
        w!(i4) = 0.0;
        dcopy(n3, mem, p(w_b, i4), 0, p(w_b, i4), 1);
        w!(if_ - 1 + n) = 0.0;
    }
    dscal(n, -one, mem, p(w_b, if_), 1);
    let ic = if_ + n;
    let id = ic + meq * n;
    if meq > 0 {
        // Recover the matrix C from the upper part of A.
        for i in 1..=meq {
            dcopy(n, mem, p(a_b, i + a_dim1), la, p(w_b, ic - 1 + i), meq);
        }
        // Recover the vector d from the upper part of b.
        dcopy(meq, mem, p(b_b, 1), 1, p(w_b, id), 1);
        dscal(meq, -one, mem, p(w_b, id), 1);
    }
    let ig = id + meq;
    if mineq > 0 {
        // Recover the matrix G from the lower part of A.
        for i in 1..=mineq {
            dcopy(n, mem, p(a_b, meq + i + a_dim1), la, p(w_b, ig - 1 + i), m1);
        }
    }
    // Augment G by +I and −I (rows for infinite bounds left zero).
    let ip = ig + mineq;
    for i in 1..=n {
        w!(ip - 1 + i) = 0.0;
        dcopy(n, mem, p(w_b, ip - 1 + i), 0, p(w_b, ip - 1 + i), m1);
    }
    let stride = m1 + 1;
    for i in 1..=n {
        if !xl!(i).is_infinite() {
            w!((ip - stride) + i * stride) = 1.0;
        }
    }
    let im = ip + n;
    for i in 1..=n {
        w!(im - 1 + i) = 0.0;
        dcopy(n, mem, p(w_b, im - 1 + i), 0, p(w_b, im - 1 + i), m1);
    }
    for i in 1..=n {
        if !xu!(i).is_infinite() {
            w!((im - stride) + i * stride) = -1.0;
        }
    }
    let ih = ig + m1 * n;
    if mineq > 0 {
        // Recover h from the lower part of b.
        dcopy(mineq, mem, p(b_b, meq + 1), 1, p(w_b, ih), 1);
        dscal(mineq, -one, mem, p(w_b, ih), 1);
    }
    // Augment h by xl and −xu.
    let il = ih + mineq;
    let iu = il + n;
    for i in 1..=n {
        w!((il - 1) + i) = if xl!(i).is_infinite() { 0.0 } else { xl!(i) };
        w!((iu - 1) + i) = if xu!(i).is_infinite() { 0.0 } else { -xu!(i) };
    }
    let iw = iu + n;
    let mut xnorm = 0.0;
    lsei(mem, jw, p(w_b, ic), p(w_b, id), p(w_b, ie), p(w_b, if_), p(w_b, ig), p(w_b, ih), max2(1, meq), meq, n, n, m1, m1, n, p(x_b, 1), &mut xnorm, p(w_b, iw), jw_p, mode);
    if *mode == 1 {
        // Restore the Lagrange multipliers.
        dcopy(m, mem, p(w_b, iw), 1, p(y_b, 1), 1);
        dcopy(n3, mem, p(w_b, iw + m), 1, p(y_b, m + 1), 1);
        dcopy(n3, mem, p(w_b, iw + m + n), 1, p(y_b, m + n3 + 1), 1);
        // Within the bounds despite roundoff (nlopt's 2010 change: nlopt keeps bounds strictly).
        for i in 1..=n {
            if x!(i) < xl!(i) {
                x!(i) = xl!(i);
            } else if x!(i) > xu!(i) {
                x!(i) = xu!(i);
            }
        }
    }
}

/// Updates the LDLᵀ factors a (n(n+1)/2, columnwise) by the rank-one matrix σ z zᵀ (Fletcher &
/// Powell's method); z is overwritten, w is working space when σ < 0. (slsqp.c: ldl_)
fn ldl(mem: &mut [f64], n: i32, a: usize, z: usize, sigma: f64, w: usize) {
    let one = 1.0;
    let four = 4.0;
    let epmach = 2.22e-16;
    let w_b = base(w, 1);
    let z_b = base(z, 1);
    let a_b = base(a, 1);
    macro_rules! w {
        ($i:expr) => {
            mem[at(w_b, $i)]
        };
    }
    macro_rules! z {
        ($i:expr) => {
            mem[at(z_b, $i)]
        };
    }
    macro_rules! a {
        ($i:expr) => {
            mem[at(a_b, $i)]
        };
    }
    if sigma == 0.0 {
        return;
    }
    let mut ij = 1;
    let mut t = one / sigma;
    if !(sigma > 0.0) {
        // Prepare the negative update.
        for i in 1..=n {
            w!(i) = z!(i);
        }
        for i in 1..=n {
            let v = w!(i);
            t += v * v / a!(ij);
            for j in i + 1..=n {
                ij += 1;
                w!(j) -= v * a!(ij);
            }
            ij += 1;
        }
        if t >= 0.0 {
            t = epmach / sigma;
        }
        for i in 1..=n {
            let j = n + 1 - i;
            ij -= i;
            let u = w!(j);
            w!(j) = t;
            t -= u * u / a!(ij);
        }
    }
    // Here the updating begins.
    for i in 1..=n {
        let v = z!(i);
        let delta = v / a!(ij);
        let tp = if sigma < 0.0 { w!(i) } else { t + delta * v };
        let alpha = tp / t;
        a!(ij) = alpha * a!(ij);
        if i == n {
            return;
        }
        let beta = delta / tp;
        if alpha > four {
            let gamma = t / tp;
            for j in i + 1..=n {
                ij += 1;
                let u = a!(ij);
                a!(ij) = gamma * u + beta * z!(j);
                z!(j) -= v * u;
            }
        } else {
            for j in i + 1..=n {
                ij += 1;
                z!(j) -= v * a!(ij);
                a!(ij) += beta * z!(j);
            }
        }
        ij += 1;
        t = tp;
    }
}

// ---------------------------------------------------------------------------------------------
// The SQP iteration (slsqp.c: slsqpb_, slsqp).

/// What slsqpb keeps from one call to the next (slsqp.c: slsqpb_state). `x0` is the index of its
/// last iterate in the memory, set once slsqp has run.
#[derive(Debug, Clone, Copy, Default)]
struct SlsqpbState {
    t: f64,
    f0: f64,
    h1: f64,
    h2: f64,
    h3: f64,
    h4: f64,
    n1: i32,
    n2: i32,
    n3: i32,
    t0: f64,
    gs: f64,
    tol: f64,
    line: i32,
    alpha: f64,
    iexact: i32,
    incons: i32,
    ireset: i32,
    itermx: i32,
    x0: Option<usize>,
}

/// slsqpb's labels (the targets of its gotos).
#[derive(Clone, Copy)]
enum Label {
    Initialize,       // L100
    ResetHessian,     // L110
    Iterate,          // L130
    LineSearchStep,   // L190
    LineSearchTest,   // L200
    ExactLineSearch,  // L210
    Evaluated,        // L220
    Converged,        // L240
    RelaxedConverged, // L255
    UpdateHessian,    // L260
}

/// The body of SLSQP, called again after each evaluation it asks for (`mode` −2/−1: the objective
/// and gradient at x; 1: the objective at x) until it reports an end (0: converged, 2–9:
/// failures, see slsqp). (slsqp.c: slsqpb_)
#[allow(clippy::too_many_arguments)]
fn slsqpb(mem: &mut [f64], jw: &mut [i32], m: i32, meq: i32, la: i32, n: i32, x: usize, xl: usize, xu: usize, f: f64, c: usize, g: usize, a: usize, acc: &mut f64, iter: &mut i32, mode: &mut i32, r: usize, l: usize, x0: usize, mu: usize, s: usize, u: usize, v: usize, w: usize, iw: usize, state: &mut SlsqpbState) {
    let one = 1.0;
    let alfmin = 0.1;
    let hun = 100.0;
    let ten = 10.0;
    let two = 2.0;
    let SlsqpbState { mut t, mut f0, mut h1, mut h2, mut h3, mut h4, mut n1, mut n2, mut n3, mut t0, mut gs, mut tol, mut line, mut alpha, mut iexact, mut incons, mut ireset, mut itermx, .. } = *state;
    let mu_b = base(mu, 1);
    let c_b = base(c, 1);
    let v_b = base(v, 1);
    let u_b = base(u, 1);
    let s_b = base(s, 1);
    let x0_b = base(x0, 1);
    let l_b = base(l, 1);
    let r_b = base(r, 1);
    let a_dim1 = la;
    let a_b = base(a, 1 + a_dim1);
    let g_b = base(g, 1);
    let xu_b = base(xu, 1);
    let xl_b = base(xl, 1);
    let x_b = base(x, 1);
    let w_b = base(w, 1);
    let p = |b: isize, k: i32| at(b, k);
    macro_rules! mu {
        ($i:expr) => {
            mem[at(mu_b, $i)]
        };
    }
    macro_rules! c {
        ($i:expr) => {
            mem[at(c_b, $i)]
        };
    }
    macro_rules! v {
        ($i:expr) => {
            mem[at(v_b, $i)]
        };
    }
    macro_rules! u {
        ($i:expr) => {
            mem[at(u_b, $i)]
        };
    }
    macro_rules! s {
        ($i:expr) => {
            mem[at(s_b, $i)]
        };
    }
    macro_rules! l {
        ($i:expr) => {
            mem[at(l_b, $i)]
        };
    }
    macro_rules! r {
        ($i:expr) => {
            mem[at(r_b, $i)]
        };
    }
    macro_rules! a {
        ($i:expr) => {
            mem[at(a_b, $i)]
        };
    }
    macro_rules! g {
        ($i:expr) => {
            mem[at(g_b, $i)]
        };
    }
    macro_rules! x {
        ($i:expr) => {
            mem[at(x_b, $i)]
        };
    }
    macro_rules! xl {
        ($i:expr) => {
            mem[at(xl_b, $i)]
        };
    }
    macro_rules! xu {
        ($i:expr) => {
            mem[at(xu_b, $i)]
        };
    }
    let mut label = match *mode {
        -1 => Label::UpdateHessian,
        0 => Label::Initialize,
        _ => Label::Evaluated,
    };
    'l330: loop {
        match label {
            Label::Initialize => {
                itermx = *iter;
                iexact = if *acc >= 0.0 { 0 } else { 1 };
                *acc = acc.abs();
                tol = ten * *acc;
                *iter = 0;
                ireset = 0;
                n1 = n + 1;
                n2 = n1 * n / 2;
                n3 = n2 + 1;
                s!(1) = 0.0;
                mu!(1) = 0.0;
                dcopy(n, mem, p(s_b, 1), 0, p(s_b, 1), 1);
                dcopy(m, mem, p(mu_b, 1), 0, p(mu_b, 1), 1);
                label = Label::ResetHessian;
            }
            Label::ResetHessian => {
                ireset += 1;
                if ireset > 5 {
                    label = Label::RelaxedConverged;
                    continue;
                }
                l!(1) = 0.0;
                dcopy(n2, mem, p(l_b, 1), 0, p(l_b, 1), 1);
                let mut j = 1;
                for i in 1..=n {
                    l!(j) = one;
                    j = j + n1 - i;
                }
                label = Label::Iterate;
            }
            Label::Iterate => {
                // Main iteration: search direction, step length, LDLᵀ update.
                *iter += 1;
                *mode = 9;
                if *iter > itermx && itermx > 0 {
                    // (nlopt's 2010 change: no limit when itermx ≤ 0)
                    break 'l330;
                }
                // The search direction as the solution of the QP subproblem.
                dcopy(n, mem, p(xl_b, 1), 1, p(u_b, 1), 1);
                dcopy(n, mem, p(xu_b, 1), 1, p(v_b, 1), 1);
                daxpy(n, -one, mem, p(x_b, 1), 1, p(u_b, 1), 1);
                daxpy(n, -one, mem, p(x_b, 1), 1, p(v_b, 1), 1);
                h4 = one;
                lsq(mem, jw, m, meq, n, n3, la, p(l_b, 1), p(g_b, 1), a, p(c_b, 1), p(u_b, 1), p(v_b, 1), p(s_b, 1), p(r_b, 1), p(w_b, 1), iw, mode);
                // The augmented problem for an inconsistent linearization.
                if *mode == 6 && n == meq {
                    *mode = 4;
                }
                if *mode == 4 {
                    for j in 1..=m {
                        if j <= meq {
                            a!(j + n1 * a_dim1) = -c!(j);
                        } else {
                            a!(j + n1 * a_dim1) = max2(-c!(j), 0.0);
                        }
                    }
                    s!(1) = 0.0;
                    dcopy(n, mem, p(s_b, 1), 0, p(s_b, 1), 1);
                    h3 = 0.0;
                    g!(n1) = 0.0;
                    l!(n3) = hun;
                    s!(n1) = one;
                    u!(n1) = 0.0;
                    v!(n1) = one;
                    incons = 0;
                    loop {
                        lsq(mem, jw, m, meq, n1, n3, la, p(l_b, 1), p(g_b, 1), a, p(c_b, 1), p(u_b, 1), p(v_b, 1), p(s_b, 1), p(r_b, 1), p(w_b, 1), iw, mode);
                        h4 = one - s!(n1);
                        if *mode == 4 {
                            l!(n3) = ten * l!(n3);
                            incons += 1;
                            if incons > 5 {
                                break 'l330;
                            }
                            continue;
                        } else if *mode != 1 {
                            break 'l330;
                        }
                        break;
                    }
                } else if *mode != 1 {
                    break 'l330;
                }
                // Update the multipliers for the L1 test.
                for i in 1..=n {
                    v!(i) = g!(i) - ddot(m, mem, p(a_b, i * a_dim1 + 1), 1, p(r_b, 1), 1);
                }
                f0 = f;
                dcopy(n, mem, p(x_b, 1), 1, p(x0_b, 1), 1);
                gs = ddot(n, mem, p(g_b, 1), 1, p(s_b, 1), 1);
                h1 = gs.abs();
                h2 = 0.0;
                for j in 1..=m {
                    h3 = if j <= meq { c!(j) } else { 0.0 };
                    h2 += max2(-c!(j), h3);
                    h3 = r!(j).abs();
                    mu!(j) = max2(h3, (mu!(j) + h3) / two);
                    h1 += h3 * c!(j).abs();
                }
                // Check convergence.
                *mode = 0;
                if h1 < *acc && h2 < *acc {
                    break 'l330;
                }
                h1 = 0.0;
                for j in 1..=m {
                    h3 = if j <= meq { c!(j) } else { 0.0 };
                    h1 += mu!(j) * max2(-c!(j), h3);
                }
                t0 = f + h1;
                h3 = gs - h1 * h4;
                *mode = 8;
                if h3 >= 0.0 {
                    label = Label::ResetHessian;
                    continue;
                }
                // Line search with an L1 test function.
                line = 0;
                alpha = one;
                label = if iexact == 1 { Label::ExactLineSearch } else { Label::LineSearchStep };
            }
            Label::LineSearchStep => {
                // Inexact line search.
                line += 1;
                h3 = alpha * h3;
                dscal(n, alpha, mem, p(s_b, 1), 1);
                dcopy(n, mem, p(x0_b, 1), 1, p(x_b, 1), 1);
                daxpy(n, one, mem, p(s_b, 1), 1, p(x_b, 1), 1);
                // (nlopt's 2010 change: roundoff must not push x past the bounds)
                for i in 1..=n {
                    if x!(i) < xl!(i) {
                        x!(i) = xl!(i);
                    } else if x!(i) > xu!(i) {
                        x!(i) = xu!(i);
                    }
                }
                // (nlopt's mode −2: the first step asks for the gradient with the objective, as it
                // is usually accepted, saving the mode −1 call after it)
                *mode = if line == 1 { -2 } else { 1 };
                break 'l330;
            }
            Label::LineSearchTest => {
                if h1.is_finite() {
                    if h1 <= h3 / ten || line > 10 {
                        label = Label::Converged;
                        continue;
                    }
                    alpha = max2(h3 / (two * (h3 - h1)), alfmin);
                } else {
                    alpha = max2(alpha * 0.5, alfmin);
                }
                label = Label::LineSearchStep;
            }
            Label::ExactLineSearch => {
                // (nlopt has no exact line search: mode 9, a failure, without saving the state)
                *mode = 9;
                return;
            }
            Label::Evaluated => {
                // The functions at the current x.
                t = f;
                for j in 1..=m {
                    h1 = if j <= meq { c!(j) } else { 0.0 };
                    t += mu!(j) * max2(-c!(j), h1);
                }
                h1 = t - t0;
                label = if iexact == 0 { Label::LineSearchTest } else { Label::ExactLineSearch };
            }
            Label::Converged => {
                h3 = 0.0;
                for j in 1..=m {
                    h1 = if j <= meq { c!(j) } else { 0.0 };
                    h3 += max2(-c!(j), h1);
                }
                *mode = if ((f - f0).abs() < *acc || dnrm2(n, mem, p(s_b, 1), 1) < *acc) && h3 < *acc { 0 } else { -1 };
                break 'l330;
            }
            Label::RelaxedConverged => {
                // Relaxed convergence test after a positive directional derivative.
                *mode = if ((f - f0).abs() < tol || dnrm2(n, mem, p(s_b, 1), 1) < tol) && h3 < tol { 0 } else { 8 };
                break 'l330;
            }
            Label::UpdateHessian => {
                // The gradient at the current x: update the factors by the modified BFGS formula.
                for i in 1..=n {
                    u!(i) = g!(i) - ddot(m, mem, p(a_b, i * a_dim1 + 1), 1, p(r_b, 1), 1) - v!(i);
                }
                // Lᵀ s
                let mut k = 0;
                for i in 1..=n {
                    h1 = 0.0;
                    k += 1;
                    for j in i + 1..=n {
                        k += 1;
                        h1 += l!(k) * s!(j);
                    }
                    v!(i) = s!(i) + h1;
                }
                // D Lᵀ s
                k = 1;
                for i in 1..=n {
                    v!(i) = l!(k) * v!(i);
                    k = k + n1 - i;
                }
                // L D Lᵀ s
                for i in (1..=n).rev() {
                    h1 = 0.0;
                    k = i;
                    for j in 1..i {
                        h1 += l!(k) * v!(j);
                        k = k + n - j;
                    }
                    v!(i) += h1;
                }
                h1 = ddot(n, mem, p(s_b, 1), 1, p(u_b, 1), 1);
                h2 = ddot(n, mem, p(s_b, 1), 1, p(v_b, 1), 1);
                h3 = h2 * 0.2;
                if h1 < h3 {
                    h4 = (h2 - h3) / (h2 - h1);
                    h1 = h3;
                    dscal(n, h4, mem, p(u_b, 1), 1);
                    daxpy(n, one - h4, mem, p(v_b, 1), 1, p(u_b, 1), 1);
                }
                ldl(mem, n, p(l_b, 1), p(u_b, 1), one / h1, p(v_b, 1));
                ldl(mem, n, p(l_b, 1), p(v_b, 1), -one / h2, p(u_b, 1));
                label = Label::Iterate;
            }
        }
    }
    *state = SlsqpbState { t, f0, h1, h2, h3, h4, n1, n2, n3, t0, gs, tol, line, alpha, iexact, incons, ireset, itermx, x0: state.x0 };
}

/// The head of SLSQP: checks the work array's length and lays out slsqpb's arrays in it.
/// `mode` 2: more equality constraints than variables; 3: LSQ iteration limit; 4: inconsistent
/// constraints; 5, 6: singular E, C in LSQ; 7: rank defect in HFTI; 8: positive directional
/// derivative in the line search; 9: iteration limit; ≥ 10: work arrays too small.
/// (slsqp.c: slsqp)
#[allow(clippy::too_many_arguments)]
fn slsqp(mem: &mut [f64], jw: &mut [i32], m: i32, meq: i32, la: i32, n: i32, x: usize, xl: usize, xu: usize, f: f64, c: usize, g: usize, a: usize, acc: &mut f64, iter: &mut i32, mode: &mut i32, w: usize, l_w: i32, jw_p: usize, l_jw: i32, state: &mut SlsqpbState) {
    let n1 = n + 1;
    let mineq = m - meq + n1 + n1;
    let il = (n1 * 3 + m) * (n1 + 1) + (n1 - meq + 1) * (mineq + 2) + (mineq << 1) + (n1 + mineq) * (n1 - meq) + (meq << 1) + n1 * n / 2 + (m << 1) + n * 3 + (n1 << 2) + 1;
    let im = max2(mineq, n1 - meq);
    if l_w < il || l_jw < im {
        *mode = max2(10, il) * 1000;
        *mode += max2(10, im);
        return;
    }
    // The initial addresses in w (one-based).
    let im = 1;
    let il = im + la;
    let ix = il + n1 * n / 2 + 1;
    let ir = ix + n;
    let is = ir + n + n + la;
    let iu = is + n1;
    let iv = iu + n1;
    let iw = iv + n1;
    let w_at = |k: i32| w + (k - 1) as usize;
    slsqpb(mem, jw, m, meq, la, n, x, xl, xu, f, c, g, a, acc, iter, mode, w_at(ir), w_at(il), w_at(ix), w_at(im), w_at(is), w_at(iu), w_at(iv), w_at(iw), jw_p, state);
    state.x0 = Some(w_at(ix));
}

#[cfg(test)]
mod tests {
    //! Each test's expected values come from nlopt 2.11.0 itself (its Python module, LD_SLSQP,
    //! CadQuery's stopping settings), objective for objective the same arithmetic as here. The
    //! points evaluated were compared one by one with nlopt's: identical to the last bit, with
    //! the same choice of evaluating the gradient or not, in every test.
    use super::*;

    const INF: f64 = f64::INFINITY;

    /// CadQuery's settings (sketch.py: ftol_abs 0, ftol_rel 0, xtol_rel 1e-9, xtol_abs 1e-12).
    fn cadquery(n: usize, max_evaluations: usize) -> Stopping {
        Stopping { ftol_rel: 0.0, ftol_abs: 0.0, xtol_rel: 1e-9, xtol_abs: vec![1e-12; n], max_evaluations }
    }

    type Objective = fn(&[f64], Option<&mut [f64]>) -> f64;

    /// The outcome, and the points evaluated (each with whether the gradient was asked for).
    fn run(objective: Objective, x0: &[f64], lower: &[f64], upper: &[f64], max_evaluations: usize) -> (Outcome, Vec<(bool, Vec<f64>)>) {
        let mut trace = Vec::new();
        let outcome = minimize(x0, lower, upper, |x, gradient| {
            trace.push((gradient.is_some(), x.to_vec()));
            objective(x, gradient)
        }, &cadquery(x0.len(), max_evaluations));
        (outcome, trace)
    }

    fn assert_close(found: &[f64], expected: &[f64]) {
        assert_eq!(found.len(), expected.len());
        for (f, e) in found.iter().zip(expected) {
            assert!(f == e || (f - e).abs() <= 1e-12 * e.abs(), "{found:?} vs nlopt's {expected:?}");
        }
    }

    fn check(outcome: &Outcome, evaluations: usize, result: Result, x: &[f64], minimum: f64) {
        assert_eq!((outcome.evaluations, outcome.result), (evaluations, result), "{outcome:?}");
        assert_close(&outcome.x, x);
        assert_close(&[outcome.minimum], &[minimum]);
    }

    fn rosenbrock(x: &[f64], gradient: Option<&mut [f64]>) -> f64 {
        let (a, b) = (x[0], x[1]);
        if let Some(g) = gradient {
            g[0] = -2.0 * (1.0 - a) - 400.0 * a * (b - a * a);
            g[1] = 200.0 * (b - a * a);
        }
        (1.0 - a).powf(2.0) + 100.0 * (b - a * a).powf(2.0)
    }

    /// Rosenbrock's function from (−1.2, 1), unbounded. nlopt: 61 evaluations, XTOL_REACHED,
    /// x = (1.000000000000167, 1.0000000000003388), minimum 3.0267804072425e-26.
    #[test]
    fn rosenbrock_unbounded() {
        let (outcome, _) = run(rosenbrock, &[-1.2, 1.0], &[-INF, -INF], &[INF, INF], 0);
        check(&outcome, 61, Result::XtolReached, &[1.000000000000167, 1.0000000000003388], 3.0267804072425e-26);
    }

    /// Rosenbrock's function from (−1.2, 1) with x₀ ≤ 0.5, a bound active at the minimum. nlopt:
    /// 43 evaluations, XTOL_REACHED, x = (0.5, 0.25), minimum 0.25.
    #[test]
    fn rosenbrock_bounded() {
        let (outcome, _) = run(rosenbrock, &[-1.2, 1.0], &[-INF, -INF], &[0.5, INF], 0);
        check(&outcome, 43, Result::XtolReached, &[0.5, 0.25], 0.25);
    }

    /// Rosenbrock's function stopped after 10 evaluations. nlopt: 10 evaluations,
    /// MAXEVAL_REACHED, best x = (−0.4887764829211292, 0.18144781457729364), minimum
    /// 2.5465589322049.
    #[test]
    fn rosenbrock_evaluation_limit() {
        let (outcome, _) = run(rosenbrock, &[-1.2, 1.0], &[-INF, -INF], &[INF, INF], 10);
        check(&outcome, 10, Result::MaxevalReached, &[-0.4887764829211292, 0.18144781457729364], 2.5465589322049);
    }

    /// A chain of seven: the sum of squares of x_i − c_i and ½ (x_{i+1} − x_i²), with x₂ ≥ 0
    /// where c₂ = −2 pulls it below.
    fn chain(x: &[f64], gradient: Option<&mut [f64]>) -> f64 {
        const TARGETS: [f64; 7] = [1.0, -0.5, -2.0, 0.7, 1.2, -0.3, 0.4];
        let n = x.len();
        let mut f = 0.0;
        let mut g = vec![0.0; n];
        for i in 0..n {
            let r = x[i] - TARGETS[i];
            f += r * r;
            g[i] += 2.0 * r;
        }
        for i in 0..n - 1 {
            let r = 0.5 * (x[i + 1] - x[i] * x[i]);
            f += r * r;
            g[i + 1] += 2.0 * r * 0.5;
            g[i] += 2.0 * r * 0.5 * (-2.0 * x[i]);
        }
        if let Some(gradient) = gradient {
            gradient.copy_from_slice(&g);
        }
        f
    }

    /// nlopt: 42 evaluations, XTOL_REACHED, x = (0.7145665030370411, −0.28829445614082605, 0,
    /// 0.656919187851867, 0.800382664544451, −0.12753860216944587, 0.3232532189860444), minimum
    /// 4.797732343827071.
    #[test]
    fn chain_with_an_active_lower_bound() {
        let mut lower = [-INF; 7];
        lower[2] = 0.0;
        let (outcome, _) = run(chain, &[0.5; 7], &lower, &[INF; 7], 0);
        let x = [0.7145665030370411, -0.28829445614082605, 0.0, 0.656919187851867, 0.800382664544451, -0.12753860216944587, 0.3232532189860444];
        check(&outcome, 42, Result::XtolReached, &x, 4.797732343827071);
    }

    /// (x₀ + 1)² + (x₁ − 2)² + ½ x₀ x₁ + (x₀⁴ + x₁⁴) / 10, whose unbounded minimum has x₀ < 0,
    /// from (0, 0.3) on the bound x₀ ≥ 0 (and x₁ ≤ 1.5, not active at the end).
    fn corner(x: &[f64], gradient: Option<&mut [f64]>) -> f64 {
        let (a, b) = (x[0], x[1]);
        if let Some(g) = gradient {
            g[0] = 2.0 * (a + 1.0) + 0.5 * b + 0.4 * a.powf(3.0);
            g[1] = 2.0 * (b - 2.0) + 0.5 * a + 0.4 * b.powf(3.0);
        }
        (a + 1.0).powf(2.0) + (b - 2.0).powf(2.0) + 0.5 * a * b + 0.1 * (a.powf(4.0) + b.powf(4.0))
    }

    /// nlopt: 6 evaluations, XTOL_REACHED, x = (0, 1.423318338918933), minimum
    /// 1.7429625209960962; the points evaluated (each with the gradient) were (0, 0.3), (0, 1.5),
    /// (0, 1.3876765083440308), (4.126655166089014e-15, 1.4222716978709813),
    /// (2.4979481222136942e-15, 1.4233328033949326), (0, 1.423318338918933).
    #[test]
    fn started_on_a_bound() {
        let (outcome, trace) = run(corner, &[0.0, 0.3], &[0.0, -INF], &[INF, 1.5], 0);
        check(&outcome, 6, Result::XtolReached, &[0.0, 1.423318338918933], 1.7429625209960962);
        let points = [[0.0, 0.3], [0.0, 1.5], [0.0, 1.3876765083440308], [4.126655166089014e-15, 1.4222716978709813], [2.4979481222136942e-15, 1.4233328033949326], [0.0, 1.423318338918933]];
        assert_eq!(trace, points.iter().map(|p| (true, p.to_vec())).collect::<Vec<_>>());
    }

    /// (x₀ + 1)² + (x₁ − 2)² from (0, 1) with x₀ ≥ 0. nlopt: 4 evaluations, XTOL_REACHED,
    /// x = (0, 2), minimum 1; evaluated (0, 1) and (0, 3) with the gradient, (0, 2) without (a
    /// second line search step), then (0, 2) again with it (mode −1 after mode 1). It ends
    /// with mode 8 (no descent direction left), passing nlopt's relaxed x test.
    #[test]
    fn against_a_bound() {
        fn against_bound(x: &[f64], gradient: Option<&mut [f64]>) -> f64 {
            if let Some(g) = gradient {
                g[0] = 2.0 * (x[0] + 1.0);
                g[1] = 2.0 * (x[1] - 2.0);
            }
            (x[0] + 1.0).powf(2.0) + (x[1] - 2.0).powf(2.0)
        }
        let (outcome, trace) = run(against_bound, &[0.0, 1.0], &[0.0, -INF], &[INF, INF], 0);
        check(&outcome, 4, Result::XtolReached, &[0.0, 2.0], 1.0);
        assert_eq!(trace, vec![(true, vec![0.0, 1.0]), (true, vec![0.0, 3.0]), (false, vec![0.0, 2.0]), (true, vec![0.0, 2.0])]);
    }

    /// (x₀ − 1)² + x₁² from its minimum (1, 0): no descent direction, mode 8 at once. nlopt: 1
    /// evaluation, XTOL_REACHED, x = (1, 0), minimum 0.
    #[test]
    fn started_at_the_minimum() {
        fn at_minimum(x: &[f64], gradient: Option<&mut [f64]>) -> f64 {
            if let Some(g) = gradient {
                g[0] = 2.0 * (x[0] - 1.0);
                g[1] = 2.0 * x[1];
            }
            (x[0] - 1.0).powf(2.0) + x[1].powf(2.0)
        }
        let (outcome, _) = run(at_minimum, &[1.0, 0.0], &[-INF, -INF], &[INF, INF], 0);
        check(&outcome, 1, Result::XtolReached, &[1.0, 0.0], 0.0);
    }

    /// x₀ + x₁ = 1 and x₂ x₃ = 2 as two squared residuals in four variables: minima along
    /// a surface, as in an underdetermined sketch. nlopt: 17 evaluations, XTOL_REACHED,
    /// x = (0.20000000000003693, 0.800000000000037, 1.4142135623729186, 1.4142135623729186),
    /// minimum 2.5440261294550552e-25.
    #[test]
    fn underdetermined() {
        fn residuals(x: &[f64], gradient: Option<&mut [f64]>) -> f64 {
            let r0 = x[0] + x[1] - 1.0;
            let r1 = x[2] * x[3] - 2.0;
            if let Some(g) = gradient {
                g[0] = 2.0 * r0;
                g[1] = 2.0 * r0;
                g[2] = 2.0 * r1 * x[3];
                g[3] = 2.0 * r1 * x[2];
            }
            r0 * r0 + r1 * r1
        }
        let (outcome, _) = run(residuals, &[0.3, 0.9, 0.5, 0.5], &[-INF; 4], &[INF; 4], 0);
        check(&outcome, 17, Result::XtolReached, &[0.20000000000003693, 0.800000000000037, 1.4142135623729186, 1.4142135623729186], 2.5440261294550552e-25);
    }

    /// x − ln x + (x − 3)⁴ / 10 from 8, infinite (gradient untouched) for x ≤ 0, where the first
    /// step lands. nlopt: 12 evaluations, XTOL_REACHED, x = 1.9349850848618124, minimum
    /// 1.4035393083996515; every evaluation asks for the gradient, as an infinite value leaves
    /// nlopt wanting it, and 1.640625 is evaluated twice.
    #[test]
    fn infinite_values() {
        fn barrier(x: &[f64], gradient: Option<&mut [f64]>) -> f64 {
            if x[0] <= 0.0 {
                return INF;
            }
            if let Some(g) = gradient {
                g[0] = 1.0 - 1.0 / x[0] + 0.4 * (x[0] - 3.0).powf(3.0);
            }
            x[0] - x[0].ln() + 0.1 * (x[0] - 3.0).powf(4.0)
        }
        let (outcome, trace) = run(barrier, &[8.0], &[-INF], &[INF], 0);
        check(&outcome, 12, Result::XtolReached, &[1.9349850848618124], 1.4035393083996515);
        let points = [8.0, -42.875, -17.4375, 1.640625, 1.640625, 1.7164987740459585, 1.891254034960864, 1.927719512042081, 1.9347161511313826, 1.9349833923476116, 1.934985084466003, 1.9349850848618124];
        assert_eq!(trace, points.iter().map(|&p| (true, vec![p])).collect::<Vec<_>>());
    }

    /// Σ (x_i − c_i)² + x_i x_{i+1} + x_i⁴ / 20 in eight variables, all within [−0.6, 0.6], seven
    /// of the bounds active at the end, from 0. nlopt: 8 evaluations, XTOL_REACHED, x below,
    /// minimum 2.8183394944800955.
    #[test]
    fn mostly_on_bounds() {
        fn boxed(x: &[f64], gradient: Option<&mut [f64]>) -> f64 {
            const TARGETS: [f64; 8] = [1.0, -1.0, 0.5, 2.0, -2.0, 0.1, -0.4, 0.8];
            let n = x.len();
            let mut f = 0.0;
            let mut g = vec![0.0; n];
            for i in 0..n {
                f += (x[i] - TARGETS[i]).powf(2.0) + 0.05 * x[i].powf(4.0);
                g[i] += 2.0 * (x[i] - TARGETS[i]) + 0.2 * x[i].powf(3.0);
            }
            for i in 0..n - 1 {
                f += x[i] * x[i + 1];
                g[i] += x[i + 1];
                g[i + 1] += x[i];
            }
            if let Some(gradient) = gradient {
                gradient.copy_from_slice(&g);
            }
            f
        }
        let (outcome, _) = run(boxed, &[0.0; 8], &[-0.6; 8], &[0.6; 8], 0);
        let x = [0.5999999999999991, -0.5999999999999978, 0.48835331287253886, 0.5999999999999958, -0.5999999999999954, 0.5999999999999999, -0.5999999999999991, 0.599999999999999];
        check(&outcome, 8, Result::XtolReached, &x, 2.8183394944800955);
    }

    /// A start outside the bounds is refused, as nlopt_optimize refuses it.
    #[test]
    fn start_outside_the_bounds() {
        let (outcome, trace) = run(rosenbrock, &[-1.2, 1.0], &[-1.0, -INF], &[INF, INF], 0);
        assert_eq!((outcome.result, outcome.evaluations, trace.len()), (Result::InvalidArgs, 0, 0));
    }
}
