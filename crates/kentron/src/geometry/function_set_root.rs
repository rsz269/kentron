//! A root of a system of equations near a start, within a box: a bounded local minimisation
//! of `½‖F‖²` by Newton (Gauss–Newton) steps.
//!
//! Each step solves `J Δ = −F` in the least-squares, minimum-norm sense from a singular value
//! decomposition of the Jacobian (one-sided Jacobi; singular values below a relative threshold
//! dropped), so that square, singular, tall and wide systems are all handled alike; it falls
//! back to steepest descent where that is not a descent direction. Steps are cut to the box
//! and to half its extent, and searched along: accepted when `½‖F‖²` drops enough (Armijo),
//! extended while it keeps dropping when it drops slowly, otherwise cut back by safeguarded
//! quadratic interpolation. Variables on a face of the box with the gradient pushing outwards
//! are held there and the step is solved for the others, so that the iterate slides along the
//! box (Nocedal & Wright, *Numerical Optimization*, ch. 10–11 and 16; Dennis & Schnabel,
//! ch. 6).

/// A system of `nb_equations()` functions of `nb_variables()` variables with its Jacobian.
pub(crate) trait FunctionSetWithDerivatives {
    fn nb_variables(&self) -> usize;
    fn nb_equations(&self) -> usize;
    /// `F(x)` into `f` and the Jacobian into `d`, row-major
    /// (`d[i·nb_variables() + j] = ∂F_i/∂x_j`); `false` where it cannot be evaluated.
    fn values(&mut self, x: &[f64], f: &mut [f64], d: &mut [f64]) -> bool;
    /// `F(x)` alone.
    fn value(&mut self, x: &[f64], f: &mut [f64]) -> bool {
        let mut d = vec![0.0; self.nb_equations() * self.nb_variables()];
        self.values(x, f, &mut d)
    }
    /// Notification that the last evaluation was at the solution.
    fn get_state_number(&mut self) {}
}

/// The solver, reusable for several starts.
pub(crate) struct FunctionSetRoot {
    n: usize,
    iteration_limit: usize,
    tolerance: Vec<f64>,
    done: bool,
    root: Vec<f64>,
    iterations: usize,
}

/// Values at least this large stand for infinity.
const HUGE: f64 = 0.5 * crate::tolerance::INFINITE;

/// A squared norm this small is zero.
const NEGLIGIBLE_SQUARED: f64 = 1e-300;

/// The state at an iterate.
#[derive(Clone)]
struct Point {
    x: Vec<f64>,
    f: Vec<f64>,
    /// The Jacobian, row-major.
    j: Vec<f64>,
    /// `½‖F‖²`.
    phi: f64,
    /// `Jᵀ F`.
    gradient: Vec<f64>,
}

impl FunctionSetRoot {
    pub(crate) fn new(n: usize, iteration_limit: usize) -> Self {
        FunctionSetRoot { n, iteration_limit, tolerance: vec![0.0; n], done: false, root: vec![0.0; n], iterations: 0 }
    }

    /// Sets the per-variable step tolerances (the first `n` entries of `tolerance`).
    pub(crate) fn set_tolerance(&mut self, tolerance: &[f64]) {
        self.tolerance.copy_from_slice(&tolerance[..self.n]);
    }

    pub(crate) fn is_done(&self) -> bool {
        self.done
    }

    /// The last point reached (the solution when done).
    pub(crate) fn root(&self) -> &[f64] {
        &self.root
    }

    pub(crate) fn nb_iterations(&self) -> usize {
        self.iterations
    }

    /// Solves from `start` (clamped into the box `inf ≤ x ≤ sup`).
    pub(crate) fn perform(&mut self, f: &mut dyn FunctionSetWithDerivatives, start: &[f64], inf: &[f64], sup: &[f64]) {
        let n = self.n;
        let m = f.nb_equations();
        assert!(m > 0 && f.nb_variables() == n && start.len() == n && inf.len() == n && sup.len() == n, "a system of the solver's size");
        self.done = false;
        self.iterations = 0;
        let lo: Vec<f64> = (0..n).map(|i| inf[i].min(sup[i])).collect();
        let hi: Vec<f64> = (0..n).map(|i| inf[i].max(sup[i])).collect();
        let x0: Vec<f64> = (0..n).map(|i| start[i].clamp(lo[i], hi[i])).collect();
        self.root.copy_from_slice(&x0);
        let last_evaluated: std::cell::RefCell<Vec<f64>> = std::cell::RefCell::new(Vec::new());
        let eval = |f: &mut dyn FunctionSetWithDerivatives, x: &[f64]| -> Option<Point> {
            let mut fv = vec![0.0; m];
            let mut j = vec![0.0; m * n];
            *last_evaluated.borrow_mut() = x.to_vec();
            if !f.values(x, &mut fv, &mut j) || fv.iter().any(|v| !v.is_finite() || v.abs() >= HUGE) {
                return None;
            }
            let gradient: Vec<f64> = (0..n).map(|c| (0..m).map(|r| j[r * n + c] * fv[r]).sum()).collect();
            if gradient.iter().any(|g| !g.is_finite()) {
                return None;
            }
            let phi = 0.5 * fv.iter().map(|v| v * v).sum::<f64>();
            Some(Point { x: x.to_vec(), f: fv, j, phi, gradient })
        };
        let Some(mut here) = eval(f, &x0) else { return };
        let stationary = |p: &Point, free: &[bool]| {
            2.0 * p.phi <= NEGLIGIBLE_SQUARED || p.gradient.iter().zip(free).filter(|(_, &fr)| fr).map(|(g, _)| g * g).sum::<f64>() <= NEGLIGIBLE_SQUARED
        };
        if stationary(&here, &vec![true; n]) {
            self.done = true;
            f.get_state_number();
            return;
        }
        let extent: Vec<f64> = (0..n).map(|i| (hi[i] - lo[i]).min(1e100)).collect();
        let mut slow = 0;
        for iteration in 1..=self.iteration_limit {
            self.iterations = iteration;
            // Variables held on a face of the box: the gradient pushing out of it.
            let on_lo = |i: usize, p: &Point| p.x[i] <= lo[i];
            let on_hi = |i: usize, p: &Point| p.x[i] >= hi[i];
            let mut free: Vec<bool> = (0..n).map(|i| !((on_lo(i, &here) && here.gradient[i] > 0.0) || (on_hi(i, &here) && here.gradient[i] < 0.0) || lo[i] == hi[i])).collect();
            if !free.iter().any(|&b| b) || stationary(&here, &free) {
                break;
            }
            // The Newton direction for the free variables, held ones added where it points out.
            let mut direction = Vec::new();
            for _ in 0..=n {
                let d = newton_step(&here, &free, m, n).filter(|d| {
                    let slope: f64 = d.iter().zip(&here.gradient).map(|(a, b)| a * b).sum();
                    slope < 0.0 && d.iter().all(|v| v.is_finite())
                });
                let d = match d {
                    Some(d) => d,
                    None => (0..n).map(|i| if free[i] { -here.gradient[i] } else { 0.0 }).collect(),
                };
                // A free variable on a face the step leaves by: held too.
                let blocked: Vec<usize> = (0..n).filter(|&i| free[i] && ((on_lo(i, &here) && d[i] < 0.0) || (on_hi(i, &here) && d[i] > 0.0))).collect();
                direction = d;
                if blocked.is_empty() {
                    break;
                }
                for i in blocked {
                    free[i] = false;
                }
                if !free.iter().any(|&b| b) {
                    direction = vec![0.0; n];
                    break;
                }
            }
            if direction.iter().all(|&v| v == 0.0) {
                break;
            }
            // Steps of at most half the box's extent.
            let mut scale: f64 = 1.0;
            for i in 0..n {
                if direction[i] != 0.0 && extent[i] > 0.0 {
                    scale = scale.min(0.5 * extent[i] / direction[i].abs());
                }
            }
            for v in &mut direction {
                *v *= scale;
            }
            // The longest step inside the box.
            let mut t_box = f64::INFINITY;
            for i in 0..n {
                if direction[i] > 0.0 {
                    t_box = t_box.min((hi[i] - here.x[i]) / direction[i]);
                } else if direction[i] < 0.0 {
                    t_box = t_box.min((lo[i] - here.x[i]) / direction[i]);
                }
            }
            let at = |p: &Point, t: f64| -> Vec<f64> { (0..n).map(|i| (p.x[i] + t * direction[i]).clamp(lo[i], hi[i])).collect() };
            let slope: f64 = direction.iter().zip(&here.gradient).map(|(a, b)| a * b).sum();
            let Some(next) = line_search(&mut |x: &[f64]| eval(f, x), &here, &at, slope, t_box) else {
                if self.root.iter().zip(&here.x).any(|(a, b)| a != b) {
                    self.root.copy_from_slice(&here.x);
                }
                return;
            };
            let Some((next, full)) = next else {
                // No decrease along the step: a (constrained) minimum, within rounding.
                break;
            };
            let moved_within = (0..n).all(|i| (next.x[i] - here.x[i]).abs() <= self.tolerance[i]);
            let drop = here.phi - next.phi;
            if drop <= 1e-6 * here.phi {
                slow += 1;
            } else {
                slow = 0;
            }
            here = next;
            self.root.copy_from_slice(&here.x);
            if stationary(&here, &vec![true; n]) || (moved_within && full) || slow >= 3 {
                self.done = true;
                break;
            }
            if iteration == self.iteration_limit {
                return;
            }
        }
        self.root.copy_from_slice(&here.x);
        self.done = true;
        let last = last_evaluated.borrow().clone();
        if last != here.x {
            if eval(f, &here.x).is_none() {
                self.done = false;
                return;
            }
        }
        f.get_state_number();
    }
}

/// The least-squares, minimum-norm solution of `J_free Δ = −F` (zero for held variables), by
/// the singular value decomposition of the free columns.
fn newton_step(p: &Point, free: &[bool], m: usize, n: usize) -> Option<Vec<f64>> {
    let columns: Vec<usize> = (0..n).filter(|&i| free[i]).collect();
    let k = columns.len();
    if k == 0 {
        return None;
    }
    // A = J restricted to the free columns, m × k, column-major for the Jacobi rotations.
    let mut a: Vec<Vec<f64>> = columns.iter().map(|&c| (0..m).map(|r| p.j[r * n + c]).collect()).collect();
    let mut v: Vec<Vec<f64>> = (0..k).map(|i| (0..k).map(|j| f64::from(i == j)).collect()).collect();
    // One-sided Jacobi: rotate column pairs until they are orthogonal.
    for _sweep in 0..60 {
        let mut rotated = false;
        for i in 0..k {
            for j in i + 1..k {
                let (alpha, beta, gamma) = (dot(&a[i], &a[i]), dot(&a[j], &a[j]), dot(&a[i], &a[j]));
                if gamma == 0.0 || gamma.abs() <= 1e-15 * (alpha * beta).sqrt() {
                    continue;
                }
                rotated = true;
                let zeta = (beta - alpha) / (2.0 * gamma);
                let t = zeta.signum() / (zeta.abs() + (1.0 + zeta * zeta).sqrt());
                let t = if zeta == 0.0 { 1.0 } else { t };
                let c = 1.0 / (1.0 + t * t).sqrt();
                let s = c * t;
                for r in 0..m {
                    let (x, y) = (a[i][r], a[j][r]);
                    a[i][r] = c * x - s * y;
                    a[j][r] = s * x + c * y;
                }
                for r in 0..k {
                    let (x, y) = (v[i][r], v[j][r]);
                    v[i][r] = c * x - s * y;
                    v[j][r] = s * x + c * y;
                }
            }
        }
        if !rotated {
            break;
        }
    }
    // A V = U Σ: column i of A is σ_i u_i. Δ = −V Σ⁺ Uᵀ F.
    let sigma: Vec<f64> = a.iter().map(|col| dot(col, col).sqrt()).collect();
    let largest = sigma.iter().copied().fold(0.0, f64::max);
    if !(largest > 0.0) || !largest.is_finite() {
        return None;
    }
    let cut = largest * 1e-12 * (m.max(k) as f64);
    let mut step = vec![0.0; n];
    for i in 0..k {
        if sigma[i] <= cut {
            continue;
        }
        // uᵢᵀF / σᵢ = (A_i · F) / σᵢ².
        let coefficient = -dot(&a[i], &p.f) / (sigma[i] * sigma[i]);
        for (r, &c) in columns.iter().enumerate() {
            step[c] += coefficient * v[i][r];
        }
    }
    Some(step)
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

/// A step along `at(here, t)`, `0 < t ≤ min(1, t_box)` first, `φ` decreasing: `Ok(Some((point,
/// full)))` with `full` when the whole first step was taken (or more), `Ok(None)` when no
/// decrease is found, `Err` on an evaluation failure (as `None` of the outer option).
#[allow(clippy::type_complexity)]
fn line_search(eval: &mut dyn FnMut(&[f64]) -> Option<Point>, here: &Point, at: &dyn Fn(&Point, f64) -> Vec<f64>, slope: f64, t_box: f64) -> Option<Option<(Point, bool)>> {
    const ARMIJO: f64 = 1e-4;
    let t_first = t_box.min(1.0);
    if !(t_first > 0.0) {
        return Some(None);
    }
    let first = eval(&at(here, t_first))?;
    if first.phi <= here.phi + ARMIJO * t_first * slope {
        // Dropping slowly (towards a multiple root, along a valley): on while it drops.
        let (mut best, mut t) = (first, t_first);
        if best.phi > 0.25 * here.phi {
            for _ in 0..20 {
                if t >= t_box {
                    break;
                }
                let t_next = (2.0 * t).min(t_box);
                let Some(p) = eval(&at(here, t_next)) else { break };
                if p.phi >= best.phi {
                    break;
                }
                (best, t) = (p, t_next);
            }
        }
        return Some(Some((best, true)));
    }
    // Back by safeguarded quadratic interpolation.
    let (mut t, mut phi_t) = (t_first, first.phi);
    let mut best: Option<Point> = (first.phi < here.phi).then_some(first);
    for _ in 0..40 {
        let denominator = 2.0 * (phi_t - here.phi - slope * t);
        let mut t_new = if denominator > 0.0 { -slope * t * t / denominator } else { 0.5 * t };
        t_new = t_new.clamp(0.1 * t, 0.5 * t);
        if t_new <= 1e-20 * t_first.max(1.0) {
            break;
        }
        let p = eval(&at(here, t_new))?;
        if p.x == here.x {
            break;
        }
        t = t_new;
        phi_t = p.phi;
        if p.phi <= here.phi + ARMIJO * t * slope {
            return Some(Some((p, false)));
        }
        if p.phi < here.phi && best.as_ref().is_none_or(|b| p.phi < b.phi) {
            best = Some(p);
        }
    }
    Some(best.map(|b| (b, false)))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A system from a closure computing F and J.
    struct System<G: FnMut(&[f64], &mut [f64], &mut [f64]) -> bool> {
        n: usize,
        m: usize,
        g: G,
        evaluations: usize,
        states: usize,
        last: Vec<f64>,
    }

    impl<G: FnMut(&[f64], &mut [f64], &mut [f64]) -> bool> FunctionSetWithDerivatives for System<G> {
        fn nb_variables(&self) -> usize {
            self.n
        }
        fn nb_equations(&self) -> usize {
            self.m
        }
        fn values(&mut self, x: &[f64], f: &mut [f64], d: &mut [f64]) -> bool {
            self.evaluations += 1;
            self.last = x.to_vec();
            (self.g)(x, f, d)
        }
        fn get_state_number(&mut self) {
            self.states += 1;
        }
    }

    fn system<G: FnMut(&[f64], &mut [f64], &mut [f64]) -> bool>(n: usize, m: usize, g: G) -> System<G> {
        System { n, m, g, evaluations: 0, states: 0, last: Vec::new() }
    }

    const BIG: f64 = 1e100;

    #[test]
    fn linear_system() {
        let mut s = system(3, 3, |x, f, d| {
            let a = [[2.0, 1.0, -1.0], [-3.0, -1.0, 2.0], [-2.0, 1.0, 2.0]];
            let b = [8.0, -11.0, -3.0];
            for r in 0..3 {
                f[r] = (0..3).map(|c| a[r][c] * x[c]).sum::<f64>() - b[r];
                for c in 0..3 {
                    d[r * 3 + c] = a[r][c];
                }
            }
            true
        });
        let mut solver = FunctionSetRoot::new(3, 100);
        solver.set_tolerance(&[1e-10; 3]);
        solver.perform(&mut s, &[0.0; 3], &[-BIG; 3], &[BIG; 3]);
        assert!(solver.is_done());
        assert!(solver.nb_iterations() <= 2, "{}", solver.nb_iterations());
        for (x, e) in solver.root().iter().zip([2.0, 3.0, -1.0]) {
            assert!((x - e).abs() < 1e-12);
        }
        assert!(s.states >= 1);
        assert_eq!(s.last, solver.root());
    }

    /// The line (t, t, t)·(1, 0.5, 0.25) + (0, 0.1, 0) meets the unit sphere.
    fn sphere_and_line(x: &[f64], f: &mut [f64], d: &mut [f64]) -> bool {
        // Unknowns (u, v, w): sphere angles u, v and line parameter w.
        let (u, v, w) = (x[0], x[1], x[2]);
        let p = [v.cos() * u.cos(), v.cos() * u.sin(), v.sin()];
        let du = [-v.cos() * u.sin(), v.cos() * u.cos(), 0.0];
        let dv = [-v.sin() * u.cos(), -v.sin() * u.sin(), v.cos()];
        let dir = [1.0, 0.5, 0.25];
        let l = [w * dir[0], 0.1 + w * dir[1], w * dir[2]];
        for r in 0..3 {
            f[r] = p[r] - l[r];
            d[r * 3] = du[r];
            d[r * 3 + 1] = dv[r];
            d[r * 3 + 2] = -dir[r];
        }
        true
    }

    #[test]
    fn line_meets_sphere() {
        let mut s = system(3, 3, sphere_and_line);
        let mut solver = FunctionSetRoot::new(3, 100);
        solver.set_tolerance(&[1e-10; 3]);
        let pi = std::f64::consts::PI;
        solver.perform(&mut s, &[0.3, 0.1, 0.8], &[0.0, -pi / 2.0, -10.0], &[2.0 * pi, pi / 2.0, 10.0]);
        assert!(solver.is_done());
        let r = solver.root();
        let mut f = [0.0; 3];
        let mut d = [0.0; 9];
        sphere_and_line(r, &mut f, &mut d);
        assert!(f.iter().all(|v| v.abs() < 1e-14), "{f:?}");
        // The exact w: |w·dir + (0, 0.1, 0)| = 1.
        let (a, b, c): (f64, f64, f64) = (1.0 + 0.25 + 0.0625, 2.0 * 0.1 * 0.5, 0.01 - 1.0);
        let w = (-b + (b * b - 4.0 * a * c).sqrt()) / (2.0 * a);
        assert!((r[2] - w).abs() < 1e-12);
    }

    #[test]
    fn root_on_a_bound_from_outside() {
        let mut s = system(2, 2, |x, f, d| {
            f[0] = x[0] - 1.0;
            f[1] = x[1] * x[1] - 0.25;
            d.copy_from_slice(&[1.0, 0.0, 0.0, 2.0 * x[1]]);
            true
        });
        let mut solver = FunctionSetRoot::new(2, 100);
        solver.set_tolerance(&[1e-12; 2]);
        solver.perform(&mut s, &[3.0, 2.0], &[0.0, 0.0], &[1.0, 1.0]);
        assert!(solver.is_done());
        assert_eq!(solver.root()[0], 1.0);
        assert!((solver.root()[1] - 0.5).abs() < 1e-12);
    }

    #[test]
    fn least_squares_along_a_bound() {
        let c = [1.7, -0.4, 0.3];
        let mut s = system(3, 3, move |x, f, d| {
            for i in 0..3 {
                f[i] = x[i] - c[i];
                for j in 0..3 {
                    d[i * 3 + j] = f64::from(i == j);
                }
            }
            true
        });
        let mut solver = FunctionSetRoot::new(3, 100);
        solver.set_tolerance(&[1e-12; 3]);
        solver.perform(&mut s, &[0.5, 0.5, 0.5], &[0.0; 3], &[1.0; 3]);
        assert!(solver.is_done());
        assert_eq!(solver.root(), [1.0, 0.0, 0.3]);
    }

    #[test]
    fn singular_jacobian_at_the_start() {
        // x² = 1, y = x: J singular at x = 0.
        let mut s = system(2, 2, |x, f, d| {
            f[0] = x[0] * x[0] - 1.0;
            f[1] = x[1] - x[0];
            d.copy_from_slice(&[2.0 * x[0], 0.0, -1.0, 1.0]);
            true
        });
        let mut solver = FunctionSetRoot::new(2, 100);
        solver.set_tolerance(&[1e-12; 2]);
        solver.perform(&mut s, &[0.0, 0.3], &[-5.0, -5.0], &[5.0, 5.0]);
        assert!(solver.is_done());
        let r = solver.root();
        assert!((r[0].abs() - 1.0).abs() < 1e-12 && (r[1] - r[0]).abs() < 1e-12, "{r:?}");
    }

    #[test]
    fn already_a_root() {
        let mut s = system(2, 2, |x, f, d| {
            f[0] = x[0] - 0.5;
            f[1] = x[1] + 0.25;
            d.copy_from_slice(&[1.0, 0.0, 0.0, 1.0]);
            true
        });
        let mut solver = FunctionSetRoot::new(2, 100);
        solver.perform(&mut s, &[0.5, -0.25], &[-1.0; 2], &[1.0; 2]);
        assert!(solver.is_done());
        assert_eq!(solver.nb_iterations(), 0);
        assert_eq!(s.evaluations, 1);
        assert_eq!(s.states, 1);
    }

    #[test]
    fn wide_and_tall_systems() {
        // Wide: x + y + z = 1, x − y = 0: any point of the line.
        let mut wide = system(3, 2, |x, f, d| {
            f[0] = x[0] + x[1] + x[2] - 1.0;
            f[1] = x[0] - x[1];
            d.copy_from_slice(&[1.0, 1.0, 1.0, 1.0, -1.0, 0.0]);
            true
        });
        let mut solver = FunctionSetRoot::new(3, 100);
        solver.set_tolerance(&[1e-12; 3]);
        solver.perform(&mut wide, &[0.0; 3], &[-BIG; 3], &[BIG; 3]);
        assert!(solver.is_done());
        let r = solver.root().to_vec();
        assert!((r[0] + r[1] + r[2] - 1.0).abs() < 1e-14 && (r[0] - r[1]).abs() < 1e-14);
        // Tall, consistent: x = 2, y = −1, x + y = 1.
        let mut tall = system(2, 3, |x, f, d| {
            f[0] = x[0] - 2.0;
            f[1] = x[1] + 1.0;
            f[2] = x[0] + x[1] - 1.0;
            d.copy_from_slice(&[1.0, 0.0, 0.0, 1.0, 1.0, 1.0]);
            true
        });
        let mut solver = FunctionSetRoot::new(2, 100);
        solver.set_tolerance(&[1e-12; 2]);
        solver.perform(&mut tall, &[0.0; 2], &[-BIG; 2], &[BIG; 2]);
        assert!(solver.is_done());
        assert!((solver.root()[0] - 2.0).abs() < 1e-14 && (solver.root()[1] + 1.0).abs() < 1e-14);
    }

    #[test]
    fn rosenbrock_and_repeated_starts() {
        let mut s = system(2, 2, |x, f, d| {
            f[0] = 10.0 * (x[1] - x[0] * x[0]);
            f[1] = 1.0 - x[0];
            d.copy_from_slice(&[-20.0 * x[0], 10.0, -1.0, 0.0]);
            true
        });
        let mut solver = FunctionSetRoot::new(2, 100);
        solver.set_tolerance(&[1e-12; 2]);
        for start in [[-1.2, 1.0], [3.0, -2.0], [0.0, 0.0]] {
            solver.perform(&mut s, &start, &[-5.0; 2], &[5.0; 2]);
            assert!(solver.is_done());
            assert!((solver.root()[0] - 1.0).abs() < 1e-12 && (solver.root()[1] - 1.0).abs() < 1e-12, "{:?}", solver.root());
        }
        // A root on a corner of the box.
        let mut corner = system(2, 2, |x, f, d| {
            f[0] = x[0] - 1.0;
            f[1] = x[1] + 2.0;
            d.copy_from_slice(&[1.0, 0.0, 0.0, 1.0]);
            true
        });
        solver.perform(&mut corner, &[0.1, 0.1], &[-2.0, -2.0], &[1.0, 1.0]);
        assert!(solver.is_done());
        assert_eq!(solver.root(), [1.0, -2.0]);
    }

    #[test]
    fn failing_at_the_start() {
        let mut s = system(2, 2, |_, _, _| false);
        let mut solver = FunctionSetRoot::new(2, 100);
        solver.perform(&mut s, &[0.0; 2], &[-1.0; 2], &[1.0; 2]);
        assert!(!solver.is_done());
        assert_eq!(s.states, 0);
    }
}
