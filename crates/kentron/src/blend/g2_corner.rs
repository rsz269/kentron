//! The curvature-continuous corner of three set-back G2 fillets, by its support function.
//!
//! At a square corner whose three edges are rounded by the order-2 blend section (E. Hartmann,
//! "Parametric Gⁿ blending of curves and surfaces", The Visual Computer 17, 2001) and set back
//! `s` from the corner, the corner is closed by a convex patch given by its support function.
//! In the corner's frame (`u` along the three edges from the corner, the solid where all three
//! are positive) every boundary point near the corner has its inward normal `ν` in the positive
//! octant, and the solid near the corner is the set of points `x` with `x·ν ≥ g(ν)` for those
//! normals: `g` is its (lower) support function, concave and of degree one, and the boundary
//! point with normal `ν` is `∇g(ν)`. The faces are its vertices `ν = e_i` (`g = 0`), the fillets
//! its edges (`ν_i = 0`), the patch the octant's inside.
//!
//! `g` is a sum of support functions (a Minkowski sum of convex bodies, so convex):
//!
//! - `r·p(ν_j, ν_l)` for each edge `i`, where `p` is the support function of the unit section's
//!   region (the quadrant with its corner cut by the section, contacts `(1, 0)` and `(0, 1)`):
//!   the planar region of fillet `i`'s section in the plane square to its edge. On edge `i` of
//!   the normals it gives fillet `i` its section; on the other two it adds `r` along `e_i` to
//!   where fillet `i`'s neighbours' rulings start.
//! - `(s − 2r)·Φ(ν)` with `Φ = ⅓ Σ p(p(ν_a, ν_b), ν_c)` over the three turns of `(a, b, c)`:
//!   zero on every edge of the normals with gradient `e_i` on edge `i`, so it moves each
//!   fillet's start from `2r` to `s` without changing any section; near a face's normal it is
//!   `p` itself, so on each face the flat ends along the unit section scaled by `s − r` with its
//!   apex where the face's two fillets' contact lines would cross, `r` in from both edges.
//!
//! Its radii of curvature (the Hessian of `g` on the sphere) are the sums of the summands'.
//! Each section's radius grows without bound towards its contacts (it starts as straight as
//! the leg), so towards fillet `i` the patch's radius across the edge grows without bound and
//! along it tends to the section's: the patch meets the fillet with its curvature, and the
//! faces likewise (curvature continuous throughout, G2).
//!
//! For `s ≥ 2r` every summand is concave. For `r < s < 2r` the second is subtracted; then `g`
//! stays concave because `Σ p(ν_j, ν_l) − Φ` is concave on the octant: the radii of `Φ` are
//! nowhere larger than those of the three sections' sum, equal only in the limit at a face's
//! normal, where both are the section's. So the patch is convex for every `s > r`, and
//! degenerates at `s = r`, where the flats' boundaries shrink to the points where the contact
//! lines cross. (Checked on a dense grid of normals in the tests below.)
//!
//! Convex, the patch lies within the solid's faces (above every tangent plane, of the faces
//! too), and within the corner's cube of side `s`: each summand keeps `u_i` at most its share
//! of `s`. Its points are `∇g`; its height along the corner's diagonal over a point of the
//! plane square to it is the maximum over the normals of `g(ν) − x·ν`, found by Newton's
//! method (the objective concave).

use crate::error::{Error, Result};

/// Value and first three derivatives of a function of one variable.
#[derive(Clone, Copy, Debug)]
struct Jet {
    v: f64,
    d1: f64,
    d2: f64,
    d3: f64,
}

impl Jet {
    fn var(t: f64) -> Jet {
        Jet { v: t, d1: 1.0, d2: 0.0, d3: 0.0 }
    }

    fn constant(x: f64) -> Jet {
        Jet { v: x, d1: 0.0, d2: 0.0, d3: 0.0 }
    }

    fn add(self, o: Jet) -> Jet {
        Jet { v: self.v + o.v, d1: self.d1 + o.d1, d2: self.d2 + o.d2, d3: self.d3 + o.d3 }
    }

    fn sub(self, o: Jet) -> Jet {
        Jet { v: self.v - o.v, d1: self.d1 - o.d1, d2: self.d2 - o.d2, d3: self.d3 - o.d3 }
    }

    fn scale(self, k: f64) -> Jet {
        Jet { v: self.v * k, d1: self.d1 * k, d2: self.d2 * k, d3: self.d3 * k }
    }

    /// Leibniz's rule.
    fn mul(self, o: Jet) -> Jet {
        Jet {
            v: self.v * o.v,
            d1: self.d1 * o.v + self.v * o.d1,
            d2: self.d2 * o.v + 2.0 * self.d1 * o.d1 + self.v * o.d2,
            d3: self.d3 * o.v + 3.0 * self.d2 * o.d1 + 3.0 * self.d1 * o.d2 + self.v * o.d3,
        }
    }

    /// `1/f`: `−f'/f²`, `(2f'² − f f'')/f³`, `(−6f'³ + 6f f' f'' − f² f''')/f⁴`.
    fn recip(self) -> Jet {
        let (f, a, b, c) = (self.v, self.d1, self.d2, self.d3);
        let g = 1.0 / f;
        Jet { v: g, d1: -a * g * g, d2: (2.0 * a * a * g - b) * g * g, d3: (-6.0 * a * a * a * g * g + 6.0 * a * b * g - c) * g * g }
    }

    fn div(self, o: Jet) -> Jet {
        self.mul(o.recip())
    }
}

/// The order-2 blend section of thumbweight `w` on unit legs along the axes from the origin,
/// from `(1, 0)` at `t = 0` to `(0, 1)` at `t = 1`: `(α, β)` with `α = f·(1 − p₁)`,
/// `β = (1 − f)·p₂`, `f = t̄³/(t̄³ + t³)`, `p₁ = t/(1 − w + wt)`, `p₂ = t(1 − w)/(1 − wt)`.
#[derive(Clone, Debug)]
pub(crate) struct UnitSection {
    pub(crate) w: f64,
    /// Where the section's tangent turns square to the normals `(ρ, 1)`, `ρ` from 0 to 1 (the
    /// first half): `(√ρ, t)` at parameters spread as squares towards 0, where `ρ` grows as
    /// `t²`. Newton's method starts from it.
    guide: Vec<[f64; 2]>,
}

impl UnitSection {
    /// The section of thumbweight `w`.
    pub(crate) fn new(w: f64) -> UnitSection {
        let mut section = UnitSection { w, guide: Vec::new() };
        let n = 256;
        section.guide = (0..=n)
            .map(|k| {
                let t = 0.5 * (k as f64 / n as f64).powi(2);
                let [x, y] = section.jets(t);
                let rho = if k == 0 { 0.0 } else { (-y.d1 / x.d1).clamp(0.0, 1.0) };
                [rho.sqrt(), t]
            })
            .collect();
        section
    }

    fn jets(&self, t: f64) -> [Jet; 2] {
        let w = self.w;
        let t = Jet::var(t);
        let one = Jet::constant(1.0);
        let s = one.sub(t);
        let (a, b) = (s.mul(s).mul(s), t.mul(t).mul(t));
        // (`f`, `1 − f`, `1 − p₁` and `p₂` each as a quotient of its own, without the
        // cancellation of `1 − f` near `t = 0`, where `β` is of the order of `t³`.)
        let sum = a.add(b);
        let (f, g) = (a.div(sum), b.div(sum));
        let q1 = s.scale(1.0 - w).div(Jet::constant(1.0 - w).add(t.scale(w)));
        let p2 = t.scale(1.0 - w).div(one.sub(t.scale(w)));
        [f.mul(q1), g.mul(p2)]
    }

    /// The point at `t`.
    pub(crate) fn point(&self, t: f64) -> [f64; 2] {
        let [a, b] = self.jets(t);
        [a.v, b.v]
    }

    /// The share of the triangle on the legs that the section takes off, `∫ αβ' − α'β dt`
    /// (twice the area between the legs and the section), by Gauss-Legendre on 64 panels.
    pub(crate) fn area_share(&self) -> f64 {
        const X: [f64; 4] = [-0.861_136_311_594_052_6, -0.339_981_043_584_856_3, 0.339_981_043_584_856_3, 0.861_136_311_594_052_6];
        const W: [f64; 4] = [0.347_854_845_137_453_9, 0.652_145_154_862_546_1, 0.652_145_154_862_546_1, 0.347_854_845_137_453_9];
        let panels = 64;
        let h = 1.0 / panels as f64;
        let mut sum = 0.0;
        for k in 0..panels {
            let mid = (k as f64 + 0.5) * h;
            for (x, wt) in X.iter().zip(W) {
                let [a, b] = self.jets(mid + 0.5 * h * x);
                sum += wt * 0.5 * h * (a.v * b.d1 - a.d1 * b.v);
            }
        }
        sum
    }

    /// The section at the equal-area thumbweight of a square edge: the `w` at which it takes
    /// off what the quarter circle tangent at the same contacts does, `2 − π/2` of the
    /// triangle. The share falls as `w` rises (to 0 as `w` reaches 1), so `w` is bisected on a
    /// bracket stepped down from 0.
    pub(crate) fn equal_area_square() -> UnitSection {
        let target = 2.0 - std::f64::consts::FRAC_PI_2;
        let share = |w: f64| UnitSection { w, guide: Vec::new() }.area_share();
        let (mut lo, mut hi) = (0.0f64, 1.0 - 1e-9);
        while share(lo) < target {
            hi = lo;
            lo = 2.0 * lo - 1.0;
        }
        for _ in 0..200 {
            let mid = 0.5 * (lo + hi);
            if share(mid) > target {
                lo = mid;
            } else {
                hi = mid;
            }
            if hi - lo <= 1e-16 * (1.0 + lo.abs()) {
                break;
            }
        }
        UnitSection::new(0.5 * (lo + hi))
    }

    /// The support function of the section's region at the inward normal `(a, b)` (both at
    /// least 0): the least of `a·x + b·y` over it, with its gradient (the point of the section
    /// where it is least) and its Hessian `−c'c'ᵀ/(c''·n)` at that point's parameter. On an
    /// axis the point is the contact there and the Hessian is left zero (it is unbounded).
    fn support(&self, a: f64, b: f64) -> Support {
        if !(a > 0.0) || !(b > 0.0) {
            let point = if a > 0.0 { [0.0, 1.0] } else { [1.0, 0.0] };
            return Support { value: a.max(0.0) * point[0] + b.max(0.0) * point[1], point, hessian: [[0.0; 2]; 2] };
        }
        // (The section is its own mirror in the diagonal, `β(1 − t) = α(t)`: with `a > b` the
        // least lies past the middle, found as the mirror's at `(b, a)`, so that its parameter
        // is always near 0 rather than 1, where it is resolved finer.)
        if a > b {
            let m = self.support(b, a);
            let h = m.hessian;
            return Support { value: m.value, point: [m.point[1], m.point[0]], hessian: [[h[1][1], h[1][0]], [h[0][1], h[0][0]]] };
        }
        // The least is where the tangent is square to the normal: a·α' + b·β' rises through 0
        // (the section is convex). Newton's method kept within a shrinking bracket.
        let (mut lo, mut hi) = (0.0f64, 0.5f64 + 1e-12);
        let mut t = {
            let root = (a / b).sqrt();
            let k = self.guide.partition_point(|g| g[0] < root).clamp(1, self.guide.len().max(2) - 1);
            match (self.guide.get(k - 1), self.guide.get(k)) {
                (Some(p), Some(q)) if q[0] > p[0] => p[1] + (q[1] - p[1]) * (root - p[0]) / (q[0] - p[0]),
                _ => 0.25,
            }
        };
        for _ in 0..200 {
            let [x, y] = self.jets(t);
            let h = a * x.d1 + b * y.d1;
            if h < 0.0 {
                lo = t;
            } else {
                hi = t;
            }
            let dh = a * x.d2 + b * y.d2;
            let mut next = t - h / dh;
            if !(next > lo && next < hi) {
                next = 0.5 * (lo + hi);
            }
            let done = (next - t).abs() <= 4.0 * f64::EPSILON * t.max(1e-300) || hi - lo <= f64::EPSILON;
            t = next;
            if done {
                break;
            }
        }
        let [x, y] = self.jets(t);
        let c1 = [x.d1, y.d1];
        let k = a * x.d2 + b * y.d2;
        Support {
            value: a * x.v + b * y.v,
            point: [x.v, y.v],
            hessian: [[-c1[0] * c1[0] / k, -c1[0] * c1[1] / k], [-c1[1] * c1[0] / k, -c1[1] * c1[1] / k]],
        }
    }
}

/// A support function's value, gradient and Hessian at a normal of the plane.
#[derive(Clone, Copy, Debug)]
struct Support {
    value: f64,
    point: [f64; 2],
    hessian: [[f64; 2]; 2],
}

/// A support function's value, gradient and Hessian at a normal of space.
#[derive(Clone, Copy, Debug)]
pub(crate) struct Support3 {
    pub(crate) value: f64,
    pub(crate) point: [f64; 3],
    pub(crate) hessian: [[f64; 3]; 3],
}

impl Support3 {
    fn zero() -> Support3 {
        Support3 { value: 0.0, point: [0.0; 3], hessian: [[0.0; 3]; 3] }
    }

    fn add_scaled(&mut self, k: f64, o: &Support3) {
        self.value += k * o.value;
        for i in 0..3 {
            self.point[i] += k * o.point[i];
            for j in 0..3 {
                self.hessian[i][j] += k * o.hessian[i][j];
            }
        }
    }
}

/// The corner's support function `g = r·Σ p(ν_j, ν_l) + (s − 2r)·Φ` (see the module's notes),
/// in the corner's frame: coordinates along its three edges, normals in the positive octant.
#[derive(Clone, Debug)]
pub(crate) struct CornerSupport {
    pub(crate) section: UnitSection,
    pub(crate) radius: f64,
    pub(crate) setback: f64,
}

impl CornerSupport {
    /// The corner of fillets of tangent offset `radius` at `section`, set back `setback`:
    /// refused unless `setback > radius > 0` (at `setback = radius` the patch degenerates).
    pub(crate) fn new(section: UnitSection, radius: f64, setback: f64) -> Result<CornerSupport> {
        if !(radius > 0.0) || !(setback > radius) || !setback.is_finite() {
            return Err(Error::InvalidInput(format!("a G2 set-back corner needs a setback {setback} greater than the radius {radius} > 0")));
        }
        Ok(CornerSupport { section, radius, setback })
    }

    /// `Σ p(ν_j, ν_l)`: the three sections' regions in their planes.
    pub(crate) fn sections(&self, n: [f64; 3]) -> Support3 {
        let mut out = Support3::zero();
        for i in 0..3 {
            let (j, l) = ((i + 1) % 3, (i + 2) % 3);
            let s = self.section.support(n[j], n[l]);
            out.value += s.value;
            out.point[j] += s.point[0];
            out.point[l] += s.point[1];
            let at = [j, l];
            for p in 0..2 {
                for q in 0..2 {
                    out.hessian[at[p]][at[q]] += s.hessian[p][q];
                }
            }
        }
        out
    }

    /// `Φ = ⅓ Σ p(p(ν_a, ν_b), ν_c)` over the three turns of the axes, by the chain rule.
    pub(crate) fn set_back(&self, n: [f64; 3]) -> Support3 {
        let mut out = Support3::zero();
        for i in 0..3 {
            let (a, b, c) = (i, (i + 1) % 3, (i + 2) % 3);
            let inner = self.section.support(n[a], n[b]);
            let outer = self.section.support(inner.value, n[c]);
            out.value += outer.value / 3.0;
            out.point[a] += outer.point[0] * inner.point[0] / 3.0;
            out.point[b] += outer.point[0] * inner.point[1] / 3.0;
            out.point[c] += outer.point[1] / 3.0;
            // The Jacobian of (inner, n_c) by n, column by column.
            let column = |k: usize| -> [f64; 2] {
                if k == a {
                    [inner.point[0], 0.0]
                } else if k == b {
                    [inner.point[1], 0.0]
                } else {
                    [0.0, 1.0]
                }
            };
            for p in 0..3 {
                for q in 0..3 {
                    let (cp, cq) = (column(p), column(q));
                    let mut sum = 0.0;
                    for x in 0..2 {
                        for y in 0..2 {
                            sum += cp[x] * outer.hessian[x][y] * cq[y];
                        }
                    }
                    out.hessian[p][q] += sum / 3.0;
                }
            }
            let ab = [a, b];
            for x in 0..2 {
                for y in 0..2 {
                    out.hessian[ab[x]][ab[y]] += outer.point[0] * inner.hessian[x][y] / 3.0;
                }
            }
        }
        out
    }

    /// `g` at the inward normal `n` (all components positive), with its gradient (the
    /// boundary point whose normal is `n`) and Hessian.
    pub(crate) fn at(&self, n: [f64; 3]) -> Support3 {
        let mut out = Support3::zero();
        out.add_scaled(self.radius, &self.sections(n));
        out.add_scaled(self.setback - 2.0 * self.radius, &self.set_back(n));
        out
    }

    /// The patch's height above the point `x` of the plane through the corner square to its
    /// diagonal `(1, 1, 1)/√3`, with `x` given in the plane's basis `(1, −1, 0)/√2`,
    /// `(1, 1, −2)/√6`: `√3·max (g(ν) − x·ν)` over the normals with `Σν = 1`, and the normal
    /// where it is reached. `None` where the maximum is not inside the octant: off the
    /// patch, over a fillet or a face.
    pub(crate) fn height(&self, x: [f64; 2]) -> Option<(f64, [f64; 3])> {
        self.height_from(x, None)
    }

    /// [`CornerSupport::height`], Newton's method started at the normal `start` (with
    /// `Σν = 1`) if one is given: a neighbouring point's.
    pub(crate) fn height_from(&self, x: [f64; 2], start: Option<[f64; 3]>) -> Option<(f64, [f64; 3])> {
        let p = plane_point(x);
        let objective = |n: [f64; 3]| self.at(n).value - dot(p, n);
        let nu = |q: [f64; 2]| [q[0], q[1], 1.0 - q[0] - q[1]];
        let (u, v) = ([1.0, 0.0, -1.0], [0.0, 1.0, -1.0]);
        let mut q = match start {
            Some(n) if n.iter().all(|&c| c > 1e-9) => [n[0], n[1]],
            _ => [1.0 / 3.0, 1.0 / 3.0],
        };
        let mut converged = false;
        let mut g = self.at(nu(q));
        for _ in 0..200 {
            let n = nu(q);
            let d = [g.point[0] - p[0], g.point[1] - p[1], g.point[2] - p[2]];
            let gradient = [dot(d, u), dot(d, v)];
            let h = |a: [f64; 3], b: [f64; 3]| (0..3).map(|i| (0..3).map(|j| a[i] * g.hessian[i][j] * b[j]).sum::<f64>()).sum::<f64>();
            let m = [[h(u, u), h(u, v)], [h(v, u), h(v, v)]];
            let det = m[0][0] * m[1][1] - m[0][1] * m[1][0];
            if !(det > 0.0) || !det.is_finite() {
                return None;
            }
            let step = [-(m[1][1] * gradient[0] - m[0][1] * gradient[1]) / det, -(m[0][0] * gradient[1] - m[1][0] * gradient[0]) / det];
            let here = g.value - dot(p, n);
            let mut t = 1.0;
            let mut moved = false;
            for _ in 0..80 {
                let next = [q[0] + t * step[0], q[1] + t * step[1]];
                let m = nu(next);
                if m.iter().all(|&c| c > 0.0) {
                    let trial = self.at(m);
                    if trial.value - dot(p, m) >= here - 1e-15 * here.abs().max(1.0) {
                        q = next;
                        g = trial;
                        moved = true;
                        break;
                    }
                }
                t *= 0.5;
            }
            let size = t * (step[0].abs() + step[1].abs());
            if !moved || size <= 1e-15 * (1.0 + q[0].abs() + q[1].abs()) {
                converged = moved || size <= 1e-15;
                break;
            }
            // (Run off to the octant's edge: the maximum is on it, over a fillet or a face.)
            if nu(q).iter().any(|&c| c < 1e-14) {
                return None;
            }
        }
        let n = nu(q);
        if !converged && n.iter().any(|&c| c < 1e-12) {
            return None;
        }
        Some((3f64.sqrt() * objective(n), n))
    }
}

/// The volume of the solid within the corner's cube of side `s` (a reference for tests, from
/// the support function alone): by the divergence theorem about the corner, the faces through
/// it add nothing, each wall `u_i = s` adds `s/3` times its area (the cube's side less the
/// fillet's section, `s² − ½·share·r²`), and the patch `−⅓∫ g dA` over its normals, the area
/// element the product of the radii, `det M dω`. The octant's normals are cut into six
/// triangles, each from a corner `e_k` to the middle of an edge beside it and the centre,
/// mapped from a square by Duffy's map from `e_k`, crowded as fourth powers towards `e_k`
/// and towards the octant's edge (where the radii grow without bound).
#[cfg(test)]
pub(crate) fn volume_in_cube(c: &CornerSupport, panels: usize) -> f64 {
    // (As fourth powers: the area element grows faster than the inverse square root of the
    // distance to the edge, and squares left the sum converging as its cube root.)
    let (px, py) = (4, 4);
    const X: [f64; 4] = [-0.861_136_311_594_052_6, -0.339_981_043_584_856_3, 0.339_981_043_584_856_3, 0.861_136_311_594_052_6];
    const W: [f64; 4] = [0.347_854_845_137_453_9, 0.652_145_154_862_546_1, 0.652_145_154_862_546_1, 0.347_854_845_137_453_9];
    let nodes: Vec<(f64, f64)> = (0..panels).flat_map(|k| X.iter().zip(W).map(move |(x, w)| ((k as f64 + 0.5 + 0.5 * x) / panels as f64, w / panels as f64 * 0.5))).collect();
    let centre = [1.0 / 3.0; 3];
    let mut integral = 0.0;
    for k in 0..3 {
        let mut corner = [0.0; 3];
        corner[k] = 1.0;
        for other in [(k + 1) % 3, (k + 2) % 3] {
            let mut middle = [0.0; 3];
            middle[k] = 0.5;
            middle[other] = 0.5;
            // ν = corner + x((1 − y)(middle − corner) + y(centre − corner)), x = ξ², y = η².
            let area2 = {
                let (a, b): ([f64; 3], [f64; 3]) = (std::array::from_fn(|q| middle[q] - corner[q]), std::array::from_fn(|q| centre[q] - corner[q]));
                // The (ν_0, ν_1) parallelogram's area: da db = |det| x dx dy.
                (a[0] * b[1] - a[1] * b[0]).abs()
            };
            for &(xi, wx) in &nodes {
                let x = xi.powi(px);
                for &(eta, wy) in &nodes {
                    let y = eta.powi(py);
                    let n: [f64; 3] = std::array::from_fn(|q| corner[q] + x * ((1.0 - y) * (middle[q] - corner[q]) + y * (centre[q] - corner[q])));
                    if n.iter().any(|&v| !(v > 0.0)) {
                        continue;
                    }
                    let g = c.at(n);
                    let len2 = dot(n, n);
                    let len = len2.sqrt();
                    let unit = n.map(|v| v / len);
                    let a = if unit[0].abs() < 0.9 { [1.0, 0.0, 0.0] } else { [0.0, 1.0, 0.0] };
                    let d = dot(a, unit);
                    let t1 = {
                        let t = [a[0] - d * unit[0], a[1] - d * unit[1], a[2] - d * unit[2]];
                        let l = dot(t, t).sqrt();
                        t.map(|v| v / l)
                    };
                    let t2 = [unit[1] * t1[2] - unit[2] * t1[1], unit[2] * t1[0] - unit[0] * t1[2], unit[0] * t1[1] - unit[1] * t1[0]];
                    let q = |p: [f64; 3], r: [f64; 3]| -(0..3).map(|i| (0..3).map(|j| p[i] * g.hessian[i][j] * r[j]).sum::<f64>()).sum::<f64>();
                    let det = q(t1, t1) * q(t2, t2) - q(t1, t2) * q(t2, t1);
                    integral += wx * wy * (px as f64 * xi.powi(px - 1)) * (py as f64 * eta.powi(py - 1)) * area2 * x * g.value * det / len2;
                }
            }
        }
    }
    let (r, s) = (c.radius, c.setback);
    let section = 0.5 * c.section.area_share() * r * r;
    s * (s * s - section) - integral / 3.0
}

impl CornerSupport {
    /// The patch's principal curvatures (largest first) at its point of inward normal `n`:
    /// the reciprocals of its principal radii, the eigenvalues of minus `g`'s Hessian on the
    /// plane square to `n`.
    #[cfg(test)]
    pub(crate) fn curvatures(&self, n: [f64; 3]) -> [f64; 2] {
        let len = dot(n, n).sqrt();
        let unit = n.map(|v| v / len);
        let h = self.at(unit).hessian;
        let a = if unit[0].abs() < 0.9 { [1.0, 0.0, 0.0] } else { [0.0, 1.0, 0.0] };
        let d = dot(a, unit);
        let t1 = {
            let t = [a[0] - d * unit[0], a[1] - d * unit[1], a[2] - d * unit[2]];
            let l = dot(t, t).sqrt();
            t.map(|v| v / l)
        };
        let t2 = [unit[1] * t1[2] - unit[2] * t1[1], unit[2] * t1[0] - unit[0] * t1[2], unit[0] * t1[1] - unit[1] * t1[0]];
        let q = |p: [f64; 3], r: [f64; 3]| -(0..3).map(|i| (0..3).map(|j| p[i] * h[i][j] * r[j]).sum::<f64>()).sum::<f64>();
        let (m00, m01, m11) = (q(t1, t1), 0.5 * (q(t1, t2) + q(t2, t1)), q(t2, t2));
        let (mean, disc) = (0.5 * (m00 + m11), (0.25 * (m00 - m11).powi(2) + m01 * m01).sqrt());
        [1.0 / (mean - disc), 1.0 / (mean + disc)]
    }
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// The basis of the plane square to the corner's diagonal, in the corner's frame.
pub(crate) const PLANE: [[f64; 3]; 2] = [
    [std::f64::consts::FRAC_1_SQRT_2, -std::f64::consts::FRAC_1_SQRT_2, 0.0],
    [0.408_248_290_463_863, 0.408_248_290_463_863, -0.816_496_580_927_726],
];

/// The point of the plane with coordinates `x` in [`PLANE`], in the corner's frame.
pub(crate) fn plane_point(x: [f64; 2]) -> [f64; 3] {
    std::array::from_fn(|k| x[0] * PLANE[0][k] + x[1] * PLANE[1][k])
}

/// A point of the corner's frame projected onto the plane, in [`PLANE`]'s coordinates.
pub(crate) fn onto_plane(u: [f64; 3]) -> [f64; 2] {
    [dot(u, PLANE[0]), dot(u, PLANE[1])]
}

/// A point's height along the corner's diagonal.
pub(crate) fn diagonal(u: [f64; 3]) -> f64 {
    (u[0] + u[1] + u[2]) / 3f64.sqrt()
}

/// The patch as the graph of a B-spline function over a rectangle of the plane square to the
/// corner's diagonal (coordinates in [`PLANE`]), of degree [`Graph::DEGREE`] with even spans:
/// its height at the knot-span grid's points is `heights` (control values, one per basis
/// function), so that the surface `corner + x·PLANE + h·diagonal` is a B-spline surface whose
/// parameters are the plane's own coordinates.
#[derive(Clone, Debug)]
pub(crate) struct Graph {
    pub(crate) lo: [f64; 2],
    pub(crate) hi: [f64; 2],
    pub(crate) spans: usize,
    /// `heights[i][j]`: the control value of basis function `i` across, `j` up.
    pub(crate) heights: Vec<Vec<f64>>,
    /// The most the graph strays from the patch at the points it was checked at (none of
    /// them fitted): along its six boundary curves and inside.
    pub(crate) deviation: f64,
}

impl Graph {
    pub(crate) const DEGREE: usize = 5;

    /// The clamped knot vector of one direction, as distinct knots and multiplicities.
    pub(crate) fn knots(&self, k: usize) -> (Vec<f64>, Vec<usize>) {
        let p = Self::DEGREE;
        let knots: Vec<f64> = (0..=self.spans).map(|i| self.lo[k] + (self.hi[k] - self.lo[k]) * i as f64 / self.spans as f64).collect();
        let mut mults = vec![1; self.spans + 1];
        mults[0] = p + 1;
        mults[self.spans] = p + 1;
        (knots, mults)
    }

    /// The Greville abscissae of one direction: where a B-spline with these knots puts the
    /// control point of each basis function for the identity function.
    pub(crate) fn greville(&self, k: usize) -> Vec<f64> {
        let full = full_knots(self.lo[k], self.hi[k], self.spans);
        let p = Self::DEGREE;
        (0..self.spans + p).map(|i| full[i + 1..=i + p].iter().sum::<f64>() / p as f64).collect()
    }

    pub(crate) fn height(&self, x: [f64; 2]) -> f64 {
        let (iu, bu) = basis(self.lo[0], self.hi[0], self.spans, x[0]);
        let (iv, bv) = basis(self.lo[1], self.hi[1], self.spans, x[1]);
        let mut sum = 0.0;
        for (a, wa) in bu.iter().enumerate() {
            for (b, wb) in bv.iter().enumerate() {
                sum += wa * wb * self.heights[iu + a][iv + b];
            }
        }
        sum
    }
}

/// The full clamped knot vector of `spans` even spans from `lo` to `hi`.
fn full_knots(lo: f64, hi: f64, spans: usize) -> Vec<f64> {
    let p = Graph::DEGREE;
    let mut k = vec![lo; p + 1];
    k.extend((1..spans).map(|i| lo + (hi - lo) * i as f64 / spans as f64));
    k.extend(vec![hi; p + 1]);
    k
}

/// The first basis function not zero at `x`, and the values of the `DEGREE + 1` not zero
/// there (Cox-de Boor), on the even clamped knots.
fn basis(lo: f64, hi: f64, spans: usize, x: f64) -> (usize, Vec<f64>) {
    let p = Graph::DEGREE;
    let knots = full_knots(lo, hi, spans);
    let x = x.clamp(lo, hi);
    let span = (((x - lo) / (hi - lo) * spans as f64).floor() as usize).min(spans - 1) + p;
    let mut b = vec![0.0; p + 1];
    b[0] = 1.0;
    let (mut left, mut right) = (vec![0.0; p + 1], vec![0.0; p + 1]);
    for j in 1..=p {
        left[j] = x - knots[span + 1 - j];
        right[j] = knots[span + j] - x;
        let mut saved = 0.0;
        for r in 0..j {
            let t = b[r] / (right[r + 1] + left[j - r]);
            b[r] = saved + right[r + 1] * t;
            saved = left[j - r] * t;
        }
        b[j] = saved;
    }
    (span - p, b)
}

/// The values and first two derivatives of the `DEGREE + 1` basis functions not zero at `x`,
/// from the first of them: the derivatives by the recurrence
/// `N'ᵢ,ₚ = p·(Nᵢ,ₚ₋₁/(tᵢ₊ₚ − tᵢ) − Nᵢ₊₁,ₚ₋₁/(tᵢ₊ₚ₊₁ − tᵢ₊₁))`, applied to the functions of
/// one and two degrees lower (a term over a zero interval left out).
fn basis_jets(lo: f64, hi: f64, spans: usize, x: f64) -> (usize, [Vec<f64>; 3]) {
    let p = Graph::DEGREE;
    let knots = full_knots(lo, hi, spans);
    let x = x.clamp(lo, hi);
    let span = (((x - lo) / (hi - lo) * spans as f64).floor() as usize).min(spans - 1) + p;
    // The functions of degree q not zero at x: N_{span-q..=span, q}.
    let of_degree = |q: usize| -> Vec<f64> {
        let mut b = vec![0.0; q + 1];
        b[0] = 1.0;
        let (mut left, mut right) = (vec![0.0; q + 1], vec![0.0; q + 1]);
        for j in 1..=q {
            left[j] = x - knots[span + 1 - j];
            right[j] = knots[span + j] - x;
            let mut saved = 0.0;
            for r in 0..j {
                let t = b[r] / (right[r + 1] + left[j - r]);
                b[r] = saved + right[r + 1] * t;
                saved = left[j - r] * t;
            }
            b[j] = saved;
        }
        b
    };
    // From the values (or derivatives) `lower` of the functions of degree q − 1 not zero at
    // x, those of degree q: element m of each is the function span − q + m.
    let raise = |q: usize, lower: &[f64]| -> Vec<f64> {
        (0..=q)
            .map(|m| {
                let i = span - q + m;
                let term = |k: usize, value: f64| {
                    let d = knots[k + q] - knots[k];
                    if d > 0.0 { value / d } else { 0.0 }
                };
                let a = if m >= 1 { term(i, lower[m - 1]) } else { 0.0 };
                let c = if m < q { term(i + 1, lower[m]) } else { 0.0 };
                q as f64 * (a - c)
            })
            .collect()
    };
    let (values, below, below2) = (of_degree(p), of_degree(p - 1), of_degree(p - 2));
    let slopes = raise(p, &below);
    let bends = raise(p, &raise(p - 1, &below2));
    (span - p, [values, slopes, bends])
}

impl CornerSupport {
    /// The patch's six boundary curves, each sampled at `n + 1` parameters: fillet `i`'s end
    /// section `s` along edge `i`, then the flat's boundary on face `i` (the section scaled by
    /// `s − r` with its apex `r` in from both edges), in the corner's frame.
    pub(crate) fn boundary(&self, n: usize) -> Vec<Vec<[f64; 3]>> {
        let (r, s) = (self.radius, self.setback);
        let mut out = Vec::new();
        for i in 0..3 {
            let (j, l) = ((i + 1) % 3, (i + 2) % 3);
            let (mut seam, mut flat) = (Vec::new(), Vec::new());
            for k in 0..=n {
                let p = self.section.point(k as f64 / n as f64);
                let mut u = [0.0; 3];
                u[i] = s;
                u[j] = r * p[0];
                u[l] = r * p[1];
                seam.push(u);
                let mut v = [0.0; 3];
                v[j] = r + (s - r) * p[0];
                v[l] = r + (s - r) * p[1];
                flat.push(v);
            }
            out.push(seam);
            out.push(flat);
        }
        out
    }

    /// The patch's height, slopes and second derivatives over the plane at
    /// [`CornerSupport::boundary`]'s points (of the faces' and fillets' there, which it
    /// shares): with `ν` the inward normal, `Tₐ = PLANE[a] + hₐ·diagonal` the tangents and
    /// `II` the second fundamental form (the section's curvature across fillet `i`'s edge,
    /// none along it; none on a face), `hₐ = −ν·PLANE[a]/ν·diagonal` and
    /// `hₐᵦ = II(Tₐ, Tᵦ)/ν·diagonal`.
    fn boundary_jets(&self, n: usize) -> Vec<([f64; 2], f64, [f64; 2], [f64; 3])> {
        let r = self.radius;
        let curves = self.boundary(n);
        let mut out = Vec::new();
        let up = [1.0 / 3f64.sqrt(); 3];
        for (c, curve) in curves.iter().enumerate() {
            let (i, seam) = (c / 2, c % 2 == 0);
            let (j, l) = ((i + 1) % 3, (i + 2) % 3);
            for (k, u) in curve.iter().enumerate() {
                let (normal, bend) = if seam {
                    let [x, y] = self.section.jets(k as f64 / n as f64);
                    let speed = x.d1.hypot(y.d1);
                    let curvature = (x.d1 * y.d2 - y.d1 * x.d2).abs() / speed.powi(3) / r;
                    let mut tangent = [0.0; 3];
                    tangent[j] = x.d1 / speed;
                    tangent[l] = y.d1 / speed;
                    let mut normal = [0.0; 3];
                    normal[j] = y.d1 / speed;
                    normal[l] = -x.d1 / speed;
                    if normal[j] < 0.0 || (normal[j] == 0.0 && normal[l] < 0.0) {
                        normal = normal.map(|v| -v);
                    }
                    (normal, Some((curvature, tangent)))
                } else {
                    let mut normal = [0.0; 3];
                    normal[i] = 1.0;
                    (normal, None)
                };
                let along = dot(normal, up);
                let slopes = [-dot(normal, PLANE[0]) / along, -dot(normal, PLANE[1]) / along];
                let tangents: [[f64; 3]; 2] = std::array::from_fn(|a| std::array::from_fn(|q| PLANE[a][q] + slopes[a] * up[q]));
                let second = match bend {
                    Some((curvature, t)) => {
                        let (a, b) = (dot(tangents[0], t), dot(tangents[1], t));
                        [curvature * a * a / along, curvature * a * b / along, curvature * b * b / along]
                    }
                    None => [0.0; 3],
                };
                out.push((onto_plane(*u), diagonal(*u), slopes, second));
            }
        }
        out
    }

    /// The patch fitted as a [`Graph`] of `spans` spans each way over its boundary's box (a
    /// fortieth larger each side): least squares through the patch's heights on a grid two
    /// to a span inside it, and along each of its boundary curves at 16 points a span through
    /// the heights, slopes and second derivatives of the fillets and faces there (so that the
    /// patch's boundary, where the fit has no points beyond, bends as they do). A slight
    /// weight on the control values' second differences keeps the system definite (the
    /// functions no point holds, outside the patch, taken from their neighbours).
    pub(crate) fn graph(&self, spans: usize) -> Result<Graph> {
        let p = Graph::DEGREE;
        let ncp = spans + p;
        let edges = self.boundary(32 * spans);
        let (mut lo, mut hi) = ([f64::INFINITY; 2], [f64::NEG_INFINITY; 2]);
        for u in edges.iter().flatten() {
            let x = onto_plane(*u);
            for k in 0..2 {
                lo[k] = lo[k].min(x[k]);
                hi[k] = hi[k].max(x[k]);
            }
        }
        let margin = (hi[0] - lo[0]).max(hi[1] - lo[1]) / 40.0;
        for k in 0..2 {
            lo[k] -= margin;
            hi[k] += margin;
        }
        let grid = 2 * spans;
        let at = |k: usize, i: usize| lo[k] + (hi[k] - lo[k]) * (i as f64 + 0.5) / grid as f64;
        let lines: Vec<usize> = (0..grid).collect();
        let rows: Vec<Vec<([f64; 2], f64)>> = crate::parallel::map(&lines, |&i| {
            let mut last = None;
            (0..grid).filter_map(|j| {
                let x = [at(0, i), at(1, j)];
                let found = self.height_from(x, last);
                last = found.map(|(_, n)| n);
                found.map(|(h, _)| (x, h))
            }).collect()
        });
        // The normal equations, banded: the basis functions (i, j) and (i', j') meet only
        // with |i − i'| and |j − j'| at most the degree.
        let n = ncp * ncp;
        let band = p * ncp + p;
        let mut a = vec![vec![0.0; band + 1]; n];
        let mut rhs = vec![0.0; n];
        let mut fit = |x: [f64; 2], terms: &[(usize, usize, f64)], value: f64, a: &mut Vec<Vec<f64>>| {
            // `terms`: (derivatives across, up, weight) summed against `value`.
            let (iu, bu) = basis_jets(lo[0], hi[0], spans, x[0]);
            let (iv, bv) = basis_jets(lo[1], hi[1], spans, x[1]);
            let row: Vec<(usize, f64)> = (0..=p)
                .flat_map(|q| (0..=p).map(move |m| (q, m)))
                .map(|(q, m)| ((iu + q) * ncp + iv + m, terms.iter().map(|&(du, dv, w)| w * bu[du][q] * bv[dv][m]).sum::<f64>()))
                .collect();
            for &(k, wk) in &row {
                rhs[k] += wk * value;
                for &(l, wl) in &row {
                    if k <= l {
                        a[k][l - k] += wk * wl;
                    }
                }
            }
        };
        for (x, h) in rows.into_iter().flatten() {
            fit(x, &[(0, 0, 1.0)], h, &mut a);
        }
        // Along the boundary, the height and its slopes and second derivatives (the patch
        // shares the fillets' and faces'), each weighed by a span's width to its order: a
        // slope off by e over a span puts the height off by about e times the width.
        let width = ((hi[0] - lo[0]) / spans as f64).max((hi[1] - lo[1]) / spans as f64);
        for (x, h, g, b) in self.boundary_jets(16 * spans) {
            fit(x, &[(0, 0, 1.0)], h, &mut a);
            fit(x, &[(1, 0, width)], width * g[0], &mut a);
            fit(x, &[(0, 1, width)], width * g[1], &mut a);
            fit(x, &[(2, 0, width * width)], width * width * b[0], &mut a);
            fit(x, &[(1, 1, width * width)], width * width * b[1], &mut a);
            fit(x, &[(0, 2, width * width)], width * width * b[2], &mut a);
        }
        let add = |k: usize, l: usize, v: f64, a: &mut Vec<Vec<f64>>| {
            let (k, l) = if k <= l { (k, l) } else { (l, k) };
            a[k][l - k] += v;
        };
        let weight = 1e-9;
        for i in 0..ncp {
            for j in 0..ncp {
                for (di, dj) in [(1, 0), (0, 1)] {
                    if i + 2 * di < ncp && j + 2 * dj < ncp {
                        let ks = [i * ncp + j, (i + di) * ncp + j + dj, (i + 2 * di) * ncp + j + 2 * dj];
                        let w = [1.0, -2.0, 1.0];
                        for x in 0..3 {
                            for y in 0..3 {
                                if ks[x] <= ks[y] {
                                    add(ks[x], ks[y], weight * w[x] * w[y], &mut a);
                                }
                            }
                        }
                    }
                }
            }
        }
        let values = banded_cholesky_solve(a, rhs, band).ok_or_else(|| Error::Failed("a G2 set-back corner's patch could not be fitted".into()))?;
        let heights: Vec<Vec<f64>> = (0..ncp).map(|i| values[i * ncp..(i + 1) * ncp].to_vec()).collect();
        let mut graph = Graph { lo, hi, spans, heights, deviation: 0.0 };
        // Checked between the fitted points: along the boundary curves, and inside on a grid
        // offset from the fitting one.
        let mut deviation: f64 = 0.0;
        for curve in &edges {
            for u in curve.iter().skip(1).step_by(2) {
                deviation = deviation.max((graph.height(onto_plane(*u)) - diagonal(*u)).abs());
            }
        }
        let check = 3 * spans / 2 + 1;
        let lines: Vec<usize> = (0..check).collect();
        let inside: Vec<f64> = crate::parallel::map(&lines, |&i| {
            let mut last = None;
            (0..check).filter_map(|j| {
                let x = [lo[0] + (hi[0] - lo[0]) * (i as f64 + 0.31) / check as f64, lo[1] + (hi[1] - lo[1]) * (j as f64 + 0.73) / check as f64];
                let found = self.height_from(x, last);
                last = found.map(|(_, n)| n);
                found.map(|(h, _)| (graph.height(x) - h).abs())
            }).fold(0.0, f64::max)
        });
        deviation = inside.into_iter().fold(deviation, f64::max);
        graph.deviation = deviation;
        Ok(graph)
    }
}

/// Solves the symmetric positive definite system whose upper band (`a[k][d]` the entry at
/// row `k`, column `k + d`) is given, by Cholesky's factorization in the band, row by row:
/// `L` lower with `L Lᵀ = A`, each row's entries in the band in one slice, so that each entry
/// is a dot product of two rows' slices.
fn banded_cholesky_solve(a: Vec<Vec<f64>>, mut b: Vec<f64>, band: usize) -> Option<Vec<f64>> {
    let n = b.len();
    let w = band + 1;
    // l[k * w + d] = L[k][k − band + d] (d = band the diagonal).
    let mut l = vec![0.0; n * w];
    for k in 0..n {
        let first = k.saturating_sub(band);
        for j in first..=k {
            // A[k][j] = A[j][k], stored at a[j][k − j].
            let mut sum = a[j][k - j];
            let from = first.max(j.saturating_sub(band));
            let (row_k, row_j) = (&l[k * w + (from + band - k)..k * w + (j + band - k)], &l[j * w + (from + band - j)..j * w + band]);
            sum -= row_k.iter().zip(row_j).map(|(x, y)| x * y).sum::<f64>();
            if j == k {
                if !(sum > 0.0) {
                    return None;
                }
                l[k * w + band] = sum.sqrt();
            } else {
                l[k * w + (j + band - k)] = sum / l[j * w + band];
            }
        }
    }
    for k in 0..n {
        let first = k.saturating_sub(band);
        let sum: f64 = (first..k).map(|m| l[k * w + (m + band - k)] * b[m]).sum();
        b[k] = (b[k] - sum) / l[k * w + band];
    }
    for k in (0..n).rev() {
        let last = (k + band).min(n - 1);
        let sum: f64 = (k + 1..=last).map(|j| l[j * w + (k + band - j)] * b[j]).sum();
        b[k] = (b[k] - sum) / l[k * w + band];
    }
    Some(b)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn corner(s: f64) -> CornerSupport {
        CornerSupport::new(UnitSection::equal_area_square(), 1.0, s).unwrap()
    }

    /// The 2x2 radii of curvature (minus the Hessian on the sphere's tangent plane at `n`).
    fn radii(h: &[[f64; 3]; 3], n: [f64; 3]) -> [[f64; 2]; 2] {
        let len = dot(n, n).sqrt();
        let n = n.map(|c| c / len);
        let a = if n[0].abs() < 0.9 { [1.0, 0.0, 0.0] } else { [0.0, 1.0, 0.0] };
        let d = dot(a, n);
        let t1 = {
            let t = [a[0] - d * n[0], a[1] - d * n[1], a[2] - d * n[2]];
            let l = dot(t, t).sqrt();
            t.map(|c| c / l)
        };
        let t2 = [n[1] * t1[2] - n[2] * t1[1], n[2] * t1[0] - n[0] * t1[2], n[0] * t1[1] - n[1] * t1[0]];
        let q = |x: [f64; 3], y: [f64; 3]| -(0..3).map(|i| (0..3).map(|j| x[i] * h[i][j] * y[j]).sum::<f64>()).sum::<f64>();
        [[q(t1, t1), q(t1, t2)], [q(t2, t1), q(t2, t2)]]
    }

    /// The least `λ` with `det(A − λB) = 0`, for `B` positive definite.
    fn least_ratio(a: [[f64; 2]; 2], b: [[f64; 2]; 2]) -> f64 {
        let qa = b[0][0] * b[1][1] - b[0][1] * b[1][0];
        let qb = -(a[0][0] * b[1][1] + a[1][1] * b[0][0] - a[0][1] * b[1][0] - a[1][0] * b[0][1]);
        let qc = a[0][0] * a[1][1] - a[0][1] * a[1][0];
        let disc = (qb * qb - 4.0 * qa * qc).max(0.0).sqrt();
        ((-qb - disc) / (2.0 * qa)).min((-qb + disc) / (2.0 * qa))
    }

    /// The equal-area section takes off `2 − π/2` of its triangle (the quarter circle's
    /// share), here by the midpoint rule on the definition with central differences; its
    /// thumbweight is about −0.525, and its middle's radius about 1.389 times the offset.
    #[test]
    fn the_equal_area_section_takes_off_the_circle_s_share() {
        let s = UnitSection::equal_area_square();
        let n = 200_000;
        let h = 1e-7;
        let mut share = 0.0;
        for k in 0..n {
            let t = (k as f64 + 0.5) / n as f64;
            let (p, a, b) = (s.point(t), s.point(t - h), s.point(t + h));
            share += (p[0] * (b[1] - a[1]) - (b[0] - a[0]) * p[1]) / (2.0 * h) / n as f64;
        }
        assert!((share - (2.0 - std::f64::consts::FRAC_PI_2)).abs() < 1e-9, "{share}");
        assert!((s.w + 0.525).abs() < 1e-3, "{}", s.w);
        let [x, y] = s.jets(0.5);
        let radius = (x.d1 * x.d1 + y.d1 * y.d1).powf(1.5) / (x.d1 * y.d2 - y.d1 * x.d2).abs();
        assert!((radius - 1.389).abs() < 1e-3, "{radius}");
    }

    /// The support function's gradient is its point, and its Hessian the gradient's
    /// derivative, against central differences.
    #[test]
    fn the_support_functions_derivatives_are_their_differences() {
        let c = corner(2.6);
        for n in [[0.2, 0.3, 0.5], [0.05, 0.9, 0.05], [0.6, 0.01, 0.39]] {
            let g = c.at(n);
            let h = 1e-6;
            for k in 0..3 {
                let mut a = n;
                let mut b = n;
                a[k] -= h;
                b[k] += h;
                let (ga, gb) = (c.at(a), c.at(b));
                let slope = (gb.value - ga.value) / (2.0 * h);
                assert!((slope - g.point[k]).abs() < 1e-8, "{n:?} {k}: {slope} {}", g.point[k]);
                for j in 0..3 {
                    let d = (gb.point[j] - ga.point[j]) / (2.0 * h);
                    assert!((d - g.hessian[j][k]).abs() < 1e-5 * (1.0 + d.abs()), "{n:?} {j} {k}: {d} {}", g.hessian[j][k]);
                }
            }
            // Of degree one: g(n) = n·∇g(n).
            assert!((dot(n, g.point) - g.value).abs() < 1e-12, "{n:?}");
        }
    }

    /// The bound behind every set-back above the radius: the radii of `Φ` are nowhere larger
    /// than the sections' sum (the least generalized eigenvalue at least 1), on a grid of
    /// normals spreading by powers of ten towards the octant's edges and corners; and both are
    /// positive definite there.
    #[test]
    fn the_set_back_term_never_outbends_the_sections() {
        let c = corner(3.0);
        let mut least = f64::INFINITY;
        let mut check = |n: [f64; 3]| {
            let (k, d) = (radii(&c.sections(n).hessian, n), radii(&c.set_back(n).hessian, n));
            for m in [k, d] {
                assert!(m[0][0] > 0.0 && m[0][0] * m[1][1] - m[0][1] * m[1][0] > 0.0, "{n:?} {m:?}");
            }
            least = least.min(least_ratio(k, d));
        };
        let steps = 60;
        for i in 0..=steps {
            for j in 0..=steps {
                let (a, b) = (10f64.powf(-9.0 * i as f64 / steps as f64), 10f64.powf(-9.0 * j as f64 / steps as f64));
                check([a, b, 1.0]);
                check([a, 1.0, b]);
            }
        }
        let k = 120;
        for i in 1..k {
            for j in 1..k - i {
                check([i as f64, j as f64, (k - i - j) as f64].map(|x| x / k as f64));
            }
        }
        assert!(least >= 1.0 && least < 1.01, "{least}");
    }

    /// The fitted graph converges on the patch as its spans shrink, checked between the
    /// fitted points: for a corner set back 3.5 radii, by 1e-3 of the radius at 40 spans
    /// across its five radii, 1e-4 at 120.
    #[test]
    fn the_graph_converges_on_the_patch() {
        let c = corner(3.5);
        let (coarse, fine) = (c.graph(40).unwrap(), c.graph(120).unwrap());
        assert!(coarse.deviation < 2e-3 && fine.deviation < 1e-4 && fine.deviation < coarse.deviation / 8.0, "{} {}", coarse.deviation, fine.deviation);
        // On the patch's boundary curves themselves, as on the patch.
        for curve in c.boundary(50) {
            for u in curve {
                assert!((fine.height(onto_plane(u)) - diagonal(u)).abs() < 1e-4);
            }
        }
    }

    /// Curvature continuity in the limit: towards fillet i's edge of the normals the patch's
    /// curvatures tend to the section's (across the edge) and none (along it); towards a
    /// face's normal, to none. The differences shrink with the normal's distance from the
    /// edge or corner, to a thousandth of the section's at a millionth.
    #[test]
    fn the_patch_s_curvature_tends_to_the_fillets_and_the_faces() {
        for s in [1.3, 2.9] {
            let c = corner(s);
            for t in [0.15, 0.3, 0.5] {
                let [x, y] = c.section.jets(t);
                let normal = [y.d1, -x.d1];
                let normal = if normal[0] < 0.0 { [-normal[0], -normal[1]] } else { normal };
                let section = (x.d1 * y.d2 - y.d1 * x.d2).abs() / (x.d1 * x.d1 + y.d1 * y.d1).powf(1.5);
                let mut last = f64::INFINITY;
                for e in [1e-2, 1e-4, 1e-6, 1e-8] {
                    let k = c.curvatures([e, normal[0], normal[1]]);
                    let off = (k[0] - section).abs().max(k[1].abs());
                    assert!(off < last, "{s} {t} {e}: {k:?} {section}");
                    last = off;
                }
                assert!(last < 1e-3 * section, "{s} {t}: {last}");
                let mut last = f64::INFINITY;
                for e in [1e-2, 1e-4, 1e-6, 1e-8] {
                    let k = c.curvatures([1.0, e * normal[0], e * normal[1]]);
                    assert!(k[0] < last, "{s} {t} {e}: {k:?}");
                    last = k[0];
                }
                assert!(last < 1e-3, "{s} {t}: {last}");
            }
        }
    }

    /// The patch's limits: at the normals of fillet i's edge it reaches fillet i's section `s`
    /// along its edge; at a face's normal, the flat's boundary is the section scaled by `s − r`
    /// with its apex `r` in from both edges; everywhere it lies within the cube of side `s`
    /// and inside the faces.
    #[test]
    fn the_patch_ends_on_the_sections_and_the_flats() {
        for s in [1.1, 1.7, 2.5, 3.5] {
            let c = corner(s);
            for t in [0.1, 0.3, 0.5, 0.8] {
                // Towards edge 0 of the normals (n0 → 0) at the section's normal at t.
                let [x, y] = c.section.jets(t);
                let normal = [y.d1, -x.d1];
                let normal = if normal[0] < 0.0 { [-normal[0], -normal[1]] } else { normal };
                let edge = c.at([1e-24, normal[0], normal[1]]).point;
                let p = c.section.point(t);
                assert!((edge[0] - s).abs() < 1e-5 && (edge[1] - p[0]).abs() < 1e-5 && (edge[2] - p[1]).abs() < 1e-5, "{s} {t}: {edge:?}");
                // Towards face 0 (n1, n2 → 0) in the same direction.
                let flat = c.at([1.0, 1e-24 * normal[0], 1e-24 * normal[1]]).point;
                assert!(flat[0].abs() < 1e-5, "{flat:?}");
                assert!((flat[1] - (1.0 + (s - 1.0) * p[0])).abs() < 1e-5 && (flat[2] - (1.0 + (s - 1.0) * p[1])).abs() < 1e-5, "{s} {t}: {flat:?} {p:?}");
            }
            let k = 60;
            for i in 1..k {
                for j in 1..k - i {
                    let u = c.at([i as f64, j as f64, (k - i - j) as f64].map(|x| x / k as f64)).point;
                    assert!(u.iter().all(|&v| v > 0.0 && v < s), "{s}: {u:?}");
                }
            }
        }
    }
}
