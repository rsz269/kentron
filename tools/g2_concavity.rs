//! A computer-assisted check, in rigorous interval arithmetic, of the concavity lemma behind
//! kentron's curvature-continuous set-back corners (see docs/g2-setback-corner-proof.md).
//!
//! Build and run (no dependencies; the proof's step is `scaled`, a few seconds):
//!
//!     rustc -O --edition 2021 tools/g2_concavity.rs -o /tmp/g2_concavity
//!     /tmp/g2_concavity scaled 0 1 0 1        # Lemma 5: N positive definite on [0, 1]²
//!
//! It works on the slice c = 1 of the octant of normals, where (by the proof's Lemma 2) minus
//! the Hessian of H(a, b) = Σp − Φ is the 2 × 2 matrix
//!
//!     M = κ(a,b)(1 − f'(X)/3) q qᵀ + F(a)(1 − f'(Z/b)/3) e_a e_aᵀ + F(b)(1 − f'(Y/a)/3) e_b e_bᵀ
//!         − ⅓ [F(X) ∇X ∇Xᵀ + κ(Y,a) v vᵀ + κ(Z,b) v' v'ᵀ]
//!
//! with f(s) = p(s, 1), F = −f'', κ(x, y) = F(x/y)/y³, X = p(a, b), Y = f(b), Z = f(a),
//! q = (b, −a), v = (Y, −a f'(b)), v' = (b f'(a), −Z), ∇X = (f'(a/b), f(a/b) − (a/b) f'(a/b)).
//!
//! Modes:
//!
//! - `scaled B0 B1 R0 R1 [N MAXDEPTH THREADS CELLS]`: the proof's computation (Lemma 5). With
//!   b = β³, a/b = ρ³, proves N = β²·diag(ρ, 1)·M·diag(ρ, 1), in its form free of singular
//!   factors (Lemma 4), positive definite over β in [B0, B1], ρ in [R0, R1], the boundary
//!   included; [0, 1]² covers the lemma. The section's functions come from a table of
//!   enclosures over cells of its parameter t (the table's construction certifies Lemma 1's
//!   sign conditions, or stops).
//! - `check`: in floating point, M against central differences of H itself (H computed from
//!   the section's support function alone): Lemma 2's formula, cross-checked.
//! - `ncheck`: N's singularity-free form against β²·diag(ρ,1)·M·diag(ρ,1) from M's plain form,
//!   at a few points: Lemma 4's algebra, cross-checked.
//! - `certify B0 B1 R0 R1 [N MAXDEPTH THREADS]`: an independent check of the inside, not part
//!   of the proof: M itself (its plain form, mean-value enclosures over boxes of b and a/b)
//!   certified positive definite over b in [B0, B1], a/b in [R0, R1].
//! - `counterexample`: certifies, at one point, that the per-turn term Σp − p(p(a, b), c) is
//!   not concave (the decomposition that would have made the lemma immediate fails).
//!
//! The thumbweight w is taken as the interval W_LO..W_HI, which holds the value kentron's
//! `UnitSection::equal_area_square` computes (−0.5252124045775235) with room to spare.

use std::ops::{Add, Div, Mul, Neg, Sub};

const W_LO: f64 = -0.525_212_404_577_6;
const W_HI: f64 = -0.525_212_404_577_4;

/// A closed interval of reals, its ends rounded outwards after every operation (each
/// operation's float result lies within half an ulp of the exact one; one ulp more each side
/// holds it).
#[derive(Clone, Copy, Debug, PartialEq)]
struct I {
    lo: f64,
    hi: f64,
}

fn down(x: f64) -> f64 {
    x.next_down()
}
fn up(x: f64) -> f64 {
    x.next_up()
}

impl I {
    fn new(lo: f64, hi: f64) -> I {
        assert!(lo <= hi, "{lo} > {hi}");
        I { lo, hi }
    }
    fn point(x: f64) -> I {
        I { lo: x, hi: x }
    }
    fn hull(self, o: I) -> I {
        I { lo: self.lo.min(o.lo), hi: self.hi.max(o.hi) }
    }
    fn contains_zero(self) -> bool {
        self.lo <= 0.0 && self.hi >= 0.0
    }
    fn sqr(self) -> I {
        if self.lo.is_nan() || self.hi.is_nan() {
            return ALL;
        }
        let (a, b) = (self.lo * self.lo, self.hi * self.hi);
        if self.contains_zero() {
            I { lo: 0.0, hi: up(a.max(b)) }
        } else {
            I { lo: down(a.min(b)), hi: up(a.max(b)) }
        }
    }
    fn powi(self, n: u32) -> I {
        let mut r = I::point(1.0);
        for _ in 0..n {
            r = r * self;
        }
        r
    }
    fn mid(self) -> f64 {
        0.5 * (self.lo + self.hi)
    }
}

/// The whole line: what an operation gives that cannot bound its result (a division by an
/// interval holding 0, or infinities meeting); no box is certified through it.
const ALL: I = I { lo: f64::NEG_INFINITY, hi: f64::INFINITY };

fn checked(lo: f64, hi: f64) -> I {
    if lo.is_nan() || hi.is_nan() { ALL } else { I { lo: down(lo), hi: up(hi) } }
}

impl Add for I {
    type Output = I;
    fn add(self, o: I) -> I {
        checked(self.lo + o.lo, self.hi + o.hi)
    }
}
impl Sub for I {
    type Output = I;
    fn sub(self, o: I) -> I {
        checked(self.lo - o.hi, self.hi - o.lo)
    }
}
impl Neg for I {
    type Output = I;
    fn neg(self) -> I {
        I { lo: -self.hi, hi: -self.lo }
    }
}
impl Mul for I {
    type Output = I;
    fn mul(self, o: I) -> I {
        let c = [self.lo * o.lo, self.lo * o.hi, self.hi * o.lo, self.hi * o.hi];
        if c.iter().any(|x| x.is_nan()) {
            return ALL;
        }
        let lo = c.iter().cloned().fold(f64::INFINITY, f64::min);
        let hi = c.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        // (A product of nonnegatives is nonnegative: its rounded lower end is kept at 0.)
        let lo = if self.lo >= 0.0 && o.lo >= 0.0 { down(lo).max(0.0) } else { down(lo) };
        I { lo, hi: up(hi) }
    }
}
impl Div for I {
    type Output = I;
    fn div(self, o: I) -> I {
        if o.contains_zero() || o.lo.is_nan() || o.hi.is_nan() {
            return ALL;
        }
        let c = [self.lo / o.lo, self.lo / o.hi, self.hi / o.lo, self.hi / o.hi];
        if c.iter().any(|x| x.is_nan()) {
            return ALL;
        }
        let lo = c.iter().cloned().fold(f64::INFINITY, f64::min);
        let hi = c.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        // (A nonnegative over a positive is nonnegative.)
        let lo = if self.lo >= 0.0 && o.lo > 0.0 { down(lo).max(0.0) } else { down(lo) };
        I { lo, hi: up(hi) }
    }
}
impl Mul<f64> for I {
    type Output = I;
    fn mul(self, k: f64) -> I {
        self * I::point(k)
    }
}
impl Add<f64> for I {
    type Output = I;
    fn add(self, k: f64) -> I {
        self + I::point(k)
    }
}
impl Sub<I> for f64 {
    type Output = I;
    fn sub(self, o: I) -> I {
        I::point(self) - o
    }
}

/// A function's value and first three derivatives at a point (all intervals).
#[derive(Clone, Copy, Debug)]
struct J {
    v: I,
    d1: I,
    d2: I,
    d3: I,
}

impl J {
    fn var(t: I) -> J {
        J { v: t, d1: I::point(1.0), d2: I::point(0.0), d3: I::point(0.0) }
    }
    fn c(x: I) -> J {
        J { v: x, d1: I::point(0.0), d2: I::point(0.0), d3: I::point(0.0) }
    }
    fn add(self, o: J) -> J {
        J { v: self.v + o.v, d1: self.d1 + o.d1, d2: self.d2 + o.d2, d3: self.d3 + o.d3 }
    }
    fn sub(self, o: J) -> J {
        J { v: self.v - o.v, d1: self.d1 - o.d1, d2: self.d2 - o.d2, d3: self.d3 - o.d3 }
    }
    /// Leibniz's rule.
    fn mul(self, o: J) -> J {
        J {
            v: self.v * o.v,
            d1: self.d1 * o.v + self.v * o.d1,
            d2: self.d2 * o.v + self.d1 * o.d1 * 2.0 + self.v * o.d2,
            d3: self.d3 * o.v + self.d2 * o.d1 * 3.0 + self.d1 * o.d2 * 3.0 + self.v * o.d3,
        }
    }
    /// 1/f: its derivatives −f'/f², (2f'² − f f'')/f³ and (−6f'³ + 6f f' f'' − f² f''')/f⁴.
    fn recip(self) -> J {
        let g = I::point(1.0) / self.v;
        let g2 = g * g;
        let (a, b, c) = (self.d1, self.d2, self.d3);
        J {
            v: g,
            d1: -(a * g2),
            d2: (a.sqr() * g * 2.0 - b) * g2,
            d3: (-(a.powi(3) * g2 * 6.0) + a * b * g * 6.0 - c) * g2,
        }
    }
    fn div(self, o: J) -> J {
        self.mul(o.recip())
    }
}

/// The unit section, α = (1−w)(1−t)⁴/(Q(t)(1−w+wt)), β = (1−w)t⁴/(Q(t)(1−wt)), Q = 1 − 3t + 3t²
/// (the order-2 blend section f(1 − p₁), (1 − f)p₂ with f = (1−t)³/((1−t)³ + t³) written over
/// one denominator): α and β with two derivatives at t.
fn section(t: I, w: I) -> (J, J) {
    if t.lo == t.hi {
        return section_naive(t, w);
    }
    // The mean-value form: each derivative at the middle, plus the next one over the whole
    // interval times the distance from the middle (tighter than the plain extension, whose
    // dependency widens a narrow interval of t many times over); intersected with the plain
    // one.
    let m = t.mid();
    let (am, bm) = section_naive(I::point(m), w);
    let (aw, bw) = section_naive(t, w);
    let d = t - I::point(m);
    let mv = |mid: J, wide: J| J {
        v: meet(mid.v + wide.d1 * d, wide.v),
        d1: meet(mid.d1 + wide.d2 * d, wide.d1),
        d2: meet(mid.d2 + wide.d3 * d, wide.d2),
        d3: wide.d3,
    };
    (mv(am, aw), mv(bm, bw))
}

/// The intersection of two intervals that both hold a value.
fn meet(a: I, b: I) -> I {
    let (lo, hi) = (a.lo.max(b.lo), a.hi.min(b.hi));
    if lo <= hi { I { lo, hi } } else { a }
}

fn section_naive(t: I, w: I) -> (J, J) {
    let t = J::var(t);
    let one = J::c(I::point(1.0));
    let s = one.sub(t);
    let k = J::c(1.0 - w);
    let q = one.sub(t.mul(J::c(I::point(3.0)))).add(t.mul(t).mul(J::c(I::point(3.0))));
    let s4 = s.mul(s).mul(s).mul(s);
    let t4 = t.mul(t).mul(t).mul(t);
    let alpha = k.mul(s4).div(q.mul(k.add(t.mul(J::c(w)))));
    let beta = k.mul(t4).div(q.mul(one.sub(t.mul(J::c(w)))));
    (alpha, beta)
}

/// h(t; s) = s·α'(t) + β'(t), increasing in t, zero where the section's tangent is square to the
/// normal (s, 1).
fn slope(t: I, s: I, w: I) -> I {
    let (a, b) = section(t, w);
    s * a.d1 + b.d1
}

/// An interval holding t*(s), the parameter where the section's tangent is square to (s, 1),
/// for 0 < s ≤ 1 (so t* ≤ ½), and so for every s in [s.lo, s.hi] (t* rises with s): found by
/// float bisection, each end then certified by the sign of h.
fn root(s: I, w: I) -> I {
    // (The guess by float bisection on h evaluated at the interval's middle; only the
    // certification below is rigorous.)
    let wm = I::point(w.mid());
    let find = |s: f64| {
        // Newton's method on h kept within a shrinking bracket.
        let (mut lo, mut hi, mut t) = (0.0f64, 0.5f64 + 1e-9, 0.25f64);
        for _ in 0..100 {
            let (a, b) = section_naive(I::point(t), wm);
            let h = (I::point(s) * a.d1 + b.d1).mid();
            if h < 0.0 {
                lo = t;
            } else {
                hi = t;
            }
            let dh = (I::point(s) * a.d2 + b.d2).mid();
            let mut next = t - h / dh;
            if !(next > lo && next < hi) {
                next = 0.5 * (lo + hi);
            }
            if (next - t).abs() <= 1e-16 * t.max(1e-300) || hi - lo <= 1e-17 {
                return next;
            }
            t = next;
        }
        t
    };
    let mut lo = find(s.lo);
    let mut step = 1e-15f64.max(lo * 1e-14);
    while !(slope(I::point(lo), I::point(s.lo), w).hi < 0.0) {
        lo -= step;
        step *= 2.0;
        assert!(lo > 0.0, "no certified lower root for {s:?}");
    }
    let mut hi = find(s.hi);
    let mut step = 1e-15f64.max(hi * 1e-14);
    // (h is increasing in t and falls as s rises: the root for s.hi bounds every root above.)
    while !(slope(I::point(hi), I::point(s.hi), w).lo > 0.0) {
        hi += step;
        step *= 2.0;
        assert!(hi < 1.0, "no certified upper root for {s:?}");
    }
    I::new(lo, hi)
}

/// f(s) = p(s, 1), f'(s) and F(s) = −f''(s) for 0 < s (an interval), by the section at t*(s):
/// f = sα + β (increasing in s: bounded by its values at the ends), f' = α (falling in s),
/// F = −α'/S' with S = −β'/α', that is F = α'³/(β''α' − β'α''). For s > 1, by the section's
/// mirror: f(s) = s f(1/s), f'(s) = f(1/s) − f'(1/s)/s, F(s) = F(1/s)/s³.
#[derive(Clone, Copy, Debug)]
struct Fs {
    f: I,
    fp: I,
    big_f: I,
}

fn fs(s: I, w: I) -> Fs {
    if !(s.lo > 0.0) || !s.hi.is_finite() {
        return Fs { f: ALL, fp: ALL, big_f: ALL };
    }
    if s.hi <= 1.0 {
        let (r1, r2) = (root(I::point(s.lo), w), root(I::point(s.hi), w));
        let at = |r: I, x: f64| {
            let (a, b) = section(r, w);
            I::point(x) * a.v + b.v
        };
        let f = I { lo: at(r1, s.lo).lo, hi: at(r2, s.hi).hi };
        let (a1, _) = section(r1, w);
        let (a2, _) = section(r2, w);
        let fp = I { lo: a2.v.lo, hi: a1.v.hi };
        let (a, b) = section(r1.hull(r2), w);
        let big_f = a.d1.powi(3) / (b.d2 * a.d1 - b.d1 * a.d2);
        Fs { f, fp, big_f }
    } else if s.lo >= 1.0 {
        // (1/s ≤ 1 exactly; its rounded interval is cut back to that.)
        let q = I::point(1.0) / s;
        let inv = I { lo: q.lo, hi: q.hi.min(1.0) };
        let m = fs(inv, w);
        Fs { f: s * m.f, fp: m.f - m.fp * inv, big_f: m.big_f * inv.powi(3) }
    } else {
        let (a, b) = (fs(I::new(s.lo, 1.0), w), fs(I::new(1.0, s.hi), w));
        Fs { f: a.f.hull(b.f), fp: a.fp.hull(b.fp), big_f: a.big_f.hull(b.big_f) }
    }
}

/// A quantity over a box of (b, r): an interval holding its values and intervals holding its
/// two partial derivatives (forward differentiation, each rule applied in intervals).
#[derive(Clone, Copy, Debug)]
struct D {
    v: I,
    g: [I; 2],
}

impl D {
    fn c(x: I) -> D {
        D { v: x, g: [I::point(0.0); 2] }
    }
    fn add(self, o: D) -> D {
        D { v: self.v + o.v, g: [self.g[0] + o.g[0], self.g[1] + o.g[1]] }
    }
    fn sub(self, o: D) -> D {
        D { v: self.v - o.v, g: [self.g[0] - o.g[0], self.g[1] - o.g[1]] }
    }
    fn neg(self) -> D {
        D { v: -self.v, g: [-self.g[0], -self.g[1]] }
    }
    fn mul(self, o: D) -> D {
        D { v: self.v * o.v, g: [self.g[0] * o.v + self.v * o.g[0], self.g[1] * o.v + self.v * o.g[1]] }
    }
    fn scale(self, k: I) -> D {
        D { v: self.v * k, g: [self.g[0] * k, self.g[1] * k] }
    }
    fn recip(self) -> D {
        let inv = I::point(1.0) / self.v;
        let k = -(inv * inv);
        D { v: inv, g: [self.g[0] * k, self.g[1] * k] }
    }
    fn div(self, o: D) -> D {
        self.mul(o.recip())
    }
    fn sqr(self) -> D {
        D { v: self.v.sqr(), g: [self.g[0] * self.v * 2.0, self.g[1] * self.v * 2.0] }
    }
    fn cube(self) -> D {
        self.sqr().mul(self)
    }
    fn hull(self, o: D) -> D {
        D { v: self.v.hull(o.v), g: [self.g[0].hull(o.g[0]), self.g[1].hull(o.g[1])] }
    }
}

/// f, f' and F (as `D`) of an argument given as a `D`: values as [`fs`] bounds them over the
/// argument's interval, derivatives by the chain rule with f' = α, f'' = −F and
/// F' = dF/ds = (dF/dt)/(dS/dt) over the parameters of the interval, where with
/// K = β''α' − β'α'': F = α'³/K, dS/dt = −K/α'², dF/dt = (3α'²α''K − α'³K')/K²,
/// K' = β'''α' − β'α'''.
fn fs_d(s: D, w: I) -> (D, D, D) {
    let all = D { v: ALL, g: [ALL; 2] };
    if !(s.v.lo > 0.0) || !s.v.hi.is_finite() {
        return (all, all, all);
    }
    if s.v.hi <= 1.0 {
        let base = fs(s.v, w);
        let t = root(I::point(s.v.lo), w).hull(root(I::point(s.v.hi), w));
        let (a, b) = section(t, w);
        let k = b.d2 * a.d1 - b.d1 * a.d2;
        let k1 = b.d3 * a.d1 - b.d1 * a.d3;
        let a1 = a.d1;
        let dfdt = (a1.sqr() * a.d2 * k * 3.0 - a1.powi(3) * k1) / k.sqr();
        let dsdt = -(k / a1.sqr());
        let fprime = dfdt / dsdt;
        let f = D { v: base.f, g: [s.g[0] * base.fp, s.g[1] * base.fp] };
        let fp = D { v: base.fp, g: [-(s.g[0] * base.big_f), -(s.g[1] * base.big_f)] };
        let big_f = D { v: base.big_f, g: [s.g[0] * fprime, s.g[1] * fprime] };
        (f, fp, big_f)
    } else if s.v.lo >= 1.0 {
        // f(s) = s f(1/s), f'(s) = f(1/s) − f'(1/s)/s, F(s) = F(1/s)/s³ (1/s ≤ 1 exactly: its
        // rounded interval cut back to that).
        let mut inv = s.recip();
        inv.v = I { lo: inv.v.lo, hi: inv.v.hi.min(1.0) };
        let (f, fp, big_f) = fs_d(inv, w);
        (s.mul(f), f.sub(fp.mul(inv)), big_f.mul(inv.cube()))
    } else {
        let below = D { v: I { lo: s.v.lo, hi: 1.0 }, g: s.g };
        let above = D { v: I { lo: 1.0, hi: s.v.hi }, g: s.g };
        let (x, y) = (fs_d(below, w), fs_d(above, w));
        (x.0.hull(y.0), x.1.hull(y.1), x.2.hull(y.2))
    }
}

/// M's entries (aa, ab, bb) as `D` over the box b in `b`, r = a/b in `r` (0 < r ≤ 1,
/// 0 < b ≤ 1), with their partial derivatives in b and r.
fn matrix_d(b: D, r: D, w: I) -> [D; 3] {
    let a = r.mul(b);
    let (fa, fpa, big_fa) = fs_d(a, w);
    let (fb, fpb, big_fb) = fs_d(b, w);
    let (fr, fpr, big_fr) = fs_d(r, w);
    let x = b.mul(fr);
    let (_, fpx, big_fx) = fs_d(x, w);
    let (y, z) = (fb, fa);
    let (_, fpya, big_fya) = fs_d(y.div(a), w);
    let (_, fpzb, big_fzb) = fs_d(z.div(b), w);
    let kab = big_fr.div(b.cube());
    let kya = big_fya.div(a.cube());
    let kzb = big_fzb.div(b.cube());
    let third = I::point(1.0) / I::point(3.0);
    let one = D::c(I::point(1.0));
    let (xa, xb) = (fpr, fr.sub(r.mul(fpr)));
    let q = [b, a.neg()];
    let v = [y, a.mul(fpb).neg()];
    let vp = [b.mul(fpa), z.neg()];
    let outer = |u: [D; 2], k: D| [k.mul(u[0].sqr()), k.mul(u[0]).mul(u[1]), k.mul(u[1].sqr())];
    let c1 = kab.mul(one.sub(fpx.scale(third)));
    let pq = outer(q, c1);
    let ea = big_fa.mul(one.sub(fpzb.scale(third)));
    let eb = big_fb.mul(one.sub(fpya.scale(third)));
    let gx = outer([xa, xb], big_fx.scale(third));
    let gv = outer(v, kya.scale(third));
    let gvp = outer(vp, kzb.scale(third));
    [
        pq[0].add(ea).sub(gx[0]).sub(gv[0]).sub(gvp[0]),
        pq[1].sub(gx[1]).sub(gv[1]).sub(gvp[1]),
        pq[2].add(eb).sub(gx[2]).sub(gv[2]).sub(gvp[2]),
    ]
}

/// M's entries over the box by the mean-value form: their values at the box's middle (a thin
/// evaluation) plus their derivatives over the box times the box's spread about the middle,
/// intersected with the plain evaluation over the box.
fn matrix(b: I, r: I, w: I) -> [I; 3] {
    let var = |x: I, k: usize| {
        let mut g = [I::point(0.0); 2];
        g[k] = I::point(1.0);
        D { v: x, g }
    };
    let wide = matrix_d(var(b, 0), var(r, 1), w);
    if b.lo == b.hi && r.lo == r.hi {
        return [wide[0].v, wide[1].v, wide[2].v];
    }
    let (bm, rm) = (b.mid(), r.mid());
    let mid = matrix_d(var(I::point(bm), 0), var(I::point(rm), 1), w);
    let (db, dr) = (b - I::point(bm), r - I::point(rm));
    std::array::from_fn(|k| meet(mid[k].v + wide[k].g[0] * db + wide[k].g[1] * dr, wide[k].v))
}

/// Whether M is certified positive definite over the box: aa > 0 and aa·bb − ab² > 0, after
/// scaling the rows and columns by d = (√|aa|, √|bb|)⁻¹ at the box's middle (a congruence,
/// which keeps definiteness and narrows the determinant's interval).
fn certified(b: I, r: I, w: I) -> bool {
    let m = matrix(b, r, w);
    if !(m.iter().all(|e| e.lo.is_finite() && e.hi.is_finite()) && m[0].lo > 0.0 && m[2].lo > 0.0) {
        return false;
    }
    let (da, db) = (1.0 / m[0].mid().sqrt(), 1.0 / m[2].mid().sqrt());
    let (aa, ab, bb) = (m[0] * (da * da), m[1] * (da * db), m[2] * (db * db));
    let det = aa * bb - ab.sqr();
    det.lo.is_finite() && det.lo > 0.0
}

/// Box (b, r) in log coordinates, split until certified (depth at most `max`).
fn cover(b: (f64, f64), r: (f64, f64), w: I, depth: usize, max: usize, stats: &mut (usize, usize, Vec<((f64, f64), (f64, f64))>)) {
    let bi = I::new(b.0.exp(), if b.1 == 0.0 { 1.0 } else { b.1.exp() });
    // (The box's ends as floats; adjacent boxes share them, the region's are its own.)
    let ri = I::new(r.0.exp(), if r.1 == 0.0 { 1.0 } else { r.1.exp().min(1.0) });
    if certified(bi, ri, w) {
        stats.0 += 1;
        return;
    }
    if depth >= max {
        stats.1 += 1;
        stats.2.push((b, r));
        return;
    }
    let (bm, rm) = (0.5 * (b.0 + b.1), 0.5 * (r.0 + r.1));
    for bb in [(b.0, bm), (bm, b.1)] {
        for rr in [(r.0, rm), (rm, r.1)] {
            cover(bb, rr, w, depth + 1, max, stats);
        }
    }
}

// ---- Floating-point counterparts, for `check` and `counterexample`. ----

fn f_float(s: f64, w: f64) -> (f64, f64, f64) {
    let m = fs(I::point(s), I::new(w, w));
    (m.f.mid(), m.fp.mid(), m.big_f.mid())
}

/// p(x, y) = y f(x/y) for x, y ≥ 0.
fn p_float(x: f64, y: f64, w: f64) -> f64 {
    if x <= 0.0 || y <= 0.0 {
        return 0.0;
    }
    y * f_float(x / y, w).0
}

/// H(a, b) on the slice c = 1, from p alone.
fn h_float(a: f64, b: f64, w: f64) -> f64 {
    let x = p_float(a, b, w);
    let (fa, fb) = (p_float(a, 1.0, w), p_float(b, 1.0, w));
    x + fa + fb - (p_float(x, 1.0, w) + p_float(fb, a, w) + p_float(fa, b, w)) / 3.0
}

/// The per-turn term G = Σp − p(p(a, b), c), on the slice c = 1.
fn g_float(a: f64, b: f64, w: f64) -> f64 {
    let x = p_float(a, b, w);
    x + p_float(a, 1.0, w) + p_float(b, 1.0, w) - p_float(x, 1.0, w)
}

fn hessian_fd(h: &dyn Fn(f64, f64) -> f64, a: f64, b: f64, e: f64) -> [f64; 3] {
    let (ea, eb) = (e * a, e * b);
    let aa = (h(a + ea, b) - 2.0 * h(a, b) + h(a - ea, b)) / (ea * ea);
    let bb = (h(a, b + eb) - 2.0 * h(a, b) + h(a, b - eb)) / (eb * eb);
    let ab = (h(a + ea, b + eb) - h(a + ea, b - eb) - h(a - ea, b + eb) + h(a - ea, b - eb)) / (4.0 * ea * eb);
    [-aa, -ab, -bb]
}

// ---- The scaled matrix N on the closed rectangle (β, ρ) ∈ [0, 1]², b = β³, r = a/b = ρ³. ----

/// The cube root of a positive interval (each end's root moved out by four ulps, more than
/// the error of `f64::cbrt`).
fn cbrt(x: I) -> I {
    if !(x.lo >= 0.0) || !x.hi.is_finite() {
        return ALL;
    }
    let widen = |v: f64, k: i32, f: fn(f64) -> f64| (0..k).fold(v, |v, _| f(v));
    I { lo: widen(x.lo.cbrt(), 4, down).max(0.0), hi: widen(x.hi.cbrt(), 4, up) }
}

/// The ingredient functions of the section at its parameter t (orders above those used left
/// to `ALL`): with A = c/(Q(c + wt)), B = c/(Q(1 − wt)), c = 1 − w, A₁ = 4A − (1−t)A',
/// B₁ = 4B + tB', σ = B₁/((1−t)³A₁) (so S = t³σ), N₁ = 1 − 3t + (4−w)t² − (1−w)t³,
/// γ = σN₁/(Q(c + wt)) − B, S₁ = 3σ + tσ', E = (1−t)³A₁/S₁:
///
/// s = t³σ, G = g(s)/s^{4/3} = γ/σ^{4/3}, G₁ = g'(s)/s^{1/3} = N₁/(Q(c + wt)σ^{1/3}),
/// Fs = F(s)s^{2/3} = Eσ^{2/3}, Fp = F'(s)s^{5/3} = (tE' − 2E)σ^{5/3}/S₁.
#[derive(Clone, Copy, Debug)]
struct Ing {
    s: I,
    g: I,
    g1: I,
    fs: I,
    fp: I,
    /// The least of A₁, B₁ and S₁ over the cell (all must be positive: α' < 0 < β' and S
    /// rising).
    least: f64,
}

fn shift(j: J) -> J {
    J { v: j.d1, d1: j.d2, d2: j.d3, d3: ALL }
}

fn ingredients(t: I, w: I) -> Ing {
    let tj = J::var(t);
    let one = J::c(I::point(1.0));
    let c = J::c(1.0 - w);
    let wj = J::c(w);
    let three = J::c(I::point(3.0));
    let q = one.sub(tj.mul(three)).add(tj.mul(tj).mul(three));
    let l1 = c.add(tj.mul(wj));
    let l2 = one.sub(tj.mul(wj));
    let a = c.div(q.mul(l1));
    let b = c.div(q.mul(l2));
    let s1t = one.sub(tj);
    let a1 = a.mul(J::c(I::point(4.0))).sub(s1t.mul(shift(a)));
    let b1 = b.mul(J::c(I::point(4.0))).add(tj.mul(shift(b)));
    let cube = s1t.mul(s1t).mul(s1t);
    let sigma = b1.div(cube.mul(a1));
    let n1 = one.sub(tj.mul(three)).add(tj.mul(tj).mul(J::c(4.0 - w))).sub(tj.mul(tj).mul(tj).mul(c));
    let gamma = sigma.mul(n1).div(q.mul(l1)).sub(b);
    let big_s1 = sigma.mul(three).add(tj.mul(shift(sigma)));
    let e = cube.mul(a1).div(big_s1);
    let s3 = cbrt(sigma.v);
    let s = t.powi(3) * sigma.v;
    Ing {
        s,
        g: gamma.v / (s3.powi(4)),
        g1: n1.v / (q.v * l1.v * s3),
        fs: e.v * s3.sqr(),
        fp: (t * e.d1 - e.v * 2.0) * s3.powi(5) / big_s1.v,
        least: a1.v.lo.min(b1.v.lo).min(big_s1.v.lo),
    }
}

/// The ingredients tabulated over the cells [t_k, t_k+1] of an even grid of t from 0 to `TMAX`:
/// each cell's s from the (rising) S at its ends, the rest over the cell.
struct Table {
    nodes: Vec<I>,
    cells: Vec<Ing>,
}

const TMAX: f64 = 0.75;

impl Table {
    fn new(n: usize, w: I) -> Table {
        let ts: Vec<f64> = (0..=n).map(|k| TMAX * k as f64 / n as f64).collect();
        let nodes: Vec<I> = ts.iter().map(|&t| if t == 0.0 { I::point(0.0) } else { ingredients(I::point(t), w).s }).collect();
        let cells: Vec<Ing> = (0..n)
            .map(|k| {
                let mut c = ingredients(I::new(ts[k], ts[k + 1]), w);
                c.s = I { lo: nodes[k].lo, hi: nodes[k + 1].hi };
                c
            })
            .collect();
        // (Lemma 1's conditions, on t in [0, TMAX] (TMAX ≥ ½, so by the mirror on all of
        // [0, 1]): A₁, B₁ > 0, so α' < 0 < β', and S₁ > 0, so S rises within each cell; and the
        // nodes' intervals of S ordered. Else the lookup below would be unsound.)
        for (k, c) in cells.iter().enumerate() {
            assert!(c.least > 0.0, "A₁, B₁ or S₁ not certified positive on cell {k}");
        }
        for k in 1..=n {
            assert!(nodes[k].lo > nodes[k - 1].hi, "S not certified rising at node {k}");
        }
        Table { nodes, cells }
    }

    /// The ingredients over an interval of s (0 ≤ s): the hull of the cells whose s may lie in it.
    fn at(&self, s: I) -> Ing {
        let all = Ing { s: ALL, g: ALL, g1: ALL, fs: ALL, fp: ALL, least: f64::NAN };
        if !(s.lo >= 0.0) || !(s.hi <= self.nodes[self.nodes.len() - 1].lo) {
            return all;
        }
        // Cells k with nodes[k+1].hi ≥ s.lo and nodes[k].lo ≤ s.hi.
        let first = self.nodes.partition_point(|n| n.hi < s.lo).saturating_sub(1);
        let last = self.nodes.partition_point(|n| n.lo <= s.hi).min(self.cells.len());
        let mut out = self.cells[first];
        for c in &self.cells[first..last.max(first + 1)] {
            out = Ing { s: out.s.hull(c.s), g: out.g.hull(c.g), g1: out.g1.hull(c.g1), fs: out.fs.hull(c.fs), fp: out.fp.hull(c.fp), least: out.least.min(c.least) };
        }
        out
    }
}

/// N = β²·diag(ρ, 1)·M·diag(ρ, 1) over the box (β, ρ), in its form free of singular factors (see
/// the proof's Lemma 5): entries (aa, ab, bb).
fn scaled(beta: I, rho: I, tab: &Table) -> [I; 3] {
    let third = I::point(1.0) / I::point(3.0);
    let one = I::point(1.0);
    let (r, b) = (rho.powi(3), beta.powi(3));
    let a = r * b;
    let (ir, ia, ib) = (tab.at(r), tab.at(a), tab.at(b));
    let xi3 = one - rho * ir.g;
    let xi = cbrt(xi3);
    let x = b * r * xi3;
    let ix = tab.at(x);
    let za = one - rho * beta * ia.g;
    let z = r * za;
    let iz = tab.at(z);
    let fpz = one - rho * cbrt(za) * iz.g1;
    let yb = one - beta * ib.g;
    let u = r / yb;
    let iu = tab.at(u);
    let u43 = (rho / cbrt(yb)).powi(4);
    let fpu = u43 * (iu.g1 - iu.g);
    let eta = tab.at(r.hull(u));
    let zeta43 = I::point(1.0).hull((one / cbrt(yb)).powi(4));
    let mut k5 = zeta43 * (eta.fs * 3.0 + eta.fp) * ib.g / yb;
    // (Where β is not small, also the difference quotient itself: K₅ = ((1−βG_b)^{−7/3}Fs(u) −
    // Fs(r))/β, and the value lies in both.)
    if beta.lo > 0.0 {
        k5 = meet(k5, ((one / cbrt(yb)).powi(7) * iu.fs - ir.fs) / beta);
    }
    let theta = tab.at(z.hull(r));
    let zeta53 = I::point(1.0).hull((one / cbrt(za)).powi(5));
    let mut k6 = -(theta.fp * zeta53 * ia.g);
    // (Likewise K₆ = (Fs(z)(1−ρβG_a)^{−2/3} − Fs(r))/(ρβ).)
    if (rho * beta).lo > 0.0 {
        k6 = meet(k6, (iz.fs / cbrt(za).sqr() - ir.fs) / (rho * beta));
    }
    let (g1r_gr, oneg1r) = (ir.g1 - ir.g, one - rho * ir.g1);
    let t1 = ir.fs * ix.g1 * xi * third;
    let t4 = ix.fs / xi.sqr() * third;
    let ybg1 = one - beta * ib.g1;
    let (zag1, zag) = (one - rho * beta * ia.g1, za);
    let aa = t1 * rho + ia.fs * (one - fpz * third) - t4 * oneg1r.sqr() - k5 * yb.sqr() * third
        - ir.fs * (ib.g * -2.0 + beta * ib.g.sqr()) * third
        - k6 * rho * zag1.sqr() * third
        - ir.fs * (rho * ia.g1 * -2.0 + rho.sqr() * beta * ia.g1.sqr()) * third;
    let ab = -(t1 * rho.powi(3)) - t4 * oneg1r * rho.powi(3) * g1r_gr
        + k5 * rho.sqr() * yb * ybg1 * third
        - ir.fs * rho.sqr() * (ib.g + ib.g1 - beta * ib.g * ib.g1) * third
        + k6 * rho.powi(3) * zag1 * zag * third
        - ir.fs * (rho.powi(3) * (ia.g1 + ia.g) - rho.powi(4) * beta * ia.g1 * ia.g) * third;
    let bb = t1 * rho.powi(5) + ib.fs * (one - fpu * third) - t4 * rho.powi(6) * g1r_gr.sqr()
        - k5 * rho.powi(4) * ybg1.sqr() * third
        - ir.fs * rho.powi(4) * (ib.g1 * -2.0 + beta * ib.g1.sqr()) * third
        - k6 * rho.powi(5) * zag.sqr() * third
        - ir.fs * (rho.powi(5) * ia.g * -2.0 + rho.powi(6) * beta * ia.g.sqr()) * third;
    [aa, ab, bb]
}

fn pd(m: [I; 3]) -> bool {
    if !(m.iter().all(|e| e.lo.is_finite() && e.hi.is_finite()) && m[0].lo > 0.0 && m[2].lo > 0.0) {
        return false;
    }
    let (da, db) = (1.0 / m[0].mid().sqrt(), 1.0 / m[2].mid().sqrt());
    let (aa, ab, bb) = (m[0] * (da * da), m[1] * (da * db), m[2] * (db * db));
    let det = aa * bb - ab.sqr();
    det.lo.is_finite() && det.lo > 0.0
}

/// The box (β, ρ) split until N is certified positive definite on each piece.
fn cover_scaled(beta: (f64, f64), rho: (f64, f64), tab: &Table, depth: usize, max: usize, stats: &mut (usize, usize, Vec<((f64, f64), (f64, f64))>)) {
    if pd(scaled(I::new(beta.0, beta.1), I::new(rho.0, rho.1), tab)) {
        stats.0 += 1;
        return;
    }
    if depth >= max {
        stats.1 += 1;
        stats.2.push((beta, rho));
        return;
    }
    let (bm, rm) = (0.5 * (beta.0 + beta.1), 0.5 * (rho.0 + rho.1));
    for bb in [(beta.0, bm), (bm, beta.1)] {
        for rr in [(rho.0, rm), (rm, rho.1)] {
            cover_scaled(bb, rr, tab, depth + 1, max, stats);
        }
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let w = I::new(W_LO, W_HI);
    let wf = -0.5252124045775235;
    match args.get(1).map(|s| s.as_str()) {
        Some("check") => {
            for (a, b) in [(0.3, 0.6), (0.05, 0.7), (0.5, 0.5), (0.9, 1.0), (0.01, 0.02), (0.2, 0.9)] {
                let m = matrix(I::point(b), I::point(a / b), I::new(wf, wf));
                let fd = hessian_fd(&|x, y| h_float(x, y, wf), a, b, 1e-4);
                println!("a {a} b {b}: formula ({:.6}, {:.6}, {:.6})  differences ({:.6}, {:.6}, {:.6})", m[0].mid(), m[1].mid(), m[2].mid(), fd[0], fd[1], fd[2]);
            }
        }
        Some("certify") => {
            // certify B0 B1 R0 R1 [N MAXDEPTH THREADS]: b in [B0, B1], r = a/b in [R0, R1], an
            // N × N grid of boxes even in log b and log r, each split as needed.
            let v: Vec<f64> = args[2..6].iter().map(|x| x.parse().unwrap()).collect();
            let n: usize = args.get(6).map(|s| s.parse().unwrap()).unwrap_or(16);
            let max: usize = args.get(7).map(|s| s.parse().unwrap()).unwrap_or(10);
            let threads: usize = args.get(8).map(|s| s.parse().unwrap()).unwrap_or(8);
            let (b0, b1, r0, r1) = (v[0].ln(), v[1].ln(), v[2].ln(), v[3].ln());
            let edge = |lo: f64, hi: f64, k: usize| if k == 0 { lo } else if k == n { hi } else { lo + (hi - lo) * k as f64 / n as f64 };
            let next = std::sync::atomic::AtomicUsize::new(0);
            let results = std::sync::Mutex::new((0usize, 0usize, Vec::new()));
            std::thread::scope(|scope| {
                for _ in 0..threads {
                    scope.spawn(|| loop {
                        let k = next.fetch_add(1, std::sync::atomic::Ordering::SeqCst);
                        if k >= n * n {
                            break;
                        }
                        let (i, j) = (k / n, k % n);
                        let mut stats = (0usize, 0usize, Vec::new());
                        cover((edge(b0, b1, i), edge(b0, b1, i + 1)), (edge(r0, r1, j), edge(r0, r1, j + 1)), w, 0, max, &mut stats);
                        let mut all = results.lock().unwrap();
                        all.0 += stats.0;
                        all.1 += stats.1;
                        all.2.extend(stats.2);
                    });
                }
            });
            let stats = results.into_inner().unwrap();
            println!("b in [{}, {}], a/b in [{}, {}]: {} boxes certified, {} not", v[0], v[1], v[2], v[3], stats.0, stats.1);
            for (b, r) in stats.2.iter().take(20) {
                println!("  not certified: b in [{:.4e}, {:.4e}], r in [{:.4e}, {:.4e}]", b.0.exp(), b.1.exp(), r.0.exp(), r.1.exp());
            }
            if stats.1 > 0 {
                std::process::exit(1);
            }
        }
        Some("ncheck") => {
            // N by its singularity-free form against β²·diag(ρ,1)·M·diag(ρ,1) from M's plain form.
            let tab = Table::new(1 << 16, w);
            for (beta, rho) in [(0.8, 0.7), (0.5, 0.9), (0.3, 0.4), (0.9, 0.2), (0.6, 0.6), (0.99, 0.99)] {
                let (bb, rr) = (beta * beta * beta, rho * rho * rho);
                let m = matrix(I::point(bb), I::point(rr), w);
                let k = beta * beta;
                let want = [m[0].mid() * k * rho * rho, m[1].mid() * k * rho, m[2].mid() * k];
                let n = scaled(I::point(beta), I::point(rho), &tab);
                println!("β {beta} ρ {rho}: N {:?}\n   want {want:?}", n.map(|x| (x.lo, x.hi)));
            }
        }
        Some("scaled") => {
            // scaled B0 B1 R0 R1 [N MAXDEPTH THREADS CELLS]: certify N positive definite over the
            // rectangle β in [B0, B1], ρ in [R0, R1] (all of [0, 1]² for the lemma).
            let v: Vec<f64> = args[2..6].iter().map(|x| x.parse().unwrap()).collect();
            let n: usize = args.get(6).map(|s| s.parse().unwrap()).unwrap_or(64);
            let max: usize = args.get(7).map(|s| s.parse().unwrap()).unwrap_or(12);
            let threads: usize = args.get(8).map(|s| s.parse().unwrap()).unwrap_or(8);
            let cells: usize = args.get(9).map(|s| s.parse().unwrap()).unwrap_or(1 << 18);
            let tab = Table::new(cells, w);
            let edge = |lo: f64, hi: f64, k: usize| if k == 0 { lo } else if k == n { hi } else { lo + (hi - lo) * k as f64 / n as f64 };
            let next = std::sync::atomic::AtomicUsize::new(0);
            let results = std::sync::Mutex::new((0usize, 0usize, Vec::new()));
            std::thread::scope(|scope| {
                for _ in 0..threads {
                    scope.spawn(|| loop {
                        let k = next.fetch_add(1, std::sync::atomic::Ordering::SeqCst);
                        if k >= n * n {
                            break;
                        }
                        let (i, j) = (k / n, k % n);
                        let mut stats = (0usize, 0usize, Vec::new());
                        cover_scaled((edge(v[0], v[1], i), edge(v[0], v[1], i + 1)), (edge(v[2], v[3], j), edge(v[2], v[3], j + 1)), &tab, 0, max, &mut stats);
                        let mut all = results.lock().unwrap();
                        all.0 += stats.0;
                        all.1 += stats.1;
                        all.2.extend(stats.2);
                    });
                }
            });
            let stats = results.into_inner().unwrap();
            println!("β in [{}, {}], ρ in [{}, {}]: {} boxes certified, {} not", v[0], v[1], v[2], v[3], stats.0, stats.1);
            for (b, r) in stats.2.iter().take(20) {
                println!("  not certified: β in [{:.6}, {:.6}], ρ in [{:.6}, {:.6}]", b.0, b.1, r.0, r.1);
            }
            if stats.1 > 0 {
                std::process::exit(1);
            }
        }
        Some("counterexample") => {
            // Σp − p(p(a,b),c) at (a, b, c) = (27, 2, 1): its Hessian on the slice has a
            // positive eigenvalue (not concave). Certified by the interval matrix of the
            // per-turn term: G's minus-Hessian is that of H with ⅓ → 1 on the c-turn and no
            // other turns; computed here by differences in float as evidence, then by the
            // interval formula below.
            let (a, b) = (27.0, 2.0);
            let fd = hessian_fd(&|x, y| g_float(x, y, wf), a, b, 1e-4);
            let det = fd[0] * fd[2] - fd[1] * fd[1];
            println!("differences of G at (27, 2, 1): minus-Hessian ({:.6}, {:.6}, {:.6}), determinant {:.6e}", fd[0], fd[1], fd[2], det);
            let m = per_turn(I::point(b), I::point(a / b), w);
            let det = m[0] * m[2] - m[1].sqr();
            println!("interval minus-Hessian of G: aa {:?} ab {:?} bb {:?}; determinant {:?}", m[0], m[1], m[2], det);
            assert!(det.hi < 0.0, "not certified");
            println!("certified: the determinant is negative, G is not concave there");
        }
        _ => eprintln!("usage: g2_concavity scaled B0 B1 R0 R1 [N MAXDEPTH THREADS CELLS] | check | ncheck | certify B0 B1 R0 R1 [N MAXDEPTH THREADS] | counterexample"),
    }
}

/// Minus the Hessian of the per-turn term G = p(a,b) + f(a) + f(b) − f(X) on the slice c = 1:
/// κ(a,b)(1 − f'(X)) q qᵀ + F(a) e_a e_aᵀ + F(b) e_b e_bᵀ − F(X) ∇X ∇Xᵀ (for a, b ≤ 1 or not:
/// `fs` handles arguments above 1 by the mirror).
fn per_turn(b: I, r: I, w: I) -> [I; 3] {
    let a = r * b;
    let (fa, fb) = (fs(a, w), fs(b, w));
    let fr = if r.hi <= 1.0 { fs(r, w) } else { fs(r, w) };
    let x = b * fr.f;
    let fx = fs(x, w);
    let kab = fr.big_f / b.powi(3);
    let (xa, xb) = (fr.fp, fr.f - r * fr.fp);
    let c1 = kab * (1.0 - fx.fp);
    [
        c1 * b.sqr() + fa.big_f - fx.big_f * xa.sqr(),
        -(c1 * a * b) - fx.big_f * xa * xb,
        c1 * a.sqr() + fb.big_f - fx.big_f * xb.sqr(),
    ]
}
