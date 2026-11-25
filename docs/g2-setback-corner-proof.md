# G2 set-back corners: the concavity lemma

`blend::g2_corner` closes a square corner whose three edges carry curvature-continuous fillets
of tangent offset `r`, set back `s` from the corner, with the convex patch whose (lower) support
function on the octant of inward normals ν = (a, b, c) is

    g = r·Σp + (s − 2r)·Φ,   Σp = p(b, c) + p(c, a) + p(a, b),   Φ = ⅓ Σ_cyc p(p(a, b), c),

where `p` is the support function of the unit section's region (see the module's notes). For
`s ≥ 2r` both summands are concave. For `r < s < 2r` the patch is convex if and only if `g` is
concave, and

    g = (2r − s)·(Σp − Φ) + (s − r)·Σp,

with both coefficients positive, so it is enough that

> **Theorem.** H = Σp − Φ is concave on the closed positive octant, for the section kentron
> uses: the order-2 blend section at its equal-area thumbweight, w in
> [−0.5252124045776, −0.5252124045774] (kentron computes w = −0.5252124045775235). Its
> Hessian is negative definite on the planes square to ν, at every ν inside the octant.

Then `g` is concave for every `s ≥ r`, strictly inside the octant (the patch strictly convex).
At `s = r` the patch degenerates only at the octant's corners, where the faces' flats shrink to
points.

The proof reduces the claim to the positive definiteness of a 2 × 2 matrix of one-variable
functions (Lemmas 1 to 4, by hand), and proves that in interval arithmetic over a closed square
(Lemma 5, by `tools/g2_concavity.rs`). What the machine proves, and what it relies on, is
listed at the end.

## Lemma 1: the section and the function f

With `c = 1 − w`, `Q(t) = 1 − 3t + 3t²`, `A = c/(Q(c + wt))`, `B = c/(Q(1 − wt))`, the order-2
blend section (Hartmann's `f(1 − p₁)`, `(1 − f)p₂` with `f = (1−t)³/((1−t)³ + t³)` over the one
denominator `(1−t)³ + t³ = Q`) is

    α(t) = (1−t)⁴ A(t),   β(t) = t⁴ B(t),   t ∈ [0, 1],

from `(1, 0)` to `(0, 1)`, its own mirror: `α(1 − t) = β(t)`. Its derivatives are

    α' = −(1−t)³ A₁,  A₁ = 4A − (1−t)A',     β' = t³ B₁,  B₁ = 4B + tB'.

*Certified conditions* (by the tool on t ∈ [0, ¾], hence by the mirror on all of [0, 1]):
`A₁ > 0`, `B₁ > 0`, so `α' < 0 < β'` inside (0, 1); and, with `S(t) = −β'/α' = t³σ(t)`,
`σ = B₁/((1−t)³A₁)`, `S' = t² S₁` with `S₁ = 3σ + tσ' > 0`: S rises from 0 to ∞.

`p(x, y) = min_t xα(t) + yβ(t)` (x, y ≥ 0) is the support function of the region (what
`UnitSection::support` computes). For the normal (s, 1), `s > 0`,
`sα'(t) + β'(t) = α'(t)(s − S(t))` is negative before `t* = S⁻¹(s)` and positive after, so the
minimum is at `t*` alone. With `f(s) = p(s, 1)`:

    p(x, y) = y f(x/y),   f(s) = s f(1/s)  (the mirror),
    f(s) = sα(t*) + β(t*),   f'(s) = α(t*),   F(s) := −f''(s) = −α'(t*)/S'(t*) > 0

(the envelope rule twice). In t:

    F = (1−t)³A₁ / (t² S₁),   g(s) := s − f(s) = t⁴γ,  γ = σN₁/(Q(c + wt)) − B,
    g'(s) = 1 − α = t N₁/(Q(c + wt)),   N₁ = 1 − 3t + (4 − w)t² − (1 − w)t³,

the last because `Q(c + wt) − c(1−t)⁴ = t·N₁` exactly. Since `σ(0) > 0`, the *scaled
ingredients*

    G(s) = g(s)/s^{4/3} = γ/σ^{4/3},          G₁(s) = g'(s)/s^{1/3} = N₁/(Q(c + wt)σ^{1/3}),
    Fs(s) = F(s)s^{2/3} = Eσ^{2/3},           Fp(s) = F'(s)s^{5/3} = (tE' − 2E)σ^{5/3}/S₁,

with `E = (1−t)³A₁/S₁`, are bounded and continuous on [0, S(¾)], including s = 0 (t = 0).
(`F = E/t²`, so `dF/dt = (tE' − 2E)/t³` and `F' = (dF/dt)/S'`.)

## Lemma 2: the Hessian on the slice c = 1

On `c = 1`, `H(a, b, 1) = p(a, b) + f(a) + f(b) − ⅓[f(X) + p(f(b), a) + p(f(a), b)]`, with
`X = p(a, b)` (the three terms of Φ: `p(p(a,b),1) = f(X)`, `p(p(b,1), a)`, `p(p(1,a), b)`).

The Hessian of `p` in its two arguments is rank one: `∇²p(x, y) = −κ(x, y)(y, −x)(y, −x)ᵀ`
with `κ(x, y) = F(x/y)/y³` (from `p = y f(x/y)`). The chain rule on each term then gives,
for the Hessian in (a, b),

    −∇²H = M = κ(a,b)(1 − f'(X)/3) q qᵀ
              + F(a)(1 − f'(Z/b)/3) e_a e_aᵀ + F(b)(1 − f'(Y/a)/3) e_b e_bᵀ
              − ⅓[F(X) ∇X∇Xᵀ + κ(Y,a) v vᵀ + κ(Z,b) v'v'ᵀ],

with `Y = f(b)`, `Z = f(a)`, `q = (b, −a)`, `v = (Y, −a f'(b))`, `v' = (b f'(a), −Z)`,
`∇X = (f'(a/b), f(a/b) − (a/b) f'(a/b))`. (For instance `∇²[p(f(b), a)] = −κ(Y,a) v vᵀ
− f'(Y/a) F(b) e_b e_bᵀ`, the second term from `f''(b)` in the chain rule.)

Checked against central differences of H (`g2_concavity check`) and against kentron's own
Hessians of Σp and Φ (`g2_corner::tests::the_slice_hessian_is_the_proof_s_matrix`).

## Lemma 3: reduction to a triangle

H is continuous on the closed octant and twice differentiable inside it. A twice
differentiable function is concave on the open convex cone if and only if its Hessian is
negative semidefinite everywhere, and concavity passes to the closure by continuity. H is of
degree one, so ν spans its Hessian's kernel, and the Hessian is negative semidefinite on R³
if and only if it is on any plane missing ν, such as `{dc = 0}` when c > 0. That restriction is
the Hessian of `(a, b) ↦ H(a, b, c)`.

H is symmetric under every permutation of (a, b, c): Σp is, p being symmetric, and each term
`p(p(a, b), c)` is symmetric in a and b, so their mean over the three turns is symmetric under
all permutations. Permuting to make c the largest coordinate and scaling it to 1, it is
enough that **M is positive definite for 0 < a ≤ b ≤ 1**.

## Lemma 4: a form of M free of singular factors

Put `b = β³`, `a/b = r = ρ³` (β, ρ ∈ (0, 1]) and `N = β²·diag(ρ, 1)·M·diag(ρ, 1)`, a
congruence: N is positive definite exactly when M is.

M is not fit for evaluation near the triangle's sides: its terms grow like `b⁻¹` towards the
corner a = b = 0 (where they cancel to `b^{−2/3}`) and like `a^{−2/3}` towards a = 0. The
cancellation is exact algebra. Write `f = s − g`; then `κ(b, a) = κ(a, b)` (by `F(1/s) = s³F(s)`)
and `v, v' → q` as `g → 0`. Split

    κ(Y,a) v vᵀ = [κ(Y,a) − κ(b,a)] v vᵀ + κ(a,b) (v vᵀ − q qᵀ) + κ(a,b) q qᵀ,

and likewise the v' term. The two `κ q qᵀ` parts cancel the first term's `−⅔κ q qᵀ`, leaving
`(κ/3) g'(X) q qᵀ`. Each remaining difference is written by the mean value theorem:

    F(Y/a) − F(b/a) = Ψ(u) − Ψ(r) = Ψ'(η)(u − r),   Ψ(u) = u³F(u) = F(1/u),
        Ψ'(η) = η^{4/3}(3Fs(η) + Fp(η)),  u = a/Y = r/(1 − βG(b)),  η between r and u;
    F(Z/b) − F(r) = F'(θ)(z − r),   z = Z/b = r(1 − ρβG(a)),  θ between z and r,

and every other quantity through the scaled ingredients, e.g. `g(b)/b = βG(b)`,
`f'(b) = 1 − βG₁(b)`, `X = b r (1 − ρG(r))`, `X^{1/3} = βρξ` with `ξ = (1 − ρG(r))^{1/3}`,
`f'(Y/a) = u^{4/3}(G₁(u) − G(u))`. All powers of β and ρ then cancel to non-negative ones.
With subscripts naming the ingredients' arguments (`G_b = G(b)`, `Fs_r = Fs(r)`, ...),
`yb = 1 − βG_b`, `za = 1 − ρβG_a`,

    K₅ = ζ^{4/3}(3Fs_η + Fp_η) G_b / yb,   ζ = η/r ∈ [1, yb⁻¹]
    K₆ = −Fp_θ ζ'^{−5/3} G_a,             ζ' = θ/r ∈ [za, 1]
    T₁ = Fs_r G₁_X ξ/3,   T₄ = Fs_X/(3ξ²)

    N_aa = T₁ρ + Fs_a(1 − f'(z)/3) − T₄(1 − ρG₁_r)² − ⅓K₅ yb² + ⅓Fs_r(2G_b − βG_b²)
           − ⅓K₆ρ(1 − ρβG₁_a)² + ⅓Fs_r(2ρG₁_a − ρ²βG₁_a²)
    N_ab = −T₁ρ³ − T₄ρ³(1 − ρG₁_r)(G₁_r − G_r) + ⅓K₅ρ² yb(1 − βG₁_b)
           − ⅓Fs_r ρ²(G_b + G₁_b − βG_bG₁_b) + ⅓K₆ρ³(1 − ρβG₁_a) za
           − ⅓Fs_r(ρ³(G₁_a + G_a) − ρ⁴βG₁_aG_a)
    N_bb = T₁ρ⁵ + Fs_b(1 − f'(Y/a)/3) − T₄ρ⁶(G₁_r − G_r)² − ⅓K₅ρ⁴(1 − βG₁_b)²
           + ⅓Fs_r ρ⁴(2G₁_b − βG₁_b²) − ⅓K₆ρ⁵ za² + ⅓Fs_r(2ρ⁵G_a − ρ⁶βG_a²)

with `f'(z) = 1 − ρ za^{1/3} G₁_z`. Every function here is bounded and continuous on the
closed square `(β, ρ) ∈ [0, 1]²` (its arguments stay below 2, inside the table's range; the
tool checks that on every box), so N extends continuously
to the closed square, and it is enough that **N is positive definite on [0, 1]²**. (Where β or
ρβ is not small, K₅ and K₆ are also the plain difference quotients
`((yb^{−7/3})Fs_u − Fs_r)/β` and `(Fs_z za^{−2/3} − Fs_r)/(ρβ)`; the tool intersects both
enclosures.) Cross-checked against β²·diag(ρ,1)·M·diag(ρ,1) from M's plain form
(`g2_concavity ncheck`).

At the sides the limits show why the inequality holds there. As ρ → 0 (a → 0 with b
fixed), `N_ab → 0`, `N_bb → Fs(b)` and
`N_aa → Fs(0)[⅓ − (yb^{−1/3} − yb²)/(3β) + (2G_b − βG_b²)/3]` (K₅ tending to
`Fs(0)(yb^{−7/3} − 1)/β`), about 0.028 at β = 1. At the corner β = ρ = 0 this becomes
`N_aa = (4/81)G₀(3 − G₀) ≈ 0.032` with `G₀ = G(0) ≈ 0.234` (using `Fs(0) = (4/9)G₀`), and
`N_bb = Fs(0) ≈ 0.104`. In the unscaled radii this is the face normal, where Φ's radii and
the sections' sum agree to first order: their leading terms cancel exactly, and what is left
is positive.

## Lemma 5 (computer-assisted): N is positive definite on [0, 1]²

`g2_concavity scaled 0 1 0 1`:

- **The table.** The ingredients G, G₁, Fs and Fp, and `s = S(t)`, are enclosed over each
  cell of an even grid of 2¹⁸ cells of t in [0, ¾], by interval arithmetic (natural
  extension of the formulas of Lemma 1 on the cell). Each cell's range of s is taken from S at
  its two ends; S rising (`S₁ > 0` on the cell) makes this sound. The construction asserts
  `A₁, B₁, S₁ > 0` on every cell and the nodes' intervals of S strictly ordered (Lemma 1's
  conditions); it stops otherwise.
- **The lookup.** An ingredient over an interval of s is the hull of the cells whose s-range
  meets it (t* is monotone in s).
- **The boxes.** The square is cut into 64 × 64 boxes (each split in four, recursively, if it
  fails; none needs it). On each box N's entries are enclosed from the formulas of Lemma 4.
  The tool then certifies `N_aa > 0`, `N_bb > 0` and `det N > 0`, the latter after a diagonal
  congruence scaling, which keeps definiteness.
- **The arithmetic.** IEEE binary64 with every result's ends moved one ulp outwards (round to
  nearest is within half an ulp), and the cube root's ends four ulps (see below). An
  operation that cannot bound its result (a division by an interval holding 0, NaN) gives
  the whole line, and no box is certified through it. w is the interval
  [−0.5252124045776, −0.5252124045774].

Result: **4096 boxes certified, 0 not** (about 4 s, most of it the table). The margins
are wide: in samples along β = 0, 0.3, 0.6 and 1, N's smallest diagonal entry (0.028 to
0.032) is reached along ρ = 0, and `det N` stays above 3·10⁻³.

## What is proved, and what it rests on

- By hand: Lemmas 1 to 4 (the algebra and the reduction). The formulas are cross-checked
  numerically (`check`, `ncheck`, the kentron unit test), not proved by those checks.
- By machine: Lemma 1's sign conditions on [0, ¾], and Lemma 5.
- Trusted: binary64 arithmetic as IEEE specifies it, `f64::cbrt` within four ulps of the true
  root (the platform's libm; Rust does not guarantee its accuracy), and the compiler. Kentron's
  w lies in the interval used; the theorem covers every w in it.
- **Not proved: other thumbweights.** `fillet_setback_g2` uses only the equal-area one. A
  numerical survey (w from −2 to 0.9) found the same inequality, but that is evidence only.
  The property does depend on the section. With `p(a, b) = a + b − (a^q + b^q)^{1/q}` (the
  region of an L^q unit ball's quarter), H appears concave for q = 2 (the circle, borderline)
  but not for q = 3 (finite differences only, not certified). So no proof can rest on
  concavity and symmetry of p alone.
- **Decompositions that fail.** The per-turn split `Σp − Φ = ⅓ Σ_cyc [Σp − p(p(a, b), c)]`
  would make the theorem immediate if each bracket were concave. It is not:
  `g2_concavity counterexample` certifies that at (a, b, c) = (27, 2, 1) the bracket's Hessian
  on the slice c = 1 has determinant in [−3.84985·10⁻⁷, −3.84985·10⁻⁷] < 0 (so it is
  indefinite). Neither is `p(a,c) + p(b,c) − p(p(a,b),c)`, which differs from the bracket by
  the concave `p(a,b)`. A Cauchy–Schwarz decoupling of the rank-one terms also fails
  (numerically, by orders of magnitude near the octant's edges). The cancellation at the
  face normals is what makes the proof need the whole matrix.

## Reproducing

    rustc -O --edition 2021 tools/g2_concavity.rs -o /tmp/g2_concavity
    /tmp/g2_concavity scaled 0 1 0 1      # Lemma 5 (and Lemma 1's conditions)
    /tmp/g2_concavity check               # Lemma 2 against differences of H
    /tmp/g2_concavity ncheck              # Lemma 4 against Lemma 2
    /tmp/g2_concavity counterexample      # the per-turn split fails
    /tmp/g2_concavity certify 0.01 1 0.01 1 40 10 12   # optional: M itself on b, a/b ≥ 0.01
                                                       # (75,262 boxes; 3 min on 12 threads)
