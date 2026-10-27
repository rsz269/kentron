//! B-spline curves as Bézier segments and surfaces as grids of Bézier patches, for fast
//! evaluation: made once per curve or surface (by inserting knots until every break has the
//! degree's multiplicity) and shared by its clones, they spare each evaluation the knot
//! search and the basis recurrences.

use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, OnceLock};

use super::bspline::KnotVector;
use crate::math::{Point, Vector};

/// What is made from a curve or surface on first use; clones share it (until changed).
/// (Made from the start, so that every evaluation of a curve or surface takes the same
/// way and agrees with the others to the last bit.)
pub(super) struct Shared<T>(Arc<OnceLock<Option<T>>>);

/// The Bézier patches of a surface.
pub(super) type PatchCache = Shared<Patches>;
/// The Bézier segments of a curve.
pub(super) type SegmentCache = Shared<Segments>;

impl<T> Clone for Shared<T> {
    fn clone(&self) -> Self {
        Shared(self.0.clone())
    }
}

impl<T> Default for Shared<T> {
    fn default() -> Self {
        Shared(Arc::new(OnceLock::new()))
    }
}

impl<T> PartialEq for Shared<T> {
    /// (What the curve or surface is, not the cache, makes them equal.)
    fn eq(&self, _: &Self) -> bool {
        true
    }
}

impl<T> std::fmt::Debug for Shared<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str("Shared")
    }
}

impl<T> Shared<T> {
    /// What `make` makes the first time (`None`: too much to keep).
    pub(super) fn get(&self, make: impl FnOnce() -> Option<T>) -> Option<&T> {
        self.0.get_or_init(make).as_ref()
    }
}

/// Which surfaces share a patch cache (clones of one surface, unchanged): a key for other
/// caches of what is computed from a surface, holding on to it so that it stays unique.
#[derive(Clone)]
pub(crate) struct SurfaceId(Arc<OnceLock<Option<Patches>>>);

impl PartialEq for SurfaceId {
    fn eq(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.0, &other.0)
    }
}

impl PatchCache {
    pub(super) fn id(&self) -> SurfaceId {
        SurfaceId(self.0.clone())
    }
}

/// Bézier patches: their breaks along u and v and the control points `(w P, w)` of each,
/// patch by patch.
pub(super) struct Patches {
    u_breaks: Vec<f64>,
    v_breaks: Vec<f64>,
    u_degree: usize,
    v_degree: usize,
    /// Patch `(a, b)` holds `(u_degree + 1) × (v_degree + 1)` points from
    /// `(a · v_patches + b) · (u_degree + 1) · (v_degree + 1)`, row by row along u.
    points: Vec<[f64; 4]>,
    /// Weights other than 1 (else the sums skip them and the quotient rule).
    rational: bool,
    /// The patches along u and v found last, tried first.
    hints: Hint<2>,
}

/// The pieces found last by a lookup among breaks (see [`locate`]): where the next
/// evaluation most often is, Newton steps and marching moving little between them. On a
/// cache line of its own, so that its updates do not take from other threads the line of
/// the breaks they read.
#[repr(align(64))]
pub(super) struct Hint<const N: usize>([AtomicUsize; N]);

impl<const N: usize> Default for Hint<N> {
    fn default() -> Self {
        Hint(std::array::from_fn(|_| AtomicUsize::new(0)))
    }
}

/// More control points than this are not kept (the surface is evaluated as before): up to
/// 64 MB of them, as a rolling-ball tool's surfaces of a few hundred sections need.
const MOST_POINTS: usize = 1 << 21;

impl Patches {
    /// The patches of the surface with these knots and control points `(w P, w)`
    /// (`net[i][j]`, `i` along u); `None` if they would be too many.
    pub(super) fn new(u: &KnotVector, v: &KnotVector, net: &[Vec<[f64; 4]>]) -> Option<Patches> {
        let (nu, nv) = (net.len(), net[0].len());
        let (u_flat, v_flat) = (flat_knots(u), flat_knots(v));
        let (u_breaks, v_breaks) = (breaks(&u_flat, u.degree, nu), breaks(&v_flat, v.degree, nv));
        let (pu, pv) = (u.degree, v.degree);
        let (npu, npv) = (u_breaks.len() - 1, v_breaks.len() - 1);
        if npu * (pu + 1) * npv * (pv + 1) > MOST_POINTS {
            return None;
        }
        // Along u, each control point a row of the net.
        let rows: Vec<f64> = net.iter().flat_map(|row| row.iter().flatten().copied()).collect();
        let rows = split(&u_flat, u.degree, nu, &rows, 4 * nv);
        // Then along v, all rows at once (each control point along v a column of the rows:
        // the same sums, number by number, as row by row).
        let r = rows.len() / (4 * nv);
        let mut columns: Vec<f64> = Vec::with_capacity(rows.len());
        for j in 0..nv {
            for row in rows.chunks(4 * nv) {
                columns.extend_from_slice(&row[4 * j..4 * j + 4]);
            }
        }
        let grid = split(&v_flat, v.degree, nv, &columns, 4 * r);
        // Patch by patch: row a (pu + 1) + i, column b (pv + 1) + j at (column r + row) 4.
        let mut points = Vec::with_capacity(grid.len() / 4);
        for a in 0..npu {
            for b in 0..npv {
                for i in 0..=pu {
                    for j in 0..=pv {
                        let at = ((b * (pv + 1) + j) * r + a * (pu + 1) + i) * 4;
                        points.push([grid[at], grid[at + 1], grid[at + 2], grid[at + 3]]);
                    }
                }
            }
        }
        let rational = points.iter().any(|q| q[3] != 1.0);
        Some(Patches { u_breaks, v_breaks, u_degree: pu, v_degree: pv, points, rational, hints: Hint::default() })
    }

    /// The parameter bounds `(u0, u1, v0, v1)`.
    pub(super) fn bounds(&self) -> (f64, f64, f64, f64) {
        let (u, v) = (&self.u_breaks, &self.v_breaks);
        (u[0], u[u.len() - 1], v[0], v[v.len() - 1])
    }

    /// The point and its derivatives up to `order` (at most 2) at `(u, v)`; `None` above
    /// degree 25, and beyond the bounds by more than an end patch's length (where a short end
    /// patch's polynomial, far from its own range, magnifies rounding: the B-spline's knot
    /// window is evaluated instead).
    pub(super) fn evaluate(&self, u: f64, v: f64, order: usize) -> Option<Derivatives> {
        if !(near(&self.u_breaks, u) && near(&self.v_breaks, v)) {
            return None;
        }
        match order {
            0 => self.evaluate_to::<0>(u, v),
            1 => self.evaluate_to::<1>(u, v),
            _ => self.evaluate_to::<2>(u, v),
        }
    }

    /// The point and its derivatives `S(k, l)` of orders `k + l <= 3` at `(u, v)` as
    /// coordinates (`[k][l]`, the others zero); `None` as for [`Self::evaluate`].
    pub(super) fn evaluate3(&self, u: f64, v: f64) -> Option<[[[f64; 3]; 4]; 4]> {
        let (u0, u1, v0, v1) = self.bounds();
        if !(u0..=u1).contains(&u) || !(v0..=v1).contains(&v) {
            return None;
        }
        match self.u_degree.max(self.v_degree) {
            0..=3 => Some(self.sums3::<4>(u, v)),
            4..=8 => Some(self.sums3::<9>(u, v)),
            9..=25 => Some(self.sums3::<26>(u, v)),
            _ => None,
        }
    }

    /// [`Self::evaluate3`] with arrays of `N` (more than the degrees).
    fn sums3<const N: usize>(&self, u: f64, v: f64) -> [[[f64; 3]; 4]; 4] {
        const ORDER: usize = 3;
        let (pu, pv) = (self.u_degree, self.v_degree);
        let (a, s, ku) = locate(&self.u_breaks, u, &self.hints.0[0]);
        let (b, t, kv) = locate(&self.v_breaks, v, &self.hints.0[1]);
        let bu = bernstein3::<N>(pu, s, ku);
        let bv = bernstein3::<N>(pv, t, kv);
        let npv = self.v_breaks.len() - 1;
        let patch = &self.points[(a * npv + b) * (pu + 1) * (pv + 1)..][..(pu + 1) * (pv + 1)];
        // Homogeneous sums h[k][l] = Σ B_i^(k) B_j^(l) (w P, w) for k + l <= 3, along v first.
        let mut h = [[[0.0f64; 4]; 4]; 4];
        for i in 0..=pu {
            let row = &patch[i * (pv + 1)..][..pv + 1];
            let mut r = [[0.0f64; 4]; 4];
            for j in 0..=pv {
                let q = row[j];
                for l in 0..=ORDER {
                    let f = bv[l][j];
                    for c in 0..4 {
                        r[l][c] += f * q[c];
                    }
                }
            }
            for k in 0..=ORDER {
                let f = bu[k][i];
                for l in 0..=ORDER - k {
                    for c in 0..4 {
                        h[k][l][c] += f * r[l][c];
                    }
                }
            }
        }
        let mut out = [[[0.0f64; 3]; 4]; 4];
        if !self.rational {
            for k in 0..=ORDER {
                for l in 0..=ORDER - k {
                    out[k][l] = [h[k][l][0], h[k][l][1], h[k][l][2]];
                }
            }
            return out;
        }
        // The quotient rule (The NURBS Book, A4.4).
        let binom = |n: usize, r: usize| (1..=r).fold(1.0, |acc, i| acc * (n - r + i) as f64 / i as f64);
        for k in 0..=ORDER {
            for l in 0..=ORDER - k {
                let mut c = [h[k][l][0], h[k][l][1], h[k][l][2]];
                for j in 1..=l {
                    let f = binom(l, j) * h[0][j][3];
                    for d in 0..3 {
                        c[d] -= f * out[k][l - j][d];
                    }
                }
                for i in 1..=k {
                    let f = binom(k, i) * h[i][0][3];
                    for d in 0..3 {
                        c[d] -= f * out[k - i][l][d];
                    }
                    for j in 1..=l {
                        let g = binom(k, i) * binom(l, j) * h[i][j][3];
                        for d in 0..3 {
                            c[d] -= g * out[k - i][l - j][d];
                        }
                    }
                }
                for d in 0..3 {
                    c[d] /= h[0][0][3];
                }
                out[k][l] = c;
            }
        }
        out
    }

    /// [`Self::evaluate`] at `(u, v)` for each `u` of `us` (best ascending), into `out`, bit
    /// for bit the same: each patch row is summed along v once for all the `u`s in the
    /// patch, then each `u` sums the rows. `false` (and `out` untouched) where `v` is out of
    /// bounds or the degree above 25; a `u` out of bounds gets `None`.
    pub(super) fn evaluate_along_u(&self, us: &[f64], v: f64, order: usize, out: &mut Vec<Option<Derivatives>>) -> bool {
        let (_, _, v0, v1) = self.bounds();
        if !(v0..=v1).contains(&v) || self.u_degree.max(self.v_degree) > 25 {
            return false;
        }
        macro_rules! sized {
            ($nu:expr, $nv:expr, $exact:literal) => {
                match (order, self.rational) {
                    (0, false) => self.along_u::<{ $nu }, { $nv }, 0, $exact, 3>(us, v, out),
                    (0, true) => self.along_u::<{ $nu }, { $nv }, 0, $exact, 4>(us, v, out),
                    (1, false) => self.along_u::<{ $nu }, { $nv }, 1, $exact, 3>(us, v, out),
                    (1, true) => self.along_u::<{ $nu }, { $nv }, 1, $exact, 4>(us, v, out),
                    (_, false) => self.along_u::<{ $nu }, { $nv }, 2, $exact, 3>(us, v, out),
                    (_, true) => self.along_u::<{ $nu }, { $nv }, 2, $exact, 4>(us, v, out),
                }
            };
        }
        // (The usual low degrees with sizes known to the compiler, as for single points.)
        macro_rules! exact {
            ($(($pu:literal, $pv:literal)),*) => {
                match (self.u_degree, self.v_degree) {
                    $(($pu, $pv) => {
                        sized!($pu + 1, $pv + 1, true);
                        return true;
                    })*
                    _ => {}
                }
            };
        }
        exact!((1, 1), (1, 2), (1, 3), (2, 1), (2, 2), (2, 3), (3, 1), (3, 2), (3, 3), (1, 4), (1, 7));
        match self.u_degree.max(self.v_degree) {
            0..=3 => sized!(4, 4, false),
            4..=8 => sized!(9, 9, false),
            _ => sized!(26, 26, false),
        }
        true
    }

    /// With arrays of `NU` and `NV` (the degrees plus one if `EXACT`, else more than them).
    fn along_u<const NU: usize, const NV: usize, const ORDER: usize, const EXACT: bool, const C: usize>(&self, us: &[f64], v: f64, out: &mut Vec<Option<Derivatives>>) {
        let (pu, pv) = if EXACT { (NU - 1, NV - 1) } else { (self.u_degree, self.v_degree) };
        let (u0, u1, _, _) = self.bounds();
        let (b, t, kv) = locate(&self.v_breaks, v, &self.hints.0[1]);
        let bv = bernstein::<NV, ORDER>(pv, t, kv);
        let npv = self.v_breaks.len() - 1;
        // The rows of the patch in use summed along v, per derivative order in v (as `sums`
        // forms them, over the same ranges).
        let mut rows = [[[0.0f64; 4]; 3]; NU];
        let mut summed = usize::MAX;
        out.reserve(us.len());
        for &u in us {
            if !(u0..=u1).contains(&u) {
                out.push(None);
                continue;
            }
            let (a, s, ku) = locate(&self.u_breaks, u, &self.hints.0[0]);
            if a != summed {
                summed = a;
                let patch = &self.points[(a * npv + b) * (pu + 1) * (pv + 1)..][..(pu + 1) * (pv + 1)];
                for i in 0..NU {
                    if i > pu {
                        break;
                    }
                    let row = &patch[i * (pv + 1)..][..pv + 1];
                    let mut r = [[0.0f64; 4]; 3];
                    for j in 0..NV {
                        if j > pv {
                            break;
                        }
                        let q = row[j];
                        for l in 0..ORDER + 1 {
                            let f = bv[l][j];
                            for c in 0..C {
                                r[l][c] += f * q[c];
                            }
                        }
                    }
                    rows[i] = r;
                }
            }
            let bu = bernstein::<NU, ORDER>(pu, s, ku);
            let mut h = [[[0.0f64; 4]; 3]; 3];
            for i in 0..NU {
                if i > pu {
                    break;
                }
                let r = &rows[i];
                for k in 0..ORDER + 1 {
                    let f = bu[k][i];
                    for l in 0..ORDER + 1 {
                        if k + l <= ORDER {
                            for c in 0..C {
                                h[k][l][c] += f * r[l][c];
                            }
                        }
                    }
                }
            }
            out.push(Some(if C == 3 {
                let vector = |c: [f64; 4]| Vector::new(c[0], c[1], c[2]);
                Derivatives {
                    point: Point::new(h[0][0][0], h[0][0][1], h[0][0][2]),
                    du: vector(h[1][0]),
                    dv: vector(h[0][1]),
                    duu: vector(h[2][0]),
                    dvv: vector(h[0][2]),
                    duv: vector(h[1][1]),
                }
            } else {
                project(&h, ORDER)
            }));
        }
    }

    fn evaluate_to<const ORDER: usize>(&self, u: f64, v: f64) -> Option<Derivatives> {
        if self.rational { self.evaluate_with::<ORDER, 4>(u, v) } else { self.evaluate_with::<ORDER, 3>(u, v) }
    }

    /// With `C` components of the control points: 4 (`w P` and `w`), or 3 when the weights
    /// are all 1.
    fn evaluate_with<const ORDER: usize, const C: usize>(&self, u: f64, v: f64) -> Option<Derivatives> {
        // The usual low degrees with sizes known to the compiler, others in larger arrays.
        macro_rules! exact {
            ($(($pu:literal, $pv:literal)),*) => {
                match (self.u_degree, self.v_degree) {
                    $(($pu, $pv) => return Some(self.sums::<{ $pu + 1 }, { $pv + 1 }, ORDER, true, C>(u, v)),)*
                    _ => {}
                }
            };
        }
        exact!((1, 1), (1, 2), (1, 3), (2, 1), (2, 2), (2, 3), (3, 1), (3, 2), (3, 3), (1, 4), (1, 7), (1, 8), (3, 4), (4, 3), (4, 4), (3, 8), (6, 6), (8, 6), (6, 8));
        match self.u_degree.max(self.v_degree) {
            0..=8 => Some(self.sums::<9, 9, ORDER, false, C>(u, v)),
            9..=25 => Some(self.sums::<26, 26, ORDER, false, C>(u, v)),
            _ => None,
        }
    }

    /// The evaluation with arrays of `NU` and `NV` (the degrees plus one if `EXACT`, else
    /// more than them).
    #[inline(always)]
    fn sums<const NU: usize, const NV: usize, const ORDER: usize, const EXACT: bool, const C: usize>(&self, u: f64, v: f64) -> Derivatives {
        let (pu, pv) = if EXACT { (NU - 1, NV - 1) } else { (self.u_degree, self.v_degree) };
        let (a, s, ku) = locate(&self.u_breaks, u, &self.hints.0[0]);
        let (b, t, kv) = locate(&self.v_breaks, v, &self.hints.0[1]);
        let bu = bernstein::<NU, ORDER>(pu, s, ku);
        let bv = bernstein::<NV, ORDER>(pv, t, kv);
        let npv = self.v_breaks.len() - 1;
        let patch = &self.points[(a * npv + b) * (pu + 1) * (pv + 1)..][..(pu + 1) * (pv + 1)];
        // Homogeneous sums h[k][l] = Σ B_i^(k) B_j^(l) (w P, w) for k + l <= ORDER.
        // (Every loop over a range known to the compiler, cut short at the degree: with the
        // exact sizes they unroll, the sums kept in registers, in the same order.)
        let mut h = [[[0.0f64; 4]; 3]; 3];
        for i in 0..NU {
            if i > pu {
                break;
            }
            let row = &patch[i * (pv + 1)..][..pv + 1];
            // Along v first: row sums for each derivative order in v.
            let mut r = [[0.0f64; 4]; 3];
            for j in 0..NV {
                if j > pv {
                    break;
                }
                let q = row[j];
                for l in 0..ORDER + 1 {
                    let f = bv[l][j];
                    for c in 0..C {
                        r[l][c] += f * q[c];
                    }
                }
            }
            for k in 0..ORDER + 1 {
                let f = bu[k][i];
                for l in 0..ORDER + 1 {
                    if k + l <= ORDER {
                        for c in 0..C {
                            h[k][l][c] += f * r[l][c];
                        }
                    }
                }
            }
        }
        if C == 3 {
            // Weights of 1: the sums are the point and derivatives.
            let vector = |c: [f64; 4]| Vector::new(c[0], c[1], c[2]);
            return Derivatives {
                point: Point::new(h[0][0][0], h[0][0][1], h[0][0][2]),
                du: vector(h[1][0]),
                dv: vector(h[0][1]),
                duu: vector(h[2][0]),
                dvv: vector(h[0][2]),
                duv: vector(h[1][1]),
            };
        }
        project(&h, ORDER)
    }
}

/// Bézier segments of a curve: their breaks and the control points `(w P, w)` of each
/// (`degree + 1` per segment); 2D points have `z = 0`.
pub(super) struct Segments {
    breaks: Vec<f64>,
    degree: usize,
    points: Vec<[f64; 4]>,
    /// Weights other than 1.
    rational: bool,
    /// For a periodic curve, its period's range (parameters outside are wrapped into it).
    period: Option<(f64, f64)>,
    /// The segment found last, tried first.
    hint: Hint<1>,
}

impl Segments {
    /// The segments of the curve with these flat knots and control points `(w P, w)` (as many
    /// as the flat knots need, a periodic curve's repeated); `None` if too many.
    pub(super) fn new(flat: &[f64], degree: usize, points: &[[f64; 4]], periodic: bool) -> Option<Segments> {
        let count = points.len();
        let breaks = breaks(flat, degree, count);
        if (breaks.len() - 1) * (degree + 1) > MOST_POINTS {
            return None;
        }
        let flat_points: Vec<f64> = points.iter().flatten().copied().collect();
        let split = split(flat, degree, count, &flat_points, 4);
        let points: Vec<[f64; 4]> = split.chunks(4).map(|c| [c[0], c[1], c[2], c[3]]).collect();
        let rational = points.iter().any(|q| q[3] != 1.0);
        let period = periodic.then(|| (breaks[0], breaks[breaks.len() - 1]));
        Some(Segments { breaks, degree, points, rational, period, hint: Hint::default() })
    }

    /// The point and its derivatives up to `order` (at most 2) at `t` as coordinates;
    /// `None` above degree 25, and outside the range of a curve that is not periodic (as for
    /// patches).
    pub(super) fn evaluate(&self, t: f64, order: usize) -> Option<[[f64; 3]; 3]> {
        let (lo, hi) = (self.breaks[0], self.breaks[self.breaks.len() - 1]);
        let t = match self.period {
            Some(_) if !(lo..=hi).contains(&t) => lo + (t - lo).rem_euclid(hi - lo),
            None if !(lo..=hi).contains(&t) => return None,
            _ => t,
        };
        match order {
            0 => self.evaluate_to::<0>(t),
            1 => self.evaluate_to::<1>(t),
            _ => self.evaluate_to::<2>(t),
        }
    }

    fn evaluate_to<const ORDER: usize>(&self, t: f64) -> Option<[[f64; 3]; 3]> {
        if self.rational { self.evaluate_with::<ORDER, 4>(t) } else { self.evaluate_with::<ORDER, 3>(t) }
    }

    fn evaluate_with<const ORDER: usize, const C: usize>(&self, t: f64) -> Option<[[f64; 3]; 3]> {
        Some(match self.degree {
            1 => self.sums::<2, ORDER, true, C>(t),
            2 => self.sums::<3, ORDER, true, C>(t),
            3 => self.sums::<4, ORDER, true, C>(t),
            0..=8 => self.sums::<9, ORDER, false, C>(t),
            9..=25 => self.sums::<26, ORDER, false, C>(t),
            _ => return None,
        })
    }

    #[inline(always)]
    fn sums<const N: usize, const ORDER: usize, const EXACT: bool, const C: usize>(&self, t: f64) -> [[f64; 3]; 3] {
        let p = if EXACT { N - 1 } else { self.degree };
        let (k, s, scale) = locate(&self.breaks, t, &self.hint.0[0]);
        let b = bernstein::<N, ORDER>(p, s, scale);
        let segment = &self.points[k * (p + 1)..][..p + 1];
        // (Loops over ranges known to the compiler, as in the patches' sums.)
        let mut h = [[0.0f64; 4]; 3];
        for i in 0..N {
            if i > p {
                break;
            }
            let q = segment[i];
            for l in 0..ORDER + 1 {
                let f = b[l][i];
                for c in 0..C {
                    h[l][c] += f * q[c];
                }
            }
        }
        if C == 3 {
            // Weights of 1: the sums are the point and derivatives.
            return [[h[0][0], h[0][1], h[0][2]], [h[1][0], h[1][1], h[1][2]], [h[2][0], h[2][1], h[2][2]]];
        }
        // From homogeneous sums to the point and derivatives (the quotient rule).
        let w = h[0][3];
        let mut out = [[0.0f64; 3]; 3];
        for c in 0..3 {
            out[0][c] = h[0][c] / w;
        }
        if ORDER >= 1 {
            for c in 0..3 {
                out[1][c] = (h[1][c] - h[1][3] * out[0][c]) / w;
            }
        }
        if ORDER >= 2 {
            for c in 0..3 {
                out[2][c] = (h[2][c] - 2.0 * h[1][3] * out[1][c] - h[2][3] * out[0][c]) / w;
            }
        }
        out
    }
}

/// A point of a surface and its derivatives (those not asked for are zero).
pub(super) struct Derivatives {
    pub point: Point,
    pub du: Vector,
    pub dv: Vector,
    pub duu: Vector,
    pub dvv: Vector,
    pub duv: Vector,
}

/// From homogeneous sums to the point and derivatives (the quotient rule).
fn project(a: &[[[f64; 4]; 3]; 3], order: usize) -> Derivatives {
    let w = a[0][0][3];
    let s = [a[0][0][0] / w, a[0][0][1] / w, a[0][0][2] / w];
    let vector = |c: [f64; 3]| Vector::new(c[0], c[1], c[2]);
    let zero = Vector::new(0.0, 0.0, 0.0);
    let point = Point::new(s[0], s[1], s[2]);
    if order == 0 {
        return Derivatives { point, du: zero, dv: zero, duu: zero, dvv: zero, duv: zero };
    }
    let first = |k: usize, l: usize| std::array::from_fn::<f64, 3, _>(|c| (a[k][l][c] - a[k][l][3] * s[c]) / w);
    let (s_u, s_v) = (first(1, 0), first(0, 1));
    if order == 1 {
        return Derivatives { point, du: vector(s_u), dv: vector(s_v), duu: zero, dvv: zero, duv: zero };
    }
    let s_uu: [f64; 3] = std::array::from_fn(|c| (a[2][0][c] - 2.0 * a[1][0][3] * s_u[c] - a[2][0][3] * s[c]) / w);
    let s_vv: [f64; 3] = std::array::from_fn(|c| (a[0][2][c] - 2.0 * a[0][1][3] * s_v[c] - a[0][2][3] * s[c]) / w);
    let s_uv: [f64; 3] = std::array::from_fn(|c| (a[1][1][c] - a[1][0][3] * s_v[c] - a[0][1][3] * s_u[c] - a[1][1][3] * s[c]) / w);
    Derivatives { point, du: vector(s_u), dv: vector(s_v), duu: vector(s_uu), dvv: vector(s_vv), duv: vector(s_uv) }
}

/// Whether `t` is within the breaks, or beyond them by at most the end patch's length (where
/// its polynomial, the B-spline's own there, is still near its range). Not a NaN.
#[inline(always)]
fn near(breaks: &[f64], t: f64) -> bool {
    let (lo, hi) = (breaks[0], breaks[breaks.len() - 1]);
    t >= lo - (breaks[1] - lo) && t <= hi + (hi - breaks[breaks.len() - 2])
}

/// The patch of `t` among the breaks (the last whose start is at most `t`, an end one
/// outside), `t`'s place in it (0 to 1 inside) and the derivative of that place in `t`.
/// The patch `hint` holds is tried first, and the one found kept there.
#[inline(always)]
fn locate(breaks: &[f64], t: f64, hint: &AtomicUsize) -> (usize, f64, f64) {
    let last = breaks.len() - 2;
    // (The hinted patch where it is the one the search would find: breaks strictly increase,
    // so that is the patch with breaks[k] <= t < breaks[k + 1], the end ones reaching out. A
    // NaN fails both tests and is searched for, as before. Relaxed: a patch another thread
    // found is as good a guess.)
    let h = hint.load(Ordering::Relaxed);
    let k = if h <= last && (h == 0 || breaks[h] <= t) && (h == last || t < breaks[h + 1]) {
        h
    } else {
        let k = breaks[1..=last].partition_point(|&b| b <= t);
        hint.store(k, Ordering::Relaxed);
        k
    };
    let (a, b) = (breaks[k], breaks[k + 1]);
    let scale = 1.0 / (b - a);
    (k, (t - a) * scale, scale)
}

/// The Bernstein polynomials of degree `p` (below `N`) at `s` and, up to `ORDER`, their
/// derivatives in the parameter whose patch place is `s` (`scale` = ds/dt).
#[inline(always)]
fn bernstein<const N: usize, const ORDER: usize>(p: usize, s: f64, scale: f64) -> [[f64; N]; 3] {
    let r = 1.0 - s;
    let mut out = [[0.0f64; N]; 3];
    // Degree k from degree k - 1, from the top down; the lower degrees the derivatives need
    // kept in out[1] (p - 1) and out[2] (p - 2) on the way. (Over ranges known to the
    // compiler, cut short at the degree, so that the small sizes unroll and keep the arrays
    // in registers: the same operations in the same order. Not the largest, for degrees up
    // to 25, where the masked steps would outnumber the others.)
    let unrolled = N <= 9;
    let mut b = [0.0f64; N];
    b[0] = 1.0;
    for k in 1..N {
        if k > p {
            break;
        }
        if ORDER >= 2 && k + 1 == p {
            out[2] = b;
        }
        if ORDER >= 1 && k == p {
            out[1] = b;
        }
        b[k] = s * b[k - 1];
        if unrolled {
            for i in (1..N).rev() {
                if i < k {
                    b[i] = r * b[i] + s * b[i - 1];
                }
            }
        } else {
            for i in (1..k).rev() {
                b[i] = r * b[i] + s * b[i - 1];
            }
        }
        b[0] *= r;
    }
    if ORDER >= 2 {
        // d²/ds² B_i^p = p (p - 1) (B_{i-2}^{p-2} - 2 B_{i-1}^{p-2} + B_i^{p-2}).
        let low = out[2];
        let f = if p >= 2 { (p * (p - 1)) as f64 * scale * scale } else { 0.0 };
        for i in 0..N {
            if i > p {
                break;
            }
            let at = |j: usize| if j + 2 <= p { low[j] } else { 0.0 };
            let (a, b) = (if i >= 2 { at(i - 2) } else { 0.0 }, if i >= 1 { at(i - 1) } else { 0.0 });
            out[2][i] = f * (a - 2.0 * b + at(i));
        }
    }
    if ORDER >= 1 {
        // d/ds B_i^p = p (B_{i-1}^{p-1} - B_i^{p-1}).
        let low = out[1];
        let f = p as f64 * scale;
        for i in 0..N {
            if i > p {
                break;
            }
            let left = if i >= 1 { low[i - 1] } else { 0.0 };
            let right = if i < p { low[i] } else { 0.0 };
            out[1][i] = f * (left - right);
        }
    }
    out[0] = b;
    out
}

/// The Bernstein polynomials of degree `p` (below `N`) at `s` and their derivatives up to the
/// third in the parameter whose patch place is `s` (`scale` = ds/dt).
fn bernstein3<const N: usize>(p: usize, s: f64, scale: f64) -> [[f64; N]; 4] {
    let r = 1.0 - s;
    let mut out = [[0.0f64; N]; 4];
    // Degree k from degree k - 1; the lower degrees p - m the derivatives need kept in
    // out[m] on the way.
    let mut b = [0.0f64; N];
    b[0] = 1.0;
    for k in 1..=p {
        for m in 1..=3 {
            if k + m == p + 1 {
                out[m] = b;
            }
        }
        b[k] = s * b[k - 1];
        for i in (1..k).rev() {
            b[i] = r * b[i] + s * b[i - 1];
        }
        b[0] *= r;
    }
    // d^m/ds^m B_i^p = p (p - 1) ... (p - m + 1) Σ_n (-1)^n C(m, n) B_{i-m+n}^{p-m}.
    for m in 1..=3 {
        if m > p {
            out[m] = [0.0; N];
            continue;
        }
        let low = out[m];
        let f = (0..m).map(|t| (p - t) as f64 * scale).product::<f64>();
        let binom = [[1.0, 0.0, 0.0, 0.0], [1.0, 1.0, 0.0, 0.0], [1.0, 2.0, 1.0, 0.0], [1.0, 3.0, 3.0, 1.0]];
        for i in 0..=p {
            let mut sum = 0.0;
            for n in 0..=m {
                if i + n >= m && i + n - m <= p - m {
                    let sign = if n % 2 == 0 { 1.0 } else { -1.0 };
                    sum += sign * binom[m][n] * low[i + n - m];
                }
            }
            out[m][i] = f * sum;
        }
    }
    out[0] = b;
    out
}

/// The distinct knots over the parameter range of `count` poles (from flat knot `degree` to
/// flat knot `count`).
fn breaks(flat: &[f64], degree: usize, count: usize) -> Vec<f64> {
    let (lo, hi) = (flat[degree], flat[count]);
    let mut out: Vec<f64> = flat[degree..=count].to_vec();
    out.dedup();
    debug_assert!(out[0] == lo && out[out.len() - 1] == hi);
    out
}

fn flat_knots(kv: &KnotVector) -> Vec<f64> {
    kv.knots.iter().zip(&kv.multiplicities).flat_map(|(&k, &m)| std::iter::repeat_n(k, m)).collect()
}

/// The control points (each `stride` numbers) of the Bézier segments of the B-spline with
/// these flat knots and `count` control points, segment by segment (`degree + 1` each): knots
/// inserted until each break has the degree's multiplicity (Boehm's algorithm).
fn split(flat: &[f64], degree: usize, count: usize, points: &[f64], stride: usize) -> Vec<f64> {
    let p = degree;
    let breaks = breaks(flat, degree, count);
    let hi = breaks[breaks.len() - 1];
    // The knots and points as gap buffers: those before the place of insertion in order,
    // those after it reversed (each point's numbers kept in order), the place moving on
    // with the breaks, so that each insertion costs the degree, not the length.
    let mut knots: Vec<f64> = Vec::with_capacity(flat.len() + breaks.len() * p);
    let mut knots_after: Vec<f64> = flat[..count + p + 1].iter().rev().copied().collect();
    let mut head: Vec<f64> = Vec::with_capacity(points.len() + breaks.len() * p * stride);
    let mut after: Vec<f64> = Vec::with_capacity(points.len() + stride);
    for point in points.chunks(stride).rev() {
        after.extend_from_slice(point);
    }
    let mut fresh = vec![0.0; p * stride];
    for &t in &breaks {
        // Every knot up to t before the place.
        while knots_after.last().is_some_and(|&k| k <= t) {
            knots.push(knots_after.pop().unwrap());
        }
        let have = knots.iter().rev().take_while(|&&k| k == t).count();
        for _ in have..p {
            // The span with flat[k] <= t <= flat[k + 1] and flat[k] < flat[k + 1].
            let k = if t < hi { knots.len() - 1 } else { knots.len() - knots.iter().rev().take_while(|&&k| k == t).count() - 1 };
            let knot = |i: usize| if i < knots.len() { knots[i] } else { knots_after[knots_after.len() - 1 - (i - knots.len())] };
            // Points 0..=k before the place.
            while head.len() < (k + 1) * stride {
                let from = after.len() - stride;
                head.extend_from_slice(&after[from..]);
                after.truncate(from);
            }
            for (n, i) in (k + 1 - p..=k).enumerate() {
                let alpha = (t - knot(i)) / (knot(i + p) - knot(i));
                for c in 0..stride {
                    fresh[n * stride + c] = alpha * head[i * stride + c] + (1.0 - alpha) * head[(i - 1) * stride + c];
                }
            }
            // Points k + 1 - p..k replaced by the p new ones, point k kept after them.
            after.extend_from_slice(&head[k * stride..(k + 1) * stride]);
            head.truncate((k + 1 - p) * stride);
            head.extend_from_slice(&fresh);
            knots.push(t);
        }
    }
    let flat: Vec<f64> = knots.into_iter().chain(knots_after.into_iter().rev()).collect();
    for point in after.rchunks(stride) {
        head.extend_from_slice(point);
    }
    let points = head;
    let n = points.len() / stride;
    // Segment by segment: the points of the span starting at each break.
    let mut out = Vec::with_capacity((breaks.len() - 1) * (p + 1) * stride);
    for w in breaks.windows(2) {
        let span = flat[..n].partition_point(|&x| x <= w[0]) - 1;
        out.extend_from_slice(&points[(span - p) * stride..(span + 1) * stride]);
    }
    out
}
