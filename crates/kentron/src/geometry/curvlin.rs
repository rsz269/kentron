//! OCCT's `Approx_CurvlinFunc`: a curve reparametrized by its normalized arc length
//! `s` in `[0, 1]`, ported for a 3D curve (OCCT's `myCase == 1`).
//!
//! The curve-on-surface cases (`myCase` 2 and 3: one or two 2D curves on surfaces, through
//! `Adaptor3d_CurveOnSurface`, with `EvalCase2`, `EvalCase3` and `EvalCurOnSur`) are not
//! ported; kentron has no `Adaptor3d_CurveOnSurface` port for them to run on.
//!
//! `Init` tabulates the arc length at ten equal parameter steps in each C3 span
//! (`myUi_1`, `mySi_1`, the lengths normalized by the whole). `GetUParameter` finds the
//! table entry below `s`, interpolates a guess for the parameter by a cubic through four
//! table points (one of them replaced by the previous answer when that falls between them)
//! and solves for the exact parameter with `GCPnts_AbscissaPoint` from that entry.

use std::cell::Cell;

use super::Continuity;
use super::abscissa::{AbscissaPoint, length_between};
use super::adaptor::CurveAdaptor;
use crate::{Error, Result};

/// `cubic(X, Xi, Yi)`: Newton's divided-difference cubic through four points, at `x`.
fn cubic(x: f64, xi: &[f64; 4], yi: &[f64; 4]) -> f64 {
    let i1 = (yi[0] - yi[1]) / (xi[0] - xi[1]);
    let i2 = (yi[1] - yi[2]) / (xi[1] - xi[2]);
    let i3 = (yi[2] - yi[3]) / (xi[2] - xi[3]);

    let i21 = (i1 - i2) / (xi[0] - xi[2]);
    let i22 = (i2 - i3) / (xi[1] - xi[3]);

    let i31 = (i21 - i22) / (xi[0] - xi[3]);

    yi[0] + (x - xi[0]) * (i1 + (x - xi[1]) * (i21 + (x - xi[2]) * i31))
}

/// `findfourpoints`: four table points around interval `n_interval`, with `(prev_s, prev_u)`
/// put in place of one of them when it falls strictly between two. An error where OCCT
/// raises `Standard_ConstructionError` (fewer than three table intervals).
fn find_four_points(mut n_interval: usize, si: &[f64], ui: &[f64], prev_s: f64, prev_u: f64) -> Result<([f64; 4], [f64; 4])> {
    let nb_int = si.len() - 1;
    if nb_int < 3 {
        return Err(Error::InvalidInput("Approx_CurvlinFunc::GetUParameter: fewer than three intervals".into()));
    }
    if n_interval < 1 {
        n_interval = 1;
    } else if n_interval > nb_int - 2 {
        n_interval = nb_int - 2;
    }
    let mut xi = [0.0; 4];
    let mut yi = [0.0; 4];
    for i in 0..4 {
        xi[i] = si[n_interval - 1 + i];
        yi[i] = ui[n_interval - 1 + i];
    }
    // Try to insert (S, U).
    for i in 0..3 {
        if xi[i] < prev_s && prev_s < xi[i + 1] {
            for j in 0..i {
                xi[j] = xi[j + 1];
                yi[j] = yi[j + 1];
            }
            xi[i] = prev_s;
            yi[i] = prev_u;
            break;
        }
    }
    Ok((xi, yi))
}

/// The curve held: the one given, or a trimmed copy after [`CurvlinFunc::trim`].
enum Held<'a> {
    Borrowed(&'a dyn CurveAdaptor),
    Owned(Box<dyn CurveAdaptor>),
}

impl Held<'_> {
    fn get(&self) -> &dyn CurveAdaptor {
        match self {
            Held::Borrowed(c) => *c,
            Held::Owned(c) => c.as_ref(),
        }
    }
}

/// `Approx_CurvlinFunc` for a 3D curve.
#[doc(alias = "Approx_CurvlinFunc")]
pub(crate) struct CurvlinFunc<'a> {
    /// `myC3D`.
    curve: Held<'a>,
    first_s: f64,
    last_s: f64,
    /// `myFirstU1`, `myLastU1`: the parameter range at construction.
    first_u1: f64,
    last_u1: f64,
    /// `myLength`: the length of the curve at construction (kept by `Trim`).
    length: f64,
    /// `myTolLen`.
    tol_len: f64,
    /// `myPrevS`, `myPrevU`: the last answer of `GetUParameter`, which seeds the next.
    prev_s: Cell<f64>,
    prev_u: Cell<f64>,
    /// `myUi_1`, `mySi_1`.
    ui: Vec<f64>,
    si: Vec<f64>,
}

impl<'a> CurvlinFunc<'a> {
    /// `Approx_CurvlinFunc(C, Tol)`: lengths are computed to `tolerance`.
    pub(crate) fn new(curve: &'a dyn CurveAdaptor, tolerance: f64) -> Self {
        let mut function = CurvlinFunc {
            curve: Held::Borrowed(curve),
            first_s: 0.0,
            last_s: 1.0,
            first_u1: 0.0,
            last_u1: 0.0,
            length: 0.0,
            tol_len: tolerance,
            prev_s: Cell::new(0.0),
            prev_u: Cell::new(0.0),
            ui: Vec::new(),
            si: Vec::new(),
        };
        function.init();
        function
    }

    /// `Init()`.
    fn init(&mut self) {
        let curve = self.curve.get();
        let (si, ui) = self.init_tables(curve);
        self.first_u1 = curve.first_parameter();
        self.last_u1 = curve.last_parameter();
        self.si = si;
        self.ui = ui;
        self.compute_length();
    }

    /// `Init(C, Si, Ui)`: the table of normalized arc lengths at ten equal steps of each C3
    /// span ("cut at the curve's discontinuities"). The parameters are accumulated step by
    /// step from the first, as OCCT does.
    fn init_tables(&self, curve: &dyn CurveAdaptor) -> (Vec<f64>, Vec<f64>) {
        let first_u = curve.first_parameter();
        let last_u = curve.last_parameter();
        const NB_INT: usize = 10;
        let intervals = curve.intervals(Continuity::C3);
        let nb_int_c3 = intervals.len() - 1;
        let disc = if nb_int_c3 > 1 { intervals } else { vec![first_u, last_u] };

        let mut ui = vec![0.0; nb_int_c3 * NB_INT + 1];
        let mut si = vec![0.0; nb_int_c3 * NB_INT + 1];
        ui[0] = first_u;
        si[0] = 0.0;
        let mut i = 1;
        for j in 1..=nb_int_c3 {
            let step = (disc[j] - disc[j - 1]) / NB_INT as f64;
            for _ in 1..=NB_INT {
                ui[i] = ui[i - 1] + step;
                si[i] = si[i - 1] + self.length_of(curve, ui[i - 1], ui[i]);
                i += 1;
            }
        }
        let len = si[si.len() - 1];
        for s in &mut si {
            *s /= len;
        }
        self.prev_s.set(self.first_s);
        self.prev_u.set(first_u);
        (si, ui)
    }

    /// `SetTol(Tol)`.
    pub(crate) fn set_tolerance(&mut self, tolerance: f64) {
        self.tol_len = tolerance;
    }

    /// `FirstParameter()`.
    pub(crate) fn first_parameter(&self) -> f64 {
        self.first_s
    }

    /// `LastParameter()`.
    pub(crate) fn last_parameter(&self) -> f64 {
        self.last_s
    }

    /// The curve as now trimmed.
    pub(crate) fn curve(&self) -> &dyn CurveAdaptor {
        self.curve.get()
    }

    /// `NbIntervals(S)`: those of the curve.
    pub(crate) fn nb_intervals(&self, continuity: Continuity) -> usize {
        self.curve.get().intervals(continuity).len() - 1
    }

    /// `Intervals(T, S)`: the curve's intervals, as arc-length parameters.
    pub(crate) fn intervals(&self, continuity: Continuity) -> Vec<f64> {
        self.curve.get().intervals(continuity).into_iter().map(|u| self.get_s_parameter(u)).collect()
    }

    /// `Trim(First, Last, Tol)`: restricts the curve to the arc-length range `[first,
    /// last]` (the curve brought back to its original range first). Nothing changes for a
    /// range shorter than `tolerance`. An error where OCCT raises `Standard_OutOfRange` (a
    /// range outside `[0, 1]`) or where `GetUParameter` fails.
    #[doc(alias = "Approx_CurvlinFunc::Trim")]
    pub(crate) fn trim(&mut self, first: f64, last: f64, tolerance: f64) -> Result<()> {
        if first < 0.0 || last > 1.0 {
            return Err(Error::InvalidInput("Approx_CurvlinFunc::Trim: range outside [0, 1]".into()));
        }
        if last - first < tolerance {
            return Ok(());
        }
        self.curve = Held::Owned(self.curve.get().trimmed(self.first_u1, self.last_u1, tolerance));
        let first_u = self.get_u_parameter_on(self.curve.get(), first)?;
        let last_u = self.get_u_parameter_on(self.curve.get(), last)?;
        self.curve = Held::Owned(self.curve.get().trimmed(first_u, last_u, tolerance));
        self.first_s = first;
        self.last_s = last;
        Ok(())
    }

    /// `Length()`: the length of the whole curve.
    fn compute_length(&mut self) {
        let curve = self.curve.get();
        let (first_u, last_u) = (curve.first_parameter(), curve.last_parameter());
        self.length = self.length_of(curve, first_u, last_u);
    }

    /// `Length(C, FirstU, LastU)`: `GCPnts_AbscissaPoint::Length` to the tolerance.
    fn length_of(&self, curve: &dyn CurveAdaptor, first_u: f64, last_u: f64) -> f64 {
        length_between(curve, first_u, last_u, Some(self.tol_len))
    }

    /// `GetLength()`.
    pub(crate) fn get_length(&self) -> f64 {
        self.length
    }

    /// `GetSParameter(U)`: the arc-length parameter of the curve's parameter `u`.
    #[doc(alias = "Approx_CurvlinFunc::GetSParameter")]
    pub(crate) fn get_s_parameter(&self, u: f64) -> f64 {
        let curve = self.curve.get();
        let origin = curve.first_parameter();
        self.first_s + self.length_of(curve, origin, u) / self.length
    }

    /// `GetUParameter(C, S, 1)` on the curve held.
    pub(crate) fn get_u_parameter(&self, s: f64) -> Result<f64> {
        self.get_u_parameter_on(self.curve.get(), s)
    }

    /// `GetUParameter(C, S, 1)`: the parameter of `curve` at the arc-length parameter `s`,
    /// from the table of the curve at construction. An error where OCCT raises
    /// `Standard_ConstructionError` (`s` outside `[0, 1]`) or `StdFail_NotDone` (the
    /// solution not found).
    #[doc(alias = "Approx_CurvlinFunc::GetUParameter")]
    pub(crate) fn get_u_parameter_on(&self, curve: &dyn CurveAdaptor, s: f64) -> Result<f64> {
        if !(0.0..=1.0).contains(&s) {
            return Err(Error::InvalidInput("Approx_CurvlinFunc::GetUParameter: s outside [0, 1]".into()));
        }
        let (init_u, init_s) = (&self.ui, &self.si);
        let length = self.length;
        let nb_int = init_u.len() - 1;
        let n_interval = if s == 1.0 {
            nb_int - 1
        } else {
            let mut i = 0;
            while i < nb_int {
                if init_s[i] <= s && s < init_s[i + 1] {
                    break;
                }
                i += 1;
            }
            i
        };
        if s == init_s[n_interval] {
            return Ok(init_u[n_interval]);
        }
        if s == init_s[n_interval + 1] {
            return Ok(init_u[n_interval + 1]);
        }

        let base = init_u[n_interval];
        let delta_s = (s - init_s[n_interval]) * length;

        // An initial point.
        let (xi, yi) = find_four_points(n_interval, init_s, init_u, self.prev_s.get(), self.prev_u.get())?;
        let u_guess = cubic(s, &xi, &yi);

        let u = AbscissaPoint::with_guess_and_tolerance(curve, delta_s, base, u_guess, self.tol_len)
            .parameter()
            .ok_or_else(|| Error::Failed("Approx_CurvlinFunc::GetUParameter: GCPnts_AbscissaPoint not done".into()))?;
        self.prev_s.set(s);
        self.prev_u.set(u);
        Ok(u)
    }

    /// `EvalCase1(S, Order, Result)`: the point (order 0) or the first or second derivative
    /// with respect to the arc-length parameter at `s`; `None` for a higher order, where OCCT
    /// returns false with zeros.
    #[doc(alias = "Approx_CurvlinFunc::EvalCase1")]
    pub(crate) fn eval_case1(&self, s: f64, order: u32) -> Result<Option<[f64; 3]>> {
        let curve = self.curve.get();
        let u = self.get_u_parameter_on(curve, s)?;
        Ok(match order {
            0 => {
                let c = curve.value(u);
                Some([c.x, c.y, c.z])
            }
            1 => {
                let (_, dc_du) = curve.d1(u);
                let mag = dc_du.length();
                let du_ds = self.length / mag;
                let dc_ds = dc_du * du_ds;
                Some([dc_ds.x, dc_ds.y, dc_ds.z])
            }
            2 => {
                let (_, dc_du, d2c_du2) = curve.d2(u);
                let mag = dc_du.length();
                let du_ds = self.length / mag;
                let d2u_ds2 = -self.length * dc_du.dot(d2c_du2) * du_ds / (mag * mag * mag);
                let d2c_ds2 = d2c_du2 * du_ds * du_ds + dc_du * d2u_ds2;
                Some([d2c_ds2.x, d2c_ds2.y, d2c_ds2.z])
            }
            _ => None,
        })
    }
}
