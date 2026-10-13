//! Default tolerances for geometric comparisons.

/// Distance below which two points are considered coincident (OCCT `Precision::Confusion`).
#[doc(alias = "Precision::Confusion")]
pub const LINEAR: f64 = 1e-7;

/// `LINEAR` squared, for comparing squared distances.
#[doc(alias = "Precision::SquareConfusion")]
pub const LINEAR_SQUARED: f64 = LINEAR * LINEAR;

/// Angle (radians) below which two directions are considered parallel (OCCT
/// `Precision::Angular`).
#[doc(alias = "Precision::Angular")]
pub const ANGULAR: f64 = 1e-12;

/// Tolerance for comparing curve and surface parameters (OCCT `Precision::PConfusion`).
#[doc(alias = "Precision::PConfusion")]
pub const PARAMETRIC: f64 = LINEAR * 0.01;

/// A value treated as infinite (OCCT `Precision::Infinite`).
#[doc(alias = "Precision::Infinite")]
pub const INFINITE: f64 = 2e100;

/// The smallest meaningful length: anything at or below it is treated as zero, e.g. when
/// deciding whether two points are distinct enough to define a direction (OCCT
/// `gp::Resolution`).
#[doc(alias = "gp::Resolution")]
pub const RESOLUTION: f64 = f64::MIN_POSITIVE;

thread_local! {
    /// Whether operations on this thread run with their looser tolerances (see
    /// [`loosened`]).
    static LOOSENED: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
}

/// Whether operations run with their looser tolerances: a retry for shapes whose points
/// meant to coincide are only nearly so (data written to few digits).
pub(crate) fn loosened() -> bool {
    LOOSENED.with(|l| l.get())
}

/// `f` run with the looser tolerances on (or off), as said; the setting restored after.
pub(crate) fn with_loosened<R>(loose: bool, f: impl FnOnce() -> R) -> R {
    let before = LOOSENED.with(|l| l.replace(loose));
    struct Restore(bool);
    impl Drop for Restore {
        fn drop(&mut self) {
            LOOSENED.with(|l| l.set(self.0));
        }
    }
    let _restore = Restore(before);
    f()
}
