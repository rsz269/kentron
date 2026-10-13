//! Assertion helpers mirroring the GoogleTest macros used by OCCT's GTests.
//!
//! `EXPECT_*` and `ASSERT_*` both map to panicking assertions: the first failure ends the test.

/// `EXPECT_NEAR(a, b, tol)`.
#[macro_export]
macro_rules! assert_near {
    ($a:expr, $b:expr, $tol:expr $(,)?) => {
        $crate::assert_near!($a, $b, $tol, "")
    };
    ($a:expr, $b:expr, $tol:expr, $($msg:tt)+) => {{
        let (a, b, tol): (f64, f64, f64) = ($a, $b, $tol);
        assert!(
            (a - b).abs() <= tol,
            "assert_near failed: |{} - {}| = {} > {}: {}",
            a, b, (a - b).abs(), tol, format_args!($($msg)+)
        );
    }};
}

/// `EXPECT_DOUBLE_EQ(a, b)`: equal within 4 ULPs, as in GoogleTest.
#[macro_export]
macro_rules! assert_double_eq {
    ($a:expr, $b:expr $(,)?) => {
        $crate::assert_double_eq!($a, $b, "")
    };
    ($a:expr, $b:expr, $($msg:tt)+) => {{
        let (a, b): (f64, f64) = ($a, $b);
        assert!(
            $crate::gtest::ulps_eq(a, b, 4),
            "assert_double_eq failed: {} != {}: {}",
            a, b, format_args!($($msg)+)
        );
    }};
}

/// `EXPECT_THROW(stmt, Exception)` / `EXPECT_ANY_THROW(stmt)`: OCCT exceptions are panics.
///
/// A `todo!()` panic does not count: unimplemented shims must not make these tests pass.
#[macro_export]
macro_rules! assert_panics {
    ($e:expr $(,)?) => {{
        let r = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            let _ = $e;
        }));
        match r {
            Ok(_) => panic!("expected a panic from `{}`", stringify!($e)),
            Err(payload) => $crate::gtest::reject_todo(payload),
        }
    }};
}

/// Re-raises a panic caught by [`assert_panics!`] if it came from `todo!()`.
pub fn reject_todo(payload: Box<dyn std::any::Any + Send>) {
    let msg = payload
        .downcast_ref::<String>()
        .map(String::as_str)
        .or_else(|| payload.downcast_ref::<&str>().copied())
        .unwrap_or("");
    if msg.starts_with("not yet implemented") {
        std::panic::resume_unwind(payload);
    }
}

/// Equality within `max_ulps` units in the last place (GoogleTest's `AlmostEquals`).
pub fn ulps_eq(a: f64, b: f64, max_ulps: u64) -> bool {
    if a.is_nan() || b.is_nan() {
        return false;
    }
    fn biased(x: f64) -> u64 {
        let bits = x.to_bits();
        const SIGN: u64 = 1 << 63;
        if bits & SIGN != 0 { !bits + 1 } else { bits | SIGN }
    }
    let (x, y) = (biased(a), biased(b));
    x.max(y) - x.min(y) <= max_ulps
}
