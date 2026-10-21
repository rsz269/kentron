//! Support code for tests converted from OCCT.
//!
//! * [`draw`]: a DRAW-like session used by tests generated from OCCT's Tcl scripts
//!   (`tests/draw_*`, produced by `tools/draw2rs.py`).
//! * [`gtest`]: assertion macros for the hand-translated GTests (`tests/gtest_*`).
//! * [`cases`]: the runner of test cases written in YAML (`cases/`, see `docs/test-cases.md`).

pub mod cases;
pub mod draw;
pub mod gtest;
pub mod limit;

/// The allocator of the test binaries (which link this library): mimalloc, as pykentron's
/// `fast` profile has it, a few percent quicker than the system's on these tests.
#[global_allocator]
static ALLOCATOR: mimalloc::MiMalloc = mimalloc::MiMalloc;

/// Whether tests expected to fail (`expected_failures.txt`, `expect: not_done`) are left out:
/// with `KENTRON_SKIP_EXPECTED_FAILURES` set (to anything but `0`), for a quicker run, which
/// does not notice one starting to pass.
pub fn skip_expected_failures() -> bool {
    std::env::var("KENTRON_SKIP_EXPECTED_FAILURES").is_ok_and(|v| v != "0")
}
