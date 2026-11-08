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

/// The tests to run, if only some: with `KENTRON_TEST_LIST` naming a file of test names (one a
/// line), as `tools/run_tests.py --only` gives the custom harnesses (their runner takes a
/// single filter); the others are left out, as filtered out.
pub fn test_list() -> Option<std::collections::HashSet<String>> {
    let file = std::env::var_os("KENTRON_TEST_LIST")?;
    let text = std::fs::read_to_string(&file).unwrap_or_else(|e| panic!("KENTRON_TEST_LIST {file:?}: {e}"));
    Some(text.lines().map(|l| l.trim().to_string()).filter(|l| !l.is_empty()).collect())
}
