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
