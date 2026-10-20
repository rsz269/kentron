//! Runs the YAML test cases in `cases/` (see `docs/test-cases.md`), each as a test named
//! `<file>::<case>`. A case marked `expect: not_done` passes while it fails and is reported
//! when it starts to pass; `skip` cases are ignored. A case taking more than
//! `KENTRON_CASE_TIMEOUT` seconds of processor time (40 by default: the slowest passing case
//! takes under 20 s), or six times that in wall time, fails; with `KENTRON_CASE_TIMES` set, each case's time is reported. With
//! `KENTRON_CASE_ORDER` naming a file of such times (`<seconds> <case>` lines, as
//! `tools/run_tests.py` keeps them), the slowest cases start first, so that the run does
//! not end waiting for a few long ones started last. With
//! `KENTRON_SKIP_EXPECTED_FAILURES` set, `expect: not_done` cases are ignored (a quicker run,
//! which does not notice one starting to pass); with `KENTRON_EXPECTED_REASONS` set, each
//! such case's failure is reported.

use std::path::Path;

use libtest_mimic::{Arguments, Failed, Trial};
use occt_tests::cases;

fn main() {
    let args = Arguments::from_args();
    // The runner turns panics (kentron's `todo!`s) into failures with their message.
    std::panic::set_hook(Box::new(cases::quiet_panic_hook));
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("cases");
    let skip_expected = occt_tests::skip_expected_failures();
    let mut trials = Vec::new();
    for file in cases::files(&root) {
        let prefix = file.strip_prefix(&root).unwrap_or(&file).with_extension("").to_string_lossy().replace('/', "::");
        let loaded = match cases::load(&file) {
            Ok(c) => c,
            Err(e) => {
                trials.push(Trial::test(prefix.clone(), move || Err(Failed::from(e))));
                continue;
            }
        };
        for case in loaded {
            let name = format!("{prefix}::{}", case.name);
            let ignored = case.skip.is_some() || (skip_expected && case.not_done);
            let timed = name.clone();
            let trial = Trial::test(name, move || {
                let not_done = case.not_done;
                // (With `KENTRON_CASE_TIMES` set, each case's time is reported.)
                let started = std::time::Instant::now();
                let outcome = run_with_timeout(case);
                if std::env::var_os("KENTRON_CASE_TIMES").is_some() {
                    eprintln!("case time {:.3}s {timed}", started.elapsed().as_secs_f64());
                }
                // (With `KENTRON_EXPECTED_REASONS` set, why each expected failure fails.)
                if let (Err(e), true) = (&outcome, not_done) {
                    if std::env::var_os("KENTRON_EXPECTED_REASONS").is_some() {
                        eprintln!("expected failure {timed}: {}", e.lines().next().unwrap_or(""));
                    }
                }
                match (outcome, not_done) {
                    (Ok(()), false) | (Err(_), true) => Ok(()),
                    (Err(e), false) => Err(Failed::from(e)),
                    (Ok(()), true) => Err(Failed::from("passes, but is marked `expect: not_done`")),
                }
            });
            trials.push(trial.with_ignored_flag(ignored));
        }
    }
    if let Some(times) = std::env::var_os("KENTRON_CASE_ORDER").and_then(|f| std::fs::read_to_string(f).ok()) {
        let seconds: std::collections::HashMap<&str, f64> = times
            .lines()
            .filter_map(|l| {
                let (t, name) = l.split_once(' ')?;
                Some((name.trim(), t.trim_end_matches('s').parse().ok()?))
            })
            .collect();
        // (A stable sort: cases without a time keep their order, after the timed ones.)
        trials.sort_by(|a, b| seconds.get(b.name()).unwrap_or(&0.0).total_cmp(seconds.get(a.name()).unwrap_or(&0.0)));
    }
    libtest_mimic::run(&args, trials).exit();
}

/// Runs a case on a thread of its own with its time limit (see [`occt_tests::limit`]).
fn run_with_timeout(case: cases::Case) -> Result<(), String> {
    // (An explicit KENTRON_CASE_TIMEOUT overrides a case's own `timeout`.)
    let limit = std::env::var("KENTRON_CASE_TIMEOUT").ok().and_then(|s| s.parse().ok()).or(case.timeout).unwrap_or(40.0);
    occt_tests::limit::run_limited(limit, move || cases::run(&case).map_err(|e| e.to_string()))?
}
