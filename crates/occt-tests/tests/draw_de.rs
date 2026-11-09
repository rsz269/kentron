//! Runs OCCT's data-exchange tests (`de/`, copied from OCCT's `tests/de` by
//! `tools/de_import.py`): each STEP or IGES file read, written and read again, both readings'
//! statistics compared with OCCT's reference data (`occt_tests::de`), each as a test named
//! `de::<grid>::<case>`. A case listed in `expected_failures.txt` (as `tests/de/<grid>/<case>`)
//! passes while it fails and is reported when it starts to pass. The cases whose file is not
//! in OCCT's public dataset are not copied (`de/excluded.txt`). Each has the YAML cases' limit
//! on its processor time (`KENTRON_CASE_TIMEOUT`, 40 s). With `KENTRON_SKIP_EXPECTED_FAILURES`
//! set, the listed cases are ignored; with `KENTRON_EXPECTED_REASONS` set, why each fails.

use std::path::Path;

use libtest_mimic::{Arguments, Failed, Trial};

fn main() {
    let args = Arguments::from_args();
    std::panic::set_hook(Box::new(occt_tests::cases::quiet_panic_hook));
    let crate_dir = Path::new(env!("CARGO_MANIFEST_DIR"));
    let root = crate_dir.join("de");
    let expected: std::collections::HashSet<String> = std::fs::read_to_string(crate_dir.join("expected_failures.txt"))
        .unwrap_or_default()
        .lines()
        .map(|l| l.trim().to_string())
        .filter(|l| !l.is_empty() && !l.starts_with('#'))
        .collect();
    let skip_expected = occt_tests::skip_expected_failures();
    let limit: f64 = std::env::var("KENTRON_CASE_TIMEOUT").ok().and_then(|s| s.parse().ok()).unwrap_or(40.0);
    let mut trials = Vec::new();
    let mut grids: Vec<_> = std::fs::read_dir(&root).into_iter().flatten().flatten().map(|e| e.path()).filter(|p| p.is_dir()).collect();
    grids.sort();
    for grid in grids {
        let mut cases: Vec<_> = std::fs::read_dir(&grid).into_iter().flatten().flatten().map(|e| e.path()).filter(|p| p.is_file()).collect();
        cases.sort();
        for path in cases {
            let (g, c) = (grid.file_name().unwrap().to_string_lossy().to_string(), path.file_name().unwrap().to_string_lossy().to_string());
            let listed = expected.contains(&format!("tests/de/{g}/{c}"));
            let name = format!("de::{g}::{c}");
            let timed = name.clone();
            let trial = Trial::test(name, move || {
                let outcome = occt_tests::limit::run_limited(limit, move || {
                    let case = occt_tests::de::load(&path)?;
                    match occt_tests::de::run(&case)? {
                        differences if differences.is_empty() => Ok(()),
                        differences => Err(differences.iter().map(|(r, c)| format!("reference: {r}\n  kentron: {c}")).collect::<Vec<_>>().join("\n")),
                    }
                });
                let outcome = outcome.and_then(|o| o);
                if listed && std::env::var_os("KENTRON_EXPECTED_REASONS").is_some() {
                    if let Err(e) = &outcome {
                        eprintln!("expected failure {timed}: {}", e.lines().next().unwrap_or(""));
                    }
                }
                match (outcome, listed) {
                    (Ok(()), false) | (Err(_), true) => Ok(()),
                    (Err(e), false) => Err(Failed::from(e)),
                    (Ok(()), true) => Err(Failed::from("passes, but is listed in expected_failures.txt")),
                }
            });
            trials.push(trial.with_ignored_flag(skip_expected && listed));
        }
    }
    if let Some(list) = occt_tests::test_list() {
        trials.retain(|t| list.contains(t.name()));
    }
    libtest_mimic::run(&args, trials).exit();
}
