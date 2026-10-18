//! Runs OCCT's DRAW scripts that are not converted to Rust (`scripts/`, copied by
//! `tools/draw2rs.py`) on the Tcl interpreter (`occt_tests::draw::tcl`), each as a test
//! named `<group>::<grid>::<case>`. A script listed in `expected_failures.txt` passes while
//! it fails and is reported when it starts to pass. Scripts needing a data file that is
//! missing, or declaring a known failure of OCCT's (a TODO for this platform), are ignored.
//! Each has the YAML cases' limit on its processor time (`KENTRON_CASE_TIMEOUT`, 40 s).
//! With `KENTRON_SKIP_EXPECTED_FAILURES` set, the listed scripts are ignored; with
//! `KENTRON_CASE_TIMES` set, each script's time is reported.

use std::path::{Path, PathBuf};

use libtest_mimic::{Arguments, Failed, Trial};
use occt_tests::draw::tcl::{run_script, Outcome};

fn main() {
    let args = Arguments::from_args();
    std::panic::set_hook(Box::new(occt_tests::cases::quiet_panic_hook));
    let crate_dir = Path::new(env!("CARGO_MANIFEST_DIR"));
    let root = crate_dir.join("scripts");
    let expected: std::collections::HashSet<String> = std::fs::read_to_string(crate_dir.join("expected_failures.txt"))
        .unwrap_or_default()
        .lines()
        .map(|l| l.trim().to_string())
        .filter(|l| !l.is_empty() && !l.starts_with('#'))
        .collect();
    let data_file = regex::Regex::new(r"locate_data_file\s+([^\]\s}]+)").expect("a pattern");
    let skip_expected = occt_tests::skip_expected_failures();
    let limit: f64 = std::env::var("KENTRON_CASE_TIMEOUT").ok().and_then(|s| s.parse().ok()).unwrap_or(40.0);
    let mut trials = Vec::new();
    for (group, grid, case) in scripts(&root) {
        let text = std::fs::read_to_string(root.join(&group).join(&grid).join(&case)).unwrap_or_default();
        let missing = data_file.captures_iter(&text).map(|c| c[1].to_string()).find(|f| !f.contains('$') && occt_tests::draw::find_data_file(f).is_none());
        let known = occt_tests::draw::tcl::known_failure(&text);
        let listed = expected.contains(&format!("tests/{group}/{grid}/{case}"));
        let name = format!("{group}::{grid}::{case}");
        let ignored = missing.is_some() || known.is_some() || (skip_expected && listed);
        let root = root.clone();
        let timed = name.clone();
        let trial = Trial::test(name, move || {
            // (With `KENTRON_CASE_TIMES` set, each script's time is reported, as the cases'.)
            let started = std::time::Instant::now();
            let outcome = occt_tests::limit::run_limited(limit, move || run_script(&root, &group, &grid, &case));
            if std::env::var_os("KENTRON_CASE_TIMES").is_some() {
                eprintln!("case time {:.3}s {timed}", started.elapsed().as_secs_f64());
            }
            match (outcome, listed) {
                (Ok(Outcome::Passed | Outcome::Skipped(_) | Outcome::KnownFailure(_)), false) => Ok(()),
                (Ok(Outcome::Failed(_)) | Err(_), true) => Ok(()),
                (Ok(Outcome::Failed(lines)), false) => Err(Failed::from(lines.join("\n"))),
                (Err(e), false) => Err(Failed::from(e)),
                (Ok(_), true) => Err(Failed::from("passes, but is listed in expected_failures.txt")),
            }
        });
        trials.push(trial.with_ignored_flag(ignored));
    }
    libtest_mimic::run(&args, trials).exit();
}

/// The scripts under `root`: `<group>/<grid>/<case>`, without the groups' and grids' own
/// files.
fn scripts(root: &Path) -> Vec<(String, String, String)> {
    let mut out = Vec::new();
    let dirs = |p: &Path| -> Vec<PathBuf> {
        let mut v: Vec<PathBuf> = std::fs::read_dir(p).into_iter().flatten().flatten().map(|e| e.path()).filter(|p| p.is_dir()).collect();
        v.sort();
        v
    };
    for group in dirs(root) {
        for grid in dirs(&group) {
            let mut cases: Vec<PathBuf> = std::fs::read_dir(&grid).into_iter().flatten().flatten().map(|e| e.path()).filter(|p| p.is_file()).collect();
            cases.sort();
            for case in cases {
                let name = case.file_name().map(|n| n.to_string_lossy().to_string()).unwrap_or_default();
                if ["begin", "end", "parse.rules"].contains(&name.as_str()) {
                    continue;
                }
                let part = |p: &Path| p.file_name().map(|n| n.to_string_lossy().to_string()).unwrap_or_default();
                out.push((part(&group), part(&grid), name));
            }
        }
    }
    out
}
