//! Runs the YAML test cases in `cases/` (see `docs/test-cases.md`), each as a test named
//! `<file>::<case>`. A case marked `expect: not_done` passes while it fails and is reported
//! when it starts to pass; `skip` cases are ignored. A case taking more than
//! `KENTRON_CASE_TIMEOUT` seconds of processor time (40 by default: the slowest passing case
//! takes under 20 s), or three times that in wall time, fails; with `KENTRON_CASE_TIMES` set, each case's time is reported. With
//! `KENTRON_CASE_ORDER` naming a file of such times (`<seconds> <case>` lines, as
//! `tools/run_tests.py` keeps them), the slowest cases start first, so that the run does
//! not end waiting for a few long ones started last.

use std::path::Path;

use libtest_mimic::{Arguments, Failed, Trial};
use occt_tests::cases;

fn main() {
    let args = Arguments::from_args();
    // The runner turns panics (kentron's `todo!`s) into failures with their message.
    std::panic::set_hook(Box::new(cases::quiet_panic_hook));
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("cases");
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
            let ignored = case.skip.is_some();
            let timed = name.clone();
            let trial = Trial::test(name, move || {
                let not_done = case.not_done;
                // (With `KENTRON_CASE_TIMES` set, each case's time is reported.)
                let started = std::time::Instant::now();
                let outcome = run_with_timeout(case);
                if std::env::var_os("KENTRON_CASE_TIMES").is_some() {
                    eprintln!("case time {:.3}s {timed}", started.elapsed().as_secs_f64());
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

/// Runs a case on a thread of its own, failing it if it runs too long (the thread is left
/// to finish or to end with the process): by the processor time of the thread running it,
/// which a machine loaded by other tests does not stretch as it stretches wall time. Its operations' parallel steps run on a pool of
/// its own, of one thread (cases already run side by side): a case left running would
/// otherwise hold the shared pool's threads and starve the others.
fn run_with_timeout(case: cases::Case) -> Result<(), String> {
    // (An explicit KENTRON_CASE_TIMEOUT overrides a case's own `timeout`.)
    let limit = std::env::var("KENTRON_CASE_TIMEOUT").ok().and_then(|s| s.parse().ok()).or(case.timeout).unwrap_or(40.0);
    let (tx, rx) = std::sync::mpsc::channel();
    let (clock_tx, clock_rx) = std::sync::mpsc::channel();
    let handle = std::thread::spawn(move || {
        let pool = rayon::ThreadPoolBuilder::new().num_threads(1).build().expect("a thread pool");
        let _ = tx.send(pool.install(|| {
            let _ = clock_tx.send(ThreadClock::current());
            cases::run(&case).map_err(|e| e.to_string())
        }));
    });
    let started = std::time::Instant::now();
    let clock = clock_rx.recv_timeout(std::time::Duration::from_secs(10)).ok().flatten();
    loop {
        match rx.recv_timeout(std::time::Duration::from_millis(250)) {
            Ok(result) => return result,
            Err(std::sync::mpsc::RecvTimeoutError::Timeout) => {
                let used = clock.as_ref().and_then(ThreadClock::seconds).unwrap_or_else(|| started.elapsed().as_secs_f64());
                if used > limit {
                    return Err(format!("timed out after {limit} s of processor time"));
                }
                if started.elapsed().as_secs_f64() > 3.0 * limit {
                    return Err(format!("timed out after {} s", 3.0 * limit));
                }
            }
            Err(std::sync::mpsc::RecvTimeoutError::Disconnected) => {
                return match handle.join() {
                    Err(panic) => std::panic::resume_unwind(panic),
                    Ok(()) => Err("the case ended without a result".into()),
                };
            }
        }
    }
}

/// The processor time clock of a thread (where the platform has one).
struct ThreadClock(#[cfg(unix)] libc::clockid_t);

impl ThreadClock {
    /// The calling thread's clock.
    fn current() -> Option<ThreadClock> {
        #[cfg(unix)]
        {
            let mut id: libc::clockid_t = 0;
            // SAFETY: pthread_self is the calling thread; id is a valid place to write to.
            (unsafe { libc::pthread_getcpuclockid(libc::pthread_self(), &mut id) } == 0).then_some(ThreadClock(id))
        }
        #[cfg(not(unix))]
        None
    }

    /// The thread's processor time so far; `None` once it is gone.
    fn seconds(&self) -> Option<f64> {
        #[cfg(unix)]
        {
            let mut t = libc::timespec { tv_sec: 0, tv_nsec: 0 };
            // SAFETY: the clock id came from pthread_getcpuclockid; t is a valid timespec.
            (unsafe { libc::clock_gettime(self.0, &mut t) } == 0).then(|| t.tv_sec as f64 + t.tv_nsec as f64 * 1e-9)
        }
        #[cfg(not(unix))]
        None
    }
}
