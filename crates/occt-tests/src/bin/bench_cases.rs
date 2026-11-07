//! Times YAML cases' operations and checks for `tools/bench_cases.py`: for each file given, each
//! case `REPEATS` times (default 3), a line of JSON with the least of each step's process
//! CPU time and wall time, in nanoseconds.
//!
//! `bench_cases FILE...` (`ONLY`: a file of case names, one a line, to run alone)

use std::time::Instant;

use occt_tests::cases;

fn cpu_nanos() -> u64 {
    let mut t = libc::timespec { tv_sec: 0, tv_nsec: 0 };
    // SAFETY: `t` is a valid timespec for the call to write.
    unsafe { libc::clock_gettime(libc::CLOCK_PROCESS_CPUTIME_ID, &mut t) };
    t.tv_sec as u64 * 1_000_000_000 + t.tv_nsec as u64
}

fn main() {
    let repeats: usize = std::env::var("REPEATS").ok().and_then(|s| s.parse().ok()).unwrap_or(3);
    let only: Option<Vec<String>> = std::env::var("ONLY").ok().map(|f| std::fs::read_to_string(f).unwrap().lines().map(str::to_string).collect());
    let start = Instant::now();
    let clock = || [cpu_nanos(), start.elapsed().as_nanos() as u64];
    for file in std::env::args().skip(1) {
        let stem = std::path::Path::new(&file).file_stem().unwrap().to_string_lossy().to_string();
        for case in cases::load(std::path::Path::new(&file)).unwrap() {
            let name = format!("{stem}::{}", case.name);
            if only.as_ref().is_some_and(|o| !o.contains(&name)) {
                continue;
            }
            let mut best: Option<Vec<(String, [u64; 2])>> = None;
            let mut failed = None;
            for _ in 0..repeats {
                match cases::timed(&case, clock) {
                    Ok(times) => match &mut best {
                        None => best = Some(times),
                        Some(b) => b.iter_mut().zip(times).for_each(|(b, t)| {
                            b.1 = [b.1[0].min(t.1[0]), b.1[1].min(t.1[1])];
                        }),
                    },
                    Err(e) => failed = Some(e),
                }
            }
            match (best, failed) {
                (Some(steps), None) => {
                    let steps: Vec<String> = steps.iter().map(|(n, t)| format!("[\"{n}\",{},{}]", t[0], t[1])).collect();
                    println!("{{\"case\":\"{name}\",\"steps\":[{}]}}", steps.join(","));
                }
                (_, e) => println!("{{\"case\":\"{name}\",\"error\":{:?}}}", e.unwrap_or_default()),
            }
        }
    }
}
