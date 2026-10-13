//! Independent steps of an operation run on a thread pool (with the `parallel` feature),
//! their results in the order of their inputs, so that what follows is the same either way.

/// How long items are worked through here before the rest go to the pool: most steps (of
/// operations on a few faces) take less than the pool costs, in waking its threads and in
/// the clock speed the other cores' work takes from this one.
#[cfg(feature = "parallel")]
const ALONE: std::time::Duration = std::time::Duration::from_millis(2);

/// About how long each task given to the pool takes (items grouped to that).
#[cfg(feature = "parallel")]
const TASK: std::time::Duration = std::time::Duration::from_micros(50);

/// The pool: half the processors' threads (about a thread per core, and sparing the
/// slower cores of mixed processors, which keep a whole step waiting), or as many as
/// `RAYON_NUM_THREADS` says.
#[cfg(feature = "parallel")]
fn pool() -> &'static rayon::ThreadPool {
    static POOL: std::sync::OnceLock<rayon::ThreadPool> = std::sync::OnceLock::new();
    POOL.get_or_init(|| {
        let threads = std::env::var("RAYON_NUM_THREADS")
            .ok()
            .and_then(|s| s.parse().ok())
            .unwrap_or_else(|| std::thread::available_parallelism().map_or(1, |n| (n.get() / 2).max(1)));
        rayon::ThreadPoolBuilder::new().num_threads(threads).thread_name(|i| format!("kentron-{i}")).build().expect("a thread pool")
    })
}

/// `f` of each item, in order.
#[cfg(feature = "parallel")]
pub(crate) fn map<T: Sync, R: Send>(items: &[T], f: impl Fn(&T) -> R + Sync + Send) -> Vec<R> {
    use rayon::prelude::*;
    let start = std::time::Instant::now();
    let mut out = Vec::with_capacity(items.len());
    for (i, item) in items.iter().enumerate() {
        out.push(f(item));
        let spent = start.elapsed();
        if spent > ALONE && i + 1 < items.len() {
            // The rest in groups taking about TASK each, judging by the items so far: on
            // the pool this runs in (a caller's own, as a test's), else on ours.
            let each = spent / (i as u32 + 1);
            let group = (TASK.as_nanos() / each.as_nanos().max(1)).max(1) as usize;
            // (On the pool's threads with the caller's tolerances.)
            let loose = crate::tolerance::loosened();
            let rest = || items[i + 1..].par_iter().with_min_len(group).map(|x| crate::tolerance::with_loosened(loose, || f(x))).collect::<Vec<R>>();
            out.extend(if rayon::current_thread_index().is_some() { rest() } else { pool().install(rest) });
            break;
        }
    }
    out
}

/// `f` of each item, in order.
#[cfg(not(feature = "parallel"))]
pub(crate) fn map<T: Sync, R: Send>(items: &[T], f: impl Fn(&T) -> R + Sync + Send) -> Vec<R> {
    items.iter().map(f).collect()
}
