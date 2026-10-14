//! Running a test with a limit on its time: on a thread of its own, failing it if it runs
//! too long (the thread is left to finish or to end with the process). The limit is on the
//! processor time of the thread running it, which a machine loaded by other tests does not
//! stretch as it stretches wall time (and six times that in wall time, a backstop for a
//! machine loaded by all the test binaries at once). Its operations'
//! parallel steps run on a pool of its own, of one thread (tests already run side by side):
//! a test left running would otherwise hold the shared pool's threads and starve the others.

/// `f`'s result, or an error if it takes more than `limit` seconds of processor time (or
/// six times that in wall time). A panic in `f` is resumed here.
pub fn run_limited<T: Send + 'static>(limit: f64, f: impl FnOnce() -> T + Send + 'static) -> Result<T, String> {
    let (tx, rx) = std::sync::mpsc::channel();
    let (clock_tx, clock_rx) = std::sync::mpsc::channel();
    let handle = std::thread::spawn(move || {
        let pool = rayon::ThreadPoolBuilder::new().num_threads(1).build().expect("a thread pool");
        let _ = tx.send(pool.install(|| {
            let _ = clock_tx.send(ThreadClock::current());
            f()
        }));
    });
    let started = std::time::Instant::now();
    let clock = clock_rx.recv_timeout(std::time::Duration::from_secs(10)).ok().flatten();
    loop {
        match rx.recv_timeout(std::time::Duration::from_millis(250)) {
            Ok(result) => return Ok(result),
            Err(std::sync::mpsc::RecvTimeoutError::Timeout) => {
                let used = clock.as_ref().and_then(ThreadClock::seconds).unwrap_or_else(|| started.elapsed().as_secs_f64());
                if used > limit {
                    return Err(format!("timed out after {limit} s of processor time"));
                }
                if started.elapsed().as_secs_f64() > 6.0 * limit {
                    return Err(format!("timed out after {} s", 6.0 * limit));
                }
            }
            Err(std::sync::mpsc::RecvTimeoutError::Disconnected) => {
                return match handle.join() {
                    Err(panic) => std::panic::resume_unwind(panic),
                    Ok(()) => Err("the test ended without a result".into()),
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
