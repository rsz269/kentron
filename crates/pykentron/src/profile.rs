//! Profiles trading memory for speed: which allocator kentron's memory comes from, and how
//! much the model's garbage collection frees.
//!
//! | profile | allocator | collection |
//! |---|---|---|
//! | `low-memory` | the system's | also forgets values derived from removed nodes |
//! | `balanced` | the system's | removes unreachable nodes, keeps derived values |
//! | `fast` (default) | mimalloc | as `balanced` |
//!
//! (Built without the `mimalloc` feature, `balanced` is the default.)
//!
//! The allocator is chosen once, at the first allocation (memory must be freed by the
//! allocator that gave it), from `PYKENTRON_PROFILE` in the environment: the profile is set
//! before `pykentron` is imported, or later with [`set_profile`] between the two profiles that
//! share the allocator in use.

use std::ffi::{c_char, CStr};
use std::sync::atomic::{AtomicU8, Ordering};

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

/// The environment variable naming the profile to start with.
const VARIABLE: &CStr = c"PYKENTRON_PROFILE";

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Profile {
    LowMemory = 1,
    Balanced = 2,
    Fast = 3,
}

impl Profile {
    const ALL: [Profile; 3] = [Profile::LowMemory, Profile::Balanced, Profile::Fast];

    fn name(self) -> &'static str {
        match self {
            Profile::LowMemory => "low-memory",
            Profile::Balanced => "balanced",
            Profile::Fast => "fast",
        }
    }

    fn named(name: &[u8]) -> Option<Profile> {
        Profile::ALL.into_iter().find(|p| p.name().as_bytes() == name)
    }

    /// Whether the profile's memory comes from mimalloc.
    fn uses_mimalloc(self) -> bool {
        self == Profile::Fast
    }

    /// Whether garbage collection also forgets the values derived from removed nodes.
    pub(crate) fn forgets_derived(self) -> bool {
        self == Profile::LowMemory
    }
}

/// The profile `PYKENTRON_PROFILE` names: `None` if it is not set, `Err` with its value if it
/// names none. (Read without allocating: the allocator is chosen from it.)
fn from_environment() -> Option<Result<Profile, &'static [u8]>> {
    unsafe extern "C" {
        fn getenv(name: *const c_char) -> *const c_char;
    }
    // SAFETY: getenv takes a NUL-terminated name and returns NULL or a NUL-terminated string
    // that lives as long as the environment entry (which pykentron never changes).
    let value = unsafe { getenv(VARIABLE.as_ptr()) };
    if value.is_null() {
        return None;
    }
    let value = unsafe { CStr::from_ptr(value) }.to_bytes();
    Some(Profile::named(value).ok_or(value))
}

/// The profile in use, 0 until it is first asked for.
static CURRENT: AtomicU8 = AtomicU8::new(0);

/// The profile used where the environment names none: `fast` where built with mimalloc.
const DEFAULT: Profile = if cfg!(feature = "mimalloc") { Profile::Fast } else { Profile::Balanced };

/// The profile in use: the environment's at first ([`DEFAULT`] if it names none).
pub(crate) fn current() -> Profile {
    match CURRENT.load(Ordering::Relaxed) {
        0 => {
            let profile = match from_environment() {
                Some(Ok(p)) if !p.uses_mimalloc() || cfg!(feature = "mimalloc") => p,
                Some(Ok(_)) => Profile::Balanced,
                _ => DEFAULT,
            };
            CURRENT.store(profile as u8, Ordering::Relaxed);
            profile
        }
        n => Profile::ALL[n as usize - 1],
    }
}

/// Fails if `PYKENTRON_PROFILE` names no profile pykentron has (checked when it is imported).
pub(crate) fn check_environment() -> PyResult<()> {
    match from_environment() {
        Some(Err(value)) => Err(PyValueError::new_err(format!(
            "PYKENTRON_PROFILE={}: not a profile (low-memory, balanced, fast)",
            String::from_utf8_lossy(value)
        ))),
        Some(Ok(p)) if p.uses_mimalloc() && !cfg!(feature = "mimalloc") => {
            Err(PyValueError::new_err("PYKENTRON_PROFILE=fast: this pykentron is built without mimalloc (its `mimalloc` feature)"))
        }
        _ => Ok(()),
    }
}

/// The profile in use: "low-memory", "balanced" or "fast".
#[pyfunction]
pub(crate) fn profile() -> &'static str {
    current().name()
}

/// Switches to the profile `name`: "low-memory" or "balanced" (which use the system's
/// allocator), or "fast" (mimalloc). The allocator is chosen at startup, so a switch to a
/// profile with another one fails: set `PYKENTRON_PROFILE` before importing pykentron instead.
#[pyfunction]
pub(crate) fn set_profile(name: &str) -> PyResult<()> {
    let profile = Profile::named(name.as_bytes())
        .ok_or_else(|| PyValueError::new_err(format!("{name}: not a profile (low-memory, balanced, fast)")))?;
    let current = current();
    if profile.uses_mimalloc() != current.uses_mimalloc() {
        return Err(PyValueError::new_err(format!(
            "{name}: needs another allocator than {} (it is chosen at startup): set PYKENTRON_PROFILE={name} before importing pykentron",
            current.name()
        )));
    }
    CURRENT.store(profile as u8, Ordering::Relaxed);
    Ok(())
}

#[cfg(feature = "mimalloc")]
mod allocator {
    use std::alloc::{GlobalAlloc, Layout, System};
    use std::sync::atomic::{AtomicU8, Ordering};

    use mimalloc::MiMalloc;

    /// mimalloc or the system's allocator, as the profile at the first allocation says.
    struct Chosen;

    /// 0 until chosen, then 1 for the system's allocator, 2 for mimalloc.
    static CHOSEN: AtomicU8 = AtomicU8::new(0);

    fn mimalloc() -> bool {
        match CHOSEN.load(Ordering::Relaxed) {
            0 => {
                // (Threads allocating first at once all choose the same.)
                let chosen = match super::from_environment() {
                    Some(Ok(p)) => p.uses_mimalloc(),
                    // (A name of no profile fails the import: the default's allocator meanwhile.)
                    _ => super::DEFAULT.uses_mimalloc(),
                };
                CHOSEN.store(if chosen { 2 } else { 1 }, Ordering::Relaxed);
                chosen
            }
            n => n == 2,
        }
    }

    // SAFETY: every call goes to one allocator, the same for the whole process.
    unsafe impl GlobalAlloc for Chosen {
        unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
            if mimalloc() { MiMalloc.alloc(layout) } else { System.alloc(layout) }
        }

        unsafe fn alloc_zeroed(&self, layout: Layout) -> *mut u8 {
            if mimalloc() { MiMalloc.alloc_zeroed(layout) } else { System.alloc_zeroed(layout) }
        }

        unsafe fn dealloc(&self, ptr: *mut u8, layout: Layout) {
            if mimalloc() { MiMalloc.dealloc(ptr, layout) } else { System.dealloc(ptr, layout) }
        }

        unsafe fn realloc(&self, ptr: *mut u8, layout: Layout, new_size: usize) -> *mut u8 {
            if mimalloc() { MiMalloc.realloc(ptr, layout, new_size) } else { System.realloc(ptr, layout, new_size) }
        }
    }

    #[global_allocator]
    static ALLOCATOR: Chosen = Chosen;
}
