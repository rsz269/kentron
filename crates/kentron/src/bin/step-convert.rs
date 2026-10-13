//! Converts STEP files between text (`.step`, `.stp`) and binary (`.binstep`), value for
//! value the same.
//!
//! ```text
//! step-convert INPUT [OUTPUT]
//! ```
//!
//! The input's form is told from its first bytes; the output's from its extension
//! (`.binstep`: binary, else text). Without an output: the input's name with `.binstep` for
//! text input, `.step` for binary. The result is read back and checked against the input.

use std::path::{Path, PathBuf};
use std::process::ExitCode;

use kentron::io::step::{binstep, convert, exchange_of};

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    if args.is_empty() || args.len() > 2 || args.iter().any(|a| a == "-h" || a == "--help") {
        eprintln!("usage: step-convert INPUT [OUTPUT]   (.binstep <-> .step/.stp, by OUTPUT's extension)");
        return ExitCode::from(2);
    }
    match run(Path::new(&args[0]), args.get(1).map(PathBuf::from)) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("step-convert: {e}");
            ExitCode::FAILURE
        }
    }
}

fn run(input: &Path, output: Option<PathBuf>) -> Result<(), String> {
    let bytes = std::fs::read(input).map_err(|e| format!("{}: {e}", input.display()))?;
    let from_binary = binstep::is_binstep(&bytes);
    let output = output.unwrap_or_else(|| input.with_extension(if from_binary { "step" } else { "binstep" }));
    let to_binary = output.extension().is_some_and(|e| e.eq_ignore_ascii_case("binstep"));
    let start = std::time::Instant::now();
    let converted = convert(&bytes, to_binary).map_err(|e| format!("{}: {e}", input.display()))?;
    let took = start.elapsed();
    // The same values both ways (as text, the written form being one per value).
    let (before, after) = (exchange_of(&bytes), exchange_of(&converted));
    match (before, after) {
        (Ok(a), Ok(b)) if a.to_text() == b.to_text() => {}
        _ => return Err(format!("{}: the result does not read back the same", output.display())),
    }
    std::fs::write(&output, &converted).map_err(|e| format!("{}: {e}", output.display()))?;
    let form = |binary: bool| if binary { "binstep" } else { "STEP text" };
    println!(
        "{} ({}, {} bytes) -> {} ({}, {} bytes) in {:.3} s",
        input.display(),
        form(from_binary),
        bytes.len(),
        output.display(),
        form(to_binary),
        converted.len(),
        took.as_secs_f64()
    );
    Ok(())
}
