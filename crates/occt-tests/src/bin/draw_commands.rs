//! Runs DRAW commands on kentron's interpreter, for tools that need the shape a DRAW test
//! builds (`tools/gap_figures.py`): one command a line on standard input, its words split
//! at spaces. What the commands print goes to standard output; a failing command ends the
//! run with its error.
//!
//! `draw_commands < FILE` (a `save NAME FILE` among the commands writes a shape as BRep)

use std::io::Read;

use occt_tests::draw::{dispatch::dispatch, Draw};

fn main() {
    let mut text = String::new();
    std::io::stdin().read_to_string(&mut text).unwrap();
    let mut draw = Draw::new();
    for line in text.lines().map(str::trim).filter(|line| !line.is_empty() && !line.starts_with('#')) {
        let words: Vec<&str> = line.split_whitespace().collect();
        match dispatch(&mut draw, words[0], &words[1..]) {
            Some(Ok(output)) => print!("{output}"),
            Some(Err(error)) => {
                eprintln!("{line}: {error}");
                std::process::exit(1);
            }
            None => {
                eprintln!("{line}: not a command");
                std::process::exit(1);
            }
        }
    }
}
