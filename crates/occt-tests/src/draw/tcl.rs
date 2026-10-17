//! Running OCCT's DRAW scripts as they are, Tcl and all, on a Tcl interpreter (Molt), its
//! commands those of [`Draw`] (by [`dispatch`](super::dispatch::dispatch)) and the Tcl
//! commands Molt lacks (`regexp`, `regsub`, `lsearch`, `lsort`, `format`, `source`, `file`,
//! DRAW's `repeat`, ...).
//!
//! A script runs as OCCT's test harness runs it: after its group's and grid's `begin` files,
//! before their `end` files, with the harness's variables set. What it prints, and the
//! results of the DRAW commands it runs (which DRAW echoes), make its log; it fails where a
//! line of the log matches a `FAILED` rule of OCCT's `parse.rules` (an error, an exception,
//! a failure), or where a command fails outside a `catch`.

use std::path::{Path, PathBuf};

use molt_forked::prelude::*;
use molt_forked::{molt_err, molt_ok};

use super::{find_data_file, Draw, DrawError};

/// What a script run gives.
#[derive(Debug, Clone, PartialEq)]
pub enum Outcome {
    Passed,
    /// The failures: lines of the log matching a `FAILED` rule, or the error that ended it.
    Failed(Vec<String>),
    /// A data file the script needs is missing.
    Skipped(String),
    /// The script declares an OCCT TODO (a known failure of OCCT's) that applies on Linux.
    KnownFailure(String),
}

/// The state of a script run: DRAW's, and the log.
pub struct Session {
    pub draw: Draw,
    pub log: Vec<String>,
    /// Directories `source` looks in (after the path as given).
    pub source_dirs: Vec<PathBuf>,
    /// A data file found missing.
    pub missing: Option<String>,
    /// The OCCT TODOs the script declares (`puts "TODO id platforms: message"`) that apply.
    pub todos: Vec<String>,
    pub exited: bool,
}

/// Commands that only show or time things (DRAW's viewers and the like): they do nothing.
const IGNORED: &[&str] = &[
    "pload", "cpulimit", "checkview", "smallview", "fit", "2dfit", "donly", "display", "erase", "clear", "axo", "top", "front",
    "left", "right", "bottom", "back", "xwd", "don", "isos", "protect", "unprotect", "dchrono", "chrono", "dlog", "decho",
    "dtracelevel", "view", "wait", "meminfo", "checkcolor", "trinfo_display", "dversion", "dsetsignal", "set_display",
    "autodisplay", "repaint", "dflush", "dbreak", "datadir", "vinit", "pu", "pd", "pl", "pr", "mu", "md", "mu4", "md4", "2dpu",
    "2dpd", "2dpl", "2dpr", "2dmu", "2dmd", "2dzoom", "zoom", "wzoom", "vclear", "vdisplay", "vfit", "vtop",
    "vaxo", "vfront", "vright", "vleft", "vback", "vbottom", "vremove", "verase", "vzoom", "vrotate", "vsetcolor", "vaspects",
    "vcamera", "vviewparams", "vrenderparams", "vlight", "vdrawparray", "vtrihedron", "vzbufftrihedron", "vfps", "vglinfo",
    "vcaps", "vsetmaterial", "vsettransparency", "vshowfaceboundary", "vdefaults", "vclose", "vraytrace", "vupdate", "vredraw",
    "vexport", "vdump", "dtyp", "disp", "pick", "checkdisplay", "tclean_display", "checktrinfo_display", "chkfi", "dhistory",
    "dfreemem", "dmalloc", "dputs", "dsetenv", "testgrid", "help",
    // (Boolean operations' history, and fillets' tolerances: kentron keeps no history and
    // uses its own tolerances.)
    "setfillhistory", "tolblend",
];

/// The interpreter's commands: Tcl's (Molt's, and ours where it lacks them), then DRAW's.
fn execute(name: &str, interp: &mut Interp<Session>, argv: &[Value]) -> MoltResult {
    if let Some(r) = extra(name, interp, argv) {
        return r;
    }
    if let Some(r) = molt_forked::__private::execute_builtin(name, interp, argv) {
        return r;
    }
    if let Some(r) = interp.try_execute_proc(name, argv) {
        return r;
    }
    if IGNORED.contains(&name) {
        return molt_ok!();
    }
    // (As in DRAW, a DRAW variable is a Tcl variable holding its own name: an argument that
    // is a Tcl variable, as one linked by `upvar`, holding a DRAW variable's name stands
    // for that.)
    let resolved: Vec<String> = argv[1..]
        .iter()
        .map(|v| {
            let name = v.as_str();
            if interp.context_mut().draw.has(name) {
                return name.to_string();
            }
            match interp.var(v) {
                Ok(value) if interp.context_mut().draw.has(value.as_str()) => value.as_str().to_string(),
                _ => name.to_string(),
            }
        })
        .collect();
    let args: Vec<&str> = resolved.iter().map(String::as_str).collect();
    let ran = run_draw_command(interp, name, &args);
    // The DRAW variables made are Tcl variables (global) too.
    let created = std::mem::take(&mut interp.context_mut().draw.created);
    if !created.is_empty() {
        // (At the top level, where scripts mostly are, asked of the scope itself: listing
        // every variable for each command made is slow with thousands of them.)
        let globals = if interp.scope_level() == 0 { Vec::new() } else { interp.vars_in_global_scope() };
        for n in created {
            let exists = if interp.scope_level() == 0 { interp.var_exists(&Value::from(n.as_str())) } else { globals.iter().any(|g| g.as_str() == n) };
            if !exists {
                let set = Value::from(vec![Value::from("set"), Value::from(n.as_str()), Value::from(n.as_str())]);
                let _ = interp.eval_at_scope(0, &set);
            }
        }
    }
    ran
}

/// Runs a DRAW command (of the dispatcher's) on the session, logging what it reports.
fn run_draw_command(interp: &mut Interp<Session>, name: &str, args: &[&str]) -> MoltResult {
    let session = interp.context_mut();
    // (A panic is an exception, as a C++ one in DRAW: an error a script may catch.)
    let draw = &mut session.draw;
    let ran = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| super::dispatch::dispatch(draw, name, args)));
    let ran = match ran {
        Ok(r) => r,
        Err(p) => {
            let what = p.downcast_ref::<String>().cloned().or_else(|| p.downcast_ref::<&str>().map(|s| s.to_string())).unwrap_or_default();
            return molt_err!(format!("An exception was caught: {what}"));
        }
    };
    // (`checkshape` reports a faulty shape, as DRAW does, rather than failing: the script
    // may look at what it says, and the log's rules tell it a failure.)
    let ran = match ran {
        Some(Err(DrawError::Failed(m))) if name == "checkshape" && m.starts_with("Faulty shape ") => Some(Ok(faulty_report(args, &m))),
        other => other,
    };
    match ran {
        Some(Ok(out)) => {
            // (DRAW echoes what a command reports.)
            if !out.is_empty() {
                session.log.extend(out.lines().map(str::to_string));
            }
            Ok(Value::from(out))
        }
        Some(Err(e)) => {
            // (A missing data file skips the script, as OCCT's harness does.)
            let message = e.to_string();
            if message.contains("data file") && message.contains("not found") {
                session.missing.get_or_insert(message.clone());
            }
            molt_err!(match e {
                DrawError::NotImplemented { .. } => message,
                DrawError::Failed(m) => format!("Error: {m}"),
            })
        }
        None => molt_err!(format!("invalid command name \"{name}\"")),
    }
}

/// What DRAW's `checkshape` says of a faulty shape, from the problems in `message` (the
/// dispatcher's failure): each in a faulty shape of its own, or, given a prefix for their
/// names, how many there are (without the word OCCT's log rules fail on).
fn faulty_report(args: &[&str], message: &str) -> String {
    let problems: Vec<&str> = message.split_once(": ").map_or(message, |(_, p)| p).split("; ").collect();
    let names: Vec<&&str> = args.iter().filter(|a| !a.starts_with('-')).collect();
    match names[..] {
        // (Counted, as DRAW does, the same each time the shape is checked.)
        [shape, prefix, ..] => format!(" -- The Shape {shape} has problems :\n*** Shapes with problems : {} (named {prefix}_*)\n", problems.len()),
        _ => {
            let lines: String = problems.iter().enumerate().map(|(i, p)| format!("On Shape faulty_{} :\n{p}\n", i + 1)).collect();
            format!("{lines}Faulty shapes in variables faulty_1 to faulty_{} \n", problems.len())
        }
    }
}

fn command_type(name: &str, interp: &Interp<Session>) -> Option<CommandKind> {
    if molt_forked::__private::is_builtin(name, interp.standard_library()) {
        return Some(CommandKind::Native);
    }
    None
}

/// The Tcl commands Molt lacks, and those it has that a DRAW session runs otherwise.
fn extra(name: &str, interp: &mut Interp<Session>, argv: &[Value]) -> Option<MoltResult> {
    let args: Vec<String> = argv[1..].iter().map(|v| v.as_str().to_string()).collect();
    Some(match name {
        "unset" => {
            // (Unsetting a DRAW variable deletes it, as in DRAW; a name that is not a Tcl
            // variable is no error then.)
            let quiet = args.iter().any(|a| a == "-nocomplain");
            let mut words = vec![argv[0].clone()];
            if !quiet {
                words.push(Value::from("-nocomplain"));
            }
            for a in &args {
                if !a.starts_with('-') && interp.context_mut().draw.has(a) {
                    let _ = interp.context_mut().draw.take_var(a);
                } else if !quiet && !a.starts_with('-') && interp.var(&Value::from(a.as_str())).is_err() {
                    return Some(molt_err!(format!("can't unset \"{a}\": no such variable")));
                }
                words.push(Value::from(a.as_str()));
            }
            return molt_forked::__private::execute_builtin("unset", interp, &words);
        }
        "puts" => {
            // `puts ?-nonewline? ?channel? string`: a TODO declares a known failure.
            let text = args.last().cloned().unwrap_or_default();
            let session = interp.context_mut();
            if let Some(rest) = text.trim_start().strip_prefix("TODO ") {
                let platforms = rest.split(':').next().unwrap_or("").split_whitespace().skip(1).collect::<Vec<_>>();
                if platforms.iter().any(|p| ["ALL", "All", "Linux", "Debian", "Ubuntu"].contains(&p.trim_matches(','))) {
                    session.todos.push(text.clone());
                }
                return Some(molt_ok!());
            }
            session.log.extend(text.lines().map(str::to_string));
            molt_ok!()
        }
        "info" => {
            // (Molt's, with Tcl's unique abbreviations of its subcommands, and `info script`.)
            let Some(sub) = args.first() else { return None };
            let full = ["exists", "script", "commands", "procs", "vars", "globals", "locals", "args", "body", "default", "level", "complete", "hostname", "patchlevel", "tclversion"]
                .into_iter()
                .filter(|f| f.starts_with(sub.as_str()))
                .collect::<Vec<_>>();
            match full[..] {
                ["script"] => molt_ok!(),
                ["hostname"] => molt_ok!("localhost"),
                ["patchlevel"] => molt_ok!("8.6.13"),
                ["tclversion"] => molt_ok!("8.6"),
                [one] if one != sub.as_str() => {
                    let mut words: Vec<Value> = argv.to_vec();
                    words[1] = Value::from(one);
                    return molt_forked::__private::execute_builtin("info", interp, &words);
                }
                _ => return None,
            }
        }
        "whatis" => {
            // DRAW's `whatis name ...`: what each variable is.
            let session = interp.context_mut();
            let mut out = Vec::new();
            for n in &args {
                use super::DrawObject as O;
                let what = match session.draw.get(n) {
                    Some(O::Shape(s)) => {
                        let s = *s;
                        let orientation = format!("{:?}", s.orientation).to_uppercase();
                        let closed = if shape_closed(&session.draw, s) { " Closed" } else { "" };
                        format!("is a shape {} {orientation} Orientable{closed}", format!("{:?}", s.kind()).to_uppercase())
                    }
                    Some(O::Curve(_)) => "is a 3d curve".into(),
                    Some(O::Curve2(_)) => "is a 2d curve".into(),
                    Some(O::Surface(_)) => "is a surface".into(),
                    Some(O::Point(_)) => "is a point".into(),
                    Some(O::Point2(_)) => "is a 2d point".into(),
                    Some(O::Number(_)) => "is a numeric variable".into(),
                    Some(_) => "is a variable".into(),
                    None => "is not a DRAW variable".into(),
                };
                out.push(format!("{n} {what}"));
            }
            let text = out.join("\n");
            session.log.extend(out);
            molt_ok!(text)
        }
        "dump" => {
            // DRAW's `dump name ...`: points and curves as OCCT prints them (in part).
            let session = interp.context_mut();
            let mut text = String::new();
            for n in &args {
                use super::DrawObject as O;
                let Some(object) = session.draw.get(n) else { continue };
                text += &format!("\n\n*********** Dump of {n} *************\n");
                text += &match object {
                    O::Point(p) => format!("Point : {:.15e}, {:.15e}, {:.15e}\n", p.x, p.y, p.z),
                    O::Point2(p) => format!("Point 2d : {:.15e}, {:.15e}\n", p.x, p.y),
                    O::Curve(c) => describe_curve(c),
                    O::Surface(kentron::geometry::surface::Surface::Plane(p)) => {
                        let f = p.frame();
                        let xyz = |v: [f64; 3]| format!("{}, {}, {} ", v[0], v[1], v[2]);
                        let d = |d: kentron::math::Direction| xyz([d.x(), d.y(), d.z()]);
                        let o = f.origin();
                        format!("Plane\n  Origin :{}\n  Axis   :{}\n  XAxis  :{}\n  YAxis  :{}\n\n", xyz([o.x, o.y, o.z]), d(f.z()), d(f.x()), d(f.y()))
                    }
                    // (A number without a newline after it, as OCCT prints it.)
                    O::Number(x) => format!("{x}"),
                    O::Curve2(kentron::geometry::curve2::Curve2::BSpline(b)) => {
                        // (OCCT's layout, with its missing newline before the poles.)
                        let mut t = format!("BSplineCurve\n  Degree {}, {} Poles, {}  KnotsPoles :\n", b.degree(), b.poles().len(), b.knots().len());
                        for (i, p) in b.poles().iter().enumerate() {
                            t += &format!("   {} : {}, {} \n", i + 1, p.x, p.y);
                        }
                        t += "Knots :\n";
                        for (i, (k, m)) in b.knots().iter().zip(b.multiplicities()).enumerate() {
                            t += &format!("   {} : {} {}\n", i + 1, k, m);
                        }
                        t
                    }
                    other => format!("{other:?}\n"),
                };
            }
            molt_ok!(text)
        }
        "boptions" => {
            let text = "Fuzzy value : 0\nRunParallel : 0\nNonDestructive : 0\nGlue : 0\nUseOBB : 0\n";
            molt_ok!(text)
        }
        "directory" => {
            // DRAW's `directory ?pattern?`: the variables' names matching the pattern.
            let pattern = args.first().cloned().unwrap_or_else(|| "*".into());
            let mut names: Vec<String> = interp.context_mut().draw.variable_names().into_iter().filter(|n| glob_match(&pattern, n)).collect();
            names.sort();
            molt_ok!(Value::from(names.into_iter().map(Value::from).collect::<Vec<_>>()))
        }
        "exit" => {
            interp.context_mut().exited = true;
            Err(Exception::molt_return(Value::empty()))
        }
        "locate_data_file" => match args.first().and_then(|n| find_data_file(n).map(|p| (n, p))) {
            Some((_, path)) => molt_ok!(path),
            None => {
                let n = args.first().cloned().unwrap_or_default();
                interp.context_mut().missing.get_or_insert(format!("data file {n} not found"));
                molt_err!(format!("Tcl Exception: File {n} could not be found"))
            }
        },
        "repeat" => {
            // DRAW's `repeat count body`.
            let [count, body] = &args[..] else { return Some(molt_err!("wrong # args: should be \"repeat count body\"")) };
            let Ok(n) = count.trim().parse::<i64>() else { return Some(molt_err!(format!("repeat: {count} is not a count"))) };
            let body = Value::from(body.as_str());
            for _ in 0..n {
                match interp.eval_value(&body) {
                    Ok(_) => {}
                    Err(e) if e.code() == ResultCode::Break => break,
                    Err(e) if e.code() == ResultCode::Continue => continue,
                    Err(e) => return Some(Err(e)),
                }
            }
            molt_ok!()
        }
        "isdraw" => {
            let known = args.first().is_some_and(|n| interp.context_mut().draw.get(n).is_some());
            molt_ok!(Value::from(known))
        }
        "source" => {
            let Some(file) = args.last() else { return Some(molt_err!("wrong # args: should be \"source fileName\"")) };
            let session = interp.context_mut();
            let candidates = std::iter::once(PathBuf::from(file)).chain(session.source_dirs.iter().map(|d| d.join(file)));
            let Some(text) = candidates.into_iter().find_map(|p| std::fs::read_to_string(p).ok()) else {
                return Some(molt_err!(format!("couldn't read file \"{file}\"")));
            };
            interp.eval(&text)
        }
        "regexp" => regexp(interp, &args),
        "regsub" => regsub(interp, &args),
        "lsearch" => lsearch(&args),
        "lset" => lset(interp, &args),
        "lsort" => lsort(&args),
        "format" => format(&args),
        "file" => file(&args),
        "checkplatform" => molt_ok!(Value::from(if args.iter().any(|a| a == "-linux") { 1 } else { 0 })),
        _ => return None,
    })
}

/// A Tcl regular expression in Rust's syntax (`\y` a word boundary, `\m` and `\M` its ends).
fn rust_regex(pattern: &str, nocase: bool, line: bool) -> Result<regex::Regex, Exception> {
    let p = pattern.replace("\\y", "\\b").replace("\\m", "\\b").replace("\\M", "\\b");
    regex::RegexBuilder::new(&p)
        .case_insensitive(nocase)
        .multi_line(line)
        .dot_matches_new_line(!line)
        .build()
        .map_err(|e| Exception::molt_err(Value::from(format!("couldn't compile regular expression pattern: {e}"))))
}

/// Leading switches (up to `--`), and the rest.
fn switches(args: &[String]) -> (Vec<&str>, &[String]) {
    let mut k = 0;
    let mut out = Vec::new();
    while k < args.len() && args[k].starts_with('-') && args.len() - k > 1 {
        if args[k] == "--" {
            k += 1;
            break;
        }
        out.push(args[k].as_str());
        k += 1;
    }
    (out, &args[k..])
}

/// `regexp ?-nocase? ?-all? ?-inline? ?-line? exp string ?matchVar? ?subMatchVar ...?`
fn regexp(interp: &mut Interp<Session>, args: &[String]) -> MoltResult {
    let (sw, rest) = switches(args);
    let [exp, string, vars @ ..] = rest else { return molt_err!("wrong # args: should be \"regexp ?switches? exp string ?matchVar? ?subMatchVar ...?\"") };
    let re = rust_regex(exp, sw.contains(&"-nocase"), sw.contains(&"-line"))?;
    let groups = |c: &regex::Captures| -> Vec<Value> { (0..c.len()).map(|i| Value::from(c.get(i).map_or("", |m| m.as_str()))).collect() };
    if sw.contains(&"-inline") {
        let mut out = Vec::new();
        for c in re.captures_iter(string) {
            out.extend(groups(&c));
            if !sw.contains(&"-all") {
                break;
            }
        }
        return molt_ok!(Value::from(out));
    }
    let all: Vec<regex::Captures> = if sw.contains(&"-all") { re.captures_iter(string).collect() } else { re.captures(string).into_iter().collect() };
    // (No match leaves the variables as they were, as in Tcl.)
    if let Some(last) = all.last() {
        let g = groups(last);
        for (k, var) in vars.iter().enumerate() {
            interp.set_var(&Value::from(var.as_str()), g.get(k).cloned().unwrap_or_else(Value::empty))?;
        }
    }
    molt_ok!(Value::from(all.len() as MoltInt))
}

/// `regsub ?-all? ?-nocase? exp string subSpec ?varName?`
fn regsub(interp: &mut Interp<Session>, args: &[String]) -> MoltResult {
    let (sw, rest) = switches(args);
    let [exp, string, spec, var @ ..] = rest else { return molt_err!("wrong # args: should be \"regsub ?switches? exp string subSpec ?varName?\"") };
    let re = rust_regex(exp, sw.contains(&"-nocase"), sw.contains(&"-line"))?;
    // `&` and `\0` the match, `\1`... its groups.
    let mut replacement = String::new();
    let mut chars = spec.chars().peekable();
    while let Some(c) = chars.next() {
        match c {
            '&' => replacement.push_str("${0}"),
            '\\' => match chars.next() {
                Some(d) if d.is_ascii_digit() => replacement.push_str(&format!("${{{d}}}")),
                Some(o) => replacement.push(o),
                None => replacement.push('\\'),
            },
            '$' => replacement.push_str("$$"),
            c => replacement.push(c),
        }
    }
    let count = if sw.contains(&"-all") { re.find_iter(string).count() } else { usize::from(re.is_match(string)) };
    let result = if sw.contains(&"-all") { re.replace_all(string, replacement.as_str()) } else { re.replace(string, replacement.as_str()) };
    match var.first() {
        Some(v) => {
            interp.set_var(&Value::from(v.as_str()), Value::from(result.to_string()))?;
            molt_ok!(Value::from(count as MoltInt))
        }
        None => molt_ok!(Value::from(result.to_string())),
    }
}

/// Tcl's `string match`: `*`, `?`, `[chars]` and `\x`.
pub fn glob_match(pattern: &str, text: &str) -> bool {
    fn go(p: &[char], t: &[char]) -> bool {
        match p.first() {
            None => t.is_empty(),
            Some('*') => (0..=t.len()).any(|k| go(&p[1..], &t[k..])),
            Some('?') => !t.is_empty() && go(&p[1..], &t[1..]),
            Some('[') => {
                let Some(end) = p.iter().position(|&c| c == ']') else { return !t.is_empty() && t[0] == '[' && go(&p[1..], &t[1..]) };
                let Some(&c) = t.first() else { return false };
                let set = &p[1..end];
                let mut ok = false;
                let mut k = 0;
                while k < set.len() {
                    if k + 2 < set.len() && set[k + 1] == '-' {
                        ok |= set[k] <= c && c <= set[k + 2];
                        k += 3;
                    } else {
                        ok |= set[k] == c;
                        k += 1;
                    }
                }
                ok && go(&p[end + 1..], &t[1..])
            }
            Some('\\') if p.len() > 1 => !t.is_empty() && t[0] == p[1] && go(&p[2..], &t[1..]),
            Some(&c) => !t.is_empty() && t[0] == c && go(&p[1..], &t[1..]),
        }
    }
    let (p, t): (Vec<char>, Vec<char>) = (pattern.chars().collect(), text.chars().collect());
    go(&p, &t)
}

/// `lsearch ?-exact|-glob|-regexp? ?-all? ?-inline? ?-not? ?-nocase? list pattern`
fn lsearch(args: &[String]) -> MoltResult {
    let (sw, rest) = switches(args);
    let [list, pattern] = rest else { return molt_err!("wrong # args: should be \"lsearch ?options? list pattern\"") };
    let items = Value::from(list.as_str()).as_list()?;
    let nocase = sw.contains(&"-nocase");
    let re = if sw.contains(&"-regexp") { Some(rust_regex(pattern, nocase, false)?) } else { None };
    let fold = |s: &str| if nocase { s.to_lowercase() } else { s.to_string() };
    let matches = |s: &str| -> bool {
        let m = if let Some(re) = &re {
            re.is_match(s)
        } else if sw.contains(&"-exact") {
            fold(s) == fold(pattern)
        } else {
            glob_match(&fold(pattern), &fold(s))
        };
        m != sw.contains(&"-not")
    };
    let found: Vec<usize> = items.iter().enumerate().filter(|(_, v)| matches(v.as_str())).map(|(k, _)| k).collect();
    let inline = sw.contains(&"-inline");
    if sw.contains(&"-all") {
        let out: Vec<Value> = found.iter().map(|&k| if inline { items[k].clone() } else { Value::from(k as MoltInt) }).collect();
        return molt_ok!(Value::from(out));
    }
    match found.first() {
        Some(&k) => molt_ok!(if inline { items[k].clone() } else { Value::from(k as MoltInt) }),
        None => molt_ok!(if inline { Value::empty() } else { Value::from(-1 as MoltInt) }),
    }
}

/// `lset varName ?index ...? value`: the variable's list with the element at the indices
/// (one in each nested list; `end`, `end-N`) replaced by `value`.
fn lset(interp: &mut Interp<Session>, args: &[String]) -> MoltResult {
    let [name, rest @ ..] = args else { return molt_err!("wrong # args: should be \"lset listVar ?index? ?index ...? value\"") };
    let Some((value, indices)) = rest.split_last() else { return molt_err!("wrong # args: should be \"lset listVar ?index? ?index ...? value\"") };
    // (Several indices given as one list, as Tcl allows.)
    let indices: Vec<String> = match indices {
        [one] => Value::from(one.as_str()).as_list()?.iter().map(|v| v.as_str().to_string()).collect(),
        many => many.to_vec(),
    };
    let var = Value::from(name.as_str());
    let list = interp.var(&var)?;
    fn set_in(list: &Value, indices: &[String], value: &str) -> MoltResult {
        let Some((first, rest)) = indices.split_first() else { return molt_ok!(Value::from(value)) };
        let mut items: Vec<Value> = list.as_list()?.to_vec();
        let n = items.len() as i64;
        let k = match first.trim() {
            "end" => n - 1,
            t if t.starts_with("end-") => n - 1 - t[4..].parse::<i64>().map_err(|_| Exception::molt_err(Value::from(format!("bad index \"{t}\""))))?,
            t => t.parse::<i64>().map_err(|_| Exception::molt_err(Value::from(format!("bad index \"{t}\""))))?,
        };
        if k < 0 || k > n || (k == n && !rest.is_empty()) {
            return molt_err!("list index out of range");
        }
        let inner = if k < n { set_in(&items[k as usize], rest, value)? } else { Value::from(value) };
        if k == n {
            items.push(inner);
        } else {
            items[k as usize] = inner;
        }
        molt_ok!(Value::from(items))
    }
    let new = set_in(&list, &indices, value)?;
    interp.set_var(&var, new.clone())?;
    molt_ok!(new)
}

/// `lsort ?-ascii|-dictionary|-integer|-real? ?-increasing|-decreasing? ?-unique? list`
fn lsort(args: &[String]) -> MoltResult {
    let (sw, rest) = switches(args);
    let [list] = rest else { return molt_err!("wrong # args: should be \"lsort ?options? list\"") };
    let mut items: Vec<Value> = Value::from(list.as_str()).as_list()?.to_vec();
    let numeric = sw.contains(&"-integer") || sw.contains(&"-real");
    let mut failed = None;
    items.sort_by(|a, b| {
        if numeric {
            match (a.as_float(), b.as_float()) {
                (Ok(x), Ok(y)) => x.total_cmp(&y),
                _ => {
                    failed = Some(format!("expected number but got \"{}\"", a.as_str()));
                    std::cmp::Ordering::Equal
                }
            }
        } else if sw.contains(&"-dictionary") {
            a.as_str().to_lowercase().cmp(&b.as_str().to_lowercase())
        } else {
            a.as_str().cmp(b.as_str())
        }
    });
    if let Some(f) = failed {
        return molt_err!(f);
    }
    if sw.contains(&"-decreasing") {
        items.reverse();
    }
    if sw.contains(&"-unique") {
        items.dedup_by(|a, b| a.as_str() == b.as_str());
    }
    molt_ok!(Value::from(items))
}

/// `format formatString ?arg ...?`: `%d %i %u %x %o %c %s %f %e %g` (and their capitals),
/// with flags, width and precision.
fn format(args: &[String]) -> MoltResult {
    let Some((spec, rest)) = args.split_first() else { return molt_err!("wrong # args: should be \"format formatString ?arg ...?\"") };
    let mut out = String::new();
    let mut next = rest.iter();
    let chars: Vec<char> = spec.chars().collect();
    let mut k = 0;
    while k < chars.len() {
        if chars[k] != '%' {
            out.push(chars[k]);
            k += 1;
            continue;
        }
        k += 1;
        if chars.get(k) == Some(&'%') {
            out.push('%');
            k += 1;
            continue;
        }
        let mut flags = String::new();
        while k < chars.len() && "-+ 0#".contains(chars[k]) {
            flags.push(chars[k]);
            k += 1;
        }
        let mut width = String::new();
        while k < chars.len() && chars[k].is_ascii_digit() {
            width.push(chars[k]);
            k += 1;
        }
        let mut precision: Option<usize> = None;
        if chars.get(k) == Some(&'.') {
            k += 1;
            let mut p = String::new();
            while k < chars.len() && chars[k].is_ascii_digit() {
                p.push(chars[k]);
                k += 1;
            }
            precision = Some(p.parse().unwrap_or(0));
        }
        while k < chars.len() && "hlLqjzt".contains(chars[k]) {
            k += 1;
        }
        let Some(&conv) = chars.get(k) else { return molt_err!("format string ended in middle of field specifier") };
        k += 1;
        let Some(arg) = next.next() else { return molt_err!("not enough arguments for all format specifiers") };
        let value = Value::from(arg.as_str());
        let body = match conv {
            'd' | 'i' | 'u' => {
                let n = value.as_float().map_err(|_| Exception::molt_err(Value::from(format!("expected integer but got \"{arg}\""))))?;
                let n = n.trunc() as i64;
                if flags.contains('+') && n >= 0 { format!("+{n}") } else { n.to_string() }
            }
            'x' | 'X' | 'o' => {
                let n = value.as_int()?;
                match conv {
                    'x' => format!("{n:x}"),
                    'X' => format!("{n:X}"),
                    _ => format!("{n:o}"),
                }
            }
            'c' => char::from_u32(value.as_int()? as u32).map(String::from).unwrap_or_default(),
            's' => match precision {
                Some(p) => arg.chars().take(p).collect(),
                None => arg.clone(),
            },
            'f' | 'F' | 'e' | 'E' | 'g' | 'G' => {
                let x = value.as_float()?;
                let p = precision.unwrap_or(6);
                let s = match conv {
                    'f' | 'F' => format!("{x:.p$}"),
                    'e' | 'E' => c_exponent(x, p, conv == 'E'),
                    _ => c_general(x, if p == 0 { 1 } else { p }, flags.contains('#'), conv == 'G'),
                };
                if flags.contains('+') && x >= 0.0 { format!("+{s}") } else { s }
            }
            other => return molt_err!(format!("bad field specifier \"{other}\"")),
        };
        let width: usize = width.parse().unwrap_or(0);
        let padded = if body.chars().count() >= width {
            body
        } else if flags.contains('-') {
            format!("{body:<width$}")
        } else if flags.contains('0') && !matches!(conv, 's' | 'c') {
            let (sign, digits) = if body.starts_with(['-', '+']) { body.split_at(1) } else { ("", body.as_str()) };
            format!("{sign}{digits:0>w$}", w = width - sign.len())
        } else {
            format!("{body:>width$}")
        };
        out.push_str(&padded);
    }
    molt_ok!(Value::from(out))
}

/// C's `%e`: a mantissa with `precision` decimals and a signed exponent of two digits or more.
fn c_exponent(x: f64, precision: usize, upper: bool) -> String {
    let s = format!("{x:.precision$e}");
    let (mantissa, exponent) = s.split_once('e').unwrap_or((&s, "0"));
    let e: i32 = exponent.parse().unwrap_or(0);
    let sign = if e < 0 { '-' } else { '+' };
    let out = format!("{mantissa}e{sign}{:02}", e.abs());
    if upper { out.to_uppercase() } else { out }
}

/// C's `%g`: `precision` significant digits, as `%e` or `%f` by the exponent, without
/// trailing zeros unless `alternate`.
fn c_general(x: f64, precision: usize, alternate: bool, upper: bool) -> String {
    if x == 0.0 {
        return "0".into();
    }
    let e = x.abs().log10().floor() as i32;
    let trim = |s: String| -> String {
        if alternate || !s.contains('.') {
            return s;
        }
        let (m, rest) = match s.find(['e', 'E']) {
            Some(k) => (s[..k].to_string(), s[k..].to_string()),
            None => (s.clone(), String::new()),
        };
        let m = m.trim_end_matches('0').trim_end_matches('.').to_string();
        m + &rest
    };
    if e < -4 || e >= precision as i32 {
        trim(c_exponent(x, precision - 1, upper))
    } else {
        let decimals = (precision as i32 - 1 - e).max(0) as usize;
        trim(format!("{x:.decimals$}"))
    }
}

/// `file exists|isfile|isdirectory|delete|join|dirname|tail|rootname|extension|size|mkdir|normalize ...`
fn file(args: &[String]) -> MoltResult {
    let Some((sub, rest)) = args.split_first() else { return molt_err!("wrong # args: should be \"file subcommand ?arg ...?\"") };
    let path = || rest.last().map(PathBuf::from).unwrap_or_default();
    match sub.as_str() {
        "exists" => molt_ok!(Value::from(path().exists())),
        "isfile" => molt_ok!(Value::from(path().is_file())),
        "isdirectory" => molt_ok!(Value::from(path().is_dir())),
        "delete" => {
            for p in rest.iter().filter(|a| !a.starts_with('-')) {
                let _ = std::fs::remove_file(p);
            }
            molt_ok!()
        }
        "mkdir" => {
            for p in rest {
                let _ = std::fs::create_dir_all(p);
            }
            molt_ok!()
        }
        "join" => molt_ok!(Value::from(rest.iter().collect::<PathBuf>().to_string_lossy().to_string())),
        "dirname" => molt_ok!(Value::from(path().parent().map(|p| p.to_string_lossy().to_string()).unwrap_or_else(|| ".".into()))),
        "tail" => molt_ok!(Value::from(path().file_name().map(|p| p.to_string_lossy().to_string()).unwrap_or_default())),
        "rootname" => {
            let p = path();
            molt_ok!(Value::from(p.with_extension("").to_string_lossy().to_string()))
        }
        "extension" => molt_ok!(Value::from(path().extension().map(|e| format!(".{}", e.to_string_lossy())).unwrap_or_default())),
        "size" => match std::fs::metadata(path()) {
            Ok(m) => molt_ok!(Value::from(m.len() as MoltInt)),
            Err(e) => molt_err!(format!("could not read \"{}\": {e}", path().display())),
        },
        "normalize" => molt_ok!(Value::from(std::fs::canonicalize(path()).unwrap_or_else(|_| path()).to_string_lossy().to_string())),
        other => molt_err!(format!("file {other} is not supported")),
    }
}

/// The rules telling a failure in the log: OCCT's `parse.rules` (`FAILED` and `IGNORE`, in
/// order), the group's and grid's before the global ones.
fn rules(tests: &Path, group: &str, grid: &str) -> Vec<(bool, regex::Regex)> {
    let mut out = Vec::new();
    for file in [tests.join(group).join(grid).join("parse.rules"), tests.join(group).join("parse.rules"), tests.join("parse.rules")] {
        let Ok(text) = std::fs::read_to_string(file) else { continue };
        for line in text.lines() {
            let line = line.trim();
            let Some((status, rest)) = line.split_once(char::is_whitespace) else { continue };
            let failed = match status {
                "FAILED" => true,
                "IGNORE" | "SKIPPED" | "OK" | "BAD" => false,
                _ => continue,
            };
            let rest = rest.trim();
            let Some(body) = rest.strip_prefix('/') else { continue };
            let Some(end) = body.rfind('/') else { continue };
            if let Ok(re) = rust_regex(&body[..end], false, false) {
                out.push((failed, re));
            }
        }
    }
    // (A faulty shape fails a script in every group, not only those whose rules say so:
    // kentron's results are to be valid wherever they are made.)
    out.push((true, regex::Regex::new(r"\bFaulty\b").expect("a pattern")));
    out
}

/// The OCCT TODO a script declares that applies on this platform (`puts "TODO id platforms:
/// message"`): a known failure of OCCT's.
pub fn known_failure(script: &str) -> Option<String> {
    let todo = regex::Regex::new(r#"^\s*puts\s+"(TODO\s+\S+\s+([^:"]+):[^"]*)""#).expect("a pattern");
    // (Only a TODO outside braces holds whatever runs: one in a block, as under an `if` on
    // the grid's variables, is told when the script puts it.)
    let mut depth = 0i32;
    for line in script.lines() {
        if depth == 0 {
            if let Some(c) = todo.captures(line) {
                if c[2].split([' ', ',']).any(|p| ["ALL", "All", "Linux", "Debian", "Ubuntu"].contains(&p)) {
                    return Some(c[1].to_string());
                }
            }
        }
        depth += line.matches('{').count() as i32 - line.matches('}').count() as i32;
    }
    None
}

/// The file of the script `<group>/<grid>/<case>` under `tests`: in the grid's directory, or
/// in the one its `cases.list` names (a grid running another's scripts with its own begin
/// and end files).
fn case_file(tests: &Path, group: &str, grid: &str, case: &str) -> std::path::PathBuf {
    let dir = tests.join(group).join(grid);
    let file = dir.join(case);
    if file.is_file() {
        return file;
    }
    match std::fs::read_to_string(dir.join("cases.list")).ok().and_then(|t| t.lines().map(str::trim).find(|l| !l.is_empty()).map(String::from)) {
        Some(alias) => dir.join(alias).join(case),
        None => file,
    }
}

/// Runs OCCT's script `tests/<group>/<grid>/<case>` from the checkout `tests` (its `tests`
/// directory) as OCCT's harness runs it.
pub fn run_script(tests: &Path, group: &str, grid: &str, case: &str) -> Outcome {
    // (A TODO that applies here makes it a known failure of OCCT's: not run, as the
    // converted scripts with one are not.)
    let script = case_file(tests, group, grid, case);
    if let Some(todo) = std::fs::read_to_string(&script).ok().and_then(|t| known_failure(&t)) {
        return Outcome::KnownFailure(todo);
    }
    let session = Session {
        draw: Draw::new(),
        log: Vec::new(),
        source_dirs: vec![tests.join(group).join(grid), tests.join(group), tests.to_path_buf()],
        missing: None,
        todos: Vec::new(),
        exited: false,
    };
    let set = CommandSet::new(execute, command_type, &[], &[]);
    let mut interp = InterpBuilder::new(session, set).standard_library(StandardLibrary::Full).environment(true).build();
    let temp = std::env::temp_dir();
    for (var, value) in [
        ("groupname", group.to_string()),
        ("gridname", grid.to_string()),
        ("casename", case.to_string()),
        ("dirname", tests.to_string_lossy().to_string()),
        ("imagedir", temp.to_string_lossy().to_string()),
        ("test_image", case.to_string()),
        ("logdir", temp.to_string_lossy().to_string()),
    ] {
        let _ = interp.set_scalar(var, Value::from(value));
    }
    let files = [
        tests.join(group).join("begin"),
        tests.join(group).join(grid).join("begin"),
        script,
        tests.join(group).join(grid).join("end"),
        tests.join(group).join("end"),
    ];
    let mut errors = Vec::new();
    for f in files {
        let Ok(text) = std::fs::read_to_string(&f) else { continue };
        if let Err(e) = interp.eval(&text) {
            if e.code() != ResultCode::Return {
                errors.push(e.value().as_str().to_string());
            }
        }
        let session = interp.context_mut();
        if let Some(m) = &session.missing {
            return Outcome::Skipped(m.clone());
        }
        if session.exited || !errors.is_empty() {
            break;
        }
    }
    let rules = rules(tests, group, grid);
    let session = interp.context_mut();
    // (For debugging: the whole log.)
    if std::env::var_os("KENTRON_DRAW_LOG").is_some() {
        for line in &session.log {
            eprintln!("{line}");
        }
        for e in &errors {
            eprintln!("ERROR {e}");
        }
    }
    if let Some(todo) = session.todos.first() {
        return Outcome::KnownFailure(todo.clone());
    }
    let mut failures = errors;
    for line in session.log.iter() {
        if let Some((failed, _)) = rules.iter().find(|(_, re)| re.is_match(line)) {
            if *failed {
                failures.push(line.clone());
            }
        }
    }
    if failures.is_empty() { Outcome::Passed } else { Outcome::Failed(failures) }
}

/// A curve as DRAW's `dump` begins describing it.
fn describe_curve(c: &kentron::geometry::curve::Curve) -> String {
    use kentron::geometry::curve::Curve;
    match c {
        Curve::Trimmed(t) => {
            let (a, b) = t.range();
            format!("Trimmed curve\nParameters : {a} {b}\nBasis curve :\n{}", describe_curve(t.basis()))
        }
        Curve::Line(_) => "Line\n".into(),
        Curve::Circle(_) => "Circle\n".into(),
        Curve::Ellipse(_) => "Ellipse\n".into(),
        Curve::Hyperbola(_) => "Hyperbola\n".into(),
        Curve::Parabola(_) => "Parabola\n".into(),
        Curve::Bezier(_) => "BezierCurve\n".into(),
        Curve::BSpline(b) => {
            // (OCCT's layout: degree, poles (with weights if rational), knots.)
            let rational = b.is_rational();
            let mut t = format!("BSplineCurve{}\n  Degree {}, {} Poles, {}  Knots\nPoles :\n", if rational { " rational" } else { "" }, b.degree(), b.poles().len(), b.knots().len());
            for (i, p) in b.poles().iter().enumerate() {
                t += &format!("   {} : {}, {}, {} ", i + 1, p.x, p.y, p.z);
                if rational {
                    t += &format!(" {}", b.weight(i));
                }
                t += "\n";
            }
            t += "Knots :\n";
            for (i, (k, m)) in b.knots().iter().zip(b.multiplicities()).enumerate() {
                t += &format!("   {} :  {} {}\n", i + 1, k, m);
            }
            t
        }
        Curve::Offset(_) => "OffsetCurve\n".into(),
    }
}

/// Whether a shape is closed, as OCCT flags it: an edge from a vertex back to it, a wire
/// round a loop, a shell enclosing a volume; a solid only as flagged by the command making
/// it (OCCT's builders differ: a loft's is, a prism's is not).
fn shape_closed(draw: &super::Draw, s: kentron::topology::Shape) -> bool {
    use kentron::topology::{NodeKind, ShellId, WireId};
    let model = &draw.model;
    match s.kind() {
        NodeKind::Edge => model.is_closed_edge(s),
        NodeKind::Wire => WireId::try_from(s.node).is_ok_and(|w| model.is_closed_wire(w)),
        NodeKind::Shell => ShellId::try_from(s.node).is_ok_and(|sh| model.is_closed_shell(sh)),
        NodeKind::Solid => draw.closed_solids.contains(&s.node),
        _ => false,
    }
}
