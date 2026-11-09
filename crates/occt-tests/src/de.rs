//! OCCT's data-exchange tests (`tests/de`, copied to `de/` by `tools/de_import.py`): a STEP
//! or IGES file read, written and read again, statistics of both readings compared with
//! OCCT's reference data, as the grids' `end` script compares them.
//!
//! Each case names its file and holds the reference (`ref_data`): one line for each
//! statistic, the second reading's value first and the first reading's in parentheses. The
//! lines are made here as the script makes them, from kentron's readings, and compared by
//! its rules: CHECKSHAPE (faulty wires, faces, shells and solids) no more than OCCT's;
//! TOLERANCE (largest and mean) no more than OCCT's by 1% (and 1e-4); every other line equal.
//! DATA and TPSTAT, the numbers of warnings and failures OCCT's translators report, are
//! OCCT's own and not compared.
//!
//! The document statistics (LABELS, PROPS, COLORS, LAYERS) are those of OCCT's document of
//! the file (its XCAF labels, as XStat counts them), worked out from kentron's assembly
//! tree: a label for each product, each of an assembly's components and each face or edge
//! of a part with a name, colour or layer of its own.

pub mod colors;

use std::path::Path;

use kentron::io::step::Node;
use kentron::topology::{Model, NodeKind, Shape};
use kentron::math::Transform;

use crate::draw::{Draw, DrawObject};

/// A case: the file it reads and OCCT's reference lines.
#[derive(Debug, Clone)]
pub struct Case {
    pub file: String,
    pub reference: Vec<String>,
    /// `set ProductMode OFF`: the file read without its product structure.
    pub product_mode_off: bool,
}

/// The case at `path` (a copy of an OCCT `tests/de` case script).
pub fn load(path: &Path) -> Result<Case, String> {
    let text = std::fs::read_to_string(path).map_err(|e| format!("{}: {e}", path.display()))?;
    let file = regex::Regex::new(r"set filename\s+(\S+)").expect("a pattern").captures(&text).map(|c| c[1].to_string()).ok_or("no filename")?;
    let mut reference = Vec::new();
    let mut inside = false;
    for line in text.lines() {
        if line.trim_start().starts_with("set ref_data") {
            inside = true;
            continue;
        }
        if inside {
            if line.trim() == "}" {
                break;
            }
            if !line.trim().is_empty() {
                reference.push(line.trim_end().to_string());
            }
        }
    }
    let product_mode_off = text.lines().any(|l| l.trim() == "set ProductMode OFF");
    Ok(Case { file, reference, product_mode_off })
}

/// The statistics of a reading, as the `end` script gathers them.
#[derive(Debug, Clone, Default)]
pub struct Stats {
    /// Faulty wires, faces, shells and solids (`checkshape res ff`).
    pub faulty: [usize; 4],
    /// Solids, shells and faces (`nbshapes`).
    pub counts: [usize; 3],
    /// Solids, shells, faces and free wires (`statshape`).
    pub statshape: [usize; 4],
    /// The largest and mean tolerance (`tolerance`).
    pub tolerance: (f64, f64),
    /// Labels at levels 0, 1 and 2, all, with names, with colours, with layers (`XStat`).
    pub labels: [usize; 7],
    /// Centroid, volume and area properties.
    pub props: [usize; 3],
    /// The colours' names and the layers' names, each sorted.
    pub colors: Vec<String>,
    pub layers: Vec<String>,
}

/// A file read: its assembly tree and its shape (all its products' shapes placed, as one).
pub struct Reading {
    pub nodes: Vec<Node>,
    pub shape: Shape,
}

/// The file at `path` read as STEP or IGES (by its extension).
pub fn read(model: &mut Model, path: &str) -> Result<Reading, String> {
    let bytes = std::fs::read(path).map_err(|e| format!("{path}: {e}"))?;
    read_bytes(model, &bytes, is_iges(path))
}

fn is_iges(path: &str) -> bool {
    let lower = path.to_ascii_lowercase();
    lower.ends_with(".igs") || lower.ends_with(".iges")
}

fn read_bytes(model: &mut Model, bytes: &[u8], iges: bool) -> Result<Reading, String> {
    if iges {
        let shape = kentron::io::iges::read(model, &mut &bytes[..]).map_err(|e| e.to_string())?;
        return Ok(Reading { nodes: vec![Node { shape: Some(shape), ..Default::default() }], shape });
    }
    let nodes = kentron::io::step::read_assembly(model, &mut &bytes[..]).map_err(|e| e.to_string())?;
    let mut shapes = Vec::new();
    fn collect(node: &Node, placement: Transform, out: &mut Vec<Shape>) {
        let here = placement * node.placement;
        if let Some(s) = node.shape {
            out.push(s.moved(&here));
        }
        for child in &node.children {
            collect(child, here, out);
        }
    }
    for root in &nodes {
        collect(root, Transform::IDENTITY, &mut shapes);
    }
    let shape = if shapes.len() == 1 {
        shapes[0]
    } else {
        let compound = model.new_empty(NodeKind::Compound);
        for s in shapes {
            model.add(compound, s).map_err(|e| e.to_string())?;
        }
        compound
    };
    Ok(Reading { nodes, shape })
}

/// The statistics of `reading` (its shape held in `d`'s model).
pub fn stats(d: &mut Draw, reading: &Reading) -> Stats {
    let mut s = Stats::default();
    d.set("res", DrawObject::Shape(reading.shape));
    let run = |d: &mut Draw, cmd: &str| -> String {
        match crate::draw::dispatch::dispatch(d, cmd, &["res"]) {
            Some(Ok(out)) => out,
            _ => String::new(),
        }
    };
    let grab = |text: &str, pattern: &str| -> usize {
        regex::Regex::new(pattern).expect("a pattern").captures(text).and_then(|c| c[1].parse().ok()).unwrap_or(0)
    };
    let nb = run(d, "nbshapes");
    s.counts = [grab(&nb, r"SOLID +: +(\d+)"), grab(&nb, r"SHELL +: +(\d+)"), grab(&nb, r"FACE +: +(\d+)")];
    let st = run(d, "statshape");
    s.statshape = [grab(&st, r"(\d+)\t+ +SOLID"), grab(&st, r"(\d+)\t+ +SHELL"), grab(&st, r"(\d+)\t+ +FACE"), grab(&st, r"(\d+)\t+ +WIRE \(Free\)")];
    let tol = run(d, "tolerance");
    if let Some(c) = regex::Regex::new(r"Tolerance +MAX=([-0-9.+eE]+) +AVG=([-0-9.+eE]+)").expect("a pattern").captures(&tol) {
        s.tolerance = (c[1].parse().unwrap_or(0.0), c[2].parse().unwrap_or(0.0));
    }
    // (Each faulty sub-shape once, by its kind.)
    let report = kentron::validity::check(&d.model, reading.shape);
    let mut seen = std::collections::HashSet::new();
    for p in report.errors() {
        if !seen.insert(p.node) {
            continue;
        }
        let slot = match p.node.kind {
            NodeKind::Wire => 0,
            NodeKind::Face => 1,
            NodeKind::Shell => 2,
            NodeKind::Solid => 3,
            _ => continue,
        };
        s.faulty[slot] += 1;
    }
    document(&reading.nodes, &mut s);
    s
}

/// The document statistics: OCCT's labels of the file's products, as XStat counts them.
fn document(nodes: &[Node], s: &mut Stats) {
    // Each product once (its first occurrence), in the order met.
    let mut products: Vec<&Node> = Vec::new();
    let mut seen = std::collections::HashSet::new();
    fn gather<'a>(node: &'a Node, products: &mut Vec<&'a Node>, seen: &mut std::collections::HashSet<usize>) {
        let fresh = match node.definition {
            Some(d) => seen.insert(d),
            None => true,
        };
        if fresh {
            products.push(node);
        }
        for child in &node.children {
            gather(child, products, seen);
        }
    }
    for root in nodes {
        gather(root, &mut products, &mut seen);
    }
    let mut levels = [0usize; 3];
    let (mut named, mut coloured, mut layered) = (0, 0, 0);
    let mut colors: Vec<[f64; 3]> = Vec::new();
    let mut layers: Vec<String> = Vec::new();
    let mut add_color = |c: [f64; 3]| {
        if !colors.contains(&c) {
            colors.push(c);
        }
    };
    for product in &products {
        levels[0] += 1;
        named += !product.name.is_empty() as usize;
        if let Some(c) = product.color {
            coloured += 1;
            add_color(c);
        }
        if !product.children.is_empty() {
            for child in &product.children {
                levels[1] += 1;
                named += !child.name.is_empty() as usize;
            }
        } else {
            for sub in &product.subshapes {
                levels[1] += 1;
                named += sub.name.is_some() as usize;
                if let Some(c) = sub.color {
                    coloured += 1;
                    add_color([c[0], c[1], c[2]]);
                }
                if let Some(l) = &sub.layer {
                    layered += 1;
                    if !layers.contains(l) {
                        layers.push(l.clone());
                    }
                }
            }
        }
    }
    let total = levels.iter().sum();
    s.labels = [levels[0], levels[1], levels[2], total, named, coloured, layered];
    let mut names: Vec<String> = colors.iter().map(|&c| color_name(c).to_string()).collect();
    names.sort();
    s.colors = names;
    layers.sort();
    s.layers = layers;
}

/// OCCT's name of a colour (Quantity_Color::Name): the nearest named one, in sRGB.
pub fn color_name(rgb: [f64; 3]) -> &'static str {
    let c = [rgb[0] as f32, rgb[1] as f32, rgb[2] as f32];
    let mut best = (f32::MAX, colors::COLORS[0].0);
    for &(name, v) in colors::COLORS {
        let d = (c[0] - v[0]).powi(2) + (c[1] - v[1]).powi(2) + (c[2] - v[2]).powi(2);
        if d < best.0 {
            best = (d, name);
            if d == 0.0 {
                break;
            }
        }
    }
    best.1
}

/// A Tcl list of `items`, as `lsort` prints it (an item with spaces, or empty, in braces).
fn tcl_list(items: &[String]) -> String {
    items
        .iter()
        .map(|i| if i.is_empty() || i.chars().any(|c| c.is_whitespace() || "{}\"\\[]$;".contains(c)) { format!("{{{i}}}") } else { i.clone() })
        .collect::<Vec<_>>()
        .join(" ")
}

/// The result lines (as the `end` script formats them) of the second reading `b` and the
/// first `a`, but DATA and TPSTAT.
pub fn lines(b: &Stats, a: &Stats) -> Vec<String> {
    let f = &b.faulty;
    let g = &a.faulty;
    vec![
        format!("CHECKSHAPE  : Wires    = {}  ( {} )  Faces    = {}  ( {} )  Shells   = {}  ( {} )   Solids   = {} ( {} )", f[0], g[0], f[1], g[1], f[2], g[2], f[3], g[3]),
        format!("NBSHAPES    : Solid    = {}  ( {} )  Shell    = {}  ( {} )  Face     = {}  ( {} ) ", b.counts[0], a.counts[0], b.counts[1], a.counts[1], b.counts[2], a.counts[2]),
        format!(
            "STATSHAPE   : Solid    = {}  ( {} )  Shell    = {}  ( {} )  Face     = {}  ( {} )   FreeWire = {}  ( {} ) ",
            b.statshape[0], a.statshape[0], b.statshape[1], a.statshape[1], b.statshape[2], a.statshape[2], b.statshape[3], a.statshape[3]
        ),
        format!("TOLERANCE   : MaxTol   = {}  ( {} )  AvgTol   =  {}  (  {} )", b.tolerance.0, a.tolerance.0, b.tolerance.1, a.tolerance.1),
        format!(
            "LABELS      : N0Labels = {}  ( {} )  N1Labels = {}  ( {} )  N2Labels = {}  ( {} )   TotalLabels = {}  ( {} )   NameLabels = {}  ( {} )   ColorLabels = {}  ( {} )   LayerLabels = {}  ( {} )",
            b.labels[0], a.labels[0], b.labels[1], a.labels[1], b.labels[2], a.labels[2], b.labels[3], a.labels[3], b.labels[4], a.labels[4], b.labels[5], a.labels[5], b.labels[6], a.labels[6]
        ),
        format!("PROPS       : Centroid = {}  ( {} )  Volume   = {}  ( {} )  Area     = {}  ( {} )", b.props[0], a.props[0], b.props[1], a.props[1], b.props[2], a.props[2]),
        format!("NCOLORS     : NColors  = {}  ( {} )", b.colors.len(), a.colors.len()),
        format!("COLORS      : Colors   = {}  ( {} )", tcl_list(&b.colors), tcl_list(&a.colors)),
        format!("NLAYERS     : NLayers  = {}  ( {} )", b.layers.len(), a.layers.len()),
        format!("LAYERS      : Layers   = {}  ( {} )", tcl_list(&b.layers), tcl_list(&a.layers)),
    ]
}

/// The numbers in a line (integers or reals), in order.
fn numbers(line: &str) -> Vec<f64> {
    line.split_whitespace().filter_map(|w| w.parse::<f64>().ok()).collect()
}

/// Whether the line `current` meets the reference line `reference` (of the same statistic),
/// by the `end` script's rules.
pub fn meets(reference: &str, current: &str) -> bool {
    let (r, c) = (numbers(reference), numbers(current));
    if reference.starts_with("CHECKSHAPE") {
        return r.len() == c.len() && c.iter().zip(&r).all(|(c, r)| c <= r);
    }
    if reference.starts_with("TOLERANCE") {
        return r.len() == c.len() && c.iter().zip(&r).all(|(c, r)| !(c - r > 1e-4 && c - r > 0.01 * r));
    }
    // (Lines compared as strings; whitespace at their ends aside.)
    reference.trim_end() == current.trim_end()
}

/// Runs a case: the file read, written and read again; the differences from the reference
/// (each as the reference line and kentron's), or the failure reading or writing.
pub fn run(case: &Case) -> Result<Vec<(String, String)>, String> {
    if case.product_mode_off {
        return Err("read.step.product.mode OFF is not supported".into());
    }
    let path = crate::draw::find_data_file(&case.file).ok_or_else(|| format!("test data not found: {}", case.file))?;
    let mut d = Draw::new();
    let first = read(&mut d.model, &path).map_err(|e| format!("First - file was not read: {e}"))?;
    let a = stats(&mut d, &first);
    if is_iges(&path) {
        return Err("Here is writing problem: kentron writes no IGES".into());
    }
    let mut bytes = Vec::new();
    kentron::io::step::write_assembly(&d.model, &first.nodes, &mut bytes).map_err(|e| format!("First - file was not written: {e}"))?;
    let second = read_bytes(&mut d.model, &bytes, false).map_err(|e| format!("Second - file was not read: {e}"))?;
    let b = stats(&mut d, &second);
    let current = lines(&b, &a);
    let mut differences = Vec::new();
    for reference in &case.reference {
        let key = &reference[..reference.find(':').unwrap_or(0)];
        if key.trim().is_empty() || key.starts_with("DATA") || key.starts_with("TPSTAT") {
            continue;
        }
        let Some(line) = current.iter().find(|l| l.starts_with(key)) else { continue };
        if !meets(reference, line) {
            differences.push((reference.clone(), line.clone()));
        }
    }
    Ok(differences)
}
