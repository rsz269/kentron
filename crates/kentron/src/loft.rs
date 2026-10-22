//! Lofting: a shell or solid passing through a series of sections.
//!
//! The sections are wires (or vertices at either end, where the loft closes to a point).
//! Each face of the result runs between corresponding edges of consecutive sections:
//!
//! 1. The sections are made compatible ([`Matching`]): wires with fewer edges are split (at
//!    the fractions of their length where the others have vertices, or where rays from the
//!    sections' centres through the others' vertices meet them, as OCCT matches closed
//!    sections), and each wire is turned to run like the previous one and, if closed, to
//!    start at the vertex matching its start. Each piece remembers the section edge it
//!    comes from, for [`Loft::generated`].
//! 2. Corresponding edges are sampled at the same fractions of their parameter ranges and
//!    interpolated by B-splines sharing their knots (lines stay exact, of degree 1).
//! 3. Across the sections, the poles are interpolated by one polynomial span (smooth loft,
//!    of degree one less than the number of sections, as OCCT's approximation gives) or
//!    joined pairwise by straight lines (ruled loft), giving one B-spline surface per edge
//!    (per edge and pair of sections when ruled). When the first section is also the last,
//!    the loft closes on itself.
//! 4. The section curves and the curves across them are the surfaces' isoparametric
//!    curves, so faces meet exactly along shared edges. A solid is closed by planar faces
//!    on the end sections.
//!
//! Between two sections, as OCCT does, corresponding lines in one plane make a planar face
//! and corresponding coaxial arcs (at the same angles) a conical, cylindrical or planar
//! one; those faces keep the sections' own lines and circles as edges.

use crate::error::{Error, Result};
use crate::geometry::curve::{BSplineCurve, Curve};
use crate::geometry::curve2::{BSplineCurve2, Curve2};
use crate::geometry::fit;
use crate::geometry::surface::{BSplineSurface, Cone, Cylinder, KnotVector, Plane, Surface};
use crate::geometry::ParametricSurface;
use crate::geometry::ParametricCurve;
use crate::math::{Direction, Frame, Point, Point2, Transform2, Transformable, Transformable2, Vector, Vector2};
use crate::properties;
use crate::topology::{EdgeId, Model, NodeId, NodeKind, Orientation, Shape, VertexId};

/// How [`loft`] joins its sections.
#[doc(alias = "BRepOffsetAPI_ThruSections")]
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct LoftOptions {
    /// Close the loft with faces on its end sections into a solid (planar, or fillings of
    /// sections in no one plane).
    pub solid: bool,
    /// Join consecutive sections by straight lines instead of a smooth surface through all.
    pub ruled: bool,
    /// How the sections' edges are matched up.
    pub matching: Matching,
    /// Where a smooth loft passes through each section, along the way across them.
    pub spacing: Spacing,
}

/// Where a smooth [`loft`] passes through its sections, along its parameter across them
/// (which shapes the surface between them: through three or more sections, it bulges
/// differently).
#[doc(alias = "Approx_ParametrizationType")]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Spacing {
    /// Spaced as the sections are apart (by the mean distance between their points), as
    /// OCCT's lofts are by default.
    #[default]
    ByDistance,
    /// By the square roots of those distances.
    ByRootDistance,
    /// Evenly, whatever the distances (as CadQuery's lofts ask for).
    Even,
}

/// How [`loft`] matches up the edges of its sections (splitting some, so that all have as
/// many).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Matching {
    /// Sections with different numbers of edges are split at the fractions of their length
    /// where the others have vertices; closed ones start at the vertex nearest to the
    /// previous section's start.
    #[default]
    ByLength,
    /// As OCCT matches them: closed sections by rays from their centres, through each
    /// vertex of the neighbouring section, splitting the edge they meet unless they pass
    /// near a vertex (OCCT's polar method), when the numbers of edges differ (a vertex
    /// section counting as one) or all sections are smooth; open ones by length then.
    ByRays,
    /// The sections as they are: their edges matched in order, none turned or restarted
    /// (split by length only if their numbers of edges differ).
    AsGiven,
}

/// A loft, with the shapes generated from its sections' edges and vertices.
#[derive(Debug, Clone, PartialEq)]
pub struct Loft {
    pub shape: Shape,
    generated: Vec<(NodeId, Vec<Shape>)>,
}

impl Loft {
    /// The shapes lofted from `input`: from an edge of a section, the faces running from
    /// its pieces (across all sections); from a vertex, the edges across the sections from
    /// it (all of them from a vertex section). Empty if `input` is not part of a section.
    pub fn generated(&self, input: Shape) -> Vec<Shape> {
        self.generated.iter().filter(|(node, _)| *node == input.node).flat_map(|(_, s)| s.clone()).collect()
    }
}

/// The shell (or solid, with `options.solid`) through `sections`: wires, or vertices as the
/// first or last section.
#[doc(alias = "BRepOffsetAPI_ThruSections")]
#[doc(alias = "thrusections")]
pub fn loft(model: &mut Model, sections: &[Shape], options: &LoftOptions) -> Result<Loft> {
    loft_traced(model, sections, options, &[], false, false)
}

/// [`loft`] of sections closing on themselves (the last the first again), periodic across
/// them: smooth where it closes, as a tool round a closed path is. (A loft closing on itself
/// otherwise passes through the first section again as through any other, as OCCT's does.)
pub(crate) fn loft_periodic(model: &mut Model, sections: &[Shape], options: &LoftOptions) -> Result<Loft> {
    loft_traced(model, sections, options, &[], false, true)
}

/// [`loft`], its faces through circle arcs that all sections have at the same place exact:
/// rational across them, their homogeneous poles interpolated along the loft, the arcs the
/// section edges as they are (where the arcs' B-spline forms match; else as [`loft`] has
/// them, sampled).
pub(crate) fn loft_exact_arcs(model: &mut Model, sections: &[Shape], options: &LoftOptions) -> Result<Loft> {
    loft_traced(model, sections, options, &[], true, false)
}

/// [`loft`], a smooth one smooth only between the sections at `breaks` (their indices): a face
/// per edge and stretch between breaks, meeting the next at the break's edges.
pub(crate) fn loft_breaking(model: &mut Model, sections: &[Shape], options: &LoftOptions, breaks: &[usize]) -> Result<Shape> {
    loft_traced(model, sections, options, breaks, false, false).map(|l| l.shape)
}

/// `sections` (wires, or vertices first and last) made to match as a loft matches them:
/// wires of as many edges, split where `matching` says, turned and started alike. Vertex
/// sections as they are.
pub(crate) fn compatible_sections(model: &mut Model, sections: &[Shape], matching: Matching) -> Result<Vec<Shape>> {
    let mut marks = Marks::default();
    let mut parsed = sections.iter().map(|&s| Section::read(model, s, &mut marks)).collect::<Result<Vec<_>>>()?;
    let closed = parsed.iter().find_map(|s| if let Section::Wire(w) = s { Some(w.closed) } else { None }).unwrap_or(false);
    if parsed.iter().any(|s| matches!(s, Section::Wire(w) if w.closed != closed)) {
        return Err(Error::InvalidInput("sections must be all closed or all open".into()));
    }
    make_compatible(&mut parsed, closed, false, matching, &mut marks);
    let mut out = Vec::with_capacity(sections.len());
    for (section, &original) in parsed.iter().zip(sections) {
        match section {
            Section::Point(..) => out.push(original),
            Section::Wire(w) => {
                let mut edges = Vec::with_capacity(w.pieces.len());
                for piece in &w.pieces {
                    let (curve, first, last) = piece.forward();
                    edges.push(crate::make::edge_from_curve_range(model, curve, first, last)?);
                }
                out.push(crate::make::wire_within(model, &edges, 1e3 * crate::tolerance::LINEAR)?);
            }
        }
    }
    Ok(out)
}

fn loft_traced(model: &mut Model, sections: &[Shape], options: &LoftOptions, breaks: &[usize], exact_arcs: bool, periodic: bool) -> Result<Loft> {
    if sections.len() < 2 {
        return Err(Error::InvalidInput("a loft needs two sections or more".into()));
    }
    let mut marks = Marks::default();
    let mut parsed = Vec::with_capacity(sections.len());
    for (k, &s) in sections.iter().enumerate() {
        let section = Section::read(model, s, &mut marks)?;
        if matches!(section, Section::Point(..)) && k != 0 && k != sections.len() - 1 {
            return Err(Error::InvalidInput("only the first and last sections of a loft can be vertices".into()));
        }
        parsed.push(section);
    }
    let closed = match parsed.iter().find_map(|s| if let Section::Wire(w) = s { Some(w.closed) } else { None }) {
        Some(c) => c,
        None => return Err(Error::InvalidInput("a loft needs a wire among its sections".into())),
    };
    if parsed.iter().any(|s| matches!(s, Section::Wire(w) if w.closed != closed)) {
        return Err(Error::InvalidInput("a loft's wires must be all closed or all open".into()));
    }
    if options.solid && !closed {
        return Err(Error::InvalidInput("a solid loft needs closed sections".into()));
    }
    // The same wire first and last: the loft closes on itself, like a ring.
    let looped = sections.len() > 2 && sections[0] == sections[sections.len() - 1];
    if looped && options.solid {
        return Err(Error::InvalidInput("a loft closing on itself has no end sections to make a solid".into()));
    }
    make_compatible(&mut parsed, closed, looped, options.matching, &mut marks);
    let grid = Grid::sample(&parsed, options.spacing, exact_arcs);
    let built = build(model, &parsed, &grid, closed, looped, periodic && looped, options, breaks)?;
    // What each section's edges and vertices generated: the faces of their pieces' columns
    // and the edges across at their boundaries.
    let mut generated: Vec<(NodeId, Vec<Shape>)> = Vec::new();
    let mut add = |node: NodeId, shapes: &[Shape]| {
        let at = match generated.iter().position(|(n, _)| *n == node) {
            Some(at) => at,
            None => {
                generated.push((node, Vec::new()));
                generated.len() - 1
            }
        };
        for s in shapes {
            if !generated[at].1.contains(s) {
                generated[at].1.push(*s);
            }
        }
    };
    for section in &parsed {
        match section {
            Section::Point(_, Some(vertex)) => add(*vertex, &built.across.concat()),
            Section::Point(_, None) => {}
            Section::Wire(w) => {
                for (i, piece) in w.pieces.iter().enumerate() {
                    if let Some(edge) = piece.edge {
                        add(edge, &built.columns[i]);
                    }
                    if let Some(vertex) = marks.sources[piece.start] {
                        add(vertex, &built.across[i]);
                    }
                    if !w.closed && i == w.pieces.len() - 1 {
                        if let Some(vertex) = marks.sources[piece.end] {
                            add(vertex, &built.across[i + 1]);
                        }
                    }
                }
            }
        }
    }
    Ok(Loft { shape: built.shape, generated })
}

/// One edge of a section (or part of one): its curve over a range, in the direction of the
/// wire.
#[derive(Clone)]
struct Piece {
    curve: Curve,
    first: f64,
    last: f64,
    /// The vertices at its start and end, as [`Marks`] number them.
    start: usize,
    end: usize,
    /// The section's edge it is (part of).
    edge: Option<NodeId>,
}

impl Piece {
    /// The point at fraction `s` of the range.
    fn at(&self, s: f64) -> Point {
        self.curve.point(self.first + (self.last - self.first) * s)
    }

    fn straight(&self) -> bool {
        matches!(self.curve, Curve::Line(_))
    }

    fn reversed(&self) -> Piece {
        Piece { curve: self.curve.clone(), first: self.last, last: self.first, start: self.end, end: self.start, edge: self.edge }
    }

    /// The pieces before and after fraction `s`, meeting at vertex `mark`.
    fn split(&self, s: f64, mark: usize) -> (Piece, Piece) {
        let t = self.first + (self.last - self.first) * s;
        (Piece { last: t, end: mark, ..self.clone() }, Piece { first: t, start: mark, ..self.clone() })
    }

    /// The curve with its range running forwards (the curve reversed if need be).
    fn forward(&self) -> (Curve, f64, f64) {
        if self.first <= self.last {
            return (self.curve.clone(), self.first, self.last);
        }
        let mut curve = self.curve.clone();
        let (first, last) = (curve.reversed_parameter(self.first), curve.reversed_parameter(self.last));
        curve.reverse();
        (curve, first, last)
    }

    fn length(&self) -> f64 {
        (0..16).map(|i| self.at(i as f64 / 16.0).distance(self.at((i + 1) as f64 / 16.0))).sum()
    }
}

/// The sections' vertices, numbered, each with the vertex of the input it is (none where
/// edges were split).
#[derive(Default)]
struct Marks {
    sources: Vec<Option<NodeId>>,
}

impl Marks {
    fn add(&mut self, source: Option<NodeId>) -> usize {
        self.sources.push(source);
        self.sources.len() - 1
    }

    /// The number of input vertex `source` (the same for each use of it).
    fn of(&mut self, source: NodeId) -> usize {
        match self.sources.iter().position(|s| *s == Some(source)) {
            Some(mark) => mark,
            None => self.add(Some(source)),
        }
    }
}

#[derive(Clone)]
struct WireSection {
    pieces: Vec<Piece>,
    closed: bool,
}

#[derive(Clone)]
enum Section {
    Wire(WireSection),
    /// A vertex section: its point, and the input vertex.
    Point(Point, Option<NodeId>),
}

impl Section {
    fn read(model: &Model, shape: Shape, marks: &mut Marks) -> Result<Section> {
        if shape.kind() == NodeKind::Vertex {
            return Ok(Section::Point(model.point(shape), Some(shape.node)));
        }
        let mut pieces = Vec::new();
        for edge in model.explore(shape, NodeKind::Edge) {
            let Some((curve, first, last)) = model.edge_curve(edge) else { continue };
            let (start, end) = model.edge_vertices(edge);
            let piece = Piece { curve, first, last, start: 0, end: 0, edge: Some(edge.node) };
            let mut piece = if edge.orientation == Orientation::Reversed { piece.reversed() } else { piece };
            (piece.start, piece.end) = (marks.of(start.node), marks.of(end.node));
            pieces.push(piece);
        }
        if pieces.is_empty() {
            return match model.explore(shape, NodeKind::Vertex).next() {
                Some(v) => Ok(Section::Point(model.point(v), Some(v.node))),
                None => Err(Error::InvalidInput("an empty loft section".into())),
            };
        }
        let size = pieces.iter().map(Piece::length).sum::<f64>();
        let tolerance = 1e-6 * size.max(1.0);
        let mut pieces = chain(pieces, tolerance)?;
        let count = pieces.len();
        let gap = pieces[0].at(0.0).distance(pieces[count - 1].at(1.0));
        let closed = gap <= tolerance;
        if closed {
            pieces[count - 1].end = pieces[0].start;
        }
        Ok(Section::Wire(WireSection { pieces, closed }))
    }
}

/// The pieces in the order they connect, each turned to run along the chain (a wire used
/// reversed lists its edges backwards).
fn chain(mut rest: Vec<Piece>, tolerance: f64) -> Result<Vec<Piece>> {
    let mut chain = vec![rest.remove(0)];
    while !rest.is_empty() {
        let (head, tail) = (chain[0].at(0.0), chain[chain.len() - 1].at(1.0));
        let found = rest.iter().enumerate().find_map(|(i, p)| {
            let (a, b) = (p.at(0.0), p.at(1.0));
            if a.distance(tail) <= tolerance {
                Some((i, false, false))
            } else if b.distance(tail) <= tolerance {
                Some((i, true, false))
            } else if b.distance(head) <= tolerance {
                Some((i, false, true))
            } else if a.distance(head) <= tolerance {
                Some((i, true, true))
            } else {
                None
            }
        });
        let Some((i, reverse, before)) = found else {
            return Err(Error::InvalidInput("a loft section's edges do not connect".into()));
        };
        let piece = rest.remove(i);
        let piece = if reverse { piece.reversed() } else { piece };
        if before {
            chain.insert(0, piece);
        } else {
            chain.push(piece);
        }
    }
    Ok(chain)
}

fn wire(s: &Section) -> &WireSection {
    match s {
        Section::Wire(w) => w,
        Section::Point(..) => unreachable!(),
    }
}

fn wire_mut(s: &mut Section) -> &mut WireSection {
    match s {
        Section::Wire(w) => w,
        Section::Point(..) => unreachable!(),
    }
}

/// Splits wires to equal edge counts, and orients and rotates each wire to follow the
/// previous one.
fn make_compatible(sections: &mut Vec<Section>, closed: bool, looped: bool, matching: Matching, marks: &mut Marks) {
    let wires: Vec<usize> = (0..sections.len()).filter(|&k| matches!(sections[k], Section::Wire(_))).collect();
    let counts: Vec<usize> = wires.iter().map(|&k| wire(&sections[k]).pieces.len()).collect();
    let split = match matching {
        Matching::ByLength => counts.iter().any(|&c| c != counts[0]),
        Matching::AsGiven => {
            if counts.iter().any(|&c| c != counts[0]) {
                let all = all_fractions(sections, &wires);
                for &k in &wires {
                    split_at_fractions(wire_mut(&mut sections[k]), &all, marks);
                }
            }
            return;
        }
        Matching::ByRays => {
            // As OCCT: by length (or rays) unless all sections have as many edges (a vertex
            // section one) and some is not smooth.
            let edges = |s: &Section| match s {
                Section::Wire(w) => w.pieces.len(),
                Section::Point(..) => 1,
            };
            let differ = sections.iter().any(|s| edges(s) != edges(&sections[0]));
            let report = differ || wires.iter().all(|&k| smooth(wire(&sections[k])));
            if closed && report && !looped {
                let saved = sections.clone();
                if match_by_rays(sections, &wires, marks) {
                    return;
                }
                *sections = saved;
                counts.iter().any(|&c| c != counts[0])
            } else {
                report && !closed && counts.iter().any(|&c| c > 1)
            }
        }
    };
    // Orientation first (it changes where the vertices fall along the wire).
    for pair in wires.windows(2) {
        let (prev, this) = (wire(&sections[pair[0]]), wire(&sections[pair[1]]));
        let flip = if closed {
            let (a, b) = (normal(prev), normal(this));
            // (Where the planes are near square to each other, the normals do not tell: the
            // way round whose vertices, from the best start, lie nearest the previous one's.)
            if a.dot(b).abs() < 0.5 * a.length() * b.length() && prev.pieces.len() == this.pieces.len() {
                let reversed: Vec<Piece> = this.pieces.iter().rev().map(Piece::reversed).collect();
                gap_from_best_start(prev, &reversed) < gap_from_best_start(prev, &this.pieces)
            } else {
                a.dot(b) < 0.0
            }
        } else {
            let (a0, a1) = (prev.pieces[0].at(0.0), prev.pieces[prev.pieces.len() - 1].at(1.0));
            let (b0, b1) = (this.pieces[0].at(0.0), this.pieces[this.pieces.len() - 1].at(1.0));
            a0.distance(b1) + a1.distance(b0) < a0.distance(b0) + a1.distance(b1)
        };
        if flip {
            let w = wire_mut(&mut sections[pair[1]]);
            w.pieces = w.pieces.iter().rev().map(Piece::reversed).collect();
        }
    }
    // Equal edge counts: split each wire where the others have vertices (by length).
    if split {
        let all = all_fractions(sections, &wires);
        for &k in &wires {
            split_at_fractions(wire_mut(&mut sections[k]), &all, marks);
        }
    }
    // A closed wire starts at the vertex nearest to the previous wire's start; or, where the
    // two lie in planes near square to each other, where its vertices in order lie nearest
    // the previous wire's in all (of an outline turned aside, the vertex nearest the previous
    // start may be its far end).
    if closed {
        for pair in wires.windows(2) {
            let prev = wire(&sections[pair[0]]).clone();
            let aside = {
                let (a, b) = (normal(&prev), normal(wire(&sections[pair[1]])));
                a.dot(b).abs() < 0.5 * a.length() * b.length()
            };
            let w = wire_mut(&mut sections[pair[1]]);
            let n = w.pieces.len();
            let gap = |i: usize| -> f64 {
                if aside && prev.pieces.len() == n {
                    (0..n).map(|k| prev.pieces[k].at(0.0).distance(w.pieces[(k + i) % n].at(0.0))).sum()
                } else {
                    w.pieces[i].at(0.0).distance(prev.pieces[0].at(0.0))
                }
            };
            if let Some(i) = (0..n).min_by(|&a, &b| gap(a).total_cmp(&gap(b))) {
                w.pieces.rotate_left(i);
            }
        }
    }
}

/// The least sum of distances between `prev`'s vertices and `pieces`' taken in order from
/// some start.
fn gap_from_best_start(prev: &WireSection, pieces: &[Piece]) -> f64 {
    let n = pieces.len();
    (0..n)
        .map(|i| (0..n).map(|k| prev.pieces[k].at(0.0).distance(pieces[(k + i) % n].at(0.0))).sum::<f64>())
        .fold(f64::INFINITY, f64::min)
}

/// The fractions of their length where the wires have vertices, of all wires.
fn all_fractions(sections: &[Section], wires: &[usize]) -> Vec<f64> {
    let fractions = |w: &WireSection| -> Vec<f64> {
        let lengths: Vec<f64> = w.pieces.iter().map(Piece::length).collect();
        let total: f64 = lengths.iter().sum();
        let mut at = 0.0;
        lengths[..lengths.len() - 1].iter().map(|l| {
            at += l;
            at / total
        }).collect()
    };
    let mut all: Vec<f64> = wires.iter().flat_map(|&k| fractions(wire(&sections[k]))).collect();
    all.sort_by(f64::total_cmp);
    all.dedup_by(|a, b| (*a - *b).abs() < 1e-9);
    all
}

/// Whether the wire's edges meet with the same tangent vectors (each end too, if closed).
fn smooth(w: &WireSection) -> bool {
    let n = w.pieces.len();
    let joints = if w.closed { n } else { n - 1 };
    (0..joints).all(|i| {
        let (a, b) = (&w.pieces[i], &w.pieces[(i + 1) % n]);
        let tangent = |p: &Piece, t: f64| p.curve.d1(t).1 * (p.last - p.first).signum();
        let (ta, tb) = (tangent(a, a.last), tangent(b, b.first));
        (ta - tb).length() <= 1e-9 * ta.length().max(tb.length())
    })
}

/// Splits `w`'s pieces so that it has a vertex at each of `fractions` of its length (those
/// not already at a vertex).
fn split_at_fractions(w: &mut WireSection, fractions: &[f64], marks: &mut Marks) {
    let lengths: Vec<f64> = w.pieces.iter().map(Piece::length).collect();
    let total: f64 = lengths.iter().sum();
    let mut out = Vec::new();
    let mut start = 0.0;
    for (piece, length) in w.pieces.iter().zip(&lengths) {
        let end = start + length / total;
        let mut rest = piece.clone();
        let mut rest_from = start;
        for &f in fractions.iter().filter(|&&f| f > start + 1e-9 && f < end - 1e-9) {
            // The fraction of what remains of the piece, by parameter (proportional to
            // length for lines and circles).
            let s = (f - rest_from) / (end - rest_from);
            let (a, b) = rest.split(s, marks.add(None));
            out.push(a);
            rest = b;
            rest_from = f;
        }
        out.push(rest);
        start = end;
    }
    w.pieces = out;
}

/// Matches closed sections as OCCT does (`BRepFill_CompatibleWires::SameNumberByPolarMethod`)
/// by rays from each section's centre: through each vertex of the next section (carried into
/// this one's plane), from the last section down to the first, then through each vertex
/// still without a match on the next section, back up. A ray splits the edge it meets, unless
/// it passes within a tenth of the edge from one of its vertices not matched yet. Each wire
/// is then turned and started to follow the first. False if the matching fails.
fn match_by_rays(sections: &mut [Section], wires: &[usize], marks: &mut Marks) -> bool {
    let (first, last) = (wires[0], wires[wires.len() - 1]);
    let mut planes = vec![(Point::ORIGIN, Vector::ZERO); sections.len()];
    for &k in wires {
        let w = wire(&sections[k]);
        let Some(axis) = Direction::try_from_vector(normal(w)) else { return false };
        planes[k] = (centre(w), Vector::from(axis));
    }
    // Down: each vertex of section k matched on section k - 1 (`down[v]`), in turn.
    let mut down: Vec<(usize, Vec<usize>)> = Vec::new();
    for k in (first + 1..=last).rev() {
        for v in vertex_sequence(wire(&sections[k])) {
            down.retain(|(key, _)| *key != v);
            down.push((v, Vec::new()));
            let towards = carry(vertex_point(wire(&sections[k]), v), planes[k], planes[k - 1]);
            if towards.distance(planes[k - 1].0) > crate::tolerance::LINEAR {
                let taken = |m: usize| down.iter().any(|(_, list)| list.contains(&m));
                let Some(hit) = meet(wire_mut(&mut sections[k - 1]), planes[k - 1].0, towards, &taken, marks) else { return false };
                down.last_mut().expect("pushed").1.push(hit);
            }
        }
    }
    // The vertex above `v`: the one matched down onto it.
    let above = |v: usize| down.iter().find(|(_, list)| list.contains(&v)).map(|(key, _)| *key);
    // Up: chains of matching vertices from each vertex of the first section, completed on
    // the sections where they have none.
    let mut chains: Vec<Vec<usize>> = Vec::new();
    for v in vertex_sequence(wire(&sections[first])) {
        let mut chain = vec![v];
        let mut at = v;
        while let Some(next) = above(at) {
            chain.push(next);
            if next == v || chain.len() >= sections.len() {
                break;
            }
            at = next;
        }
        chains.push(chain);
    }
    for k in first..last {
        for v in vertex_sequence(wire(&sections[k])) {
            let Some(c) = chains.iter().position(|chain| chain.contains(&v)) else { return false };
            let next = wire(&sections[k + 1]);
            if chains[c].iter().any(|&m| has_vertex(next, m)) {
                continue;
            }
            let towards = carry(vertex_point(wire(&sections[k]), v), planes[k], planes[k + 1]);
            if towards.distance(planes[k + 1].0) > crate::tolerance::LINEAR {
                let taken = |m: usize| chains.iter().any(|chain| chain.contains(&m));
                let Some(hit) = meet(wire_mut(&mut sections[k + 1]), planes[k + 1].0, towards, &taken, marks) else { return false };
                chains[c].push(hit);
            }
        }
    }
    // Each wire turned and started to follow the first: its first edge between the matches of
    // the first wire's first edge's ends, the one closest in direction (carried back down).
    let (head, head_from) = (wire(&sections[first]).pieces[0].clone(), planes[first].0);
    let towards_head = (head.at(0.9) - head_from).normalized();
    for k in first + 1..=last {
        let rank = k - first;
        let matched = |v: usize| chains.iter().find(|chain| chain[0] == v).and_then(|chain| chain.get(rank).copied());
        let (Some(v1), Some(v2)) = (matched(head.start), matched(head.end)) else { return false };
        let w = wire(&sections[k]);
        let closeness = |p: Point| {
            let p = (first + 1..=k).rev().fold(p, |p, j| carry(p, planes[j], planes[j - 1]));
            (p - head_from).normalized().dot(towards_head)
        };
        let mut best: Option<(usize, bool, f64)> = None;
        for (j, piece) in w.pieces.iter().enumerate() {
            for (backwards, from, to, along) in [(false, piece.start, piece.end, 0.9), (true, piece.end, piece.start, 0.1)] {
                if from == v1 && to == v2 && best.is_none_or(|b| closeness(piece.at(along)) > b.2) {
                    best = Some((j, backwards, closeness(piece.at(along))));
                }
            }
        }
        let Some((j, backwards, _)) = best else { return false };
        let w = wire_mut(&mut sections[k]);
        if backwards {
            w.pieces = w.pieces.iter().rev().map(Piece::reversed).collect();
            let j = w.pieces.len() - 1 - j;
            w.pieces.rotate_left(j);
        } else {
            w.pieces.rotate_left(j);
        }
    }
    let count = wire(&sections[first]).pieces.len();
    wires.iter().all(|&k| wire(&sections[k]).pieces.len() == count)
}

/// The distinct vertices of a wire, in the order OCCT lists them (`SeqOfVertices`): along the
/// wire, each edge's in the direction of its curve.
fn vertex_sequence(w: &WireSection) -> Vec<usize> {
    let mut sequence = Vec::new();
    for piece in &w.pieces {
        let ends = if piece.first <= piece.last { [piece.start, piece.end] } else { [piece.end, piece.start] };
        for v in ends {
            if !sequence.contains(&v) {
                sequence.push(v);
            }
        }
    }
    sequence
}

fn has_vertex(w: &WireSection, v: usize) -> bool {
    w.pieces.iter().any(|p| p.start == v || p.end == v)
}

fn vertex_point(w: &WireSection, v: usize) -> Point {
    let piece = w.pieces.iter().find(|p| p.start == v || p.end == v).expect("a vertex of the wire");
    piece.at(if piece.start == v { 0.0 } else { 1.0 })
}

/// The centre of a closed wire, as OCCT takes it for matching: a lone circle's or ellipse's
/// centre, otherwise the centre of mass of its length.
fn centre(w: &WireSection) -> Point {
    match (w.pieces.len(), &w.pieces[0].curve) {
        (1, Curve::Circle(c)) => return c.center(),
        (1, Curve::Ellipse(e)) => return e.center(),
        _ => {}
    }
    let (mut sum, mut length) = (Vector::ZERO, 0.0);
    for piece in &w.pieces {
        for i in 0..32 {
            let (a, b) = (piece.at(i as f64 / 32.0), piece.at((i + 1) as f64 / 32.0));
            let l = a.distance(b);
            sum = sum + ((a - Point::ORIGIN) + (b - Point::ORIGIN)) * (0.5 * l);
            length += l;
        }
    }
    Point::ORIGIN + sum / length
}

/// `p` carried from one section's plane (centre and normal) to another's: moved by the
/// offset between the centres and, if the normals differ, turned about the second centre by
/// the angle between them (OCCT's `Transform` in BRepFill_CompatibleWires).
fn carry(p: Point, (from, from_axis): (Point, Vector), (to, to_axis): (Point, Vector)) -> Point {
    let moved = p + (to - from);
    let angle = from_axis.angle(to_axis);
    if angle <= 1e-4 || angle >= std::f64::consts::PI - 1e-4 {
        return moved;
    }
    let offset = to - from;
    let a1 = if offset.dot(from_axis) < -1e-7 { -from_axis } else { from_axis };
    let a2 = if offset.dot(to_axis) < 1e-7 { -to_axis } else { to_axis };
    let Some(about) = Direction::try_from_vector(a1.cross(a2)) else { return moved };
    let sign_test = a2 * offset.dot(a1) - a1 * offset.dot(a2);
    let (alpha, beta) = (sign_test.dot(a1), sign_test.dot(a2));
    let sign = if alpha * beta > 0.0 && alpha.abs() > 1e-4 && beta.abs() > 1e-4 { -1.0 } else { 1.0 };
    let turn = crate::math::Transform::rotation(crate::math::Axis::new(to, about), sign * a1.angle(a2));
    moved.transformed(&turn)
}

/// Where the ray from `from` through `towards` meets wire `w` (the meeting nearest to
/// `towards`), as a vertex of the wire: the vertex it passes through, or near (within a
/// tenth of the edge's range) unless `taken` holds for that one, or else a new vertex
/// splitting the edge there (OCCT's `EdgeIntersectOnWire`).
fn meet(w: &mut WireSection, from: Point, towards: Point, taken: &dyn Fn(usize) -> bool, marks: &mut Marks) -> Option<usize> {
    let samples: Vec<Point> = w.pieces.iter().flat_map(|p| (0..=16).map(move |i| p.at(i as f64 / 16.0))).collect();
    let (lo, hi) = samples.iter().fold((samples[0], samples[0]), |(lo, hi), p| {
        (Point::new(lo.x.min(p.x), lo.y.min(p.y), lo.z.min(p.z)), Point::new(hi.x.max(p.x), hi.y.max(p.y), hi.z.max(p.z)))
    });
    let reach = lo.distance(hi);
    let direction = Direction::try_from_vector(towards - from)?;
    let ray_distance = |p: Point| {
        let t = (p - from).dot(Vector::from(direction)).clamp(0.0, reach);
        p.distance(from + Vector::from(direction) * t)
    };
    // The local minima of the distance along each piece: (distance, piece, fraction).
    let mut minima: Vec<(f64, usize, f64)> = Vec::new();
    const N: usize = 64;
    for (j, piece) in w.pieces.iter().enumerate() {
        let d: Vec<f64> = (0..=N).map(|i| ray_distance(piece.at(i as f64 / N as f64))).collect();
        for i in 0..=N {
            if (i > 0 && d[i - 1] < d[i]) || (i < N && d[i + 1] < d[i]) {
                continue;
            }
            // Golden-section search between the neighbouring samples.
            let (mut a, mut b) = (i.saturating_sub(1) as f64 / N as f64, (i + 1).min(N) as f64 / N as f64);
            let g = 0.5 * (5f64.sqrt() - 1.0);
            for _ in 0..60 {
                let (x1, x2) = (b - g * (b - a), a + g * (b - a));
                if ray_distance(piece.at(x1)) < ray_distance(piece.at(x2)) {
                    b = x2;
                } else {
                    a = x1;
                }
            }
            let s = [0.5 * (a + b), i as f64 / N as f64].into_iter().min_by(|&x, &y| ray_distance(piece.at(x)).total_cmp(&ray_distance(piece.at(y))))?;
            minima.push((ray_distance(piece.at(s)), j, s));
        }
    }
    let least = minima.iter().map(|m| m.0).fold(f64::INFINITY, f64::min);
    let &(_, j, s) = minima
        .iter()
        .filter(|m| m.0 <= least + crate::tolerance::LINEAR)
        .min_by(|x, y| w.pieces[x.1].at(x.2).distance(towards).total_cmp(&w.pieces[y.1].at(y.2).distance(towards)))?;
    let piece = &w.pieces[j];
    let point = piece.at(s);
    // Through a vertex.
    for (end, v) in [(0.0, piece.start), (1.0, piece.end)] {
        if point.distance(piece.at(end)) <= crate::tolerance::LINEAR {
            return Some(v);
        }
    }
    // Near one not matched yet.
    let near = if s < 0.1 { Some(piece.start) } else if s > 0.9 { Some(piece.end) } else { None };
    if let Some(v) = near.filter(|&v| !taken(v)) {
        return Some(v);
    }
    let mark = marks.add(None);
    let (a, b) = piece.split(s, mark);
    w.pieces.splice(j..=j, [a, b]);
    Some(mark)
}

/// The Newell normal of a closed wire (its enclosed area's vector).
fn normal(w: &WireSection) -> Vector {
    let points: Vec<Point> = w.pieces.iter().flat_map(|p| (0..8).map(move |i| p.at(i as f64 / 8.0))).collect();
    let n = points.len();
    (0..n).fold(Vector::ZERO, |acc, i| acc + (points[i] - Point::ORIGIN).cross(points[(i + 1) % n] - Point::ORIGIN) * 0.5)
}

/// Per edge: the sections' curves as B-splines with shared knots, the parameters along the
/// edge and across the sections.
struct Grid {
    /// `curves[i][k]`: edge `i` of section `k`, over [0, 1].
    curves: Vec<Vec<BSplineCurve>>,
    /// Whether edge `i`'s curves are its sections' arcs exactly (rational), not samples.
    exact: Vec<bool>,
    /// The sections' positions across the loft, from 0 to 1.
    across: Vec<f64>,
}

impl Grid {
    fn sample(sections: &[Section], spacing: Spacing, exact_arcs: bool) -> Grid {
        let edges = sections.iter().find_map(|s| if let Section::Wire(w) = s { Some(w.pieces.len()) } else { None }).expect("a wire");
        let mut curves: Vec<Vec<BSplineCurve>> = Vec::with_capacity(edges);
        let mut exact = vec![false; edges];
        // (The rows the spacing across is measured on: the sampled ones also where the
        // surface has exact rows, which have far fewer poles.)
        let mut measured: Vec<Option<Vec<BSplineCurve>>> = Vec::with_capacity(edges);
        for i in 0..edges {
            let exact_row = exact_arcs.then(|| exact_rows(sections, i)).flatten();
            let straight = sections.iter().all(|s| match s {
                Section::Wire(w) => w.pieces[i].straight(),
                Section::Point(..) => true,
            });
            // Each section's edge sampled 64 times a turn (at least 64 times), in powers of two
            // (as one winding round many times, a helix, needs to be told apart turn by turn):
            // the same edge the same way whatever it is lofted to (faces lofted from it apart
            // meet along one curve); the coarser ones' knots then refined to the finest's (a
            // fit's knots lie among those of one twice as fine).
            let samples = |s: &Section| -> usize {
                let Section::Wire(w) = s else { return 1 };
                if straight {
                    return 1;
                }
                const PROBE: usize = 1024;
                let points: Vec<Point> = (0..=PROBE).map(|j| w.pieces[i].at(j as f64 / PROBE as f64)).collect();
                let turn: f64 = points.windows(3).map(|q| (q[1] - q[0]).angle(q[2] - q[1])).filter(|a| a.is_finite()).sum();
                let mut n = 64;
                while (n as f64) < turn / std::f64::consts::TAU * 64.0 && n < 8192 {
                    n *= 2;
                }
                n
            };
            let counts: Vec<usize> = sections.iter().map(samples).collect();
            let finest = counts.iter().copied().max().expect("a section");
            let row = |(s, &n): (&Section, &usize)| -> BSplineCurve {
                let n = if matches!(s, Section::Point(..)) { finest } else { n };
                let params: Vec<f64> = (0..=n).map(|j| j as f64 / n as f64).collect();
                let points: Vec<Point> = match s {
                    Section::Wire(w) => params.iter().map(|&t| w.pieces[i].at(t)).collect(),
                    Section::Point(p, _) => vec![*p; params.len()],
                };
                let mut curve = fit::interpolate_with_parameters(&points, &params).expect("increasing parameters");
                if n < finest {
                    let params: Vec<f64> = (0..=finest).map(|j| j as f64 / finest as f64).collect();
                    let fine = fit::interpolate_with_parameters(&params.iter().map(|&t| Point::new(t, 0.0, 0.0)).collect::<Vec<_>>(), &params).expect("increasing parameters");
                    let knots: Vec<(f64, usize)> = fine.knots().iter().map(|&k| (k, 1)).collect();
                    curve.insert_knots(&knots);
                }
                curve
            };
            let sampled: Vec<BSplineCurve> = sections.iter().zip(&counts).map(row).collect();
            match exact_row {
                Some(rows) => {
                    curves.push(rows);
                    measured.push(Some(sampled));
                    exact[i] = true;
                }
                None => {
                    curves.push(sampled);
                    measured.push(None);
                }
            }
        }
        // Across: by the mean distance between consecutive sections' points (or as `spacing`
        // says otherwise).
        let mut across = vec![0.0];
        for k in 1..sections.len() {
            let mut sum = 0.0;
            let mut count = 0;
            for (row, sampled) in curves.iter().zip(&measured) {
                let row = sampled.as_ref().unwrap_or(row);
                for j in 0..row[k].poles().len() {
                    sum += row[k].poles()[j].distance(row[k - 1].poles()[j]);
                    count += 1;
                }
            }
            let distance = (sum / count as f64).max(1e-9);
            let step = match spacing {
                Spacing::ByDistance => distance,
                Spacing::ByRootDistance => distance.sqrt(),
                Spacing::Even => 1.0,
            };
            across.push(across[k - 1] + step);
        }
        let total = across[across.len() - 1];
        for v in &mut across {
            *v /= total;
        }
        Grid { curves, exact, across }
    }

    /// The surface of edge `i` over sections `from..=to`: its poles interpolated across the
    /// sections, with the curves across them.
    /// With `periodic` (the rows closing on themselves, the last the first again): periodic
    /// across, through the rows but the last, smooth where it closes.
    fn surface(&self, i: usize, from: usize, to: usize, periodic: bool) -> (Surface, Vec<BSplineCurve>) {
        let rows = &self.curves[i][from..=to];
        let across: Vec<f64> = self.across[from..=to].to_vec();
        let count = rows[0].poles().len();
        let through = |points: &[Point]| -> BSplineCurve {
            if periodic {
                let m = points.len() - 1;
                if let Ok(c) = fit::interpolate_periodic_with_parameters(&points[..m], &across[..m], across[m] - across[0]) {
                    return c;
                }
            }
            polynomial_through(points, &across)
        };
        let knots_across = |c: &BSplineCurve| {
            if c.is_periodic() {
                KnotVector::periodic(c.knots().to_vec(), c.multiplicities().to_vec(), c.degree())
            } else {
                KnotVector::new(c.knots().to_vec(), c.multiplicities().to_vec(), c.degree())
            }
        };
        if rows.iter().any(|r| r.is_rational()) {
            // Rational rows: their homogeneous poles (w P, and w) interpolated alike.
            let weight = |r: &BSplineCurve, j: usize| r.weights().map_or(1.0, |w| w[j]);
            let columns: Vec<(BSplineCurve, BSplineCurve)> = (0..count)
                .map(|j| {
                    let points: Vec<Point> = rows.iter().map(|r| Point::ORIGIN + (r.poles()[j] - Point::ORIGIN) * weight(r, j)).collect();
                    let weights: Vec<Point> = rows.iter().map(|r| Point::new(weight(r, j), 0.0, 0.0)).collect();
                    (through(&points), through(&weights))
                })
                .collect();
            let projected = |(h, w): &(BSplineCurve, BSplineCurve)| -> (Vec<Point>, Vec<f64>) {
                let ws: Vec<f64> = w.poles().iter().map(|p| p.x).collect();
                (h.poles().iter().zip(&ws).map(|(&p, &w)| Point::ORIGIN + (p - Point::ORIGIN) * (1.0 / w)).collect(), ws)
            };
            let (poles, weights): (Vec<Vec<Point>>, Vec<Vec<f64>>) = columns.iter().map(projected).unzip();
            let u = KnotVector::new(rows[0].knots().to_vec(), rows[0].multiplicities().to_vec(), rows[0].degree());
            let (h0, _) = &columns[0];
            let v = knots_across(h0);
            let surface = Surface::BSpline(BSplineSurface::new(poles, u, v).with_weights(weights));
            // The curves across at u = 0 and u = 1 (rational where their weights vary).
            let ends = [0, count - 1].map(|j| {
                let (p, w) = projected(&columns[j]);
                let curve = if h0.is_periodic() {
                    BSplineCurve::new_periodic(p, h0.knots().to_vec(), h0.multiplicities().to_vec(), h0.degree())
                } else {
                    BSplineCurve::new(p, h0.knots().to_vec(), h0.multiplicities().to_vec(), h0.degree())
                };
                if w.iter().all(|&x| x == 1.0) { curve } else { curve.with_weights(w) }
            });
            return (surface, ends.to_vec());
        }
        let columns: Vec<BSplineCurve> = (0..count)
            .map(|j| {
                let points: Vec<Point> = rows.iter().map(|r| r.poles()[j]).collect();
                through(&points)
            })
            .collect();
        let poles: Vec<Vec<Point>> = columns.iter().map(|c| c.poles().to_vec()).collect();
        let u = KnotVector::new(rows[0].knots().to_vec(), rows[0].multiplicities().to_vec(), rows[0].degree());
        let v = knots_across(&columns[0]);
        let surface = Surface::BSpline(BSplineSurface::new(poles, u, v));
        // The curves across at u = 0 and u = 1 interpolate the end points of the rows.
        let ends = [0, count - 1].map(|j| columns[j].clone());
        (surface, ends.to_vec())
    }
}

/// Edge `i` of each section as its circle arc exactly, a rational B-spline over [0, 1] from
/// the piece's start to its end, where every section's is an arc and their B-spline forms
/// have the same degree and knots; else `None`.
fn exact_rows(sections: &[Section], i: usize) -> Option<Vec<BSplineCurve>> {
    let mut rows = Vec::with_capacity(sections.len());
    for s in sections {
        let Section::Wire(w) = s else { return None };
        let piece = &w.pieces[i];
        if !matches!(piece.curve, Curve::Circle(_)) {
            return None;
        }
        let (lo, hi) = (piece.first.min(piece.last), piece.first.max(piece.last));
        let trimmed = Curve::Trimmed(crate::geometry::curve::TrimmedCurve::new(piece.curve.clone(), lo, hi));
        let mut b = crate::geometry::convert::curve_to_bspline(&trimmed)?.clamped();
        if piece.first > piece.last {
            b.reverse();
        }
        let (k0, k1) = (b.knots()[0], b.knots()[b.knots().len() - 1]);
        if !(k1 > k0) {
            return None;
        }
        let knots: Vec<f64> = b.knots().iter().map(|k| (k - k0) / (k1 - k0)).collect();
        let row = BSplineCurve::new(b.poles().to_vec(), knots, b.multiplicities().to_vec(), b.degree());
        rows.push(match b.weights() {
            Some(w) => row.with_weights(w.to_vec()),
            None => row,
        });
    }
    let first = rows.first()?;
    let alike = |r: &BSplineCurve| {
        r.degree() == first.degree() && r.multiplicities() == first.multiplicities() && r.knots().len() == first.knots().len() && r.knots().iter().zip(first.knots()).all(|(a, b)| (a - b).abs() <= 1e-12)
    };
    rows.iter().all(alike).then_some(rows)
}

/// The polynomial curve (one Bezier span, of degree one less than the number of points)
/// through `points` at `params`, as OCCT's loft approximation makes across its sections. With
/// more than nine points, a cubic spline.
fn polynomial_through(points: &[Point], params: &[f64]) -> BSplineCurve {
    let n = points.len();
    if n > 9 {
        return fit::interpolate_with_parameters(points, params).expect("increasing parameters");
    }
    let degree = n - 1;
    let (t0, t1) = (params[0], params[n - 1]);
    // Bernstein polynomials of the degree at each parameter (mapped to [0, 1]).
    let bernstein = |t: f64| -> Vec<f64> {
        let s = (t - t0) / (t1 - t0);
        let mut b = vec![0.0; n];
        b[0] = 1.0;
        for d in 1..=degree {
            for k in (0..=d).rev() {
                let left = if k > 0 { b[k - 1] * s } else { 0.0 };
                b[k] = b[k] * (1.0 - s) + left;
            }
        }
        b
    };
    let matrix: Vec<Vec<f64>> = params.iter().map(|&t| bernstein(t)).collect();
    let rhs: Vec<Vector> = points.iter().map(|&p| p - Point::ORIGIN).collect();
    let poles = solve(matrix, rhs).into_iter().map(|v| Point::ORIGIN + v).collect();
    BSplineCurve::new(poles, vec![t0, t1], vec![degree + 1, degree + 1], degree)
}

/// The solution of the linear system `matrix · x = rhs` (Gaussian elimination with partial
/// pivoting), for vector unknowns.
fn solve(mut matrix: Vec<Vec<f64>>, mut rhs: Vec<Vector>) -> Vec<Vector> {
    let n = rhs.len();
    for col in 0..n {
        let pivot = (col..n).max_by(|&a, &b| matrix[a][col].abs().total_cmp(&matrix[b][col].abs())).expect("rows");
        matrix.swap(col, pivot);
        rhs.swap(col, pivot);
        for row in col + 1..n {
            let f = matrix[row][col] / matrix[col][col];
            for k in col..n {
                matrix[row][k] -= f * matrix[col][k];
            }
            rhs[row] = rhs[row] - rhs[col] * f;
        }
    }
    let mut x = vec![Vector::ZERO; n];
    for row in (0..n).rev() {
        let mut sum = rhs[row];
        for k in row + 1..n {
            sum = sum - x[k] * matrix[row][k];
        }
        x[row] = sum / matrix[row][row];
    }
    x
}

/// A 2D segment from `a` to `b` over parameters `[t0, t1]`.
fn segment(a: Point2, b: Point2, t0: f64, t1: f64) -> Curve2 {
    Curve2::BSpline(BSplineCurve2::new(vec![a, b], vec![t0, t1], vec![2, 2], 1))
}

/// A built loft: the shape, the faces of each column (running from the sections' edges at an
/// index) and the edges across at each boundary between columns.
struct Built {
    shape: Shape,
    columns: Vec<Vec<Shape>>,
    across: Vec<Vec<Shape>>,
}

fn build(model: &mut Model, sections: &[Section], grid: &Grid, closed: bool, looped: bool, periodic: bool, options: &LoftOptions, breaks: &[usize]) -> Result<Built> {
    let edges = grid.curves.len();
    let n = sections.len();
    // The section whose vertices and edges stand for section `k` (the first for the last
    // when the loft closes on itself).
    let same = |k: usize| if looped && k == n - 1 { 0 } else { k };
    let spans: Vec<(usize, usize)> = if options.ruled {
        (0..n - 1).map(|k| (k, k + 1)).collect()
    } else {
        let mut stops: Vec<usize> = std::iter::once(0).chain(breaks.iter().copied().filter(|&b| b > 0 && b < n - 1)).chain(std::iter::once(n - 1)).collect();
        stops.dedup();
        stops.windows(2).map(|w| (w[0], w[1])).collect()
    };
    // Vertices: per section and boundary between edges (one for a point section).
    let boundaries = if closed { edges } else { edges + 1 };
    let vertices: Vec<Vec<VertexId>> = (0..n)
        .map(|k| match &sections[k] {
            Section::Point(p, _) => vec![model.add_vertex(*p, crate::tolerance::LINEAR); boundaries],
            Section::Wire(_) => (0..boundaries)
                .map(|b| {
                    let p = if b < edges { grid.curves[b][k].point(0.0) } else { grid.curves[edges - 1][k].point(1.0) };
                    model.add_vertex(p, crate::tolerance::LINEAR)
                })
                .collect(),
        })
        .collect();
    let vertex = |k: usize, b: usize| vertices[same(k)][b % boundaries];
    // Section edges where faces start or end: (edge, section) -> edge id (degenerate at a
    // point section).
    // Analytic faces between two sections (their section edges keep the exact curves).
    let analytic: Vec<Option<Surface>> = (0..edges)
        .map(|i| match (n, looped, &sections[0], &sections[n - 1]) {
            (2, false, Section::Wire(a), Section::Wire(b)) => ruled_surface(&a.pieces[i], &b.pieces[i]),
            _ => None,
        })
        .collect();
    // (Each with the range of its curve that the face's parameter from 0 to 1 across runs
    // over.)
    let mut section_edges: std::collections::HashMap<(usize, usize), (EdgeId, (f64, f64))> = std::collections::HashMap::new();
    let mut section_edge = |model: &mut Model, i: usize, k: usize| -> Result<(EdgeId, (f64, f64))> {
        let k = same(k);
        if let Some(&e) = section_edges.get(&(i, k)) {
            return Ok(e);
        }
        let id = model.add_edge(vertex(k, i), vertex(k, i + 1), crate::tolerance::LINEAR)?;
        let mut range = (0.0, 1.0);
        if let Section::Wire(w) = &sections[k] {
            let (curve, first, last) = w.pieces[i].forward();
            if analytic[i].is_some() {
                model.set_edge_curve(id.into(), curve, crate::tolerance::LINEAR)?;
                model.set_edge_range(id.into(), first, last)?;
                range = (first, last);
            } else if matches!(curve, Curve::Line(_) | Curve::Circle(_) | Curve::Ellipse(_)) && !grid.exact[i] {
                // A line or conic stays one (as its section's edge, meeting the faces beyond
                // it on the same curve), over the same fraction of its range as the face's
                // boundary across: within the tolerance of that boundary from it.
                let grid_curve = &grid.curves[i][k];
                let off = (0..=64)
                    .map(|j| j as f64 / 64.0)
                    .map(|t| grid_curve.point(t).distance(curve.point(first + (last - first) * t)))
                    .fold(0.0, f64::max);
                model.set_edge_curve(id.into(), curve, crate::tolerance::LINEAR.max(2.0 * off))?;
                model.set_edge_range(id.into(), first, last)?;
                range = (first, last);
            } else {
                model.set_edge_curve(id.into(), Curve::BSpline(grid.curves[i][k].clone()), crate::tolerance::LINEAR)?;
            }
        }
        section_edges.insert((i, k), (id, range));
        Ok((id, range))
    };
    let mut faces = Vec::new();
    let mut columns: Vec<Vec<Shape>> = vec![Vec::new(); edges];
    let mut across: Vec<Vec<Shape>> = vec![Vec::new(); boundaries];
    for &(from, to) in &spans {
        let (v0, v1) = (grid.across[from], grid.across[to]);
        // The edges across, at each boundary (shared by the faces on either side).
        let mut across_edges: Vec<Option<EdgeId>> = vec![None; boundaries];
        let mut surfaces = Vec::with_capacity(edges);
        for i in 0..edges {
            // (Periodic across where asked, a loft closing on itself in one span.)
            let (surface, ends) = grid.surface(i, from, to, periodic && from == 0 && to == n - 1);
            for (b, curve) in [(i, &ends[0]), (i + 1, &ends[1])] {
                let slot = b % boundaries;
                if across_edges[slot].is_none() {
                    let id = model.add_edge(vertex(from, b), vertex(to, b), crate::tolerance::LINEAR)?;
                    // Over the sections' positions, rescaled to the span.
                    model.set_edge_curve(id.into(), Curve::BSpline(curve.clone()), crate::tolerance::LINEAR)?;
                    model.set_edge_range(id.into(), v0, v1)?;
                    across_edges[slot] = Some(id);
                    across[slot].push(id.into());
                }
            }
            surfaces.push(surface);
        }
        for (i, surface) in surfaces.into_iter().enumerate() {
            let face = model.new_empty(NodeKind::Face);
            let (bottom, bottom_range) = section_edge(model, i, from)?;
            let (top, top_range) = section_edge(model, i, to)?;
            let (left, right) = (across_edges[i % boundaries].expect("edge"), across_edges[(i + 1) % boundaries].expect("edge"));
            let uv = |u: f64, v: f64| Point2::new(u, v);
            let (surface, mut parts) = match &analytic[i] {
                Some(analytic) => (analytic.clone(), analytic_parts(model, analytic, [bottom, right, top, left], (v0, v1))?),
                None => (
                    surface,
                    vec![
                        (bottom, Orientation::Forward, segment(uv(0.0, v0), uv(1.0, v0), bottom_range.0, bottom_range.1), bottom_range),
                        (right, Orientation::Forward, segment(uv(1.0, v0), uv(1.0, v1), v0, v1), (v0, v1)),
                        (top, Orientation::Reversed, segment(uv(0.0, v1), uv(1.0, v1), top_range.0, top_range.1), top_range),
                        (left, Orientation::Reversed, segment(uv(0.0, v0), uv(0.0, v1), v0, v1), (v0, v1)),
                    ],
                ),
            };
            model.set_face_surface(face, surface, crate::tolerance::LINEAR)?;
            // A loop turning clockwise in the surface's parameters (an analytic surface facing
            // the other way): the face is used reversed, its loop run backwards.
            let reversed = loop_area(&parts) < 0.0;
            if reversed {
                parts.reverse();
                for part in &mut parts {
                    part.1 = part.1.reversed();
                }
            }
            let mut coedges = Vec::with_capacity(4);
            for (edge, orientation, pcurve, (a, b)) in parts {
                let coedge = model.add_coedge(edge, orientation)?;
                model.set_pcurve_with_range(coedge, pcurve, a, b)?;
                coedges.push(coedge);
            }
            let wire = model.add_wire(&coedges)?;
            model.add(face, wire.into())?;
            let face = if reversed { face.reversed() } else { face };
            faces.push(face);
            columns[i].push(face);
        }
    }
    let shell = model.new_empty(NodeKind::Shell);
    for &face in &faces {
        model.add(shell, face)?;
    }
    if !options.solid {
        return Ok(Built { shape: shell, columns, across });
    }
    // Caps on the end sections, using their edges the other way round.
    for (k, outward) in [(0, Orientation::Reversed), (n - 1, Orientation::Forward)] {
        let Section::Wire(w) = &sections[k] else { continue };
        // The cap's boundary runs counter-clockwise about its plane's normal. (A section not
        // in one plane, as a wire projected onto a curved face: capped by the filling of
        // least bending energy, on its edges.)
        let order: Vec<usize> = if outward == Orientation::Reversed { (0..edges).rev().collect() } else { (0..edges).collect() };
        let plane = match plane_of(w, if outward == Orientation::Reversed { -1.0 } else { 1.0 }) {
            Ok(plane) => plane,
            Err(Error::InvalidInput(m)) if m.contains("not planar") => {
                let mut runs = Vec::with_capacity(edges);
                for &i in &order {
                    let (edge, _) = section_edge(model, i, k)?;
                    let (curve, first, last) = model.edge_curve(Shape::from(edge)).expect("a section edge has a curve");
                    runs.push((edge, outward, curve, first, last));
                }
                let cap = crate::fill::face_on_loop(model, &runs, &Default::default())?;
                model.add(shell, cap)?;
                continue;
            }
            Err(e) => return Err(e),
        };
        let face = model.new_empty(NodeKind::Face);
        model.set_face_surface(face, Surface::Plane(plane), crate::tolerance::LINEAR)?;
        let mut coedges = Vec::with_capacity(edges);
        for i in order {
            let (edge, _) = section_edge(model, i, k)?;
            let coedge = model.add_coedge(edge, outward)?;
            let (curve, first, last) = model.edge_curve(Shape::from(edge)).expect("a section edge has a curve");
            let pcurve = crate::geometry::project::analytic_pcurve(&curve, &Surface::Plane(plane))
                .ok_or_else(|| Error::InvalidInput("a loft's end section is not planar".into()))?;
            model.set_pcurve_with_range(coedge, pcurve, first, last)?;
            coedges.push(coedge);
        }
        let wire = model.add_wire(&coedges)?;
        model.add(face, wire.into())?;
        model.add(shell, face)?;
    }
    let solid = model.new_empty(NodeKind::Solid);
    // (Its sign from the faces' terms, kept in the model for the checks to come.)
    let volume = properties::closed_shell_volume(model, shell);
    model.add(solid, if volume < 0.0 { shell.reversed() } else { shell })?;
    Ok(Built { shape: solid, columns, across })
}

/// The analytic surface ruled between two corresponding pieces, if they make one: lines in
/// one plane, or arcs about one axis at the same angles (a cone, a cylinder, or a flat ring).
fn ruled_surface(a: &Piece, b: &Piece) -> Option<Surface> {
    const TOL: f64 = crate::tolerance::LINEAR;
    let ends = [a.at(0.0), a.at(1.0), b.at(0.0), b.at(1.0)];
    match (&a.curve, &b.curve) {
        (Curve::Line(_), Curve::Line(_)) => {
            let [a0, a1, b0, b1] = ends;
            let span = (a1 - a0).length().max((b1 - b0).length());
            let normal = [(b0 - a0), (b1 - a0)]
                .into_iter()
                .map(|d| (a1 - a0).cross(d))
                .max_by(|x, y| x.length().total_cmp(&y.length()))?;
            if normal.length() <= 1e-9 * span * span {
                return None;
            }
            let plane = Plane::new(a0, Direction::try_from_vector(normal)?);
            ends.iter().all(|&p| plane.distance(p) <= TOL).then_some(Surface::Plane(plane))
        }
        (Curve::Circle(ca), Curve::Circle(cb)) => {
            let axis = |c: &crate::geometry::curve::Circle| Vector::from(c.frame().x()).cross(Vector::from(c.frame().y()));
            let z = axis(ca);
            let offset = cb.center() - ca.center();
            let height = offset.dot(z);
            if axis(cb).cross(z).length() > 1e-9 || (offset - z * height).length() > TOL {
                return None;
            }
            // Corresponding points at the same angles about the axis.
            let radial = |p: Point, c: Point| {
                let d = p - c;
                Direction::try_from_vector(d - z * d.dot(z))
            };
            let aligned = [0.0, 0.5, 1.0].iter().all(|&s| match (radial(a.at(s), ca.center()), radial(b.at(s), cb.center())) {
                (Some(p), Some(q)) => Vector::from(p).dot(Vector::from(q)) > 1.0 - 1e-12,
                _ => false,
            });
            if !aligned {
                return None;
            }
            let x = radial(a.at(0.0), ca.center())?;
            let (ra, rb) = (ca.radius(), cb.radius());
            if height.abs() <= TOL {
                if (ra - rb).abs() <= TOL {
                    return None;
                }
                return Some(Surface::Plane(Plane::new(ca.center(), Direction::try_from_vector(z)?)));
            }
            // The axis towards the second arc.
            let up = Direction::try_from_vector(z * height.signum())?;
            let frame = Frame::new(ca.center(), up, x);
            if (ra - rb).abs() <= TOL {
                return Some(Surface::Cylinder(Cylinder::new(frame, ra)));
            }
            Some(Surface::Cone(Cone::new(frame, (rb - ra).atan2(height.abs()), ra)))
        }
        _ => None,
    }
}

/// The coedges (edge, orientation, pcurve, range) of an analytic face between sections,
/// from its `[bottom, right, top, left]` edges: the section edges' pcurves in closed form,
/// the edges across straight in the surface's parameters.
#[allow(clippy::type_complexity)]
fn analytic_parts(
    model: &Model,
    surface: &Surface,
    [bottom, right, top, left]: [EdgeId; 4],
    (v0, v1): (f64, f64),
) -> Result<Vec<(EdgeId, Orientation, Curve2, (f64, f64))>> {
    let pcurve = |edge: EdgeId| -> Result<(Curve2, f64, f64)> {
        let (curve, first, last) = model.edge_curve(Shape::from(edge)).expect("a section edge has a curve");
        let p = crate::geometry::project::analytic_pcurve(&curve, surface)
            .ok_or_else(|| Error::InvalidInput("no pcurve on a loft's analytic face".into()))?;
        Ok((p, first, last))
    };
    let (pb, fb, lb) = pcurve(bottom)?;
    let (mut pt, ft, lt) = pcurve(top)?;
    let (b0, b1) = (pb.point(fb), pb.point(lb));
    // The top in the same period as the bottom (their ends at the same angles).
    if let Some(period) = surface.u_period() {
        let du = ((b0.x - pt.point(ft).x) / period).round() * period;
        if du != 0.0 {
            pt = pt.transformed(&Transform2::translation(Vector2::new(du, 0.0)));
        }
    }
    let (t0, t1) = (pt.point(ft), pt.point(lt));
    Ok(vec![
        (bottom, Orientation::Forward, pb, (fb, lb)),
        (right, Orientation::Forward, segment(b1, t1, v0, v1), (v0, v1)),
        (top, Orientation::Reversed, pt, (ft, lt)),
        (left, Orientation::Reversed, segment(b0, t0, v0, v1), (v0, v1)),
    ])
}

/// The signed area of a loop of coedges in parameter space (positive counter-clockwise).
fn loop_area(parts: &[(EdgeId, Orientation, Curve2, (f64, f64))]) -> f64 {
    let mut points = Vec::new();
    for (_, orientation, pcurve, (a, b)) in parts {
        let (a, b) = if *orientation == Orientation::Forward { (*a, *b) } else { (*b, *a) };
        for j in 0..16 {
            points.push(pcurve.point(a + (b - a) * j as f64 / 16.0));
        }
    }
    let n = points.len();
    (0..n).map(|i| points[i].x * points[(i + 1) % n].y - points[(i + 1) % n].x * points[i].y).sum::<f64>() / 2.0
}

/// The plane of a closed planar wire, normal along its Newell normal times `sign`. Fails if
/// the wire is not planar.
fn plane_of(w: &WireSection, sign: f64) -> Result<Plane> {
    let n = normal(w) * sign;
    let direction = Direction::try_from_vector(n).ok_or_else(|| Error::InvalidInput("a loft's end section encloses no area".into()))?;
    let points: Vec<Point> = w.pieces.iter().flat_map(|p| (0..8).map(move |i| p.at(i as f64 / 8.0))).collect();
    let center = Point::ORIGIN + points.iter().fold(Vector::ZERO, |acc, p| acc + (*p - Point::ORIGIN)) / points.len() as f64;
    let plane = Plane::new(center, direction);
    let size = points.iter().map(|p| p.distance(center)).fold(0.0, f64::max);
    if points.iter().any(|&p| plane.distance(p) > 1e-7 * size.max(1.0)) {
        return Err(Error::InvalidInput("a loft's end section is not planar".into()));
    }
    Ok(plane)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::curve::Circle;
    use crate::make;
    use crate::math::Frame;

    fn circle(model: &mut Model, z: f64, radius: f64, normal: Direction) -> Shape {
        let frame = Frame::new(Point::new(0.0, 0.0, z), normal, Direction::X);
        let edge = make::edge_from_curve(model, Curve::Circle(Circle::new(frame, radius))).unwrap();
        make::wire(model, &[edge]).unwrap()
    }

    fn polygon(model: &mut Model, points: &[[f64; 3]]) -> Shape {
        let p: Vec<Point> = points.iter().map(|c| Point::new(c[0], c[1], c[2])).collect();
        let edges: Vec<Shape> = (0..p.len()).map(|i| make::edge_between(model, p[i], p[(i + 1) % p.len()]).unwrap()).collect();
        make::wire(model, &edges).unwrap()
    }

    /// A smooth loft through circles of radius 5, 6 and 5 at heights 0, 1 and 10: through
    /// them evenly (as CadQuery lofts), as OCCT's (volume 1011.593); by their distances
    /// apart, it bulges further between the two near ones.
    #[test]
    fn three_sections_spaced_evenly_or_by_distance() {
        let mut model = Model::new();
        let sections: Vec<Shape> = [(0.0, 5.0), (1.0, 6.0), (10.0, 5.0)].iter().map(|&(z, r)| circle(&mut model, z, r, Direction::Z)).collect();
        let volume = |model: &mut Model, spacing: Spacing| {
            let options = LoftOptions { solid: true, spacing, ..Default::default() };
            let solid = loft(model, &sections, &options).unwrap().shape;
            properties::volume(model, solid).mass
        };
        let even = volume(&mut model, Spacing::Even);
        assert!((even - 1011.593).abs() < 1e-3, "{even}");
        assert!(volume(&mut model, Spacing::ByDistance) > even + 100.0);
    }

    /// A loft from a circle up to a rectangle through a smaller circle: the bottom circle's
    /// pieces (split to match the rectangle's sides) stay arcs of it, the rectangle's sides
    /// lines, as the faces beyond (a cylinder the loft stands on) meet them; the faces between
    /// hold to them within the edges' tolerances (the middle section is no edge: the faces run
    /// on through it).
    #[test]
    fn a_loft_keeps_its_sections_lines_and_arcs() {
        let mut model = Model::new();
        let sections = [
            circle(&mut model, 0.0, 3.0, Direction::Z),
            circle(&mut model, 2.0, 2.5, Direction::Z),
            polygon(&mut model, &[[-2.0, -1.0, 4.0], [2.0, -1.0, 4.0], [2.0, 1.0, 4.0], [-2.0, 1.0, 4.0]]),
        ];
        let options = LoftOptions { solid: true, ..Default::default() };
        let solid = loft(&mut model, &sections, &options).unwrap().shape;
        let report = crate::validity::check(&model, solid);
        assert!(report.is_valid(), "{report:?}");
        let mut arcs = 0;
        let mut lines = 0;
        for e in model.unique_subshapes(solid, NodeKind::Edge) {
            let Some((curve, first, last)) = model.edge_curve(e) else { continue };
            let (a, b) = (curve.point(first), curve.point(last));
            if (a.z - b.z).abs() > 1e-9 {
                continue;
            }
            match curve {
                Curve::Circle(_) => arcs += 1,
                Curve::Line(_) => lines += 1,
                _ => panic!("a section edge at z {} made {curve:?}", a.z),
            }
        }
        assert_eq!((arcs, lines), (4, 4));
    }

    /// A strip 0.1 high ruled between a helix winding forty times round a unit cylinder
    /// (pitch 0.4) and the same 0.1 higher: its area is 0.1 times the cylinder's
    /// circumference forty times, its turns told apart.
    #[test]
    fn a_strip_ruled_between_helices_winding_forty_times() {
        use crate::geometry::fit;
        let mut model = Model::new();
        let (turns, pitch) = (40.0, 0.4);
        let helix = |model: &mut Model, dz: f64| {
            let n = 32 * turns as usize;
            let params: Vec<f64> = (0..=n).map(|i| i as f64 / n as f64).collect();
            let points: Vec<Point> = params.iter().map(|&t| {
                let a = std::f64::consts::TAU * turns * t;
                Point::new(a.cos(), a.sin(), pitch * turns * t + dz)
            }).collect();
            let edge = make::edge_from_curve(model, Curve::BSpline(fit::interpolate_with_parameters(&points, &params).unwrap())).unwrap();
            make::wire(model, &[edge]).unwrap()
        };
        let (low, high) = (helix(&mut model, 0.0), helix(&mut model, 0.1));
        let strip = loft(&mut model, &[low, high], &LoftOptions { ruled: true, ..Default::default() }).unwrap().shape;
        let area = properties::area(&model, strip).mass;
        let expected = 0.1 * std::f64::consts::TAU * turns;
        assert!((area - expected).abs() < 1e-4 * expected, "area {area} not {expected}");
    }

    #[test]
    fn circles_running_opposite_ways_make_a_frustum() {
        let mut model = Model::new();
        let bottom = circle(&mut model, 0.0, 10.0, Direction::Z);
        let top = circle(&mut model, 40.0, 6.0, Direction::NEG_Z);
        let shell = loft(&mut model, &[bottom, top], &LoftOptions::default()).unwrap().shape;
        assert_eq!(shell.kind(), NodeKind::Shell);
        let slant = (16.0f64 + 1600.0).sqrt();
        let area = properties::area(&model, shell).mass;
        assert!((area - std::f64::consts::PI * 16.0 * slant).abs() < 1e-3 * area, "area {area}");
        assert!(crate::validity::check(&model, shell).errors().next().is_none());
    }

    #[test]
    fn a_ruled_solid_between_squares_of_different_edge_counts() {
        let mut model = Model::new();
        let bottom = polygon(&mut model, &[[0., 0., 0.], [2., 0., 0.], [2., 2., 0.], [0., 2., 0.]]);
        // The same square at height 1, with its first side in two edges.
        let top = polygon(&mut model, &[[0., 0., 1.], [1., 0., 1.], [2., 0., 1.], [2., 2., 1.], [0., 2., 1.]]);
        let solid = loft(&mut model, &[bottom, top], &LoftOptions { solid: true, ruled: true, ..Default::default() }).unwrap().shape;
        assert_eq!(solid.kind(), NodeKind::Solid);
        assert_eq!(model.unique_subshapes(solid, NodeKind::Face).len(), 5 + 2);
        let volume = properties::volume(&model, solid).mass;
        assert!((volume - 4.0).abs() < 1e-9, "volume {volume}");
        assert!(crate::validity::check(&model, solid).errors().next().is_none());
    }

    #[test]
    fn a_loft_to_a_vertex_is_a_cone() {
        let mut model = Model::new();
        let base = circle(&mut model, 0.0, 10.0, Direction::Z);
        let apex = make::vertex(&mut model, Point::new(0.0, 0.0, 40.0));
        let solid = loft(&mut model, &[base, apex], &LoftOptions { solid: true, ..Default::default() }).unwrap().shape;
        let volume = properties::volume(&model, solid).mass;
        let cone = std::f64::consts::PI * 100.0 * 40.0 / 3.0;
        assert!((volume - cone).abs() < 1e-4 * cone, "volume {volume}");
    }

    #[test]
    fn a_square_matched_to_a_circle_by_rays_keeps_its_history() {
        let mut model = Model::new();
        let square = polygon(&mut model, &[[-10., -10., 0.], [10., -10., 0.], [10., 10., 0.], [-10., 10., 0.]]);
        let circle = circle(&mut model, 20.0, 10.0, Direction::Z);
        let options = LoftOptions { ruled: true, matching: Matching::ByRays, ..Default::default() };
        let lofted = loft(&mut model, &[square, circle], &options).unwrap();
        // The rays from the circle's centre through the square's corners split the circle,
        // and the one through the circle's vertex (at x = 10) the square's side there.
        assert_eq!(model.unique_subshapes(lofted.shape, NodeKind::Face).len(), 5);
        let sides: Vec<usize> = model.unique_subshapes(square, NodeKind::Edge).into_iter().map(|e| lofted.generated(e).len()).collect();
        assert_eq!(sides, [1, 2, 1, 1]);
        let arc = model.unique_subshapes(circle, NodeKind::Edge)[0];
        assert_eq!(lofted.generated(arc).len(), 5);
        for corner in model.unique_subshapes(square, NodeKind::Vertex) {
            let across = lofted.generated(corner);
            assert_eq!(across.len(), 1);
            // From the corner straight up and in to the circle.
            let length = properties::length(&model, across[0]).mass;
            assert!((length - (400.0 + (200f64.sqrt() - 10.0).powi(2)).sqrt()).abs() < 1e-6, "length {length}");
        }
    }

    /// A loft from a square to one with a corner lifted (a section in no one plane): the end
    /// capped by a filling on its edges, the solid closed.
    #[test]
    fn a_section_in_no_plane_is_capped_by_a_filling() {
        let mut model = Model::new();
        let low = polygon(&mut model, &[[0.0, 0.0, 0.0], [10.0, 0.0, 0.0], [10.0, 10.0, 0.0], [0.0, 10.0, 0.0]]);
        let high = polygon(&mut model, &[[0.0, 0.0, 10.0], [10.0, 0.0, 10.0], [10.0, 10.0, 12.0], [0.0, 10.0, 10.0]]);
        let built = loft(&mut model, &[low, high], &LoftOptions { solid: true, ruled: true, ..Default::default() }).unwrap();
        let report = crate::validity::check(&model, built.shape);
        assert!(report.is_valid(), "{report:?}");
        // (Between 1000, the box, and 1000 plus the lifted corner's pyramid over the top.)
        let volume = crate::properties::volume(&model, built.shape).mass;
        assert!(volume > 1000.0 && volume < 1000.0 + 100.0 * 2.0 / 2.0, "{volume}");
    }
}
