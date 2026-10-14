//! Boolean operations and sections on shapes other than solids: faces and shells, wires and
//! edges, vertices, as OCCT defines them.
//!
//! An operand's dimension is that of its highest-dimensional parts (3 with solids, 2 with
//! faces, 1 with edges, 0 with vertices). The operands are split by each other as for
//! solids, then each part (a region of a face, a piece of an edge, a vertex) is placed with
//! respect to the other operand: inside a solid of it, on it, or outside it.
//!
//! * An intersection keeps the parts of the lower-dimensional operand inside or on the other
//!   (of the object if both have the same dimension).
//! * A difference keeps the parts of the object outside the tool; it is defined when the tool's
//!   dimension is at least the object's.
//! * A union keeps all parts, those of the tool on the object once; it is defined for
//!   operands of the same dimension.
//! * A section is the edges where the operands meet: curves where their faces cross, and
//!   pieces of edges of either lying on the other.
//!
//! Faces of a shell argument stay together in shells, edges of a wire argument in wires.

use crate::hashing::{FastMap as HashMap, FastSet as HashSet};

use crate::error::{Error, Result};
use crate::math::Point;
use crate::topology::{EdgeId, Model, NodeKind, Orientation, Shape};

use super::build::{Kept, Output};
use super::geom::{self, on_tol};
use super::graph::{Graph, Location};
use super::split::{self, ArcOnFace, Piece, Region};
use super::Operation;

/// Where a part of one operand lies with respect to the other.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Place {
    /// Inside one of its solids.
    In,
    /// On one of its faces, edges or vertices.
    On,
    Out,
}

impl Graph {
    /// True if both operands are solids and nothing else.
    pub fn only_solids(&self) -> bool {
        self.dimension == [3, 3]
            && self.faces.iter().all(|f| f.solid)
            && self.edges.iter().all(|e| !e.free)
            && self.free_vertices.is_empty()
    }

    /// Where `p` lies with respect to operand `k`.
    fn place(&self, p: Point, k: usize) -> Place {
        for face in self.index().faces_containing(k, p).into_iter().map(|g| &self.faces[g]) {
            if geom::signed_distance(&face.surface, p).abs() <= on_tol()
                && face.locate(p, &self.edges, &self.vertices) != Location::Outside
            {
                return Place::On;
            }
        }
        for e in self.index().edges_containing(k, p).into_iter().map(|i| &self.edges[i]).filter(|e| e.free) {
            if let Some(curve) = &e.curve {
                if geom::project(curve, e.first, e.last, p).1 <= on_tol() {
                    return Place::On;
                }
            }
        }
        if self.free_vertices.iter().any(|&(o, v)| o == k && self.vertices.points[v].distance(p) <= on_tol()) {
            return Place::On;
        }
        let solid_faces: Vec<usize> = (0..self.faces.len()).filter(|&i| self.faces[i].operand == k && self.faces[i].solid).collect();
        let inside = !solid_faces.is_empty()
            && if self.open[k] { self.behind(&solid_faces, p) } else { self.contains(k, &solid_faces, p) };
        if inside { Place::In } else { Place::Out }
    }

    /// The edges of operand `k` passing through `p`.
    fn edges_along(&self, p: Point, k: usize) -> Vec<usize> {
        self.index()
            .edges_containing(k, p)
            .into_iter()
            .filter(|&i| {
                let e = &self.edges[i];
                e.curve.as_ref().is_some_and(|c| geom::project(c, e.first, e.last, p).1 <= on_tol())
            })
            .collect()
    }

    /// A point in the middle of a piece.
    fn piece_middle(&self, piece: &Piece) -> Point {
        self.edges[piece.edge].point(0.5 * (piece.first + piece.last), &self.vertices)
    }
}

thread_local! {
    /// Whether faces keep the arcs splitting them that bound nothing, as slits.
    static KEEP_SLITS: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
}

/// `f` with the faces of operations on shapes other than solids keeping the arcs on them
/// that bound nothing as slits (loops along the arc and back), for later operations to
/// split there: shapes fused one after another as if all at once.
pub(super) fn with_slits<T>(f: impl FnOnce() -> T) -> T {
    let before = KEEP_SLITS.with(|k| k.replace(true));
    let result = f();
    KEEP_SLITS.with(|k| k.set(before));
    result
}

/// Splits both operands by each other.
fn intersect(graph: &mut Graph) -> (Vec<Piece>, Vec<Vec<usize>>) {
    graph.intersect_edges_with_faces();
    graph.intersect_free_edges();
    graph.intersect_edges_on_faces();
    graph.share_vertices();
    graph.intersect_faces();
    graph.pieces()
}

/// The boolean operation on operands not both solids (already added to `graph`).
pub(super) fn run(model: &mut Model, mut graph: Graph, operation: Operation) -> Result<Shape> {
    let dims = graph.dimension;
    match operation {
        Operation::Union if dims[0] != dims[1] => {
            return Err(Error::InvalidInput(format!(
                "the union of a {}-dimensional and a {}-dimensional shape is not defined",
                dims[0], dims[1]
            )))
        }
        Operation::Difference if dims[1] < dims[0] => {
            return Err(Error::InvalidInput(format!(
                "the difference of a {}-dimensional shape and a {}-dimensional one is not defined",
                dims[0], dims[1]
            )))
        }
        _ => {}
    }
    let lowest = dims[0].min(dims[1]);
    // `shared`: of a wire or shell of the tool that a union keeps whole (see below).
    let keeps = |k: usize, place: Place, shared: bool| match operation {
        Operation::Intersection if dims[0] == dims[1] => k == 0 && place != Place::Out,
        Operation::Intersection => dims[k] == lowest && place != Place::Out,
        Operation::Difference => k == 0 && place == Place::Out,
        Operation::Union => k == 0 || place == Place::Out || shared,
    };
    // Parts of solids are kept only by operations on solids alone.
    for k in 0..2 {
        let wanted = match operation {
            Operation::Intersection => dims[k] == lowest,
            Operation::Difference => k == 0,
            Operation::Union => true,
        };
        if wanted && graph.faces.iter().any(|f| f.operand == k && f.solid) {
            todo!("boolean operations on solids mixed with other shapes");
        }
    }
    let (mut pieces, of_edge) = intersect(&mut graph);

    // Where the parts lie: regions of faces, pieces of edges.
    let mut arcs_of: Vec<Vec<ArcOnFace>> = Vec::with_capacity(graph.faces.len());
    for fi in 0..graph.faces.len() {
        let free = !graph.faces[fi].solid;
        arcs_of.push(if free { graph.arcs_on_face(fi, &mut pieces, &of_edge, None) } else { Vec::new() });
    }
    let mut regions: Vec<(usize, Region, Place)> = Vec::new();
    for (fi, arcs) in arcs_of.iter().enumerate() {
        let face = &graph.faces[fi];
        if face.solid {
            continue;
        }
        let scale = face.bbox.0.distance(face.bbox.1);
        for region in split::split_face(arcs, &pieces, &face.polylines, scale, KEEP_SLITS.with(|k| k.get())) {
            let p = geom::surface_point(&face.surface, region.interior);
            let place = graph.place(p, 1 - face.operand);
            regions.push((fi, region, place));
        }
    }
    let mut piece_places: HashMap<usize, Place> = HashMap::default();
    for (ei, e) in graph.edges.iter().enumerate() {
        if let (Some(k), true) = (e.operand, e.free) {
            for &pi in &of_edge[ei] {
                piece_places.insert(pi, graph.place(graph.piece_middle(&pieces[pi]), 1 - k));
            }
        }
    }
    // A union keeps each wire or shell of the tool whole, its parts on the object shared
    // with it, unless it lies entirely on the object (as OCCT, which avoids duplicating
    // containers).
    let mut partly_out: HashSet<usize> = HashSet::default();
    for (fi, _, place) in &regions {
        let face = &graph.faces[*fi];
        if let (1, Place::Out, Some(c)) = (face.operand, place, face.container) {
            partly_out.insert(c);
        }
    }
    for (&pi, place) in &piece_places {
        let e = &graph.edges[pieces[pi].edge];
        if let (Some(1), Place::Out, Some(c)) = (e.operand, place, e.container) {
            partly_out.insert(c);
        }
    }
    let shared = |k: usize, container: Option<usize>| k == 1 && container.is_some_and(|c| partly_out.contains(&c));
    // A wire or shell of the object lying entirely on the tool is avoided too when the tool's
    // are kept whole (they hold its parts).
    let mut partly_off: HashSet<usize> = HashSet::default();
    for (fi, _, place) in &regions {
        let face = &graph.faces[*fi];
        if let (0, Some(c)) = (face.operand, face.container) {
            if *place != Place::On {
                partly_off.insert(c);
            }
        }
    }
    for (&pi, place) in &piece_places {
        let e = &graph.edges[pieces[pi].edge];
        if let (Some(0), Some(c)) = (e.operand, e.container) {
            if *place != Place::On {
                partly_off.insert(c);
            }
        }
    }
    let avoided = |k: usize, container: Option<usize>| {
        operation == Operation::Union && k == 0 && !partly_out.is_empty() && container.is_some_and(|c| !partly_off.contains(&c))
    };

    let mut kept: Vec<Kept> = Vec::new();
    for (fi, region, place) in regions {
        let face = &graph.faces[fi];
        if keeps(face.operand, place, shared(face.operand, face.container)) && !avoided(face.operand, face.container) {
            kept.push(Kept { face: fi, region, flip: false, far: false });
        }
    }
    // Pieces of edges. Along an edge of the other operand, an edge is split only
    // at the operands' own vertices: consecutive pieces running along the same edge are
    // merged across other points (where the curves cross while within tolerance of each
    // other).
    let original: HashSet<usize> = graph.edges.iter().filter(|e| e.operand.is_some()).flat_map(|e| [e.start, e.end]).collect();
    let mut kept_pieces: Vec<usize> = Vec::new();
    for (ei, e) in graph.edges.iter().enumerate() {
        let (Some(k), true) = (e.operand, e.free) else { continue };
        let mut previous_on: Option<usize> = None;
        let mut previous_along: Vec<usize> = Vec::new();
        for &pi in &of_edge[ei] {
            let place = piece_places[&pi];
            if !keeps(k, place, shared(k, e.container)) || avoided(k, e.container) {
                previous_on = None;
                previous_along.clear();
                continue;
            }
            // (Not where the piece is shared with a wire of the tool kept whole.)
            let merging = place == Place::On
                && (operation == Operation::Intersection || (operation == Operation::Union && partly_out.is_empty()));
            let along = if merging { graph.edges_along(graph.piece_middle(&pieces[pi]), 1 - k) } else { Vec::new() };
            match previous_on {
                Some(last)
                    if merging
                        && !original.contains(&pieces[pi].start)
                        && along.iter().any(|e| previous_along.contains(e)) =>
                {
                    let (a, b) = (&pieces[kept_pieces[last]], &pieces[pi]);
                    let merged = Piece { edge: ei, first: a.first, last: b.last, start: a.start, end: b.end, rep: pieces.len() };
                    pieces.push(merged);
                    kept_pieces[last] = pieces.len() - 1;
                }
                _ => {
                    kept_pieces.push(pi);
                    previous_on = merging.then_some(kept_pieces.len() - 1);
                }
            }
            previous_along = along;
        }
    }
    // Vertices.
    let kept_vertices: Vec<usize> = graph
        .free_vertices
        .iter()
        .filter(|&&(k, v)| keeps(k, graph.place(graph.vertices.points[v], 1 - k), false))
        .map(|&(_, v)| v)
        .collect();

    let mut out = Output::default();
    let result = model.new_empty(NodeKind::Compound);
    let faces = out.faces(model, &graph, &pieces, &arcs_of, &kept)?;
    let containers: Vec<Option<usize>> = kept.iter().map(|k| graph.faces[k.face].container).collect();
    for group in groups(faces.len(), &out.neighbours(), |i, j| containers[i].is_some() && containers[i] == containers[j]) {
        if containers[group[0]].is_none() {
            model.add(result, faces[group[0]])?;
            continue;
        }
        let shell = model.new_empty(NodeKind::Shell);
        for i in group {
            model.add(shell, faces[i])?;
        }
        model.add(result, shell)?;
    }
    let edges: Vec<Shape> = kept_pieces.iter().map(|&pi| out.edge(model, &graph, &pieces, pi)).collect::<Result<_>>()?;
    let wire_of: Vec<Option<usize>> = kept_pieces.iter().map(|&pi| graph.edges[pieces[pi].edge].container).collect();
    let mut touching = Vec::new();
    for i in 0..kept_pieces.len() {
        for j in i + 1..kept_pieces.len() {
            let (a, b) = (&pieces[kept_pieces[i]], &pieces[kept_pieces[j]]);
            if [a.start, a.end].iter().any(|v| *v == b.start || *v == b.end) {
                touching.push((i, j));
            }
        }
    }
    for group in groups(edges.len(), &touching, |i, j| wire_of[i].is_some() && wire_of[i] == wire_of[j]) {
        if wire_of[group[0]].is_none() {
            model.add(result, edges[group[0]])?;
            continue;
        }
        // Each edge used as its piece runs in the input wire.
        let coedges = group
            .iter()
            .map(|&i| {
                let id = EdgeId::try_from(edges[i].node).expect("an edge");
                let reversed = graph.edges[pieces[kept_pieces[i]].edge].reversed;
                let orientation = if reversed { edges[i].orientation.reversed() } else { edges[i].orientation };
                model.add_coedge(id, orientation)
            })
            .collect::<Result<Vec<_>>>()?;
        let wire = model.add_wire(&coedges)?;
        model.set_coedge_order(wire, &coedges)?;
        model.add(result, wire.into())?;
    }
    for v in kept_vertices {
        let id = out.vertex(model, &graph, v);
        model.add(result, id.into())?;
    }
    out.set_tolerances(model, &graph, &pieces)?;
    Ok(result)
}

/// The section of two operands (already added to `graph`): a compound of the edges where
/// they meet, and of the points where they touch elsewhere.
pub(super) fn section(model: &mut Model, mut graph: Graph) -> Result<Shape> {
    let (pieces, of_edge) = intersect(&mut graph);
    let mut taken: HashSet<usize> = HashSet::default();
    let mut chosen: Vec<usize> = Vec::new();
    for (ei, e) in graph.edges.iter().enumerate() {
        for &pi in &of_edge[ei] {
            let piece = &pieces[pi];
            let meets = match e.operand {
                None => true,
                Some(k) => e.curve.is_some() && graph.place(graph.piece_middle(piece), 1 - k) == Place::On,
            };
            if meets && taken.insert(piece.rep) {
                chosen.push(pi);
            }
        }
    }
    let mut out = Output::default();
    let result = model.new_empty(NodeKind::Compound);
    let mut ends: HashSet<usize> = HashSet::default();
    for pi in chosen {
        let edge = out.edge(model, &graph, &pieces, pi)?;
        model.add(result, Shape { orientation: Orientation::Forward, ..edge })?;
        ends.extend([pieces[pi].start, pieces[pi].end]);
    }
    // Vertices of either operand on the other, away from the section's edges.
    let mut vertices: Vec<(usize, usize)> = graph.free_vertices.clone();
    for e in &graph.edges {
        if let Some(k) = e.operand {
            vertices.extend([(k, e.start), (k, e.end)]);
        }
    }
    for &(k, v) in &vertices {
        if !ends.contains(&v) && graph.place(graph.vertices.points[v], 1 - k) == Place::On {
            ends.insert(v);
            let id = out.vertex(model, &graph, v);
            model.add(result, id.into())?;
        }
    }
    out.set_tolerances(model, &graph, &pieces)?;
    Ok(result)
}

/// The groups of `0..count` connected by `links` where `joins` allows, each in increasing
/// order, ordered by their first member.
fn groups(count: usize, links: &[(usize, usize)], joins: impl Fn(usize, usize) -> bool) -> Vec<Vec<usize>> {
    let mut parent: Vec<usize> = (0..count).collect();
    fn root(parent: &mut [usize], mut i: usize) -> usize {
        while parent[i] != i {
            parent[i] = parent[parent[i]];
            i = parent[i];
        }
        i
    }
    for &(i, j) in links {
        if joins(i, j) {
            let (a, b) = (root(&mut parent, i), root(&mut parent, j));
            parent[a.max(b)] = a.min(b);
        }
    }
    let mut out: Vec<Vec<usize>> = Vec::new();
    let mut index_of: crate::hashing::FastMap<usize, usize> = crate::hashing::FastMap::default();
    for i in 0..count {
        let r = root(&mut parent, i);
        let g = *index_of.entry(r).or_insert_with(|| {
            out.push(Vec::new());
            out.len() - 1
        });
        out[g].push(i);
    }
    out
}
