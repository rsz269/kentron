//! A general fuse of shapes without solids in one graph: every shape split by all the others
//! at once (where [`super::cells`] otherwise splits them one after another, each step an
//! operation of its own).
//!
//! Each part of a shape (a region of a face, a piece of an edge, a vertex) is a cell, inside
//! its own shape and those of the others it lies on, of at least its dimension. Faces are
//! split where they meet faces, not by the edges of other shapes lying on them; parts of
//! edges and vertices that coincide are one cell.

use crate::hashing::FastMap as HashMap;

use crate::error::Result;
use crate::topology::{Model, NodeKind, Shape};

use super::build::{Kept, Output};
use super::cells::Cell;
use super::geom::{self, on_tol};
use super::graph::{Graph, Location};
use super::split;

/// The cells of `shapes` (as [`super::cells`]), found in one graph; `None` where that is not
/// done: with solids, empty shapes, or faces of two shapes lying on one another.
pub(super) fn cells(model: &mut Model, shapes: &[Shape]) -> Result<Option<Vec<Cell>>> {
    let n = shapes.len();
    if n < 2 || shapes.iter().any(|&s| model.explore(s, NodeKind::Solid).next().is_some() || super::is_empty(model, s)) {
        return Ok(None);
    }
    let mut graph = Graph::with_operands(n);
    for (k, &shape) in shapes.iter().enumerate() {
        graph.add_operand(model, shape, k, &Default::default())?;
    }
    graph.intersect_edges_with_faces();
    graph.intersect_free_edges();
    graph.intersect_edges_on_faces();
    graph.share_vertices();
    graph.intersect_faces();
    graph.intersect_sections();
    graph.share_vertices();
    let (mut pieces, of_edge) = graph.pieces();
    // (Faces are not split by other shapes' edges lying on them.)
    let free: Vec<bool> = graph.edges.iter().map(|e| e.free).collect();
    for face in &mut graph.faces {
        face.on_edges.retain(|&e| !free[e]);
    }
    let dimension = graph.dimension.clone();

    // The regions of the faces, each inside its own shape only (on another's face, the
    // cells of both would be one: not done here).
    let mut arcs_of = Vec::with_capacity(graph.faces.len());
    for fi in 0..graph.faces.len() {
        arcs_of.push(graph.arcs_on_face(fi, &mut pieces, &of_edge, None));
    }
    let mut kept: Vec<Kept> = Vec::new();
    let mut region_marks: Vec<usize> = Vec::new();
    for (fi, arcs) in arcs_of.iter().enumerate() {
        let face = &graph.faces[fi];
        let scale = face.bbox.0.distance(face.bbox.1);
        for region in split::split_face(arcs, &pieces, &face.polylines, scale, false) {
            let p = geom::surface_point(&face.surface, region.interior);
            let on_other = graph.index().faces_containing_others(face.operand, p).into_iter().any(|g| {
                let other = &graph.faces[g];
                geom::signed_distance(&other.surface, p).abs() <= on_tol() && other.locate(p, &graph.edges, &graph.vertices) != Location::Outside
            });
            if on_other {
                return Ok(None);
            }
            region_marks.push(face.operand);
            kept.push(Kept { face: fi, region, flip: false, far: false });
        }
    }

    // The shapes (other than `own`) a point of a part of dimension `d` lies on.
    let marks_at = |graph: &Graph, p: crate::math::Point, own: usize, d: u8| -> Vec<usize> {
        let mut candidates: Vec<usize> = graph.index().faces_containing_others(own, p).into_iter().map(|g| graph.faces[g].operand).collect();
        candidates.extend(graph.index().edges_containing_others(own, p).into_iter().filter_map(|e| graph.edges[e].operand));
        candidates.extend(graph.free_vertices.iter().filter(|&&(k, v)| k != own && graph.vertices.points[v].distance(p) <= on_tol()).map(|&(k, _)| k));
        candidates.sort_unstable();
        candidates.dedup();
        candidates.into_iter().filter(|&j| dimension[j] >= d && graph.is_on(p, j)).collect()
    };

    // Pieces of edges not on faces: coinciding ones (of shapes lying along one another)
    // one cell, inside all of them.
    let mut piece_cells: Vec<(usize, Vec<usize>)> = Vec::new();
    let mut by_rep: HashMap<usize, usize> = HashMap::default();
    for (ei, e) in graph.edges.iter().enumerate() {
        let (Some(k), true) = (e.operand, e.free) else { continue };
        for &pi in &of_edge[ei] {
            let mut marks = marks_at(&graph, graph.piece_point(&pieces[pi]), k, 1);
            marks.push(k);
            match by_rep.get(&pieces[pi].rep) {
                Some(&c) => piece_cells[c].1.extend(marks),
                None => {
                    by_rep.insert(pieces[pi].rep, piece_cells.len());
                    piece_cells.push((pi, marks));
                }
            }
        }
    }
    // Vertices, likewise.
    let mut vertex_cells: Vec<(usize, Vec<usize>)> = Vec::new();
    let mut by_vertex: HashMap<usize, usize> = HashMap::default();
    for &(k, v) in &graph.free_vertices {
        let mut marks = marks_at(&graph, graph.vertices.points[v], k, 0);
        marks.push(k);
        match by_vertex.get(&v) {
            Some(&c) => vertex_cells[c].1.extend(marks),
            None => {
                by_vertex.insert(v, vertex_cells.len());
                vertex_cells.push((v, marks));
            }
        }
    }

    let mut out = Output::default();
    let faces = out.faces(model, &graph, &pieces, &arcs_of, &kept)?;
    let inside = |marks: &[usize]| {
        let mut flags = vec![false; n];
        for &k in marks {
            flags[k] = true;
        }
        flags
    };
    let mut cells: Vec<Cell> = faces.into_iter().zip(&region_marks).map(|(shape, &k)| Cell { inside: inside(&[k]), shape }).collect();
    for (pi, marks) in piece_cells {
        let shape = out.edge(model, &graph, &pieces, pi)?;
        cells.push(Cell { inside: inside(&marks), shape });
    }
    for (v, marks) in vertex_cells {
        let shape = Shape::from(out.vertex(model, &graph, v));
        cells.push(Cell { inside: inside(&marks), shape });
    }
    out.set_tolerances(model, &graph, &pieces)?;
    // In the order of shapes added one after another: by the shapes each is inside, the
    // first first, inside one before outside it.
    cells.sort_by(|a, b| b.inside.cmp(&a.inside));
    Ok(Some(cells))
}
