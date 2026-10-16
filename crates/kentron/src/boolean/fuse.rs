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
/// done: with solids and other shapes together, or empty shapes, or where faces of shapes
/// (not solids) lying on one another are not split alike there.
pub(super) fn cells(model: &mut Model, shapes: &[Shape]) -> Result<Option<Vec<Cell>>> {
    let n = shapes.len();
    if n < 2 || shapes.iter().any(|&s| super::is_empty(model, s)) {
        return Ok(None);
    }
    let solid = |s: Shape| model.explore(s, NodeKind::Solid).next().is_some();
    if shapes.iter().all(|&s| super::graph::made_of_solids(model, s)) {
        return solid_cells(model, shapes);
    }
    if shapes.iter().any(|&s| solid(s)) {
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

    // The regions of the faces, each inside its own shape and those whose faces it lies on:
    // a region lying on faces of other shapes is one cell with theirs, kept from the first
    // of the shapes (each face split by the other's edges lying on it, alike there).
    let mut arcs_of = Vec::with_capacity(graph.faces.len());
    for fi in 0..graph.faces.len() {
        arcs_of.push(graph.arcs_on_face(fi, &mut pieces, &of_edge, None));
    }
    let mut kept: Vec<Kept> = Vec::new();
    let mut region_marks: Vec<Vec<usize>> = Vec::new();
    let mut shared_regions = false;
    for (fi, arcs) in arcs_of.iter().enumerate() {
        let face = &graph.faces[fi];
        let scale = face.bbox.0.distance(face.bbox.1);
        for region in split::split_face(arcs, &pieces, &face.polylines, scale, false) {
            let p = geom::surface_point(&face.surface, region.interior);
            let mut on: Vec<usize> = graph
                .index()
                .faces_containing_others(face.operand, p)
                .into_iter()
                .filter(|&g| {
                    let other = &graph.faces[g];
                    geom::signed_distance(&other.surface, p).abs() <= on_tol() && other.locate(p, &graph.edges, &graph.vertices) != Location::Outside
                })
                .map(|g| graph.faces[g].operand)
                .collect();
            if !on.is_empty() {
                shared_regions = true;
            }
            if on.iter().any(|&k| k < face.operand) {
                continue;
            }
            on.push(face.operand);
            on.sort_unstable();
            on.dedup();
            region_marks.push(on);
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
    // Where regions lie on one another, each shape's faces made up of the cells inside it,
    // else (faces split differently where they meet) none of them.
    if shared_regions {
        let areas: Vec<f64> = faces.iter().map(|&f| crate::properties::area(model, f).mass).collect();
        for (k, &shape) in shapes.iter().enumerate() {
            if dimension[k] != 2 {
                continue;
            }
            let whole = crate::properties::area(model, shape).mass;
            let parts: f64 = region_marks.iter().zip(&areas).filter(|(m, _)| m.contains(&k)).map(|(_, a)| a).sum();
            if (parts - whole).abs() > 1e-6 * whole.abs() + 1e-12 {
                return Ok(None);
            }
        }
    }
    let mut cells: Vec<Cell> = faces.into_iter().zip(&region_marks).map(|(shape, marks)| Cell { inside: inside(marks), shape }).collect();
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

/// The cells of shapes made of solids: every face split by all the other shapes in one
/// graph, the pieces of space the parts enclose, each inside the shapes a point of it is in
/// (those inside none, enclosed between the shapes, left out). `None` with a shape whose
/// shell is open (a half-space), or where the cells do not make up the shapes (faces left
/// unsplit where they meet).
fn solid_cells(model: &mut Model, shapes: &[Shape]) -> Result<Option<Vec<Cell>>> {
    let n = shapes.len();
    let mut graph = Graph::with_operands(n);
    for (k, &shape) in shapes.iter().enumerate() {
        graph.add_operand(model, shape, k, &Default::default())?;
    }
    if graph.open.iter().any(|&o| o) {
        return Ok(None);
    }
    graph.intersect_edges_with_faces();
    graph.intersect_edges_on_faces();
    graph.share_vertices();
    graph.intersect_faces();
    graph.intersect_sections();
    graph.share_vertices();
    let (mut pieces, of_edge) = graph.pieces();
    let mut arcs_of = Vec::with_capacity(graph.faces.len());
    for fi in 0..graph.faces.len() {
        arcs_of.push(graph.arcs_on_face(fi, &mut pieces, &of_edge, None));
    }
    // Every region of every face (those of faces lying on one another made one when the
    // solids are gathered).
    let mut kept: Vec<Kept> = Vec::new();
    for (fi, arcs) in arcs_of.iter().enumerate() {
        let face = &graph.faces[fi];
        let scale = face.bbox.0.distance(face.bbox.1);
        for region in split::split_face(arcs, &pieces, &face.polylines, scale, false) {
            kept.push(Kept { face: fi, region, flip: false, far: false });
        }
    }
    let mut out = Output::default();
    let faces = out.faces(model, &graph, &pieces, &arcs_of, &kept)?;
    out.set_tolerances(model, &graph, &pieces)?;
    // (Each face held by the box of the face it is a region of.)
    let known: crate::hashing::FastMap<crate::topology::NodeId, crate::math::BoundingBox> =
        faces.iter().zip(&kept).map(|(f, k)| (f.node, crate::math::BoundingBox::new(graph.faces[k.face].bbox.0, graph.faces[k.face].bbox.1))).collect();
    // Where a point lies: in the shapes whose boxes hold it and which contain it.
    let mut boxes = vec![crate::math::BoundingBox::EMPTY; n];
    for face in &graph.faces {
        boxes[face.operand].add_box(&crate::math::BoundingBox::new(face.bbox.0, face.bbox.1));
    }
    // A step behind a face as a shell holds it, from the middle of its region (the face built
    // faces out of its shape: what lies behind it is behind the region where the shell holds
    // it as built).
    let built: crate::hashing::FastMap<crate::topology::NodeId, usize> = faces.iter().enumerate().map(|(i, f)| (f.node, i)).collect();
    let mut extent = crate::math::BoundingBox::EMPTY;
    for b in &boxes {
        extent.add_box(b);
    }
    let step = (1e-5 * (extent.max - extent.min).length()).max(10.0 * geom::on_tol());
    let marks_at = |_: &Model, held: Shape| -> Option<Vec<bool>> {
        let i = *built.get(&held.node)?;
        let face = &graph.faces[kept[i].face];
        let on = geom::surface_point(&face.surface, kept[i].region.interior);
        let behind = held.orientation == faces[i].orientation;
        let p = on + face.outward_normal(on) * if behind { -step } else { step };
        Some((0..n).map(|k| boxes[k].contains(p) && graph.contains(k, graph.faces_of(k), p)).collect())
    };
    let enclosed = super::volume::enclosed_solids_marked(model, &faces, &known, &marks_at)?;
    // The pieces inside some shape are the cells (those inside none, enclosed between the
    // shapes, left out).
    let mut cells: Vec<Cell> = enclosed.into_iter().filter(|(_, inside)| inside.iter().any(|&b| b)).map(|(shape, inside)| Cell { inside, shape }).collect();
    // Each shape made up of its cells, else (faces left unsplit where shapes meet) none of
    // them.
    let volume = |model: &Model, s: Shape| crate::properties::solid_volume(model, s).map_or_else(|| crate::properties::volume(model, s).mass, |v| v.0);
    let slack = 1e-9 * (extent.max - extent.min).length().powi(3);
    let cell_volumes: Vec<f64> = cells.iter().map(|c| volume(model, c.shape)).collect();
    for (k, &shape) in shapes.iter().enumerate() {
        let whole = volume(model, shape);
        let parts: f64 = cells.iter().zip(&cell_volumes).filter(|(c, _)| c.inside[k]).map(|(_, v)| v).sum();
        // (Within what integrating the faces either way may differ by.)
        if (parts - whole).abs() > 1e-4 * whole.abs() + slack || cell_volumes.iter().any(|&v| v <= 0.0) {
            return Ok(None);
        }
    }
    cells.sort_by(|a, b| b.inside.cmp(&a.inside));
    Ok(Some(cells))
}
