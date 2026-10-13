//! TopExp_Test.cxx

use kentron::topology::{Model, NodeKind};

use crate::{first, make_box};

#[test]
fn map_shapes_faces() {
    let mut model = Model::new();
    let solid = make_box(&mut model, 10.0, 20.0, 30.0);
    assert_eq!(model.unique_subshapes(solid, NodeKind::Face).len(), 6);
}

#[test]
fn map_shapes_edges() {
    let mut model = Model::new();
    let solid = make_box(&mut model, 10.0, 20.0, 30.0);
    assert_eq!(model.unique_subshapes(solid, NodeKind::Edge).len(), 12);
}

#[test]
fn map_shapes_vertices() {
    let mut model = Model::new();
    let solid = make_box(&mut model, 10.0, 20.0, 30.0);
    assert_eq!(model.unique_subshapes(solid, NodeKind::Vertex).len(), 8);
}

#[test]
fn first_vertex_last_vertex() {
    let mut model = Model::new();
    let solid = make_box(&mut model, 10.0, 20.0, 30.0);
    let edge = first(&model, solid, NodeKind::Edge);
    let (start, end) = model.edge_vertices(edge);
    assert!(!start.is_same(&end), "a box edge has distinct end vertices");
}

#[test]
fn vertices() {
    let mut model = Model::new();
    let solid = make_box(&mut model, 10.0, 20.0, 30.0);
    let edge = first(&model, solid, NodeKind::Edge);
    let (start, end) = model.edge_vertices(edge);
    assert_eq!(start.kind(), NodeKind::Vertex);
    assert_eq!(end.kind(), NodeKind::Vertex);
}

#[test]
fn common_vertex() {
    let mut model = Model::new();
    let solid = make_box(&mut model, 10.0, 20.0, 30.0);
    let edges = model.unique_subshapes(solid, NodeKind::Edge);
    assert!(edges.len() >= 2);

    let found = edges.iter().enumerate().any(|(i, e1)| {
        edges[i + 1..].iter().any(|e2| model.common_vertex(*e1, *e2).is_some())
    });
    assert!(found, "a box has adjacent edges");
}

/// Issue #1194: exploring the edges of each face, inside a face exploration, terminates.
#[test]
fn explorer_nested_face_edge_terminates() {
    let mut model = Model::new();
    let solid = make_box(&mut model, 10.0, 20.0, 30.0);

    let mut total_edges = 0;
    let mut faces = 0;
    for face in model.explore(solid, NodeKind::Face) {
        faces += 1;
        let edges = model.explore(face, NodeKind::Edge).take(100).count();
        assert_eq!(edges, 4, "each box face has 4 edges");
        total_edges += edges;
    }
    assert_eq!(faces, 6);
    assert_eq!(total_edges, 24, "6 faces * 4 edges");
}

/// Reusing an explorer for each face (here: a fresh one per face) terminates.
#[test]
fn explorer_re_init_terminates() {
    let mut model = Model::new();
    let solid = make_box(&mut model, 10.0, 20.0, 30.0);
    for face in model.explore(solid, NodeKind::Face) {
        assert_eq!(model.explore(face, NodeKind::Edge).take(100).count(), 4);
    }
}

#[test]
fn map_shapes_and_ancestors() {
    let mut model = Model::new();
    let solid = make_box(&mut model, 10.0, 20.0, 30.0);
    let edge_faces = model.subshape_ancestors(solid, NodeKind::Edge, NodeKind::Face);
    assert_eq!(edge_faces.len(), 12);
    for (_, faces) in &edge_faces {
        assert_eq!(faces.len(), 2, "each box edge is shared by 2 faces");
    }
}
