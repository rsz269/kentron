//! BRepGraph_Polygon_Test.cxx: meshes stored on faces and coedges.
//!
//! Not ported: TopoDS round trips (`*_Roundtrip_*`), `Polygon3D_Captured_WhenPresent`
//! (compares with the source TopoDS shape) and `VertexPointRepresentations_*` (vertex
//! parameters on curves and surfaces, which kentron does not store).

use kentron::geometry::ParametricCurve;
use std::collections::HashMap;

use kentron::math::Point2;
use kentron::meshing;
use kentron::primitives::Cylinder;
use kentron::topology::{CoEdge, Face, Model, NodeKind};

use crate::graph_support::box_model;

#[test]
fn meshing_a_box_meshes_every_face() {
    let (mut model, solid) = box_model();
    meshing::mesh(&mut model, solid, 0.5).unwrap();
    for (face, _) in model.table::<Face>().iter() {
        let mesh = model.face_mesh(face).expect("face without mesh");
        assert!(!mesh.triangles().is_empty());
    }
}

#[test]
fn meshing_a_box_gives_coedges_mesh_polylines() {
    let (mut model, solid) = box_model();
    meshing::mesh(&mut model, solid, 0.5).unwrap();
    let with_polyline = model
        .table::<CoEdge>()
        .iter()
        .filter(|&(c, _)| model.coedge_mesh_polyline(c).is_some())
        .count();
    assert!(with_polyline > 0);
    for (_, coedge) in model.table::<CoEdge>().iter() {
        assert!(model.get(coedge.face).is_some(), "coedge without face");
    }
}

#[test]
fn box_coedges_have_non_origin_uv_ends() {
    let (model, _) = box_model();
    let non_origin = model
        .table::<CoEdge>()
        .iter()
        .filter(|&(c, _)| {
            let (curve, first, last) = model.effective_pcurve(c).expect("coedge without pcurve");
            [first, last].iter().any(|&t| curve.point(t).distance(Point2::ORIGIN) > kentron::tolerance::LINEAR)
        })
        .count();
    assert!(non_origin > 0);
}

#[test]
fn meshed_cylinder_seam_has_two_mesh_polylines_on_one_face() {
    let mut model = Model::new();
    let cylinder = Cylinder::new(5.0, 10.0).build(&mut model).unwrap();
    meshing::mesh(&mut model, cylinder, 0.1).unwrap();
    let found = model.explore(cylinder, NodeKind::Edge).any(|edge| {
        let mut per_face: HashMap<_, usize> = HashMap::new();
        for (c, coedge) in model.table::<CoEdge>().iter() {
            if coedge.edge == edge.node && model.coedge_mesh_polyline(c).is_some() {
                *per_face.entry(coedge.face).or_default() += 1;
            }
        }
        per_face.values().any(|&n| n >= 2)
    });
    assert!(found, "no seam with two mesh polylines on one face");
}
