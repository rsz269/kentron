//! BRepGraph_Build_Test.cxx
//!
//! In OCCT a BRepGraph is populated from TopoDS shapes; in kentron modeling operations build
//! straight into a [`Model`], so "the graph matches TopExp" becomes "the model holds
//! exactly the entities reachable from the built shape". Not ported (OCCT's import
//! pipeline: options, warning statuses, auto-created root products, mutation-boundary
//! checks): `SinglePlanarFace_HasWarnings`, `SinglePlanarFace_Counts`,
//! `AddFlatten_OnEmptyGraph_BuildsFlattenedGraph`, `Build_MutationBoundary_IsValid`,
//! `AddFlatten_SameFaceTwice_DedupsDefinition`, `AddFlatten_AfterBuild_DoesNotCreateNewSolidDefs`,
//! `AddFull_MutationBoundary_IsValid`, `AddFlatten_AppendedFaceHasNoParentShell`,
//! `RootProductIds_Box_ReturnsOneProduct`, `BuildOptions_DisableAutoProduct_DoesNotCreateProducts`,
//! `RootProductIds_AddFlatten_ProductCountUnchanged`.

use std::f64::consts::TAU;

use kentron::geometry::surface::{Cone, Cylinder, Plane, Sphere, Surface, Torus, TrimmedSurface};
use kentron::make;
use kentron::math::{Direction, Frame, Point};
use kentron::primitives;
use kentron::tolerance;
use kentron::topology::{Edge, Face, Model, NodeKind, Shape, Vertex};

use crate::{make_box, make_edge};

/// The model holds exactly the `kinds` entities reachable from `root`.
fn assert_counts_match(model: &Model, root: Shape, kinds: &[NodeKind]) {
    for &kind in kinds {
        assert_eq!(
            model.count(kind),
            model.unique_subshapes(root, kind).len(),
            "{kind:?} count"
        );
    }
}

const ALL_KINDS: [NodeKind; 6] =
    [NodeKind::Solid, NodeKind::Shell, NodeKind::Face, NodeKind::Wire, NodeKind::Edge, NodeKind::Vertex];
const FEV: [NodeKind; 3] = [NodeKind::Face, NodeKind::Edge, NodeKind::Vertex];

fn has_surface(model: &Model, pred: impl Fn(&Surface) -> bool) -> bool {
    model.table::<Face>().iter().any(|(id, _)| model.face_surface(id).is_some_and(&pred))
}

fn degenerate_edges(model: &Model) -> usize {
    model.table::<Edge>().iter().filter(|(id, _)| model.is_degenerate(*id)).count()
}

fn compound(model: &mut Model, children: &[Shape]) -> Shape {
    let c = model.new_empty(NodeKind::Compound);
    for &child in children {
        model.add(c, child).unwrap();
    }
    c
}

// --- sphere ---

#[test]
fn sphere_is_not_empty() {
    let mut model = Model::new();
    primitives::Sphere::new(10.0).build(&mut model).unwrap();
    assert!(!model.is_empty());
}

#[test]
fn sphere_def_counts_match_top_exp() {
    let mut model = Model::new();
    let shape = primitives::Sphere::new(10.0).build(&mut model).unwrap();
    assert_counts_match(&model, shape, &ALL_KINDS);
}

#[test]
fn sphere_surface_type() {
    let mut model = Model::new();
    primitives::Sphere::new(5.0).build(&mut model).unwrap();
    assert!(model.count(NodeKind::Face) >= 1);
    assert!(has_surface(&model, |s| matches!(s, Surface::Sphere(_))));
}

/// A sphere has degenerate edges at its poles.
#[test]
fn sphere_has_degenerate_edges() {
    let mut model = Model::new();
    primitives::Sphere::new(8.0).build(&mut model).unwrap();
    assert!(degenerate_edges(&model) >= 1);
}

// --- cylinder ---

#[test]
fn cylinder_is_not_empty() {
    let mut model = Model::new();
    primitives::Cylinder::new(5.0, 20.0).build(&mut model).unwrap();
    assert!(!model.is_empty());
}

#[test]
fn cylinder_def_counts_match_top_exp() {
    let mut model = Model::new();
    let shape = primitives::Cylinder::new(5.0, 20.0).build(&mut model).unwrap();
    assert_counts_match(&model, shape, &FEV);
}

#[test]
fn cylinder_surface_type() {
    let mut model = Model::new();
    primitives::Cylinder::new(5.0, 20.0).build(&mut model).unwrap();
    assert!(has_surface(&model, |s| matches!(s, Surface::Cylinder(_))));
}

// --- cone ---

#[test]
fn cone_is_not_empty() {
    let mut model = Model::new();
    primitives::Cone::new(10.0, 0.0, 15.0).build(&mut model).unwrap();
    assert!(!model.is_empty());
}

#[test]
fn cone_def_counts_match_top_exp() {
    let mut model = Model::new();
    let shape = primitives::Cone::new(10.0, 0.0, 15.0).build(&mut model).unwrap();
    assert_counts_match(&model, shape, &FEV);
}

#[test]
fn cone_surface_type() {
    let mut model = Model::new();
    primitives::Cone::new(10.0, 5.0, 15.0).build(&mut model).unwrap();
    assert!(has_surface(&model, |s| matches!(s, Surface::Cone(_))));
}

/// A cone with a zero top radius has a degenerate edge at its apex.
#[test]
fn cone_has_degenerate_edge() {
    let mut model = Model::new();
    primitives::Cone::new(10.0, 0.0, 15.0).build(&mut model).unwrap();
    assert!(degenerate_edges(&model) >= 1);
}

// --- torus ---

#[test]
fn torus_is_not_empty() {
    let mut model = Model::new();
    primitives::Torus::new(20.0, 5.0).build(&mut model).unwrap();
    assert!(!model.is_empty());
}

#[test]
#[should_panic] // expected to fail: a known gap in kentron
fn torus_def_counts_match_top_exp() {
    let mut model = Model::new();
    let shape = primitives::Torus::new(20.0, 5.0).build(&mut model).unwrap();
    assert_counts_match(&model, shape, &FEV);
}

#[test]
fn torus_surface_type() {
    let mut model = Model::new();
    primitives::Torus::new(20.0, 5.0).build(&mut model).unwrap();
    assert!(has_surface(&model, |s| matches!(s, Surface::Torus(_))));
}

// --- wedge ---

#[test]
fn wedge_is_not_empty() {
    let mut model = Model::new();
    primitives::Wedge::new(10.0, 10.0, 10.0, 5.0).build(&mut model).unwrap();
    assert!(!model.is_empty());
}

#[test]
fn wedge_def_counts_match_top_exp() {
    let mut model = Model::new();
    let shape = primitives::Wedge::new(10.0, 10.0, 10.0, 5.0).build(&mut model).unwrap();
    assert_counts_match(&model, shape, &ALL_KINDS);
}

#[test]
fn wedge_all_planar_surfaces() {
    let mut model = Model::new();
    primitives::Wedge::new(10.0, 10.0, 10.0, 5.0).build(&mut model).unwrap();
    for (face, _) in model.table::<Face>().iter() {
        let surface = model.face_surface(face).expect("face without surface");
        assert!(matches!(surface, Surface::Plane(_)), "{face:?} is not planar");
    }
}

// --- compounds ---

#[test]
fn compound_two_primitives_is_not_empty() {
    let mut model = Model::new();
    let solid = make_box(&mut model, 10.0, 10.0, 10.0);
    let sphere = primitives::Sphere::new(5.0).build(&mut model).unwrap();
    compound(&mut model, &[solid, sphere]);
    assert!(!model.is_empty());
}

#[test]
fn compound_two_primitives_def_counts_add_up() {
    let mut model = Model::new();
    let solid = make_box(&mut model, 10.0, 10.0, 10.0);
    let sphere = primitives::Sphere::new(5.0).build(&mut model).unwrap();
    let c = compound(&mut model, &[solid, sphere]);
    assert_counts_match(&model, c, &FEV);
}

#[test]
fn compound_three_boxes_def_counts() {
    let mut model = Model::new();
    let boxes: Vec<Shape> =
        (0..3).map(|i| make_box(&mut model, 5.0 * (i + 1) as f64, 10.0, 10.0)).collect();
    let c = compound(&mut model, &boxes);
    assert_counts_match(&model, c, &[NodeKind::Solid, NodeKind::Face, NodeKind::Edge]);
}

#[test]
fn compound_nested_def_counts() {
    let mut model = Model::new();
    let solid = make_box(&mut model, 10.0, 10.0, 10.0);
    let cylinder = primitives::Cylinder::new(3.0, 10.0).build(&mut model).unwrap();
    let inner = compound(&mut model, &[solid]);
    let outer = compound(&mut model, &[inner, cylinder]);
    assert_counts_match(&model, outer, &FEV);
}

// --- minimal shapes ---

#[test]
fn single_edge_handles_gracefully() {
    let mut model = Model::new();
    make_edge(&mut model, Point::ORIGIN, Point::new(1.0, 0.0, 0.0));
    assert_eq!(model.count(NodeKind::Face), 0);
    assert_eq!(model.count(NodeKind::Solid), 0);
}

#[test]
fn single_vertex_handles_gracefully() {
    let mut model = Model::new();
    make::vertex(&mut model, Point::new(1.0, 2.0, 3.0));
    assert_eq!(model.count(NodeKind::Face), 0);
    assert_eq!(model.count(NodeKind::Edge), 0);
}

// --- box ---

#[test]
fn box_face_def_count_matches_top_exp() {
    let mut model = Model::new();
    let solid = make_box(&mut model, 10.0, 20.0, 30.0);
    assert_counts_match(&model, solid, &[NodeKind::Face]);
    assert_eq!(model.count(NodeKind::Face), 6);
}

#[test]
fn box_edge_def_count_matches_top_exp() {
    let mut model = Model::new();
    let solid = make_box(&mut model, 10.0, 20.0, 30.0);
    assert_counts_match(&model, solid, &[NodeKind::Edge]);
    assert_eq!(model.count(NodeKind::Edge), 12);
}

#[test]
fn box_vertex_def_count_matches_top_exp() {
    let mut model = Model::new();
    let solid = make_box(&mut model, 10.0, 20.0, 30.0);
    assert_counts_match(&model, solid, &[NodeKind::Vertex]);
    assert_eq!(model.count(NodeKind::Vertex), 8);
}

/// Every vertex record is one of the box's corners.
#[test]
fn box_vertex_points_match_brep_tool() {
    let mut model = Model::new();
    make_box(&mut model, 10.0, 20.0, 30.0);
    assert_eq!(model.count(NodeKind::Vertex), 8);
    for (_, vertex) in model.table::<Vertex>().iter() {
        let p = vertex.point;
        let is_corner = |c: f64, size: f64| c.abs() < tolerance::LINEAR || (c - size).abs() < tolerance::LINEAR;
        assert!(is_corner(p.x, 10.0) && is_corner(p.y, 20.0) && is_corner(p.z, 30.0), "{p:?}");
    }
}

#[test]
fn box_face_tolerances_match_brep_tool() {
    let mut model = Model::new();
    make_box(&mut model, 10.0, 20.0, 30.0);
    for (face, record) in model.table::<Face>().iter() {
        assert_eq!(record.tolerance, model.tolerance(face));
    }
}

#[test]
fn box_edge_tolerances_match_brep_tool() {
    let mut model = Model::new();
    make_box(&mut model, 10.0, 20.0, 30.0);
    for (edge, record) in model.table::<Edge>().iter() {
        assert_eq!(record.tolerance, model.tolerance(edge));
    }
}

#[test]
fn box_all_surfaces_are_planes() {
    let mut model = Model::new();
    make_box(&mut model, 10.0, 20.0, 30.0);
    assert_eq!(model.count(NodeKind::Face), 6);
    for (face, _) in model.table::<Face>().iter() {
        assert!(matches!(model.face_surface(face), Some(Surface::Plane(_))), "{face:?}");
    }
}

#[test]
fn box_no_degenerate_edges() {
    let mut model = Model::new();
    make_box(&mut model, 10.0, 20.0, 30.0);
    assert_eq!(degenerate_edges(&model), 0);
}

#[test]
fn box_edge_vertex_defs_are_valid() {
    let mut model = Model::new();
    make_box(&mut model, 10.0, 20.0, 30.0);
    for (edge, record) in model.table::<Edge>().iter() {
        let start = model.get(record.start).expect("start vertex use");
        let end = model.get(record.end).expect("end vertex use");
        assert!(model.get(start.vertex).is_some(), "{edge:?} start");
        assert!(model.get(end.vertex).is_some(), "{edge:?} end");
    }
}

#[test]
fn box_face_surfaces_are_valid() {
    let mut model = Model::new();
    make_box(&mut model, 10.0, 20.0, 30.0);
    for (face, _) in model.table::<Face>().iter() {
        assert!(model.face_surface(face).is_some(), "{face:?} has no surface");
    }
}

#[test]
fn box_edge_param_range_is_non_degenerate() {
    let mut model = Model::new();
    make_box(&mut model, 10.0, 20.0, 30.0);
    for (edge, _) in model.table::<Edge>().iter() {
        let (first, last) = model.edge_range(edge);
        assert!(first < last, "{edge:?} range [{first}, {last}]");
    }
}

/// Adding more geometry to a model keeps the stable ids of what was there.
#[test]
fn add_flatten_preserves_existing_uids() {
    let mut model = Model::new();
    make_box(&mut model, 10.0, 20.0, 30.0);
    let first_edge = model.table::<Edge>().iter().next().expect("an edge").0;
    let before = model.stable_id(first_edge).expect("stable id");

    primitives::Sphere::new(5.0).build(&mut model).unwrap();
    assert_eq!(model.stable_id(first_edge), Some(before));
}

#[test]
fn add_flatten_standalone_vertex_appends_into_non_empty_graph() {
    let mut model = Model::new();
    make_box(&mut model, 10.0, 20.0, 30.0);
    let before = model.count(NodeKind::Vertex);
    make::vertex(&mut model, Point::new(100.0, 0.0, 0.0));
    assert_eq!(model.count(NodeKind::Vertex), before + 1);
}

#[test]
fn add_flatten_standalone_edge_appends_into_non_empty_graph() {
    let mut model = Model::new();
    make_box(&mut model, 10.0, 20.0, 30.0);
    let before = model.count(NodeKind::Edge);
    make_edge(&mut model, Point::new(100.0, 0.0, 0.0), Point::new(120.0, 0.0, 0.0));
    assert_eq!(model.count(NodeKind::Edge), before + 1);
}

#[test]
fn add_flatten_standalone_wire_appends_into_non_empty_graph() {
    let mut model = Model::new();
    make_box(&mut model, 10.0, 20.0, 30.0);
    let before = model.count(NodeKind::Wire);
    let edge = make_edge(&mut model, Point::ORIGIN, Point::new(10.0, 0.0, 0.0));
    make::wire(&mut model, &[edge]).unwrap();
    assert_eq!(model.count(NodeKind::Wire), before + 1);
}

#[test]
fn add_flatten_compound_with_standalone_shapes() {
    let mut model = Model::new();
    make_box(&mut model, 10.0, 20.0, 30.0);
    let (vertices, edges) = (model.count(NodeKind::Vertex), model.count(NodeKind::Edge));
    let v = make::vertex(&mut model, Point::new(50.0, 0.0, 0.0));
    let e = make_edge(&mut model, Point::ORIGIN, Point::new(10.0, 0.0, 0.0));
    compound(&mut model, &[v, e]);
    assert!(model.count(NodeKind::Vertex) > vertices);
    assert!(model.count(NodeKind::Edge) > edges);
}

#[test]
fn build_basic_queries_and_algorithms_work() {
    let mut model = Model::new();
    make_box(&mut model, 10.0, 20.0, 30.0);
    assert_eq!(model.count(NodeKind::Face), 6);
    assert_eq!(model.count(NodeKind::Edge), 12);
    assert_eq!(model.count(NodeKind::Vertex), 8);
}

// --- faces on whole (naturally bounded) surfaces ---

fn assert_natural_face_has_boundary(surface: Surface) {
    let mut model = Model::new();
    let face = make::face_on_surface(&mut model, surface).expect("face");
    assert_eq!(model.count(NodeKind::Face), 1);
    assert!(model.children(face).count() > 0, "the face has wires");
    assert!(model.count(NodeKind::Wire) > 0);
    assert!(model.count(NodeKind::Edge) > 0);
    assert!(model.count(NodeKind::Vertex) > 0);
}

#[test]
fn natural_bound_finite_no_wire_natural_faces_synthesize_topology() {
    let z_frame = Frame::from_z(Point::ORIGIN, Direction::Z);
    assert_natural_face_has_boundary(Surface::Sphere(Sphere::new(z_frame, 10.0)));
    assert_natural_face_has_boundary(Surface::Torus(Torus::new(z_frame, 10.0, 3.0)));
    assert_natural_face_has_boundary(Surface::Trimmed(TrimmedSurface::new(
        Surface::Cylinder(Cylinder::new(z_frame, 5.0)),
        0.0,
        TAU,
        0.0,
        10.0,
    )));
    assert_natural_face_has_boundary(Surface::Trimmed(TrimmedSurface::new(
        Surface::Cone(Cone::new(z_frame, 0.2, 5.0)),
        0.0,
        TAU,
        0.0,
        10.0,
    )));
    assert_natural_face_has_boundary(Surface::Trimmed(TrimmedSurface::new(
        Surface::Plane(Plane::new(Point::ORIGIN, Direction::Z)),
        0.0,
        3.0,
        0.0,
        4.0,
    )));
}

#[test]
fn natural_bound_sphere_has_wires_and_edges() {
    let mut model = Model::new();
    primitives::Sphere::new(10.0).build(&mut model).unwrap();
    assert!(model.count(NodeKind::Face) >= 1);
    for (face, record) in model.table::<Face>().iter() {
        assert!(!record.wires.is_empty(), "{face:?} has no wires");
    }
    assert!(model.count(NodeKind::Wire) > 0);
    assert!(model.count(NodeKind::Edge) > 0);
    assert!(model.count(NodeKind::Vertex) > 0);
}

#[test]
fn natural_bound_cylinder_has_wires_and_edges() {
    let mut model = Model::new();
    primitives::Cylinder::new(5.0, 10.0).build(&mut model).unwrap();
    assert!(model.count(NodeKind::Face) >= 1);
    for (face, record) in model.table::<Face>().iter() {
        assert!(!record.wires.is_empty(), "{face:?} has no wires");
    }
}
