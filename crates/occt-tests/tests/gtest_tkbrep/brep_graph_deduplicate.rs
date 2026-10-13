//! BRepGraph_Deduplicate_Test.cxx: sharing equal geometry and merging identical topology.
//!
//! Not ported: tests adding duplicate pcurves through coedges outside any wire
//! (`MultiplePCurveDups_*`, `PCurveDup_*`, `NoPCurveDuplicates_*`, `*PCurve*_Nullify*`,
//! `AfterDedup_AllInlinePCurvesHaveCurve2d`, `PCurveEntryCount_*`, the disabled
//! `PCurveDedup_*`), TopoDS round trips (`RoundTrip_*`, `Build_SharedTFace_*`),
//! `ParallelBuild_*`, and the history-flag restoration tests (kentron never toggles the
//! history's enabled flag). Tests that check history opt in with `Model::enable_history`.

use std::collections::HashSet;

use kentron::math::Point;
use kentron::primitives::{Cone, Cuboid, Cylinder, Sphere};
use kentron::make;
use kentron::topology::{
    ChangeKind, CoEdge, DeduplicateOptions, DeduplicateReport, Edge, Face, Model, NodeKind,
    Shape,
};
use kentron::tolerance;
use kentron::validity;

use crate::graph_support::{box_model, compound};

fn merge() -> DeduplicateOptions {
    DeduplicateOptions { merge_entities: true, ..Default::default() }
}

fn analyze() -> DeduplicateOptions {
    DeduplicateOptions { analyze_only: true, ..Default::default() }
}

/// A model holding a compound of `copies` imported from `source`.
fn imported(source: &Model, shapes: &[Shape]) -> Model {
    let mut model = Model::new();
    let copies: Vec<Shape> = shapes.iter().map(|&s| model.import(source, s).unwrap()).collect();
    compound(&mut model, &copies);
    model
}

/// The first face of a primitive built in `source`.
fn first_face(source: &Model, solid: Shape) -> Shape {
    crate::first(source, solid, NodeKind::Face)
}

/// A model with a compound of `n` copies of the first face of a 10 x 20 x 30 box.
fn copied_box_faces(n: usize) -> Model {
    let (source, solid) = box_model();
    let face = first_face(&source, solid);
    imported(&source, &vec![face; n])
}

/// Two copies of a box face plus the first face of a cylinder.
fn mixed_compound() -> Model {
    let (mut source, solid) = box_model();
    let face = first_face(&source, solid);
    let cylinder = Cylinder::new(5.0, 15.0).build(&mut source).unwrap();
    let cylinder_face = first_face(&source, cylinder);
    imported(&source, &[face, face, cylinder_face])
}

/// ((copy, copy), copy) of a box face.
fn nested_compound() -> Model {
    let (source, solid) = box_model();
    let face = first_face(&source, solid);
    let mut model = Model::new();
    let a = model.import(&source, face).unwrap();
    let b = model.import(&source, face).unwrap();
    let inner = compound(&mut model, &[a, b]);
    let c = model.import(&source, face).unwrap();
    compound(&mut model, &[inner, c]);
    model
}

fn three_distinct_primitives() -> Model {
    let mut model = Model::new();
    let a = Cuboid::new(10.0, 20.0, 30.0).build(&mut model).unwrap();
    let b = Sphere::new(5.0).build(&mut model).unwrap();
    let c = Cone::new(3.0, 1.0, 8.0).build(&mut model).unwrap();
    compound(&mut model, &[a, b, c]);
    model
}

/// `n` independently built 10 x 20 x 30 boxes in a compound.
fn identical_boxes(n: usize) -> Model {
    let mut model = Model::new();
    let boxes: Vec<Shape> = (0..n).map(|_| crate::make_box(&mut model, 10.0, 20.0, 30.0)).collect();
    compound(&mut model, &boxes);
    model
}

/// A box edge and two copies of it.
fn three_copied_edges() -> Model {
    let (source, solid) = box_model();
    let edge = crate::first(&source, solid, NodeKind::Edge);
    imported(&source, &[edge, edge, edge])
}

fn face_touching_boxes() -> Model {
    let mut model = Model::new();
    let a = Cuboid::from_corners(Point::ORIGIN, Point::new(10.0, 10.0, 10.0)).build(&mut model).unwrap();
    let b = Cuboid::from_corners(Point::new(0.0, 0.0, 10.0), Point::new(10.0, 10.0, 20.0))
        .build(&mut model)
        .unwrap();
    compound(&mut model, &[a, b]);
    model
}

/// Distinct surface slots used by the live faces.
fn distinct_face_surfaces(model: &Model) -> usize {
    let slots: HashSet<_> = model
        .table::<Face>()
        .iter()
        .filter(|(_, f)| model.get(f.surface).is_some())
        .map(|(_, f)| f.surface)
        .collect();
    slots.len()
}

/// Distinct curve slots used by the live edges.
fn distinct_edge_curves(model: &Model) -> usize {
    let slots: HashSet<_> = model
        .table::<Edge>()
        .iter()
        .filter(|(_, e)| model.get(e.curve).is_some())
        .map(|(_, e)| e.curve)
        .collect();
    slots.len()
}

fn events_named(model: &Model, operation: &str) -> usize {
    model.history().unwrap().events().iter().filter(|e| e.operation == operation).count()
}

// ---- geometry sharing ----

#[test]
fn analyze_only_does_not_rewrite() {
    let mut model = copied_box_faces(2);
    assert_eq!(distinct_face_surfaces(&model), 2);
    let report = model.deduplicate_with(analyze());
    assert_eq!(report.shared_surfaces, 0);
    assert_eq!(distinct_face_surfaces(&model), 2);
}

#[test]
fn analyze_only_reports_canonical_candidates() {
    let mut model = copied_box_faces(2);
    let report = model.deduplicate_with(analyze());
    assert_eq!(report.distinct_surfaces, 1);
    assert_eq!(report.distinct_curves, 4);
    assert_eq!(report.shared_surfaces, 0);
    assert_eq!(report.shared_curves, 0);
}

/// OCCT records one history event per canonicalized surface or curve; kentron does not
/// record geometry sharing (see `DeduplicateOptions::record_history`).
#[test]
fn canonicalize_surfaces_rewrites_without_history() {
    let mut model = copied_box_faces(2);
    model.enable_history();
    let history_before = model.history().unwrap().len();
    let report = model.deduplicate();
    assert_eq!(report.shared_surfaces, 1);
    assert_eq!(distinct_face_surfaces(&model), 1);
    assert_eq!(report.history_events, 0);
    assert_eq!(model.history().unwrap().len(), history_before);
}

#[test]
fn canonicalize_curves_rewrites_and_reduces_unique() {
    let mut model = copied_box_faces(2);
    assert_eq!(distinct_edge_curves(&model), 8);
    let report = model.deduplicate();
    assert_eq!(report.shared_curves, 4);
    assert_eq!(distinct_edge_curves(&model), 4);
}

#[test]
fn history_off_does_not_add_history() {
    let mut model = copied_box_faces(2);
    model.enable_history();
    let report = model.deduplicate_with(DeduplicateOptions {
        record_history: false,
        merge_entities: true,
        ..Default::default()
    });
    assert_eq!(report.shared_surfaces, 1);
    assert_eq!(report.shared_curves, 4);
    assert!(report.merge_applied);
    assert_eq!(model.history().unwrap().len(), 0);
}

#[test]
fn default_options_perform_works() {
    let mut model = copied_box_faces(2);
    let report = model.deduplicate();
    assert_eq!(report.distinct_surfaces, 1);
    assert_eq!(report.shared_surfaces, 1);
    assert_eq!(report.shared_curves, 4);
    assert!(!report.merge_applied);
}

#[test]
fn single_face_no_surface_rewrite() {
    let report = copied_box_faces(1).deduplicate();
    assert_eq!(report.shared_surfaces, 0);
}

#[test]
fn empty_model_returns_empty_report() {
    let report = Model::new().deduplicate_with(merge());
    assert_eq!(report, DeduplicateReport::default());
}

#[test]
fn idempotent_second_run_no_rewrites() {
    let mut model = copied_box_faces(2);
    let first = model.deduplicate();
    assert_eq!(first.shared_surfaces, 1);
    assert_eq!(first.shared_curves, 4);
    let (surfaces, curves) = (distinct_face_surfaces(&model), distinct_edge_curves(&model));
    let second = model.deduplicate();
    assert_eq!(second.shared_surfaces, 0);
    assert_eq!(second.shared_curves, 0);
    assert_eq!(distinct_face_surfaces(&model), surfaces);
    assert_eq!(distinct_edge_curves(&model), curves);
}

#[test]
fn full_box_all_surfaces_unique() {
    let (mut model, _) = box_model();
    let report = model.deduplicate();
    assert_eq!(report.distinct_surfaces, 6);
    assert_eq!(report.shared_surfaces, 0);
}

#[test]
fn result_counters_consistency() {
    let mut model = copied_box_faces(2);
    assert_eq!(model.count(NodeKind::Face), 2);
    assert_eq!(model.count(NodeKind::Edge), 8);
    assert_eq!(model.table::<CoEdge>().iter().filter(|(_, c)| model.get(c.pcurve).is_some()).count(), 8);
    let report = model.deduplicate_with(analyze());
    assert_eq!(report.distinct_surfaces, 1);
    assert_eq!(report.distinct_curves, 4);
}

#[test]
fn empty_compound_no_rewrites() {
    let mut model = Model::new();
    compound(&mut model, &[]);
    let report = model.deduplicate();
    assert_eq!(report.shared_surfaces, 0);
    assert_eq!(report.shared_curves, 0);
}

#[test]
fn multiple_copies_n_way_dedup() {
    for n in [3, 4, 10] {
        let mut model = copied_box_faces(n);
        assert_eq!(model.count(NodeKind::Face), n);
        let report = model.deduplicate();
        assert_eq!(report.distinct_surfaces, 1);
        assert_eq!(report.shared_surfaces, n - 1);
        assert_eq!(distinct_face_surfaces(&model), 1);
    }
}

#[test]
fn mixed_geometry_only_identical_deduped() {
    let mut model = mixed_compound();
    assert_eq!(model.count(NodeKind::Face), 3);
    assert_eq!(model.count(NodeKind::Edge), 11);
    let report = model.deduplicate();
    assert_eq!(report.distinct_surfaces, 2);
    assert_eq!(report.distinct_curves, 7);
    assert_eq!(report.shared_surfaces, 1);
    assert_eq!(report.shared_curves, 4);
    assert_eq!(distinct_face_surfaces(&model), 2);
    assert_eq!(distinct_edge_curves(&model), 7);
}

#[test]
fn idempotent_mixed_compound() {
    let mut model = mixed_compound();
    let first = model.deduplicate();
    assert_eq!((first.shared_surfaces, first.shared_curves), (1, 4));
    let (surfaces, curves) = (distinct_face_surfaces(&model), distinct_edge_curves(&model));
    model.deduplicate();
    assert_eq!(distinct_face_surfaces(&model), surfaces);
    assert_eq!(distinct_edge_curves(&model), curves);
}

#[test]
fn default_options_have_expected_values() {
    let options = DeduplicateOptions::default();
    assert!(!options.analyze_only);
    assert!(options.record_history);
    assert!(!options.merge_entities);
    assert_eq!(options.angular_tolerance, tolerance::ANGULAR);
    assert_eq!(options.linear_tolerance, tolerance::LINEAR);
}

#[test]
fn nested_compound_all_copies_deduped() {
    let mut model = nested_compound();
    assert_eq!(model.count(NodeKind::Face), 3);
    assert_eq!(model.count(NodeKind::Edge), 12);
    let report = model.deduplicate();
    assert_eq!(report.distinct_surfaces, 1);
    assert_eq!(report.distinct_curves, 4);
    assert_eq!(report.shared_surfaces, 2);
    assert_eq!(report.shared_curves, 8);
    assert_eq!(distinct_face_surfaces(&model), 1);
    assert_eq!(distinct_edge_curves(&model), 4);
}

#[test]
fn three_distinct_primitives_minimal_dedup() {
    let mut model = three_distinct_primitives();
    assert_eq!(model.count(NodeKind::Face), 10);
    assert_eq!(model.count(NodeKind::Edge), 18);
    let report = model.deduplicate();
    // The box's bottom face and the cone's base lie on the same plane.
    assert_eq!(report.distinct_surfaces, 9);
    assert_eq!(report.shared_surfaces, 1);
    assert_eq!(distinct_face_surfaces(&model), 9);
    assert!(distinct_edge_curves(&model) <= 18);
}

#[test]
fn two_identical_boxes_surfaces_and_curves_deduped() {
    let mut model = identical_boxes(2);
    assert_eq!(model.count(NodeKind::Face), 12);
    assert_eq!(model.count(NodeKind::Edge), 24);
    let report = model.deduplicate();
    assert_eq!(report.distinct_surfaces, 6);
    assert_eq!(report.distinct_curves, 12);
    assert_eq!(report.shared_surfaces, 6);
    assert_eq!(report.shared_curves, 12);
    assert_eq!(distinct_face_surfaces(&model), 6);
    assert_eq!(distinct_edge_curves(&model), 12);
    let second = model.deduplicate();
    assert_eq!((second.shared_surfaces, second.shared_curves), (0, 0));
    assert_eq!(distinct_face_surfaces(&model), 6);
    assert_eq!(distinct_edge_curves(&model), 12);
}

#[test]
fn curve_rewrite_count_matches_duplicate_edge_curves() {
    let mut model = copied_box_faces(2);
    assert_eq!(model.count(NodeKind::Edge), 8);
    assert_eq!(model.deduplicate_with(analyze()).distinct_curves, 4);
    assert_eq!(model.deduplicate().shared_curves, 4);
}

#[test]
fn after_dedup_all_faces_and_edges_keep_geometry() {
    let mut model = copied_box_faces(4);
    model.deduplicate();
    for (face, _) in model.table::<Face>().iter() {
        assert!(model.surface(face).is_some(), "{face:?} lost its surface");
    }
    for (edge, _) in model.table::<Edge>().iter() {
        if !model.is_degenerate(edge) {
            assert!(model.edge_curve(edge).is_some(), "{edge:?} lost its curve");
        }
    }
}

#[test]
fn counts_unchanged_after_dedup() {
    let mut model = copied_box_faces(2);
    let coedges = model.count(NodeKind::CoEdge);
    model.deduplicate();
    assert_eq!(model.count(NodeKind::Face), 2);
    assert_eq!(model.count(NodeKind::Edge), 8);
    assert_eq!(model.count(NodeKind::CoEdge), coedges);
}

#[test]
fn two_copied_sphere_faces_deduped() {
    let mut source = Model::new();
    let sphere = Sphere::new(10.0).build(&mut source).unwrap();
    let face = first_face(&source, sphere);
    let mut model = imported(&source, &[face, face]);
    assert_eq!(model.count(NodeKind::Face), 2);
    assert_eq!(model.count(NodeKind::Edge), 6);
    let report = model.deduplicate();
    assert_eq!(report.distinct_surfaces, 1);
    assert_eq!(report.shared_surfaces, 1);
    let faces: Vec<_> = model.table::<Face>().iter().map(|(_, f)| f.surface).collect();
    assert_eq!(faces[0], faces[1]);
}

#[test]
fn two_copied_cylinder_faces_deduped() {
    let mut source = Model::new();
    let cylinder = Cylinder::new(5.0, 20.0).build(&mut source).unwrap();
    let face = first_face(&source, cylinder);
    let mut model = imported(&source, &[face, face]);
    assert_eq!(model.count(NodeKind::Face), 2);
    assert_eq!(model.count(NodeKind::Edge), 6);
    let report = model.deduplicate();
    assert_eq!(report.distinct_surfaces, 1);
    assert_eq!(report.distinct_curves, 3);
    assert_eq!(report.shared_surfaces, 1);
    assert_eq!(report.shared_curves, 3);
}

#[test]
fn different_sized_cylinders_not_deduped() {
    let mut source = Model::new();
    let a = Cylinder::new(5.0, 20.0).build(&mut source).unwrap();
    let b = Cylinder::new(10.0, 20.0).build(&mut source).unwrap();
    let fa = first_face(&source, a);
    let fb = first_face(&source, b);
    let mut model = imported(&source, &[fa, fb]);
    assert_eq!(model.count(NodeKind::Edge), 6);
    let report = model.deduplicate();
    assert_eq!(report.distinct_surfaces, 2);
    assert_eq!(report.distinct_curves, 6);
    assert_eq!(report.shared_surfaces, 0);
    assert_eq!(report.shared_curves, 0);
    assert_eq!(report.history_events, 0);
}

#[test]
fn shared_geometry_is_one_slot() {
    let mut model = copied_box_faces(2);
    let surfaces: Vec<_> = model.table::<Face>().iter().map(|(_, f)| f.surface).collect();
    assert_ne!(surfaces[0], surfaces[1]);
    model.deduplicate();
    let faces: Vec<_> = model.table::<Face>().iter().map(|(id, f)| (id, f.surface)).collect();
    assert_eq!(faces[0].1, faces[1].1);
    assert!(model.surface(faces[0].0).is_some());
    assert!(model.surface(faces[1].0).is_some());
    assert_eq!(distinct_edge_curves(&model), 4);
}

#[test]
fn analyze_only_leaves_geometry_slots_alone() {
    let mut model = copied_box_faces(2);
    let surfaces: Vec<_> = model.table::<Face>().iter().map(|(_, f)| f.surface).collect();
    let curves: Vec<_> = model.table::<Edge>().iter().map(|(_, e)| e.curve).collect();
    model.deduplicate_with(analyze());
    let surfaces_after: Vec<_> = model.table::<Face>().iter().map(|(_, f)| f.surface).collect();
    let curves_after: Vec<_> = model.table::<Edge>().iter().map(|(_, e)| e.curve).collect();
    assert_eq!(surfaces_after, surfaces);
    assert_eq!(curves_after, curves);
}

// ---- entity merging ----

#[test]
fn merge_vertices_shared_vertices_reduced() {
    let mut model = copied_box_faces(2);
    let vertices = model.count(NodeKind::Vertex);
    let report = model.deduplicate_with(merge());
    assert!(report.merged(NodeKind::Vertex) > 0);
    assert!(report.merge_applied);
    assert!(model.count(NodeKind::Vertex) < vertices);
}

#[test]
fn merge_edges_shared_edges_reduced() {
    let report = copied_box_faces(2).deduplicate_with(merge());
    assert!(report.merged(NodeKind::Edge) > 0);
}

#[test]
fn merge_faces_identical_faces_merged() {
    let mut model = copied_box_faces(2);
    let report = model.deduplicate_with(merge());
    assert_eq!(report.merged(NodeKind::Wire), 1);
    assert_eq!(report.merged(NodeKind::Face), 1);
    model.compact();
    assert_eq!(model.count(NodeKind::Face), 1);
    assert!(validity::check_model_structure(&model).is_valid());
}

#[test]
fn no_merge_by_default() {
    let report = copied_box_faces(2).deduplicate();
    for kind in [NodeKind::Vertex, NodeKind::Edge, NodeKind::Wire, NodeKind::Face] {
        assert_eq!(report.merged(kind), 0);
    }
    assert!(!report.merge_applied);
}

#[test]
fn analyze_only_merge_counts_only() {
    let mut model = copied_box_faces(2);
    let vertices = model.count(NodeKind::Vertex);
    let report = model.deduplicate_with(DeduplicateOptions { analyze_only: true, ..merge() });
    assert!(report.merged(NodeKind::Vertex) > 0);
    assert_eq!(model.count(NodeKind::Vertex), vertices);
    assert!(!report.merge_applied);
}

#[test]
fn history_records_merge_events() {
    let mut model = copied_box_faces(2);
    model.enable_history();
    let report = model.deduplicate_with(merge());
    let merged: usize = NodeKind::ALL.iter().map(|&k| report.merged(k)).sum();
    assert!(merged > 0);
    assert_eq!(report.history_events, merged);
    assert_eq!(model.history().unwrap().len(), merged);
    for event in model.history().unwrap().events() {
        assert_eq!(event.kind, ChangeKind::Replaced);
        assert_eq!(event.mapping.len(), 1);
        let (old, images) = &event.mapping[0];
        assert_eq!(images.len(), 1);
        assert_eq!(old.kind, images[0].kind);
        assert_eq!(model.history().unwrap().find_derived(*old), images.clone());
    }
    assert_eq!(events_named(&model, "deduplicate"), merged);
}

#[test]
fn merging_without_enabled_history_records_nothing() {
    let mut model = copied_box_faces(2);
    let report = model.deduplicate_with(merge());
    assert!(report.merge_applied);
    assert_eq!(report.history_events, 0);
    assert!(model.history().is_none());
}

#[test]
fn after_merge_validate_no_issues() {
    let mut model = copied_box_faces(2);
    model.deduplicate_with(merge());
    model.compact();
    assert!(validity::check_model_structure(&model).is_valid());
}

#[test]
fn merge_faces_same_surface_different_wires_not_merged() {
    let mut model = Model::new();
    let rectangle = |model: &mut Model, x0: f64, x1: f64, y1: f64| {
        let p = [
            Point::new(x0, 0.0, 0.0),
            Point::new(x1, 0.0, 0.0),
            Point::new(x1, y1, 0.0),
            Point::new(x0, y1, 0.0),
        ];
        let edges: Vec<Shape> =
            (0..4).map(|i| crate::make_edge(model, p[i], p[(i + 1) % 4])).collect();
        let wire = make::wire(model, &edges).unwrap();
        make::planar_face(model, wire).unwrap()
    };
    let a = rectangle(&mut model, 0.0, 10.0, 20.0);
    let b = rectangle(&mut model, 30.0, 40.0, 10.0);
    compound(&mut model, &[a, b]);
    let report = model.deduplicate_with(merge());
    let surfaces: Vec<_> = model.table::<Face>().iter().map(|(_, f)| f.surface).collect();
    assert_eq!(surfaces[0], surfaces[1]);
    assert_eq!(model.count(NodeKind::Face), 2);
    assert_eq!(report.merged(NodeKind::Face), 0);
}

fn coedges_with_pcurve(model: &Model) -> usize {
    model.table::<CoEdge>().iter().filter(|(_, c)| model.get(c.pcurve).is_some()).count()
}

fn seam_coedges(model: &Model) -> usize {
    model.table::<CoEdge>().iter().filter(|&(id, _)| model.is_seam(id)).count()
}

#[test]
fn sphere_seam_pcurves_survive_dedup_compact() {
    let mut model = Model::new();
    let sphere = Sphere::new(10.0).build(&mut model).unwrap();
    let sphere_sid = model.node_stable_id(sphere.node).unwrap();
    assert_eq!(coedges_with_pcurve(&model), 4);
    assert_eq!(seam_coedges(&model), 2);
    model.deduplicate_with(merge());
    assert_eq!(coedges_with_pcurve(&model), 4);
    model.compact();
    assert_eq!(coedges_with_pcurve(&model), 4);
    assert_eq!(seam_coedges(&model), 2);
    assert!(validity::check_model_structure(&model).is_valid());
    let sphere = Shape::from(model.find_stable_node(sphere_sid).unwrap());
    assert!(validity::check(&model, sphere).is_valid());
}

#[test]
fn three_copied_cylinder_faces_pcurves_survive_full_pipeline() {
    let mut source = Model::new();
    let cylinder = Cylinder::new(5.0, 20.0).build(&mut source).unwrap();
    let face = first_face(&source, cylinder);
    let mut model = imported(&source, &[face, face, face]);
    assert_eq!(model.count(NodeKind::Face), 3);
    assert_eq!(model.count(NodeKind::Edge), 9);
    assert_eq!(coedges_with_pcurve(&model), 12);
    model.deduplicate_with(merge());
    let report = model.compact();
    assert_eq!(report.removed(NodeKind::Edge), 6);
    assert_eq!(coedges_with_pcurve(&model), 4);
    assert!(validity::check_model_structure(&model).is_valid());
    for (face, _) in model.table::<Face>().iter() {
        assert!(validity::check(&model, face.into()).is_valid());
    }
}

fn assert_three_boxes_merged(report: &DeduplicateReport, model: &mut Model) {
    assert_eq!(report.merged(NodeKind::Vertex), 16);
    assert_eq!(report.merged(NodeKind::Edge), 24);
    assert_eq!(report.merged(NodeKind::Wire), 12);
    assert_eq!(report.merged(NodeKind::Face), 12);
    let compacted = model.compact();
    assert_eq!(compacted.removed(NodeKind::Face), 12);
    assert_eq!(compacted.removed(NodeKind::Edge), 24);
    for (_, face) in model.table::<Face>().iter() {
        assert_eq!(face.wires.len(), 1);
    }
    assert!(validity::check_model_structure(model).is_valid());
}

#[test]
fn three_copied_boxes_dedup_merge_compact_faces_keep_outer_wires() {
    let mut model = identical_boxes(3);
    assert_eq!(model.count(NodeKind::Face), 18);
    assert_eq!(model.count(NodeKind::Solid), 3);
    let report = model.deduplicate_with(merge());
    assert_three_boxes_merged(&report, &mut model);
    for (solid, _) in model.table::<kentron::topology::Solid>().iter() {
        assert!(validity::check(&model, solid.into()).is_valid());
    }
}

#[test]
fn two_boxes_merge_and_stay_valid() {
    let mut model = Model::new();
    let a = Cuboid::new(10.0, 10.0, 10.0).build(&mut model).unwrap();
    let b = Cuboid::new(10.0, 10.0, 10.0).build(&mut model).unwrap();
    compound(&mut model, &[a, b]);
    assert_eq!(model.count(NodeKind::Solid), 2);
    model.deduplicate_with(merge());
    model.compact();
    assert!(validity::check_model_structure(&model).is_valid());
    for (solid, _) in model.table::<kentron::topology::Solid>().iter() {
        assert!(validity::check(&model, solid.into()).is_valid());
    }
}

#[test]
fn edge_merge_three_copied_edges_validates() {
    let mut model = three_copied_edges();
    assert_eq!(model.count(NodeKind::Edge), 3);
    let report = model.deduplicate_with(merge());
    assert_eq!(report.merged(NodeKind::Edge), 2);
    model.compact();
    assert!(validity::check_model_structure(&model).is_valid());
}

#[test]
fn two_face_touching_boxes_edges_with_two_faces() {
    let mut model = face_touching_boxes();
    model.deduplicate_with(merge());
    let with_two_faces = model
        .table::<Edge>()
        .iter()
        .filter(|&(edge, _)| {
            let faces: HashSet<_> = model.ancestors_of::<Face>(edge).into_iter().collect();
            faces.len() == 2
        })
        .count();
    assert_eq!(with_two_faces, 16);
    model.compact();
    assert!(validity::check_model_structure(&model).is_valid());
}

#[test]
fn two_spheres_different_size_not_merged() {
    let mut model = Model::new();
    let a = Sphere::new(5.0).build(&mut model).unwrap();
    let b = Sphere::new(10.0).build(&mut model).unwrap();
    compound(&mut model, &[a, b]);
    let report = model.deduplicate_with(merge());
    assert_eq!(report.merged(NodeKind::Face), 0);
    assert_eq!(report.merged(NodeKind::Edge), 0);
    assert_eq!(model.count(NodeKind::Face), 2);
}

#[test]
fn sphere_self_loop_edges_survive_merge() {
    let mut model = Model::new();
    let sphere = Sphere::new(10.0).build(&mut model).unwrap();
    let sid = model.node_stable_id(sphere.node).unwrap();
    model.deduplicate_with(merge());
    model.compact();
    assert!(validity::check_model_structure(&model).is_valid());
    let sphere = Shape::from(model.find_stable_node(sid).unwrap());
    assert!(validity::check(&model, sphere).is_valid());
}
