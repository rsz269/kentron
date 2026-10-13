//! BRepGraph_ParentExplorer_Test.cxx (`Model::explore_up`)
//!
//! OCCT auto-creates a root product (with an occurrence) for every imported shape, so its
//! upward walks from a box face end at that product; kentron only has the products it is
//! asked to build, so the walk ends at the solid unless the test builds an assembly. Not
//! ported: the link-kind / reference introspection (`*_ExposeChildAndRef`,
//! `*_IsStructural`, `OccurrenceParent_ExposeOccurrenceRef`),
//! `ShapeRootDirectParent_IsOccurrenceNotProductShortcut` (occurrence-wrapped part roots)
//! and `OccurrenceParents_SetRefOccurrenceDefRejectsSharing` (internal setter guard).

use occt_tests::assert_near;
use kentron::math::{Axis, Point, Transform, Transformable, Vector};
use kentron::tolerance;
use kentron::topology::{CoEdge, Face, Model, NodeId, NodeKind, Shape, Shell, Solid};

use crate::graph_support::{box_model, compound, first};

fn translation(x: f64, y: f64, z: f64) -> Transform {
    Transform::translation(Vector::new(x, y, z))
}

fn rotation_z(angle: f64) -> Transform {
    Transform::rotation(Axis::Z, angle)
}

/// Placements compared by where they send a probe point.
fn equivalent(a: &Transform, b: &Transform) -> bool {
    let probe = Point::new(3.0, 5.0, 7.0);
    probe.transformed(a).distance(probe.transformed(b)) <= tolerance::LINEAR
}

fn nodes(explorer: impl Iterator<Item = Shape>) -> Vec<NodeId> {
    explorer.map(|s| s.node).collect()
}

#[test]
fn face_parents_all_count_and_order() {
    let (model, _) = box_model();
    let face = first::<Face>(&model);
    let expected = [NodeId::from(first::<Shell>(&model)), first::<Solid>(&model).into()];
    assert_eq!(nodes(model.explore_up(face)), expected);
}

#[test]
fn face_parents_typed_solid_one_result() {
    let (model, _) = box_model();
    let found = nodes(model.explore_up(first::<Face>(&model)).of_kind(NodeKind::Solid));
    assert_eq!(found, [NodeId::from(first::<Solid>(&model))]);
}

#[test]
fn face_parents_direct_parents_stops_at_immediate_shell() {
    let (model, _) = box_model();
    let found = nodes(model.explore_up(first::<Face>(&model)).direct_only());
    assert_eq!(found, [NodeId::from(first::<Shell>(&model))]);
}

/// The box as a part instanced in an assembly.
fn box_in_assembly() -> (Model, Shape) {
    let (mut model, solid) = box_model();
    let part = model.new_part(solid).unwrap();
    let assembly = model.new_assembly();
    model.add_instance(assembly, part, Transform::IDENTITY).unwrap();
    (model, solid)
}

#[test]
fn avoid_kind_solid_prunes_products() {
    let (model, _) = box_in_assembly();
    let found = model
        .explore_up(first::<Face>(&model))
        .of_kind(NodeKind::Product)
        .skip_above(NodeKind::Solid)
        .count();
    assert_eq!(found, 0);
}

#[test]
fn avoid_kind_emit_boundary_returns_solid_instead_of_products() {
    let (model, _) = box_in_assembly();
    let found = nodes(
        model.explore_up(first::<Face>(&model)).of_kind(NodeKind::Product).stop_at(NodeKind::Solid),
    );
    assert_eq!(found, [NodeId::from(first::<Solid>(&model))]);
}

#[test]
fn avoid_kind_same_as_target_is_ignored() {
    let (model, _) = box_in_assembly();
    let found = nodes(
        model.explore_up(first::<Face>(&model)).of_kind(NodeKind::Solid).skip_above(NodeKind::Solid),
    );
    assert_eq!(found, [NodeId::from(first::<Solid>(&model))]);
}

#[test]
fn all_parents_avoid_solid_prunes_products() {
    let (model, _) = box_in_assembly();
    let found = nodes(model.explore_up(first::<Face>(&model)).skip_above(NodeKind::Solid));
    assert_eq!(found, [NodeId::from(first::<Shell>(&model))]);
}

#[test]
fn all_parents_avoid_solid_emit_boundary_returns_shell_and_solid() {
    let (model, _) = box_in_assembly();
    let found = nodes(model.explore_up(first::<Face>(&model)).stop_at(NodeKind::Solid));
    assert_eq!(found, [NodeId::from(first::<Shell>(&model)), first::<Solid>(&model).into()]);
}

/// One part instanced twice: walking up from its solid reaches the part twice, in two
/// different contexts.
#[test]
fn shared_product_product_parents_keep_distinct_contexts() {
    let (mut model, solid) = box_model();
    let part = model.new_part(solid).unwrap();
    let assembly = model.new_assembly();
    model.add_instance(assembly, part, translation(10.0, 0.0, 0.0)).unwrap();
    model.add_instance(assembly, part, translation(25.0, 0.0, 0.0)).unwrap();

    let mut explorer = model.explore_up(solid.node).of_kind(NodeKind::Product);
    let mut leaf_placements = Vec::new();
    while let Some(product) = explorer.next() {
        if product.node == part {
            leaf_placements.push(explorer.leaf_placement());
        }
    }
    assert_eq!(leaf_placements.len(), 2);
    assert!(!equivalent(&leaf_placements[0], &leaf_placements[1]));
}

/// Six nested assemblies, each placing its child by +level along X.
#[test]
fn deep_product_occurrence_chain_composes_to_top_product() {
    let (mut model, solid) = box_model();
    let mut child = model.new_part(solid).unwrap();
    let mut expected_x = 0.0;
    for level in 1..=6 {
        let parent = model.new_assembly();
        model.add_instance(parent, child, translation(level as f64, 0.0, 0.0)).unwrap();
        expected_x += level as f64;
        child = parent;
    }
    let top = child;

    let mut explorer = model.explore_up(solid.node).of_kind(NodeKind::Product);
    let mut found_top = false;
    while let Some(product) = explorer.next() {
        if product.node == top {
            found_top = true;
            assert_near!(explorer.leaf_placement().translation_part().x, expected_x, tolerance::LINEAR);
        }
    }
    assert!(found_top);
}

#[test]
fn deep_product_occurrence_chain_reaches_nested_occurrences() {
    let mut model = Model::new();
    let part = model.new_assembly();
    let sub = model.new_assembly();
    let top = model.new_assembly();
    let part_instance = model.add_instance(sub, part, Transform::IDENTITY).unwrap();
    let sub_instance = model.add_instance(top, sub, Transform::IDENTITY).unwrap();

    let found = nodes(model.explore_up(part).of_kind(NodeKind::Instance));
    assert_eq!(found.len(), 2);
    assert!(found.contains(&part_instance.into()));
    assert!(found.contains(&sub_instance.into()));
}

#[test]
fn deep_compound_chain_reaches_top_compound() {
    let (mut model, solid) = box_model();
    let mut child = solid;
    for _ in 0..8 {
        child = compound(&mut model, &[child]);
    }
    let found = nodes(model.explore_up(solid.node).of_kind(NodeKind::Compound));
    assert_eq!(found.len(), 8);
    assert!(found.contains(&child.node));
}

/// Compound placements accumulate on the way up.
#[test]
fn vertex_child_of_deep_compound_chain_reaches_top_compound() {
    let mut model = Model::new();
    let vertex = model.add_vertex(Point::new(1.0, 2.0, 3.0), 1.0e-7);
    let mut child = Shape::from(vertex);
    let mut expected_x = 0.0;
    for level in 1..=4 {
        let step = translation(level as f64, 0.0, 0.0);
        child = compound(&mut model, &[child.with_placement(step)]);
        expected_x += level as f64;
    }
    let top = child.node;

    let mut explorer = model.explore_up(vertex).of_kind(NodeKind::Compound);
    let (mut count, mut found_top) = (0, false);
    while let Some(c) = explorer.next() {
        count += 1;
        if c.node == top {
            found_top = true;
            assert_near!(explorer.leaf_placement().translation_part().x, expected_x, tolerance::LINEAR);
        }
    }
    assert_eq!(count, 4);
    assert!(found_top);
    assert!(kentron::validity::check_model(&model).is_valid());
}

/// Placements through compounds, a placed part root, and two assembly levels.
#[test]
fn mixed_product_occurrence_and_compound_chain_preserves_locations() {
    let (mut model, solid) = box_model();
    let comp1_loc = translation(2.0, 0.0, 0.0);
    let comp2_loc = rotation_z(0.25);
    let comp3_loc = translation(0.0, 3.0, 0.0);
    let inst1_loc = translation(0.0, 0.0, 5.0);
    let inst2_loc = rotation_z(-0.5);

    let comp1 = compound(&mut model, &[solid.with_placement(comp1_loc)]);
    let comp2 = compound(&mut model, &[comp1.with_placement(comp2_loc)]);
    let comp3 = compound(&mut model, &[comp2.with_placement(comp3_loc)]);
    let part = model.new_part(comp3).unwrap();
    let mid = model.new_assembly();
    let top = model.new_assembly();
    let inst1 = model.add_instance(mid, part, inst1_loc).unwrap();
    let inst2 = model.add_instance(top, mid, inst2_loc).unwrap();

    let to_part = inst2_loc * inst1_loc;
    let to_comp3 = to_part;
    let to_comp2 = to_comp3 * comp3_loc;
    let to_comp1 = to_comp2 * comp2_loc;
    let to_solid = to_comp1 * comp1_loc;

    let expected: [(NodeId, Transform); 8] = [
        (top.into(), Transform::IDENTITY),
        (inst2.into(), inst2_loc),
        (mid.into(), inst2_loc),
        (inst1.into(), to_part),
        (part.into(), to_part),
        (comp3.node, to_comp3),
        (comp2.node, to_comp2),
        (comp1.node, to_comp1),
    ];
    let mut seen = [false; 8];
    let mut explorer = model.explore_up(solid.node);
    while let Some(ancestor) = explorer.next() {
        if !equivalent(&explorer.leaf_placement(), &to_solid) {
            continue;
        }
        if let Some(i) = expected.iter().position(|(n, _)| *n == ancestor.node) {
            assert!(equivalent(&ancestor.placement, &expected[i].1), "placement of {:?}", ancestor.node);
            seen[i] = true;
        }
    }
    assert!(seen.iter().all(|&s| s), "missing ancestors: {seen:?}");
}

#[test]
fn target_product_avoid_compound_stops_at_compound() {
    let (mut model, solid) = box_model();
    let c = compound(&mut model, &[solid]);
    let part = model.new_part(c).unwrap();
    let assembly = model.new_assembly();
    model.add_instance(assembly, part, Transform::IDENTITY).unwrap();

    let found = model
        .explore_up(solid.node)
        .of_kind(NodeKind::Product)
        .skip_above(NodeKind::Compound)
        .count();
    assert_eq!(found, 0);
}

#[test]
fn co_edge_parents_immediate_wire_is_visible() {
    let (model, _) = box_model();
    let face = first::<Face>(&model);
    let wire = model[model[face].wires[0]].wire;
    let coedge = model[wire].coedges[0];
    let found = nodes(model.explore_up(coedge).take(2));
    assert_eq!(found, [NodeId::from(wire), face.into()]);
}

#[test]
fn co_edge_child_of_compound_immediate_compound_is_visible() {
    let (mut model, _) = box_model();
    let coedge = first::<CoEdge>(&model);
    let wire = model[coedge].wire;
    let c = compound(&mut model, &[coedge.into()]);
    let found = nodes(model.explore_up(coedge).direct_only());
    assert!(found.contains(&wire.into()));
    assert!(found.contains(&c.node));
}

#[test]
fn product_root_has_no_parents() {
    let (mut model, _) = box_model();
    let root = model.new_assembly();
    assert_eq!(model.explore_up(root).count(), 0);
}
