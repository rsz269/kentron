//! BRepGraph_Assembly_Test.cxx: products (parts and assemblies) and their instances.
//!
//! OCCT creates a root product (and an occurrence) for every shape added to a graph;
//! kentron creates products only on request (`Model::new_part`). Document roots are not
//! stored either: `Model::root_products` returns the products that are not instantiated.
//!
//! Not ported: auto-product and import-option tests (`Build_*AutoCreatesRootProduct`,
//! `Add_*`), the `ShapesView_*` reconstruction tests, `parentOccurrence` disambiguation
//! (`LinkProducts_ParentOccurrenceMustMatchParentProduct` and the nested-occurrence
//! filters), mutation-guard tests (`MutProduct_RAII`, `MutInvalidAssemblyDefs_*`),
//! `NodeId_Helpers`, `UID_IsAssembly`.

use kentron::math::{Transform, Vector};
use kentron::topology::{Model, NodeKind, Product, ProductId, Shape};

use crate::graph_support::box_model;

/// A box and a part holding it.
fn part_model() -> (Model, Shape, ProductId) {
    let (mut model, solid) = box_model();
    let part = model.new_part(solid).unwrap();
    (model, solid, part)
}

fn translation(x: f64, y: f64, z: f64) -> Transform {
    Transform::translation(Vector::new(x, y, z))
}

fn near(a: f64, b: f64) -> bool {
    (a - b).abs() <= kentron::tolerance::LINEAR
}

#[test]
fn new_part_is_part() {
    let (model, _, part) = part_model();
    assert!(model[part].is_part());
    assert!(!model[part].is_assembly());
    assert_eq!(model.root_products(), vec![part]);
}

#[test]
fn new_part_with_invalid_shape_fails() {
    let (mut model, solid, part) = part_model();
    assert!(model.new_part(part).is_err());
    model.remove_node(solid.node);
    assert!(model.new_part(solid).is_err());
}

#[test]
fn new_empty_product_is_neither_part_nor_assembly() {
    let mut model = Model::new();
    let product = model.new_assembly();
    assert!(!model[product].is_assembly());
    assert!(!model[product].is_part());
}

#[test]
fn add_instance_links_products() {
    let (mut model, _, part) = part_model();
    let assembly = model.new_assembly();
    let instance = model.add_instance(assembly, part, translation(100.0, 0.0, 0.0)).unwrap();
    assert_eq!(model[instance].target, part);
    let (parent, placement) = model.instance_placement(instance).unwrap();
    assert_eq!(parent, assembly);
    assert!(near(placement.translation_part().x, 100.0));
    assert_eq!(model[assembly].children.len(), 1);
    assert_eq!(model[model[assembly].children[0]].instance, instance);
    assert!(model[assembly].is_assembly());
}

#[test]
fn shared_part_multiple_instances() {
    let (mut model, _, part) = part_model();
    let assembly = model.new_assembly();
    let a = model.add_instance(assembly, part, translation(100.0, 0.0, 0.0)).unwrap();
    let b = model.add_instance(assembly, part, translation(200.0, 0.0, 0.0)).unwrap();
    assert_ne!(a, b);
    assert_eq!(model[a].target, model[b].target);
    assert_eq!(model[assembly].children.len(), 2);
    assert_eq!(model.instances_of(part), vec![a, b]);
}

#[test]
fn root_products_follow_instancing() {
    let (mut model, _, part) = part_model();
    assert_eq!(model.root_products(), vec![part]);
    let assembly = model.new_assembly();
    model.add_instance(assembly, part, Transform::IDENTITY).unwrap();
    assert_eq!(model.root_products(), vec![assembly]);
}

#[test]
fn remove_instance_updates_parent() {
    let (mut model, _, part) = part_model();
    let assembly = model.new_assembly();
    let instance = model.add_instance(assembly, part, Transform::IDENTITY).unwrap();
    let child = model[assembly].children[0];
    model.remove_instance(instance).unwrap();
    assert!(model.get(instance).is_none());
    assert!(model[assembly].children.is_empty());
    assert!(model.get(child).is_none());
}

#[test]
fn remove_assembly_removes_its_instances() {
    let (mut model, _, part) = part_model();
    let assembly = model.new_assembly();
    let a = model.add_instance(assembly, part, Transform::IDENTITY).unwrap();
    let b = model.add_instance(assembly, part, Transform::IDENTITY).unwrap();
    model.remove_node(assembly);
    assert!(model.get(assembly).is_none());
    assert!(model.get(a).is_none());
    assert!(model.get(b).is_none());
    assert!(model.get(part).is_some());
}

#[test]
fn remove_part_removes_instances_of_it() {
    let (mut model, _, part) = part_model();
    let assembly = model.new_assembly();
    let instance = model.add_instance(assembly, part, Transform::IDENTITY).unwrap();
    model.remove_node(part);
    assert!(model.get(part).is_none());
    assert!(model.get(instance).is_none());
    assert!(model[assembly].children.is_empty());
}

#[test]
fn remove_nested_instance_keeps_products() {
    let (mut model, _, leaf) = part_model();
    let middle = model.new_assembly();
    let top = model.new_assembly();
    let upper = model.add_instance(top, middle, translation(10.0, 0.0, 0.0)).unwrap();
    let lower = model.add_instance(middle, leaf, translation(0.0, 20.0, 0.0)).unwrap();
    model.remove_instance(upper).unwrap();
    assert!(model.get(upper).is_none());
    assert!(model.get(lower).is_some());
    for product in [leaf, middle, top] {
        assert!(model.get(product).is_some());
    }
}

#[test]
fn touching_a_product_bumps_its_revision() {
    let (mut model, _, part) = part_model();
    model.touch(part);
    assert!(model.revision(part).unwrap().own > 0);
}

#[test]
fn move_instance_sets_placement() {
    let (mut model, _, part) = part_model();
    let assembly = model.new_assembly();
    let instance = model.add_instance(assembly, part, Transform::IDENTITY).unwrap();
    model.move_instance(instance, &translation(50.0, 0.0, 0.0)).unwrap();
    let (_, placement) = model.instance_placement(instance).unwrap();
    assert!(near(placement.translation_part().x, 50.0));
}

#[test]
fn instance_placements_are_local() {
    let (mut model, _, part) = part_model();
    let sub = model.new_assembly();
    let root = model.new_assembly();
    let upper = model.add_instance(root, sub, translation(0.0, 200.0, 0.0)).unwrap();
    let lower = model.add_instance(sub, part, translation(100.0, 0.0, 0.0)).unwrap();
    let (_, p) = model.instance_placement(lower).unwrap();
    assert!(near(p.translation_part().x, 100.0) && near(p.translation_part().y, 0.0));
    let (_, p) = model.instance_placement(upper).unwrap();
    assert!(near(p.translation_part().x, 0.0) && near(p.translation_part().y, 200.0));
}

#[test]
fn three_level_nesting_placements() {
    let (mut model, _, leaf) = part_model();
    let middle = model.new_assembly();
    let root = model.new_assembly();
    let top = model.new_assembly();
    let a = model.add_instance(top, root, translation(0.0, 0.0, 3.0)).unwrap();
    let b = model.add_instance(root, middle, translation(0.0, 2.0, 0.0)).unwrap();
    let c = model.add_instance(middle, leaf, translation(1.0, 0.0, 0.0)).unwrap();
    assert!(near(model.instance_placement(c).unwrap().1.translation_part().x, 1.0));
    assert!(near(model.instance_placement(b).unwrap().1.translation_part().y, 2.0));
    assert!(near(model.instance_placement(a).unwrap().1.translation_part().z, 3.0));
    assert_eq!(model.root_products(), vec![top]);
    // The placed leaf in the top assembly accumulates all three placements.
    let placed: Vec<Shape> = model.explore(Shape::from(top), NodeKind::Solid).collect();
    assert_eq!(placed.len(), 1);
    let t = placed[0].placement.translation_part();
    assert!(near(t.x, 1.0) && near(t.y, 2.0) && near(t.z, 3.0));
}

#[test]
fn shared_part_distinct_paths_give_distinct_placements() {
    let (mut model, _, part) = part_model();
    let assembly = model.new_assembly();
    model.add_instance(assembly, part, translation(100.0, 0.0, 0.0)).unwrap();
    model.add_instance(assembly, part, translation(0.0, 200.0, 0.0)).unwrap();
    let placed: Vec<Shape> = model.explore(Shape::from(assembly), NodeKind::Solid).collect();
    assert_eq!(placed.len(), 2);
    // The same solid (OCCT's IsPartner), placed twice.
    assert_eq!(placed[0].node, placed[1].node);
    assert!(near(placed[0].placement.translation_part().x, 100.0));
    assert!(near(placed[1].placement.translation_part().y, 200.0));
}

#[test]
fn node_count_includes_products_and_instances() {
    let (mut model, _, part) = part_model();
    let total = |m: &Model| NodeKind::ALL.iter().map(|&k| m.count(k)).sum::<usize>();
    let before = total(&model);
    let assembly = model.new_assembly();
    model.add_instance(assembly, part, Transform::IDENTITY).unwrap();
    assert_eq!(total(&model), before + 2);
}

#[test]
fn product_and_instance_counts() {
    let (mut model, _, part) = part_model();
    let assembly = model.new_assembly();
    model.add_instance(assembly, part, Transform::IDENTITY).unwrap();
    model.add_instance(assembly, part, Transform::IDENTITY).unwrap();
    assert_eq!(model.table::<Product>().iter().count(), 2);
    assert_eq!(model.count(NodeKind::Instance), 2);
}

#[test]
fn add_instance_invalid_products_fail() {
    let (mut model, _, part) = part_model();
    let missing = ProductId::from_index(999);
    assert!(model.add_instance(missing, part, Transform::IDENTITY).is_err());
    let assembly = model.new_assembly();
    assert!(model.add_instance(assembly, missing, Transform::IDENTITY).is_err());
}

#[test]
fn add_instance_into_part_or_self_fails() {
    let (mut model, _, part) = part_model();
    assert!(model.add_instance(part, part, Transform::IDENTITY).is_err());
    let assembly = model.new_assembly();
    assert!(model.add_instance(assembly, assembly, Transform::IDENTITY).is_err());
}

#[test]
fn add_instance_cycle_fails() {
    let mut model = Model::new();
    let a = model.new_assembly();
    let b = model.new_assembly();
    model.add_instance(a, b, Transform::IDENTITY).unwrap();
    assert!(model.add_instance(b, a, Transform::IDENTITY).is_err());
}

#[test]
fn removed_instance_does_not_change_roots() {
    let (mut model, _, part) = part_model();
    let assembly = model.new_assembly();
    let instance = model.add_instance(assembly, part, Transform::IDENTITY).unwrap();
    assert_eq!(model.root_products(), vec![assembly]);
    model.remove_instance(instance).unwrap();
    assert!(model.root_products().contains(&assembly));
}

#[test]
fn add_instance_with_removed_product_fails() {
    let (mut model, _, part) = part_model();
    let assembly = model.new_assembly();
    model.remove_node(assembly);
    assert!(model.add_instance(assembly, part, Transform::IDENTITY).is_err());
    let other = model.new_assembly();
    model.remove_node(part);
    assert!(model.add_instance(other, part, Transform::IDENTITY).is_err());
}

#[test]
fn root_products_after_repeated_removals() {
    let mut model = Model::new();
    let a = model.new_assembly();
    let b = model.new_assembly();
    let c = model.new_assembly();
    assert_eq!(model.root_products().len(), 3);
    model.remove_node(b);
    let roots = model.root_products();
    assert_eq!(roots.len(), 2);
    assert!(roots.contains(&a) && roots.contains(&c) && !roots.contains(&b));
    let d = model.new_assembly();
    model.remove_node(a);
    model.remove_node(d);
    assert_eq!(model.root_products(), vec![c]);
}
