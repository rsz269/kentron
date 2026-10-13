//! BRepGraph_BatchOps_Test.cxx: all-or-nothing batch editing.

use kentron::geometry::surface::{Plane, Surface};
use kentron::math::{Direction, Point, Transform};
use kentron::topology::{
    CompSolid, Compound, CompoundChild, FaceUse, Model, NodeKind, Orientation, Product, Shape,
    Shell, ShellUse, Solid, SolidUse,
};
use kentron::validity;

use crate::graph_support::id_of;

/// A face on the XY plane bounded by a one-edge wire (curveless edge).
fn simple_face(model: &mut Model) -> Shape {
    let v0 = model.add_vertex(Point::ORIGIN, 1.0e-7);
    let v1 = model.add_vertex(Point::new(10.0, 0.0, 0.0), 1.0e-7);
    let edge = model.add_edge(v0, v1, 1.0e-7).unwrap();
    let coedge = model.add_coedge(edge, Orientation::Forward).unwrap();
    let wire = model.add_wire(&[coedge]).unwrap();
    let face = model.new_empty(NodeKind::Face);
    model.set_face_surface(face, Surface::Plane(Plane::new(Point::ORIGIN, Direction::Z)), 1.0e-7).unwrap();
    model.add(face, wire.into()).unwrap();
    face
}

fn shell_of(model: &mut Model, faces: &[Shape]) -> Shape {
    let shell = model.new_empty(NodeKind::Shell);
    model.add_all(shell, faces).unwrap();
    shell
}

fn solid_of(model: &mut Model, shells: &[Shape]) -> Shape {
    let solid = model.new_empty(NodeKind::Solid);
    model.add_all(solid, shells).unwrap();
    solid
}

#[test]
fn shell_ops_append_batch_all_forward() {
    let mut model = Model::new();
    let faces: Vec<Shape> = (0..3).map(|_| simple_face(&mut model)).collect();
    let shell = shell_of(&mut model, &faces);
    let uses = &model[id_of::<Shell>(shell)].faces;
    assert_eq!(uses.len(), 3);
    assert!(uses.iter().all(|&u| model.get(u).is_some()));
}

#[test]
fn shell_ops_append_batch_with_orientations() {
    let mut model = Model::new();
    let f0 = simple_face(&mut model);
    let f1 = simple_face(&mut model);
    let shell = shell_of(&mut model, &[f0, f1.reversed()]);
    let uses = &model[id_of::<Shell>(shell)].faces;
    assert_eq!(uses.len(), 2);
    assert_eq!(model[uses[0]].orientation, Orientation::Forward);
    assert_eq!(model[uses[1]].orientation, Orientation::Reversed);
}

#[test]
fn shell_ops_append_batch_invalid_face_returns_empty() {
    let mut model = Model::new();
    let f0 = simple_face(&mut model);
    let shell = model.new_empty(NodeKind::Shell);
    let missing = Shape::from(kentron::topology::FaceId::from_index(999));
    assert!(model.add_all(shell, &[f0, missing]).is_err());
    assert!(model[id_of::<Shell>(shell)].faces.is_empty(), "nothing is added");
}

#[test]
fn shell_ops_append_batch_invalid_shell_returns_empty() {
    let mut model = Model::new();
    let f0 = simple_face(&mut model);
    let missing = Shape::from(kentron::topology::ShellId::from_index(999));
    assert!(model.add_all(missing, &[f0]).is_err());
}

#[test]
fn shell_ops_remove_faces_all_removed() {
    let mut model = Model::new();
    let faces = [simple_face(&mut model), simple_face(&mut model)];
    let shell = shell_of(&mut model, &faces);
    let uses = model[id_of::<Shell>(shell)].faces.clone();
    model.remove_uses(&uses).unwrap();
    assert!(model[id_of::<Shell>(shell)].faces.is_empty());
}

#[test]
fn shell_ops_remove_faces_invalid_ref_returns_false() {
    let mut model = Model::new();
    let f0 = simple_face(&mut model);
    let shell = shell_of(&mut model, &[f0]);
    let valid = model[id_of::<Shell>(shell)].faces[0];
    assert!(model.remove_uses(&[valid, kentron::topology::Id::<FaceUse>::from_index(999)]).is_err());
    assert_eq!(model[id_of::<Shell>(shell)].faces.len(), 1, "nothing is removed");
}

#[test]
fn solid_ops_append_batch_all_forward() {
    let mut model = Model::new();
    let s0 = { let f = simple_face(&mut model); shell_of(&mut model, &[f]) };
    let s1 = { let f = simple_face(&mut model); shell_of(&mut model, &[f]) };
    let solid = solid_of(&mut model, &[s0, s1]);
    assert_eq!(model[id_of::<Solid>(solid)].shells.len(), 2);
}

#[test]
fn solid_ops_remove_shells_all_removed() {
    let mut model = Model::new();
    let f = simple_face(&mut model);
    let s0 = shell_of(&mut model, &[f]);
    let solid = solid_of(&mut model, &[s0]);
    let uses: Vec<kentron::topology::Id<ShellUse>> = model[id_of::<Solid>(solid)].shells.clone();
    model.remove_uses(&uses).unwrap();
    assert!(model[id_of::<Solid>(solid)].shells.is_empty());
}

#[test]
fn compound_ops_append_batch_all_forward() {
    let mut model = Model::new();
    let faces = [simple_face(&mut model), simple_face(&mut model)];
    let c = model.new_empty(NodeKind::Compound);
    model.add_all(c, &faces).unwrap();
    assert_eq!(model[id_of::<Compound>(c)].children.len(), 2);
}

#[test]
fn compound_ops_remove_children_all_removed() {
    let mut model = Model::new();
    let f = simple_face(&mut model);
    let c = model.new_empty(NodeKind::Compound);
    model.add(c, f).unwrap();
    let uses: Vec<kentron::topology::Id<CompoundChild>> = model[id_of::<Compound>(c)].children.clone();
    model.remove_uses(&uses).unwrap();
    assert!(model[id_of::<Compound>(c)].children.is_empty());
}

#[test]
fn compound_ops_replace_child() {
    let mut model = Model::new();
    let f0 = simple_face(&mut model);
    let f1 = simple_face(&mut model);
    let c = model.new_empty(NodeKind::Compound);
    model.add(c, f0).unwrap();
    model.replace_child(c, f0, f1).unwrap();
    let children: Vec<Shape> = model.children(c).collect();
    assert_eq!(children.len(), 1);
    assert_eq!(children[0].node, f1.node);
}

fn two_solids(model: &mut Model) -> (Shape, Shape) {
    let f0 = simple_face(model);
    let f1 = simple_face(model);
    let s0 = shell_of(model, &[f0]);
    let s1 = shell_of(model, &[f1]);
    (solid_of(model, &[s0]), solid_of(model, &[s1]))
}

#[test]
fn comp_solid_ops_append_batch_all_forward() {
    let mut model = Model::new();
    let (a, b) = two_solids(&mut model);
    let cs = model.new_empty(NodeKind::CompSolid);
    model.add_all(cs, &[a, b]).unwrap();
    assert_eq!(model[id_of::<CompSolid>(cs)].solids.len(), 2);
}

#[test]
fn comp_solid_ops_remove_solids_all_removed() {
    let mut model = Model::new();
    let (a, _) = two_solids(&mut model);
    let cs = model.new_empty(NodeKind::CompSolid);
    model.add(cs, a).unwrap();
    let uses: Vec<kentron::topology::Id<SolidUse>> = model[id_of::<CompSolid>(cs)].solids.clone();
    model.remove_uses(&uses).unwrap();
    assert!(model[id_of::<CompSolid>(cs)].solids.is_empty());
}

#[test]
fn comp_solid_ops_replace_solid() {
    let mut model = Model::new();
    let (a, b) = two_solids(&mut model);
    let cs = model.new_empty(NodeKind::CompSolid);
    model.add(cs, a).unwrap();
    model.replace_child(cs, a, b).unwrap();
    assert_eq!(model.children(cs).map(|s| s.node).collect::<Vec<_>>(), [b.node]);
}

#[test]
fn product_ops_append_batch() {
    let mut model = Model::new();
    let (a, b) = two_solids(&mut model);
    let part_a = model.new_part(a).unwrap();
    let part_b = model.new_part(b).unwrap();
    let assembly = model.new_assembly();
    let instances = model
        .add_instances(assembly, &[(part_a, Transform::IDENTITY), (part_b, Transform::IDENTITY)])
        .unwrap();
    assert_eq!(instances.len(), 2);
    assert!(instances.iter().all(|&i| model.get(i).is_some()));
}

#[test]
fn product_ops_remove_occurrences_all_removed() {
    let mut model = Model::new();
    let (a, _) = two_solids(&mut model);
    let part = model.new_part(a).unwrap();
    let assembly = model.new_assembly();
    let instance = model.add_instance(assembly, part, Transform::IDENTITY).unwrap();
    model.remove_instance(instance).unwrap();
    let product: &Product = &model[assembly];
    assert!(product.children.is_empty());
}

#[test]
fn validation_after_batch_operations() {
    let mut model = Model::new();
    let faces: Vec<Shape> = (0..3).map(|_| simple_face(&mut model)).collect();
    let shell = shell_of(&mut model, &faces);
    let uses = model[id_of::<Shell>(shell)].faces.clone();
    assert_eq!(uses.len(), 3);
    model.remove_uses(&uses[..2]).unwrap();
    assert_eq!(model[id_of::<Shell>(shell)].faces.len(), 1);
    assert!(validity::check_model_structure(&model).is_valid());
}
