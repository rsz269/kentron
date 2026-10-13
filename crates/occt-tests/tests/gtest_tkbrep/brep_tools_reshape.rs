//! BRepTools_ReShape_Test.cxx (BRepTools_ReShape is `Substitution`)

use kentron::math::{Direction, Frame, Point};
use kentron::primitives::Cuboid;
use kentron::topology::{Model, NodeKind, Shape, Substitution};
use kentron::make;

fn vertex(model: &mut Model, x: f64) -> Shape {
    make::vertex(model, Point::new(x, 0.0, 0.0))
}

fn edge(model: &mut Model, v1: Shape, v2: Shape) -> Shape {
    make::edge(model, v1, v2).unwrap()
}

fn compound(model: &mut Model, children: &[Shape]) -> Shape {
    let c = model.new_empty(NodeKind::Compound);
    for &child in children {
        model.add(c, child).unwrap();
    }
    c
}

fn same(a: Option<Shape>, b: Shape) -> bool {
    a.is_some_and(|a| a.is_same(&b))
}

/// `get` gives the direct replacement (one hop); `resolve` follows the chain to its end.
#[test]
fn value_leaf_follows_chain_to_leaf() {
    let mut model = Model::new();
    let [a, b, c] = [0.0, 1.0, 2.0].map(|x| vertex(&mut model, x));
    let mut subst = Substitution::new();
    subst.replace(a, b);
    subst.replace(b, c);

    assert!(same(subst.get(a), b), "get gives only the direct replacement");
    assert!(same(subst.resolve(a), c), "resolve walks the whole chain");
    assert!(same(subst.resolve(b), c));
    assert!(same(subst.resolve(c), c), "the end of a chain is a fixpoint");
}

#[test]
fn value_leaf_unrecorded_shape_is_identity() {
    let mut model = Model::new();
    let a = vertex(&mut model, 0.0);
    assert!(same(Substitution::new().resolve(a), a));
}

/// A chain ending in a removal resolves to nothing.
#[test]
fn value_leaf_chain_ending_in_remove_returns_null() {
    let mut model = Model::new();
    let [a, b] = [0.0, 1.0].map(|x| vertex(&mut model, x));
    let mut subst = Substitution::new();
    subst.replace(a, b);
    subst.remove(b);
    assert!(subst.resolve(a).is_none());
}

/// a -> b, b -> a: both records are kept, and applying terminates.
#[test]
fn replace_direct_cycle_handled_by_apply() {
    let mut model = Model::new();
    let [a, b] = [0.0, 1.0].map(|x| vertex(&mut model, x));
    let mut subst = Substitution::new();
    subst.replace(a, b);
    subst.replace(b, a);

    assert!(same(subst.get(a), b));
    assert!(same(subst.get(b), a));
    assert!(subst.apply(&mut model, a).is_some(), "apply must terminate with a shape");
}

/// a -> b -> c -> a: accepted when recorded, broken when applied.
#[test]
fn replace_longer_cycle_handled_by_apply() {
    let mut model = Model::new();
    let [a, b, c] = [0.0, 1.0, 2.0].map(|x| vertex(&mut model, x));
    let mut subst = Substitution::new();
    subst.replace(a, b);
    subst.replace(b, c);
    subst.replace(c, a);

    assert!(same(subst.get(a), b));
    assert!(same(subst.get(b), c));
    assert!(same(subst.get(c), a));
    assert!(subst.apply(&mut model, a).is_some(), "apply must terminate with a shape");
}

/// Replacing a shape by a compound containing it must not recurse forever.
#[test]
fn apply_handles_structural_containment_without_crash() {
    let mut model = Model::new();
    let [v1, v2, v3] = [0.0, 1.0, 2.0].map(|x| vertex(&mut model, x));
    let e = edge(&mut model, v1, v2);
    let container = compound(&mut model, &[e, v3]);

    let mut subst = Substitution::new();
    subst.replace(e, container); // e -> compound{e, v3}
    assert!(subst.apply(&mut model, e).is_some());
}

/// A sub-shape reached along two paths (two edges sharing a vertex) is replaced on both.
#[test]
fn apply_diamond_sharing_is_processed_correctly() {
    let mut model = Model::new();
    let [v1, v2, v3, v2_new] = [0.0, 1.0, 2.0, 1.5].map(|x| vertex(&mut model, x));
    let e1 = edge(&mut model, v1, v2);
    let e2 = edge(&mut model, v2, v3);
    let pair = compound(&mut model, &[e1, e2]);

    let mut subst = Substitution::new();
    subst.replace(v2, v2_new);
    let result = subst.apply(&mut model, pair).expect("result");

    let mut replaced = 0;
    for v in model.explore(result, NodeKind::Vertex) {
        assert_ne!(v.node, v2.node, "v2 must be replaced on every path");
        if v.node == v2_new.node {
            replaced += 1;
        }
    }
    assert!(replaced >= 2, "the replacement must appear on both paths");
}

#[test]
fn clear_drops_all_bindings() {
    let mut model = Model::new();
    let [a, b] = [0.0, 1.0].map(|x| vertex(&mut model, x));
    let mut subst = Substitution::new();
    subst.replace(a, b);
    assert!(same(subst.get(a), b));

    subst.clear();
    assert!(same(subst.get(a), a));
    assert!(same(subst.resolve(a), a));
}

/// Replacing the same shape again overwrites the earlier record.
#[test]
fn replace_last_wins() {
    let mut model = Model::new();
    let [a, b, c] = [0.0, 1.0, 2.0].map(|x| vertex(&mut model, x));
    let mut subst = Substitution::new();
    subst.replace(a, b);
    subst.replace(a, c);
    assert!(same(subst.get(a), c));
    assert!(same(subst.resolve(a), c));
}

#[test]
fn replace_self_is_no_op() {
    let mut model = Model::new();
    let a = vertex(&mut model, 0.0);
    let mut subst = Substitution::new();
    subst.replace(a, a);
    assert!(!subst.is_recorded(a));
    assert!(same(subst.get(a), a));
}

/// A removal drops the sub-shape from the rebuilt parent.
#[test]
fn apply_remove_drops_sub_shape_from_parent() {
    let mut model = Model::new();
    let [v1, v2, v3] = [0.0, 1.0, 2.0].map(|x| vertex(&mut model, x));
    let parent = compound(&mut model, &[v1, v2, v3]);

    let mut subst = Substitution::new();
    subst.remove(v2);
    let result = subst.apply(&mut model, parent).expect("result");

    let children: Vec<_> = model.children(result).collect();
    assert!(children.iter().all(|c| c.node != v2.node), "v2 should be removed");
    assert_eq!(children.len(), 2);
}

/// An edge shared by two parents is replaced in both.
#[test]
fn apply_shared_edge_across_two_parents() {
    let mut model = Model::new();
    let [v1, v2, v3, v_new] = [0.0, 1.0, 2.0, 0.5].map(|x| vertex(&mut model, x));
    let e1 = edge(&mut model, v1, v2);
    let e2 = edge(&mut model, v2, v3);
    let e1_new = edge(&mut model, v_new, v2);
    let parent_a = compound(&mut model, &[e1, v3]);
    let parent_b = compound(&mut model, &[e1, e2]);

    let mut subst = Substitution::new();
    subst.replace(e1, e1_new);
    let result_a = subst.apply(&mut model, parent_a).expect("result a");
    let result_b = subst.apply(&mut model, parent_b).expect("result b");

    for result in [result_a, result_b] {
        let edges: Vec<_> = model.explore(result, NodeKind::Edge).collect();
        assert!(edges.iter().all(|e| e.node != e1.node));
        assert!(edges.iter().any(|e| e.node == e1_new.node));
    }
}

/// Applying without records returns the shape itself.
#[test]
fn apply_no_bindings_is_identity() {
    let mut model = Model::new();
    let [v1, v2] = [0.0, 1.0].map(|x| vertex(&mut model, x));
    let e = edge(&mut model, v1, v2);
    let result = Substitution::new().apply(&mut model, e).expect("result");
    assert!(result.is_same(&e));
}

/// A replacement containing the replaced shape two levels down must not recurse forever.
#[test]
fn apply_deep_structural_containment_without_crash() {
    let mut model = Model::new();
    let [v1, v2, extra] = [0.0, 1.0, 3.0].map(|x| vertex(&mut model, x));
    let e = edge(&mut model, v1, v2);
    let inner = compound(&mut model, &[e]);
    let outer = compound(&mut model, &[inner, extra]);

    let mut subst = Substitution::new();
    subst.replace(e, outer); // e -> outer{inner{e}, extra}
    assert!(subst.apply(&mut model, e).is_some());
}

#[test]
fn value_leaf_considers_location_mode() {
    let mut model = Model::new();
    let [a, b] = [0.0, 1.0].map(|x| vertex(&mut model, x));
    let mut subst = Substitution::new().matching_placement();
    subst.replace(a, b);
    assert!(same(subst.resolve(a), b));
}

/// From tests/bugs/moddata_3/bug29662: replacing a sub-compound keeps the untouched
/// branch's solids and adds the replacement's.
#[test]
fn mod_data_bug_29662_replace_compound() {
    let mut model = Model::new();
    let box_at = |model: &mut Model, x: f64, z: f64| {
        Cuboid::new(10.0, 10.0, 10.0)
            .with_frame(Frame::new(Point::new(x, 0.0, z), Direction::Z, Direction::X))
            .build(model)
            .unwrap()
    };
    let box1 = box_at(&mut model, 0.0, 0.0);
    let box2 = box_at(&mut model, 10.0, 0.0);
    let box3 = box_at(&mut model, 0.0, 10.0);
    let box4 = box_at(&mut model, 10.0, 10.0);
    let box5 = box_at(&mut model, 20.0, 0.0);

    let c1 = compound(&mut model, &[box1, box2]);
    let c2 = compound(&mut model, &[box3, box4]);
    let c = compound(&mut model, &[c1, c2]);
    let c3 = compound(&mut model, &[box5]);

    let mut subst = Substitution::new();
    subst.replace(c2, c3);
    let result = subst.apply_until(&mut model, c, NodeKind::Compound).expect("result");
    assert_eq!(model.explore(result, NodeKind::Solid).count(), 3);
}
