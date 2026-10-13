//! BVH_Tree_Test.cxx (BVH_Tree: `math::Bvh`)
//!
//! Not ported: the tests building and editing trees node by node (`AddLeafNode`,
//! `AddInnerNode`, `MinMaxPoints`, `Clear`, `Reserve`, `SetOuterInner`, `Level`,
//! `AddLeafNodeWithBox`, `AddInnerNodeWithBox`, `EstimateSAH`, `NbPrimitives`,
//! `SinglePrimitive`, `Float3DTree`, `Tree2D`, `Tree4D`, `MultipleLeaves`, `DeepTree`,
//! `ModifyPrimitiveIndices`, `ModifyMinMaxPoints`, `ChangeChild`, `NodeInfoBuffer`,
//! `MinMaxPointBuffers`): a `Bvh` is built from boxes and read through `root`, `node_count`,
//! `depth`, `leaf_sizes` and `sah_cost` (tested with the builders).

use kentron::math::Bvh;

#[test]
fn default_constructor() {
    let tree = Bvh::default();
    assert_eq!(tree.node_count(), 0);
    assert_eq!(tree.depth(), 0);
}
