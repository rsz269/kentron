//! BRepPrimAPI_MakePrism_Test.cxx

use kentron::math::{Point, Vector};
use kentron::topology::Model;
use kentron::{make, sweep};

/// OCC31294: asking what a shape that is not part of the profile generated must give an
/// empty answer (it used to crash).
#[test]
fn occ31294_generated_list_for_non_base_shape() {
    let mut model = Model::new();
    let vertex = make::vertex(&mut model, Point::ORIGIN);
    let dummy = make::vertex(&mut model, Point::ORIGIN);
    let prism = sweep::extrude(&mut model, vertex, Vector::new(0.0, 0.0, 1.0)).unwrap();

    // the profile vertex generated the prism's edge
    assert_eq!(prism.generated(vertex).len(), 1);
    // a vertex outside the profile generated nothing
    assert_eq!(prism.generated(dummy).len(), 0);
}
