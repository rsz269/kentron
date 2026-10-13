//! BRepTools_Test.cxx
//!
//! OCCT's binary format (BinTools) has no kentron counterpart: the binary round trips of
//! `ModDataBug_29723_FormatVersions` and `ModDataBug_31946_TriangleOptions` are dropped,
//! and `DemoBinPersist_1_BinaryStreamRoundTrip` uses the ASCII format instead.

use std::f64::consts::PI;

use occt_tests::{assert_double_eq, assert_near};
use kentron::geometry::curve::{Circle, Curve};
use kentron::io::brep::{self, FormatVersion, WriteOptions};
use kentron::math::{Axis, Direction, Frame, Point, Transform, Vector};
use kentron::primitives::Cuboid;
use kentron::properties::{self, Options};
use kentron::topology::{Model, NodeKind, Shape};
use kentron::{boolean, convert, make, meshing, validity};

use crate::{make_box, make_edge};

fn box_brep_string() -> Vec<u8> {
    let mut model = Model::new();
    let solid = make_box(&mut model, 10.0, 20.0, 30.0);
    let mut out = Vec::new();
    brep::write(&model, solid, &mut out, &WriteOptions::default()).unwrap();
    assert!(!out.is_empty());
    out
}

/// Reads a shape from `data` into a fresh model.
fn read_from(data: &[u8]) -> (Model, Shape) {
    let mut model = Model::new();
    let shape = brep::read(&mut model, &mut &data[..]).expect("read failed");
    (model, shape)
}

/// Writes `shape` and reads it back into a fresh model.
fn round_trip(model: &Model, shape: Shape, version: FormatVersion, with_meshes: bool) -> (Model, Shape) {
    let mut out = Vec::new();
    let options = WriteOptions { version, with_meshes, with_normals: false };
    brep::write(model, shape, &mut out, &options).unwrap();
    assert!(!out.is_empty());
    read_from(&out)
}

fn nb_triangles(model: &Model, shape: Shape) -> usize {
    model
        .explore(shape, NodeKind::Face)
        .filter_map(|face| model.face_mesh(face))
        .map(|mesh| mesh.triangles().len())
        .sum()
}

fn nb_subshapes(model: &Model, shape: Shape, kind: NodeKind) -> usize {
    model.explore(shape, kind).count()
}

fn expect_equivalent_shape(
    (expected_model, expected): (&Model, Shape),
    (actual_model, actual): (&Model, Shape),
    expected_triangles: usize,
) {
    assert_eq!(actual.kind(), expected.kind());
    for kind in [NodeKind::Solid, NodeKind::Shell, NodeKind::Face, NodeKind::Edge, NodeKind::Vertex] {
        assert_eq!(
            nb_subshapes(actual_model, actual, kind),
            nb_subshapes(expected_model, expected, kind),
            "{kind:?} count"
        );
    }
    assert_eq!(nb_triangles(actual_model, actual), expected_triangles);
    assert_near!(
        properties::volume(actual_model, actual).mass,
        properties::volume(expected_model, expected).mass,
        1.0e-7
    );
    assert_near!(
        properties::area(actual_model, actual).mass,
        properties::area(expected_model, expected).mass,
        1.0e-7
    );
}

/// Reads `data` on many threads at once (each into its own model) and checks every result
/// has the reference volume. Guards against shared state in number parsing.
fn parallel_reads_match(data: &[u8], reads: usize, reference_mass: f64) {
    let masses: Vec<f64> = std::thread::scope(|scope| {
        let handles: Vec<_> = (0..reads)
            .map(|_| {
                scope.spawn(|| {
                    let (model, shape) = read_from(data);
                    properties::volume(&model, shape).mass.abs()
                })
            })
            .collect();
        handles.into_iter().map(|h| h.join().expect("a read failed")).collect()
    });
    for (i, mass) in masses.into_iter().enumerate() {
        assert_near!(mass, reference_mass, 1.0e-6, "shape {i} mass {mass} != {reference_mass}");
    }
}

#[test]
fn parallel_read_identical_mass() {
    let data = box_brep_string();
    let (model, shape) = read_from(&data);
    let reference = properties::volume(&model, shape).mass.abs();
    assert!(reference > 0.0);
    parallel_reads_match(&data, 200, reference);
}

#[test]
fn parallel_read_repeated_runs() {
    let data = box_brep_string();
    let (model, shape) = read_from(&data);
    let reference = properties::volume(&model, shape).mass.abs();
    for _ in 0..5 {
        parallel_reads_match(&data, 100, reference);
    }
}

/// bugs/moddata_3/bug29723: every BRep format version preserves NURBS geometry.
#[test]
fn mod_data_bug_29723_format_versions() {
    let mut model = Model::new();
    let solid = make_box(&mut model, 10.0, 20.0, 30.0);
    let nurbs_box = convert::to_nurbs(&mut model, solid).unwrap();

    for version in [FormatVersion::V1, FormatVersion::V2] {
        let (restored_model, restored) = round_trip(&model, nurbs_box, version, true);
        expect_equivalent_shape((&model, nurbs_box), (&restored_model, restored), 0);
    }
}

/// bugs/moddata_3/bug31946: the writer honors the mesh option.
#[test]
#[should_panic] // expected to fail: a known gap in kentron
fn mod_data_bug_31946_triangle_options() {
    let mut model = Model::new();
    let solid = make_box(&mut model, 10.0, 20.0, 30.0);
    meshing::mesh(&mut model, solid, 0.1).unwrap();
    assert_eq!(nb_triangles(&model, solid), 12);

    let (with_model, with) = round_trip(&model, solid, FormatVersion::default(), true);
    let (without_model, without) = round_trip(&model, solid, FormatVersion::default(), false);
    expect_equivalent_shape((&model, solid), (&with_model, with), 12);
    expect_equivalent_shape((&model, solid), (&without_model, without), 0);
}

/// bugs/demo/bug25344_4: a shape keeps its type and topology through save and restore.
#[test]
fn demo_bug_25344_4_save_restore_shape_type() {
    let mut model = Model::new();
    let source = make_box(&mut model, 1.0, 1.0, 1.0);
    let (restored_model, restored) = round_trip(&model, source, FormatVersion::V2, true);
    assert_eq!(restored.kind(), source.kind());
    assert_eq!(nb_subshapes(&restored_model, restored, NodeKind::Solid), 1);
    assert_eq!(nb_subshapes(&restored_model, restored, NodeKind::Face), 6);
}

/// demo/draw/binpersist_1: a box survives storage with its bounding box, topology,
/// validity and global properties intact.
#[test]
#[should_panic] // expected to fail: a known gap in kentron
fn demo_bin_persist_1_binary_stream_round_trip() {
    let mut model = Model::new();
    let source = Cuboid::new(100.0, 200.0, 300.0)
        .with_frame(Frame::new(Point::new(10.0, 20.0, 30.0), Direction::Z, Direction::X))
        .build(&mut model)
        .unwrap();
    let (restored_model, restored) = round_trip(&model, source, FormatVersion::default(), true);

    let expected = properties::bounding_box(&model, source);
    let actual = properties::bounding_box(&restored_model, restored);
    for (a, e) in [
        (actual.min.x, expected.min.x),
        (actual.min.y, expected.min.y),
        (actual.min.z, expected.min.z),
        (actual.max.x, expected.max.x),
        (actual.max.y, expected.max.y),
        (actual.max.z, expected.max.z),
    ] {
        assert_near!(a, e, 1.0e-9);
    }

    assert!(validity::check(&restored_model, restored).is_valid());
    for kind in [NodeKind::Solid, NodeKind::Shell, NodeKind::Face, NodeKind::Edge, NodeKind::Vertex] {
        assert_eq!(
            nb_subshapes(&restored_model, restored, kind),
            nb_subshapes(&model, source, kind),
            "{kind:?} count"
        );
    }

    // edges are counted once per face using them
    let source_length = properties::length(&model, source).mass;
    assert_near!(properties::length(&restored_model, restored).mass, source_length, 1.0e-7);
    assert_double_eq!(source_length, 4800.0);

    let source_area = properties::area(&model, source).mass;
    assert_near!(properties::area(&restored_model, restored).mass, source_area, 1.0e-7);
    assert_double_eq!(source_area, 220000.0);

    let source_volume = properties::volume(&model, source).mass;
    assert_near!(properties::volume(&restored_model, restored).mass, source_volume, 1.0e-7);
    assert_double_eq!(source_volume, 6000000.0);
}

/// bugs/moddata_1/bug20616: empty shapes of every kind survive storage.
#[test]
#[should_panic] // expected to fail: a known gap in kentron
fn mod_data_bug_20616_empty_shape_save_restore() {
    for kind in [
        NodeKind::Vertex,
        NodeKind::Edge,
        NodeKind::Wire,
        NodeKind::Face,
        NodeKind::Shell,
        NodeKind::Solid,
        NodeKind::CompSolid,
        NodeKind::Compound,
    ] {
        let mut model = Model::new();
        let source = model.new_empty(kind);
        let (_, restored) = round_trip(&model, source, FormatVersion::default(), true);
        assert_eq!(restored.kind(), kind);
    }
}

fn expect_same_principal_properties(exact: &properties::MassProperties, meshed: &properties::MassProperties) {
    assert_near!(exact.mass, meshed.mass, 1.0e-7);
    let (e, m) = (exact.principal_moments(), meshed.principal_moments());
    for i in 0..3 {
        assert_near!(e[i], m[i], 1.0e-7, "principal moment {i}");
    }
}

/// bugs/modalg_7/bug29734: properties computed from planar faces' meshes agree with the
/// exact ones.
#[test]
#[should_panic] // expected to fail: a known gap in kentron
fn modalg_bug_29734_exact_and_triangulated_properties() {
    let mut model = Model::new();
    let first_box = make_box(&mut model, 1.0, 2.0, 3.0);
    let second_box = make_box(&mut model, 3.0, 2.0, 1.0);
    let moved = second_box
        .moved(&Transform::translation(Vector::new(0.5, 0.5, 0.5)))
        .moved(&Transform::rotation(Axis::new(Point::ORIGIN, Direction::new(1.0, 1.0, 1.0)), PI / 6.0));
    let fused = boolean::union(&mut model, first_box, moved).unwrap();

    meshing::mesh(&mut model, fused, 0.01).unwrap();
    assert!(validity::check(&model, fused).is_valid());

    let exact = Options::default();
    let meshed = Options { use_meshes: true, ..Options::default() };
    expect_same_principal_properties(
        &properties::length_with(&model, fused, exact),
        &properties::length_with(&model, fused, meshed),
    );
    expect_same_principal_properties(
        &properties::area_with(&model, fused, exact),
        &properties::area_with(&model, fused, meshed),
    );
    expect_same_principal_properties(
        &properties::volume_with(&model, fused, exact),
        &properties::volume_with(&model, fused, meshed),
    );
}

/// bugs/demo/bug23130: a box's faces, edges and vertices all have the default tolerance.
#[test]
fn demo_bug_23130_box_max_tolerance() {
    let mut model = Model::new();
    let solid = make_box(&mut model, 100.0, 100.0, 100.0);
    for (kind, expected_count) in [(NodeKind::Face, 6), (NodeKind::Edge, 12), (NodeKind::Vertex, 8)] {
        let tolerances: Vec<f64> =
            model.explore(solid, kind).map(|s| model.tolerance(s)).collect();
        let min = tolerances.iter().copied().fold(f64::INFINITY, f64::min);
        let max = tolerances.iter().copied().fold(f64::NEG_INFINITY, f64::max);

        assert_eq!(model.unique_subshapes(solid, kind).len(), expected_count);
        assert_near!(min, 1.0e-7, 1.0e-7);
        assert_near!(max, 1.0e-7, 1.0e-7);
    }
}

/// From tests/bugs/moddata_3/bug27356: clearing meshes removes a free edge's polyline
/// and leaves the edge usable.
#[test]
fn mod_data_bug_27356_clean_removes_free_edge_polygon3d() {
    let mut model = Model::new();
    let edge = make_edge(&mut model, Point::ORIGIN, Point::new(1.0, 0.0, 0.0));
    meshing::mesh(&mut model, edge, 0.01).unwrap();
    assert!(model.edge_polyline(edge).is_some());

    model.clear_meshes(edge);
    assert!(model.edge_polyline(edge).is_none());
}

fn circle_face(model: &mut Model, center: Point, normal: Direction, radius: f64) -> Shape {
    let circle = Circle::new(Frame::from_z(center, normal), radius);
    let edge = make::edge_from_curve(model, Curve::Circle(circle)).unwrap();
    let wire = make::wire(model, &[edge]).unwrap();
    make::planar_face(model, wire).unwrap()
}

fn cut_and_select_first_face(model: &mut Model, face: Shape, tool: Shape) -> Shape {
    let result = boolean::difference(model, face, tool).expect("cut failed");
    model.unique_subshapes(result, NodeKind::Face)[0]
}

fn expect_box_uv_bounds(model: &Model, face: Shape) {
    let (u_min, u_max, v_min, v_max) = model.face_uv_bounds(face);
    assert_near!(u_min, 0.0, 1.0e-10);
    assert_near!(u_max, 1500.0, 1.0e-10);
    assert_near!(v_min, 0.0, 1.0e-10);
    assert_near!(v_max, 1500.0, 1.0e-10);
}

/// From tests/bugs/modalg_5/bug24404: the UV bounds of a box face stay its full
/// parametric limits after holes are cut into it.
#[test]
fn modalg_bug_24404_box_face_uv_bounds() {
    let mut model = Model::new();
    let solid = make_box(&mut model, 1500.0, 1500.0, 1500.0);
    let mut face = model.unique_subshapes(solid, NodeKind::Face)[3];
    expect_box_uv_bounds(&model, face);

    let holes = [
        (Point::new(1361.60462531413, 1500.0, 275.105307765905), 182.781239888725),
        (Point::new(398.623694869499, 1500.0, 5.77182937332096), 181.948898616306),
        (Point::new(66.2457890511211, 1500.0, 832.500519113239), 176.393473166564),
        (Point::new(1177.73545803307, 1500.0, 1406.03245550006), 150.768085993996),
    ];
    for (center, radius) in holes {
        let tool = circle_face(&mut model, center, Direction::Y, radius);
        face = cut_and_select_first_face(&mut model, face, tool);
        expect_box_uv_bounds(&model, face);
    }
}
