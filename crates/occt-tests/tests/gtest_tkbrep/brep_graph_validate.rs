//! BRepGraph_Validate_Test.cxx: model consistency checks.
//!
//! OCCT's `Lightweight` / `Audit` modes are `validity::check_model_structure` /
//! `check_model`. Corruptions are made by editing records directly (their fields are
//! public), as OCCT does through mutation guards.
//!
//! Not ported: `DocumentRootReferencedByOccurrence_Detected` (kentron derives root products
//! instead of storing them), `Audit_DetectsSharedOwnedUse` (slots shared by several owners
//! are legal: deduplication makes them), `RemoveNodeMaintainsActiveCount` and
//! `Audit_ValidatesCoEdgeUIDs*` (internal counters and UID tables).

use kentron::geometry::ParametricCurve;
use std::f64::consts::PI;

use kentron::geometry::curve::{Circle, Curve, Line};
use kentron::geometry::surface::{Cylinder, Plane, Sphere, Surface};
use kentron::math::{Direction, Frame, Point, Transform, Vector};
use kentron::topology::{
    CoEdge, Edge, EdgeId, Face, FaceId, Id, Instance, Model, NodeId, NodeKind, Orientation,
    ShellId, SolidId, VertexId, Wire, WireId,
};
use kentron::validity::{self, ProblemKind, Severity};

use crate::graph_support::{box_model, compound, first, id_of};

fn audit(model: &Model) -> validity::Report {
    validity::check_model(model)
}

fn light(model: &Model) -> validity::Report {
    validity::check_model_structure(model)
}

/// A vertex that is added and removed again.
fn removed_vertex(model: &mut Model) -> VertexId {
    let v = model.add_vertex(Point::new(100.0, 0.0, 0.0), 1.0e-7);
    model.remove_node(v);
    v
}

#[test]
fn clean_model_no_problems() {
    let (model, _) = box_model();
    for report in [audit(&model), light(&model)] {
        assert!(report.is_valid());
        assert_eq!(report.errors().count(), 0);
        assert_eq!(report.warnings().count(), 0);
    }
}

#[test]
fn after_geometry_deduplicate_no_problems() {
    let (source, solid) = box_model();
    let face = crate::first(&source, solid, NodeKind::Face);
    let mut model = Model::new();
    let a = model.import(&source, face).unwrap();
    let b = model.import(&source, face).unwrap();
    compound(&mut model, &[a, b]);
    model.deduplicate();
    let report = audit(&model);
    assert!(report.is_valid());
    assert_eq!(report.errors().count(), 0);
}

#[test]
fn detects_removed_node_reference() {
    let (mut model, _) = box_model();
    let gone = removed_vertex(&mut model);
    let edge = first::<Edge>(&model);
    let start = model[edge].start;
    model[start].vertex = gone;
    for report in [light(&model), audit(&model)] {
        assert!(!report.is_valid());
        assert!(report.has(ProblemKind::DanglingReference, Severity::Error));
    }
}

#[test]
fn wire_connectivity_disconnected_edges_is_a_warning() {
    let (mut model, _) = box_model();
    let (wire, _) = model.table::<Wire>().iter().find(|(_, w)| w.coedges.len() >= 2).unwrap();
    let edge = model[model[wire].coedges[0]].edge;
    let (start, end) = model.edge_end_vertices(edge);
    let (other, _) = model
        .table::<kentron::topology::Vertex>()
        .iter()
        .find(|&(v, _)| v != start && v != end)
        .unwrap();
    let end_use = model[edge].end;
    model[end_use].vertex = other;
    let report = audit(&model);
    assert!(!report.has(ProblemKind::DisconnectedWire, Severity::Error));
    assert!(report.has(ProblemKind::DisconnectedWire, Severity::Warning));
}

#[test]
fn orphan_curve_slot_is_an_error() {
    let (mut model, _) = box_model();
    let edge = first::<Edge>(&model);
    assert!(model.get(model[edge].curve).is_some());
    model[edge].curve = Id::INVALID;
    let report = audit(&model);
    assert!(!report.is_valid());
    assert!(report.has(ProblemKind::OrphanSlot, Severity::Error));
}

#[test]
fn after_split_edge_audit_stays_clean() {
    let (mut model, _) = box_model();
    let edges = model.slot_count(NodeKind::Edge);
    let (edge, _) = model
        .table::<Edge>()
        .iter()
        .find(|&(e, _)| !model.is_degenerate(e) && model.edge_curve(e).is_some())
        .unwrap();
    let (first, last) = model.edge_range(edge);
    let mid = 0.5 * (first + last);
    let (curve, _, _) = model.edge_curve(edge).unwrap();
    let vertex = model.add_vertex(curve.point(mid), model.tolerance(edge));
    let (a, b) = model.split_edge(edge, vertex, mid).unwrap();
    assert_eq!(model.slot_count(NodeKind::Edge), edges + 2);
    assert!(model.get(a).is_some() && model.get(b).is_some());
    assert!(model.get(edge).is_none());
    assert!(audit(&model).is_valid());
}

#[test]
fn coedge_face_out_of_bounds_detected() {
    let (mut model, _) = box_model();
    let coedge = first::<CoEdge>(&model);
    model[coedge].face = FaceId::from_index(model.slot_count(NodeKind::Face) as u32 + 999);
    assert!(!light(&model).is_valid());
    assert!(!audit(&model).is_valid());
}

#[test]
fn assembly_valid_product_no_problems() {
    let (mut model, solid) = box_model();
    let part = model.new_part(solid).unwrap();
    let assembly = model.new_assembly();
    model.add_instance(assembly, part, Transform::IDENTITY).unwrap();
    model.add_instance(assembly, part, Transform::translation(Vector::new(20.0, 0.0, 0.0))).unwrap();
    assert_eq!(model[assembly].children.len(), 2);
    assert!(audit(&model).is_valid());
}

fn instance_model() -> (Model, Id<Instance>) {
    let (mut model, solid) = box_model();
    let part = model.new_part(solid).unwrap();
    let assembly = model.new_assembly();
    let instance = model.add_instance(assembly, part, Transform::IDENTITY).unwrap();
    (model, instance)
}

#[test]
fn instance_target_out_of_bounds_detected() {
    for kind in [NodeKind::Solid, NodeKind::Product] {
        let (mut model, instance) = instance_model();
        model[instance].target = NodeId { kind, index: model.slot_count(kind) as u32 + 999 };
        let report = audit(&model);
        assert!(!report.is_valid());
        assert!(report.has(ProblemKind::DanglingReference, Severity::Error));
    }
}

#[test]
fn instance_targeting_an_instance_detected() {
    let (mut model, instance) = instance_model();
    model[instance].target = instance.into();
    let report = audit(&model);
    assert!(!report.is_valid());
    assert!(report.errors().next().is_some());
}

#[test]
fn removed_vertex_reference_found_by_both_modes() {
    let (mut model, _) = box_model();
    let gone = removed_vertex(&mut model);
    let start = model[first::<Edge>(&model)].start;
    model[start].vertex = gone;
    assert!(!light(&model).is_valid());
    let report = audit(&model);
    assert!(!report.is_valid());
    assert!(report.has(ProblemKind::DanglingReference, Severity::Error));
}

#[test]
fn face_with_surface_but_no_wires_is_a_warning() {
    let (mut model, _) = box_model();
    let face = first::<Face>(&model);
    let wire_use = *model[face].wires.last().unwrap();
    model.remove_use(wire_use);
    let report = audit(&model);
    assert!(report.is_valid());
    let warning = report.problems.first().expect("expected a warning");
    assert_eq!(warning.severity, Severity::Warning);
    assert_eq!(warning.kind, ProblemKind::FaceWithoutBoundary);
}

// ---- models built from low-level records ----

fn line_edge(model: &mut Model, a: VertexId, b: VertexId) -> EdgeId {
    let (pa, pb) = (model[a].point, model[b].point);
    let edge = model.add_edge(a, b, 1.0e-7).unwrap();
    let direction = Direction::from_vector(pb - pa);
    model.set_edge_curve(edge.into(), Curve::Line(Line::new(pa, direction)), 1.0e-7).unwrap();
    model.set_edge_range(edge.into(), 0.0, pa.distance(pb)).unwrap();
    edge
}

fn curve_edge(model: &mut Model, a: VertexId, b: VertexId, curve: Curve, last: f64) -> EdgeId {
    let edge = model.add_edge(a, b, 1.0e-7).unwrap();
    model.set_edge_curve(edge.into(), curve, 1.0e-7).unwrap();
    model.set_edge_range(edge.into(), 0.0, last).unwrap();
    edge
}

fn wire_of(model: &mut Model, coedges: &[(EdgeId, Orientation)]) -> WireId {
    let ids: Vec<_> = coedges.iter().map(|&(e, o)| model.add_coedge(e, o).unwrap()).collect();
    model.add_wire(&ids).unwrap()
}

fn face_of(model: &mut Model, surface: Surface, wire: WireId) -> FaceId {
    let face = model.new_empty(NodeKind::Face);
    model.set_face_surface(face, surface, 1.0e-7).unwrap();
    model.add(face, wire.into()).unwrap();
    id_of(face)
}

fn solid_of(model: &mut Model, faces: &[FaceId]) -> (ShellId, SolidId) {
    let shell = model.new_empty(NodeKind::Shell);
    for &f in faces {
        model.add(shell, f.into()).unwrap();
    }
    let solid = model.new_empty(NodeKind::Solid);
    model.add(solid, shell).unwrap();
    (id_of(shell), id_of(solid))
}

fn plane(origin: Point, normal: Direction) -> Surface {
    Surface::Plane(Plane::new(origin, normal))
}

#[test]
fn synthetic_box_audit_clean() {
    use Orientation::{Forward as F, Reversed as R};
    let mut m = Model::new();
    let p = [
        (0.0, 0.0, 0.0),
        (10.0, 0.0, 0.0),
        (10.0, 20.0, 0.0),
        (0.0, 20.0, 0.0),
        (0.0, 0.0, 30.0),
        (10.0, 0.0, 30.0),
        (10.0, 20.0, 30.0),
        (0.0, 20.0, 30.0),
    ];
    let v: Vec<VertexId> = p.iter().map(|&(x, y, z)| m.add_vertex(Point::new(x, y, z), 1.0e-7)).collect();
    let pairs = [(0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4), (0, 4), (1, 5), (2, 6), (3, 7)];
    let e: Vec<EdgeId> = pairs.iter().map(|&(a, b)| line_edge(&mut m, v[a], v[b])).collect();
    let wires = [
        [(e[0], F), (e[1], F), (e[2], F), (e[3], F)],
        [(e[4], F), (e[5], F), (e[6], F), (e[7], F)],
        [(e[0], F), (e[9], F), (e[4], R), (e[8], R)],
        [(e[1], F), (e[10], F), (e[5], R), (e[9], R)],
        [(e[2], F), (e[11], F), (e[6], R), (e[10], R)],
        [(e[3], F), (e[8], F), (e[7], R), (e[11], R)],
    ];
    let planes = [
        plane(Point::ORIGIN, Direction::NEG_Z),
        plane(Point::new(0.0, 0.0, 30.0), Direction::Z),
        plane(Point::ORIGIN, Direction::NEG_Y),
        plane(Point::new(10.0, 0.0, 0.0), Direction::X),
        plane(Point::new(10.0, 20.0, 0.0), Direction::Y),
        plane(Point::new(0.0, 20.0, 0.0), Direction::NEG_X),
    ];
    let faces: Vec<FaceId> = wires
        .iter()
        .zip(planes)
        .map(|(w, s)| {
            let wire = wire_of(&mut m, w);
            face_of(&mut m, s, wire)
        })
        .collect();
    solid_of(&mut m, &faces);
    let report = audit(&m);
    assert!(report.is_valid());
    assert_eq!(report.errors().count(), 0);
}

#[test]
fn synthetic_cylinder_audit_clean() {
    use Orientation::{Forward as F, Reversed as R};
    let mut m = Model::new();
    let bottom = m.add_vertex(Point::new(5.0, 0.0, 0.0), 1.0e-7);
    let top = m.add_vertex(Point::new(5.0, 0.0, 15.0), 1.0e-7);
    let circle = |z: f64| Curve::Circle(Circle::new(Frame::from_z(Point::new(0.0, 0.0, z), Direction::Z), 5.0));
    let e_bottom = curve_edge(&mut m, bottom, bottom, circle(0.0), 2.0 * PI);
    let e_top = curve_edge(&mut m, top, top, circle(15.0), 2.0 * PI);
    let seam = line_edge(&mut m, bottom, top);
    let w_bottom = wire_of(&mut m, &[(e_bottom, F), (seam, R)]);
    let w_top = wire_of(&mut m, &[(e_top, F), (seam, F)]);
    let w_side = wire_of(&mut m, &[(seam, F), (seam, R)]);
    let f_bottom = face_of(&mut m, plane(Point::ORIGIN, Direction::NEG_Z), w_bottom);
    let f_top = face_of(&mut m, plane(Point::new(0.0, 0.0, 15.0), Direction::Z), w_top);
    let side = Surface::Cylinder(Cylinder::new(Frame::WORLD, 5.0));
    let f_side = face_of(&mut m, side, w_side);
    solid_of(&mut m, &[f_bottom, f_top, f_side]);
    assert!(audit(&m).is_valid());
}

#[test]
fn synthetic_sphere_audit_clean() {
    use Orientation::{Forward as F, Reversed as R};
    let mut m = Model::new();
    let v = m.add_vertex(Point::new(8.0, 0.0, 0.0), 1.0e-7);
    let circle = Curve::Circle(Circle::new(Frame::from_z(Point::ORIGIN, Direction::Y), 8.0));
    let seam = curve_edge(&mut m, v, v, circle, 2.0 * PI);
    let wire = wire_of(&mut m, &[(seam, F), (seam, R)]);
    let face = face_of(&mut m, Surface::Sphere(Sphere::new(Frame::WORLD, 8.0)), wire);
    solid_of(&mut m, &[face]);
    assert!(audit(&m).is_valid());
}

/// A triangle face in its own solid, with corners `a`, `b`, `c`.
fn triangle_solid(m: &mut Model, a: Point, b: Point, c: Point) -> SolidId {
    let v: Vec<VertexId> = [a, b, c].iter().map(|&p| m.add_vertex(p, 1.0e-7)).collect();
    let e: Vec<EdgeId> = (0..3).map(|i| line_edge(m, v[i], v[(i + 1) % 3])).collect();
    let wire = wire_of(m, &[(e[0], Orientation::Forward), (e[1], Orientation::Forward), (e[2], Orientation::Forward)]);
    let face = face_of(m, plane(Point::ORIGIN, Direction::Z), wire);
    solid_of(m, &[face]).1
}

#[test]
fn synthetic_compound_audit_clean() {
    let mut m = Model::new();
    let a = triangle_solid(&mut m, Point::ORIGIN, Point::new(10.0, 0.0, 0.0), Point::new(5.0, 10.0, 0.0));
    let b = triangle_solid(
        &mut m,
        Point::new(20.0, 0.0, 0.0),
        Point::new(30.0, 0.0, 0.0),
        Point::new(25.0, 10.0, 0.0),
    );
    compound(&mut m, &[a.into(), b.into()]);
    assert!(audit(&m).is_valid());
}

#[test]
fn synthetic_assembly_audit_clean() {
    let mut m = Model::new();
    let v: Vec<VertexId> = [(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)]
        .iter()
        .map(|&(x, y)| m.add_vertex(Point::new(x, y, 0.0), 1.0e-7))
        .collect();
    let e: Vec<(EdgeId, Orientation)> =
        (0..4).map(|i| (line_edge(&mut m, v[i], v[(i + 1) % 4]), Orientation::Forward)).collect();
    let wire = wire_of(&mut m, &e);
    let face = face_of(&mut m, plane(Point::ORIGIN, Direction::Z), wire);
    let (_, solid) = solid_of(&mut m, &[face]);
    let part = m.new_part(solid).unwrap();
    let assembly = m.new_assembly();
    m.add_instance(assembly, part, Transform::IDENTITY).unwrap();
    m.add_instance(assembly, part, Transform::translation(Vector::new(20.0, 0.0, 0.0))).unwrap();
    assert!(audit(&m).is_valid());
}
