"""pykentron against CadQuery: the same constructions measured by both.

Run with a Python that has CadQuery installed (any test needs it, to compare with):

    python -m pytest crates/pykentron/tests      # or: python crates/pykentron/tests/test_against_cadquery.py

Each case is a function of the module (`cq`): it returns a shape, workplane or number, and
both libraries' answers are compared (volume, area, bounding box, centre, counts).
"""

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "python"))
import pykentron  # noqa: E402

try:
    import pytest

    cadquery = pytest.importorskip("cadquery")
    parametrize = pytest.mark.parametrize
except ImportError:  # run as a script (python test_against_cadquery.py)
    import cadquery

    def parametrize(name, values, ids=None):
        return lambda f: f


def box(cq):
    return cq.Solid.makeBox(10, 20, 30, cq.Vector(1, 2, 3))


def cylinder_along_y(cq):
    return cq.Solid.makeCylinder(5, 15, (1, 0, 0), (0, 1, 0))


def cone(cq):
    return cq.Solid.makeCone(5, 2, 10, (0, 0, 0), (1, 0, 0))


def workplane_box(cq):
    return cq.Workplane("XY").box(10, 4, 3).translate((1, 2, 3)).val()


def fused_boxes(cq):
    a = cq.Solid.makeBox(10, 10, 10)
    return a.fuse(cq.Solid.makeBox(10, 10, 10, cq.Vector(5, 5, 5)))


def cut_bores(cq):
    a = cq.Solid.makeBox(20, 6, 20)
    return a.cut(*(cq.Solid.makeCylinder(1.7, 10, (x, -2, z), (0, 1, 0)) for x in (5, 15) for z in (5, 15)))


def intersected(cq):
    return cq.Solid.makeBox(10, 10, 10).intersect(cq.Solid.makeCylinder(6, 20, (5, 5, -5)))


def polyline_prism_xz(cq):
    return cq.Workplane("XZ").polyline([(0, 0), (10, 0), (8, 5), (2, 5)]).close().extrude(4).val()


def polyline_prism_custom_plane(cq):
    plane = cq.Plane((1, 2, 3), (0, 1, 0), (1, 0, 0))
    return cq.Workplane(plane).polyline([(0, 0), (6, 0), (6, 2), (3, 4), (0, 2)]).close().extrude(3).val()


def hexagon(cq):
    return cq.Workplane("XY").polygon(6, 10).extrude(2).val()


def face_with_hole(cq):
    outer = cq.Wire.makePolygon([(0, 0, 0), (20, 0, 0), (20, 10, 0), (0, 10, 0)], close=True)
    inner = cq.Wire.makePolygon([(5, 3, 0), (8, 3, 0), (8, 6, 0), (5, 6, 0)], close=True)
    face = cq.Face.makeFromWires(outer, [inner])
    return cq.Solid.extrudeLinear(face, cq.Vector(0, 0, 4))


def loft(cq):
    w1 = cq.Wire.makePolygon([(0, 0, 0), (10, 0, 0), (10, 10, 0), (0, 10, 0)], close=True)
    w2 = cq.Wire.makePolygon([(2, 2, 8), (8, 2, 8), (8, 8, 8), (2, 8, 8)], close=True)
    return cq.Solid.makeLoft([w1, w2], ruled=True)


def rotated_mirrored(cq):
    a = cq.Solid.makeBox(10, 2, 3).rotate((0, 0, 0), (0, 0, 1), 30).translate((5, 0, 0))
    return a.mirror("YZ")


def revolved(cq):
    return cq.Workplane("XZ").polyline([(2, 0), (4, 0), (4, 3), (2, 3)]).close().revolve(360, (0, 0, 0), (0, 1, 0)).val()


def located(cq):
    return cq.Solid.makeBox(1, 2, 3).moved(cq.Location(cq.Vector(1, 2, 3), cq.Vector(0, 0, 1), 90))


def euler_location(cq):
    return cq.Solid.makeBox(1, 2, 3).moved(cq.Location(1, 2, 3, 10, 20, 30))


def offset_outline(cq):
    w = cq.Wire.makePolygon([(0, 0, 0), (10, 0, 0), (10, 6, 0), (0, 6, 0)], close=True)
    return w.offset2D(1.5, "arc")[0]


def smoothed_helix_spline(cq):
    # (What Workplane.parametricCurve makes: a smoothed approximation.)
    at = lambda t: cq.Vector(5 * math.cos(6 * math.pi * t), 5 * math.sin(6 * math.pi * t), 4 * t)
    return cq.Edge.makeSplineApprox([at(i / 400) for i in range(401)], tol=1e-6, smoothing=(1, 1, 1), maxDeg=6)


def chamfered_arc_plate(cq):
    # (The arc's circle comes back into the plate near its front corners.)
    plate = cq.Workplane("XY").lineTo(100, 0).lineTo(100, 76.2).threePointArc((50, 101.6), (0, 76.2)).close().extrude(10).val()
    arc = [e for e in plate.Edges() if e.geomType() == "CIRCLE" and e.Center().z > 5]
    return plate.chamfer(3, None, arc)


def asymmetric_chamfer_of_an_extrusion(cq):
    # (The first distance along the edge's first face: of an extrusion, its side, at both
    # ends.)
    rod = cq.Workplane("XY").circle(1.5).extrude(10).val()
    rims = [e for e in rod.Edges() if e.geomType() == "CIRCLE"]
    return rod.chamfer(1.0, 0.5, rims)


def mirrored_outline(cq):
    # (Two halves of an outline joined along the mirror line: both wind one way, and the
    # line out and back between them is dropped.)
    half = cq.Wire.makePolygon([(0, 0, 0), (0, 6, 0), (10, 6, 0), (22, 0, 0)], close=True)
    wire = cq.Wire.combine([half, half.mirror("XZ")])[0]
    return cq.Solid.extrudeLinear(cq.Face.makeFromWires(wire), cq.Vector(0, 0, 4))


def fillet_with_bore_edges(cq):
    # (The bore's vertical edges, where its cylinder meets itself, are left out.)
    block = cq.Solid.makeBox(20, 20, 10).cut(cq.Solid.makeCylinder(4, 30, cq.Vector(10, 10, -10)))
    vertical = [e for e in block.Edges() if abs(e.startPoint().z - e.endPoint().z) > 1]
    return block.fillet(2, vertical)


def loft_turned_aside(cq):
    # (Outlines in planes square to each other: matched by their vertices, not by normals.)
    a = cq.Wire.makePolygon([(0, 0, 0), (10, 0, 0), (12, 6, 0), (0, 8, 0)], close=True)
    b = cq.Wire.makePolygon([(30, 0, 30), (30, 0, 20), (30, 6, 18), (30, 8, 30)], close=True)
    return cq.Solid.makeLoft([a, b])


def sweep_round_hexagon(cq):
    # (A closed spine with corners, the profile mitred at each, as CadQuery's default.)
    r = 3 / math.sin(math.pi / 3)
    pts = [(r * math.cos(i * math.pi / 3 + math.pi / 6), r * math.sin(i * math.pi / 3 + math.pi / 6), 0) for i in range(7)]
    spine = cq.Wire.makePolygon(pts)
    profile = cq.Face.makeFromWires(cq.Wire.makeCircle(0.5, cq.Vector(3, 0, 0), cq.Vector(0, 1, 0)))
    return cq.Solid.sweep(profile, spine, True, False, None, "right")


def frenet_sweep_round_helix(cq):
    helix = cq.Wire.makeHelix(pitch=0.4, height=2.4, radius=0.5)
    square = cq.Wire.makePolygon([(0.35, 0, 0.1), (0.5, 0, 0.05), (0.5, 0, 0.35), (0.35, 0, 0.3)], close=True)
    return cq.Solid.sweep(cq.Face.makeFromWires(square), helix, True, True)


def tapered_split_square(cq):
    # (Sides split in line: the taper's offset joins them.)
    pts = [(0, 0, 0), (5, 0, 0), (10, 0, 0), (10, 5, 0), (10, 10, 0), (5, 10, 0), (0, 10, 0), (0, 5, 0)]
    face = cq.Face.makeFromWires(cq.Wire.makePolygon(pts, close=True))
    return cq.Solid.extrudeLinear(face, cq.Vector(0, 0, 1), taper=25)


def rotated_then_located(cq):
    # (A rotation moves the geometry; a location after it places that, not the original.)
    turned = cq.Solid.makeBox(4, 2, 1).translate((1, 0, 0)).rotate((0, 0, 0), (0, 0, 1), 90)
    return turned.located(cq.Location(cq.Vector(10, 0, 0)))


def spline_surface_through_a_grid(cq):
    at = lambda x, y: cq.Vector(x, y, 0.3 * x * y - 0.1 * x * x + math.sin(0.7 * y))
    grid = [[at(i * 0.2, j * 0.25) for j in range(20)] for i in range(30)]
    return cq.Face.makeSplineApprox(grid, tol=1e-3, minDeg=3, maxDeg=8)


def plane_of_negative_size(cq):
    return cq.Face.makePlane(-2.0, 3.0)


def subclassed_wire(cq):
    # (Constructors make the class they are called on, as scripts subclass to add methods.)
    class Outline(cq.Wire):
        def face(self):
            return cq.Face.makeFromWires(self)

    return Outline.makePolygon([(0, 0, 0), (3, 0, 0), (3, 2, 0)], close=True).face()


def array_of_vectors(cq):
    # (numpy keeps vectors whole, as objects: one per point.)
    import numpy

    points = numpy.array([cq.Vector(0, 0, 0), cq.Vector(4, 0, 0), cq.Vector(4, 3, 0)])
    return cq.Face.makeFromWires(cq.Wire.makePolygon(list(points[[0, 1, 2]]) * (points.size == 3), close=True))


def sweep_through_three_circles(cq):
    # (Sections blended smoothly, not linearly between neighbours; the ends' caps facing out.)
    spine = cq.Wire.assembleEdges([cq.Edge.makeLine(cq.Vector(-10, 0, 0), cq.Vector(10, 0, 0))])
    circles = [cq.Wire.makeCircle(r, cq.Vector(x, 0, 0), cq.Vector(1, 0, 0)) for x, r in ((-10, 2), (0, 1), (10, 2))]
    return cq.Solid.sweep_multi(circles, spine)


def sweep_from_circle_to_square(cq):
    # (Sections of different numbers of edges, matched as a loft matches them.)
    spine = cq.Wire.assembleEdges([cq.Edge.makeLine(cq.Vector(-10, 0, 0), cq.Vector(10, 0, 0))])
    circle = cq.Wire.makeCircle(1, cq.Vector(-10, 0, 0), cq.Vector(1, 0, 0))
    square = cq.Wire.makePolygon([(10, -1, -1), (10, 1, -1), (10, 1, 1), (10, -1, 1)], close=True)
    return cq.Solid.sweep_multi([circle, square], spine)


def fuse_of_many(cq):
    # (Fused pairwise, then pairs of pairs: as one union of them all.)
    cubes = [cq.Solid.makeBox(2, 2, 2, cq.Vector(1.5 * i, (i % 3) * 0.5, 0)) for i in range(25)]
    return cubes[0].fuse(*cubes[1:])


def rounded_pocket_rim(cq):
    # (The rim's corners are concave on the top: the fillets meet there in mitres.)
    block = cq.Solid.makeBox(20, 20, 6).cut(cq.Solid.makeBox(10, 8, 4, cq.Vector(5, 6, 3)))
    rim = [e for e in block.Edges() if abs(e.Center().z - 6) < 1e-9 and 1 < e.Center().x < 19 and 1 < e.Center().y < 19]
    return block.fillet(1, rim)


def fillet_of_a_mitre(cq):
    # (Where fillets of a box's two sides meet, their mitre goes on smoothly from the box's
    # edge between those sides (filleted with it), and the fillet narrows to nothing where
    # both reach the top.)
    block = cq.Solid.makeBox(10, 10, 10)
    top = [e for e in block.Edges() if abs(e.Center().z - 10) < 1e-9 and (abs(e.Center().x - 10) < 1e-9 or abs(e.Center().y - 10) < 1e-9)]
    block = block.fillet(3, top)
    return block.fillet(1, [e for e in block.Edges() if e.geomType() == "ELLIPSE"])


def rounded_recess_cut(cq):
    # (A tapered recess with rounded corners and bottom: the bottom's fillets go on from the
    # straight sides round the corners' fillets; the recess cut from a block.)
    def rect(w, d, z):
        return cq.Wire.makePolygon([cq.Vector(x, y, z) for x, y in ((-w / 2, -d / 2), (w / 2, -d / 2), (w / 2, d / 2), (-w / 2, d / 2))], close=True)

    recess = cq.Solid.makeLoft([rect(16, 56, 3), rect(30, 70, 12)], True)
    recess = recess.fillet(6, [e for e in recess.Edges() if abs(e.Center().x) > 1 and abs(e.Center().y) > 1 and e.Center().z > 3.5])
    recess = recess.fillet(4, [e for e in recess.Edges() if abs(e.Center().z - 3) < 1e-6])
    return cq.Solid.makeBox(40, 80, 12, cq.Vector(-20, -40, 0)).cut(recess)


def filled_square_through_a_point(cq):
    # (A plate spanning a square, bent through a point off its plane: the least bending.)
    corners = [cq.Vector(0, 0, 0), cq.Vector(0, 10, 0), cq.Vector(0, 10, 10), cq.Vector(0, 0, 10)]
    edges = cq.Wire.makePolygon(corners, close=True).Edges()
    return cq.Face.makeNSidedSurface(edges, [cq.Vector(5, 5, 5)])


CASES = [
    box, cylinder_along_y, cone, workplane_box, fused_boxes, cut_bores, intersected,
    polyline_prism_xz, polyline_prism_custom_plane, hexagon, face_with_hole, loft,
    rotated_mirrored, revolved, located, euler_location, offset_outline, smoothed_helix_spline,
    chamfered_arc_plate, mirrored_outline, fillet_with_bore_edges,
    loft_turned_aside, sweep_round_hexagon, frenet_sweep_round_helix, tapered_split_square,
    rotated_then_located, spline_surface_through_a_grid, plane_of_negative_size, subclassed_wire,
    array_of_vectors, sweep_through_three_circles, sweep_from_circle_to_square,
    fuse_of_many, rounded_pocket_rim, fillet_of_a_mitre, rounded_recess_cut, filled_square_through_a_point,
    asymmetric_chamfer_of_an_extrusion,
]


def measures(shape):
    if hasattr(shape, "val"):
        shape = shape.val()
    bb = shape.BoundingBox()
    kind = shape.ShapeType() if callable(shape.ShapeType) else shape.ShapeType
    m = {"bbox": (bb.xmin, bb.ymin, bb.zmin, bb.xmax, bb.ymax, bb.zmax), "solids": len(shape.Solids())}
    if kind in ("Solid", "Compound", "CompSolid"):
        m["volume"] = shape.Volume()
    if kind in ("Wire", "Edge"):
        m["length"] = shape.Length() if hasattr(shape, "Length") else None
    else:
        m["area"] = shape.Area()
    c = shape.Center()
    m["center"] = (c.x, c.y, c.z)
    return m


@parametrize("case", CASES, ids=[c.__name__ for c in CASES])
def test_same_as_cadquery(case):
    ours, theirs = measures(case(pykentron)), measures(case(cadquery))
    for key, expected in theirs.items():
        got = ours[key]
        if isinstance(expected, tuple):
            # Bounding boxes: CadQuery's are enlarged by tolerances; allow 1e-3.
            assert all(math.isclose(a, b, rel_tol=1e-4, abs_tol=2e-3) for a, b in zip(got, expected)), f"{key}: {got} != {expected}"
        else:
            assert math.isclose(got, expected, rel_tol=1e-4, abs_tol=1e-6), f"{key}: {got} != {expected}"


if __name__ == "__main__":
    import traceback

    failed = 0
    for case in CASES:
        try:
            test_same_as_cadquery(case)
            print(f"ok      {case.__name__}")
        except Exception as e:
            failed += 1
            print(f"FAILED  {case.__name__}: {type(e).__name__}: {e}")
    print(f"{len(CASES) - failed} passed, {failed} failed")
    sys.exit(1 if failed else 0)
