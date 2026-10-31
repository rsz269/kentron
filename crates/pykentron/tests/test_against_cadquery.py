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


def rounded_corners_between_arcs_and_lines(cq):
    # (A D-shaped outline: its corners between the arc and the lines rounded along the arc,
    # the one between the two lines as the others.)
    V = cq.Vector
    edges = [cq.Edge.makeLine(V(0, 0, 0), V(20, 0, 0)), cq.Edge.makeLine(V(20, 0, 0), V(20, 10, 0)),
             cq.Edge.makeThreePointArc(V(20, 10, 0), V(10, 18, 0), V(0, 10, 0)), cq.Edge.makeLine(V(0, 10, 0), V(0, 0, 0))]
    face = cq.Face.makeFromWires(cq.Wire.assembleEdges(edges))
    return face.fillet2D(2.0, face.Vertices())


def loft_through_three_unevenly_spaced_circles(cq):
    # (CadQuery's loft, as a workplane's, passes through its sections evenly, however far
    # apart they are.)
    import importlib

    V = cq.Vector
    circles = [cq.Wire.makeCircle(r, V(0, 0, z), V(0, 0, 1)) for z, r in ((0, 5), (1, 6), (10, 5))]
    return importlib.import_module(cq.__name__ + ".occ_impl.shapes").loft(circles, cap=True)


def vertices_of_a_face_placed_twice(cq):
    # (One face placed at two points: told apart by their hashes, as a sketch tells apart
    # what it selects, the vertices of both (a small cube at each).)
    face = cq.Face.makePlane(1, 2)
    faces = [face.moved(cq.Location(cq.Vector(x, 0, 0))) for x in (-3, 3)]
    vertices = {hash(v): v for f in faces for v in f.Vertices()}.values()
    return cq.Compound.makeCompound([cq.Solid.makeBox(0.1, 0.1, 0.1, v.Center()) for v in vertices])


def text_in_fonts_not_installed_and_cff_ones(cq):
    # (A family not installed set as OCCT sets it, in a generic one; a font with CFF
    # outlines, Nimbus Sans, standing as high as OCCT sets it, aligned at the top.)
    texts = [cq.Compound.makeText("Hag", 10, 1, font=f, halign="left", valign="top").translate(cq.Vector(0, 20 * i, 0))
             for i, f in enumerate(("Arial", "Serif", "Courier", "Nimbus Sans"))]
    return cq.Compound.makeCompound(texts)


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


def reeded_prisms_under_their_lids(cq):
    # (A prism ringed by round ribs, a rounded lid standing on it, its sides at the ribs'
    # crests: each rib's arc meets a side a hair either way of tangent, crossing it twice
    # close together, and the slivers between are faces of the union too, closed.)
    def outline_point(w, d, r, s):
        fx, fy = w / 2 - r, d / 2 - r
        flats = [((0, d / 2), (1, 0), fx), ((w / 2, fy), (0, -1), 2 * fy), ((fx, -d / 2), (-1, 0), 2 * fx), ((-w / 2, -fy), (0, 1), 2 * fy), ((-fx, d / 2), (1, 0), fx)]
        arcs = [((fx, fy), 90), ((fx, -fy), 0), ((-fx, -fy), -90), ((-fx, fy), 180)]
        for i, ((x, y), (dx, dy), length) in enumerate(flats):
            if s <= length or i == 4:
                return (x + dx * s, y + dy * s), (-dy, dx)
            s -= length
            (cx, cy), a0 = arcs[i]
            if s <= math.pi * r / 2:
                a = math.radians(a0) - s / r
                return (cx + r * math.cos(a), cy + r * math.sin(a)), (math.cos(a), math.sin(a))
            s -= math.pi * r / 2

    def module(w, d, z0, lip, z1, rib=0.8):
        r = min(12.0, min(w, d) / 2 - 0.5)
        n = 4 * round((2 * (w + d) - 8 * r + 2 * math.pi * r) / 12)
        step = (2 * (w + d) - 8 * r + 2 * math.pi * r) / n
        grooves = [cq.Vector(*outline_point(w, d, r, k * step)[0]) for k in range(n)]
        edges = []
        for k in range(n):
            (x, y), (ax, ay) = outline_point(w, d, r, (k + 0.5) * step)
            edges.append(cq.Edge.makeThreePointArc(grooves[k], cq.Vector(x + ax * rib, y + ay * rib), grooves[(k + 1) % n]))
        reeds = cq.Solid.extrudeLinear(cq.Face.makeFromWires(cq.Wire.assembleEdges(edges)), cq.Vector(0, 0, lip - z0)).translate(cq.Vector(0, 0, z0))
        lid = cq.Workplane("XY").workplane(offset=lip).rect(w + 2 * rib, d + 2 * rib).extrude(z1 - lip).edges("|Z").fillet(r + rib).val()
        return reeds.fuse(lid)

    return cq.Compound.makeCompound([module(48.7, 48.7, 14.9, 31.7, 38.7), module(99.5, 48.7, 14.9, 100.0, 107.0).translate(cq.Vector(0, 80, 0))])


def cut_of_overlapping_solids(cq):
    # (A compound's solids that overlap are split by one another, what they share a solid of
    # its own; the solid apart from them is left as it is.)
    boxes = [cq.Solid.makeBox(10, 10, 10, cq.Vector(x, 0, 0)) for x in (0, 8, 30)]
    return cq.Compound.makeCompound(boxes).cut(cq.Solid.makeCylinder(2, 20, cq.Vector(4, 5, -5)))


def hollowed_turned_box(cq):
    # (A placed solid hollowed to a closed shell: its own shell, placed, bounds the cavity.)
    return cq.Solid.makeBox(12, 4.5, 26.2).rotate(cq.Vector(), cq.Vector(0, 0, 1), 30).hollow([], 0.1)


def corners_of_a_placed_cutter_rounded(cq):
    # (A disc notched by a square placed across its rim, the square's corners in the disc then
    # rounded, as a sketch rounds the corners it tagged: found among the cut face's corners,
    # though the cut made them anew where it baked the square's placement.)
    disc = cq.Face.makeFromWires(cq.Wire.makeCircle(3, cq.Vector(), cq.Vector(0, 0, 1)))
    square = cq.Face.makePlane(1, 1).moved(cq.Location(cq.Vector(2.8, 0, 0)))
    face = disc.cut(square).Faces()[0]
    corners = [v for v in square.Vertices() if v in face.Vertices()]
    assert len(corners) == 2
    return cq.Solid.extrudeLinear(face.fillet2D(0.2, corners), cq.Vector(0, 0, 1))


CASES = [
    box, cylinder_along_y, cone, workplane_box, fused_boxes, cut_bores, intersected,
    polyline_prism_xz, polyline_prism_custom_plane, hexagon, face_with_hole, loft,
    rotated_mirrored, revolved, located, euler_location, offset_outline, smoothed_helix_spline,
    chamfered_arc_plate, mirrored_outline, fillet_with_bore_edges,
    loft_turned_aside, sweep_round_hexagon, frenet_sweep_round_helix, tapered_split_square,
    rotated_then_located, spline_surface_through_a_grid, plane_of_negative_size, subclassed_wire,
    array_of_vectors, sweep_through_three_circles, sweep_from_circle_to_square,
    fuse_of_many, rounded_pocket_rim, fillet_of_a_mitre, rounded_recess_cut, filled_square_through_a_point,
    asymmetric_chamfer_of_an_extrusion, rounded_corners_between_arcs_and_lines,
    loft_through_three_unevenly_spaced_circles, vertices_of_a_face_placed_twice,
    text_in_fonts_not_installed_and_cff_ones, reeded_prisms_under_their_lids,
    cut_of_overlapping_solids, hollowed_turned_box, corners_of_a_placed_cutter_rounded,
]


def measures(shape):
    if hasattr(shape, "val"):
        shape = shape.val()
    bb = shape.BoundingBox()
    kind = shape.ShapeType() if callable(shape.ShapeType) else shape.ShapeType
    m = {"bbox": (bb.xmin, bb.ymin, bb.zmin, bb.xmax, bb.ymax, bb.zmax), "solids": len(shape.Solids()), "valid": shape.isValid()}
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


def _gltf_scene(path):
    """A glTF file's nodes (name, translation, rotation), materials' colours, and the volume
    its meshes enclose, placed."""
    import json
    import struct

    raw = Path(path).read_bytes()
    if raw[:4] == b"glTF":
        n = struct.unpack_from("<I", raw, 12)[0]
        doc = json.loads(raw[20:20 + n])
        data = raw[20 + n + 8:]
    else:
        doc = json.loads(raw)
        data = (Path(path).parent / doc["buffers"][0]["uri"]).read_bytes()

    def accessor(i):
        a = doc["accessors"][i]
        v = doc["bufferViews"][a["bufferView"]]
        size = 3 if a["type"] == "VEC3" else 1
        values = struct.unpack_from("<%d%s" % (size * a["count"], {5126: "f", 5125: "I", 5123: "H"}[a["componentType"]]), data, v["byteOffset"] + a.get("byteOffset", 0))
        return [values[k:k + 3] for k in range(0, len(values), 3)] if size == 3 else list(values)

    def placement(node):
        x, y, z, w = node.get("rotation", [0, 0, 0, 1])
        t = node.get("translation", [0, 0, 0])
        r = [[1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)], [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)], [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)]]
        return [r[i] + [t[i]] for i in range(3)]

    def times(a, b):
        return [[sum(a[i][k] * b[k][j] for k in range(3)) for j in range(3)] + [sum(a[i][k] * b[k][3] for k in range(3)) + a[i][3]] for i in range(3)]

    volume, nodes = 0.0, []

    def walk(i, m):
        nonlocal volume
        node = doc["nodes"][i]
        m = times(m, placement(node))
        nodes.append((node.get("name"), node.get("translation", [0, 0, 0]), node.get("rotation", [0, 0, 0, 1])))
        for primitive in doc["meshes"][node["mesh"]]["primitives"] if "mesh" in node else []:
            points = [[sum(m[r][k] * p[k] for k in range(3)) + m[r][3] for r in range(3)] for p in accessor(primitive["attributes"]["POSITION"])]
            indices = accessor(primitive["indices"])
            for k in range(0, len(indices), 3):
                a, b, c = (points[j] for j in indices[k:k + 3])
                volume += (a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) + a[2] * (b[0] * c[1] - b[1] * c[0])) / 6
        for child in node.get("children", []):
            walk(child, m)

    walk(doc["scenes"][0]["nodes"][0], [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0]])
    return nodes, [m["pbrMetallicRoughness"]["baseColorFactor"] for m in doc.get("materials", [])], volume


@parametrize("suffix", [".gltf", ".glb"])
def test_gltf_assembly_as_cadquery(tmp_path, suffix):
    # (A placed assembly, coloured parts and a placed sub-assembly: the same nodes, placed
    # alike, turned Z up to Y up, the same materials, and the meshes enclosing the parts.)
    def write(cq, path):
        sub = cq.Assembly(name="sub", loc=cq.Location((0, 30, 0), (0, 0, 1), 45)).add(cq.Solid.makeSphere(4), name="ball", color=cq.Color("green"))
        assy = cq.Assembly(name="top", loc=cq.Location((1, 2, 3)))
        assy.add(cq.Solid.makeBox(10, 20, 30), name="box", color=cq.Color("red"))
        assy.add(cq.Solid.makeCylinder(3, 5), name="cyl", loc=cq.Location((20, 0, 0), (1, 0, 0), 30), color=cq.Color(0, 0, 1, 0.5))
        assy.add(sub)
        assy.export(str(path))
        return _gltf_scene(path)

    ours = write(pykentron, tmp_path / f"ours{suffix}")
    theirs = write(cadquery, tmp_path / f"theirs{suffix}")
    assert [n[0] for n in ours[0]] == [n[0] for n in theirs[0]]
    for (_, t1, r1), (_, t2, r2) in zip(ours[0], theirs[0]):
        assert all(math.isclose(a, b, abs_tol=1e-9) for a, b in zip(t1 + r1, t2 + r2)), (t1, r1, t2, r2)
    assert ours[1] == theirs[1]
    # (Meshes of the same parts, finer or coarser.)
    assert math.isclose(ours[2], theirs[2], rel_tol=1e-3), (ours[2], theirs[2])


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
    import tempfile

    for suffix in (".gltf", ".glb"):
        with tempfile.TemporaryDirectory() as folder:
            try:
                test_gltf_assembly_as_cadquery(Path(folder), suffix)
                print(f"ok      gltf_assembly{suffix}")
            except Exception as e:
                failed += 1
                print(f"FAILED  gltf_assembly{suffix}: {type(e).__name__}: {e}")
    print(f"{len(CASES) + 2 - failed} passed, {failed} failed")
    sys.exit(1 if failed else 0)
