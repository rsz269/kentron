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


CASES = [
    box, cylinder_along_y, cone, workplane_box, fused_boxes, cut_bores, intersected,
    polyline_prism_xz, polyline_prism_custom_plane, hexagon, face_with_hole, loft,
    rotated_mirrored, revolved, located, euler_location, offset_outline,
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
