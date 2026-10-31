"""DXF export (CadQuery's `DxfDocument` and `exportDXF`), with ezdxf: lines, circles and
their arcs, ellipses and their arcs, and B-splines (other curves as B-splines), each from
kentron's edge curves as CadQuery takes them from OCCT's; or other curves approximated by
circular arcs, or B-splines by cubic ones (`approx`)."""

from __future__ import annotations

import math
from typing import Any, Dict, Iterable, List, Literal, Optional, Protocol, Tuple, Union, runtime_checkable

import ezdxf
from ezdxf import units, zoom
from ezdxf.entities import factory

from ... import _kentron as native
from ...geom import Plane, Vector
from ..shapes import Compound, Edge, Face, Shape, compound

RAD2DEG = 360.0 / (2 * math.pi)

ApproxOptions = Literal["spline", "arc"]
DxfEntityAttributes = Tuple[Literal["ARC", "CIRCLE", "ELLIPSE", "LINE", "SPLINE"], Dict[str, Any]]


@runtime_checkable
class WorkplaneLike(Protocol):
    @property
    def plane(self) -> Plane:
        ...

    def __iter__(self) -> Iterable[Shape]:
        ...


def _signed_angle(a: Vector, b: Vector, ref: Vector) -> float:
    """The angle from `a` to `b`, positive if their cross product points along `ref`."""
    angle = a.getAngle(b)
    return angle if a.cross(b).dot(ref) >= 0 else -angle


class DxfDocument:
    """A DXF document of CadQuery objects (an ezdxf document, as `document`)."""

    CURVE_TOLERANCE = 1e-9

    def __init__(self, dxfversion: str = "AC1027", setup: Union[bool, List[str]] = False, doc_units: int = units.MM, *, metadata: Union[Dict[str, str], None] = None, approx: Optional[ApproxOptions] = None, tolerance: float = 1e-3):
        self._DISPATCH_MAP = {
            "LINE": self._dxf_line,
            "CIRCLE": self._dxf_circle,
            "ELLIPSE": self._dxf_ellipse,
        }
        self.approx = approx
        self.tolerance = tolerance
        self.document = ezdxf.new(dxfversion=dxfversion, setup=setup, units=doc_units)
        self.msp = self.document.modelspace()
        doc_metadata = self.document.ezdxf_metadata()
        for key, value in (metadata or {}).items():
            doc_metadata[key] = value

    def add_layer(self, name: str, *, color: int = 7, linetype: str = "CONTINUOUS") -> "DxfDocument":
        self.document.layers.add(name, color=color, linetype=linetype)
        return self

    def add_shape(self, shape: Union[WorkplaneLike, Shape], layer: str = "") -> "DxfDocument":
        """Adds a shape's (or a workplane's shapes', in its plane's coordinates) edges."""
        if isinstance(shape, WorkplaneLike):
            shape_ = compound(*shape.__iter__()).transformShape(shape.plane.fG)
        else:
            shape_ = shape
        general = {"layer": layer} if layer else {}
        edges = self._ordered_edges(shape_)
        for edge in edges:
            kind = edge.geomType()
            # (With approx "arc", every curve but lines and circles as circular arcs within
            # the tolerance; with "spline", B-splines of degree above 3 as cubic ones.)
            if self.approx == "arc" and kind not in ("LINE", "CIRCLE"):
                for entity_type, attributes in self._dxf_arcs(edge, self.tolerance):
                    self.msp.add_entity(factory.new(entity_type, dxfattribs={**attributes, **general}))
                continue
            converter = self._DISPATCH_MAP.get(kind)
            if converter:
                entity_type, attributes = converter(edge)
                self.msp.add_entity(factory.new(entity_type, dxfattribs={**attributes, **general}))
            else:
                if self.approx == "spline" and native.edge_bspline(edge.wrapped)[0] > 3:
                    points = [edge.positionAt(i / 64) for i in range(65)]
                    edge = Edge.makeSplineApprox(points, self.tolerance, minDeg=1, maxDeg=3)
                _, attributes = self._dxf_spline(edge, None)
                self.msp.add_spline(dxfattribs=general).apply_construction_tool(ezdxf.math.BSpline(**attributes))
        return self

    @staticmethod
    def _dxf_arcs(edge: Edge, tolerance: float) -> List[DxfEntityAttributes]:
        """`edge` (in the XY plane) as circular arcs (lines where straight), each through its
        piece's ends and middle and within `tolerance` of it, pieces halved until they are."""
        point = lambda t: edge.positionAt(t)
        out: List[DxfEntityAttributes] = []

        def piece(a: float, b: float, depth: int) -> None:
            p0, pm, p1 = point(a), point(0.5 * (a + b)), point(b)
            probes = [point(a + (b - a) * f) for f in (0.25, 0.75)]
            # The circle through the three (in XY): its centre where the chords' bisectors meet.
            ax, ay, bx, by, cx, cy = p0.x, p0.y, pm.x, pm.y, p1.x, p1.y
            d = 2 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by))
            if abs(d) <= 1e-12 * max(1.0, (p1 - p0).Length ** 2):
                if all(_off_line(q, p0, p1) <= tolerance for q in probes) or depth > 20:
                    out.append(("LINE", {"start": p0.toTuple(), "end": p1.toTuple()}))
                    return
            else:
                ux = ((ax * ax + ay * ay) * (by - cy) + (bx * bx + by * by) * (cy - ay) + (cx * cx + cy * cy) * (ay - by)) / d
                uy = ((ax * ax + ay * ay) * (cx - bx) + (bx * bx + by * by) * (ax - cx) + (cx * cx + cy * cy) * (bx - ax)) / d
                r = math.hypot(ax - ux, ay - uy)
                if all(abs(math.hypot(q.x - ux, q.y - uy) - r) <= tolerance for q in probes) or depth > 20:
                    angle = lambda p: math.degrees(math.atan2(p.y - uy, p.x - ux))
                    # (DXF arcs run counter-clockwise: from the end where the turn is one.)
                    ccw = (bx - ax) * (cy - ay) - (by - ay) * (cx - ax) > 0
                    start, end = (angle(p0), angle(p1)) if ccw else (angle(p1), angle(p0))
                    out.append(("ARC", {"center": (ux, uy, p0.z), "radius": r, "start_angle": start, "end_angle": end}))
                    return
            piece(a, 0.5 * (a + b), depth + 1)
            piece(0.5 * (a + b), b, depth + 1)

        piece(0.0, 1.0, 0)
        return out

    @staticmethod
    def _ordered_edges(s: Shape) -> List[Edge]:
        rv: List[Edge] = []
        for w in s.Wires():
            rv.extend(w)
        if isinstance(s, Compound):
            rv.extend(e for e in s if isinstance(e, Edge))
        elif isinstance(s, Edge):
            rv.append(s)
        return rv

    @staticmethod
    def _dxf_line(edge: Edge) -> DxfEntityAttributes:
        return ("LINE", {"start": edge.startPoint().toTuple(), "end": edge.endPoint().toTuple()})

    @staticmethod
    def _dxf_circle(edge: Edge) -> DxfEntityAttributes:
        _, centre, _, ydir, zdir, radius, _ = native.edge_conic(edge.wrapped)
        first, last = edge.bounds()
        ydir, zdir = Vector(ydir), Vector(zdir)
        phi = _signed_angle(ydir, Vector(0, 1, 0), zdir)
        if zdir.z > 0:
            a1, a2 = RAD2DEG * (first - phi), RAD2DEG * (last - phi)
        else:
            a1, a2 = -RAD2DEG * (last - phi) + 180, -RAD2DEG * (first - phi) + 180
        if edge.IsClosed():
            return ("CIRCLE", {"center": tuple(centre), "radius": radius})
        return ("ARC", {"center": tuple(centre), "radius": radius, "start_angle": a1, "end_angle": a2})

    @staticmethod
    def _dxf_ellipse(edge: Edge) -> DxfEntityAttributes:
        _, centre, xdir, _, zdir, major, minor = native.edge_conic(edge.wrapped)
        first, last = edge.bounds()
        if zdir[2] > 0:
            start, end = first, last
        else:
            # (The arc on the ellipse reversed, as CadQuery takes it.)
            start, end = 2 * math.pi - last, 2 * math.pi - first
        return ("ELLIPSE", {
            "center": tuple(centre),
            "major_axis": tuple(major * c for c in xdir),
            "ratio": minor / major,
            "start_param": start,
            "end_param": end,
        })

    @classmethod
    def _dxf_spline(cls, edge: Edge, plane: Optional[Plane]) -> DxfEntityAttributes:
        degree, poles, knots, weights, periodic = native.edge_bspline(edge.wrapped)
        return ("SPLINE", {"control_points": [tuple(p) for p in poles], "order": degree + 1, "knots": list(knots), "weights": weights})


def _off_line(q: Vector, a: Vector, b: Vector) -> float:
    """The distance of `q` from the line through `a` and `b` (in XY)."""
    dx, dy = b.x - a.x, b.y - a.y
    length = math.hypot(dx, dy)
    if length == 0:
        return math.hypot(q.x - a.x, q.y - a.y)
    return abs((q.x - a.x) * dy - (q.y - a.y) * dx) / length


def exportDXF(w, fname: str, approx: Optional[ApproxOptions] = None, tolerance: float = 1e-3, *, doc_units: int = units.MM) -> None:
    """Writes a workplane's (or shapes') edges to a DXF file."""
    dxf = DxfDocument(approx=approx, tolerance=tolerance, doc_units=doc_units)
    if isinstance(w, (WorkplaneLike, Shape)):
        dxf.add_shape(w)
    else:
        for s in w:
            dxf.add_shape(s)
    zoom.extents(dxf.msp)
    dxf.document.saveas(fname)
