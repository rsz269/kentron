"""DXF import (CadQuery's): a DXF file's layers' curves joined into wires, the wires into
faces, with ezdxf."""

from __future__ import annotations

import math
from collections import OrderedDict
from typing import List

from ... import _kentron as native
from ...geom import Vector
from ..shapes import Edge, Face, Shape, _edges_to_wires, sortWiresByBuildOrder

RAD2DEG = 360.0 / (2 * math.pi)


def _dxf_line(el):
    try:
        return (Edge.makeLine(Vector(el.dxf.start.xyz), Vector(el.dxf.end.xyz)),)
    except Exception:
        return ()


def _dxf_circle(el):
    try:
        return (Edge.makeCircle(el.dxf.radius, Vector(el.dxf.center.xyz)),)
    except Exception:
        return ()


def _dxf_arc(el):
    try:
        return (Edge.makeCircle(el.dxf.radius, Vector(el.dxf.center.xyz), angle1=el.dxf.start_angle, angle2=el.dxf.end_angle),)
    except Exception:
        return ()


def _dxf_polyline(el):
    rv = (DXF_CONVERTERS[e.dxf.dxftype](e) for e in el.virtual_entities())
    return (e[0] for e in rv if e)


def _dxf_spline(el):
    try:
        knots = OrderedDict()
        for k in el.knots:
            knots[k] = knots.get(k, 0) + 1
        weights = list(el.weights) if el.weights else None
        poles = [tuple(p) for p in el.control_points]
        handle = native.make_bspline(poles, weights, list(knots), list(knots.values()), el.dxf.degree, bool(el.closed))
        return (Shape.cast(handle),)
    except Exception:
        return ()


def _dxf_ellipse(el):
    try:
        return (Edge.makeEllipse(
            el.dxf.major_axis.magnitude,
            el.minor_axis.magnitude,
            pnt=Vector(el.dxf.center.xyz),
            dir=el.dxf.extrusion.xyz,
            xdir=Vector(el.dxf.major_axis.xyz),
            angle1=el.dxf.start_param * RAD2DEG,
            angle2=el.dxf.end_param * RAD2DEG,
        ),)
    except Exception:
        return ()


DXF_CONVERTERS = {
    "LINE": _dxf_line,
    "CIRCLE": _dxf_circle,
    "ARC": _dxf_arc,
    "POLYLINE": _dxf_polyline,
    "LWPOLYLINE": _dxf_polyline,
    "SPLINE": _dxf_spline,
    "ELLIPSE": _dxf_ellipse,
}


def _dxf_convert(elements, tol: float) -> List[Shape]:
    edges = []
    for el in elements:
        conv = DXF_CONVERTERS.get(el.dxf.dxftype)
        if conv:
            edges.extend(conv(el))
    return _edges_to_wires(edges, tol) if edges else []


def _importDXF(filename: str, tol: float = 1e-6, exclude: List[str] = [], include: List[str] = []) -> List[Face]:
    """The faces bounded by the wires of the file's (selected) layers' curves."""
    if exclude and include:
        raise ValueError("you may specify either 'include' or 'exclude' but not both")
    # (Imported here, on the first file read: importing the shim, whose sketches import this
    # module, would otherwise wait a tenth of a second for ezdxf.)
    import ezdxf

    dxf = ezdxf.readfile(filename)
    faces = []
    layers = dxf.modelspace().groupby(dxfattrib="layer")
    names = set(name.lower() for name in layers.keys())
    if include:
        selected = names & set(name.lower() for name in include)
    elif exclude:
        selected = names - set(name.lower() for name in exclude)
    else:
        selected = names
    if not selected:
        raise ValueError("no DXF layers selected")
    for name, layer in layers.items():
        if name.lower() in selected:
            wires = _dxf_convert(layer, tol)
            for wire_set in sortWiresByBuildOrder(wires):
                faces.append(Face.makeFromWires(wire_set[0], wire_set[1:]))
    return faces
