"""Shapes drawn to SVG as CadQuery draws them: edges projected along a view direction,
visible ones solid and hidden ones dashed, in CadQuery's page template.

Hidden lines are found against the shape's triangulation (a depth test, as a z-buffer
does): a point of an edge is hidden where a triangle lies in front of it. Outlines of curved
faces (silhouettes) are not drawn.
"""

from __future__ import annotations

import math
from typing import Dict, List, Tuple

from ... import _kentron as native
from ...geom import Vector, _default_x
from ...shapes import Shape

DISCRETIZATION_TOLERANCE = 1e-3

SVG_TEMPLATE = """<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<svg
   xmlns:svg="http://www.w3.org/2000/svg"
   xmlns="http://www.w3.org/2000/svg"
   width="%(width)s"
   height="%(height)s"
>
    <g transform="scale(%(unitScale)s, -%(unitScale)s)   translate(%(xTranslate)s,%(yTranslate)s)" stroke-width="%(strokeWidth)s"  fill="none">
       <!-- hidden lines -->
       <g  stroke="rgb(%(hiddenColor)s)" fill="none" stroke-dasharray="%(strokeWidth)s,%(strokeWidth)s" >
%(hiddenContent)s
       </g>
       <!-- solid lines -->
       <g  stroke="rgb(%(strokeColor)s)" fill="none">
%(visibleContent)s
       </g>
    </g>
    %(axesIndicator)s
</svg>
"""

# The axes indicator (CadQuery's).
AXES_TEMPLATE = """<g transform="translate(20,%(textboxY)s)" stroke="rgb(0,0,255)">
        <line x1="30" y1="-30" x2="75" y2="-33" stroke-width="3" stroke="#000000" />
         <text x="80" y="-30" style="stroke:#000000">X </text>

        <line x1="30" y1="-30" x2="30" y2="-75" stroke-width="3" stroke="#000000" />
         <text x="25" y="-85" style="stroke:#000000">Y </text>

        <line x1="30" y1="-30" x2="58" y2="-15" stroke-width="3" stroke="#000000" />
         <text x="65" y="-5" style="stroke:#000000">Z </text>
        <!--
            <line x1="0" y1="0" x2="%(unitScale)s" y2="0" stroke-width="3" />
            <text x="0" y="20" style="stroke:#000000">1  %(uom)s </text>
        -->
    </g>"""

PATHTEMPLATE = '\t\t\t<path d="%s" />\n'


class UNITS:
    MM = "mm"
    IN = "in"


def guessUnitOfMeasure(shape) -> str:
    """Millimetres or inches, guessed from the shape's size (as CadQuery guesses)."""
    bb = shape.BoundingBox()
    dims = [bb.xlen, bb.ylen, bb.zlen]
    if max(dims) > 10:
        return UNITS.MM
    if min(dims) < 0.1:
        return UNITS.IN
    if sum(dims) < 10:
        return UNITS.IN
    return UNITS.MM


class _DepthTest:
    """The shape's triangles seen along the view, bucketed on a grid of the view plane."""

    def __init__(self, shape: Shape, x: Vector, y: Vector, z: Vector, size: float):
        self.cells: Dict[Tuple[int, int], list] = {}
        self.tolerance = 1e-3 * size
        self.slack = 2 * self.tolerance + 1e-9 * size
        self.step = size / 64 if size > 0 else 1.0
        if not shape.Faces():
            return
        points, triangles = native.triangulate(shape.wrapped, self.tolerance, 0.2)
        view = [(x.dot(p), y.dot(p), z.dot(p)) for p in map(Vector, points)]
        for a, b, c in triangles:
            t = (view[a], view[b], view[c])
            u0, u1 = min(p[0] for p in t), max(p[0] for p in t)
            v0, v1 = min(p[1] for p in t), max(p[1] for p in t)
            for i in range(self._cell(u0), self._cell(u1) + 1):
                for j in range(self._cell(v0), self._cell(v1) + 1):
                    self.cells.setdefault((i, j), []).append(t)

    def _cell(self, s: float) -> int:
        return int(math.floor(s / self.step))

    def hidden(self, u: float, v: float, depth: float) -> bool:
        """A triangle lies in front of the point (u, v) at `depth` (towards the viewer)."""
        for (ax, ay, az), (bx, by, bz), (cx, cy, cz) in self.cells.get((self._cell(u), self._cell(v)), ()):
            det = (by - ay) * (cx - ax) - (bx - ax) * (cy - ay)
            if abs(det) < 1e-300:
                continue
            l1 = ((by - cy) * (u - cx) + (cx - bx) * (v - cy)) / det
            l2 = ((cy - ay) * (u - cx) + (ax - cx) * (v - cy)) / det
            l3 = 1 - l1 - l2
            eps = 1e-9
            if l1 < -eps or l2 < -eps or l3 < -eps:
                continue
            if l1 * az + l2 * bz + l3 * cz > depth + self.slack:
                return True
        return False


def _projected_runs(shape: Shape, direction) -> Tuple[List[list], List[list]]:
    """The edges projected along `direction` (towards the viewer), split into visible and
    hidden runs of 2D points."""
    z = Vector(direction).normalized()
    x = _default_x(z)
    y = z.cross(x)
    size = shape.BoundingBox().DiagonalLength
    depth = _DepthTest(shape, x, y, z, size)
    visible, hidden = [], []
    for e in shape.Edges():
        u0, u1 = e.bounds()
        # (Lines sampled too: visibility can change along them.)
        n = 24 if e.geomType() == "LINE" else 48
        pts = [Vector(d[0]) for d in e._derivatives([u0 + (u1 - u0) * i / n for i in range(n + 1)])]
        run, run_hidden = [], None
        for p in pts:
            q = (x.dot(p), y.dot(p))
            h = depth.hidden(q[0], q[1], z.dot(p))
            if run_hidden is not None and h != run_hidden:
                (hidden if run_hidden else visible).append(run)
                run = [run[-1]]
            run.append(q)
            run_hidden = h
        if len(run) > 1:
            (hidden if run_hidden else visible).append(run)
    return visible, hidden


def _path(run) -> str:
    (x0, y0), rest = run[0], run[1:]
    return f"M{x0},{y0} " + "".join(f"L{x},{y} " for x, y in rest)


def getSVG(shape, opts=None) -> str:
    """The shape drawn as SVG text (CadQuery's options: width, height, marginLeft,
    marginTop, projectionDir, showAxes, strokeWidth, strokeColor, hiddenColor, showHidden;
    focus, for a perspective view, is not supported)."""
    d = {
        "width": 800,
        "height": 240,
        "marginLeft": 200,
        "marginTop": 20,
        "projectionDir": (-1.75, 1.1, 5),
        "showAxes": True,
        "strokeWidth": -1.0,  # (-1: from the scale)
        "strokeColor": (0, 0, 0),
        "hiddenColor": (160, 160, 160),
        "showHidden": True,
        "focus": None,
    }
    if opts:
        d.update(opts)
    uom = guessUnitOfMeasure(shape)
    width = float(d["width"]) if d["width"] is not None else None
    height = float(d["height"]) if d["height"] is not None else None
    marginLeft = float(d["marginLeft"])
    marginTop = float(d["marginTop"])
    projectionDir = tuple(d["projectionDir"])
    showAxes = bool(d["showAxes"])
    strokeWidth = float(d["strokeWidth"])
    strokeColor = tuple(d["strokeColor"])
    hiddenColor = tuple(d["hiddenColor"])
    showHidden = bool(d["showHidden"])
    if d.get("focus"):
        raise NotImplementedError("perspective SVG views")

    visible, hidden = _projected_runs(shape, projectionDir)
    points = [p for run in visible + hidden for p in run] or [(0.0, 0.0)]
    xmin, xmax = min(p[0] for p in points), max(p[0] for p in points)
    ymin, ymax = min(p[1] for p in points), max(p[1] for p in points)
    xlen, ylen = max(xmax - xmin, 1e-9), max(ymax - ymin, 1e-9)

    if width is None or height is None:
        if width is None:
            width = (height - 2.0 * marginTop) * (xlen / ylen) + 2.0 * marginLeft
        else:
            height = (width - 2.0 * marginLeft) * (ylen / xlen) + 2.0 * marginTop
        unitScale = (width - 2.0 * marginLeft) / xlen
    else:
        unitScale = min(width / xlen * 0.75, height / ylen * 0.75)
    xTranslate, yTranslate = (0 - xmin) + marginLeft / unitScale, (0 - ymax) - marginTop / unitScale
    if strokeWidth == -1.0:
        strokeWidth = 1.0 / unitScale

    hiddenContent = "".join(PATHTEMPLATE % _path(r) for r in hidden) if showHidden else ""
    visibleContent = "".join(PATHTEMPLATE % _path(r) for r in visible)
    if showAxes and projectionDir == (-1.75, 1.1, 5):
        axesIndicator = AXES_TEMPLATE % {"unitScale": str(unitScale), "textboxY": str(height - 30), "uom": str(uom)}
    else:
        axesIndicator = ""
    return SVG_TEMPLATE % {
        "unitScale": str(unitScale),
        "strokeWidth": str(strokeWidth),
        "strokeColor": ",".join(str(c) for c in strokeColor),
        "hiddenColor": ",".join(str(c) for c in hiddenColor),
        "hiddenContent": hiddenContent,
        "visibleContent": visibleContent,
        "xTranslate": str(xTranslate),
        "yTranslate": str(yTranslate),
        "width": str(width),
        "height": str(height),
        "textboxY": str(height - 30),
        "uom": str(uom),
        "axesIndicator": axesIndicator,
    }


def exportSVG(shape, fileName: str, opts=None) -> None:
    """Writes a workplane's (first) shape as SVG."""
    with open(fileName, "w") as f:
        f.write(getSVG(shape.val() if hasattr(shape, "val") else shape, opts))
