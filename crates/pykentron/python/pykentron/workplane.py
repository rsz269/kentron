"""The fluent `Workplane` (the part of CadQuery's that 2D-sketch-then-extrude code uses).

A workplane holds a plane, a stack of objects (shapes, points) and pending wires drawn on it;
operations return new workplanes whose parent is the one they were called on.
"""

from __future__ import annotations

import math
from typing import List, Optional, Sequence

from .geom import Location, Plane, Vector
from .shapes import Compound, Edge, Face, Shape, Solid, Wire
from . import selectors


class Workplane:
    """CadQuery's `Workplane`: a plane with a stack of objects."""

    def __init__(self, inPlane="XY", origin=(0, 0, 0), obj=None):
        if isinstance(inPlane, Plane):
            plane = inPlane
        elif isinstance(inPlane, str):
            plane = Plane.named(inPlane, origin)
        elif isinstance(inPlane, (Vector, Location, Shape)):
            # An object for the stack, on the XY plane (as CadQuery's `Workplane(obj)`).
            obj = inPlane
            plane = Plane.named("XY", origin)
        else:
            raise ValueError(f"a workplane needs a plane, got {inPlane!r}")
        if isinstance(inPlane, Plane) and tuple(origin) != (0, 0, 0):
            plane = Plane(origin, plane.xDir, plane.zDir)
        self.plane = plane
        self.objects: list = [obj] if obj is not None else []
        self.parent: Optional[Workplane] = None
        # Wires being drawn: edges of the current one, and finished (closed) wires.
        self._edges: List[Edge] = []
        self._wires: List[Wire] = []

    def _child(self, objects, edges=None, wires=None) -> "Workplane":
        w = Workplane.__new__(Workplane)
        w.plane = self.plane
        w.objects = list(objects)
        w.parent = self
        w._edges = list(self._edges if edges is None else edges)
        w._wires = list(self._wires if wires is None else wires)
        return w

    # --- the stack ------------------------------------------------------------------------

    def val(self):
        return self.objects[0] if self.objects else None

    def vals(self) -> list:
        return list(self.objects)

    def first(self) -> "Workplane":
        return self._child(self.objects[:1])

    def last(self) -> "Workplane":
        return self._child(self.objects[-1:])

    def size(self) -> int:
        return len(self.objects)

    def all(self) -> list:
        return [self._child([o]) for o in self.objects]

    def add(self, obj) -> "Workplane":
        if isinstance(obj, Workplane):
            new = obj.objects
        elif isinstance(obj, (list, tuple)):
            new = list(obj)
        else:
            new = [obj]
        return self._child(self.objects + new)

    def end(self, n: int = 1) -> "Workplane":
        w = self
        for _ in range(n):
            w = w.parent
        return w

    def __iter__(self):
        return iter(self.objects)

    def _solid(self) -> Optional[Shape]:
        """The solid this chain works on: the nearest on the stack or a parent's."""
        w = self
        while w is not None:
            for o in w.objects:
                if isinstance(o, (Solid, Compound)):
                    return o
            w = w.parent
        return None

    def findSolid(self) -> Shape:
        s = self._solid()
        if s is None:
            raise ValueError("no solid on the stack")
        return s

    # --- drawing --------------------------------------------------------------------------

    def _world(self, p) -> Vector:
        return self.plane.toWorldCoords(Vector(p))

    def _centres(self) -> list:
        pts = [o for o in self.objects if isinstance(o, Vector)]
        return pts or [self.plane.origin]

    def moveTo(self, x: float = 0, y: float = 0) -> "Workplane":
        return self._child([self._world((x, y))], edges=[])

    def polyline(self, listOfXYTuple: Sequence, forConstruction: bool = False, includeCurrent: bool = False) -> "Workplane":
        pts = [self._world(p) for p in listOfXYTuple]
        if includeCurrent and self.objects and isinstance(self.objects[-1], Vector):
            pts.insert(0, self.objects[-1])
        edges = [Edge.makeLine(a, b) for a, b in zip(pts, pts[1:]) if (a - b).Length > 1e-9]
        return self._child([pts[-1]], edges=self._edges + edges)

    def lineTo(self, x: float, y: float, forConstruction: bool = False) -> "Workplane":
        start = self.objects[-1] if self.objects and isinstance(self.objects[-1], Vector) else self.plane.origin
        end = self._world((x, y))
        return self._child([end], edges=self._edges + [Edge.makeLine(start, end)])

    def line(self, xDist: float, yDist: float, forConstruction: bool = False) -> "Workplane":
        start = self.objects[-1] if self.objects and isinstance(self.objects[-1], Vector) else self.plane.origin
        local = self.plane.toLocalCoords(start)
        return self.lineTo(local.x + xDist, local.y + yDist)

    def hLine(self, distance: float, forConstruction: bool = False) -> "Workplane":
        return self.line(distance, 0)

    def vLine(self, distance: float, forConstruction: bool = False) -> "Workplane":
        return self.line(0, distance)

    def threePointArc(self, point1, point2, forConstruction: bool = False) -> "Workplane":
        start = self.objects[-1] if self.objects and isinstance(self.objects[-1], Vector) else self.plane.origin
        mid, end = self._world(point1), self._world(point2)
        return self._child([end], edges=self._edges + [Edge.makeThreePointArc(start, mid, end)])

    def close(self) -> "Workplane":
        edges = list(self._edges)
        if not edges:
            raise ValueError("nothing to close")
        start = edges[0].startPoint()
        end = edges[-1].endPoint()
        if (start - end).Length > 1e-9:
            edges.append(Edge.makeLine(end, start))
        wire = Wire.assembleEdges(edges)
        return self._child([wire], edges=[], wires=self._wires + [wire])

    def wire(self, forConstruction: bool = False) -> "Workplane":
        wire = Wire.assembleEdges(self._edges)
        return self._child([wire], edges=[], wires=self._wires + [wire])

    def rect(self, xLen: float, yLen: float, centered=True, forConstruction: bool = False) -> "Workplane":
        cx, cy = (centered, centered) if isinstance(centered, bool) else centered
        wires = []
        for c in self._centres():
            local = self.plane.toLocalCoords(c)
            x0 = local.x - (xLen / 2 if cx else 0)
            y0 = local.y - (yLen / 2 if cy else 0)
            pts = [self._world(p) for p in ((x0, y0), (x0 + xLen, y0), (x0 + xLen, y0 + yLen), (x0, y0 + yLen))]
            wires.append(Wire.makePolygon(pts, close=True))
        return self._child(wires, wires=self._wires + wires)

    def circle(self, radius: float, forConstruction: bool = False) -> "Workplane":
        wires = [Wire.makeCircle(radius, c, self.plane.zDir) for c in self._centres()]
        return self._child(wires, wires=self._wires + wires)

    def polygon(self, nSides: int, diameter: float, forConstruction: bool = False, circumscribed: bool = False) -> "Workplane":
        # The circle through the corners has the given diameter (through the sides' midpoints
        # if circumscribed); the first corner on the plane's X axis.
        r = diameter / 2
        if circumscribed:
            r /= math.cos(math.pi / nSides)
        wires = []
        for c in self._centres():
            local = self.plane.toLocalCoords(c)
            pts = [
                self._world((local.x + r * math.cos(2 * math.pi * i / nSides), local.y + r * math.sin(2 * math.pi * i / nSides)))
                for i in range(nSides)
            ]
            wires.append(Wire.makePolygon(pts, close=True))
        return self._child(wires, wires=self._wires + wires)

    def pushPoints(self, pntList) -> "Workplane":
        return self._child([self._world(p) if not isinstance(p, Location) else p._apply_point((0, 0, 0)) for p in pntList])

    def center(self, x: float, y: float) -> "Workplane":
        w = self._child(self.objects)
        w.plane = Plane(self.plane.toWorldCoords((x, y)), self.plane.xDir, self.plane.zDir)
        return w

    def transformed(self, rotate=(0, 0, 0), offset=(0, 0, 0)) -> "Workplane":
        w = self._child(self.objects)
        p = self.plane.rotated(rotate)
        w.plane = Plane(self.plane.toWorldCoords(offset), p.xDir, p.zDir)
        return w

    def workplane(self, offset: float = 0.0, invert: bool = False, centerOption: str = "ProjectedOrigin", origin=None) -> "Workplane":
        w = self._child([])
        z = -self.plane.zDir if invert else self.plane.zDir
        base = Vector(origin) if origin is not None else self.plane.origin
        w.plane = Plane(base + self.plane.zDir * offset, self.plane.xDir, z)
        return w

    # --- solids ---------------------------------------------------------------------------

    def _faces(self) -> List[Face]:
        """Faces of the pending wires: each outer wire with the wires nested inside it."""
        wires = list(self._wires)
        if not wires:
            raise ValueError("no pending wires")
        faces_of = [Face.makeFromWires(w) for w in wires]
        areas = [f.Area() for f in faces_of]
        order = sorted(range(len(wires)), key=lambda i: -areas[i])
        result, used = [], set()
        for i in order:
            if i in used:
                continue
            used.add(i)
            box = faces_of[i].BoundingBox()
            inner = []
            for j in order:
                if j in used:
                    continue
                if box.isInside(faces_of[j].BoundingBox()) and faces_of[i].intersect(faces_of[j]).Area() > 0.999 * areas[j]:
                    inner.append(wires[j])
                    used.add(j)
            result.append(Face.makeFromWires(wires[i], inner))
        return result

    def _combine(self, new: Shape, combine, clean: bool = True) -> "Workplane":
        base = self._solid()
        if combine == "cut":
            if base is None:
                raise ValueError("nothing to cut from")
            new = base.cut(new)
        elif combine and base is not None:
            new = base.fuse(new)
        if clean:
            new = new.clean()
        return self._child([_single(new)], edges=[], wires=[])

    def extrude(self, until, combine=True, clean: bool = True, both: bool = False, taper=None) -> "Workplane":
        if taper:
            raise NotImplementedError("tapered extrusions")
        d = float(until)
        normal = self.plane.zDir
        solids = []
        for face in self._faces():
            if both:
                face = face.translate(normal * -d)
                solids.append(Solid.extrudeLinear(face, normal * (2 * d)))
            else:
                solids.append(Solid.extrudeLinear(face, normal * d))
        new = solids[0] if len(solids) == 1 else Compound.makeCompound(solids)
        return self._combine(new, combine, clean)

    def cutBlind(self, until, clean: bool = True, both: bool = False, taper=None) -> "Workplane":
        return self.extrude(until, combine="cut", clean=clean, both=both, taper=taper)

    def revolve(self, angleDegrees: float = 360.0, axisStart=None, axisEnd=None, combine=True, clean: bool = True) -> "Workplane":
        start = self._world(axisStart if axisStart is not None else (0, 0))
        end = self._world(axisEnd if axisEnd is not None else (0, 1))
        solids = [Solid.revolve(face, angleDegrees, start, end) for face in self._faces()]
        new = solids[0] if len(solids) == 1 else Compound.makeCompound(solids)
        return self._combine(new, combine, clean)

    def box(self, length: float, width: float, height: float, centered=True, combine=True, clean: bool = True) -> "Workplane":
        cx, cy, cz = (centered,) * 3 if isinstance(centered, bool) else centered
        boxes = []
        for c in self._centres():
            local = self.plane.toLocalCoords(c)
            corner = (
                local.x - (length / 2 if cx else 0),
                local.y - (width / 2 if cy else 0),
                local.z - (height / 2 if cz else 0),
            )
            box = Solid.makeBox(length, width, height)
            boxes.append(box.moved(Location(self.plane)).translate(self.plane.xDir * corner[0] + self.plane.yDir * corner[1] + self.plane.zDir * corner[2]))
        new = boxes[0] if len(boxes) == 1 else Compound.makeCompound(boxes)
        return self._combine(new, combine, clean)

    def cylinder(self, height: float, radius: float, direct=None, angle: float = 360, centered=True, combine=True, clean: bool = True) -> "Workplane":
        cx, cy, cz = (centered,) * 3 if isinstance(centered, bool) else centered
        shapes = []
        for c in self._centres():
            local = self.plane.toLocalCoords(c)
            base = self._world((local.x + (0 if cx else radius), local.y + (0 if cy else radius), local.z - (height / 2 if cz else 0)))
            shapes.append(Solid.makeCylinder(radius, height, base, direct or self.plane.zDir, angle))
        new = shapes[0] if len(shapes) == 1 else Compound.makeCompound(shapes)
        return self._combine(new, combine, clean)

    def union(self, other, clean: bool = True) -> "Workplane":
        other = other.val() if isinstance(other, Workplane) else other
        base = self._solid()
        result = base.fuse(other) if base is not None else other
        return self._child([result.clean() if clean else result])

    def cut(self, other, clean: bool = True) -> "Workplane":
        other = other.val() if isinstance(other, Workplane) else other
        result = self.findSolid().cut(other)
        return self._child([result.clean() if clean else result])

    def intersect(self, other, clean: bool = True) -> "Workplane":
        other = other.val() if isinstance(other, Workplane) else other
        result = self.findSolid().intersect(other)
        return self._child([result.clean() if clean else result])

    def translate(self, vec) -> "Workplane":
        return self._child([o.translate(vec) if isinstance(o, Shape) else o + Vector(vec) for o in self.objects])

    def rotate(self, axisStartPoint, axisEndPoint, angleDegrees: float) -> "Workplane":
        return self._child([o.rotate(axisStartPoint, axisEndPoint, angleDegrees) for o in self.objects])

    def mirror(self, mirrorPlane="XY", basePointVector=(0, 0, 0), union: bool = False) -> "Workplane":
        mirrored = [o.mirror(mirrorPlane, basePointVector) for o in self.objects]
        if union:
            return self._child([o.fuse(m) for o, m in zip(self.objects, mirrored)])
        return self._child(mirrored)

    def fillet(self, radius: float) -> "Workplane":
        solid = self.findSolid()
        edges = [o for o in self.objects if isinstance(o, Edge)]
        if not edges:
            raise ValueError("fillet needs edges selected")
        return self._child([Shape.cast(_native().fillet(solid.wrapped, [e.wrapped for e in edges], radius))])

    # --- selection ------------------------------------------------------------------------

    def _select(self, kind: str, selector) -> "Workplane":
        found, seen = [], set()
        for o in self.objects:
            if isinstance(o, Shape):
                for s in getattr(o, kind)():
                    if s.hashCode() not in seen or not any(s.isSame(f) for f in found):
                        seen.add(s.hashCode())
                        found.append(s)
        if selector is not None:
            found = selectors.select(found, selector)
        return self._child(found)

    def solids(self, selector=None) -> "Workplane":
        return self._select("Solids", selector)

    def faces(self, selector=None) -> "Workplane":
        return self._select("Faces", selector)

    def edges(self, selector=None) -> "Workplane":
        return self._select("Edges", selector)

    def wires(self, selector=None) -> "Workplane":
        return self._select("Wires", selector)

    def vertices(self, selector=None) -> "Workplane":
        return self._select("Vertices", selector)

    def shells(self, selector=None) -> "Workplane":
        return self._select("Shells", selector)

    def toPending(self) -> "Workplane":
        wires = [o for o in self.objects if isinstance(o, Wire)]
        return self._child(self.objects, wires=self._wires + wires)

    def __repr__(self) -> str:
        return f"<Workplane objects={self.objects!r}>"


def _single(shape: Shape) -> Shape:
    """A compound of one solid as that solid (as CadQuery leaves a workplane's result)."""
    if isinstance(shape, Compound):
        solids = shape.Solids()
        if len(solids) == 1:
            return solids[0]
    return shape


def _native():
    from . import _kentron

    return _kentron
