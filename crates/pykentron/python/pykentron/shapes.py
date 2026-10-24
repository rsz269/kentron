"""Shapes (CadQuery's `Shape` and its kinds), over kentron's shape handles.

Shapes are values: operations return new shapes, and the shapes they are made from stay as
they were (they share what they did not change, in kentron's model). `move` and `locate`
are the exceptions, as in CadQuery: they re-place the shape in place.

Curves are read as CadQuery reads them through OCCT's adaptors: an edge by its own curve and
parameter range, whatever the edge's orientation; a wire by its edges in order, each run
along the wire, their ranges laid end to end from 0.
"""

from __future__ import annotations

import math
import os
import sys
import tempfile
import warnings
from typing import Any, Callable, Dict, Iterable, List, Optional, Sequence, Tuple

from . import _kentron as native
from .geom import BoundBox, Location, Matrix, Plane, Vector, _compose, _default_x, _invert, _new, _rotation_rows, _rows_identity, _vector


def _rows_translation(v) -> list:
    v = Vector(v)
    rows = _rows_identity()
    rows[0][3], rows[1][3], rows[2][3] = v.x, v.y, v.z
    return rows


def _rows_rotation(start, end, degrees) -> list:
    start, end = Vector(start), Vector(end)
    rows = _rotation_rows(end - start, math.radians(degrees))
    # About the axis through `start`: p -> R (p - s) + s.
    for i in range(3):
        rows[i][3] = start[i] - sum(rows[i][k] * start[k] for k in range(3))
    return rows


def _mirror_rows(normal, base) -> list:
    n = Vector(normal).normalized()
    b = Vector(base)
    rows = []
    for i in range(3):
        row = [(1.0 if i == j else 0.0) - 2 * n[i] * n[j] for j in range(3)]
        row.append(2 * n[i] * b.dot(n))
        rows.append(row)
    return rows


_PLANE_NORMALS = {"XY": (0, 0, 1), "YX": (0, 0, 1), "XZ": (0, 1, 0), "ZX": (0, 1, 0), "YZ": (1, 0, 0), "ZY": (1, 0, 0)}


def _selector(selector):
    """A selector object for `selector` (a string in CadQuery's syntax, or one already)."""
    if not isinstance(selector, str):
        return selector
    module = sys.modules.get("cadquery.selectors")
    if module is None and "cadquery" in sys.modules:
        try:
            from cadquery import selectors as module  # (the shim's: CadQuery's own)
        except ImportError:
            module = None
    if module is None:
        from . import selectors as module
    return module.StringSyntaxSelector(selector)


def _tuple(v) -> Tuple[float, float, float]:
    v = _vector(v)
    return (v.x, v.y, v.z)


class Shape:
    """A shape: a vertex, edge, wire, face, shell, solid or compound."""

    def __init__(self, obj):
        self.wrapped = obj
        self.label = ""
        self.forConstruction = False

    @staticmethod
    def cast(obj, forConstruction: bool = False) -> "Shape":
        return _cast_as(_KINDS.get(obj.kind(), Shape), obj, forConstruction)

    def _new(self, handle) -> "Shape":
        return Shape.cast(handle)

    def _copy_with(self, handle) -> "Shape":
        """The shape as `handle`, keeping this one's construction flag and label (as
        CadQuery's copies do, for placements)."""
        s = Shape.cast(handle, self.forConstruction)
        s.label = self.label
        return s

    # --- identity -------------------------------------------------------------------------

    def ShapeType(self) -> str:
        return self.wrapped.kind()

    def geomType(self) -> str:
        return self.wrapped.kind().upper()

    def isNull(self) -> bool:
        return self.wrapped is None or self.wrapped.IsNull()

    def isSame(self, other: "Shape") -> bool:
        return self.wrapped.is_same(other.wrapped)

    def isEqual(self, other: "Shape") -> bool:
        return self.wrapped.is_equal(other.wrapped)

    def hashCode(self) -> int:
        return self.wrapped.node_hash()

    def __hash__(self) -> int:
        return self.hashCode()

    def __eq__(self, other) -> bool:
        return isinstance(other, Shape) and self.isSame(other)

    def __repr__(self) -> str:
        return f"<{self.wrapped.kind()}>"

    def __getstate__(self):
        # A BRep image, as CadQuery's (BinTools); the handle itself is the model's.
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "shape.brep")
            native.write_brep(self.wrapped, path)
            with open(path, "rb") as f:
                return (f.read(), self.forConstruction)

    def __setstate__(self, state):
        data, forConstruction = state
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "shape.brep")
            with open(path, "wb") as f:
                f.write(data)
            self.wrapped = native.read_brep(path)
        self.label = ""
        self.forConstruction = forConstruction

    # --- placement ------------------------------------------------------------------------

    # (OCCT's translate, rotate and transformShape move a copy's geometry, which its location
    # (`location`, `located`, `locate`) then leaves alone; pykentron places the shape
    # instead, and keeps that placement as `_baked`: the location is what is placed beyond
    # it.)

    def _transformed_copy(self, rows) -> "Shape":
        s = self._new(self.wrapped.moved(rows))
        s._baked = [list(r) for r in s.wrapped.placement()]
        return s

    def _baked_rows(self):
        return self.__dict__.get("_baked") or _rows_identity()

    def _with_baked(self, s: "Shape") -> "Shape":
        if "_baked" in self.__dict__:
            s._baked = self._baked
        return s

    def translate(self, vector) -> "Shape":
        return self._transformed_copy(_rows_translation(vector))

    def rotate(self, startVector, endVector, angleDegrees: float) -> "Shape":
        return self._transformed_copy(_rows_rotation(startVector, endVector, angleDegrees))

    def mirror(self, mirrorPlane="XY", basePointVector=(0, 0, 0)) -> "Shape":
        normal = _PLANE_NORMALS[mirrorPlane] if isinstance(mirrorPlane, str) else mirrorPlane
        # A reflection: the geometry is moved (placements are rotations).
        return self._new(native.transformed(self.wrapped, _mirror_rows(normal, basePointVector)))

    def scale(self, factor: float) -> "Shape":
        rows = [[factor if i == j else 0.0 for j in range(3)] + [0.0] for i in range(3)]
        return self._new(native.transformed(self.wrapped, rows))

    def transformShape(self, m) -> "Shape":
        rows = _matrix_rows(m)
        if _is_rigid(rows):
            return self._transformed_copy(rows)
        return self._new(native.transformed(self.wrapped, rows))

    def transformGeometry(self, m) -> "Shape":
        return self._new(native.transformed(self.wrapped, _matrix_rows(m)))

    def location(self) -> Location:
        placement = [list(r) for r in self.wrapped.placement()]
        if "_baked" not in self.__dict__:
            return Location(placement)
        return Location(_compose(placement, _invert(self._baked)))

    def locate(self, loc: Location) -> "Shape":
        self.wrapped = self.wrapped.located(_compose(loc.rows, self._baked_rows()))
        self.__dict__.pop("_boxes", None)
        return self

    def located(self, loc: Location) -> "Shape":
        return self._with_baked(self._copy_with(self.wrapped.located(_compose(loc.rows, self._baked_rows()))))

    def move(self, *args, **kwargs) -> "Shape":
        loc = _location_of(*args, **kwargs)
        self.wrapped = self.wrapped.moved(loc.rows)
        self.__dict__.pop("_boxes", None)
        return self

    def moved(self, *args, **kwargs) -> "Shape":
        # Another shape: copies at its locations (CadQuery's `moved(s.toLocs())`).
        if len(args) == 1 and not kwargs and isinstance(args[0], Shape):
            return self.moved(args[0].toLocs())
        if kwargs:
            return self._with_baked(self._copy_with(self.wrapped.moved(_location_of(*args, **kwargs).rows)))
        # Several locations: a compound of copies.
        if len(args) == 1 and isinstance(args[0], (list, tuple)) and args[0] and isinstance(args[0][0], (Location, list, tuple, Vector)) and not _is_vector_like(args[0]):
            args = tuple(args[0])
        if len(args) > 1 and all(isinstance(a, (Location, Vector)) or _is_vector_like(a) for a in args) and not _is_scalars(args):
            return Compound.makeCompound([self.moved(a) for a in args])
        loc = _location_of(*args)
        return self._with_baked(self._copy_with(self.wrapped.moved(loc.rows)))

    def copy(self, mesh: bool = False) -> "Shape":
        return self._with_baked(self._copy_with(self.wrapped))

    def reverse(self) -> "Shape":
        """The shape with its orientation reversed."""
        return self._new(self.wrapped.reversed())

    # --- booleans -------------------------------------------------------------------------

    def fuse(self, *toFuse: "Shape", glue: bool = False, tol: Optional[float] = None) -> "Shape":
        # (Pairwise, then pairs of pairs: each step a union of shapes of like size, not each
        # shape added in turn to a whole that grows with every one.)
        shapes = [self.wrapped] + [other.wrapped for other in toFuse]
        while len(shapes) > 1:
            shapes = [native.union(shapes[i], shapes[i + 1]) if i + 1 < len(shapes) else shapes[i] for i in range(0, len(shapes), 2)]
        return Shape.cast(shapes[0])

    def cut(self, *toCut: "Shape", tol: Optional[float] = None) -> "Shape":
        result = self.wrapped
        for other in toCut:
            result = _cut_solid_by_solid(result, other.wrapped)
        return Shape.cast(result)

    def intersect(self, *toIntersect: "Shape", tol: Optional[float] = None) -> "Shape":
        # (With several, as OCCT's common: what lies within any of them, their union.)
        if not toIntersect:
            return self
        tools = toIntersect[0] if len(toIntersect) == 1 else Shape.fuse(*toIntersect)
        return Shape.cast(native.intersection(self.wrapped, tools.wrapped))

    def split(self, *splitters: "Shape") -> "Shape":
        """This shape cut into pieces by `splitters` (BRepAlgoAPI_Splitter): a compound of
        them; the piece itself where it is not cut (as OCCT gives it)."""
        pieces = _split_pieces(self, splitters)
        return pieces[0] if len(pieces) == 1 else Compound.makeCompound(pieces)

    def __add__(self, other):
        from .occ_impl.shapes import fuse

        return fuse(self, other)

    def __sub__(self, other):
        from .occ_impl.shapes import cut

        return cut(self, other)

    def __mul__(self, other):
        from .occ_impl.shapes import intersect

        return intersect(self, other)

    def __truediv__(self, other):
        from .occ_impl.shapes import split

        return split(self, other)

    def _members(self) -> set:
        return set(self) if isinstance(self, Compound) else {self}

    def __and__(self, other: "Shape") -> "Compound":
        """The shapes in both selections."""
        return Compound.makeCompound(self._members() & other._members())

    def __or__(self, other: "Shape") -> "Compound":
        """The shapes in either selection."""
        return Compound.makeCompound(self._members() | other._members())

    def __mod__(self, other: "Shape") -> "Compound":
        """The shapes of this selection not in `other`."""
        return Compound.makeCompound(self._members() - other._members())

    def clean(self) -> "Shape":
        """Faces on the same surfaces, and edges on the same curves, merged."""
        return self._new(native.unify_same_domain(self.wrapped))

    def fix(self) -> "Shape":
        return self

    # --- queries --------------------------------------------------------------------------

    def isValid(self) -> bool:
        return native.is_valid(self.wrapped)

    def Volume(self, tolerance: Optional[float] = None) -> float:
        return native.volume(self.wrapped)[0]

    def Area(self) -> float:
        return native.area(self.wrapped)[0]

    def Length(self) -> float:
        return native.length(self.wrapped)[0]

    def Center(self) -> Vector:
        return Shape.centerOfMass(self)

    @staticmethod
    def _dimension(obj: "Shape") -> int:
        """The dimension CadQuery measures a shape in: a compound's first leaf's (a solid's
        if empty)."""
        kind = obj.ShapeType()
        while kind == "Compound":
            children = list(obj)
            if not children:
                return 3
            obj = children[0]
            kind = obj.ShapeType()
        return {"Vertex": 0, "Edge": 1, "Wire": 1, "Face": 2, "Shell": 2}.get(kind, 3)

    @staticmethod
    def computeMass(obj: "Shape", tol: Optional[float] = None) -> float:
        """The length, area or volume of `obj`, by its dimension."""
        d = Shape._dimension(obj)
        if d == 0:
            return 0.0
        return (native.length, native.area, native.volume)[d - 1](obj.wrapped)[0]

    @staticmethod
    def centerOfMass(obj: "Shape") -> Vector:
        kind = obj.wrapped.kind()
        if kind == "Vertex":
            return _new(*native.vertex_point(obj.wrapped))
        d = Shape._dimension(obj)
        return _new(*(native.length, native.area, native.volume)[max(d, 1) - 1](obj.wrapped)[1])

    @staticmethod
    def matrixOfInertia(obj: "Shape") -> List[List[float]]:
        """The matrix of inertia about the centre of mass (unit density)."""
        d = Shape._dimension(obj)
        if d == 0:
            raise NotImplementedError("the inertia of a vertex")
        return [list(r) for r in native.inertia(obj.wrapped, d)]

    def CenterOfBoundBox(self, tolerance: Optional[float] = None) -> Vector:
        return self.BoundingBox().center

    def BoundingBox(self, tolerance: Optional[float] = None, optimal: bool = True) -> BoundBox:
        # Once tessellated, a shape is bounded by its meshes, widened by their deviation from
        # the surfaces, as CadQuery's (OCCT's) boxes are.
        # Shapes do not change, so the box is computed once per shape (and mesh state).
        meshed = getattr(self, "_meshed", False)
        cache = self.__dict__.setdefault("_boxes", {})
        if meshed not in cache:
            cache[meshed] = native.bounding_box(self.wrapped, meshed)
        (x0, y0, z0), (x1, y1, z1) = cache[meshed]
        return BoundBox(x0, y0, z0, x1, y1, z1)

    def distance(self, other: "Shape") -> float:
        """The smallest distance between the shapes (zero where they touch, cross, or one
        lies inside a solid of the other)."""
        return native.distance(self.wrapped, other.wrapped)[0]

    def distances(self, *others: "Shape"):
        for o in others:
            yield self.distance(o)

    def facesIntersectedByLine(self, point, axis, tol: float = 1e-4, direction: Optional[str] = None) -> List["Face"]:
        """The faces the line through `point` along `axis` meets (only those ahead of the
        point, "AlongAxis", or behind it, "Opposite", if given), nearest first (a face met
        twice counted at its nearest meeting)."""
        if direction not in (None, "AlongAxis", "Opposite"):
            raise ValueError("Invalid direction specification.\nValid specification are 'AlongAxis' and 'Opposite'.")
        p, a = Vector(point), Vector(axis).normalized()
        box = self.BoundingBox()
        reach = box.DiagonalLength + (p - box.center).Length + 1.0
        line = native.make_line((p - a * reach).toTuple(), (p + a * reach).toTuple())
        found = []
        for f in self.Faces():
            d, _, q = native.distance(line, f.wrapped)
            if d > tol:
                continue
            offset = Vector(q) - p
            along = offset.dot(a)
            if direction == "AlongAxis" and not (along > 0 and offset.Length > tol):
                continue
            if direction == "Opposite" and not (along < 0 and offset.Length > tol):
                continue
            found.append((offset.Length, f))
        found.sort(key=lambda t: t[0])
        return [f for _, f in found]

    def tessellate(self, tolerance: float, angularTolerance: float = 0.1):
        """Points (Vectors) and triangles (index triples), counter-clockwise seen from outside:
        each face's points of its own (as CadQuery's)."""
        vertices, triangles = [], []
        for face in self.Faces():
            points, tris = native.triangulate(face.wrapped, tolerance, angularTolerance)
            if not tris:
                continue
            offset = len(vertices)
            vertices.extend(Vector(*p) for p in points)
            triangles.extend(tuple(i + offset for i in t) for t in tris)
        self._meshed = True
        return vertices, triangles

    def mesh(self, tolerance: float, angularTolerance: float = 0.1) -> None:
        """Meshes the shape's faces (where they are not yet)."""
        native.triangulate(self.wrapped, tolerance, angularTolerance)
        self._meshed = True

    # --- sub-shapes -----------------------------------------------------------------------

    def _sub(self, kind: str) -> list:
        cls = _KINDS[kind]
        return [_cast_as(cls, h) for h in native.subshapes(self.wrapped, kind, kind == "Edge")]

    def Vertices(self) -> List["Vertex"]:
        return self._sub("Vertex")

    def Edges(self) -> List["Edge"]:
        """The edges (but those collapsed to a point)."""
        return self._sub("Edge")

    def Wires(self) -> List["Wire"]:
        return self._sub("Wire")

    def Faces(self) -> List["Face"]:
        return self._sub("Face")

    def Shells(self) -> List["Shell"]:
        return self._sub("Shell")

    def Solids(self) -> List["Solid"]:
        return self._sub("Solid")

    def CompSolids(self) -> list:
        return self._sub("CompSolid")

    def Compounds(self) -> list:
        return self._sub("Compound")

    def __iter__(self):
        return (Shape.cast(h) for h in native.children(self.wrapped))

    def _entities(self, topo_type: str) -> list:
        """The native shapes of the sub-shapes of a type (CadQuery's, by its type name)."""
        return [h for h in native.subshapes(self.wrapped, topo_type)]

    def _entitiesFrom(self, child_type: str, parent_type: str) -> Dict["Shape", List["Shape"]]:
        """Each sub-shape of `child_type` with the sub-shapes of `parent_type` holding it."""
        out: Dict[Shape, List[Shape]] = {}
        for parent in self._sub(parent_type):
            for child in parent._sub(child_type):
                out.setdefault(child, []).append(parent)
        return out

    def ancestors(self, ctx: "Shape", kind: str) -> "Compound":
        """The shapes of `kind` within `ctx` holding this one."""
        mine = self.ShapeType()
        return Compound.makeCompound([p for p in ctx._sub(kind) if any(c.isSame(self) for c in p._sub(mine))])

    def siblings(self, ctx: "Shape", kind: str, level=1) -> "Compound":
        """The shapes within `ctx` sharing sub-shapes of `kind` with this one (and, with
        `level` above 1, theirs in turn)."""
        return Compound.makeCompound(_siblings([self], ctx, kind, level))

    def _filter(self, selector, objs: Iterable["Shape"]) -> "Shape":
        selected = _selector(selector).filter(list(objs)) if selector else list(objs)
        return selected[0] if len(selected) == 1 else Compound.makeCompound(selected)

    def _filter_single(self, selector, objs: Iterable["Shape"]) -> "Shape":
        selected = _selector(selector).filter(list(objs)) if selector else list(objs)
        if not selected:
            raise ValueError("No Shape available for selection")
        return selected[0]

    def vertices(self, selector=None) -> "Shape":
        return self._filter(selector, self._sub("Vertex"))

    def edges(self, selector=None) -> "Shape":
        return self._filter(selector, self._sub("Edge"))

    def wires(self, selector=None) -> "Shape":
        return self._filter(selector, self._sub("Wire"))

    def faces(self, selector=None) -> "Shape":
        return self._filter(selector, self._sub("Face"))

    def shells(self, selector=None) -> "Shape":
        return self._filter(selector, self._sub("Shell"))

    def solids(self, selector=None) -> "Shape":
        return self._filter(selector, self._sub("Solid"))

    def vertex(self, selector=None) -> "Vertex":
        return self._filter_single(selector, self._sub("Vertex"))

    def edge(self, selector=None) -> "Edge":
        return self._filter_single(selector, self._sub("Edge"))

    def wire(self, selector=None) -> "Wire":
        return self._filter_single(selector, self._sub("Wire"))

    def face(self, selector=None) -> "Face":
        return self._filter_single(selector, self._sub("Face"))

    def shell(self, selector=None) -> "Shell":
        return self._filter_single(selector, self._sub("Shell"))

    def solid(self, selector=None) -> "Solid":
        return self._filter_single(selector, self._sub("Solid"))

    def filter(self, f: Callable[["Shape"], bool]) -> "Compound":
        return Compound.makeCompound([s for s in self if f(s)])

    def sort(self, key: Callable[["Shape"], Any]) -> "Compound":
        return Compound.makeCompound(sorted(self, key=key))

    def __getitem__(self, item):
        children = list(self)
        if isinstance(item, slice):
            return Compound.makeCompound(children[item])
        return children[item]

    def size(self) -> int:
        return len(list(self))

    def toLocs(self) -> List[Location]:
        """The shape as locations: a face's or edge's at its middle, else its centre's."""
        t = self.ShapeType()
        if t == "Compound":
            return [loc for el in self for loc in el.toLocs()]
        if t == "Face":
            u0, u1, v0, v1 = self.uvBounds()
            return [self.locationAt((u0 + u1) / 2, (v0 + v1) / 2)]
        if t == "Edge":
            u0, u1 = self.bounds()
            return [self.locationAt((u0 + u1) / 2, mode="parameter")]
        return [Location(self.Center())]

    @staticmethod
    def CombinedCenter(objects) -> "Vector":
        """The centre of mass of several shapes together."""
        objects = list(objects)
        masses = [Shape.computeMass(o) for o in objects]
        total = sum(masses)
        acc = Vector(0, 0, 0)
        for o, m in zip(objects, masses):
            acc = acc + Shape.centerOfMass(o) * m
        return acc * (1.0 / total)

    @staticmethod
    def CombinedCenterOfBoundBox(objects) -> "Vector":
        """The mean of several shapes' bounding-box centres."""
        objects = list(objects)
        acc = Vector(0, 0, 0)
        for o in objects:
            acc = acc + o.BoundingBox().center
        return acc * (1.0 / len(objects))

    def Closed(self) -> bool:
        """Closed: an edge or wire ending where it starts, a shell whose every edge is
        shared, a solid."""
        kind = self.ShapeType()
        if kind in ("Edge", "Wire"):
            return self.IsClosed()
        if kind == "Shell":
            counts: Dict[Shape, int] = {}
            for f in self.Faces():
                for e in f.Edges():
                    counts[e] = counts.get(e, 0) + 1
            return bool(counts) and all(c >= 2 for c in counts.values())
        return kind in ("Solid", "CompSolid")

    # --- files ----------------------------------------------------------------------------

    def exportBrep(self, f) -> bool:
        """Writes a BRep file (to a path, or a binary file object)."""
        if hasattr(f, "write"):
            f.write(self.__getstate__()[0])
            return True
        from .exporters import _in_a_directory

        if not _in_a_directory(str(f)):
            return False
        native.write_brep(self.wrapped, str(f))
        return True

    @classmethod
    def importBrep(cls, f) -> "Shape":
        """Reads a BRep file (from a path, or a binary file object)."""
        if hasattr(f, "read"):
            s = Shape.__new__(Shape)
            s.__setstate__((f.read(), False))
            return Shape.cast(s.wrapped)
        return Shape.cast(native.read_brep(str(f)))

    def exportStep(self, fileName: str, unit: str = "MM", outputUnit=None, **kwargs):
        from . import exporters

        exporters.export(self, fileName, "STEP", unit=unit, outputUnit=outputUnit)

    def exportStl(self, fileName: str, tolerance: float = 1e-3, angularTolerance: float = 0.1, ascii: bool = False, relative: bool = True, **kwargs):
        from . import exporters

        if not exporters._in_a_directory(str(fileName)):
            return False
        exporters.export(self, fileName, "STL", tolerance=tolerance, angularTolerance=angularTolerance, opt={"ascii": ascii}, relative=relative)
        return True

    def export(self, fname: str, tolerance: float = 0.1, angularTolerance: float = 0.1, unit: str = "MM", outputUnit=None, opt=None) -> None:
        from . import exporters

        exporters.export(self, fname, tolerance=tolerance, angularTolerance=angularTolerance, opt=opt)


def _siblings(shapes, ctx: Shape, kind: str, level) -> set:
    if not isinstance(level, int):
        out = set()
        for lvl in level:
            out |= _siblings(shapes, ctx, kind, lvl)
        return out
    types = {s.ShapeType() for s in shapes}
    candidates = [c for t in types for c in ctx._sub(t)]
    exclude = set()

    def step(current, lvl):
        exclude.update(current)
        found = set()
        for s in current:
            shared = set(s._sub(kind))
            found.update(c for c in candidates if c not in exclude and shared & set(c._sub(kind)))
        return found if lvl == 1 else step(found, lvl - 1)

    return step(list(shapes), level)


def _cut_solid_by_solid(shape, tool):
    """`shape` with `tool` cut away; a compound of solids in groups of solids that overlap,
    a solid alone kept as it is where its box is clear of the tool's. (Kentron, given the
    compound whole, puts all its faces through its general fuse, the solids split by one
    another first; that only changes solids that overlap, which OCCT splits too: what two of
    them share comes out a solid of its own. Solids of different groups do not overlap.)"""
    if shape.kind() != "Compound":
        return native.difference(shape, tool)
    solids = native.children(shape)
    if len(solids) < 2 or any(s.kind() != "Solid" for s in solids):
        return native.difference(shape, tool)
    gap = 1e-6
    boxes = [native.bounding_box(solid, False) for solid in solids]
    apart = lambda a, b: any(a[0][k] > b[1][k] + gap or b[0][k] > a[1][k] + gap for k in range(3))

    def overlap(i, j):
        """Whether solids i and j share more than a sliver `gap` thick across the larger box."""
        common = native.intersection(solids[i], solids[j])
        volume = sum(native.volume(c)[0] for c in (native.children(common) if common.kind() == "Compound" else [common]))
        size = max(boxes[k][1][d] - boxes[k][0][d] for k in (i, j) for d in range(3))
        return volume > gap * size * size

    # The groups (of the solids' indices): the boxes swept along x, those that meet checked.
    group = list(range(len(solids)))
    def root(i):
        while group[i] != i:
            group[i] = group[group[i]]
            i = group[i]
        return i
    open_ = []
    for i in sorted(range(len(solids)), key=lambda i: boxes[i][0][0]):
        open_ = [j for j in open_ if boxes[j][1][0] + gap >= boxes[i][0][0]]
        for j in open_:
            if root(i) != root(j) and not apart(boxes[i], boxes[j]) and overlap(i, j):
                group[root(i)] = root(j)
        open_.append(i)
    groups = {}
    for i in range(len(solids)):
        groups.setdefault(root(i), []).append(i)
    if len(groups) == 1:
        return native.difference(shape, tool)
    box = native.bounding_box(tool, False)
    out = []
    for members in groups.values():
        if len(members) == 1 and apart(boxes[members[0]], box):
            out.append(solids[members[0]])
            continue
        part = solids[members[0]] if len(members) == 1 else native.make_compound([solids[i] for i in members])
        cut = native.difference(part, tool)
        out.extend(native.children(cut) if cut.kind() == "Compound" else [cut])
    return native.make_compound(out)


def _split_pieces(shape: Shape, splitters) -> List[Shape]:
    """The pieces of `shape` cut by `splitters`: the cells inside it, of its dimension."""
    tools = [t for s in splitters for t in (list(s) if s.ShapeType() == "Compound" else [s])]
    if shape.ShapeType() == "Edge" and all(t.ShapeType() == "Vertex" for t in tools):
        return _split_edge(shape, tools)
    dim = Shape._dimension(shape)
    pieces = []
    for leaf in (list(shape) if shape.ShapeType() == "Compound" else [shape]):
        for inside, h in native.cells([leaf.wrapped] + [t.wrapped for t in tools]):
            cell = Shape.cast(h)
            if inside[0] and Shape._dimension(cell) == dim:
                pieces.append(cell)
    return pieces


def _split_edge(edge: "Edge", vertices) -> List[Shape]:
    """The edge cut at the vertices on it."""
    first, last = edge.bounds()
    points = [v.Center() for v in vertices]
    us = edge.params(points)
    cuts = sorted(u for u, p in zip(us, points) if (edge.positionAt(u, "parameter") - p).Length <= 1e-6 and first + 1e-9 < u < last - 1e-9)
    if not cuts:
        return [edge]
    bounds = [first] + cuts + [last]
    return [edge.trim(a, b) for a, b in zip(bounds, bounds[1:]) if b - a > 1e-12]


def _is_scalars(args) -> bool:
    return all(isinstance(a, (int, float)) for a in args)


def _is_vector_like(a) -> bool:
    return isinstance(a, Vector) or (isinstance(a, (tuple, list)) and len(a) in (2, 3) and _is_scalars(a))


def _location_of(*args, **kwargs) -> Location:
    if kwargs:
        # CadQuery's (x, y, z, rx, ry, rz) form.
        return Location(*args, **kwargs)
    if len(args) == 1 and isinstance(args[0], Location):
        return args[0]
    if len(args) == 1 and _is_vector_like(args[0]):
        return Location(Vector(args[0]))
    # CadQuery's (x, y, z, rx, ry, rz) given in order.
    if args and len(args) <= 6 and _is_scalars(args):
        return Location(**dict(zip(("x", "y", "z", "rx", "ry", "rz"), args)))
    return Location(*args)


def _matrix_rows(m) -> list:
    rows = getattr(m, "rows", m)
    return [list(map(float, r))[:4] for r in rows[:3]]


def _is_rigid(rows) -> bool:
    """A rotation and translation (a placement), no reflection or scaling."""
    m = [r[:3] for r in rows]
    for i in range(3):
        for j in range(3):
            d = sum(m[k][i] * m[k][j] for k in range(3))
            if abs(d - (1.0 if i == j else 0.0)) > 1e-9:
                return False
    det = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
    return det > 0


class Vertex(Shape):
    def __init__(self, obj, forConstruction: bool = False):
        super().__init__(obj)
        self.forConstruction = forConstruction

    @property
    def X(self) -> float:
        return native.vertex_point(self.wrapped)[0]

    @property
    def Y(self) -> float:
        return native.vertex_point(self.wrapped)[1]

    @property
    def Z(self) -> float:
        return native.vertex_point(self.wrapped)[2]

    def toTuple(self):
        return native.vertex_point(self.wrapped)

    def Center(self) -> Vector:
        return Vector(native.vertex_point(self.wrapped))

    @classmethod
    def makeVertex(cls, x: float, y: float, z: float) -> "Vertex":
        return Shape.cast(native.make_vertex((x, y, z)))


class Mixin1D:
    """Edges and wires as curves (CadQuery's `Mixin1D`)."""

    def _info(self):
        return native.curve_info(self.wrapped)

    def bounds(self) -> Tuple[float, float]:
        """The parameter range."""
        first, last, _, _, _ = self._info()
        return first, last

    _bounds = bounds

    def _derivatives(self, us) -> list:
        return native.curve_derivatives(self.wrapped, [float(u) for u in us])

    def startPoint(self) -> Vector:
        return self._ends()[0]

    def endPoint(self) -> Vector:
        return self._ends()[1]

    def _ends(self) -> Tuple[Vector, Vector]:
        """The first and last points along the edge (as oriented) or wire."""
        if self.ShapeType() == "Edge":
            a, b = native.edge_points(self.wrapped, [0.0, 1.0])
            return _new(*a), _new(*b)
        edges = native.ordered_edges(self.wrapped)
        return _new(*native.edge_points(edges[0], [0.0])[0]), _new(*native.edge_points(edges[-1], [1.0])[0])

    def paramAt(self, d) -> float:
        """The parameter at fraction `d` of the length (or of the point nearest `d`, a
        Vector)."""
        if isinstance(d, Vector):
            return native.curve_parameters_of_points(self.wrapped, [d.toTuple()])[0]
        return native.curve_parameters_at_lengths(self.wrapped, [float(d)])[0]

    def params(self, pts, tol: float = 1e-6) -> List[float]:
        return native.curve_parameters_of_points(self.wrapped, [_tuple(p) for p in pts])

    def paramsLength(self, locations) -> List[float]:
        return native.curve_parameters_at_lengths(self.wrapped, [float(d) for d in locations])

    def _param(self, d: float, mode: str) -> float:
        return self.paramAt(d) if mode == "length" else d

    def tangentAt(self, locationParam: float = 0.5, mode: str = "length") -> Vector:
        return Vector(self._derivatives([self._param(locationParam, mode)])[0][1]).normalized()

    def tangents(self, locations, mode: str = "length") -> List[Vector]:
        us = self.paramsLength(locations) if mode == "length" else list(locations)
        return [Vector(d[1]).normalized() for d in self._derivatives(us)]

    def normal(self) -> Vector:
        """The normal of the plane the curve lies in."""
        n = native.curve_plane_normal(self.wrapped)
        if n is None:
            raise ValueError("Normal not defined")
        n = Vector(n)
        if self.geomType() not in ("CIRCLE", "ELLIPSE"):
            # A plane through the points has no side: its normal's largest component made
            # positive (as OCCT's plane finder tends to give it).
            if max(n, key=abs) < 0:
                n = -n
        return n

    def Center(self) -> Vector:
        return Vector(native.length(self.wrapped)[1])

    def Length(self) -> float:
        return native.curve_length(self.wrapped)

    def radius(self) -> float:
        """The radius of the (first) edge's circle."""
        edges = [self] if self.ShapeType() == "Edge" else list(self)
        circle = native.edge_circle(edges[0].wrapped)
        if circle is None:
            raise ValueError("Shape could not be reduced to a circle")
        return circle[0]

    def IsClosed(self) -> bool:
        """Closed topologically (its ends one vertex, as OCCT's BRep_Tool::IsClosed), or with
        its ends at one point."""
        if native.is_closed(self.wrapped):
            return True
        a, b = self._ends()
        return (a - b).Length <= 1e-7

    def positionAt(self, d: float, mode: str = "length") -> Vector:
        return Vector(self._derivatives([self._param(d, mode)])[0][0])

    def positions(self, ds, mode: str = "length") -> List[Vector]:
        us = self.paramsLength(ds) if mode == "length" else list(ds)
        return [Vector(d[0]) for d in self._derivatives(us)]

    def sample(self, n) -> Tuple[List[Vector], List[float]]:
        """Points spread evenly along the curve (`n` of them, or as many as keep within
        deflection `n`), and their parameters."""
        first, last, closed, _, _ = self._info()
        if not isinstance(n, int):
            count = max(2, int(math.ceil(self.Length() / max(math.sqrt(8 * n * max(self.Length(), 1e-9)), 1e-9))) + 1)
        else:
            count = n + 1 if closed else n
        fractions = [i / (count - 1) for i in range(count)]
        us = self.paramsLength(fractions)
        if closed:
            us = us[:-1]
        return [Vector(d[0]) for d in self._derivatives(us)], us

    def locationAt(self, d: float, mode: str = "length", frame: str = "frenet", planar: bool = False) -> Location:
        """The location at `d`: Z along the tangent, X along the curve's normal (or, if
        `planar`, Z up and X across the tangent)."""
        p, d1, d2 = self._derivatives([self._param(d, mode)])[0]
        t = Vector(d1).normalized()
        if planar:
            z = Vector(0, 0, 1)
            return Location(Plane(Vector(p), z.cross(t).normalized(), z))
        n = Vector(d2) - t * Vector(d2).dot(t)
        x = n.normalized() if n.Length > 1e-9 * max(Vector(d2).Length, 1.0) else _default_x(t)
        return Location(Plane(Vector(p), x, t))

    def locations(self, ds, mode: str = "length", frame: str = "frenet", planar: bool = False) -> List[Location]:
        return [self.locationAt(d, mode, frame, planar) for d in ds]

    def curvatureAt(self, d: float, mode: str = "length", resolution: float = 1e-6) -> float:
        _, d1, d2 = self._derivatives([self._param(d, mode)])[0]
        d1, d2 = Vector(d1), Vector(d2)
        speed = d1.Length
        return d1.cross(d2).Length / speed ** 3 if speed > 0 else 0.0

    def curvatures(self, ds, mode: str = "length", resolution: float = 1e-6) -> List[float]:
        return [self.curvatureAt(d, mode, resolution) for d in ds]

    def project(self, face: "Face", d, closest: bool = True):
        """The curve projected onto the face's surface along `d` (the projection nearest it):
        a wire."""
        return Shape.cast(native.project_along(self.wrapped, face.wrapped, Vector(d).toTuple()))

    Closed = IsClosed


# (pykentron's earlier name.)
_OneDimensional = Mixin1D


def _ordered_edges(edges):
    return edges


def _rotated_about(v: Vector, axis: Vector, angle: float) -> Vector:
    return v.transform(Matrix(_rotation_rows(axis, angle)))


class Edge(Shape, Mixin1D):
    def geomType(self) -> str:
        return native.curve_type(self.wrapped)

    def _geomAdaptor(self) -> "Edge":
        return self

    def close(self) -> "Shape":
        return self if self.IsClosed() else Wire.assembleEdges([self]).close()

    def arcCenter(self) -> Vector:
        g = self.geomType()
        if g == "CIRCLE":
            return Vector(native.edge_circle(self.wrapped)[1])
        if g == "ELLIPSE":
            # Halfway between points half a period apart.
            u = self.bounds()[0]
            a, b = self._derivatives([u, u + math.pi])
            return (Vector(a[0]) + Vector(b[0])) * 0.5
        raise ValueError(f"{g} has no arc center")

    def trim(self, u0: float, u1: float) -> "Edge":
        """The edge's curve over `(u0, u1)`."""
        return Shape.cast(native.edge_trimmed(self.wrapped, u0, u1))

    def radius(self) -> float:
        circle = native.edge_circle(self.wrapped)
        if circle is None:
            raise ValueError("Shape could not be reduced to a circle")
        return circle[0]

    @classmethod
    def makeLine(cls, v1, v2) -> "Edge":
        return Shape.cast(native.make_line(_tuple(v1), _tuple(v2)))

    @classmethod
    def makeCircle(cls, radius: float, pnt=(0, 0, 0), dir=(0, 0, 1), angle1: float = 360, angle2: float = 360, orientation: bool = True) -> "Edge":
        pnt, dir = Vector(pnt), Vector(dir)
        if angle1 == angle2:
            return Shape.cast(native.make_circle(radius, pnt.toTuple(), dir.toTuple(), 0.0, 360.0))
        def make(r, centre, normal, x, a1, a2):
            return native.make_circle(r, centre, normal, a1, a2, x)

        return cls._arc(make, (radius,), pnt, dir, _default_x(dir.normalized()), angle1, angle2, orientation)

    @staticmethod
    def _arc(make, radii, pnt: Vector, dir: Vector, x: Vector, angle1: float, angle2: float, sense: bool) -> "Edge":
        """An arc of a conic (GC_MakeArcOf...): counter-clockwise from `angle1` to `angle2`
        (degrees about `dir` from `x`), run backwards unless `sense` (on the conic reversed:
        its normal reversed, angles negated)."""
        if not sense:
            dir, angle1, angle2 = -dir, -angle2, -angle1
        a1 = angle1 % 360.0
        a2 = a1 + (angle2 - angle1) % 360.0
        if a2 <= a1:
            a2 += 360.0
        return Shape.cast(make(*radii, pnt.toTuple(), dir.toTuple(), x.toTuple(), a1, a2))

    @classmethod
    def makeEllipse(cls, x_radius: float, y_radius: float, pnt=(0, 0, 0), dir=(0, 0, 1), xdir=(1, 0, 0), angle1: float = 360.0, angle2: float = 360.0, sense: int = 1) -> "Edge":
        pnt, dir, xdir = Vector(pnt), Vector(dir).normalized(), Vector(xdir)
        xdir = (xdir - dir * xdir.dot(dir)).normalized()
        correction = 0.0
        major, minor = x_radius, y_radius
        if y_radius > x_radius:
            # The major axis along y: the ellipse turned a quarter about its normal.
            major, minor, correction = y_radius, x_radius, 90.0
            xdir = dir.cross(xdir)

        def make(a, b, centre, normal, x, a1, a2):
            return native.make_ellipse(a, b, centre, normal, x, a1, a2)

        if angle1 == angle2:
            return Shape.cast(native.make_ellipse(major, minor, pnt.toTuple(), dir.toTuple(), xdir.toTuple(), 0.0, 360.0))
        return cls._arc(make, (major, minor), pnt, dir, xdir, angle1 - correction, angle2 - correction, sense == 1)

    @classmethod
    def makeSpline(cls, listOfVector, tangents=None, periodic: bool = False, parameters=None, scale: bool = True, tol: float = 1e-6) -> "Edge":
        """The spline through the points, as OCCT interpolates it."""
        pts = [_tuple(v) for v in listOfVector]
        if parameters is not None:
            parameters = [float(p) for p in parameters]
            if len(parameters) != len(pts) + periodic:
                raise ValueError(
                    "There must be one parameter for each interpolation point "
                    "(plus one if periodic), or none specified. Parameter count: "
                    f"{len(parameters)}, point count: {len(pts)}"
                )
        tgts = None
        if tangents:
            tangents = list(tangents)
            if not (len(tangents) == 2 and len(pts) != 2) and len(tangents) != len(pts):
                raise ValueError(
                    "There must be one tangent for each interpolation point, "
                    "or just two end point tangents. Tangent count: "
                    f"{len(tangents)}, point count: {len(pts)}"
                )
            tgts = [None if t is None else _tuple(t) for t in tangents]
        return Shape.cast(native.make_spline(pts, tgts, parameters, periodic, scale, tol))

    @classmethod
    def makeSplineApprox(cls, listOfVector, tol: float = 1e-3, smoothing=None, minDeg: int = 1, maxDeg: int = 6) -> "Edge":
        smoothing = [float(w) for w in smoothing] if smoothing else None
        return Shape.cast(native.make_spline_approx([_tuple(v) for v in listOfVector], tol, minDeg, maxDeg, smoothing))

    @classmethod
    def makeThreePointArc(cls, v1, v2, v3) -> "Edge":
        return Shape.cast(native.make_three_point_arc(_tuple(v1), _tuple(v2), _tuple(v3)))

    @classmethod
    def makeTangentArc(cls, v1, v2, v3) -> "Edge":
        """The arc from `v1`, leaving along `v2`, to `v3`."""
        p1, t, p3 = Vector(v1), Vector(v2), Vector(v3)
        d = p3 - p1
        n = t.cross(d)
        if n.Length <= 1e-12 * t.Length * d.Length:
            raise ValueError("a tangent arc along its chord")
        n = n.normalized()
        u = n.cross(t).normalized()  # towards the centre
        r = d.dot(d) / (2 * u.dot(d))
        c = p1 + u * r
        a, b = p1 - c, p3 - c
        # The angle from a to b turning the way the tangent leaves (about n).
        sweep = math.atan2(a.cross(b).dot(n), a.dot(b)) % (2 * math.pi)
        mid = c + _rotated_about(a, n, sweep / 2)
        return cls.makeThreePointArc(p1, mid, p3)

    @classmethod
    def makeBezier(cls, points) -> "Edge":
        return Shape.cast(native.make_bezier([_tuple(p) for p in points]))


class Wire(Shape, Mixin1D):
    def _nbEdges(self) -> int:
        return len(native.ordered_edges(self.wrapped))

    def _geomAdaptor(self) -> "Wire":
        return self

    def __iter__(self):
        """The edges in order along the wire."""
        return (Shape.cast(h) for h in native.ordered_edges(self.wrapped))

    def Vertices(self) -> List[Vertex]:
        """The vertices in order along the wire."""
        out = []
        edges = list(self)
        for e in edges:
            out.append(_start_vertex(e))
        if edges:
            out.append(_start_vertex(edges[-1].reverse()))
        if self.IsClosed() and len(out) > 1:
            out = out[:-1]
        return out

    @classmethod
    def assembleEdges(cls, listOfEdges: Iterable[Edge]) -> "Wire":
        """The wire of the edges (in any order and orientation, as OCCT's wire builder takes
        them: chained end to end)."""
        edges = [e for el in listOfEdges for e in (list(el) if el.ShapeType() == "Wire" else [el])]
        try:
            return Shape.cast(native.make_wire([e.wrapped for e in edges]))
        except ValueError:
            chains = _edges_to_wires(edges, 1e-7)
            if len(chains) != 1:
                raise
            return chains[0]

    @classmethod
    def combine(cls, listOfWires, tol: float = 1e-9) -> List["Wire"]:
        """The wires and edges joined into wires, end to end."""
        edges = Compound.makeCompound(listOfWires).Edges()
        return _edges_to_wires(edges, tol)

    @classmethod
    def makePolygon(cls, listOfVertices, forConstruction: bool = False, close: bool = False) -> "Wire":
        pts = [_tuple(p) for p in listOfVertices]
        return Shape.cast(native.make_polygon(pts, close), forConstruction)

    @classmethod
    def makeCircle(cls, radius: float, center, normal) -> "Wire":
        return cls.assembleEdges([Edge.makeCircle(radius, center, normal)])

    @classmethod
    def makeEllipse(cls, x_radius: float, y_radius: float, center, normal, xDir, angle1: float = 360.0, angle2: float = 360.0, rotation_angle: float = 0.0, closed: bool = True) -> "Wire":
        e = Edge.makeEllipse(x_radius, y_radius, center, normal, xDir, angle1, angle2)
        if angle1 != angle2 and closed:
            w = cls.assembleEdges([e, Edge.makeLine(e.endPoint(), e.startPoint())])
        else:
            w = cls.assembleEdges([e])
        if rotation_angle != 0.0:
            w = w.rotate(center, Vector(center) + Vector(normal), rotation_angle)
        return w

    @classmethod
    def makeHelix(cls, pitch: float, height: float, radius: float, center=(0, 0, 0), dir=(0, 0, 1), angle: float = 360.0, lefthand: bool = False) -> "Wire":
        """A helix on the cylinder of `radius` about `dir` through `center` (on a cone of
        half-angle `angle` degrees if not 360), rising `pitch` a turn, to `height`.

        As OCCT's (a line in the surface's parameters, its 3D curve a spline), parametrised by
        the length along that line: a spline through points with their derivatives."""
        c, z = Vector(center), Vector(dir).normalized()
        x = _default_x(z)
        y = z.cross(x)
        cone = angle != 360.0
        half = math.radians(angle) if cone else 0.0
        sign = -1.0 if lefthand else 1.0
        speed = math.hypot(2 * math.pi, pitch)
        turns = height / pitch
        total = turns * speed
        steps = max(8, int(math.ceil(abs(turns) * 32)))

        def at(s):
            # The surface's (u, v) along the line, and their rates.
            u, v = sign * 2 * math.pi * s / speed, pitch * s / speed
            du, dv = sign * 2 * math.pi / speed, pitch / speed
            if cone:
                r, h = radius + v * math.sin(half), v * math.cos(half)
                dr, dh = dv * math.sin(half), dv * math.cos(half)
            else:
                r, h, dr, dh = radius, v, 0.0, dv
            radial = x * math.cos(u) + y * math.sin(u)
            around = x * -math.sin(u) + y * math.cos(u)
            return c + radial * r + z * h, radial * dr + around * (r * du) + z * dh

        params = [total * i / steps for i in range(steps + 1)]
        samples = [at(s) for s in params]
        e = native.make_spline([p.toTuple() for p, _ in samples], [d.toTuple() for _, d in samples], params, False, False, 1e-9)
        return cls.assembleEdges([Shape.cast(e)])

    def close(self) -> "Wire":
        if self.IsClosed():
            return self
        return Wire.combine((self, Edge.makeLine(self.endPoint(), self.startPoint())))[0]

    def stitch(self, other: "Wire") -> "Wire":
        return Wire.assembleEdges(list(self) + list(other))

    def offset2D(self, d: float, kind: str = "arc") -> List["Wire"]:
        result = Shape.cast(native.offset_wire(self.wrapped, d, kind))
        if result.wrapped.kind() == "Wire":
            return [result]
        return result.Wires()

    def fillet2D(self, radius: float, vertices: Iterable[Vertex]) -> "Wire":
        """The wire with its corners at `vertices` rounded (between straight edges)."""
        return _round_corners(self, radius, list(vertices), lambda a, b, p: _fillet_corner(a, b, p, radius))

    def chamfer2D(self, d: float, vertices: Iterable[Vertex]) -> "Wire":
        """The wire with its corners at `vertices` cut back by `d` along both edges."""
        return _round_corners(self, d, list(vertices), lambda a, b, p: _chamfer_corner(a, b, p, d))

    def fillet(self, radius: float, vertices=None) -> "Wire":
        """Corners rounded to `radius` (at `vertices`, else all), between straight edges."""
        return _round_corners(self, radius, vertices, lambda a, b, p: _fillet_corner(a, b, p, radius))


def _start_vertex(edge: Edge) -> Vertex:
    """The vertex an edge (as oriented) starts at."""
    start = Vector(native.edge_points(edge.wrapped, [0.0])[0])
    vs = edge._sub("Vertex")
    return min(vs, key=lambda v: (v.Center() - start).Length)


def _edges_to_wires(edges: List[Edge], tol: float) -> List[Wire]:
    """The edges joined into wires, end to end (ShapeAnalysis_FreeBounds)."""
    remaining = list(edges)
    wires = []
    tol = max(tol, 1e-7)
    while remaining:
        chain = [remaining.pop(0)]
        grown = True
        while grown:
            grown = False
            head, tail = chain[0].startPoint(), chain[-1].endPoint()
            for k, e in enumerate(remaining):
                s, t = e.startPoint(), e.endPoint()
                if (s - tail).Length <= tol:
                    chain.append(remaining.pop(k))
                elif (t - tail).Length <= tol:
                    chain.append(remaining.pop(k).reverse())
                elif (t - head).Length <= tol:
                    chain.insert(0, remaining.pop(k))
                elif (s - head).Length <= tol:
                    chain.insert(0, remaining.pop(k).reverse())
                else:
                    continue
                grown = True
                break
        wires.append(Shape.cast(native.make_wire([e.wrapped for e in _slits_closed(_loops_one_way(chain, tol), tol)], tol)))
    return wires


def _slits_closed(chain: List[Edge], tol: float) -> List[Edge]:
    """`chain` less each edge followed by the same edge back (a slit, as where two mirrored
    halves joined both run along the mirror line): the two bound nothing, and would leave
    faces back to back inside a solid made from them."""

    def back(a: Edge, b: Edge) -> bool:
        return (
            (a.startPoint() - b.endPoint()).Length <= tol
            and (a.endPoint() - b.startPoint()).Length <= tol
            and (a.positionAt(0.5) - b.positionAt(0.5)).Length <= tol
        )

    out: List[Edge] = []
    for e in chain:
        if out and back(out[-1], e):
            out.pop()
        else:
            out.append(e)
    # (Round the end of a closed chain too.)
    while len(out) > 2 and back(out[-1], out[0]):
        out = out[1:-1]
    return out if len(out) >= 1 else chain


def _lines_joined(wire: "Wire") -> "Wire":
    """A closed wire with its straight edges in line joined (across its start too), its
    curved edges kept as they are."""
    if not wire.IsClosed():
        return wire
    edges = [Shape.cast(e) for e in native.ordered_edges(wire.wrapped)]
    n = len(edges)
    if n < 2:
        return wire

    def continues(a, b) -> bool:
        # (Two lines, the second going on straight from the first.)
        if a.geomType() != "LINE" or b.geomType() != "LINE":
            return False
        da, db = a.endPoint() - a.startPoint(), b.endPoint() - b.startPoint()
        return da.cross(db).Length <= 1e-9 * da.Length * db.Length and da.dot(db) > 0

    # From an edge the one before does not go on into, runs of lines in line as one each.
    start = next((i for i in range(n) if not continues(edges[i - 1], edges[i])), None)
    if start is None:
        return wire
    order = edges[start:] + edges[:start]
    joined, i = [], 0
    while i < n:
        j = i
        while j + 1 < n and continues(order[j], order[j + 1]):
            j += 1
        joined.append(order[i] if j == i else Edge.makeLine(order[i].startPoint(), order[j].endPoint()))
        i = j + 1
    if len(joined) == n:
        return wire
    return Wire.assembleEdges(joined)


def _box_within(inner, outer, tol: float = 1e-6) -> bool:
    """Whether box `inner` lies within box `outer` (up to `tol`)."""
    return (inner.xmin >= outer.xmin - tol and inner.ymin >= outer.ymin - tol and inner.zmin >= outer.zmin - tol
            and inner.xmax <= outer.xmax + tol and inner.ymax <= outer.ymax + tol and inner.zmax <= outer.zmax + tol)


def _has_corners(wire: "Wire") -> bool:
    """Whether the wire's edges meet at an angle somewhere (at its start too, if closed)."""
    edges = wire.Edges()
    pairs = list(zip(edges, edges[1:])) + ([(edges[-1], edges[0])] if wire.IsClosed() and len(edges) > 1 else [])
    for a, b in pairs:
        ta, tb = a.tangentAt(1.0), b.tangentAt(0.0)
        if ta.cross(tb).Length > 1e-6 or ta.dot(tb) < 0:
            return True
    return False


def _loops_one_way(chain: List[Edge], tol: float) -> List[Edge]:
    """`chain` with the loops after its first (where it passes its start again, as two
    mirrored halves joined do) turned to wind as the first does: they then add to the area
    it bounds rather than taking it away."""
    start = chain[0].startPoint()
    cuts = [i + 1 for i, e in enumerate(chain) if (e.endPoint() - start).Length <= tol]
    if len(cuts) < 2:
        return chain
    loops = [chain[a:b] for a, b in zip([0] + cuts[:-1], cuts)] + ([chain[cuts[-1]:]] if cuts[-1] < len(chain) else [])

    def winding(loop):
        # (Newell's: twice the vector area of a polygon through points along the edges.)
        pts = [e.positionAt(t, mode="parameter") for e in loop for t in (0.0, 0.25, 0.5, 0.75)]
        total = Vector(0, 0, 0)
        for a, b in zip(pts, pts[1:] + pts[:1]):
            total = total + a.cross(b)
        return total

    first = winding(loops[0])
    out = list(loops[0])
    for loop in loops[1:]:
        if winding(loop).dot(first) < 0:
            loop = [e.reverse() for e in reversed(loop)]
        out.extend(loop)
    return out


def _line_corner(a: Edge, b: Edge):
    """The straight edges meeting at a corner: (corner point, direction back along `a`,
    direction on along `b`); None unless both are straight."""
    if a.geomType() != "LINE" or b.geomType() != "LINE":
        return None
    p = a.endPoint()
    return p, (a.startPoint() - p).normalized(), (b.endPoint() - p).normalized()


def _fillet_corner(a: Edge, b: Edge, corner, radius: float):
    """`a` and `b` trimmed back and the arc of `radius` between them."""
    p, u, v = corner
    half = u.getAngle(v) / 2
    if half <= 1e-9 or abs(half - math.pi / 2) <= 1e-12:
        return None
    back = radius / math.tan(half)
    if back > (a.startPoint() - p).Length + 1e-9 or back > (b.endPoint() - p).Length + 1e-9:
        raise ValueError("the fillet is too large for its edges")
    pa, pb = p + u * back, p + v * back
    bisector = (u + v).normalized()
    mid = p + bisector * (radius / math.sin(half) - radius)
    return pa, Edge.makeThreePointArc(pa, mid, pb), pb


def _chamfer_corner(a: Edge, b: Edge, corner, d: float):
    p, u, v = corner
    pa, pb = p + u * d, p + v * d
    return pa, Edge.makeLine(pa, pb), pb


def _round_corners(wire: Wire, size: float, vertices, corner_edge) -> Wire:
    """`wire`'s corners at `vertices` (all if none) replaced by `corner_edge`'s edges."""
    edges = list(wire)
    closed = wire.IsClosed()
    wanted = [Vector(v.Center()) for v in vertices] if vertices else None
    n = len(edges)
    starts = [e.startPoint() for e in edges]
    ends = [e.endPoint() for e in edges]
    inserts = {}
    for i in range(n if closed else n - 1):
        j = (i + 1) % n
        p = ends[i]
        if wanted is not None and not any((p - w).Length <= 1e-6 for w in wanted):
            continue
        corner = _line_corner(Edge.makeLine(starts[i], ends[i]) if edges[i].geomType() == "LINE" else edges[i], Edge.makeLine(starts[j], ends[j]) if edges[j].geomType() == "LINE" else edges[j])
        if corner is None:
            if wanted is not None:
                raise NotImplementedError("2D fillets and chamfers between curved edges")
            continue
        made = corner_edge(edges[i], edges[j], corner)
        if made is None:
            continue
        pa, e, pb = made
        ends[i], starts[j] = pa, pb
        inserts[i] = e
    out = []
    for i in range(n):
        if (starts[i] - ends[i]).Length > 1e-9:
            out.append(Edge.makeLine(starts[i], ends[i]) if edges[i].geomType() == "LINE" else edges[i])
        if i in inserts:
            out.append(inserts[i])
    return Wire.assembleEdges(out)


class Face(Shape):
    def geomType(self) -> str:
        return native.surface_type(self.wrapped)

    def _geomAdaptor(self) -> "Face":
        return self

    def uvBounds(self) -> Tuple[float, float, float, float]:
        return native.face_uv_bounds(self.wrapped)

    _uvBounds = uvBounds

    def paramAt(self, pt) -> Tuple[float, float]:
        return native.face_parameters(self.wrapped, [_tuple(pt)])[0]

    def params(self, pts, tol: float = 1e-9) -> Tuple[List[float], List[float]]:
        uvs = native.face_parameters(self.wrapped, [_tuple(p) for p in pts])
        return [u for u, _ in uvs], [v for _, v in uvs]

    def _at(self, uvs) -> list:
        return native.face_points(self.wrapped, [(float(u), float(v)) for u, v in uvs])

    def positionAt(self, u: float, v: float) -> Vector:
        return Vector(self._at([(u, v)])[0][0])

    def positions(self, uvs) -> List[Vector]:
        return [Vector(r[0]) for r in self._at(uvs)]

    def locationAt(self, u: float, v: float) -> Location:
        p, _, du, dv = self._at([(u, v)])[0]
        du, dv = Vector(du), Vector(dv)
        return Location(Plane(Vector(p), du, du.cross(dv)))

    def project(self, other: "Face", d) -> "Face":
        """The face projected onto `other`'s surface along `d`: a face there bounded by its
        wires' projections (as CadQuery's `Face.project`)."""
        return Shape.cast(native.project_along(self.wrapped, other.wrapped, Vector(d).toTuple()))

    def normalAt(self, *args, **kwargs):
        """The unit normal (as the face is used) at the point nearest `locationVector` (the
        middle of the parameter range if none); or, given parameters `u, v`, the normal and
        the point there."""
        if len(args) == 2 and all(isinstance(a, (int, float)) for a in args):
            p, n, _, _ = self._at([args])[0]
            return Vector(n), Vector(p)
        location = args[0] if args else kwargs.get("locationVector")
        if location is None:
            u0, u1, v0, v1 = self.uvBounds()
            uv = (0.5 * (u0 + u1), 0.5 * (v0 + v1))
        else:
            uv = self.paramAt(location)
        return Vector(self._at([uv])[0][1])

    def normals(self, us, vs) -> Tuple[List[Vector], List[Vector]]:
        rows = self._at(list(zip(us, vs)))
        return [Vector(r[1]) for r in rows], [Vector(r[0]) for r in rows]

    def tangentAt(self, u: float, v: float) -> Tuple[Vector, Vector, Vector]:
        p, _, du, dv = self._at([(u, v)])[0]
        return Vector(du).normalized(), Vector(dv).normalized(), Vector(p)

    def Center(self) -> Vector:
        return Vector(native.area(self.wrapped)[1])

    def outerWire(self) -> Wire:
        return self.Wires()[0]

    def innerWires(self) -> List[Wire]:
        return self.Wires()[1:]

    @classmethod
    def makeFromWires(cls, outerWire: Wire, innerWires: Sequence[Wire] = ()) -> "Shape":
        """The face bounded by `outerWire` with `innerWires` as holes. As CadQuery's (whose face
        fixing splits it so), "inner" wires clear of the outer one bound faces of their own
        (with the wires inside them as holes), the result then a compound of faces."""
        # (A wire given again, the same shape, bounds nothing more: as where a script queues
        # its pending wire twice.)
        innerWires = [w for i, w in enumerate(innerWires) if not w.isSame(outerWire) and not any(w.isSame(v) for v in innerWires[:i])]
        if innerWires and not outerWire.IsClosed():
            raise ValueError("Cannot build face(s): outer wire is not closed")
        for w in innerWires:
            if not w.IsClosed():
                raise ValueError("Cannot build face(s): inner wire is not closed")

        # A wire is a hole of a group when it lies within the group's outer wire: a point on
        # it at no distance from the face the outer wire bounds (bounding boxes overlap for
        # wires side by side, rotated rectangles for one).
        bounded = {}

        # (Only where the outer wire's box holds the wire's: most wires side by side are
        # told apart by their boxes, without a face.)
        points = {}

        def within(outer, w) -> bool:
            o, b = outer.BoundingBox(), w.BoundingBox()
            tol = 1e-7 * (1.0 + max(o.DiagonalLength, 1.0))
            if b.xmin < o.xmin - tol or b.ymin < o.ymin - tol or b.zmin < o.zmin - tol or b.xmax > o.xmax + tol or b.ymax > o.ymax + tol or b.zmax > o.zmax + tol:
                return False
            key = id(outer)
            if key not in bounded:
                bounded[key] = Shape.cast(native.make_face(outer.wrapped, []))
            face = bounded[key]
            if id(w) not in points:
                points[id(w)] = w.positionAt(0.5)
            p = points[id(w)]
            return face.distance(Vertex.makeVertex(p.x, p.y, p.z)) <= tol

        groups = [[outerWire]]
        for w in innerWires:
            for g in groups:
                if within(g[0], w):
                    g.append(w)
                    break
            else:
                # (A wire around groups' outer wires, given after them, is the outer one: those
                # become its holes, their own holes islands bounding faces of their own.)
                around = [g for g in groups if within(w, g[0])]
                if around:
                    first = groups.index(around[0])
                    groups = [g for g in groups if g not in around]
                    groups.insert(first, [w] + [g[0] for g in around])
                    groups.extend([h] for g in around for h in g[1:])
                else:
                    groups.append([w])
        faces = [Shape.cast(native.make_face(g[0].wrapped, [w.wrapped for w in g[1:]])) for g in groups]
        return faces[0] if len(faces) == 1 else Compound.makeCompound(faces)

    @classmethod
    def makePlane(cls, length=None, width=None, basePnt=(0, 0, 0), dir=(0, 0, 1)) -> "Face":
        """The plane face through `basePnt` normal to `dir`, `width` along its X and
        `length` along its Y (as large as kentron takes for "infinite" if not given)."""
        if length and width:
            u, v = width / 2, length / 2
        else:
            u = v = 1e4
        # (Of negative sizes, as of positive: CadQuery's takes them so.)
        u, v = abs(u), abs(v)
        return Shape.cast(native.make_plane_face(_tuple(basePnt), _tuple(dir), -u, u, -v, v))

    @classmethod
    def makeRuledSurface(cls, edgeOrWire1, edgeOrWire2) -> "Shape":
        w1, w2 = (Wire.assembleEdges([e]) if e.ShapeType() == "Edge" else e for e in (edgeOrWire1, edgeOrWire2))
        shell = Shape.cast(native.make_loft([w1.wrapped, w2.wrapped], False, True))
        faces = shell.Faces()
        return faces[0] if len(faces) == 1 and edgeOrWire1.ShapeType() == "Edge" else shell

    @classmethod
    def makeNSidedSurface(cls, edges, constraints=(), continuity=None, degree: int = 3, *args, **kwargs) -> "Face":
        """The face spanning the closed loop of `edges` (edges or wires) through the points of
        `constraints`, of least bending energy (OCCT's filling: `degree` the energy's order;
        the other settings, its approximation's, are not needed)."""
        points = []
        for c in constraints:
            if isinstance(c, Shape):
                raise NotImplementedError("filling surfaces constrained by curves within them")
            points.append(_tuple(c))
        return Shape.cast(native.make_filling([e.wrapped for e in edges], points, int(degree)))

    @classmethod
    def makeSplineApprox(cls, points, tol: float = 1e-2, smoothing=None, minDeg: int = 1, maxDeg: int = 3) -> "Face":
        """The spline face near the grid of points (rows along u), within `tol`. (Smoothing
        weights are not used: the surface passes within the tolerance as without them.)"""
        grid = [[_tuple(p) for p in row] for row in points]
        return Shape.cast(native.make_spline_surface_approx(grid, tol, minDeg, maxDeg))

    def fillet2D(self, radius: float, vertices: Iterable[Vertex]) -> "Face":
        """The face with its corners at `vertices` rounded (between straight edges)."""
        return self._corners(vertices, lambda a, b, p: _fillet_corner(a, b, p, radius), radius)

    def chamfer2D(self, d: float, vertices: Iterable[Vertex]) -> "Face":
        """The face with its corners at `vertices` cut back by `d` along both edges."""
        return self._corners(vertices, lambda a, b, p: _chamfer_corner(a, b, p, d))

    def _corners(self, vertices, corner_edge, radius=None) -> "Face":
        vertices = list(vertices)
        wires = []
        for w in self.Wires():
            mine = [v for v in vertices if any(v.isSame(x) for x in w.Vertices())]
            if not mine:
                wires.append(w)
                continue
            try:
                wires.append(_round_corners(w, 0.0, mine, corner_edge))
            except NotImplementedError:
                # (Between curved edges: rounded by the kernel, along the curves.)
                if radius is None:
                    raise
                wires.append(Shape.cast(native.fillet_wire_2d(w.wrapped, radius, [v.wrapped for v in mine])))
        return Face.makeFromWires(wires[0], wires[1:])

    def thicken(self, thickness: float) -> "Solid":
        """The solid between the face and the face offset by `thickness` along its normal."""
        return _thickened(self, thickness)

    def trim(self, *args) -> "Face":
        """The face's surface over `(u0, u1) × (v0, v1)`."""
        if len(args) in (4, 5) and all(isinstance(a, (int, float)) for a in args):
            u0, u1, v0, v1 = args[:4]
            return Shape.cast(native.face_trimmed(self.wrapped, u0, u1, v0, v1))
        raise NotImplementedError("trimming faces by wires or parameter polygons")

    def isoline(self, param: float, direction: str = "v") -> Edge:
        raise NotImplementedError("isolines")

    def addHole(self, *inner) -> "Face":
        wires = [w if w.ShapeType() == "Wire" else Wire.assembleEdges([w]) for w in inner]
        if self.geomType() != "PLANE":
            return Shape.cast(native.face_with_holes(self.wrapped, [w.wrapped for w in wires]))
        return Face.makeFromWires(self.outerWire(), self.innerWires() + wires)

    def toPln(self):
        """The plane of a planar face as OCCT's (`gp_Pln`), on its surface's frame."""
        from .gp import gp_Ax3, gp_Pln

        if self.geomType() != "PLANE":
            raise ValueError("not a planar face")
        p, n, du, _ = self._at([(0.0, 0.0)])[0]
        return gp_Pln(gp_Ax3(Vector(p).toPnt(), Vector(n).toDir(), Vector(du).toDir()))


class Shell(Shape):
    @classmethod
    def makeShell(cls, listOfFaces: Iterable[Face]) -> "Shape":
        return Shape.cast(native.sew([f.wrapped for f in listOfFaces]))


class Mixin3D:
    def isInside(self, point, tolerance: float = 1e-6) -> bool:
        """True if the point is inside the solid or on its boundary (closer to it than
        `tolerance`; a point at the tolerance, up to rounding, is outside, as in OCCT)."""
        if native.classify_point(self.wrapped, _tuple(point)) != "out":
            return True
        return native.distance(native.make_vertex(_tuple(point)), self.wrapped)[0] < tolerance * (1.0 - 1e-9)

    def fillet(self, radius: float, edgeList) -> "Shape":
        """The solid with `edgeList`'s edges rounded to `radius`."""
        edges = [e.wrapped for e in edgeList]
        return Shape.cast(native.fillet(self.wrapped, edges, radius))

    def chamfer(self, length: float, length2=None, edgeList=()) -> "Shape":
        """The solid with `edgeList`'s edges bevelled: cut back by `length` along each edge's
        first face (of this shape's faces, in order) and by `length2` (`length` if not given)
        along the other."""
        faces = self.Faces()
        spec = []
        for e in edgeList:
            first = next((f for f in faces if any(e.isSame(x) for x in f.Edges())), None)
            spec.append((e.wrapped, length, length2 or length, first.wrapped if first is not None else None))
        return Shape.cast(native.chamfer(self.wrapped, spec))

    def hollow(self, faceList, thickness: float, tolerance: float = 0.0001, kind: str = "arc") -> "Shape":
        """The solid made a shell of `thickness` (outwards if positive), open where
        `faceList`'s faces are (closed, a cavity inside, if there are none)."""
        faces = [f.wrapped for f in (faceList or [])]
        if faces:
            return Shape.cast(native.thick_solid(self.wrapped, faces, thickness, kind))
        offset = Shape.cast(native.offset_shape(self.wrapped, thickness, kind))
        mine, theirs = self.Shells()[0], offset.Shells()[0]
        outer, inner = (theirs, mine) if thickness > 0 else (mine, theirs)
        return Shape.cast(native.make_solid([outer.wrapped, inner.wrapped.reversed()]))

    def dprism(self, basis, profiles, depth=None, taper=0, upToFace=None, thruAll=True, additive=True):
        """A prism feature (BRepFeat_MakeDPrism): each profile (faces, or wires grouped into
        faces) extruded along its normal up to `upToFace` (its plane), by `depth`, or through
        everything both ways, and fused with the solid (or cut from it unless `additive`).
        Tapered prisms are not supported."""
        if taper:
            raise NotImplementedError("tapered prism features")
        profiles = list(profiles)
        if profiles and isinstance(profiles[0], Wire):
            from .occ_impl.shapes import sortWiresByBuildOrder

            profiles = [Face.makeFromWires(g[0], g[1:]) for g in sortWiresByBuildOrder(profiles)]
        result = self
        reach = 2.0 * self.BoundingBox().DiagonalLength
        for face in profiles:
            n = face.normalAt()
            if upToFace is not None and upToFace.geomType() == "PLANE":
                # Far along the normal, then only what lies on the profile's side of the
                # face's plane: the prism ends on the plane, at whatever slant it has.
                origin, normal = upToFace.Center(), upToFace.normalAt()
                size = 4.0 * (reach + face.BoundingBox().DiagonalLength + (origin - face.Center()).Length)
                towards = n if (origin - face.Center()).dot(n) >= 0 else -n
                far = Shape.cast(native.extrude(face.wrapped, (towards * size).toTuple()))
                side = 1.0 if (face.Center() - origin).dot(normal) >= 0 else -1.0
                slab = Shape.cast(native.extrude(Face.makePlane(size, size, origin, normal).wrapped, (normal * (side * size)).toTuple()))
                prism = far.intersect(slab).wrapped
            elif upToFace is not None:
                prism = native.extrude(face.wrapped, (n * (upToFace.Center() - face.Center()).dot(n)).toTuple())
            elif thruAll or depth is None:
                d = n * (reach + face.BoundingBox().DiagonalLength)
                prism = native.extrude(face.moved(Location(-d)).wrapped, (d * 2).toTuple())
            else:
                prism = native.extrude(face.wrapped, (n * depth).toTuple())
            prism = Shape.cast(prism)
            result = result.fuse(prism) if additive else result.cut(prism)
        return result


class Solid(Shape, Mixin3D):
    @staticmethod
    def isSolid(obj) -> bool:
        return hasattr(obj, "ShapeType") and (obj.ShapeType() == "Solid" or (obj.ShapeType() == "Compound" and len(obj.Solids()) > 0))

    @classmethod
    def makeSolid(cls, shell: Shell) -> "Solid":
        return Shape.cast(native.make_solid([shell.wrapped]))

    @classmethod
    def makeBox(cls, length: float, width: float, height: float, pnt=(0, 0, 0), dir=(0, 0, 1)) -> "Solid":
        return Shape.cast(native.make_box(length, width, height, _tuple(pnt), _tuple(dir)))

    @classmethod
    def makeCylinder(cls, radius: float, height: float, pnt=(0, 0, 0), dir=(0, 0, 1), angleDegrees: float = 360) -> "Solid":
        return Shape.cast(native.make_cylinder(radius, height, _tuple(pnt), _tuple(dir), angleDegrees))

    @classmethod
    def makeCone(cls, radius1: float, radius2: float, height: float, pnt=(0, 0, 0), dir=(0, 0, 1), angleDegrees: float = 360) -> "Solid":
        return Shape.cast(native.make_cone(radius1, radius2, height, _tuple(pnt), _tuple(dir), angleDegrees))

    @classmethod
    def makeSphere(cls, radius: float, pnt=(0, 0, 0), dir=(0, 0, 1), angleDegrees1: float = 0, angleDegrees2: float = 90, angleDegrees3: float = 360) -> "Solid":
        return Shape.cast(native.make_sphere(radius, _tuple(pnt), _tuple(dir), angleDegrees1, angleDegrees2, angleDegrees3))

    @classmethod
    def makeTorus(cls, radius1: float, radius2: float, pnt=(0, 0, 0), dir=(0, 0, 1), angleDegrees1: float = 0, angleDegrees2: float = 360) -> "Solid":
        return Shape.cast(native.make_torus(radius1, radius2, _tuple(pnt), _tuple(dir)))

    @classmethod
    def makeWedge(cls, dx: float, dy: float, dz: float, xmin: float, zmin: float, xmax: float, zmax: float, pnt=(0, 0, 0), dir=(0, 0, 1)) -> "Solid":
        return Shape.cast(native.make_wedge(dx, dy, dz, xmin, zmin, xmax, zmax, _tuple(pnt), _tuple(dir)))

    @classmethod
    def makeLoft(cls, listOfWire: Sequence[Wire], ruled: bool = False) -> "Solid":
        if len(listOfWire) < 2:
            raise ValueError("More than one wire is required")
        return Shape.cast(native.make_loft([w.wrapped for w in listOfWire], True, ruled))

    @classmethod
    def extrudeLinear(cls, *args, taper: float = 0) -> "Solid":
        # (face, vector[, taper]) or (outerWire, innerWires, vector[, taper]).
        if isinstance(args[0], Face):
            face, vec = args[0], args[1]
            rest = args[2:]
        else:
            outer, inners, vec = args[0], args[1], args[2]
            rest = args[3:]
            face = None
        taper = rest[0] if rest else taper
        if taper:
            # Drafted: the outer wire to its offset (narrower for a positive taper, along
            # the extrusion), the sides ruled between them.
            outer = face.outerWire() if face is not None else outer
            v = Vector(vec)
            tops = outer.offset2D(-v.Length * math.tan(math.radians(taper)), "intersection")
            if len(tops) == 1 and len(tops[0].Edges()) != len(outer.Edges()):
                # (The offset joins straight edges in line: the outline's joined too.)
                outer = _lines_joined(outer)
            if len(tops) != 1 or len(tops[0].Edges()) != len(outer.Edges()):
                # The draft closes some of the outline's arcs: as OCCT's drafted prism.
                if face is None:
                    face = Face.makeFromWires(outer, inners)
                return Shape.cast(native.extrude_drafted(face.wrapped, _tuple(v), math.radians(taper)))
            return cls.makeLoft([outer, tops[0].translate(v)], ruled=True)
        if face is None:
            face = Face.makeFromWires(outer, inners)
        return Shape.cast(native.extrude(face.wrapped, _tuple(vec)))

    @classmethod
    def extrudeLinearWithRotation(cls, *args) -> "Solid":
        """A twisted prism: the profile extruded along `vecNormal` while turning
        `angleDegrees` about it (through `vecCenter`), as a loft through turned copies."""
        if isinstance(args[0], Face):
            outer, inners = args[0].outerWire(), args[0].innerWires()
            center, normal, angle = args[1:4]
        else:
            outer, inners, center, normal, angle = args[:5]
        center, normal = Vector(center), Vector(normal)
        steps = max(2, int(math.ceil(abs(angle) / 15.0)))

        def twisted(w: Wire) -> "Shape":
            sections = [w.rotate(center, center + normal, angle * i / steps).translate(normal * (i / steps)) for i in range(steps + 1)]
            # (Each vertex to its own turned copy: a loft matching starts by nearness would
            # skip round a small profile far from the axis, its copies moved more than its size.)
            return Shape.cast(native.make_loft([s.wrapped for s in sections], True, False, "chordal", "as_given"))

        result = twisted(outer)
        for w in inners:
            result = result.cut(twisted(w))
        return result

    @classmethod
    def revolve(cls, *args) -> "Solid":
        # (face, angle, axisStart, axisEnd) or (outerWire, innerWires, angle, axisStart, axisEnd).
        if isinstance(args[0], Face):
            face, (angleDegrees, axisStart, axisEnd) = args[0], args[1:4]
        else:
            face = Face.makeFromWires(args[0], args[1])
            angleDegrees, axisStart, axisEnd = args[2:5]
        start, end = Vector(axisStart), Vector(axisEnd)
        return Shape.cast(native.revolve(face.wrapped, start.toTuple(), (end - start).toTuple(), angleDegrees))

    @staticmethod
    def _toWire(p) -> Wire:
        return Wire.assembleEdges([p]) if p.ShapeType() == "Edge" else p

    @classmethod
    def sweep(cls, *args, makeSolid: bool = True, isFrenet: bool = False, mode=None, transitionMode: str = "transformed") -> "Shape":
        """The profile (a face, or an outer wire and inner wires) swept along `path`."""
        if isinstance(args[0], Face):
            face, path, rest = args[0], args[1], args[2:]
            outer, inners = face.outerWire(), face.innerWires()
        else:
            outer, inners, path, rest = args[0], list(args[1]), args[2], args[3:]
        options = list(rest) + [None] * 4
        makeSolid = options[0] if rest else makeSolid
        if mode is None and len(rest) > 2:
            mode = rest[2]
        if mode is not None:
            raise NotImplementedError("sweeps guided by an auxiliary spine or binormal")
        spine = cls._toWire(path)
        transition = options[3] if len(rest) > 3 and options[3] else transitionMode
        isFrenet = bool(options[1]) if len(rest) > 1 else isFrenet
        # (Only where the spine has corners: elsewhere every transition is the same sweep.)
        corner = transition if transition in ("right", "round") and _has_corners(spine) else "transformed"
        if makeSolid:
            # (Holes of one face swept at once where the inner wires lie within the outer one;
            # otherwise as CadQuery: each wire swept, the inner ones cut from the outer.)
            box = outer.BoundingBox()
            if all(_box_within(w.BoundingBox(), box) for w in inners):
                profile = Face.makeFromWires(outer, inners)
                return Shape.cast(native.pipe(profile.wrapped, spine.wrapped, corner=corner, frenet=isFrenet))
            solids = [Shape.cast(native.pipe(Face.makeFromWires(w).wrapped, spine.wrapped, corner=corner, frenet=isFrenet)) for w in [outer] + inners]
            return solids[0].cut(*solids[1:])
        shapes = [Shape.cast(native.pipe(w.wrapped, spine.wrapped, corner=corner, frenet=isFrenet)) for w in [outer] + inners]
        return shapes[0] if len(shapes) == 1 else Compound.makeCompound(shapes)

    @classmethod
    def sweep_multi(cls, profiles, path, makeSolid: bool = True, isFrenet: bool = False, mode=None) -> "Shape":
        if mode is not None:
            raise NotImplementedError("sweeps guided by an auxiliary spine or binormal")
        wires = [p if p.ShapeType() == "Wire" else p.outerWire() for p in profiles]
        return Shape.cast(native.pipe_through([w.wrapped for w in wires], cls._toWire(path).wrapped, makeSolid))

    def outerShell(self) -> Shell:
        return self.Shells()[0]

    def innerShells(self) -> List[Shell]:
        return self.Shells()[1:]

    def addCavity(self, *shells) -> "Solid":
        inner = [(s if s.ShapeType() == "Shell" else s.outerShell()).wrapped.reversed() for s in shells]
        return Shape.cast(native.make_solid([s.wrapped for s in self.Shells()] + inner))


class CompSolid(Shape, Mixin3D):
    pass


class Compound(Shape, Mixin3D):
    @classmethod
    def makeCompound(cls, listOfShapes: Iterable[Shape]) -> "Compound":
        return Shape.cast(native.make_compound([s.wrapped for s in listOfShapes]))

    @staticmethod
    def _makeCompound(listOfShapes) -> Any:
        return native.make_compound(list(listOfShapes))

    def remove(self, *shape: Shape) -> "Compound":
        kept = [c for c in self if not any(c.isSame(s) for s in shape)]
        self.wrapped = native.make_compound([c.wrapped for c in kept])
        return self

    @classmethod
    def makeText(
        cls,
        text: str,
        size: float,
        height: float,
        font: str = "Arial",
        fontPath: Optional[str] = None,
        kind: str = "regular",
        halign: str = "center",
        valign: str = "center",
        position: Plane = Plane.XY(),
    ) -> "Shape":
        """Text set in the font (`fontPath`, else the system's font best matching `font` and
        `kind`) at em height `size`, aligned as OCCT aligns it, extruded by `height` (flat
        faces if 0), placed on `position`."""
        from .fonts import font_file

        path = font_file(font, fontPath, kind)
        flat = Shape.cast(native.text(text, float(size), path, halign, valign))
        if height != 0 and flat.Faces():
            normal = flat.Faces()[0].normalAt() * height
            rv = cls(native.extrude(flat.wrapped, normal.toTuple())).transformShape(position.rG)
        else:
            rv = flat.transformShape(position.rG)
        return rv

    def __bool__(self) -> bool:
        return len(native.children(self.wrapped)) > 0

    def fuse(self, *toFuse: Shape, glue: bool = False, tol: Optional[float] = None) -> "Shape":
        args = list(self) + list(toFuse)
        if len(args) <= 1:
            return args[0] if args else self
        # (A compound of solids and its tools: one union, which takes the compound's solids
        # clear of the tools as they are and fuses only those the tools meet, where the
        # pairwise tree would fuse the compound's own pieces with one another again.)
        if toFuse and all(c.wrapped.kind() == "Solid" for c in self) and all(t.wrapped.kind() in ("Solid", "Compound") for t in toFuse):
            tools = toFuse[0] if len(toFuse) == 1 else Shape.fuse(*toFuse)
            return Shape.cast(native.union(self.wrapped, tools.wrapped))
        return Shape.fuse(args[0], *args[1:])

    def ancestors(self, ctx: Shape, kind: str) -> "Compound":
        return Compound.makeCompound([a for s in self for a in s.ancestors(ctx, kind)])

    def siblings(self, ctx: Shape, kind: str, level=1) -> "Compound":
        return Compound.makeCompound(_siblings(list(self), ctx, kind, level))


def _cast_as(cls, handle, forConstruction: bool = False) -> Shape:
    """A shape of class `cls` for `handle` (as `Shape.cast`, its kind known)."""
    s = _object_new(cls)
    s.wrapped, s.label, s.forConstruction = handle, "", forConstruction
    return s


_object_new = object.__new__

_KINDS = {
    "Vertex": Vertex,
    "Edge": Edge,
    "Wire": Wire,
    "Face": Face,
    "Shell": Shell,
    "Solid": Solid,
    "CompSolid": CompSolid,
    "Compound": Compound,
}


def _made_as(name: str, made):
    """The constructor classmethod `made`, returning the class it is called on where that
    is a subclass of what it makes (as CadQuery's, which make `cls(...)`: scripts subclass
    `Wire` and the like to add methods)."""

    def wrapper(cls, *args, **kwargs):
        shape = made.__get__(None, cls)(*args, **kwargs)
        if isinstance(shape, Shape) and cls is not type(shape) and issubclass(cls, type(shape)):
            shape.__class__ = cls
        return shape

    wrapper.__name__, wrapper.__qualname__, wrapper.__doc__ = name, made.__func__.__qualname__, made.__func__.__doc__
    return classmethod(wrapper)


for _cls in _KINDS.values():
    for _name, _attr in list(vars(_cls).items()):
        if isinstance(_attr, classmethod) and (_name.startswith("make") or _name in ("assembleEdges", "combine", "extrudeLinear", "revolve", "sweep", "loft")):
            setattr(_cls, _name, _made_as(_name, _attr))
del _cls, _name, _attr


def _thickened(face: "Face", thickness: float) -> "Solid":
    """The solid between `face` and its offset by `thickness` along its normal (what
    `Face.thicken` makes, kept apart from it: extensions replace that method)."""
    if face.geomType() == "PLANE":
        return Shape.cast(native.extrude(face.wrapped, (face.normalAt() * thickness).toTuple()))
    # The face, its offset, and a ruled wall between each edge of its boundary and that
    # edge's offset (none at a seam, where the face meets itself), sewn.
    moved = Shape.cast(native.offset_shape(face.wrapped, thickness, "intersection"))

    def boundary(f):
        uses = [e for w in f.Wires() for e in w.Edges()]
        return [e for i, e in enumerate(uses) if not any(e.isSame(x) for j, x in enumerate(uses) if j != i)]

    theirs = boundary(moved)
    middles = [x.positionAt(0.5) for x in theirs]
    walls = []
    for e in boundary(face):
        start, end = e.startPoint(), e.endPoint()
        shifted = lambda p: p + face.normalAt(p) * thickness
        a, b = shifted(start), shifted(end)
        target = shifted(e.positionAt(0.5))
        g = theirs[min(range(len(theirs)), key=lambda k: (middles[k] - target).Length)]
        # (Run the same way as `e`.)
        if (g.startPoint() - a).Length + (g.endPoint() - b).Length > (g.endPoint() - a).Length + (g.startPoint() - b).Length:
            g = Shape.cast(g.wrapped.Reversed())
        walls.append(native.make_loft([native.make_wire([e.wrapped]), native.make_wire([g.wrapped])], False, True))
    shell = Shape.cast(native.sew([face.wrapped, moved.wrapped] + walls, 1e-6))
    return Shape.cast(native.make_solid([s.wrapped for s in shell.Shells()]))
