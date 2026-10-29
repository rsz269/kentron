"""CadQuery's `occ_impl.shapes` on pykentron: its shape classes (from `pykentron.shapes`) and
module functions (CadQuery's free-function API, `cadquery.func`).

The functions take and return shapes as CadQuery's do: results of one piece are that piece
(a shell of one open face that face, a compound of one shape that shape), several a
compound. `history` and `name` record an operation's history (CadQuery's `History`)
for extrusions, lofts, sweeps and offsets, found from where the result's faces and edges
lie (kentron's operations do not keep one); elsewhere they are accepted and ignored.
"""

from __future__ import annotations

import math
from typing import Iterable, List, Literal, Optional, Sequence, Tuple, Union

from .. import _kentron as native
from ..geom import Location, Plane, Vector, VectorLike  # noqa: F401
from ..shapes import *  # noqa: F401,F403
from ..shapes import Compound, Edge, Face, Shape, Shell, Solid, Vertex, Wire, _cut_solid_by_solid, _edges_to_wires, _split_pieces, _tuple

Shapes = Literal["Vertex", "Edge", "Wire", "Face", "Shell", "Solid", "CompSolid", "Compound"]
Real = Union[int, float]

TOLERANCE = 1e-6


def sortWiresByBuildOrder(wireList: List[Wire]) -> List[List[Wire]]:
    """The wires grouped into faces: each outer wire, then the wires of its holes (as
    CadQuery's)."""
    if len(wireList) < 2:
        return [wireList]
    faces = Face.makeFromWires(wireList[0], wireList[1:])
    return [[face.outerWire()] + face.innerWires() for face in faces.Faces()]


def wiresToFaces(wireList: List[Wire]) -> List[Face]:
    return Face.makeFromWires(wireList[0], wireList[1:]).Faces()


def edgesToWires(edges: Iterable[Edge], tol: float = 1e-6) -> List[Wire]:
    """The edges joined into wires, end to end."""
    return _edges_to_wires(list(edges), tol)


# --- OCP-level helpers ---------------------------------------------------------------------
# (CadQuery's take and return OCCT shapes, `Shape.wrapped`: here pykentron's handles.)


def downcast(obj):
    """The shape as its own kind: pykentron's handles always are."""
    return obj


def fix(obj):
    """The shape repaired: pykentron's are kept as built."""
    return obj


_KINDS = {"compound": 0, "compsolid": 1, "solid": 2, "shell": 3, "face": 4, "wire": 5, "edge": 6, "vertex": 7}


def shapetype(obj) -> int:
    """The kind of shape, as OCP's TopAbs number."""
    return _KINDS[obj.kind().lower()]


# --- utilities ---------------------------------------------------------------------------


def _types(ts) -> tuple:
    return ts if isinstance(ts, tuple) else (ts,)


def _get(s: Shape, ts) -> Iterable[Shape]:
    """`s` if of one of the types `ts`, or its members if a compound of them; else an
    error."""
    types = _types(ts)
    t = s.ShapeType()
    if t in types:
        yield s
    elif t == "Compound":
        for el in s:
            if el.ShapeType() in types:
                yield el
            else:
                raise ValueError(f"Required type(s): {types}; encountered {el.ShapeType()}")
    else:
        raise ValueError(f"Required type(s): {types}; encountered {t}")


def _get_one(s: Shape, ts) -> Shape:
    types = _types(ts)
    t = s.ShapeType()
    if t in types:
        return s
    if t == "Compound":
        for el in s:
            if el.ShapeType() in types:
                return el
            raise ValueError(f"Required type(s): {types}, encountered {el.ShapeType()}")
    raise ValueError(f"Required type(s): {types}; encountered {t}")


def _get_one_wire(s: Shape) -> Wire:
    rv = _get_one(s, ("Wire", "Edge"))
    return rv if isinstance(rv, Wire) else Wire.assembleEdges([rv])


def _get_wires(s: Shape) -> Iterable[Shape]:
    t = s.ShapeType()
    if t == "Wire":
        yield s
    elif t == "Edge":
        yield Wire.assembleEdges([s])
    elif t == "Compound":
        for el in s:
            yield from _get_wires(el)
    else:
        raise ValueError(f"Required type(s): Edge, Wire; encountered {t}")


def _get_edges(*shapes: Shape) -> Iterable[Shape]:
    for s in shapes:
        t = s.ShapeType()
        if t == "Edge":
            yield s
        elif t == "Wire":
            yield from _get_edges(s.edges())
        elif t == "Compound":
            for el in s:
                yield from _get_edges(el)
        else:
            raise ValueError(f"Required type(s): Edge, Wire; encountered {t}")


def _get_faces(*shapes: Shape) -> Iterable[Face]:
    for s in shapes:
        t = s.ShapeType()
        if t == "Face":
            yield s
        elif t in ("Edge", "Wire"):
            yield face(s)
        elif t == "Compound":
            yield from _get_faces(*s)
        else:
            raise ValueError(f"Required type(s): Edge, Wire, Face; encountered {t}")


def _get_wire_lists(s: Sequence[Shape]) -> List[list]:
    """Matching wires (or vertices at the ends) of each section, for sweeps and lofts."""
    wire_lists: List[list] = []
    last = len(s) - 1
    for i, el in enumerate(s):
        if i == 0:
            try:
                wire_lists = [[w] for w in _get_wires(el)]
            except ValueError:
                wire_lists = [[v] for v in el.Vertices()]
            if not wire_lists:
                break
        elif i == last:
            try:
                for wl, w in zip(wire_lists, _get_wires(el)):
                    wl.append(w)
            except ValueError:
                for wl, v in zip(wire_lists, el.Vertices()):
                    wl.append(v)
        else:
            for wl, w in zip(wire_lists, _get_wires(el)):
                wl.append(w)
    return wire_lists


def _get_face_lists(s: Sequence[Shape]) -> List[list]:
    face_lists: List[list] = []
    last = len(s) - 1
    for i, el in enumerate(s):
        if i == 0:
            face_lists = [[f] for f in el.Faces()]
            if not face_lists and not el.Edges():
                face_lists = [[v] for v in el.Vertices()]
            if not face_lists:
                break
        elif i == last:
            faces = el.Faces()
            if len(faces) == len(face_lists):
                for fl, f in zip(face_lists, faces):
                    fl.append(f)
            else:
                for fl, v in zip(face_lists, el.Vertices()):
                    fl.append(v)
        else:
            for fl, f in zip(face_lists, el.Faces()):
                fl.append(f)
    if any(isinstance(el[0], Vertex) and isinstance(el[1], Vertex) for el in face_lists):
        return []
    return face_lists


def _normalize(s: Shape) -> Shape:
    """A shell of one open face as the face, a compound of one shape as the shape."""
    t = s.ShapeType()
    if t == "Shell":
        faces = s.Faces()
        if len(faces) == 1 and not _closed_shell(s):
            return faces[0]
    elif t == "Compound":
        objs = list(s)
        if len(objs) == 1:
            return objs[0]
    return s


def _closed_shell(s: Shape) -> bool:
    # Every edge shared by two faces.
    counts = {}
    for f in s.Faces():
        for e in f.Edges():
            counts[e] = counts.get(e, 0) + 1
    return bool(counts) and all(c >= 2 for c in counts.values())


def _compound_or_shape(s) -> Shape:
    """Shapes (or one) as one shape, normalized, or a compound of them."""
    if isinstance(s, Shape):
        return _normalize(s)
    s = list(s)
    if len(s) == 1:
        return _normalize(s[0] if isinstance(s[0], Shape) else Shape.cast(s[0]))
    return Compound.makeCompound([_normalize(el if isinstance(el, Shape) else Shape.cast(el)) for el in s])


def _shape_to_faces_shells(s: Shape) -> Tuple[List[Face], List[Shell]]:
    faces, shells = [], []
    t = s.ShapeType()
    members = list(s) if t == "Compound" else [s]
    for el in members:
        kind = el.ShapeType()
        if kind == "Face":
            faces.append(el)
        elif kind == "Shell":
            shells.append(el)
        else:
            raise ValueError(f"Found unexpected type {kind}")
    return faces, shells


def _sequence(args) -> bool:
    return len(args) == 1 and isinstance(args[0], (list, tuple)) or (len(args) == 1 and not isinstance(args[0], Shape) and hasattr(args[0], "__iter__"))


# --- construction ------------------------------------------------------------------------


def wire(*s) -> Wire:
    """The wire through edges (and wires' edges)."""
    if _sequence(s):
        s = tuple(s[0])
    return Wire.assembleEdges([e for el in s for e in _get_edges(el)])


def face(*s) -> Face:
    """The planar face bounded by wires (or edges), the inner ones its holes."""
    if _sequence(s):
        s = tuple(s[0])
    wires = [w for el in s for w in _get_wires(el)]
    if not wires:
        raise ValueError("Face construction failed")
    # The largest wire bounds the face; those within it are its holes.
    wires.sort(key=lambda w: -w.BoundingBox().DiagonalLength)
    rv = Face.makeFromWires(wires[0], wires[1:])
    return _get_one(rv, "Face")


def shell(*s, tol: float = 1e-6, manifold: bool = True, ctx=None, history=None, name=None) -> Shape:
    """The faces sewn into shells (a face of its own a shell of one face); with `ctx`, the
    shapes the faces come from (CadQuery's local sewing), sewn with them (each face once).
    Each face's image in the history is the face the sewing made of it."""
    if _sequence(s):
        s = tuple(s[0])
    faces = [f for el in s for f in _get(el, "Face")]
    context = [] if not ctx else (ctx.Faces() if isinstance(ctx, Shape) else [f for c in ctx for f in c.Faces()])
    faces = faces + [f for f in context if not any(f.isSame(g) for g in faces)]
    if history is not None:
        native.record_begin()
    sewed = None
    try:
        sewed = Shape.cast(native.sew([f.wrapped for f in faces], tol))
    finally:
        pairs = native.record_end(sewed.wrapped) if history is not None and sewed is not None else []
        if history is not None and sewed is None:
            native.record_end(faces[0].wrapped)
    rv_faces, rv_shells = _shape_to_faces_shells(sewed)
    shells = rv_shells + [Shape.cast(native.make_shell([f.wrapped])) for f in rv_faces]
    result = _compound_or_shape(shells)
    if not manifold and len(shells) > 1:
        # (Non-manifold sewing, as OCCT's: all the faces one shell, where more than two may
        # meet at an edge. kentron's sewing joins two faces an edge at most: the shells it
        # made, their faces put into one.)
        result = Shape.cast(native.make_shell([f.wrapped for sh in shells for f in sh.Faces()]))
    if history is not None:
        made = {}
        for source, image in pairs:
            made.setdefault(_node(source), []).append(_node(image))
        placed = {_node(x.wrapped): x for x in result.Faces()}
        op = Op(name)
        history.append(op, name)
        for f in faces:
            op._tracked.add(f)
            images = [placed[n] for n in made.get(_node(f.wrapped), []) if n in placed]
            if images:
                op._images[f] = _normalize(compound(images))
                op._modified[f] = _normalize(compound(images))
    return result


def solid(*s, inner=None, tol: float = 1e-6, history=None, name=None) -> Shape:
    """The solid bounded by shells (or by the shell sewn from faces); given a sequence of
    faces, and `inner` ones, the solid with cavities."""
    if len(s) >= 1 and isinstance(s[0], (list, tuple)):
        if len(s) > 1 and inner is None:
            inner = s[1]
        outer = _get_one(shell(*s[0], tol=tol), "Shell")
        shells = [outer.wrapped]
        if inner:
            shells += [sh.wrapped.reversed() for sh in _get(shell(*inner, tol=tol), "Shell")]
        return Shape.cast(native.make_solid(shells))
    found = [f for el in s for f in _get(el, ("Shell", "Face"))]
    shells = [el for el in found if isinstance(el, Shell)]
    if not shells:
        shells = [_get_one(shell(*[el for el in found if isinstance(el, Face)], tol=tol), "Shell")]
    return _compound_or_shape([Shape.cast(native.make_solid([sh.wrapped])) for sh in shells])


def compound(*s) -> Compound:
    if _sequence(s):
        s = tuple(s[0])
    return Compound.makeCompound(s)


def vertex(*args) -> Vertex:
    p = Vector(*args) if len(args) == 3 else Vector(args[0])
    return Shape.cast(native.make_vertex(p.toTuple()))


def segment(p1, p2) -> Edge:
    return Edge.makeLine(p1, p2)


def polyline(*pts) -> Wire:
    return Wire.makePolygon(pts)


def polygon(*pts) -> Wire:
    return Wire.makePolygon(pts, close=True)


def rect(w: float, h: float) -> Wire:
    return polygon((-w / 2, -h / 2, 0), (w / 2, -h / 2, 0), (w / 2, h / 2, 0), (-w / 2, h / 2, 0))


def spline(*args, tgts=None, params=None, tol: float = 1e-6, periodic: bool = False, scale: bool = True) -> Edge:
    """The spline through points (given in turn, or as a sequence, then with tangents and
    parameters), as OCCT interpolates it."""
    if len(args) >= 1 and isinstance(args[0], (list, tuple)) and args[0] and isinstance(args[0][0], (list, tuple, Vector)):
        pts = args[0]
        if len(args) > 1:
            tgts = args[1]
        if len(args) > 2:
            params = args[2]
    else:
        pts = args
    return Edge.makeSpline([Vector(p) for p in pts], tangents=[Vector(t) for t in tgts] if tgts is not None else None, periodic=periodic, parameters=params, scale=scale, tol=tol)


def circle(r: float) -> Edge:
    return Edge.makeCircle(r)


def ellipse(r1: float, r2: float) -> Edge:
    return Edge.makeEllipse(r1, r2)


def plane(*args) -> Face:
    """The plane face `w` × `l` about the origin (or "infinite": as large as kentron takes)."""
    if args:
        w, l = args
        return Shape.cast(native.make_plane_face((0, 0, 0), (0, 0, 1), -w / 2, w / 2, -l / 2, l / 2))
    return Face.makePlane()


def box(w: float, l: float, h: float) -> Solid:
    return Solid.makeBox(w, l, h, (-w / 2, -l / 2, 0))


def cylinder(d: float, h: float) -> Solid:
    return Solid.makeCylinder(d / 2, h)


def sphere(d: float) -> Solid:
    return Solid.makeSphere(d / 2, angleDegrees1=-90, angleDegrees2=90)


def torus(d1: float, d2: float) -> Solid:
    return Solid.makeTorus(d1 / 2, d2 / 2)


def cone(*args) -> Solid:
    d1, d2, h = (args[0], 0.0, args[1]) if len(args) == 2 else args
    return Solid.makeCone(d1 / 2, d2 / 2, h)


def text(txt: str, size: Real, *args, **kwargs) -> Shape:
    """Flat text (`text(txt, size, font, path, kind, halign, valign)`), text along a spine
    (`text(txt, size, spine, planar=False, ...)`: each glyph stood on the spine at its
    position along it, upright or lying in the plane with `planar`), or text on a base face
    (`text(txt, size, spine, base, ...)`: the glyphs along the spine projected onto it), as
    CadQuery's."""
    spine = kwargs.pop("spine", None)
    base = kwargs.pop("base", None)
    args = list(args)
    if spine is None and args and isinstance(args[0], Shape):
        spine = args.pop(0)
    if spine is None:
        return _flat_text(txt, size, *args, **kwargs)
    if base is None and args and isinstance(args[0], Shape):
        base = args.pop(0)
    if base is not None:
        face = _get_one(base, "Face")
        glyphs = text(txt, size, spine, False, *args, **kwargs)
        return _normalize(compound([f.project(face, f.normalAt()) for f in glyphs.Faces()]))
    planar = kwargs.pop("planar", args.pop(0) if args and isinstance(args[0], bool) else False)
    wire = _get_one_wire(spine)
    length = wire.Length()
    placed = []
    for f in _flat_text(txt, size, *args, **kwargs).Faces():
        x = f.BoundingBox().center.x
        placed.append(f.moved(-x).moved(rx=-90 if planar else 0, ry=-90).moved(wire.locationAt(x / length)))
    return _normalize(compound(placed))


def _flat_text(
    txt: str,
    size: Real,
    font: str = "Arial",
    path: Optional[str] = None,
    kind: str = "regular",
    halign: str = "center",
    valign: str = "center",
) -> Shape:
    """The glyphs' faces, fused where they overlap."""
    from ..fonts import font_file

    flat = Shape.cast(native.text(txt, float(size), font_file(font, path, kind), halign, valign))
    return clean(compound(flat.Faces()).fuse())


# --- operations --------------------------------------------------------------------------


def setThreads(n: int) -> None:
    """(kentron's operations use their own thread pool.)"""


def fuse(s1: Shape, s2: Shape, *shapes: Shape, tol: float = 0.0, glue=None, history=None, name=None) -> Shape:
    def run():
        result = s1
        for s in (s2,) + shapes:
            result = Shape.cast(native.union(result.wrapped, s.wrapped))
        return _compound_or_shape(result)
    return _recorded(history, name, [s1, s2, *shapes], run)


def cut(s1: Shape, s2: Shape, tol: float = 0.0, glue=None, history=None, name=None) -> Shape:
    return _recorded(history, name, [s1, s2], lambda: _compound_or_shape(Shape.cast(_cut_solid_by_solid(s1.wrapped, s2.wrapped))))


def intersect(s1: Shape, s2: Shape, tol: float = 0.0, glue=None, history=None, name=None) -> Shape:
    return _recorded(history, name, [s1, s2], lambda: _compound_or_shape(Shape.cast(native.intersection(s1.wrapped, s2.wrapped))))


def split(s1: Shape, s2: Shape, tol: float = 0.0, history=None, name=None) -> Shape:
    """`s1` cut into pieces by `s2`."""
    return _compound_or_shape(_split_pieces(s1, [s2]))


def imprint(*shapes: Shape, tol: float = 0.0, glue="full", history=None, name=None) -> Shape:
    """The shapes split by one another (OCCT's general fuse): all the pieces."""
    return _recorded(history, name, list(shapes), lambda: _compound_or_shape([Shape.cast(h) for _, h in native.cells([s.wrapped for s in shapes])]), images=True)


def clean(s: Shape) -> Shape:
    return _compound_or_shape(s.clean())


def fill(s: Shape, constraints=()) -> Shape:
    """The face spanning `s`'s edges (a closed loop), through `constraints`: points, and
    curves (their edges, which the face passes along, sampled as points it passes through),
    as OCCT's filling."""
    points = []
    for c in constraints:
        if isinstance(c, Shape):
            for e in _get_edges(c):
                points.extend(_tuple(e.positionAt(k / 16)) for k in range(17))
        else:
            points.append(_tuple(Vector(c)))
    edges = [e.wrapped for e in _get_edges(s)]
    return _compound_or_shape(Shape.cast(native.make_filling(edges, points, 3, 15)))


def cap(s: Shape, ctx: Shape, constraints=()) -> Shape:
    raise NotImplementedError("filling surfaces")


def fillet(s: Shape, edges: Shape, r: float, history=None, name=None) -> Shape:
    target = _get_one(s, ("Shell", "Solid"))
    run = lambda: _compound_or_shape(Shape.cast(native.fillet(target.wrapped, [e.wrapped for e in _get_edges(edges.edges())], r)))
    return _recorded(history, name, [target, edges], run)


def chamfer(s: Shape, edges: Shape, d: float, history=None, name=None) -> Shape:
    """`s`'s (a shell's or solid's) edges bevelled by `d` along both faces."""
    target = _get_one(s, ("Shell", "Solid"))
    spec = [(e.wrapped, d, d, None) for e in _get_edges(edges.edges())]
    return _recorded(history, name, [target, edges], lambda: _compound_or_shape(Shape.cast(native.chamfer(target.wrapped, spec))))


def extrude(s: Shape, d, both: bool = False, history=None, name=None) -> Shape:
    """`s` (vertices, edges, wires or faces) swept along `d` (both ways if `both`)."""
    d = Vector(d)
    results, generated, ends = [], {}, []
    for el in _get(s, ("Vertex", "Edge", "Wire", "Face")):
        profile, by = (el.moved(-d), d * 2) if both else (el, d)
        if history is None:
            results.append(Shape.cast(native.extrude(profile.wrapped, by.toTuple())))
            continue
        # (With its history as kentron's extrusion keeps it.)
        shape, swept, firsts_lasts, _ = native.extrude_recorded(profile.wrapped, by.toTuple())
        results.append(Shape.cast(shape))
        for source, made in swept:
            generated.setdefault(_node(source), []).extend(Shape.cast(m) for m in made)
        ends.extend((_node(source), Shape.cast(a), Shape.cast(b)) for source, a, b in firsts_lasts)
    result = _compound_or_shape(results)
    if history is not None:
        _swept_op(history, name, [s], s, generated, ends)
    return result


def _swept_op(history, name, inputs, profile: Shape, generated, ends) -> None:
    """A step of `history` from a sweep's own (kentron's) history: what each profile edge
    and vertex swept (`generated`, by node), and each profile face's first and last faces
    (`ends`, by node)."""
    op = Op(name)
    history.append(op, name)
    for shape in inputs:
        op._tracked.update(shape.Faces())
        op._tracked.update(shape.Edges())
        op._tracked.update(shape.Vertices())
    for sub in list(profile.Edges()) + list(profile.Vertices()):
        if _node(sub.wrapped) in generated:
            op._generated[sub] = _normalize(compound(generated[_node(sub.wrapped)]))
    for f in profile.Faces():
        for key, first, last in ends:
            if key == _node(f.wrapped):
                op._first[f], op._last[f] = first, last
                op._first_shape = op._first_shape | first
                op._last_shape = op._last_shape | last


def revolve(s: Shape, p, d, a: float = 360, history=None, name=None) -> Shape:
    results = [Shape.cast(native.revolve(el.wrapped, _tuple(p), _tuple(d), a)) for el in _get(s, ("Vertex", "Edge", "Wire", "Face"))]
    return _compound_or_shape(results)


def offset(s: Shape, t: float, cap: bool = True, both: bool = False, tol: float = 1e-6, history=None, name=None) -> Shape:
    """Faces or shells thickened by `t` along their normals: a closed shell into the solid
    between it and its offset (joined where they meet, as CadQuery's offsets are), the faces
    of others one by one."""
    def one(t):
        out = []
        for el in _get(s, ("Face", "Shell")):
            if el.ShapeType() == "Shell" and el.Closed():
                solid = Shape.cast(native.make_solid([el.wrapped]))
                out.append(Shape.cast(native.thick_solid(solid.wrapped, [], t, "intersection")))
                continue
            for f in el.Faces():
                out.append(f.thicken(t))
        return out

    if history is not None and all(f.geomType() == "PLANE" for f in s.Faces()) and not any(el.ShapeType() == "Shell" and el.Closed() for el in _get(s, ("Face", "Shell"))):
        return _offset_planes_recorded(s, t, both, history, name)
    if both:
        pos, neg = one(t), one(-t)
        both_ = [a.fuse(b) for a, b in zip(pos, neg)]
        result = both_[0] if len(both_) == 1 else Compound.makeCompound(both_)
    else:
        result = _compound_or_shape(one(t))
    if history is not None:
        offsets = (1, -1) if both else (1,)
        _record(history, name, [s], result, moved=[(f, [f.normalAt() * (k * t) for k in offsets]) for f in s.Faces()])
    return result


def _offset_planes_recorded(s: Shape, t: float, both: bool, history, name) -> Shape:
    """`offset` of planar faces, with its history as kentron's extrusions keep it: each face
    extruded along its normal by `t` (and by `-t`, fused, with `both`), the face's offsets
    the extrusions' end faces, its edges' the copies of its edges there."""
    results, offsets = [], {}
    for f in s.Faces():
        made = {}
        parts = []
        for k in ((1, -1) if both else (1,)):
            shape, _, ends, copies = native.extrude_recorded(f.wrapped, (f.normalAt() * (k * t)).toTuple())
            parts.append(Shape.cast(shape))
            for source, _, last in ends:
                made.setdefault(_node(source), []).append(last)
            for source, copy in copies:
                made.setdefault(_node(source), []).append(copy)
        if len(parts) == 2:
            native.record_begin()
            fused = None
            try:
                fused = parts[0].fuse(parts[1])
            finally:
                pairs = native.record_end(fused.wrapped if fused is not None else parts[0].wrapped)
            made = _through(pairs, made)
            parts = [fused]
        results.append(parts[0])
        for k, v in made.items():
            offsets.setdefault(k, []).extend(v)
    result = results[0] if len(results) == 1 else Compound.makeCompound(results)
    placed = {_node(x.wrapped): x for x in list(result.Faces()) + list(result.Edges())}
    op = Op(name)
    history.append(op, name)
    op._tracked.update(s.Faces())
    op._tracked.update(s.Edges())
    op._tracked.update(s.Vertices())
    for sub in list(s.Faces()) + list(s.Edges()):
        images = [placed[_node(h)] for h in offsets.get(_node(sub.wrapped), []) if _node(h) in placed]
        if images:
            op._generated[sub] = _normalize(compound(images))
    return result


def offset2D(s: Shape, t: float, ctx=None, kind: str = "arc", closed: bool = True) -> Shape:
    """Edges or wires offset in their plane; open ones (not `closed`) on one side only."""
    results = []
    for w in _get_wires(s):
        results.extend(w.offset2D(t, kind, not closed))
    return _compound_or_shape(results)


def chamfer2D(s: Shape, verts: Shape, d: float) -> Shape:
    return _compound_or_shape(_get_one(s, "Face").chamfer2D(d, verts.vertices() if verts.ShapeType() != "Vertex" else [verts]))


def fillet2D(s: Shape, verts: Shape, r: float) -> Shape:
    vs = list(verts.Vertices())
    return _compound_or_shape(_get_one(s, "Face").fillet2D(r, vs))


def sweep(s, path: Shape, aux=None, cap: bool = False, transition: str = "transformed", history=None, name=None) -> Shape:
    """Edges, wires or faces swept along `path` (several sections, of matching wires or
    faces, blended along it)."""
    if history is not None and not isinstance(s, (list, tuple)) and aux is None and transition == "transformed" and s.Faces():
        # (One profile's faces: with the history kentron's pipe keeps.)
        spine = _get_one_wire(path)
        results, generated, ends = [], {}, []
        for f in s.Faces():
            shape, swept, firsts_lasts = native.pipe_recorded(f.wrapped, spine.wrapped)
            results.append(Shape.cast(shape))
            for source, made in swept:
                generated.setdefault(_node(source), []).extend(Shape.cast(m) for m in made)
            ends.extend((_node(source), Shape.cast(a), Shape.cast(b)) for source, a, b in firsts_lasts)
        result = _compound_or_shape(results)
        _swept_op(history, name, [s, path], s, generated, ends)
        return result
    result = _sweep(s, path, aux, cap, transition)
    if history is not None:
        # (Sections blended along the spine are made alike, split to match, before they are
        # swept: their history is found from where the result's faces lie.)
        if isinstance(s, (list, tuple)):
            first, last = s[0], s[-1]
            _record(history, name, [*s, path], result, firsts=[(f, f) for f in first.Faces()], lasts=[(f, f) for f in last.Faces()], sections=first.Edges())
        else:
            _record(history, name, [s, path], result, firsts=[(f, f) for f in s.Faces()], sections=s.Edges(), last_rest=s)
    return result


def _sweep(s, path: Shape, aux, cap: bool, transition: str) -> Shape:
    if aux is not None:
        raise NotImplementedError("sweeps guided by an auxiliary spine")
    spine = _get_one_wire(path)
    if isinstance(s, (list, tuple)):
        lists = [[f for f in el.Faces()] for el in s]
        if all(lists) and all(len(l) == len(lists[0]) for l in lists):
            results = []
            for k in range(len(lists[0])):
                sections = [l[k] for l in lists]
                outer = Shape.cast(native.pipe_through([f.outerWire().wrapped for f in sections], spine.wrapped, True))
                for j in range(len(sections[0].innerWires())):
                    hole = Shape.cast(native.pipe_through([f.innerWires()[j].wrapped for f in sections], spine.wrapped, True))
                    outer = outer.cut(hole)
                results.append(outer)
            return _compound_or_shape(results)
        wire_lists = [list(_get_wires(el)) for el in s]
        results = [Shape.cast(native.pipe_through([wl[k].wrapped for wl in wire_lists], spine.wrapped, cap)) for k in range(len(wire_lists[0]))]
        return _compound_or_shape(results)
    faces = s.Faces()
    if faces:
        return _compound_or_shape([Shape.cast(native.pipe(f.wrapped, spine.wrapped)) for f in faces])
    results = []
    for w in _get_wires(s):
        profile = Face.makeFromWires(w) if cap else w
        results.append(Shape.cast(native.pipe(profile.wrapped, spine.wrapped)))
    return _compound_or_shape(results)


def loft(*s, cap: bool = False, ruled: bool = False, continuity: str = "C2", parametrization: str = "uniform", degree: int = 3, compat: bool = True, smoothing: bool = False, weights=(1, 1, 1), history=None, name=None) -> Shape:
    """A loft through sections: faces (solids, their holes lofted too), or wires and edges
    (vertices at the ends), capped into solids if `cap`; passing through them spaced as
    `parametrization` says (by default evenly, as CadQuery's: `Solid.makeLoft` spaces them
    by their distances apart, as OCCT does by default)."""
    if len(s) == 1 and isinstance(s[0], (list, tuple)):
        s = tuple(s[0])
    if history is not None:
        return _loft_recorded(s, cap, ruled, parametrization, history, name)
    results = []
    face_lists = _get_face_lists(s)
    if face_lists:
        for el in face_lists:
            outer = [f.outerWire() if isinstance(f, Face) else f for f in el]
            result = Shape.cast(native.make_loft([w.wrapped for w in outer], True, ruled, parametrization))
            if not any(isinstance(f, Vertex) for f in el):
                for j in range(len(el[0].innerWires())):
                    hole = Shape.cast(native.make_loft([f.innerWires()[j].wrapped for f in el], True, ruled, parametrization))
                    result = result.cut(hole)
            results.append(result)
    else:
        for el in _get_wire_lists(s):
            results.append(Shape.cast(native.make_loft([w.wrapped for w in el], cap, ruled, parametrization)))
    return _compound_or_shape(results)


def _through(pairs, current):
    """`current` (node to handles) carried through a recorded boolean's `pairs`: each handle
    replaced by the faces and edges the boolean made from it (none where it is gone)."""
    made = {}
    for source, image in pairs:
        made.setdefault(_node(source), []).append(image)
    return {k: [m for h in v for m in made.get(_node(h), [])] for k, v in current.items()}


def _cut_recorded(shape: Shape, tool: Shape, current):
    """`shape` with `tool` cut away, and `current` carried through the cut."""
    native.record_begin()
    result = None
    try:
        result = shape.cut(tool)
    finally:
        pairs = native.record_end(result.wrapped if result is not None else shape.wrapped)
    return result, _through(pairs, current)


def _loft_recorded(s, cap: bool, ruled: bool, parametrization: str, history, name) -> Shape:
    """`loft` with its history as kentron's lofts keep it: what each section edge and vertex
    lofted and the caps on the first and last sections, carried through the cuts making
    the holes of faces lofted with holes."""
    results, current = [], {}
    face_lists = _get_face_lists(s)
    lists = face_lists if face_lists else _get_wire_lists(s)
    for el in lists:
        if face_lists:
            outer = [f.outerWire() if isinstance(f, Face) else f for f in el]
            shape, made, (first, last) = native.make_loft_recorded([w.wrapped for w in outer], True, ruled, parametrization)
        else:
            shape, made, (first, last) = native.make_loft_recorded([w.wrapped for w in el], cap, ruled, parametrization)
        here = {_node(src): list(images) for src, images in made}
        here[("first",)] = [first] if first is not None else []
        here[("last",)] = [last] if last is not None else []
        result = Shape.cast(shape)
        if face_lists and not any(isinstance(f, Vertex) for f in el):
            for j in range(len(el[0].innerWires())):
                hole, hole_made, _ = native.make_loft_recorded([f.innerWires()[j].wrapped for f in el], True, ruled, parametrization)
                for src, images in hole_made:
                    here.setdefault(_node(src), []).extend(images)
                result, here = _cut_recorded(result, Shape.cast(hole), here)
        results.append(result)
        for k, v in here.items():
            current.setdefault(k, []).extend(v)
    result = _compound_or_shape(results)
    placed = {_node(x.wrapped): x for x in list(result.Faces()) + list(result.Edges())}
    def found(handles):
        return [placed[_node(h)] for h in handles if _node(h) in placed]
    op = Op(name)
    history.append(op, name)
    for section in s:
        op._tracked.update(section.Faces())
        op._tracked.update(section.Edges())
        op._tracked.update(section.Vertices())
        for sub in list(section.Edges()) + list(section.Vertices()):
            images = found(current.get(_node(sub.wrapped), []))
            if images:
                op._generated[sub] = _normalize(compound(images))
    for key, shapes, per, section in ((("first",), "_first_shape", op._first, s[0]), (("last",), "_last_shape", op._last, s[-1])):
        caps = found(current.get(key, []))
        if caps:
            setattr(op, shapes, _normalize(compound(caps)))
            for f in section.Faces():
                per[f] = _normalize(compound(caps))
    return result


def project(s: Shape, base: Shape, *args, **kwargs) -> Shape:
    raise NotImplementedError("projecting curves onto shapes")


def hollow(s: Shape, *args, tol: float = 1e-3, kind: str = "intersection", history=None, name=None) -> Shape:
    """`s` made a shell of thickness `t`, open where `faces` are (closed, a cavity
    inside, without)."""
    if args and (args[0] is None or isinstance(args[0], Shape)):
        faces, t = args[0], args[1]
    else:
        faces, t = None, args[0]
    target = s.solid()
    return _compound_or_shape(target.hollow(faces.Faces() if faces else None, t, tol, kind))


def prism(ctx: Shape, base, faces: Shape, t, *args, history=None, name=None, **kwargs) -> Shape:
    """A prismatic feature, as CadQuery's (OCCT's BRepFeat_MakePrism and MakeDPrism): each of
    `faces` swept, fused with `ctx` (or, not `additive`, cut from it), as far as `t` (a
    length), up to its face (a shape), from and up to a pair of faces, or (None) through
    all; along a direction, `prism(ctx, base, faces, t, dir, additive)`, or along each
    face's normal, tapered by an angle in degrees, `prism(ctx, base, faces, t, angle,
    additive)`. (The base face, OCCT's sketch face, is not needed.)"""
    if history is not None:
        raise NotImplementedError("histories of operations")
    rest = list(args)
    direction, angle = None, 0.0
    if rest and not isinstance(rest[0], bool) and isinstance(rest[0], (int, float)):
        angle = float(rest.pop(0))
    elif rest and not isinstance(rest[0], bool):
        direction = Vector(rest.pop(0))
    if "angle" in kwargs:
        angle = float(kwargs.pop("angle"))
    if "dir" in kwargs:
        direction = Vector(kwargs.pop("dir"))
    additive = bool(rest.pop(0)) if rest else bool(kwargs.pop("additive", True))
    extent = {}
    if isinstance(t, Shape):
        extent["until"] = t.face().wrapped
    elif isinstance(t, tuple):
        extent["start"], extent["until"] = t[0].face().wrapped, t[1].face().wrapped
    elif t is not None:
        extent["length"] = float(t)
    result = ctx
    for f in _get_faces(faces):
        along = direction if direction is not None else f.normalAt()
        result = Shape.cast(native.feature_prism(result.wrapped, f.wrapped, _tuple(along), math.radians(angle), additive, **extent))
    return _compound_or_shape(result)


def draft(ctx: Shape, base: Shape, faces: Shape, *args, history=None, name=None) -> Shape:
    """`ctx` with `faces` drafted by an angle (degrees) about where they meet `base`'s
    plane, against its normal or a given direction: `draft(ctx, base, faces, angle)` or
    `draft(ctx, base, faces, dir, angle)`, as CadQuery's (OCCT's draft angle)."""
    if history is not None:
        raise NotImplementedError("histories of operations")
    if len(args) == 1:
        direction, angle = None, args[0]
    elif len(args) == 2:
        direction, angle = args
    else:
        raise TypeError("draft(ctx, base, faces, [dir,] angle)")
    base_face = base.face() if hasattr(base, "face") else base
    normal = base_face.normalAt()
    pull = Vector(direction) if direction is not None else normal
    drafted = list(_get_faces(faces))
    # (OCCT drafts planes, cylinders and cones: other faces cannot be drafted.)
    for f in drafted:
        if f.geomType() not in ("PLANE", "CYLINDER", "CONE"):
            raise ValueError(f"Face {f} cannot be used in a draft operation.")
    try:
        result = native.draft_faces(ctx.wrapped, [f.wrapped for f in drafted], _tuple(pull), math.radians(angle), _tuple(base_face.Center()), _tuple(normal))
    except ValueError as e:
        raise ValueError(f"the faces cannot be drafted: {e}") from e
    return _compound_or_shape(Shape.cast(result))


def check(s: Shape, results=None, tol: Optional[float] = None, verbose: bool = True) -> bool:
    """True if `s` is valid; its problems added to `results` (if given)."""
    problems = native.problems(s.wrapped)
    if results is not None:
        results.clear()
        results.extend(([s], p) for p in problems)
    return not problems


# --- operation history ----------------------------------------------------------------


class Op:
    """One step of an operation history (CadQuery's): the shapes it tracked, deleted,
    modified and generated, and its first and last shapes (bottom and top of a sweep)."""

    def __init__(self, name: Optional[str] = None) -> None:
        self._name = name if name else None
        self._tracked = set()
        self._deleted = []
        self._modified = {}
        self._generated = {}
        self._images = {}
        self._first = {}
        self._last = {}
        self._first_shape = compound()
        self._last_shape = compound()

    def _get(self, d, k: Shape) -> Shape:
        if k.ShapeType() == "Compound":
            tmp = []
            for el in k:
                val = d[el]
                if val.ShapeType() == "Compound":
                    tmp.extend(val)
                else:
                    tmp.append(val)
            return _normalize(compound(tmp))
        return _normalize(d[k])

    def modified(self, s: Optional[Shape] = None) -> Shape:
        """Shapes modified from s."""
        return self._get(self._modified, s) if s else _normalize(compound(*self._modified.values()))

    def generated(self, s: Optional[Shape] = None) -> Shape:
        """Shapes generated from s."""
        return self._get(self._generated, s) if s else _normalize(compound(*self._generated.values()))

    def deleted(self) -> Shape:
        """Deleted shapes."""
        return _normalize(compound(self._deleted))

    def images(self, s: Shape) -> Shape:
        """Images of s."""
        return self._get(self._images, s)

    def first(self, s: Optional[Shape] = None) -> Shape:
        """First shape (e.g. bottom face) or first shape generated from s."""
        return self._get(self._first, s) if s else _normalize(self._first_shape)

    def last(self, s: Optional[Shape] = None) -> Shape:
        """Last shape (e.g. top face) or last shape generated from s."""
        return self._get(self._last, s) if s else _normalize(self._last_shape)


def _combine_hist_dict(d1, *ds) -> None:
    """History dicts combined into `d1`: a key in both maps to the compound of both values."""
    for d in ds:
        for k in d1.keys() & d.keys():
            d1[k] |= d[k]
        for k in d.keys() - d1.keys():
            d1[k] = d[k]


class History:
    """An operation history: its steps in order, and by name."""

    def __init__(self) -> None:
        self.ops = []
        self.opDict = dict()

    def __getitem__(self, ix):
        return self.opDict[ix] if isinstance(ix, str) else self.ops[ix]

    def pop(self) -> Op:
        return self.ops.pop()

    def append(self, op: Op, name: Optional[str] = None) -> None:
        self.ops.append(op)
        if name:
            self.opDict[name] = op


def _alike(a: Shape, b: Shape, tol: float = 1e-6) -> bool:
    """Whether two faces or edges lie in one place (the same kind, size and centre)."""
    if a.ShapeType() != b.ShapeType():
        return False
    measure = (lambda s: s.Area()) if a.ShapeType() == "Face" else (lambda s: s.Length())
    size = max(measure(a), measure(b))
    reach = tol * (1.0 + a.BoundingBox().DiagonalLength)
    return abs(measure(a) - measure(b)) <= tol * (1.0 + size) and (a.Center() - b.Center()).Length <= reach


def _record(history, name, inputs, result: Shape, firsts=(), lasts=(), sections=(), moved=(), last_rest=None) -> None:
    """A step of `history` for `result`, made from `inputs` (all their faces, edges and
    vertices tracked), found from where the result's faces and edges lie: its first and last
    faces, those alike `firsts` and `lasts` (each also recorded as the first and last of
    the input it came from); each edge of `sections` (the profile's) generating the result's
    faces running from an edge alike it, other than the first and last; and for each
    `(input, offsets)` of `moved`, the input face (and its edges) generating the result's
    faces (edges) alike it moved by each offset; with `last_rest` (the profile), the faces
    left over (neither first nor generated) its last (a sweep's end, moved along the path)."""
    if history is None:
        return
    op = Op(name)
    history.append(op, name)
    for shape in inputs:
        op._tracked.update(shape.Faces())
        op._tracked.update(shape.Edges())
        op._tracked.update(shape.Vertices())
    faces = result.Faces()
    def found(target):
        return [f for f in faces if _alike(f, target)]
    ends = []
    for which, targets, shapes, per in (("first", firsts, "_first_shape", op._first), ("last", lasts, "_last_shape", op._last)):
        for source, target in targets:
            for f in found(target):
                setattr(op, shapes, getattr(op, shapes) | f)
                per[source] = f if source not in per else per[source] | f
                ends.append(f)
    generated = []
    for e in sections:
        sides = [f for f in faces if f not in ends and any(_alike(x, e) for x in f.Edges())]
        if sides:
            op._generated[e] = _normalize(compound(sides))
            generated.extend(sides)
    if last_rest is not None:
        rest = [f for f in faces if f not in ends and f not in generated]
        if rest:
            op._last_shape = op._last_shape | _normalize(compound(rest))
            op._last[last_rest] = _normalize(compound(rest))
    for source, offsets in moved:
        made = [f for off in offsets for f in faces if _alike(f, source.moved(Location(off)))]
        if made:
            op._generated[source] = _normalize(compound(made))
        edges = result.Edges()
        for e in source.Edges():
            images = [x for off in offsets for x in edges if _alike(x, e.moved(Location(off)))]
            if images:
                op._generated[e] = _normalize(compound(images))


def _lies_on(piece: Shape, face: Shape, tol: float = 1e-6) -> bool:
    """Whether a result's face lies on an input's face: of its kind, its centre on it."""
    c = piece.Center()
    return piece.geomType() == face.geomType() and face.distance(Vertex.makeVertex(c.x, c.y, c.z)) <= tol * (1.0 + face.BoundingBox().DiagonalLength)


_UNPLACED = [[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0]]


def _node(s) -> int:
    """A handle's node, whatever its placement (the key kentron's history records)."""
    return s.located(_UNPLACED).node_hash()


def _recorded(history, name, operands, run, images: bool = False, deletes: bool = True) -> Shape:
    """`run()`'s result, a step of `history` recorded for it as kentron's booleans made it
    (each result face and edge from the operand faces and edges it came from): each operand
    face or edge modified into the result's ones made from it (where those are not itself,
    kept as it is), each face with none deleted; with `images` (an imprint), each face's
    images, itself included."""
    if history is None:
        return run()
    native.record_begin()
    result = None
    try:
        result = run()
    finally:
        pairs = native.record_end(result.wrapped if result is not None else operands[0].wrapped)
    made = {}
    for source, image in pairs:
        made.setdefault(_node(source), []).append(_node(image))
    placed = {_node(x.wrapped): x for x in list(result.Faces()) + list(result.Edges())}
    op = Op(name)
    history.append(op, name)
    for operand in operands:
        op._tracked.update(operand.Faces())
        op._tracked.update(operand.Edges())
        op._tracked.update(operand.Vertices())
        for sub in list(operand.Faces()) + list(operand.Edges()):
            key = _node(sub.wrapped)
            found = [placed[n] for n in dict.fromkeys(made.get(key, [])) if n in placed]
            # (A face lying on one of the other operand's is kept once, as that one: the
            # result's faces lying on it are its images too.)
            if not found and sub.ShapeType() == "Face":
                found = [r for r in result.Faces() if _lies_on(r, sub)]
            if not found:
                if sub.ShapeType() == "Face" and deletes:
                    op._deleted.append(sub)
                continue
            if images and sub.ShapeType() == "Face":
                op._images[sub] = _normalize(compound(found))
            # (Of its own kind, modified into them; others, as an edge's fillet face, it
            # generated.)
            same = [x for x in found if x.ShapeType() == sub.ShapeType()]
            other = [x for x in found if x.ShapeType() != sub.ShapeType()]
            if same and not (len(same) == 1 and _node(same[0].wrapped) == key):
                op._modified[sub] = _normalize(compound(same))
            if other:
                op._generated[sub] = _normalize(compound(other))
    return result


def isSubshape(s1: Shape, s2: Shape) -> bool:
    """True if `s1` is a sub-shape of `s2`."""
    return any(s1.isSame(x) for x in s2._sub(s1.ShapeType()))


def closest(s1: Shape, s2: Shape) -> Tuple[Vector, Vector]:
    """The closest points of two shapes."""
    _, a, b = native.distance(s1.wrapped, s2.wrapped)
    return Vector(a), Vector(b)


def edgeOn(*args, **kwargs) -> Shape:
    raise NotImplementedError("edges on faces from parameter-space curves")


def wireOn(*args, **kwargs) -> Shape:
    raise NotImplementedError("wires on faces from parameter-space curves")


def faceOn(*args, **kwargs) -> Shape:
    raise NotImplementedError("faces on faces from parameter-space curves")
