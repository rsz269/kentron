"""CadQuery's `occ_impl.shapes` on pykentron: its shape classes (from `pykentron.shapes`) and
module functions (CadQuery's free-function API, `cadquery.func`).

The functions take and return shapes as CadQuery's do: results of one piece are that piece
(a shell of one open face that face, a compound of one shape that shape), several a
compound. `history` and `name` (CadQuery's operation history) are accepted and ignored.
"""

from __future__ import annotations

import math
from typing import Iterable, List, Literal, Optional, Sequence, Tuple, Union

from .. import _kentron as native
from ..geom import Location, Plane, Vector, VectorLike  # noqa: F401
from ..shapes import *  # noqa: F401,F403
from ..shapes import Compound, Edge, Face, Shape, Shell, Solid, Vertex, Wire, _edges_to_wires, _split_pieces, _tuple

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
    """The faces sewn into shells (a face of its own a shell of one face)."""
    if _sequence(s):
        s = tuple(s[0])
    faces = [f for el in s for f in _get(el, "Face")]
    if ctx:
        faces += ctx.Faces() if isinstance(ctx, Shape) else [f for c in ctx for f in c.Faces()]
    sewed = Shape.cast(native.sew([f.wrapped for f in faces], tol))
    rv_faces, rv_shells = _shape_to_faces_shells(sewed)
    shells = rv_shells + [Shape.cast(native.make_shell([f.wrapped])) for f in rv_faces]
    return _compound_or_shape(shells)


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
    result = s1
    for s in (s2,) + shapes:
        result = Shape.cast(native.union(result.wrapped, s.wrapped))
    return _compound_or_shape(result)


def cut(s1: Shape, s2: Shape, tol: float = 0.0, glue=None, history=None, name=None) -> Shape:
    return _compound_or_shape(Shape.cast(native.difference(s1.wrapped, s2.wrapped)))


def intersect(s1: Shape, s2: Shape, tol: float = 0.0, glue=None, history=None, name=None) -> Shape:
    return _compound_or_shape(Shape.cast(native.intersection(s1.wrapped, s2.wrapped)))


def split(s1: Shape, s2: Shape, tol: float = 0.0, history=None, name=None) -> Shape:
    """`s1` cut into pieces by `s2`."""
    return _compound_or_shape(_split_pieces(s1, [s2]))


def imprint(*shapes: Shape, tol: float = 0.0, glue="full", history=None, name=None) -> Shape:
    """The shapes split by one another (OCCT's general fuse): all the pieces."""
    return _compound_or_shape([Shape.cast(h) for _, h in native.cells([s.wrapped for s in shapes])])


def clean(s: Shape) -> Shape:
    return _compound_or_shape(s.clean())


def fill(s: Shape, constraints=()) -> Shape:
    raise NotImplementedError("filling surfaces")


def cap(s: Shape, ctx: Shape, constraints=()) -> Shape:
    raise NotImplementedError("filling surfaces")


def fillet(s: Shape, edges: Shape, r: float, history=None, name=None) -> Shape:
    target = _get_one(s, ("Shell", "Solid"))
    return _compound_or_shape(Shape.cast(native.fillet(target.wrapped, [e.wrapped for e in _get_edges(edges.edges())], r)))


def chamfer(s: Shape, edges: Shape, d: float, history=None, name=None) -> Shape:
    """`s`'s (a shell's or solid's) edges bevelled by `d` along both faces."""
    target = _get_one(s, ("Shell", "Solid"))
    spec = [(e.wrapped, d, d, None) for e in _get_edges(edges.edges())]
    return _compound_or_shape(Shape.cast(native.chamfer(target.wrapped, spec)))


def extrude(s: Shape, d, both: bool = False, history=None, name=None) -> Shape:
    """`s` (vertices, edges, wires or faces) swept along `d` (both ways if `both`)."""
    d = Vector(d)
    results = []
    for el in _get(s, ("Vertex", "Edge", "Wire", "Face")):
        if both:
            results.append(Shape.cast(native.extrude(el.moved(-d).wrapped, (d * 2).toTuple())))
        else:
            results.append(Shape.cast(native.extrude(el.wrapped, d.toTuple())))
    return _compound_or_shape(results)


def revolve(s: Shape, p, d, a: float = 360, history=None, name=None) -> Shape:
    results = [Shape.cast(native.revolve(el.wrapped, _tuple(p), _tuple(d), a)) for el in _get(s, ("Vertex", "Edge", "Wire", "Face"))]
    return _compound_or_shape(results)


def offset(s: Shape, t: float, cap: bool = True, both: bool = False, tol: float = 1e-6, history=None, name=None) -> Shape:
    """Faces or shells thickened by `t` along their normals (planar faces)."""
    def one(t):
        out = []
        for el in _get(s, ("Face", "Shell")):
            for f in el.Faces():
                out.append(f.thicken(t))
        return out

    if both:
        pos, neg = one(t), one(-t)
        both_ = [a.fuse(b) for a, b in zip(pos, neg)]
        return both_[0] if len(both_) == 1 else Compound.makeCompound(both_)
    return _compound_or_shape(one(t))


def offset2D(s: Shape, t: float, ctx=None, kind: str = "arc", closed: bool = True) -> Shape:
    results = []
    for w in _get_wires(s):
        results.extend(w.offset2D(t, kind))
    return _compound_or_shape(results)


def chamfer2D(s: Shape, verts: Shape, d: float) -> Shape:
    return _compound_or_shape(_get_one(s, "Face").chamfer2D(d, verts.vertices() if verts.ShapeType() != "Vertex" else [verts]))


def fillet2D(s: Shape, verts: Shape, r: float) -> Shape:
    vs = list(verts.Vertices())
    return _compound_or_shape(_get_one(s, "Face").fillet2D(r, vs))


def sweep(s, path: Shape, aux=None, cap: bool = False, transition: str = "transformed", history=None, name=None) -> Shape:
    """Edges, wires or faces swept along `path` (several sections, of matching wires or
    faces, blended along it)."""
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
    (vertices at the ends), capped into solids if `cap`."""
    if len(s) == 1 and isinstance(s[0], (list, tuple)):
        s = tuple(s[0])
    results = []
    face_lists = _get_face_lists(s)
    if face_lists:
        for el in face_lists:
            outer = [f.outerWire() if isinstance(f, Face) else f for f in el]
            result = Shape.cast(native.make_loft([w.wrapped for w in outer], True, ruled))
            if not any(isinstance(f, Vertex) for f in el):
                for j in range(len(el[0].innerWires())):
                    hole = Shape.cast(native.make_loft([f.innerWires()[j].wrapped for f in el], True, ruled))
                    result = result.cut(hole)
            results.append(result)
    else:
        for el in _get_wire_lists(s):
            results.append(Shape.cast(native.make_loft([w.wrapped for w in el], cap, ruled)))
    return _compound_or_shape(results)


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


def prism(*args, **kwargs) -> Shape:
    raise NotImplementedError("prismatic features")


def draft(*args, **kwargs) -> Shape:
    raise NotImplementedError("drafts")


def check(s: Shape, results=None, tol: Optional[float] = None, verbose: bool = True) -> bool:
    """True if `s` is valid; its problems added to `results` (if given)."""
    problems = native.problems(s.wrapped)
    if results is not None:
        results.clear()
        results.extend(([s], p) for p in problems)
    return not problems


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
