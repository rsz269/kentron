"""OCP.BRepBuilderAPI: copying and transforming shapes, as pykentron provides them; the
other names are placeholders."""

import copy as _copy


class BRepBuilderAPI_Copy:
    """A copy of a shape. pykentron's geometry is immutable, so the copy is a new handle to
    the same shape (which `Location(loc)` can then place without moving the original)."""

    def __init__(self, shape=None, copyGeom: bool = True, copyMesh: bool = False):
        self._shape = None if shape is None else _copy.copy(shape)

    def Perform(self, shape, copyGeom: bool = True, copyMesh: bool = False) -> None:
        self._shape = _copy.copy(shape)

    def IsDone(self) -> bool:
        return self._shape is not None

    def Shape(self):
        return self._shape


class BRepBuilderAPI_Transform:
    """A shape transformed by a `gp_Trsf`: moved, if a rigid motion, else rebuilt."""

    def __init__(self, shape=None, trsf=None, copy: bool = False, copyMesh: bool = False):
        self._trsf = trsf
        self._shape = None
        if shape is not None:
            self.Perform(shape, copy)

    def Perform(self, shape, copy: bool = False) -> None:
        from pykentron import _kentron

        self._shape = _kentron.transformed(shape, self._trsf.rows)

    def IsDone(self) -> bool:
        return self._shape is not None

    def Shape(self):
        return self._shape

    def ModifiedShape(self, shape):
        return self._shape

class BRepBuilderAPI_Sewing:
    """Faces sewn into shells where their boundaries meet within the tolerance."""

    def __init__(self, tolerance: float = 1e-6, option1: bool = True, option2: bool = True, option3: bool = True, option4: bool = False):
        self._tolerance = tolerance
        self._shapes = []
        self._result = None

    def Add(self, shape) -> None:
        self._shapes.append(shape)

    def Perform(self, *args) -> None:
        from pykentron import _kentron

        self._result = _kentron.sew(self._shapes, self._tolerance)

    def SewedShape(self):
        if self._result is None:
            self.Perform()
        return self._result


class BRepBuilderAPI_MakeFace:
    """A face: from a closed planar wire, or a face again (`Add` cutting holes in it)."""

    def __init__(self, shape=None, *args):
        from pykentron import _kentron

        self._face = None
        if shape is None:
            return
        # (TopAbs numbers: 4 a face, 5 a wire.)
        kind = int(shape.ShapeType()) if hasattr(shape, "ShapeType") else None
        if kind == 4:
            self._face = shape
        elif kind == 5:
            self._face = _kentron.make_face(shape, [])
        else:
            raise NotImplementedError("OCP.BRepBuilderAPI.BRepBuilderAPI_MakeFace from surfaces is not provided by pykentron")

    def Add(self, wire) -> None:
        from pykentron import _kentron
        from pykentron.shapes import Shape

        # (A wire that does not close bounds no hole: OCCT adds it, the face then invalid;
        # here it is left out.)
        if not Shape.cast(wire).IsClosed():
            return
        self._face = _kentron.face_with_holes(self._face, [wire])

    def IsDone(self) -> bool:
        return self._face is not None

    def Face(self):
        return self._face

    Shape = Face


class BRepBuilderAPI_MakeEdge:
    """An edge: here, from an edge already made (what pykentron's curves of intersections
    are, `GeomAPI_IntSS.Line`), or between two points."""

    def __init__(self, *args):
        from pykentron import _kentron

        if len(args) == 1 and hasattr(args[0], "ShapeType"):
            self._edge = args[0]
        elif len(args) == 2 and all(hasattr(a, "X") for a in args):
            self._edge = _kentron.make_line(*[(a.X(), a.Y(), a.Z()) for a in args])
        else:
            raise NotImplementedError("OCP.BRepBuilderAPI.BRepBuilderAPI_MakeEdge from curves is not provided by pykentron")

    def IsDone(self) -> bool:
        return True

    def Edge(self):
        return self._edge

    Shape = Edge


class BRepBuilderAPI_MakeWire:
    """A wire of edges added end to end (not done where they do not join)."""

    def __init__(self, *edges):
        self._edges = []
        for e in edges:
            self.Add(e)

    def Add(self, shapes) -> None:
        from pykentron.shapes import Shape

        if hasattr(shapes, "ShapeType"):
            shapes = [shapes]
        for s in shapes:
            self._edges.extend(e.wrapped for e in Shape.cast(s).Edges())

    def _wire(self):
        from pykentron import _kentron

        try:
            return _kentron.make_wire(self._edges)
        except ValueError:
            return None

    def IsDone(self) -> bool:
        return bool(self._edges) and self._wire() is not None

    def Wire(self):
        return self._wire()

    Shape = Wire


from ._placeholder import fallback as _fallback  # noqa: E402

__getattr__ = _fallback("BRepBuilderAPI")
