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

from ._placeholder import fallback as _fallback  # noqa: E402

__getattr__ = _fallback("BRepBuilderAPI")
