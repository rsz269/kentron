"""OCP.BRepOffset: offsets of shapes (`BRepOffset_MakeOffset`, as CadQuery extensions use it
to thicken faces); the other names are placeholders."""

from ._placeholder import fallback as _fallback

BRepOffset_Skin = 0
BRepOffset_Pipe = 1
BRepOffset_RectoVerso = 2


class BRepOffset_MakeOffset:
    """A shape offset by a distance: thickened into a solid (`Thickening`), else its faces
    moved (a shell's or solid's)."""

    def __init__(self, *args, **kwargs):
        self._shape = None
        self._result = None
        if args or kwargs:
            self.Initialize(*args, **kwargs)

    def Initialize(self, S, Offset: float, Tol: float = 1e-7, Mode=BRepOffset_Skin, Intersection: bool = False, SelfInter: bool = False, Join=None, Thickening: bool = False, RemoveIntEdges: bool = False):
        self._shape, self._offset, self._thickening = S, Offset, Thickening

    def MakeOffsetShape(self) -> None:
        from pykentron import _kentron
        from pykentron.shapes import Shape, _thickened

        shape = Shape.cast(self._shape)
        if self._thickening:
            faces = shape.Faces()
            if len(faces) != 1:
                raise NotImplementedError("OCP.BRepOffset.BRepOffset_MakeOffset thickening several faces is not provided by pykentron")
            self._result = _thickened(faces[0], self._offset).wrapped
        else:
            self._result = _kentron.offset_shape(self._shape, self._offset, "intersection")

    def MakeThickSolid(self) -> None:
        self.MakeOffsetShape()

    def IsDone(self) -> bool:
        return self._result is not None

    def Shape(self):
        return self._result


__getattr__ = _fallback("BRepOffset")
