"""OCP.GeomAPI: intersections of surfaces (`GeomAPI_IntSS`, here of faces: what
`BRepAdaptor_Surface(face).Surface().Surface()` stands for); the other names are placeholders."""

from ._placeholder import fallback as _fallback


class GeomAPI_IntSS:
    """The curves where two surfaces meet: here the faces standing for them (as
    `BRepAdaptor_Surface` gives them), their section's edges."""

    def __init__(self, *args):
        self._lines = []
        if args:
            self.Perform(*args)

    def Perform(self, a, b, tolerance: float = 1e-7) -> None:
        from pykentron import _kentron
        from pykentron.shapes import Shape

        a, b = (getattr(x, "_face", x) for x in (a, b))
        self._lines = [e.wrapped for e in Shape.cast(_kentron.section(a, b)).Edges()]

    def IsDone(self) -> bool:
        return True

    def NbLines(self) -> int:
        return len(self._lines)

    def Line(self, index: int):
        return self._lines[index - 1]


__getattr__ = _fallback("GeomAPI")
