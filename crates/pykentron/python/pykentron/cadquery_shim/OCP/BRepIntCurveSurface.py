"""OCP.BRepIntCurveSurface: where a line crosses a shape's faces; the other names are
placeholders."""

from ._placeholder import fallback as _fallback


class BRepIntCurveSurface_Inter:
    """The points where a line (`gp_Lin`) crosses the faces of a shape, in the line's order:
    iterated with `More`, `Next`; `Pnt`, `Face` and `W` (the line's parameter) at each."""

    def __init__(self):
        self._hits = []
        self._at = 0

    def Init(self, shape, line, tolerance: float) -> None:
        from pykentron import _kentron

        p, d = line.Location(), line.Direction()
        self._hits = _kentron.line_intersections(shape, (p.X(), p.Y(), p.Z()), (d.X(), d.Y(), d.Z()), max(tolerance, 1e-7))
        self._at = 0

    def More(self) -> bool:
        return self._at < len(self._hits)

    def Next(self) -> None:
        self._at += 1

    def Pnt(self):
        from pykentron.gp import gp_Pnt

        return gp_Pnt(*self._hits[self._at][0])

    def Face(self):
        return self._hits[self._at][1]

    def W(self) -> float:
        return self._hits[self._at][2]

    def IsDone(self) -> bool:
        return True


__getattr__ = _fallback("BRepIntCurveSurface")
