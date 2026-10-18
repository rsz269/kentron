"""OCP.gce: constructions of basic geometry (`gce_MakeLin`); the other names are
placeholders."""

from ._placeholder import fallback as _fallback


class gce_MakeLin:
    """A line through a point along a direction (or through two points)."""

    def __init__(self, a, b=None):
        from pykentron.gp import gp_Dir, gp_Lin, gp_Pnt

        if isinstance(b, gp_Pnt):
            b = gp_Dir(b.X() - a.X(), b.Y() - a.Y(), b.Z() - a.Z())
        self._line = gp_Lin(a, b) if b is not None else gp_Lin(a)

    def IsDone(self) -> bool:
        return True

    def Value(self):
        return self._line

    Operator = Value


__getattr__ = _fallback("gce")
