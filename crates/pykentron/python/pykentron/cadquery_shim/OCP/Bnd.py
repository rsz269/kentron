class Bnd_Box:
    """An axis-aligned box, grown by `BRepBndLib.Add*_s`."""

    def __init__(self):
        self._box = None

    def IsVoid(self) -> bool:
        return self._box is None

    def Get(self):
        if self._box is None:
            raise ValueError("an empty box")
        (x0, y0, z0), (x1, y1, z1) = self._box
        return x0, y0, z0, x1, y1, z1

    def _add(self, lo, hi):
        if self._box is None:
            self._box = (tuple(lo), tuple(hi))
        else:
            (a, b) = self._box
            self._box = (tuple(map(min, a, lo)), tuple(map(max, b, hi)))

from ._placeholder import fallback as _fallback  # noqa: E402

__getattr__ = _fallback("Bnd")
