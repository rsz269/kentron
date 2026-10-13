"""Shape kinds, as OCP's TopAbs."""

TopAbs_COMPOUND = 0
TopAbs_COMPSOLID = 1
TopAbs_SOLID = 2
TopAbs_SHELL = 3
TopAbs_FACE = 4
TopAbs_WIRE = 5
TopAbs_EDGE = 6
TopAbs_VERTEX = 7
TopAbs_SHAPE = 8

from ._placeholder import fallback as _fallback  # noqa: E402

__getattr__ = _fallback("TopAbs")
