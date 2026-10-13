from pykentron import _kentron as _native


class BRepBndLib:
    @staticmethod
    def AddOptimal_s(shape, box, useTriangulation: bool = True, useShapeTolerance: bool = False):
        lo, hi = _native.bounding_box(shape)
        box._add(lo, hi)

    @staticmethod
    def Add_s(shape, box, useTriangulation: bool = True):
        lo, hi = _native.bounding_box(shape)
        box._add(lo, hi)

from ._placeholder import fallback as _fallback  # noqa: E402

__getattr__ = _fallback("BRepBndLib")
