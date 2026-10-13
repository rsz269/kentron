"""Placeholders: pykentron does not provide OCP.TopoDS."""

from ._placeholder import placeholders

globals().update(placeholders("TopoDS", "TopoDS TopoDS_Shape TopoDS_Solid TopoDS_Face"))

from ._placeholder import fallback as _fallback  # noqa: E402

__getattr__ = _fallback("TopoDS")
