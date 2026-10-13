"""Placeholders: pykentron does not provide OCP.XCAFDoc."""

from ._placeholder import placeholders

globals().update(placeholders("XCAFDoc", "XCAFDoc_DocumentTool XCAFDoc_ColorType XCAFDoc_ColorTool XCAFDoc_ShapeTool"))

from ._placeholder import fallback as _fallback  # noqa: E402

__getattr__ = _fallback("XCAFDoc")
