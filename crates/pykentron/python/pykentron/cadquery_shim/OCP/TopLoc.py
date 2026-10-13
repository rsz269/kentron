"""Placeholders: pykentron does not provide OCP.TopLoc."""

from ._placeholder import placeholders

globals().update(placeholders("TopLoc", "TopLoc_Location"))

from ._placeholder import fallback as _fallback  # noqa: E402

__getattr__ = _fallback("TopLoc")
