"""Placeholders: pykentron does not provide OCP.TopExp."""

from ._placeholder import placeholders

globals().update(placeholders("TopExp", "TopExp_Explorer TopExp"))

from ._placeholder import fallback as _fallback  # noqa: E402

__getattr__ = _fallback("TopExp")
