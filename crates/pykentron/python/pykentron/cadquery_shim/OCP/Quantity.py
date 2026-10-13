"""Placeholders: pykentron does not provide OCP.Quantity."""

from ._placeholder import placeholders

globals().update(placeholders("Quantity", "Quantity_Color"))

from ._placeholder import fallback as _fallback  # noqa: E402

__getattr__ = _fallback("Quantity")
