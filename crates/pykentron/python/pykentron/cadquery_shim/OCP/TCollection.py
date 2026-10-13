"""Placeholders: pykentron does not provide OCP.TCollection."""

from ._placeholder import placeholders

globals().update(placeholders("TCollection", "TCollection_ExtendedString TCollection_AsciiString"))

from ._placeholder import fallback as _fallback  # noqa: E402

__getattr__ = _fallback("TCollection")
