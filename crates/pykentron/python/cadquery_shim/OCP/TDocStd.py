"""Placeholders: pykentron does not provide OCP.TDocStd."""

from ._placeholder import placeholders

globals().update(placeholders("TDocStd", "TDocStd_Document"))

from ._placeholder import fallback as _fallback  # noqa: E402

__getattr__ = _fallback("TDocStd")
