"""Placeholders: pykentron does not provide OCP.TDF."""

from ._placeholder import placeholders

globals().update(placeholders("TDF", "TDF_Label TDF_LabelSequence"))

from ._placeholder import fallback as _fallback  # noqa: E402

__getattr__ = _fallback("TDF")
