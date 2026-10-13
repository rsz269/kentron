"""Placeholders: pykentron does not provide OCP.IFSelect."""

from ._placeholder import placeholders

globals().update(placeholders("IFSelect", "IFSelect_ReturnStatus"))

IFSelect_RetVoid, IFSelect_RetDone, IFSelect_RetError, IFSelect_RetFail, IFSelect_RetStop = range(5)

from ._placeholder import fallback as _fallback  # noqa: E402

__getattr__ = _fallback("IFSelect")
