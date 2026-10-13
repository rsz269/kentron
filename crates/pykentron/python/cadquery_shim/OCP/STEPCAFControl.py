"""Placeholders: pykentron does not provide OCP.STEPCAFControl."""

from ._placeholder import placeholders

globals().update(placeholders("STEPCAFControl", "STEPCAFControl_Reader STEPCAFControl_Writer"))

from ._placeholder import fallback as _fallback  # noqa: E402

__getattr__ = _fallback("STEPCAFControl")
