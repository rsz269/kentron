"""OCP.ShapeFix: fixing faces (`ShapeFix_Face`), which pykentron's need no fixing for; the
other names are placeholders."""

from ._placeholder import fallback as _fallback


class ShapeFix_Face:
    """A face fixed: kept as it is."""

    def __init__(self, face=None):
        self._face = face

    def Init(self, face) -> None:
        self._face = face

    def FixOrientation(self) -> bool:
        return False

    def Perform(self) -> bool:
        return False

    def Face(self):
        return self._face


__getattr__ = _fallback("ShapeFix")
