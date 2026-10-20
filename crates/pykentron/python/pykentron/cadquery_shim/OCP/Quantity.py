"""OCP.Quantity on pykentron: the names of colours (Quantity_NameOfColor), as a colour's
`Name()` gives them; Quantity_Color is not provided."""

from ._placeholder import placeholders

globals().update(placeholders("Quantity", "Quantity_Color"))

from ._placeholder import fallback as _fallback  # noqa: E402


class Quantity_NameOfColor:
    """Names of colours, by their linear components (those `Color.GetRGB().Name()` picks
    from: the primaries, secondaries and greys)."""

    _components = {
        "Quantity_NOC_BLACK": (0.0, 0.0, 0.0),
        "Quantity_NOC_WHITE": (1.0, 1.0, 1.0),
        "Quantity_NOC_RED": (1.0, 0.0, 0.0),
        "Quantity_NOC_GREEN": (0.0, 1.0, 0.0),
        "Quantity_NOC_BLUE1": (0.0, 0.0, 1.0),
        "Quantity_NOC_YELLOW": (1.0, 1.0, 0.0),
        "Quantity_NOC_CYAN1": (0.0, 1.0, 1.0),
        "Quantity_NOC_MAGENTA1": (1.0, 0.0, 1.0),
        "Quantity_NOC_GRAY50": (0.214041, 0.214041, 0.214041),
    }

    @classmethod
    def _nearest(cls, rgb):
        name = min(cls._components, key=lambda n: sum((a - b) ** 2 for a, b in zip(cls._components[n], rgb)))
        return getattr(cls, name)


for _name in Quantity_NameOfColor._components:
    setattr(Quantity_NameOfColor, _name, _name)
Quantity_NameOfColor.Quantity_NOC_BLUE = Quantity_NameOfColor.Quantity_NOC_BLUE1
Quantity_NameOfColor.Quantity_NOC_CYAN = Quantity_NameOfColor.Quantity_NOC_CYAN1
Quantity_NameOfColor.Quantity_NOC_MAGENTA = Quantity_NameOfColor.Quantity_NOC_MAGENTA1
Quantity_NOC_RED = Quantity_NameOfColor.Quantity_NOC_RED

__getattr__ = _fallback("Quantity")
