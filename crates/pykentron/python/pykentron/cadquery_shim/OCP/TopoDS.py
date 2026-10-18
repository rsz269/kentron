"""OCP.TopoDS: only the empty shape (`TopoDS_Shape()`, which code written for OCP sets as a
shape's `wrapped` to make it null); the other names are placeholders."""

from ._placeholder import placeholders

globals().update(placeholders("TopoDS", "TopoDS TopoDS_Solid TopoDS_Face"))


class TopoDS_Shape:
    """The empty shape."""

    def IsNull(self) -> bool:
        return True

    def kind(self) -> str:
        raise ValueError("an empty shape")

from ._placeholder import fallback as _fallback  # noqa: E402

__getattr__ = _fallback("TopoDS")
