"""OCP.TopTools: sequences of shapes (`TopTools_HSequenceOfShape`); the other names are
placeholders."""

from ._placeholder import fallback as _fallback


class TopTools_HSequenceOfShape(list):
    """A sequence of shapes."""

    def Append(self, shape) -> None:
        self.append(shape)

    def Length(self) -> int:
        return len(self)

    def Size(self) -> int:
        return len(self)

    def Value(self, index: int):
        return self[index - 1]


TopTools_ListOfShape = TopTools_HSequenceOfShape

__getattr__ = _fallback("TopTools")
