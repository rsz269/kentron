"""OCP.ShapeUpgrade: merging faces on the same surfaces and edges on the same curves, as
pykentron's `unify_same_domain`; the other names are placeholders."""


class ShapeUpgrade_UnifySameDomain:
    """Faces on the same surfaces, and edges on the same curves, merged. (pykentron always
    merges both, and never keeps internal edges.)"""

    def __init__(self, shape=None, unifyEdges: bool = True, unifyFaces: bool = True, concatBSplines: bool = False):
        self._shape = shape
        self._result = None

    def Initialize(self, shape, unifyEdges: bool = True, unifyFaces: bool = True, concatBSplines: bool = False) -> None:
        self._shape = shape
        self._result = None

    def AllowInternalEdges(self, allow: bool) -> None:
        pass

    def SetLinearTolerance(self, tol: float) -> None:
        pass

    def SetAngularTolerance(self, tol: float) -> None:
        pass

    def Build(self) -> None:
        from pykentron import _kentron

        self._result = _kentron.unify_same_domain(self._shape)

    def Shape(self):
        if self._result is None:
            self.Build()
        return self._result

from ._placeholder import fallback as _fallback  # noqa: E402

__getattr__ = _fallback("ShapeUpgrade")
