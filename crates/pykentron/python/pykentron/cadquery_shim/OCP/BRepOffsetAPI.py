"""OCP.BRepOffsetAPI: filling surfaces (`BRepOffsetAPI_MakeFilling`); the other names are
placeholders."""

from ._placeholder import fallback as _fallback


class BRepOffsetAPI_MakeFilling:
    """The face spanning a closed loop of edges, through points added within it: of least
    bending energy of order `Degree` (the approximation's settings are not needed)."""

    def __init__(self, Degree: int = 3, NbPtsOnCur: int = 15, NbIter: int = 2, Anisotropie: bool = False, Tol2d: float = 1e-5, Tol3d: float = 1e-4, TolAng: float = 1e-2, TolCurv: float = 0.1, MaxDeg: int = 8, MaxSegments: int = 9):
        self._degree = Degree
        self._edges = []
        self._points = []
        self._face = None

    def Add(self, constraint, *args) -> int:
        from pykentron.gp import gp_Pnt

        if isinstance(constraint, gp_Pnt):
            self._points.append((constraint.X(), constraint.Y(), constraint.Z()))
            return len(self._points)
        # (A support face (TopAbs 4) given: the surface to keep tangent to there.)
        if args and hasattr(args[0], "ShapeType") and int(args[0].ShapeType()) == 4:
            raise NotImplementedError("OCP.BRepOffsetAPI.BRepOffsetAPI_MakeFilling on support faces is not provided by pykentron")
        self._edges.append(constraint)
        return len(self._edges)

    def Build(self) -> None:
        from pykentron import _kentron

        self._face = _kentron.make_filling(self._edges, self._points, self._degree)

    def IsDone(self) -> bool:
        return self._face is not None

    def Shape(self):
        if self._face is None:
            self.Build()
        return self._face


__getattr__ = _fallback("BRepOffsetAPI")
