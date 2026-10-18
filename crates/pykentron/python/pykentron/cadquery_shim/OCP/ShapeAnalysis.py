"""OCP.ShapeAnalysis: edges joined into wires (`ShapeAnalysis_FreeBounds`); the other names
are placeholders."""

from ._placeholder import fallback as _fallback


class ShapeAnalysis_FreeBounds:
    """Free boundaries: here, edges joined end to end into wires."""

    @staticmethod
    def ConnectEdgesToWires_s(edges, tolerance: float, shared: bool, wires) -> None:
        from pykentron.shapes import Shape, _edges_to_wires

        for w in _edges_to_wires([Shape.cast(e) for e in edges], tolerance):
            wires.Append(w.wrapped)

    ConnectEdgesToWires = ConnectEdgesToWires_s


__getattr__ = _fallback("ShapeAnalysis")
