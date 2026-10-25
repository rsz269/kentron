"""OCP.BRepProj: a wire projected onto a shape, along a direction or from a point."""

from ._placeholder import fallback as _fallback


class BRepProj_Projection:
    """The wires along which a wire, projected onto a shape along a direction (`gp_Dir`) or
    from a point (`gp_Pnt`, conically), meets it: the section of the shape by the surface the
    wire sweeps, joined into wires. Iterated with `More`, `Current` and `Next`."""

    def __init__(self, wire, shape, along):
        from pykentron import _kentron
        from pykentron.gp import gp_Dir
        from pykentron.shapes import Shape, _edges_to_wires

        source, target = Shape.cast(wire), Shape.cast(shape)
        a, b = source.BoundingBox(), target.BoundingBox()
        lo = [min(x, y) for x, y in zip((a.xmin, a.ymin, a.zmin), (b.xmin, b.ymin, b.zmin))]
        hi = [max(x, y) for x, y in zip((a.xmax, a.ymax, a.zmax), (b.xmax, b.ymax, b.zmax))]
        # (Far enough each way to pass through all of both.)
        reach = 2.0 * (1.0 + sum((h - l) ** 2 for l, h in zip(lo, hi)) ** 0.5)
        if isinstance(along, gp_Dir):
            d = (along.X(), along.Y(), along.Z())
            start = _kentron.transformed(source.wrapped, [[1, 0, 0, -reach * d[0]], [0, 1, 0, -reach * d[1]], [0, 0, 1, -reach * d[2]]])
            sweep = _kentron.extrude(start, tuple(2 * reach * x for x in d))
        else:
            # Rays from the point through the wire: the wire scaled about the point, from
            # near it out past the shape.
            c = (along.X(), along.Y(), along.Z())
            nearest = min(sum((p - q) ** 2 for p, q in zip(v.toTuple(), c)) ** 0.5 for v in source.Vertices()) if source.Vertices() else 1.0
            # (Only out past the shape's farthest corner from the point: rays carried far
            # beyond it make long surfaces whose small meetings with it are missed.)
            corners = [(x, y, z) for x in (b.xmin, b.xmax) for y in (b.ymin, b.ymax) for z in (b.zmin, b.zmax)]
            farthest = max(sum((p - q) ** 2 for p, q in zip(k, c)) ** 0.5 for k in corners)
            far = 1.1 * farthest / max(nearest, 1e-9)
            # (And in from well short of its nearest face: surfaces meeting it near their
            # narrow ends are missed there.)
            from pykentron.shapes import Compound, Vertex
            closest = Vertex.makeVertex(*c).distance(Compound.makeCompound(target.Faces())) if target.Faces() else 0.0
            widest = max(sum((p - q) ** 2 for p, q in zip(v.toTuple(), c)) ** 0.5 for v in source.Vertices()) if source.Vertices() else 1.0
            near = max(0.5 * closest / max(widest, 1e-9), 1e-3 * far)
            far = max(far, 1.5 * near)

            def scaled(shape, s):
                return _kentron.transformed(shape, [[s, 0, 0, c[0] * (1 - s)], [0, s, 0, c[1] * (1 - s)], [0, 0, s, c[2] * (1 - s)]])

            # (Ruled between each edge's two copies.)
            sweep = _kentron.make_compound([
                _kentron.make_loft([_kentron.make_wire([scaled(e.wrapped, near)]), _kentron.make_wire([scaled(e.wrapped, far)])], False, True)
                for e in source.Edges()
            ])
        edges = Shape.cast(_kentron.section(sweep, target.wrapped)).Edges()
        self._wires = [w.wrapped for w in _edges_to_wires(edges, 1e-6)] if edges else []
        self._at = 0

    def IsDone(self) -> bool:
        return True

    def Init(self) -> None:
        self._at = 0

    def More(self) -> bool:
        return self._at < len(self._wires)

    def Current(self):
        return self._wires[self._at]

    def Next(self) -> None:
        self._at += 1


__getattr__ = _fallback("BRepProj")
