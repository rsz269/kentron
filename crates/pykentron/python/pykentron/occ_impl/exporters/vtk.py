"""CadQuery's `occ_impl.exporters.vtk` on pykentron, without VTK: a shape's mesh as VTK XML
PolyData text (`toString`), its edges and its faces apart."""

from __future__ import annotations

from typing import List, Sequence, Tuple

from ... import _kentron as native
from ...shapes import Shape


def _polydata(points: Sequence[Sequence[float]], lines: Sequence[Sequence[int]], polys: Sequence[Sequence[int]]) -> str:
    """VTK's XML PolyData (ASCII) of `points` and the `lines` (polylines) and `polys`
    (triangles) through them."""

    def cells(tag: str, items: Sequence[Sequence[int]]) -> str:
        offsets, at = [], 0
        for c in items:
            at += len(c)
            offsets.append(at)
        connectivity = " ".join(str(i) for c in items for i in c)
        return (
            f"<{tag}>"
            f'<DataArray type="Int64" Name="connectivity" format="ascii">{connectivity}</DataArray>'
            f'<DataArray type="Int64" Name="offsets" format="ascii">{" ".join(map(str, offsets))}</DataArray>'
            f"</{tag}>"
        )

    coords = " ".join(f"{x:.9g}" for p in points for x in p)
    return (
        '<?xml version="1.0"?>\n'
        '<VTKFile type="PolyData" version="0.1" byte_order="LittleEndian" header_type="UInt32">'
        "<PolyData>"
        f'<Piece NumberOfPoints="{len(points)}" NumberOfVerts="0" NumberOfLines="{len(lines)}" NumberOfStrips="0" NumberOfPolys="{len(polys)}">'
        f'<Points><DataArray type="Float32" Name="Points" NumberOfComponents="3" format="ascii">{coords}</DataArray></Points>'
        f"{cells('Lines', lines)}{cells('Polys', polys)}"
        "</Piece></PolyData></VTKFile>\n"
    )


def _edges(shape: Shape, tolerance: float) -> Tuple[List[Tuple[float, float, float]], List[List[int]]]:
    """Each edge of `shape` as a polyline, finer where it turns."""
    points, lines = [], []
    for e in shape.Edges():
        if e.geomType() == "LINE":
            n = 1
        else:
            # (Chords within the tolerance of an arc of the edge's length, at most 256.)
            n = max(4, min(256, int(e.Length() / max(8.0 * tolerance, 1e-9) ** 0.5) + 1))
        start = len(points)
        points.extend(e.positionAt(i / n).toTuple() for i in range(n + 1))
        lines.append(list(range(start, start + n + 1)))
    return points, lines


def toString(shape: Shape, tolerance: float = 1e-3, angularTolerance: float = 0.1) -> Tuple[str, str]:
    """The shape's edges and its faces, each as VTK XML PolyData text."""
    edge_points, lines = _edges(shape, tolerance)
    face_points, triangles = [], []
    if shape.Faces():
        face_points, triangles = native.triangulate(shape.wrapped, tolerance, angularTolerance)
    return _polydata(edge_points, lines, []), _polydata(face_points, [], [list(t) for t in triangles])
