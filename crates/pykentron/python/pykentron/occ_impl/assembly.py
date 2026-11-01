"""CadQuery's `occ_impl.assembly` on pykentron: colours and materials of assembly nodes, the
assembly as JSON for VTK.js (`toJSON`) and its solids imprinted on one another (`imprint`)."""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple

from ..geom import Color  # noqa: F401


class Material:
    """A named material with a density (CadQuery's, OCCT's XCAFDoc_Material: "Default",
    7.85 g/cm^3 unless given)."""

    def __init__(self, name: Optional[str] = None, **kwargs):
        # (Visual properties: none yet, as CadQuery's.)
        self.wrapped_vis = _VisualMaterial()
        self._data = (
            name or "Default",
            kwargs.get("description", ""),
            float(kwargs.get("density", 7.85)),
            kwargs.get("densityUnit", "g/cm^3"),
        )

    @property
    def name(self) -> str:
        return self._data[0]

    @property
    def description(self) -> str:
        return self._data[1]

    @property
    def density(self) -> float:
        return self._data[2]

    @property
    def densityUnit(self) -> str:
        return self._data[3]

    def toTuple(self) -> Tuple[str, str, float, str]:
        return self._data

    def __hash__(self):
        return hash(self._data)

    def __eq__(self, other):
        return isinstance(other, Material) and self._data == other._data

    def __getstate__(self):
        return self._data

    def __setstate__(self, data):
        self._data = tuple(data)
        self.wrapped_vis = _VisualMaterial()

    def __repr__(self) -> str:
        return f"Material{self._data}"


class _VisualMaterial:
    """A material's visual properties (OCCT's XCAFDoc_VisMaterial): empty."""

    def IsEmpty(self) -> bool:
        return True


def toJSON(assy, color=(1.0, 1.0, 1.0, 1.0), edgecolor=(0.0, 0.0, 0.0, 1.0), tolerance: float = 1e-3) -> List[Dict[str, Any]]:
    """For each part of the assembly, its edges and its faces (VTK XML PolyData text), with
    their colours and the part's position and orientation (radians), as CadQuery's."""
    from math import radians

    from .exporters.vtk import toString

    rv = []
    for shape, _, loc, col in assy:
        edges, faces = toString(shape, tolerance)
        trans, rot = loc.toTuple()
        orientation = tuple(radians(r) for r in rot)
        rv.append({"shape": edges, "color": edgecolor, "position": trans, "orientation": orientation})
        rv.append({"shape": faces, "color": col.toTuple() if col else color, "position": trans, "orientation": orientation})
    return rv


def imprint(assy) -> Tuple[Any, Dict[Any, Tuple[str, ...]]]:
    """The assembly's solids, placed, split where they meet and sharing the faces they touch
    along (OCCT's BOPAlgo_MakeConnected), and for each solid of the result the names of the
    parts it comes from."""
    from .. import _kentron as native
    from ..shapes import Shape, Vertex
    from .shapes import _compound_or_shape

    named = []
    for obj, name, loc, _ in assy:
        for s in obj.moved(loc).Solids():
            named.append((s, name))
    if not named:
        raise ValueError("an assembly without solids")
    cells = [Shape.cast(h) for _, h in native.cells([s.wrapped for s, _ in named])]
    result = _compound_or_shape(cells)
    origins: Dict[Any, Tuple[str, ...]] = {}
    for piece in result.Solids():
        # (The parts holding its centre; where none does (a piece not convex), the nearest.)
        centre = piece.Center()
        names = tuple(name for s, name in named if s.isInside(centre))
        if not names:
            point = Vertex.makeVertex(*centre.toTuple())
            names = (min(named, key=lambda n: n[0].distance(point))[1],)
        origins[piece] = names
    return result, origins
