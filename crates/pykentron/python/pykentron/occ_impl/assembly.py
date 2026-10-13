"""CadQuery's `occ_impl.assembly` on pykentron: colours and materials of assembly nodes."""

from __future__ import annotations

from typing import Optional, Tuple

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
