"""pykentron: kentron, the CAD kernel, from Python, with a CadQuery-compatible API.

Scripts written for CadQuery run on kentron by importing `pykentron` in its place (or through
the `cadquery` shim package next to this one, which re-exports it).

Memory against speed: `PYKENTRON_PROFILE` (`low-memory`, `balanced`, the default, or `fast`)
set before the import chooses a profile; `profile()` and `set_profile()` read and change it.
"""

from . import _kentron
from ._kentron import profile, set_profile
from .geom import BoundBox, Color, Location, Plane, Vector
from .shapes import Compound, CompSolid, Edge, Face, Mixin3D, Shape, Shell, Solid, Vertex, Wire
from .workplane import Workplane
from .assembly import Assembly
from . import exporters, importers, selectors

__version__ = _kentron.__version__

# CadQuery's older name for Workplane, kept by it.
CQ = Workplane

__all__ = [
    "Assembly", "BoundBox", "Color", "Compound", "CompSolid", "Edge", "Face", "Location", "Mixin3D",
    "CQ", "Plane", "Shape", "Shell", "Solid", "Vector", "Vertex", "Wire", "Workplane",
    "exporters", "importers", "selectors", "profile", "set_profile",
]
