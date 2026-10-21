"""pykentron: kentron, the CAD kernel, from Python, with a CadQuery-compatible API.

Scripts written for CadQuery run on kentron by importing `pykentron` in its place (or through
the `cadquery` shim package next to this one, which re-exports it).

Memory against speed: `PYKENTRON_PROFILE` (`low-memory` or `fast`, the default) set before the
import chooses a profile (its allocator); `profile()` reads it.
"""

import os as _os


def _check_cpu() -> None:
    """The native module is built for x86-64-v2 (SSE4.2, SSSE3, POPCNT): on a processor
    without them (or an emulator hiding them), a clear error, not an illegal instruction."""
    # (os.uname, not platform, whose import alone is a fifth of this module's.)
    if not hasattr(_os, "uname") or _os.uname().machine.lower() not in ("x86_64", "amd64"):
        return
    try:
        with open("/proc/cpuinfo") as f:
            flags = next((line.split(":", 1)[1].split() for line in f if line.startswith("flags")), None)
    except OSError:
        return
    missing = [n for n in ("sse4_2", "ssse3", "popcnt") if flags is not None and n not in flags]
    if missing:
        raise ImportError(f"pykentron needs an x86-64-v2 processor (this one lacks {', '.join(missing)}); build it from source with RUSTFLAGS= for any x86-64")


_check_cpu()

from . import _kentron  # noqa: E402
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
