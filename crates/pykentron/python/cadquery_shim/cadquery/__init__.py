"""`cadquery` on pykentron: put this directory first on PYTHONPATH to run CadQuery scripts on
kentron.

CadQuery's own pure-Python layers (Workplane, Sketch, selectors, cqgi, free functions,
assemblies: see VENDORED.txt) run here unchanged on its `occ_impl` layer, which comes from
pykentron: each `cadquery.occ_impl.*` module without a file here takes its names from the
matching pykentron module, and names pykentron does not have yet are placeholders that raise
NotImplementedError when used (so that code importing more than it uses runs)."""

import importlib
import importlib.abc
import importlib.machinery
import importlib.util
import sys
import types as _pytypes  # (as the package's own `types` module takes the name)

# Where the names of each CadQuery module are found in pykentron, in order.
_SOURCES = {
    "cadquery.occ_impl.shapes": ["pykentron.shapes", "pykentron"],
    "cadquery.occ_impl.geom": ["pykentron.geom", "pykentron"],
    "cadquery.occ_impl.assembly": ["pykentron.assembly"],
    "cadquery.occ_impl.exporters": ["pykentron.exporters"],
    "cadquery.occ_impl.importers": ["pykentron.importers"],
}


def _placeholder(module: str, name: str):
    """A name pykentron does not provide: a class whose use (calling or instantiating it)
    raises NotImplementedError."""

    def fail(*args, **kwargs):
        raise NotImplementedError(f"{module}.{name} is not provided by pykentron")

    return type(name, (), {"__init__": fail, "__new__": lambda cls, *a, **k: fail(), "__class_getitem__": classmethod(lambda cls, item: cls)})


def _sources(module: str) -> list:
    """The pykentron modules a CadQuery module's names come from, in order: pykentron's own
    `occ_impl` modules first (CadQuery's names and signatures)."""
    sources = list(_SOURCES.get(module, []))
    if module.startswith("cadquery.occ_impl"):
        sources.insert(0, "pykentron" + module[len("cadquery"):])
    return sources


def _lookup(module: str, name: str):
    for source in _sources(module):
        try:
            mod = importlib.import_module(source)
        except ImportError:
            continue
        if hasattr(mod, name):
            return getattr(mod, name)
    return _placeholder(module, name)


class _ShimModule(_pytypes.ModuleType):
    """A CadQuery module pykentron has no file for: names resolved by `_lookup`."""

    def __getattr__(self, name):
        if name == "__all__":
            # (For star imports: the public names of the first pykentron module it draws on,
            # as CadQuery's module would give its own.)
            names = []
            for source in _sources(self.__name__):
                try:
                    mod = importlib.import_module(source)
                except ImportError:
                    continue
                public = getattr(mod, "__all__", None) or [n for n in vars(mod) if not n.startswith("_")]
                names = [n for n in public if not isinstance(getattr(mod, n, None), _pytypes.ModuleType)]
                break
            self.__all__ = names
            return names
        if name.startswith("__"):
            raise AttributeError(name)
        value = _lookup(self.__name__, name)
        setattr(self, name, value)
        return value


class _ShimFinder(importlib.abc.MetaPathFinder, importlib.abc.Loader):
    def find_spec(self, fullname, path, target=None):
        if not fullname.startswith("cadquery."):
            return None
        parent = sys.modules.get(fullname.rsplit(".", 1)[0])
        search = getattr(parent, "__path__", None)
        if search and importlib.machinery.PathFinder.find_spec(fullname, search) is not None:
            return None
        # (A package, so that its own submodules import too.)
        return importlib.util.spec_from_loader(fullname, self, is_package=True)

    def create_module(self, spec):
        module = _ShimModule(spec.name)
        module.__path__ = []
        return module

    def exec_module(self, module):
        pass


if not any(isinstance(f, _ShimFinder) for f in sys.meta_path):
    sys.meta_path.insert(0, _ShimFinder())


from pykentron import __version__  # noqa: E402,F401

# As CadQuery's own `cadquery/__init__.py`.
from .occ_impl.geom import Plane, BoundBox, Vector, Matrix, Location  # noqa: E402
from .occ_impl.shapes import Shape, Vertex, Edge, Face, Wire, Solid, Shell, Compound, sortWiresByBuildOrder  # noqa: E402
from .occ_impl import exporters  # noqa: E402
from .occ_impl import importers  # noqa: E402
from .selectors import (  # noqa: E402
    NearestToPointSelector,
    ParallelDirSelector,
    DirectionSelector,
    PerpendicularDirSelector,
    TypeSelector,
    DirectionMinMaxSelector,
    StringSyntaxSelector,
    Selector,
)
from .sketch import Sketch  # noqa: E402
from .cq import CQ, Workplane  # noqa: E402
from .assembly import Assembly, Color, Constraint, Material  # noqa: E402
from .types import UnitLiterals  # noqa: E402
from . import selectors  # noqa: E402
from . import plugins  # noqa: E402

__all__ = [
    "CQ", "Workplane", "Assembly", "Color", "Constraint", "Material", "plugins", "selectors", "Plane",
    "BoundBox", "Matrix", "Vector", "Location", "sortWiresByBuildOrder", "Shape", "Vertex", "Edge", "Wire",
    "Face", "Solid", "Shell", "Compound", "exporters", "importers", "NearestToPointSelector",
    "ParallelDirSelector", "DirectionSelector", "PerpendicularDirSelector", "TypeSelector",
    "DirectionMinMaxSelector", "StringSyntaxSelector", "Selector", "Sketch", "UnitLiterals",
]
