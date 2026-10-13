"""The parts of OCP (OpenCascade's Python bindings) that CadQuery scripts call directly on
shapes, backed by pykentron. Other names import, but raise NotImplementedError when used."""

import importlib.abc
import importlib.machinery
import importlib.util
import sys
import types

from ._placeholder import Placeholder


class _PlaceholderModule(types.ModuleType):
    """An OCP module pykentron does not provide: every name in it a placeholder."""

    def __getattr__(self, name):
        if name.startswith("__"):
            raise AttributeError(name)
        cls = type(name, (Placeholder,), {"_name": f"{self.__name__[4:]}.{name}"})
        setattr(self, name, cls)
        return cls


class _PlaceholderFinder(importlib.abc.MetaPathFinder, importlib.abc.Loader):
    """Imports OCP modules with no file here as placeholder modules (scripts and tests import
    many OCP names they use only in some paths)."""

    def find_spec(self, fullname, path, target=None):
        if not fullname.startswith("OCP.") or fullname.count(".") != 1:
            return None
        if importlib.machinery.PathFinder.find_spec(fullname, __path__) is not None:
            return None
        return importlib.util.spec_from_loader(fullname, self)

    def create_module(self, spec):
        return _PlaceholderModule(spec.name)

    def exec_module(self, module):
        pass


if not any(isinstance(f, _PlaceholderFinder) for f in sys.meta_path):
    sys.meta_path.append(_PlaceholderFinder())


def __getattr__(name):
    # (`OCP.X` without importing it first, as OCP's own package allows.)
    if name.startswith("__"):
        raise AttributeError(name)
    return importlib.import_module(f"OCP.{name}")
