"""CadQuery's `vis` (rendering with VTK) where VTK (and numba, scipy, casadi) is installed: this module is then CadQuery's
own (`_vis_vtk`, vendored), as code that patches its names (tests faking VTK's windows) needs;
else `show` does nothing and its other names are placeholders.

CadQuery's own is loaded when first used, not on import: scripts import `show` from here at
the top and often never call it, and loading VTK, numba (with LLVM) and casadi costs 0.2 s and
some 140 MB. Until then `show` and `show_object` load it when called, and any other name read
or set here loads it and is read or set on it (so patching a name patches CadQuery's)."""

import importlib
import importlib.util
import sys
import types

# (With VTK and what CadQuery's NURBS module needs, numba, scipy and casadi: its `vis` imports it.)
if all(importlib.util.find_spec(m) is not None for m in ("vtkmodules", "numba", "scipy", "casadi")):

    def _load():
        """CadQuery's `vis`, loaded, and put in this module's place."""
        real = importlib.import_module("cadquery._vis_vtk")
        sys.modules[__name__] = real
        import cadquery

        cadquery.vis = real
        return real

    def show(*args, **kwargs):
        return _load().show(*args, **kwargs)

    def show_object(*args, **kwargs):
        return _load().show_object(*args, **kwargs)

    class _Lazy(types.ModuleType):
        def __getattr__(self, name):
            if name.startswith("__"):
                raise AttributeError(name)
            return getattr(_load(), name)

        def __setattr__(self, name, value):
            if name.startswith("__"):
                super().__setattr__(name, value)
            else:
                setattr(_load(), name, value)

        def __delattr__(self, name):
            delattr(_load(), name)

    sys.modules[__name__].__class__ = _Lazy
else:
    from pykentron.vis import show, show_object  # noqa: F401

    def __getattr__(name):
        """Other names of cadquery.vis: placeholders (see the `cadquery` shim)."""
        if name.startswith("__"):
            raise AttributeError(name)
        import cadquery

        return cadquery._placeholder("cadquery.vis", name)
