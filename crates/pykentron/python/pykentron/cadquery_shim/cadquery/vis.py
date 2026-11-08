"""CadQuery's `vis` (rendering with VTK) where VTK (and numba, scipy, casadi) is installed: this module is then CadQuery's
own (`_vis_vtk`, vendored), as code that patches its names (tests faking VTK's windows) needs;
else `show` does nothing and its other names are placeholders."""

import importlib
import importlib.util
import sys

# (With VTK and what CadQuery's NURBS module needs, numba, scipy and casadi: its `vis` imports it.)
if all(importlib.util.find_spec(m) is not None for m in ("vtkmodules", "numba", "scipy", "casadi")):
    sys.modules[__name__] = importlib.import_module("cadquery._vis_vtk")
else:
    from pykentron.vis import show, show_object  # noqa: F401

    def __getattr__(name):
        """Other names of cadquery.vis: placeholders (see the `cadquery` shim)."""
        if name.startswith("__"):
            raise AttributeError(name)
        import cadquery

        return cadquery._placeholder("cadquery.vis", name)
