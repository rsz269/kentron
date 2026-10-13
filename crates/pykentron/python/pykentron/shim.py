"""CadQuery scripts on pykentron, unchanged: makes `import cadquery` (and the few `OCP` names
scripts use directly) resolve to pykentron's CadQuery layer.

    import pykentron.shim          # first, before anything imports cadquery
    import cadquery as cq          # pykentron underneath

or, for a script left as it is: `pykentron-run script.py [args...]`.

The shim's `cadquery` package is CadQuery's own pure-Python layers (Workplane, Sketch,
selectors, cqgi, assemblies; copied from CadQuery, Apache-2.0) on pykentron's shapes; its
`OCP` package provides a few of OCP's names and placeholders for the rest. Nothing named
`cadquery` or `OCP` is installed: a real CadQuery installed beside pykentron is unaffected
unless the shim is activated.
"""

import os
import runpy
import sys

SHIM = os.path.join(os.path.dirname(os.path.abspath(__file__)), "cadquery_shim")


def activate() -> None:
    """Puts the shim first on `sys.path`. Fails if the real `cadquery` or `OCP` was imported
    already (the two cannot be mixed in one process)."""
    for name in ("cadquery", "OCP"):
        module = sys.modules.get(name)
        if module is not None and not (getattr(module, "__file__", None) or "").startswith(SHIM):
            raise ImportError(f"`{name}` was imported before pykentron.shim: import pykentron.shim first")
    if SHIM not in sys.path:
        sys.path.insert(0, SHIM)


activate()


def main() -> None:
    """`pykentron-run SCRIPT [ARGS...]`: runs a CadQuery script on pykentron."""
    if len(sys.argv) < 2:
        sys.exit("usage: pykentron-run SCRIPT [ARGS...]")
    script = os.path.abspath(sys.argv[1])
    sys.argv = sys.argv[1:]
    # The script's directory is importable, as when Python runs it (after the shim).
    sys.path.insert(1, os.path.dirname(script))
    runpy.run_path(script, run_name="__main__")
