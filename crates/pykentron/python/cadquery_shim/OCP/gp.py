"""OCP.gp: OCCT's basic geometry classes, those pykentron provides in pure Python
(`pykentron.gp`); the others are placeholders."""

from pykentron.gp import *  # noqa: F401,F403

from ._placeholder import fallback as _fallback  # noqa: E402

__getattr__ = _fallback("gp")
