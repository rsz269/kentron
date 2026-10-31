"""OCP.BOPAlgo: the statuses a shape's check reports (`check`, pykentron's); the others are
placeholders."""

from pykentron.occ_impl.shapes import BOPAlgo_CheckStatus  # noqa: F401

from ._placeholder import fallback as _fallback  # noqa: E402

__getattr__ = _fallback("BOPAlgo")
