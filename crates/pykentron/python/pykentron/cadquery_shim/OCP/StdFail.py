"""OCP.StdFail: OCCT's failures of algorithms not done. kentron's failures surface as
ValueError, so these are ValueError here: code catching them catches kentron's."""

StdFail_NotDone = ValueError
StdFail_InfiniteSolutions = ValueError
StdFail_UndefinedDerivative = ValueError
StdFail_UndefinedValue = ValueError

from ._placeholder import fallback as _fallback  # noqa: E402

__getattr__ = _fallback("StdFail")
