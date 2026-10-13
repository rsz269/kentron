"""OCP.Standard: OCCT's exceptions. kentron's failures surface as ValueError, so OCCT's failure
exceptions are ValueError here: code catching them catches kentron's."""

Standard_Failure = ValueError
Standard_ConstructionError = ValueError
Standard_DomainError = ValueError
Standard_NoSuchObject = ValueError
Standard_OutOfRange = ValueError

from ._placeholder import fallback as _fallback  # noqa: E402

__getattr__ = _fallback("Standard")
