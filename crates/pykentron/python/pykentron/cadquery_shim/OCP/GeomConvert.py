"""OCP.GeomConvert: curves as B-splines (`GeomConvert.CurveToBSplineCurve_s`), here of edges:
what `Edge._geomAdaptor().Curve().Curve()` stands for. The other names are placeholders."""

from ._placeholder import fallback as _fallback


class GeomConvert:
    @staticmethod
    def CurveToBSplineCurve_s(curve, *args):
        """The edge with its curve as a (rational) B-spline."""
        from pykentron import _kentron

        return _kentron.to_nurbs(getattr(curve, "wrapped", curve))


__getattr__ = _fallback("GeomConvert")
