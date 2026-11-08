"""OCP.Geom's B-spline curves and surfaces as data (poles, knots, degrees), answering OCCT's
queries about them: what CadQuery's `occ_impl.nurbs` makes of its B-splines and reads back
(the edges and faces themselves are kentron's). Other names are placeholders."""

from ._placeholder import fallback as _fallback


class _Pnt:
    """A pole, as OCCT's gp_Pnt answers for its coordinates."""

    def __init__(self, xyz):
        self._xyz = tuple(float(c) for c in xyz)

    def X(self):
        return self._xyz[0]

    def Y(self):
        return self._xyz[1]

    def Z(self):
        return self._xyz[2]

    def Coord(self):
        return self._xyz


class Geom_BSplineCurve:
    """A B-spline curve: poles (x, y, z), knots and their multiplicities, degree, periodic."""

    def __init__(self, poles, knots, multiplicities, degree, periodic=False):
        self._poles = [tuple(p) for p in poles]
        self._knots = [float(k) for k in knots]
        self._mults = [int(m) for m in multiplicities]
        self._degree = int(degree)
        self._periodic = bool(periodic)

    def IsPeriodic(self):
        return self._periodic

    def Degree(self):
        return self._degree

    def NbPoles(self):
        return len(self._poles)

    def Pole(self, i):
        return _Pnt(self._poles[i - 1])

    def Poles(self):
        return [_Pnt(p) for p in self._poles]

    def NbKnots(self):
        return len(self._knots)

    def Knots(self):
        return list(self._knots)

    def Multiplicities(self):
        return list(self._mults)


class Geom_BSplineSurface:
    """A B-spline surface: poles (rows along u), knots and multiplicities along u and v,
    degrees, periodic along each."""

    def __init__(self, poles, uknots, vknots, umults, vmults, udegree, vdegree, uperiodic=False, vperiodic=False):
        self._poles = [[tuple(p) for p in row] for row in poles]
        self._uknots, self._vknots = [float(k) for k in uknots], [float(k) for k in vknots]
        self._umults, self._vmults = [int(m) for m in umults], [int(m) for m in vmults]
        self._udegree, self._vdegree = int(udegree), int(vdegree)
        self._uperiodic, self._vperiodic = bool(uperiodic), bool(vperiodic)

    def IsUPeriodic(self):
        return self._uperiodic

    def IsVPeriodic(self):
        return self._vperiodic

    def UDegree(self):
        return self._udegree

    def VDegree(self):
        return self._vdegree

    def NbUPoles(self):
        return len(self._poles)

    def NbVPoles(self):
        return len(self._poles[0]) if self._poles else 0

    def Pole(self, i, j):
        return _Pnt(self._poles[i - 1][j - 1])

    def UKnots(self):
        return list(self._uknots)

    def VKnots(self):
        return list(self._vknots)

    def UMultiplicities(self):
        return list(self._umults)

    def VMultiplicities(self):
        return list(self._vmults)


__getattr__ = _fallback("Geom")
