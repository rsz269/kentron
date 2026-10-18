"""OCCT's basic geometry classes (`gp`), the part of them CadQuery's code uses, in pure
Python: points, vectors, directions, axes, planes, lines, rigid transformations and
quaternions, with OCCT's method names (the `OCP.gp` the `cadquery` shim provides).
"""

from __future__ import annotations

import math
from typing import Tuple


def _xyz(args) -> Tuple[float, float, float]:
    if len(args) == 3:
        return tuple(float(a) for a in args)
    if len(args) == 1:
        a = args[0]
        if hasattr(a, "XYZ"):
            a = a.XYZ()
        if isinstance(a, gp_XYZ):
            return a._c
        return tuple(float(c) for c in a)
    if not args:
        return (0.0, 0.0, 0.0)
    raise TypeError(f"expected coordinates, got {args!r}")


class gp_XYZ:
    """A coordinate triple."""

    __slots__ = ("_c",)

    def __init__(self, *args):
        self._c = _xyz(args)

    def X(self) -> float:
        return self._c[0]

    def Y(self) -> float:
        return self._c[1]

    def Z(self) -> float:
        return self._c[2]

    def Coord(self, i: int = None):
        return self._c if i is None else self._c[i - 1]

    def SetCoord(self, *args) -> None:
        if len(args) == 2:
            c = list(self._c)
            c[args[0] - 1] = float(args[1])
            self._c = tuple(c)
        else:
            self._c = _xyz(args)

    def SetX(self, x: float) -> None:
        self.SetCoord(1, x)

    def SetY(self, y: float) -> None:
        self.SetCoord(2, y)

    def SetZ(self, z: float) -> None:
        self.SetCoord(3, z)

    def XYZ(self) -> "gp_XYZ":
        return gp_XYZ(*self._c)

    def Modulus(self) -> float:
        return math.sqrt(sum(c * c for c in self._c))

    def Added(self, o) -> "gp_XYZ":
        return gp_XYZ(*(a + b for a, b in zip(self._c, _xyz((o,)))))

    def Subtracted(self, o) -> "gp_XYZ":
        return gp_XYZ(*(a - b for a, b in zip(self._c, _xyz((o,)))))

    def Multiplied(self, s: float) -> "gp_XYZ":
        return gp_XYZ(*(a * s for a in self._c))

    __add__ = Added
    __sub__ = Subtracted
    __mul__ = Multiplied

    def IsEqual(self, other, tolerance: float) -> bool:
        return all(abs(a - b) <= tolerance for a, b in zip(self._c, _xyz((other,))))

    def __repr__(self) -> str:
        return f"{type(self).__name__}{self._c}"


def _dot(a, b) -> float:
    return sum(x * y for x, y in zip(a, b))


def _cross(a, b):
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])


def _norm(a) -> float:
    return math.sqrt(_dot(a, a))


def _angle(a, b) -> float:
    return math.atan2(_norm(_cross(a, b)), _dot(a, b))


class gp_Pnt(gp_XYZ):
    """A point."""

    __slots__ = ()

    def Distance(self, p) -> float:
        return _norm(tuple(a - b for a, b in zip(self._c, _xyz((p,)))))

    def SquareDistance(self, p) -> float:
        return self.Distance(p) ** 2

    def Translated(self, *args) -> "gp_Pnt":
        v = _xyz(args) if len(args) == 1 else _xyz((gp_Vec(*args),))
        return gp_Pnt(*(a + b for a, b in zip(self._c, v)))

    def Transformed(self, t: "gp_Trsf") -> "gp_Pnt":
        return gp_Pnt(*t._apply(self._c))

    def Transform(self, t: "gp_Trsf") -> None:
        self._c = t._apply(self._c)

    def IsEqual(self, other, tolerance: float) -> bool:
        return self.Distance(other) <= tolerance


class gp_Vec(gp_XYZ):
    """A vector (made from coordinates, a direction, or two points: from the first to the
    second)."""

    __slots__ = ()

    def __init__(self, *args):
        if len(args) == 2:
            a, b = _xyz((args[0],)), _xyz((args[1],))
            self._c = tuple(y - x for x, y in zip(a, b))
        else:
            self._c = _xyz(args)

    def Magnitude(self) -> float:
        return _norm(self._c)

    def SquareMagnitude(self) -> float:
        return _dot(self._c, self._c)

    def Normalized(self) -> "gp_Vec":
        m = self.Magnitude()
        if m == 0:
            raise ValueError("gp_Vec::Normalized: a null vector")
        return gp_Vec(*(c / m for c in self._c))

    def Normalize(self) -> None:
        self._c = self.Normalized()._c

    def Dot(self, o) -> float:
        return _dot(self._c, _xyz((o,)))

    def Crossed(self, o) -> "gp_Vec":
        return gp_Vec(*_cross(self._c, _xyz((o,))))

    def Added(self, o) -> "gp_Vec":
        return gp_Vec(*(a + b for a, b in zip(self._c, _xyz((o,)))))

    def Subtracted(self, o) -> "gp_Vec":
        return gp_Vec(*(a - b for a, b in zip(self._c, _xyz((o,)))))

    def Multiplied(self, s: float) -> "gp_Vec":
        return gp_Vec(*(a * s for a in self._c))

    def Divided(self, s: float) -> "gp_Vec":
        return gp_Vec(*(a / s for a in self._c))

    def Reversed(self) -> "gp_Vec":
        return gp_Vec(*(-a for a in self._c))

    __add__ = Added
    __sub__ = Subtracted
    __mul__ = Multiplied
    __rmul__ = Multiplied
    __truediv__ = Divided
    __neg__ = Reversed

    def Angle(self, o) -> float:
        return _angle(self._c, _xyz((o,)))

    def AngleWithRef(self, o, ref) -> float:
        """The angle to `o`, positive if their cross product points along `ref`."""
        o = _xyz((o,))
        a = _angle(self._c, o)
        return a if _dot(_cross(self._c, o), _xyz((ref,))) >= 0 else -a

    def IsEqual(self, other, linear: float, angular: float) -> bool:
        o = _xyz((other,))
        a, b = _norm(self._c), _norm(o)
        if abs(a - b) > linear:
            return False
        return a <= linear or b <= linear or _angle(self._c, o) <= angular

    def IsParallel(self, other, angular: float) -> bool:
        a = _angle(self._c, _xyz((other,)))
        return a <= angular or math.pi - a <= angular

    def Transformed(self, t: "gp_Trsf") -> "gp_Vec":
        return gp_Vec(*t._apply_vector(self._c))

    def Rotated(self, axis: "gp_Ax1", angle: float) -> "gp_Vec":
        """The vector turned by `angle` (radians) about the direction of `axis`."""
        t = gp_Trsf()
        t.SetRotation(axis, angle)
        return type(self)(*t._apply_vector(self._c))


class gp_Dir(gp_Vec):
    """A unit vector."""

    __slots__ = ()

    def __init__(self, *args):
        c = _xyz(args)
        m = _norm(c)
        if m <= 1e-300:
            raise ValueError("gp_Dir: a null vector")
        self._c = tuple(x / m for x in c)

    def Crossed(self, o) -> "gp_Dir":
        return gp_Dir(*_cross(self._c, _xyz((o,))))

    def Reversed(self) -> "gp_Dir":
        return gp_Dir(*(-a for a in self._c))

    __neg__ = Reversed

    def IsOpposite(self, other, angular: float) -> bool:
        return math.pi - _angle(self._c, _xyz((other,))) <= angular

    def IsNormal(self, other, angular: float) -> bool:
        return abs(math.pi / 2 - _angle(self._c, _xyz((other,)))) <= angular

    def Transformed(self, t: "gp_Trsf") -> "gp_Dir":
        return gp_Dir(*t._apply_vector(self._c))


def _default_x(n) -> gp_Dir:
    """The X direction OCCT gives a frame with only its main direction (gp_Ax2)."""
    a, b, c = n
    aa, ab, ac = abs(a), abs(b), abs(c)
    if ab <= aa and ab <= ac:
        d = (-c, 0.0, a) if aa > ac else (c, 0.0, -a)
    elif aa <= ab and aa <= ac:
        d = (0.0, -c, b) if ab > ac else (0.0, c, -b)
    else:
        d = (-b, a, 0.0) if aa > ab else (b, -a, 0.0)
    return gp_Dir(*d)


class gp_Ax1:
    """An axis: a point and a direction."""

    def __init__(self, p=None, d=None):
        self._p = gp_Pnt(p) if p is not None else gp_Pnt()
        self._d = gp_Dir(d) if d is not None else gp_Dir(0, 0, 1)

    def Location(self) -> gp_Pnt:
        return gp_Pnt(self._p)

    def Direction(self) -> gp_Dir:
        return gp_Dir(self._d)

    def SetLocation(self, p) -> None:
        self._p = gp_Pnt(p)

    def SetDirection(self, d) -> None:
        self._d = gp_Dir(d)


class gp_Ax2(gp_Ax1):
    """A right-handed frame: a point, a main direction and an X direction (OCCT's choice if
    not given)."""

    def __init__(self, p=None, n=None, x=None):
        super().__init__(p, n)
        self._x = self._square(x)

    def _square(self, x) -> gp_Dir:
        if x is None:
            return _default_x(self._d._c)
        x = _xyz((x,))
        z = self._d._c
        return gp_Dir(*_cross(_cross(z, x), z))

    def XDirection(self) -> gp_Dir:
        return gp_Dir(self._x)

    def YDirection(self) -> gp_Dir:
        return gp_Dir(*_cross(self._d._c, self._x._c))

    def Axis(self) -> gp_Ax1:
        return gp_Ax1(self._p, self._d)

    def SetDirection(self, d) -> None:
        self._d = gp_Dir(d)
        self._x = self._square(self._x)

    def SetXDirection(self, x) -> None:
        self._x = self._square(x)


class gp_Ax3(gp_Ax2):
    """A frame (as gp_Ax2, OCCT's possibly left-handed one kept right-handed here)."""


class gp_Pln:
    """A plane: through a point with a normal, or on a frame."""

    def __init__(self, *args):
        if len(args) == 1 and isinstance(args[0], gp_Ax2):
            a = args[0]
            self._ax = gp_Ax3(a.Location(), a.Direction(), a.XDirection())
        elif len(args) == 2:
            self._ax = gp_Ax3(args[0], args[1])
        else:
            self._ax = gp_Ax3()

    def Location(self) -> gp_Pnt:
        return self._ax.Location()

    def Axis(self) -> gp_Ax1:
        return self._ax.Axis()

    def Position(self) -> gp_Ax3:
        return self._ax

    def XAxis(self) -> gp_Ax1:
        return gp_Ax1(self._ax.Location(), self._ax.XDirection())

    def YAxis(self) -> gp_Ax1:
        return gp_Ax1(self._ax.Location(), self._ax.YDirection())

    def Distance(self, p) -> float:
        d = tuple(a - b for a, b in zip(_xyz((p,)), self._ax._p._c))
        return abs(_dot(d, self._ax._d._c))

    def Contains(self, *args) -> bool:
        """A point within `tolerance` of the plane, or a line (`linear`, `angular`
        tolerances) in it."""
        if isinstance(args[0], gp_Lin):
            line, linear, angular = args
            return self.Distance(line.Location()) <= linear and abs(math.pi / 2 - line.Direction().Angle(self._ax._d)) <= angular
        return self.Distance(args[0]) <= args[1]


class gp_Lin:
    """A line: a point and a direction."""

    def __init__(self, *args):
        if len(args) == 1:
            self._ax = gp_Ax1(args[0].Location(), args[0].Direction())
        else:
            self._ax = gp_Ax1(*args)

    def Location(self) -> gp_Pnt:
        return self._ax.Location()

    def Direction(self) -> gp_Dir:
        return self._ax.Direction()

    def Position(self) -> gp_Ax1:
        return self._ax

    def Distance(self, p) -> float:
        d = tuple(a - b for a, b in zip(_xyz((p,)), self._ax._p._c))
        return _norm(_cross(d, self._ax._d._c))


class gp_EulerSequence:
    gp_Extrinsic_XYZ = "Extrinsic_XYZ"
    gp_Intrinsic_XYZ = "Intrinsic_XYZ"
    gp_Extrinsic_ZYX = "Extrinsic_ZYX"
    gp_Intrinsic_ZYX = "Intrinsic_ZYX"


gp_Extrinsic_XYZ = gp_EulerSequence.gp_Extrinsic_XYZ
gp_Intrinsic_XYZ = gp_EulerSequence.gp_Intrinsic_XYZ
gp_Extrinsic_ZYX = gp_EulerSequence.gp_Extrinsic_ZYX
gp_Intrinsic_ZYX = gp_EulerSequence.gp_Intrinsic_ZYX


def _axis_quaternion(axis, angle) -> Tuple[float, float, float, float]:
    x, y, z = gp_Dir(*axis)._c
    s = math.sin(angle / 2)
    return (x * s, y * s, z * s, math.cos(angle / 2))


def _qmul(a, b):
    ax, ay, az, aw = a
    bx, by, bz, bw = b
    return (
        aw * bx + ax * bw + ay * bz - az * by,
        aw * by - ax * bz + ay * bw + az * bx,
        aw * bz + ax * by - ay * bx + az * bw,
        aw * bw - ax * bx - ay * by - az * bz,
    )


# The sequences' rotations: axes in the order applied, extrinsic (about fixed axes) or not.
_SEQUENCES = {
    "Extrinsic_XYZ": ("XYZ", True),
    "Intrinsic_XYZ": ("XYZ", False),
    "Extrinsic_ZYX": ("ZYX", True),
    "Intrinsic_ZYX": ("ZYX", False),
}
_AXES = {"X": (1, 0, 0), "Y": (0, 1, 0), "Z": (0, 0, 1)}


class gp_Quaternion:
    """A rotation as a unit quaternion (x, y, z, w)."""

    def __init__(self, *args):
        if len(args) == 4:
            self._q = tuple(float(a) for a in args)
        elif len(args) == 2:
            self._q = _axis_quaternion(_xyz((args[0],)), float(args[1]))
        else:
            self._q = (0.0, 0.0, 0.0, 1.0)

    def X(self) -> float:
        return self._q[0]

    def Y(self) -> float:
        return self._q[1]

    def Z(self) -> float:
        return self._q[2]

    def W(self) -> float:
        return self._q[3]

    def Set(self, x: float, y: float, z: float, w: float) -> None:
        self._q = (float(x), float(y), float(z), float(w))

    def SetVectorAndAngle(self, axis, angle: float) -> None:
        self._q = _axis_quaternion(_xyz((axis,)), angle)

    def SetEulerAngles(self, order, a: float, b: float, c: float) -> None:
        axes, extrinsic = _SEQUENCES[order]
        q = (0.0, 0.0, 0.0, 1.0)
        for axis, angle in zip(axes, (a, b, c)):
            r = _axis_quaternion(_AXES[axis], angle)
            q = _qmul(r, q) if extrinsic else _qmul(q, r)
        self._q = q

    def GetEulerAngles(self, order) -> Tuple[float, float, float]:
        m = self.GetMatrix()
        axes, extrinsic = _SEQUENCES[order]
        if (axes, extrinsic) == ("XYZ", True) or (axes, extrinsic) == ("ZYX", False):
            # R = Rz Ry Rx: the X, Y, Z angles (reversed for the intrinsic Z-Y-X order).
            ry = math.asin(max(-1.0, min(1.0, -m[2][0])))
            if abs(math.cos(ry)) > 1e-12:
                rx, rz = math.atan2(m[2][1], m[2][2]), math.atan2(m[1][0], m[0][0])
            else:
                rx, rz = math.atan2(-m[1][2], m[1][1]), 0.0
            return (rx, ry, rz) if extrinsic else (rz, ry, rx)
        # R = Rx Ry Rz.
        ry = math.asin(max(-1.0, min(1.0, m[0][2])))
        if abs(math.cos(ry)) > 1e-12:
            rx, rz = math.atan2(-m[1][2], m[2][2]), math.atan2(-m[0][1], m[0][0])
        else:
            rx, rz = math.atan2(m[2][1], m[1][1]), 0.0
        return (rx, ry, rz) if not extrinsic else (rz, ry, rx)

    def GetMatrix(self):
        x, y, z, w = self.Normalized()._q
        return [
            [1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)],
            [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
            [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)],
        ]

    def GetRotationAngle(self) -> float:
        """The angle, from -pi to pi: negative where w is (as OCCT's)."""
        x, y, z, w = self._q
        s = math.sqrt(x * x + y * y + z * z)
        return 2 * math.atan2(-s, -w) if w < 0 else 2 * math.atan2(s, w)

    def GetVectorAndAngle(self):
        x, y, z, _ = self._q
        s = math.sqrt(x * x + y * y + z * z)
        if s <= 1e-290:
            return gp_Vec(0, 0, 1), 0.0
        return gp_Vec(x / s, y / s, z / s), self.GetRotationAngle()

    def Normalized(self) -> "gp_Quaternion":
        n = math.sqrt(sum(c * c for c in self._q))
        return gp_Quaternion(*(c / n for c in self._q)) if n else gp_Quaternion()

    def Multiply(self, v):
        """`v` (a gp_Vec) rotated."""
        m = self.GetMatrix()
        c = _xyz((v,))
        return gp_Vec(*(_dot(r, c) for r in m))

    def __mul__(self, other: "gp_Quaternion") -> "gp_Quaternion":
        return gp_Quaternion(*_qmul(self._q, other._q))

    @staticmethod
    def _from_matrix(m) -> "gp_Quaternion":
        t = m[0][0] + m[1][1] + m[2][2]
        if t > 0:
            s = math.sqrt(t + 1.0) * 2
            q = ((m[2][1] - m[1][2]) / s, (m[0][2] - m[2][0]) / s, (m[1][0] - m[0][1]) / s, 0.25 * s)
        elif m[0][0] > m[1][1] and m[0][0] > m[2][2]:
            s = math.sqrt(1.0 + m[0][0] - m[1][1] - m[2][2]) * 2
            q = (0.25 * s, (m[0][1] + m[1][0]) / s, (m[0][2] + m[2][0]) / s, (m[2][1] - m[1][2]) / s)
        elif m[1][1] > m[2][2]:
            s = math.sqrt(1.0 + m[1][1] - m[0][0] - m[2][2]) * 2
            q = ((m[0][1] + m[1][0]) / s, 0.25 * s, (m[1][2] + m[2][1]) / s, (m[0][2] - m[2][0]) / s)
        else:
            s = math.sqrt(1.0 + m[2][2] - m[0][0] - m[1][1]) * 2
            q = ((m[0][2] + m[2][0]) / s, (m[1][2] + m[2][1]) / s, 0.25 * s, (m[1][0] - m[0][1]) / s)
        return gp_Quaternion(*q)


def _identity():
    return [[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0]]


class gp_Trsf:
    """A rigid motion, possibly with a uniform scale or a reflection: the rows of a 3×4
    matrix `p -> M p + t`."""

    def __init__(self, rows=None):
        self.rows = [list(map(float, r)) for r in rows] if rows is not None else _identity()

    def _apply(self, p):
        r = self.rows
        return tuple(r[i][0] * p[0] + r[i][1] * p[1] + r[i][2] * p[2] + r[i][3] for i in range(3))

    def _apply_vector(self, v):
        r = self.rows
        return tuple(r[i][0] * v[0] + r[i][1] * v[1] + r[i][2] * v[2] for i in range(3))

    def Value(self, row: int, col: int) -> float:
        return self.rows[row - 1][col - 1]

    def SetValues(self, *a) -> None:
        self.rows = [list(map(float, a[0:4])), list(map(float, a[4:8])), list(map(float, a[8:12]))]

    def TranslationPart(self) -> gp_XYZ:
        return gp_XYZ(*(r[3] for r in self.rows))

    def SetTranslationPart(self, v) -> None:
        for r, c in zip(self.rows, _xyz((v,))):
            r[3] = c

    def SetTranslation(self, *args) -> None:
        v = _xyz((gp_Vec(*args),)) if len(args) == 2 else _xyz(args)
        self.rows = _identity()
        self.SetTranslationPart(v)

    def _set_linear(self, m) -> None:
        for i in range(3):
            self.rows[i][:3] = list(m[i])

    def SetRotation(self, *args) -> None:
        """A rotation: by a quaternion (keeping the translation), or about an axis (gp_Ax1) by
        an angle (replacing the transformation)."""
        if len(args) == 1:
            self._set_linear(args[0].GetMatrix())
            return
        axis, angle = args
        m = gp_Quaternion(axis.Direction(), angle).GetMatrix()
        o = axis.Location()._c
        self.rows = [list(m[i]) + [o[i] - _dot(m[i], o)] for i in range(3)]

    def SetRotationPart(self, q: gp_Quaternion) -> None:
        self._set_linear(q.GetMatrix())

    def GetRotation(self) -> gp_Quaternion:
        s = self.ScaleFactor()
        m = [[c / s for c in r[:3]] for r in self.rows]
        return gp_Quaternion._from_matrix(m)

    def ScaleFactor(self) -> float:
        m = [r[:3] for r in self.rows]
        det = _dot(m[0], _cross(m[1], m[2]))
        return math.copysign(abs(det) ** (1.0 / 3.0), det) if det else 1.0

    def SetScale(self, p, s: float) -> None:
        o = _xyz((p,))
        self.rows = [[s if i == j else 0.0 for j in range(3)] + [o[i] * (1 - s)] for i in range(3)]

    def SetMirror(self, a) -> None:
        """The reflection in a point, an axis (gp_Ax1: a half turn about it) or a plane
        (gp_Ax2: in the plane normal to its direction)."""
        if isinstance(a, gp_Ax2):
            n, o = a.Direction()._c, a.Location()._c
            m = [[(1.0 if i == j else 0.0) - 2 * n[i] * n[j] for j in range(3)] for i in range(3)]
        elif isinstance(a, gp_Ax1):
            n, o = a.Direction()._c, a.Location()._c
            m = [[2 * n[i] * n[j] - (1.0 if i == j else 0.0) for j in range(3)] for i in range(3)]
        else:
            o = _xyz((a,))
            m = [[-1.0 if i == j else 0.0 for j in range(3)] for i in range(3)]
        self.rows = [list(m[i]) + [o[i] - _dot(m[i], o)] for i in range(3)]

    def SetTransformation(self, to, from_=None) -> None:
        """The change of coordinates from frame `from_` (the global one if not given) to
        frame `to`'s coordinates."""

        def frame_rows(ax):
            x, y, z, o = ax.XDirection()._c, ax.YDirection()._c, ax.Direction()._c, ax.Location()._c
            return [[x[i], y[i], z[i], o[i]] for i in range(3)]

        to_global = gp_Trsf(frame_rows(to))
        result = to_global.Inverted()
        if from_ is not None:
            result = result * gp_Trsf(frame_rows(from_))
        self.rows = result.rows

    def Multiplied(self, other: "gp_Trsf") -> "gp_Trsf":
        """`self * other` (`other` applied first)."""
        a, b = self.rows, other.rows
        out = []
        for i in range(3):
            row = [sum(a[i][k] * b[k][j] for k in range(3)) for j in range(4)]
            row[3] += a[i][3]
            out.append(row)
        return gp_Trsf(out)

    __mul__ = Multiplied

    def Multiply(self, other: "gp_Trsf") -> None:
        self.rows = self.Multiplied(other).rows

    def PreMultiply(self, other: "gp_Trsf") -> None:
        self.rows = other.Multiplied(self).rows

    def Inverted(self) -> "gp_Trsf":
        (a, b, c), (d, e, f), (g, h, i) = (r[:3] for r in self.rows)
        det = a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)
        inv = [
            [(e * i - f * h) / det, (c * h - b * i) / det, (b * f - c * e) / det],
            [(f * g - d * i) / det, (a * i - c * g) / det, (c * d - a * f) / det],
            [(d * h - e * g) / det, (b * g - a * h) / det, (a * e - b * d) / det],
        ]
        t = [r[3] for r in self.rows]
        return gp_Trsf([inv[r] + [-_dot(inv[r], t)] for r in range(3)])

    def Invert(self) -> None:
        self.rows = self.Inverted().rows

    def Powered(self, n: int) -> "gp_Trsf":
        base = self if n >= 0 else self.Inverted()
        out = gp_Trsf()
        for _ in range(abs(n)):
            out = out * base
        return out

    def __repr__(self) -> str:
        return f"gp_Trsf({self.rows})"


class gp:
    """OCCT's `gp` package functions."""

    @staticmethod
    def Origin_s() -> gp_Pnt:
        return gp_Pnt()

    @staticmethod
    def OX_s() -> gp_Ax1:
        return gp_Ax1(gp_Pnt(), gp_Dir(1, 0, 0))

    @staticmethod
    def OY_s() -> gp_Ax1:
        return gp_Ax1(gp_Pnt(), gp_Dir(0, 1, 0))

    @staticmethod
    def OZ_s() -> gp_Ax1:
        return gp_Ax1(gp_Pnt(), gp_Dir(0, 0, 1))

    @staticmethod
    def DX_s() -> gp_Dir:
        return gp_Dir(1, 0, 0)

    @staticmethod
    def DY_s() -> gp_Dir:
        return gp_Dir(0, 1, 0)

    @staticmethod
    def DZ_s() -> gp_Dir:
        return gp_Dir(0, 0, 1)

    @staticmethod
    def XOY_s() -> gp_Ax2:
        return gp_Ax2(gp_Pnt(), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0))

    @staticmethod
    def YOZ_s() -> gp_Ax2:
        return gp_Ax2(gp_Pnt(), gp_Dir(1, 0, 0), gp_Dir(0, 1, 0))

    @staticmethod
    def ZOX_s() -> gp_Ax2:
        return gp_Ax2(gp_Pnt(), gp_Dir(0, 1, 0), gp_Dir(0, 0, 1))

    @staticmethod
    def Resolution_s() -> float:
        return 1e-290
