"""Points, vectors, placements, planes and boxes, as CadQuery has them.

Placements (`Location`) are 3x4 matrices: rows `[m0, m1, m2, t]` of `p -> M p + t`, the form
kentron takes them in.
"""

from __future__ import annotations

import math
from typing import Iterable, Sequence, Tuple, Union

Tuple3 = Tuple[float, float, float]


class Vector:
    """A 3D point or direction (CadQuery's `Vector`)."""

    __slots__ = ("x", "y", "z")

    def __init__(self, *args):
        if len(args) == 3:
            x, y, z = args
        elif len(args) == 2:
            (x, y), z = args, 0.0
        elif len(args) == 1:
            a = args[0]
            if isinstance(a, Vector):
                x, y, z = a.x, a.y, a.z
            elif hasattr(a, "XYZ"):
                # (OCCT's points, vectors and directions: `gp`'s.)
                x, y, z = a.X(), a.Y(), a.Z()
            elif len(a) == 3:
                x, y, z = a
            elif len(a) == 2:
                (x, y), z = a, 0.0
            else:
                raise TypeError(f"expected 2 or 3 coordinates, got {a!r}")
        elif not args:
            x = y = z = 0.0
        else:
            raise TypeError(f"expected 0 to 3 arguments, got {len(args)}")
        self.x, self.y, self.z = float(x), float(y), float(z)

    # CadQuery's capitalised accessors.
    @property
    def X(self) -> float:
        return self.x

    @property
    def Y(self) -> float:
        return self.y

    @property
    def Z(self) -> float:
        return self.z

    @property
    def Length(self) -> float:
        return math.sqrt(self.x * self.x + self.y * self.y + self.z * self.z)

    @property
    def wrapped(self):
        return self

    def toTuple(self) -> Tuple3:
        return (self.x, self.y, self.z)

    def __iter__(self):
        return iter((self.x, self.y, self.z))

    def __getitem__(self, i: int) -> float:
        return (self.x, self.y, self.z)[i]

    def __len__(self) -> int:
        return 3

    def add(self, v) -> "Vector":
        v = Vector(v)
        return Vector(self.x + v.x, self.y + v.y, self.z + v.z)

    def sub(self, v) -> "Vector":
        v = Vector(v)
        return Vector(self.x - v.x, self.y - v.y, self.z - v.z)

    def multiply(self, s: float) -> "Vector":
        return Vector(self.x * s, self.y * s, self.z * s)

    def __add__(self, v):
        return self.add(v)

    __radd__ = __add__

    def __sub__(self, v):
        return self.sub(v)

    def __rsub__(self, v):
        return Vector(v).sub(self)

    def __mul__(self, s):
        return self.multiply(s)

    __rmul__ = __mul__

    def __truediv__(self, s):
        return Vector(self.x / s, self.y / s, self.z / s)

    def __neg__(self):
        return Vector(-self.x, -self.y, -self.z)

    def __abs__(self):
        return self.Length

    def dot(self, v) -> float:
        v = Vector(v)
        return self.x * v.x + self.y * v.y + self.z * v.z

    def cross(self, v) -> "Vector":
        v = Vector(v)
        return Vector(self.y * v.z - self.z * v.y, self.z * v.x - self.x * v.z, self.x * v.y - self.y * v.x)

    def normalized(self) -> "Vector":
        length = self.Length
        if length == 0:
            raise ValueError("cannot normalize a zero vector")
        return self / length

    def Center(self) -> "Vector":
        """The vector itself (as shapes have a centre)."""
        return self

    def getAngle(self, v) -> float:
        v = Vector(v)
        c = self.dot(v) / (self.Length * v.Length)
        return math.acos(max(-1.0, min(1.0, c)))

    def getSignedAngle(self, v) -> float:
        # Signed about -Z, as CadQuery (gp_Vec::AngleWithRef with -Z): positive clockwise.
        v = Vector(v)
        angle = self.getAngle(v)
        return angle if self.cross(v).z < 0 else -angle

    def distanceToLine(self):
        raise NotImplementedError("Have not needed this yet, but OCCT supports it!")

    def distanceToPlane(self):
        raise NotImplementedError("Have not needed this yet, but OCCT supports it!")

    def projectToLine(self, line) -> "Vector":
        line = Vector(line)
        return line * (self.dot(line) / line.dot(line))

    def projectToPlane(self, plane: "Plane") -> "Vector":
        n = plane.zDir
        return self - n * ((self - plane.origin).dot(n) / n.dot(n))

    def distance(self, v) -> float:  # not in CadQuery, but handy
        return (self - Vector(v)).Length

    def toPnt(self):
        """The point as OCCT's (`gp_Pnt`)."""
        from .gp import gp_Pnt

        return gp_Pnt(self.x, self.y, self.z)

    def toDir(self):
        """The direction as OCCT's (`gp_Dir`)."""
        from .gp import gp_Dir

        return gp_Dir(self.x, self.y, self.z)

    def transform(self, t) -> "Vector":
        """The point moved by `t` (a `Matrix` or `Location`)."""
        return Vector(*_apply(t.rows, self))

    def __eq__(self, other) -> bool:
        # As gp_Vec::IsEqual: lengths within 1e-5 and, unless one is that short, directions.
        if not isinstance(other, Vector):
            return False
        a, b = self.Length, other.Length
        if abs(a - b) > 1e-5:
            return False
        return a <= 1e-5 or b <= 1e-5 or self.getAngle(other) <= 1e-5

    def __getstate__(self):
        return (self.x, self.y, self.z)

    def __setstate__(self, state):
        self.x, self.y, self.z = state

    def __hash__(self):
        return hash((round(self.x, 6), round(self.y, 6), round(self.z, 6)))

    def __repr__(self) -> str:
        return f"Vector: ({self.x}, {self.y}, {self.z})"

    __str__ = __repr__


VectorLike = Union[Vector, Tuple3, Tuple[float, float]]


def _rows_identity():
    return [[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0]]


def _rotation_rows(axis, angle_rad) -> list:
    """The rows of a rotation about `axis` through the origin."""
    a = Vector(axis).normalized()
    c, s = math.cos(angle_rad), math.sin(angle_rad)
    t = 1 - c
    x, y, z = a.x, a.y, a.z
    return [
        [t * x * x + c, t * x * y - s * z, t * x * z + s * y, 0.0],
        [t * x * y + s * z, t * y * y + c, t * y * z - s * x, 0.0],
        [t * x * z - s * y, t * y * z + s * x, t * z * z + c, 0.0],
    ]


def _compose(a, b):
    """The rows of `a * b` (b applied first)."""
    out = []
    for i in range(3):
        row = []
        for j in range(4):
            v = sum(a[i][k] * b[k][j] for k in range(3))
            if j == 3:
                v += a[i][3]
            row.append(v)
        out.append(row)
    return out


def _invert(m):
    """The inverse of an affine map (rows)."""
    (a, b, c), (d, e, f), (g, h, i) = (r[:3] for r in m)
    det = a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)
    if det == 0:
        raise ValueError("the matrix is singular")
    inv = [
        [(e * i - f * h) / det, (c * h - b * i) / det, (b * f - c * e) / det],
        [(f * g - d * i) / det, (a * i - c * g) / det, (c * d - a * f) / det],
        [(d * h - e * g) / det, (b * g - a * h) / det, (a * e - b * d) / det],
    ]
    t = [m[0][3], m[1][3], m[2][3]]
    return [inv[r] + [-sum(inv[r][k] * t[k] for k in range(3))] for r in range(3)]


def _apply(m, p) -> Tuple3:
    return tuple(m[i][0] * p[0] + m[i][1] * p[1] + m[i][2] * p[2] + m[i][3] for i in range(3))


class Matrix:
    """A 3D affine transformation, as a 4x4 matrix whose last row is [0, 0, 0, 1]
    (CadQuery's `Matrix`).

    Made from nothing (the identity), or from rows: three or four of four values each.
    """

    def __init__(self, matrix=None):
        if matrix is None:
            self.rows = _rows_identity()
        elif isinstance(matrix, (Matrix, Location)):
            self.rows = [r[:] for r in matrix.rows]
        elif isinstance(matrix, (list, tuple)):
            if not (len(matrix) in (3, 4) and all(isinstance(r, (list, tuple)) and len(r) == 4 for r in matrix)):
                raise TypeError(f"Matrix constructor requires 2d list of 4x3 or 4x4, but got: {matrix!r}")
            if len(matrix) == 4 and tuple(matrix[3]) != (0, 0, 0, 1):
                raise ValueError(f"Expected the last row to be [0,0,0,1], but got: {matrix[3]!r}")
            self.rows = [[float(e) for e in r] for r in matrix[:3]]
        else:
            raise TypeError(f"Invalid param to matrix constructor: {matrix}")

    @property
    def wrapped(self):
        return self

    def rotateX(self, angle: float):
        """Rotates about X by `angle` (radians), before this transformation."""
        self.rows = _compose(self.rows, _rotation_rows((1, 0, 0), angle))

    def rotateY(self, angle: float):
        self.rows = _compose(self.rows, _rotation_rows((0, 1, 0), angle))

    def rotateZ(self, angle: float):
        self.rows = _compose(self.rows, _rotation_rows((0, 0, 1), angle))

    def inverse(self) -> "Matrix":
        return Matrix(_invert(self.rows))

    def multiply(self, other):
        """`other` (a `Vector`) transformed, or the transformation `self * other`."""
        if isinstance(other, Vector):
            return other.transform(self)
        return Matrix(_compose(self.rows, other.rows))

    def transposed_list(self) -> Sequence[float]:
        data = [r[:] for r in self.rows] + [[0.0, 0.0, 0.0, 1.0]]
        return [data[j][i] for i in range(4) for j in range(4)]

    def __getitem__(self, rc) -> float:
        if not isinstance(rc, tuple) or len(rc) != 2:
            raise IndexError("Matrix subscript must provide (row, column)")
        r, c = rc
        if not (0 <= r <= 3 and 0 <= c <= 3):
            raise IndexError(f"Out of bounds access into 4x4 matrix: {rc!r}")
        return self.rows[r][c] if r < 3 else [0.0, 0.0, 0.0, 1.0][c]

    def __repr__(self) -> str:
        t = self.transposed_list()
        body = ",\n        ".join(str(t[i::4]) for i in range(4))
        return f"Matrix([{body}])"

    def __getstate__(self):
        return [r[:] for r in self.rows]

    def __setstate__(self, state):
        self.rows = [list(r) for r in state]


class Location:
    """A placement: a rotation and translation (CadQuery's `Location`)."""

    def __init__(self, *args, **kwargs):
        rows = _rows_identity()
        if kwargs:
            # Keywords: CadQuery's (x, y, z, rx, ry, rz) form.
            names = ("x", "y", "z", "rx", "ry", "rz")
            values = dict(zip(names, args))
            values.update(kwargs)
            unknown = set(values) - set(names)
            if unknown:
                raise TypeError(f"unexpected arguments {sorted(unknown)}")
            args = tuple(values.get(n, 0) for n in names)
        if not args:
            pass
        elif len(args) == 1 and (isinstance(args[0], (Location, Matrix)) or type(args[0]).__name__ == "gp_Trsf"):
            rows = [r[:] for r in args[0].rows]
        elif len(args) == 1 and isinstance(args[0], Plane):
            p = args[0]
            rows = _plane_rows(p.origin, p.xDir, p.yDir, p.zDir)
        elif len(args) == 1 and isinstance(args[0], list) and len(args[0]) == 3 and len(args[0][0]) == 4:
            rows = [list(map(float, r)) for r in args[0]]
        elif len(args) == 1:
            v = Vector(args[0])
            rows[0][3], rows[1][3], rows[2][3] = v.x, v.y, v.z
        elif len(args) == 2 and isinstance(args[0], Plane):
            p, v = args[0], Vector(args[1])
            rows = _plane_rows(v, p.xDir, p.yDir, p.zDir)
        elif len(args) == 2:
            # Translation and three extrinsic X, Y, Z rotation angles (degrees).
            rows = _euler_rows(*args[1])
            v = Vector(args[0])
            rows[0][3], rows[1][3], rows[2][3] = v.x, v.y, v.z
        elif len(args) == 3 and not all(isinstance(a, (int, float)) for a in args):
            t, ax, angle = args
            rows = _rotation_rows(ax, math.radians(angle))
            v = Vector(t)
            rows[0][3], rows[1][3], rows[2][3] = v.x, v.y, v.z
        else:
            x, y, z, rx, ry, rz = (list(args) + [0.0] * 6)[:6]
            rows = _euler_rows(rx, ry, rz)
            rows[0][3], rows[1][3], rows[2][3] = float(x), float(y), float(z)
        self.rows = rows

    @property
    def wrapped(self):
        return self

    @property
    def inverse(self) -> "Location":
        return Location(_invert(self.rows))

    def __mul__(self, other: "Location") -> "Location":
        return Location(_compose(self.rows, other.rows))

    def __pow__(self, n: int) -> "Location":
        result = Location()
        base = self if n >= 0 else self.inverse
        for _ in range(abs(n)):
            result = result * base
        return result

    def IsEqual(self, other: "Location", tolerance: float = 1e-9) -> bool:
        """The same placement (as OCCT's TopLoc_Location, read from `loc.wrapped`)."""
        return all(abs(a - b) <= tolerance for r, q in zip(self.rows, other.rows) for a, b in zip(r, q))

    def Transformation(self):
        """The placement as OCCT's transformation (`gp_Trsf`): CadQuery reads it from
        `loc.wrapped`, which is the location itself here."""
        from .gp import gp_Trsf

        return gp_Trsf(self.rows)

    def _apply_point(self, p) -> Vector:
        return Vector(*_apply(self.rows, Vector(p)))

    def __getstate__(self):
        return [r[:] for r in self.rows]

    def __setstate__(self, state):
        self.rows = [list(r) for r in state]

    def toTuple(self):
        r = self.rows
        # Extrinsic X, Y, Z angles of R = Rz Ry Rx.
        ry = math.asin(max(-1.0, min(1.0, -r[2][0])))
        if abs(math.cos(ry)) > 1e-12:
            rx = math.atan2(r[2][1], r[2][2])
            rz = math.atan2(r[1][0], r[0][0])
        else:
            rx = math.atan2(-r[1][2], r[1][1])
            rz = 0.0
        return (r[0][3], r[1][3], r[2][3]), (math.degrees(rx), math.degrees(ry), math.degrees(rz))

    @property
    def plane(self) -> "Plane":
        r = self.rows
        return Plane((r[0][3], r[1][3], r[2][3]), (r[0][0], r[1][0], r[2][0]), (r[0][2], r[1][2], r[2][2]))

    def __repr__(self) -> str:
        t, a = self.toTuple()
        return f"Location({t}, {a})"


def _euler_rows(rx, ry, rz):
    """Extrinsic rotations about X, then Y, then Z (degrees)."""
    rows = _rows_identity()
    for axis, angle in (((1, 0, 0), rx), ((0, 1, 0), ry), ((0, 0, 1), rz)):
        if angle:
            rows = _compose(_rotation_rows(axis, math.radians(angle)), rows)
    return rows


def _plane_rows(origin, x, y, z):
    o, x, y, z = Vector(origin), Vector(x), Vector(y), Vector(z)
    return [[x.x, y.x, z.x, o.x], [x.y, y.y, z.y, o.y], [x.z, y.z, z.z, o.z]]


_NAMED_PLANES = {
    # name: (xDir, normal)
    "XY": ((1, 0, 0), (0, 0, 1)),
    "YZ": ((0, 1, 0), (1, 0, 0)),
    "ZX": ((0, 0, 1), (0, 1, 0)),
    "XZ": ((1, 0, 0), (0, -1, 0)),
    "YX": ((0, 1, 0), (0, 0, -1)),
    "ZY": ((0, 0, 1), (-1, 0, 0)),
    "front": ((1, 0, 0), (0, 0, 1)),
    "back": ((-1, 0, 0), (0, 0, -1)),
    "left": ((0, 0, 1), (-1, 0, 0)),
    "right": ((0, 0, -1), (1, 0, 0)),
    "top": ((1, 0, 0), (0, 1, 0)),
    "bottom": ((1, 0, 0), (0, -1, 0)),
}


def _default_x(n: Vector) -> Vector:
    """The x direction OCCT gives a coordinate system with only its normal (gp_Ax2)."""
    a, b, c = n
    aa, ab, ac = abs(a), abs(b), abs(c)
    if ab <= aa and ab <= ac:
        d = (-c, 0.0, a) if aa > ac else (c, 0.0, -a)
    elif aa <= ab and aa <= ac:
        d = (0.0, -c, b) if ab > ac else (0.0, c, -b)
    else:
        d = (-b, a, 0.0) if aa > ab else (b, -a, 0.0)
    return Vector(d).normalized()


def _named(name: str):
    def make(cls, origin=(0, 0, 0), xDir=None):
        plane = cls.named(name, origin)
        if xDir is not None:
            plane._setPlaneDir(xDir)
        return plane

    make.__name__ = name
    make.__doc__ = f"The {name} plane through `origin`."
    return classmethod(make)


class Plane:
    """A 2D coordinate system in space (CadQuery's `Plane`)."""

    _eq_tolerance_origin = 1e-6
    _eq_tolerance_dot = 1e-6

    def __init__(self, origin=(0, 0, 0), xDir=None, normal=(0, 0, 1)):
        if isinstance(origin, Location):
            p = origin.plane
            origin, xDir, normal = p.origin, p.xDir, p.zDir
        z = Vector(normal)
        if z.Length == 0.0:
            raise ValueError("normal should be non null")
        self.zDir = z.normalized()
        if xDir is None:
            xDir = _default_x(self.zDir)
        elif Vector(xDir).Length == 0.0:
            raise ValueError("xDir should be non null")
        self._setPlaneDir(xDir)
        self.origin = Vector(origin)

    def _setPlaneDir(self, xDir) -> None:
        x = Vector(xDir)
        # (Square to the normal, as the coordinate system it stands for.)
        x = x - self.zDir * x.dot(self.zDir)
        if x.Length == 0.0:
            raise ValueError("xDir should not be along the normal")
        self.xDir = x.normalized()
        self.yDir = self.zDir.cross(self.xDir).normalized()

    @classmethod
    def named(cls, stdName: str, origin=(0, 0, 0)) -> "Plane":
        try:
            x, n = _NAMED_PLANES[stdName]
        except KeyError:
            raise ValueError(f"Supported names are {list(_NAMED_PLANES)}")
        return cls(origin, x, n)

    XY = _named("XY")
    YZ = _named("YZ")
    ZX = _named("ZX")
    XZ = _named("XZ")
    YX = _named("YX")
    ZY = _named("ZY")
    front = _named("front")
    back = _named("back")
    left = _named("left")
    right = _named("right")
    top = _named("top")
    bottom = _named("bottom")

    @property
    def origin(self) -> Vector:
        return self._origin

    @origin.setter
    def origin(self, value) -> None:
        self._origin = Vector(value)

    @property
    def location(self) -> Location:
        return Location(self)

    @property
    def rG(self) -> Matrix:
        """Local to global coordinates."""
        return Matrix(_plane_rows(self.origin, self.xDir, self.yDir, self.zDir))

    @property
    def fG(self) -> Matrix:
        """Global to local coordinates."""
        return self.rG.inverse()

    def toWorldCoords(self, tuplePoint) -> Vector:
        return Vector(tuplePoint).transform(self.rG)

    def toLocalCoords(self, obj):
        """A point, or a shape, in this plane's coordinates."""
        from .shapes import Shape

        if isinstance(obj, Vector):
            return obj.transform(self.fG)
        if isinstance(obj, Shape):
            return obj.transformShape(self.fG)
        raise ValueError(f"Don't know how to convert type {type(obj)} to local coordinates")

    def setOrigin2d(self, x: float, y: float) -> None:
        self.origin = self.toWorldCoords((x, y))

    def rotated(self, rotate=(0, 0, 0)) -> "Plane":
        """A copy turned about its own X, Y and Z axes by the angles (degrees); the turns
        compose as CadQuery's (Rx Ry Rz)."""
        rx, ry, rz = (math.radians(a) for a in Vector(rotate))
        m = _rows_identity()
        for axis, angle in ((self.xDir, rx), (self.yDir, ry), (self.zDir, rz)):
            m = _compose(m, _rotation_rows(axis, angle))
        t = Matrix(m)
        return Plane(self.origin, self.xDir.transform(t), self.zDir.transform(t))

    def mirrorInPlane(self, listOfShapes, axis: str = "X"):
        """The shapes turned half about this plane's X (or Y) axis."""
        if axis == "X":
            a = self.xDir
        elif axis == "Y":
            a = self.yDir
        else:
            raise NotImplementedError
        rows = _rotation_rows(a, math.pi)
        o = self.origin
        for i in range(3):
            rows[i][3] = o[i] - sum(rows[i][k] * o[k] for k in range(3))
        return [s.transformShape(Matrix(rows)) for s in listOfShapes]

    def toPln(self):
        """The plane as OCCT's (`gp_Pln`)."""
        from .gp import gp_Ax3, gp_Pln

        return gp_Pln(gp_Ax3(self.origin.toPnt(), self.zDir.toDir(), self.xDir.toDir()))

    def __eq__(self, other) -> bool:
        cls = type(self)
        return (
            isinstance(other, Plane)
            and abs(self.origin - other.origin) < cls._eq_tolerance_origin
            and abs(self.zDir.dot(other.zDir) - 1) < cls._eq_tolerance_dot
            and abs(self.xDir.dot(other.xDir) - 1) < cls._eq_tolerance_dot
        )

    def __ne__(self, other) -> bool:
        return not self.__eq__(other)

    __hash__ = None

    def __getstate__(self):
        return (self.xDir, self.yDir, self.zDir, self._origin)

    def __setstate__(self, data):
        self.xDir, self.yDir, self.zDir, self.origin = data

    def __repr__(self) -> str:
        return f"Plane(origin={self.origin.toTuple()}, xDir={self.xDir.toTuple()}, normal={self.zDir.toTuple()})"


TOL = 1e-2  # (CadQuery's default gap for boxes grown by `add`)


class BoundBox:
    """An axis-aligned box (CadQuery's `BoundBox`): the extremes, grown by a gap (as OCCT's
    Bnd_Box: a box's gap is the largest of those it was made from)."""

    def __init__(self, xmin, ymin, zmin, xmax, ymax, zmax, gap: float = 0.0):
        self._raw = (xmin, ymin, zmin, xmax, ymax, zmax)
        self.gap = gap
        self.xmin, self.ymin, self.zmin = xmin - gap, ymin - gap, zmin - gap
        self.xmax, self.ymax, self.zmax = xmax + gap, ymax + gap, zmax + gap
        self.xlen = self.xmax - self.xmin
        self.ylen = self.ymax - self.ymin
        self.zlen = self.zmax - self.zmin
        self.center = Vector((self.xmin + self.xmax) / 2, (self.ymin + self.ymax) / 2, (self.zmin + self.zmax) / 2)
        self.DiagonalLength = math.sqrt(self.xlen ** 2 + self.ylen ** 2 + self.zlen ** 2)

    def add(self, obj, tol: float = None) -> "BoundBox":
        """A box holding this and `obj` (a point or a box)."""
        # (As CadQuery's Bnd_Box: `tol` is the gap of the empty box this one is copied into,
        # so the copy keeps this box's gap instead.)
        if isinstance(obj, BoundBox):
            o, gap = obj._raw, max(self.gap, obj.gap)
        else:
            p = Vector(obj)
            o, gap = (p.x, p.y, p.z, p.x, p.y, p.z), self.gap
        s = self._raw
        return BoundBox(*(min(s[i], o[i]) for i in range(3)), *(max(s[i], o[i]) for i in range(3, 6)), gap=gap)

    def enlarge(self, tol: float) -> "BoundBox":
        """The box with its gap at least `tol` (as Bnd_Box::Enlarge)."""
        return BoundBox(*self._raw, gap=max(self.gap, abs(tol)))

    def isInside(self, b2: "BoundBox") -> bool:
        """True if `b2` lies strictly inside this box."""
        return (
            b2.xmin > self.xmin and b2.ymin > self.ymin and b2.zmin > self.zmin
            and b2.xmax < self.xmax and b2.ymax < self.ymax and b2.zmax < self.zmax
        )

    @classmethod
    def _fromTopoDS(cls, shape, tol: float = None, optimal: bool = True) -> "BoundBox":
        """The box of a native shape (CadQuery's, of an OCCT shape)."""
        from . import _kentron as native

        (x0, y0, z0), (x1, y1, z1) = native.bounding_box(shape)
        return cls(x0, y0, z0, x1, y1, z1)

    @staticmethod
    def findOutsideBox2D(bb1: "BoundBox", bb2: "BoundBox"):
        if bb1.xmin < bb2.xmin and bb1.xmax > bb2.xmax and bb1.ymin < bb2.ymin and bb1.ymax > bb2.ymax:
            return bb1
        if bb2.xmin < bb1.xmin and bb2.xmax > bb1.xmax and bb2.ymin < bb1.ymin and bb2.ymax > bb1.ymax:
            return bb2
        return None

    def __repr__(self) -> str:
        return f"BoundBox(({self.xmin}, {self.ymin}, {self.zmin}), ({self.xmax}, {self.ymax}, {self.zmax}))"


def _to_linear(c: float) -> float:
    """An sRGB component as a linear one."""
    return c / 12.92 if c <= 0.04045 else ((c + 0.055) / 1.055) ** 2.4


def _to_srgb(c: float) -> float:
    """A linear component as an sRGB one."""
    return c * 12.92 if c <= 0.0031308 else 1.055 * c ** (1 / 2.4) - 0.055


class Color:
    """An RGBA colour (CadQuery's `Color`), components from 0 to 1: by name, hex string
    ("#RGB", "#RGBA", "#RRGGBB", "#RRGGBBAA"), or components (sRGB, unless `srgb` is false:
    linear). Read as OCCT's from `wrapped`: `GetRGB()` gives the linear components."""

    _NAMED = {
        "red": (1, 0, 0), "green": (0, 1, 0), "blue": (0, 0, 1), "yellow": (1, 1, 0),
        "white": (1, 1, 1), "black": (0, 0, 0), "gray": (0.5, 0.5, 0.5), "grey": (0.5, 0.5, 0.5),
        "orange": (1, 0.65, 0), "cyan": (0, 1, 1), "magenta": (1, 0, 1),
    }

    def __init__(self, *args, **kwargs):
        if not args:
            r, g, b, a = 1.0, 1.0, 0.0, 1.0
        elif len(args) == 1 and isinstance(args[0], Color):
            r, g, b, a = args[0].toTuple()
        elif len(args) == 1 and isinstance(args[0], str):
            r, g, b, a = self._parse(args[0])
        elif len(args) in (3, 4, 5):
            r, g, b = (float(c) for c in args[:3])
            a = float(args[3]) if len(args) > 3 else float(kwargs.get("a") or 1.0)
            if len(args) == 5 and not args[4]:
                r, g, b = _to_srgb(r), _to_srgb(g), _to_srgb(b)
        else:
            raise ValueError(f"Unsupported arguments: {args}, {kwargs}")
        self._rgba = (float(r), float(g), float(b), float(a))

    @classmethod
    def _parse(cls, text: str):
        name = text.lower()
        if name in cls._NAMED:
            return (*cls._NAMED[name], 1.0)
        if text.startswith("#") and len(text) in (4, 5, 7, 9):
            digits = text[1:]
            step = 1 if len(digits) in (3, 4) else 2
            values = [int(digits[i:i + step] * (2 // step), 16) / 255.0 for i in range(0, len(digits), step)]
            if len(values) == 3:
                values.append(1.0)
            return tuple(values)
        raise ValueError(f"Unknown color name: {text}")

    def toTuple(self):
        return self._rgba

    def GetRGB(self) -> "_RGB":
        """The linear red, green and blue components (OCCT's Quantity_Color, as CadQuery
        reads them from `color.wrapped`)."""
        return _RGB(*(_to_linear(c) for c in self._rgba[:3]))

    def Alpha(self) -> float:
        return self._rgba[3]

    GetAlpha = Alpha

    @property
    def wrapped(self):
        return self

    def __eq__(self, other):
        return isinstance(other, Color) and self._rgba == other._rgba

    def __hash__(self):
        return hash(self._rgba)

    def __getstate__(self):
        return self._rgba

    def __setstate__(self, state):
        self._rgba = tuple(state)

    def __repr__(self) -> str:
        return f"Color{self._rgba}"


class _RGB:
    """A colour's components, read as OCCT's Quantity_Color's."""

    def __init__(self, r: float, g: float, b: float):
        self._rgb = (r, g, b)

    def Red(self) -> float:
        return self._rgb[0]

    def Green(self) -> float:
        return self._rgb[1]

    def Blue(self) -> float:
        return self._rgb[2]
