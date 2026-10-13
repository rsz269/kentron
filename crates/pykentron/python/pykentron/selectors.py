"""String selectors on lists of shapes (the axis-based subset of CadQuery's).

* `|X`: straight edges along the axis, planar faces facing along it (either way).
* `#X`: straight edges across the axis, planar faces facing across it.
* `+X`, `-X`: straight edges and planar faces running or facing along the axis, or against it.
* `>X`, `<X`: the shapes farthest along the axis, or farthest against it (by centre).
* `%TYPE`: shapes of a geometry type (`%LINE`, `%CIRCLE`, `%PLANE`, `%CYLINDER`, ...).

Axes are `X`, `Y`, `Z`, or a vector given as `(x, y, z)`.
"""

from __future__ import annotations

from .geom import Vector

_AXES = {"X": Vector(1, 0, 0), "Y": Vector(0, 1, 0), "Z": Vector(0, 0, 1)}
TOLERANCE = 1e-4


def _axis(text: str) -> Vector:
    text = text.strip()
    if text.upper() in _AXES:
        return _AXES[text.upper()]
    if text.startswith("("):
        return Vector(*map(float, text.strip("()").split(","))).normalized()
    raise ValueError(f"unknown selector axis {text!r}")


def _direction(shape):
    """A straight edge's direction, or a planar face's normal; None otherwise."""
    kind = shape.wrapped.kind()
    if kind == "Edge":
        if shape.geomType() != "LINE":
            return None
        return (shape.endPoint() - shape.startPoint()).normalized()
    if kind == "Face":
        if shape.geomType() != "PLANE":
            return None
        return shape.normalAt()
    return None


def select(shapes: list, selector) -> list:
    if not isinstance(selector, str):
        return selector.filter(shapes)
    s = selector.strip()
    if s.startswith("%"):
        kind = s[1:].strip().upper()
        return [shape for shape in shapes if shape.geomType() == kind]
    op, axis = s[0], _axis(s[1:])
    if op in "<>":
        if not shapes:
            return []
        key = [shape.Center().dot(axis) for shape in shapes]
        best = max(key) if op == ">" else min(key)
        return [shape for shape, k in zip(shapes, key) if abs(k - best) < TOLERANCE]
    found = []
    for shape in shapes:
        d = _direction(shape)
        if d is None:
            continue
        # The edge's direction or the face's normal, against the axis.
        along = abs(d.dot(axis))
        if op == "|":
            ok = along > 1 - TOLERANCE
        elif op == "#":
            ok = along < TOLERANCE
        elif op == "+":
            ok = d.dot(axis) > 1 - TOLERANCE
        elif op == "-":
            ok = d.dot(axis) < -1 + TOLERANCE
        else:
            raise ValueError(f"unknown selector {selector!r}")
        if ok:
            found.append(shape)
    return found
