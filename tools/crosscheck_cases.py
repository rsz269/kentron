#!/usr/bin/env python3
"""Checks YAML test cases (docs/test-cases.md) on references other than kentron.

Each case's definition is interpreted by every leg, and each of its checks measured there:

  smlib  NVIDIA's USD BRep (the SMLib kernel): `pip install usd-brep`, in a Python 3.12 of its own
  occt   OCCT, through OCP (cadquery-ocp)
  brute  no kernel: the shapes as point-membership predicates written from the case's
         definition, with volume, centroid and bounds from quasi-random samples, area from
         lines crossing the boundary (Cauchy-Crofton), and solids counted by labelling a grid

The legs share nothing but the case's text. A shape the brute leg has no rule for takes a
predicate from the case's `brute:` map (ignored by kentron's runner), from a name to a numpy
expression over x, y, z and the names bound so far, called as `s(x, y, z)`:

  brute: { result: "s(x,y,z) & ~((x<1) & (y<1) & ((x-1)**2 + (y-1)**2 > 1))" }

Writing that expression from the shape's definition, not from a kernel's output, is what
keeps the leg independent. Each case runs in a process of its own: SMLib can crash or hang.

  tools/crosscheck_cases.py --smlib-python P --occt-python P [--only SUBSTRING] [-j N] [--json OUT] FILE...
"""
import argparse
import concurrent.futures as futures
import json
import math
import os
import re
import subprocess
import sys
import time

import yaml

SMLIB_DEFAULT_TOLERANCE = 2.0e-4  # What usd-brep's own tests allow on a volume.


# --- Reading a case --------------------------------------------------------------------------

def load_cases(path):
    with open(path) as f:
        return [d for d in yaml.safe_load_all(f) if d and "case" in d]


class Values:
    """The parts of a case's text every leg reads the same way: numbers, points, frames."""

    FUNCTIONS = {n: getattr(math, n) for n in "sqrt sin cos tan asin acos atan exp log atan2 pow".split()}
    FUNCTIONS.update(abs=abs, min=min, max=max, pi=math.pi)

    def __init__(self, params):
        self.params = {}
        for k, v in (params or {}).items():
            self.params[k] = self.number(v)

    def number(self, y):
        if isinstance(y, (int, float)):
            return float(y)
        text = str(y).replace("^", "**")
        return float(eval(text, {"__builtins__": {}}, {**self.FUNCTIONS, **self.params}))

    def point(self, y):
        return tuple(self.number(c) for c in y)

    def direction(self, y):
        if isinstance(y, str):
            sign = -1.0 if y.startswith("-") else 1.0
            return tuple(sign * float(y.lstrip("-") == a) for a in "xyz")
        v = self.point(y)
        n = math.sqrt(sum(c * c for c in v))
        return tuple(c / n for c in v)

    def frame(self, y):
        """(origin, z, x), x None when left to the kernel; None for the world frame."""
        if y is None:
            return None
        return (self.point(y.get("origin", [0, 0, 0])), self.direction(y.get("z", "z")),
                self.direction(y["x"]) if "x" in y else None)


class Selection:
    def __init__(self, parent, kind, index):
        self.parent, self.kind, self.index = parent, kind, index


class Interpreter:
    """Runs a case on a kernel (a leg's class), and measures its checks."""

    def __init__(self, kernel, case):
        self.kernel, self.case = kernel, case
        self.values = Values(case.get("params"))
        self.names = {}
        self.results = []

    def shape(self, y):
        value = self.names[y] if isinstance(y, str) else self.op(y)
        if isinstance(value, Selection):
            return self.kernel.subshapes(value.parent, value.kind)[value.index]
        return value

    def shapes(self, y):
        return [self.shape(i) for i in y] if isinstance(y, list) else [self.shape(y)]

    def op(self, y):
        if isinstance(y, str):
            return self.names[y]
        (name, args), = y.items()
        method = getattr(self, "op_" + name, None)
        if method is None:
            raise NotImplementedError("operation " + name)
        return method(args if args is not None else {})

    def run(self):
        for step in self.case.get("build", []):
            (key, value), = step.items()
            if key == "check":
                self.check(value)
            else:
                self.bind(key, value)
        for check in self.case.get("checks", []):
            self.check(check)
        return self.results

    def bind(self, name, value):
        try:
            self.names[name] = self.op(value)
        except NotImplementedError:
            self.names[name] = None  # (Whatever depends on it goes unmeasured on this leg.)
        except (AttributeError, TypeError, KeyError, RuntimeError):
            if not any(v is None for v in self.names.values()):
                raise
            self.names[name] = None

    def check(self, c):
        record = {k: v for k, v in c.items()}
        kind = next((k for k in ("valid", "volume", "area", "length", "count", "kind", "centroid", "closed",
                                 "bounds", "inertia", "distance", "contains", "exists", "faults", "tolerance") if k in c), None)
        if kind is None:
            return
        record["measure"] = kind
        try:
            target = self.names.get(c[kind]) if isinstance(c[kind], str) else None
            if isinstance(c[kind], str) and target is None:
                raise NotImplementedError("no value")
            if isinstance(target, Selection):
                sub = self.kernel.subshapes(target.parent, target.kind)[target.index]
                record["value"] = (list(self.kernel.sub_centroid(sub, target.kind)) if kind == "centroid"
                                   else self.kernel.sub_measure(sub, target.kind))
                self.results.append(record)
                return
            if kind == "distance":
                a, b = (self.names[n] for n in c[kind])
                if a is None or b is None:
                    raise NotImplementedError("no value")
                record["value"] = self.kernel.distance(a, b)
            elif kind == "contains":
                record["value"] = [self.kernel.contains(target, self.values.point(p)) for p in c["points"]]
            elif kind == "count":
                record["value"] = self.kernel.counts(target)
            elif kind in ("exists", "faults", "tolerance", "closed", "kind"):
                record["value"] = None
            else:
                record["value"] = getattr(self.kernel, kind)(target)
        except NotImplementedError as e:
            record["unavailable"] = str(e)
        self.results.append(record)

    # Operations, in the terms of docs/test-cases.md.

    def op_box(self, a):
        return self.kernel.box(self.values.point(a["size"]), self.values.point(a.get("corner", [0, 0, 0])),
                               self.values.frame(a.get("frame")))

    def _angle(self, a):
        return self.values.number(a["angle"]) if "angle" in a else None

    def op_cylinder(self, a):
        n = self.values.number
        return self.kernel.cylinder(n(a["radius"]), n(a["height"]), self._angle(a), self.values.frame(a.get("frame")))

    def op_cone(self, a):
        n = self.values.number
        return self.kernel.cone(n(a["r1"]), n(a["r2"]), n(a["height"]), self._angle(a), self.values.frame(a.get("frame")))

    def op_sphere(self, a):
        latitudes = [self.values.number(v) for v in a["latitudes"]] if "latitudes" in a else None
        return self.kernel.sphere(self.values.number(a["radius"]), latitudes, self._angle(a), self.values.frame(a.get("frame")))

    def op_torus(self, a):
        n = self.values.number
        return self.kernel.torus(n(a["major"]), n(a["minor"]), self._angle(a), self.values.frame(a.get("frame")))

    def op_wedge(self, a):
        size = self.values.point(a["size"])
        top = [self.values.number(v) for v in a["top"]] if "top" in a else [0, 0, self.values.number(a["top_dx"]), size[2]]
        return self.kernel.wedge(size, top)

    def op_translate(self, a):
        return self.kernel.translate(self.shape(a["shape"]), self.values.point(a["by"]))

    def op_rotate(self, a):
        axis = a["axis"]
        return self.kernel.rotate(self.shape(a["shape"]), self.values.point(axis.get("origin", [0, 0, 0])),
                                  self.values.direction(axis["direction"]), self.values.number(a["angle"]))

    def op_scale(self, a):
        return self.kernel.scale(self.shape(a["shape"]), self.values.point(a.get("center", [0, 0, 0])),
                                 self.values.number(a["factor"]))

    def op_mirror(self, a):
        return self.kernel.mirror(self.shape(a["shape"]), self.values.point(a.get("origin", [0, 0, 0])),
                                  self.values.direction(a["normal"]))

    def _boolean(self, op, a):
        fold = lambda shapes: shapes[0] if len(shapes) == 1 else self.kernel.compound(shapes)
        return getattr(self.kernel, op)(fold(self.shapes(a["object"])), fold(self.shapes(a["tool"])))

    def op_union(self, a):
        return self._boolean("union", a)

    def op_difference(self, a):
        return self._boolean("difference", a)

    def op_intersection(self, a):
        return self._boolean("intersection", a)

    def op_section(self, a):
        return self._boolean("section", a)

    def op_compound(self, a):
        return self.kernel.compound(self.shapes(a["of"]))

    def op_select(self, a):
        parent = self.names[a["of"]]
        (kind,) = [k for k in ("vertex", "edge", "face", "solid") if k in a]
        target = self.values.point(a[kind])
        measure = next((self.values.number(a[m]) for m in ("length", "area", "volume") if m in a), None)
        tolerance = 1e-5 * (1 + math.sqrt(sum(c * c for c in target)))
        fits = []
        for index, sub in enumerate(self.kernel.subshapes(parent, kind)):
            c = self.kernel.sub_centroid(sub, kind)
            if max(abs(p - q) for p, q in zip(c, target)) > tolerance:
                continue
            if measure is not None and abs(self.kernel.sub_measure(sub, kind) - measure) > 1e-5 * (1 + abs(measure)):
                continue
            fits.append(index)
        if len(fits) != 1:
            raise ValueError(f"select: {len(fits)} {kind}s at {target}")
        return Selection(parent, kind, fits[0])

    def op_fillet(self, a):
        n = self.values.number
        edges = a["edges"]
        if isinstance(edges, (str, dict)):
            edges = [edges]
        pairs = []
        if isinstance(edges[0], dict) and "law" in edges[0]:
            (u0, r0), (u1, r1) = edges[0]["law"]
            return self.kernel.fillet_law(self.shape(a["shape"]), self.names[edges[0]["edge"]], n(r0), n(r1))
        for e in edges:
            if isinstance(e, dict) and "edge" in e:
                pairs.append((self.names[e["edge"]], n(e["radius"])))
            else:
                pairs.append((self.names[e], n(a["radius"])))
        return self.kernel.fillet(self.shape(a["shape"]), pairs)

    def op_thick_solid(self, a):
        faces = a.get("faces", [])
        faces = [faces] if isinstance(faces, str) else faces
        return self.kernel.thick_solid(self.shape(a["shape"]), [self.names[f] for f in faces],
                                       self.values.number(a["thickness"]), a.get("join", "arc"))

    def op_offset_shape(self, a):
        return self.kernel.offset_shape(self.shape(a["shape"]), self.values.number(a["distance"]), a.get("join", "arc"))

    def op_polyline(self, a):
        return ("polygon", [self.values.point(p) for p in a["points"]])

    def op_circle(self, a):
        return ("circle", self.values.point(a.get("center", [0, 0, 0])), self.values.number(a["radius"]),
                self.values.direction(a.get("normal", "z")))

    def op_bspline(self, a):
        return ("bspline", [self.values.point(p) for p in a["poles"]], [self.values.number(k[0]) for k in a["knots"]],
                [int(k[1]) for k in a["knots"]], [self.values.number(v) for v in a.get("weights", [1] * len(a["poles"]))],
                int(a["degree"]))

    def op_helix(self, a):
        n = self.values.number
        return self.kernel.helix(n(a["start_radius"]), n(a.get("end_radius", a["start_radius"])), n(a["height"]),
                                 n(a["turns"]), bool(a.get("left_handed", False)))

    def op_offset_wire(self, a):
        return self.kernel.offset_wire(self.shape(a["shape"]), self.values.number(a["distance"]))

    def op_bezier(self, a):
        return ("bezier", [self.values.point(p) for p in a["poles"]])

    def op_loft(self, a):
        return self.kernel.loft([self.names[s] for s in a["sections"]], a.get("solid", False), a.get("ruled", False))

    def op_to_nurbs(self, a):
        return self.kernel.to_nurbs(self.shape(a["shape"]))

    def op_vertex(self, a):
        return self.kernel.vertex(self.values.point(a["at"]))

    def op_edge(self, a):
        if "curve" in a:
            curve = self.names[a["curve"]]
            if "range" in a:
                t0, t1 = (self.values.number(t) for t in a["range"])
                curve = ("arc", curve[1], curve[2], curve[3], t0, t1)
            return self.kernel.curve_edge(curve)
        return self.kernel.segment(self.values.point(a["from"]), self.values.point(a["to"]))

    def op_wire(self, a):
        edges = a["edges"] if isinstance(a["edges"], list) else [a["edges"]]
        return self.names[edges[0]] if len(edges) == 1 else self.kernel.wire([self.names[e] for e in edges])

    def op_face(self, a):
        return self.kernel.face(self.names[a["wire"]])

    def op_pipe(self, a):
        return self.kernel.pipe(self.shape(a["profile"]), self.shape(a["spine"]), a.get("corner"))

    def op_chamfer(self, a):
        edges = a["edges"] if isinstance(a["edges"], list) else [a["edges"]]
        return self.kernel.chamfer(self.shape(a["shape"]), [self.names[e] for e in edges], self.values.number(a["distance"]))

    def op_remove_faces(self, a):
        faces = a["faces"] if isinstance(a["faces"], list) else [a["faces"]]
        return self.kernel.remove_faces(self.shape(a["shape"]), [self.names[f] for f in faces])

    def op_extrude(self, a):
        if "angle" in a:
            return self.kernel.extrude_drafted(self.shape(a["profile"]), self.values.point(a["by"]), self.values.number(a["angle"]))
        return self.kernel.extrude(self.shape(a["profile"]), self.values.point(a["by"]))

    def op_revolve(self, a):
        axis = a["axis"]
        return self.kernel.revolve(self.shape(a["profile"]), self.values.point(axis.get("origin", [0, 0, 0])),
                                   self.values.direction(axis["direction"]), self.values.number(a.get("angle", 360)))


# --- SMLib -----------------------------------------------------------------------------------

class Section(list):
    """SMLib's section: curves, and the points where the operands only touch."""

    def __init__(self, curves, points):
        super().__init__(curves)
        self.points = [tuple(p) for p in points]


class Described:
    """Curves and edges kernels make from a description (polygon, circle): the shared part."""

    def curve_edge(self, description):
        return description

    def vertex(self, p):
        return ("point", p)

    def segment(self, p, q):
        return ("segment", p, q)

    def wire(self, edges):
        return ("wire", edges)

    def helix(self, r0, r1, height, turns, left_handed):
        raise NotImplementedError("helix")

    def offset_wire(self, s, distance):
        raise NotImplementedError("offset_wire")


class Smlib(Described):
    """usd_brep's operations consume their operands, and transforms work in place: every
    operation here works on copies, so that a name keeps its value."""

    def __init__(self):
        import usd_brep
        self.sm = usd_brep
        self.profiles = {}  # What each planar face was made from, by id.
        self.keep = []

    def copy(self, s):
        return self.sm.brep_copy(s)

    def planar(self, *shapes):
        return all(id(s) in self.profiles for s in shapes)

    def _planar_result(self, s):
        self.profiles[id(s)] = None
        self.keep.append(s)
        return s

    def _curves(self, profile):
        if profile[0] == "polygon":
            return [self.sm.create_line_segment(p, q) for p, q in zip(profile[1], profile[1][1:])]
        if profile[0] == "segment":
            return [self.sm.create_line_segment(profile[1], profile[2])]
        if profile[0] == "wire":
            return [c for e in profile[1] for c in self._curves(e)]
        if profile[0] == "curves":
            return [c for p in profile[1] for c in self._curves(p)]
        if profile[0] == "arc":
            return [self.sm.create_arc(profile[1], profile[2], math.degrees(profile[4]), math.degrees(profile[5]))]
        if profile[0] == "bspline":
            _, poles, knots, multiplicities, weights, degree = profile
            return [self.sm.create_canonical_curve(poles, knots, multiplicities, weights, degree=degree,
                                                   form=self.sm.BSplineCurveForm.CIRCULAR_ARC)]
        if profile[0] == "helix":
            return [profile[1]]
        if profile[0] == "bezier":
            return [self.sm.create_curve(profile[1])]  # (Four points make a cubic Bezier curve of them.)
        return [self.sm.create_circle(profile[1], profile[2])]

    def _swept(self, profile):
        """The curves of a profile, and whether they enclose an area to cap."""
        if isinstance(profile, tuple):
            return self._curves(profile), False
        return self._curves(self.profiles[id(profile)]), True

    def place(self, s, frame):
        if frame is None:
            return s
        origin, z, x = frame
        axis = (-z[1], z[0], 0.0)
        norm = math.hypot(axis[0], axis[1])
        if norm > 1e-12:
            self.sm.rotate(s, (0, 0, 0), axis, math.degrees(math.atan2(norm, z[2])))
        elif z[2] < 0:
            self.sm.rotate(s, (0, 0, 0), (1, 0, 0), 180.0)
        if x is not None and (z[0], z[1]) != (0.0, 0.0):
            raise NotImplementedError("a frame's x off a Z-turned frame")
        if x is not None:
            self.sm.rotate(s, (0, 0, 0), (0, 0, 1), math.degrees(math.atan2(x[1], x[0])) * (1 if z[2] > 0 else -1))
        self.sm.translate(s, origin)
        return s

    def box(self, size, corner, frame):
        return self.place(self.sm.create_box(corner, *size), frame)

    def cylinder(self, r, h, angle, frame):
        if angle is None or angle >= 360:
            s = self.sm.create_cylinder((0, 0, 0), r, h)
        else:
            s = self.sm.create_partial_cylinder((0, 0, 0), r, h, 0.0, angle)
        return self.place(s, frame)

    def cone(self, r1, r2, h, angle, frame):
        if angle is None or angle >= 360:
            s = self.sm.create_cone((0, 0, 0), r1, r2, h)
        else:
            s = self.sm.create_partial_cone((0, 0, 0), r1, r2, h, 0.0, angle)
        return self.place(s, frame)

    def sphere(self, r, latitudes, angle, frame):
        if latitudes is not None:
            raise NotImplementedError("a sphere between latitudes")
        if angle is None or angle >= 360:
            s = self.sm.create_sphere((0, 0, 0), r)
        else:
            s = self.sm.create_partial_sphere((0, 0, 0), r, 0.0, angle)
        return self.place(s, frame)

    def torus(self, major, minor, angle, frame):
        if angle is None or angle >= 360:
            s = self.sm.create_torus((0, 0, 0), major, minor)
        else:
            s = self.sm.create_partial_torus((0, 0, 0), major, minor, 0.0, angle)
        return self.place(s, frame)

    def wedge(self, size, top):
        # (A pyramid only: its height along Y, as a wedge's, its apex over the base's centre.)
        if size[0] != size[2] or top != [size[0] / 2, size[2] / 2, size[0] / 2, size[2] / 2]:
            raise NotImplementedError("wedge")
        s = self.sm.create_pyramid((0, 0, 0), size[0], size[1])
        self.sm.rotate(s, (0, 0, 0), (1, 0, 0), -90.0)
        self.sm.translate(s, (size[0] / 2, 0, size[2] / 2))
        return s

    def translate(self, s, v):
        s = self.copy(s)
        self.sm.translate(s, v)
        return s

    def rotate(self, s, origin, axis, angle):
        s = self.copy(s)
        self.sm.rotate(s, origin, axis, angle)
        return s

    def scale(self, s, center, factor):
        s = self.copy(s)
        self.sm.translate(s, tuple(-c for c in center))
        self.sm.scale(s, (factor, factor, factor))
        self.sm.translate(s, center)
        return s

    def mirror(self, s, origin, normal):
        axes = [i for i, c in enumerate(normal) if abs(c) > 1e-12]
        if len(axes) != 1:
            raise NotImplementedError("a mirror plane off the axes")
        s = self.copy(s)
        self.sm.translate(s, tuple(-c for c in origin))
        self.sm.scale(s, tuple(-1.0 if i == axes[0] else 1.0 for i in range(3)))
        self.sm.translate(s, origin)
        return s

    def _boolean_2d(self, a, b, operation):
        result = self._planar_result(self.sm.boolean_2d(self.copy(a), self.copy(b), operation=operation))
        self.profiles[id(result)] = ("curves", [self.profiles[id(a)], self.profiles[id(b)]])
        return result

    def union(self, a, b):
        if self.planar(a, b):
            return self._boolean_2d(a, b, self.sm.BooleanOp2D.UNION)
        return self.sm.boolean_union(self.copy(a), self.copy(b))

    def difference(self, a, b):
        if self.planar(a, b):
            return self._boolean_2d(a, b, self.sm.BooleanOp2D.DIFFERENCE)
        return self.sm.boolean_difference(self.copy(a), self.copy(b))

    def intersection(self, a, b):
        if self.planar(a, b):
            return self._boolean_2d(a, b, self.sm.BooleanOp2D.INTERSECTION)
        return self.sm.boolean_intersection(self.copy(a), self.copy(b))

    def section(self, a, b):
        if isinstance(a, tuple) and isinstance(b, tuple):
            points = [p for c in self._curves(a) for d in self._curves(b) for p in self.sm.intersect_curves(c, d)[0]]
            return Section([], points)
        if isinstance(b, tuple):
            a, b = b, a
        if isinstance(a, tuple):
            found = [(self.sm.intersect_curve_face(c, b) if isinstance(b, self.sm.Face) else self.sm.intersect_curve_brep(c, b))
                     for c in self._curves(a)]
            return Section([c for f in found for c in f[0]], [p for f in found for p in f[1]])
        curves, points = self.sm.intersect_breps(self.copy(a), self.copy(b))
        return Section(curves, points)

    def compound(self, shapes):
        raise NotImplementedError("compound")

    def subshapes(self, s, kind):
        return {"vertex": s.vertices, "edge": s.edges, "face": s.faces}[kind]()

    def sub_centroid(self, sub, kind):
        if kind == "vertex":
            return sub.point()
        if kind == "face":
            return sub.centroid(relative_accuracy=1e-7)
        # An edge's centre of length, from a dense polygon (usd_brep gives none).
        t0, t1 = sub.parameter_range()
        count = 2000
        points = [sub.point_at(t0 + (t1 - t0) * i / count) for i in range(count + 1)]
        weights = [math.dist(p, q) for p, q in zip(points, points[1:])]
        total = sum(weights)
        return tuple(sum(w * (p[i] + q[i]) / 2 for w, p, q in zip(weights, points, points[1:])) / total for i in range(3))

    def sub_measure(self, sub, kind):
        return sub.length() if kind == "edge" else self._accurately(lambda a: sub.area(relative_accuracy=a))

    def fillet(self, s, pairs):
        radii = {r for _, r in pairs}
        s = self.copy(s)
        edges = s.edges()
        chosen = [edges[sel.index] for sel, _ in pairs]
        if len(radii) == 1:
            return self.sm.fillet_edges(s, chosen, radius=radii.pop())
        return self.sm.fillet_edges_per_edge(s, chosen, [r for _, r in pairs],
                                             [self.sm.FilletXSect.CIRCULAR] * len(pairs))

    def helix(self, r0, r1, height, turns, left_handed):
        return ("helix", self.sm.create_helix((0, 0, 0), height, r0, r1, turns, right_handed=not left_handed))

    def fillet_law(self, s, sel, r0, r1):
        s = self.copy(s)
        return self.sm.variable_radius_fillet(s, [s.edges()[sel.index]], start_radius=r0, end_radius=r1)

    def pipe(self, profile, spine, corner=None):
        # (SMLib's sweep along a path mitres its corners: `corner: mitred`.)
        description = self.profiles[id(profile)]
        if spine[0] == "bezier" and description[0] == "circle":
            return self.sm.pipe_sweep(self._curves(spine)[0], radius=description[2])
        return self.sm.sweep_along_planar_path(self._curves(description), self._curves(spine), move_profile=True,
                                               cap_ends=self.sm.CapEnds.BOTH)

    def loft(self, sections, solid, ruled):
        curves = [self._curves(s) for s in sections]
        return self.sm.create_skin_primitive([len(c) for c in curves], [c for cs in curves for c in cs],
                                             degree=1 if ruled else 3, cap_ends=solid)

    def to_nurbs(self, s):
        s = self.copy(s)
        self.sm.turn_to_nurbs(s)
        return s

    def extrude_drafted(self, profile, by, angle):
        return self.sm.taper_extrude(self._curves(self.profiles[id(profile)]), by[2], angle)

    def chamfer(self, s, edges, distance):
        # (SMLib's is the whole edges' only, its radius the chamfer's width.)
        if len(edges) != len(s.edges()):
            raise NotImplementedError("a chamfer of some edges")
        return self.sm.chamfer_fillet(self.copy(s), distance * math.sqrt(2))

    def remove_faces(self, s, faces):
        # (SMLib finds the fillets itself: the case's faces are to be all of them.)
        s = self.copy(s)
        self.sm.remove_fillet(s)
        return s

    def thick_solid(self, s, faces, thickness, join):
        s = self.copy(s)
        all_faces = s.faces()
        return self.sm.shell_brep(s, thickness, create_solid=False, faces_to_shell=[all_faces[f.index] for f in faces])

    def offset_shape(self, s, distance, join):
        return self.sm.offset_brep(self.copy(s), distance)

    def face(self, profile):
        if profile[0] == "circle" and abs(profile[3][2]) < 1 - 1e-12:
            s = self.place(self.sm.create_planar_circle((0, 0, 0), profile[2]), (profile[1], profile[3], None))
        else:
            s = self.sm.create_planar_faces(self._curves(profile))
        self.profiles[id(s)] = profile
        self.keep.append(s)
        return s

    def extrude(self, profile, by):
        length = math.sqrt(sum(c * c for c in by))
        curves, caps = self._swept(profile)
        return self.sm.linear_sweep(curves, tuple(c / length for c in by), length, cap_ends=caps)

    def revolve(self, profile, origin, axis, angle):
        curves, caps = self._swept(profile)
        return self.sm.rotational_sweep(curves, origin, axis, angle, cap_ends=caps)

    def _accurately(self, measure):
        # (A tight accuracy can fail, as it does on a whole sphere: the next is taken.)
        for accuracy in (1e-8, 1e-7, 1e-6, 1e-5, 1e-4):
            try:
                return measure(accuracy)
            except RuntimeError:
                continue
        raise RuntimeError("no accuracy succeeded")

    def volume(self, s):
        return self._accurately(lambda a: s.volume(relative_accuracy=a))

    def area(self, s):
        if self.planar(s):
            return math.fsum(self._accurately(lambda a: f.area(relative_accuracy=a)) for f in s.faces())
        return self._accurately(lambda a: s.area(relative_accuracy=a))

    def length(self, s):
        if isinstance(s, list):
            return sum(c.length() for c in s)
        if isinstance(s, tuple):
            return sum(c.length() for c in self._curves(s))
        return sum(e.length() for e in s.edges())

    def centroid(self, s):
        return list(self._accurately(lambda a: s.mass_properties(relative_accuracy=a))["centroid"])

    def bounds(self, s):
        lo, hi = s.bounding_box(tight=True)
        return [list(lo), list(hi)]

    def valid(self, s):
        return bool(s.is_manifold_solid())

    def counts(self, s):
        if isinstance(s, Section):
            return {"vertices": len(s.points)} if not len(s) else {}
        return {"faces": s.face_count(), "edges": s.edge_count(), "vertices": s.vertex_count()}

    def distance(self, a, b):
        if isinstance(b, tuple) and b[0] == "point":
            a, b = b, a
        if isinstance(a, tuple) and a[0] == "point":
            if isinstance(b, Section):
                return min(math.dist(a[1], p) for p in b.points)
            if isinstance(b, tuple):
                return min(self.sm.curve_closest_point(c, a[1])[2] for c in self._curves(b))
            if self.planar(b):
                return min(f.closest_point(a[1])[1] for f in b.faces())
            return self.sm.brep_closest_point(b, a[1])[1]
        return self.sm.brep_distance(a, b)

    def contains(self, s, p):
        raise NotImplementedError("contains")

    def inertia(self, s):
        return list(self._accurately(lambda a: s.mass_properties(relative_accuracy=a))["moments_of_inertia"])


# --- OCCT ------------------------------------------------------------------------------------

class Occt(Described):
    def __init__(self):
        import OCP
        from OCP import (BRep, BRepAlgoAPI, BRepBndLib, BRepBuilderAPI, BRepCheck, BRepClass3d, BRepExtrema,
                         BRepFilletAPI, BRepGProp, BRepOffset, BRepOffsetAPI, BRepPrimAPI, Bnd, GeomAbs, GProp,
                         TopAbs, TopExp, TopoDS, TopTools, gp)
        self.m = dict(locals())

    def __getattr__(self, name):
        return self.m[name]

    def pnt(self, p):
        return self.gp.gp_Pnt(*p)

    def ax2(self, frame):
        if frame is None:
            return self.gp.gp_Ax2()
        origin, z, x = frame
        if x is None:
            return self.gp.gp_Ax2(self.pnt(origin), self.gp.gp_Dir(*z))
        return self.gp.gp_Ax2(self.pnt(origin), self.gp.gp_Dir(*z), self.gp.gp_Dir(*x))

    def box(self, size, corner, frame):
        s = self.BRepPrimAPI.BRepPrimAPI_MakeBox(self.pnt(corner), *size).Shape()
        if frame is None:
            return s
        t = self.gp.gp_Trsf()
        t.SetTransformation(self.gp.gp_Ax3(self.ax2(frame)), self.gp.gp_Ax3())
        return self._moved(s, t)

    def cylinder(self, r, h, angle, frame):
        args = (self.ax2(frame), r, h) + ((math.radians(angle),) if angle is not None else ())
        return self.BRepPrimAPI.BRepPrimAPI_MakeCylinder(*args).Shape()

    def cone(self, r1, r2, h, angle, frame):
        args = (self.ax2(frame), r1, r2, h) + ((math.radians(angle),) if angle is not None else ())
        return self.BRepPrimAPI.BRepPrimAPI_MakeCone(*args).Shape()

    def sphere(self, r, latitudes, angle, frame):
        args = [self.ax2(frame), r]
        if latitudes is not None:
            args += [math.radians(v) for v in latitudes]
        if angle is not None:
            args.append(math.radians(angle))
        return self.BRepPrimAPI.BRepPrimAPI_MakeSphere(*args).Shape()

    def torus(self, major, minor, angle, frame):
        args = (self.ax2(frame), major, minor) + ((math.radians(angle),) if angle is not None else ())
        return self.BRepPrimAPI.BRepPrimAPI_MakeTorus(*args).Shape()

    def wedge(self, size, top):
        return self.BRepPrimAPI.BRepPrimAPI_MakeWedge(*size, *top).Shape()

    def _moved(self, s, t):
        return self.BRepBuilderAPI.BRepBuilderAPI_Transform(s, t, True).Shape()

    def translate(self, s, v):
        t = self.gp.gp_Trsf()
        t.SetTranslation(self.gp.gp_Vec(*v))
        return self._moved(s, t)

    def rotate(self, s, origin, axis, angle):
        t = self.gp.gp_Trsf()
        t.SetRotation(self.gp.gp_Ax1(self.pnt(origin), self.gp.gp_Dir(*axis)), math.radians(angle))
        return self._moved(s, t)

    def scale(self, s, center, factor):
        t = self.gp.gp_Trsf()
        t.SetScale(self.pnt(center), factor)
        return self._moved(s, t)

    def mirror(self, s, origin, normal):
        t = self.gp.gp_Trsf()
        t.SetMirror(self.gp.gp_Ax2(self.pnt(origin), self.gp.gp_Dir(*normal)))
        return self._moved(s, t)

    def _boolean(self, cls, a, b):
        op = cls(a, b)
        if not op.IsDone():
            raise RuntimeError("boolean failed")
        return op.Shape()

    def union(self, a, b):
        return self._boolean(self.BRepAlgoAPI.BRepAlgoAPI_Fuse, a, b)

    def difference(self, a, b):
        return self._boolean(self.BRepAlgoAPI.BRepAlgoAPI_Cut, a, b)

    def intersection(self, a, b):
        return self._boolean(self.BRepAlgoAPI.BRepAlgoAPI_Common, a, b)

    def section(self, a, b):
        return self._boolean(self.BRepAlgoAPI.BRepAlgoAPI_Section, a, b)

    def compound(self, shapes):
        c = self.TopoDS.TopoDS_Compound()
        builder = self.BRep.BRep_Builder()
        builder.MakeCompound(c)
        for s in shapes:
            builder.Add(c, s)
        return c

    def _map(self, s, kind):
        m = self.TopTools.TopTools_IndexedMapOfShape()
        self.TopExp.TopExp.MapShapes_s(s, kind, m)
        return [m.FindKey(i) for i in range(1, m.Extent() + 1)]

    def subshapes(self, s, kind):
        kinds = self.TopAbs.TopAbs_ShapeEnum
        return self._map(s, {"vertex": kinds.TopAbs_VERTEX, "edge": kinds.TopAbs_EDGE, "face": kinds.TopAbs_FACE,
                             "solid": kinds.TopAbs_SOLID}[kind])

    def _props(self, s, dimension):
        p = self.GProp.GProp_GProps()
        g = self.BRepGProp.BRepGProp
        {1: g.LinearProperties_s, 2: g.SurfaceProperties_s, 3: g.VolumeProperties_s}[dimension](s, p)
        return p

    def sub_centroid(self, sub, kind):
        if kind == "vertex":
            p = self.BRep.BRep_Tool.Pnt_s(self.TopoDS.TopoDS.Vertex(sub))
        else:
            p = self._props(sub, {"edge": 1, "face": 2, "solid": 3}[kind]).CentreOfMass()
        return (p.X(), p.Y(), p.Z())

    def sub_measure(self, sub, kind):
        return self._props(sub, {"edge": 1, "face": 2, "solid": 3}[kind]).Mass()

    def fillet(self, s, pairs):
        maker = self.BRepFilletAPI.BRepFilletAPI_MakeFillet(s)
        edges = self.subshapes(s, "edge")
        for sel, r in pairs:
            maker.Add(r, self.TopoDS.TopoDS.Edge(edges[sel.index]))
        maker.Build()
        if not maker.IsDone():
            raise RuntimeError("fillet failed")
        return maker.Shape()

    def vertex(self, p):
        return self.BRepBuilderAPI.BRepBuilderAPI_MakeVertex(self.pnt(p)).Vertex()

    def segment(self, p, q):
        return self.BRepBuilderAPI.BRepBuilderAPI_MakeEdge(self.pnt(p), self.pnt(q)).Edge()

    def offset_wire(self, s, distance):
        wire = self._outer_wire(self.face(s)) if isinstance(s, tuple) else s
        maker = self.BRepOffsetAPI.BRepOffsetAPI_MakeOffset(wire, self.GeomAbs.GeomAbs_Arc)
        maker.Perform(distance)
        return maker.Shape()

    def _outer_wire(self, face):
        from OCP.BRepTools import BRepTools
        return BRepTools.OuterWire_s(self.TopoDS.TopoDS.Face(face))

    def curve_edge(self, d):
        if d[0] == "bspline":
            from OCP.Geom import Geom_BSplineCurve
            from OCP.TColgp import TColgp_Array1OfPnt
            from OCP.TColStd import TColStd_Array1OfInteger, TColStd_Array1OfReal
            _, poles, knots, multiplicities, weights, degree = d
            P, W = TColgp_Array1OfPnt(1, len(poles)), TColStd_Array1OfReal(1, len(poles))
            K, M = TColStd_Array1OfReal(1, len(knots)), TColStd_Array1OfInteger(1, len(knots))
            for i, (p, w) in enumerate(zip(poles, weights), 1):
                P.SetValue(i, self.pnt(p))
                W.SetValue(i, w)
            for i, (k, m) in enumerate(zip(knots, multiplicities), 1):
                K.SetValue(i, k)
                M.SetValue(i, m)
            return self.BRepBuilderAPI.BRepBuilderAPI_MakeEdge(Geom_BSplineCurve(P, W, K, M, degree)).Edge()
        if d[0] == "bezier":
            from OCP.Geom import Geom_BezierCurve
            from OCP.TColgp import TColgp_Array1OfPnt
            poles = TColgp_Array1OfPnt(1, len(d[1]))
            for i, p in enumerate(d[1]):
                poles.SetValue(i + 1, self.pnt(p))
            return self.BRepBuilderAPI.BRepBuilderAPI_MakeEdge(Geom_BezierCurve(poles)).Edge()
        circle = self.gp.gp_Circ(self.gp.gp_Ax2(self.pnt(d[1]), self.gp.gp_Dir(*d[3])), d[2])
        if d[0] == "arc":
            return self.BRepBuilderAPI.BRepBuilderAPI_MakeEdge(circle, d[4], d[5]).Edge()
        return self.BRepBuilderAPI.BRepBuilderAPI_MakeEdge(circle).Edge()

    def wire(self, edges):
        maker = self.BRepBuilderAPI.BRepBuilderAPI_MakeWire()
        for e in edges:
            maker.Add(e)
        return maker.Wire()

    def _as_wire(self, shape):
        if isinstance(shape, tuple) and shape[0] == "polygon":
            maker = self.BRepBuilderAPI.BRepBuilderAPI_MakePolygon()
            closed = shape[1][0] == shape[1][-1]
            for p in shape[1][:-1] if closed else shape[1]:
                maker.Add(self.pnt(p))
            if closed:
                maker.Close()
            return maker.Wire()
        if isinstance(shape, self.TopoDS.TopoDS_Wire):
            return shape
        if shape.ShapeType() == self.TopAbs.TopAbs_ShapeEnum.TopAbs_WIRE:
            return self.TopoDS.TopoDS.Wire(shape)
        return self.BRepBuilderAPI.BRepBuilderAPI_MakeWire(shape).Wire()

    def loft(self, sections, solid, ruled):
        maker = self.BRepOffsetAPI.BRepOffsetAPI_ThruSections(solid, ruled)
        for s in sections:
            maker.AddWire(self._as_wire(s))
        maker.Build()
        return maker.Shape()

    def to_nurbs(self, s):
        return self.BRepBuilderAPI.BRepBuilderAPI_NurbsConvert(s).Shape()

    def extrude_drafted(self, profile, by, angle):
        from OCP.LocOpe import LocOpe_DPrism
        # (Its height runs along the leaning sides: the result's is that times the cosine.)
        length = math.sqrt(sum(c * c for c in by)) / math.cos(math.radians(angle))
        return LocOpe_DPrism(self.TopoDS.TopoDS.Face(profile), length, math.radians(angle)).Shape()

    def fillet_law(self, s, sel, r0, r1):
        maker = self.BRepFilletAPI.BRepFilletAPI_MakeFillet(s)
        maker.Add(r0, r1, self.TopoDS.TopoDS.Edge(self.subshapes(s, "edge")[sel.index]))
        maker.Build()
        if not maker.IsDone():
            raise RuntimeError("fillet failed")
        return maker.Shape()

    def pipe(self, profile, spine, corner=None):
        if corner == "mitred":
            from OCP.BRepBuilderAPI import BRepBuilderAPI_TransitionMode
            from OCP.BRepTools import BRepTools
            maker = self.BRepOffsetAPI.BRepOffsetAPI_MakePipeShell(self._as_wire(spine))
            maker.SetTransitionMode(BRepBuilderAPI_TransitionMode.BRepBuilderAPI_RightCorner)
            maker.Add(BRepTools.OuterWire_s(self.TopoDS.TopoDS.Face(profile)))
            maker.Build()
            maker.MakeSolid()
            return maker.Shape()
        return self.BRepOffsetAPI.BRepOffsetAPI_MakePipe(self._as_wire(spine), profile).Shape()

    def _join(self, join):
        return self.GeomAbs.GeomAbs_Intersection if join == "intersection" else self.GeomAbs.GeomAbs_Arc

    def chamfer(self, s, edges, distance):
        maker = self.BRepFilletAPI.BRepFilletAPI_MakeChamfer(s)
        all_edges = self.subshapes(s, "edge")
        for sel in edges:
            maker.Add(distance, self.TopoDS.TopoDS.Edge(all_edges[sel.index]))
        maker.Build()
        if not maker.IsDone():
            raise RuntimeError("chamfer failed")
        return maker.Shape()

    def remove_faces(self, s, faces):
        from OCP.BRepAlgoAPI import BRepAlgoAPI_Defeaturing
        maker = BRepAlgoAPI_Defeaturing()
        maker.SetShape(s)
        all_faces = self.subshapes(s, "face")
        for f in faces:
            maker.AddFaceToRemove(all_faces[f.index])
        maker.Build()
        if not maker.IsDone():
            raise RuntimeError("defeaturing failed")
        return maker.Shape()

    def thick_solid(self, s, faces, thickness, join):
        if not faces:
            # (OCCT's maker then only offsets: the shell between the two is made here.)
            offset = self.offset_shape(s, thickness, join)
            return self.difference(offset, s) if thickness > 0 else self.difference(s, offset)
        all_faces = self.subshapes(s, "face")
        closing = self.TopTools.TopTools_ListOfShape()
        for f in faces:
            closing.Append(all_faces[f.index])
        maker = self.BRepOffsetAPI.BRepOffsetAPI_MakeThickSolid()
        maker.MakeThickSolidByJoin(s, closing, thickness, 1e-3, self.BRepOffset.BRepOffset_Skin, False, False,
                                   self._join(join))
        maker.Build()
        return maker.Shape()

    def offset_shape(self, s, distance, join):
        maker = self.BRepOffsetAPI.BRepOffsetAPI_MakeOffsetShape()
        maker.PerformByJoin(s, distance, 1e-3, self.BRepOffset.BRepOffset_Skin, False, False, self._join(join))
        return maker.Shape()

    def face(self, profile):
        if not isinstance(profile, tuple):
            wire = self._as_wire(profile)
        elif profile[0] == "polygon":
            w = self.BRepBuilderAPI.BRepBuilderAPI_MakePolygon()
            for p in profile[1][:-1]:
                w.Add(self.pnt(p))
            w.Close()
            wire = w.Wire()
        else:
            wire = self.BRepBuilderAPI.BRepBuilderAPI_MakeWire(self.curve_edge(profile)).Wire()
        return self.BRepBuilderAPI.BRepBuilderAPI_MakeFace(wire).Shape()

    def extrude(self, profile, by):
        profile = self._as_wire(profile) if isinstance(profile, tuple) else profile
        return self.BRepPrimAPI.BRepPrimAPI_MakePrism(profile, self.gp.gp_Vec(*by)).Shape()

    def revolve(self, profile, origin, axis, angle):
        ax = self.gp.gp_Ax1(self.pnt(origin), self.gp.gp_Dir(*axis))
        return self.BRepPrimAPI.BRepPrimAPI_MakeRevol(profile, ax, math.radians(angle)).Shape()

    def volume(self, s):
        return self._props(s, 3).Mass()  # (Of an empty result, 0.)

    def area(self, s):
        return self._props(s, 2).Mass()

    def length(self, s):
        if isinstance(s, tuple):  # A polygon, not yet a shape.
            return sum(math.dist(p, q) for p, q in zip(s[1], s[1][1:]))
        return self._props(s, 1).Mass()

    def centroid(self, s):
        p = self._props(s, 3).CentreOfMass()
        return [p.X(), p.Y(), p.Z()]

    def bounds(self, s):
        b = self.Bnd.Bnd_Box()
        self.BRepBndLib.BRepBndLib.AddOptimal_s(s, b, False, False)
        v = b.Get()
        return [list(v[:3]), list(v[3:])]

    def valid(self, s):
        return bool(self.BRepCheck.BRepCheck_Analyzer(s).IsValid())

    def counts(self, s):
        return {p: len(self.subshapes(s, k)) for k, p in (("vertex", "vertices"), ("edge", "edges"), ("face", "faces"), ("solid", "solids"))}

    def distance(self, a, b):
        return self.BRepExtrema.BRepExtrema_DistShapeShape(a, b).Value()

    def inertia(self, s):
        m = self._props(s, 3).MatrixOfInertia()
        return [m.Value(1, 1), m.Value(2, 2), m.Value(3, 3)]

    def contains(self, s, p):
        c = self.BRepClass3d.BRepClass3d_SolidClassifier(s, self.pnt(p), 1e-7)
        return c.State() == self.TopAbs.TopAbs_State.TopAbs_IN


# --- Brute force -----------------------------------------------------------------------------

class Region:
    """A set of points: a vectorised membership test, and a box holding it."""

    def __init__(self, test, lo, hi, planar=False):
        self.test, self.lo, self.hi, self.planar = test, lo, hi, planar


class Brute(Described):
    SAMPLES = 1 << 22  # Per shift; independent random shifts give an honest error estimate.
    SHIFTS = 8
    LINES = 1 << 16
    GRID = 200

    def __init__(self):
        import numpy
        self.np = numpy

    def vertex(self, p):
        return Region(None, self.np.array(p, float), self.np.array(p, float))

    def segment(self, p, q):
        raise NotImplementedError("a segment")

    def fillet_law(self, s, sel, r0, r1):
        raise NotImplementedError("fillet needs a `brute:` predicate")

    def pipe(self, profile, spine, corner=None):
        raise NotImplementedError("pipe needs a `brute:` predicate")

    def loft(self, sections, solid, ruled):
        raise NotImplementedError("loft needs a `brute:` predicate")

    def to_nurbs(self, s):
        return s  # (The same set of points.)

    def extrude_drafted(self, profile, by, angle):
        raise NotImplementedError("a drafted extrusion needs a `brute:` predicate")

    def _frame_axes(self, frame):
        np = self.np
        if frame is None:
            return np.zeros(3), np.eye(3)
        origin, z, x = frame
        z = np.array(z)
        if x is None:
            # (The axes around Z matter to partial shapes only: gp_Ax2's choice is OCCT's.)
            helper = np.array([1.0, 0, 0]) if abs(z[0]) < 0.9 else np.array([0, 1.0, 0])
            x = helper - z * (helper @ z)
            x /= np.linalg.norm(x)
        else:
            x = np.array(x)
        return np.array(origin), np.stack([x, np.cross(z, x), z])

    def _placed(self, region, frame):
        """The region, defined in the frame, in world coordinates."""
        np = self.np
        origin, axes = self._frame_axes(frame)
        return self._mapped(region, lambda P: (P - origin) @ axes.T, lambda P: P @ axes + origin)

    def _mapped(self, region, to_local, to_world):
        np = self.np
        lo, hi = region.lo, region.hi
        corners = np.array([[a, b, c] for a in (lo[0], hi[0]) for b in (lo[1], hi[1]) for c in (lo[2], hi[2])])
        world = to_world(corners)
        return Region(lambda P: region.test(to_local(P)), world.min(axis=0), world.max(axis=0))

    def _sector(self, angle):
        np = self.np
        if angle is None or angle >= 360:
            return lambda P: np.ones(len(P), bool)
        return lambda P: np.mod(np.arctan2(P[:, 1], P[:, 0]), 2 * np.pi) <= math.radians(angle)

    def box(self, size, corner, frame):
        np = self.np
        lo, hi = np.array(corner), np.array(corner) + np.array(size)
        return self._placed(Region(lambda P: ((P >= lo) & (P <= hi)).all(axis=1), lo, hi), frame)

    def cylinder(self, r, h, angle, frame):
        return self.cone(r, r, h, angle, frame)

    def cone(self, r1, r2, h, angle, frame):
        np = self.np
        sector = self._sector(angle)
        rmax = max(r1, r2)

        def test(P):
            radius = r1 + (r2 - r1) * P[:, 2] / h
            return (P[:, 2] >= 0) & (P[:, 2] <= h) & (P[:, 0] ** 2 + P[:, 1] ** 2 <= radius ** 2) & sector(P)
        return self._placed(Region(test, np.array([-rmax, -rmax, 0.0]), np.array([rmax, rmax, h])), frame)

    def sphere(self, r, latitudes, angle, frame):
        np = self.np
        sector = self._sector(angle)
        low, high = (math.radians(v) for v in (latitudes or (-90, 90)))

        def test(P):
            d = np.sqrt((P * P).sum(axis=1))
            with np.errstate(invalid="ignore", divide="ignore"):
                latitude = np.arcsin(np.clip(P[:, 2] / d, -1, 1))
            return (d <= r) & (latitude >= low) & (latitude <= high) & sector(P)
        return self._placed(Region(test, np.full(3, -r), np.full(3, r)), frame)

    def torus(self, major, minor, angle, frame):
        np = self.np
        sector = self._sector(angle)
        extent = major + minor

        def test(P):
            ring = np.sqrt(P[:, 0] ** 2 + P[:, 1] ** 2) - major
            return (ring ** 2 + P[:, 2] ** 2 <= minor ** 2) & sector(P)
        return self._placed(Region(test, np.array([-extent, -extent, -minor]), np.array([extent, extent, minor])), frame)

    def wedge(self, size, top):
        np = self.np
        dx, dy, dz = size
        xmin, zmin, xmax, zmax = top

        def test(P):
            t = P[:, 1] / dy
            return ((P[:, 1] >= 0) & (P[:, 1] <= dy) & (P[:, 0] >= xmin * t) & (P[:, 0] <= dx + (xmax - dx) * t)
                    & (P[:, 2] >= zmin * t) & (P[:, 2] <= dz + (zmax - dz) * t))
        return Region(test, np.zeros(3), np.array(size, float))

    def translate(self, s, v):
        np = self.np
        v = np.array(v)
        return self._mapped(s, lambda P: P - v, lambda P: P + v)

    def rotate(self, s, origin, axis, angle):
        np = self.np
        k, o, a = np.array(axis), np.array(origin), math.radians(angle)
        K = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
        R = np.eye(3) + math.sin(a) * K + (1 - math.cos(a)) * K @ K
        return self._mapped(s, lambda P: (P - o) @ R + o, lambda P: (P - o) @ R.T + o)

    def scale(self, s, center, factor):
        np = self.np
        c = np.array(center)
        return self._mapped(s, lambda P: (P - c) / factor + c, lambda P: (P - c) * factor + c)

    def mirror(self, s, origin, normal):
        np = self.np
        o, n = np.array(origin), np.array(normal)
        flip = lambda P: P - 2 * ((P - o) @ n)[:, None] * n
        return self._mapped(s, flip, flip)

    def union(self, a, b):
        np = self.np
        return Region(lambda P: a.test(P) | b.test(P), np.minimum(a.lo, b.lo), np.maximum(a.hi, b.hi), a.planar)

    compound = None  # Set below: a list of shapes is their union.

    def difference(self, a, b):
        return Region(lambda P: a.test(P) & ~b.test(P), a.lo, a.hi, a.planar)

    def intersection(self, a, b):
        np = self.np
        return Region(lambda P: a.test(P) & b.test(P), np.maximum(a.lo, b.lo), np.minimum(a.hi, b.hi), a.planar)

    def section(self, a, b):
        raise NotImplementedError("section")

    def fillet(self, s, pairs):
        raise NotImplementedError("fillet needs a `brute:` predicate")

    def face(self, profile):
        """A planar region at z = 0 (or the polygon's z), tested on x and y alone."""
        np = self.np
        if profile[0] not in ("circle", "polygon") or (profile[0] == "circle" and abs(profile[3][2]) < 1 - 1e-12):
            raise NotImplementedError("a face needing a `brute:` predicate")
        if profile[0] == "circle":
            c, r = np.array(profile[1]), profile[2]
            return Region(lambda P: (P[:, 0] - c[0]) ** 2 + (P[:, 1] - c[1]) ** 2 <= r * r,
                          c - np.array([r, r, 0]), c + np.array([r, r, 0]), True)
        corners = np.array(profile[1], float)

        def test(P):
            inside = np.zeros(len(P), bool)
            for a, b in zip(corners, corners[1:]):
                crosses = (a[1] > P[:, 1]) != (b[1] > P[:, 1])
                with np.errstate(divide="ignore", invalid="ignore"):
                    x = a[0] + (P[:, 1] - a[1]) * (b[0] - a[0]) / (b[1] - a[1])
                inside ^= crosses & (P[:, 0] < x)
            return inside
        return Region(test, corners.min(axis=0), corners.max(axis=0), True)

    def extrude(self, profile, by):
        np = self.np
        if isinstance(profile, tuple) or not by[2]:
            raise NotImplementedError("an extrusion needing a `brute:` predicate")
        v, z0 = np.array(by), profile.lo[2]

        def test(P):
            t = (P[:, 2] - z0) / v[2]
            return profile.test(P - t[:, None] * v) & (t >= 0) & (t <= 1)
        return Region(test, np.minimum(profile.lo, profile.lo + v), np.maximum(profile.hi, profile.hi + v))

    def revolve(self, profile, origin, axis, angle):
        raise NotImplementedError("revolve needs a `brute:` predicate")

    def subshapes(self, s, kind):
        raise NotImplementedError("sub-shapes")

    def _points(self, region, shift):
        """A quasi-random sequence (the R3 additive recurrence) over the region's box, in chunks."""
        np = self.np
        g = 1.2207440846057596  # The real root of x^4 = x + 1.
        alpha = np.array([1 / g, 1 / g ** 2, 1 / g ** 3])
        offset = np.random.default_rng(1000 + shift).random(3)  # (A Cranley-Patterson rotation.)
        chunk = 1 << 20
        for start in range(0, self.SAMPLES, chunk):
            n = (start + np.arange(1, chunk + 1, dtype=np.float64))[:, None]
            yield region.lo + (region.hi - region.lo) * np.modf(offset + n * alpha)[0]

    def _volume(self, region):
        np = self.np
        box = float(np.prod(region.hi - region.lo))
        if box <= 0 or not (region.hi > region.lo).all():
            return 0.0, 0.0, np.zeros(3)
        estimates, moments = [], []
        for shift in range(self.SHIFTS):
            inside, moment = 0, np.zeros(3)
            for P in self._points(region, shift):
                t = region.test(P)
                inside += int(t.sum())
                moment += P[t].sum(axis=0)
            estimates.append(box * inside / self.SAMPLES)
            moments.append(moment / max(inside, 1))
        mean = float(np.mean(estimates))
        self.centroid_error = 4 * float(np.std(moments, axis=0).max()) / math.sqrt(self.SHIFTS)
        return mean, float(np.std(estimates) / math.sqrt(self.SHIFTS)), np.mean(moments, axis=0)

    def _planar_area(self, s):
        np = self.np
        if not (s.hi[:2] > s.lo[:2]).all():
            return {"value": 0.0, "error": 0.0}
        estimates = []
        box = float(np.prod((s.hi - s.lo)[:2]))
        flat = Region(s.test, s.lo.copy(), s.hi.copy())
        flat.lo[2], flat.hi[2] = s.lo[2], s.lo[2]
        for shift in range(self.SHIFTS):
            inside = sum(int(s.test(P).sum()) for P in self._points(flat, shift))
            estimates.append(box * inside / self.SAMPLES)
        mean = float(np.mean(estimates))
        return {"value": mean, "error": 4 * float(np.std(estimates)) / math.sqrt(self.SHIFTS)}

    def volume(self, s):
        v, e, _ = self._volume(s)
        return {"value": v, "error": 4 * e}

    def centroid(self, s):
        centre = self._volume(s)[2]
        return {"value": [float(c) for c in centre], "error": self.centroid_error}

    def area(self, s):
        """Cauchy-Crofton: isotropic lines through a sphere of radius R cross a surface of area A
        on average A / (2 pi R^2) times."""
        np = self.np
        if s.planar:
            return self._planar_area(s)
        rng = np.random.default_rng(1)
        centre = (s.lo + s.hi) / 2
        R = float(np.linalg.norm(s.hi - s.lo)) / 2 * 1.001
        steps = 2048
        t = (np.arange(steps) + 0.5) / steps * 2 - 1
        batches = []
        for _ in range(16):
            n = self.LINES // 16
            d = rng.normal(size=(n, 3))
            d /= np.linalg.norm(d, axis=1, keepdims=True)
            q = rng.normal(size=(n, 3))
            q -= d * (q * d).sum(axis=1, keepdims=True)
            q /= np.linalg.norm(q, axis=1, keepdims=True)
            rho = R * np.sqrt(rng.random(n))
            half = np.sqrt(R * R - rho * rho)
            crossings = np.zeros(n)
            for start in range(0, n, 512):
                sl = slice(start, start + 512)
                P = (centre + q[sl, None, :] * rho[sl, None, None]
                     + d[sl, None, :] * half[sl, None, None] * t[None, :, None])
                inside = s.test(P.reshape(-1, 3)).reshape(-1, steps)
                crossings[sl] = (inside[:, 1:] != inside[:, :-1]).sum(axis=1) + inside[:, 0] + inside[:, -1]
            batches.append(crossings.mean() * 2 * math.pi * R * R)
        mean = float(np.mean(batches))
        return {"value": mean, "error": 3 * float(np.std(batches) / 4) + mean * 2e-3}

    def _grid(self, s):
        np = self.np
        pad = (s.hi - s.lo) * 0.02 + 1e-9
        axes = [np.linspace(s.lo[i] - pad[i], s.hi[i] + pad[i], self.GRID) for i in range(3)]
        X, Y, Z = np.meshgrid(*axes, indexing="ij")
        return s.test(np.stack([X.ravel(), Y.ravel(), Z.ravel()], axis=1)).reshape(X.shape), axes

    def counts(self, s):
        from scipy import ndimage
        grid, _ = self._grid(s)
        solids = ndimage.label(grid)[1]
        cavities = ndimage.label(~grid)[1] - 1
        return {"solids": solids, "cavities": cavities}

    def bounds(self, s):
        np = self.np
        grid, axes = self._grid(s)
        lo, hi = [], []
        for i in range(3):
            filled = np.nonzero(grid.any(axis=tuple(j for j in range(3) if j != i)))[0]
            lo.append(float(axes[i][filled[0]]))
            hi.append(float(axes[i][filled[-1]]))
        step = [float(a[1] - a[0]) for a in axes]
        return {"value": [lo, hi], "error": max(step) * 1.5}

    def valid(self, s):
        raise NotImplementedError("validity")

    def length(self, s):
        raise NotImplementedError("length")

    def _surface_points(self, s):
        if not isinstance(s, Region):
            raise NotImplementedError("the distance to a curve")
        return self._region_surface_points(s)

    def _region_surface_points(self, s):
        """The grid's cells of a region next to cells outside it, as points."""
        from scipy import ndimage
        np = self.np
        if s.test is None:
            return s.lo[None, :], 0.0
        grid, axes = self._grid(s)
        rim = grid & ~ndimage.binary_erosion(grid)
        i, j, k = np.nonzero(rim)
        step = max(float(a[1] - a[0]) for a in axes)
        return np.stack([axes[0][i], axes[1][j], axes[2][k]], axis=1), step

    def distance(self, a, b):
        from scipy.spatial import cKDTree
        (P, h1), (Q, h2) = self._surface_points(a), self._surface_points(b)
        if not len(P) or not len(Q):
            raise NotImplementedError("a flat shape")
        d = float(cKDTree(Q).query(P)[0].min())
        return {"value": d, "error": 2 * (h1 + h2)}

    def inertia(self, s):
        np = self.np
        box = float(np.prod(s.hi - s.lo))
        rows = []
        for shift in range(self.SHIFTS):
            n, first, second = 0, np.zeros(3), np.zeros(3)
            for P in self._points(s, shift):
                Q = P[s.test(P)]
                n += len(Q)
                first += Q.sum(axis=0)
                second += (Q * Q).sum(axis=0)
            mean, square, mass = first / n, second / n, box * n / self.SAMPLES
            variance = square - mean ** 2
            rows.append([mass * (variance[1] + variance[2]), mass * (variance[0] + variance[2]), mass * (variance[0] + variance[1])])
        rows = np.array(rows)
        return {"value": [float(c) for c in rows.mean(axis=0)], "error": 4 * float(rows.std(axis=0).max()) / math.sqrt(self.SHIFTS)}

    def contains(self, s, p):
        return bool(s.test(self.np.array([p]))[0])


def brute_compound(self, shapes):
    region = shapes[0]
    for s in shapes[1:]:
        region = self.union(region, s)
    return region


Brute.compound = brute_compound


class BruteInterpreter(Interpreter):
    """Takes a name's predicate from the case's `brute:` map, where it has one."""

    def bind(self, name, value):
        override = (self.case.get("brute") or {}).get(name)
        if override is None:
            try:
                self.names[name] = self.op(value)
            except (NotImplementedError, AttributeError, TypeError, KeyError):
                self.names[name] = None  # (Whatever depends on it goes unmeasured here.)
            return
        np = self.kernel.np
        known = {k: v for k, v in self.names.items() if isinstance(v, Region)}
        if not known and not (isinstance(override, dict) and "box" in override):
            raise ValueError(f"{name}: a predicate with nothing to bound it needs a box")
        callable_names = {k: (lambda r: lambda x, y, z: r.test(np.stack([x, y, z], axis=1)))(v) for k, v in known.items()}
        constants = {**Values.FUNCTIONS, **self.values.params, "np": np}
        code = compile(override if isinstance(override, str) else override["test"], name, "eval")

        def test(P):
            return eval(code, {"__builtins__": {}}, {**constants, **callable_names, "x": P[:, 0], "y": P[:, 1], "z": P[:, 2]})
        if isinstance(override, dict) and "box" in override:
            lo, hi = (np.array(self.values.point(p)) for p in override["box"])
        else:
            lo = np.min([r.lo for r in known.values()], axis=0)
            hi = np.max([r.hi for r in known.values()], axis=0)
        self.names[name] = Region(test, lo, hi)

    def op_select(self, a):
        return None  # The predicate of what is done with a selection is the case's own.

    def op_fillet(self, a):
        raise NotImplementedError("fillet needs a `brute:` predicate")

    def op_thick_solid(self, a):
        raise NotImplementedError("thick_solid needs a `brute:` predicate")

    def op_offset_shape(self, a):
        raise NotImplementedError("offset_shape needs a `brute:` predicate")

    def op_chamfer(self, a):
        raise NotImplementedError("chamfer needs a `brute:` predicate")

    def op_remove_faces(self, a):
        raise NotImplementedError("remove_faces needs a `brute:` predicate")


# --- Running one case, in a process of its own ------------------------------------------------

def run_one(leg, path, name):
    case = next(c for c in load_cases(path) if c["case"] == name)
    kernel = {"smlib": Smlib, "occt": Occt, "brute": Brute}[leg]()
    interpreter = (BruteInterpreter if leg == "brute" else Interpreter)(kernel, case)
    print("\n@@RESULT@@" + json.dumps(interpreter.run(), default=float), flush=True)


def spawn(python, leg, path, name, timeout):
    started = time.time()
    # (From a directory of the tool's own: nothing beside the case files is imported.)
    command = [python, os.path.abspath(__file__), "--one", leg, path, name]
    try:
        done = subprocess.run(command, capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        return {"error": f"timeout after {timeout} s"}
    if "@@RESULT@@" in done.stdout:
        return {"checks": json.loads(done.stdout.split("@@RESULT@@")[1]), "seconds": round(time.time() - started, 1)}
    tail = done.stderr.strip().splitlines()
    return {"error": tail[-1][:200] if tail else f"exit code {done.returncode}"}


def expected_of(check, values):
    for key in ("eq", "lt", "le", "gt", "ge"):
        if key in check:
            return key, check[key]
    return None, None


def agrees(check, record, values, leg):
    """Whether a leg's measure fits the check; None when the leg cannot say."""
    if record.get("unavailable") or record.get("value") is None:
        return None
    measure = record["measure"]
    value, error = record["value"], 0.0
    if isinstance(value, dict) and "value" in value:
        value, error = value["value"], value["error"]
    relation, expected = expected_of(check, values)
    if measure == "valid":
        return bool(value)
    if measure == "count":
        if not isinstance(expected, dict):
            return None
        shared = [k for k in expected if k in value]
        if not shared:
            return None
        return all(value[k] == expected[k] for k in shared) if relation == "eq" else None
    if measure == "contains":
        return all(v == check["inside"] for v in value)
    if expected is None:
        return None
    if measure == "inertia":
        wanted = [values.number(v) for v in expected]
        scale = max(abs(v) for v in wanted)
        tolerance = max(error, check.get("rel", 1e-6) * scale, SMLIB_DEFAULT_TOLERANCE * scale if leg == "smlib" else 0)
        return all(abs(a - b) <= tolerance for a, b in zip(value, wanted))
    if measure in ("bounds", "centroid"):
        flat = lambda v: [c for row in v for c in (row if isinstance(row, (list, tuple)) else [row])]
        tolerance = max(error, check.get("abs", 1e-6), 1e-4 if leg == "smlib" else 0)
        want = flat([values.point(r) if isinstance(r, list) else values.number(r) for r in expected])
        return all(abs(a - b) <= tolerance for a, b in zip(flat(value), want))
    expected = values.number(expected)
    if relation != "eq":
        return {"lt": value < expected, "le": value <= expected, "gt": value > expected, "ge": value >= expected}[relation]
    rel, absolute = check.get("rel"), check.get("abs", 0.0)
    tolerance = (1e-6 if rel is None and not absolute else (rel or 0.0)) * abs(expected) + absolute
    if leg == "smlib":
        tolerance = max(tolerance, SMLIB_DEFAULT_TOLERANCE * abs(expected))
    return abs(value - expected) <= tolerance + error


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "--one":
        run_one(*sys.argv[2:5])
        return
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("files", nargs="+")
    parser.add_argument("--smlib-python", required=True)
    parser.add_argument("--occt-python", required=True)
    parser.add_argument("--only", default="")
    parser.add_argument("--legs", default="smlib,occt,brute")
    parser.add_argument("-j", type=int, default=max(2, (os.cpu_count() or 4) // 2))
    parser.add_argument("--timeout", type=float, default=900)
    parser.add_argument("--json")
    parser.add_argument("--quiet", action="store_true", help="only the cases with a disagreement")
    args = parser.parse_args()

    legs = args.legs.split(",")
    python = {"smlib": args.smlib_python, "occt": args.occt_python, "brute": args.occt_python}
    jobs = []
    for path in args.files:
        for case in load_cases(path):
            if args.only in case["case"]:
                jobs.append((path, case))
    results = {}
    with futures.ThreadPoolExecutor(args.j) as pool:
        running = {pool.submit(spawn, python[leg], leg, path, case["case"], args.timeout): (case["case"], leg)
                   for path, case in jobs for leg in legs}
        for f in futures.as_completed(running):
            name, leg = running[f]
            results.setdefault(name, {})[leg] = f.result()

    verdicts = {"agree": 0, "disagree": 0}
    for path, case in jobs:
        name = case["case"]
        values = Values(case.get("params"))
        lines, bad = [], False
        for leg in legs:
            outcome = results[name][leg]
            if "error" in outcome:
                lines.append(f"    {leg}: ERROR {outcome['error']}")
                if leg != "smlib":
                    bad = True
                continue
            for record in outcome["checks"]:
                fits = agrees(record, record, values, leg)
                shown = record.get("value", record.get("unavailable"))
                relation, expected = expected_of(record, values)
                mark = {True: "ok ", False: "DIFF", None: "-  "}[fits]
                if fits is False:
                    bad = True
                lines.append(f"    {leg:6}{mark} {record['measure']:9} {shown}  (case: {relation} {expected})")
        verdicts["disagree" if bad else "agree"] += 1
        if bad or not args.quiet:
            print(("! " if bad else "  ") + name)
            print("\n".join(lines))
    print(f"{len(jobs)} cases: {verdicts['agree']} with no disagreement, {verdicts['disagree']} with one")
    if args.json:
        with open(args.json, "w") as f:
            json.dump(results, f, indent=1, default=float)


if __name__ == "__main__":
    main()
