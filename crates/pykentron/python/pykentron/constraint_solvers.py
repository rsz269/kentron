"""CadQuery's assembly and sketch constraint solvers, solved by kentron.

CadQuery states its constraint problems in `cadquery.occ_impl.solver` (assemblies, solved by
IPOPT through casadi) and `cadquery.occ_impl.sketch_solver` (sketches, solved by nlopt's
SLSQP). The shim's copies of those modules end by handing their solver classes to
`assembly_solver` and `sketch_solver` here, which put in their place classes taking the same
arguments and giving the same kind of result, the problems solved by kentron (least squares,
with exact derivatives) without casadi, nlopt or numpy.

`PYKENTRON_SOLVER=cadquery` keeps CadQuery's own solvers (casadi and nlopt then needed), to
compare the two on the same geometry.
"""

import os
from math import pi


def _kentron_solves() -> bool:
    return os.environ.get("PYKENTRON_SOLVER", "kentron").strip().lower() != "cadquery"


def _xyz(p):
    return [p.X(), p.Y(), p.Z()]


def _assembly_numbers(ks, markers, kind, val):
    """One of CadQuery's simple assembly constraints (`ConstraintSpec.toPODs`) as the native
    solver takes it: (kind, parts, numbers); `None` for "Fixed" (the part is locked)."""
    from OCP.gp import gp_Quaternion, gp_Extrinsic_XYZ

    parts = list(ks)
    if kind == "Point":
        m1, m2 = markers
        return kind, parts, _xyz(m1) + _xyz(m2) + [0.0 if val is None else float(val)]
    if kind == "Axis":
        m1, m2 = markers
        val = pi if val is None else val
        mode = 0 if val == 0 else 1 if val == pi else 2
        return kind, parts, _xyz(m1) + _xyz(m2) + [float(mode), float(val)]
    if kind == "PointInPlane":
        m1, m2 = markers
        axis = m2.Axis()
        return kind, parts, _xyz(m1) + _xyz(axis.Location()) + _xyz(axis.Direction()) + [0.0 if val is None else float(val)]
    if kind == "PointOnLine":
        m1, m2 = markers
        return kind, parts, _xyz(m1) + _xyz(m2.Location()) + _xyz(m2.Direction()) + [0.0 if val is None else float(val)]
    if kind == "Fixed":
        return None
    if kind == "FixedPoint":
        (m1,) = markers
        return kind, parts, _xyz(m1) + [float(v) for v in val]
    if kind == "FixedAxis":
        (m1,) = markers
        return kind, parts, _xyz(m1) + [float(v) for v in val]
    if kind == "FixedRotation":
        q = gp_Quaternion()
        q.SetEulerAngles(gp_Extrinsic_XYZ, *val)
        return kind, parts, [q.W(), q.X(), q.Y(), q.Z()]
    raise ValueError(f"Unknown constraint kind {kind}")


def assembly_solver(original):
    """The assembly solver class to use: kentron's, or CadQuery's `original`."""
    if not _kentron_solves():
        return original

    class ConstraintSolver:
        """CadQuery's assembly constraint problem (`cadquery.occ_impl.solver`), solved by
        kentron: the parts' placements making the constraints' costs least."""

        def __init__(self, entities, constraints, locked=[], scale=1):
            self.entities = list(entities)
            self.constraints = constraints
            self.locked = list(locked)
            self.scale = scale
            self.ne = len(self.entities)
            self.nc = len(self.constraints)

        _locToDOF6 = staticmethod(original._locToDOF6)

        @staticmethod
        def _build_transform(t, r):
            from OCP.gp import gp_Trsf, gp_Quaternion, gp_Vec

            rv = gp_Trsf()
            a, b, c = r
            m = a ** 2 + b ** 2 + c ** 2
            rv.SetRotation(gp_Quaternion(2 * a / (m + 1), 2 * b / (m + 1), 2 * c / (m + 1), (1 - m) / (m + 1)))
            rv.SetTranslationPart(gp_Vec(*t))
            return rv

        def solve(self, verbosity: int = 0):
            from . import _kentron
            from cadquery.occ_impl.geom import Location

            starts = [self._locToDOF6(loc) for loc in self.entities]
            locked = [i in self.locked for i in range(self.ne)]
            numbers = []
            for ks, (markers, kind, val) in self.constraints:
                n = _assembly_numbers(ks, markers, kind, val)
                if n is not None:
                    numbers.append(n)
            placements, (cost, steps, evaluations, converged), optimal, regularizations = _kentron.solve_assembly(starts, locked, numbers, float(self.scale))
            locs = [Location(self._build_transform(t, r)) for t, r in placements]
            result = {
                "success": converged,
                "return_status": "Solve_Succeeded" if optimal else "Solved_To_Acceptable_Level" if converged else "Maximum_Iterations_Exceeded",
                "iter_count": steps,
                # (Unconstrained least squares: no primal infeasibility.)
                "iterations": {"inf_pr": [0.0], "obj": [cost]},
                "cost": cost,
                "evaluations": evaluations,
                # (Each iteration's Hessian regularization, IPOPT's "rg" column.)
                "regularizations": regularizations,
                "solver": "kentron",
            }
            return locs, result

    return ConstraintSolver


def _sketch_params(kind, val):
    if kind in ("Fixed", "Coincident"):
        return []
    if kind == "FixedPoint":
        return [None if val is None else float(val)]
    if kind == "Distance":
        v1, v2, d = val
        return [None if v1 is None else float(v1), None if v2 is None else float(v2), float(d)]
    if kind == "Orientation":
        return [float(v) for v in val]
    return [float(val)]


def _arc_point(x, val):
    """An arc's point at fraction `val` of its sweep (`None`: its centre), as a tuple."""
    from math import sin, cos

    if val is None:
        return (x[0], x[1])
    a = x[3] + val * x[4]
    return (x[0] + x[2] * sin(a), x[1] + x[2] * cos(a))


def _arc_first(x):
    return _arc_point(x, 0.0)


def _arc_last(x):
    return _arc_point(x, 1.0)


def sketch_solver(original, arc_first, arc_last, arc_point):
    """The sketch solver class and the arc point functions `Sketch.solve` builds edges with:
    kentron's (points as tuples, without numpy), or CadQuery's `original` and its own."""
    if not _kentron_solves():
        return original, arc_first, arc_last, arc_point

    class SketchConstraintSolver:
        """CadQuery's sketch constraint problem (`cadquery.occ_impl.sketch_solver`), solved by
        kentron: the lines' and arcs' variables making the constraints' costs least."""

        def __init__(self, entities, constraints, geoms):
            self.entities = list(entities)
            self.constraints = list(constraints)
            self.geoms = list(geoms)
            self.ne = len(self.entities)
            self.nc = len(self.constraints)
            self.ixs = [0]
            for e in self.entities:
                self.ixs.append(self.ixs[-1] + len(e))

        def solve(self):
            from . import _kentron

            entities = [(g, [float(v) for v in e]) for g, e in zip(self.geoms, self.entities)]
            constraints = [(kind, e1, e2, _sketch_params(kind, val)) for (e1, e2), kind, val in self.constraints]
            x, (cost, steps, evaluations, converged) = _kentron.solve_sketch(entities, constraints)
            status = {
                "entities": self.entities,
                "cost": cost,
                "iters": evaluations,
                # (nlopt's codes: 4 XTOL_REACHED, where CadQuery's solves end; 5 MAXEVAL_REACHED.)
                "status": 4 if converged else 5,
                "solver": "kentron",
            }
            return [tuple(v) for v in x], status

    return SketchConstraintSolver, _arc_first, _arc_last, _arc_point
