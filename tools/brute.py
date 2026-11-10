#!/usr/bin/env python3
"""Brute-force references for checking a kernel's results (kentron's or OCCT's).

When kentron and CadQuery disagree, neither is the reference: these compute the answer by
slow, simple means from the operands, so that either result can be judged. All read BRep
files and use OCCT (through OCP, as in CadQuery's environment) only for what is safe to
trust: point classification against a solid's faces and fine meshing. Run with the corpus
venv's Python (it has OCP, numpy and scipy).

Library (import tools/brute.py, or run it):

- `load(path)`, `triangles(shape, deflection)`: a shape and its mesh, face placements
  applied, as an (n, 3, 3) array.
- `Classifier(shape).states(points)`: 'I', 'O' or 'B' for each point.
- `voxelize(shape, lo, hi, h)`: an inside mask on a grid of cubes of side `h`, by the
  parity of crossings along z of a fine mesh (vectorised; millions of voxels in a minute,
  where classifying each takes hours).
- `opening(mask, r, h)`, `closing(mask, r, h)`: the union of balls of radius `r` inside the
  solid (rounding its convex edges, as a fillet does) and its complement's (concave ones).
- `mesh_volume(shape, deflection)`: the volume and area of a fine mesh (thin solids too).
- `monte_carlo_volume(shape, n)`: the volume by random points in the box, with its error.
- `boundary_area(a, b, state, deflection)`: the area of a's boundary lying `state` ('I' or
  'O') with respect to b, from a's mesh (triangles classified at their centroids). The
  area of a union is boundary_area(a, b, 'O') + boundary_area(b, a, 'O'); of an
  intersection, both with 'I'; of a difference a - b, boundary_area(a, b, 'O') +
  boundary_area(b, a, 'I').
- `faces_against(result, operands, op)`: each face of a boolean's result, sampled on a
  lattice of its parameters and classified against the operands: faces where the result
  should not have boundary are flagged (a union's faces inside the other operand, etc.).
- `inside_both(a, b, centre, radius, n)`: points near `centre` inside both solids (do
  shapes said to touch overlap?).
- `polytope_volume(halfspaces)`: the exact volume of {x : A x + b <= 0}.
- `wire_offset(wire, distance, n)`: a planar wire's offset with round joins (and ends), as
  the boundary of its points within `distance` (shapely's buffer of the wire sampled `n`
  times an edge): around an open wire, outwards (inwards for a negative distance) from the
  region a closed one bounds. Its length and the number of loops.

CLI:

    brute.py volume SHAPE.brep [N]                 Monte Carlo volume
    brute.py meshed SHAPE.brep [DEFLECTION]        meshed volume and area, finer and finer
    brute.py area union|common|cut A.brep B.brep   boolean's area from the operands' meshes
    brute.py faces union|common|cut R.brep A.brep B.brep   misplaced faces of a result
        (in a union, faces inside an operand or off both; and so on)
    brute.py wire WIRE.brep FRACTION              a wire's offset (round joins) by FRACTION of
        its length (as OCCT's offset tests take it), both ways for a closed wire: length, loops
    brute.py round open|close PRE.brep R.brep... --box=x0,y0,z0,x1,y1,z1 --r=R [--h=H] [--compare=BOX]
        a fillet's rolling-ball reference: the voxels where each result differs from the
        opening (convex edges) or closing (concave) of PRE, the solid before the fillet.
        Only meaningful where every edge compared near is filleted with radius R: limit the
        comparison with --compare=x0,y0,z0,x1,y1,z1 (a box within --box). The box should
        reach R + 2H past it (its sides are a false boundary; voxels within R + H of them
        are not compared).
"""
import sys

import numpy as np


def load(path):
    from OCP.BRep import BRep_Builder
    from OCP.BRepTools import BRepTools
    from OCP.TopoDS import TopoDS_Shape

    shape = TopoDS_Shape()
    if not BRepTools.Read_s(shape, str(path), BRep_Builder()):
        raise ValueError(f"cannot read {path}")
    return shape


def faces(shape):
    from OCP.TopAbs import TopAbs_FACE
    from OCP.TopExp import TopExp_Explorer
    from OCP.TopoDS import TopoDS

    out, e = [], TopExp_Explorer(shape, TopAbs_FACE)
    while e.More():
        out.append(TopoDS.Face_s(e.Current()))
        e.Next()
    return out


def triangles(shape, deflection=1e-3, angular=0.1):
    """The shape's mesh as an (n, 3, 3) array of triangles, faces' placements applied."""
    from OCP.BRep import BRep_Tool
    from OCP.BRepMesh import BRepMesh_IncrementalMesh
    from OCP.TopAbs import TopAbs_REVERSED
    from OCP.TopLoc import TopLoc_Location

    BRepMesh_IncrementalMesh(shape, deflection, False, angular, True)
    out = []
    for f in faces(shape):
        loc = TopLoc_Location()
        tri = BRep_Tool.Triangulation_s(f, loc)
        if tri is None:
            raise ValueError("a face was not meshed")
        t = loc.Transformation()
        nodes = [tri.Node(i).Transformed(t) for i in range(1, tri.NbNodes() + 1)]
        p = np.array([[q.X(), q.Y(), q.Z()] for q in nodes])
        idx = np.array([tri.Triangle(k).Get() for k in range(1, tri.NbTriangles() + 1)]) - 1
        # (Wound as the face is used: a reversed face's triangles turned, so that they all
        # face out of a solid.)
        if f.Orientation() == TopAbs_REVERSED:
            idx = idx[:, ::-1]
        out.append(p[idx])
    return np.concatenate(out) if out else np.zeros((0, 3, 3))


def bounding_box(shape):
    from OCP.Bnd import Bnd_Box
    from OCP.BRepBndLib import BRepBndLib

    box = Bnd_Box()
    BRepBndLib.Add_s(shape, box)
    v = box.Get()
    return np.array(v[:3]), np.array(v[3:])


class Classifier:
    """Points against a solid: 'I' inside, 'O' outside, 'B' on its boundary (within tol)."""

    def __init__(self, shape, tol=1e-7):
        from OCP.BRepClass3d import BRepClass3d_SolidClassifier

        self.c, self.tol = BRepClass3d_SolidClassifier(shape), tol

    def state(self, p):
        from OCP.gp import gp_Pnt

        self.c.Perform(gp_Pnt(*map(float, p)), self.tol)
        return "IOB"[int(self.c.State())] if int(self.c.State()) < 3 else "?"

    def states(self, points):
        return [self.state(p) for p in points]


def axes(lo, hi, h):
    """Voxel centres along each axis (nudged off round values, so that rays miss mesh edges)."""
    lo, hi = np.asarray(lo, float), np.asarray(hi, float)
    n = np.maximum(np.round((hi - lo) / h).astype(int), 1)
    return [lo[i] + (np.arange(n[i]) + 0.5) * h + 1.2345678e-7 * (i + 1) for i in range(3)]


def voxelize(shape, lo, hi, h, deflection=None):
    """Inside mask on the grid `axes(lo, hi, h)`, by ray parity along z through a mesh."""
    tris = triangles(shape, deflection or min(1e-3, 0.05 * h))
    xs, ys, zs = axes(lo, hi, h)
    hits = {}
    for (x0, y0, z0), (x1, y1, z1), (x2, y2, z2) in tris:
        i0, i1 = np.searchsorted(xs, min(x0, x1, x2)), np.searchsorted(xs, max(x0, x1, x2))
        j0, j1 = np.searchsorted(ys, min(y0, y1, y2)), np.searchsorted(ys, max(y0, y1, y2))
        if i0 >= i1 or j0 >= j1:
            continue
        d = (y1 - y2) * (x0 - x2) + (x2 - x1) * (y0 - y2)
        if d == 0:
            continue
        gx, gy = np.meshgrid(xs[i0:i1], ys[j0:j1], indexing="ij")
        a = ((y1 - y2) * (gx - x2) + (x2 - x1) * (gy - y2)) / d
        b = ((y2 - y0) * (gx - x2) + (x0 - x2) * (gy - y2)) / d
        c = 1 - a - b
        m = (a >= 0) & (b >= 0) & (c >= 0)
        z = a * z0 + b * z1 + c * z2
        for ii, jj, zz in zip(*np.nonzero(m), z[m]):
            hits.setdefault((i0 + ii, j0 + jj), []).append(zz)
    mask = np.zeros((len(xs), len(ys), len(zs)), bool)
    for (i, j), found in hits.items():
        found = np.sort(found)
        # (A ray through an edge two triangles share counts it once.)
        found = found[np.concatenate([[True], np.diff(found) > 1e-9])]
        mask[i, j] = np.searchsorted(found, zs) % 2 == 1
    return mask


def opening(mask, r, h):
    """The union of the balls of radius r within the solid: its convex edges rounded."""
    from scipy.ndimage import distance_transform_edt as edt

    centres = edt(mask, sampling=h) >= r
    return edt(~centres, sampling=h) <= r


def closing(mask, r, h):
    """The complement of the union of balls of radius r outside it: concave edges rounded."""
    return ~opening(~mask, r, h)


def mesh_volume(shape, deflection=1e-4):
    """The volume a fine mesh of the shape encloses (signed tetrahedra from the origin), and
    its area; good for thin solids, where random points would hardly ever land inside. The
    mesh lies within `deflection` of the faces: the error is about that times the area."""
    t = triangles(shape, deflection)
    volume = float(np.einsum("ij,ij->i", t[:, 0], np.cross(t[:, 1], t[:, 2])).sum() / 6.0)
    area = float(0.5 * np.linalg.norm(np.cross(t[:, 1] - t[:, 0], t[:, 2] - t[:, 0]), axis=1).sum())
    return volume, area


def monte_carlo_volume(shape, n=200_000, seed=1):
    lo, hi = bounding_box(shape)
    lo, hi = lo - 0.01 * (hi - lo), hi + 0.01 * (hi - lo)
    rng = np.random.default_rng(seed)
    c = Classifier(shape)
    inside = sum(c.state(lo + rng.random(3) * (hi - lo)) == "I" for _ in range(n))
    box = float(np.prod(hi - lo))
    f = inside / n
    return box * f, box * np.sqrt(f * (1 - f) / n)


def subdivided(tris, size):
    """Triangles split in four until no edge is longer than `size`."""
    out = []
    while len(tris):
        longest = np.max(np.linalg.norm(tris - np.roll(tris, 1, axis=1), axis=2), axis=1)
        out.append(tris[longest <= size])
        t = tris[longest > size]
        a, b, c = t[:, 0], t[:, 1], t[:, 2]
        ab, bc, ca = (a + b) / 2, (b + c) / 2, (c + a) / 2
        tris = np.concatenate([np.stack(x, axis=1) for x in ((a, ab, ca), (ab, b, bc), (ca, bc, c), (ab, bc, ca))]) if len(t) else t
    return np.concatenate(out)


def boundary_area(a, b, state, deflection=1e-4, near=None, size=None):
    """Area of a's boundary lying `state` ('I' or 'O') against b, from a's mesh.

    Triangles near b (within its box, or `near`, a box (lo, hi)) are split to edges no
    longer than `size` (by default a five-hundredth of b's box across) and classified at
    their centroids; the others are outside b. The error is about `size` times the length
    of where b crosses a's boundary."""
    tris = triangles(a, deflection)
    blo, bhi = near if near is not None else bounding_box(b)
    size = size or 2e-3 * float(np.linalg.norm(bhi - blo))
    tlo, thi = tris.min(axis=1), tris.max(axis=1)
    hit = np.all(thi >= blo, axis=1) & np.all(tlo <= bhi, axis=1)
    tris = np.concatenate([tris[~hit], subdivided(tris[hit], size)])
    areas = 0.5 * np.linalg.norm(np.cross(tris[:, 1] - tris[:, 0], tris[:, 2] - tris[:, 0]), axis=1)
    tlo, thi = tris.min(axis=1), tris.max(axis=1)
    candidates = np.all(thi >= blo, axis=1) & np.all(tlo <= bhi, axis=1)
    states = np.full(len(tris), "O")
    c = Classifier(b)
    for k in np.nonzero(candidates)[0]:
        states[k] = c.state(tris[k].mean(axis=0))
    return float(areas[states == state].sum())


def boolean_area(op, a, b, deflection=1e-4):
    if op == "union":
        return boundary_area(a, b, "O", deflection) + boundary_area(b, a, "O", deflection)
    if op == "common":
        return boundary_area(a, b, "I", deflection) + boundary_area(b, a, "I", deflection)
    if op == "cut":
        return boundary_area(a, b, "O", deflection) + boundary_area(b, a, "I", deflection)
    raise ValueError(op)


# Where a boolean's result may have boundary: (state against a, state against b) of its
# faces' points, for each operation.
ALLOWED = {
    "union": {("B", "O"), ("O", "B"), ("B", "B")},
    "common": {("B", "I"), ("I", "B"), ("B", "B")},
    "cut": {("B", "O"), ("I", "B"), ("B", "B")},
}


def faces_against(result, a, b, op, n=8):
    """For each face of `result`: its area, the counts of (state against a, against b) of
    points on an n by n lattice of its parameters within it, and whether a quarter or more
    of them lie where the operation leaves no boundary."""
    from OCP.BRepAdaptor import BRepAdaptor_Surface
    from OCP.BRepGProp import BRepGProp
    from OCP.BRepTopAdaptor import BRepTopAdaptor_FClass2d
    from OCP.GProp import GProp_GProps
    from OCP.gp import gp_Pnt2d

    ca, cb = Classifier(a, 1e-6), Classifier(b, 1e-6)
    out = []
    for f in faces(result):
        g = GProp_GProps()
        BRepGProp.SurfaceProperties_s(f, g)
        s, fc = BRepAdaptor_Surface(f), BRepTopAdaptor_FClass2d(f, 1e-9)
        u0, u1, v0, v1 = s.FirstUParameter(), s.LastUParameter(), s.FirstVParameter(), s.LastVParameter()
        counts = {}
        for i in range(1, n):
            for j in range(1, n):
                uv = gp_Pnt2d(u0 + (u1 - u0) * i / n, v0 + (v1 - v0) * j / n)
                if int(fc.Perform(uv)) != 0:
                    continue
                p = s.Value(uv.X(), uv.Y())
                key = (ca.state((p.X(), p.Y(), p.Z())), cb.state((p.X(), p.Y(), p.Z())))
                counts[key] = counts.get(key, 0) + 1
        # (Misplaced where a quarter of its points are: a point or two classified off a
        # free-form face, by its tolerance, is not.)
        total = sum(counts.values())
        bad = total > 0 and sum(v for k, v in counts.items() if k not in ALLOWED[op]) >= 0.25 * total
        out.append((g.Mass(), counts, bad))
    return out


def inside_both(a, b, centre, radius, n=2000, seed=1):
    """Random points within `radius` of `centre` inside both solids."""
    rng = np.random.default_rng(seed)
    ca, cb = Classifier(a), Classifier(b)
    found = []
    for _ in range(n):
        d = rng.normal(size=3)
        p = np.asarray(centre) + d / np.linalg.norm(d) * radius * rng.random() ** (1 / 3)
        if ca.state(p) == "I" and cb.state(p) == "I":
            found.append(p)
    return found


def polytope_volume(halfspaces):
    """Exact volume of {x : A x + b <= 0}, rows [A | b] (bounded, with an interior)."""
    from scipy.optimize import linprog
    from scipy.spatial import ConvexHull, HalfspaceIntersection

    hs = np.asarray(halfspaces, float)
    norms = np.linalg.norm(hs[:, :3], axis=1)
    # The Chebyshev centre, an interior point.
    res = linprog([0, 0, 0, -1], A_ub=np.hstack([hs[:, :3], norms[:, None]]), b_ub=-hs[:, 3], bounds=[(None, None)] * 3 + [(0, None)])
    return ConvexHull(HalfspaceIntersection(hs, res.x[:3]).intersections).volume


def wire_points(wire, n=2000):
    """The wire's points in order, `n` along each edge, as 2D coordinates in the plane that
    fits them best; whether it is closed; and the plane (a point and two axes: a 3D point p
    is at (p - point) @ axes.T)."""
    from OCP.BRepAdaptor import BRepAdaptor_Curve
    from OCP.BRepTools import BRepTools_WireExplorer
    from OCP.TopAbs import TopAbs_REVERSED, TopAbs_WIRE
    from OCP.TopExp import TopExp_Explorer
    from OCP.TopoDS import TopoDS

    e = TopExp_Explorer(wire, TopAbs_WIRE)
    w = TopoDS.Wire_s(e.Current())
    pts = []
    it = BRepTools_WireExplorer(w)
    while it.More():
        edge = it.Current()
        c = BRepAdaptor_Curve(edge)
        ts = np.linspace(c.FirstParameter(), c.LastParameter(), n)
        if edge.Orientation() == TopAbs_REVERSED:
            ts = ts[::-1]
        q = [c.Value(t) for t in ts]
        pts.extend([(p.X(), p.Y(), p.Z()) for p in q][0 if not pts else 1:])
        it.Next()
    pts = np.array(pts)
    closed = np.linalg.norm(pts[0] - pts[-1]) < 1e-6 * np.ptp(pts, axis=0).max()
    # The plane: the points' best fit.
    c = pts.mean(axis=0)
    _, _, vt = np.linalg.svd(pts - c, full_matrices=False)
    uv = (pts - c) @ vt[:2].T
    return uv, closed, (c, vt[:2])


def wire_offset_loops(wire, distance, n=2000):
    """The loops bounding a planar wire's offset with round joins: each an array of points,
    as 2D coordinates in the plane wire_points gives."""
    from shapely.geometry import LineString, Polygon

    uv, closed, _ = wire_points(wire, n)
    base = Polygon(uv) if closed else LineString(uv)
    if not closed:
        distance = abs(distance)
    region = base.buffer(distance, quad_segs=512)
    polys = list(region.geoms) if hasattr(region, "geoms") else ([region] if not region.is_empty else [])
    return [np.array(r.coords) for p in polys for r in [p.exterior, *p.interiors]]


def wire_offset(wire, distance, n=2000):
    """The length and number of loops of a planar wire's offset with round joins."""
    loops = wire_offset_loops(wire, distance, n)
    return sum(np.linalg.norm(np.diff(loop, axis=0), axis=1).sum() for loop in loops), len(loops)


def _round(argv):
    kind, pre, results = argv[0], argv[1], [a for a in argv[2:] if not a.startswith("--")]
    opts = dict(a[2:].split("=", 1) for a in argv if a.startswith("--") and "=" in a)
    box = np.array([float(x) for x in opts["box"].split(",")])
    r, h = float(opts["r"]), float(opts.get("h", 0.04))
    from scipy.ndimage import distance_transform_edt as edt

    lo, hi = box[:3], box[3:]
    reference = (opening if kind == "open" else closing)(voxelize(load(pre), lo, hi, h), r, h)
    xs, ys, zs = axes(lo, hi, h)
    X, Y, Z = np.meshgrid(xs, ys, zs, indexing="ij")
    # (Away from the box's sides by r and a voxel, and from the reference's boundary by more
    # than a voxel and a half, where voxelization alone differs.)
    m = (X > lo[0] + r + h) & (X < hi[0] - r - h) & (Y > lo[1] + r + h) & (Y < hi[1] - r - h) & (Z > lo[2] + r + h) & (Z < hi[2] - r - h)
    if "compare" in opts:
        c = np.array([float(x) for x in opts["compare"].split(",")])
        m &= (X > c[0]) & (Y > c[1]) & (Z > c[2]) & (X < c[3]) & (Y < c[4]) & (Z < c[5])
    m &= np.where(reference, edt(reference, sampling=h), edt(~reference, sampling=h)) >= 1.5 * h
    print(f"compared {int(m.sum())} voxels of {h}")
    for path in results:
        v = voxelize(load(path), lo, hi, h)
        extra, missing = (v & ~reference & m).sum() * h**3, (~v & reference & m).sum() * h**3
        wrong = (v != reference) & m
        where = ""
        if wrong.any():
            idx = np.argwhere(wrong)
            where = " around " + ", ".join(f"{a[int(k)]:.3f}" for a, k in zip((xs, ys, zs), idx.mean(axis=0)))
        print(f"{path}: wrong {int(wrong.sum())} voxels, material extra {extra:.4g}, missing {missing:.4g}{where}")


def main(argv):
    if not argv or argv[0] in ("-h", "--help"):
        print(__doc__)
        return
    cmd, rest = argv[0], argv[1:]
    if cmd == "meshed":
        for d in (1e-3, 1e-4, float(rest[1]) if len(rest) > 1 else 2e-5):
            v, a = mesh_volume(load(rest[0]), d)
            print(f"deflection {d:g}: volume {v:.6g} area {a:.6g}")
    elif cmd == "volume":
        v, e = monte_carlo_volume(load(rest[0]), int(rest[1]) if len(rest) > 1 else 200_000)
        print(f"volume {v:.6g} ± {e:.2g}")
    elif cmd == "area":
        print(f"area {boolean_area(rest[0], load(rest[1]), load(rest[2])):.6g}")
    elif cmd == "faces":
        op, r, a, b = rest[0], load(rest[1]), load(rest[2]), load(rest[3])
        rows = faces_against(r, a, b, op)
        for i, (area, counts, bad) in enumerate(rows):
            if bad:
                print(f"face {i}: area {area:.5g} states {counts}  <-- misplaced")
        print(f"{sum(b for *_, b in rows)} of {len(rows)} faces misplaced, area {sum(a for a, _, b in rows if b):.5g}")
    elif cmd == "round":
        _round(rest)
    elif cmd == "wire":
        from OCP.GProp import GProp_GProps
        from OCP.BRepGProp import BRepGProp

        w = load(rest[0])
        props = GProp_GProps()
        BRepGProp.LinearProperties_s(w, props)
        d = props.Mass() * float(rest[1])
        for sign in ((1, -1) if wire_points(w, 2)[1] else (1,)):
            length, loops = wire_offset(w, sign * d)
            print(f"distance {sign * d:.6g}: length {length:.6g}, {loops} loop(s)")
    else:
        raise SystemExit(f"unknown command {cmd}")


if __name__ == "__main__":
    main(sys.argv[1:])
