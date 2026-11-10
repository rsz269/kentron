"""Straight docking bars: one per module, glued to the skin where the module goes.

A module's clip has a straight T-slot open at both ends and slides onto the bar from its open end with almost no
force (a printed flexure tongue in the clip roof clicks into a dimple on the bar), instead of being pressed onto a curved rail
(which needed about 2,000 N). The closed end has a stop. Each bar is a chord of the old rail ellipse: its ends sit on
the ellipse and its middle sits a few mm inside it; the foot is cut to the skin so the bar stands on it, and a
conforming flange gives the glue some area.
"""
import functools
import math

import cadquery as cq
import numpy as np

import build
import layout

R = build.R
CLIP_LENGTH = 56          # the clip's length along the bar
TRAVEL = 12               # slide distance each way
BAR_LENGTH = CLIP_LENGTH + 2 * TRAVEL
FLANGE_W, FLANGE_T = 30.0, 1.5   # wide flange: a heavy module tries to peel the bar's inner edge, and the peel force falls as the flange widens
PEG_R, PEG_BARB_R, PEG_BARB_L = 1.5, 2.2, 1.4   # locating peg: shaft radius, barb radius at its shoulder, barb length (the barb snaps through a hole in the skin)
PEG_HOLE_R = 1.6                                 # hole in the skin, 0.1 mm clear of the shaft
PEG_X, PEG_Y = -9.0, 26.0                        # pegs sit near the flange's inner edge (where the peel starts), this far along the bar either side of centre
STOP_H = 2.5              # stop block rises this far above the cap
DIMPLES = (-8.0, 0.0, 8.0)   # detent positions along the bar, mm from the centre


@functools.lru_cache(maxsize=None)
def frame(angle: float) -> dict:
    """Bar frame at rail angle `angle`: centre, direction, outward normal, heading and the heights along the bar."""
    a, b, _ = build.surface(R["position"])
    pt = lambda t: np.array([a * math.cos(t), b * math.sin(t)])
    t0 = math.radians(angle)
    lo, hi = 0.0, 1.0
    for _ in range(60):
        mid = (lo + hi) / 2
        if np.linalg.norm(pt(t0 + mid) - pt(t0 - mid)) < BAR_LENGTH:
            lo = mid
        else:
            hi = mid
    A, B = pt(t0 - lo), pt(t0 + lo)
    M = (A + B) / 2
    d = (B - A) / np.linalg.norm(B - A)
    n = np.array([d[1], -d[0]])
    if n @ M < 0:
        n = -n
    zmax = max(build.z_of(*(M + y * d + x * n)) for y in np.linspace(-BAR_LENGTH / 2, BAR_LENGTH / 2, 17) for x in (-5.0, 0.0, 5.0))
    stem_top = zmax + R["stem_height"]
    return {"M": M, "d": d, "n": n, "heading": math.degrees(math.atan2(n[1], n[0])), "zmax": zmax,
            "stem_top": stem_top, "cap_top": stem_top + R["cap_height"]}


def peg_points(angle: float) -> list[np.ndarray]:
    """Plan positions of the bar's two locating pegs."""
    f = frame(angle)
    return [f["M"] + PEG_X * f["n"] + s * PEG_Y * f["d"] for s in (-1, 1)]


def peg_axis(px: float, py: float):
    """Point on the skin's top surface, the outward unit normal there, and the skin's thickness along that normal."""
    z = build.z_of(px, py)
    eps = 0.5
    gx = (build.z_of(px + eps, py) - build.z_of(px - eps, py)) / (2 * eps)
    gy = (build.z_of(px, py + eps) - build.z_of(px, py - eps)) / (2 * eps)
    n = np.array([-gx, -gy, 1.0])
    n /= np.linalg.norm(n)
    return np.array([px, py, z]), n, build.Y["thickness"] * n[2]


def _place(shape: cq.Workplane, f: dict) -> cq.Workplane:
    return shape.rotate((0, 0, 0), (0, 0, 1), f["heading"]).translate((float(f["M"][0]), float(f["M"][1]), 0))


def _box(w: float, length: float, z0: float, z1: float, cy: float = 0.0) -> cq.Workplane:
    return cq.Workplane("XY").box(w, length, z1 - z0, centered=(True, True, False)).translate((0, cy, z0))


def bar(angle: float, sink: float = 0.0) -> cq.Workplane:
    """One docking bar in place. sink > 0 lets the foot overlap the skin by that much (for the formed model's boolean);
    use 0 for the printed part, whose base is the skin surface."""
    f = frame(angle)
    z_low = f["zmax"] - 10
    stem = _box(R["stem_width"], BAR_LENGTH, z_low, f["stem_top"])
    cap = _box(R["cap_width"], BAR_LENGTH, f["stem_top"], f["cap_top"])
    stop = _box(8, 3, f["stem_top"], f["cap_top"] + STOP_H, cy=BAR_LENGTH / 2 - 1.5)
    body = stem.union(cap).union(stop)
    for y in DIMPLES:   # shallow dimples the clip's printed detent rib clicks into (kept inside the clip's travel: it stops 9 mm out)
        body = body.cut(cq.Workplane("XY").workplane(offset=f["cap_top"] - 0.5).center(0, y).circle(0.9).extrude(1.5))
    body = _place(body, f).cut(build.cone(-sink))
    flange_plan = _place(_box(FLANGE_W, BAR_LENGTH, -100, 100), f)
    slab = build.cone(FLANGE_T).cut(build.cone(-sink))
    out = build.safe_union(body, flange_plan.intersect(slab))
    for px, py in peg_points(angle):   # two pegs along the skin's normal: shaft through the skin, then a barb that snaps under it
        P, n, t_n = peg_axis(px, py)
        shaft = cq.Solid.makeCylinder(PEG_R, 0.6 + t_n, cq.Vector(*(P + 0.6 * n)), cq.Vector(*(-n)))
        barb = cq.Solid.makeCone(PEG_BARB_R, 1.2, PEG_BARB_L, cq.Vector(*(P - t_n * n)), cq.Vector(*(-n)))
        out = build.safe_union(out, cq.Workplane("XY").add(shaft))
        out = build.safe_union(out, cq.Workplane("XY").add(barb))
    return out


def bars_formed():
    """All bars of the active kit, as they sit on the formed collar (feet overlap the skin 0.3 mm), or None."""
    shapes = [bar(a, sink=0.3) for a in layout.angles().values()]
    if not shapes:
        return None
    out = shapes[0]
    for s in shapes[1:]:
        out = out.union(s)
    return out
