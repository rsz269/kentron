"""Parametric festival yoke: sloped shoulder collar with a module rail and an LED lip.
Run: uv run --no-project --with cadquery build.py"""
import math
import tomllib
from pathlib import Path

import cadquery as cq

import sizing

ROOT = Path(__file__).parent
cfg = tomllib.loads((ROOT / "config.toml").read_text())
SIZE, W = sizing.resolve(cfg)
OUT = ROOT / "out" / SIZE
Y, R, G, P, J = cfg["yoke"], cfg["rail"], cfg["glow"], cfg["print"], cfg["joint"]
PAD, D = cfg["pad"], cfg["diffuser"]

AI, BI = W["neck_opening_a"], W["neck_opening_b"]
AO, BO = W["shoulder_half_width"], W["shoulder_half_depth"]
DROP = W["shoulder_drop"]


def surface(u: float) -> tuple[float, float, float]:
    """Ellipse on the collar surface; u=0 neck edge, u=1 outer edge. Returns (a, b, z)."""
    return AI + (AO - AI) * u, BI + (BO - BI) * u, -DROP * u


def z_of(x: float, y: float) -> float:
    """Height of the collar's top surface above the plan point (x, y)."""
    lo, hi = 0.0, 1.0
    for _ in range(40):
        u = (lo + hi) / 2
        a, bb = AI + (AO - AI) * u, BI + (BO - BI) * u
        if (x / a) ** 2 + (y / bb) ** 2 > 1:
            lo = u
        else:
            hi = u
    return -DROP * (lo + hi) / 2


def ring(u: float, inner: float, outer: float, z0: float, h: float) -> cq.Workplane:
    """Elliptical wall centred on surface(u), extending `inner`/`outer` mm either side, from z0 up h."""
    a, b, _ = surface(u)
    big = cq.Workplane("XY").workplane(offset=z0).ellipse(a + outer, b + outer).extrude(h)
    small = cq.Workplane("XY").workplane(offset=z0 - 1).ellipse(a - inner, b - inner).extrude(h + 2)
    return big.cut(small)


def cone(dz: float) -> cq.Workplane:
    """Solid frustum under the collar surface, shifted down by dz."""
    return (
        cq.Workplane("XY").workplane(offset=dz).ellipse(AI, BI)
        .workplane(offset=-DROP).ellipse(AO, BO).loft()
    )


def skin_blank() -> cq.Workplane:
    """Collar skin with no slots, rounding or gap: the bare shell, used to develop the flat patterns."""
    t = Y["thickness"]
    s = cone(0).cut(cone(-t))
    k = (AO - AI) * t / DROP + 1
    return s.cut(cq.Workplane("XY").workplane(offset=-DROP - 5).ellipse(AI + k, BI + k).extrude(DROP + 10))


def skin() -> cq.Workplane:
    """Collar skin with rounded edges. Rounding is attempted at decreasing radii until the finished skin is a valid
    solid (the kernel is fragile on some sizes); with none valid it is left sharp and a warning is printed."""
    for r in (Y["round"], Y["round"] * 0.8, Y["round"] * 0.6, Y["round"] * 0.4):
        s = _skin(r)
        if s.val().isValid() and s.val().Volume() > 0:
            return s
    print("WARN edge rounding skipped for this size (kernel could not produce a valid rounded skin)")
    return _skin(0)


def slot_blocked(deg: float) -> bool:
    """True if a lightening slot at this angle would sit under a docking bar's flange (the bar needs solid skin to glue to)."""
    import dock
    import layout
    th = math.radians(deg)
    a, b, _ = surface(0.42)
    cx, cy = math.cos(th) * a, math.sin(th) * b
    half = max((AO - AI) * 0.30 - Y["slot_width"], 0) / 2
    ux, uy = math.cos(th), math.sin(th)
    pts = [(cx + s * half * ux, cy + s * half * uy) for s in (-1, 0, 1)]
    for angle in layout.angles().values():
        f = dock.frame(angle)
        M, d = f["M"], f["d"]
        for px, py in pts:
            t = max(-dock.BAR_LENGTH / 2, min(dock.BAR_LENGTH / 2, (px - M[0]) * d[0] + (py - M[1]) * d[1]))
            dist = math.hypot(px - (M[0] + t * d[0]), py - (M[1] + t * d[1]))
            if dist < dock.FLANGE_W / 2 + Y["slot_width"] / 2 + 1.0:
                return True
    return False


def _skin(r: float) -> cq.Workplane:
    t = Y["thickness"]
    s = cone(0).cut(cone(-t))
    k = (AO - AI) * t / DROP + 1
    s = s.cut(cq.Workplane("XY").workplane(offset=-DROP - 5).ellipse(AI + k, BI + k).extrude(DROP + 10))
    s = safe_cut(s, gap())   # the front gap goes in before the edges are rounded, so its edges are rounded too
    if r:
        s = s.edges().fillet(r)  # round the neck, gap and shoulder edges before the slots are cut
    # lightening slots, skipping the front gap (60-120 deg)
    length = (AO - AI) * 0.30
    for deg in range(15, 360, 30):  # offset so no slot sits on a seam
        if 60 <= deg <= 120 or slot_blocked(deg):
            continue
        th = math.radians(deg)
        a, b, _ = surface(0.42)
        slot = (
            cq.Workplane("XY").workplane(offset=-DROP - 5).center(math.cos(th) * a, math.sin(th) * b)
            .slot2D(length, Y["slot_width"], deg).extrude(DROP + 10)
        )
        s = s.cut(slot)
    s = cut_peg_holes(s)
    return s


def cut_peg_holes(s: cq.Workplane) -> cq.Workplane:
    """Holes through the skin for each docking bar's locating pegs (the bars' feet are glued over them)."""
    import dock
    import layout
    for angle in layout.angles().values():
        for px, py in dock.peg_points(angle):
            P, n, t_n = dock.peg_axis(float(px), float(py))
            hole = cq.Workplane("XY").add(cq.Solid.makeCylinder(dock.PEG_HOLE_R, t_n + 6, cq.Vector(*(P + 3 * n)), cq.Vector(*(-n))))
            try:
                trial = s.cut(hole)
                if trial.val().isValid() and 0 < trial.val().Volume() <= s.val().Volume():
                    s = trial
            except Exception:
                pass
    return s


def rail():
    """The docking bars of the active kit (one straight bar per module), in the formed collar; None when the kit has no modules."""
    import dock
    return dock.bars_formed()


BEAD_LEAD, BEAD_STEP = 0.2, 0.35                      # how far the two bead steps stand out from the hook's inner face
SNAP_DEPTH, SNAP_HEIGHT, SNAP_DROP = 0.5, 1.4, 1.2   # groove in the lip's outer face: depth, height, distance below the lip top
LIP_U = G["position"]


def lip_top() -> float:
    return surface(LIP_U)[2] + G["lip_height"]


def snap_z() -> float:
    """Bottom of the snap groove."""
    return lip_top() - SNAP_DROP - SNAP_HEIGHT


def lip(groove: bool = True) -> cq.Workplane:
    w = G["lip_width"]
    body = ring(LIP_U, w / 2, w / 2, surface(LIP_U)[2] - 1.5, G["lip_height"] + 1.5)  # sunk 1.5 mm into the skin, standing on it
    if not groove:
        return body
    return body.cut(ring(LIP_U, 0.5 - w / 2, w / 2 + 0.1, snap_z(), SNAP_HEIGHT))  # groove the diffuser bead snaps into


def led_channel() -> float:
    """Narrowest width between the rail cap and the lip, mm."""
    a0, b0, _ = surface(R["position"])
    al, bl, _ = surface(LIP_U)
    w = G["lip_width"]
    return min((al - w / 2) - (a0 + R["cap_width"] / 2), (bl - w / 2) - (b0 + R["cap_width"] / 2))


def strip_width() -> float:
    """Widest whole-mm LED strip (max 8) that lies in the channel beside the cover's inner wall with 0.5 mm to spare."""
    return min(8.0, math.floor(led_channel() - FLANGE_T - 1.5))


def pads() -> cq.Workplane:
    """Skin-conforming seats under the collar for stick-on foam."""
    t = Y["thickness"]
    band = cone(-t + 0.3).cut(cone(-t - PAD["thickness"]))  # overlaps the skin by 0.3 mm so the union is clean
    a, b, _ = surface(PAD["position"])
    out = None
    for deg in PAD["angles"]:
        th = math.radians(deg)
        cx, cy = a * math.cos(th), b * math.sin(th)
        heading = math.degrees(math.atan2(a * math.sin(th), b * math.cos(th)))
        oval = (
            cq.Workplane("XY").workplane(offset=-DROP - 8).center(cx, cy)
            .ellipse(PAD["length"] / 2, PAD["width"] / 2).extrude(DROP + 16)
            .rotate((cx, cy, 0), (cx, cy, 1), heading)
        )
        piece = band.intersect(oval)
        out = piece if out is None else out.union(piece)
    return out


FLANGE_T = 1.0          # the cover's inner wall, standing on the skin


def diffuser() -> cq.Workplane:
    """Thin cover over the LED channel between rail and lip. A hook on the outer edge grips the lip and a bead on it
    snaps into the groove in the lip."""
    t, gp = D["thickness"], D["gap"]
    a0, b0, zs = surface(R["position"])
    al, bl, _ = surface(LIP_U)
    w = G["lip_width"]
    rail_top = zs + R["stem_height"] + R["cap_height"]
    zA, zB = rail_top + gp + t, lip_top() + gp + t
    aA, bA = a0 + R["cap_width"] / 2, b0 + R["cap_width"] / 2
    aB, bB = al + w / 2 + 1.5, bl + w / 2 + 1.5
    dz = zA - zB
    ka, kb = (aB - aA) * t / dz + 0.5, (bB - bA) * t / dz + 0.5

    def frustum(shift):
        return (
            cq.Workplane("XY").workplane(offset=zA + shift).ellipse(aA - ka, bA - kb)
            .workplane(offset=-dz).ellipse(aB, bB).loft()
        )

    band = frustum(0).cut(frustum(-t))
    band = band.cut(cq.Workplane("XY").workplane(offset=zB - 5).ellipse(aA, bA).extrude(dz + 10))
    face = w / 2 + 0.1                                     # hook's inner face stands 0.1 mm off the lip's outer face
    hook = ring(LIP_U, -face, w / 2 + 1.5, zB - t - D["hook"], D["hook"] + t)
    # thin wall on the inner edge stands on the skin and holds the inner edge up (there is no continuous rail any more)
    flange = ring(R["position"], -(R["cap_width"] / 2), R["cap_width"] / 2 + FLANGE_T, zs - 10, zA - (zs - 10)).cut(cone(0.15))
    # two-step bead on the hook's inner face: narrow lead-in below, full step above with a flat top that locks in the groove
    z0 = snap_z() + 0.1
    lead = ring(LIP_U, -face + BEAD_LEAD, face + 0.2, z0, 0.5)
    step = ring(LIP_U, -face + BEAD_STEP, face + 0.2, z0 + 0.5, 0.7)
    return safe_cut(band.union(hook).union(flange).union(lead).union(step), gap())


def safe_cut(shape: cq.Workplane, tool: cq.Workplane) -> cq.Workplane:
    """shape.cut(tool), checked: the CAD kernel sometimes returns an inverted or empty solid for particular dimensions.
    Retries with the tool nudged by hundredths of a millimetre, and refuses to return a result it cannot trust."""
    base = shape.val().Volume()
    for dx, dy in ((0, 0), (0.013, 0.007), (-0.011, 0.017), (0.021, -0.013), (-0.017, -0.019)):
        try:
            out = shape.cut(tool.translate((dx, dy, 0))) if (dx or dy) else shape.cut(tool)
            v = out.val().Volume()
        except Exception:
            continue
        if out.val().isValid() and 0 < v <= base * 1.0001:
            return out
    raise RuntimeError("boolean cut failed on every retry: change a dimension by 0.1 mm and try again")


def drop_junk(shape: cq.Workplane) -> cq.Workplane:
    """Remove zero-volume or inverted solids that failed booleans can leave in a result: an inverted solid reads as
    'everything except itself' in the next intersection and shows up as a huge false overlap."""
    solids = [s for s in shape.solids().vals() if s.Volume() > 0.5]
    if len(solids) == len(shape.solids().vals()):
        return shape
    if len(solids) == 1:
        return cq.Workplane("XY").newObject([solids[0]])
    return cq.Workplane("XY").newObject([cq.Compound.makeCompound(solids)])


def safe_union(a: cq.Workplane, b: cq.Workplane) -> cq.Workplane:
    """a.union(b), checked like safe_cut: the result must be valid and no smaller than the larger input nor bigger than both added."""
    va, vb = a.val().Volume(), b.val().Volume()
    for dx, dy, dz in ((0, 0, 0), (0.013, 0.007, 0.004), (-0.011, 0.017, -0.006), (0.021, -0.013, 0.009), (-0.017, -0.019, -0.011)):
        try:
            out = a.union(b.translate((dx, dy, dz))) if (dx or dy or dz) else a.union(b)
            v = out.val().Volume()
        except Exception:
            continue
        if out.val().isValid() and max(va, vb) * 0.999 <= v <= (va + vb) * 1.0001:
            return drop_junk(out)
    # the kernel cannot fuse these two: keep them as one compound (they overlap by tenths of a millimetre at most), so
    # clash and balance checks still work; nothing is printed from this combined shape
    return cq.Workplane("XY").newObject([cq.Compound.makeCompound([a.val(), b.val()])])


def overlap_volume(a: cq.Workplane, b: cq.Workplane) -> float:
    """Volume a and b share, mm3. b is taken apart into its solids first: a compound whose solids overlap each other (the
    fallback when the kernel cannot fuse the formed collar) gives garbage in a single boolean."""
    total = 0.0
    for s in b.solids().vals():
        c = a.intersect(cq.Workplane("XY").newObject([s]))
        if c.vals():
            total += sum(x.Volume() for x in c.solids().vals())
    return total


def gap() -> cq.Workplane:
    return cq.Workplane("XY").box(Y["front_gap"], BO * 3, DROP * 4).translate((0, BO, 0))


def components() -> dict[str, cq.Workplane]:
    g = gap()
    comps = {"skin": skin(), "lip": safe_cut(lip(), g)}   # the foam seats are engraved marks in the flat skin, not printed pads
    bars = rail()
    if bars is not None:
        # the bars sit clear of the front gap; if a bar's flange does reach it, trim it, and if the kernel cannot, keep it whole
        if bars.intersect(g).vals() and sum(s.Volume() for s in bars.intersect(g).solids().vals()) > 1.0:
            try:
                comps["rail"] = safe_cut(bars, g)
            except RuntimeError:
                print("WARN a docking bar reaches the front gap and could not be trimmed: move that module's angle")
                comps["rail"] = bars
        else:
            comps["rail"] = bars
    return comps


def yoke(with_bars: bool = True) -> cq.Workplane:
    c = components()
    out = safe_union(c["skin"], c["lip"])
    return drop_junk(safe_union(out, c["rail"])) if with_bars and "rail" in c else out


QUARTERS = {"fr": (1, 1), "br": (1, -1), "fl": (-1, 1), "bl": (-1, -1)}


def quarter(shape: cq.Workplane, sx: int, sy: int) -> cq.Workplane:
    big = 2000
    return shape.intersect(cq.Workplane("XY").box(big, big, big).translate((sx * big / 2, sy * big / 2, 0)))


def dovetail(c: float = 0.0) -> list[tuple[float, float]]:
    """Dovetail tab in seam coordinates (s along the seam, d across it into the receiving part)."""
    r, t, n = J["root_width"] / 2 + c, J["tip_width"] / 2 + c, J["length"] + c
    return [(-r, -3), (r, -3), (r, 0), (t, n), (-t, n), (-r, 0)]


def prism(pts: list[tuple[float, float]]) -> cq.Workplane:
    return cq.Workplane("XY").workplane(offset=-DROP - 6).polyline(pts).close().extrude(DROP + 12)


def parts() -> dict[str, cq.Workplane]:
    """The four print quarters with dovetail joints: each seam has tabs on one part, matching pockets on the other."""
    y = yoke(with_bars=False)   # the seam joints are a property of the skin; the bars are separate parts
    q = {name: quarter(y, sx, sy) for name, (sx, sy) in QUARTERS.items()}

    # (tab part, pocket part, seam point at u -> (x, y), local (s, d) -> global (x, y))
    seams = [
        ("fr", "br", J["positions_side"], lambda u: (surface(u)[0], 0), lambda p, s, d: (p[0] + s, -d)),
        ("fl", "bl", J["positions_side"], lambda u: (-surface(u)[0], 0), lambda p, s, d: (p[0] - s, -d)),
        ("br", "bl", J["positions_back"], lambda u: (0, -surface(u)[1]), lambda p, s, d: (-d, p[1] + s)),
    ]
    for tab_part, pocket_part, us, at, to_global in seams:
        for u in us:
            p = at(u)
            tab = [to_global(p, s, d) for s, d in dovetail()]
            pocket = [to_global(p, s, d) for s, d in dovetail(J["clearance"])]
            q[tab_part] = q[tab_part].union(y.intersect(prism(tab)))
            q[pocket_part] = q[pocket_part].cut(prism(pocket))
    return q


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    y_full = yoke()
    y = yoke(with_bars=False)
    q = parts()
    bed, total = P["bed"], 0.0
    print(f"size {SIZE}")
    for line in sizing.fit_report(W, Y["front_gap"]):
        print("  " + line)
    print(f"{'part':<10}{'bbox mm (x y z)':<24}{'solid g':>9}  fits bed")
    for name, shp in q.items():
        bb = shp.val().BoundingBox()
        dims = sorted([bb.xlen, bb.ylen])
        fits = dims[0] <= min(bed[:2]) and dims[1] <= max(bed[:2]) and bb.zlen <= bed[2]
        g = shp.val().Volume() / 1000 * P["density_g_cc"]
        total += g
        cq.exporters.export(shp, str(OUT / f"yoke_{name}.stl"))
        print(f"yoke_{name:<5}{bb.xlen:6.0f}{bb.ylen:6.0f}{bb.zlen:6.0f}      {g:>9.1f}  {'ok' if fits else 'TOO BIG'}")
    cq.exporters.export(y_full, str(OUT / "assembly.stl"))

    # LED diffuser cover, printed in the same four quarters
    dif = diffuser()
    dtotal = 0.0
    for name, (sx, sy) in QUARTERS.items():
        dq = quarter(dif, sx, sy)
        bb = dq.val().BoundingBox()
        dims = sorted([bb.xlen, bb.ylen])
        fits = dims[0] <= min(bed[:2]) and dims[1] <= max(bed[:2]) and bb.zlen <= bed[2]
        g = dq.val().Volume() / 1000 * P["density_g_cc"]
        dtotal += g
        total += g
        cq.exporters.export(dq, str(OUT / f"diffuser_{name}.stl"))
        print(f"diffuser_{name:<2}{bb.xlen:6.0f}{bb.ylen:6.0f}{bb.zlen:6.0f}      {g:>9.1f}  {'ok' if fits else 'TOO BIG'}")
    touch = dif.intersect(y)
    touch_mm3 = sum(s.Volume() for s in touch.solids().vals()) if touch.vals() else 0.0

    # joint check: tabs must not overlap their pockets, and the glued-up whole should match the unsplit yoke
    overlap = sum(q[a].intersect(q[b]).val().Volume() for a, b in (("fr", "br"), ("fl", "bl"), ("br", "bl")) if q[a].intersect(q[b]).vals())
    missing = y.val().Volume() - sum(s.val().Volume() for s in q.values())
    loose = [n for n, s in q.items() if len(s.solids().vals()) != 1]
    need = max(max(s.val().BoundingBox().xlen, s.val().BoundingBox().ylen) for s in list(q.values()) + [quarter(dif, sx, sy) for sx, sy in QUARTERS.values()])
    print(f"\nSmallest square bed that fits every part: {need:.0f} mm (configured {min(bed[:2])} mm)")
    print(f"Total solid-plastic upper bound: {total:.0f} g (skin is thin, real print is lighter)")
    print(f"Joint check: part-to-part overlap {overlap:.2f} mm3 (want 0), clearance gaps {missing:.0f} mm3, loose crumbs: {loose or "none"}")
    channel = led_channel()
    print(f"LED channel between rail and lip: {channel:.1f} mm wide at its narrowest, takes a strip up to {strip_width():.0f} mm wide"
          + ("" if strip_width() >= 8 else "  <- narrower than the usual 8 mm strip: buy a thinner one"))
    plain = safe_cut(lip(groove=False), gap()).intersect(dif)
    flex = sum(s.Volume() for s in plain.solids().vals()) if plain.vals() else 0.0
    print(f"Snap: bead overlaps a groove-less lip by {flex:.0f} mm3 (that is what flexes on install); {SNAP_DEPTH} mm deep groove locks it; install by pressing down")
    print(f"Diffuser: {dtotal:.0f} g of the total, overlap with yoke {touch_mm3:.2f} mm3 (want 0)")
    print(f"Highest point: {y_full.val().BoundingBox().zmax:.1f} mm vs neck edge. Nothing above the shoulders.")


if __name__ == "__main__":
    main()
