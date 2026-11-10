# SPDX-FileCopyrightText: 2026 Jamey Warren
# SPDX-License-Identifier: MIT

"""
Fork-yoke — wraparound bracket that follows the earcup (v0.4 sleek-arm pass).

Two arms sweep from a top swivel hub down AROUND the round cup to two pivot eyes
on its ±X sides. Each arm follows a quarter-ellipse (semi-axes a = eye x, b = hub
z) so it mimics the cup's circular outline with clearance — wider over the top than
at the sides so the cup can tilt without striking the bracket (cf. the Extreme
Isolation yoke). Each eye takes an M3 shoulder screw into the cup's pivot boss (the
tilt joint, ±tilt_range); the swivel hub's vertical bore mates the slider above.

SLEEK-ARM FORM (v0.4): the arms are no longer flat tapered bars. Each arm is now a
single LOFTED tube — a rounded-rectangle cross-section (width tapers 9→6, thickness
6, with ~2 mm corner arcs) lofted through stations placed perpendicular to the
quarter-ellipse tangent. A plane cut across an arm shows ARCS, not a sharp
rectangle: no 90° edges run along the arm, so it reads as a soft organic form
rather than a blocky bar. The rounding is in the swept 2D section (the construction),
NOT a 3D edge fillet — `.fillet()` is unusable on this OCC build once the part has
cuts/unions (see docs/cadquery-build-notes.md). `revolve` is dead; `loft` works.

Local frame: pivot axis at z=0 (eyes at ±yoke_pivot_centres/2, 0, 0); the hub is at
z=yoke_fork_height. In assembly the frame is lifted so z=0 lands on the cup's pivot
bosses (cup mid-height).

All dimensions are ESTIMATES flagged in params.py.
"""

import math

import cadquery as cq
from cadquery import Plane, Vector
from params import P


def _ellipse_pt(sign, t):
    """Point on the quarter-ellipse at parameter t in [0,1] (the arm path).
    t=0 → eye at (sign*a, 0); t=1 → hub at (0, b). Returns (x, z)."""
    a = P.yoke_pivot_centres / 2
    b = P.yoke_fork_height
    ang = math.radians(90.0 * t)
    return (sign * a * math.cos(ang), b * math.sin(ang))


def _ellipse_tangent(sign, t):
    """Unit tangent (dx, dz) of the arm path at t, pointing eye→hub. Used to set
    each loft section's plane perpendicular to the path."""
    a = P.yoke_pivot_centres / 2
    b = P.yoke_fork_height
    ang = math.radians(90.0 * t)
    dx = -sign * a * math.sin(ang)
    dz = b * math.cos(ang)
    n = math.hypot(dx, dz)
    return (dx / n, dz / n)


def _rounded_rect_wire(plane, width, thick, r):
    """A rounded-rectangle WIRE on `plane`: local-x span = width, local-y span =
    thick, corner radius r. Built explicitly from 4 line segments + 4 tangent arcs
    (robust: no fillet helper, no findSolid — 2D fillets fail on this OCC build,
    but a hand-built arc wire lofts cleanly). Returns a cq Wire."""
    hw = width / 2.0
    ht = thick / 2.0
    r = min(r, hw - 1e-3, ht - 1e-3)
    sx = hw - r                              # half-length of the flats (x)
    sy = ht - r                              # half-length of the flats (y)
    k = r * 0.70710678                       # arc midpoint offset (45°)
    wp = cq.Workplane(plane).moveTo(-sx, -ht)
    wp = wp.lineTo(sx, -ht)                                            # bottom edge
    wp = wp.threePointArc((hw - r + k, -ht + r - k), (hw, -sy))        # BR corner
    wp = wp.lineTo(hw, sy)                                             # right edge
    wp = wp.threePointArc((hw - r + k, ht - r + k), (sx, ht))          # TR corner
    wp = wp.lineTo(-sx, ht)                                            # top edge
    wp = wp.threePointArc((-hw + r - k, ht - r + k), (-hw, sy))        # TL corner
    wp = wp.lineTo(-hw, -sy)                                           # left edge
    wp = wp.threePointArc((-hw + r - k, -ht + r - k), (-sx, -ht))      # BL corner
    wp = wp.close()
    return wp.val()


def _make_arm(sign, n_stations=12):
    """Loft one sleek arm through rounded-rect sections placed perpendicular to the
    quarter-ellipse tangent at each station. Width tapers yoke_arm_width →
    yoke_arm_hub_width; thickness = yoke_arm_thickness. The corner radius is capped
    by the local half-width so the slim hub end never self-intersects. Returns a
    Workplane carrying one solid (loft → a single continuous organic tube)."""
    arm_w = P.yoke_arm_width
    arm_w_hub = P.yoke_arm_hub_width
    arm_t = P.yoke_arm_thickness
    wires = []
    for i in range(n_stations + 1):
        t = i / n_stations
        cx, cz = _ellipse_pt(sign, t)
        tx, tz = _ellipse_tangent(sign, t)
        width = arm_w + (arm_w_hub - arm_w) * t
        r = min(P.yoke_arm_corner_radius, width / 2 - 1e-3, arm_t / 2 - 1e-3)
        # section plane ⟂ the tangent: normal along the tangent (in XZ), local-x
        # in the XZ plane (the WIDTH axis), local-y along global Y (the THICKNESS).
        pl = Plane(origin=Vector(cx, 0.0, cz),
                   xDir=Vector(-tz, 0.0, tx),
                   normal=Vector(tx, 0.0, tz))
        wires.append(_rounded_rect_wire(pl, width, arm_t, r))
    return cq.Workplane(obj=cq.Solid.makeLoft(wires, ruled=False))


def make_yoke() -> cq.Workplane:
    a = P.yoke_pivot_centres / 2             # 49 — eye x = ellipse semi-axis (sides)
    b = P.yoke_fork_height                   # 55 — hub z = ellipse semi-axis (top)
    arm_t = P.yoke_arm_thickness
    hub_z = P.yoke_fork_height

    yoke = None
    for sign in (+1, -1):
        x_eye = sign * a
        # eye: cylinder axis X around the pivot hole
        eye = (
            cq.Workplane("YZ")
            .workplane(offset=x_eye - arm_t / 2)
            .center(0, 0)
            .circle(P.yoke_pivot_eye_diameter / 2)
            .extrude(arm_t)
        )
        # The arm WRAPS the round cup: it follows a quarter-ellipse from the eye
        # (a, 0) at the cup's side up and over to the hub (0, b) at the top,
        # mimicking the cup's circular outline (cf. the Extreme Isolation bracket).
        # The ellipse is taller than wide (b > a), so the gap to the cup grows from
        # ~4 mm at the sides to ~(b − cup_r) at the top — the cup needs that extra
        # top room to tilt in/out without striking the bracket. The arm is a SLEEK
        # lofted tube (rounded section, no sharp edges along it) — see _make_arm.
        # The loft's first/last sections sit at the eye/hub centres, so each end
        # overlaps DEEP into the eye / hub cylinder and the union fuses to ONE solid
        # (a tangent kiss would leave two disjoint solids on this OCC build).
        arm = _make_arm(sign)
        piece = eye.union(arm)
        yoke = piece if yoke is None else yoke.union(piece)

    # Junction MOUNT-BOSS for the bought adjustment post — an ISO 7379 Ø6×M5 SHOULDER SCREW. The
    # arms tie into a short boss at the apex; an M5 HEAT-SET goes in the boss, the screw's M5 thread
    # screws into it, and its ground Ø6 SHOULDER rises as the post (the slider slides + swivels on
    # it — a smooth metal bearing, the plastic barrel is sacrificial; no rod machining). Separating
    # the post lets THIS fork print flat. The boss OD (Ø13) is wider than the slider barrel (Ø12), so
    # the shoulder seats on the boss top face and the boss is the BOTTOM END-STOP — a loosened slider
    # bottoms on it and can't ride off (the screw HEAD stops it the other way). M5 bore cut w/ pivots.
    bd = P.yoke_socket_boss_diameter
    boss_top = hub_z + 4                          # boss top = where the shoulder seats / exposed post starts
    boss_bot = boss_top - P.yoke_rod_mount_depth
    hub = (
        cq.Workplane("XY").workplane(offset=boss_bot)
        .circle(bd / 2).extrude(boss_top - boss_bot)
    )
    yoke = yoke.union(hub)

    # bores: pivot holes (axis X) through each eye. (The adjustment post is SOLID —
    # no bore; it slides in the slider and the slider thumbscrew locks it.) Cut LAST,
    # after every union, so the booleans stay on clean geometry.
    for sign in (+1, -1):
        x = sign * a
        bore = (
            cq.Workplane("YZ")
            .workplane(offset=x - (arm_t / 2 + 1))
            .center(0, 0)
            .circle(P.yoke_pivot_hole_diameter / 2)
            .extrude(arm_t + 2)
        )
        yoke = yoke.cut(bore)

    # M5 HEAT-SET bore in the boss for the shoulder screw's M5 thread (recessed below the boss top
    # so the screw's Ø6 shoulder seats on the full-diameter boss face, clear of the thread undercut).
    insert = (
        cq.Workplane("XY").workplane(offset=hub_z + 4 + 1)
        .circle(P.m5_insert_hole_diameter / 2)
        .extrude(-(P.yoke_rod_thread_length + 2))
    )
    yoke = yoke.cut(insert)

    # (Over-rotation stop slot REMOVED 2026-06-26 — the cup now rotates freely in the
    # yoke, Grado-style. The eye is a clean bored cylinder; no notch weakening it.)

    # NB no 3D edge fillet here: the arm's roundness is in the lofted 2D section, and
    # `.fillet()` fails on this OCC build once the part has cuts/unions (the eyes,
    # hub, bores, and stop slots above). See docs/cadquery-build-notes.md.

    # Tilt clearance: the wraparound arms follow an ellipse that clears the cup by
    # ~4 mm at the sides and ~10 mm over the top, so the cup tilts in/out without
    # striking the bracket. VERIFIED IN-CAD by gate.py (pivot-tilt-clearance):
    # rotating the cup through the full ±tilt_range adds <1% to the 0° cup∩yoke
    # bearing overlap. A test print should still confirm the real friction/feel.
    return yoke


if __name__ == "__main__":
    cq.exporters.export(make_yoke(), "output/yoke.stl")
    print("wrote output/yoke.stl")
