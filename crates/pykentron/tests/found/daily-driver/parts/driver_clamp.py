# SPDX-FileCopyrightText: 2026 Jamey Warren
# SPDX-License-Identifier: MIT

"""
Driver clamp ring — retains the driver from the BACK (3-bolt, the maker's prototype).

A 3-ear ring that slips over the back of the driver: the driver's rear nests SLIGHTLY
into a RECESS in the ring, the recess wall captures the driver's outer edge, and the
recess floor (an inner shoulder) bears on the back of the frame rim — pressing the
driver forward into the baffle's recess. The magnet protrudes back through the open
centre. The 3 ears bolt to standoff bosses on the baffle BACK. Independent of the
baffle→frame mount.

Form pass: the ear plate (floor + 3 mounting ears) is built as ONE 2D footprint so
the post→ring junctions and the ear pads come out BLENDED, not square — a sketch-
level (2D) fillet, because this OCC build won't fillet a 3D part once it has
cuts/unions (round BEFORE cut). The plate perimeter then gets a small roundover.

Local frame: the shoulder face (bears on the driver rear rim) is at z=0; the RECESS
opens +Z (the driver nests in from +Z); the floor + ears run -Z. All driver-fit dims
are REF / driver-pending (params.py).
"""

import math

import cadquery as cq
from params import P


def _safe_fillet(wp: cq.Workplane, selector: str, r: float, label: str) -> cq.Workplane:
    """Fillet selected edges, but only KEEP it if the result stays a valid solid.
    This OCC build will accept a fillet radius that silently invalidates the solid
    (no exception) — so we verify and fall back to the unfilleted shape. Warn, never
    mask: a skipped roundover is a cosmetic loss, an invalid solid is a broken part."""
    if r <= 0:
        return wp
    try:
        out = wp.edges(selector).fillet(r)
        if out.val().isValid():
            return out
        print(f"  [warn] driver_clamp: {label} r={r} made the solid invalid — skipped.")
    except Exception as e:  # noqa: BLE001
        print(f"  [warn] driver_clamp: {label} skipped ({e}).")
    return wp


def make_driver_clamp() -> cq.Workplane:
    ir = P.driver_clamp_inner_diameter / 2                       # open centre (magnet clears)
    rec_r = (P.driver_od + P.driver_clamp_recess_clearance) / 2  # recess: the driver nests in
    orr = rec_r + P.driver_clamp_wall                            # outer recess wall OD
    floor_th = P.driver_clamp_floor_thickness                    # shoulder floor
    rec_d = P.driver_clamp_recess_depth                          # how far the driver nests
    bcr = P.driver_clamp_bolt_circle / 2
    ear_r = P.driver_clamp_ear_diameter / 2
    n = P.driver_clamp_count
    fr = P.driver_clamp_fillet
    er = P.driver_clamp_edge_round

    angs = [math.radians(i * 360.0 / n) for i in range(n)]       # 0 / 120 / 240

    # --- Ear plate: shoulder floor/ring [ir..orr] + 3 mounting ears (pad + necked
    #     post to the ring), z(-floor_th .. 0). Built from clean 3D primitives, then
    #     the post→ring / post→pad corners are BLENDED in 3D on the still-uncut plate
    #     (OCC takes a fillet on these isolated junction edges; it refuses once holes
    #     are cut — so round BEFORE cut). The post is narrower than the pad so every
    #     junction is a real reflex corner that will take a fillet. ---
    bar_len = bcr - orr + 2 * ear_r          # ring → past the pad centre
    bar_w = P.driver_clamp_post_width        # < ear pad dia → reflex (filletable) corners
    clamp = (
        cq.Workplane("XY").workplane(offset=-floor_th)
        .circle(orr).circle(ir).extrude(floor_th)            # shoulder floor / ring
    )
    for a in angs:
        cx, cy = bcr * math.cos(a), bcr * math.sin(a)
        post = (
            cq.Workplane("XY").workplane(offset=-floor_th)
            .transformed(rotate=(0, 0, math.degrees(a)))
            .center((orr + bcr) / 2, 0).rect(bar_len, bar_w).extrude(floor_th)
        )
        pad = (
            cq.Workplane("XY").workplane(offset=-floor_th)
            .center(cx, cy).circle(ear_r).extrude(floor_th)
        )
        clamp = clamp.union(post).union(pad)

    # Blend the post↔ring / post↔pad corners (the vertical junction edges), then
    # roundover the perimeter rim (the 90° top/bottom corners). Both validated.
    clamp = _safe_fillet(clamp, "|Z", fr, "junction fillet")
    clamp = _safe_fillet(clamp, ">Z or <Z", er, "perimeter roundover")

    # Recess wall [rec_r..orr] rising +Z to rec_d — captures the driver's outer edge
    # as it nests into the recess (the floor [ir..rec_r] bears on the rear rim).
    wall = cq.Workplane("XY").circle(orr).circle(rec_r).extrude(rec_d)
    clamp = clamp.union(wall)

    # M3 clearance holes through each ear (axial) — cut LAST.
    for a in angs:
        ex, ey = bcr * math.cos(a), bcr * math.sin(a)
        hole = (
            cq.Workplane("XY").workplane(offset=-floor_th - 0.5).center(ex, ey)
            .circle(P.m3_clearance_hole / 2).extrude(floor_th + 1.0)
        )
        clamp = clamp.cut(hole)

    return clamp


if __name__ == "__main__":
    cq.exporters.export(make_driver_clamp(), "output/driver_clamp.stl")
    print("wrote output/driver_clamp.stl")
