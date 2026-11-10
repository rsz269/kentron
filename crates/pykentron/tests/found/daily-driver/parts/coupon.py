# SPDX-FileCopyrightText: 2026 Jamey Warren
# SPDX-License-Identifier: MIT

"""
Fit coupons — small printable QA pieces that LOCK a toleranced interface against
real hardware before committing a full cup/baffle print.

Two independent coupons, each derived ENTIRELY from the real interface params (so
a coupon can never drift from the part it validates — change params.py and the
coupon regenerates with it):

  driver_coupon — reproduces the baffle's BACK driver interface: the shallow seat +
    locating collar (the measured driver registers into them) plus the 3 standoff
    bosses at the clamp bolt circle, so the REAL driver_clamp ring bolts on and you can
    confirm the seat/collar fit + rear-rim capture / standoff height. Frame matches
    baffle.py: BACK face z=0, FRONT z=baffle_thickness; recess on the back, bosses
    run −z. A central recess puck + 3 spokes to the bosses (trims the full Ø77
    baffle plate to just the bits under test).

  pad_coupon — a FULL grip ring at cup_outer_diameter (the dia the Tier-1 pad's
    mount-skirt grips) with the DT770-style retention lip on top. Full ring, not an
    arc: the question is "does the body grip ~91.4?" — a hoop-interference question
    that only a full-circumference ring answers honestly (an arc lets the elastic
    skirt splay and under-reads the grip). Open ring (no floor) — the pad slides
    over it from below and hooks behind the lip.

All FIT dims come from params; only coupon scaffolding (coupon_*) is coupon-local.
"""

import math

import cadquery as cq
from params import P


def make_driver_coupon() -> cq.Workplane:
    t = P.baffle_thickness                          # real plate thickness → recess depth + lamina are true
    rec_r = P.driver_recess_diameter / 2            # back recess the driver frame drops into
    ap_r = P.driver_aperture / 2                    # acoustic aperture (driver fires through)
    puck_r = rec_r + P.coupon_driver_puck_margin    # central puck OD = recess + 2× margin
    cbr = P.driver_clamp_bolt_circle / 2            # standoff bolt circle (the clamp ears land here)
    boss_r = P.insert_boss_diameter / 2
    boss_h = P.driver_clamp_standoff                # = driver_body_depth − driver_recess_depth (the protrusion)

    # 1. Central plate puck (clean), z 0..t — hosts the recess + aperture.
    coupon = cq.Workplane("XY").circle(puck_r).extrude(t)

    # 2. Three spokes out to the clamp bolt circle, full plate thickness so each
    #    standoff's M3 insert bore has backing material above it (insert_boss_depth
    #    reaches into the plate, not just the 2 mm boss).
    for i in range(P.driver_clamp_count):
        a = i * 360.0 / P.driver_clamp_count        # 0 / 120 / 240
        r0, r1 = puck_r - 1.0, cbr + boss_r + 1.0   # overlap the puck; reach past the boss
        spoke = (
            cq.Workplane("XY")
            .transformed(rotate=(0, 0, a))
            .center((r0 + r1) / 2, 0.0)
            .rect(r1 - r0, P.coupon_driver_tab_width)
            .extrude(t)
        )
        coupon = coupon.union(spoke)

    # 3. Standoff bosses on the BACK (−z), height = the driver's behind-baffle
    #    protrusion — identical to baffle.py so the captured stack-up is the real one.
    for i in range(P.driver_clamp_count):
        a = math.radians(i * 360.0 / P.driver_clamp_count)
        cx, cy = cbr * math.cos(a), cbr * math.sin(a)
        boss = (
            cq.Workplane("XY").workplane(offset=-boss_h).center(cx, cy)
            .circle(boss_r).extrude(boss_h)
        )
        coupon = coupon.union(boss)

    # 3b. Driver locating COLLAR on the BACK (z = -collar_height .. 0) — same short
    #     wall as the baffle, so the coupon checks the real seat + collar location.
    collar = (
        cq.Workplane("XY").workplane(offset=-P.driver_collar_height)
        .circle(P.driver_recess_diameter / 2 + P.driver_collar_wall)
        .circle(P.driver_recess_diameter / 2).extrude(P.driver_collar_height)
    )
    coupon = coupon.union(collar)

    # 4. Cuts LAST (this OCC build rounds/cuts in that order): aperture, back recess,
    #    then the M3 insert bores through each standoff. Mirrors baffle.py exactly.
    aperture = cq.Workplane("XY").workplane(offset=-0.5).circle(ap_r).extrude(t + 1.0)
    coupon = coupon.cut(aperture)
    recess = (
        cq.Workplane("XY").workplane(offset=-0.5)
        .circle(rec_r).extrude(P.driver_recess_depth + 0.5)
    )
    coupon = coupon.cut(recess)
    for i in range(P.driver_clamp_count):
        a = math.radians(i * 360.0 / P.driver_clamp_count)
        cx, cy = cbr * math.cos(a), cbr * math.sin(a)
        bore = (
            cq.Workplane("XY").workplane(offset=-boss_h - 0.5).center(cx, cy)
            .circle(P.m3_insert_hole_diameter / 2).extrude(P.insert_boss_depth + 0.5)
        )
        coupon = coupon.cut(bore)

    return coupon


def make_pad_coupon() -> cq.Workplane:
    R = P.cup_outer_diameter / 2          # the grip surface — the dia the pad skirt grips
    wall = P.coupon_pad_ring_wall
    H = P.coupon_pad_ring_height
    inner = R - wall
    lip_ext = P.pad_lip_extension         # DT770-style lip sticks OUT this far radially
    lip_th = P.pad_lip_thickness

    # Grip wall: open annular ring (no floor), z 0..H. The pad slides over the
    # R-radius outer face from below.
    ring = cq.Workplane("XY").circle(R).circle(inner).extrude(H)

    # Retention lip at the rim (top), protruding OUT to R+lip_ext — the pad skirt's
    # return hooks under its z=(H−lip_th) underside. Clean annular extrude, unioned.
    lip = (
        cq.Workplane("XY").workplane(offset=H - lip_th)
        .circle(R + lip_ext).circle(inner).extrude(lip_th)
    )
    ring = ring.union(lip)

    # Best-effort brim roundover (eases the pad on, matches the cup lip feel). This
    # OCC build declines fillets on already-cut/complex parts — warn, don't mask.
    try:
        ring = ring.edges(">Z").fillet(P.pad_lip_round)
    except Exception as e:  # noqa: BLE001
        print(f"  [warn] pad_coupon: lip roundover skipped ({e}).")

    return ring


if __name__ == "__main__":
    cq.exporters.export(make_driver_coupon(), "output/driver_coupon.stl")
    cq.exporters.export(make_pad_coupon(), "output/pad_coupon.stl")
    print("wrote output/driver_coupon.stl + output/pad_coupon.stl")
