# SPDX-FileCopyrightText: 2026 Jamey Warren
# SPDX-License-Identifier: MIT

"""
Baffle plate — front-mount driver-mount plate with an integral guard (v0.4).

Local frame: BACK face at z=0, FRONT face at z=baffle_thickness. The driver mounts
cup-side (recess on the BACK) and fires forward through the aperture; the guard sits
recessed below the front face. In assembly the whole part is lifted to baffle_seat_z
so its front is flush with the cup rim. Pad retention is a lip on the CUP's outer rim
(DT770-style), NOT on the baffle — the baffle front is a clean plate.

v0.4 (maker pass): the plate is no longer a uniform slab. It is FULL thickness only in
the central driver/guard HUB (that depth is dome-gated — see baffle_thickness); the
outer RING is THINNER (front recessed) to shed bulk. The front venting opens from a few
tiny holes to large ARC-SLOTS between the 3 clamp-standoff SECTORS (the solid support),
backed by a glued acoustic paper/mesh that sits in a shallow FRONT depression — that
layer sets the back→front resistance (grade MEASUREMENT-GATED). All recesses open UPWARD
so the back-face-down print is unchanged (guard spokes still bridge fine).

All dimensions are ESTIMATES flagged in params.py.
"""

import math

import cadquery as cq
from params import P


def _arc_sector(ri, ro, a_center_deg, a_half_deg, z0, dz):
    """An annular SECTOR (radii ri..ro, centred on a_center_deg, ±a_half_deg) extruded
    dz from z0. Two radial edges + an outer and inner arc — a clean wedge cutter."""
    a0, a1, am = a_center_deg - a_half_deg, a_center_deg + a_half_deg, a_center_deg

    def p(rad, deg):
        a = math.radians(deg)
        return (rad * math.cos(a), rad * math.sin(a))

    return (cq.Workplane("XY").workplane(offset=z0)
            .moveTo(*p(ri, a0)).lineTo(*p(ro, a0))
            .threePointArc(p(ro, am), p(ro, a1))
            .lineTo(*p(ri, a1))
            .threePointArc(p(ri, am), p(ri, a0))
            .close().extrude(dz))


def make_baffle() -> cq.Workplane:
    r = P.baffle_outer_diameter / 2
    t = P.baffle_thickness                       # HUB (full) thickness — dome-gated
    ring_t = P.baffle_ring_thickness             # outer RING thickness (front recessed)
    hub_r = P.baffle_hub_radius                  # full-thickness out to here (driver + guard)
    ap_r = P.driver_aperture / 2

    # Aperture SHAPE hook — only round is authored today (see params + DESIGN-LOG).
    if P.driver_aperture_shape != "round":
        raise NotImplementedError(
            f"driver_aperture_shape={P.driver_aperture_shape!r}: only 'round' is "
            "built today. Author the non-round aperture/recess/guard before enabling."
        )

    # 1. STEPPED plate: a full disc at the RING thickness + a HUB standing up to full t.
    #    BACK face stays flat at z=0 (the driver-mount side + the print bed); the outer ring's
    #    FRONT is recessed by (t − ring_t) so the plate reads light instead of a thick slab.
    baffle = cq.Workplane("XY").circle(r).extrude(ring_t)
    baffle = baffle.union(cq.Workplane("XY").circle(hub_r).extrude(t))

    # 2. Acoustic aperture — through-hole (full t).
    aperture = cq.Workplane("XY").workplane(offset=-0.5).circle(ap_r).extrude(t + 1.0)
    baffle = baffle.cut(aperture)

    # 3. Driver SEAT on the BACK: a shallow pocket (driver_recess_diameter) cut
    #    driver_recess_depth up from the back face. The driver frame rim registers
    #    into it and seats on the ledge at z=driver_recess_depth; the narrower
    #    aperture carries on to the front. Shallow (1 mm) on purpose — the dome then
    #    stays low, clear of the front guard; the COLLAR below adds lateral location.
    recess = (
        cq.Workplane("XY")
        .workplane(offset=-0.5)
        .circle(P.driver_recess_diameter / 2)
        .extrude(P.driver_recess_depth + 0.5)
    )
    baffle = baffle.cut(recess)

    # 3b. Driver locating COLLAR — a short wall around the driver on the BACK,
    #    continuing the seat wall proud of the back face (z = -collar_height .. 0).
    collar = (
        cq.Workplane("XY")
        .workplane(offset=-P.driver_collar_height)
        .circle(P.driver_recess_diameter / 2 + P.driver_collar_wall)
        .circle(P.driver_recess_diameter / 2)
        .extrude(P.driver_collar_height)
    )
    baffle = baffle.union(collar)

    # 4. Integral driver GUARD across the aperture — concentric RINGS tied by radial
    #    SPOKES (a classic driver grille). It lives in the HUB's front lamina
    #    (z = recess_depth .. t). DOME CLEARANCE: the guard floor sits guard_dome_clearance
    #    ABOVE the dome's forward-most (excursed) position so the diaphragm never touches
    #    it. The lamina is thin, so the rib is auto-thinned to fit and the build WARNS the
    #    true clearances. (driver_dome_proud/excursion are still estimates — MEASURE.)
    lamina = t - P.driver_recess_depth
    dome_static = P.driver_recess_depth + P.driver_dome_proud       # cone peak at REST
    dome_dynamic = dome_static + P.driver_dome_excursion            # forward-most IN PLAY
    g_bot = dome_dynamic + P.guard_dome_clearance                  # guard floor
    g_th = max(0.8, min(P.guard_thickness, t - g_bot))            # fit under the front face
    g_top = g_bot + g_th
    pad_setback = t - g_top
    if pad_setback < P.guard_setback - 1e-6:
        how = "BLOWN — guard would sit past the front face" if pad_setback < -1e-6 else "tight"
        print(
            f"  [warn] baffle: dome clearance budget {how}. seat {P.driver_recess_depth}, dome "
            f"{P.driver_dome_proud}+excursion {P.driver_dome_excursion} → dynamic dome z{dome_dynamic:.1f}; "
            f"guard z{g_bot:.1f}–{g_top:.1f}, pad setback {pad_setback:.2f} (want {P.guard_setback}). "
            f"Both dome figures are estimates — MEASURE; shallower seat / thinner guard / deeper baffle if blown."
        )
    w = P.guard_member_width
    hub_grille_r = P.guard_hub_diameter / 2
    spoke_len = 2 * (ap_r + 1.0)
    guard = cq.Workplane("XY").workplane(offset=g_bot).circle(hub_grille_r).extrude(g_th)  # central hub
    for k in range(P.guard_ring_count):                         # concentric rings
        rk = hub_grille_r + (ap_r - hub_grille_r) * (k + 1) / (P.guard_ring_count + 1)
        ring = (
            cq.Workplane("XY").workplane(offset=g_bot)
            .circle(rk + w / 2).circle(rk - w / 2).extrude(g_th)
        )
        guard = guard.union(ring)
    for i in range(P.guard_spoke_count):                        # radial spokes
        ang = i * 360.0 / P.guard_spoke_count
        spoke = (
            cq.Workplane("XY").workplane(offset=g_bot)
            .transformed(rotate=(0, 0, ang)).rect(spoke_len, w).extrude(g_th)
        )
        guard = guard.union(spoke)
    baffle = baffle.union(guard)

    # 5. Four M3 clearance holes on the bolt circle (diagonals), counterbored from the
    #    FRONT. The screws live in the thinned RING, so the counterbore is referenced to
    #    the RING front (z=ring_t), heads sinking below it (hidden under the pad).
    bcr = P.baffle_screw_radius
    for i in range(P.baffle_screw_count):
        a = math.radians(45 + i * 360 / P.baffle_screw_count)
        cx, cy = bcr * math.cos(a), bcr * math.sin(a)
        through = (
            cq.Workplane("XY").workplane(offset=-0.5)
            .center(cx, cy).circle(P.m3_clearance_hole / 2).extrude(ring_t + 1.0)
        )
        cbore = (
            cq.Workplane("XY").workplane(offset=ring_t - P.baffle_counterbore_depth)
            .center(cx, cy).circle(P.baffle_counterbore_diameter / 2).extrude(P.baffle_counterbore_depth + 0.5)
        )
        baffle = baffle.cut(through).cut(cbore)

    # 6. Front venting = a SERIES OF HOLES in N "hot-dog" zones, each covered by a glued ARC STRIP of
    #    acoustic paper sitting in a shallow front DEPRESSION (maker: holes are easier to glue over
    #    consistently than open slots, and a few arc strips waste far less paper than one big annulus).
    #    The strip zones sit BETWEEN the 4 mounting screws (centred 0/90/180/270 for strip_count=4) so
    #    the screw bosses keep their strength; the holes auto-SKIP the 3 clamp standoffs (0/120/240) on
    #    the back. The paper (not the hole size) sets the back→front resistance — grade measurement-gated.
    #    All cut UPWARD-friendly for the back-down print (depression + holes open toward the front).
    n = P.baffle_vent_strip_count
    strip_half = P.baffle_vent_strip_half
    vin, vout = P.baffle_vent_inner_r, P.baffle_vent_outer_r
    r_mid = (vin + vout) / 2
    rec_d = P.baffle_paper_recess_depth
    hd = P.baffle_vent_hole_diameter
    # clamp-standoff footprints to dodge (back-side bosses at the clamp bolt circle)
    so_r = P.driver_clamp_bolt_circle / 2
    so_keepout = P.insert_boss_diameter / 2 + hd / 2 + 1.0
    standoffs = [(so_r * math.cos(math.radians(k * 360.0 / P.driver_clamp_count)),
                  so_r * math.sin(math.radians(k * 360.0 / P.driver_clamp_count)))
                 for k in range(P.driver_clamp_count)]
    n_skipped = 0
    for s in range(n):
        zc = s * 360.0 / n                                     # strip centre (0/90/180/270 for n=4)
        # 6a. paper-strip DEPRESSION — a "hot-dog" arc pocket (a touch wider/longer than the holes so
        #     the strip has a glue seat all round).
        dep = _arc_sector(vin - 0.8, vout + 0.8, zc, strip_half + 3.0, ring_t - rec_d, rec_d + 0.5)
        baffle = baffle.cut(dep)
        # 6b. the SERIES OF HOLES along the strip arc (count from arc length / pitch; centred). Skip
        #     any hole sitting over a clamp standoff.
        arc_len = 2 * math.radians(strip_half) * r_mid
        nh = max(1, int(arc_len // P.baffle_vent_hole_pitch) + 1)
        for j in range(nh):
            frac = 0.0 if nh == 1 else (j / (nh - 1) - 0.5)    # −0.5 .. 0.5 across the strip
            a = math.radians(zc + frac * 2 * strip_half)
            hx, hy = r_mid * math.cos(a), r_mid * math.sin(a)
            if any(math.hypot(hx - sx, hy - sy) < so_keepout for sx, sy in standoffs):
                n_skipped += 1
                continue                                       # dodge a clamp standoff
            hole = (cq.Workplane("XY").workplane(offset=-0.5)
                    .center(hx, hy).circle(hd / 2).extrude(ring_t + 1.0))
            baffle = baffle.cut(hole)
    if n_skipped:
        print(f"  [info] baffle: {n_skipped} vent hole(s) skipped where a strip crosses a clamp standoff.")

    # 7. Driver-clamp STANDOFFS — 3 bosses on the BACK face (z=0) at the clamp bolt circle, each with
    #    an M3 heat-set bore, at 0/120/240. The vent holes dodge them (above). Standoff ≈
    #    driver_body_depth − driver_recess_depth.
    cbr = P.driver_clamp_bolt_circle / 2
    boss_h = P.driver_clamp_standoff
    for i in range(P.driver_clamp_count):
        a = math.radians(i * 360.0 / P.driver_clamp_count)     # 0 / 120 / 240
        cx, cy = cbr * math.cos(a), cbr * math.sin(a)
        boss = (cq.Workplane("XY").workplane(offset=-boss_h).center(cx, cy)
                .circle(P.insert_boss_diameter / 2).extrude(boss_h))   # on the BACK (−z)
        baffle = baffle.union(boss)
        bore = (cq.Workplane("XY").workplane(offset=-boss_h - 0.5).center(cx, cy)
                .circle(P.m3_insert_hole_diameter / 2).extrude(P.insert_boss_depth + 0.5))
        baffle = baffle.cut(bore)

    return baffle


if __name__ == "__main__":
    cq.exporters.export(make_baffle(), "output/baffle.stl")
    print("wrote output/baffle.stl")
