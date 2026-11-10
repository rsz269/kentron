# SPDX-FileCopyrightText: 2026 Jamey Warren
# SPDX-License-Identifier: MIT

"""
Cup shell — the heart of the Daily Driver (v0.3 engineering pass).

Build order: acoustic void → shell → clean concentric grille (decoupled from the
bosses) → 4 wall-blended baffle bosses (front-facing heat-set bore) → 2 external
yoke pivot bosses on the side walls (radial heat-set bore for the fork shoulder
screw). Every dimension is an ESTIMATE flagged in params.py — measured values
overwrite cleanly.

Reconciled to design-spec.md v0.3: the baffle bosses moved OFF the grille and
into the perimeter wall, so the grille's outer ring is now its own radius
(P.grille_outer_ring_radius), no longer pinned to the boss circle.
"""

import math

import cadquery as cq
from params import P


def make_cup() -> cq.Workplane:
    od = P.cup_outer_diameter
    total_h = P.cup_total_height

    # 1. Solid blank — the front (+Z) stays a CYLINDER (pad seat + void + pivot bosses); the rear
    #    cup_dome_height flows into a CONVEX DOMED back (DT880/Denon family) that bulges from the OD at
    #    the dome top inward to a FLAT grille face of cup_back_face_radius at z=0 (the grille / closed-back
    #    ports sit on that flat). Lofted stacked circles — a sin profile so the dome meets the cylinder
    #    tangent-vertical (smooth) and bulges convex below. Replaces the old cylinder + single roundover.
    dome_h = P.cup_dome_height
    r_back = P.cup_back_face_radius
    nseg = 16                                          # dense sampling → smooth dome under a RULED loft
    wires = []
    for i in range(nseg + 1):
        t = i / nseg
        r = r_back + (od / 2 - r_back) * math.sin(t * math.pi / 2)   # convex; vertical tangent at the top
        wires.append(cq.Wire.makeCircle(r, cq.Vector(0, 0, dome_h * t), cq.Vector(0, 0, 1)))
    wires.append(cq.Wire.makeCircle(od / 2, cq.Vector(0, 0, total_h), cq.Vector(0, 0, 1)))   # cylinder to the front
    # RULED (straight between sections) — a smooth/spline loft overshoots and bulges the dome way out.
    cup = cq.Workplane(obj=cq.Solid.makeLoft(wires, ruled=True))

    # 2. Hollow the acoustic void from the front (+Z), leaving side walls of
    #    wall_thickness and a thicker closed back (cup_back_thickness) — the grille
    #    substrate the chamfer lives in. Explicit cut (not shell) so the back band
    #    can be thicker than the side wall.
    void_r = P.cup_interior_diameter / 2
    void = (
        cq.Workplane("XY").workplane(offset=P.cup_back_thickness)
        .circle(void_r).extrude(total_h)          # up through the open front (and beyond)
    )
    cup = cup.cut(void)

    # 2b. Soften the back-FACE edge: a small round where the convex dome meets the flat grille face
    #     (the dome itself does the main rounding now). Done HERE while the bottom is still a clean disc —
    #     OCC on this build declines fillets once the grille/bosses/flange complicate it.
    cup = cup.edges("<Z").fillet(P.cup_back_round)

    # 3. Rear STRUCTURAL GRILLE (Stage 1b) — a rigid TRIANGULAR ×3 lattice (three
    #    opposing bar layers at 0/60/120°) that carries the protection + stiffness,
    #    with the LOGO rings + dot riding FLUSH on top (single colour, co-planar). This
    #    inverts the old logo-as-structure grille: the mesh is the structure, the logo
    #    is decoration. Built as a "keep" union of solid members; the COMPLEMENT
    #    (zone disc − keep) is cut from the closed back so the members stay and the
    #    gaps open to the driver. Cut BEFORE the bosses so nothing slices a boss.
    r_out = P.grille_outer_ring_radius
    hub_r = P.grille_hub_diameter / 2
    zone_r = r_out + P.grille_outer_ring_width / 2  # outer edge of the outer ring = zone

    z0 = -1.0
    cut_h = P.cup_back_thickness + 2.0          # pierce the full (thicker) back band

    def _disc(radius):
        return cq.Workplane("XY").workplane(offset=z0).circle(radius).extrude(cut_h)

    def _ring(rc, rw):
        return (cq.Workplane("XY").workplane(offset=z0)
                .circle(rc + rw / 2).circle(max(rc - rw / 2, 0.01)).extrude(cut_h))

    zone = _disc(zone_r)

    # Structural triangular lattice FIRST (the mesh behind the logo). Parallel bars of
    # grille_lattice_member_width at grille_lattice_pitch, centred on the hub, in three
    # opposing layers. Bars run past the zone; the complement-cut trims them to the disc.
    mw = P.grille_lattice_member_width
    pitch = P.grille_lattice_pitch
    nbar = int(math.ceil(zone_r / pitch)) + 1
    keep = None
    for a in P.grille_lattice_angles:
        for k in range(-nbar, nbar + 1):
            bar = (cq.Workplane("XY").workplane(offset=z0)
                   .transformed(rotate=(0, 0, a))
                   .center(0, k * pitch)
                   .rect(2 * zone_r + pitch, mw)
                   .extrude(cut_h))
            keep = bar if keep is None else keep.union(bar)

    # LOGO on top — two concentric rings + the centre dot, flush in the same plane.
    keep = keep.union(_ring(P.grille_inner_ring_radius, P.grille_inner_ring_width))
    keep = keep.union(_ring(r_out, P.grille_outer_ring_width))
    keep = keep.union(_disc(hub_r))             # centre dot (logo)

    # OPEN grille (default) or CLOSED back + pluggable tuning ports — ONE toggle, two
    # variants (the closed-back conversion). The grille solids above are built either way
    # (cheap, discarded when closed); only the back treatment differs here.
    if P.cup_open_back:
        cup = cup.cut(zone.cut(keep))                       # open the grille gaps
    else:
        pcr = P.cup_port_circle_diameter / 2                # CLOSED: solid back + tuning ports
        for i in range(P.cup_port_count):
            a = math.radians(i * 360.0 / P.cup_port_count)
            px, py = pcr * math.cos(a), pcr * math.sin(a)
            port = (cq.Workplane("XY").workplane(offset=z0)
                    .center(px, py).circle(P.cup_port_diameter / 2).extrude(cut_h))
            cup = cup.cut(port)

    # 3c. DAMPING retaining RING — a thin ring on the interior back floor that locates a
    #     felt / open-cell disc over the grille (light rear damping; the felt is a soft good,
    #     BOM). Sits inside the baffle-boss circle (r35) so it never fouls a boss; embedded
    #     0.5 mm into the floor so it fuses to the solid lattice members in open-back mode.
    dfr = P.damping_felt_diameter / 2
    ring = (cq.Workplane("XY").workplane(offset=P.cup_interior_floor_z - 0.5)
            .circle(dfr + P.damping_ring_wall).circle(dfr)
            .extrude(P.damping_ring_height + 0.5))
    cup = cup.union(ring)

    # 4. Baffle-mounting bosses — BUTTRESSED columns (the maker flagged the bare
    #    columns as fragile / snap-off-able). Each is a column at the bolt circle
    #    (r=baffle_screw_radius), from the joint-lap top (baffle_boss_floor_z,
    #    frame-only) up to the baffle underside, PLUS a wider base FLARE that merges
    #    the boss into the cup wall over a much wider arc than the bare column's thin
    #    lens. Built solid (column + flare) THEN bored from the front-facing top —
    #    round-before-cut, so the flare adds real support material and the heat-set
    #    bore stays clear (a base fillet here yields an invalid solid — the boss
    #    floors in open cavity with no host floor to round into).
    boss_points = [
        (
            P.baffle_screw_radius * math.cos(math.radians(45 + i * 360 / P.baffle_screw_count)),
            P.baffle_screw_radius * math.sin(math.radians(45 + i * 360 / P.baffle_screw_count)),
        )
        for i in range(P.baffle_screw_count)
    ]
    fz, bh = P.baffle_boss_floor_z, P.baffle_boss_height
    for bx, by in boss_points:
        col = (cq.Workplane("XY").workplane(offset=fz).center(bx, by)
               .circle(P.baffle_boss_diameter / 2).extrude(bh))
        flare = (cq.Workplane("XY").workplane(offset=fz).center(bx, by)
                 .circle(P.baffle_boss_flare_diameter / 2).extrude(P.baffle_boss_flare_height))
        cup = cup.union(col).union(flare)
    for bx, by in boss_points:  # bores last (round-before-cut)
        bore = (cq.Workplane("XY").workplane(offset=fz + bh).center(bx, by)
                .circle(P.m3_insert_hole_diameter / 2).extrude(-P.insert_boss_depth))
        cup = cup.cut(bore)

    # 5. Yoke pivot bosses — two external bosses at 0/180 on the cup side walls,
    #    at mid-height, each a radial cylinder spanning the wall to an outer seat so
    #    it fully houses an M3 heat-set insert bored from the outside. The fork's
    #    shoulder screw threads into it. With the pad-driven 90 mm cup the wall is
    #    6 mm, so the boss (span 9) is ~4 mm proud and its inner end now stops IN the
    #    wall — no lug into the cavity (the old thin-wall lug is gone).
    r_out_boss = P.pivot_boss_outer_radius              # 49 (yoke_pivot_centres/2)
    r_in_boss = r_out_boss - P.pivot_boss_through_span  # 40 — inside the 6 mm wall
    span = P.pivot_boss_through_span
    zc = P.pivot_boss_z
    for sign in (+1, -1):
        # boss body: cylinder axis ±X, dia pivot_boss_diameter
        body = (
            cq.Workplane("YZ")
            .workplane(offset=sign * r_in_boss)
            .center(0, zc)
            .circle(P.pivot_boss_diameter / 2)
            .extrude(sign * span)
        )
        cup = cup.union(body)
        # heat-set bore from the outer face inward (insert installed from outside)
        bore = (
            cq.Workplane("YZ")
            .workplane(offset=sign * r_out_boss)
            .center(0, zc)
            .circle(P.m3_insert_hole_diameter / 2)
            .extrude(-sign * P.insert_boss_depth)
        )
        cup = cup.cut(bore)

    # 6. Edge treatment: the back-outer comfort/print break is now the chamfer in
    #    step 2b (the form pass "set the outer profile"), so the old no-op outer-
    #    wall fillet — which only ever warned on the bare cylinder — is retired.

    # 7. Earpad RETAINING FLANGE (DT770-style) — a thin brim at the cup's front
    #    OUTER edge that extends the perimeter OUTWARD (radially), so the earpad's
    #    skirt wraps over it and hooks BEHIND it. It sticks OUT toward the perimeter,
    #    NOT up toward the head, so the baffle stays flush (not recessed). The brim
    #    sits at the front edge (top flush with the rim); below it the wall steps
    #    back in, giving the pad skirt an undercut to grip. The brim edges are
    #    ROUNDED (soft-form) so it eases the pad over and feels good in the hand —
    #    done on the clean disc before the union (OCC won't fillet it after).
    flange_ir = od / 2 - 1.0                          # overlap the wall → solid union
    flange_or = od / 2 + P.pad_lip_extension          # brim sticks OUT to here
    flange = (
        cq.Workplane("XY")
        .workplane(offset=total_h - P.pad_lip_thickness)
        .circle(flange_or).circle(flange_ir)
        .extrude(P.pad_lip_thickness)
    )
    try:
        flange = flange.edges().fillet(P.pad_lip_round)   # round the whole brim
    except Exception as e:  # noqa: BLE001 — report, don't mask; the flange still stands
        print(f"  [warn] cup: pad-flange roundover skipped ({e}).")
    cup = cup.union(flange)

    # 8. Cable exit — a hole through the −Y wall (the cup's BOTTOM when worn: T_cup
    #    maps cup −Y → global −Z) so the driver cable leaves the cup. At the pivot
    #    mid-height (depth), clear of the ±X pivot bosses. A clean through-cut, last.
    cable = cq.Solid.makeCylinder(
        P.cable_exit_diameter / 2, od / 2 + 2.0,
        cq.Vector(0, -(od / 2 + 2.0), P.pivot_boss_z), cq.Vector(0, 1, 0))
    cup = cup.cut(cq.Workplane(obj=cable))

    return cup


if __name__ == "__main__":
    cq.exporters.export(make_cup(), "output/cup.stl")
    print("wrote output/cup.stl")
