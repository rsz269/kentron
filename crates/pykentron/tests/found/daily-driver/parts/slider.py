# SPDX-FileCopyrightText: 2026 Jamey Warren
# SPDX-License-Identifier: MIT

"""
Slider — OFFSET-OUTER clamp, flat rounded LOZENGE (v0.9, Beyerdynamic end-cap form).

Layout (maker's call): the metal band rides INSIDE (head-side); the rod + tube ride
OUTSIDE it. The clamp BODY is a flat rounded LOZENGE — a stadium outline (rounded ends),
beveled toward the barrel — like the Beyer headband end-cap, sleek and low-profile. The
post-bore TUBE sits on its OUTER wall (the rod for the yoke post), with enough material
around it for the side M3 THUMBSCREW that friction-locks the post.

    head ──  COVER │ BAND │ LOZENGE (+RIB) │ TUBE + POST  ── outside

The band's two prongs drop into a RECESS in the lozenge's inner face; a short registering
RIB enters the bow's open channel; the COVER plate (parts/headband_clamp.py) spans only
the grip region so the band SWEEPS OUT cleanly above it. Two M3 screws through the
prong-tip holes into inserts in the lozenge.

Frame: z=0 at the barrel MID; lozenge centred there; bolts at hole_z. ESTIMATES (params).
Rounded outlines are hand-built arc wires (this OCC build's 2D/3D fillets fail) and lofted.
"""

import math

import cadquery as cq
from params import P


def _rounded_rect_wire(plane, width, thick, r):
    """A rounded-rectangle WIRE on `plane` (local-x = width, local-y = thick, corner r).
    Four lines + four tangent arcs — robust where fillet helpers fail on this OCC build."""
    hw = width / 2.0
    ht = thick / 2.0
    r = min(r, hw - 1e-3, ht - 1e-3)
    sx, sy = hw - r, ht - r
    k = r * 0.70710678
    wp = (cq.Workplane(plane).moveTo(-sx, -ht)
          .lineTo(sx, -ht).threePointArc((hw - r + k, -ht + r - k), (hw, -sy))
          .lineTo(hw, sy).threePointArc((hw - r + k, ht - r + k), (sx, ht))
          .lineTo(-sx, ht).threePointArc((-hw + r - k, ht - r + k), (-hw, sy))
          .lineTo(-hw, -sy).threePointArc((-hw + r - k, -ht + r - k), (-sx, -ht))
          .close())
    return wp.val()


def _lozenge_wire(y, w, h, r):
    """Stadium outline in the X-Z plane at depth Y=y (centred on the barrel mid)."""
    plane = cq.Plane(origin=(0, y, 0), xDir=(1, 0, 0), normal=(0, 1, 0))
    return _rounded_rect_wire(plane, w, h, r)


def make_slider() -> cq.Workplane:
    R = P.slider_collar_diameter / 2
    h = P.slider_collar_height
    bore_r = (P.yoke_post_diameter + P.slider_post_clearance) / 2
    LX = P.slider_clamp_width
    LZ = P.slider_clamp_height
    rr = P.slider_clamp_corner_r
    bev = P.slider_clamp_bevel
    bevh = P.slider_clamp_bevel_head
    pd = P.slider_clamp_standoff
    bt = P.bow_thickness
    s = P.bow_endtab_hole_spacing / 2
    z_lo = P.slider_clamp_z_lo
    ovl = 4.0
    y_in = -R - pd                                       # lozenge inner (head-side) face
    y_out = -R + ovl                                     # lozenge outer face (buried in the barrel)

    # BARREL = post-bore TUBE, centred z=0 — the rod housing on the OUTER wall. The end rims
    # are CHAMFERED by construction (this OCC build's edge .fillet silently fails): loft a
    # smaller circle at each end face out to full R one rim-width in, giving a 45° eased rim
    # top and bottom. Mid-wall stays full R (the gate's collar-wall reads the nominal wall).
    rim = P.slider_collar_rim_round
    collar = cq.Workplane(obj=cq.Solid.makeLoft([
        cq.Wire.makeCircle(R - rim, cq.Vector(0, 0, -h / 2), cq.Vector(0, 0, 1)),
        cq.Wire.makeCircle(R, cq.Vector(0, 0, -h / 2 + rim), cq.Vector(0, 0, 1)),
        cq.Wire.makeCircle(R, cq.Vector(0, 0, h / 2 - rim), cq.Vector(0, 0, 1)),
        cq.Wire.makeCircle(R - rim, cq.Vector(0, 0, h / 2), cq.Vector(0, 0, 1))]))

    # CLAMP LOZENGE — flat rounded stadium, a shallow LENS in section: it eases to BOTH
    # faces from a widest mid-band, so neither face presents a proud square lip. Three arc
    # wires (lofts cleanly where fillets fail):
    #   • HEAD side (-Y, y_in): an INSET eased lip (bevh) so the perimeter RECEDES from the
    #     temple — the contact face is a smaller pillow, the widest section sits bevh behind
    #     it. This is the ergonomic head-relief (was a full-size square edge pointing at skin).
    #   • MID band (y_in+bevh): the full stadium — widest, set back from the head.
    #   • OUTER (+Y, y_out, buried in the barrel): the existing cosmetic bevel (the sleek edge).
    pill = cq.Solid.makeLoft([
        _lozenge_wire(y_in, LX - 2 * bevh, LZ - 2 * bevh, rr - bevh),
        _lozenge_wire(y_in + bevh, LX, LZ, rr),
        _lozenge_wire(y_out, LX - 2 * bev, LZ - 2 * bev, rr - bev)])
    collar = collar.union(cq.Workplane(obj=pill))

    # GUSSETS — fair the tube into the lozenge so it reads grown-in, not stuck-on. A small
    # triangular web on each side of the vertical tube↔lozenge junction; each web buries one
    # vertex deep in the lozenge AND one in the tube wall so the union stays one solid.
    xj = math.sqrt(max(R * R - y_out * y_out, 1.0))      # where the tube crosses the outer face
    g = P.slider_tube_gusset
    gz = P.slider_tube_gusset_z / 2
    for sx in (+1, -1):
        tri = [(sx * (xj + g), y_out - 1.5),             # outboard along the face, into the lozenge
               (sx * (xj - 2.0), y_out - 1.5),           # corner, buried in both bodies
               (sx * (R - 2.0), y_out + 3.0)]            # up into the tube wall (clear of the bore)
        web = (cq.Workplane("XY").workplane(offset=-gz)
               .polyline(tri).close().extrude(2 * gz))
        collar = collar.union(web)

    # RECESS — shallow pocket in the lozenge's inner face the band's prongs register into.
    rec_d = bt + 0.4
    rec_top = P.slider_clamp_hole_z + 5.0
    collar = collar.cut(
        cq.Workplane("XY").workplane(offset=(z_lo + rec_top) / 2)
        .center(0, y_in + rec_d / 2).box(P.bow_width + 1.0, rec_d + 0.02, rec_top - z_lo))
    seat_y = y_in + rec_d

    # Anti-rotation RIB into the bow's open channel; +Y end buried in the lozenge so it fuses.
    channel_w = P.bow_width - 2 * P.bow_rail_width
    rib_far = seat_y + 1.5
    rib_near = seat_y - P.slider_clamp_rib_depth
    collar = collar.union(
        cq.Workplane("XY")
        .workplane(offset=P.slider_clamp_rib_z - P.slider_clamp_rib_height / 2)
        .center(0, (rib_far + rib_near) / 2)
        .box(channel_w - 1.0, rib_far - rib_near, P.slider_clamp_rib_height,
             centered=(True, True, False)))

    # POST BORE through the barrel.
    collar = collar.cut(
        cq.Workplane("XY").workplane(offset=-h / 2 - 1).circle(bore_r).extrude(h + 2))

    # Eased post-bore MOUTHS — a shallow countersink cone at each end gives the post a lead-in
    # and a clean printed edge (constructed; the kernel won't fillet the bore rim).
    bc = P.slider_bore_chamfer
    for zc, dz in ((h / 2, -1.0), (-h / 2, 1.0)):
        cone = cq.Solid.makeCone(bore_r + bc, bore_r, bc,
                                 cq.Vector(0, 0, zc), cq.Vector(0, 0, dz))
        collar = collar.cut(cq.Workplane(obj=cone))

    # Two M3 inserts in the lozenge (axis +Y from the recess floor) at the prong-hole pitch.
    for x in (+s, -s):
        bore = cq.Solid.makeCylinder(
            P.m3_insert_hole_diameter / 2, P.insert_boss_depth,
            cq.Vector(x, seat_y, P.slider_clamp_hole_z), cq.Vector(0, 1, 0))
        collar = collar.cut(cq.Workplane(obj=bore))

    # THUMBSCREW boss + CAPTIVE-SHOE POCKET on the +Y OUTBOARD face, CENTRED on the barrel mid
    # (boss_z, on the x=0 centreline → keeps the slider L/R symmetric). The 8-32 knob threads a
    # heat-set in the boss and its tip presses a conformal SHOE (parts/slider_shoe.py) that
    # cradles the post — the metal never touches the printed post, so the bearing isn't gouged.
    # In the worn pose local +Y → global +X (straight out the side of the head) = the natural
    # two-finger reach with the phones ON.
    boss_h = P.slider_thumbscrew_boss_proud
    bovl = 4.0
    bz = P.slider_thumbscrew_boss_z
    boss = cq.Solid.makeCylinder(P.slider_thumbscrew_boss / 2, boss_h + bovl,
                                 cq.Vector(0, R - bovl, bz), cq.Vector(0, 1, 0))
    collar = collar.union(cq.Workplane(obj=boss))

    # SHOE POCKET — a rectangular slot in the barrel wall that hosts the pressure shoe. It opens
    # to the bore on −Y (the saddle reaches the post) and meets the insert on +Y (the screw tip
    # enters to push the shoe). The shoe drops in through the bore at assembly; the post traps it.
    pc = P.slider_shoe_clearance
    pkt_w = P.slider_shoe_width + 2 * pc
    pkt_h = P.slider_shoe_height + 2 * pc
    p_lo = P.yoke_post_diameter / 2 - 0.5            # ~3.5: opens into the bore so the saddle reaches the post
    p_hi = R                                          # barrel OD: the insert/screw sits just outboard of here
    pocket = cq.Solid.makeBox(pkt_w, p_hi - p_lo, pkt_h,
                              cq.Vector(-pkt_w / 2, p_lo, bz - pkt_h / 2))
    collar = collar.cut(cq.Workplane(obj=pocket))

    # 8-32 heat-set in the boss, ABOVE the pocket (so the thread engages solid material); the
    # screw tip protrudes past it into the pocket to press the shoe. Screw NEVER reaches the post.
    ins = cq.Solid.makeCylinder(P.slider_thumbscrew_insert_hole / 2, boss_h,
                                cq.Vector(0, R + boss_h, bz), cq.Vector(0, -1, 0))
    collar = collar.cut(cq.Workplane(obj=ins))
    # Countersink lead-in at the insert mouth (eases the heat-set start).
    sbc = P.slider_boss_chamfer
    icone = cq.Solid.makeCone(P.slider_thumbscrew_insert_hole / 2 + sbc,
                              P.slider_thumbscrew_insert_hole / 2, sbc,
                              cq.Vector(0, R + boss_h, bz), cq.Vector(0, -1, 0))
    collar = collar.cut(cq.Workplane(obj=icone))

    # FINGER SCALLOPS — a shallow concave channel down each ±X END of the lozenge gives the
    # hand a defined pinch to slide the block on the post (the height-adjust motion). A tall
    # vertical cylinder set just outboard of the end face, dipping `depth` in; centred in the
    # lozenge's Y depth, outboard of the inserts (x=±13) and recess (x=±17) so the gate is
    # untouched. Cut LAST, on the finished solid.
    sc_r = P.slider_grip_scallop_r
    sc_d = P.slider_grip_scallop_depth
    y_mid = (y_in + y_out) / 2
    for sx in (+1, -1):
        dish = cq.Solid.makeCylinder(
            sc_r, LZ + 8, cq.Vector(sx * (LX / 2 + sc_r - sc_d), y_mid, -(LZ + 8) / 2),
            cq.Vector(0, 0, 1))
        collar = collar.cut(cq.Workplane(obj=dish))

    # GUARANTEE L/R SYMMETRY → the SAME print serves both ears (no mirrored second part). The
    # part is symmetric BY CONSTRUCTION (every feature is centred or in a ±x pair, thumbscrew
    # on x=0) EXCEPT the two grip-scallop cuts, which this OCC kernel resolves a hair unevenly
    # against the lofted lozenge ends (~0.8 mm). So fold the finished part onto its own YZ
    # mirror: keep the +X half (half-space cut) and union its mirror. This drives the residual
    # to exactly 0. (intersect(part, mirror(part)) silently no-ops on this build; keeping a
    # half-space half then mirror-unioning is the reliable route.)
    half = collar.intersect(cq.Workplane("XY").transformed(offset=(200, 0, 0)).box(400, 400, 400))
    collar = half.union(half.mirror("YZ"))
    return collar


if __name__ == "__main__":
    cq.exporters.export(make_slider(), "output/slider.stl")
    print("wrote output/slider.stl")
