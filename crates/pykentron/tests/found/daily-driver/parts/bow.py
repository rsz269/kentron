# SPDX-FileCopyrightText: 2026 Jamey Warren
# SPDX-License-Identifier: MIT

"""
Head bow — REFERENCE BODY ONLY (v0.3 engineering pass).

This is NOT a printed part. The real bow is the BOUGHT Beyerdynamic Metal Head
Bow (or a DIY 1095 spring-steel blank to the same geometry). This module models a
reference body so the assembly reads correctly and the slider channel + end
interface have something to mate. It exports to STEP only (build.py excludes it
from the printed-parts STL set).

Every dimension is an ESTIMATE/REF (params.py) — TBD from the measured Beyer part.

Built from EXTRUDED annular sectors (this OCP build's `revolve` is unusable):
an arc band curving over the top (+Z) in the XZ plane, width along Y. Radius,
developed length, and width are MEASURED off the real Beyer band; the arc derives
from them (params helper). Two end-tab mounting holes per end model where the
slider mechanism fastens (their exact dia/pitch are TBD from the real part).
"""

import math

import cadquery as cq
from params import P


def _arc_band(ri, ro, deg_half, y_width):
    """Annular sector (ri..ro) centred on +Z (90°), ±deg_half, extruded ±Y."""
    a0, a1, am = 90 - deg_half, 90 + deg_half, 90

    def p(r, a):
        a = math.radians(a)
        return (r * math.cos(a), r * math.sin(a))

    Ao, Mo, Bo = p(ro, a0), p(ro, am), p(ro, a1)
    Bi, Mi, Ai = p(ri, a1), p(ri, am), p(ri, a0)
    return (
        cq.Workplane("XZ")
        .moveTo(*Ao).threePointArc(Mo, Bo)
        .lineTo(*Bi).threePointArc(Mi, Ai)
        .close()
        .extrude(y_width / 2, both=True)
    )


def _radial_cutter(profile, a_deg, R, depth=6.0):
    """A polygonal prism that pierces the thin band. `profile` is (s, y) points in a
    frame tangent to the arc at a_deg (s = along-arc, y = band width), extruded
    ±depth radially so it fully cuts the band despite the arc's curvature.
    """
    a = math.radians(a_deg)
    radial = (math.cos(a), 0.0, math.sin(a))
    tang = (-math.sin(a), 0.0, math.cos(a))
    origin = (R * math.cos(a), 0.0, R * math.sin(a))
    plane = cq.Plane(origin=origin, xDir=tang, normal=radial)
    return cq.Workplane(plane).polyline(profile).close().extrude(depth, both=True)


def make_bow(radius: float = None, arc_degrees: float = None) -> cq.Workplane:
    """Reference band. Defaults to the MEASURED at-rest geometry (params); the
    assembly passes the flexed (worn) radius/arc so the same strap can be posed
    opened-out on a head. Width + developed length are invariant under flex.
    """
    R = P.bow_radius if radius is None else radius
    th = P.bow_thickness
    arc = P.bow_arc_degrees if arc_degrees is None else arc_degrees
    half_arc = arc / 2

    # 1. Band: thin arc (thickness th, radial) × width (along Y), over the top.
    W = P.bow_width
    band = _arc_band(R - th / 2, R + th / 2, half_arc, W)

    # 2. Open the band like the real metal bow (maker's flat-lay photo): two outer RAILS
    #    run the FULL length and out to the ends as two separate PRONGS — the space
    #    BETWEEN them is HOLLOW all the way to the tips (a big open gap), one screw hole
    #    near each prong tip, and a single central X braces the middle. There is NO solid
    #    end tab. So: (a) cut the between-rails space fully open over the spans outside the
    #    central X (right out to the ends), then (b) cut the X's four triangular voids.
    if P.bow_pattern_enabled:
        a_start, a_end = 90 - half_arc, 90 + half_arc
        yin = W / 2 - P.bow_rail_width                           # rail inner edge (y)
        sw = P.bow_strut_width
        half_pat = math.degrees((P.bow_pattern_length / 2) / R)  # central X half-span (length-based)
        x_lo, x_hi = 90 - half_pat, 90 + half_pat

        # (a) hollow between the rails over the open spans — ALL THE WAY OUT to the ends
        #     (only the central X is left), so each end is two separate prongs with a big
        #     gap. Cut in short segments so the straight cutter follows the arc.
        for o_lo, o_hi in ((a_start, x_lo), (x_hi, a_end)):
            if o_hi - o_lo < 0.5:
                continue
            nseg = max(1, int(math.ceil((o_hi - o_lo) / 7.0)))
            dseg = (o_hi - o_lo) / nseg
            for i in range(nseg):
                a_c = o_lo + (i + 0.5) * dseg
                hs = math.radians(dseg) * R / 2 + 0.3            # half segment + overlap
                rect = [(-hs, yin), (hs, yin), (hs, -yin), (-hs, -yin)]
                band = band.cut(_radial_cutter(rect, a_c, R))

        # (b) central X bracing: n bays of four triangular voids leave crossing struts.
        n = P.bow_pattern_bays
        ang_bay = (x_hi - x_lo) / n
        for i in range(n):
            a_c = x_lo + (i + 0.5) * ang_bay
            hs = math.radians(ang_bay) * R / 2                   # half bay arc-length (mm)
            voids = [
                [(-hs + sw, yin), (hs - sw, yin), (0.0, sw)],        # top
                [(-hs + sw, -yin), (hs - sw, -yin), (0.0, -sw)],     # bottom
                [(-hs, yin - sw), (-hs, -yin + sw), (-sw, 0.0)],     # left
                [(hs, yin - sw), (hs, -yin + sw), (sw, 0.0)],        # right
            ]
            for v in voids:
                band = band.cut(_radial_cutter(v, a_c, R))

    # 3. PRONG-TIP mounting holes — one near each prong tip (at the two RAIL CENTRES,
    #    pitch bow_endtab_hole_spacing), set in from the tip, drilled radially through the
    #    band thickness. The band bolts to the slider clamp with two M3 screws per end; the
    #    rib registers in the big gap between the prongs. Layout from the maker's flat-lay.
    hole_r = P.bow_endtab_hole_diameter / 2
    s = P.bow_endtab_hole_spacing / 2                    # half-pitch across width (y)
    inset = math.degrees(P.bow_endtab_hole_inset / R)    # holes set in from the prong TIP
    for end_sign in (+1, -1):
        a = math.radians(90 + end_sign * (half_arc - inset))
        radial = cq.Vector(math.cos(a), 0, math.sin(a))
        base0 = cq.Vector(R * math.cos(a), 0, R * math.sin(a)) - radial * 2
        for y in (+s, -s):
            drill = cq.Solid.makeCylinder(hole_r, 4.0, base0 + cq.Vector(0, y, 0), radial)
            band = band.cut(cq.Workplane(obj=drill))

    # 4. PRONG-TIP ROUNDING — clip the two outer corners of each prong tip (a 45° chamfer per
    #    corner, via the radial cutter) so the rails end rounded, not square. The clip (tipr)
    #    is shorter than the hole inset, so the tip hole stays clear. Reference body, but it's
    #    what the manual/renders show. Only when the band is open into separate prongs.
    if P.bow_pattern_enabled:
        tipr = P.bow_prong_tip_r
        a_lo, a_hi = 90 - half_arc, 90 + half_arc
        yo = W / 2                              # prong OUTER edge
        yi = W / 2 - P.bow_rail_width           # prong INNER edge
        for a_deg, sgn_in in ((a_lo, +1.0), (a_hi, -1.0)):   # +s points toward the band centre
            for sgn_y in (+1.0, -1.0):                       # the two prongs (±Y rails)
                yO, yI = sgn_y * yo, sgn_y * yi
                band = band.cut(_radial_cutter(                # outer corner
                    [(0, yO), (sgn_in * tipr, yO), (0, yO - sgn_y * tipr)], a_deg, R))
                band = band.cut(_radial_cutter(                # inner corner
                    [(0, yI), (sgn_in * tipr, yI), (0, yI + sgn_y * tipr)], a_deg, R))

    return band


if __name__ == "__main__":
    cq.exporters.export(make_bow(), "output/bow.step")
    print("wrote output/bow.step  (REFERENCE BODY ONLY)")
