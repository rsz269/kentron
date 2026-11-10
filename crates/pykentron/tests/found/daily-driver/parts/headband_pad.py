# SPDX-FileCopyrightText: 2026 Jamey Warren
# SPDX-License-Identifier: MIT

"""
Headband pad — crown cushion that WRAPS the bow (rough-draft mockup).

A soft comfort cushion (foam, or printed TPU) that runs the full central arc of the
bow — BETWEEN the two end tabs (which stay bare for the slider screws) — and wraps
AROUND the steel band: a thick cushion on the head side, walls up both edges, and a
lip over the top, with the band nested in a channel through the middle. This reads as
a real wrapped headband pad rather than the old 80° crown sliver on the inner face.
Form, retention, and material are still TBD; every dimension is an ESTIMATE (params).

Built from extruded annular sectors (this OCC build's `revolve` is unusable). Apex at
the top (+Z) in the XZ plane, width along Y, so it takes the bow's assembly transform
directly. The band channel is cut as a narrower annular sector, open at the pad's two
arc ends where the band continues out to the tabs.
"""

import math

import cadquery as cq
from params import P


def _arc_band(rin, rout, half_deg, y_width):
    """Annular sector (rin..rout) centred on +Z (90°), ±half_deg, extruded ±Y."""
    b0, bm, b1 = 90 - half_deg, 90, 90 + half_deg

    def p(r, a):
        a = math.radians(a)
        return (r * math.cos(a), r * math.sin(a))

    Ao, Mo, Bo = p(rout, b0), p(rout, bm), p(rout, b1)
    Bi, Mi, Ai = p(rin, b1), p(rin, bm), p(rin, b0)
    return (
        cq.Workplane("XZ")
        .moveTo(*Ao).threePointArc(Mo, Bo)
        .lineTo(*Bi).threePointArc(Mi, Ai)
        .close()
        .extrude(y_width / 2, both=True)
    )


def make_headband_pad(radius: float = None, arc_degrees: float = None) -> cq.Workplane:
    """Crown cushion wrapping the bow. Defaults to the at-rest bow geometry; the
    assembly passes the flexed (worn) radius + arc so the pad hugs the worn band.
    The pad spans the FULL band arc minus the end tabs (bare for the slider screws).
    """
    R = P.bow_radius if radius is None else radius
    full_arc = P.bow_arc_degrees if arc_degrees is None else arc_degrees
    bt = P.bow_thickness
    clr = P.headband_pad_channel_clearance

    # Span the central arc BETWEEN the end tabs (leave the screw tabs bare).
    tab_ang = math.degrees(P.bow_endtab_length / R) + 4.0        # tab + small margin
    half = max(10.0, full_arc / 2 - tab_ang)

    # Cushion cross-section WRAPS the band: head-side cushion (thickness) below, a lip
    # over the top (wrap) above, walls up both edges (width > band width).
    ri = R - bt / 2 - P.headband_pad_thickness                  # head-side (inner) cushion face
    ro = R + bt / 2 + P.headband_pad_wrap                       # over-the-top (outer) cushion face
    pad = _arc_band(ri, ro, half, P.headband_pad_width)

    # Band channel: the bow nests in; cut slightly longer (half + 1°) so it's OPEN at
    # the pad's two arc ends, where the band continues out to the tabs.
    ci = R - bt / 2 - clr
    co = R + bt / 2 + clr
    channel = _arc_band(ci, co, half + 1.0, P.bow_width + 2 * clr)
    pad = pad.cut(channel)

    # Leather-cushion PLEATS: shallow rounded transverse grooves quilting the head-side
    # face (à la a Beyerdynamic pad — but NO snap buttons, kept generic). Each is a Y-axis
    # cylinder grazing the inner face at ri; spread evenly across the pad arc.
    gr, gd = P.headband_pad_pleat_radius, P.headband_pad_pleat_depth
    rc = ri - gr + gd                                            # centre just below the inner face
    n_pl = P.headband_pad_pleats
    for k in range(n_pl):
        a = math.radians(90 - half + (k + 0.5) * (2 * half) / n_pl)
        seam = cq.Solid.makeCylinder(
            gr, P.headband_pad_width + 4.0,
            cq.Vector(rc * math.cos(a), -(P.headband_pad_width / 2 + 2.0), rc * math.sin(a)),
            cq.Vector(0, 1, 0))
        pad = pad.cut(cq.Workplane(obj=seam))
    return pad


if __name__ == "__main__":
    cq.exporters.export(make_headband_pad(), "output/headband_pad.stl")
    print("wrote output/headband_pad.stl")
