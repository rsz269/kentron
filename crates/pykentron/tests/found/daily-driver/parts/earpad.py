# SPDX-FileCopyrightText: 2026 Jamey Warren
# SPDX-License-Identifier: MIT

"""
Earpad — a generalised round cushion MOCKUP (reference only, NOT printed).

The Daily Driver doesn't ship an earpad: builders fit their own round pad from the
broad Beyerdynamic range or aftermarket makers (Brainwavz, Dekoni, …). This module
models a representative round cushion so the assembly and the website views read
like a finished headphone, and so the cup's pad-grip OD / lip have something to seat.

Modelled as a torus cushion with a flat mounting base (a puffy ring with an ear
opening) — built from the makeTorus PRIMITIVE, not a revolve (revolve is unusable on
this OCC build). Local frame: the flat mounting base is at z=0 (seats on the cup
front rim); the cushion rises +Z toward the head, the ear opening down the centre.
All dimensions are a generic ESTIMATE (params.py) — measure/refine for a real pad.
"""

import cadquery as cq
from params import P


def make_earpad(depth: float = None) -> cq.Workplane:
    od = P.earpad_outer_diameter
    idd = P.earpad_inner_diameter
    rm = (od + idd) / 4.0          # torus mean radius (ring centre)
    tr = (od - idd) / 4.0          # tube radius → spans ID..OD (radial); height set below
    bf = P.earpad_base_flat
    d = P.earpad_depth if depth is None else depth   # relaxed (acoustic) depth, or a WORN/compressed one

    # Torus, axis +Z. The bare torus would be 2·tr tall; the pad DEPTH is `d` (earpad_depth for the
    # relaxed/acoustic height, or earpad_worn_depth for the compressed worn-fit view), so Z-scale the
    # solid to that height (transformGeometry — revolve/loft of a true profile is dead on this OCC
    # build). z ∈ [-tr·sz, tr·sz].
    sz = d / (2.0 * tr)
    torus = cq.Solid.makeTorus(rm, tr).transformGeometry(
        cq.Matrix([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, sz, 0]]))
    pad = cq.Workplane(obj=torus).translate((0, 0, tr * sz))   # base tangents z=0

    # Flatten the mounting base so it seats on the cup front rim, then drop the
    # flattened face to z=0 (cushion above, ear opening down the centre).
    pad = pad.cut(cq.Workplane("XY").workplane(offset=-1.0).circle(od).extrude(bf + 1.0))
    return pad.translate((0, 0, -bf))


if __name__ == "__main__":
    cq.exporters.export(make_earpad(), "output/earpad.step")
    print("wrote output/earpad.step  (REFERENCE MOCKUP — bring your own pad)")
