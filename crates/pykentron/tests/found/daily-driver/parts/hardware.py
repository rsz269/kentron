# SPDX-FileCopyrightText: 2026 Jamey Warren
# SPDX-License-Identifier: MIT

"""
Fastener hardware — M3 shoulder screw + heat-set insert (the yoke↔cup pivot).

WHY THIS EXISTS: the gate's pivot checks used to be pure parameter arithmetic
against a placeholder. This module emits the REAL fastener geometry so the fit is
validated against actual parts (does the shoulder span the eye, does the screw
reach the insert, is there wall around the installed insert, does the cup clear
the yoke through the full tilt range).

GEOMETRY SOURCE OF TRUTH = clean cadquery PRIMITIVES (`*_envelope`), always
available on the core deps. The gate and the STL build never need anything beyond
`requirements.txt`. If cq_warehouse IS installed (requirements-dev.txt, optional)
its accurate library geometry is used for the ASSEMBLY VISUALISATION only —
guarded, with the primitive as the fallback, and at the SAME envelope dimensions
so a clearance result never depends on which path ran.

VERIFIED against cq_warehouse 0.8.0 (git-only, dormant since 2023):
  - It has NO `ShoulderScrew` class — confirmed by introspection — so the shoulder
    screw is composed here from primitives (stacked cylinders: head + smooth
    shoulder + thread). Its dims are ESTIMATE / caliper-pending in params.py.
  - `HeatSetNut("M3-0.5-Standard", "McMaster-Carr")` builds and measures OD 4.70 /
    length 5.70 — those are the REF insert dims in params.py.
"""

import math

import cadquery as cq
from params import P


# ---- Primitive envelopes (core deps; the gate + build use these) -------------

def shoulder_screw_envelope() -> cq.Workplane:
    """M3 shoulder screw as stacked cylinders, axis +Z.

    Datum z=0 is the shoulder↔thread boundary (the face the eye/boss stack starts
    at). Thread runs DOWN (-Z) into the insert; the smooth shoulder (the pivot
    bearing the yoke eye rides on) and the head run UP (+Z). cq_warehouse has no
    shoulder-screw class, so this primitive IS the screw — see module docstring.
    """
    thread_l = P.shoulder_screw_thread_length
    shoulder_l = P.shoulder_screw_shoulder_length
    head_h = P.shoulder_screw_head_height
    screw = (
        cq.Workplane("XY")
        .circle(P.shoulder_screw_thread_diameter / 2).extrude(-thread_l)        # thread, -Z
        .faces(">Z").workplane()
        .circle(P.shoulder_screw_shoulder_diameter / 2).extrude(shoulder_l)     # shoulder, +Z
        .faces(">Z").workplane()
        .circle(P.shoulder_screw_head_diameter / 2).extrude(head_h)             # head, +Z
    )
    return screw


def heatset_insert_envelope() -> cq.Workplane:
    """M3 heat-set insert as a plain cylinder (OD x length), axis +Z, base at z=0.

    Envelope of the INSTALLED brass OD (4.70 mm), not the undersized thermal bore —
    so wall checks measure material around the seated insert (the conservative case).
    """
    return (
        cq.Workplane("XY")
        .circle(P.heatset_insert_diameter / 2)
        .extrude(P.heatset_insert_length)
    )


# ---- Accurate (optional, cq_warehouse) for assembly visualisation only -------

def make_heatset_insert() -> cq.Workplane:
    """Accurate HeatSetNut if cq_warehouse is present, else the primitive envelope.

    Same OD/length either way, so swapping paths never changes a clearance result.
    """
    try:
        from cq_warehouse.fastener import HeatSetNut
        nut = HeatSetNut(size="M3-0.5-Standard", fastener_type="McMaster-Carr", simple=True)
        return cq.Workplane(obj=cq.Solid(nut.wrapped))
    except Exception:  # noqa: BLE001 — not installed / version drift → primitive
        return heatset_insert_envelope()


def make_shoulder_screw() -> cq.Workplane:
    """The shoulder screw. Always the primitive composition (no library class)."""
    return shoulder_screw_envelope()


def make_thumbscrew() -> cq.Workplane:
    """8-32 knurled-knob thumbscrew MOCKUP — knurled hand-grip head + threaded shaft, axis +Z.

    Datum z=0 is the shaft TIP (the face that presses the post); the shaft runs +Z
    out through the slider boss to the head on top. Assembly-viz only — the slider's
    heat-set insert + clearance bore are the real interface (this just shows the lock).
    """
    shaft = (cq.Workplane("XY")
             .circle(P.slider_thumbscrew_diameter / 2).extrude(P.thumbscrew_shaft_length))
    head = (cq.Workplane("XY").workplane(offset=P.thumbscrew_shaft_length)
            .circle(P.thumbscrew_head_diameter / 2).extrude(P.thumbscrew_head_height))
    ts = shaft.union(head)
    # Knurl hint: a ring of shallow flutes down the head rim so it reads as a
    # hand-grip thumbscrew, not a plain screw. Cosmetic, best-effort (cuts are safe).
    try:
        hr = P.thumbscrew_head_diameter / 2
        for i in range(16):
            a = math.radians(i * 360.0 / 16)
            flute = cq.Solid.makeCylinder(
                0.5, P.thumbscrew_head_height + 1.0,
                cq.Vector(hr * math.cos(a), hr * math.sin(a), P.thumbscrew_shaft_length - 0.5),
                cq.Vector(0, 0, 1))
            ts = ts.cut(cq.Workplane(obj=flute))
    except Exception:  # noqa: BLE001 — knurl is cosmetic; never fail the viz
        pass
    return ts
