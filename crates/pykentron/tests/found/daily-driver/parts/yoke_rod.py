# SPDX-FileCopyrightText: 2026 Jamey Warren
# SPDX-License-Identifier: MIT

"""
Yoke adjustment post — a bought ISO 7379 SHOULDER SCREW. REFERENCE body (not printed).

The post is an off-the-shelf socket shoulder screw — Ø6 mm ground shoulder × M5 thread, 18-8
stainless, 50 mm shoulder — chosen so NOTHING is machined on the rod:
  • the ground Ø6 SHOULDER is the post — the printed/Delrin slider barrel SLIDES + SWIVELS on it
    (f9 ground ≈ Ø5.96–5.99 → ~0.4 mm clearance in the Ø6.4 bore; the plastic is the wear surface);
  • the M5 THREAD screws into an M5 heat-set in the printed fork boss (the easy step);
  • the HEAD (≈Ø10.2 × 4.5) is the built-in TOP STOP — wider than the bore, so a loosened slider
    can't slide off (the fork boss is the bottom stop).
This single fastener replaces the old tapped-rod + epoxy socket + separate stop-knob.

Built head-up from z=0 = the shoulder seat (= the fork boss top): thread below (−Z, into the
fork), the Ø6 shoulder up, the head on top.
"""

import cadquery as cq
from params import P


def make_yoke_rod() -> cq.Workplane:
    ds = P.yoke_post_diameter                 # Ø6 shoulder (the bearing)
    ls = P.yoke_post_length                   # 50 mm shoulder length (exposed post)
    dh = P.yoke_rod_head_diameter
    kh = P.yoke_rod_head_height
    lt = P.yoke_rod_thread_length
    dt = 5.0                                  # M5 thread

    screw = cq.Workplane("XY").circle(dt / 2).extrude(-lt)               # M5 thread, into the fork
    screw = screw.union(cq.Workplane("XY").circle(ds / 2).extrude(ls))   # Ø6 shoulder = the post
    screw = screw.union(cq.Workplane("XY").workplane(offset=ls)
                        .circle(dh / 2).extrude(kh))                     # head = top stop
    # Hex socket hint in the head top (it's a socket shoulder screw; cosmetic).
    try:
        screw = screw.faces(">Z").workplane().polygon(6, 3.0).cutBlind(-2.5)
    except Exception:  # noqa: BLE001 — cosmetic; never fail the reference body
        pass
    return screw


if __name__ == "__main__":
    cq.exporters.export(make_yoke_rod(), "output/yoke_rod.step")
    print("wrote output/yoke_rod.step")
