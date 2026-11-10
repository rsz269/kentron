# SPDX-FileCopyrightText: 2026 Jamey Warren
# SPDX-License-Identifier: MIT

"""
Vent plug — a press-fit plug for the CLOSED-BACK cup's tuning ports.

The closed-back variant (params.cup_open_back = False) replaces the open grille with a
ring of tuning PORTS. Press one of these plugs into a port to close it: plug N of M and
the rear openness becomes a MEASURABLE, REVERSIBLE tuning knob (no reprint, no tools).
Print a handful per build. Frame: shank up from z=0, a flange cap on top (won't push
through the port; doubles as a finger/pliers pull-grip). Sized off the same port params
as the cup, so the fit can't drift from the hole it plugs.
"""

import cadquery as cq
from params import P


def make_vent_plug() -> cq.Workplane:
    shank_r = P.cup_port_diameter / 2 - P.vent_plug_clearance   # press-fit into the port
    shank_h = P.cup_back_thickness                              # spans the back-band thickness
    flange_r = P.cup_port_diameter / 2 + P.vent_plug_flange     # head lip + pull grip
    flange_h = 1.5

    plug = cq.Workplane("XY").circle(shank_r).extrude(shank_h)
    cap = (cq.Workplane("XY").workplane(offset=shank_h)
           .circle(flange_r).extrude(flange_h))
    return plug.union(cap)


if __name__ == "__main__":
    cq.exporters.export(make_vent_plug(), "output/vent_plug.stl")
    print("wrote output/vent_plug.stl")
