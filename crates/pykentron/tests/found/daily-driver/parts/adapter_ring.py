# SPDX-FileCopyrightText: 2026 Jamey Warren
# SPDX-License-Identifier: MIT

"""
Step-down driver adapter ring — accessory, "design big, adapt down".

A printed ring that drops into a baffle built for a LARGER driver and presents a
seat + aperture for a SMALLER one, so one baffle can host several driver sizes
with no reprint — a real driver-testing workflow. This is a WORKED EXAMPLE
(50 mm-class host → the 40 mm reference driver); it does NOT change the reference
build, whose driver stays driver_od=42. Its params are the adapter_* set.

Frame: BACK at z=0, FRONT at z=adapter_height. The small driver drops into the
back pocket and seats on the front ledge; sound exits the front aperture. The ring
itself sits in the host baffle's driver recess.

ACOUSTIC HONESTY: a step-down ring lengthens/steps the front cavity — it is NOT
acoustically neutral. Ring variants are REW-loop items (see DESIGN-LOG), not a
free swap. All dims are ESTIMATE pending measured drivers.
"""

import cadquery as cq
from params import P


def make_adapter_ring() -> cq.Workplane:
    host_r = P.adapter_host_diameter / 2          # OD: fits the host baffle recess
    target_r = P.adapter_target_driver_od / 2     # back pocket for the smaller driver
    ap_r = P.adapter_target_aperture / 2          # front aperture (derives from target od)
    h = P.adapter_height
    seat = P.adapter_seat_thickness               # front seat floor the driver rests on

    ring = cq.Workplane("XY").circle(host_r).extrude(h)              # outer disc

    # Back pocket: the smaller driver frame drops in from behind and seats on the
    # ledge at z = h - seat (the front lamina the aperture pierces).
    pocket = (
        cq.Workplane("XY").workplane(offset=-0.5)
        .circle(target_r).extrude(h - seat + 0.5)
    )
    ring = ring.cut(pocket)

    # Front aperture through the seat lamina.
    aperture = (
        cq.Workplane("XY").workplane(offset=h - seat - 0.5)
        .circle(ap_r).extrude(seat + 1.0)
    )
    ring = ring.cut(aperture)

    return ring


if __name__ == "__main__":
    cq.exporters.export(make_adapter_ring(), "output/adapter_ring.stl")
    print("wrote output/adapter_ring.stl")
