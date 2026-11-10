# SPDX-FileCopyrightText: 2026 Jamey Warren
# SPDX-License-Identifier: MIT

"""
Driver — a representative 40 mm dynamic driver MOCKUP (reference only, NOT printed).

Shown in the assembly so the driver ↔ baffle ↔ clamp-ring fit reads. Real geometry is
the bought part; all dims here are REF / driver-pending (params.py).

Local frame: the mounting-flange FRONT face is at z=0 (the datum that seats on the
baffle's back recess ledge); the diaphragm fires +Z. The frame body + magnet run -Z.
"""

import cadquery as cq
from params import P


def make_driver() -> cq.Workplane:
    od = P.driver_od
    bd = P.driver_body_depth

    drv = cq.Workplane("XY").circle(od / 2).extrude(-bd)           # basket z(-bd..0); rim = the mount
    # diaphragm dome — a low cone firing +Z
    dome = cq.Solid.makeCone(
        P.driver_diaphragm_diameter / 2, P.driver_diaphragm_diameter * 0.15,
        P.driver_dome_proud, cq.Vector(0, 0, 0), cq.Vector(0, 0, 1))
    drv = drv.union(cq.Workplane(obj=dome))
    # rear magnet
    drv = drv.union(
        cq.Workplane("XY").workplane(offset=-bd)
        .circle(P.driver_magnet_diameter / 2).extrude(-P.driver_magnet_depth)
    )
    return drv


if __name__ == "__main__":
    cq.exporters.export(make_driver(), "output/driver.step")
    print("wrote output/driver.step")
