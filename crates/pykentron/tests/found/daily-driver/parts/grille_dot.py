# SPDX-FileCopyrightText: 2026 Jamey Warren
# SPDX-License-Identifier: MIT

"""
Grille accent dot — the warm-orange center cap of the rear grille.

The makerphones mark's dot is its only accent (the rings are never filled). On a
single-extruder FDM build you can't print one part two-tone, so the accent is a
SEPARATE small part printed in the warm orange and pressed into the grille hub —
modular, on-brand, serviceable. Modelled at the cup's grille face (z=0) standing
proud toward the outboard (-Z); placed + coloured in the assembly. ESTIMATE dims.
"""

import cadquery as cq
from params import P


def make_grille_dot() -> cq.Workplane:
    return (
        cq.Workplane("XY")
        .circle(P.grille_dot_diameter / 2)
        .extrude(-P.grille_dot_proud)        # proud of the grille face, outboard
    )


if __name__ == "__main__":
    cq.exporters.export(make_grille_dot(), "output/grille_dot.stl")
    print("wrote output/grille_dot.stl")
