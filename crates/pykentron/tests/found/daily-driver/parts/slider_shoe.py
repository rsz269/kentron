# SPDX-FileCopyrightText: 2026 Jamey Warren
# SPDX-License-Identifier: MIT

"""
Slider pressure SHOE — the conformal pad the slider thumbscrew presses against the post.

The height lock must NOT drive a metal screw point onto the round Ø6 PRINTED PETG post: a
point load ploughs a divot and, over repeated adjustment, destroys the very surface that has
to keep sliding + swivelling cleanly. So the screw presses THIS small pad instead, and the
pad's concave face cradles the post over an AREA — the lock holds firm (HP1000-style) without
marring the post. This mirrors the aftermarket aluminium Grado rod-block's silicone 'slider
pill' (screw → pill → rod). Print in PETG (or cut from a Delrin blank); ×2 (one per ear).

Frame: the shoe is built centred at the origin — a small block with a concave cylindrical
SADDLE (axis Z, matching the vertical post) cut into its −Y face. The assembly poses it in
the slider's wall pocket with the saddle coaxial with the post. It drops into the pocket
through the bore during assembly; the post then traps it.
"""

import cadquery as cq
from params import P


def make_slider_shoe() -> cq.Workplane:
    sw = P.slider_shoe_width          # X (along the lozenge long axis)
    st = P.slider_shoe_thickness      # Y (radial: screw face → saddle)
    sh = P.slider_shoe_height         # Z (up the post)
    sr = P.slider_shoe_saddle_r
    sd = P.slider_shoe_saddle_depth

    # Block centred at the origin; the +Y face takes the screw tip, the −Y face the saddle.
    shoe = cq.Workplane("XY").box(sw, st, sh)

    # Concave SADDLE on the −Y face: cut a Z-axis cylinder positioned so it dishes the face by
    # `sd`. Radius sr (≈ post radius) → the cradle conforms to the Ø6 post over an area.
    cyl = cq.Solid.makeCylinder(
        sr, sh + 4, cq.Vector(0, -st / 2 - sr + sd, -(sh + 4) / 2), cq.Vector(0, 0, 1))
    shoe = shoe.cut(cq.Workplane(obj=cyl))
    return shoe


# Where the shoe sits in the slider frame: saddle coaxial with the post (Z axis), so the
# concave cradles the post's +Y surface. The assembly imports this to pose the shoe.
def shoe_offset_y() -> float:
    """+Y translation that puts the shoe's saddle axis on the post axis (y=0)."""
    return P.slider_shoe_thickness / 2 + P.slider_shoe_saddle_r - P.slider_shoe_saddle_depth


if __name__ == "__main__":
    cq.exporters.export(make_slider_shoe(), "output/slider_shoe.stl")
    print("wrote output/slider_shoe.stl")
