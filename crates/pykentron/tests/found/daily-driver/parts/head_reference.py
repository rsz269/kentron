# SPDX-FileCopyrightText: 2026 Jamey Warren
# SPDX-License-Identifier: MIT

"""
Reference HEAD — a stylised, KU100-flavoured dummy head for worn-fit context. VIZ ONLY.

Not a printed part. It exists so the 3D viewer can show the headphone ON a head: an egg ovoid
with a hint of a FACE (nose + brow), a NECK stub, and — the point — a simplified PINNA (ear) on
each side at the ear position, so you can see the earcup landing AROUND the ear. Deliberately
abstract (translucent, not photoreal), like a Neumann KU100. Built from scaled-sphere ellipsoids
(`transformGeometry`, since revolve is dead on this OCC build) unioned together; scales uniformly
by ear-to-ear breadth so S/M/L share a shape. The assembly drops the ears at the cup level.

Frame: +X = right ear, +Y = FACE (front), +Z = up; centred at the origin (the assembly lifts it
by head_ref_z so the ears land at the cups).
"""

import cadquery as cq
from params import P


def _ovoid(half, center=(0.0, 0.0, 0.0)) -> cq.Workplane:
    """Ellipsoid as a scaled unit sphere. `half` = (ax, ay, az) semi-axes; `center` = (cx,cy,cz)."""
    ax, ay, az = half
    cx, cy, cz = center
    sph = cq.Solid.makeSphere(1.0, angleDegrees1=-90, angleDegrees2=90, angleDegrees3=360)
    m = cq.Matrix([[ax, 0, 0, cx], [0, ay, 0, cy], [0, 0, az, cz]])
    return cq.Workplane(obj=sph.transformGeometry(m))


def make_head_reference(ear_half: float = None) -> cq.Workplane:
    """A stylised reference head. `ear_half` = ear-to-ear HALF breadth; the MEDIUM head is scaled
    UNIFORMLY by ear_half / medium so S/M/L share a shape. Defaults to the medium head."""
    eh = P.head_ref_ear_half if ear_half is None else ear_half
    s = eh / P.head_ref_ear_half
    HX, HY, HZ = eh, P.head_ref_depth_half * s, P.head_ref_height_half * s
    ear_z = -P.head_ref_z          # local z that lands the ears at the cup level (global ≈ 0)

    # 1. Head: an egg ovoid (round cranium, narrower jaw via a slight downward taper cut later).
    head = _ovoid((HX, HY, HZ))

    # 2. NECK stub (KU100-style) — a SHORT cylinder dropping from under the jaw, ending just below the
    #    chin (a bust cut). Keeps the head off the viewer "ground": the contact shadow parks at the
    #    visible bottom, so a short neck puts the ground just under the chin instead of mid-face.
    neck = (cq.Workplane("XY").workplane(offset=-HZ - 10 * s)
            .circle(33 * s).extrude(26 * s).translate((0, -6 * s, 0)))
    head = head.union(neck)

    # 3. FACE hints on the +Y front (subtle person-likeness, not photoreal):
    #    a NOSE ridge on the centreline + a BROW above it.
    nose = _ovoid((9 * s, 18 * s, 24 * s), center=(0, HY * 0.92, ear_z + 14 * s))
    brow = _ovoid((34 * s, 12 * s, 8 * s), center=(0, HY * 0.84, ear_z + 40 * s))
    head = head.union(nose).union(brow)

    # 4. EARS — a simplified PINNA on each side: a vertical flattened ovoid proud of the head at the
    #    ear position. Sits slightly BACK (−Y) and at the cup level. Union first, then dish the CONCHA
    #    + a canal dimple (round/union BEFORE cut, per this OCC build).
    ear_y = -12 * s
    for sx in (+1, -1):
        pinna = _ovoid((10 * s, 17 * s, 30 * s), center=(sx * (HX - 1 * s), ear_y, ear_z))
        head = head.union(pinna)
    for sx in (+1, -1):
        concha = _ovoid((7 * s, 11 * s, 15 * s), center=(sx * (HX + 6 * s), ear_y, ear_z))
        head = head.cut(concha)
        canal = cq.Solid.makeCylinder(
            3.5 * s, 16 * s, cq.Vector(sx * (HX + 8 * s), ear_y, ear_z), cq.Vector(-sx, 0, 0))
        head = head.cut(cq.Workplane(obj=canal))

    return head


if __name__ == "__main__":
    cq.exporters.export(make_head_reference(), "output/head_reference.step")
    print("wrote output/head_reference.step")
