"""
Parametric SFR strongback — ZLP V3.

Built by revolving a closed half-section profile 360° around the Z axis,
then cutting the central bore and the small bolt/instrument holes.

All geometry dimensions are required — no defaults. This follows the same
convention as create_reactor_vessel, keeping reactor-specific values in
the assembly file rather than the component builder.

Profile is computed from six intuitive dimensions:

                        ← flange_radius →
    ┌─────────────────────────────┐  ─┐
    │         top flange          │   │ total_height
    │                              \  │
    │                               \ │  ← tapered edge
    │                          ──────┘ │  (flange_radius → skirt_outer_radius)
    │                          │    ─┐ │  ← taper_bottom_z
    │                          │     │ │
    └──────────────────────────┘     │ ┘
       ← skirt_inner_radius          │   skirt_height
       ←── skirt_outer_radius ───────┘

profile_pts can be supplied as an advanced override for non-standard shapes.

Single public function:
    create_strongback()  — returns a single cq.Workplane solid
"""

from __future__ import annotations
import math
import cadquery as cq

from profile_from_straight_connections import create_profile_from_straight_connections
from utils import revolve_profile


def _cut_vertical_cylinder(
    solid: cq.Workplane,
    radius: float,
    z_bottom: float,
    z_top: float,
    x: float = 0.0,
    y: float = 0.0,
) -> cq.Workplane:
    """Cut a vertical cylinder through a solid."""
    h = z_top - z_bottom
    cutter = (
        cq.Workplane("XY")
        .workplane(offset=z_bottom)
        .circle(radius)
        .extrude(h)
        .translate((x, y, 0))
    )
    return solid.cut(cutter)


def create_strongback(
    # ── Geometry (required) ───────────────────────────────────────────────
    total_height:       float,   # overall axial height [m]
    flange_radius:      float,   # outer radius of flat top flange [m]
    skirt_outer_radius: float,   # outer radius of lower skirt [m]
    skirt_inner_radius: float,   # inner radius of lower skirt [m]
    skirt_height:       float,   # axial height of lower skirt [m]
    taper_bottom_z:     float,   # z where tapered edge meets skirt outer surface [m]

    # ── Bore and holes (required) ─────────────────────────────────────────
    bore_radius:            float,   # central bore radius [m]
    small_hole_radius:      float,   # small hole radius [m]
    small_hole_count:       int,
    small_hole_placement_r: float,   # radial distance of hole centres from axis [m]

    # ── Global position ───────────────────────────────────────────────────
    z_bottom: float = 0.0,

    # ── Advanced override ─────────────────────────────────────────────────
    profile_pts: list[tuple[float, float]] | None = None,
) -> cq.Workplane:
    """
    Build the strongback by revolving a half-section profile 360° around Z,
    then cutting the central bore and the small bolt/instrument holes.

    By default the profile is computed from the six geometry dimensions.
    Supply profile_pts directly only if you need a non-standard cross-section
    (all dimension parameters are then ignored).

    Heights (profile_z_top, small_hole_z_bottom) are always derived from
    whichever profile is used, so no extra height parameters are needed.

    Args:
        total_height:           overall axial height [m]
        flange_radius:          outer radius of flat top flange [m]
        skirt_outer_radius:     outer radius of lower skirt [m]
        skirt_inner_radius:     inner radius of lower skirt [m]
        skirt_height:           axial height of lower skirt section [m]
        taper_bottom_z:         z where tapered edge meets skirt outer surface [m]
        bore_radius:            radius of central bore [m]
        small_hole_radius:      radius of small bolt/instrument holes [m]
        small_hole_count:       number of small holes
        small_hole_placement_r: radial distance of small hole centres from axis [m]
        z_bottom:               translate whole solid so base sits at this z [m]
        profile_pts:            advanced override — (r, z) half-section points.
                                If provided, all dimension parameters are ignored.

    Returns:
        cq.Workplane — single solid with holes cut
    """
    # ── Build or use profile points ───────────────────────────────────────
    if profile_pts is not None:
        pts = profile_pts
    else:
        pts = [
            (0.0,               total_height),
            (flange_radius,     total_height),
            (skirt_outer_radius, taper_bottom_z),
            (skirt_outer_radius, 0.0),
            (skirt_inner_radius, 0.0),
            (skirt_inner_radius, skirt_height),
            (0.0,               skirt_height),
        ]

    # ── Derive heights from the profile ───────────────────────────────────
    profile_z_top = max(z for _, z in pts)

    inner_z_values = [z for r, z in pts if r == 0.0 and z > 0.0]
    if not inner_z_values:
        raise ValueError(
            "Cannot derive small_hole_z_bottom from profile_pts: "
            "no points with r=0 and z>0 found."
        )
    small_hole_z_bottom = min(inner_z_values)

    # ── 1. Revolve profile ────────────────────────────────────────────────
    profile = create_profile_from_straight_connections(pts, plane="XZ", closed=True)
    solid = revolve_profile(profile, angle=360, axis="Z")

    # ── 2. Central bore (full height) ─────────────────────────────────────
    solid = _cut_vertical_cylinder(
        solid,
        radius   = bore_radius,
        z_bottom = 0.0,
        z_top    = profile_z_top,
    )

    # ── 3. Small holes (through upper flange only) ────────────────────────
    for i in range(small_hole_count):
        angle = 2 * math.pi * i / small_hole_count
        solid = _cut_vertical_cylinder(
            solid,
            radius   = small_hole_radius,
            z_bottom = small_hole_z_bottom,
            z_top    = profile_z_top,
            x        = small_hole_placement_r * math.cos(angle),
            y        = small_hole_placement_r * math.sin(angle),
        )

    # ── 4. Translate ──────────────────────────────────────────────────────
    if z_bottom != 0.0:
        solid = solid.translate((0, 0, z_bottom))

    return solid


if __name__ == "__main__":
    from ocp_vscode import show
    sb = create_strongback(
        total_height        = 1.242,
        flange_radius       = 2.684,
        skirt_outer_radius  = 3.030,
        skirt_inner_radius  = 2.243,
        skirt_height        = 0.436,
        taper_bottom_z      = 0.356,
        bore_radius             = 0.303,
        small_hole_radius       = 0.0755,
        small_hole_count        = 6,
        small_hole_placement_r  = 0.900,
    )
    show(sb)



"""
Parametric SFR strongback — ZLP V3.

Built by revolving a closed half-section profile 360° around the Z axis,
then cutting the central bore and the Ø151 bolt/instrument holes.

Profile points (r, z) in metres from ESFR-SMART Section A-A:
    (0,      1.242)
    (2.684,  1.242)
    (3.030,  0.356)
    (3.030,  0.000)
    (2.243,  0.000)
    (2.243,  0.436)
    (0,      0.436)

Holes (from top-view drawing):
    - Central bore  : Ø606  mm, through full height
    - Small holes   : Ø151  mm × 6, equally spaced on r ≈ 900 mm placement radius,
                      cut through the top section (z = 0.436 → 1.242)

Single public function:
    create_strongback()  — returns a single cq.Workplane solid
"""
"""
from __future__ import annotations
import math
import cadquery as cq

from profile_from_straight_connections import create_profile_from_straight_connections
from utils import revolve_profile


# Default profile points (r, z) in metres
_DEFAULT_PROFILE = [
    (0.000, 1.242),
    (2.684, 1.242),
    (3.030, 0.356),
    (3.030, 0.000),
    (2.243, 0.000),
    (2.243, 0.436),
    (0.000, 0.436),
]


def _cut_vertical_cylinder(
    solid: cq.Workplane,
    radius: float,
    z_bottom: float,
    z_top: float,
    x: float = 0.0,
    y: float = 0.0,
) -> cq.Workplane:
    # Cut a vertical cylinder through a solid.
    h = z_top - z_bottom
    cutter = (
        cq.Workplane("XY")
        .workplane(offset=z_bottom)
        .circle(radius)
        .extrude(h)
        .translate((x, y, 0))
    )
    return solid.cut(cutter)


def create_strongback(
    profile_pts: list[tuple[float, float]] | None = None,

    # Central bore
    bore_radius: float          = 0.303,   # Ø606 / 2  [m]

    # Small holes (Ø151, through top section only)
    small_hole_radius: float    = 0.0755,  # Ø151 / 2  [m]
    small_hole_count: int       = 6,
    small_hole_placement_r: float = 0.900, # placement radius [m]
    small_hole_z_bottom: float  = 0.436,   # bottom of top section
    small_hole_z_top: float     = 1.242,   # top face

    z_bottom: float             = 0.0,
) -> cq.Workplane:
    """
    # Build the strongback by revolving a closed half-section profile 360° around Z,
    # then cutting the central bore and the small bolt/instrument holes.

    # Args:
    #     profile_pts:           (r, z) points for the half-section. Defaults to ESFR-SMART.
    #     bore_radius:           radius of central bore [m] (Ø606 → 0.303)
    #     small_hole_radius:     radius of small holes [m] (Ø151 → 0.0755)
    #     small_hole_count:      number of small holes (default 6)
    #     small_hole_placement_r: radial distance of small hole centres from Z axis [m]
    #     small_hole_z_bottom:   z of bottom face of small holes [m]
    #     small_hole_z_top:      z of top face of small holes [m]
    #     z_bottom:              translate whole solid so base sits at this z [m]

    # Returns:
    #     cq.Workplane — single solid with holes cut
"""
    pts = profile_pts or _DEFAULT_PROFILE

    # ── 1. Revolve profile ───────────────────────────────────────────────
    profile = create_profile_from_straight_connections(pts, plane="XZ", closed=True)
    solid = revolve_profile(profile, angle=360, axis="Z")

    # ── 2. Central bore ──────────────────────────────────────────────────
    solid = _cut_vertical_cylinder(
        solid,
        radius   = bore_radius,
        z_bottom = 0.0,
        z_top    = 1.242,
    )

    # ── 3. Small holes ───────────────────────────────────────────────────
    for i in range(small_hole_count):
        angle = 2 * math.pi * i / small_hole_count
        hx = small_hole_placement_r * math.cos(angle)
        hy = small_hole_placement_r * math.sin(angle)
        solid = _cut_vertical_cylinder(
            solid,
            radius   = small_hole_radius,
            z_bottom = small_hole_z_bottom,
            z_top    = small_hole_z_top,
            x        = hx,
            y        = hy,
        )

    # ── 4. Translate ─────────────────────────────────────────────────────
    if z_bottom != 0.0:
        solid = solid.translate((0, 0, z_bottom))

    return solid


# ---------------------------------------------------------------------------
# Example
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    from ocp_vscode import show
    sb = create_strongback()
    show(sb)
"""