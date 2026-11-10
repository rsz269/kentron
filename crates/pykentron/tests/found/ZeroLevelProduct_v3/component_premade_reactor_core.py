"""
Parametric SFR reactor core � ZLP V3.

Single public function:
    create_reactor_core()  � solid cylinder or regular n-sided prism

The cross-section is controlled by:
    radius      circumscribed radius (vertex-to-centre), same convention as
                cq.Sketch().regularPolygon()
    n_sides     number of polygon sides; None or 0 ? circular cylinder

Returns a single cq.Workplane solid (same pattern as create_reactor_vessel).
"""

from __future__ import annotations

import cadquery as cq


def create_reactor_core(
    radius: float,
    height: float,
    z_bottom: float = 0.0,
    n_sides: int | None = None,
) -> cq.Workplane:
    """
    Build a reactor core solid (cylinder or regular polygon prism).

    Args:
        radius:    circumscribed radius [m]
        height:    axial height [m]
        z_bottom:  z coordinate of the bottom face [m] (default 0)
        n_sides:   number of polygon sides (>= 3), or None / 0 for a cylinder

    Returns:
        cq.Workplane solid with bottom face at z_bottom
    """
    if radius <= 0:
        raise ValueError(f"radius must be > 0, got {radius}")
    if height <= 0:
        raise ValueError(f"height must be > 0, got {height}")

    if n_sides:
        if int(n_sides) < 3:
            raise ValueError(f"n_sides must be >= 3, got {n_sides}")
        sketch = cq.Sketch().regularPolygon(radius, int(n_sides))
        solid = cq.Workplane("XY").placeSketch(sketch).extrude(height)
    else:
        solid = cq.Workplane("XY").circle(radius).extrude(height)

    return solid.translate((0, 0, z_bottom))


# ---------------------------------------------------------------------------
# Examples
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    from ocp_vscode import show, set_defaults
    import time
    set_defaults(reset_camera=True)

    # 1 — circular cylinder (no n_sides)
    core_cylinder = create_reactor_core(
        radius   = 1.5,    # [m]
        height   = 2.0,    # [m]
        z_bottom = 0.0,
    )
    show(core_cylinder, names=["cylinder"])
    time.sleep(2)

    # 2 — hexagonal prism (6 sides, typical SFR lattice)
    core_hex = create_reactor_core(
        radius   = 1.5,
        height   = 2.0,
        z_bottom = 0.0,
        n_sides  = 6,
    )
    show(core_hex, names=["hexagonal"])
    time.sleep(2)

    # 3 — square prism, offset in z
    core_square = create_reactor_core(
        radius   = 1.0,
        height   = 3.0,
        z_bottom = -1.5,
        n_sides  = 4,
    )
    show(core_square, names=["square"])
